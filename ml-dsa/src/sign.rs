//! FIPS 204 §6.2 / §5.2 — signing (Algorithms 7 and 2).
//!
//! `sign_internal` is the rejection-sampling loop (Fiat–Shamir with Aborts). The
//! **number** of iterations is secret-dependent and leaks — this is the documented,
//! accepted leakage analysed in Part 2; the loop *body* must not leak beyond that.
//! The accept/reject decisions use norm comparisons that Part 2 will make constant
//! time; the faithful baseline keeps them straightforward.

use rand_core::{CryptoRng, RngCore};

use crate::error::{Error, Result};
use crate::expand::{expand_a, expand_mask};
use crate::hash::H;
use crate::ntt::ntt;
use crate::ntt_arith::{matrix_vector_ntt, scalar_vector_ntt};
use crate::params::ParameterSet;
use crate::sample::sample_in_ball;
use crate::serdes::{sig_encode, sk_decode};
use crate::vecops::{
    add_vec, center_vec, count_ones, high_bits_vec, inf_norm, inv_ntt_vec, low_bits_vec,
    make_hint_vec, neg_vec, ntt_vec, sub_vec,
};

/// FIPS 204, Algorithm 7 — ML-DSA.Sign_internal: signature for a formatted message `M'`.
pub fn sign_internal<P: ParameterSet, const K: usize, const L: usize>(
    sk: &[u8],
    m_prime: &[u8],
    rnd: &[u8; 32],
) -> Vec<u8> {
    sign_internal_traced::<P, K, L>(sk, m_prime, rnd).0
}

/// Like [`sign_internal`], but also returns the number of rejection-loop iterations
/// (attempts until acceptance). That count is the documented, accepted leakage of
/// Fiat–Shamir with Aborts — the signal the Part 2 side-channel chapter studies.
pub fn sign_internal_traced<P: ParameterSet, const K: usize, const L: usize>(
    sk: &[u8],
    m_prime: &[u8],
    rnd: &[u8; 32],
) -> (Vec<u8>, u32) {
    let (rho, k_seed, tr, s1, s2, t0) =
        sk_decode::<P, K, L>(sk).expect("valid secret key");
    let s1_hat = ntt_vec(&s1); // line 2
    let s2_hat = ntt_vec(&s2); // line 3
    let t0_hat = ntt_vec(&t0); // line 4
    let a_hat = expand_a::<K, L>(&rho); // line 5

    // line 6: μ ← H(BytesToBits(tr) || M', 64)  == H(tr || m_prime, 64)
    let mut hm = H::init();
    hm.absorb(&tr);
    hm.absorb(m_prime);
    let mut mu = [0u8; 64];
    hm.finalize().squeeze(&mut mu);

    // line 7: ρ'' ← H(K || rnd || μ, 64)
    let mut hr = H::init();
    hr.absorb(&k_seed);
    hr.absorb(rnd);
    hr.absorb(&mu);
    let mut rho_pp = [0u8; 64];
    hr.finalize().squeeze(&mut rho_pp);

    let mut kappa = 0u32; // line 8
    let mut attempts = 0u32; // rejection-loop iteration counter (not part of the spec)
    loop {
        attempts += 1;
        // line 11-13
        let y = expand_mask::<P, L>(&rho_pp, kappa);
        let w = inv_ntt_vec(&matrix_vector_ntt(&a_hat, &ntt_vec(&y)));
        let w1 = high_bits_vec::<P, K>(&w);

        // line 15-17: commitment hash, challenge
        let mut hc = H::init();
        hc.absorb(&mu);
        hc.absorb(&crate::serdes::w1_encode::<P, K>(&w1));
        let c_tilde = hc.finalize().squeeze_vec(P::C_TILDE_BYTES);
        let c = sample_in_ball::<P>(&c_tilde);
        let c_hat = ntt(&c);

        // line 18-21
        let cs1 = inv_ntt_vec(&scalar_vector_ntt(&c_hat, &s1_hat));
        let cs2 = inv_ntt_vec(&scalar_vector_ntt(&c_hat, &s2_hat));
        let z = add_vec(&y, &cs1);
        let r0 = low_bits_vec::<P, K>(&sub_vec(&w, &cs2));

        // line 23: first validity check
        if inf_norm(&z) >= P::GAMMA1 - P::BETA || inf_norm(&r0) >= P::GAMMA2 - P::BETA {
            kappa += L as u32;
            continue;
        }

        // line 25-26: hint
        let ct0 = inv_ntt_vec(&scalar_vector_ntt(&c_hat, &t0_hat));
        let r_arg = add_vec(&sub_vec(&w, &cs2), &ct0); // w − cs2 + ct0
        let h = make_hint_vec::<P, K>(&neg_vec(&ct0), &r_arg);

        // line 28: second validity check
        if inf_norm(&ct0) >= P::GAMMA2 || count_ones(&h) > P::OMEGA {
            kappa += L as u32;
            continue;
        }

        // line 33: σ ← sigEncode(c~, z mod± q, h)
        return (sig_encode::<P, K, L>(&c_tilde, &center_vec(&z), &h), attempts);
    }
}

/// Format `M'` for the pure (non-prehash) variant:
/// `IntegerToBytes(0,1) || IntegerToBytes(|ctx|,1) || ctx || M`.
pub(crate) fn format_m_prime(ctx: &[u8], m: &[u8]) -> Vec<u8> {
    let mut mp = Vec::with_capacity(2 + ctx.len() + m.len());
    mp.push(0);
    mp.push(ctx.len() as u8);
    mp.extend_from_slice(ctx);
    mp.extend_from_slice(m);
    mp
}

/// FIPS 204, Algorithm 2 — ML-DSA.Sign (hedged): draws `rnd` from the injected RNG.
pub fn sign<P: ParameterSet, const K: usize, const L: usize, R: CryptoRng + RngCore>(
    sk: &[u8],
    m: &[u8],
    ctx: &[u8],
    rng: &mut R,
) -> Result<Vec<u8>> {
    if ctx.len() > 255 {
        return Err(Error::ContextTooLong);
    }
    let mut rnd = [0u8; 32];
    rng.fill_bytes(&mut rnd);
    Ok(sign_internal::<P, K, L>(sk, &format_m_prime(ctx, m), &rnd))
}

/// Deterministic variant (FIPS 204 §3.4): `rnd = {0}^32`. Not recommended where
/// side-channel attacks are a concern.
pub fn sign_deterministic<P: ParameterSet, const K: usize, const L: usize>(
    sk: &[u8],
    m: &[u8],
    ctx: &[u8],
) -> Result<Vec<u8>> {
    Ok(sign_deterministic_traced::<P, K, L>(sk, m, ctx)?.0)
}

/// Deterministic signing that also returns the rejection-loop iteration count.
pub fn sign_deterministic_traced<P: ParameterSet, const K: usize, const L: usize>(
    sk: &[u8],
    m: &[u8],
    ctx: &[u8],
) -> Result<(Vec<u8>, u32)> {
    if ctx.len() > 255 {
        return Err(Error::ContextTooLong);
    }
    Ok(sign_internal_traced::<P, K, L>(sk, &format_m_prime(ctx, m), &[0u8; 32]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::keygen::key_gen_internal;
    use crate::params::MlDsa65;

    const K: usize = MlDsa65::K;
    const L: usize = MlDsa65::L;

    #[test]
    fn sign_shapes_and_deterministic() {
        let (_pk, sk) = key_gen_internal::<MlDsa65, K, L>(&[0x11u8; 32]);
        let sig = sign_deterministic::<MlDsa65, K, L>(&sk, b"hello", b"").unwrap();
        assert_eq!(sig.len(), MlDsa65::SIG_BYTES);
        // deterministic variant is reproducible
        assert_eq!(sig, sign_deterministic::<MlDsa65, K, L>(&sk, b"hello", b"").unwrap());
    }

    #[test]
    fn ctx_too_long_is_rejected() {
        let (_pk, sk) = key_gen_internal::<MlDsa65, K, L>(&[0x22u8; 32]);
        let long = vec![0u8; 256];
        assert!(matches!(
            sign_deterministic::<MlDsa65, K, L>(&sk, b"m", &long),
            Err(Error::ContextTooLong)
        ));
    }
}

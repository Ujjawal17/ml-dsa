#![allow(clippy::needless_range_loop)]

use crate::expand::expand_a;
use crate::hash::H;
use crate::ntt::ntt;
use crate::ntt_arith::{matrix_vector_ntt, scalar_vector_ntt};
use crate::params::{ParameterSet, D, N};
use crate::poly::PolyVec;
use crate::sample::sample_in_ball;
use crate::serdes::{pk_decode, sig_decode, w1_encode};
use crate::sign::format_m_prime;
use crate::vecops::{inf_norm, inv_ntt_vec, ntt_vec, sub_vec, use_hint_vec};

/// FIPS 204, Algorithm 8 — ML-DSA.Verify_internal.
pub fn verify_internal<P: ParameterSet, const K: usize, const L: usize>(
    pk: &[u8],
    m_prime: &[u8],
    sig: &[u8],
) -> bool {
    // line 1-3: decode and if bad length or malformed hint (⊥) returns invalid.
    let (rho, t1) = match pk_decode::<P, K>(pk) {
        Ok(v) => v,
        Err(_) => return false,
    };
    let (c_tilde, z, h) = match sig_decode::<P, K, L>(sig) {
        Ok(v) => v,
        Err(_) => return false,
    };

    let a_hat = expand_a::<K, L>(&rho); // line 5

    // line 6
    let mut th = H::init();
    th.absorb(pk);
    let mut tr = [0u8; 64];
    th.finalize().squeeze(&mut tr);

    // line 7
    let mut hm = H::init();
    hm.absorb(&tr);
    hm.absorb(m_prime);
    let mut mu = [0u8; 64];
    hm.finalize().squeeze(&mut mu);

    let c = sample_in_ball::<P>(&c_tilde); // line 8

    // line 9
    let mut t1_scaled = PolyVec::<K>::zero();
    for i in 0..K {
        for j in 0..N {
            t1_scaled.v[i].coeffs[j] = t1.v[i].coeffs[j] << D;
        }
    }
    let az = matrix_vector_ntt(&a_hat, &ntt_vec(&z));
    let ct1 = scalar_vector_ntt(&ntt(&c), &ntt_vec(&t1_scaled));
    let w_prime = sub_vec(&inv_ntt_vec(&az), &inv_ntt_vec(&ct1));

    let w1_prime = use_hint_vec::<P, K>(&h, &w_prime); // line 10

    // line 12
    let mut hc = H::init();
    hc.absorb(&mu);
    hc.absorb(&w1_encode::<P, K>(&w1_prime));
    let c_tilde_prime = hc.finalize().squeeze_vec(P::C_TILDE_BYTES);

    // line 13
    inf_norm(&z) < P::GAMMA1 - P::BETA && c_tilde == c_tilde_prime
}

/// Improved-path [verify_internal]: identical structure and results, using the division-free NTT components.
pub fn verify_internal_fast<P: ParameterSet, const K: usize, const L: usize>(
    pk: &[u8],
    m_prime: &[u8],
    sig: &[u8],
) -> bool {
    use crate::ntt::ntt_fast;
    use crate::ntt_arith::{matrix_vector_ntt_fast, scalar_vector_ntt_fast};
    use crate::vecops::{inv_ntt_vec_fast, ntt_vec_fast};

    let (rho, t1) = match pk_decode::<P, K>(pk) {
        Ok(v) => v,
        Err(_) => return false,
    };
    let (c_tilde, z, h) = match sig_decode::<P, K, L>(sig) {
        Ok(v) => v,
        Err(_) => return false,
    };

    let a_hat = expand_a::<K, L>(&rho);

    let mut th = H::init();
    th.absorb(pk);
    let mut tr = [0u8; 64];
    th.finalize().squeeze(&mut tr);

    let mut hm = H::init();
    hm.absorb(&tr);
    hm.absorb(m_prime);
    let mut mu = [0u8; 64];
    hm.finalize().squeeze(&mut mu);

    let c = sample_in_ball::<P>(&c_tilde);

    let mut t1_scaled = PolyVec::<K>::zero();
    for i in 0..K {
        for j in 0..N {
            t1_scaled.v[i].coeffs[j] = t1.v[i].coeffs[j] << D;
        }
    }
    let az = matrix_vector_ntt_fast(&a_hat, &ntt_vec_fast(&z));
    let ct1 = scalar_vector_ntt_fast(&ntt_fast(&c), &ntt_vec_fast(&t1_scaled));
    let w_prime = sub_vec(&inv_ntt_vec_fast(&az), &inv_ntt_vec_fast(&ct1));

    let w1_prime = use_hint_vec::<P, K>(&h, &w_prime);

    let mut hc = H::init();
    hc.absorb(&mu);
    hc.absorb(&w1_encode::<P, K>(&w1_prime));
    let c_tilde_prime = hc.finalize().squeeze_vec(P::C_TILDE_BYTES);

    inf_norm(&z) < P::GAMMA1 - P::BETA && c_tilde == c_tilde_prime
}

/// Improved-path [verify].
pub fn verify_fast<P: ParameterSet, const K: usize, const L: usize>(
    pk: &[u8],
    m: &[u8],
    sig: &[u8],
    ctx: &[u8],
) -> bool {
    if ctx.len() > 255 {
        return false;
    }
    verify_internal_fast::<P, K, L>(pk, &format_m_prime(ctx, m), sig)
}

/// FIPS 204, Algorithm 3 — ML-DSA.Verify.
pub fn verify<P: ParameterSet, const K: usize, const L: usize>(
    pk: &[u8],
    m: &[u8],
    sig: &[u8],
    ctx: &[u8],
) -> bool {
    if ctx.len() > 255 {
        return false;
    }
    verify_internal::<P, K, L>(pk, &format_m_prime(ctx, m), sig)
}

/// FIPS 204, Algorithm 5 — HashML-DSA.Verify (pre-hash). Builds M' same as Algorithm 4 does, then goes to Verify_internal.
pub fn hash_verify<P: ParameterSet, const K: usize, const L: usize>(
    pk: &[u8],
    m: &[u8],
    sig: &[u8],
    ctx: &[u8],
    ph: crate::prehash::PreHash,
) -> bool {
    if ctx.len() > 255 {
        return false;
    }
    verify_internal::<P, K, L>(pk, &crate::sign::format_m_prime_prehash(ctx, m, ph), sig)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::keygen::key_gen_internal;
    use crate::params::MlDsa65;
    use crate::sign::sign_deterministic;

    const K: usize = MlDsa65::K;
    const L: usize = MlDsa65::L;

    #[test]
    fn sign_then_verify_round_trip() {
        let (pk, sk) = key_gen_internal::<MlDsa65, K, L>(&[0x55u8; 32]);
        let sig =
            sign_deterministic::<MlDsa65, K, L>(&sk, b"the quick brown fox", b"ctx").unwrap();
        assert!(verify::<MlDsa65, K, L>(&pk, b"the quick brown fox", &sig, b"ctx"));
    }

    #[test]
    fn verify_rejects_tampered_signature() {
        let (pk, sk) = key_gen_internal::<MlDsa65, K, L>(&[0x56u8; 32]);
        let mut sig = sign_deterministic::<MlDsa65, K, L>(&sk, b"msg", b"").unwrap();
        sig[100] ^= 1; // flip a bit
        assert!(!verify::<MlDsa65, K, L>(&pk, b"msg", &sig, b""));
    }

    #[test]
    fn verify_rejects_wrong_context_and_message() {
        let (pk, sk) = key_gen_internal::<MlDsa65, K, L>(&[0x57u8; 32]);
        let sig = sign_deterministic::<MlDsa65, K, L>(&sk, b"msg", b"ctxA").unwrap();
        assert!(!verify::<MlDsa65, K, L>(&pk, b"msg", &sig, b"ctxB"));
        assert!(!verify::<MlDsa65, K, L>(&pk, b"other", &sig, b"ctxA"));
    }

    #[test]
    fn verify_rejects_bad_lengths() {
        let (pk, _sk) = key_gen_internal::<MlDsa65, K, L>(&[0x58u8; 32]);
        assert!(!verify::<MlDsa65, K, L>(&pk, b"m", &[0u8; 10], b"")); // short sig
        assert!(!verify::<MlDsa65, K, L>(
            &[0u8; 10],
            b"m",
            &[0u8; MlDsa65::SIG_BYTES],
            b""
        )); // short pk
    }
}

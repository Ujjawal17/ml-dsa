//! FIPS 204 §6.3 / §5.3 — verification (Algorithms 8 and 3).
//!
//! Verify operates entirely on public data and adversary-controlled bytes, so it
//! length-checks `pk`/`σ` (returning false on mismatch) and rejects a malformed hint
//! (`sigDecode` → `⊥`). It reconstructs `w' = Az − c·t1·2^d` and checks the
//! commitment hash matches.
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
    // line 1-3: decode; bad length or malformed hint (⊥) ⇒ invalid.
    let (rho, t1) = match pk_decode::<P, K>(pk) {
        Ok(v) => v,
        Err(_) => return false,
    };
    let (c_tilde, z, h) = match sig_decode::<P, K, L>(sig) {
        Ok(v) => v,
        Err(_) => return false,
    };

    let a_hat = expand_a::<K, L>(&rho); // line 5

    // line 6: tr ← H(pk, 64)
    let mut th = H::init();
    th.absorb(pk);
    let mut tr = [0u8; 64];
    th.finalize().squeeze(&mut tr);

    // line 7: μ ← H(BytesToBits(tr) || M', 64)
    let mut hm = H::init();
    hm.absorb(&tr);
    hm.absorb(m_prime);
    let mut mu = [0u8; 64];
    hm.finalize().squeeze(&mut mu);

    let c = sample_in_ball::<P>(&c_tilde); // line 8

    // line 9: w' = NTT^-1(Â ∘ NTT(z) − NTT(c) ∘ NTT(t1·2^d))
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

    // line 12: c~' ← H(μ || w1Encode(w1'), λ/4)
    let mut hc = H::init();
    hc.absorb(&mu);
    hc.absorb(&w1_encode::<P, K>(&w1_prime));
    let c_tilde_prime = hc.finalize().squeeze_vec(P::C_TILDE_BYTES);

    // line 13: ‖z‖∞ < γ1 − β  AND  c~ = c~'
    inf_norm(&z) < P::GAMMA1 - P::BETA && c_tilde == c_tilde_prime
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

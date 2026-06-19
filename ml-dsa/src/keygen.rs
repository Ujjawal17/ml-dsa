//! FIPS 204 §6.1 / §5.1 — key generation (Algorithms 6 and 1).
//!
//! `key_gen_internal` is the CAVP-tested surface: deterministic in its 32-byte seed
//! `ξ`. `key_gen` is the public wrapper that draws `ξ` from an injected RNG.
use rand_core::{CryptoRng, RngCore};

use crate::expand::{expand_a, expand_s};
use crate::hash::H;
use crate::ntt_arith::matrix_vector_ntt;
use crate::params::{K, L};
use crate::serdes::{pk_encode, sk_encode};
use crate::vecops::{add_vec, inv_ntt_vec, ntt_vec, power2round_vec};

/// FIPS 204, Algorithm 6 — ML-DSA.KeyGen_internal: encoded `(pk, sk)` from seed `ξ`.
pub fn key_gen_internal(xi: &[u8; 32]) -> (Vec<u8>, Vec<u8>) {
    // line 1: (ρ, ρ', K) ← H(ξ || IntegerToBytes(k,1) || IntegerToBytes(l,1), 128)
    let mut h = H::init();
    h.absorb(xi);
    h.absorb(&[K as u8, L as u8]);
    let seed = h.finalize().squeeze_vec(128);
    let mut rho = [0u8; 32];
    rho.copy_from_slice(&seed[..32]);
    let mut rho_prime = [0u8; 64];
    rho_prime.copy_from_slice(&seed[32..96]);
    let mut k_seed = [0u8; 32];
    k_seed.copy_from_slice(&seed[96..128]);

    let a_hat = expand_a(&rho); // line 3
    let (s1, s2) = expand_s(&rho_prime); // line 4

    // line 5: t = NTT^-1(Â ∘ NTT(s1)) + s2
    let s1_hat = ntt_vec(&s1);
    let as1 = matrix_vector_ntt(&a_hat, &s1_hat);
    let t = add_vec(&inv_ntt_vec(&as1), &s2);

    let (t1, t0) = power2round_vec(&t); // line 6

    let pk = pk_encode(&rho, &t1); // line 8

    let mut th = H::init(); // line 9: tr ← H(pk, 64)
    th.absorb(&pk);
    let mut tr = [0u8; 64];
    th.finalize().squeeze(&mut tr);

    let sk = sk_encode(&rho, &k_seed, &tr, &s1, &s2, &t0); // line 10
    (pk, sk)
}

/// FIPS 204, Algorithm 1 — ML-DSA.KeyGen: draw `ξ` from the injected RNG, then expand.
pub fn key_gen<R: CryptoRng + RngCore>(rng: &mut R) -> (Vec<u8>, Vec<u8>) {
    let mut xi = [0u8; 32];
    rng.fill_bytes(&mut xi);
    key_gen_internal(&xi)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{ETA, PK_BYTES, SK_BYTES};
    use crate::serdes::{pk_decode, sk_decode};

    #[test]
    fn key_gen_internal_shapes_and_self_consistency() {
        let xi = [0x42u8; 32];
        let (pk, sk) = key_gen_internal(&xi);
        assert_eq!(pk.len(), PK_BYTES);
        assert_eq!(sk.len(), SK_BYTES);
        // deterministic in ξ
        assert_eq!(key_gen_internal(&xi), (pk.clone(), sk.clone()));
        // pk/sk decode, and the ρ stored in both must match.
        let (rho_pk, _t1) = pk_decode(&pk).unwrap();
        let (rho_sk, _k, tr, s1, s2, _t0) = sk_decode(&sk).unwrap();
        assert_eq!(rho_pk, rho_sk);
        // tr in sk must equal H(pk, 64).
        let mut th = H::init();
        th.absorb(&pk);
        let mut tr_expected = [0u8; 64];
        th.finalize().squeeze(&mut tr_expected);
        assert_eq!(tr, tr_expected);
        // secrets are in [-η, η].
        for p in s1.v.iter().chain(s2.v.iter()) {
            assert!(p.coeffs.iter().all(|&x| (-ETA..=ETA).contains(&x)));
        }
    }
}

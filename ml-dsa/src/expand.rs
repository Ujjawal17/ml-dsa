//! FIPS 204 §7.3 — ExpandA / ExpandS / ExpandMask (Algorithms 32–34).
//!
//! `ExpandA` samples the public matrix directly in the NTT domain (`Â`). `ExpandS`
//! samples the secrets `s1, s2`. `ExpandMask` derives the per-signature masking
//! vector `y` from a secret-derived seed (no rejection — `BitUnpack` always succeeds
//! because `γ1` is a power of two), so it is straightforwardly constant-time.
#![allow(clippy::needless_range_loop)]

use crate::encoding::{bit_unpack, bitlen, integer_to_bytes};
use crate::hash::H;
use crate::params::ParameterSet;
use crate::poly::{PolyMatNTT, PolyVec};
use crate::sample::{rej_bounded_poly, rej_ntt_poly};

/// FIPS 204, Algorithm 32 — ExpandA: the `K×L` matrix `Â` in the NTT domain.
pub fn expand_a<const K: usize, const L: usize>(rho: &[u8; 32]) -> PolyMatNTT<K, L> {
    let mut a = PolyMatNTT::<K, L>::zero();
    for r in 0..K {
        for s in 0..L {
            let mut seed = [0u8; 34];
            seed[..32].copy_from_slice(rho);
            seed[32] = s as u8; // IntegerToBytes(s, 1)
            seed[33] = r as u8; // IntegerToBytes(r, 1)
            a.rows[r].v[s] = rej_ntt_poly(&seed);
        }
    }
    a
}

/// FIPS 204, Algorithm 33 — ExpandS: secret vectors `s1 ∈ R^L`, `s2 ∈ R^K` in `[−η, η]`.
pub fn expand_s<P: ParameterSet, const K: usize, const L: usize>(
    rho: &[u8; 64],
) -> (PolyVec<L>, PolyVec<K>) {
    let mut s1 = PolyVec::<L>::zero();
    let mut s2 = PolyVec::<K>::zero();
    for r in 0..L {
        s1.v[r] = rej_bounded_poly::<P>(&bounded_seed(rho, r as u64));
    }
    for r in 0..K {
        s2.v[r] = rej_bounded_poly::<P>(&bounded_seed(rho, (r + L) as u64));
    }
    (s1, s2)
}

/// Build the 66-byte RejBoundedPoly seed `ρ || IntegerToBytes(nonce, 2)`.
fn bounded_seed(rho: &[u8; 64], nonce: u64) -> [u8; 66] {
    let mut seed = [0u8; 66];
    seed[..64].copy_from_slice(rho);
    let idx = integer_to_bytes(nonce, 2);
    seed[64] = idx[0];
    seed[65] = idx[1];
    seed
}

/// FIPS 204, Algorithm 34 — ExpandMask: masking vector `y ∈ R^L`, coeffs in `(−γ1, γ1]`.
pub fn expand_mask<P: ParameterSet, const L: usize>(rho: &[u8; 64], mu: u32) -> PolyVec<L> {
    let c = 1 + bitlen(P::GAMMA1 as u32 - 1) as usize; // 20 (or 18 for ML-DSA-44)
    let mut y = PolyVec::<L>::zero();
    for r in 0..L {
        let mut h = H::init();
        h.absorb(rho);
        h.absorb(&integer_to_bytes((mu + r as u32) as u64, 2)); // ρ' = ρ || (μ+r)
        let v = h.finalize().squeeze_vec(32 * c);
        y.v[r] = bit_unpack(&v, (P::GAMMA1 - 1) as u32, P::GAMMA1 as u32);
    }
    y
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{MlDsa65, ParameterSet, Q};

    const K: usize = MlDsa65::K;
    const L: usize = MlDsa65::L;

    #[test]
    fn expand_a_range_and_deterministic() {
        let rho = [5u8; 32];
        let a = expand_a::<K, L>(&rho);
        assert!(a.rows[0].v[0].coeffs.iter().all(|&x| (0..Q).contains(&x)));
        assert!(a.rows[K - 1].v[L - 1].coeffs.iter().all(|&x| (0..Q).contains(&x)));
        assert_eq!(a.rows[2].v[1].coeffs, expand_a::<K, L>(&rho).rows[2].v[1].coeffs);
    }

    #[test]
    fn expand_s_in_range() {
        let rho = [1u8; 64];
        let (s1, s2) = expand_s::<MlDsa65, K, L>(&rho);
        assert_eq!(s1.v.len(), L);
        assert_eq!(s2.v.len(), K);
        for p in s1.v.iter().chain(s2.v.iter()) {
            assert!(p.coeffs.iter().all(|&x| (-MlDsa65::ETA..=MlDsa65::ETA).contains(&x)));
        }
    }

    #[test]
    fn expand_mask_in_range() {
        let rho = [2u8; 64];
        let y = expand_mask::<MlDsa65, L>(&rho, 0);
        assert_eq!(y.v.len(), L);
        for p in y.v.iter() {
            assert!(p
                .coeffs
                .iter()
                .all(|&x| (-(MlDsa65::GAMMA1 - 1)..=MlDsa65::GAMMA1).contains(&x)));
        }
    }
}

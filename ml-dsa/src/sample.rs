//! FIPS 204 §7.3 — pseudorandom sampling (Algorithms 29–31).
//!
//! `RejNTTPoly`/`RejBoundedPoly` reject on **public** seed-derived bytes, so their
//! data-dependent loop counts are acceptable (the conventional position). The inner
//! `while j > i` loop of `SampleInBall` leaks only its *count* per signing call —
//! the documented accepted leakage (Appendix C bounds it).
#![allow(clippy::needless_range_loop)]

use crate::encoding::{bytes_to_bits, coeff_from_half_byte, coeff_from_three_bytes};
use crate::hash::{G, H};
use crate::params::{ParameterSet, N};
use crate::poly::{Poly, PolyNTT};

/// FIPS 204, Algorithm 29 — SampleInBall: `c ∈ R`, `τ` nonzero coeffs in `{−1, 1}`.
/// `rho` is the commitment hash `c~` (`λ/4` bytes).
pub fn sample_in_ball<P: ParameterSet>(rho: &[u8]) -> Poly {
    let tau = P::TAU;
    let mut c = Poly::zero();
    let mut reader = {
        let mut h = H::init();
        h.absorb(rho);
        h.finalize()
    };
    let mut s = [0u8; 8];
    reader.squeeze(&mut s); // line 4
    let hbits = bytes_to_bits(&s); // line 5: 64-bit sign string
    let mut byte = [0u8; 1];
    for i in (N - tau)..N {
        reader.squeeze(&mut byte); // line 7
        let mut j = byte[0] as usize;
        while j > i {
            // line 8: rejection sampling of a position in {0,…,i}
            reader.squeeze(&mut byte);
            j = byte[0] as usize;
        }
        c.coeffs[i] = c.coeffs[j]; // line 11
        c.coeffs[j] = if hbits[i + tau - N] == 1 { -1 } else { 1 }; // line 12: (-1)^h
    }
    c
}

/// FIPS 204, Algorithm 30 — RejNTTPoly: `â ∈ T_q` from a public 34-byte seed.
pub fn rej_ntt_poly(rho: &[u8; 34]) -> PolyNTT {
    let mut a = PolyNTT::zero();
    let mut reader = {
        let mut g = G::init();
        g.absorb(rho);
        g.finalize()
    };
    let mut j = 0usize;
    let mut s = [0u8; 3];
    while j < N {
        reader.squeeze(&mut s);
        if let Some(coeff) = coeff_from_three_bytes(s[0], s[1], s[2]) {
            a.coeffs[j] = coeff;
            j += 1;
        }
    }
    a
}

/// FIPS 204, Algorithm 31 — RejBoundedPoly: `a ∈ R`, coeffs in `[−η, η]`, 66-byte seed.
pub fn rej_bounded_poly<P: ParameterSet>(rho: &[u8; 66]) -> Poly {
    let mut a = Poly::zero();
    let mut reader = {
        let mut h = H::init();
        h.absorb(rho);
        h.finalize()
    };
    let mut j = 0usize;
    let mut z = [0u8; 1];
    while j < N {
        reader.squeeze(&mut z);
        let z0 = coeff_from_half_byte::<P>(z[0] % 16); // low nibble
        let z1 = coeff_from_half_byte::<P>(z[0] / 16); // high nibble
        if let Some(v) = z0 {
            a.coeffs[j] = v;
            j += 1;
        }
        if let Some(v) = z1 {
            if j < N {
                a.coeffs[j] = v;
                j += 1;
            }
        }
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{MlDsa44, MlDsa65, Q};

    #[test]
    fn sample_in_ball_weight_and_values() {
        let rho = [7u8; 48];
        let c = sample_in_ball::<MlDsa65>(&rho);
        let nonzero = c.coeffs.iter().filter(|&&x| x != 0).count();
        assert_eq!(nonzero, MlDsa65::TAU, "Hamming weight must be τ");
        assert!(c.coeffs.iter().all(|&x| (-1..=1).contains(&x)));
        // determinism
        assert_eq!(c.coeffs, sample_in_ball::<MlDsa65>(&rho).coeffs);
        // a different τ gives a different Hamming weight from the same seed
        let c44 = sample_in_ball::<MlDsa44>(&rho[..32]);
        assert_eq!(c44.coeffs.iter().filter(|&&x| x != 0).count(), MlDsa44::TAU);
    }

    #[test]
    fn rej_ntt_poly_range_and_deterministic() {
        let rho = [9u8; 34];
        let a = rej_ntt_poly(&rho);
        assert!(a.coeffs.iter().all(|&x| (0..Q).contains(&x)));
        assert_eq!(a.coeffs, rej_ntt_poly(&rho).coeffs);
    }

    #[test]
    fn rej_bounded_poly_in_range() {
        let rho = [3u8; 66];
        let a = rej_bounded_poly::<MlDsa65>(&rho);
        assert!(a.coeffs.iter().all(|&x| (-MlDsa65::ETA..=MlDsa65::ETA).contains(&x)));
        assert_eq!(a.coeffs, rej_bounded_poly::<MlDsa65>(&rho).coeffs);
        // η = 2 branch
        let a44 = rej_bounded_poly::<MlDsa44>(&rho);
        assert!(a44.coeffs.iter().all(|&x| (-MlDsa44::ETA..=MlDsa44::ETA).contains(&x)));
    }
}

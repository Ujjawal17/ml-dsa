//! FIPS 204 §7.4 — rounding and hints (Algorithms 35–40), plus componentwise
//! `Poly` variants.
//!
//! All inputs are reduced with `r mod q` to `[0, q)` first (spec line 1 of each),
//! then the centred `mod±` is taken where the spec requires it. The `Decompose`
//! special case `r⁺ − r₀ = q − 1` (lines 3–5) is the single most common
//! reference-implementation bug; it has a dedicated boundary test below.
//!
//! Indexed loops in the `Poly` variants mirror "for i from 0 to 255".
#![allow(clippy::needless_range_loop)]

use crate::field::mod_pm;
use crate::params::{D, GAMMA2, N, Q};
use crate::poly::Poly;

/// `2^d`.
const TWO_POW_D: i32 = 1 << D; // 8192
/// `2·gamma2`.
const TWO_GAMMA2: i32 = 2 * GAMMA2; // 523776

/// FIPS 204, Algorithm 35 — Power2Round. Returns `(r1, r0)` with `r ≡ r1·2^d + r0`.
pub fn power2round(r: i32) -> (i32, i32) {
    let rp = r.rem_euclid(Q); // r mod q
    let r0 = mod_pm(rp, TWO_POW_D); // r+ mod± 2^d
    ((rp - r0) / TWO_POW_D, r0)
}

/// FIPS 204, Algorithm 36 — Decompose. Returns `(r1, r0)` with `r ≡ r1·(2·gamma2) + r0`.
pub fn decompose(r: i32) -> (i32, i32) {
    let rp = r.rem_euclid(Q); // r mod q
    let mut r0 = mod_pm(rp, TWO_GAMMA2); // r+ mod± 2γ2
    let r1;
    if rp - r0 == Q - 1 {
        // special case: r1 wraps to 0 and r0 absorbs the -1 (lines 3-5)
        r1 = 0;
        r0 -= 1;
    } else {
        r1 = (rp - r0) / TWO_GAMMA2;
    }
    (r1, r0)
}

/// FIPS 204, Algorithm 37 — HighBits: the `r1` part of `Decompose(r)`.
pub fn high_bits(r: i32) -> i32 {
    decompose(r).0
}

/// FIPS 204, Algorithm 38 — LowBits: the `r0` part of `Decompose(r)`.
pub fn low_bits(r: i32) -> i32 {
    decompose(r).1
}

/// FIPS 204, Algorithm 39 — MakeHint: does adding `z` change the high bits of `r`?
pub fn make_hint(z: i32, r: i32) -> bool {
    high_bits(r) != high_bits(r + z)
}

/// FIPS 204, Algorithm 40 — UseHint: high bits of `r`, adjusted by hint `h`.
pub fn use_hint(h: bool, r: i32) -> i32 {
    let m = (Q - 1) / TWO_GAMMA2; // 16 for ML-DSA-65
    let (r1, r0) = decompose(r);
    if h && r0 > 0 {
        (r1 + 1).rem_euclid(m)
    } else if h && r0 <= 0 {
        (r1 - 1).rem_euclid(m)
    } else {
        r1
    }
}

// --- Componentwise `Poly` variants (FIPS 204 §7.4 applies these coefficientwise) ---

/// Power2Round applied to every coefficient; returns `(t1, t0)`.
pub fn power2round_poly(r: &Poly) -> (Poly, Poly) {
    let mut t1 = Poly::zero();
    let mut t0 = Poly::zero();
    for i in 0..N {
        let (a, b) = power2round(r.coeffs[i]);
        t1.coeffs[i] = a;
        t0.coeffs[i] = b;
    }
    (t1, t0)
}

/// HighBits applied to every coefficient.
pub fn high_bits_poly(r: &Poly) -> Poly {
    let mut out = Poly::zero();
    for i in 0..N {
        out.coeffs[i] = high_bits(r.coeffs[i]);
    }
    out
}

/// LowBits applied to every coefficient.
pub fn low_bits_poly(r: &Poly) -> Poly {
    let mut out = Poly::zero();
    for i in 0..N {
        out.coeffs[i] = low_bits(r.coeffs[i]);
    }
    out
}

/// MakeHint applied componentwise; hint coefficients are `0`/`1`.
pub fn make_hint_poly(z: &Poly, r: &Poly) -> Poly {
    let mut h = Poly::zero();
    for i in 0..N {
        h.coeffs[i] = make_hint(z.coeffs[i], r.coeffs[i]) as i32;
    }
    h
}

/// UseHint applied componentwise; `h` coefficients are interpreted as `0`/`1`.
pub fn use_hint_poly(h: &Poly, r: &Poly) -> Poly {
    let mut out = Poly::zero();
    for i in 0..N {
        out.coeffs[i] = use_hint(h.coeffs[i] != 0, r.coeffs[i]);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    struct XorShift(u64);
    impl XorShift {
        fn next_u32(&mut self) -> u32 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            (x >> 32) as u32
        }
        fn coeff(&mut self) -> i32 {
            (self.next_u32() % Q as u32) as i32
        }
    }

    #[test]
    fn decompose_identity_and_range() {
        let mut rng = XorShift(0xa1b2_c3d4_e5f6_0718);
        let m = (Q - 1) / TWO_GAMMA2;
        for _ in 0..5000 {
            let r = rng.coeff();
            let (r1, r0) = decompose(r);
            // r ≡ r1·2γ2 + r0 (mod q)
            assert_eq!((r1 as i64 * TWO_GAMMA2 as i64 + r0 as i64).rem_euclid(Q as i64) as i32, r);
            // r1 in [0, (q-1)/2γ2)
            assert!(r1 >= 0 && r1 < m);
        }
    }

    #[test]
    fn decompose_q_minus_1_special_case() {
        // The notorious edge: r⁺ − r₀ = q − 1 must give r1 = 0, r0 = -1.
        assert_eq!(decompose(Q - 1), (0, -1));
    }

    #[test]
    fn power2round_identity() {
        let mut rng = XorShift(0x0f0e_0d0c_0b0a_0908);
        for _ in 0..5000 {
            let r = rng.coeff();
            let (r1, r0) = power2round(r);
            // r⁺ = r1·2^d + r0 exactly.
            assert_eq!(r1 * TWO_POW_D + r0, r);
            assert!(r0 > -TWO_POW_D / 2 && r0 <= TWO_POW_D / 2);
        }
    }

    #[test]
    fn hint_correctness_invariant() {
        // For |z| ≤ γ2: UseHint(MakeHint(z, r), r + z) = HighBits(r).
        let mut rng = XorShift(0xfeed_face_1234_5678);
        for _ in 0..5000 {
            let r = rng.coeff();
            let z = (rng.next_u32() % (2 * GAMMA2 as u32 + 1)) as i32 - GAMMA2;
            let h = make_hint(z, r);
            assert_eq!(use_hint(h, r + z), high_bits(r));
        }
    }
}

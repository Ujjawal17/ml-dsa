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
use crate::params::{ParameterSet, D, N, Q};
use crate::poly::Poly;

/// `2^d`.
const TWO_POW_D: i32 = 1 << D; // 8192

/// FIPS 204, Algorithm 35 — Power2Round. Returns `(r1, r0)` with `r ≡ r1·2^d + r0`.
pub fn power2round(r: i32) -> (i32, i32) {
    let rp = r.rem_euclid(Q); // r mod q
    let r0 = mod_pm(rp, TWO_POW_D); // r+ mod± 2^d
    ((rp - r0) / TWO_POW_D, r0)
}

/// FIPS 204, Algorithm 36 — Decompose. Returns `(r1, r0)` with `r ≡ r1·(2·gamma2) + r0`.
pub fn decompose<P: ParameterSet>(r: i32) -> (i32, i32) {
    let two_gamma2 = 2 * P::GAMMA2;
    let rp = r.rem_euclid(Q); // r mod q
    let mut r0 = mod_pm(rp, two_gamma2); // r+ mod± 2γ2
    let r1;
    if rp - r0 == Q - 1 {
        // special case: r1 wraps to 0 and r0 absorbs the -1 (lines 3-5)
        r1 = 0;
        r0 -= 1;
    } else {
        r1 = (rp - r0) / two_gamma2;
    }
    (r1, r0)
}

/// FIPS 204, Algorithm 37 — HighBits: the `r1` part of `Decompose(r)`.
pub fn high_bits<P: ParameterSet>(r: i32) -> i32 {
    decompose::<P>(r).0
}

/// FIPS 204, Algorithm 38 — LowBits: the `r0` part of `Decompose(r)`.
pub fn low_bits<P: ParameterSet>(r: i32) -> i32 {
    decompose::<P>(r).1
}

/// FIPS 204, Algorithm 39 — MakeHint: does adding `z` change the high bits of `r`?
pub fn make_hint<P: ParameterSet>(z: i32, r: i32) -> bool {
    high_bits::<P>(r) != high_bits::<P>(r + z)
}

/// FIPS 204, Algorithm 40 — UseHint: high bits of `r`, adjusted by hint `h`.
pub fn use_hint<P: ParameterSet>(h: bool, r: i32) -> i32 {
    let m = (Q - 1) / (2 * P::GAMMA2); // 16 for ML-DSA-65/87, 44 for ML-DSA-44
    let (r1, r0) = decompose::<P>(r);
    if h && r0 > 0 {
        (r1 + 1).rem_euclid(m)
    } else if h && r0 <= 0 {
        (r1 - 1).rem_euclid(m)
    } else {
        r1
    }
}

// --- Improved path: branchless (constant-time) variants ---
//
// The baseline mirrors the spec's `if`/`else` (the `q−1` special case, the centred
// reduction) and reduces with division. These variants compute the same values with
// no secret-dependent branch and no division; the `q−1` case is folded into the
// final masked subtract. Exhaustive tests below pin them to the baseline over the
// entire input range `[0, q)` for every γ2 in use.

/// Branchless [`decompose`]: value-equal for all inputs (exhaustively tested).
///
/// `r1 ≈ round(r⁺ / 2γ2)` is computed by fixed-point multiply — the constants are
/// `1025 ≈ 2^22/(2γ2/2^7)` for `γ2 = (q−1)/32` and `11275 ≈ 2^24/(2γ2/2^7)` for
/// `γ2 = (q−1)/88` (the two classes FIPS 204 uses) — then `r0 = r⁺ − r1·2γ2` is
/// re-centred with a masked subtract of `q`, which also absorbs the spec's
/// `r⁺ − r0 = q − 1` special case (lines 3–5 of Algorithm 36).
pub fn decompose_ct<P: ParameterSet>(r: i32) -> (i32, i32) {
    use crate::ct::gt_mask;
    use crate::field::to_canonical;
    let a = to_canonical(r); // r mod q, branchless
    let mut a1 = (a + 127) >> 7;
    if P::GAMMA2 == (Q - 1) / 32 {
        a1 = (a1 * 1025 + (1 << 21)) >> 22;
        a1 &= 15;
    } else if P::GAMMA2 == (Q - 1) / 88 {
        a1 = (a1 * 11275 + (1 << 23)) >> 24;
        a1 ^= ((43 - a1) >> 31) & a1;
    } else {
        unreachable!("FIPS 204 defines only gamma2 = (q-1)/32 and (q-1)/88");
    }
    let mut a0 = a - a1 * 2 * P::GAMMA2;
    a0 -= gt_mask(a0, (Q - 1) / 2) & Q;
    (a1, a0)
}

/// Branchless [`high_bits`].
pub fn high_bits_ct<P: ParameterSet>(r: i32) -> i32 {
    decompose_ct::<P>(r).0
}

/// Branchless [`low_bits`].
pub fn low_bits_ct<P: ParameterSet>(r: i32) -> i32 {
    decompose_ct::<P>(r).1
}

/// Branchless [`make_hint`]: `0`/`1` instead of `bool`, computed by masked
/// comparison of the two high-bits values.
pub fn make_hint_ct<P: ParameterSet>(z: i32, r: i32) -> i32 {
    use crate::ct::ne_bit;
    ne_bit(high_bits_ct::<P>(r), high_bits_ct::<P>(r + z))
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
pub fn high_bits_poly<P: ParameterSet>(r: &Poly) -> Poly {
    let mut out = Poly::zero();
    for i in 0..N {
        out.coeffs[i] = high_bits::<P>(r.coeffs[i]);
    }
    out
}

/// LowBits applied to every coefficient.
pub fn low_bits_poly<P: ParameterSet>(r: &Poly) -> Poly {
    let mut out = Poly::zero();
    for i in 0..N {
        out.coeffs[i] = low_bits::<P>(r.coeffs[i]);
    }
    out
}

/// MakeHint applied componentwise; hint coefficients are `0`/`1`.
pub fn make_hint_poly<P: ParameterSet>(z: &Poly, r: &Poly) -> Poly {
    let mut h = Poly::zero();
    for i in 0..N {
        h.coeffs[i] = make_hint::<P>(z.coeffs[i], r.coeffs[i]) as i32;
    }
    h
}

/// UseHint applied componentwise; `h` coefficients are interpreted as `0`/`1`.
pub fn use_hint_poly<P: ParameterSet>(h: &Poly, r: &Poly) -> Poly {
    let mut out = Poly::zero();
    for i in 0..N {
        out.coeffs[i] = use_hint::<P>(h.coeffs[i] != 0, r.coeffs[i]);
    }
    out
}

/// Branchless [`high_bits_poly`].
pub fn high_bits_poly_ct<P: ParameterSet>(r: &Poly) -> Poly {
    let mut out = Poly::zero();
    for i in 0..N {
        out.coeffs[i] = high_bits_ct::<P>(r.coeffs[i]);
    }
    out
}

/// Branchless [`low_bits_poly`].
pub fn low_bits_poly_ct<P: ParameterSet>(r: &Poly) -> Poly {
    let mut out = Poly::zero();
    for i in 0..N {
        out.coeffs[i] = low_bits_ct::<P>(r.coeffs[i]);
    }
    out
}

/// Branchless [`make_hint_poly`].
pub fn make_hint_poly_ct<P: ParameterSet>(z: &Poly, r: &Poly) -> Poly {
    let mut h = Poly::zero();
    for i in 0..N {
        h.coeffs[i] = make_hint_ct::<P>(z.coeffs[i], r.coeffs[i]);
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{MlDsa44, MlDsa65, MlDsa87};

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

    fn decompose_identity_and_range_for<P: ParameterSet>() {
        let mut rng = XorShift(0xa1b2_c3d4_e5f6_0718);
        let two_gamma2 = 2 * P::GAMMA2;
        let m = (Q - 1) / two_gamma2;
        for _ in 0..5000 {
            let r = rng.coeff();
            let (r1, r0) = decompose::<P>(r);
            // r ≡ r1·2γ2 + r0 (mod q)
            assert_eq!((r1 as i64 * two_gamma2 as i64 + r0 as i64).rem_euclid(Q as i64) as i32, r);
            // r1 in [0, (q-1)/2γ2)
            assert!(r1 >= 0 && r1 < m);
        }
    }

    #[test]
    fn decompose_identity_and_range() {
        decompose_identity_and_range_for::<MlDsa44>();
        decompose_identity_and_range_for::<MlDsa65>();
        decompose_identity_and_range_for::<MlDsa87>();
    }

    #[test]
    fn decompose_q_minus_1_special_case() {
        // The notorious edge: r⁺ − r₀ = q − 1 must give r1 = 0, r0 = -1.
        assert_eq!(decompose::<MlDsa65>(Q - 1), (0, -1));
        assert_eq!(decompose::<MlDsa44>(Q - 1), (0, -1));
        assert_eq!(decompose::<MlDsa87>(Q - 1), (0, -1));
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

    fn hint_correctness_invariant_for<P: ParameterSet>() {
        // For |z| ≤ γ2: UseHint(MakeHint(z, r), r + z) = HighBits(r).
        let mut rng = XorShift(0xfeed_face_1234_5678);
        for _ in 0..5000 {
            let r = rng.coeff();
            let z = (rng.next_u32() % (2 * P::GAMMA2 as u32 + 1)) as i32 - P::GAMMA2;
            let h = make_hint::<P>(z, r);
            assert_eq!(use_hint::<P>(h, r + z), high_bits::<P>(r));
        }
    }

    #[test]
    fn hint_correctness_invariant() {
        hint_correctness_invariant_for::<MlDsa44>();
        hint_correctness_invariant_for::<MlDsa65>();
        hint_correctness_invariant_for::<MlDsa87>();
    }

    /// The branchless decompose must be value-equal to the baseline over the
    /// ENTIRE input domain [0, q) — exhaustive, not sampled, because the fixed-
    /// point constants only "round correctly" if they do so for every input.
    /// MlDsa44 covers γ2 = (q−1)/88; MlDsa65 covers (q−1)/32 (shared with 87).
    fn decompose_ct_exhaustive_for<P: ParameterSet>() {
        for r in 0..Q {
            assert_eq!(decompose_ct::<P>(r), decompose::<P>(r), "r = {r}");
        }
    }

    #[test]
    fn decompose_ct_exhaustive_gamma2_32() {
        decompose_ct_exhaustive_for::<MlDsa65>();
    }

    #[test]
    fn decompose_ct_exhaustive_gamma2_88() {
        decompose_ct_exhaustive_for::<MlDsa44>();
    }

    /// Branchless decompose must also handle non-canonical inputs like the
    /// baseline does (both reduce mod q first).
    #[test]
    fn decompose_ct_non_canonical_inputs() {
        for r in [-1, -Q, Q, Q + 1, 2 * Q - 1, i32::MAX, i32::MIN] {
            assert_eq!(decompose_ct::<MlDsa65>(r), decompose::<MlDsa65>(r), "r = {r}");
            assert_eq!(decompose_ct::<MlDsa44>(r), decompose::<MlDsa44>(r), "r = {r}");
        }
    }

    #[test]
    fn make_hint_ct_equals_baseline() {
        let mut rng = XorShift(0x3141_5926_5358_9793);
        for _ in 0..100_000 {
            let r = rng.coeff();
            let z = rng.coeff() - Q / 2;
            assert_eq!(
                make_hint_ct::<MlDsa65>(z, r),
                make_hint::<MlDsa65>(z, r) as i32,
                "z = {z}, r = {r}"
            );
            assert_eq!(
                make_hint_ct::<MlDsa44>(z, r),
                make_hint::<MlDsa44>(z, r) as i32,
                "z = {z}, r = {r}"
            );
        }
    }
}

//! `Z_q` modular arithmetic primitives.
//! Baseline uses plain hardware division (`rem_euclid`), while the optimized path uses division-free Montgomery reduction.

use crate::ct::caddq;
use crate::params::Q;

/// Reduce a wide value into the canonical range `[0, q)`.
#[inline]
pub fn reduce_q(x: i64) -> i32 {
    x.rem_euclid(Q as i64) as i32
}

//Improved-path reducers (division-free)

/// q^{-1} mod 2^32
const QINV: u32 = 58_728_449;

/// R^2 mod q where R = 2^32
pub(crate) const R2: i64 = ((1u128 << 64) % Q as u128) as i64;

/// FIPS 204, Algorithm 49 — MontgomeryReduce
#[inline]
pub fn montgomery_reduce(a: i64) -> i32 {
    debug_assert!(a.unsigned_abs() < (1u64 << 31) * Q as u64);
    let t = (a as u32).wrapping_mul(QINV) as i32; // a·q^{-1} mod± 2^32
    ((a - t as i64 * Q as i64) >> 32) as i32
}

///for any i32, returns r ≡ a (mod q) with |r| < q, via one shift and one multiply.
#[inline]
pub(crate) fn reduce32(a: i32) -> i32 {
    let t = ((a as i64 + (1 << 22)) >> 23) as i32; // ≈ round(a / 2^23)
    (a as i64 - t as i64 * Q as i64) as i32
}

/// Branchless canonicalization to [0, q) without division and branch.
#[inline]
pub(crate) fn to_canonical(a: i32) -> i32 {
    caddq(reduce32(a))
}

///the representative of r in (-alpha/2, alpha/2].
#[inline]
pub fn mod_pm(r: i32, alpha: i32) -> i32 {
    let mut x = r.rem_euclid(alpha);
    if x > alpha / 2 {
        x -= alpha;
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Q;

    #[test]
    fn reduce_q_handles_negatives_and_overflow() {
        assert_eq!(reduce_q(-1), Q - 1);
        assert_eq!(reduce_q(Q as i64), 0);
        assert_eq!(reduce_q(2 * Q as i64 + 7), 7);
        assert_eq!(reduce_q(-(Q as i64) - 3), Q - 3);
    }

    #[test]
    fn mod_pm_is_centred() {
        // alpha = 8 -> representatives in (-4, 4].
        assert_eq!(mod_pm(0, 8), 0);
        assert_eq!(mod_pm(4, 8), 4);
        assert_eq!(mod_pm(5, 8), -3);
        assert_eq!(mod_pm(7, 8), -1);
        assert_eq!(mod_pm(8, 8), 0);
    }

    struct XorShift(u64);
    impl XorShift {
        fn next_u64(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
    }

    #[test]
    fn qinv_is_inverse_of_q_mod_2_32() {
        assert_eq!((Q as u32).wrapping_mul(QINV), 1);
        // R2 = 2^64 mod q, in [0, q).
        assert!((0..Q as i64).contains(&R2));
        assert_eq!(R2 as i128, (1i128 << 64) % Q as i128);
    }

    #[test]
    fn montgomery_reduce_is_congruent_and_bounded() {
        let bound = (1i64 << 31) * Q as i64;
        let mut rng = XorShift(0x0dd0_51ac_e997_1234);
        let check = |a: i64| {
            let r = montgomery_reduce(a);
            assert!(r.unsigned_abs() < Q as u32, "|r| < q for a = {a}");
            // r · 2^32 ≡ a (mod q)
            assert_eq!(((r as i128) * (1i128 << 32) - a as i128).rem_euclid(Q as i128), 0);
        };
        for a in [0, 1, -1, Q as i64, -(Q as i64), bound - 1, -(bound - 1)] {
            check(a);
        }
        for _ in 0..100_000 {
            check(rng.next_u64() as i64 % bound);
        }
    }

    #[test]
    fn to_canonical_equals_rem_euclid_on_full_i32_range() {
        let mut rng = XorShift(0x5eed_5eed_5eed_5eed);
        let check = |a: i32| {
            assert_eq!(to_canonical(a), a.rem_euclid(Q), "a = {a}");
        };
        for a in [0, 1, -1, Q - 1, Q, Q + 1, -Q, -(Q - 1), i32::MAX, i32::MIN, i32::MIN + 1] {
            check(a);
        }
        for _ in 0..1_000_000 {
            check(rng.next_u64() as i32);
        }
    }
}

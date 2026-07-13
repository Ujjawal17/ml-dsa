//! Constant-time helpers. Secret-dependent code must use these rather than a data
//! `if`, so that running time does not depend on secret values — the property the
//! Part 2 side-channel chapter measures. These are branchless: on `i32`, `a >> 31`
//! is all-ones iff `a` is negative, which gives a mask with no branch.
//!
//! (Allowed dead-code in Phase 0: `field.rs` in Phase 1 is the first consumer.)
#![allow(dead_code)]

use crate::params::Q;

/// Conditional add of `q`: returns `a + q` if `a < 0`, else `a`. Normalizes a
/// possibly-negative centred value toward `[0, q)` without a data-dependent branch.
#[inline(always)]
pub(crate) fn caddq(a: i32) -> i32 {
    a + ((a >> 31) & Q)
}

/// Conditional subtract of `q`: returns `a - q` if `a >= q`, else `a`.
#[inline(always)]
pub(crate) fn csubq(a: i32) -> i32 {
    let b = a - Q;
    b + ((b >> 31) & Q)
}

/// All-ones mask iff `a > b` (signed), else zero — without a branch.
/// Valid when `b - a` does not overflow `i32` (all uses here are values in
/// `(-q, q)`-scale ranges, far from the `i32` limits).
#[inline(always)]
pub(crate) fn gt_mask(a: i32, b: i32) -> i32 {
    (b - a) >> 31
}

/// `1` iff `a != b`, else `0` — without a branch: `x | -x` has its sign bit set
/// exactly when `x != 0`.
#[inline(always)]
pub(crate) fn ne_bit(a: i32, b: i32) -> i32 {
    let x = a ^ b;
    ((x | x.wrapping_neg()) >> 31) & 1
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Q;

    #[test]
    fn caddq_normalizes_negatives() {
        assert_eq!(caddq(-1), Q - 1);
        assert_eq!(caddq(0), 0);
        assert_eq!(caddq(5), 5);
        assert_eq!(caddq(-(Q - 1)), 1);
    }

    #[test]
    fn csubq_reduces_at_q() {
        assert_eq!(csubq(Q), 0);
        assert_eq!(csubq(Q + 5), 5);
        assert_eq!(csubq(Q - 1), Q - 1);
        assert_eq!(csubq(0), 0);
    }
}

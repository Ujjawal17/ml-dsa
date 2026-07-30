//! Constant-time helpers. Secret-dependent code use these rather than a data if, so that running time does not depend on secret values. 
//! These are branchless on i32, a >> 31 is all-ones iff a is negative, which gives a mask with no branch.

#![allow(dead_code)] //to avoid unused function warnings

use crate::params::Q;

/// conitional add of q and normalize towards range [0,q) by adding Q
#[inline(always)]
pub(crate) fn caddq(a: i32) -> i32 {
    a + ((a >> 31) & Q) 
}

/// conditional subtract of q
#[inline(always)]
pub(crate) fn csubq(a: i32) -> i32 {
    let b = a - Q;
    b + ((b >> 31) & Q)
}

/// All-ones mask iff a > b (signed), else zero
#[inline(always)]
pub(crate) fn gt_mask(a: i32, b: i32) -> i32 {
    (b - a) >> 31
}

/// 1 iff a != b, else 0
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

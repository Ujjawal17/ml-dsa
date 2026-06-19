//! `Z_q` arithmetic primitives for the faithful baseline.
//!
//! The baseline works in the **canonical range `[0, q)`** and reduces with plain
//! `mod q` (FIPS 204's NTT pseudocode uses `mod q` directly). Montgomery / lazy
//! reduction is the Part 2 optimization (FIPS 204 Appendix A, Algorithm 49), not
//! used here. The centred representative `mod±` is provided for the few places the
//! spec calls for it explicitly (Power2Round, Decompose).

use crate::params::Q;

/// Reduce a wide value into the canonical range `[0, q)`.
#[inline]
pub fn reduce_q(x: i64) -> i32 {
    x.rem_euclid(Q as i64) as i32
}

/// Centred modulo `mod± alpha`: the representative of `r` in `(-alpha/2, alpha/2]`.
/// `alpha` is assumed positive and even (it is `2^d` or `2*gamma2` in this crate).
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
}

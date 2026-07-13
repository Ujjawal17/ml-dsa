//! Vector lifts of the poly-level operations, shared by keygen / sign / verify.
//! FIPS 204 applies NTT, rounding, and hints componentwise over `R^k` / `R^l`
//! (§7.4, §7.6); these helpers express that lifting once.
#![allow(clippy::needless_range_loop)]

use crate::field::{mod_pm, reduce_q};
use crate::ntt::{inv_ntt, ntt};
use crate::params::{ParameterSet, N, Q};
use crate::poly::{PolyVec, PolyVecNTT};
use crate::rounding::{
    high_bits_poly, low_bits_poly, make_hint_poly, power2round_poly, use_hint_poly,
};

/// NTT each polynomial of a vector.
pub fn ntt_vec<const M: usize>(v: &PolyVec<M>) -> PolyVecNTT<M> {
    let mut out = PolyVecNTT::<M>::zero();
    for i in 0..M {
        out.v[i] = ntt(&v.v[i]);
    }
    out
}

/// Inverse-NTT each polynomial of a vector.
pub fn inv_ntt_vec<const M: usize>(v: &PolyVecNTT<M>) -> PolyVec<M> {
    let mut out = PolyVec::<M>::zero();
    for i in 0..M {
        out.v[i] = inv_ntt(&v.v[i]);
    }
    out
}

/// [`ntt_vec`] on the improved (division-free) component.
pub fn ntt_vec_fast<const M: usize>(v: &PolyVec<M>) -> PolyVecNTT<M> {
    let mut out = PolyVecNTT::<M>::zero();
    for i in 0..M {
        out.v[i] = crate::ntt::ntt_fast(&v.v[i]);
    }
    out
}

/// [`inv_ntt_vec`] on the improved (division-free) component.
pub fn inv_ntt_vec_fast<const M: usize>(v: &PolyVecNTT<M>) -> PolyVec<M> {
    let mut out = PolyVec::<M>::zero();
    for i in 0..M {
        out.v[i] = crate::ntt::inv_ntt_fast(&v.v[i]);
    }
    out
}

/// Componentwise `a + b` in `R^M` (reduced to `[0, q)`).
pub fn add_vec<const M: usize>(a: &PolyVec<M>, b: &PolyVec<M>) -> PolyVec<M> {
    let mut out = PolyVec::<M>::zero();
    for i in 0..M {
        for j in 0..N {
            out.v[i].coeffs[j] = reduce_q(a.v[i].coeffs[j] as i64 + b.v[i].coeffs[j] as i64);
        }
    }
    out
}

/// Componentwise `a - b` in `R^M` (reduced to `[0, q)`).
pub fn sub_vec<const M: usize>(a: &PolyVec<M>, b: &PolyVec<M>) -> PolyVec<M> {
    let mut out = PolyVec::<M>::zero();
    for i in 0..M {
        for j in 0..N {
            out.v[i].coeffs[j] = reduce_q(a.v[i].coeffs[j] as i64 - b.v[i].coeffs[j] as i64);
        }
    }
    out
}

/// Componentwise negation `-a` in `R^M` (reduced to `[0, q)`).
pub fn neg_vec<const M: usize>(a: &PolyVec<M>) -> PolyVec<M> {
    let mut out = PolyVec::<M>::zero();
    for i in 0..M {
        for j in 0..N {
            out.v[i].coeffs[j] = reduce_q(-(a.v[i].coeffs[j] as i64));
        }
    }
    out
}

/// Map each coefficient to its centred representative `mod± q` (in `(-q/2, q/2]`).
pub fn center_vec<const M: usize>(a: &PolyVec<M>) -> PolyVec<M> {
    let mut out = PolyVec::<M>::zero();
    for i in 0..M {
        for j in 0..N {
            out.v[i].coeffs[j] = mod_pm(a.v[i].coeffs[j], Q);
        }
    }
    out
}

/// Power2Round over a vector; returns `(t1, t0)`.
pub fn power2round_vec<const M: usize>(v: &PolyVec<M>) -> (PolyVec<M>, PolyVec<M>) {
    let mut t1 = PolyVec::<M>::zero();
    let mut t0 = PolyVec::<M>::zero();
    for i in 0..M {
        let (a, b) = power2round_poly(&v.v[i]);
        t1.v[i] = a;
        t0.v[i] = b;
    }
    (t1, t0)
}

/// HighBits over a vector.
pub fn high_bits_vec<P: ParameterSet, const M: usize>(v: &PolyVec<M>) -> PolyVec<M> {
    let mut out = PolyVec::<M>::zero();
    for i in 0..M {
        out.v[i] = high_bits_poly::<P>(&v.v[i]);
    }
    out
}

/// LowBits over a vector.
pub fn low_bits_vec<P: ParameterSet, const M: usize>(v: &PolyVec<M>) -> PolyVec<M> {
    let mut out = PolyVec::<M>::zero();
    for i in 0..M {
        out.v[i] = low_bits_poly::<P>(&v.v[i]);
    }
    out
}

/// MakeHint over vectors (componentwise), giving a 0/1 hint vector.
pub fn make_hint_vec<P: ParameterSet, const M: usize>(
    z: &PolyVec<M>,
    r: &PolyVec<M>,
) -> PolyVec<M> {
    let mut out = PolyVec::<M>::zero();
    for i in 0..M {
        out.v[i] = make_hint_poly::<P>(&z.v[i], &r.v[i]);
    }
    out
}

/// UseHint over vectors (componentwise).
pub fn use_hint_vec<P: ParameterSet, const M: usize>(
    h: &PolyVec<M>,
    r: &PolyVec<M>,
) -> PolyVec<M> {
    let mut out = PolyVec::<M>::zero();
    for i in 0..M {
        out.v[i] = use_hint_poly::<P>(&h.v[i], &r.v[i]);
    }
    out
}

/// Centred infinity norm `‖v‖∞` (each coefficient reduced `mod± q` first).
pub fn inf_norm<const M: usize>(v: &PolyVec<M>) -> i32 {
    let mut m = 0;
    for i in 0..M {
        for j in 0..N {
            let a = mod_pm(v.v[i].coeffs[j], Q).abs();
            if a > m {
                m = a;
            }
        }
    }
    m
}

/// Number of nonzero (`1`) coefficients across a hint vector.
pub fn count_ones<const M: usize>(v: &PolyVec<M>) -> usize {
    let mut c = 0;
    for i in 0..M {
        for j in 0..N {
            if v.v[i].coeffs[j] != 0 {
                c += 1;
            }
        }
    }
    c
}

// --- Improved path: branchless variants ---

/// Branchless `‖v‖∞ >= bound`: value-equal to `inf_norm(v) >= bound` (tested),
/// with no per-coefficient branch — a violation mask is OR-accumulated and only
/// the *aggregate* answer becomes a branch at the caller, which is exactly the
/// accept/reject decision whose count is the documented rejection-loop leakage.
pub fn exceeds_bound_ct<const M: usize>(v: &PolyVec<M>, bound: i32) -> bool {
    use crate::ct::gt_mask;
    use crate::field::to_canonical;
    let half = (Q - 1) / 2;
    let mut violation = 0i32;
    for i in 0..M {
        for j in 0..N {
            let x = to_canonical(v.v[i].coeffs[j]); // [0, q)
            // |x mod± q| = x if x <= (q-1)/2, else q - x — selected by mask.
            let m = gt_mask(x, half);
            let ax = (x & !m) | ((Q - x) & m);
            violation |= gt_mask(ax, bound - 1); // ax >= bound
        }
    }
    violation != 0
}

/// Branchless [`count_ones`]. Contract: coefficients must be `0`/`1` (true for
/// every hint vector by construction), so counting is plain summation.
pub fn count_ones_ct<const M: usize>(v: &PolyVec<M>) -> usize {
    let mut c = 0usize;
    for i in 0..M {
        for j in 0..N {
            debug_assert!(v.v[i].coeffs[j] == 0 || v.v[i].coeffs[j] == 1);
            c += v.v[i].coeffs[j] as usize;
        }
    }
    c
}

/// Branchless HighBits over a vector.
pub fn high_bits_vec_ct<P: ParameterSet, const M: usize>(v: &PolyVec<M>) -> PolyVec<M> {
    let mut out = PolyVec::<M>::zero();
    for i in 0..M {
        out.v[i] = crate::rounding::high_bits_poly_ct::<P>(&v.v[i]);
    }
    out
}

/// Branchless LowBits over a vector.
pub fn low_bits_vec_ct<P: ParameterSet, const M: usize>(v: &PolyVec<M>) -> PolyVec<M> {
    let mut out = PolyVec::<M>::zero();
    for i in 0..M {
        out.v[i] = crate::rounding::low_bits_poly_ct::<P>(&v.v[i]);
    }
    out
}

/// Branchless MakeHint over vectors.
pub fn make_hint_vec_ct<P: ParameterSet, const M: usize>(
    z: &PolyVec<M>,
    r: &PolyVec<M>,
) -> PolyVec<M> {
    let mut out = PolyVec::<M>::zero();
    for i in 0..M {
        out.v[i] = crate::rounding::make_hint_poly_ct::<P>(&z.v[i], &r.v[i]);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{MlDsa44, MlDsa65};

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
    }

    fn random_vec(rng: &mut XorShift) -> PolyVec<3> {
        let mut v = PolyVec::<3>::zero();
        for p in v.v.iter_mut() {
            for c in p.coeffs.iter_mut() {
                *c = (rng.next_u32() % Q as u32) as i32;
            }
        }
        v
    }

    /// `exceeds_bound_ct` must agree with `inf_norm >= bound` everywhere,
    /// including exactly at the boundary (the off-by-one a KAT would miss).
    #[test]
    fn exceeds_bound_ct_equals_inf_norm_check() {
        let mut rng = XorShift(0x600d_b075_1234_5678);
        for _ in 0..500 {
            let v = random_vec(&mut rng);
            let norm = inf_norm(&v);
            for bound in [1, norm - 1, norm, norm + 1, Q / 2] {
                if bound < 1 {
                    continue;
                }
                assert_eq!(exceeds_bound_ct(&v, bound), norm >= bound, "bound = {bound}");
            }
        }
        // boundary coefficient exactly at (q-1)/2 (max centred magnitude)
        let mut v = PolyVec::<3>::zero();
        v.v[1].coeffs[7] = (Q - 1) / 2;
        assert!(exceeds_bound_ct(&v, (Q - 1) / 2));
        assert!(!exceeds_bound_ct(&v, (Q - 1) / 2 + 1));
    }

    /// The branchless vector lifts must be value-equal to the baseline lifts.
    #[test]
    fn ct_vec_lifts_equal_baseline() {
        let mut rng = XorShift(0xc0de_c0de_c0de_c0de);
        for _ in 0..50 {
            let v = random_vec(&mut rng);
            let w = random_vec(&mut rng);
            assert_eq!(high_bits_vec_ct::<MlDsa65, 3>(&v).v, high_bits_vec::<MlDsa65, 3>(&v).v);
            assert_eq!(high_bits_vec_ct::<MlDsa44, 3>(&v).v, high_bits_vec::<MlDsa44, 3>(&v).v);
            assert_eq!(low_bits_vec_ct::<MlDsa65, 3>(&v).v, low_bits_vec::<MlDsa65, 3>(&v).v);
            let base = make_hint_vec::<MlDsa65, 3>(&v, &w);
            assert_eq!(make_hint_vec_ct::<MlDsa65, 3>(&v, &w).v, base.v);
            assert_eq!(count_ones_ct(&base), count_ones(&base));
        }
    }
}

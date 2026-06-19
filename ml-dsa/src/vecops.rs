//! Vector lifts of the poly-level operations, shared by keygen / sign / verify.
//! FIPS 204 applies NTT, rounding, and hints componentwise over `R^k` / `R^l`
//! (§7.4, §7.6); these helpers express that lifting once.
#![allow(clippy::needless_range_loop)]

use crate::field::{mod_pm, reduce_q};
use crate::ntt::{inv_ntt, ntt};
use crate::params::{N, Q};
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
pub fn high_bits_vec<const M: usize>(v: &PolyVec<M>) -> PolyVec<M> {
    let mut out = PolyVec::<M>::zero();
    for i in 0..M {
        out.v[i] = high_bits_poly(&v.v[i]);
    }
    out
}

/// LowBits over a vector.
pub fn low_bits_vec<const M: usize>(v: &PolyVec<M>) -> PolyVec<M> {
    let mut out = PolyVec::<M>::zero();
    for i in 0..M {
        out.v[i] = low_bits_poly(&v.v[i]);
    }
    out
}

/// MakeHint over vectors (componentwise), giving a 0/1 hint vector.
pub fn make_hint_vec<const M: usize>(z: &PolyVec<M>, r: &PolyVec<M>) -> PolyVec<M> {
    let mut out = PolyVec::<M>::zero();
    for i in 0..M {
        out.v[i] = make_hint_poly(&z.v[i], &r.v[i]);
    }
    out
}

/// UseHint over vectors (componentwise).
pub fn use_hint_vec<const M: usize>(h: &PolyVec<M>, r: &PolyVec<M>) -> PolyVec<M> {
    let mut out = PolyVec::<M>::zero();
    for i in 0..M {
        out.v[i] = use_hint_poly(&h.v[i], &r.v[i]);
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

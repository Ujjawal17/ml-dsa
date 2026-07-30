//! Shared leakage-analysis
pub mod fault;

use ml_dsa::field::montgomery_reduce;
use ml_dsa::params::MlDsa65;
use ml_dsa::rounding::{decompose, decompose_ct, make_hint, make_hint_ct};

/// Coefficients per measurement batch.
pub const BATCH: usize = 256;

/// Baseline (branchy) Decompose over a batch
#[inline(never)]
pub fn decompose_baseline_batch(coeffs: &[i32]) -> i32 {
    let mut acc = 0i32;
    for &r in coeffs {
        let (r1, r0) = decompose::<MlDsa65>(r);
        acc = acc.wrapping_add(r1).wrapping_add(r0);
    }
    acc
}

/// Branchless Decompose over a batch.
#[inline(never)]
pub fn decompose_ct_batch(coeffs: &[i32]) -> i32 {
    let mut acc = 0i32;
    for &r in coeffs {
        let (r1, r0) = decompose_ct::<MlDsa65>(r);
        acc = acc.wrapping_add(r1).wrapping_add(r0);
    }
    acc
}

/// A batch whose coefficients all trigger the baseline's q − 1 special branch.
pub fn batch_special() -> Vec<i32> {
    vec![ml_dsa::params::Q - 1; BATCH]
}

/// A batch of ordinary coefficients.
pub fn batch_common() -> Vec<i32> {
    vec![1_234_567; BATCH]
}

//MakeHint

/// Baseline (branchy) MakeHint over a batch.
#[inline(never)]
pub fn make_hint_baseline_batch(coeffs: &[i32]) -> i32 {
    let mut acc = 0i32;
    for &r in coeffs {
        acc = acc.wrapping_add(make_hint::<MlDsa65>(1, r) as i32);
    }
    acc
}

/// Branchless MakeHint over a batch.
#[inline(never)]
pub fn make_hint_ct_batch(coeffs: &[i32]) -> i32 {
    let mut acc = 0i32;
    for &r in coeffs {
        acc = acc.wrapping_add(make_hint_ct::<MlDsa65>(1, r));
    }
    acc
}

//Montgomery reduction

/// MontgomeryReduce over a batch of wide products.
#[inline(never)]
pub fn montgomery_batch(vals: &[i64]) -> i32 {
    let mut acc = 0i32;
    for &a in vals {
        acc = acc.wrapping_add(montgomery_reduce(a));
    }
    acc
}

/// A batch of "wide product" inputs for Montgomery reduction.
pub fn wide_fixed(v: i64) -> Vec<i64> {
    vec![v; BATCH]
}

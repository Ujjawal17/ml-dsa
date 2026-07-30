//! NTT / NTT^-1 in isolation: faithful (`rem_euclid`) vs improved (Montgomery + deferred reduction).

use std::hint::black_box;

use iai_callgrind::{library_benchmark, library_benchmark_group, main};
use ml_dsa::ntt::{inv_ntt, inv_ntt_fast, ntt, ntt_fast};
use ml_dsa::{Poly, PolyNTT};

/// Deterministic pseudo-random polynomial with coefficients in [0, q).
fn sample_poly() -> Poly {
    let mut p = Poly::zero();
    let mut x = 0x1234_5678_9abc_def0u64;
    for c in p.coeffs.iter_mut() {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        *c = (x % 8_380_417) as i32;
    }
    p
}

fn sample_ntt() -> PolyNTT {
    ntt(&sample_poly())
}

#[library_benchmark(setup = sample_poly)]
fn ntt_baseline(p: Poly) -> PolyNTT {
    black_box(ntt(black_box(&p)))
}

#[library_benchmark(setup = sample_poly)]
fn ntt_montgomery(p: Poly) -> PolyNTT {
    black_box(ntt_fast(black_box(&p)))
}

#[library_benchmark(setup = sample_ntt)]
fn inv_ntt_baseline(p: PolyNTT) -> Poly {
    black_box(inv_ntt(black_box(&p)))
}

#[library_benchmark(setup = sample_ntt)]
fn inv_ntt_montgomery(p: PolyNTT) -> Poly {
    black_box(inv_ntt_fast(black_box(&p)))
}

library_benchmark_group!(
    name = transform;
    benchmarks = ntt_baseline, ntt_montgomery, inv_ntt_baseline, inv_ntt_montgomery
);
main!(library_benchmark_groups = transform);

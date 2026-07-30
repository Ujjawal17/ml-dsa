//! Deterministic constant-time audit via callgrind instruction counts.
//! For each function two secret input classes are run andequal counts means constant-time, so differing counts means leak.
//!
//!   rounding   decompose   baseline vs _ct   
//!   hint       make_hint   baseline vs _ct
//!   arithmetic montgomery_reduce                
//!   sampling   sample_in_ball, rej_bounded_poly 

use std::hint::black_box;

use iai_callgrind::{library_benchmark, library_benchmark_group, main};
use ml_dsa::params::MlDsa65;
use ml_dsa::sample::{rej_bounded_poly, sample_in_ball};
use ml_dsa::Poly;
use ml_dsa_leakage::{
    batch_common, batch_special, decompose_baseline_batch, decompose_ct_batch,
    make_hint_baseline_batch, make_hint_ct_batch, montgomery_batch, wide_fixed,
};

//Decompose

#[library_benchmark(setup = batch_special)]
fn decompose_baseline_special(v: Vec<i32>) -> i32 {
    black_box(decompose_baseline_batch(black_box(&v)))
}
#[library_benchmark(setup = batch_common)]
fn decompose_baseline_common(v: Vec<i32>) -> i32 {
    black_box(decompose_baseline_batch(black_box(&v)))
}
#[library_benchmark(setup = batch_special)]
fn decompose_ct_special(v: Vec<i32>) -> i32 {
    black_box(decompose_ct_batch(black_box(&v)))
}
#[library_benchmark(setup = batch_common)]
fn decompose_ct_common(v: Vec<i32>) -> i32 {
    black_box(decompose_ct_batch(black_box(&v)))
}

library_benchmark_group!(
    name = decompose;
    benchmarks = decompose_baseline_special, decompose_baseline_common,
                 decompose_ct_special, decompose_ct_common
);

//MakeHint

#[library_benchmark(setup = batch_special)]
fn make_hint_baseline_a(v: Vec<i32>) -> i32 {
    black_box(make_hint_baseline_batch(black_box(&v)))
}
#[library_benchmark(setup = batch_common)]
fn make_hint_baseline_b(v: Vec<i32>) -> i32 {
    black_box(make_hint_baseline_batch(black_box(&v)))
}
#[library_benchmark(setup = batch_special)]
fn make_hint_ct_a(v: Vec<i32>) -> i32 {
    black_box(make_hint_ct_batch(black_box(&v)))
}
#[library_benchmark(setup = batch_common)]
fn make_hint_ct_b(v: Vec<i32>) -> i32 {
    black_box(make_hint_ct_batch(black_box(&v)))
}

library_benchmark_group!(
    name = hint;
    benchmarks = make_hint_baseline_a, make_hint_baseline_b, make_hint_ct_a, make_hint_ct_b
);

//Montgomery reduction

fn wide_a() -> Vec<i64> {
    wide_fixed(12_345)
}
fn wide_b() -> Vec<i64> {
    wide_fixed(9_876_543_210)
}

#[library_benchmark(setup = wide_a)]
fn montgomery_a(v: Vec<i64>) -> i32 {
    black_box(montgomery_batch(black_box(&v)))
}
#[library_benchmark(setup = wide_b)]
fn montgomery_b(v: Vec<i64>) -> i32 {
    black_box(montgomery_batch(black_box(&v)))
}

library_benchmark_group!(name = arithmetic; benchmarks = montgomery_a, montgomery_b);

//Rejection samplers

fn ball_seed_a() -> Vec<u8> {
    vec![7u8; 48]
}
fn ball_seed_b() -> Vec<u8> {
    vec![199u8; 48]
}
fn rej_seed_a() -> [u8; 66] {
    [3u8; 66]
}
fn rej_seed_b() -> [u8; 66] {
    [222u8; 66]
}

#[library_benchmark(setup = ball_seed_a)]
fn sample_in_ball_a(seed: Vec<u8>) -> Poly {
    black_box(sample_in_ball::<MlDsa65>(black_box(&seed)))
}
#[library_benchmark(setup = ball_seed_b)]
fn sample_in_ball_b(seed: Vec<u8>) -> Poly {
    black_box(sample_in_ball::<MlDsa65>(black_box(&seed)))
}
#[library_benchmark(setup = rej_seed_a)]
fn rej_bounded_a(seed: [u8; 66]) -> Poly {
    black_box(rej_bounded_poly::<MlDsa65>(black_box(&seed)))
}
#[library_benchmark(setup = rej_seed_b)]
fn rej_bounded_b(seed: [u8; 66]) -> Poly {
    black_box(rej_bounded_poly::<MlDsa65>(black_box(&seed)))
}

library_benchmark_group!(
    name = sampling;
    benchmarks = sample_in_ball_a, sample_in_ball_b, rej_bounded_a, rej_bounded_b
);

main!(library_benchmark_groups = decompose, hint, arithmetic, sampling);

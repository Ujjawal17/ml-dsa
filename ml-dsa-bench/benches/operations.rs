//! The three top-level operations for reference vs improved path (single call each).

use std::hint::black_box;

use iai_callgrind::{library_benchmark, library_benchmark_group, main};
use ml_dsa::ml_dsa_65;
use ml_dsa_bench::{keypair, m_prime, SEED};

//KeyGen

#[library_benchmark]
fn keygen_baseline() -> (Vec<u8>, Vec<u8>) {
    black_box(ml_dsa_65::key_gen_internal(black_box(&SEED)))
}

#[library_benchmark]
fn keygen_improved() -> (Vec<u8>, Vec<u8>) {
    black_box(ml_dsa_65::key_gen_internal_fast(black_box(&SEED)))
}

//Sign (single-shot)

fn sign_fixture() -> (Vec<u8>, Vec<u8>) {
    let (_pk, sk) = keypair();
    (sk, m_prime(b"benchmark message"))
}

#[library_benchmark(setup = sign_fixture)]
fn sign_baseline(fx: (Vec<u8>, Vec<u8>)) -> Vec<u8> {
    let (sk, mp) = fx;
    black_box(ml_dsa_65::sign_internal(black_box(&sk), black_box(&mp), black_box(&[0u8; 32])))
}

#[library_benchmark(setup = sign_fixture)]
fn sign_improved_oneshot(fx: (Vec<u8>, Vec<u8>)) -> Vec<u8> {
    let (sk, mp) = fx;
    black_box(ml_dsa_65::sign_internal_fast(black_box(&sk), black_box(&mp), black_box(&[0u8; 32])))
}

//Verify

fn verify_fixture() -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    let (pk, sk) = keypair();
    let mp = m_prime(b"benchmark message");
    let sig = ml_dsa_65::sign_internal(&sk, &mp, &[0u8; 32]);
    (pk, mp, sig)
}

#[library_benchmark(setup = verify_fixture)]
fn verify_baseline(fx: (Vec<u8>, Vec<u8>, Vec<u8>)) -> bool {
    let (pk, mp, sig) = fx;
    black_box(ml_dsa_65::verify_internal(black_box(&pk), black_box(&mp), black_box(&sig)))
}

#[library_benchmark(setup = verify_fixture)]
fn verify_improved(fx: (Vec<u8>, Vec<u8>, Vec<u8>)) -> bool {
    let (pk, mp, sig) = fx;
    black_box(ml_dsa_65::verify_internal_fast(black_box(&pk), black_box(&mp), black_box(&sig)))
}

library_benchmark_group!(
    name = operations;
    benchmarks =
        keygen_baseline, keygen_improved,
        sign_baseline, sign_improved_oneshot,
        verify_baseline, verify_improved
);
main!(library_benchmark_groups = operations);

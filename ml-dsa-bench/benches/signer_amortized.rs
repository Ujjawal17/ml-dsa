//! The centerpiece: amortized multi-signature workload (one key, N signatures).
//!
//! Baseline re-runs `skDecode + ExpandA + NTT(s1,s2,t0)` on every signature;
//! the improved path runs them once in `Signer::from_sk` and reuses `Â, ŝ1, ŝ2, t̂0`.
//! Both sign the same N messages under the same key, so they hit identical
//! rejection-loop counts — the instruction-count delta is exactly the amortization
//! (plus the Montgomery/lazy-NTT win in the shared loop body). Divide by N for the
//! per-signature figure.

use std::hint::black_box;

use iai_callgrind::{library_benchmark, library_benchmark_group, main};
use ml_dsa::ml_dsa_65::{self, Signer};
use ml_dsa_bench::{keypair, messages};

/// Number of signatures in the amortized workload.
const N: usize = 10;

fn fixture() -> (Vec<u8>, Vec<Vec<u8>>) {
    let (_pk, sk) = keypair();
    (sk, messages(N))
}

#[library_benchmark(setup = fixture)]
fn amortized_baseline(fx: (Vec<u8>, Vec<Vec<u8>>)) -> usize {
    let (sk, msgs) = fx;
    let mut acc = 0usize;
    for mp in &msgs {
        acc += black_box(ml_dsa_65::sign_internal(
            black_box(&sk),
            black_box(mp),
            black_box(&[0u8; 32]),
        ))
        .len();
    }
    black_box(acc)
}

#[library_benchmark(setup = fixture)]
fn amortized_signer(fx: (Vec<u8>, Vec<Vec<u8>>)) -> usize {
    let (sk, msgs) = fx;
    let signer = Signer::from_sk(black_box(&sk)).expect("valid secret key");
    let mut acc = 0usize;
    for mp in &msgs {
        acc += black_box(signer.sign_internal(black_box(mp), black_box(&[0u8; 32]))).len();
    }
    black_box(acc)
}

library_benchmark_group!(
    name = amortized;
    benchmarks = amortized_baseline, amortized_signer
);
main!(library_benchmark_groups = amortized);

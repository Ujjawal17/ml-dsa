//! Statistical timing leak test (dudect: Welch's t-test on real `std::time::Instant`
//! measurements). Two input classes — a **fixed** secret (Left) and a **random**
//! secret (Right) — are timed through `Decompose`; if the timing distributions
//! differ (|t| ≳ 4.5), timing depends on the secret ⇒ leak.
//!
//! **All inputs are generated up front**, then timed in a separate loop, so no RNG
//! or fill work is interleaved with the measurements — otherwise the per-class setup
//! cost (constant fill vs 256 RNG draws) leaks into the timing and produces a false
//! positive even for a constant-time target.
//!
//! Run the compiled binary directly (dudect-bencher rejects cargo's `--bench` arg):
//!   `cargo bench -p ml-dsa-leakage --no-run` then execute the `dudect_ct-*` binary.
//! Expect: `baseline` |t| vs `constant_time` |t| — the before/after leakage signal.

use std::hint::black_box;

use dudect_bencher::{ctbench_main, BenchRng, Class, CtRunner};
use ml_dsa::params::Q;
use ml_dsa_leakage::{decompose_baseline_batch, decompose_ct_batch, BATCH};
use rand::Rng;

const SAMPLES: usize = 80_000;

/// Drive `target` over pre-built inputs: Left = a fixed secret (all `q−1`, the
/// baseline's special-case branch), Right = a fresh random secret.
fn run(runner: &mut CtRunner, rng: &mut BenchRng, target: fn(&[i32]) -> i32) {
    let mut samples: Vec<(Class, Vec<i32>)> = Vec::with_capacity(SAMPLES);
    for _ in 0..SAMPLES {
        if rng.gen::<bool>() {
            samples.push((Class::Left, vec![Q - 1; BATCH]));
        } else {
            let v: Vec<i32> = (0..BATCH).map(|_| rng.gen_range(0..Q)).collect();
            samples.push((Class::Right, v));
        }
    }
    // Timing loop: only `target` runs inside run_one; no generation here.
    for (class, input) in &samples {
        runner.run_one(*class, || black_box(target(black_box(input))));
    }
}

fn baseline(runner: &mut CtRunner, rng: &mut BenchRng) {
    run(runner, rng, decompose_baseline_batch);
}

fn constant_time(runner: &mut CtRunner, rng: &mut BenchRng) {
    run(runner, rng, decompose_ct_batch);
}

ctbench_main!(baseline, constant_time);

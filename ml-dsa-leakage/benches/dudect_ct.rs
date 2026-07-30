//! Statistical timing leak test (dudect: Welch's t-test on real std::time::Instant measurements). 
//! Two input classes where a fixed secret (Left) and a random secret (Right) are timed through 
//! Decompose and if the timing distributions differ, since timing depends on the secret it causes leak.

//! To run the compiled binary:
//!   `cargo bench -p ml-dsa-leakage --no-run` then execute the `dudect_ct-*` binary.

use std::hint::black_box;

use dudect_bencher::{ctbench_main, BenchRng, Class, CtRunner};
use ml_dsa::params::Q;
use ml_dsa_leakage::{decompose_baseline_batch, decompose_ct_batch, BATCH};
use rand::Rng;

const SAMPLES: usize = 80_000;

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
    // In timing loop, only target runs inside run_one and no generation is done here.
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

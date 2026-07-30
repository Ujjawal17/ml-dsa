//! Wall-clock signing latency for ML-DSA-65
//! ONly wal clock measurement as everything else is deterministic instruction counts

//! To Run: `cargo run -p ml-dsa-bench --example latency --release`.

use std::time::Instant;

use ml_dsa::ml_dsa_65::Signer;
use ml_dsa_bench::{keypair, messages};

fn main() {
    let (_pk, sk) = keypair();
    let signer = Signer::from_sk(&sk).expect("valid secret key");
    let n = 10_000;
    let msgs = messages(n);
    let rnd = [0u8; 32];

    // Warm up (page-in, branch predictors, caches).
    for m in msgs.iter().take(100) {
        std::hint::black_box(signer.sign_internal(std::hint::black_box(m), &rnd));
    }

    let mut ns: Vec<u128> = Vec::with_capacity(n);
    for m in &msgs {
        let t = Instant::now();
        let s = signer.sign_internal(std::hint::black_box(m), &rnd);
        ns.push(t.elapsed().as_nanos());
        std::hint::black_box(&s);
    }

    ns.sort_unstable();
    let us = |x: u128| x as f64 / 1000.0;
    let median = us(ns[n / 2]);
    let mean = us(ns.iter().sum::<u128>() / n as u128);
    let p99 = us(ns[(n * 99) / 100]);
    println!("ML-DSA-65 signing latency, one key, amortized Signer, {n} signs:");
    println!("  median = {median:.1} us");
    println!("  mean   = {mean:.1} us");
    println!("  p99    = {p99:.1} us");
}

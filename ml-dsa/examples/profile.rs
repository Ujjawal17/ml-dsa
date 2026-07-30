//! Instruction-count profiling target for ML-DSA-65.
//!
//! Runs one KeyGen -> Sign -> Verify at a **fixed seed** (so the figures are
//! reproducible) for either the faithful `baseline` path (default) or the
//! `improved` path (fast/branchless arithmetic + amortized `Signer`). callgrind is
//! deterministic, so a single call per operation yields an exact profile, no
//! iteration count to tune.
//!
//! Each operation is wrapped in an #[inline(never)] bench_* function with a
//! distinctive name, so its symbol always survives release inlining and can be
//! toggled reliably. --toggle-collect collects the toggled function and every
//! function it calls, so toggling bench_sign captures all of signing while the
//! preceding key setup (which runs outside it) is excluded.
//!
//! By Isolating one operation with --collect-atstart=no and one --toggle-collect:
//!
//!   baseline (run: profile):
//!     '*bench_keygen*'   '*bench_sign*'   '*bench_verify*'
//!   improved (run: `profile improved`):
//!     '*bench_keygen_fast*'   '*bench_signer_setup*'   (one-time: ExpandA + NTTs)
//!     '*bench_signer_sign*'   (per-signature, amortized)   '*bench_verify_fast*'
//!
//! Example:
//!     cargo build --release --example profile
//!     mkdir -p results
//!     valgrind --tool=callgrind --collect-atstart=no \
//!       --toggle-collect='*bench_sign*' \
//!       --callgrind-out-file=results/cg.baseline.sign.out \
//!       target/release/examples/profile
//!     callgrind_annotate results/cg.baseline.sign.out

use std::collections::BTreeMap;
use std::hint::black_box;

use ml_dsa::ml_dsa_65::{self, Signer};

const DEFAULT_MSG: &[u8] = b"instruction-count profiling message";

/// M' for the pure variant with an empty context.
fn m_prime(msg: &[u8]) -> Vec<u8> {
    let mut mp = vec![0u8, 0u8];
    mp.extend_from_slice(msg);
    mp
}

/// Printed the deterministic rejection-loop iteration count for a batch of messages, so a low- and high-iteration case can be picked for the profiling comparison.
/// Key is fixed (seed 0x42), so only the message varies the count.
///
/// Also listed the first message found for each distinct iteration count, so the cost model total(n) can be profiled at intermediate n.
fn scan() {
    let (_pk, sk) = ml_dsa_65::key_gen_internal(&[0x42u8; 32]);
    let mut min = (u32::MAX, String::new());
    let mut max = (0u32, String::new());
    let mut first_per_n: BTreeMap<u32, String> = BTreeMap::new();
    for i in 0..200u32 {
        let msg = format!("case-{i:03}");
        let (_sig, attempts) =
            ml_dsa_65::sign_deterministic_traced(&sk, msg.as_bytes(), b"").unwrap();
        if attempts < min.0 {
            min = (attempts, msg.clone());
        }
        if attempts > max.0 {
            max = (attempts, msg.clone());
        }
        first_per_n.entry(attempts).or_insert(msg);
    }
    // Also reported the default profiling message for reference.
    let (_s, a) = ml_dsa_65::sign_deterministic_traced(&sk, DEFAULT_MSG, b"").unwrap();
    println!("default message: {a} iteration(s)");
    println!("min: {} iteration(s) at message {:?}", min.0, min.1);
    println!("max: {} iteration(s) at message {:?}", max.0, max.1);
    println!("first message per iteration count:");
    for (n, msg) in &first_per_n {
        println!("  n={n}: {msg:?}");
    }
}

//Baseline operations

#[inline(never)]
fn bench_keygen(xi: &[u8; 32]) -> (Vec<u8>, Vec<u8>) {
    ml_dsa_65::key_gen_internal(black_box(xi))
}

#[inline(never)]
fn bench_sign(sk: &[u8], mp: &[u8]) -> Vec<u8> {
    ml_dsa_65::sign_internal(black_box(sk), black_box(mp), black_box(&[0u8; 32]))
}

#[inline(never)]
fn bench_verify(pk: &[u8], mp: &[u8], sig: &[u8]) -> bool {
    ml_dsa_65::verify_internal(black_box(pk), black_box(mp), black_box(sig))
}

/// Faithful baseline: plain `mod q` arithmetic, branchy rounding, per-call ExpandA.
fn profile_baseline(msg: &[u8]) {
    let (pk, sk) = bench_keygen(&[0x42u8; 32]);
    let mp = m_prime(msg);
    let sig = bench_sign(&sk, &mp);
    black_box(bench_verify(&pk, &mp, &sig));
}

//Improved operations

#[inline(never)]
fn bench_keygen_fast(xi: &[u8; 32]) -> (Vec<u8>, Vec<u8>) {
    ml_dsa_65::key_gen_internal_fast(black_box(xi))
}

#[inline(never)]
fn bench_signer_setup(sk: &[u8]) -> Signer {
    Signer::from_sk(black_box(sk)).expect("valid secret key")
}

#[inline(never)]
fn bench_signer_sign(signer: &Signer, mp: &[u8]) -> Vec<u8> {
    signer.sign_internal(black_box(mp), black_box(&[0u8; 32]))
}

#[inline(never)]
fn bench_verify_fast(pk: &[u8], mp: &[u8], sig: &[u8]) -> bool {
    ml_dsa_65::verify_internal_fast(black_box(pk), black_box(mp), black_box(sig))
}

/// Improved path: Montgomery/lazy NTT + branchless rounding, and the amortized
/// Signer split into from_sk (one-time) and sign_internal (per-signature).
fn profile_improved(msg: &[u8]) {
    let (pk, sk) = bench_keygen_fast(&[0x42u8; 32]);
    let mp = m_prime(msg);
    let signer = bench_signer_setup(&sk);
    let sig = bench_signer_sign(&signer, &mp);
    black_box(bench_verify_fast(&pk, &mp, &sig));
}

fn main() {
    // argv[1] = mode ("baseline" | "improved" | "scan"); argv[2] = optional message.
    let args: Vec<String> = std::env::args().collect();
    let mode = args.get(1).map(String::as_str).unwrap_or("baseline");
    let msg: &[u8] = args.get(2).map(|s| s.as_bytes()).unwrap_or(DEFAULT_MSG);
    match mode {
        "scan" => scan(),
        "improved" => profile_improved(msg),
        _ => profile_baseline(msg),
    }
}

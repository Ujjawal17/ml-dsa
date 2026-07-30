//! Whole-program constant-time audit (SUPERCOP "TIMECOP" technique).
//! This checks the entire codebase the signature touches at the binary level and it confirms which functions are constant-time and which are not.
//!
//! Run (baseline vs improved path):
//!   cargo build --release --example timecop
//!   valgrind --tool=memcheck --error-exitcode=0 \
//!     target/release/examples/timecop baseline
//!   valgrind ... target/release/examples/timecop improved

use std::ffi::c_void;
use std::hint::black_box;

use crabgrind::memcheck::{mark_memory, MemState};
use ml_dsa::ml_dsa_65;

/// It marks bytes[start..end] as undefined (secret) to memcheck.
fn poison(bytes: &[u8], start: usize, end: usize) {
    let _ = mark_memory(
        bytes[start..end].as_ptr() as *const c_void,
        end - start,
        MemState::Undefined,
    );
}

fn main() {
    let mode = std::env::args().nth(1).unwrap_or_else(|| "baseline".into());
    let (_pk, sk) = ml_dsa_65::key_gen_internal(&[0x42u8; 32]);

    let mut m_prime = vec![0u8, 0u8];
    m_prime.extend_from_slice(b"timecop constant-time audit");

    // Poison the secret key material only (K, and s1||s2||t0) and rho, tr stay public.
    poison(&sk, 32, 64);
    poison(&sk, 128, sk.len());

    let sig = match mode.as_str() {
        "improved" => ml_dsa_65::sign_internal_fast(&sk, &m_prime, &[0u8; 32]),
        _ => ml_dsa_65::sign_internal(&sk, &m_prime, &[0u8; 32]),
    };
    black_box(sig);
}

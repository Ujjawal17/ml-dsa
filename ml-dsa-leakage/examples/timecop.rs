//! Whole-program constant-time audit (the ctgrind / SUPERCOP "TIMECOP" technique).
//!
//! Marks the *secret* key material as undefined to Valgrind's memcheck, then runs a
//! full `sign_internal`. memcheck propagates "undefined-ness" through the computation
//! and reports EVERY place a secret-derived value reaches:
//!   * a conditional branch  -> "Conditional jump or move depends on uninitialised
//!     value(s)"  = a control-flow (timing) leak;
//!   * a memory address       -> a secret-dependent memory access = a cache leak.
//!
//! This checks the entire codebase the signature touches at the binary level (so it
//! respects the compiler's actual lowering), instead of hand-picking functions —
//! i.e. it confirms which functions are constant-time and which are not.
//!
//! Run (baseline vs improved path):
//!   cargo build --release --example timecop
//!   valgrind --tool=memcheck --error-exitcode=0 \
//!     target/release/examples/timecop baseline
//!   valgrind ... target/release/examples/timecop improved
//!
//! Secret-key byte layout (skEncode): rho[0..32] K[32..64] tr[64..128] s1|s2|t0[128..].
//! rho and tr are public, so only [32..64] and [128..] are poisoned.

use std::ffi::c_void;
use std::hint::black_box;

use crabgrind::memcheck::{mark_memory, MemState};
use ml_dsa::ml_dsa_65;

/// Mark `bytes[start..end]` as undefined (secret) to memcheck; no-op without Valgrind.
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

    // Poison the secret key material only (K, and s1||s2||t0); rho, tr stay public.
    poison(&sk, 32, 64);
    poison(&sk, 128, sk.len());

    let sig = match mode.as_str() {
        "improved" => ml_dsa_65::sign_internal_fast(&sk, &m_prime, &[0u8; 32]),
        _ => ml_dsa_65::sign_internal(&sk, &m_prime, &[0u8; 32]),
    };
    black_box(sig);
}

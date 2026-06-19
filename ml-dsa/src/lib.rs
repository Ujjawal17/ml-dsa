#![forbid(unsafe_code)]
//! ML-DSA (FIPS 204) — a faithful, constant-time reference implementation in
//! pure safe Rust.
//!
//! Phase 0 (the "spine") fixes the shared vocabulary every later module builds
//! against: parameters, error type, core polynomial types, the SHAKE wrappers,
//! and constant-time helpers. The implementation is hardcoded to **ML-DSA-65**
//! for the first NIST KAT pass; generalization to 44/87 behind a `ParameterSet`
//! trait comes later (plan Phase 5).
//!
//! `#![forbid(unsafe_code)]` above makes "pure safe Rust" a compiler-enforced
//! property of the whole crate — a hardening claim the C reference cannot state.

pub mod encoding;
pub mod error;
pub mod expand;
pub mod field;
pub mod hash;
pub mod hint;
pub mod keygen;
pub mod ntt;
pub mod ntt_arith;
pub mod params;
pub mod poly;
pub mod rounding;
pub mod sample;
pub mod serdes;
pub mod sign;
pub mod vecops;
pub mod verify;

mod ct;

pub use error::{Error, Result};
pub use ntt::{inv_ntt, ntt};
pub use poly::{Poly, PolyMatNTT, PolyNTT, PolyVec, PolyVecNTT};

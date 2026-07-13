//! Error type for the ML-DSA API: length checks, malformed input, and the FIPS
//! "bottom" (⊥) reject cases.

use core::fmt;

/// Errors returned by encode/decode and the high-level API.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// A byte slice had the wrong length for what it was meant to encode.
    InvalidLength { expected: usize, got: usize },
    /// Input bytes were structurally malformed (e.g. an out-of-range packed coefficient).
    MalformedInput,
    /// Context string exceeded 255 bytes (FIPS 204 Algs 2/3/4/5, line 1).
    ContextTooLong,
    /// FIPS "bottom" (⊥): a decode or verify step rejected.
    Reject,
    /// Verify-after-sign failed: the freshly produced signature did not verify,
    /// indicating a fault during signing; the signature was withheld
    /// (Bruinderink–Pessl countermeasure).
    FaultDetected,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::InvalidLength { expected, got } => {
                write!(f, "invalid length: expected {expected}, got {got}")
            }
            Error::MalformedInput => write!(f, "malformed input"),
            Error::ContextTooLong => write!(f, "context string exceeds 255 bytes"),
            Error::Reject => write!(f, "rejected (bottom)"),
            Error::FaultDetected => {
                write!(f, "verify-after-sign failed: fault detected, signature withheld")
            }
        }
    }
}

impl std::error::Error for Error {}

/// Crate-wide result alias.
pub type Result<T> = core::result::Result<T, Error>;

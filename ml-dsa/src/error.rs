//! Error type for the ML-DSA API: length checks, malformed input, and the FIPS
//! (⊥) reject cases.

use core::fmt;

/// Errors returned by encode/decode and the high-level API.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    ///byte slice got the wrong length
    InvalidLength { expected: usize, got: usize },
    ///Input bytes were malformed 
    MalformedInput,
    ///Context string exceeded 255 bytes
    ContextTooLong,
    ///FIPS "bottom" (⊥): a decode or verify step rejected.
    Reject,
    ///Verify-after-sign failed: the newly generated signature did not verify, indicating a fault during signing withholding signature
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

///Crate-wide result alias.
pub type Result<T> = core::result::Result<T, Error>;

//! FIPS 204 parameters, hardcoded for **ML-DSA-65** (security category 3).
//!
//! These are the Table 1 (scheme) and Table 2 (sizes) constants. The plan is to
//! get ML-DSA-65 byte-exact against the NIST KAT first, then lift these constants
//! behind a `ParameterSet` trait (Phase 5) to also cover ML-DSA-44 and -87.

/// Polynomial degree: the ring is `Z_q[X] / (X^n + 1)`.
pub const N: usize = 256;

/// Modulus `q = 2^23 - 2^13 + 1` (an NTT-friendly prime).
pub const Q: i32 = 8_380_417;

/// Number of low-order bits dropped by Power2Round (`t = t1·2^d + t0`).
pub const D: usize = 13;

/// Number of `±1` coefficients in the challenge polynomial `c`.
pub const TAU: usize = 49;

/// Collision-strength parameter, in bits. The commitment hash `c~` is `lambda/4` bytes.
pub const LAMBDA: usize = 192;

/// Coefficient range of the masking vector `y`: `(-gamma1, gamma1]`.
pub const GAMMA1: i32 = 1 << 19; // 524288

/// Low-order rounding parameter for Decompose / HighBits / LowBits.
pub const GAMMA2: i32 = (Q - 1) / 32; // 261888

/// Rows of the public matrix `A` (A is `K x L` over `R_q`).
pub const K: usize = 6;

/// Columns of the public matrix `A`.
pub const L: usize = 5;

/// Secret coefficient bound: `s1, s2` coefficients lie in `[-eta, eta]`.
pub const ETA: i32 = 4;

/// Rejection bound `beta = tau * eta`.
pub const BETA: i32 = TAU as i32 * ETA; // 196

/// Maximum number of `1`s in the hint vector `h`.
pub const OMEGA: usize = 55;

// --- Derived byte sizes (FIPS 204 Table 2) ---

/// Seed / per-call randomness size.
pub const SEED_BYTES: usize = 32;

/// Length of the commitment hash `c~`, `= lambda/4` bytes.
pub const C_TILDE_BYTES: usize = LAMBDA / 4; // 48

/// Encoded public key length, in bytes.
pub const PK_BYTES: usize = 1952;

/// Encoded secret key length, in bytes.
pub const SK_BYTES: usize = 4032;

/// Encoded signature length, in bytes.
pub const SIG_BYTES: usize = 3309;

#[cfg(test)]
mod tests {
    use super::*;

    /// Cross-check the derived sizes against FIPS 204 Table 2 for ML-DSA-65, so a
    /// typo in a constant is caught immediately rather than as a KAT mismatch later.
    #[test]
    fn ml_dsa_65_sizes() {
        // bitlen(q-1) = 23, so t1 keeps (23 - d) = 10 bits per coefficient.
        assert_eq!(PK_BYTES, 32 + 32 * K * (23 - D));
        // c~ is lambda/4 bytes.
        assert_eq!(C_TILDE_BYTES, 48);
        // beta = tau * eta.
        assert_eq!(BETA, 196);
        assert_eq!(GAMMA2, 261_888);
    }
}

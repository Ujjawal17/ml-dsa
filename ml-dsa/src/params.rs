use crate::encoding::bitlen;

/// Polynomial degree, the ring is Z_q[X] / (X^n + 1).
pub const N: usize = 256;

/// Modulus q = 2^23 - 2^13 + 1.
pub const Q: i32 = 8_380_417;

/// Number of low-order bits dropped by Power2Round (t = t1·2^d + t0).
pub const D: usize = 13;

/// Seed / per-call randomness size.
pub const SEED_BYTES: usize = 32;

pub trait ParameterSet {
    /// ACVP name of the set.
    const NAME: &'static str;

    //Table 1
    /// Rows of the public matrix A (A is K x L over R_q).
    const K: usize;
    /// Columns of the public matrix A.
    const L: usize;
    /// Secret coefficient bound: s1, s2 coefficients lie in [-eta, eta].
    const ETA: i32;
    /// Number of ±1 coefficients in the challenge polynomial c.
    const TAU: usize;
    /// Collision-strength parameter, in bits. c~ is lambda/4 bytes.
    const LAMBDA: usize;
    /// Coefficient range of the masking vector y: (-gamma1, gamma1].
    const GAMMA1: i32;
    /// Low-order rounding parameter for Decompose / HighBits / LowBits.
    const GAMMA2: i32;
    /// Maximum number of 1s in the hint vector h.
    const OMEGA: usize;

    //Derived (Table 1)
    /// Rejection bound beta = tau * eta.
    const BETA: i32 = Self::TAU as i32 * Self::ETA;

    //Derived byte sizes (Table 2)
    /// Length of the commitment hash c~, = lambda/4 bytes.
    const C_TILDE_BYTES: usize = Self::LAMBDA / 4;
    /// Encoded public key: rho || SimpleBitPack(t1) per row, t1 at bitlen(q-1) - d = 10 bits per coefficient.
    const PK_BYTES: usize = 32 + 32 * Self::K * (bitlen((Q - 1) as u32) as usize - D);
    /// Encoded secret key: rho || K || tr, then s1, s2 at bitlen(2*eta) bits and t0 at d bits per coefficient.
    const SK_BYTES: usize = 32
        + 32
        + 64
        + 32 * ((Self::K + Self::L) * bitlen(2 * Self::ETA as u32) as usize + D * Self::K);
    /// Encoded signature: c~ || BitPack(z) || HintBitPack(h), z at 1 + bitlen(gamma1 - 1) bits per coefficient.
    const SIG_BYTES: usize = Self::LAMBDA / 4
        + 32 * Self::L * (1 + bitlen(Self::GAMMA1 as u32 - 1) as usize)
        + Self::OMEGA
        + Self::K;
}

/// ML-DSA-44 (security category 2).
pub struct MlDsa44;

impl ParameterSet for MlDsa44 {
    const NAME: &'static str = "ML-DSA-44";
    const K: usize = 4;
    const L: usize = 4;
    const ETA: i32 = 2;
    const TAU: usize = 39;
    const LAMBDA: usize = 128;
    const GAMMA1: i32 = 1 << 17; // 131072
    const GAMMA2: i32 = (Q - 1) / 88; // 95232
    const OMEGA: usize = 80;
}

/// ML-DSA-65 (security category 3).
pub struct MlDsa65;

impl ParameterSet for MlDsa65 {
    const NAME: &'static str = "ML-DSA-65";
    const K: usize = 6;
    const L: usize = 5;
    const ETA: i32 = 4;
    const TAU: usize = 49;
    const LAMBDA: usize = 192;
    const GAMMA1: i32 = 1 << 19; // 524288
    const GAMMA2: i32 = (Q - 1) / 32; // 261888
    const OMEGA: usize = 55;
}

/// ML-DSA-87 (security category 5).
pub struct MlDsa87;

impl ParameterSet for MlDsa87 {
    const NAME: &'static str = "ML-DSA-87";
    const K: usize = 8;
    const L: usize = 7;
    const ETA: i32 = 2;
    const TAU: usize = 60;
    const LAMBDA: usize = 256;
    const GAMMA1: i32 = 1 << 19; // 524288
    const GAMMA2: i32 = (Q - 1) / 32; // 261888
    const OMEGA: usize = 75;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_2_sizes_all_sets() {
        assert_eq!(MlDsa44::BETA, 78);
        assert_eq!(MlDsa44::C_TILDE_BYTES, 32);
        assert_eq!(MlDsa44::GAMMA2, 95_232);
        assert_eq!(MlDsa44::PK_BYTES, 1312);
        assert_eq!(MlDsa44::SK_BYTES, 2560);
        assert_eq!(MlDsa44::SIG_BYTES, 2420);

        assert_eq!(MlDsa65::BETA, 196);
        assert_eq!(MlDsa65::C_TILDE_BYTES, 48);
        assert_eq!(MlDsa65::GAMMA2, 261_888);
        assert_eq!(MlDsa65::PK_BYTES, 1952);
        assert_eq!(MlDsa65::SK_BYTES, 4032);
        assert_eq!(MlDsa65::SIG_BYTES, 3309);

        assert_eq!(MlDsa87::BETA, 120);
        assert_eq!(MlDsa87::C_TILDE_BYTES, 64);
        assert_eq!(MlDsa87::GAMMA2, 261_888);
        assert_eq!(MlDsa87::PK_BYTES, 2592);
        assert_eq!(MlDsa87::SK_BYTES, 4896);
        assert_eq!(MlDsa87::SIG_BYTES, 4627);
    }
}

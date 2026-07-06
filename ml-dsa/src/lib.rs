#![forbid(unsafe_code)]
//! ML-DSA (FIPS 204) — a faithful reference implementation in pure safe Rust.
//!
//! The core algorithms are generic over a [`params::ParameterSet`] (plus explicit
//! `K`/`L` const generics — a stable-Rust constraint, see `params.rs`); the
//! [`ml_dsa_44`], [`ml_dsa_65`], and [`ml_dsa_87`] modules instantiate them per
//! parameter set so callers never spell the generics out.
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

/// Instantiate the §5/§6 interfaces for one parameter set as a module of plain
/// (non-generic) functions.
macro_rules! parameter_set_api {
    ($mod_name:ident, $doc:literal, $param:ty) => {
        #[doc = $doc]
        pub mod $mod_name {
            use rand_core::{CryptoRng, RngCore};

            use crate::params::ParameterSet;
            use crate::Result;

            /// The `ParameterSet` marker type of this module.
            pub type Params = $param;
            const K: usize = <$param>::K;
            const L: usize = <$param>::L;

            /// Encoded public key length, in bytes (FIPS 204 Table 2).
            pub const PK_BYTES: usize = <$param>::PK_BYTES;
            /// Encoded secret key length, in bytes (FIPS 204 Table 2).
            pub const SK_BYTES: usize = <$param>::SK_BYTES;
            /// Encoded signature length, in bytes (FIPS 204 Table 2).
            pub const SIG_BYTES: usize = <$param>::SIG_BYTES;

            /// FIPS 204, Algorithm 1 — ML-DSA.KeyGen.
            pub fn key_gen<R: CryptoRng + RngCore>(rng: &mut R) -> (Vec<u8>, Vec<u8>) {
                crate::keygen::key_gen::<Params, K, L, R>(rng)
            }

            /// FIPS 204, Algorithm 6 — ML-DSA.KeyGen_internal (deterministic in `ξ`).
            pub fn key_gen_internal(xi: &[u8; 32]) -> (Vec<u8>, Vec<u8>) {
                crate::keygen::key_gen_internal::<Params, K, L>(xi)
            }

            /// FIPS 204, Algorithm 2 — ML-DSA.Sign (hedged).
            pub fn sign<R: CryptoRng + RngCore>(
                sk: &[u8],
                m: &[u8],
                ctx: &[u8],
                rng: &mut R,
            ) -> Result<Vec<u8>> {
                crate::sign::sign::<Params, K, L, R>(sk, m, ctx, rng)
            }

            /// Deterministic signing variant (FIPS 204 §3.4): `rnd = {0}^32`.
            pub fn sign_deterministic(sk: &[u8], m: &[u8], ctx: &[u8]) -> Result<Vec<u8>> {
                crate::sign::sign_deterministic::<Params, K, L>(sk, m, ctx)
            }

            /// Deterministic signing that also returns the rejection-loop count.
            pub fn sign_deterministic_traced(
                sk: &[u8],
                m: &[u8],
                ctx: &[u8],
            ) -> Result<(Vec<u8>, u32)> {
                crate::sign::sign_deterministic_traced::<Params, K, L>(sk, m, ctx)
            }

            /// FIPS 204, Algorithm 7 — ML-DSA.Sign_internal.
            pub fn sign_internal(sk: &[u8], m_prime: &[u8], rnd: &[u8; 32]) -> Vec<u8> {
                crate::sign::sign_internal::<Params, K, L>(sk, m_prime, rnd)
            }

            /// FIPS 204, Algorithm 3 — ML-DSA.Verify.
            pub fn verify(pk: &[u8], m: &[u8], sig: &[u8], ctx: &[u8]) -> bool {
                crate::verify::verify::<Params, K, L>(pk, m, sig, ctx)
            }

            /// FIPS 204, Algorithm 8 — ML-DSA.Verify_internal.
            pub fn verify_internal(pk: &[u8], m_prime: &[u8], sig: &[u8]) -> bool {
                crate::verify::verify_internal::<Params, K, L>(pk, m_prime, sig)
            }
        }
    };
}

parameter_set_api!(ml_dsa_44, "ML-DSA-44 (security category 2).", crate::params::MlDsa44);
parameter_set_api!(ml_dsa_65, "ML-DSA-65 (security category 3).", crate::params::MlDsa65);
parameter_set_api!(ml_dsa_87, "ML-DSA-87 (security category 5).", crate::params::MlDsa87);

#[cfg(test)]
mod tests {
    /// Full sign→verify round trip through the per-set wrappers, all three sets.
    #[test]
    fn round_trip_all_parameter_sets() {
        macro_rules! round_trip {
            ($api:ident) => {
                let (pk, sk) = crate::$api::key_gen_internal(&[0x77u8; 32]);
                assert_eq!(pk.len(), crate::$api::PK_BYTES);
                assert_eq!(sk.len(), crate::$api::SK_BYTES);
                let sig = crate::$api::sign_deterministic(&sk, b"msg", b"ctx").unwrap();
                assert_eq!(sig.len(), crate::$api::SIG_BYTES);
                assert!(crate::$api::verify(&pk, b"msg", &sig, b"ctx"));
                assert!(!crate::$api::verify(&pk, b"other", &sig, b"ctx"));
            };
        }
        round_trip!(ml_dsa_44);
        round_trip!(ml_dsa_65);
        round_trip!(ml_dsa_87);
    }
}

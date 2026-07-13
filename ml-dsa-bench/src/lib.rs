//! Shared fixtures for the ML-DSA-65 instruction-count benchmarks.
//!
//! All fixtures are deterministic (fixed seed) so the benchmarks are reproducible.
//! The reference and improved paths are byte-identical, so signing the *same*
//! messages under the *same* key makes the two hit the same rejection-loop counts —
//! the difference in instruction counts is then purely the optimization.

use ml_dsa::ml_dsa_65;

/// Fixed KeyGen seed.
pub const SEED: [u8; 32] = [0x42u8; 32];

/// Deterministic key pair `(pk, sk)`.
pub fn keypair() -> (Vec<u8>, Vec<u8>) {
    ml_dsa_65::key_gen_internal(&SEED)
}

/// `M'` for the pure variant, empty context: `IntegerToBytes(0,1)^2 || M`.
pub fn m_prime(msg: &[u8]) -> Vec<u8> {
    let mut mp = vec![0u8, 0u8];
    mp.extend_from_slice(msg);
    mp
}

/// `n` distinct formatted messages for the amortized multi-signature workload.
pub fn messages(n: usize) -> Vec<Vec<u8>> {
    (0..n).map(|i| m_prime(format!("bench-msg-{i:04}").as_bytes())).collect()
}

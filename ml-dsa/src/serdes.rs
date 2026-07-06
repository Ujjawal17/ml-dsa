//! FIPS 204 §7.2 — key and signature (de)serialization (Algorithms 22–28).
//!
//! Decoders that may see adversary-controlled bytes (`pk_decode`, `sig_decode`)
//! length-check their input first and return an error on mismatch — a security
//! requirement, not a convenience. `sk_decode` is for trusted input (the spec notes
//! it may return out-of-range values on malformed input).
//!
//! All encoded sizes are derived from the `ParameterSet` scalars; the per-set
//! totals are pinned to FIPS 204 Table 2 by the tests in `params.rs`.
#![allow(clippy::needless_range_loop)]

use crate::encoding::{bit_pack, bit_unpack, bitlen, simple_bit_pack, simple_bit_unpack};
use crate::error::{Error, Result};
use crate::hint::{hint_bit_pack, hint_bit_unpack};
use crate::params::{ParameterSet, D, Q};
use crate::poly::PolyVec;

// --- Per-polynomial encoded sizes (FIPS 204 Table 2 derivations) ---

/// Top of the `t1` coefficient range: `2^(bitlen(q-1) − d) − 1` = 1023.
const T1_TOP: u32 = (1u32 << (bitlen((Q - 1) as u32) - D as u32)) - 1;
/// Bytes per `t1` polynomial under SimpleBitPack: `32 · bitlen(T1_TOP)` = 320.
const T1_BYTES: usize = 32 * bitlen(T1_TOP) as usize;
/// `t0` is packed with `a = 2^(d-1) − 1`, `b = 2^(d-1)`.
const T0_A: u32 = (1u32 << (D - 1)) - 1; // 4095
const T0_B: u32 = 1u32 << (D - 1); // 4096
/// Bytes per `t0` polynomial: `32 · bitlen(a+b)` = 416.
const T0_BYTES: usize = 32 * bitlen(T0_A + T0_B) as usize;

/// Bytes per `s1`/`s2` polynomial under BitPack(η, η): `32 · bitlen(2η)`.
fn eta_bytes<P: ParameterSet>() -> usize {
    32 * bitlen(2 * P::ETA as u32) as usize
}

/// Bytes per `z` polynomial under BitPack(γ1−1, γ1): `32 · bitlen(2γ1−1)`.
fn z_bytes<P: ParameterSet>() -> usize {
    32 * bitlen((P::GAMMA1 - 1) as u32 + P::GAMMA1 as u32) as usize
}

/// Top of the `w1` coefficient range: `(q-1)/(2γ2) − 1` (15 or 43).
fn w1_top<P: ParameterSet>() -> u32 {
    ((Q - 1) / (2 * P::GAMMA2)) as u32 - 1
}

/// FIPS 204, Algorithm 22 — pkEncode.
pub fn pk_encode<P: ParameterSet, const K: usize>(rho: &[u8; 32], t1: &PolyVec<K>) -> Vec<u8> {
    let mut pk = Vec::with_capacity(P::PK_BYTES);
    pk.extend_from_slice(rho);
    for i in 0..K {
        pk.extend_from_slice(&simple_bit_pack(&t1.v[i], T1_TOP));
    }
    pk
}

/// FIPS 204, Algorithm 23 — pkDecode. Length-checked (untrusted input).
pub fn pk_decode<P: ParameterSet, const K: usize>(pk: &[u8]) -> Result<([u8; 32], PolyVec<K>)> {
    if pk.len() != P::PK_BYTES {
        return Err(Error::InvalidLength { expected: P::PK_BYTES, got: pk.len() });
    }
    let mut rho = [0u8; 32];
    rho.copy_from_slice(&pk[..32]);
    let mut t1 = PolyVec::<K>::zero();
    for i in 0..K {
        let start = 32 + i * T1_BYTES;
        t1.v[i] = simple_bit_unpack(&pk[start..start + T1_BYTES], T1_TOP);
    }
    Ok((rho, t1))
}

type DecodedSk<const K: usize, const L: usize> =
    ([u8; 32], [u8; 32], [u8; 64], PolyVec<L>, PolyVec<K>, PolyVec<K>);

/// FIPS 204, Algorithm 24 — skEncode.
pub fn sk_encode<P: ParameterSet, const K: usize, const L: usize>(
    rho: &[u8; 32],
    k_seed: &[u8; 32],
    tr: &[u8; 64],
    s1: &PolyVec<L>,
    s2: &PolyVec<K>,
    t0: &PolyVec<K>,
) -> Vec<u8> {
    let mut sk = Vec::with_capacity(P::SK_BYTES);
    sk.extend_from_slice(rho);
    sk.extend_from_slice(k_seed);
    sk.extend_from_slice(tr);
    for i in 0..L {
        sk.extend_from_slice(&bit_pack(&s1.v[i], P::ETA as u32, P::ETA as u32));
    }
    for i in 0..K {
        sk.extend_from_slice(&bit_pack(&s2.v[i], P::ETA as u32, P::ETA as u32));
    }
    for i in 0..K {
        sk.extend_from_slice(&bit_pack(&t0.v[i], T0_A, T0_B));
    }
    sk
}

/// FIPS 204, Algorithm 25 — skDecode. Trusted input (may be out of range if malformed).
pub fn sk_decode<P: ParameterSet, const K: usize, const L: usize>(
    sk: &[u8],
) -> Result<DecodedSk<K, L>> {
    if sk.len() != P::SK_BYTES {
        return Err(Error::InvalidLength { expected: P::SK_BYTES, got: sk.len() });
    }
    let eta_bytes = eta_bytes::<P>();
    let mut rho = [0u8; 32];
    rho.copy_from_slice(&sk[..32]);
    let mut k_seed = [0u8; 32];
    k_seed.copy_from_slice(&sk[32..64]);
    let mut tr = [0u8; 64];
    tr.copy_from_slice(&sk[64..128]);

    let mut off = 128;
    let mut s1 = PolyVec::<L>::zero();
    for i in 0..L {
        s1.v[i] = bit_unpack(&sk[off..off + eta_bytes], P::ETA as u32, P::ETA as u32);
        off += eta_bytes;
    }
    let mut s2 = PolyVec::<K>::zero();
    for i in 0..K {
        s2.v[i] = bit_unpack(&sk[off..off + eta_bytes], P::ETA as u32, P::ETA as u32);
        off += eta_bytes;
    }
    let mut t0 = PolyVec::<K>::zero();
    for i in 0..K {
        t0.v[i] = bit_unpack(&sk[off..off + T0_BYTES], T0_A, T0_B);
        off += T0_BYTES;
    }
    Ok((rho, k_seed, tr, s1, s2, t0))
}

/// FIPS 204, Algorithm 26 — sigEncode.
pub fn sig_encode<P: ParameterSet, const K: usize, const L: usize>(
    c_tilde: &[u8],
    z: &PolyVec<L>,
    h: &PolyVec<K>,
) -> Vec<u8> {
    let mut sigma = Vec::with_capacity(P::SIG_BYTES);
    sigma.extend_from_slice(c_tilde);
    for i in 0..L {
        sigma.extend_from_slice(&bit_pack(&z.v[i], (P::GAMMA1 - 1) as u32, P::GAMMA1 as u32));
    }
    sigma.extend_from_slice(&hint_bit_pack::<P, K>(h));
    sigma
}

/// FIPS 204, Algorithm 27 — sigDecode. Length-checked; returns `⊥` on a malformed hint.
pub fn sig_decode<P: ParameterSet, const K: usize, const L: usize>(
    sigma: &[u8],
) -> Result<(Vec<u8>, PolyVec<L>, PolyVec<K>)> {
    if sigma.len() != P::SIG_BYTES {
        return Err(Error::InvalidLength { expected: P::SIG_BYTES, got: sigma.len() });
    }
    let z_bytes = z_bytes::<P>();
    let c_tilde = sigma[..P::C_TILDE_BYTES].to_vec();
    let mut off = P::C_TILDE_BYTES;
    let mut z = PolyVec::<L>::zero();
    for i in 0..L {
        z.v[i] =
            bit_unpack(&sigma[off..off + z_bytes], (P::GAMMA1 - 1) as u32, P::GAMMA1 as u32);
        off += z_bytes;
    }
    let h = hint_bit_unpack::<P, K>(&sigma[off..]).ok_or(Error::Reject)?;
    Ok((c_tilde, z, h))
}

/// FIPS 204, Algorithm 28 — w1Encode.
pub fn w1_encode<P: ParameterSet, const K: usize>(w1: &PolyVec<K>) -> Vec<u8> {
    let mut out = Vec::new();
    for i in 0..K {
        out.extend_from_slice(&simple_bit_pack(&w1.v[i], w1_top::<P>()));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::{MlDsa44, MlDsa65};

    const K: usize = MlDsa65::K;
    const L: usize = MlDsa65::L;

    struct XorShift(u64);
    impl XorShift {
        fn next_u32(&mut self) -> u32 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            (x >> 32) as u32
        }
        fn range(&mut self, lo: i32, hi: i32) -> i32 {
            let span = (hi - lo + 1) as u32;
            lo + (self.next_u32() % span) as i32
        }
    }

    #[test]
    fn pk_round_trip() {
        let mut rng = XorShift(0x11);
        let rho = [0xABu8; 32];
        let mut t1 = PolyVec::<K>::zero();
        for p in t1.v.iter_mut() {
            for c in p.coeffs.iter_mut() {
                *c = rng.range(0, T1_TOP as i32);
            }
        }
        let pk = pk_encode::<MlDsa65, K>(&rho, &t1);
        assert_eq!(pk.len(), MlDsa65::PK_BYTES);
        let (rho2, t12) = pk_decode::<MlDsa65, K>(&pk).unwrap();
        assert_eq!(rho2, rho);
        for i in 0..K {
            assert_eq!(t12.v[i].coeffs, t1.v[i].coeffs);
        }
    }

    #[test]
    fn sk_round_trip() {
        let mut rng = XorShift(0x22);
        let rho = [1u8; 32];
        let k_seed = [2u8; 32];
        let tr = [3u8; 64];
        let mut s1 = PolyVec::<L>::zero();
        for p in s1.v.iter_mut() {
            for c in p.coeffs.iter_mut() {
                *c = rng.range(-MlDsa65::ETA, MlDsa65::ETA);
            }
        }
        let mut s2 = PolyVec::<K>::zero();
        for p in s2.v.iter_mut() {
            for c in p.coeffs.iter_mut() {
                *c = rng.range(-MlDsa65::ETA, MlDsa65::ETA);
            }
        }
        let mut t0 = PolyVec::<K>::zero();
        for p in t0.v.iter_mut() {
            for c in p.coeffs.iter_mut() {
                *c = rng.range(-(T0_A as i32), T0_B as i32);
            }
        }
        let sk = sk_encode::<MlDsa65, K, L>(&rho, &k_seed, &tr, &s1, &s2, &t0);
        assert_eq!(sk.len(), MlDsa65::SK_BYTES);
        let (r2, k2, tr2, s1b, s2b, t0b) = sk_decode::<MlDsa65, K, L>(&sk).unwrap();
        assert_eq!((r2, k2, tr2), (rho, k_seed, tr));
        for i in 0..L {
            assert_eq!(s1b.v[i].coeffs, s1.v[i].coeffs);
        }
        for i in 0..K {
            assert_eq!(s2b.v[i].coeffs, s2.v[i].coeffs);
            assert_eq!(t0b.v[i].coeffs, t0.v[i].coeffs);
        }
    }

    #[test]
    fn sig_round_trip() {
        let mut rng = XorShift(0x33);
        let c_tilde = vec![0x5Au8; MlDsa65::C_TILDE_BYTES];
        let mut z = PolyVec::<L>::zero();
        for p in z.v.iter_mut() {
            for c in p.coeffs.iter_mut() {
                *c = rng.range(-(MlDsa65::GAMMA1 - 1), MlDsa65::GAMMA1);
            }
        }
        let mut h = PolyVec::<K>::zero();
        h.v[0].coeffs[3] = 1;
        h.v[0].coeffs[100] = 1;
        h.v[4].coeffs[42] = 1;
        let sigma = sig_encode::<MlDsa65, K, L>(&c_tilde, &z, &h);
        assert_eq!(sigma.len(), MlDsa65::SIG_BYTES);
        let (c2, z2, h2) = sig_decode::<MlDsa65, K, L>(&sigma).unwrap();
        assert_eq!(c2, c_tilde);
        for i in 0..L {
            assert_eq!(z2.v[i].coeffs, z.v[i].coeffs);
        }
        for i in 0..K {
            assert_eq!(h2.v[i].coeffs, h.v[i].coeffs);
        }
    }

    #[test]
    fn sig_round_trip_44() {
        // ML-DSA-44 exercises the 18-bit z packing (γ1 = 2^17) and ω = 80.
        let mut rng = XorShift(0x44);
        let c_tilde = vec![0xC3u8; MlDsa44::C_TILDE_BYTES];
        let mut z = PolyVec::<{ MlDsa44::L }>::zero();
        for p in z.v.iter_mut() {
            for c in p.coeffs.iter_mut() {
                *c = rng.range(-(MlDsa44::GAMMA1 - 1), MlDsa44::GAMMA1);
            }
        }
        let mut h = PolyVec::<{ MlDsa44::K }>::zero();
        h.v[1].coeffs[7] = 1;
        h.v[3].coeffs[200] = 1;
        let sigma = sig_encode::<MlDsa44, { MlDsa44::K }, { MlDsa44::L }>(&c_tilde, &z, &h);
        assert_eq!(sigma.len(), MlDsa44::SIG_BYTES);
        let (c2, z2, h2) =
            sig_decode::<MlDsa44, { MlDsa44::K }, { MlDsa44::L }>(&sigma).unwrap();
        assert_eq!(c2, c_tilde);
        for i in 0..MlDsa44::L {
            assert_eq!(z2.v[i].coeffs, z.v[i].coeffs);
        }
        for i in 0..MlDsa44::K {
            assert_eq!(h2.v[i].coeffs, h.v[i].coeffs);
        }
    }

    #[test]
    fn decoders_reject_wrong_length() {
        assert!(matches!(
            pk_decode::<MlDsa65, K>(&[0u8; 10]),
            Err(Error::InvalidLength { .. })
        ));
        assert!(matches!(
            sk_decode::<MlDsa65, K, L>(&[0u8; 10]),
            Err(Error::InvalidLength { .. })
        ));
        assert!(matches!(
            sig_decode::<MlDsa65, K, L>(&[0u8; 10]),
            Err(Error::InvalidLength { .. })
        ));
    }
}

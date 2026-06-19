//! FIPS 204 §7.2 — hint packing (Algorithms 20, 21).
//!
//! The hint `h` is a vector of `K` polynomials with binary coefficients and at most
//! `ω` ones in total. HintBitPack stores, in its first `ω` bytes, the positions of
//! the ones (per polynomial, ascending), and in its last `K` bytes a running count.
//! HintBitUnpack reverses this and returns `⊥` (None) on any malformed encoding —
//! a security requirement (FIPS 204 §6.3): Verify must reject a bad hint.
#![allow(clippy::needless_range_loop)]

use crate::params::{K, N, OMEGA};
use crate::poly::PolyVec;

/// FIPS 204, Algorithm 20 — HintBitPack: encode binary hint `h` into `ω + K` bytes.
pub fn hint_bit_pack(h: &PolyVec<K>) -> [u8; OMEGA + K] {
    let mut y = [0u8; OMEGA + K];
    let mut index = 0usize;
    for i in 0..K {
        for j in 0..N {
            if h.v[i].coeffs[j] != 0 {
                y[index] = j as u8; // position of a nonzero coeff
                index += 1;
            }
        }
        y[OMEGA + i] = index as u8; // running count after poly i
    }
    y
}

/// FIPS 204, Algorithm 21 — HintBitUnpack: reverse HintBitPack, or `⊥` (None) if malformed.
pub fn hint_bit_unpack(y: &[u8]) -> Option<PolyVec<K>> {
    if y.len() != OMEGA + K {
        return None;
    }
    let mut h = PolyVec::<K>::zero();
    let mut index = 0usize;
    for i in 0..K {
        let end = y[OMEGA + i] as usize;
        if end < index || end > OMEGA {
            return None; // malformed: count out of order / too large
        }
        let first = index;
        while index < end {
            if index > first && y[index - 1] >= y[index] {
                return None; // positions must be strictly increasing
            }
            h.v[i].coeffs[y[index] as usize] = 1;
            index += 1;
        }
    }
    // Any leftover bytes in the first ω must be zero.
    for i in index..OMEGA {
        if y[i] != 0 {
            return None;
        }
    }
    Some(h)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hint_with_ones(positions: &[(usize, usize)]) -> PolyVec<K> {
        let mut h = PolyVec::<K>::zero();
        for &(poly, coeff) in positions {
            h.v[poly].coeffs[coeff] = 1;
        }
        h
    }

    #[test]
    fn hint_round_trip() {
        let h = hint_with_ones(&[(0, 5), (0, 10), (0, 200), (2, 1), (5, 255)]);
        let packed = hint_bit_pack(&h);
        assert_eq!(packed.len(), OMEGA + K);
        assert_eq!(hint_bit_unpack(&packed).unwrap().v, h.v);
    }

    #[test]
    fn empty_hint_round_trips() {
        let h = PolyVec::<K>::zero();
        assert_eq!(hint_bit_unpack(&hint_bit_pack(&h)).unwrap().v, h.v);
    }

    #[test]
    fn rejects_non_increasing_positions() {
        let mut y = [0u8; OMEGA + K];
        // poly 0 claims two positions but they are not strictly increasing.
        y[0] = 10;
        y[1] = 5;
        y[OMEGA] = 2; // count for poly 0
        assert!(hint_bit_unpack(&y).is_none());
    }

    #[test]
    fn rejects_count_above_omega() {
        let mut y = [0u8; OMEGA + K];
        y[OMEGA] = (OMEGA + 1) as u8; // impossible count
        assert!(hint_bit_unpack(&y).is_none());
    }

    #[test]
    fn rejects_wrong_length() {
        assert!(hint_bit_unpack(&[0u8; OMEGA + K - 1]).is_none());
    }
}

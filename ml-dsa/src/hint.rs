#![allow(clippy::needless_range_loop)]

use crate::params::{ParameterSet, N};
use crate::poly::PolyVec;

/// FIPS 204, Algorithm 20 — HintBitPack to encode/serialize binary hint h into ω + K bytes.
pub fn hint_bit_pack<P: ParameterSet, const K: usize>(h: &PolyVec<K>) -> Vec<u8> {
    let mut y = vec![0u8; P::OMEGA + K];
    let mut index = 0usize;
    for i in 0..K {
        for j in 0..N {
            if h.v[i].coeffs[j] != 0 {
                y[index] = j as u8; //position of a nonzero coeff
                index += 1;
            }
        }
        y[P::OMEGA + i] = index as u8; //running count after poly i
    }
    y 
}

/// FIPS 204, Algorithm 21 — HintBitUnpack to reverse HintBitPack, or ⊥ (None) if malformed.
pub fn hint_bit_unpack<P: ParameterSet, const K: usize>(y: &[u8]) -> Option<PolyVec<K>> {
    if y.len() != P::OMEGA + K {
        return None;
    }
    let mut h = PolyVec::<K>::zero();
    let mut index = 0usize;
    for i in 0..K {
        let end = y[P::OMEGA + i] as usize;
        if end < index || end > P::OMEGA {
            return None; //malformed: either count out of order or too large
        }
        let first = index;
        while index < end {
            if index > first && y[index - 1] >= y[index] {
                return None; //positions must be strictly increasing
            }
            h.v[i].coeffs[y[index] as usize] = 1;
            index += 1;
        }
    }
    // Any leftover bytes in the first ω must be zero.
    for i in index..P::OMEGA {
        if y[i] != 0 {
            return None;
        }
    }
    Some(h)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::MlDsa65;

    const K: usize = MlDsa65::K;
    const OMEGA: usize = MlDsa65::OMEGA;

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
        let packed = hint_bit_pack::<MlDsa65, K>(&h);
        assert_eq!(packed.len(), OMEGA + K);
        assert_eq!(hint_bit_unpack::<MlDsa65, K>(&packed).unwrap().v, h.v);
    }

    #[test]
    fn empty_hint_round_trips() {
        let h = PolyVec::<K>::zero();
        let packed = hint_bit_pack::<MlDsa65, K>(&h);
        assert_eq!(hint_bit_unpack::<MlDsa65, K>(&packed).unwrap().v, h.v);
    }

    #[test]
    fn rejects_non_increasing_positions() {
        let mut y = vec![0u8; OMEGA + K];
        // poly 0 claims two positions but they are not strictly increasing.
        y[0] = 10;
        y[1] = 5;
        y[OMEGA] = 2; // count for poly 0
        assert!(hint_bit_unpack::<MlDsa65, K>(&y).is_none());
    }

    #[test]
    fn rejects_count_above_omega() {
        let mut y = vec![0u8; OMEGA + K];
        y[OMEGA] = (OMEGA + 1) as u8; // impossible count
        assert!(hint_bit_unpack::<MlDsa65, K>(&y).is_none());
    }

    #[test]
    fn rejects_wrong_length() {
        assert!(hint_bit_unpack::<MlDsa65, K>(&[0u8; OMEGA + K - 1]).is_none());
    }
}

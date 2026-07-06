//! FIPS 204 §7.1 — data conversion (Algorithms 9–19). **Little-endian** throughout.
//!
//! Implemented literally against the pseudocode (building explicit bit vectors) so
//! it reads against the standard; the fast direct-bit-twiddling version is a Part 2
//! optimization. The pack/unpack routines run over fixed-length loops with no
//! secret-dependent branches, so they are constant-time on the secret surface.
//!
//! Indexed loops mirror "for i from 0 to 255" / "for i from 0 to α−1".
#![allow(clippy::needless_range_loop)]

use crate::params::{ParameterSet, N, Q};
use crate::poly::Poly;

/// `bitlen b` — number of bits in the binary representation of `b` (`bitlen 0 = 0`).
pub const fn bitlen(b: u32) -> u32 {
    32 - b.leading_zeros()
}

/// FIPS 204, Algorithm 9 — IntegerToBits: `x mod 2^α` as `α` bits, little-endian.
pub fn integer_to_bits(x: u64, alpha: usize) -> Vec<u8> {
    let mut y = vec![0u8; alpha];
    let mut xp = x;
    for i in 0..alpha {
        y[i] = (xp & 1) as u8;
        xp >>= 1;
    }
    y
}

/// FIPS 204, Algorithm 10 — BitsToInteger: integer value of an `α`-bit little-endian string.
pub fn bits_to_integer(y: &[u8], alpha: usize) -> u64 {
    let mut x = 0u64;
    for i in 1..=alpha {
        x = 2 * x + y[alpha - i] as u64;
    }
    x
}

/// FIPS 204, Algorithm 11 — IntegerToBytes: `x mod 256^α` as `α` bytes, little-endian.
pub fn integer_to_bytes(x: u64, alpha: usize) -> Vec<u8> {
    let mut y = vec![0u8; alpha];
    let mut xp = x;
    for i in 0..alpha {
        y[i] = (xp & 0xff) as u8;
        xp >>= 8;
    }
    y
}

/// FIPS 204, Algorithm 12 — BitsToBytes: pack a bit string into bytes, LSB-first.
pub fn bits_to_bytes(y: &[u8]) -> Vec<u8> {
    let alpha = y.len();
    let mut z = vec![0u8; alpha.div_ceil(8)];
    for i in 0..alpha {
        z[i / 8] += y[i] << (i % 8);
    }
    z
}

/// FIPS 204, Algorithm 13 — BytesToBits: unpack bytes into a bit string, LSB-first.
pub fn bytes_to_bits(z: &[u8]) -> Vec<u8> {
    let alpha = z.len();
    let mut y = vec![0u8; 8 * alpha];
    for i in 0..alpha {
        let mut zi = z[i];
        for j in 0..8 {
            y[8 * i + j] = zi & 1;
            zi >>= 1;
        }
    }
    y
}

/// FIPS 204, Algorithm 14 — CoeffFromThreeBytes: 3 bytes → an element of `{0,…,q−1}` or `⊥`.
pub fn coeff_from_three_bytes(b0: u8, b1: u8, b2: u8) -> Option<i32> {
    let b2p = (b2 & 0x7f) as i32; // clear the top bit of b2
    let z = (b2p << 16) + ((b1 as i32) << 8) + b0 as i32; // 0 ≤ z ≤ 2^23 − 1
    if z < Q {
        Some(z)
    } else {
        None // ⊥ (rejection)
    }
}

/// FIPS 204, Algorithm 15 — CoeffFromHalfByte: nibble → an element of `{−η,…,η}` or `⊥`.
/// The two branches are the spec's lines 1–2 (`η = 2`, with the `mod 5` fold) and
/// lines 3–4 (`η = 4`).
pub fn coeff_from_half_byte<P: ParameterSet>(b: u8) -> Option<i32> {
    if P::ETA == 2 && b < 15 {
        Some(2 - (b % 5) as i32)
    } else if P::ETA == 4 && b < 9 {
        Some(4 - b as i32)
    } else {
        None // ⊥ (rejection)
    }
}

/// FIPS 204, Algorithm 16 — SimpleBitPack: pack `w` (coeffs in `[0, b]`) into bytes.
pub fn simple_bit_pack(w: &Poly, b: u32) -> Vec<u8> {
    let c = bitlen(b) as usize;
    let mut bits = Vec::with_capacity(N * c);
    for i in 0..N {
        bits.extend_from_slice(&integer_to_bits(w.coeffs[i] as u64, c));
    }
    bits_to_bytes(&bits)
}

/// FIPS 204, Algorithm 17 — BitPack: pack `w` (coeffs in `[−a, b]`) into bytes.
pub fn bit_pack(w: &Poly, a: u32, b: u32) -> Vec<u8> {
    let c = bitlen(a + b) as usize;
    let mut bits = Vec::with_capacity(N * c);
    for i in 0..N {
        let v = b as i64 - w.coeffs[i] as i64; // b − w_i ∈ [0, a + b]
        bits.extend_from_slice(&integer_to_bits(v as u64, c));
    }
    bits_to_bytes(&bits)
}

/// FIPS 204, Algorithm 18 — SimpleBitUnpack: reverse SimpleBitPack.
pub fn simple_bit_unpack(v: &[u8], b: u32) -> Poly {
    let c = bitlen(b) as usize;
    let z = bytes_to_bits(v);
    let mut w = Poly::zero();
    for i in 0..N {
        w.coeffs[i] = bits_to_integer(&z[i * c..i * c + c], c) as i32;
    }
    w
}

/// FIPS 204, Algorithm 19 — BitUnpack: reverse BitPack.
pub fn bit_unpack(v: &[u8], a: u32, b: u32) -> Poly {
    let c = bitlen(a + b) as usize;
    let z = bytes_to_bits(v);
    let mut w = Poly::zero();
    for i in 0..N {
        w.coeffs[i] = b as i32 - bits_to_integer(&z[i * c..i * c + c], c) as i32;
    }
    w
}

#[cfg(test)]
mod tests {
    use super::*;

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
    }

    #[test]
    fn bitlen_known() {
        assert_eq!(bitlen(0), 0);
        assert_eq!(bitlen(1), 1);
        assert_eq!(bitlen(1023), 10);
        assert_eq!(bitlen(1024), 11);
    }

    #[test]
    fn integer_bits_round_trip() {
        for &(x, a) in &[(0u64, 8usize), (1, 8), (255, 8), (0b1011, 4), (12345, 20)] {
            let bits = integer_to_bits(x, a);
            assert_eq!(bits.len(), a);
            assert_eq!(bits_to_integer(&bits, a), x);
        }
    }

    #[test]
    fn integer_to_bytes_is_little_endian() {
        assert_eq!(integer_to_bytes(0x01_02_03, 3), vec![0x03, 0x02, 0x01]);
    }

    #[test]
    fn bits_bytes_round_trip_and_endianness() {
        // 0b0000_0001 -> bit 0 set, rest clear.
        assert_eq!(bytes_to_bits(&[1])[..8], [1, 0, 0, 0, 0, 0, 0, 0]);
        let bytes = [0xa5u8, 0x3c, 0xff, 0x00];
        let bits = bytes_to_bits(&bytes);
        assert_eq!(bits_to_bytes(&bits), bytes);
    }

    #[test]
    fn coeff_from_three_bytes_rejects_at_q() {
        // q = 0x7F_E0_01.
        assert_eq!(coeff_from_three_bytes(0x00, 0x00, 0x00), Some(0));
        assert_eq!(coeff_from_three_bytes(0x00, 0xE0, 0x7F), Some(Q - 1));
        assert_eq!(coeff_from_three_bytes(0x01, 0xE0, 0x7F), None); // == q
        // top bit of b2 is ignored.
        assert_eq!(coeff_from_three_bytes(0x00, 0x00, 0x80), Some(0));
    }

    #[test]
    fn coeff_from_half_byte_eta4() {
        use crate::params::MlDsa65;
        assert_eq!(coeff_from_half_byte::<MlDsa65>(0), Some(4));
        assert_eq!(coeff_from_half_byte::<MlDsa65>(8), Some(-4));
        assert_eq!(coeff_from_half_byte::<MlDsa65>(9), None);
        assert_eq!(coeff_from_half_byte::<MlDsa65>(15), None);
    }

    #[test]
    fn coeff_from_half_byte_eta2() {
        use crate::params::MlDsa44;
        // η = 2: 2 − (b mod 5) for b < 15; b = 15 rejects.
        assert_eq!(coeff_from_half_byte::<MlDsa44>(0), Some(2));
        assert_eq!(coeff_from_half_byte::<MlDsa44>(4), Some(-2));
        assert_eq!(coeff_from_half_byte::<MlDsa44>(5), Some(2));
        assert_eq!(coeff_from_half_byte::<MlDsa44>(14), Some(-2));
        assert_eq!(coeff_from_half_byte::<MlDsa44>(15), None);
    }

    #[test]
    fn simple_bit_pack_round_trip_and_length() {
        let mut rng = XorShift(0x2468_ace0_1357_9bdf);
        let b = 1023u32; // 10-bit coeffs (t1 range)
        let mut w = Poly::zero();
        for c in w.coeffs.iter_mut() {
            *c = (rng.next_u32() % (b + 1)) as i32;
        }
        let packed = simple_bit_pack(&w, b);
        assert_eq!(packed.len(), 32 * bitlen(b) as usize); // 320
        assert_eq!(simple_bit_unpack(&packed, b).coeffs, w.coeffs);
    }

    #[test]
    fn bit_pack_round_trip_signed() {
        let mut rng = XorShift(0x1111_2222_3333_4444);
        let (a, b) = (4u32, 4u32); // coeffs in [-4, 4] (η = 4)
        let mut w = Poly::zero();
        for c in w.coeffs.iter_mut() {
            *c = (rng.next_u32() % 9) as i32 - 4;
        }
        let packed = bit_pack(&w, a, b);
        assert_eq!(packed.len(), 32 * bitlen(a + b) as usize); // 128
        assert_eq!(bit_unpack(&packed, a, b).coeffs, w.coeffs);
    }
}

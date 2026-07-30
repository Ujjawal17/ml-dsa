use sha2::Digest as _;

/// An approved pre-hash function.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PreHash {
    Sha2_224,
    Sha2_256,
    Sha2_384,
    Sha2_512,
    Sha2_512_224,
    Sha2_512_256,
    Sha3_224,
    Sha3_256,
    Sha3_384,
    Sha3_512,
    Shake128,
    Shake256,
}

impl PreHash {
    /// Map an ACVP hashAlg string to a [PreHash], or None if unsupported.
    pub fn from_acvp(name: &str) -> Option<Self> {
        Some(match name {
            "SHA2-224" => Self::Sha2_224,
            "SHA2-256" => Self::Sha2_256,
            "SHA2-384" => Self::Sha2_384,
            "SHA2-512" => Self::Sha2_512,
            "SHA2-512/224" => Self::Sha2_512_224,
            "SHA2-512/256" => Self::Sha2_512_256,
            "SHA3-224" => Self::Sha3_224,
            "SHA3-256" => Self::Sha3_256,
            "SHA3-384" => Self::Sha3_384,
            "SHA3-512" => Self::Sha3_512,
            "SHAKE-128" => Self::Shake128,
            "SHAKE-256" => Self::Shake256,
            _ => return None,
        })
    }

    /// The 11-byte DER-encoded OID.
    pub fn oid(self) -> [u8; 11] {
        let last: u8 = match self {
            Self::Sha2_256 => 0x01,
            Self::Sha2_384 => 0x02,
            Self::Sha2_512 => 0x03,
            Self::Sha2_224 => 0x04,
            Self::Sha2_512_224 => 0x05,
            Self::Sha2_512_256 => 0x06,
            Self::Sha3_224 => 0x07,
            Self::Sha3_256 => 0x08,
            Self::Sha3_384 => 0x09,
            Self::Sha3_512 => 0x0A,
            Self::Shake128 => 0x0B,
            Self::Shake256 => 0x0C,
        };
        [0x06, 0x09, 0x60, 0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02, last]
    }

    /// PH(M), the pre-hash digest of the message.
    pub fn digest(self, m: &[u8]) -> Vec<u8> {
        use sha2::{Sha224, Sha256, Sha384, Sha512, Sha512_224, Sha512_256};
        use sha3::{Sha3_224, Sha3_256, Sha3_384, Sha3_512};
        match self {
            Self::Sha2_224 => Sha224::digest(m).to_vec(),
            Self::Sha2_256 => Sha256::digest(m).to_vec(),
            Self::Sha2_384 => Sha384::digest(m).to_vec(),
            Self::Sha2_512 => Sha512::digest(m).to_vec(),
            Self::Sha2_512_224 => Sha512_224::digest(m).to_vec(),
            Self::Sha2_512_256 => Sha512_256::digest(m).to_vec(),
            Self::Sha3_224 => Sha3_224::digest(m).to_vec(),
            Self::Sha3_256 => Sha3_256::digest(m).to_vec(),
            Self::Sha3_384 => Sha3_384::digest(m).to_vec(),
            Self::Sha3_512 => Sha3_512::digest(m).to_vec(),
            Self::Shake128 => shake(m, 32, true),
            Self::Shake256 => shake(m, 64, false),
        }
    }
}

/// SHAKE XOF squeezed to out_len bytes.
fn shake(m: &[u8], out_len: usize, is128: bool) -> Vec<u8> {
    use sha3::digest::{ExtendableOutput, Update, XofReader};
    let mut out = vec![0u8; out_len];
    if is128 {
        let mut h = sha3::Shake128::default();
        h.update(m);
        h.finalize_xof().read(&mut out);
    } else {
        let mut h = sha3::Shake256::default();
        h.update(m);
        h.finalize_xof().read(&mut out);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn oid_matches_fips204_examples() {
        // The three OIDs printed in FIPS 204 Algorithm 4.
        assert_eq!(PreHash::Sha2_256.oid()[10], 0x01);
        assert_eq!(PreHash::Sha2_512.oid()[10], 0x03);
        assert_eq!(PreHash::Shake128.oid()[10], 0x0B);
        assert_eq!(&PreHash::Sha2_256.oid()[..10], &[0x06, 0x09, 0x60, 0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02]);
    }

    #[test]
    fn digest_lengths() {
        let m = b"abc";
        assert_eq!(PreHash::Sha2_224.digest(m).len(), 28);
        assert_eq!(PreHash::Sha2_256.digest(m).len(), 32);
        assert_eq!(PreHash::Sha2_384.digest(m).len(), 48);
        assert_eq!(PreHash::Sha2_512.digest(m).len(), 64);
        assert_eq!(PreHash::Sha2_512_224.digest(m).len(), 28);
        assert_eq!(PreHash::Sha2_512_256.digest(m).len(), 32);
        assert_eq!(PreHash::Sha3_224.digest(m).len(), 28);
        assert_eq!(PreHash::Sha3_512.digest(m).len(), 64);
        assert_eq!(PreHash::Shake128.digest(m).len(), 32);
        assert_eq!(PreHash::Shake256.digest(m).len(), 64);
    }
}

//! SHAKE wrappers named after FIPS 204: `H = SHAKE256`, `G = SHAKE128`.
//!
//! The API is **incremental** (absorb, then squeeze) because rejection sampling
//! (RejNTTPoly / RejBoundedPoly / SampleInBall) consumes an amount of output that
//! is not known in advance. We wrap the `sha3` crate rather than implement Keccak;
//! the standard's properties are independent of which SHAKE is plugged in.

use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::{Shake128, Shake256};

/// `H` = SHAKE256 (FIPS 204 §3.7) — absorb side.
pub struct H {
    state: Shake256,
}

/// Squeeze side of an [`H`] instance.
pub struct HReader {
    reader: <Shake256 as ExtendableOutput>::Reader,
}

impl H {
    pub fn init() -> Self {
        Self { state: Shake256::default() }
    }

    pub fn absorb(&mut self, data: &[u8]) {
        self.state.update(data);
    }

    pub fn finalize(self) -> HReader {
        HReader { reader: self.state.finalize_xof() }
    }
}

impl HReader {
    /// Squeeze exactly `out.len()` bytes into `out`.
    pub fn squeeze(&mut self, out: &mut [u8]) {
        self.reader.read(out);
    }

    /// Squeeze `n` bytes into a freshly allocated buffer.
    pub fn squeeze_vec(&mut self, n: usize) -> Vec<u8> {
        let mut buf = vec![0u8; n];
        self.reader.read(&mut buf);
        buf
    }
}

/// `G` = SHAKE128 (FIPS 204 §3.7) — absorb side.
pub struct G {
    state: Shake128,
}

/// Squeeze side of a [`G`] instance.
pub struct GReader {
    reader: <Shake128 as ExtendableOutput>::Reader,
}

impl G {
    pub fn init() -> Self {
        Self { state: Shake128::default() }
    }

    pub fn absorb(&mut self, data: &[u8]) {
        self.state.update(data);
    }

    pub fn finalize(self) -> GReader {
        GReader { reader: self.state.finalize_xof() }
    }
}

impl GReader {
    /// Squeeze exactly `out.len()` bytes into `out`.
    pub fn squeeze(&mut self, out: &mut [u8]) {
        self.reader.read(out);
    }

    /// Squeeze `n` bytes into a freshly allocated buffer.
    pub fn squeeze_vec(&mut self, n: usize) -> Vec<u8> {
        let mut buf = vec![0u8; n];
        self.reader.read(&mut buf);
        buf
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// SHAKE256("") first 8 bytes are a well-known value — confirms wiring/endianness.
    #[test]
    fn shake256_empty_known_answer() {
        let mut r = H::init().finalize();
        let got = r.squeeze_vec(8);
        assert_eq!(got, [0x46, 0xb9, 0xdd, 0x2b, 0x0b, 0xa8, 0x8d, 0x13]);
    }

    /// Incremental absorb must equal one-shot absorb of the concatenation.
    #[test]
    fn incremental_absorb_matches_oneshot() {
        let mut a = G::init();
        a.absorb(b"hello ");
        a.absorb(b"world");
        let mut b = G::init();
        b.absorb(b"hello world");
        assert_eq!(a.finalize().squeeze_vec(32), b.finalize().squeeze_vec(32));
    }
}

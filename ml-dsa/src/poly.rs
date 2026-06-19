//! Core polynomial types over `R_q = Z_q[X]/(X^256 + 1)` and its NTT image `T_q`.
//!
//! A coefficient `Zq` is stored as a **centred `i32`**. During the NTT, values may
//! transiently exceed `[-(q-1)/2, (q-1)/2]`; they are normalized back only at
//! serialization boundaries. `Poly` (normal domain) and `PolyNTT` (NTT domain) are
//! deliberately **distinct types**, so a non-NTT polynomial cannot be
//! pointwise-multiplied by accident — the type system catches the mistake.
//!
//! This is the Phase 0 integration contract. The method surface is kept minimal on
//! purpose; later phases add operations against these frozen types.

use crate::params::N;

/// A coefficient in `Z_q`, centred representation (`i32`).
pub type Zq = i32;

/// An element of `R_q`: 256 coefficients in the normal domain.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Poly {
    pub coeffs: [Zq; N],
}

/// An element of `T_q`: 256 coefficients in the NTT domain.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct PolyNTT {
    pub coeffs: [Zq; N],
}

impl Poly {
    /// The zero polynomial.
    pub const fn zero() -> Self {
        Self { coeffs: [0; N] }
    }
}

impl Default for Poly {
    fn default() -> Self {
        Self::zero()
    }
}

impl PolyNTT {
    /// The zero polynomial (NTT domain).
    pub const fn zero() -> Self {
        Self { coeffs: [0; N] }
    }
}

impl Default for PolyNTT {
    fn default() -> Self {
        Self::zero()
    }
}

/// A length-`K` vector of polynomials (normal domain).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct PolyVec<const K: usize> {
    pub v: [Poly; K],
}

impl<const K: usize> PolyVec<K> {
    pub fn zero() -> Self {
        Self { v: [Poly::zero(); K] }
    }
}

impl<const K: usize> Default for PolyVec<K> {
    fn default() -> Self {
        Self::zero()
    }
}

/// A length-`K` vector of polynomials (NTT domain).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct PolyVecNTT<const K: usize> {
    pub v: [PolyNTT; K],
}

impl<const K: usize> PolyVecNTT<K> {
    pub fn zero() -> Self {
        Self { v: [PolyNTT::zero(); K] }
    }
}

impl<const K: usize> Default for PolyVecNTT<K> {
    fn default() -> Self {
        Self::zero()
    }
}

/// A `K x L` matrix of polynomials in the NTT domain (the expanded `A-hat`).
/// FIPS 204 keeps `A` in NTT form throughout, so there is no normal-domain matrix.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct PolyMatNTT<const K: usize, const L: usize> {
    pub rows: [PolyVecNTT<L>; K],
}

impl<const K: usize, const L: usize> PolyMatNTT<K, L> {
    pub fn zero() -> Self {
        Self { rows: [PolyVecNTT::<L>::zero(); K] }
    }
}

//! FIPS 204 §7.6 — linear algebra in the NTT domain `T_q` (faithful, plain mod q).
//!
//! Indexed loops mirror the spec pseudocode ("for i from 0 to 255"), so the
//! `needless_range_loop` lint is intentionally allowed here for readability against
//! the standard.
#![allow(clippy::needless_range_loop)]

use crate::field::reduce_q;
use crate::params::N;
use crate::poly::{PolyMatNTT, PolyNTT, PolyVecNTT};

/// FIPS 204, Algorithm 44 — AddNTT: componentwise sum in `T_q`.
pub fn add_ntt(a: &PolyNTT, b: &PolyNTT) -> PolyNTT {
    let mut c = PolyNTT::zero();
    for i in 0..N {
        c.coeffs[i] = reduce_q(a.coeffs[i] as i64 + b.coeffs[i] as i64);
    }
    c
}

/// FIPS 204, Algorithm 45 — MultiplyNTT: componentwise (pointwise) product in `T_q`.
pub fn multiply_ntt(a: &PolyNTT, b: &PolyNTT) -> PolyNTT {
    let mut c = PolyNTT::zero();
    for i in 0..N {
        c.coeffs[i] = reduce_q(a.coeffs[i] as i64 * b.coeffs[i] as i64);
    }
    c
}

/// FIPS 204, Algorithm 46 — AddVectorNTT: sum of two length-`L` vectors over `T_q`.
pub fn add_vector_ntt<const L: usize>(v: &PolyVecNTT<L>, w: &PolyVecNTT<L>) -> PolyVecNTT<L> {
    let mut u = PolyVecNTT::<L>::zero();
    for i in 0..L {
        u.v[i] = add_ntt(&v.v[i], &w.v[i]);
    }
    u
}

/// FIPS 204, Algorithm 47 — ScalarVectorNTT: scalar `c` times a length-`L` vector.
pub fn scalar_vector_ntt<const L: usize>(c: &PolyNTT, v: &PolyVecNTT<L>) -> PolyVecNTT<L> {
    let mut w = PolyVecNTT::<L>::zero();
    for i in 0..L {
        w.v[i] = multiply_ntt(c, &v.v[i]);
    }
    w
}

/// FIPS 204, Algorithm 48 — MatrixVectorNTT: `K x L` matrix times length-`L` vector.
pub fn matrix_vector_ntt<const K: usize, const L: usize>(
    m: &PolyMatNTT<K, L>,
    v: &PolyVecNTT<L>,
) -> PolyVecNTT<K> {
    let mut w = PolyVecNTT::<K>::zero();
    for i in 0..K {
        for j in 0..L {
            w.v[i] = add_ntt(&w.v[i], &multiply_ntt(&m.rows[i].v[j], &v.v[j]));
        }
    }
    w
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ntt::{inv_ntt, ntt};
    use crate::params::Q;
    use crate::poly::Poly;

    #[test]
    fn add_ntt_is_componentwise() {
        let mut a = PolyNTT::zero();
        let mut b = PolyNTT::zero();
        a.coeffs[0] = Q - 1;
        b.coeffs[0] = 5;
        assert_eq!(add_ntt(&a, &b).coeffs[0], 4); // (q-1)+5 = q+4 ≡ 4
    }

    /// MatrixVectorNTT over a 1x1 matrix must equal a single pointwise product,
    /// and inverse-transform to the schoolbook product of the two operands.
    #[test]
    fn matrix_vector_1x1_matches_pointwise() {
        let mut a = Poly::zero();
        let mut b = Poly::zero();
        a.coeffs[1] = 3;
        b.coeffs[2] = 7; // product should be 21 * X^3
        let mat = PolyMatNTT::<1, 1> { rows: [PolyVecNTT::<1> { v: [ntt(&a)] }] };
        let vec = PolyVecNTT::<1> { v: [ntt(&b)] };
        let out = matrix_vector_ntt(&mat, &vec);
        let prod = inv_ntt(&out.v[0]);
        let mut expected = Poly::zero();
        expected.coeffs[3] = 21;
        assert_eq!(prod.coeffs, expected.coeffs);
    }
}

#![allow(clippy::needless_range_loop)]

use crate::field::reduce_q;
use crate::params::N;
use crate::poly::{PolyMatNTT, PolyNTT, PolyVecNTT};

/// FIPS 204, Algorithm 44 — AddNTT to calculate the componentwise sum in T_q.
pub fn add_ntt(a: &PolyNTT, b: &PolyNTT) -> PolyNTT {
    let mut c = PolyNTT::zero();
    for i in 0..N {
        c.coeffs[i] = reduce_q(a.coeffs[i] as i64 + b.coeffs[i] as i64);
    }
    c
}

/// FIPS 204, Algorithm 45 — MultiplyNTT to calculate the componentwise (pointwise) product in T_q.
pub fn multiply_ntt(a: &PolyNTT, b: &PolyNTT) -> PolyNTT {
    let mut c = PolyNTT::zero();
    for i in 0..N {
        c.coeffs[i] = reduce_q(a.coeffs[i] as i64 * b.coeffs[i] as i64);
    }
    c   
}

/// FIPS 204, Algorithm 46 — AddVectorNTT to calculate the sum of two length-L vectors over T_q.
pub fn add_vector_ntt<const L: usize>(v: &PolyVecNTT<L>, w: &PolyVecNTT<L>) -> PolyVecNTT<L> {
    let mut u = PolyVecNTT::<L>::zero();
    for i in 0..L {
        u.v[i] = add_ntt(&v.v[i], &w.v[i]);
    }
    u
}

/// FIPS 204, Algorithm 47 — ScalarVectorNTT to calculate the product of a scalar c times a length-L vector.
pub fn scalar_vector_ntt<const L: usize>(c: &PolyNTT, v: &PolyVecNTT<L>) -> PolyVecNTT<L> {
    let mut w = PolyVecNTT::<L>::zero();
    for i in 0..L {
        w.v[i] = multiply_ntt(c, &v.v[i]);
    }
    w
}

/// FIPS 204, Algorithm 48 — MatrixVectorNTT to calculate the product of a K x L matrix times length-L vector.
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

//Improved path: division-free counterparts

/// [add_ntt] without division
pub fn add_ntt_fast(a: &PolyNTT, b: &PolyNTT) -> PolyNTT {
    use crate::ct::csubq;
    let mut c = PolyNTT::zero();
    for i in 0..N {
        c.coeffs[i] = csubq(a.coeffs[i] + b.coeffs[i]);
    }
    c
}

/// [multiply_ntt] without division
pub fn multiply_ntt_fast(a: &PolyNTT, b: &PolyNTT) -> PolyNTT {
    use crate::ct::caddq;
    use crate::field::{montgomery_reduce, R2};
    let mut c = PolyNTT::zero();
    for i in 0..N {
        let t = montgomery_reduce(a.coeffs[i] as i64 * b.coeffs[i] as i64);
        c.coeffs[i] = caddq(montgomery_reduce(t as i64 * R2));
    }
    c
}

/// [add_vector_ntt] on the fast component.
pub fn add_vector_ntt_fast<const L: usize>(
    v: &PolyVecNTT<L>,
    w: &PolyVecNTT<L>,
) -> PolyVecNTT<L> {
    let mut u = PolyVecNTT::<L>::zero();
    for i in 0..L {
        u.v[i] = add_ntt_fast(&v.v[i], &w.v[i]);
    }
    u
}

/// [scalar_vector_ntt] on the fast component.
pub fn scalar_vector_ntt_fast<const L: usize>(c: &PolyNTT, v: &PolyVecNTT<L>) -> PolyVecNTT<L> {
    let mut w = PolyVecNTT::<L>::zero();
    for i in 0..L {
        w.v[i] = multiply_ntt_fast(c, &v.v[i]);
    }
    w
}

/// [matrix_vector_ntt] on the fast components.
pub fn matrix_vector_ntt_fast<const K: usize, const L: usize>(
    m: &PolyMatNTT<K, L>,
    v: &PolyVecNTT<L>,
) -> PolyVecNTT<K> {
    let mut w = PolyVecNTT::<K>::zero();
    for i in 0..K {
        for j in 0..L {
            w.v[i] = add_ntt_fast(&w.v[i], &multiply_ntt_fast(&m.rows[i].v[j], &v.v[j]));
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

    fn random_ntt_poly(rng: &mut XorShift) -> PolyNTT {
        let mut p = PolyNTT::zero();
        for c in p.coeffs.iter_mut() {
            *c = (rng.next_u32() % Q as u32) as i32;
        }
        p
    }

    #[test]
    fn fast_arith_equals_baseline() {
        let mut rng = XorShift(0x0123_4567_89ab_cdef);
        let all_max = PolyNTT { coeffs: [Q - 1; N] };
        let mut cases: Vec<(PolyNTT, PolyNTT)> =
            vec![(all_max, all_max), (PolyNTT::zero(), all_max)];
        for _ in 0..100 {
            cases.push((random_ntt_poly(&mut rng), random_ntt_poly(&mut rng)));
        }
        for (a, b) in cases {
            assert_eq!(add_ntt_fast(&a, &b).coeffs, add_ntt(&a, &b).coeffs);
            assert_eq!(multiply_ntt_fast(&a, &b).coeffs, multiply_ntt(&a, &b).coeffs);
        }
    }

    #[test]
    fn fast_matrix_vector_equals_baseline() {
        let mut rng = XorShift(0xfeed_f00d_dead_10cc);
        for _ in 0..10 {
            let m = PolyMatNTT::<2, 3> {
                rows: [
                    PolyVecNTT::<3> {
                        v: [
                            random_ntt_poly(&mut rng),
                            random_ntt_poly(&mut rng),
                            random_ntt_poly(&mut rng),
                        ],
                    },
                    PolyVecNTT::<3> {
                        v: [
                            random_ntt_poly(&mut rng),
                            random_ntt_poly(&mut rng),
                            random_ntt_poly(&mut rng),
                        ],
                    },
                ],
            };
            let v = PolyVecNTT::<3> {
                v: [
                    random_ntt_poly(&mut rng),
                    random_ntt_poly(&mut rng),
                    random_ntt_poly(&mut rng),
                ],
            };
            let fast = matrix_vector_ntt_fast(&m, &v);
            let base = matrix_vector_ntt(&m, &v);
            for i in 0..2 {
                assert_eq!(fast.v[i].coeffs, base.v[i].coeffs);
            }
        }
    }
}

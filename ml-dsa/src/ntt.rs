//! FIPS 204 §7.5 — the Number-Theoretic Transform and its inverse.
//!
//! Faithful transcription using plain `mod q` arithmetic. The zeta convention is
//! **observable through the KAT** (the matrix `A-hat` is sampled directly in the NTT
//! domain by RejNTTPoly), so it must match the spec exactly:
//! `zetas[m] = ζ^BitRev8(m) mod q` with `ζ = 1753`. The table is computed at compile
//! time from `ζ` and `BitRev8` so it can be audited rather than trusted as a blob.
//!
//! Indexed `while`/`for` loops mirror the spec pseudocode line-for-line on purpose.

use crate::field::reduce_q;
use crate::params::{N, Q};
use crate::poly::{Poly, PolyNTT};

/// `ζ = 1753`, a 512th root of unity modulo `q` (FIPS 204 §2.5).
const ZETA: i64 = 1753;

/// `256^{-1} mod q`, applied at the end of NTT^{-1} (FIPS 204 Algorithm 42, line 21).
const F_256_INV: i64 = 8_347_681;

/// FIPS 204, Algorithm 43 — BitRev8: reverse the bits of an 8-bit integer.
pub const fn bit_rev_8(m: u8) -> u8 {
    let mut r = 0u8;
    let mut i = 0;
    while i < 8 {
        r |= ((m >> i) & 1) << (7 - i);
        i += 1;
    }
    r
}

/// `zetas[m] = ζ^BitRev8(m) mod q`, in `[0, q)` (FIPS 204 §7.5 / Appendix B).
/// Computed at compile time so the table is derivable from `ζ` and `BitRev8`.
pub const ZETAS: [i32; N] = {
    let mut z = [0i32; N];
    let mut m = 0usize;
    while m < N {
        z[m] = pow_mod(ZETA, bit_rev_8(m as u8) as u64, Q as i64) as i32;
        m += 1;
    }
    z
};

/// Compile-time modular exponentiation (used only to build `ZETAS`).
const fn pow_mod(base: i64, exp: u64, q: i64) -> i64 {
    let mut result = 1i64;
    let mut b = base % q;
    let mut e = exp;
    while e > 0 {
        if e & 1 == 1 {
            result = result * b % q;
        }
        b = b * b % q;
        e >>= 1;
    }
    result
}

/// FIPS 204, Algorithm 41 — NTT: maps `w ∈ R_q` to `ŵ ∈ T_q`.
pub fn ntt(w: &Poly) -> PolyNTT {
    let mut a = w.coeffs;
    for c in a.iter_mut() {
        *c = c.rem_euclid(Q); // line 2: copy w_j into [0, q)
    }
    let mut m = 0usize;
    let mut len = 128usize;
    while len >= 1 {
        let mut start = 0usize;
        while start < N {
            m += 1;
            let z = ZETAS[m] as i64; // line 10: z ← zetas[m]
            let mut j = start;
            while j < start + len {
                let t = reduce_q(z * a[j + len] as i64); // line 12
                a[j + len] = reduce_q(a[j] as i64 - t as i64); // line 13
                a[j] = reduce_q(a[j] as i64 + t as i64); // line 14
                j += 1;
            }
            start += 2 * len;
        }
        len /= 2;
    }
    PolyNTT { coeffs: a }
}

/// FIPS 204, Algorithm 42 — NTT^{-1}: maps `ŵ ∈ T_q` back to `w ∈ R_q`.
pub fn inv_ntt(w_hat: &PolyNTT) -> Poly {
    let mut a = w_hat.coeffs;
    for c in a.iter_mut() {
        *c = c.rem_euclid(Q);
    }
    let mut m = 256usize;
    let mut len = 1usize;
    while len < N {
        let mut start = 0usize;
        while start < N {
            m -= 1;
            let z = reduce_q(-(ZETAS[m] as i64)) as i64; // line 10: z ← -zetas[m] mod q
            let mut j = start;
            while j < start + len {
                let t = a[j]; // line 12
                a[j] = reduce_q(t as i64 + a[j + len] as i64); // line 13
                a[j + len] = reduce_q(t as i64 - a[j + len] as i64); // line 14
                a[j + len] = reduce_q(z * a[j + len] as i64); // line 15
                j += 1;
            }
            start += 2 * len;
        }
        len *= 2;
    }
    for c in a.iter_mut() {
        *c = reduce_q(F_256_INV * (*c as i64)); // lines 21-24: multiply by 256^-1
    }
    Poly { coeffs: a }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::N;

    /// Tiny deterministic PRNG so tests need no external dependency.
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
        fn coeff(&mut self) -> i32 {
            (self.next_u32() % Q as u32) as i32
        }
    }

    fn random_poly(rng: &mut XorShift) -> Poly {
        let mut p = Poly::zero();
        for c in p.coeffs.iter_mut() {
            *c = rng.coeff();
        }
        p
    }

    /// Schoolbook negacyclic multiply in R_q (X^256 = -1) — the reference the NTT
    /// path must agree with. Indexed loops mirror the textbook convolution.
    #[allow(clippy::needless_range_loop)]
    fn schoolbook(a: &Poly, b: &Poly) -> Poly {
        let mut acc = [0i64; N];
        for i in 0..N {
            for j in 0..N {
                let prod = a.coeffs[i] as i64 * b.coeffs[j] as i64;
                let k = i + j;
                if k < N {
                    acc[k] += prod;
                } else {
                    acc[k - N] -= prod; // wrap with sign flip
                }
            }
        }
        let mut out = Poly::zero();
        for i in 0..N {
            out.coeffs[i] = acc[i].rem_euclid(Q as i64) as i32;
        }
        out
    }

    #[test]
    fn ntt_inv_round_trip() {
        let mut rng = XorShift(0x1234_5678_9abc_def0);
        for _ in 0..50 {
            let p = random_poly(&mut rng);
            let back = inv_ntt(&ntt(&p));
            assert_eq!(p.coeffs, back.coeffs, "ntt∘inv_ntt must be the identity");
        }
    }

    #[test]
    fn ntt_multiply_matches_schoolbook() {
        use crate::ntt_arith::multiply_ntt;
        let mut rng = XorShift(0xdead_beef_cafe_babe);
        for _ in 0..25 {
            let a = random_poly(&mut rng);
            let b = random_poly(&mut rng);
            let via_ntt = inv_ntt(&multiply_ntt(&ntt(&a), &ntt(&b)));
            let direct = schoolbook(&a, &b);
            assert_eq!(via_ntt.coeffs, direct.coeffs, "NTT product must equal schoolbook");
        }
    }

    #[test]
    fn bit_rev_8_known_values() {
        assert_eq!(bit_rev_8(1), 128);
        assert_eq!(bit_rev_8(128), 1);
        assert_eq!(bit_rev_8(0), 0);
        assert_eq!(bit_rev_8(0b0000_0011), 0b1100_0000);
    }
}

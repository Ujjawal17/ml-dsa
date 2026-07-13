//! Differential-fault harness (Bruinderink–Pessl).
//!
//! A **mirror** of `ML-DSA.Sign_internal`, reassembled from the library's public
//! building blocks so a single-bit fault can be injected at a named intermediate
//! (`cs1 = c·s1`, after MultiplyNTT + inverse NTT) — the production crate is left
//! untouched. A fidelity test asserts the no-fault mirror is **byte-identical** to
//! the real `sign_internal`, so the mirror is a faithful stand-in, not an approximation.
//!
//! The experiment models the classic attack: a fault during signing yields a faulty
//! signature; whether that (exploitable) artifact escapes is governed by
//! **verify-after-sign** (withhold if the signature fails an internal check), reported
//! over the {deterministic, hedged} × {verify-after-sign on, off} 2×2.

use ml_dsa::expand::{expand_a, expand_mask};
use ml_dsa::hash::H;
use ml_dsa::ntt::ntt;
use ml_dsa::ntt_arith::{matrix_vector_ntt, scalar_vector_ntt};
use ml_dsa::params::ParameterSet;
use ml_dsa::sample::sample_in_ball;
use ml_dsa::serdes::{sig_encode, sk_decode, w1_encode};
use ml_dsa::vecops::{
    add_vec, center_vec, count_ones, high_bits_vec, inf_norm, inv_ntt_vec, low_bits_vec,
    make_hint_vec, neg_vec, ntt_vec, sub_vec,
};

/// A single-bit fault injected at a named signing intermediate.
#[derive(Clone, Copy, Debug)]
pub enum Fault {
    /// No fault — the mirror reproduces `sign_internal` exactly.
    None,
    /// Flip bit `bit` of coefficient `coeff` in polynomial `poly` of `cs1 = c·s1`
    /// (after MultiplyNTT + inverse NTT). Low bits keep `z` within the norm bound so
    /// the signer still accepts, while the fault propagates through verification via
    /// `A·Δ` (a full-magnitude change), so the faulty signature does not verify.
    Cs1BitFlip { poly: usize, coeff: usize, bit: u32 },
}

/// Mirror of the baseline `sign_internal`; byte-identical to it when
/// `fault == Fault::None` (asserted by the fidelity test below).
pub fn sign_internal_mirror<P: ParameterSet, const K: usize, const L: usize>(
    sk: &[u8],
    m_prime: &[u8],
    rnd: &[u8; 32],
    fault: Fault,
) -> Vec<u8> {
    let (rho, k_seed, tr, s1, s2, t0) = sk_decode::<P, K, L>(sk).expect("valid secret key");
    let s1_hat = ntt_vec(&s1);
    let s2_hat = ntt_vec(&s2);
    let t0_hat = ntt_vec(&t0);
    let a_hat = expand_a::<K, L>(&rho);

    let mut hm = H::init();
    hm.absorb(&tr);
    hm.absorb(m_prime);
    let mut mu = [0u8; 64];
    hm.finalize().squeeze(&mut mu);

    let mut hr = H::init();
    hr.absorb(&k_seed);
    hr.absorb(rnd);
    hr.absorb(&mu);
    let mut rho_pp = [0u8; 64];
    hr.finalize().squeeze(&mut rho_pp);

    let mut kappa = 0u32;
    loop {
        let y = expand_mask::<P, L>(&rho_pp, kappa);
        let w = inv_ntt_vec(&matrix_vector_ntt(&a_hat, &ntt_vec(&y)));
        let w1 = high_bits_vec::<P, K>(&w);

        let mut hc = H::init();
        hc.absorb(&mu);
        hc.absorb(&w1_encode::<P, K>(&w1));
        let c_tilde = hc.finalize().squeeze_vec(P::C_TILDE_BYTES);
        let c = sample_in_ball::<P>(&c_tilde);
        let c_hat = ntt(&c);

        let mut cs1 = inv_ntt_vec(&scalar_vector_ntt(&c_hat, &s1_hat));
        // ---- fault injection point (a plain line of code; no library hooks) ----
        if let Fault::Cs1BitFlip { poly, coeff, bit } = fault {
            cs1.v[poly].coeffs[coeff] ^= 1 << bit;
        }
        let cs2 = inv_ntt_vec(&scalar_vector_ntt(&c_hat, &s2_hat));
        let z = add_vec(&y, &cs1);
        let r0 = low_bits_vec::<P, K>(&sub_vec(&w, &cs2));

        if inf_norm(&z) >= P::GAMMA1 - P::BETA || inf_norm(&r0) >= P::GAMMA2 - P::BETA {
            kappa += L as u32;
            continue;
        }

        let ct0 = inv_ntt_vec(&scalar_vector_ntt(&c_hat, &t0_hat));
        let r_arg = add_vec(&sub_vec(&w, &cs2), &ct0);
        let h = make_hint_vec::<P, K>(&neg_vec(&ct0), &r_arg);

        if inf_norm(&ct0) >= P::GAMMA2 || count_ones(&h) > P::OMEGA {
            kappa += L as u32;
            continue;
        }

        return sig_encode::<P, K, L>(&c_tilde, &center_vec(&z), &h);
    }
}

/// Outcome of one 2×2 cell.
#[derive(Clone, Copy, Debug)]
pub struct Cell {
    /// The faulty signature differs from the unfaulted signature on the same `rnd`.
    pub faulty: bool,
    /// The faulty signature verifies under the real public key.
    pub verifies: bool,
    /// The faulty signature is released (verify-after-sign off, or on but it verified).
    pub released: bool,
    /// A faulty, exploitable artifact escaped the signer.
    pub escaped: bool,
    /// Deterministic mode: the attacker can reproduce the matching correct signature
    /// on the same `y`, giving the directly-comparable correct/faulty pair.
    pub comparable_pair: bool,
}

#[cfg(test)]
mod tests {
    use super::*;
    use ml_dsa::ml_dsa_65;
    use ml_dsa::params::MlDsa65;

    const K: usize = MlDsa65::K;
    const L: usize = MlDsa65::L;

    /// The no-fault mirror must equal the real signer byte-for-byte.
    #[test]
    fn mirror_is_byte_identical_to_real_signer() {
        for seed in 0u8..6 {
            let (_pk, sk) = ml_dsa_65::key_gen_internal(&[seed; 32]);
            for (msg, rnd) in [
                (b"".as_slice(), [0u8; 32]),
                (b"fault fidelity".as_slice(), [0x5au8; 32]),
                (b"another message".as_slice(), [0xffu8; 32]),
            ] {
                let mut m_prime = vec![0u8, 0u8];
                m_prime.extend_from_slice(msg);
                let real = ml_dsa_65::sign_internal(&sk, &m_prime, &rnd);
                let mirror =
                    sign_internal_mirror::<MlDsa65, K, L>(&sk, &m_prime, &rnd, Fault::None);
                assert_eq!(real, mirror, "seed {seed}, msg {msg:?}");
            }
        }
    }

    /// A low-bit cs1 fault keeps the signer accepting (z stays in bound) but the
    /// faulty signature must not verify (fault propagates through A·Δ).
    #[test]
    fn cs1_fault_produces_non_verifying_signature() {
        let (pk, sk) = ml_dsa_65::key_gen_internal(&[0x42u8; 32]);
        let mut m_prime = vec![0u8, 0u8];
        m_prime.extend_from_slice(b"fault effect");
        let fault = Fault::Cs1BitFlip { poly: 0, coeff: 0, bit: 3 };
        let faulty = sign_internal_mirror::<MlDsa65, K, L>(&sk, &m_prime, &[0u8; 32], fault);
        let correct = sign_internal_mirror::<MlDsa65, K, L>(&sk, &m_prime, &[0u8; 32], Fault::None);
        assert_ne!(faulty, correct, "fault must change the signature");
        assert!(!ml_dsa_65::verify_internal(&pk, &m_prime, &faulty), "faulty sig must not verify");
        assert!(ml_dsa_65::verify_internal(&pk, &m_prime, &correct), "correct sig must verify");
    }
}

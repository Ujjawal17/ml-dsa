use core::marker::PhantomData;

use rand_core::{CryptoRng, RngCore};
use zeroize::Zeroize;

use crate::error::{Error, Result};
use crate::expand::{expand_a, expand_mask};
use crate::hash::H;
use crate::ntt::ntt_fast;
use crate::ntt_arith::{matrix_vector_ntt_fast, scalar_vector_ntt_fast};
use crate::params::ParameterSet;
use crate::poly::{PolyMatNTT, PolyVecNTT};
use crate::sample::sample_in_ball;
use crate::serdes::{sig_encode, sk_decode, w1_encode};
use crate::sign::format_m_prime;
use crate::vecops::{
    add_vec, center_vec, count_ones_ct, exceeds_bound_ct, high_bits_vec_ct, inv_ntt_vec_fast,
    low_bits_vec_ct, make_hint_vec_ct, ntt_vec_fast, sub_vec,
};
use crate::verify::verify_internal_fast;

/// A prepared ML-DSA signing key, in this the key-dependent, message-independent work of Algorithm 7 (lines 1–5) done once and reused for every signature.
pub struct Signer<P: ParameterSet, const K: usize, const L: usize> {
    k_seed: [u8; 32],
    tr: [u8; 64],
    a_hat: PolyMatNTT<K, L>,
    s1_hat: PolyVecNTT<L>,
    s2_hat: PolyVecNTT<K>,
    t0_hat: PolyVecNTT<K>,
    /// Re-encoded public key and kept for verify-after-sign.
    pk: Vec<u8>,
    _params: PhantomData<P>,
}

impl<P: ParameterSet, const K: usize, const L: usize> Signer<P, K, L> {
    /// Build a prepared signer from an encoded secret key, that is skDecode + ExpandA(ρ) + NTT(s1), NTT(s2), NTT(t0), once.
    pub fn from_sk(sk: &[u8]) -> Result<Self> {
        let (rho, k_seed, tr, s1, s2, t0) = sk_decode::<P, K, L>(sk)?;
        let a_hat = expand_a::<K, L>(&rho);
        let s1_hat = ntt_vec_fast(&s1);
        let s2_hat = ntt_vec_fast(&s2);
        let t0_hat = ntt_vec_fast(&t0);

        // Reconstruct pk = pkEncode(ρ, t1) for verify-after-sign, t = A·s1 + s2 (Algorithm 6 line 5), then t1 = Power2Round(t).0.
        let as1 = matrix_vector_ntt_fast(&a_hat, &s1_hat);
        let t = add_vec(&inv_ntt_vec_fast(&as1), &s2);
        let (t1, _t0) = crate::vecops::power2round_vec(&t);
        let pk = crate::serdes::pk_encode::<P, K>(&rho, &t1);

        Ok(Self {
            k_seed,
            tr,
            a_hat,
            s1_hat,
            s2_hat,
            t0_hat,
            pk,
            _params: PhantomData,
        })
    }

    /// FIPS 204, Algorithm 7 (lines 6 onward) is byte-identical to the baseline sign_internal, with the key-dependent prefix amortized away.
    pub fn sign_internal(&self, m_prime: &[u8], rnd: &[u8; 32]) -> Vec<u8> {
        // line 6
        let mut hm = H::init();
        hm.absorb(&self.tr);
        hm.absorb(m_prime);
        let mut mu = [0u8; 64];
        hm.finalize().squeeze(&mut mu);

        // line 7
        let mut hr = H::init();
        hr.absorb(&self.k_seed);
        hr.absorb(rnd);
        hr.absorb(&mu);
        let mut rho_pp = [0u8; 64];
        hr.finalize().squeeze(&mut rho_pp);

        let mut kappa = 0u32; // line 8
        loop {
            // lines 11-13
            let y = expand_mask::<P, L>(&rho_pp, kappa);
            let w = inv_ntt_vec_fast(&matrix_vector_ntt_fast(&self.a_hat, &ntt_vec_fast(&y)));
            let w1 = high_bits_vec_ct::<P, K>(&w);

            // lines 15-17
            let mut hc = H::init();
            hc.absorb(&mu);
            hc.absorb(&w1_encode::<P, K>(&w1));
            let c_tilde = hc.finalize().squeeze_vec(P::C_TILDE_BYTES);
            let c = sample_in_ball::<P>(&c_tilde);
            let c_hat = ntt_fast(&c);

            // lines 18-21
            let cs1 = inv_ntt_vec_fast(&scalar_vector_ntt_fast(&c_hat, &self.s1_hat));
            let cs2 = inv_ntt_vec_fast(&scalar_vector_ntt_fast(&c_hat, &self.s2_hat));
            let z = add_vec(&y, &cs1);
            let r0 = low_bits_vec_ct::<P, K>(&sub_vec(&w, &cs2));

            // line 23
            if exceeds_bound_ct(&z, P::GAMMA1 - P::BETA)
                || exceeds_bound_ct(&r0, P::GAMMA2 - P::BETA)
            {
                kappa += L as u32;
                continue;
            }

            // lines 25-26
            let ct0 = inv_ntt_vec_fast(&scalar_vector_ntt_fast(&c_hat, &self.t0_hat));
            let r_arg = add_vec(&sub_vec(&w, &cs2), &ct0);
            let h = make_hint_vec_ct::<P, K>(&crate::vecops::neg_vec(&ct0), &r_arg);

            // line 28
            if exceeds_bound_ct(&ct0, P::GAMMA2) || count_ones_ct(&h) > P::OMEGA {
                kappa += L as u32;
                continue;
            }

            // line 33
            return sig_encode::<P, K, L>(&c_tilde, &center_vec(&z), &h);
        }
    }

    /// FIPS 204, Algorithm 2 — hedged signing on the prepared key.
    pub fn sign<R: CryptoRng + RngCore>(
        &self,
        m: &[u8],
        ctx: &[u8],
        rng: &mut R,
    ) -> Result<Vec<u8>> {
        if ctx.len() > 255 {
            return Err(Error::ContextTooLong);
        }
        let mut rnd = [0u8; 32];
        rng.fill_bytes(&mut rnd);
        Ok(self.sign_internal(&format_m_prime(ctx, m), &rnd))
    }

    /// Deterministic signing variant (rnd = {0}^32) on the prepared key.
    pub fn sign_deterministic(&self, m: &[u8], ctx: &[u8]) -> Result<Vec<u8>> {
        if ctx.len() > 255 {
            return Err(Error::ContextTooLong);
        }
        Ok(self.sign_internal(&format_m_prime(ctx, m), &[0u8; 32]))
    }

    /// hedged signing, then internal verification is done against this key's public key, **before** the signature is released. 
    // A verification failure means a fault occurred during signing and the (potentially exploitable) signature is withheld and an error returned instead.
    pub fn sign_verified<R: CryptoRng + RngCore>(
        &self,
        m: &[u8],
        ctx: &[u8],
        rng: &mut R,
    ) -> Result<Vec<u8>> {
        if ctx.len() > 255 {
            return Err(Error::ContextTooLong);
        }
        let mut rnd = [0u8; 32];
        rng.fill_bytes(&mut rnd);
        let m_prime = format_m_prime(ctx, m);
        let sig = self.sign_internal(&m_prime, &rnd);
        if !verify_internal_fast::<P, K, L>(&self.pk, &m_prime, &sig) {
            return Err(Error::FaultDetected);
        }
        Ok(sig)
    }

    /// Deterministic verify-after-sign (the fault-experiment 2×2 needs both).
    pub fn sign_deterministic_verified(&self, m: &[u8], ctx: &[u8]) -> Result<Vec<u8>> {
        if ctx.len() > 255 {
            return Err(Error::ContextTooLong);
        }
        let m_prime = format_m_prime(ctx, m);
        let sig = self.sign_internal(&m_prime, &[0u8; 32]);
        if !verify_internal_fast::<P, K, L>(&self.pk, &m_prime, &sig) {
            return Err(Error::FaultDetected);
        }
        Ok(sig)
    }

    /// The public key corresponding to this signing key.
    pub fn public_key(&self) -> &[u8] {
        &self.pk
    }
}

/// Zeroize the persistent secrets on drop
impl<P: ParameterSet, const K: usize, const L: usize> Drop for Signer<P, K, L> {
    fn drop(&mut self) {
        self.k_seed.zeroize();
        for p in self.s1_hat.v.iter_mut() {
            p.coeffs.zeroize();
        }
        for p in self.s2_hat.v.iter_mut() {
            p.coeffs.zeroize();
        }
        for p in self.t0_hat.v.iter_mut() {
            p.coeffs.zeroize();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::keygen::key_gen_internal;
    use crate::params::MlDsa65;

    const K: usize = MlDsa65::K;
    const L: usize = MlDsa65::L;

    #[test]
    fn signer_is_byte_identical_to_baseline_sign() {
        for seed in 0u8..8 {
            let (_pk, sk) = key_gen_internal::<MlDsa65, K, L>(&[seed; 32]);
            let signer = Signer::<MlDsa65, K, L>::from_sk(&sk).unwrap();
            for (msg, ctx, rnd) in [
                (b"message one".as_slice(), b"".as_slice(), [0u8; 32]),
                (b"message two".as_slice(), b"ctx".as_slice(), [0x5au8; 32]),
                (b"".as_slice(), b"c".as_slice(), [0xffu8; 32]),
            ] {
                let m_prime = format_m_prime(ctx, msg);
                let baseline = crate::sign::sign_internal::<MlDsa65, K, L>(&sk, &m_prime, &rnd);
                let amortized = signer.sign_internal(&m_prime, &rnd);
                assert_eq!(baseline, amortized, "seed {seed}, msg {msg:?}");
            }
        }
    }

    #[test]
    fn signer_public_key_matches_keygen() {
        let (pk, sk) = key_gen_internal::<MlDsa65, K, L>(&[0x21u8; 32]);
        let signer = Signer::<MlDsa65, K, L>::from_sk(&sk).unwrap();
        assert_eq!(signer.public_key(), &pk[..]);
    }

    #[test]
    fn sign_verified_accepts_honest_signature() {
        let (pk, sk) = key_gen_internal::<MlDsa65, K, L>(&[0x33u8; 32]);
        let signer = Signer::<MlDsa65, K, L>::from_sk(&sk).unwrap();
        let sig = signer.sign_deterministic_verified(b"msg", b"").unwrap();
        assert!(crate::verify::verify::<MlDsa65, K, L>(&pk, b"msg", &sig, b""));
    }

    #[test]
    fn from_sk_rejects_bad_length() {
        assert!(Signer::<MlDsa65, K, L>::from_sk(&[0u8; 7]).is_err());
    }
}

//! Fuzz `sig_decode` (untrusted input): must never panic on adversarial bytes,
//! must reject a malformed hint (⊥), and a successful decode must round-trip exactly
//! (a canonical signature — no non-canonical hint or out-of-range `z` slips through).
#![no_main]

use libfuzzer_sys::fuzz_target;
use ml_dsa::params::{MlDsa65, ParameterSet};
use ml_dsa::serdes::{sig_decode, sig_encode};

const K: usize = MlDsa65::K;
const L: usize = MlDsa65::L;

fuzz_target!(|data: &[u8]| {
    if let Ok((c_tilde, z, h)) = sig_decode::<MlDsa65, K, L>(data) {
        let re = sig_encode::<MlDsa65, K, L>(&c_tilde, &z, &h);
        assert_eq!(re, data, "sig decode→encode must round-trip");
    }
});

//! Fuzz `pk_decode` (untrusted input): must never panic on adversarial bytes, and a
//! successful decode must round-trip exactly (the decoder is a true inverse of the
//! encoder, so no non-canonical public key is silently accepted).
#![no_main]

use libfuzzer_sys::fuzz_target;
use ml_dsa::params::{MlDsa65, ParameterSet};
use ml_dsa::serdes::{pk_decode, pk_encode};

const K: usize = MlDsa65::K;

fuzz_target!(|data: &[u8]| {
    if let Ok((rho, t1)) = pk_decode::<MlDsa65, K>(data) {
        let re = pk_encode::<MlDsa65, K>(&rho, &t1);
        assert_eq!(re, data, "pk decode→encode must round-trip");
    }
});

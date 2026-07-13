//! Fuzz `sk_decode`: the secret key is *trusted* input (FIPS 204 notes it may yield
//! out-of-range values on malformed input), so the checked property is panic-freedom
//! only — no round-trip guarantee. The length check must still hold (no OOB reads).
#![no_main]

use libfuzzer_sys::fuzz_target;
use ml_dsa::params::{MlDsa65, ParameterSet};
use ml_dsa::serdes::sk_decode;

const K: usize = MlDsa65::K;
const L: usize = MlDsa65::L;

fuzz_target!(|data: &[u8]| {
    let _ = sk_decode::<MlDsa65, K, L>(data);
});

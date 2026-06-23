//! Annotated KeyGen -> Sign -> Verify for ML-DSA-65, with separate sub-commands.
//!
//! Full annotated demo (default message, or your own):
//!     cargo run --example roundtrip
//!     cargo run --example roundtrip -- "your message here"
//!     cargo run --example roundtrip -- demo "your message" "your context"
//!
//! Sign as a standalone step (saves the public key + signature to files):
//!     cargo run --example roundtrip -- sign "your message" ["your context"]
//!
//! Verify as a standalone step (reads those files back):
//!     cargo run --example roundtrip -- verify "your message" ["your context"]
//!
//! Signing here is deterministic (rnd = {0}^32) so results are reproducible. The
//! "hedged" variant `sign` would draw fresh randomness and give a different (still
//! valid) signature each call.

use std::fs;

use ml_dsa::keygen::key_gen_internal;
use ml_dsa::params::{PK_BYTES, SIG_BYTES, SK_BYTES};
use ml_dsa::sign::sign_deterministic_traced;
use ml_dsa::verify::verify;

const PK_FILE: &str = "roundtrip_pk.bin";
const SIG_FILE: &str = "roundtrip_sig.bin";
const DEFAULT_MSG: &str = "Hello, post-quantum world!";

/// First `n` bytes of a buffer as hex, plus the total length.
fn hex_prefix(bytes: &[u8], n: usize) -> String {
    let head: String = bytes.iter().take(n).map(|b| format!("{b:02x}")).collect();
    format!("{head}... ({} bytes total)", bytes.len())
}

/// Reproducible demo key pair (fixed seed). Real use draws a random 32-byte seed.
fn demo_keypair() -> (Vec<u8>, Vec<u8>) {
    key_gen_internal(&[0x42u8; 32])
}

/// `sign` sub-command: keygen, sign, and save the public key + signature to files
/// so `verify` can be run separately afterwards.
fn do_sign(message: &str, context: &str) {
    let (pk, sk) = demo_keypair();
    let (sig, attempts) = sign_deterministic_traced(&sk, message.as_bytes(), context.as_bytes())
        .expect("context must be <= 255 bytes");
    fs::write(PK_FILE, &pk).expect("write public key file");
    fs::write(SIG_FILE, &sig).expect("write signature file");
    println!("Signed message {message:?} (context {context:?})");
    println!("  rejection-loop iterations: {attempts}");
    println!("  public key -> {PK_FILE}   {}", hex_prefix(&pk, 8));
    println!("  signature  -> {SIG_FILE}  {}", hex_prefix(&sig, 8));
    println!("\nNow verify it:");
    println!("  cargo run --example roundtrip -- verify {message:?} {context:?}");
}

/// `verify` sub-command: read the saved public key + signature and check them
/// against the given message/context.
fn do_verify(message: &str, context: &str) {
    let pk = fs::read(PK_FILE).expect("run the `sign` sub-command first to create the files");
    let sig = fs::read(SIG_FILE).expect("run the `sign` sub-command first to create the files");
    let ok = verify(&pk, message.as_bytes(), &sig, context.as_bytes());
    println!("Verify message {message:?} (context {context:?}) against {SIG_FILE}");
    println!("  result -> {ok}");
    if ok {
        println!("  VALID: the signature matches this message + context.");
    } else {
        println!("  INVALID: wrong message/context, or the signature was changed.");
    }
}

/// `demo` sub-command (default): the full annotated round trip in one run.
fn do_demo(message: &str, context: &str) {
    println!("ML-DSA-65 round trip");
    println!("  message : {message:?}");
    println!("  context : {context:?}");

    // 1. KeyGen (FIPS 204 Algorithm 6): 32-byte seed -> (pk, sk)
    let (pk, sk) = demo_keypair();
    println!("\n1. KeyGen");
    println!("   public key : {}", hex_prefix(&pk, 8));
    println!("   secret key : {}", hex_prefix(&sk, 8));
    assert_eq!(pk.len(), PK_BYTES);
    assert_eq!(sk.len(), SK_BYTES);

    // 2. Sign (FIPS 204 Algorithm 2/7)
    let (sig, attempts) = sign_deterministic_traced(&sk, message.as_bytes(), context.as_bytes())
        .expect("context must be <= 255 bytes");
    println!("\n2. Sign");
    println!("   rejection-loop iterations : {attempts}  (re-sampled until ||z||,||r0|| small enough)");
    println!("   signature  : {}", hex_prefix(&sig, 8));
    assert_eq!(sig.len(), SIG_BYTES);

    // 3. Verify (FIPS 204 Algorithm 3/8): the honest case must accept
    println!("\n3. Verify");
    let ok = verify(&pk, message.as_bytes(), &sig, context.as_bytes());
    println!("   correct message + context -> {ok:<5}  (expected: true)");

    // 4. Negative checks: any change must be rejected
    let mut tampered = sig.clone();
    tampered[0] ^= 1; // flip a single bit
    let bad = verify(&pk, message.as_bytes(), &tampered, context.as_bytes());
    println!("   tampered signature        -> {bad:<5}  (expected: false)");
    let wrong = verify(&pk, b"a different message", &sig, context.as_bytes());
    println!("   wrong message             -> {wrong:<5}  (expected: false)");

    println!("\nTry: cargo run --example roundtrip -- sign \"my own message\"");
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let arg = |i: usize, default: &'static str| -> String {
        args.get(i).cloned().unwrap_or_else(|| default.to_string())
    };
    match args.get(1).map(String::as_str) {
        None => do_demo(DEFAULT_MSG, ""),
        Some("demo") => do_demo(&arg(2, DEFAULT_MSG), &arg(3, "")),
        Some("sign") => do_sign(&arg(2, DEFAULT_MSG), &arg(3, "")),
        Some("verify") => do_verify(&arg(2, DEFAULT_MSG), &arg(3, "")),
        // Anything else is treated as the message for the full demo.
        Some(other) => do_demo(other, &arg(2, "")),
    }
}

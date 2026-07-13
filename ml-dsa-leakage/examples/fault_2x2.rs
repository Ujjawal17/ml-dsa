//! Differential-fault 2×2: {deterministic, hedged} × {verify-after-sign on, off}.
//!
//! Injects a single-bit fault at `cs1` during signing and reports, per cell, whether
//! an exploitable faulty artifact escapes the signer. Success criterion is *escape of
//! an exploitable correct/faulty artifact*, not key recovery (cited: Bruinderink–Pessl).

use ml_dsa::ml_dsa_65;
use ml_dsa::params::{MlDsa65, ParameterSet};
use ml_dsa_leakage::fault::{sign_internal_mirror, Cell, Fault};

const K: usize = MlDsa65::K;
const L: usize = MlDsa65::L;

/// Deterministic 32-byte RNG (fixed for `deterministic`, seed-varied for `hedged`).
fn rnd_from(seed: u64) -> [u8; 32] {
    let mut out = [0u8; 32];
    let mut x = seed | 1;
    for b in out.iter_mut() {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        *b = x as u8;
    }
    out
}

fn run_cell(
    pk: &[u8],
    sk: &[u8],
    m_prime: &[u8],
    fault: Fault,
    deterministic: bool,
    verify_after_sign: bool,
) -> Cell {
    // Deterministic: rnd = 0 (reproducible, so the correct/faulty pair shares y).
    // Hedged: fresh randomness for the faulty signing call.
    let rnd = if deterministic { [0u8; 32] } else { rnd_from(0xC0FFEE) };

    let faulty = sign_internal_mirror::<MlDsa65, K, L>(sk, m_prime, &rnd, fault);
    let correct = sign_internal_mirror::<MlDsa65, K, L>(sk, m_prime, &rnd, Fault::None);

    let is_faulty = faulty != correct;
    let verifies = ml_dsa_65::verify_internal(pk, m_prime, &faulty);
    // Verify-after-sign withholds a signature that fails the internal check.
    let released = if verify_after_sign { verifies } else { true };
    let escaped = released && is_faulty;

    Cell {
        faulty: is_faulty,
        verifies,
        released,
        escaped,
        comparable_pair: deterministic,
    }
}

fn verdict(c: &Cell) -> String {
    if !c.escaped {
        "BLOCKED — faulty sig withheld by verify-after-sign".to_string()
    } else if c.comparable_pair {
        "ESCAPE — faulty sig released; correct/faulty pair on same y (exploitable)".to_string()
    } else {
        "ESCAPE — faulty sig released (hedged: pair on differing y)".to_string()
    }
}

fn main() {
    let (pk, sk) = ml_dsa_65::key_gen_internal(&[0x42u8; 32]);
    let mut m_prime = vec![0u8, 0u8];
    m_prime.extend_from_slice(b"differential fault experiment");
    // Low-bit cs1 fault: keeps z within the norm bound (signer accepts) while breaking
    // verification (fault propagates through A·Δ).
    let fault = Fault::Cs1BitFlip { poly: 0, coeff: 0, bit: 3 };

    println!("Single-bit fault at cs1[0][0] bit 3 (ML-DSA-65)\n");
    println!("{:<16}{:<52}verify-after-sign ON", "", "verify-after-sign OFF");
    for &det in &[true, false] {
        let off = run_cell(&pk, &sk, &m_prime, fault, det, false);
        let on = run_cell(&pk, &sk, &m_prime, fault, det, true);
        let label = if det { "deterministic" } else { "hedged" };
        println!("{:<16}{:<52}{}", label, verdict(&off), verdict(&on));
    }

    println!("\nPer-cell detail (faulty / verifies / released / escaped):");
    for &det in &[true, false] {
        for &vas in &[false, true] {
            let c = run_cell(&pk, &sk, &m_prime, fault, det, vas);
            println!(
                "  {:<14} v-a-s={:<5}  faulty={} verifies={} released={} escaped={}",
                if det { "deterministic" } else { "hedged" },
                vas,
                c.faulty,
                c.verifies,
                c.released,
                c.escaped
            );
        }
    }
}

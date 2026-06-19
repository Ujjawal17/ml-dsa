//! NIST ACVP known-answer-test runner for ML-DSA.
//!
//! Reads the official ACVP `internalProjection.json` vectors and checks our
//! implementation byte-for-byte. Covers ML-DSA-65 keyGen and sigGen (pure,
//! non-externalMu); sigVer is added when Verify lands.

use std::path::Path;

use ml_dsa::keygen::key_gen_internal;
use ml_dsa::sign::sign_internal;
use ml_dsa::verify::verify_internal;

fn vectors(name: &str) -> serde_json::Value {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("vectors").join(name);
    let data = std::fs::read_to_string(&path).unwrap_or_else(|_| panic!("read {name}"));
    serde_json::from_str(&data).unwrap_or_else(|_| panic!("parse {name}"))
}

fn hex_to_bytes(s: &str) -> Vec<u8> {
    (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16).expect("valid hex"))
        .collect()
}

fn hex_field(test: &serde_json::Value, key: &str) -> Vec<u8> {
    hex_to_bytes(test[key].as_str().unwrap_or(""))
}

/// `IntegerToBytes(0,1) || IntegerToBytes(|ctx|,1) || ctx || M` (pure variant).
fn format_m_prime(ctx: &[u8], m: &[u8]) -> Vec<u8> {
    let mut mp = vec![0u8, ctx.len() as u8];
    mp.extend_from_slice(ctx);
    mp.extend_from_slice(m);
    mp
}

fn run_keygen_kat() -> (usize, usize) {
    let json = vectors("keygen.json");
    let (mut pass, mut fail) = (0, 0);
    for group in json["testGroups"].as_array().unwrap() {
        if group["parameterSet"].as_str() != Some("ML-DSA-65") {
            continue;
        }
        for test in group["tests"].as_array().unwrap() {
            let mut xi = [0u8; 32];
            xi.copy_from_slice(&hex_field(test, "seed"));
            let (pk, sk) = key_gen_internal(&xi);
            if pk == hex_field(test, "pk") && sk == hex_field(test, "sk") {
                pass += 1;
            } else {
                fail += 1;
                eprintln!("  keyGen tcId {} FAILED", test["tcId"]);
            }
        }
    }
    (pass, fail)
}

fn run_siggen_kat() -> (usize, usize, usize) {
    let json = vectors("siggen.json");
    let (mut pass, mut fail, mut skipped) = (0, 0, 0);
    for group in json["testGroups"].as_array().unwrap() {
        if group["parameterSet"].as_str() != Some("ML-DSA-65") {
            continue;
        }
        let pure = group["preHash"].as_str() == Some("pure");
        let external_mu = group["externalMu"].as_bool().unwrap_or(false);
        let internal = group["signatureInterface"].as_str() == Some("internal");
        let deterministic = group["deterministic"].as_bool().unwrap_or(false);
        let tests = group["tests"].as_array().unwrap();
        if !pure || external_mu {
            skipped += tests.len(); // HashML-DSA / externalMu are out of scope here
            continue;
        }
        for test in tests {
            if test["signature"].as_str().is_none() {
                skipped += 1;
                continue;
            }
            let sk = hex_field(test, "sk");
            let msg = hex_field(test, "message");
            let rnd_vec = if deterministic { vec![0u8; 32] } else { hex_field(test, "rnd") };
            let mut rnd = [0u8; 32];
            rnd.copy_from_slice(&rnd_vec);
            let m_prime = if internal {
                msg
            } else {
                format_m_prime(&hex_field(test, "context"), &msg)
            };
            let sig = sign_internal(&sk, &m_prime, &rnd);
            if sig == hex_field(test, "signature") {
                pass += 1;
            } else {
                fail += 1;
                eprintln!("  sigGen tcId {} FAILED", test["tcId"]);
            }
        }
    }
    (pass, fail, skipped)
}

fn run_sigver_kat() -> (usize, usize, usize) {
    let json = vectors("sigver.json");
    let (mut pass, mut fail, mut skipped) = (0, 0, 0);
    for group in json["testGroups"].as_array().unwrap() {
        if group["parameterSet"].as_str() != Some("ML-DSA-65") {
            continue;
        }
        let pure = group["preHash"].as_str() == Some("pure");
        let external_mu = group["externalMu"].as_bool().unwrap_or(false);
        let internal = group["signatureInterface"].as_str() == Some("internal");
        let tests = group["tests"].as_array().unwrap();
        if !pure || external_mu {
            skipped += tests.len();
            continue;
        }
        for test in tests {
            let pk = hex_field(test, "pk");
            let sig = hex_field(test, "signature");
            let msg = hex_field(test, "message");
            let expected = test["testPassed"].as_bool().unwrap();
            let m_prime = if internal {
                msg
            } else {
                format_m_prime(&hex_field(test, "context"), &msg)
            };
            if verify_internal(&pk, &m_prime, &sig) == expected {
                pass += 1;
            } else {
                fail += 1;
                eprintln!("  sigVer tcId {} FAILED (expected {expected})", test["tcId"]);
            }
        }
    }
    (pass, fail, skipped)
}

fn main() {
    let (kp, kf) = run_keygen_kat();
    println!("ML-DSA-65 keyGen ACVP KAT: {kp} passed, {kf} failed");
    let (sp, sf, ss) = run_siggen_kat();
    println!("ML-DSA-65 sigGen ACVP KAT: {sp} passed, {sf} failed, {ss} skipped");
    let (vp, vf, vs) = run_sigver_kat();
    println!("ML-DSA-65 sigVer ACVP KAT: {vp} passed, {vf} failed, {vs} skipped");
    if kf > 0 || sf > 0 || vf > 0 {
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ml_dsa_65_keygen_kat() {
        let (pass, fail) = run_keygen_kat();
        assert!(pass > 0, "no keyGen tests found");
        assert_eq!(fail, 0, "{fail} keyGen failures");
    }

    #[test]
    fn ml_dsa_65_siggen_kat() {
        let (pass, fail, _skipped) = run_siggen_kat();
        assert!(pass > 0, "no sigGen tests found");
        assert_eq!(fail, 0, "{fail} sigGen failures");
    }

    #[test]
    fn ml_dsa_65_sigver_kat() {
        let (pass, fail, _skipped) = run_sigver_kat();
        assert!(pass > 0, "no sigVer tests found");
        assert_eq!(fail, 0, "{fail} sigVer failures");
    }
}

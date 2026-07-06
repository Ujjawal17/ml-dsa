//! NIST ACVP known-answer-test runner for ML-DSA.
//!
//! Reads the official ACVP `internalProjection.json` vectors and checks our
//! implementation byte-for-byte. Covers keyGen, sigGen, and sigVer (pure,
//! non-externalMu) for all three parameter sets: ML-DSA-44, -65, and -87.

use std::path::Path;

type KeyGenFn = fn(&[u8; 32]) -> (Vec<u8>, Vec<u8>);
type SignFn = fn(&[u8], &[u8], &[u8; 32]) -> Vec<u8>;
type VerifyFn = fn(&[u8], &[u8], &[u8]) -> bool;

/// The per-set §5/§6 entry points, so one KAT loop can drive all three sets.
struct Api {
    name: &'static str,
    key_gen_internal: KeyGenFn,
    sign_internal: SignFn,
    verify_internal: VerifyFn,
}

const APIS: [Api; 3] = [
    Api {
        name: "ML-DSA-44",
        key_gen_internal: ml_dsa::ml_dsa_44::key_gen_internal,
        sign_internal: ml_dsa::ml_dsa_44::sign_internal,
        verify_internal: ml_dsa::ml_dsa_44::verify_internal,
    },
    Api {
        name: "ML-DSA-65",
        key_gen_internal: ml_dsa::ml_dsa_65::key_gen_internal,
        sign_internal: ml_dsa::ml_dsa_65::sign_internal,
        verify_internal: ml_dsa::ml_dsa_65::verify_internal,
    },
    Api {
        name: "ML-DSA-87",
        key_gen_internal: ml_dsa::ml_dsa_87::key_gen_internal,
        sign_internal: ml_dsa::ml_dsa_87::sign_internal,
        verify_internal: ml_dsa::ml_dsa_87::verify_internal,
    },
];

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

/// Per-parameter-set pass/fail/skip tally.
#[derive(Default, Clone, Copy)]
struct Tally {
    pass: usize,
    fail: usize,
    skipped: usize,
}

fn run_keygen_kat(api: &Api) -> Tally {
    let json = vectors("keygen.json");
    let mut t = Tally::default();
    for group in json["testGroups"].as_array().unwrap() {
        if group["parameterSet"].as_str() != Some(api.name) {
            continue;
        }
        for test in group["tests"].as_array().unwrap() {
            let mut xi = [0u8; 32];
            xi.copy_from_slice(&hex_field(test, "seed"));
            let (pk, sk) = (api.key_gen_internal)(&xi);
            if pk == hex_field(test, "pk") && sk == hex_field(test, "sk") {
                t.pass += 1;
            } else {
                t.fail += 1;
                eprintln!("  {} keyGen tcId {} FAILED", api.name, test["tcId"]);
            }
        }
    }
    t
}

fn run_siggen_kat(api: &Api) -> Tally {
    let json = vectors("siggen.json");
    let mut t = Tally::default();
    for group in json["testGroups"].as_array().unwrap() {
        if group["parameterSet"].as_str() != Some(api.name) {
            continue;
        }
        let pure = group["preHash"].as_str() == Some("pure");
        let external_mu = group["externalMu"].as_bool().unwrap_or(false);
        let internal = group["signatureInterface"].as_str() == Some("internal");
        let deterministic = group["deterministic"].as_bool().unwrap_or(false);
        let tests = group["tests"].as_array().unwrap();
        if !pure || external_mu {
            t.skipped += tests.len(); // HashML-DSA / externalMu are out of scope here
            continue;
        }
        for test in tests {
            if test["signature"].as_str().is_none() {
                t.skipped += 1;
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
            let sig = (api.sign_internal)(&sk, &m_prime, &rnd);
            if sig == hex_field(test, "signature") {
                t.pass += 1;
            } else {
                t.fail += 1;
                eprintln!("  {} sigGen tcId {} FAILED", api.name, test["tcId"]);
            }
        }
    }
    t
}

fn run_sigver_kat(api: &Api) -> Tally {
    let json = vectors("sigver.json");
    let mut t = Tally::default();
    for group in json["testGroups"].as_array().unwrap() {
        if group["parameterSet"].as_str() != Some(api.name) {
            continue;
        }
        let pure = group["preHash"].as_str() == Some("pure");
        let external_mu = group["externalMu"].as_bool().unwrap_or(false);
        let internal = group["signatureInterface"].as_str() == Some("internal");
        let tests = group["tests"].as_array().unwrap();
        if !pure || external_mu {
            t.skipped += tests.len();
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
            if (api.verify_internal)(&pk, &m_prime, &sig) == expected {
                t.pass += 1;
            } else {
                t.fail += 1;
                eprintln!("  {} sigVer tcId {} FAILED (expected {expected})", api.name, test["tcId"]);
            }
        }
    }
    t
}

fn main() {
    let mut failures = 0;
    for api in &APIS {
        let k = run_keygen_kat(api);
        println!("{} keyGen ACVP KAT: {} passed, {} failed", api.name, k.pass, k.fail);
        let s = run_siggen_kat(api);
        println!(
            "{} sigGen ACVP KAT: {} passed, {} failed, {} skipped",
            api.name, s.pass, s.fail, s.skipped
        );
        let v = run_sigver_kat(api);
        println!(
            "{} sigVer ACVP KAT: {} passed, {} failed, {} skipped",
            api.name, v.pass, v.fail, v.skipped
        );
        failures += k.fail + s.fail + v.fail;
    }
    if failures > 0 {
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keygen_kat_all_sets() {
        for api in &APIS {
            let t = run_keygen_kat(api);
            assert!(t.pass > 0, "no {} keyGen tests found", api.name);
            assert_eq!(t.fail, 0, "{} {} keyGen failures", t.fail, api.name);
        }
    }

    #[test]
    fn siggen_kat_all_sets() {
        for api in &APIS {
            let t = run_siggen_kat(api);
            assert!(t.pass > 0, "no {} sigGen tests found", api.name);
            assert_eq!(t.fail, 0, "{} {} sigGen failures", t.fail, api.name);
        }
    }

    #[test]
    fn sigver_kat_all_sets() {
        for api in &APIS {
            let t = run_sigver_kat(api);
            assert!(t.pass > 0, "no {} sigVer tests found", api.name);
            assert_eq!(t.fail, 0, "{} {} sigVer failures", t.fail, api.name);
        }
    }
}

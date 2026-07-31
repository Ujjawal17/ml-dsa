//! NIST ACVP known-answer-test runner for ML-DSA.
//! It reads the official ACVP `internalProjection.json` vectors and checks the implementation byte-for-byte.

use std::path::Path;

type KeyGenFn = fn(&[u8; 32]) -> (Vec<u8>, Vec<u8>);
type SignFn = fn(&[u8], &[u8], &[u8; 32]) -> Vec<u8>;
type VerifyFn = fn(&[u8], &[u8], &[u8]) -> bool;

/// It maps entry points of a parameter set to run both reference and improved implementations under the same KAT loop as its primary correctness evidence.
struct Api {
    /// ACVP parameterSet string this row consumes.
    param_set: &'static str,
    /// Display label (to distinguishes reference from improved).
    label: &'static str,
    key_gen_internal: KeyGenFn,
    sign_internal: SignFn,
    verify_internal: VerifyFn,
}

macro_rules! api_rows {
    ($name:literal, $module:ident) => {
        [
            Api {
                param_set: $name,
                label: concat!($name, " [reference]"),
                key_gen_internal: ml_dsa::$module::key_gen_internal,
                sign_internal: ml_dsa::$module::sign_internal,
                verify_internal: ml_dsa::$module::verify_internal,
            },
            Api {
                param_set: $name,
                label: concat!($name, " [improved] "),
                key_gen_internal: ml_dsa::$module::key_gen_internal_fast,
                sign_internal: ml_dsa::$module::sign_internal_fast,
                verify_internal: ml_dsa::$module::verify_internal_fast,
            },
        ]
    };
}

const APIS: [[Api; 2]; 3] = [
    api_rows!("ML-DSA-44", ml_dsa_44),
    api_rows!("ML-DSA-65", ml_dsa_65),
    api_rows!("ML-DSA-87", ml_dsa_87),
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

fn format_m_prime(ctx: &[u8], m: &[u8]) -> Vec<u8> {
    let mut mp = vec![0u8, ctx.len() as u8];
    mp.extend_from_slice(ctx);
    mp.extend_from_slice(m);
    mp
}

/// Per-parameter-set pass/fail/skip checking.
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
        if group["parameterSet"].as_str() != Some(api.param_set) {
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
                eprintln!("  {} keyGen tcId {} FAILED", api.label, test["tcId"]);
            }
        }
    }
    t
}

fn run_siggen_kat(api: &Api) -> Tally {
    let json = vectors("siggen.json");
    let mut t = Tally::default();
    for group in json["testGroups"].as_array().unwrap() {
        if group["parameterSet"].as_str() != Some(api.param_set) {
            continue;
        }
        let prehash = group["preHash"].as_str() == Some("preHash");
        let external_mu = group["externalMu"].as_bool().unwrap_or(false);
        let internal = group["signatureInterface"].as_str() == Some("internal");
        let deterministic = group["deterministic"].as_bool().unwrap_or(false);
        let tests = group["tests"].as_array().unwrap();
        if prehash || external_mu {
            // HashML-DSA is covered by the dedicated pre-hash KATs below;
            // externalMu is out of scope. Both the pure (external) and the
            // plain internal-interface groups are exercised here.
            t.skipped += tests.len();
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
                eprintln!("  {} sigGen tcId {} FAILED", api.label, test["tcId"]);
            }
        }
    }
    t
}

fn run_sigver_kat(api: &Api) -> Tally {
    let json = vectors("sigver.json");
    let mut t = Tally::default();
    for group in json["testGroups"].as_array().unwrap() {
        if group["parameterSet"].as_str() != Some(api.param_set) {
            continue;
        }
        let prehash = group["preHash"].as_str() == Some("preHash");
        let external_mu = group["externalMu"].as_bool().unwrap_or(false);
        let internal = group["signatureInterface"].as_str() == Some("internal");
        let tests = group["tests"].as_array().unwrap();
        if prehash || external_mu {
            // HashML-DSA covered by the dedicated pre-hash KATs; externalMu out of scope.
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
                eprintln!("  {} sigVer tcId {} FAILED (expected {expected})", api.label, test["tcId"]);
            }
        }
    }
    t
}

const PARAM_SETS: [&str; 3] = ["ML-DSA-44", "ML-DSA-65", "ML-DSA-87"];

/// HashML-DSA (pre-hash, Algorithm 4) sigGen KAT for one parameter set.
fn run_prehash_siggen(param_set: &str) -> Tally {
    let json = vectors("siggen.json");
    let mut t = Tally::default();
    for group in json["testGroups"].as_array().unwrap() {
        if group["parameterSet"].as_str() != Some(param_set) {
            continue;
        }
        if group["preHash"].as_str() != Some("preHash")
            || group["externalMu"].as_bool().unwrap_or(false)
        {
            continue;
        }
        let deterministic = group["deterministic"].as_bool().unwrap_or(false);
        for test in group["tests"].as_array().unwrap() {
            if test["signature"].as_str().is_none() {
                t.skipped += 1;
                continue;
            }
            let ph = match ml_dsa::PreHash::from_acvp(test["hashAlg"].as_str().unwrap_or("")) {
                Some(p) => p,
                None => {
                    t.skipped += 1;
                    continue;
                }
            };
            let sk = hex_field(test, "sk");
            let msg = hex_field(test, "message");
            let ctx = hex_field(test, "context");
            let rnd_vec = if deterministic { vec![0u8; 32] } else { hex_field(test, "rnd") };
            let mut rnd = [0u8; 32];
            rnd.copy_from_slice(&rnd_vec);
            let sig = match param_set {
                "ML-DSA-44" => ml_dsa::ml_dsa_44::hash_sign_with_rnd(&sk, &msg, &ctx, ph, &rnd),
                "ML-DSA-65" => ml_dsa::ml_dsa_65::hash_sign_with_rnd(&sk, &msg, &ctx, ph, &rnd),
                "ML-DSA-87" => ml_dsa::ml_dsa_87::hash_sign_with_rnd(&sk, &msg, &ctx, ph, &rnd),
                _ => continue,
            }
            .expect("valid context length");
            if sig == hex_field(test, "signature") {
                t.pass += 1;
            } else {
                t.fail += 1;
                eprintln!(
                    "  {param_set} preHash sigGen tcId {} ({}) FAILED",
                    test["tcId"], test["hashAlg"]
                );
            }
        }
    }
    t
}

/// HashML-DSA (pre-hash, Algorithm 5) sigVer KAT for one parameter set.
fn run_prehash_sigver(param_set: &str) -> Tally {
    let json = vectors("sigver.json");
    let mut t = Tally::default();
    for group in json["testGroups"].as_array().unwrap() {
        if group["parameterSet"].as_str() != Some(param_set) {
            continue;
        }
        if group["preHash"].as_str() != Some("preHash")
            || group["externalMu"].as_bool().unwrap_or(false)
        {
            continue;
        }
        for test in group["tests"].as_array().unwrap() {
            let ph = match ml_dsa::PreHash::from_acvp(test["hashAlg"].as_str().unwrap_or("")) {
                Some(p) => p,
                None => {
                    t.skipped += 1;
                    continue;
                }
            };
            let pk = hex_field(test, "pk");
            let sig = hex_field(test, "signature");
            let msg = hex_field(test, "message");
            let ctx = hex_field(test, "context");
            let expected = test["testPassed"].as_bool().unwrap();
            let got = match param_set {
                "ML-DSA-44" => ml_dsa::ml_dsa_44::hash_verify(&pk, &msg, &sig, &ctx, ph),
                "ML-DSA-65" => ml_dsa::ml_dsa_65::hash_verify(&pk, &msg, &sig, &ctx, ph),
                "ML-DSA-87" => ml_dsa::ml_dsa_87::hash_verify(&pk, &msg, &sig, &ctx, ph),
                _ => continue,
            };
            if got == expected {
                t.pass += 1;
            } else {
                t.fail += 1;
                eprintln!(
                    "  {param_set} preHash sigVer tcId {} ({}) FAILED (expected {expected})",
                    test["tcId"], test["hashAlg"]
                );
            }
        }
    }
    t
}

fn main() {
    let mut failures = 0;
    for api in APIS.iter().flatten() {
        let k = run_keygen_kat(api);
        println!("{} keyGen ACVP KAT: {} passed, {} failed", api.label, k.pass, k.fail);
        let s = run_siggen_kat(api);
        println!(
            "{} sigGen ACVP KAT: {} passed, {} failed, {} skipped",
            api.label, s.pass, s.fail, s.skipped
        );
        let v = run_sigver_kat(api);
        println!(
            "{} sigVer ACVP KAT: {} passed, {} failed, {} skipped",
            api.label, v.pass, v.fail, v.skipped
        );
        failures += k.fail + s.fail + v.fail;
    }
    // HashML-DSA (pre-hash) — Algorithms 4 & 5, once per parameter set.
    for ps in PARAM_SETS {
        let s = run_prehash_siggen(ps);
        println!(
            "{ps} HashML-DSA sigGen KAT: {} passed, {} failed, {} skipped",
            s.pass, s.fail, s.skipped
        );
        let v = run_prehash_sigver(ps);
        println!(
            "{ps} HashML-DSA sigVer KAT: {} passed, {} failed, {} skipped",
            v.pass, v.fail, v.skipped
        );
        failures += s.fail + v.fail;
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
        for api in APIS.iter().flatten() {
            let t = run_keygen_kat(api);
            assert!(t.pass > 0, "no {} keyGen tests found", api.label);
            assert_eq!(t.fail, 0, "{} {} keyGen failures", t.fail, api.label);
        }
    }

    #[test]
    fn siggen_kat_all_sets() {
        for api in APIS.iter().flatten() {
            let t = run_siggen_kat(api);
            assert!(t.pass > 0, "no {} sigGen tests found", api.label);
            assert_eq!(t.fail, 0, "{} {} sigGen failures", t.fail, api.label);
        }
    }

    #[test]
    fn sigver_kat_all_sets() {
        for api in APIS.iter().flatten() {
            let t = run_sigver_kat(api);
            assert!(t.pass > 0, "no {} sigVer tests found", api.label);
            assert_eq!(t.fail, 0, "{} {} sigVer failures", t.fail, api.label);
        }
    }

    #[test]
    fn prehash_kat_all_sets() {
        for ps in PARAM_SETS {
            let s = run_prehash_siggen(ps);
            assert!(s.pass > 0, "no {ps} HashML-DSA sigGen tests found");
            assert_eq!(s.fail, 0, "{} {ps} HashML-DSA sigGen failures", s.fail);
            let v = run_prehash_sigver(ps);
            assert!(v.pass > 0, "no {ps} HashML-DSA sigVer tests found");
            assert_eq!(v.fail, 0, "{} {ps} HashML-DSA sigVer failures", v.fail);
        }
    }
}

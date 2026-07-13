# Decoder fuzzing (ML-DSA-65): pk / sig / sk decode

`cargo-fuzz` 0.13.2 (libFuzzer) on nightly, AddressSanitizer, over the adversary-facing
decoders. Crate: `fuzz/` (detached from the workspace, nightly-pinned). Corpus seeded
with one valid artifact each from the ACVP KAT vectors so the fuzzer reaches the
successful-decode (round-trip) path immediately — a random 1952/3309/4032-byte input is
astronomically unlikely to hit it otherwise.

## Properties checked

| target       | property |
|--------------|----------|
| `pk_decode`  | no panic on arbitrary bytes; on `Ok`, `pk_encode(decode) == input` (exact round-trip — no non-canonical key accepted) |
| `sig_decode` | no panic; malformed hint ⇒ ⊥; on `Ok`, `sig_encode(decode) == input` (canonical `z` and hint) |
| `sk_decode`  | no panic (trusted input — FIPS allows out-of-range values, so no round-trip guarantee) |

The round-trip assertions are enforced inside the fuzz target, so any decoder that
accepted a non-canonical encoding would trip an assertion and be reported as a crash.

## Result (40 s / target this run)

| target       | executions | crashes | corpus | notes |
|--------------|-----------:|--------:|-------:|-------|
| `pk_decode`  |   184,093  | **0**   | 2      | length-gated (needs exactly 1952 B) |
| `sig_decode` |   129,466  | **0**   | 55     | richest — the hint decoder's branches explored (cov 180) |
| `sk_decode`  |   335,972  | **0**   | 2      | length-gated (4032 B) |

No crashes, no ASan errors, no round-trip assertion failures, no artifacts written.
`#![forbid(unsafe_code)]` in the library means memory-safety holds by construction;
fuzzing additionally exercises the hand-written length checks and the hint decoder's
canonical-form validation (`HintBitUnpack → ⊥`).

## Reproduce / extend

```
cargo +nightly fuzz run pk_decode   # or sig_decode / sk_decode
```

This session's 40 s/target is a smoke run. For the thesis, the recommended budget is
**30–60 min/target** (or to a coverage plateau); the committed `fuzz/corpus/` carries
over between runs. Fuzzing is evidence of absence of shallow bugs, not a proof — stated
as such.

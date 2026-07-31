# ml-dsa

An implementation of ML-DSA (FIPS 204), the NIST post-quantum signature scheme, in
plain safe Rust. Alongside the implementation there's the benchmarking, side-channel,
and fault-analysis work that was done on top of it.

The whole crate is `#![forbid(unsafe_code)]`, so there's no unsafe anywhere.

## How it's put together

Most operations exist in two forms. There's a reference version that follows the
FIPS 204 pseudocode as literally as Rust allows (plain `mod q` arithmetic, same
branching as the spec), and an optimised version (`*_fast`, Montgomery NTT, the
reusable `Signer`) that's free to do things differently as long as it produces
exactly the same bytes as the reference for any input.

Keeping them byte-for-byte identical is the point: it means the known-answer tests
check the fast path too, not just the reference. On top of that a few of the
sensitive primitives have branchless `_ct` versions that are constant-time by
construction.

All three parameter sets (44, 65, 87) are supported through the `ParameterSet` trait
with `K`/`L` as const generics. Both signing interfaces are there — plain ML-DSA and
the HashML-DSA pre-hash variants (all twelve approved hashes). `externalMu` is the
one thing left out.

| Set | Public key | Secret key | Signature |
|-----|-----------:|-----------:|----------:|
| ML-DSA-44 | 1312 B | 2560 B | 2420 B |
| ML-DSA-65 | 1952 B | 4032 B | 3309 B |
| ML-DSA-87 | 2592 B | 4896 B | 4627 B |

## The crates

It's a Cargo workspace with four crates and a `results/` directory.

### `ml-dsa`

The library itself. Every function carries its FIPS reference in a doc comment
(`FIPS 204, Algorithm N — Name`). The modules line up with the signing pipeline:

- arithmetic: `field`, `ntt` and `ntt_arith` (reference and Montgomery), `poly`,
  `vecops`, `rounding`, `hint`
- sampling and expansion: `sample`, `expand`, `hash` (SHAKE/Keccak)
- encoding: `encoding`, `serdes` — the bit-packing of keys and signatures
- operations: `keygen`, `sign`, `verify`, plus `signer`, which caches the one-time
  key setup so you can sign many messages cheaply in a session
- interfaces: `prehash` (HashML-DSA), `params` (the trait and `MlDsa44/65/87`),
  `ct`, `error`

You call it through a per-parameter-set module — `mldsa44`, `mldsa65`, `mldsa87` —
each of which gives you `key_gen`, `sign`, `verify`, the `*_internal` and `*_fast`
variants, `hash_sign`/`hash_verify`, and the `PK_BYTES`/`SK_BYTES`/`SIG_BYTES`
constants. You always pass in your own RNG (`rand_core::CryptoRng`), which is what
makes deterministic KAT replay work.

There are two examples: `roundtrip.rs` (keygen, sign, verify) and `profile.rs`
(what the callgrind baseline runs against).

### `ml-dsa-kat`

The ACVP known-answer-test runner. The vectors are in `vectors/` (`keygen.json`,
`siggen.json`, `sigver.json`) and get run through both the reference and fast paths,
for all three sets and both interfaces.

```
cargo run -p ml-dsa-kat
```

Everything should come back green with 0 failed.

### `ml-dsa-bench`

Instruction-count benchmarks using `iai-callgrind`. Counting retired instructions
instead of wall-clock time keeps before/after comparisons exact and repeatable.

- `benches/ntt.rs` — the NTT, reference vs Montgomery
- `benches/operations.rs` — keygen, sign, verify
- `benches/signer_amortized.rs` — the reusable `Signer` across a session
- `examples/latency.rs` — the one actual wall-clock timing (thesis §7.3)

```
cargo bench -p ml-dsa-bench
```

### `ml-dsa-leakage`

The constant-time and fault harnesses.

- `benches/ct_invariance.rs` — checks instruction counts don't vary with secret
  data (a control-flow leak test under callgrind)
- `benches/dudect_ct.rs` — the statistical timing test (Welch t-test) via dudect
- `examples/timecop.rs` — whole-program constant-time audit using Valgrind client
  requests (`crabgrind`)
- `examples/cachegrind_ball.rs` — cache behaviour of `SampleInBall`
- `examples/fault_2x2.rs` — the fault experiment (signing mode against
  verify-after-sign)
- `src/fault.rs` — a copy of the signer wired up for fault injection

### `results`

The measurements the thesis refers to. The `.md` files are the write-ups —
`README.md` (baseline profile), `methods.md` (what's tested and why), `bench.md`
(before/after), `ct.md` (constant-time), `fault.md` (fault analysis), `latency.md` —
and the rest is the raw callgrind/cachegrind/dudect/timecop output behind them.

## Building and checking

```
cargo run -p ml-dsa-kat       # the KAT run
cargo test                    # unit tests, KATs, fault-fidelity tests
cargo clippy --all-targets    # kept at zero warnings
```

The toolchain is pinned in `rust-toolchain.toml` (stable 1.96.0). The library's only
dependencies are `sha3`, `sha2`, `subtle`, `zeroize`, and `rand_core`.

## License

MIT or Apache-2.0.

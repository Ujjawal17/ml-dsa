# ML-DSA (FIPS 204) — Technical Specification

A pure-Rust implementation of FIPS 204 (ML-DSA): a **faithful reference** transcription and an
**optimized, hardened** variant derived from it. Correctness is defined by **byte-exact agreement
with the NIST ACVP known-answer tests**.

## 1. Scope

- Algorithms: KeyGen, Sign (hedged + deterministic), Verify — both the external (§5) and internal
  (§6) interfaces — plus all auxiliary functions (§7, Algorithms 9–49).
- Parameter sets: ML-DSA-44, ML-DSA-65, ML-DSA-87 (§2).
- Interfaces validated: pure (non-prehash), non-`externalMu`.
- **Not covered:** HashML-DSA / pre-hash (Algorithms 4–5), AVX2/NEON vectorisation, masking
  countermeasures, `no_std`/embedded targets. SHAKE/Keccak is consumed from `sha3` (treated as a
  primitive below the line, not re-implemented).

## 2. Parameters

Common to all sets: `q = 8 380 417` (= 2²³ − 2¹³ + 1), `N = 256`, `d = 13`.

| Parameter | ML-DSA-44 | ML-DSA-65 | ML-DSA-87 |
|---|---|---|---|
| (k, l) | (4, 4) | (6, 5) | (8, 7) |
| η | 2 | 4 | 2 |
| τ | 39 | 49 | 60 |
| β = τ·η | 78 | 196 | 120 |
| γ1 | 2¹⁷ | 2¹⁹ | 2¹⁹ |
| γ2 | (q−1)/88 | (q−1)/32 | (q−1)/32 |
| ω | 80 | 55 | 75 |
| λ | 128 | 192 | 256 |
| c̃ bytes (λ/4) | 32 | 48 | 64 |
| public key bytes | 1312 | 1952 | 2592 |
| secret key bytes | 2560 | 4032 | 4896 |
| signature bytes | 2420 | 3309 | 4627 |

## 3. Core types and representation

- `Zq` = `i32`. Ring elements are held in the **canonical range `[0, q)`** (matching the spec's
  `mod q`). The centred representative (`mod± α`, in `(−α/2, α/2]`) is used **only** where the spec
  requires it: Power2Round, Decompose's `r0`, and the signing infinity-norm checks.
- `Poly` = `[Zq; 256]` (an element of `R_q = Z_q[X]/(X²⁵⁶+1)`).
- `PolyNTT` = a distinct type for the NTT image `T_q` — prevents pointwise-multiplying a non-NTT
  polynomial by accident.
- `PolyVec<const K>`, `PolyVecNTT<const K>`, `PolyMatNTT<const K, const L>`.
- Secrets that **persist beyond a single call** — the prepared signer's cached `K, ŝ1, ŝ2, t̂0`
  — are zeroised on drop. Transient stack copies inside one call are out of scope: safe Rust
  cannot scrub compiler-made moves/copies (a documented limitation, not an oversight).
- The crate is compiled under `#![forbid(unsafe_code)]` (compiler-enforced memory safety).

## 4. Algorithm inventory (FIPS 204)

- **External (§5):** 1 KeyGen · 2 Sign (hedged + deterministic) · 3 Verify.
- **Internal (§6):** 6 KeyGen_internal · 7 Sign_internal · 8 Verify_internal.
- **Data conversion (§7.1):** 9–13 IntegerToBits/Bytes, Bits/Bytes conversions · 14
  CoeffFromThreeBytes · 15 CoeffFromHalfByte · 16–19 (Simple)BitPack/Unpack.
- **Encodings (§7.2):** 20–21 HintBitPack/Unpack · 22–25 pk/sk Encode/Decode · 26–27 sig
  Encode/Decode · 28 w1Encode.
- **Sampling (§7.3):** 29 SampleInBall · 30 RejNTTPoly · 31 RejBoundedPoly · 32 ExpandA · 33
  ExpandS · 34 ExpandMask.
- **Rounding & hints (§7.4):** 35 Power2Round · 36 Decompose · 37 HighBits · 38 LowBits · 39
  MakeHint · 40 UseHint (applied componentwise over vectors).
- **NTT (§7.5–7.6):** 41 NTT · 42 NTT⁻¹ · 43 BitRev8 · 44–48 Add/Multiply/vector/matrix NTT · 49
  MontgomeryReduce.

## 5. Implementation conventions

- **Faithful baseline:** the most direct literal transcription of the FIPS 204 pseudocode, using
  **plain `mod q`** arithmetic. Montgomery multiplication (Algorithm 49 / Appendix A) is treated by
  the standard as an optimisation and therefore belongs only to the optimised variant (§7).
- **NTT zeta convention (KAT-observable):** `zetas[m] = ζ^BitRev8(m) mod q` with `ζ = 1753`. The
  table is computed at compile time from `ζ` and `BitRev8` so it is auditable rather than an opaque
  blob. NTT⁻¹ finishes by multiplying by `256⁻¹ mod q = 8 347 681`.
- **Bit/byte order:** little-endian throughout (Algorithms 9–13).
- **Message formatting (pure variant):** `M' = IntegerToBytes(0,1) ‖ IntegerToBytes(|ctx|,1) ‖ ctx ‖ M`,
  with `|ctx| ≤ 255` (else reject).
- **Hashing chain (Sign/Verify):** `μ = H(tr ‖ M', 64)`; `ρ'' = H(K ‖ rnd ‖ μ, 64)`;
  `c̃ = H(μ ‖ w1Encode(w1), λ/4)`; `tr = H(pk, 64)`. `H` = SHAKE256, `G` = SHAKE128.
- **Deterministic signing variant:** `rnd = {0}³²`. Hedged (default) draws `rnd` from an injected
  `CryptoRng`; the RNG is always injected (never `getrandom` directly) so KATs replay deterministically.
- **Decode validation:** decoders that may receive adversary-controlled bytes (`pk_decode`,
  `sig_decode`) length-check first and return an error on mismatch; a malformed hint
  (`HintBitUnpack → ⊥`) makes Verify return false.

## 6. Parameterisation (multi-set support)

- A **`ParameterSet` trait** carries the scalar parameters (η, τ, β, γ1, γ2, ω, λ, c̃-bytes, and the
  encoded key/signature sizes). `K` and `L` are threaded as **explicit `const` generics** on the
  functions (and on `PolyVec`/`PolyMat`).
- **Stable-Rust constraint:** trait associated consts cannot be used as array lengths
  (`PolyVec<{P::K}>` would require the nightly `generic_const_exprs` feature). Consequently `K`/`L`
  are explicit const generics, and the hint encoding uses a `Vec<u8>` rather than `[u8; ω+K]`.
- **Per-set behavioural differences (not just constant values):**
  - **η = 2 (ML-DSA-44, -87)** vs **η = 4 (ML-DSA-65)** — `CoeffFromHalfByte` has distinct branches.
  - **ML-DSA-44: γ1 = 2¹⁷** → `z` packs at 18 bits/coefficient (vs 20), and **γ2 = (q−1)/88**.

## 7. Optimised variant (performance)

The optimised path must produce **byte-identical output** to the reference (so it passes the same
ACVP KAT — the optimisation's correctness proof) while executing fewer instructions.

- **(a) Lazy/deferred-reduction NTT + Montgomery (Algorithm 49) + merged butterfly layers** — the
  arithmetic win, replacing the baseline's `rem_euclid` (a hardware division) per the standard's
  Appendix A.
- **(b) Amortised prepared signer** — `ExpandA(ρ)` and `NTT(s1/s2/t0)` are computed **once at
  key-load** and reused across many signatures. Per-signature output is unchanged; this removes
  `ExpandA` (a major SHAKE cost) from the per-signature path.
- **(c) SHAKE squeeze buffering** — squeeze rate-sized blocks in rejection sampling / ExpandMask.

**Correctness guardrails (the optimised path has no KAT of its own):**
1. Run the ACVP KAT **through the optimised path** as well as the reference.
2. A randomised **`optimised == reference` equivalence check** (`proptest` + `cargo-fuzz`) over many
   random inputs — required to catch lazy-reduction overflow boundaries that fixed KAT vectors
   cannot exercise.

**Measurement methodology — retired-instruction counts:**
- Metric = retired-instruction counts via **`callgrind` / `iai-callgrind`**. These are
  **deterministic and machine-independent** (identical for a given build on any machine), hence
  reproducible. Real cycles / wall-clock (`rdtsc`, `perf`) are **not** used as the optimisation
  metric: they vary across microarchitectures, are noisy, and require `unsafe`.
- Rank functions by **self-instructions** to identify the hot spot and the Amdahl ceiling (if the
  Keccak permutation accounts for N% of self-instructions, arithmetic optimisation is bounded by
  ~(100−N)%).
- Report **before/after** instruction counts per component and end-to-end on a multi-signature
  workload (where the amortised signer's benefit is visible).

## 8. Hardening

- **`#![forbid(unsafe_code)]`** — pure safe Rust; no unsafe blocks anywhere in the crate, including
  the measurement harnesses (instruction counting via `callgrind` and timing via `std::time::Instant`
  require no `unsafe`).
- **Verify-after-sign** — internally run Verify before releasing a signature; a countermeasure
  against differential fault attacks (Bruinderink–Pessl).
- **Decoder fuzzing** — `cargo-fuzz` targets over `pk_decode` / `sig_decode` / `sk_decode` (parsing
  adversary-controlled bytes), checking for panics and out-of-range outputs.

## 9. Constant-time analysis

Constant-time has two independent dimensions:
- **(i) Control-flow** — no secret-dependent branches.
- **(ii) Data-memory** — no secret-dependent memory addresses (cache-timing).

- The faithful baseline mirrors the spec's `if`/`else` (e.g. Decompose's `q−1` case, the signing
  norm checks, UseHint) and is therefore **not** constant-time; it is the "before".
- The hardened variant replaces secret-dependent branches with **branchless selection**
  (`subtle::ConditionallySelectable`, or the mask `a += (a >> 31) & q`), and the difference is
  **measured before/after**.
- **The claim is scoped to dimension (i), control-flow.** Two dimension-(ii) leaks cannot be removed
  by branch elimination and are **documented as accepted leakage** (consistent with reference
  implementations): `SampleInBall`'s secret-indexed write `c[j] = c[i]`, and the rejection-loop
  iteration count.
- **Verification (all pure-safe):**
  - `callgrind` **instruction-count invariance** under fixed-vs-random secret inputs — a
    *deterministic* control-flow leak test (a function whose instruction count depends on the secret
    contains a control-flow leak).
  - `dudect` statistical timing test (Welch's t-test) using `std::time::Instant`.
  - `cachegrind` to demonstrate the data-dependent memory-access pattern of `SampleInBall`.
- **Stated limitation:** source-level constant-time ≠ binary-level constant-time (the compiler may
  lower a branchless selection back into a conditional branch).

## 10. Fault model

- **Emulated** single-bit fault: a bit is flipped in a named signing intermediate (e.g. a
  coefficient of `c·s1` after `MultiplyNTT`), modelling the *effect* of a fault rather than a
  physical glitch.
- Success criterion = whether an **exploitable correct/faulty signature artefact escapes** (full
  key recovery is cited from Bruinderink–Pessl, not re-implemented).
- Reported as a 2×2 over **{deterministic, hedged}** signing × **{verify-after-sign on, off}**.

## 11. Dependencies

- **Runtime:** `sha3` (SHAKE128/256), `subtle` (constant-time selection), `zeroize` (clear secret
  intermediates), `rand_core` (injected `CryptoRng`).
- **Analysis / test:** `proptest`, `iai-callgrind` + `callgrind`/`cachegrind` (Valgrind),
  `criterion` (wall-clock, used only where real throughput is required), a `dudect`-style harness
  (`std::time::Instant`), `cargo-fuzz`, `serde_json` (ACVP vector parsing).

## 12. Open questions for the supervisor

1. **Parameter-set scope:** are ML-DSA-44 and ML-DSA-87 required, or is ML-DSA-65 sufficient — i.e.
   does "NIST KAT 100%" mean 100% on the implemented set(s) with ML-DSA-65 as the floor, or all
   three sets mandatory? *(This determines whether §6 parameterisation is undertaken.)*
2. **Optimisation aspect:** does a **performance** optimisation satisfy the "improve/optimise in some
   aspect" requirement?
3. **Optimisation baseline:** is an improvement over the author's **own baseline** sufficient, or must
   the optimised version be **competitive with the state of the art** (the C reference / RustCrypto
   `ml-dsa`)?
4. **Profiling granularity:** report instruction counts at the **top-level operations**
   (KeyGen/Sign/Verify) or **per function**?
5. **Optional analyses:** which to include — cache-timing demonstration (`cachegrind` on
   `SampleInBall`), rejection-count distribution, heap-allocation profiling (`dhat`), and/or a
   state-of-the-art performance comparison?
6. **Constant-time claim:** is a claim scoped to **control-flow**, with `SampleInBall`'s
   secret-indexed memory access documented as accepted leakage (as in reference implementations),
   acceptable — or is a stronger data-memory guarantee expected?

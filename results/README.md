# Baseline instruction-count profile (ML-DSA-65)

Retired-instruction counts (`Ir`) from Valgrind **callgrind**. Deterministic and
machine-independent: identical for a given build on any host, so a single call per
operation is exact (no iteration averaging needed for the counts).

## Reproduce

```
cargo build --release --example profile          # [profile.release] debug = true
mkdir -p results
BIN=target/release/examples/profile

# Each operation is an #[inline(never)] bench_* wrapper; --toggle-collect starts
# collection on entry and includes all callees, so preceding setup is excluded.
valgrind --tool=callgrind --collect-atstart=no --toggle-collect='*bench_keygen*' \
  --callgrind-out-file=results/cg.baseline.keygen.out $BIN
valgrind --tool=callgrind --collect-atstart=no --toggle-collect='*bench_sign*' \
  --callgrind-out-file=results/cg.baseline.sign.out   $BIN
valgrind --tool=callgrind --collect-atstart=no --toggle-collect='*bench_verify*' \
  --callgrind-out-file=results/cg.baseline.verify.out $BIN

callgrind_annotate --auto=no --inclusive=no results/cg.baseline.<op>.out   # flat, self
callgrind_annotate --auto=no --inclusive=yes results/cg.baseline.sign.out  # inclusive
```

Fixed seed `0x42`; empty context, deterministic signing (`rnd = 0^32`). The default
message accepts on rejection-loop iteration 1. **Sign cost is not constant** — see
"Rejection-count dependence" below; the self-breakdown table is the n=1 case.

Toolchain: rustc/cargo 1.96.0 (release + `debug = true`). Function costs are summed
across the source files they inline from (callgrind splits e.g. `keccak_p` across
`bit.rs` + `intrinsics/mod.rs`).

## Self-instruction breakdown (each op = 100%)

| subsystem            | KeyGen | Sign  | Verify |
|----------------------|-------:|------:|-------:|
| **SHAKE total**      | 34.3%  | 21.7% | 31.6%  |
|   – Keccak permute   | 28.7%  | 18.1% | 26.0%  |
|   – SHAKE plumbing   |  5.5%  |  3.6% |  5.6%  |
| NTT / poly arith     | 18.3%  | 35.7% | 34.5%  |
| encoding / serdes    | 20.9%  | 23.0% | 19.0%  |
| libc / alloc         | 19.4%  |  9.3% |  9.6%  |
| sampling / expand*   |  4.4%  |  2.2% |  3.3%  |
| rounding / field     |  0.3%  |  1.8% |  1.2%  |

`*` self-cost of the sampling control code only; its SHAKE cost is inside "SHAKE total".

Totals (Ir): KeyGen 4,954,338 · Sign 7,614,873 · Verify 4,922,544.

## Derived quantities

**Amdahl bound (SHAKE is `sha3`, below the line).** Optimizing everything except SHAKE
can remove at most the non-SHAKE fraction, so the max conceivable speedup is
`total / SHAKE = 1 / SHAKE_fraction` (a loose upper bound — non-SHAKE work cannot
actually be driven to zero):

| op          | SHAKE | max speedup (1/SHAKE) |
|-------------|------:|----------------------:|
| KeyGen      | 34.3% | 2.9×                  |
| Sign (n=1)  | 21.7% | 4.6×                  |
| Sign (n=24) |  9.6% | 10.4×                 |
| Verify      | 31.6% | 3.2×                  |

Realistic NTT-only Montgomery (assume ~2× on the NTT self-cost, n=1 Sign, NTT = 35.7%):
`1 / (1 − 0.357·0.5) ≈ 1.2×` end-to-end on Sign.

## Rejection-count dependence (Sign)

Sign is a Fiat–Shamir-with-aborts loop; total cost is linear in the accepted
iteration count `n` (fixed prefix + `n` × loop body). Measured at the two extremes
found by `profile scan` (n from 1 to 24 over 200 messages, same key):

```
total(n) = 4,636,407 + n · 2,978,286   Ir      (exact fit, 2 points)
  n=1  ->  7,614,693 Ir       n=24 -> 76,115,279 Ir
```

`fixed` = ExpandA + skDecode + NTT(s1,s2,t0) + final sig_encode.
`per_iter` = ExpandMask + A·y + SampleInBall + cs1/cs2/ct0 + rounding + w1Encode.

**Proportions shift with n** (they are *not* iteration-stable):

| subsystem      | n=1   | n=24  |
|----------------|------:|------:|
| SHAKE          | 21.7% |  9.6% |
| NTT / arith    | 35.7% | 45.2% |
| encoding       | 23.0% | 21.7% |
| loop overhead  |  6.3% | 12.6% |

As `n` grows the fixed ExpandA (SHAKE) is diluted and the per-iteration NTT dominates.

**Amortized-`Signer` benefit is n-dependent.** `ExpandA` inclusive = 1,741,814 Ir,
fixed. Its share of a signature — the per-signature saving from precomputing it in
`Signer::from_sk` — falls as the rejection count rises:

| n (iterations) |  1  |  2  |  4  |  5  | 24  |
|----------------|----:|----:|----:|----:|----:|
| Sign total (M) | 7.6 |10.6 |16.6 |19.5 |76.1 |
| ExpandA share  |22.9%|16.4%|10.5%| 8.9%| 2.3%|

ML-DSA-65's expected iteration count is ~4–5, so the *typical*-signature Signer win is
~9–11%, not the 22.9% of the lucky n=1 case; largest for low-rejection signatures.

## Squeeze-buffering decision

SHAKE plumbing (block-buffer bookkeeping around the permutation) is only ~3.6–5.6%;
`block_buffer` already caches rate-sized output blocks, so the permutation count is
near-minimal. Rate-sized squeeze buffering can recover at most that plumbing slice →
**minor**, do only if cheap. Larger non-arithmetic target: the faithful bit-vector
encoding (bit_pack / bits_to_bytes / bytes_to_bits), 19–23%, currently outside the
optimization list.

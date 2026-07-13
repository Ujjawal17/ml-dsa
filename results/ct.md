# Constant-time audit (ML-DSA-65): which functions are constant-time?

Method: **callgrind instruction-count invariance** (deterministic; counts are exact).
Each function is run over two secret input classes; **equal counts ⇒ constant-time**,
**differing counts ⇒ a control-flow / variable-loop leak**. Build: bench profile
(`opt-level = 3`, `debug = true`). Cross-checked with dudect (real timing) on
`decompose`.

## Battery result

| function (category)            | class A | class B | verdict            |
|--------------------------------|--------:|--------:|--------------------|
| `decompose` baseline (round)   |    4771 |    4771 | **constant-time**  |
| `decompose` `_ct`              |    3586 |    3586 | constant-time (leaner) |
| `make_hint` baseline (hint)    |    9501 |    9501 | **constant-time**  |
| `make_hint` `_ct`              |    4897 |    4897 | constant-time (leaner) |
| `montgomery_reduce` (arith)    |    2615 |    2615 | **constant-time**  |
| `sample_in_ball` (sampling)    |   18101 |   17748 | **NOT** (Δ 353)    |
| `rej_bounded_poly` (sampling)  |   36238 |   35246 | **NOT** (Δ 992)    |

(Batch = 256 coefficients, except the samplers which are one call on a seed.)

dudect cross-check on `decompose` (80k samples, pre-built inputs): baseline max |t| =
1.58, `_ct` max |t| = 2.84 — both < 4.5, i.e. both timing-invariant, agreeing with the
callgrind verdict.

## What this says

**Constant-time (already, at `-O3`):** the per-coefficient rounding (`decompose`,
`high_bits`/`low_bits`), hint (`make_hint`), and arithmetic (`montgomery_reduce`, and
the NTT built from it). The compiler lowers the baseline's source-level `if`s (the
`q−1` special case, `mod±`'s conditional subtract, the hint comparison) into branchless
cmov/select, so the faithful baseline is **already binary-CT for these**. The hand-
written `_ct` variants are also CT and additionally leaner (`decompose` 3586 vs 4771;
`make_hint` 4897 vs 9501 — ~2×), so the hardening costs no performance — it improves it.
Their real value is **robustness**: constant-time guaranteed at the source level,
independent of the optimizer (the source-vs-binary gap, spec §9).

**NOT constant-time (by design — documented accepted leakage):** the rejection samplers
`sample_in_ball`, `rej_bounded_poly` (and `rej_ntt_poly`), whose **iteration count is
data-dependent**. The battery correctly flags them (Δ ≠ 0), which also validates that
the method detects real leaks rather than always reporting "CT". These are the Fiat–
Shamir-with-aborts / rejection-sampling leaks accepted by all reference implementations.

## Scope of this method (and the whole-program check)

callgrind instruction counts detect **control-flow / loop-count** leaks, not pure
**data-memory (cache)** leaks. `sample_in_ball` is flagged here via its loop count; its
*other* leak — the secret-indexed write `c[j] = c[i]` — is a data-address leak that
needs **cachegrind** (or the whole-program taint audit below).

Whole-program audit (`examples/timecop.rs`, ctgrind/TIMECOP: poison the secret key,
run `sign` under memcheck, catch every secret-dependent branch and address at once) is
**built but currently blocked on this host**: Valgrind memcheck fails to start on this
Arch/glibc build ("mandatory `memcmp` redirection in `ld.so`"), and neither debuginfod
nor the archlinux debuginfod server supplies the needed `ld.so` symbols. Unblocking it
needs Arch's debug repos + `glibc-debug`. The per-function callgrind battery above is
the working, deterministic substitute.

## Cachegrind: SampleInBall's data-dependent memory access (`examples/cachegrind_ball.rs`)

`sample_in_ball` run under `cachegrind --cache-sim=yes` for two seeds (bytes 7 vs 199).

Aggregate totals (secret-dependent memory-access **counts**, the loop-count component):

| seed | I refs  | D refs (rd + wr)          | D1 misses |
|------|--------:|---------------------------|----------:|
| 7    | 347,456 | 121,525 (84,093 + 37,432) | 3,023     |
| 199  | 347,138 | 121,391 (84,018 + 37,373) | 3,023     |

Per-source-line (`cg_annotate --auto=yes --show=Dr,Dw`), the two leak components separate
cleanly:

```
                                          seed 7 (Dr/Dw)   seed 199 (Dr/Dw)
  reader.squeeze(&mut byte)  // rej. loop     0 / 11           0 / 5      <- COUNT varies
  c.coeffs[i] = c.coeffs[j]  // line 11      49 / 49          49 / 49     <- fixed count,
  c.coeffs[j] = (-1)^h       // line 12      49 / 49          49 / 49        secret ADDRESS j
```

Two distinct leaks, both documented accepted leakage:
- **Rejection-loop count** — the inner `squeeze` runs a secret-dependent number of times
  (11 vs 5), so data-access counts (and I refs, D refs) vary with the secret.
- **Secret-indexed memory access** `c[j]` — executes a **fixed** τ = 49 times regardless
  of the secret (identical 49/49 counts), so this is *not* a count leak: the leak is the
  **address** `j` (a secret-derived cache line). This is why the **D1 miss totals are
  identical (3,023)** — the 1 KB `c` array is L1-resident, so the address-dependence is
  invisible to aggregate cache statistics and shows only in the access *pattern*. On real
  hardware a cache-line-probing adversary (Prime+Probe) recovers `j`; that is the cited
  attack surface, not re-implemented here.

# Before/after instruction counts (ML-DSA-65): faithful vs improved

Retired instructions via `iai-callgrind` 0.14.0 (Valgrind callgrind), `cargo bench
-p ml-dsa-bench`. Deterministic; fixed seed `0x42`. Reference and improved paths are
byte-identical, so each pair signs the same messages under the same key and hits the
same rejection-loop counts — the delta is purely the optimization.

## NTT in isolation (pure arithmetic: `rem_euclid` vs Montgomery + deferred reduction)

| transform | baseline | improved | speedup |
|-----------|---------:|---------:|--------:|
| `ntt`     |   55,219 |   19,421 | 2.84×   |
| `inv_ntt` |   63,535 |   25,817 | 2.46×   |

## Top-level operations (single call)

| operation            | baseline  | improved  | reduction | speedup |
|----------------------|----------:|----------:|----------:|--------:|
| KeyGen               | 4,954,181 | 4,496,118 |  9.2%     | 1.10×   |
| Verify               | 4,924,648 | 3,999,662 | 18.8%     | 1.23×   |
| Sign (one-shot, n=1) | 7,615,123 | 6,940,195 |  8.9%     | 1.10×   |

One-shot Sign constructs a throwaway `Signer` (pays full key setup incl. pk
reconstruction), and this message accepts on iteration 1 (short loop → little NTT to
speed up), so the one-shot figure is the *weakest* case for the improved path — as
expected, the win is in the amortized workload below.

## Amortized multi-signature workload (one key, N = 10 signatures)

| path                        | total (10 sigs) | per signature |
|-----------------------------|----------------:|--------------:|
| baseline (`sign_internal`)  |     210,468,922 |    21,046,892 |
| improved (`Signer`)         |     128,313,615 |    12,831,362 |
| **reduction**               |   **39.0%**     |   **1.64×**   |

These 10 messages average ~5.5 rejection iterations (`(21.05M − 4.64M)/2.98M`), i.e.
realistic signatures (ML-DSA-65 expected ~4–5), not the lucky n=1 case.

## Reading

The amortized 1.64× combines **two** effects: (a) `ExpandA + skDecode + NTT(s1,s2,t0)`
run once instead of per signature (~1.74M/sig saved), and (b) Montgomery/lazy NTT in
the shared loop body. Effect (b) scales with rejection count (more iterations → more
NTTs → more Montgomery savings), so the improved path's advantage *grows* with the
loop count, and the amortized/realistic figure (1.64×) far exceeds the one-shot n=1
figure (1.10×). Verify (no rejection loop, NTT-heavy) shows the arithmetic win cleanly
at 1.23×; the isolated NTT is 2.5–2.8×, diluted at the operation level by SHAKE and
encoding per the baseline profile.

iai-callgrind stored these as the saved baseline; re-running `cargo bench` reports
deltas against them (regression guard).

# Differential fault analysis (ML-DSA-65): verify-after-sign countermeasure

Emulated single-bit fault at `cs1 = c·s1` (bit 3 of coefficient 0, polynomial 0),
injected in a **mirror** of `sign_internal` reassembled from the library's public
building blocks (`examples/fault_2x2.rs`, `src/fault.rs`). The production crate is
untouched; a fidelity test asserts the no-fault mirror is byte-identical to the real
`sign_internal`. Success criterion: escape of an exploitable correct/faulty artifact,
**not** key recovery (cited: Bruinderink–Pessl).

## Fault effect

A low-bit `cs1` flip changes one `z` coefficient by `2^3 = 8`, which stays within the
`γ1 − β` norm bound so the signer still accepts, but propagates through verification as
`A·Δ` (a full-magnitude change across the polynomial), so **the faulty signature does
not verify**. Confirmed by test: faulty ≠ correct, `verify(faulty) = false`,
`verify(correct) = true`.

## The 2×2

| signing mode  | verify-after-sign OFF                                   | verify-after-sign ON |
|---------------|--------------------------------------------------------|----------------------|
| deterministic | **ESCAPE** — faulty sig released; correct/faulty pair on the **same y** (directly exploitable) | **BLOCKED** — withheld |
| hedged        | **ESCAPE** — faulty sig released (pair on differing y)  | **BLOCKED** — withheld |

Per-cell (faulty / verifies / released / escaped):

```
  deterministic  v-a-s=off   faulty=1 verifies=0 released=1 escaped=1
  deterministic  v-a-s=on    faulty=1 verifies=0 released=0 escaped=0
  hedged         v-a-s=off   faulty=1 verifies=0 released=1 escaped=1
  hedged         v-a-s=on    faulty=1 verifies=0 released=0 escaped=0
```

## Reading

- **verify-after-sign closes the escape in every mode.** The faulty signature fails the
  internal verification, so it is withheld — `escaped = 0` in both ON cells. This is the
  countermeasure's measured value: with it off, the faulty artifact always escapes.
- **Deterministic signing is the worse case.** With fixed `rnd`, re-signing the same
  message reproduces the *correct* signature on the same `y`, handing the attacker the
  directly-comparable correct/faulty pair that the differential attack needs. Hedged
  signing still leaks the faulty artifact (v-a-s off) but on a different `y`, so the pair
  is not directly comparable — weaker for the attacker, but not a substitute for
  verify-after-sign.
- Net: the recommended configuration is **verify-after-sign on**; it neutralizes the
  fault escape regardless of deterministic vs hedged.

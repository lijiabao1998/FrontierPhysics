# PHYS001_INDEPENDENT_AUDIT (DeepSeek r1) — REVISION 2

Supersedes revision 1. Every claim that changed is marked **RETRACTED** or **CORRECTED**
with the evidence that forced the change. A verifier's own retraction is part of the
record, not something to hide.

Question under audit: is the PHYS-001 r1 `E3` failure (β=3 negative control frozen at
"mean S₂ slope ≥ 1.7", measured **1.54**) an **(A)** evaluator-threshold error,
**(B)** implementation error, **(C)** wrong physical reasoning, **(D)** finite-window
effect, or **(E)** insufficient evidence?

## 0. What revision 1 got wrong, and why it matters

| revision 1 claim | status | evidence |
|---|---|---|
| "independent value 1.4911 vs GLM's 1.54, 3% apart, **consistent with realisation noise**" | **RETRACTED** | The estimator has **zero** realisation noise: GLM's own `sf_std = 7.1e-15`, and my independent replication gives identical slopes for all 5 seeds. Noise cannot explain a 0.049 gap. Revision 1 was comparing **two different quantities**: its own number was the slope over all 234 integers in [16,249], while GLM's is the slope over GLM's 13 geometric lags. The gap is a **lag-set (sampling) difference**. |
| "implementation error (B) **excluded**" | verdict unchanged, **justification replaced** | Revision 1's justification (a tolerance) was invalid, as the reviewer said. (B) is now excluded by **exact reproduction on GLM's own terms** — see §1. |
| "the window's upper lags are **already saturated**; a second **saturation** mechanism truncates the window" | **RETRACTED** as stated | Saturation needs `k_lo·r >> 1`; the window only reaches `k_lo·r = 1.146`. The exact local exponent at `r=217` is **1.1630** and at `r=249` is **1.0723**, not ≈0. Replaced by: the window extends past the **onset of finite-band turnover** (`k_lo·r ~ 1`), so the local exponent is falling through the window. |
| "C₀ = 0.461390, stable to 5e-7" | **RETRACTED, twice.** | The script printed `0.461336` while the report used `0.461390`. The first revision then replaced both with `0.4613932125`, which **was itself wrong in the 6th decimal** — it came from `mpmath.quad(f, [1, inf])`, which mishandles the oscillatory integral and gave `F_tail = 0.481883428` where two independent routes give `0.481882378`. The validated value is **`0.4613921675492818`**, computed with the standard library only. See §4. |
| "recommended window `[4, 20]`" | **RETRACTED** — and so is its replacement | `[4,20]` violated the file's own gate (`k_hi·4 = 5.52 < 20`). The gate-satisfying window `[15,21]` was then tested and **did not achieve the gate's purpose**: it gives β=5/3 → **0.7285**, a *larger* deviation from 2/3 than GLM's window. See §3. |

## 1. (B) implementation error: excluded by exact reproduction

`reproduce_and_compare.py` runs GLM's script **verbatim** (copy hash
`b1070ea49ec142addeda3be51092eaa92624199bcf788c78013cad0c172d5a70` over LF-normalised bytes, matching the file on
`glm/PHYS-001-baseline-r1`) and independently reimplements its specification.

Derivation used for the closed form: with `u_j = Σ_k A_k cos(2πkj/N + φ_k)` and
`A_k = k^(-β/2)`, the `j`-average of the squared increment is **phase-independent**
(the cross terms vanish over a full period), so exactly

```
s2(r) = Σ_{k=3..900} k^(-β) (1 − cos(2πk r / N))        for EVERY seed
```

Measured vs closed form: agreement to `≈5e-15` (machine precision), confirming the
derivation. Independent reimplementation vs GLM's published results:

| quantity | GLM | independent | gap |
|---|---|---|---|
| β=5/3 slope | 0.6383646079612384 | 0.6383646079612274 | **1.1e-14** |
| β=3 slope | 1.5400285073651272 | 1.5400285073651140 | **1.3e-14** |
| β=5/3 per-seed std | 7.1e-15 | 0 (to ~1e-15) | — |

**Error budget: the comparison is limited only by double-precision arithmetic (~1e-14).**
There is no tolerance, no seed ensemble and no unexplained residual. (B) is excluded.

Note the claim's **scope**: this establishes that GLM's *estimator implementation* matches
its documented specification. It is not a claim about turbulence.

## 2. (C)/(D): the exact local exponent, and where the window actually sits

`d ln S2 / d ln r` of the **exact discrete mode sum** (symmetric difference in `ln r`,
`h = 0.01`; no fitting, no window):

| r | k_lo·r | β = 5/3 | β = 2 | β = 3 |
|---|---|---|---|---|
| 16 | 0.074 | 0.7993 | 1.0396 | **1.7324** |
| 32 | 0.147 | 0.6528 | 0.9596 | 1.6687 |
| 64 | 0.295 | 0.6316 | 0.9166 | 1.5717 |
| 128 | 0.589 | 0.5738 | 0.8206 | 1.4014 |
| 217 | 0.999 | 0.4615 | 0.6634 | **1.1630** |
| 249 | 1.146 | 0.4129 | 0.6001 | **1.0723** |

The reviewer's independently computed values (~1.16 at `r=217`, ~1.07 at `r=249`) are
reproduced exactly. Two consequences:

1. **Saturation is not reached, and the plateau language is retracted.** The exponent is
   *decreasing* through the window but is nowhere near 0.
2. **The frozen threshold 1.7 is the window's lower-edge local exponent.** At `r=16` the
   exact value is **1.7324** (the log approximation gives 1.7168; it is off by 0.016
   because `k_lo r = 0.074` is not yet small). So 1.7 was set at a value attainable only
   at the window's first lag and unattainable by any average over it. That is the
   proximate cause of the recorded FAIL.

**Verdict unchanged, mechanism corrected: (D) finite-window effect, with (A) as the
proximate cause. Not (B).** (C) is not wrong so much as beside the point: the β=3 `r²log`
marginal behaviour is real (§4), but the window's value is governed by the falling local
exponent, of which the log correction is one ingredient and the finite-band turnover the
other.

One revision-1 claim survives with a corrected mechanism: lowering `k_lo` (wider band)
does move the slope toward 2 — because it raises `r_* = 1/k_lo` so the window sits
further inside the log regime, **not** because a plateau was involved.

## 3. (A)/window: the proposed gate is retracted

Revision 1 proposed: a slope measures `β−1` only if `k_hi·r_a ≥ 20` and
`k_lo·r_b ≤ 0.1`. For GLM's band that requires `r_a ≥ 14.49`, `r_b ≤ 21.73` → the
integers `[15,21]`. Measured over that window:

| β | slope over [15,21] | theory | gap |
|---|---|---|---|
| 5/3 | 0.7285 | 0.6667 | **+0.0618** |
| 2 | 1.0064 | 1.0 | +0.0064 |
| 3 | 1.7233 | 2.0 | −0.2767 |

Window-perturbation spread is small (2.3e-3 … 4.5e-3), so the window is *stable* — but it
is **biased the other way**. The reason is now visible in §2: the local exponent is
**above** 2/3 at small `r` and below it at large `r`. Which side a window errs on depends
on where it sits, not on how well resolved the band is. A "resolved-band" condition alone
therefore does **not** make an OLS slope equal the asymptotic exponent. The gate is
retracted as a recommendation.

Chord-slope scan (all integers, so this is the analysis's own sampling; reported to show
the dependence, **not** to select a window):

| window | β=5/3 | β=2 | β=3 |
|---|---|---|---|
| [16,249] (GLM's span) | 0.6169 | 0.8725 | 1.4911 |
| [15,21] | 0.7285 | 1.0064 | 1.7233 |
| [16,64] | 0.6897 | 0.9690 | 1.6590 |
| [16,128] | 0.6624 | 0.9339 | 1.5942 |
| [32,249] | 0.5953 | 0.8438 | 1.4418 |
| [64,249] | 0.5627 | 0.7993 | 1.3688 |
| [16,512] | 0.5089 | 0.7295 | 1.2750 |
| [16,1024] | 0.2663 | 0.4038 | 0.7794 |

β=5/3 spans **0.27 … 0.73** — a spread of 0.46, i.e. ~10× the tolerance used for it.
**The measured slope is a property of the (band, window, lag-set) triple, not of the
spectrum.** Two further consequences:

- **Lag set matters too.** Over the same span, GLM's 13 geometric lags give 0.6384 while
  all 234 integers give 0.6169. Both are legitimate; neither is "the" exponent.
- **β=3 requires both cutoffs to be checked.** The marginal approximation
  `2 − 1/(ln(1/(k_lo r)) + 2C₀)` requires `k_lo·r ≪ 1` **and** `k_hi·r ≫ 1`.
  The explicit gates used here (`k_lo·r ≤ 0.1`, `k_hi·r ≥ 20`) give the integer
  window **[15,21]**, where both conditions hold and the OLS slope is **1.7233**.
  Requiring the approximation to exceed 1.9 gives **r < 0.0248**, but there
  `k_hi·r < 0.035`: the UV condition fails, so that approximation cannot be used.
  The exact finite-sum local exponent actually exceeds 1.9 at **r=1 (1.9743)**
  and **r=2 (1.9117)**. Those small lags approach the analytic, upper-cutoff
  `r²` regime; they do not establish marginal inertial-range scaling near 2.
  **RETRACTED:** the earlier explanation that the sole obstruction was lattice
  spacing, and its implication that no integer lag can exceed 1.9. See
  `CANONICAL_FACTS.md` P3 for the matching UV/IR table. β=3 remains a
  **reproducibility control with a closed form** in this audit.

**What should replace the retracted gate** (pre-registrable, and honest about what it can
and cannot do):
1. Report every slope together with its band, window and lag set. A window-free exponent
   is not available in this configuration.
2. Pre-register the *expected* value from the closed form and the local-exponent profile,
   not from the asymptotic exponent.
3. For β=3, pre-register the closed form and label the control as reproducibility-only.
4. A future marginal-scaling test near 2 would need a separately frozen band and
   window that satisfy both cutoffs. The present [15,21] window instead gives ≈1.72;
   the near-2 values at r=1,2 come from the analytic UV regime.

## 4. The β=3 log constant, with a supported precision

`F(x) = ∫_x^∞ q^{-3}(1−cos q) dq = ½·ln(1/x) + C₀ + O(x²)`, and by term-by-term
integration of the series for `1−cos q`,

```
C0 = F_tail + Σ_{n≥2} (−1)^(n+1) / ((2n)!·(2n−2)),   F_tail = ∫_1^∞ q^-3 (1−cos q) dq
```

computed with the **standard library only** — an earlier revision of this section said "mpmath at
40 digits", which was false: the committed script imports no `mpmath`, performs ordinary
double-precision Simpson integration, and ignores its `--dps` argument. The multiprecision check
was run interactively during development and is **not** reproduced by any committed command, so
the section no longer attributes the result to it.

- **C₀ = 0.4613921675492818** — standard library only (no third-party module, double precision).
  Computed as `F_tail` plus the rapidly converging series, where `F_tail` is a composite
  Simpson on `[1, 1000]` plus its analytic remainder `1/(2Q²) + Q⁻³sinQ − 3Q⁻⁴cosQ`.
  **`F_tail = 0.481882378019`**, validated to ~1e-12 by two independent routes: Simpson plus
  analytic tail at `Q = 1000, 2000, 5000` (stable), and `F_tail = ½ − A₃(1)` with
  `A₃(1) = ∫₁^∞ cos q / q³ dq` on a different grid to `2π·400` plus its own remainder.
  This **replaces** the earlier `0.4613932125491026`, which was wrong in the 6th decimal.
- Convergence sweep over the advertised range, `x = 1e-1 … 1e-20` (20 points):
  `F(x) − ½ln(1/x)` runs `0.4616015…` → `0.4613932…`, **max − min = 2.1e-4**.

So the value is supported to **≈2e-4** by this sweep, not to 5e-7. The sweep is the one
revision 1 advertised but never performed (it evaluated only `x = 1e-10`). All three earlier
numbers — `0.461336`, `0.461390`, and `0.4613932125` — are superseded by `0.4613921675492818`.

**Dependency note.** The first revision of this section used `mpmath` for a "40-digit"
quadrature, which was declared nowhere in the repository, so the advertised reproduction
command could not run on a clean checkout. The dependency is **removed**, not declared: the
computation is standard library only. Reviewers should treat the surviving docstring mentions
of `mpmath` as history, not as an import.

## 5. Verdict for the recorded PR (unchanged in substance, corrected in statement)

**(D) finite-window effect, with (A) as the proximate cause of the recorded FAIL. Not (B)
— now excluded by exact reproduction to ~1e-14. (C) incomplete rather than wrong.**

- **Keep the FAIL.** It is the clearest evidence in the repo that the frozen gate has
  teeth, and the number is correct for the configuration measured.
- **Correct the explanation** to: the local exponent of the exact mode sum falls from
  1.7324 to 1.0723 across the window (the window crosses the finite-band turnover onset
  `k_lo r ~ 1`); the frozen 1.7 is the lower-edge local value; and the β=3 `r²log`
  marginal form is real but is not by itself the reason the average is 1.54.
- **Do not** present the measured slope as an exponent without its window and lag set.

## 6. Residual limitations of this revision

- Still no claim about real Navier–Stokes turbulence: everything here concerns an
  estimator applied to a synthetic spectrum.
- GLM's r2 branch code was not executed; its "90 measurements, 0 violations" remains a
  documented claim. Its β=3 centre value is confirmed correct by §1.
- The reproduced agreement in §1 tests the *slope*. GLM's other estimators (the FFT
  periodogram slope, the energy budget) were not independently reimplemented; the energy
  residual being an analytic identity at 4.1e-14 was accepted as reported.
- The chord-slope scan in §3 uses all integers per window; a geometric lag set with the
  same endpoints would give slightly different numbers. This affects no conclusion: the
  spread across windows dwarfs the lag-set effect.

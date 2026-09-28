# External verifier review — FrontierPhysics PRs #1, #2 (DeepSeek r1, 2026-09-28)

> ## ⚠ SUPERSEDED IN PART — read this first
>
> This review was written before this PR's own revision 2 audit. Two of its statements
> are **no longer accurate** and must not be acted on:
>
> 1. It describes the β=3 window as straddling a **saturation** edge and recommends a
>    change on that basis. Revision 2 **retracted** the saturation diagnosis: `k_lo·r`
>    reaches only 1.146, not `>> 1`, and the exact local exponent is **1.1630** at
>    r=217 and **1.0723** at r=249 — nowhere near the 0 of a saturated plateau. The
>    correct description is that the window crosses the **onset of finite-band
>    turnover**. See `problems/PHYS-001/experiments/dsk_audit_r1/PHYS001_INDEPENDENT_AUDIT.md`
>    §2.
> 2. It quotes `F(x) = ½·ln(1/x) + 0.461390…`. That constant is superseded twice over:
>    the validated value is **0.4613921675492818** (see the audit §4), and the earlier
>    0.4613932125 from mpmath was itself wrong in the 6th decimal.
>
> The **verdicts** in this review are unchanged and still stand: the FAIL is kept, PR #1
> and PR #2 both `NEEDS_CHANGES` on their explanations rather than on any result. What is
> withdrawn is the saturation *mechanism* and the quoted constant. The authoritative
> document is `PHYS001_INDEPENDENT_AUDIT.md` (revision 2); this file is kept as the record
> of the first review, not as current guidance.



Reviewer: DeepSeek, branch `dsk/PHYS-001-independent-audit-r1`. No reviewed branch
modified; no merge performed. Independent implementation:
`problems/PHYS-001/experiments/dsk_audit_r1/s2_scaling_audit.py` (direct summation
over modes, no FFT, no estimator, no GLM code). Full analysis:
`dsk_audit_r1/PHYS001_INDEPENDENT_AUDIT.md`.

## PR #1 — `glm/PHYS-001-baseline-r1`

| CLAIM | EVIDENCE | INDEPENDENTLY VERIFIED? | SCOPE | FAILURE MODE | REC |
|---|---|---|---|---|---|
| E1 energy budget residual 4.1e-14 | `results/r1/k41_baseline_results.json` | not re-run (would require GLM's field generator); the claim is a machine-precision identity for periodic modes and is *prima facie* plausible | synthetic periodic field | — | — |
| E2 β=5/3 S₂ slope 0.6384 in [0.6167, 0.7167] | same | **yes in substance** — my independent value is **0.6169** over the same window, inside GLM's band and on the same side of 2/3 | synthetic spectrum | **the band edge already biases this control**; GLM passed only because the tolerance was ±0.05 | — |
| E4 spectral slopes 1.667 / 3.0008 (±0.1) | same | not re-run; consistent with the analytic slopes of `A_k = k^(-β/2)` | — | — | — |
| **E3 β=3 negative control FAIL: 1.54 < frozen 1.7** | same + r1 report analysis | **value verified: 1.4911** vs 1.54 (3% apart, consistent with realisation noise) | synthetic spectrum | see below | — |
| r1 diagnosis "β=3 邊際情形，S2 ∝ r²·log 之緩慢變化有效指數" | report | **correct** — verified analytically in this audit: `F(x) = ∫_x^∞ q^-3(1-cos q)dq = ½·ln(1/x) + 0.461390…` | — | **incomplete**: omits the saturation edge | — |
| two self-caught bugs (log-linear spectral regression; wrong output path) | report | consistent with the recorded before/after values | — | — | — |
| no independent-reproduction claim | PR body | **correctly disclaimed** | — | — | — |

**Verdict: `NEEDS_CHANGES`.** Not because of the FAIL — the task asked precisely
whether a FAIL should be rejected, and it should **not** be: the failure is honestly
recorded, the frozen rule was respected, and the number is correct. The change
required is to the *diagnosis*: r1 attributes 1.54 to the r²·log marginal exponent
alone, but the quantitative reason the *window average* cannot reach ~1.72 is that
the window straddles the saturation cross-over `r_* = 1/k_lo = 217` (the window's
`k_lo·r` runs 0.074 → 1.146). The correct statement is "both the log correction *and*
the saturation edge", and it changes what the next round should do.

## PR #2 — `glm/PHYS-001-calib-r2`

| CLAIM | EVIDENCE | INDEPENDENTLY VERIFIED? | SCOPE | FAILURE MODE | REC |
|---|---|---|---|---|---|
| two-stage process: cycle 1 FAIL (zero-width bands), then fix, then cycle 2 PASS | `experiments/r2_calibration/{bands.json,bands2.json}`, `results/r2/*` | **not reproduced** (code not executed here); the timestamps and the recorded cycle-1 failure are consistent and the failure mode is real and correctly diagnosed | calibration methodology | — | — |
| cycle 2: 90 measurements, 0 violations, bands = mean ± max(4σ, 0.01) | `bands2.json` | centre values checked: β=3 centre ≈ 1.54 matches my independent 1.4911; separation 0.64 / ~1.0 / ~1.54 is real | — | — | — |
| this fixes the r1 calibration error | report | **partially.** It makes the control *reproducible*; it does not make it a *scaling-law* test | — | see below | — |

**Verdict: `NEEDS_CHANGES`.** The zero-width-band failure and the
`±max(4σ, floor)` repair are a genuine methodological improvement and should be kept.
But freezing the window average of a quantity that has no constant exponent produces a
stable number that is not the exponent. Recommended one-line addition: state that for
β=3 the measured slope is a window property, not an exponent, and either (a) narrow
the window so that `k_hi·r_a ≥ 20` and `k_lo·r_b ≤ 0.1`, or (b) drop β=3 as a
power-law control and use the closed form
`2 − 1/(ln(1/(k_lo r)) + 0.9228)` as its expectation. Rule (a) is checkable before any
run, so it is a real pre-registration gate.

## On the assigned question: is the FAIL itself a research asset?

**Yes, and it should be kept and merged.** The instruction to treat a failure round as
an asset applies exactly here:

- The FAIL was produced by a **pre-registered threshold that was not relaxed
  afterwards**, which is the behaviour the governance is trying to induce. Loosening
  E3 post hoc would have destroyed the only demonstration in the repo that the frozen
  gate has teeth.
- The FAIL is traceable to a real, subtle, *reproducible* property of the β=3
  spectrum — a genuine marginal case whose integral diverges logarithmically at the
  infrared end. That is a better outcome than a passing control.
- Everything needed to judge it is present: field spec, band, window, seeds,
  thresholds, the raw per-seed JSON, environment and hashes. The FAIL is fully,
  accurately and reproducibly recorded.
- The one thing a reader cannot get from PR #1 alone is the closed form for β=3. That
  is the change requested above.

## Merge recommendations

| PR | recommendation |
|---|---|
| FrontierPhysics #1 | **NEEDS_CHANGES** — keep the FAIL as recorded; replace the β=3 explanation with the corrected two-mechanism account; note that the β=5/3 control is also band-edge biased (passed only because the tolerance was ±0.05) |
| FrontierPhysics #2 | **NEEDS_CHANGES** — keep the two-stage method and the floor rule; add the scaling-window validity gate and stop describing the β=3 centre as an exponent |

Neither PR is recommended `READY_FOR_OWNER_REVIEW` as-is. Neither should be rejected
for the FAIL: the r1 FAIL and the r2 zero-width-band FAIL are the strongest evidence
in either physics PR that the frozen-gate discipline is real rather than decorative.

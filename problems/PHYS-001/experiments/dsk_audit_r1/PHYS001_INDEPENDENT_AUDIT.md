# PHYS001_INDEPENDENT_AUDIT (DeepSeek r1)

Independent audit of the PHYS-001 r1 `E3` failure: the β=3 negative control was frozen
at "mean S₂ slope ≥ 1.7" and measured **1.54**, recorded as FAIL. Question posed to
this audit: is that failure **(A)** a wrong evaluator threshold, **(B)** an
implementation error, **(C)** wrong physical reasoning, **(D)** a finite-window
effect, or **(E)** insufficient evidence?

**Answer: D, with A as the proximate cause of the recorded FAIL. Not B. C is
incomplete rather than wrong.** Details and the arithmetic below.

No GLM code was read or reused. The implementation in
`experiments/dsk_audit_r1/s2_scaling_audit.py` evaluates the *definition*

```
S2(r) = 2 * Σ_k P(k) (1 - cos(k r)),      P(k) = k^-β
```

by direct summation over modes — no FFT, no structure-function estimator, no shared
code. A log-log least-squares slope (not log-linear; GLM's own r1 round already
caught that class of bug in its own `spectral_slope`) is fitted over GLM's frozen
window [16, 249] with GLM's band (N = 4096, mode indices j = 3..900).

## 1. Derivation, including the marginal case β = 3

With q = k·r,

```
S2(r) = 2 C r^(β-1) · J(k_lo r, k_hi r),      J(a,b) = ∫_a^b q^-β (1 - cos q) dq
```

- **1 < β < 3.** J → I(β) = ∫_0^∞ q^-β(1-cos q)dq, finite. So `S2(r) ~ r^(β-1)`,
  and β = 5/3 gives the 2/3 law. This is the ordinary case.
- **β = 3 is genuinely marginal, and the integral diverges logarithmically.** Near
  q = 0, `q^-3 (1 - cos q) = 1/(2q) - q/24 + ...`, so the infrared end diverges.
  Integrating the series term by term:

  ```
  F(x) := ∫_x^∞ q^-3 (1 - cos q) dq = (1/2)·ln(1/x) + C0 + O(x^2),
  C0 = 0.461390...      (numerically: 0.461390, stable to 5e-7 for x = 1e-1 … 1e-20)
  ```

  hence

  ```
  S2(r) = 2 C r^2 [ (1/2)·ln(1/(k_lo r)) + C0 ]        for k_lo r << 1
  ```

  So β = 3 gives **r² times a log**, and the local exponent is

  ```
  d ln S2 / d ln r = 2 - 1 / ( ln(1/(k_lo r)) + 2·C0 )
  ```

  strictly below 2, and **r-dependent**. This is exactly the "r²·log" that GLM's r1
  round named. It is correct, and this audit confirms it analytically.

- **Second, independent truncation.** Once k_lo r ≫ 1 the whole band is saturated:
  J → const and S2(r) flattens to `k_lo^-2 - k_hi^-2`. The cross-over lag is
  `r_* ≈ 1/k_lo`. For GLM's band, k_lo = 2π·3/4096 = 0.004597, so **r_* = 217.3**.

## 2. The decisive number

GLM's window is [16, 249]. Across it, `k_lo·r` runs from **0.0736 to 1.146** — i.e.
the window **straddles r_* = 217.3**. The log regime requires k_lo r ≪ 1 and only
just holds at the window's bottom edge; the top edge is already saturated.

| quantity | value |
|---|---|
| GLM reported β=3 S₂ slope | 1.54 |
| this audit's independent β=3 slope over [16,249] | **1.4911** |
| log-formula local exponent at r = 16 | 1.7168 |
| log-formula local exponent at r = 32 / 64 / 128 / 208 | 1.6477 / 1.5338 / 1.3113 / 0.9652 |
| GLM's frozen r1 threshold | 1.7 |

**1.4911 vs 1.54 agree to 3%.** GLM's measurement is reproduced by a code-independent
computation. (B) is excluded.

The local exponent from the closed form declines monotonically across the window and
its average is ≈ 1.49 — matching the measured window slope. So the 1.54 is not a
wrong exponent: it is the average of a function that has no constant exponent.

Notice also that 1.7168 — the log-formula value at the window's *first* lag — is
essentially GLM's frozen 1.7. The threshold was set at a value attainable only at the
window edge, and therefore unattainable by any *average* over the window. That is the
proximate cause of the FAIL (A), but the reason it is unattainable is a window effect
(D), so the fix is not "lower the threshold" alone.

## 3. Controls that make the diagnosis non-vacuous

| probe | result | reading |
|---|---|---|
| band widened, N = 4096 (j 1..2048) | slope 1.7137 | moves toward the β−1 = 2 limit |
| N = 16384 (j 1..4096) | slope 1.7961 | further toward 2 |
| N = 65536 (j 1..8192) | slope 1.8416 | monotone convergence to 2 |
| β = 2, GLM band | slope 0.8725 (asymptotic 1.0) | same machinery reproduces r¹ direction |
| β = 5/3, GLM band | slope 0.6169 (first half 0.6608) | reproduces the 2/3 law |
| β = 3, low-k edge raised (j_min = 20) | slope −0.0555 | slope is entirely a band-edge artefact |

The β=3 slope collapses to ≈ 0 as soon as the low-k edge moves up. That is the
signature of a saturation effect, not of a power law: with j_min = 20, k_lo r ≫ 1
across the whole window and S₂ is a plateau.

And the β = 5/3 case is instructive in the same direction: its slope is also
slightly below 2/3 (0.6169 vs 0.6667) because its window, too, reaches the
saturation edge. GLM's β = 5/3 pass only survived because the frozen tolerance was
±0.05 around 2/3, and 0.6384 (GLM) / 0.6169 (this audit) are both inside it. The
same band-edge effect is present in the *passing* control; it is simply small enough
there not to matter.

## 4. Verdict and the actual design flaw

- **(B) implementation error — excluded.** Independent definition-level computation
  reproduces the number.
- **(C) physical reasoning — not the failure.** The r²·log reasoning is correct and
  this audit verified the log divergence analytically. What the r1 explanation omits
  is the *second* edge (saturation at r_* = 1/k_lo), which is what actually makes the
  window average fall to 1.49 rather than sit near 1.72.
- **(A) threshold wrong — true but proximate.** 1.7 was set at the window-edge local
  exponent.
- **(D) finite-window effect — the substantive answer.** The window straddles r_*.
- **(E) insufficient evidence — not needed** for this question; the number is
  reproducible and its mechanism is closed-form.

**The deeper design flaw, which neither r1 nor r2 fixed:** β = 3 has *no* clean
scaling exponent in a finite band. The local exponent is a continuous, r-dependent
quantity running from 2 (small r) through the log regime to 0 (saturated). A control
that asserts a single frozen exponent for β = 3 is asserting a power law that does not
exist in the configuration being measured. GLM's r2 "two-stage calibration" is a
genuine methodological improvement and its ±max(4σ, floor) band rule is sound, but it
*calibrates around* the artefact rather than removing it: it freezes whatever the
window average happens to be, which makes the control reproducible without making it
a scaling-law test.

## 5. Concrete pre-registration rule this audit recommends

For a synthetic spectrum P(k) ∝ k^-β on a band [k_lo, k_hi] and a lag window [r_a, r_b],
a slope measurement is a measurement of the exponent β−1 **only if**

```
k_hi · r_a  ≥  20      (all band modes resolved at the smallest lag)
k_lo · r_b  ≤  0.1     (the largest lag stays well inside the log/inertial regime)
```

Under GLM's band this requires `r_b ≤ 0.1/0.004597 ≈ 21.8`, i.e. a window like
[4, 20] — much shorter than [16, 249]. With that rule, β = 5/3 and β = 2 would be
measurable cleanly and β = 3 would be *reported as having no clean exponent*, with
the closed form `2 − 1/(ln(1/(k_lo r)) + 0.9228)` quoted instead of a single number.
Both requirements can be checked before any run, so this is a real pre-registration
gate rather than a post-hoc adjustment.

## 6. Residual limitations of this audit

- The audit models the synthetic field as an exact continuum/mode sum. GLM's field is
  a finite random-phase superposition with amplitudes `A_k = k^(-β/2)`, i.e. the same
  spectrum `P(k) = k^-β`, so the comparison is like-for-like up to realisation noise;
  the 3% gap between 1.4911 and 1.54 is consistent with that and is not a discrepancy
  this audit can fully attribute without GLM's seeds.
- GLM's r2 branch (`glm/PHYS-001-calib-r2`) was read only as a report; its code was
  not executed here. The r2 claim of "90 measurements, 0 violations" is therefore
  taken as a documented claim, not independently reproduced. What *is* independently
  established is that its centre value (≈1.54 for β = 3) is the correct value for that
  band and window, so r2 is not hiding a numerical error.
- Nothing here bears on real Navier–Stokes turbulence. This is a statement about an
  estimator applied to a synthetic spectrum, which is what the control was.

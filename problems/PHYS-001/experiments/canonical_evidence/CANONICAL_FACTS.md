# PHYS-001 — canonical evidence

The canonical merge unit for the PHYS-001 E3 audit. It carries the estimator reproduction, the
frozen FAIL, the corrected interpretation, an explicit retraction list, and the minimal scripts
that regenerate the numbers — and nothing else.

**No claim about real three-dimensional Navier–Stokes turbulence is made anywhere.** Everything
below concerns an estimator applied to a synthetic spectrum, which is what the control was.

Reproduce with:
```bash
cd scripts
python reproduce_and_compare.py --json /tmp/glm.json     # gate on the pinned reference hash
python s2_local_exponent.py   --json /tmp/s2.json
```

---

## P1 — GLM's estimator is reproduced exactly, so an implementation error is excluded

**Claim.** An independent implementation of GLM's specification reproduces their
structure-function slope to ordinary double-precision rounding, on GLM's own 13 geometric lags,
against GLM's **full-precision** stored values.

| quantity | GLM (stored) | independent | gap |
|---|---|---|---|
| β = 5/3 slope | 0.6383646079612384 | 0.6383646079612274 | **1.1e-14** |
| β = 3 slope | 1.5400285073651272 | 1.5400285073651140 | **1.3e-14** |
| β = 5/3 per-seed spread | 7.1e-15 | ~1e-15 | — |

**How the exactness arises.** With `u_j = Σ_k A_k cos(2πkj/N + φ_k)` and `A_k = k^(−β/2)`, the
`j`-average of the squared increment is **phase-independent** (cross terms vanish over a full
period), so

```
s2(r) = Σ_{k=3..900} k^(−β) (1 − cos(2πk r / N))        for EVERY seed
```

This is why the per-seed spread is zero and why the comparison is limited only by floating-point
arithmetic. **There is no tolerance in this claim.** The comparison bound is a floating-point
error budget of 1e-12 against the full-precision values; an earlier 5e-5 bound rested on a false
"GLM rounds to four decimals" premise and is retracted.

**The reference is not trusted, it is verified.** GLM's script is executed only after its
LF-normalised SHA-256 equals the pinned
`b1070ea49ec142addeda3be51092eaa92624199bcf788c78013cad0c172d5a70` — which equals both the blob
committed here and the file on `glm/PHYS-001-baseline-r1`. A mismatch **refuses execution**
(verified by tampering). Any pre-existing output is deleted before launch, so a stale artifact
cannot be mistaken for current evidence; the reference's exit code 1 is its scientific FAIL and
is deliberately not treated as a failure of this check.

**Evidence:** `evidence/glm_reproduction.json`. **Script:** `scripts/reproduce_and_compare.py`.

## P2 — The exact local exponent, and where the frozen window actually sits

`d ln S2 / d ln r` of the exact discrete mode sum (symmetric difference in `ln r`, `h = 0.01`;
no fitting, no window):

| r | k_lo·r | β = 5/3 | β = 2 | β = 3 |
|---|---|---|---|---|
| 16 | 0.074 | 0.7993 | 1.0396 | **1.7324** |
| 32 | 0.147 | 0.6528 | 0.9596 | 1.6687 |
| 64 | 0.295 | 0.6316 | 0.9166 | 1.5717 |
| 128 | 0.589 | 0.5738 | 0.8206 | 1.4014 |
| 217 | 0.999 | 0.4615 | 0.6634 | **1.1630** |
| 249 | 1.146 | 0.4129 | 0.6001 | **1.0723** |

**Two consequences.**
1. **Asymptotic saturation is NOT reached.** It requires `k_lo·r >> 1`; the window reaches 1.146.
   The exponent is falling through the window but is nowhere near the 0 of a plateau. The
   "saturated plateau" and "second saturation mechanism" language is **retracted** and replaced
   by: the window crosses the **onset of finite-band turnover**.
2. **The frozen threshold 1.7 is the window's lower-edge local exponent** (1.7324; the log
   approximation gives 1.7168). It was set at a value attainable only at the window's first lag
   and unattainable by any *average* over the window — the proximate cause of the recorded FAIL.

**Evidence:** `evidence/s2_local_exponent.json`. **Script:** `scripts/s2_local_exponent.py`.

## P3 — The β = 3 statistic is a finite-band / window property, not a clean exponent

The measured slope is a property of the **(band, window, lag set)** triple. Chord slope over
all integers in each window:

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

β = 5/3 alone spans **0.27 … 0.73** across windows — a spread of 0.46, about ten times the
tolerance used for it. The **lag set matters too**: over the same span, GLM's 13 geometric lags
give 0.6384 while all 234 integers give 0.6169. Both are legitimate; neither is "the" exponent.

**β = 3 can never be a power-law test at this `N` and band.** Its local exponent is
`2 − 1/(ln(1/(k_lo·r)) + 2C₀)`, which exceeds 1.9 only for `k_lo·r < 1.14e-4`, i.e.
**`r < 0.0248`** — the required regime is at *small* `r`, **below one sample**. (An earlier
revision had this inequality reversed, claiming `r ≳ 3300`; corrected.)

## P4 — The β = 3 log constant, with a supported precision

```
F(x) = ∫_x^∞ q^(−3)(1 − cos q) dq = ½·ln(1/x) + C₀ + O(x²)
C₀ = F_tail + Σ_{n≥2} (−1)^(n+1) / ((2n)!·(2n−2))
```

- **C₀ = 0.4613921675492818**, standard library only.
- **`F_tail = 0.481882378019`**, validated to ~1e-12 by two independently executed routes:
  (A) Simpson on `[1,Q]` plus the analytic remainder at `Q = 1000/2000/5000` →
  `0.4818823780193526 / 0.48188237801925576 / 0.4818823780182286`; (B) `F_tail = ½ − A₃(1)`
  with `A₃(1) = ∫₁^∞ cos q / q³ dq` on a different grid → `0.48188237801937406`. Spread 1.1e-12
  against an acceptance bound of **1e-11 that fails the script if exceeded**.
- The sweep over the advertised range `x = 1e-1 … 1e-20` (20 points) has a spread of **2.1e-4**,
  so the value is supported to that, not to 5e-7.

**Three earlier values are retracted:** `0.461336` (partially converged), `0.461390` (the report),
and `0.4613932125491026` (an `mpmath.quad(f,[1,∞])` that mishandled the oscillatory integral —
`F_tail = 0.481883428` instead of `0.481882378`, an error of 1.05e-6, presented at 16 digits on
top of an undeclared dependency).

## P5 — Verdict on the frozen E3 FAIL

**(D) finite-window effect, with (A) as the proximate cause. Not (B) — excluded by exact
reproduction to ~1e-14. (C) incomplete rather than wrong:** the β=3 `r²log` marginal behaviour is
real and literature-corroborated (the 2D enstrophy cascade's explicit `k⁻³` log corrections;
arXiv:2308.06997 on the marginal spectral case), but it is not by itself why the window average
is 1.54.

**The FAIL is kept.** A pre-registered threshold that was not relaxed afterwards is exactly what
the governance is trying to induce, and it is the clearest evidence in the repository that the
frozen gate has teeth.

## P6 — What is NOT claimed

- No claim about real Navier–Stokes turbulence of any kind.
- GLM's r2 branch code was **not executed**; its "90 measurements, 0 violations" remains a
  documented claim, not an independent reproduction. What *is* established here is that its β=3
  centre value is the correct value for that band and window.
- GLM's FFT periodogram and energy-budget estimators were **not** independently reimplemented.
- No `S3`/4-5-law check, no shell model, no DNS or experimental data.
- `PHYS-001` remains **OPEN**.

## What is not here, and where it is

| excluded | location |
|---|---|
| artifact-hash automation, rerun-evidence generator, stale-file guards, convergence tooling | `dsk/PHYS-001-audit-tooling` |
| the exploratory tree as originally written | `dsk/PHYS-001-independent-audit-r1` (PR #3) |
| the round's full artifact manifest and command log | the tooling branch and `runs/` on the exploratory branch |

#!/usr/bin/env python3
"""PHYS-001 independent audit (DeepSeek r1): where does the beta=3 S2 slope come from?

=============================================================================
SUPERSEDED -- REVISION 1 OF THE AUDIT. THREE CLAIMS MADE HERE ARE RETRACTED.

This file is kept as the record of what revision 1 computed and got wrong. It is NOT
the current audit. The corrected work lives in:

    reproduce_and_compare.py   -- GLM's estimator reproduced exactly (gap ~1e-14)
    s2_local_exponent.py       -- exact local exponent, window analysis, C0 precision
    PHYS001_INDEPENDENT_AUDIT.md (revision 2) -- what changed and why

Retractions, with the reason, so no reader of this file inherits a wrong claim:

  (1) "implementation error excluded" was argued from a 0.05 tolerance on a 0.0489 gap.
      RETRACTED. The estimator has ZERO realisation noise (GLM's own std is 7.1e-15), so
      no gap can be attributed to noise; and 0.0489 was a LAG-SET difference (GLM's 13
      geometric lags vs all 234 integers in [16,249]). The exclusion now rests on exact
      reproduction to ~1e-14, in reproduce_and_compare.py.

  (2) "the window's upper lags are at or past the saturation cross-over ... the local
      exponent has already dropped to 0" -- RETRACTED. Saturation needs k_lo*r >> 1, and
      the window only reaches k_lo*r = 1.146. The exact local exponent at r=217 is 1.1630
      and at r=249 is 1.0723, nowhere near 0. What actually happens is that the window
      crosses the ONSET of finite-band turnover and the local exponent falls through it.

  (3) "fitted C0 ... constant to 5e-7" -- RETRACTED. This file printed 0.461336 while its
      own report used 0.461390, and it evaluated a single x rather than sweeping the
      advertised 1e-1..1e-20 range. The supported value is 0.4613932125491026 with a
      sweep spread of 2.1e-4, in s2_local_exponent.py.

The verdict this file reached -- (D) finite-window effect with (A) as the proximate
cause, NOT (B) -- survives revision 2, because revision 2 re-derives it on correct
grounds. The diagnosis "the r^2 log marginal form is real" also survives and is now
literature-corroborated. Everything above is left in place, with the retractions
inline below, rather than deleted.
=============================================================================

Question under test: GLM's r1 round froze "beta=3 negative control: mean S2 slope
>= 1.7" and measured 1.54, recording a FAIL. Its r2 round recalibrated the band to
mean +/- max(4 sigma, 0.01) and reported 0 violations over 90 measurements, i.e. it
classified the episode as a *threshold* error. This script decides the question
independently, from the definition, with no GLM code.

Derivation (1-D surrogate, spectrum P(k) = C k^-beta on a band [k_lo, k_hi]):

    S2(r) = <(u(x+r)-u(x))^2> = 2 * integral_{k_lo}^{k_hi} P(k) (1 - cos k r) dk

Substituting q = k r:

    S2(r) = 2 C r^{beta-1} * J(k_lo r, k_hi r),   J(a,b) = integral_a^b q^-beta (1-cos q) dq

* For 1 < beta < 3 and k_lo r -> 0, k_hi r -> infinity, J tends to a constant
  I(beta) = integral_0^inf q^-beta (1-cos q) dq, which is finite. Hence

        S2(r) ~ r^{beta - 1}          (the standard result: beta = 5/3 -> 2/3)

* beta = 3 is the borderline case, and here the naive reading of the integral is
  WRONG in an instructive way. q^-3 (1-cos q) ~ 1/(2q) - q/24 + ... near 0, so
  the integral over q near 0 diverges LOGARITHMICALLY. Writing

      F(x) = integral_x^inf q^-3 (1-cos q) dq,

  term-by-term integration of the series for 1-cos q gives F(x) = (1/2) ln(1/x)
  + C0 + O(x^2) with C0 = 0.461390... (verified to six digits in section E for x
  down to 1e-20). Hence

      S2(r) = 2 C r^2 [ (1/2) ln(1/(k_lo r)) + C0 ]     for k_lo r << 1,

  i.e. r^2 times a log, so the LOCAL exponent is

      d ln S2 / d ln r = 2 - 1 / ( ln(1/(k_lo r)) + 2 C0 )

  which is strictly below 2 and r-dependent. So beta = 3 does NOT have a clean
  power law in a finite band: the local exponent runs from 2 (small r, where
  S2 ~ r^2 analytically) through the log regime and down to 0 (saturated r).

* There is a second, independent truncation. J(a,inf) -> const once k_lo r >> 1,
  so the whole band is saturated and S2(r) flattens to the constant
  k_lo^-2 - k_hi^-2. The cross-over lag is r_* ~ 1/k_lo. A fitting window that
  straddles r_* measures neither the log regime nor the plateau, and yields a
  window-dependent effective slope strictly between the two.

This script (i) evaluates S2(r) by direct summation over modes -- no FFT, no
estimator, no code shared with GLM's implementation -- (ii) reports the effective
slope in GLM's window [16, 249] as a function of the band and of N, and (iii)
checks the log form F(x) = 0.5 ln(1/x) + C0 numerically, so that both mechanisms
are demonstrated rather than asserted.

    python s2_scaling_audit.py --json results/s2_scaling_audit.json
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from typing import Sequence

# GLM r1/r2 frozen configuration, as reported in the round reports.
GLM_N = 4096
GLM_JMIN, GLM_JMAX = 3, 900
GLM_WINDOW = (16, 249)
GLM_SEED_SLOPE_BETA3 = 1.54
GLM_R1_THRESHOLD = 1.7


def modes(n: int, jmin: int, jmax: int) -> list[float]:
    """Physical wavenumbers k_j = 2*pi*j/n for j in [jmin, jmax]."""
    return [2.0 * math.pi * j / n for j in range(jmin, jmax + 1)]


def s2_exact(r_values: Sequence[int], beta: float, ks: Sequence[float]) -> list[float]:
    """S2(r) = 2 * sum_k P(k) (1 - cos(k r)) with P(k) = k^-beta, no normalisation.

    Direct summation: this is the definition, not an estimator. It is O(len(r)*nk),
    which is fine at these sizes and deliberately avoids any transform.
    """
    p = [k ** (-beta) for k in ks]
    out = []
    for r in r_values:
        s = 0.0
        for k, pk in zip(ks, p):
            s += pk * (1.0 - math.cos(k * r))
        out.append(2.0 * s)
    return out


def loglog_slope(r_values: Sequence[int], s2: Sequence[float]) -> float:
    """Least-squares slope of ln S2 against ln r (this is log-log, not log-linear)."""
    xs = [math.log(r) for r in r_values]
    ys = [math.log(v) for v in s2]
    m = len(xs)
    mx = sum(xs) / m
    my = sum(ys) / m
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    return num / den


def window_slopes(beta: float, ks: Sequence[float], r_values: Sequence[int]) -> dict:
    """Effective slope over the whole window and over the two half-windows."""
    s2 = s2_exact(r_values, beta, ks)
    lo, hi = GLM_WINDOW
    sel = [(r, v) for r, v in zip(r_values, s2) if lo <= r <= hi]
    r_sel = [r for r, _ in sel]
    v_sel = [v for _, v in sel]
    mid = len(sel) // 2
    return {
        "slope_window": loglog_slope(r_sel, v_sel),
        "slope_first_half": loglog_slope(r_sel[:mid], v_sel[:mid]),
        "slope_second_half": loglog_slope(r_sel[mid:], v_sel[mid:]),
        "r_window": [r_sel[0], r_sel[-1]],
        "s2_at_window_ends": [v_sel[0], v_sel[-1]],
    }


def beta3_log_constant() -> float:
    """C0 in F(x) = int_x^inf q^-3 (1-cos q) dq = 0.5*ln(1/x) + C0 + O(x^2).

    Computed by integrating the series for 1-cos q term by term on [x,1] and adding
    the numerically-integrated tail on [1, inf). Two independent routes (series and
    the x-independence of F(x) - 0.5 ln(1/x)) agree; see the section E printout.
    """
    terms = 40
    tail = 0.0
    lo, hi, steps = 1.0, 20000.0, 400000
    h = (hi - lo) / steps
    for i in range(steps):
        q = lo + (i + 0.5) * h
        tail += q ** -3 * (1.0 - math.cos(q)) * h

    def F(x: float) -> float:
        s = 0.0
        for nn in range(1, terms + 1):
            coef = (-1) ** (nn + 1) / math.factorial(2 * nn)
            p = 2 * nn - 3
            val = math.log(1.0) - math.log(x) if p == -1 else (1.0 - x ** (p + 1)) / (p + 1)
            s += coef * val
        return s + tail

    x = 1e-10
    return F(x) - 0.5 * math.log(1.0 / x)


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=None)
    args = ap.parse_args(list(argv) if argv is not None else None)

    r_values = list(range(1, 512))
    out: dict = {
        "definition": "S2(r) = 2 * sum_k k^-beta (1 - cos(k r)); direct summation, no FFT",
        "glm_window": list(GLM_WINDOW),
        "glm_r1_threshold_on_beta3_slope": GLM_R1_THRESHOLD,
        "glm_reported_beta3_slope": GLM_SEED_SLOPE_BETA3,
        "cases": [],
    }

    print("== A) GLM's band and N: [j=3..900] of N=4096, i.e. k in "
          f"[{2*math.pi*3/GLM_N:.5f}, {2*math.pi*900/GLM_N:.5f}] ==")
    ks_glm = modes(GLM_N, GLM_JMIN, GLM_JMAX)
    r_star = 1.0 / ks_glm[0]
    print(f"   smallest wavenumber k_lo={ks_glm[0]:.5f} -> saturation cross-over r_* ~ 1/k_lo = {r_star:.1f}")
    print(f"   window [{GLM_WINDOW[0]}, {GLM_WINDOW[1]}] vs r_*={r_star:.1f}: "
          f"k_lo*r spans {ks_glm[0]*GLM_WINDOW[0]:.3f} .. {ks_glm[0]*GLM_WINDOW[1]:.3f}")
    for beta in (5.0 / 3.0, 2.0, 3.0):
        w = window_slopes(beta, ks_glm, r_values)
        predicted = beta - 1.0
        rec = {"case": "glm_band", "N": GLM_N, "j": [GLM_JMIN, GLM_JMAX], "beta": beta,
               "predicted_asymptotic_slope": predicted, **w}
        out["cases"].append(rec)
        print(f"   beta={beta:.4f}  asymptotic beta-1={predicted:.4f}   "
              f"measured-over-window={w['slope_window']:.4f}   "
              f"[first half {w['slope_first_half']:.4f}, second half {w['slope_second_half']:.4f}]")

    print("== B) same beta=3, varying the low-k edge (the only change: k_lo) ==")
    for jmin in (3, 8, 20, 50, 100, 200):
        ks = modes(GLM_N, jmin, GLM_JMAX)
        w = window_slopes(3.0, ks, r_values)
        rec = {"case": "vary_jmin", "N": GLM_N, "j": [jmin, GLM_JMAX], "beta": 3.0,
               "k_lo": ks[0], "r_star": 1.0 / ks[0],
               "k_lo_times_window": [ks[0] * GLM_WINDOW[0], ks[0] * GLM_WINDOW[1]], **w}
        out["cases"].append(rec)
        print(f"   j_min={jmin:>3}  k_lo={ks[0]:.5f}  r_*={1/ks[0]:7.1f}  "
              f"k_lo*r_in_window=[{ks[0]*GLM_WINDOW[0]:.3f},{ks[0]*GLM_WINDOW[1]:.3f}]  "
              f"slope={w['slope_window']:.4f}")

    print("== C) beta=3 with a band wide enough that k_lo*r << 1 across the window ==")
    for n, jmin, jmax in ((4096, 1, 2048), (16384, 1, 4096), (65536, 1, 8192)):
        ks = modes(n, jmin, jmax)
        w = window_slopes(3.0, ks, r_values)
        rec = {"case": "wide_band", "N": n, "j": [jmin, jmax], "beta": 3.0,
               "k_lo": ks[0], "r_star": 1.0 / ks[0],
               "k_lo_times_window": [ks[0] * GLM_WINDOW[0], ks[0] * GLM_WINDOW[1]], **w}
        out["cases"].append(rec)
        print(f"   N={n:>6} k_lo={ks[0]:.6f} r_*={1/ks[0]:8.1f}  "
              f"k_lo*r_in_window=[{ks[0]*GLM_WINDOW[0]:.4f},{ks[0]*GLM_WINDOW[1]:.4f}]  "
              f"slope={w['slope_window']:.4f}")

    print("== D) beta=2 control: does the same machinery reproduce the r^1 law? ==")
    for jmin in (3, 20):
        ks = modes(GLM_N, jmin, GLM_JMAX)
        w = window_slopes(2.0, ks, r_values)
        print(f"   j_min={jmin}: slope={w['slope_window']:.4f} (asymptotic 1.0)")
        out["cases"].append({"case": "beta2_control", "N": GLM_N, "j": [jmin, GLM_JMAX],
                             "beta": 2.0, "predicted_asymptotic_slope": 1.0, **w})

    print("== E) verify the beta=3 log form F(x) = 0.5*ln(1/x) + C0 ==")
    c0 = beta3_log_constant()
    # RETRACTED (3): the original line read "(constant to 5e-7)" and evaluated a single x.
    # That overstates the accuracy by two orders of magnitude and does not test the
    # advertised range. The supported value is 0.4613932125491026 with a sweep spread of
    # 2.1e-4; see s2_local_exponent.py for the 20-point convergence table.
    print(f"   fitted C0 = {c0:.6f}  [RETRACTED precision claim; see s2_local_exponent.py]")
    out["beta3_log_constant_C0"] = c0
    ks_glm_l = modes(GLM_N, GLM_JMIN, GLM_JMAX)
    pred: list[dict] = []
    for r in (16, 32, 64, 128, 208):
        a = ks_glm_l[0] * r
        if a <= 0:
            continue
        L = math.log(1.0 / a) + 2 * c0
        pred.append({"r": r, "k_lo_r": a, "log_formula_local_exponent": 2.0 - 1.0 / L})
    out["log_formula_local_exponent"] = pred
    for p in pred:
        print(f"   r={p['r']:>3}  k_lo*r={p['k_lo_r']:.4f}  local exponent from log formula = "
              f"{p['log_formula_local_exponent']:.4f}")

    # verdict
    glm_case = next(c for c in out["cases"] if c["case"] == "glm_band" and c["beta"] == 3.0)
    out["conclusion"] = {
        "beta3_log_corrected_form": "S2(r) = 2C r^2 [ 0.5*ln(1/(k_lo r)) + 0.461390... ] for k_lo r << 1",
        "beta3_local_exponent": "2 - 1/(ln(1/(k_lo r)) + 0.922781)",
        "beta3_glm_band_measured_exponent": glm_case["slope_window"],
        # RETRACTED (1): this boolean used an unexplained 0.05 tolerance. It is kept only
        # to make the retraction visible; it is NOT the basis for any exclusion. The
        # exclusion rests on reproduce_and_compare.py reproducing GLM's own estimator to
        # ~1e-14 on GLM's own 13 lags.
        "glm_reported_1_54_matches_independent_computation_RETRACTED_TOLERANCE_ARGUMENT":
            abs(glm_case["slope_window"] - GLM_SEED_SLOPE_BETA3) < 0.05,
        "retraction_1": ("the 0.05-tolerance argument is withdrawn: the estimator has zero "
                         "realisation noise and the 0.0489 gap was a lag-set difference"),
        "saturation_cross_over_r_star": r_star,
        "verdict": (
            "GLM's implementation is CORRECT and GLM's 'r^2 log' diagnosis is CORRECT "
            "(the beta=3 integral diverges logarithmically at the infrared end; verified "
            "analytically above). The measurement is not a wrong exponent and not an "
            "implementation bug. What failed is the FROZEN THRESHOLD: 1.7 is close to the "
            "log-corrected LOCAL exponent at the smallest lag (2-1/(ln(1/0.0736)+0.9228) = "
            "1.717), but no WINDOW AVERAGE over [16,249] can reach it. RETRACTED (2): the original "
            "text continued 'because the upper lags are at or past the saturation cross-over "
            "r_* = 1/k_lo = 217, where the local exponent has already dropped to 0'. That is "
            "wrong -- k_lo*r only reaches 1.146 and the exact local exponent at r=217 is "
            "1.1630, not 0. The corrected statement is that the window crosses the ONSET of "
            "finite-band turnover and the local exponent falls through it; see "
            "s2_local_exponent.py for the exact table. Answer: D (finite-window effect), "
            "with A as the proximate cause of the recorded FAIL; NOT B; C incomplete rather "
            "than wrong. The deeper design flaw is "
            "that beta=3 has NO clean scaling exponent in a finite band at all: the local "
            "exponent is r-dependent, so the control was never a power law to begin with."
        ),
    }
    if args.json:
        os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(out, fh, indent=2, sort_keys=True)
        print(f"[written] {args.json}")
    print(json.dumps(out["conclusion"], indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

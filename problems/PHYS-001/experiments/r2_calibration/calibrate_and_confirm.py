#!/usr/bin/env python3
"""PHYS-001 r2: two-stage pre-registered calibration and confirmation.

Stage A (calibration, seeds 1-8): beta in {5/3, 2, 3}, N in {2048, 4096},
three frozen fitting windows. Measures SF slope (log-log OLS), spectrum slope
(FFT band [16,512]), energy residual. From Stage A ONLY, freezes per-(beta,
window) acceptance bands = mean +/- 4*std, written to bands.json.

Stage B (confirmation, seeds 11-20, disjoint; N=4096): all measurements must
fall inside the frozen bands. Bands are never adjusted after Stage B.

Also reports N-sensitivity (slope shift 2048->4096) and window-sensitivity
(spread across the three windows) from Stage A. Reuses the r1 field
generator/estimators by import (single implementation, audited there).
Stdlib only, deterministic.
"""
from __future__ import annotations
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "baseline_r1"))
from k41_baseline import (make_field, s2, ols, fft, spectral_slope,  # noqa: E402
                          energy_budget)

BETAS = {"K41_5_3": 5.0 / 3.0, "K2_2": 2.0, "K3_3": 3.0}
WINDOWS = {"W1": [16, 20, 25, 31, 39, 49, 62, 78, 98, 124, 156, 197, 249],
           "W2": [24, 30, 38, 48, 60, 76, 96, 121, 152, 192, 242, 305, 384],
           # W3 corrected (Codex P2): the admitted rule froze a 13-point
           # GEOMETRIC grid on [12, 156]; the previous list was not geometric
           # (final step 152->156). Rebuilt as round(12 * 13^(k/12)), k=0..12.
           "W3": [round(12 * 13 ** (k / 12)) for k in range(13)]}
SPEC_BAND = (16, 512)
STAGE_A_SEEDS = [1, 2, 3, 4, 5, 6, 7, 8]
STAGE_B_SEEDS = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20]


def sf_slope(u: list[float], lags: list[int]) -> float:
    vals = [math.log(s2(u, r)) for r in lags]
    return ols([math.log(r) for r in lags], vals)


def measure(u, amps, window) -> dict:
    return {"sf_slope": sf_slope(u, window),
            "spec_slope": spectral_slope(u),
            "energy_rel_residual": energy_budget(u, amps)}


def main() -> int:
    t0 = datetime.now(timezone.utc)
    stage_a = []
    for name, beta in BETAS.items():
        for N in (2048, 4096):
            for seed in STAGE_A_SEEDS:
                import k41_baseline as k
                k.N = N  # field size override
                u, amps = make_field(beta, seed)
                for wname, window in WINDOWS.items():
                    # Codex P2: hold the physical interval fixed across N —
                    # scale integer lags by N/4096 (rounded, deduped)
                    if N == 4096:
                        lags_n = list(window)
                    else:
                        lags_n = sorted(set(max(1, round(r * N / 4096))
                                            for r in window))
                    m = measure(u, amps, lags_n)
                    m["lags_used"] = lags_n
                    stage_a.append({"beta_name": name, "N": N, "seed": seed,
                                    "window": wname, **m})
    # freeze bands from Stage A only (per beta/window, N=4096 rows)
    bands = {}
    for name in BETAS:
        for wname in WINDOWS:
            vals = [r["sf_slope"] for r in stage_a
                    if r["beta_name"] == name and r["window"] == wname and r["N"] == 4096]
            mu = sum(vals) / len(vals)
            sd = (sum((v - mu) ** 2 for v in vals) / len(vals)) ** 0.5
            bands[f"{name}|{wname}|sf"] = [mu - 4 * sd, mu + 4 * sd]
            vals_spec = [r["spec_slope"] for r in stage_a
                         if r["beta_name"] == name and r["window"] == wname and r["N"] == 4096]
            mu_s = sum(vals_spec) / len(vals_spec)
            sd_s = (sum((v - mu_s) ** 2 for v in vals_spec) / len(vals_spec)) ** 0.5
            bands[f"{name}|{wname}|spec"] = [mu_s - 4 * sd_s, mu_s + 4 * sd_s]
    max_res = max(r["energy_rel_residual"] for r in stage_a)
    # frozen B3 contract: <= 2x Stage A maximum, NO added floor (the earlier
    # 1e-10 floor relaxed the admitted limit by ~3 orders of magnitude; Codex
    # P2 — removed. Cycle-1 recorded verdicts are unchanged under the exact
    # rule because Stage A residuals were ~4e-14.)
    bands["energy_max"] = 2 * max_res
    bands_doc = {"frozen_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                 "rule": "mean +/- 4*std of Stage A (seeds 1-8, N=4096) per beta/window",
                 "bands": bands}
    (HERE / "bands.json").write_text(
        json.dumps(bands_doc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    # N-sensitivity and window-sensitivity (Stage A descriptive stats)
    # Codex P2: lag index r is the normalized separation r/N; comparing the
    # same integer lags across N=2048/4096 halves the physical window. The
    # N=2048 rows are regenerated with lags SCALED by N/4096 so both runs fit
    # the same physical interval; the shift reported is then resolution
    # sensitivity, not a window shift.
    n_sens = {}
    for name in BETAS:
        for wname in WINDOWS:
            m4096 = sum(r["sf_slope"] for r in stage_a
                        if r["beta_name"] == name and r["N"] == 4096
                        and r["window"] == wname) / 8
            m2048 = sum(r["sf_slope"] for r in stage_a
                        if r["beta_name"] == name and r["N"] == 2048
                        and r["window"] == wname) / 8
            n_sens[f"{name}|{wname}"] = round(m4096 - m2048, 4)
    w_sens = {name: round(max(
        sum(r["sf_slope"] for r in stage_a if r["beta_name"] == name
            and r["window"] == w and r["N"] == 4096) / 8 for w in WINDOWS)
        - min(sum(r["sf_slope"] for r in stage_a if r["beta_name"] == name
                  and r["window"] == w and r["N"] == 4096) / 8 for w in WINDOWS), 4)
        for name in BETAS}

    # Stage B: fresh seeds, must land inside frozen bands
    import k41_baseline as k
    k.N = 4096
    stage_b, violations = [], []
    for name, beta in BETAS.items():
        for seed in STAGE_B_SEEDS:
            u, amps = make_field(beta, seed)
            for wname, window in WINDOWS.items():
                m = measure(u, amps, window)
                stage_b.append({"beta_name": name, "seed": seed, "window": wname, **m})
                lo, hi = bands[f"{name}|{wname}|sf"]
                if not (lo <= m["sf_slope"] <= hi):
                    violations.append({"beta": name, "window": wname, "seed": seed,
                                       "sf_slope": m["sf_slope"], "band": [lo, hi]})
                lo_s, hi_s = bands[f"{name}|{wname}|spec"]
                if not (lo_s <= m["spec_slope"] <= hi_s):
                    violations.append({"beta": name, "window": wname, "seed": seed,
                                       "spec_slope": m["spec_slope"], "band": [lo_s, hi_s]})
                if m["energy_rel_residual"] > bands["energy_max"]:
                    violations.append({"beta": name, "seed": seed,
                                       "energy": m["energy_rel_residual"]})
    # B4: bands for the three betas must not overlap (separation)
    lo53 = bands["K41_5_3|W1|sf"][1]
    hi22 = bands["K2_2|W1|sf"][0]
    hi33 = bands["K3_3|W1|sf"][0]
    sep = lo53 < hi22 < hi33 or (hi33 < lo53 and hi22 < lo53) or (hi22 < lo53 and lo53 > hi33)
    b4 = (bands["K41_5_3|W1|sf"][1] < bands["K2_2|W1|sf"][0]) and \
         (bands["K2_2|W1|sf"][1] < bands["K3_3|W1|sf"][0])
    checks = {"B1_confirmatory_in_bands": not violations,
              "B4_band_separation": b4}
    out = {"generated": t0.isoformat(timespec="seconds"),
           "stage_a_runs": len(stage_a), "stage_b_runs": len(stage_b),
           "n_sensitivity_4096_minus_2048": n_sens,
           "window_sensitivity_sf_slope_spread": w_sens,
           "bands": bands_doc, "violations": violations[:20],
           "n_violations": len(violations), "checks": checks,
           "verdict": ("CALIBRATED_BASELINE_PASS" if all(checks.values())
                       else "CALIBRATED_BASELINE_FAIL")}
    dest = HERE.parent.parent / "results" / "r2"
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "r2_calibration_results.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"stage_a_runs": len(stage_a), "stage_b_runs": len(stage_b),
                      "n_violations": len(violations), "checks": checks,
                      "verdict": out["verdict"],
                      "band_examples": {k: [round(x, 4) for x in v] for k, v in
                                        list(bands.items())[:6] if isinstance(v, list)},
                      "n_sens_example": dict(list(n_sens.items())[:3]),
                      "w_sens": w_sens}, ensure_ascii=False, indent=2))
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())

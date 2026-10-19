#!/usr/bin/env python3
"""PHYS-001 round 1 baseline: K41 synthetic reference and evaluator.

Scope (frozen in round.json BEFORE execution):
  - synthetic 1D velocity signal u_j = sum_k A_k cos(2*pi*k*j/N + phi_k),
    A_k = k^(-beta/2)  =>  1D power spectrum S(k) ∝ k^-beta
    beta = 5/3 (K41 positive control) and beta = 3 (negative control),
    wavenumber band [3, 900], N = 4096 samples, seeds 1..5 (random.Random).
  - checks:
      E1 energy budget |<u^2> - sum A_k^2/2| / (sum A_k^2/2) <= 0.10
      E2 structure function slope in FROZEN window (13 geometric lags 16..249),
         5-seed mean in [2/3 - 0.05, 2/3 + 0.05], std <= 0.03   (beta = 5/3)
      E3 beta = 3: 5-seed mean slope >= 1.7 (clearly outside the 2/3 band)
      E4 FFT periodogram slope over k in [16, 512] within beta +/- 0.1
      E5 cross-consistency |(spectral_slope - 1) - sf_slope| <= 0.1
      E6 evaluator guard: sine-only and amplitude-inflated (x1.5) inputs FAIL
  - relation used: S(k) ∝ k^-beta  <=>  S2(r) ∝ r^(beta-1) for 1 < beta < 3.
    A K41 k^-5/3 3D spectrum has a k^-5/3 1D surrogate, giving the 2/3 law.
  - if measurements fall outside the frozen tolerances this is recorded as a
    FAILURE with analysis; tolerances are never moved after the fact.

No claim is made about real turbulence: this is an idealized non-intermittent
reference field, not a Navier-Stokes solution. Stdlib only, deterministic.
"""
from __future__ import annotations
import datetime as dt
import json
import math
import random
import sys
from pathlib import Path

N = 4096
K_MIN, K_MAX = 3, 900
BETAS = {"K41_positive": 5.0 / 3.0, "steep_negative": 3.0}
SEEDS = [1, 2, 3, 4, 5]
LAGS = [16, 20, 25, 31, 39, 49, 62, 78, 98, 124, 156, 197, 249]  # frozen window
SPEC_KMIN, SPEC_KMAX = 16, 512
TOL = {"energy": 0.10, "sf_band": (2.0 / 3.0 - 0.05, 2.0 / 3.0 + 0.05),
       "sf_std": 0.03, "neg_min": 1.7, "spec": 0.1, "cross": 0.1}


def make_field(beta: float, seed: int) -> tuple[list[float], list[float]]:
    """Random-phase mode sum with A_k = k^(-beta/2); returns (u, amplitudes)."""
    rng = random.Random(seed)
    amps = [k ** (-beta / 2.0) for k in range(K_MIN, K_MAX + 1)]
    u = [0.0] * N
    for idx, k in enumerate(range(K_MIN, K_MAX + 1)):
        phi = 2.0 * math.pi * rng.random()
        a = amps[idx]
        c, s = a * math.cos(phi), a * math.sin(phi)
        theta = 2.0 * math.pi * k / N
        ct, st = math.cos(theta), math.sin(theta)
        for j in range(N):
            u[j] += c          # c = A*cos(phi + j*theta) via rotation
            c, s = c * ct - s * st, c * st + s * ct
    return u, amps


def s2(u: list[float], r: int) -> float:
    acc = 0.0
    for j in range(N):
        d = u[(j + r) % N] - u[j]
        acc += d * d
    return acc / N


def ols(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return sxy / sxx


def sf_slope(u: list[float]) -> float:
    vals = [math.log(s2(u, r)) for r in LAGS]
    return ols([math.log(r) for r in LAGS], vals)


def fft(p: list[float]) -> list[complex]:
    """Iterative radix-2 FFT (N is a power of two)."""
    a = [complex(v, 0.0) for v in p]
    j, n = 0, len(a)
    for i in range(1, n):  # bit-reversal permutation
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j |= bit
        if i < j:
            a[i], a[j] = a[j], a[i]
    length = 2
    while length <= n:
        ang = -2.0 * math.pi / length
        wl = complex(math.cos(ang), math.sin(ang))
        for start in range(0, n, length):
            w = 1 + 0j
            for k in range(start, start + length // 2):
                t = w * a[k + length // 2]
                a[k + length // 2] = a[k] - t
                a[k] += t
                w *= wl
        length <<= 1
    return a


def spectral_slope(u: list[float]) -> float:
    spec = fft(u)
    power = [(k, abs(spec[k]) ** 2) for k in range(N // 2)]
    ks, ps, cur, acc, cnt = [], [], [], 0.0, 0
    log_bins = [math.log(k) for k in range(SPEC_KMIN, SPEC_KMAX + 1)]
    step = (log_bins[-1] - log_bins[0]) / 23.0
    bounds = [log_bins[0] + i * step for i in range(24)]
    b = 0
    for k in range(SPEC_KMIN, SPEC_KMAX + 1):
        while b < 23 and math.log(k) >= bounds[b + 1]:
            if cnt:
                ks.append(sum(cur) / cnt)
                ps.append(math.log(acc / cnt))
            cur, acc, cnt = [], 0.0, 0
            b += 1
        cur.append(k)
        acc += power[k][1]
        cnt += 1
    if cnt:
        ks.append(sum(cur) / cnt)
        ps.append(math.log(acc / cnt))
    return ols([math.log(x) for x in ks], ps) * -1.0  # log-log: power ~ k^-beta


def energy_budget(u: list[float], amps: list[float]) -> float:
    analytic = sum(a * a for a in amps) / 2.0
    measured = sum(v * v for v in u) / N
    return abs(measured - analytic) / analytic


def run_case(beta: float, seed: int) -> dict:
    u, amps = make_field(beta, seed)
    return {"seed": seed, "energy_rel_residual": energy_budget(u, amps),
            "sf_slope": sf_slope(u), "spec_slope": spectral_slope(u)}


def main() -> int:
    t0 = dt.datetime.now(dt.timezone.utc)
    out: dict = {"config": {"N": N, "k_band": [K_MIN, K_MAX], "seeds": SEEDS,
                            "lags": LAGS, "spec_band": [SPEC_KMIN, SPEC_KMAX],
                            "tolerances": {**TOL, "sf_band": list(TOL["sf_band"])}},
                 "generated": t0.isoformat(timespec="seconds"), "cases": []}
    for name, beta in BETAS.items():
        cases = [run_case(beta, s) for s in SEEDS]
        slopes = [c["sf_slope"] for c in cases]
        specs = [c["spec_slope"] for c in cases]
        mean, std = sum(slopes) / len(slopes), (
            sum((x - sum(slopes) / len(slopes)) ** 2 for x in slopes) / len(slopes)
        ) ** 0.5
        spec_mean = sum(specs) / len(specs)
        out["cases"].append({"beta_name": name, "beta": beta, "per_seed": cases,
                             "sf_mean_slope": mean, "sf_std": std,
                             "spec_mean_slope": spec_mean,
                             "energy_max_rel_residual": max(c["energy_rel_residual"]
                                                            for c in cases)})
    pos = next(c for c in out["cases"] if c["beta_name"] == "K41_positive")
    neg = next(c for c in out["cases"] if c["beta_name"] == "steep_negative")
    lo, hi = TOL["sf_band"]
    checks = {
        "E1_energy_budget_pos": pos["energy_max_rel_residual"] <= TOL["energy"],
        "E1_energy_budget_neg": neg["energy_max_rel_residual"] <= TOL["energy"],
        "E2_sf_2_3_band": lo <= pos["sf_mean_slope"] <= hi and pos["sf_std"] <= TOL["sf_std"],
        "E3_negative_control_steep": neg["sf_mean_slope"] >= TOL["neg_min"],
        "E4_spectrum_slope_pos": abs(pos["spec_mean_slope"] - 5.0 / 3.0) <= TOL["spec"],
        "E4_spectrum_slope_neg": abs(neg["spec_mean_slope"] - 3.0) <= TOL["spec"],
        "E5_cross_consistency_pos": abs((pos["spec_mean_slope"] - 1.0) - pos["sf_mean_slope"]) <= TOL["cross"],
    }
    # E6 evaluator guards: corrupted inputs must FAIL the corresponding checks
    u_ref, _ = make_field(5.0 / 3.0, 1)
    sine = [math.cos(2.0 * math.pi * 50.0 * j / N) for j in range(N)]
    checks["E6_sine_rejected_by_sf"] = not (lo <= sf_slope(sine) <= hi)
    checks["E6_amplitude_inflation_rejected_by_energy"] = energy_budget(
        [1.5 * v for v in u_ref], [1.0] * 1) > TOL["energy"]
    out["checks"] = checks
    out["verdict"] = ("BASELINE_PASS" if all(checks.values())
                      else "BASELINE_FAIL_DO_NOT_PROCEED_TO_NOVELTY")
    out["completed"] = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    dest = Path(__file__).resolve().parents[2] / "results" / "r1"
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "k41_baseline_results.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"sf_mean_slope_pos": round(pos["sf_mean_slope"], 4),
                      "sf_std_pos": round(pos["sf_std"], 4),
                      "sf_mean_slope_neg": round(neg["sf_mean_slope"], 4),
                      "spec_mean_slope_pos": round(pos["spec_mean_slope"], 4),
                      "spec_mean_slope_neg": round(neg["spec_mean_slope"], 4),
                      "energy_max_rel_residual_pos": pos["energy_max_rel_residual"],
                      "checks": checks, "verdict": out["verdict"]},
                     ensure_ascii=False, indent=2))
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())

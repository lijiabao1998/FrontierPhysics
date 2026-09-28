#!/usr/bin/env python3
"""PHYS-001 r2 cycle 2 — the ACTUAL executable that produced bands2.json and
r2_cycle2_results.json (Codex P1: previously run as an inline heredoc that was
never checked in, making the PASS unreproducible).

Semantics (Codex P1, downgraded honestly): because the synthetic field is a
full-period integer-Fourier mode sum, varying the seed changes ONLY phases;
energy, power spectrum and full-domain S2 are phase-independent up to
floating-point roundoff. Stage A/B therefore constitute a DETERMINISTIC
NUMERICAL REPRODUCIBILITY CHECK, not a statistical confirmation. A true
statistical generalization test would need to vary amplitudes/spectra or use
independent finite subwindows (deferred; see ROUND_REPORT).

Design: Stage A (seeds 21-28) -> freeze bands with rule
mean +/- max(4*std, FLOOR=0.01)  -> Stage B (seeds 31-40) must land inside.
B4: beta-band separation is checked in ALL THREE frozen windows.

Usage: python calibrate_and_confirm_cycle2.py   (writes bands2.json and
results/r2/r2_cycle2_results.json; deterministic; stdlib only)
"""
from __future__ import annotations
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "baseline_r1"))
import calibrate_and_confirm as c  # noqa: E402
import k41_baseline as k  # noqa: E402

A_SEEDS = [21, 22, 23, 24, 25, 26, 27, 28]
B_SEEDS = [31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
FLOOR = 0.01


def main() -> int:
    stage_a = []
    for name, beta in c.BETAS.items():
        k.N = 4096
        for seed in A_SEEDS:
            u, amps = c.make_field(beta, seed)
            for wname, window in c.WINDOWS.items():
                stage_a.append({"beta": name, "window": wname,
                                **c.measure(u, amps, window)})
    energy_max = 2 * max(r["energy_rel_residual"] for r in stage_a)
    bands = {"energy_max": energy_max}
    for name in c.BETAS:
        for wname in c.WINDOWS:
            for key, field in ((f"{name}|{wname}|sf", "sf_slope"),
                               (f"{name}|{wname}|spec", "spec_slope")):
                vals = [r[field] for r in stage_a
                        if r["beta"] == name and r["window"] == wname]
                mu = sum(vals) / len(vals)
                sd = (sum((v - mu) ** 2 for v in vals) / len(vals)) ** 0.5
                half = max(4 * sd, FLOOR)
                bands[key] = [mu - half, mu + half]
    bdoc = {"frozen_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "rule": "mean +/- max(4*std, 0.01); Stage A seeds 21-28, N=4096",
            "floor": FLOOR,
            "bands": {k2: [round(x, 4) for x in v] for k2, v in bands.items()}}
    (HERE / "bands2.json").write_text(
        json.dumps(bdoc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    stage_b, viol = [], []
    for name, beta in c.BETAS.items():
        k.N = 4096
        for seed in B_SEEDS:
            u, amps = c.make_field(beta, seed)
            for wname, window in c.WINDOWS.items():
                m = c.measure(u, amps, window)
                stage_b.append(m)
                lo, hi = bands[f"{name}|{wname}|sf"]
                if not (lo <= m["sf_slope"] <= hi):
                    viol.append({"beta": name, "window": wname, "seed": seed,
                                 "sf_slope": round(m["sf_slope"], 4)})
                lo, hi = bands[f"{name}|{wname}|spec"]
                if not (lo <= m["spec_slope"] <= hi):
                    viol.append({"beta": name, "window": wname, "seed": seed,
                                 "spec_slope": round(m["spec_slope"], 4)})
                # B3 (frozen contract): energy residual <= 2x Stage A maximum
                if m["energy_rel_residual"] > bands["energy_max"]:
                    viol.append({"beta": name, "window": wname, "seed": seed,
                                 "energy": m["energy_rel_residual"]})
    # B4: separation in ALL THREE windows (Codex P2)
    sep = {}
    for wname in c.WINDOWS:
        sep[wname] = (bands[f"K41_5_3|{wname}|sf"][1] < bands[f"K2_2|{wname}|sf"][0]
                      and bands[f"K2_2|{wname}|sf"][1] < bands[f"K3_3|{wname}|sf"][0])
    checks = {"B1_in_bands": not viol, "B4_separation_all_windows": all(sep.values())}
    out = {"cycle": 2, "executable": "calibrate_and_confirm_cycle2.py",
           "semantics": "deterministic numerical reproducibility check "
                        "(seeds vary phases only; not a statistical confirmation)",
           "stage_a_seeds": A_SEEDS, "stage_b_seeds": B_SEEDS,
           "stage_a_n": len(stage_a), "stage_b_n": len(stage_b),
           "n_violations": len(viol), "violations": viol[:10],
           "B4_separation_by_window": sep,
           "checks": checks,
           "verdict": ("CALIBRATED_BASELINE_PASS" if all(out_reg["checks"].values())
                       else "CALIBRATED_BASELINE_FAIL"),
           "registration_status": out_reg["registration_status"],
           "bands_doc": bdoc,
           "generated": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    dest = HERE.parent.parent / "results" / "r2"
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "r2_cycle2_results.json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("cycle", "stage_a_n", "stage_b_n",
                                          "n_violations", "registration_status", "checks", "verdict",
                                          "B4_separation_by_window")},
                     ensure_ascii=False, indent=2))
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())

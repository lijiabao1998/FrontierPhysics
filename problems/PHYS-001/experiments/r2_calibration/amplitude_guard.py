#!/usr/bin/env python3
"""Amplitude-inflation guard with teeth (Codex P2): uses the REAL make_field
amplitudes as the analytic reference. Unmodified field must PASS the frozen
energy tolerance; a 1.5x inflated field must FAIL. Writes
results/r2/amplitude_guard_check.json deterministically (no clock field)."""
from __future__ import annotations
import json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "baseline_r1"))
import k41_baseline as k

def main() -> int:
    u, amps = k.make_field(5.0 / 3.0, 1)
    res_ref = k.energy_budget(u, amps)
    res_15 = k.energy_budget([1.5 * v for v in u], amps)
    tol = k.TOL["energy"]
    guard = {"reference_residual": res_ref, "reference_pass": res_ref <= tol,
             "inflated_residual": res_15, "inflated_fail": res_15 > tol,
             "tolerance": tol, "inflation_factor": 1.5,
             "note": "guard uses the actual make_field amplitudes as the "
                     "analytic reference (Codex P2)"}
    dest = HERE.parent.parent / "results" / "r2"
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "amplitude_guard_check.json").write_text(
        json.dumps(guard, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(guard, indent=2))
    ok = guard["reference_pass"] and guard["inflated_fail"]
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())

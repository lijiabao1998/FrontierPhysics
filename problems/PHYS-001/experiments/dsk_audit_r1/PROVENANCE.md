# Provenance / environment — PHYS-001 independent audit (DeepSeek r1)

Independent implementation: `s2_scaling_audit.py`. Evaluates the **definition**

    S2(r) = 2 * sum_k k^-beta (1 - cos(k r))

by direct summation over modes. No FFT, no structure-function estimator, and no
shared code with GLM's `k41_baseline.py` or `calibrate_and_confirm.py` (neither was
read before writing this). Python 3.13.x, standard library only.

Primary source consulted: arXiv:2607.26896 (abstract level only, as in GLM r1).
GLM's r2 branch (`glm/PHYS-001-calib-r2`) was read as a **report**, not executed;
that limitation is stated in PHYS001_INDEPENDENT_AUDIT.md section 6.

Reproduce:

    python s2_scaling_audit.py --json results/s2_scaling_audit.json

Outputs `results/s2_scaling_audit.json` with every case recorded, including the
zero-perturbation controls (beta=2, beta=5/3) and the four band/window probes.

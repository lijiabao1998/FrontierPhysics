# Independent protocol review before computation

Reviewer: separate read-only `independent_review` session, 2026-10-11 UTC. No research computation or reference main was run by the reviewer.

The reviewer checked pinned governance, both protocols and the immutable source manifests; all 10 snapshots matched their SHA256 and byte counts. Historical Physics r2 residuals remain classified INHERITED_NOT_RERUN.

Before freezing, exact u[j] formulas and fixture coefficients were added. Sine is generated directly and has positive-frequency coefficient -i*N*A/2. Sign mutation uses phased mixed data; simultaneous 1.5x input/reference scaling is detected by failure of the inflated negative control. Wrong [1.0] reference is rejected because the clean control fails.

The reviewer requested independent fixed names for the top-level required checks. This was incorporated; the per-case expected name list is already generated independently from the frozen protocol. Only source-reviewed `fft` and `energy_budget` AST definitions execute, with fixture-local N; no source module top-level, original make_field or main runs.

The reviewer found the DFT sign/normalization, Parseval, mixed/sine and amplitude-poison fixtures within scope. This is independent static protocol/code review, not an independent experimental rerun. Post-run review is recorded separately.

# FROZEN — superseded as the canonical merge unit

This directory tree (on `dsk/PHYS-001-independent-audit-r1`, PR #3) is **frozen**. It is kept as
history; it is **not** the merge unit for the PHYS-001 audit evidence any more.

| unit | branch | PR |
|---|---|---|
| **canonical evidence** (the claims, with the merge gate) | `dsk/PHYS-001-canonical-evidence` | #4 |
| **audit tooling** (hash automation, rerun generator, guards, convergence tooling) | `dsk/PHYS-001-audit-tooling` | #5 |
| this tree (frozen history) | `dsk/PHYS-001-independent-audit-r1` | #3 |

## Claim hygiene applied before freezing

- `dsk_audit_r1/s2_scaling_audit.py` carries a **SUPERSEDED** banner and its three wrong claims
  (the 0.05-tolerance exclusion, the saturated-plateau mechanism, and the `C₀` precision claim)
  are retracted **in place**, so a reader of that file cannot inherit them.
- `dsk_audit_r1/PHYS001_INDEPENDENT_AUDIT.md` (revision 2) carries the retraction table and is
  the authoritative document of this tree.
- `DEEPSEEK_REVIEW_FrontierPhysics_PRs.md` carries a SUPERSEDED banner naming its two inaccurate
  statements (the saturation mechanism and the quoted constant) while keeping its verdicts.
- `runs/.../round.json` records that every artifact was regenerated with the remediated scripts,
  with hashes over git blob bytes.

The four retractions the owner required — the arbitrary 0.05 implementation-verification
threshold, the saturated-plateau claim, the invalid `[4,20]` gate (and its replacement), and the
unsupported `C₀` precision claim — are all recorded in
`problems/PHYS-001/experiments/canonical_evidence/CANONICAL_FACTS.md` on the canonical branch.

## No new functionality

This freeze adds no tooling. Further work belongs on `dsk/PHYS-001-audit-tooling`, which may
continue to receive review passes without gating the canonical evidence merge.

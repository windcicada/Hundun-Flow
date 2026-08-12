# HUNDUN-FLOW Task 11 no-96 pressure-extrema disambiguation

Status: user-approved on 2026-08-08. This amendment permanently removes
96-cubed flow simulations from the Task 11 acceptance and diagnostic schedule.
It supersedes the earlier 24/48/96 execution wording while preserving the
scientific thresholds and the product numerical path.

## Goal

Distinguish a moving global pressure-`Linf` maximum caused by geometry/stencil
phase from a genuine local immersed-pressure-row error using only 12, 24 and
48 cell simulations.

## Non-negotiable constraints

- No Task 11 command, CTest entry, systemd unit, or diagnostic runner may
  launch a 96-cubed flow solve.
- Existing 96-cubed logs remain historical evidence only; they are not rerun.
- Formal order thresholds, pressure gauge handling, PISO corrector count,
  conservation checks, and product solver rows are unchanged.
- The diagnostic uses the same production solver path and only changes
  test-side observation and surface-fixture policy.
- The fixed-surface comparison reuses the deterministic 48-level sphere STL
  for 12/24/48 grids; the control comparison regenerates the STL at each
  grid's `0.45h` target edge.

## Diagnostic design

Each diagnostic resolution records the top eight pressure-error extrema after
the existing global gauge subtraction. Each record contains the global cell,
logical cell, signed and absolute error, wall distance and wall-distance over
`h`, incident wall-link count, nearest authority-link donor count, donor
condition estimate, and deterministic donor fingerprint.

Two policies are compared:

1. `per_level_surface`: the existing production-test fixture policy;
2. `fixed_48_surface`: the same 48-level sphere triangulation for all three
   grids.

The diagnostic prints `L2`, near-wall `L2`, global `Linf`, top-eight extrema,
and the observed orders for 12→24 and 24→48. It does not replace the formal
order gate; it explains whether the gate is selecting different local rows.

The formal near-wall volume norm follows the controlling Stage 3 design: its
support has one fixed physical thickness across all refinement levels. For
the approved 12/24/48 matrix, that thickness is frozen before execution as
the coarse two-cell support,
`2*(L_ref/12) = L_ref/6`. A moving `2*h_max` support may be printed only as a
diagnostic partition and is not an acceptance row.

## Decision rule

- If fixed-surface top-eight-envelope errors become monotone and
  second-order while per-level STL remains nonmonotone, classify the failure
  as geometry/grid-phase plus max-selector sensitivity.
- If the same local row/donor family remains sub-second-order under both
  policies, continue with an A22 local residual derivation.
- If both policies are second-order in top-eight-envelope measures
  while the single global maximum still oscillates, record the global-`Linf`
  gate as a pre-asymptotic acceptance-design blocker; do not relax it silently.

## Acceptance scope

The retained Task 11 selector matrix is now 12/24/48 only for:

```text
sphere_uniform_acceptance
sphere_warped_acceptance
prism_warped_acceptance
```

The six deferred selectors remain deferred. No 96-cubed acceptance evidence
is required or produced under this amendment.

## Rollback boundary

The amendment is test-only. Reverting the diagnostic changes restores the
previous test-support and scheduling behavior without changing product code.

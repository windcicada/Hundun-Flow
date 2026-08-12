# HUNDUN-FLOW Task 11 Force-Authority Semantic Amendment (RED-S3 M2)

Status: user-approved option A on 2026-08-07 (force-path separation). This
amendment is additive to the RED-S1 signed-force authority and to the M2
force-consistency survey plan. It does not relax a threshold, add a corrector,
filter, damping or per-case tuning, and it does not change the solve rows.

## 1. Established root cause (executable evidence, 2026-08-06/07)

The RED-S3 decision-tree verdict is `C_exact` first order (1.23), and the M2
executable experiments established:

- The operator's wall-face pressure quadrature is anchored to the
  fluid-solid background grid face (the staircase surface). Even with exact
  face pressures the staircase quadrature converges only at ~1.35 order for
  the curved sphere (0.02402 -> 0.02751 N vs analytic 0.02975 N), because the
  staircase is O(h) from the true body.
- The surface quadrature (true body surface, triangle measures, wall-anchored
  authority values) converges at second order (wall-plan geometry 2.00,
  reconstruction 2.11).
- The shared row reconstruction's wall VALUE is already second order and
  already agrees with the per-link authority (the wall-value swap changes
  `C_exact` by less than 1%); the first-order component is the quadrature
  ANCHOR (face area and face value), not the wall value.
- Re-anchoring the SOLVE rows to the body surface breaks the exact Schur
  pressure operator's CG consistency: the first corrector's independent
  residual jumps to 2.2e-8 against a 1e-16 contract while the recursive CG
  residual is 2.1e-18 (a ~1e10 drift). The solve rows must keep the
  background face-flux structure.

## 2. Selected resolution: force-path separation

The operator force is the wall-row residual sum, but the wall-face pressure
value used by that sum is the shared per-link wall-anchored authority
reconstruction's extrapolated value at the body-surface patch centroid of the
link, weighted by the link's true body-surface measure vector:

```text
F_operator = sum over wall links of p_authority(centroid_link) * A_surface(link)
             + (all non-wall-face row terms, unchanged)
B_budget    = -F_operator
F_surface   = integral(sigma n_s dA), fluid-on-solid
C           = F_operator - F_surface
```

The pressure Poisson / momentum SOLVE keeps the background face-flux rows
exactly as frozen. The final force collection (`collect_final_force`) evaluates
the boundary rows in the force-authority mode, where each link's wall-face
pressure value is the per-link authority reconstruction's extrapolated value at
the body-surface patch centroid (the canonical sharp-interface boundary
pressure: ghost-cell-IBM `extrapolate_scalar`, Basilisk `embed_force`) and its
measure is the true surface measure (the same triangle-point-to-link
partition the surface quadrature uses). The operator and the surface therefore
share the same wall-value authority family while remaining independent
accumulation paths; the solve is bitwise unchanged.

Executable evidence (12/24 sphere-uniform, force rows): the unconstrained
extrapolated value yields `C_exact` order 3.30 and solved-state consistency
order 2.6; the wall-gradient-constrained value yields `C_exact` order 1.42
because its gradient-datum coupling `A_surface * boundary_coefficient * g`
carries a coherent ~0.7-order quadrature component. The open-source sharp
interface force references extrapolate the fluid pressure to the boundary and
do not pin the wall value by the gradient datum, so the unconstrained
extrapolation is the selected scientific estimator.

This replaces the M2 survey plan's row-level re-anchoring (which the RED-S3
experiments proved incompatible with the solve) with a force-only evaluation
of the same boundary row. It is the minimal amendment that keeps
`C_exact = O_h(U*) - S_h(U*)` second order without touching the solver.

## 3. Invariants preserved

- RED-S1 four-field report and sign semantics: `operator_force`,
  `budget_reaction`, `surface_traction`, `consistency`; the assembler remains
  the only conversion point and the budget is negated exactly once.
- The solve rows, the pressure operator, the corrector count and all solver
  controls are unchanged.
- No threshold, corrector, filter, damping or tuning change.
- The constant mode: the sum of the true surface-measure vectors over the
  closed body is zero, so the corrected operator force annihilates constant
  pressure exactly in the aggregate.

## 4. Verification gates (unchanged from M2, now on the force rows)

1. Exact-state wall-value oracle: the force row's wall-face value equals the
   authority wall value and stays second order.
2. `C_exact` (decision tree, on the force rows) >= 1.8.
3. `functional_selection_fast` pressure consistency >= 1.8.
4. Screen (24/48) then the formal 24/48/96 matrix.
5. Full focused regression: RED-S1 59 observations, Release 26-test set,
   Debug decomposition 5/5.

## 5. File whitelist

```text
immersed/src/ghost_stencil_plan.cpp            (per-link true surface measure)
immersed/include/hundun/immersed/ghost_stencil_plan.hpp  (private accessor)
finite_volume/src/immersed_operator.cpp        (force-authority row variant)
finite_volume/src/immersed_operator_test_access.hpp      (force rows snapshot)
finite_volume/src/immersed_boundary_authority_detail.hpp (scoped force mode)
flow/src/stage3_flow.cpp                       (force collection mode)
tests/support/stage3_mms.{hpp,cpp}             (exact-state oracle on force rows)
tests/numerical/test_laminar_ibm_order.cpp     (decision-tree selector)
```

## 6. Rollback boundary

Disposable candidate trees; the solve rows and the frozen RED-S1 evidence are
untouched; the force mode is scoped to `collect_final_force`. Any mutation
that reverts the force row to the background face-flux value must return
`C_exact` to first order (M2-M1, observed 1.231). Mutations:

| ID | Mutation | Designated failure |
|---|---|---|
| M2-M1 | force row falls back to the solve row (background face flux) | `C_exact` order ~1.23 |
| M2-M2 | wall value evaluated at the solid-side neighbour centre instead of the patch centroid | `C_exact` drops to ~1.78 (<1.8) |
| M2-M3 | wall value uses the wall-gradient-constrained authority value | `C_exact` drops to ~1.4 (constrained g-coupling) |
| M2-M4 | wall flux weighted by the background face area/normal instead of the true surface measure | `C_exact` drops below 1.8 |

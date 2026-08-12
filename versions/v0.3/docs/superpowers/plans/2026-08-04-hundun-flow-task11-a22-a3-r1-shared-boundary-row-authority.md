# Task 11 A22-A3/R1 Shared Boundary Row Authority Plan

The term "same authority" is governed by the basis-equivalence and consumer
responsibility contract in
`../specs/2026-08-04-hundun-flow-task11-a22-a3-r1-authority-equivalence-addendum.md`.
It does not require consumers to share a QR object address.

## Frozen sequence

1. Preserve the unchanged `sphere_fast` failure and exact log identities.
2. Add direct tests-only row-authority REDs before changing product behavior.
3. Factor the current A22 donor-union/constraint algebra into one private
   immutable row authority; do not duplicate its QR formulas.
4. Make the A22 mechanical row consume the factored authority and prove its
   existing polynomial/affine snapshots unchanged.
5. Route predictor current/history values and pressure gradients through the
   same row authority; A3 automatically consumes the same route.
6. Route wall-force pressure and velocity functionals through that authority
   without copying operator reaction into the surface integral.
7. Run direct Debug RED/GREEN and mutation oracles.
8. Run focused one/two/four-rank operator, PISO, transaction and force tests.
9. Run Release `sphere_probe_12`, then `sphere_fast` (`12^3/24^3`). Stop if
   either L2 or Linf order remains below 1.8.
10. Only after the fast gate is green, run the approved 24/48 development
    screen. Do not run 96 until all Task 11 focused gates are stable.

## Allowed implementation clusters

- private row authority and its tests;
- A22 adapter integration;
- predictor/A3 integration;
- force adapter integration;
- tests-only diagnostics needed to prove conditioning, exactness, rollback and
  decomposition invariance.

The clusters are one scientific repair and receive one final Task 11 verdict.

## Required direct evidence

- quadratic basis and cell-average reproduction;
- bounded close-parallel-link functional with a mutation using one independent
  link fit;
- exact `G(p+dp,g+dg)-G(p,g)` and `A du + delta G`;
- exact A3 operator/routed-face identity and wall positive zero;
- one/two/four-rank stable IDs and values;
- force uses the row authority yet remains an independent surface integral;
- failed inner/outer solve restores state and authority bitwise;
- tests-off and standalone public-header builds show no test seam or new API.

## Forbidden

No threshold relaxation, tuned distance/angle grouping, singular-value
truncation, pressure filter/high-pass, diagonal damping, extra corrector,
post-solve overwrite, case-specific branch, G1/G2/G3, Task 12, Stage 4,
private reference source access, publication or push.

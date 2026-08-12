# Task 11 A22-A3 Exact-Predictor Schur Implementation Plan

> Main-agent-led coupled numerical task. Do not delegate the derivation,
> cross-module implementation, review or acceptance.

**Goal:** Replace the rejected approximate A2 pressure path with the exact
matrix-free Schur complement of the existing A22 momentum and face-predictor
linearization.

**Authority:**
`docs/superpowers/specs/2026-08-04-hundun-flow-task11-a22-a3-exact-predictor-schur-design.md`

## Phase 1: Independent Jacobian RED

Allowed files: tests-only portions of `flow/src/stage3_flow.cpp`,
`flow/src/stage3_flow_test_access.hpp`, and `tests/mpi/test_immersed_piso.cpp`.

1. Add a wished-for `exact_predictor_schur_identity` report containing stable
   GIDs, A22 pressure residual, solved momentum response, face increments,
   divergence and candidate operator action.
2. Extend the authoritative nested equality/helper with ordinary, sign,
   duplicate and nested-size mutations.
3. Compare the literal predictor Jacobian against current A2 and bind a
   semantic RED at the elementwise difference, after structure/state-neutral
   checks pass.

## Phase 2: One private predictor-increment authority

Allowed product file: `flow/src/stage3_flow.cpp`; private finite-volume access
only if the RED proves an unavailable A22 row/gradient query.

1. Factor active-face increment evaluation from `assemble_face_predictor()`.
2. Prove the original predictor output is bitwise unchanged when the helper is
   called with its existing inputs.
3. Cover active-active, periodic, physical pressure outlet and immersed-wall
   rules. Wall flux and velocity remain positive zero.
4. Run focused 1/2/4-rank uniform/warped predictor tests.

## Phase 3: Exact Schur operator

1. Replace `A2PressureOperator` with a private `ExactPredictorSchurOperator`;
   retain no selectable A2 fallback.
2. For each apply, map scaled pressure to the complete A22 pressure residual,
   solve the three actual momentum systems with frozen inner controls, and
   evaluate the factored face increment and its conservative divergence.
3. Evaluate the wall-gradient affine source through the same path.
4. Cache the last accepted apply response only by exact revision and input
   fingerprint; stale values are rejected, never silently reused.
5. Prove constant/reference behavior, parity energy, bilinear symmetry and
   exact independent residual. If literal nonsymmetry exceeds its bound,
   switch the Stage 3 pressure composition to BiCGStab and record why.

## Phase 4: Unified correction

1. Use the A3 solve's momentum and face response directly to update trial
   velocity, face velocity and face mass flux.
2. Remove the separate A2 face correction and separate post-pressure momentum
   correction as product authorities.
3. Retain tests-only recomputation of
   `A delta_u+C_A22 delta_p+B_A22 delta_g=0` and operator-to-face divergence.
4. Cover inner/outer failure, first/second corrector rollback, lowest failing
   rank and bitwise state restoration.

## Phase 5: Direct and numerical gates

1. One-rank exact predictor identity and frozen checkerboard.
2. 1/2/4-rank identity, checkerboard, transaction and decomposition evidence.
3. Focused Debug, Release, ASan, UBSan, tests-off and public-header consumers.
4. Only after direct gates: 12/24 fast localization, 24/48 screen, and one
   final Release 24/48/96 formal matrix.
5. Main-agent complete-diff requirements review, then code-quality review,
   then exact-HEAD Task 11 acceptance.

## Prohibited scope

No public API/schema/Restart/diagnostic change, new runtime dependency,
threshold change, damping, filtering, tuning, third corrector, G1/G2/G3,
Task 12, Stage 4, private legacy source, publication or push. Preserve the
single final Task 11 commit rule and exact DCO.

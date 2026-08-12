# Task 11 A22-A2-P Predictor-Consistency Diagnostic Plan

> Main-agent-only mathematical diagnosis. No worker receives this coupled
> cross-module task.

**Goal:** Prove or disprove that the reassembled PISO predictor omits the A2
interface component of the exact pressure operator before any product repair.

**Authority:**
`docs/superpowers/specs/2026-08-04-hundun-flow-task11-a22-a2p-predictor-consistency-design.md`

**Base / candidate:** accepted Task 10 base
`0db56e463470dd1a605709ba05d8bd6a900f496b`; unchanged repository HEAD
`28f5dd541a0e3ce9ecf852e53d83981add3a5be8` with the preserved dirty Task 11
candidate.

## Task 1: Freeze evidence and write the RED

Files:

- create the A2-P evidence matrix under the Task 11 SDD directory;
- modify tests-only declarations in `flow/src/stage3_flow_test_access.hpp`;
- modify `tests/mpi/test_immersed_piso.cpp`.

Steps:

1. Map `r`, `Kp`, `Ep`, `s_g`, `F_compact` and `F_hybrid` to exact vectors,
   equality semantics, mutations and the one-rank checkerboard gate.
2. Extend the authoritative record comparison/helper first.
3. Require nonempty, correctly sized decomposition vectors and exact
   `hybrid=compact+interface`, `F_compact=r+compact`,
   `F_hybrid=F_compact+interface` identities.
4. Build the focused test and retain the intended failure caused by missing
   tests-only capture, not by environment or unrelated code.

## Task 2: Add the minimum tests-only capture

Files:

- tests-only blocks in `flow/src/stage3_flow.cpp`;
- tests-only record in `flow/src/stage3_flow_test_access.hpp`.

Steps:

1. At the exact corrector coefficient revision, form `q=S p` from the pressure
   before correction.
2. Apply `K_RC` and `K_RC+L_E` independently to the same `q`; reuse the already
   computed affine vector only as read-only evidence.
3. Store all vectors in per-volume scaling and gather them by stable global
   cell ID with the existing correction record.
4. Do not change a product field, solve input, correction, revision or
   corrector count.
5. Run the direct one-rank test and repeated-query/state-neutrality checks.

## Task 3: Diagnose and decide

1. Print signed parity, L2 and Linf for every vector for correctors one and
   two.
2. Run omission, sign and duplication mutation oracles.
3. Add constant/linear/quadratic test-only probes with analytic wall gradients
   if the checkerboard evidence selects the missing-interface hypothesis.
4. Record one of two verdicts:
   - predictor omission proved, with a bounded product-repair addendum; or
   - small repair rejected, with no product edit and an augmented-projection
     proposal.
5. Do not run large grids or acceptance matrices during this diagnosis.

## Prohibited scope

No threshold, case, solver control, corrector count, damping, filtering,
high-pass, product fallback, public interface, G1/G2/G3, Task 12, Stage 4,
private legacy source, publication or push.

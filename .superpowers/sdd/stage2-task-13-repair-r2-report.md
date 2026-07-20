# Stage 2 Task 13 Code-Quality Repair R2 Report

Date: 2026-07-20

Status: DONE

Accepted Task 13 candidate before this repair:
`d19d7372f1fd1bc3f236d131e340a6faf80ba599`

Frozen brief SHA-256:
`26b4aa9a39366d6d837b952a31ee9dee362615513e4cda5f70e34e8347674e91`

Requirements R2 SHA-256:
`dea636ea3a2f409e47390340a950b33576a7d23857c70fc5897946ea066a48ac`

Code-quality R1 SHA-256:
`a1df06ce1716339cc966f4ec20caf774d60849770b9e24baf0e0c9297e19765e`

## Scope

This repair closes only code-quality R1 Important finding I1. It preserves
the public linear API, the accepted CG and BiCGStab recurrences, successful
checkpoint order, operation-count meanings, and the 46/54 successful-path
reduction oracles. It adds no Task 14, external backend, public dependency,
or production device behavior.

## RED evidence

The tests were added before product changes.

1. CG and BiCGStab receive 63 finite right-hand-side entries equal to
   `DBL_MAX/4`. Their local scaled norm sums remain finite while the final
   globally agreed norm overflows. The expected report is
   `non_finite_value/-1`, eight counted reductions, zero matvecs, zero
   preconditioner applications, unchanged solution, and identical reports.
   The coordinator ran the 1/2/4-rank BiCGStab target matrix: all three tests
   failed at the new `lowest_failing_rank == -1` assertion with CTest exit 8.

2. Target ranks receive 63 finite entries equal to `1e154`, while other ranks
   receive `1.0`. The scaled norm remains finite, but a local dot accumulation
   overflows on rank 0 for one rank, rank 1 for two ranks, and ranks 1 and 3
   for four ranks. After the norm-only repair the coordinator reran the same
   matrix: one rank passed and the two-/four-rank cases failed at the expected
   lowest-rank assertion with CTest exit 8. This isolated the dot provenance
   defect.

## Implementation

- Added private `linear::detail::SynchronizedReductionError`, derived from the
  existing project error and carrying only `has_local_source`. It is declared
  in a private implementation header and does not alter public API.
- `VectorOps::norm()` now converts a local invalid view/value or scaled-sum
  failure into a sentinel, completes the matching counted reduction, and then
  throws the structured error. A finite local/global scaled sum followed by a
  non-finite final norm is identified as a globally derived failure with no
  local source.
- `VectorOps::dot_batch()` keeps a fixed scalar
  `first_local_source_pair`; it does not allocate. After its one existing
  reduction it reports provenance for the first non-finite result pair only.
  Thus a later pair's local invalid value cannot be falsely attached to an
  earlier pure-global overflow.
- CG `norm`/`dot` and BiCGStab `norm`/`dot`/`dot_two` catch the typed MPI
  operation error first. For a synchronized numerical-reduction error, every
  rank executes the already-required fixed checkpoint. Ranks with a local
  source submit `non_finite_value`; if no rank has a source, the synchronized
  result is `non_finite_value/-1`. Unexpected vector-operation exceptions use
  `collective_failure`.
- The structured-error path adds no successful-path collective. Global norm
  overflow remains eight reductions. First recurrence-dot overflow remains
  nineteen reductions, with one matvec; CG has one preconditioner application
  and BiCGStab has none at that point.

## Test additions

- The 1/2/4-rank BiCGStab executable now checks both solvers for global norm
  overflow and selected-rank dot accumulation overflow. It checks reason,
  lowest rank, all report fields across ranks, exact counters, no premature
  matvec, current residual policy, and unchanged solution.
- The 2/4-rank VectorOps executable has a two-pair directed test. With a
  pure-global failure in pair 0 and rank-local failure in pair 1, pair 0
  reports no local source. Reversing the pairs makes only rank 1 report a
  local source for the first failed pair. Each call remains one reduction of
  two scalars.
- Existing successful tests still assert 46 BiCGStab reductions for one full
  iteration and 54 CG reductions for two replacement iterations, including
  the no-loop-allocation probes.

## Verification

All commands used the accepted Clang 15/libc++ and OpenMPI environment. MPI
tests were run serially through CTest to avoid shared-fixture interference.

- Debug focused linear/VectorOps matrix: 14/14 passed.
- Debug complete suite after the final per-pair repair: 131/131 passed.
- Release focused matrix: 14/14 passed.
- ASan focused matrix with `ASAN_OPTIONS=detect_leaks=0`: 14/14 passed.
- UBSan focused matrix: 14/14 passed.
- Debug, Release, ASan, and UBSan builds completed successfully.
- `git diff --check` passed.
- CodeGraphF was synchronized after the final source/test changes.

## Changed files

- `linear/src/vector_ops_detail.hpp`
- `linear/src/vector_ops.cpp`
- `linear/src/conjugate_gradient.cpp`
- `linear/src/bicgstab.cpp`
- `tests/mpi/test_bicgstab.cpp`
- `tests/mpi/test_vector_ops_mpi.cpp`
- `.superpowers/sdd/stage2-task-13-repair-r2-report.md`

No private or legacy implementation was accessed. No publication or push was
performed.

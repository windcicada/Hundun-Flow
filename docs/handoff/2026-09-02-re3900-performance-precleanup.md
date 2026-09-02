# Re3900 performance and pre-cleanup reference

Date: 2026-09-02

## Purpose

This note preserves the small amount of reusable diagnostic work found before
removing old branches and worktrees.  It is not a release receipt and does not
authorize weaker numerical tolerances, fewer PISO correctors, reduced
precision, or production debug instrumentation.

## Today-change audit

The audit covered all 193 worktrees registered by Git, including detached
benchmark snapshots.  Of the dirty files, exactly one has a 2026-09-02 mtime:

- detached worktree `/home/wyf/code_dev/hd`, based on `edf752d6`, modified
  `versions/v0.4/src/app_main.cpp` at 09:44:36 +08 to print successful
  pressure--energy globalization trajectories;
- its exact diff is preserved as
  `docs/handoff/reference-patches/2026-09-02-success-globalization-debug.patch`,
  SHA-256 `818cbe9ab6571a61b0849b5743828ae5f5cc4054e2b8368f6ff2600619cdc87c`.

The only commits made today are already retained on
`codex/re3900-10d-smoke`:

- `53e4b059c24769e1e75b3d811c63a6c276288392`, V1 preregistration;
- `edf752d6b9e3e926f099f56cd8564bd269b10f61`, bounded pressure--energy
  refinement rejection, capacity 12, Evidence V7, and candidate identity V2.

No additional commit dated today was found through branch refs, reflogs, or
unreachable-object inspection.

## Other diagnostic material worth preserving

The uncommitted `/home/wyf/code_dev/hp` change adds `MPI_Pcontrol` boundaries
for stages 40/50 and candidate phases 141--143/151--153.  It is based on old
commit `f028631`, uses temporary numeric stage IDs, and does not apply cleanly
to `edf752d6`; it must not be merged as production code.  Its exact diff is
preserved as
`docs/handoff/reference-patches/2026-09-02-pcontrol-stage-markers.patch`,
SHA-256 `72bb9c2fd68d07b9daca52585cd5e33affcb3418b352f7b81e0c585e0a5aef5b`.
It can be manually adapted in a throwaway worktree together with the PMPI
wrapper `/home/wyf/code_dev/.benchmarks/c3900s/mp.c`.  Both saved patches use
zero-context hunks and should be checked/applied with `git apply
--unidiff-zero`.

The deleted instrumentation twin is preserved outside the repository as the
incremental bundle
`/home/wyf/code_dev/hundun-flow-v04-performance-reference-20260902.bundle`,
SHA-256 `71916541a25e0988e864ba8e42c78a91882da70108a3cc28b795949d640ff191`.
It contains `b7ef89cad8a7f48ad35660c95fca0e5cf3654738` and requires base
`67fb5d75f07a7363b418f1f1ba56de6d82c6b5bd`.  The bundle is reference-only:
the profiler has a four-step buffer ceiling, adds an unconditional elapsed-time
reduction in its old driver integration, and needs naming/API cleanup before
reuse.

The functional hot-path changes from that twin are already represented on the
retained mainline by `340eb2b`, `dcd41c1`, and `ba1eb00`; they are not missing
and do not need a second copy.

The minimal performance evidence set is also archived outside the repository
as `/home/wyf/code_dev/hundun-flow-v04-performance-evidence-20260902.tar.gz`,
SHA-256 `707ac58b2e28ea7ee408e30a5e7caa3d94923e4fe9dc751506d635683b2e3f2e`.
`gzip -t` passed and the archive contains the current 50-step/noisy smoke
evidence, the periodic/symmetry A/B, the September 1 PMPI/root-cause files,
the core-hotpath receipt, and the frozen historical paired statistics.

## Performance evidence

### Historical paired baseline

`formal_full20_h2c234b1_c3c22e0f_r1_03/PAIRED-05-STATS.json` records the
frozen `480x480x48 / 64-rank` comparison, or 172,800 cells/rank:

- COAST P90: 15.419435 s/step, equivalent to 4.461642 s at 50,000
  cells/rank under simple cell-count normalization;
- HUNDUN P90: 4.418378 s/step, equivalent to 1.278466 s at 50,000
  cells/rank;
- five-pair decision: `ACCEPT`.

The old problem was a `20D x 20D x piD` slip-span case.  It is a useful speed
anchor, not an equal-work measurement for the present periodic thin-span
problem.

### September 1 root-cause benchmark

`/home/wyf/code_dev/.benchmarks/c3900s/cause.md` records a
`320x320x32 / 128-rank` comparison:

- fresh BDF2 HUNDUN median 4.767824 s versus COAST 1.192407 s;
- developed HUNDUN median 3.927224 s versus COAST 1.185999 s;
- pressure stages account for 87.056% of developed HUNDUN time;
- a typical step performs 4,420 halo exchanges; four success-path control
  reductions per exchange imply 17,680 reductions that were absent from the
  old evidence counters.

The PMPI profile assigns the pressure stages 18,833 `MPI_Allreduce` and 278
`MPI_Allgatherv` calls/rank.  Payload wait is small; fine-grained global status
synchronization is the dominant MPI problem.

### Core-hotpath preservation check

`/home/wyf/code_dev/.benchmarks/pa/final.md` and run `pa/runs/x9` prove that
`ba1eb00` retained the hot-path work:

- `480x480x48 / 64 ranks`, step 2: 18.605899 s;
- two pressure solves, zero pressure--energy refinement solves;
- 74 total linear iterations;
- stages 40+50: 16.037301 s.

The current branch descends from the equivalent mainline changes.  Normalized
by cell count and pressure-like iterations, the old run costs about
1.4278 s/(million cells x iteration).

### Current 20D x 10D x 3D result

The stable same-binary run is `/home/wyf/r39m/c12/x50/evidence.jsonl`:

- candidate `edf752d6`, binary
  `e45f26c32168ccb5500793b71433e8cf95c541e905c6e64a7a0b70263d588d0d`;
- `456x256x104 / 128 ranks`, 94,848 cells/rank;
- 49 BDF2 rows have median 34.809398 s and interpolated P90 35.755134 s;
- cell-count normalization gives 18.350096 s at 50,000 cells/rank;
- median pressure--energy refinement count is six;
- step 2 uses two ordinary pressure solves plus five refinement solves,
  243 linear iterations, 6,046 recorded blocking collectives, and 6.781 s of
  recorded reduction time.

Across BDF2 rows, refinement count and step time correlate at 0.9951.  A
linear fit gives about 3.94 s per additional refinement solve.  Groups with
3/4/5/6/7 refinements have median step times of approximately
23.29/27.24/31.07/35.00/39.40 s.

The 10:38 two-step smoke is not the timing baseline.  It has identical
candidate, binary, case, CPU plan, iterations, refinements, messages, and
refills, but step 2 took 46.599512 s and recorded reductions grew from
6.781 s to 15.498 s.  This is a 50.6% same-work slowdown consistent with
machine/runtime noise.  Its numerical success remains valid.

For the cleaner current step 2, pressure work costs about
1.3120 s/(million cells x pressure-like iteration), slightly better than the
old `ba1eb00` unit-work value.  The large wall-time regression is therefore
caused primarily by pressure work volume increasing from about 65 iterations
over two solves to about 234 iterations over seven solves, not by loss of the
September 1 hot-path optimizations.

### Periodic/symmetry diagnostic A/B

The current-binary two-step A/B under `/home/wyf/r5` uses
`384x128x32 / 128 ranks`.  The two case JSON files differ only in the z-min
and z-max boundary kinds:

- periodic step 2: 7.390551 s, five refinements, 62 total linear iterations;
- symmetry step 2: 1.996615 s, zero refinements, 26 total linear iterations.

The candidate identity and executable are identical.  The case/product
fingerprints and adaptive-order IBM reconstruction counts correctly differ
because the boundary contract is part of the compiled product.  This result
strongly selects the periodic coupling path as the source of repeated
component-residual refinement, but changing the formal physical boundary is
not a valid performance fix.

## Diagnosis and next measurement

1. Keep the 50-step run as the current timing reference; do not use the noisy
   10:38 smoke as a performance claim.
2. Use the saved successful-trajectory patch in a throwaway worktree to learn
   why the thin periodic case repeatedly needs five to seven nonlinear
   refinements.
3. If finer attribution is required, adapt the saved PControl patch and PMPI
   wrapper for one short, isolated run.  Do not merge the instrumentation.
4. Optimize by reducing or reusing exact refinement work and by batching the
   remaining MG/Halo control authority.  Preserve the two PISO correctors,
   FP64 residual contracts, EOS/continuity/energy gates, and rollback
   semantics.

## Review notes for today's commits

- `edf752d6` is correctness/evidence work, not a speed optimization.  Raising
  the bound from 6 to 12 permits expensive accepted paths; it does not itself
  force extra solves.
- `app_identity_detail.hpp` uses `std::string_view` without directly including
  `<string_view>`.
- changing public report array capacity from 6 to 12 changes C++ layout; the
  release documentation should state the rebuild/ABI boundary.
- the V1 policy still lacks a wall-time gate for the new medium case.  The
  historical paired policy cannot silently be reused because the domain,
  periodic span, rank count, time step, and nonlinear work differ.

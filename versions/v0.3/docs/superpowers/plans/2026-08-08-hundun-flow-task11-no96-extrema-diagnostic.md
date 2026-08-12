# Task 11 no-96 extrema disambiguation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Permanently remove 96-cubed Task 11 solves and add a 12/24/48 top-eight extrema diagnostic comparing per-level and fixed-48 surface fixtures.

**Architecture:** Extend the test-support manufactured-case definition with a test-only surface policy and top-eight pressure-extrema detail. Keep the production solver untouched. Register one diagnostic selector and change retained acceptance schedules to 12/24/48 only.

**Tech Stack:** C++17 test support, existing HUNDUN finite-volume/immersed path, CMake/CTest, MPI rank-one Release runner, systemd user units for detached diagnostics.

## Global Constraints

- No Task 11 command, CTest entry, systemd unit, or diagnostic runner may launch a 96-cubed flow solve.
- Existing 96-cubed logs are historical evidence only and are not rerun.
- Do not change product solver rows, thresholds, PISO correctors, filters, damping, or gauge semantics.
- Use `apply_patch` for manual edits and preserve all inherited dirty files.
- Bind every detached run to exact HEAD, binary SHA-256, command, environment, log, and exit status.

---

### Task 1: Add the failing top-eight extrema oracle

**Files:**
- Modify: `tests/support/stage3_mms.hpp`
- Modify: `tests/numerical/test_laminar_ibm_order.cpp`
- Test: `test_laminar_ibm_order oracle`

**Interfaces:**
- Add `select_pressure_error_extrema(candidates, limit)` returning a
  deterministic absolute-error-descending vector with global-cell-id ties.
- Keep the existing single-extremum helper and oracle behavior unchanged.

- [x] Add the oracle assertions for empty input, limit zero, tie ordering,
  nonfinite rejection, and truncation to the requested limit.
- [x] Run the oracle and verify it fails because the new helper is absent.
- [x] Add the minimal declaration/implementation needed by the oracle.
- [x] Run the oracle again and verify it passes.

### Task 2: Add surface policy and local extrema details

**Files:**
- Modify: `tests/support/stage3_mms.hpp`
- Modify: `tests/support/stage3_mms.cpp`

**Interfaces:**
- Add `ManufacturedSurfacePolicy { per_level, fixed_48 }` to
  `ManufacturedCase`.
- Add a test-only `PressureErrorExtremumDetail` and a top-eight vector to
  `ManufacturedRunResult`.
- `run_manufactured_case` uses `1/cells` for `per_level` and `1/48` for
  `fixed_48`; all other flow setup is unchanged.

- [x] Add a test-only diagnostic flag to `ManufacturedCase` and RED assertions
  that the default policy remains `per_level`.
- [x] Collect and deterministically sort top-eight pressure extrema on the
  rank-one diagnostic path.
- [x] For each selected cell, compute wall-distance-over-`h`, incident link
  count, nearest authority-link donor count, condition estimate, and donor
  fingerprint.
- [x] Print no diagnostic detail for ordinary formal selectors unless the
  flag is enabled.

### Task 3: Register the 12/24/48 diagnostic and remove 96 schedules

**Files:**
- Modify: `tests/numerical/test_laminar_ibm_order.cpp`
- Modify: `CMakeLists.txt`

- [x] Add `pressure_extrema_screen` using grids `{12, 24, 48}` for both
  `per_level` and `fixed_48` policies.
- [x] Print `L2`, near-wall `L2`, global `Linf`, top-eight records, and the two
  observed orders for each policy.
- [x] Change retained acceptance selectors to `{12, 24, 48}` only.
- [x] Remove any new or active Task 11 schedule that contains a 96-cubed flow
  command; leave historical logs untouched.
- [x] Add a CTest entry with a bounded screen timeout and no `RUN_SERIAL`
  conflict with unrelated unit tests.

### Task 4: Verify the diagnostic in increasing cost order

**Files:**
- Create: `.superpowers/task-11-no96-extrema-disambiguation-2026-08-08.md`
- Create: `/tmp/hundun-task11-no96-extrema/manifest.txt` (runtime evidence)

- [x] Build the Release diagnostic binary and record its SHA-256.
- [x] Run the oracle and the new 12/24/48 diagnostic once per surface policy.
- [x] Verify no command line, process argument, or log contains a 96-cubed
  flow launch.
- [x] Compare top-eight cell movement, wall-distance phase, donor metadata,
  and observed orders; write the decision and any remaining blocker.
- [ ] Run retained 12/24/48 acceptance selectors only after the diagnostic
  completes.

### Task 5: Review and handoff

**Files:**
- Review: all files changed by Tasks 1–4
- Update: `docs/superpowers/plans/2026-08-07-hundun-flow-task11-acceptance-cluster-reduction.md`

- [x] Confirm the plan/spec contain no executable 96-cubed schedule.
- [x] Confirm product sources are unchanged by the diagnostic work.
- [x] Record deferred six selectors and the A22 follow-up decision.
- [ ] Keep the worktree index empty unless the user explicitly requests a
  commit; preserve all inherited evidence and logs.

### Task 6: Resolve the 12-to-24 near-wall pressure row

**Files:**
- Modify: `tests/support/stage3_mms.hpp`
- Modify: `tests/support/stage3_mms.cpp`
- Modify: `tests/numerical/test_laminar_ibm_order.cpp`
- Create: `.superpowers/task-11-near-wall-pressure-fixed-band-2026-08-08.md`

- [x] Partition the moving `0--2h` diagnostic into four half-cell bands and
  incident/non-incident cell families, including mean, RMS and centered RMS.
- [x] Prove that the immersed incident family is second-order (`1.925398`)
  while the moving non-incident support causes the aggregate `1.601264` row.
- [x] Reconcile the implementation with the controlling requirement for a
  fixed physical near-wall thickness.
- [x] Freeze the formal support at `2*(L_ref/12)=L_ref/6` without changing
  the solver, gauge, threshold or PISO path.
- [x] Verify the corrected 12-to-24 pressure order is `1.989545`.

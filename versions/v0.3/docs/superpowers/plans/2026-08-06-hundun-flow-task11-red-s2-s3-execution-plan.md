# Task 11 RED-S2/RED-S3 execution plan

Status: frozen for execution on 2026-08-06 after the RED-S1 aggregate
acceptance. This plan is scoped by the science-first execution amendment and
the Task 11 science-closure RED design. It does not accept Task 11 or Stage 3.

The 2026-08-08 no-96 amendment supersedes the final 96-cubed acceptance
wording below. All active Task 11 screens and retained acceptance selectors
use 12/24/48 only; historical 96-cubed runs are not repeated.

## Authority and inputs

- RED-S1 accepted at `29cc3e1` (aggregate record:
  `.superpowers/task-11-red-s1-aggregate-acceptance-2026-08-06.md`).
- RED-S2/RED-S3 definitions:
  `docs/superpowers/plans/2026-08-05-hundun-flow-task11-science-closure-red-design.md`
  sections 4 and 5.
- Upstream read-only cross-check (2026-08-06 user guidance): Basilisk
  `embed.h` force/flux algebra, AMReX EBFluxRegister conservation, gslib
  `findpts` global donor lookup. GPL code is algorithm-reference only; BSD
  implementations are referenceable. No source is copied.

## Facts already established (not re-investigated)

1. The decomposition RED was a test-support lookup bug, not a product solve
   failure: `MeshTopology::find_local_cell` knows only face-adjacency ghosts
   and `ImmersedDomain::links()` returns the rank-local owned subset, while
   the product resolves all row/wall donors from halo-exchanged fields with
   the fixed reach-four ghost width.
2. `QuadraticReconstruction::create` already resolves donors from global
   coordinates (no local topology lookup) and already enforces the reach-four
   cap with an explicit construction-time error. The product shared-row
   authority is therefore global by construction.
3. `tests/support/stage3_mms.cpp` now evaluates exact donor cell averages and
   the wall-functional `h` from global coordinates and gathers the global
   link->fluid-cell map, matching the product's global donor semantics.
4. The single pressure-authority chain (row reconstruction -> wall pressure
   gradient -> final momentum residual and final force) is already one
   source; no second authority is introduced.

## RED-S2 closure scope

### A. Product-side explicit donor-owner audit (no numerical change expected)

Files (read-only audit first, then edit only if a gap is proven):

```text
immersed/src/quadratic_reconstruction.cpp
immersed/src/ghost_stencil_plan.cpp
finite_volume/src/immersed_operator.cpp
finite_volume/src/immersed_boundary_authority_detail.hpp
```

Audit item: every shared-row donor must have an explicit deterministic global
identity and the row's reach must be diagnosed at construction time, not
rejected late. Current state: `QuadraticReconstruction::create` throws
"donor exceeds maximum halo reach four" at construction; row evaluation reads
the exchanged field by global id. If the audit finds any path that implicitly
assumes a donor is local (topology lookup, field-index without ghost check,
or per-link fallback), fix it with the smallest edit and a mutation RED.

### B. Direct shared-row mutation RED

Primary host: `tests/mpi/test_immersed_operator.cpp` (already contains shared
row authority probes vs the legacy per-link reconstruction and interface-row
snapshot mutations). Execute the seven product mutants one at a time in
disposable candidates and require death through designated observations:

| ID | Mutation | Designated failing oracle |
|---|---|---|
| S2-M1 | `build_boundary_row_reconstruction` keeps only first-link donors | shared-row donor union/fingerprint and value/gradient probes |
| S2-M2 | pressure/viscous affine plan uses `ghost_plan.reconstruction(link.id)` | shared-row polynomial value/gradient vs legacy |
| S2-M3 | complete A22 row omits diagonal or neighbour defect | exact polynomial residual + source-term inventory |
| S2-M4 | final boundary row omits or repeats background removal | before/background/removed/wall/after identity |
| S2-M5 | viscous wall derivative reuses stale/different authority | exact shared-row derivative with non-polynomial probe |
| S2-M6 | row execution evaluates the affine row twice | exact evaluation/write counters and final vector |
| S2-M7 | link ordering or rank decomposition changes assembly | canonical donor IDs/fingerprint, 1/2/4-rank bitwise |

### C. Regression gates

- Release: all five decomposition fast configurations plus the focused
  26-test set (already 26/26 after the test-support fix).
- Debug: all five decomposition fast configurations on 1/2/4 ranks and the
  focused contracts.
- 1-rank manufactured smoke must stay bitwise `4647cc59…`.

## RED-S3 closure scope

### D. Exact-state force-consistency decomposition selector

Files:

```text
tests/support/stage3_mms.hpp
tests/support/stage3_mms.cpp
tests/numerical/test_laminar_ibm_order.cpp
CMakeLists.txt
```

Add `exact_force_consistency_fast` (12^3 sphere-uniform and warped-prism
shared-row viscous discrimination) computing, for each force part:

```text
C_exact_h  = O_h(U*) - S_h(U*)
C_solved_h = O_h(U_h) - S_h(U_h)
C_state_h  = [O_h(U_h)-O_h(U*)] - [S_h(U_h)-S_h(U*)]
C_solved_h = C_exact_h + C_state_h   (bookkeeping closure)
```

Register as `test_laminar_ibm_exact_force_consistency_fast`,
labels `numerical;stage3;task11;fast`. Run the decision tree:

- `C_exact` first order: reconstruction/geometry/operator/surface authority
  is the root; whitelist expands to the immersed/ghost/quadratic files under
  a separate mutation RED.
- `C_exact` second order but `C_solved` first order: pressure/PISO/state
  coupling is the root; whitelist expands to `flow/src/stage3_piso_detail.hpp`
  and `flow/src/stage3_flow.cpp` under a separate mutation RED.
- pressure passes but viscous fails: shared viscous row/derivative only.

### E. 2D force-algebra cross-validation oracle (test-only)

Add a fast-layer oracle that compares the product integrated pressure/viscous
force against the Basilisk `embed_force` formula evaluated term-by-term on a
2D slice, using only product geometry/normal/weights and analytic stress.
It must not implement a replacement integrator in the product path.

### F. RED-S4 solved fast/screen

Only after D/E classify the root cause and the minimal repair is green, run
the existing `functional_selection_fast` and `functional_closure_fast`
selectors (12/24 Release) and then the 12/24/48 extrema screen. No formal
96^3 acceptance is scheduled under the no-96 amendment.

## Mutation and evidence protocol

- One disposable candidate per mutant; fresh build; serialized MPI.
- Every run binds HEAD, seed hashes, binary SHA-256, command, environment,
  log SHA-256 and exit status.
- No threshold, PISO corrector, filter, damping, per-case tuning, or
  alternative product algorithm is permitted.
- The main worktree index stays empty; product edits land only after their
  RED is observed and only through the frozen file whitelist.

## Rollback boundary

- Pre-RED-S2 immutable snapshot: `/tmp/hundun-task11-red-s1-green-before.IIsVxa`
  plus the current whole-tree hash recorded in the RED-S1 evidence.
- Test-support edits made so far (`tests/support/stage3_mms.cpp`) are
  revertible by restoring the snapshot file; the Release/ Debug evidence
  above is bound to the edited file.
- Any product edit is made only in disposable mutation candidates until its
  RED is accepted; the live product files are not touched by mutations.
- Docs-only commits record acceptances; product consolidation waits for the
  complete Task 11 candidate review.

## Completion boundary

RED-S2 is accepted when: the decomposition fast matrix is green on 1/2/4
Debug and Release; the direct shared-row RED and S2-M1..M7 are accepted; and
independent review passes. RED-S3 is accepted when the exact-state
decomposition classifies the root cause and the minimal repair passes its
RED/mutations plus the fast selectors. Neither accepts Task 11; formal
acceptance and the Task 11 verdict remain gated on RED-S4 and the full
matrix.

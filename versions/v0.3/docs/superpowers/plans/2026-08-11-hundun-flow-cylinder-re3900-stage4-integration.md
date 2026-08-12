# HUNDUN-FLOW Stage 4 Re=3900 Cylinder Integration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fuse the accepted Stage 4 tree with the periodic-cylinder work and make the frozen Re=3900 cylinder case complete its functional and full-grid startup gates.

**Architecture:** Develop from the Stage 4 governance seal, preserve the old dirty benchmark worktree, and port the already reviewed periodic IBM patch. Generalize the accepted compact `GhostedVectorHalo` to sparse active layouts, replace active-operator runtime global gathers, then use nested solve provenance to linearize the exact Schur response, select restarted GMRES for its nonsymmetric operator, and keep positive wall transport coefficients inside their frozen donor range.

**Tech Stack:** C++17, CMake, MPI-3, HUNDUN CPU reference execution, compact Buffer Halo, CG/restarted GMRES, systemd user units, Python only for maintainer-side benchmark generation and analysis.

## Global Constraints

- Parent is Stage 4 seal `033a685c90c1f9c674e93a4b82db10db4c381abe`; tested Stage 4 code is `6407cd7c591ce088db7f1dd7e296d77acd18da1c`.
- Preserve exact Schur, exactly two PISO correctors, force signs, scientific thresholds and the frozen formal cylinder case.
- Do not add damping, filtering, case-specific parameters, extra correctors, 96³, vendor solver dependencies or private source.
- Do not touch the old dirty benchmark worktree, product `main`, COAST research tree, Stage 5, remote repositories or research processes.
- Long runs start only after fast and functional gates; at most one high-memory MPI job runs at a time.

---

### Task 1: Freeze and Port the Stage 3 Benchmark Delta

**Files:**
- Modify: `.gitignore`
- Modify: `docs/numerics/stage3-contracts.md`
- Modify: `include/hundun/mesh_topology.hpp`
- Modify: `src/mesh_topology.cpp`
- Create: `src/ib_periodic_surface_window_detail.hpp`
- Modify: `src/ib_domain.cpp`
- Modify: `src/ib_ghost_stencil_plan.cpp`
- Modify: `src/ib_quadratic_reconstruction.cpp`
- Modify: `src/ib_deterministic_qr_detail.hpp`
- Modify: `src/fvm_immersed_operator.cpp`
- Modify: `src/flow_immersed.cpp`
- Modify: `src/app_immersed_flow_driver.cpp`
- Create: `benchmarks/cylinder_re3900/cases/hundun/**`, `benchmarks/cylinder_re3900/geometry/**`, and C++ periodic tests only
- External maintainer tools: `/home/wyf/code_dev/.benchmarks/cylinder-re3900-stage4/maintainer-python/**`
- External private comparison evidence: `/home/wyf/code_dev/.benchmarks/cylinder-re3900-stage4/private-comparison-snapshot/**`

**Interfaces:**
- Consumes: accepted Stage 4 source tree and the exact dirty delta in `/home/wyf/code_dev/.worktrees/hundun-flow-cylinder3900`.
- Produces: a reviewable Stage 4-based periodic IBM benchmark tree without modifying the source worktree.

- [x] **Step 1: Record both source identities and dirty path inventory.**
- [x] **Step 2: Apply the tracked product diff and mechanically copy only listed untracked benchmark assets; move Python tools/tests to the external maintainer root before configure.**
- [x] **Step 3: Run `git diff --check`, two periodic C++ REDs, the 16 external Python benchmark tests and 64-rank `--validate`.**
- [x] **Step 4: Run the affected Stage 4 source-policy, schema and IBM focused tests.**
- [x] **Step 5: Commit the fusion checkpoint with DCO after main-agent full-diff review.**

### Task 2: Generalize Compact Buffer Halo to Sparse Layouts

**Files:**
- Modify: `include/hundun/lin_ghosted_vector_halo.hpp`
- Modify: `src/lin_ghosted_vector_halo.cpp`
- Modify: `src/lin_ghosted_vector_halo_detail.hpp` only if a test seam is required
- Modify: `tests/mpi/test_ghosted_vector_halo.cpp`

**Interfaces:**
- Consumes: `StructuredDecomposition`, `MeshTopology`, `ExecutionContext`, and a `VectorLayout` ordered as owned then ghosts.
- Produces: `GhostedVectorHalo::create(decomposition, topology, context, layout)` with the existing begin/wait/exchange and counter contracts.

- [x] **Step 1: Add a 1/2/4-rank sparse-layout RED whose active IDs omit interior solid cells and whose ghost values equal a canonical-ID oracle.**
- [x] **Step 2: Verify RED fails because the sparse-layout overload is absent.**
- [x] **Step 3: Refactor plan construction to consume an explicit layout while preserving the existing overload as a wrapper.**
- [x] **Step 4: Add mutation REDs for duplicate ID, locally owned ghost, missing/inactive requested ID and collective rank disagreement.**
- [x] **Step 5: Verify old full-layout Halo tests and new sparse tests pass at 1/2/4 ranks.**
- [x] **Step 6: Commit the independently reviewable sparse-Halo change with DCO.**

### Task 3: Replace Active Runtime Global Gathers

**Files:**
- Modify: `src/flow_immersed.cpp`
- Modify: `src/flow_immersed_test_access.hpp` if existing test access needs counters
- Modify: `tests/mpi/test_active_boundary_layout.cpp`
- Modify: `tests/mpi/test_immersed_piso.cpp`

**Interfaces:**
- Consumes: sparse `GhostedVectorHalo`, active `VectorLayout`, active pressure/momentum connections.
- Produces: active operator apply and pressure-mobility diagonal exchange with peer-only runtime payloads.

- [x] **Step 1: Extend the existing small-grid gathered-reference oracle to cover the active exchange path across 1/2/4 ranks.**
- [x] **Step 2: Add a counter mutation proving runtime active apply payload is bounded exactly by ghost count.**
- [x] **Step 3: Give each active operator a `GhostedVector` workspace and sparse Halo; bind every connection to a local layout offset.**
- [x] **Step 4: Exchange three momentum diagonals through the same sparse layout before pressure mobility assembly.**
- [x] **Step 5: Remove runtime global-value buffers/counts/displacements without changing coefficient or reduction order.**
- [x] **Step 6: Run active operator, signed-force, decomposition, rollback and exactly-two-corrector focused tests.**
- [x] **Step 7: Commit the active-exchange change with DCO.**

### Task 4: Preserve Nested Linear-Solve Failure Provenance

**Files:**
- Modify: `src/flow_immersed.cpp`
- Modify: `src/app_immersed_flow_driver.cpp`
- Modify: `tests/mpi/test_immersed_piso.cpp`
- Modify: `tests/mpi/test_adaptive_time_control_diagnostics.cpp` if public reporting changes

**Interfaces:**
- Consumes: inner momentum and outer pressure `SolveReport` values.
- Produces: deterministic failure phase/component/report diagnostics without changing the solve controls.

- [x] **Step 1: Add a mutation RED that forces one exact-response momentum component to hit maximum iterations and requires phase, component, reason, iterations and residuals.**
- [x] **Step 2: Add a second RED distinguishing an outer-pressure failure from an independent-residual rejection.**
- [x] **Step 3: Implement internal attempt-local provenance and include it in rank-0 driver failure output.**
- [x] **Step 4: Verify rollback, retry identity and successful report semantics remain unchanged when no failure occurs.**
- [x] **Step 5: Run the 48³ case once and seal the exact nested failure report.** Evidence: `runs/hundun-fast48-f12d67f-20260811T080200Z` records first-corrector `pressure_outer`, 500 iterations and divergent independent residual on exact HEAD `f12d67f`.

### Task 5: Repair the Proven Solver Layer

**Files:**
- Create: `include/hundun/lin_restarted_gmres.hpp`
- Create: `src/lin_restarted_gmres.cpp`
- Modify: `src/CMakeLists.txt`
- Modify: `src/app_immersed_flow_driver.cpp`
- Create: `tests/unit/test_restarted_gmres_header_contract.cpp`
- Create: `tests/mpi/test_restarted_gmres.cpp`
- Modify: `tests/CMakeLists.txt`
- Modify: `src/flow_immersed.cpp` only if the frozen outer-failure RED exposes an operator-contract defect.
- Modify: `src/ib_quadratic_reconstruction_detail.hpp`, `src/ib_quadratic_reconstruction.cpp`, `src/ib_wall_force.cpp` and `tests/mpi/test_wall_force.cpp` only if the successful solve exposes a positive-coefficient wall reconstruction defect.

**Interfaces:**
- Consumes: Task 4 exact failure phase and mutation RED.
- Produces: a general, fixed-algorithm preconditioner or solver-contract repair that preserves the exact operator.

- [x] **Step 1: Freeze a non-normal nonsymmetric matrix RED and the exact `f12d67f` outer-pressure failure before product changes.** The first 24³ GMRES screen proved a nonlinear open-inlet matvec: recursive residual fell to `0.0711` while the true residual remained `18.94`.
- [x] **Step 2: Preserve CG for SPD momentum; do not use CG for the exact Schur because the accepted IBM oracle proves that operator is nonsymmetric.**
- [x] **Step 3: Add CPU-reference restarted right-preconditioned GMRES with two-pass modified Gram-Schmidt, Givens least squares and true-residual recomputation at every restart/claimed convergence.** The generic test covers 1/2/4 ranks, zero RHS, invalid restart and a genuine short-restart cycle.
- [x] **Step 4: Keep compact Schur/Jacobi as the single pressure preconditioner authority and switch only the driver pressure solver from BiCGStab to GMRES.** The homogeneous exact response is now `F(p,0)-F(0,0)` while the affine RHS retains `F(0,g)`; a nonzero-inlet 1/2/4-rank RED proves zero action and superposition. Wall-force `mu_eff` uses the same frozen donors with a donor-extrema bound, preventing a positive coefficient field from producing a negative quadratic wall value.
- [x] **Step 5: Run fast 12/24 diagnostics and the frozen 48³ 64-rank one-step gate.** 12³ remains structurally invalid for this geometry and is not used. The 24³/8-rank product path passed in `screens/wale-bounded-mu-gmres64-precommit-n5Um35`: two correctors, continuity `9.33e-15`, pressure residual `9.30e-13`, finite four-field force and finite WALE summary. Exact HEAD `22487ab` with one compact Jacobi application is the first mutation RED: `runs/hundun-fast48-22487ab-cWYiMB` rejects at `pressure_outer` after 500 iterations. Commit `7e240b5` replaces that non-scalable application with a fixed compact-Poisson inexact-inverse contract consumed by FGMRES and passes 24³/48³ screens plus the exact-HEAD 48³ run `runs/hundun-fast48-7e240b5-i9YapF`; the latter commits one 64-rank step with continuity `6.69e-14`, pressure residual `1.99e-11` and complete force fields in 32.29 s. The frozen 480×480×48 startup then supplies three more mutation REDs: `runs/hundun-fullgrid-inner-provenance-diagnostic-xyvdrz` proves the compact immersed operator is not safely CG-solvable (`pressure_compact_preconditioner`, 200 iterations, residual `1.0000 -> 1.2100`); rank-local IC(0) alone still fails to reduce that global mode in `runs/hundun-fullgrid-compact-icc-precommit-t23lXE`; and the two-level/inner-FGMRES candidate in `runs/hundun-fullgrid-two-level-r12-precommit-ggfnNy` reaches 243.025 GiB aggregate rank RSS because restarted FGMRES allocated its large Buffer workspace on every nested solve. The general repair keeps the exact operator, thresholds and inexact-inverse contract, uses rank-local IC(0) plus a rank-constant coarse correction, and gives the already accepted FGMRES implementation transactional persistent workspace. The allocation RED proves that a repeated same-layout solve performs zero Buffer allocations at 1/2/4 ranks and that failed layout replacement preserves the old workspace. The repaired product path passes `screens/persistent-fgmres24-precommit-yHchtF` in 7.05 s and `screens/persistent-fgmres48-precommit-FEwqdZ` in 33.45 s; exact candidate `9c18c7a` then passes the 48³ tests-on and tests-off gates with bitwise-identical product output in 33.20 s and 32.08 s. Its frozen 480×480×48 tests-off startup commits one step in 26:51.50 with continuity `7.50e-15`, pressure residual `1.78e-11`, complete four-field force diagnostics and bounded aggregate rank RSS.
- [x] **Step 6: Review the mathematical derivation, mutation sensitivity and complete diff; commit with DCO.** The final candidate keeps the exact Schur authority, uses the compact operator only as an inexact inverse, and changes no scientific threshold, corrector count or cylinder case input. The complete diff review covered public enum ordinals, MPI ownership, IC(0)/coarse algebra, nested failure provenance, persistent-workspace rollback, caller impact, allocation lifetime and private-source independence.

### Task 6: Cylinder Functional and Full-Grid Startup Acceptance

**Files:**
- Modify: `docs/benchmarks/cylinder-re3900-design.md`
- Create: `.superpowers/sdd/cylinder-re3900-stage4-acceptance.md`
- Runtime evidence: `/home/wyf/code_dev/.benchmarks/cylinder-re3900-stage4/**`

**Interfaces:**
- Consumes: frozen candidate, binary, cases and COAST one-step evidence.
- Produces: exact-HEAD ACCEPT/REJECT and permission boundary for later long statistics.

- [x] **Step 1: Run 48³/64-rank one complete step with two PISO correctors and finite force fields.**
- [x] **Step 2: Run the 480×480×48/64-rank minimum full-grid startup through one committed step; record time, RSS, solve counters, Cd and Cl.**
- [x] **Step 3: Run complete affected Debug, focused Release, small ASan/UBSan, Stage 4 regression and 1/2/4-rank decomposition gates.**
- [x] **Step 4: Freeze HEAD/tree/diff/binary/command/environment/log hashes and verify no benchmark MPI process remains.**
- [x] **Step 5: Main agent performs requirements, science, caller-impact, copyright and complete-diff review and records ACCEPT or REJECT.**

## Plan self-review

- Coverage: Stage 4 intake, periodic surface semantics, scalable exchange, nested failure diagnosis, evidence-driven solver repair, force output and full-grid startup each have one owner.
- Placeholders: Task 5 file scope is evidence-selected but bounded; it does not authorize an unspecified algorithm or threshold change.
- Type consistency: all active exchanges use the existing `VectorLayout`, `GhostedVector` and `GhostedVectorHalo` contracts; no second vector identity is introduced.

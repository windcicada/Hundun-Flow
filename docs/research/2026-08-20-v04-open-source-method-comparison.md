# HUNDUN-FLOW v0.4 Cartesian / Low-Mach Plan: Primary-Source Comparison

Date: 2026-08-20
Reviewed plan: `docs/superpowers/plans/2026-08-12-hundun-flow-v0.4-cartesian-performance-architecture.md`
Scope: OpenFOAM pressure--velocity/thermophysical sequencing, AMReX storage/halo/EB/MLMG/I/O, AMReX-Hydro and IncFlo projections, and HYPRE/PETSc linear-solver lifecycle. This review deliberately excludes AMR, body-fitted meshes, nonmatching refinement, product GPU work, SIMPLE/PIMPLE implementation, and any proposal to copy upstream source.

## Executive decision

**REVISE before implementing the complete Task 15 pressure-coupling path.** The plan's main performance architecture is sound: cold compilation, padded SoA storage, static Cartesian/IBM plans, explicit halos, one owned flux publication path, persistent operator/workspace objects, and transaction-scoped state are all worth keeping. The plan is not yet closed in four places that directly affect correctness or measured performance:

1. The current hot schedule publishes the current-step final pressure/flux before advancing enthalpy and species. Those later updates change `rho(p_abs,h,Y)`, so the terminal EOS state can disagree with the density and `drho/dp` used by the pressure equation. OpenFOAM's compressible transient ordering instead advances thermophysics before entering its PISO pressure-corrector loop; AMReX-Hydro's two projections are also placed on opposite sides of the advective state update, not implemented as two adjacent copies of one pressure correction.
2. The multigrid plan does not cleanly separate immutable coarsening topology, exact operator numeric coefficients at every level, and reusable/stale-allowed preconditioner setup. A coefficient threshold may control preconditioner setup reuse, but it may not leave the exact operator or its required coarse numeric representation ambiguous.
3. The stated one-shot freeze order has a dependency cycle: arena and communication capacity depend on graph liveness/ghost analysis, while the executable graph needs bound arena views and communication objects. It needs an explicit analyze-then-bind two-phase compiler. Numeric/hierarchy storage can be allocated at freeze, but it cannot be certified as numerically valid before initial fields and coefficients exist.
4. The performance contract counts point-to-point messages and bytes but not reductions/collectives, and the comparison uses unpaired global medians. For a pressure-based MPI code, collective count and slowest-rank wall time are first-class costs. The release protocol should use max-rank timed regions and paired HUNDUN/COAST ratios under alternating launch order.

These revisions do **not** require AMR, PIMPLE, a third pressure corrector, body-fitted grids, or case-specific product code.

## Decision summary

| Area | Decision | Plan action |
|---|---|---|
| Immutable case/geometry/boundary/operator plans | **KEEP** | Preserve cold compilation and revision identities. |
| Padded x-fast SoA and preallocated arenas | **KEEP**, minor **REVISE** | Add parallel NUMA first-touch and per-thread false-sharing tests; do not imitate AMReX's per-call allocations. |
| `rAU/rAtU/HbyA/phiHbyA` dependency model | **KEEP**, **REVISE identities** | Include density/thermo/`drho_dp`, pressure reference/gauge, flux-history and numeric-BC revisions in the pressure-face coefficient and `phiHbyA` identities. |
| Exactly two PISO pressure corrections | **KEEP count**, **REVISE schedule** | Define the thermophysical predictor and terminal EOS/continuity closure around the two corrections; do not equate this with AMReX's MAC+nodal projection pair. |
| Current-step final flux consumed by enthalpy/species after PISO | **REVISE** | It creates a coupling cycle with EOS. Choose and test an explicit second-order segregated contract; recommended v0.4 contract is described below. |
| Static IBM geometry/stencil/quadrature caches | **KEEP** | Rebuild only on geometry/topology/partition/BC-plan identity changes. |
| Calling HUNDUN ghost-cell IBM "AMReX EB" numerically | **REVISE terminology** | Borrow cache ownership/lifetime only; AMReX EB2's cut-cell metrics are a different numerical method. |
| Explicit halo begin/finish and persistent buffers/requests | **KEEP**, **REVISE coverage** | Compile the complete overlap graph including edge/corner/periodic images or prove staged face exchange fills them. `FillPatcher` is not needed for v0.4. |
| Four-layer linear lifecycle | **REVISE split** | Separate structural coarsening plan, exact numeric operator/coarse coefficients, reusable preconditioner setup, and solve workspace. |
| HYPRE Struct adapter | **KEEP**, **REVISE lifecycle** | Persist grid/stencil/matrix/vector/solver objects; refill/assemble values on coefficient change; make `Setup` reuse an explicit preconditioner policy. Keep IBM compact rows outside Struct. |
| One production freeze | **KEEP**, **REVISE implementation order** | Use registration -> analysis -> allocation/instantiation -> binding -> seal. |
| Minimal synchronous restart and VTK/Visit output | **KEEP**, minor **REVISE** | Add parent-directory durability, root header broadcast, and explicit snapshot lifetime. |
| Async output | **DEFER** | It requires either a copy or pinned state layer and often `MPI_THREAD_MULTIPLE`; neither belongs in the v0.4 hot path. |
| Advanced communication-avoiding/pipelined Krylov | **DEFER**, measure now | First record reduction counts and time; introduce only after the frozen baseline identifies reductions as material. |
| AMR FillPatch/reflux/average-down and AMReX data types | **DEFER/EXCLUDE** | They solve a problem v0.4 explicitly does not support. |

## 1. Pressure--velocity intermediates and final flux

### Public-source finding

OpenFOAM's incompressible pressure correction constructs reciprocal momentum diagonal, `HbyA`, predicted face flux, the consistent reciprocal diagonal when enabled, the pressure equation, the corrected pressure-equation flux, and the velocity correction within one pressure-correction call. The final face flux is taken from the pressure equation, and velocity is corrected from the corresponding pressure gradient. See the pinned [`incompressibleFluid/correctPressure.C`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/b9da51ab0673423aa2af6a45a72a3fbec9c66f9f/applications/modules/incompressibleFluid/correctPressure.C).

The compressible/isothermal path adds thermodynamic density and `psi`, constructs a density-weighted face coefficient, includes the full density time derivative plus the implicit pressure derivative in the pressure equation, updates thermodynamic density after the pressure solve, and corrects velocity through the same pressure path. See the pinned [`isothermalFluid/correctPressure.C`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/b9da51ab0673423aa2af6a45a72a3fbec9c66f9f/applications/modules/isothermalFluid/correctPressure.C).

OpenFOAM's pressure-corrector wrapper calls that pressure routine once per PISO correction, so corrector-local intermediates are recreated/re-evaluated; see [`isothermalFluid.C`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/b9da51ab0673423aa2af6a45a72a3fbec9c66f9f/applications/modules/isothermalFluid/isothermalFluid.C). This supports HUNDUN's requirement that corrector 2 cannot blindly reuse corrector-1 `HbyA/phiHbyA` after trial `U`, pressure, density, or flux changed.

### Decision

**KEEP** the distinct intermediate identities and the single final-flux writer. **REVISE** the dependency tuples:

- `rAU`: momentum symbolic/numeric diagonal, density used in the transient momentum term, implicit-source diagonal, constraints and numeric boundary coefficients.
- consistent reciprocal diagonal: every `rAU` dependency plus the consistent correction term.
- `HbyA`: complete momentum numeric identity, current trial `U`, current pressure-gradient subtraction convention, constraints and velocity boundary values.
- pressure face coefficient: reciprocal diagonal **and current face-density/thermo revision**; this is distinct from cell `rAU/rAtU`.
- `phiHbyA`: `HbyA`, pressure face coefficient, current trial/accepted flux-history revisions, BDF coefficients, mesh metrics, pressure/velocity BC numeric values, `p_ref`/gauge policy, and current corrector index.
- pressure storage/RHS: current `rho`, accepted density histories, `drho_dp_hY`, local absolute pressure, time coefficients, mass sources, and boundary closure.

Task 15 should mutate each newly listed dependency. A cache hit is legal only when the complete tuple matches; "same geometry" or "same momentum matrix" is insufficient.

### Clean-room boundary

OpenFOAM is GPL. HUNDUN may use the public algebraic relationship among momentum diagonal, predicted velocity, pressure equation, corrected flux and velocity, but must retain independently designed types, control flow, tests and implementation. Do not copy the cited expressions, comments, class organization, or line-by-line sequence. The fixed upstream [`COPYING`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/b9da51ab0673423aa2af6a45a72a3fbec9c66f9f/COPYING) is the license authority.

## 2. Low-Mach thermodynamics, `drho/dp`, and the two-corrector schedule

### Public-source finding

OpenFOAM's module driver orders each outer iteration as momentum predictor -> thermophysical/species/energy predictor -> pressure corrector; see pinned [`foamRun.C`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/b9da51ab0673423aa2af6a45a72a3fbec9c66f9f/applications/solvers/foamRun/foamRun.C). Its multicomponent thermophysical predictor advances species and energy and then corrects thermo before pressure correction; see [`multicomponentFluid/thermophysicalPredictor.C`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/b9da51ab0673423aa2af6a45a72a3fbec9c66f9f/applications/modules/multicomponentFluid/thermophysicalPredictor.C). The compressible pressure routine then holds the updated thermophysical variables fixed while correcting pressure-dependent density through `psi`.

AMReX-Hydro documents two projections for a typical second-order low-Mach/incompressible step, but they have different roles: a face-centered MAC projection makes half-time advective velocities satisfy the divergence constraint before flux construction, while a later cell-centered projection enforces the new-time constraint. They are not two adjacent PISO pressure corrections. See the official [Projection Methods](https://amrex-fluids.github.io/amrex-hydro/docs_html/Projections.html). IncFlo likewise uses projected MAC velocities for advection and a later nodal projection; its fixed source shows the MAC projector can retain operator/MLMG objects while updating coefficients in [`incflo_compute_MAC_projected_velocities.cpp`](https://github.com/AMReX-Fluids/incflo/blob/7307d8725c2a538f09cafbeacbfeb63e0fb11d22/src/convection/incflo_compute_MAC_projected_velocities.cpp), whereas its nodal path constructs a projector for the cell update in [`incflo_apply_nodal_projection.cpp`](https://github.com/AMReX-Fluids/incflo/blob/7307d8725c2a538f09cafbeacbfeb63e0fb11d22/src/projection/incflo_apply_nodal_projection.cpp).

### Plan defect

The current schedule is:

```text
thermo/transport -> momentum -> PISO 1 -> PISO 2/final flux
-> enthalpy -> species/passive -> commit
```

After PISO 2, enthalpy/species change `h,Y`; EOS therefore changes `rho(p_abs,h,Y)`. Unless the pressure equation already contains the exact discrete response of those later solves, the committed state can no longer simultaneously satisfy:

```text
rho = EOS(p_ref + pi, h, Y)
BDF(rho) + div(final_mass_flux) = mass_source
```

The plan currently asserts both but does not specify the missing coupled response. Two pressure corrections do not by themselves resolve this ordering cycle.

### Recommended v0.4 revision

**KEEP exactly two pressure solves and no PIMPLE outer loop**, but revise Tasks 1, 8, 14, 15 and the hot-loop schedule to define this segregated contract explicitly:

1. Form a second-order thermophysical predictor for `h,Y` before the pressure loop using only committed face-flux histories (BE on startup/restart; variable-step extrapolation thereafter) and declared sources/diffusion. This remains a consumer of the sole committed face-flux authority at known time revisions; it is not a second final-flux writer.
2. Evaluate `rho* = EOS(p_abs^k,h*,Y*)`, `drho_dp_hY`, transport and sound-speed/Mach from the predicted thermophysical state.
3. For a closed domain, perform the bounded `p_ref` mass Newton now, impose a defined mean/gauge on `pi`, and invalidate all pressure-storage/face-coefficient identities. For an open domain, keep `p_ref` under absolute-pressure boundary authority.
4. Hold `h*,Y*` fixed during the two PISO pressure corrections. On correction `k`, assemble the full BDF density defect plus the pressure increment storage `a0 * V * drho_dp_hY * delta_pi`, rather than a standalone `drho/dp` term that could double-count or omit `BDF(rho)`.
5. After each pressure solution, update/evaluate EOS with the new local absolute pressure. Corrector 2 absorbs the remaining nonlinear density/continuity defect and is the only current-step final `U`/mass-flux publisher.
6. Before commit, require four separately reported terminal defects: EOS residual, discrete continuity residual using the published face flux, closed-domain mass residual (when applicable), and pressure-gauge/boundary closure residual. Reject the transaction if any is nonfinite, violates `p_abs/T/Y` bounds, or exceeds its pre-registered tolerance.

This recommendation changes the plan sentence that current-step enthalpy/species "require the certified final flux". If the project insists that those equations consume the newly published current-step final flux, then an additional coupled response/outer iteration is mathematically required; that would contradict the current no-PIMPLE/no-extra-corrector constraint and must be designed explicitly rather than hidden in cache logic.

Add a coupled manufactured test with nonzero, independently varying pressure, temperature and composition so that mutations of thermophysical ordering, `drho/dp`, density history, or `p_ref` fail. The existing pressure-storage unit test is necessary but not sufficient.

## 3. Storage, arena and workspace lifecycle

### Public-source finding

AMReX `MultiFab` stores box-local arrays with explicit valid and ghost regions, and each `FArrayBox` owns contiguous component data. It supports an optional single-chunk allocation and selectable arenas; see official [Basics: MultiFab and Memory Allocation](https://amrex-codes.github.io/amrex/docs_html/Basics.html). AMReX's arenas are pools intended to reduce allocation overhead; see [GPU/Memory Allocation](https://amrex-codes.github.io/amrex/docs_html/GPU.html), whose arena concepts also document retained pools even though HUNDUN v0.4 is CPU-only.

### Decision

**KEEP** padded x-fast SoA, precomputed offsets, one owning allocation per NUMA owner, handle rotation, non-owning hot views and maximum-capacity solver scratch. This is an independent CPU specialization and is stricter than AMReX's general-purpose ownership.

**REVISE** Tasks 4, 6 and 18 to require:

- parallel first-touch by the exact worker/NUMA placement that will own each page;
- 64-byte alignment for field bases and cache-line separation for per-thread counters/scratch;
- capacity and overflow checks before allocation, including all Krylov bases, MG levels, halo buffers, IBM compact rows and service snapshots;
- a stable-address test after retry and after every legal numeric refill;
- separate allocation capacity from numeric validity: freeze may allocate coefficient/hierarchy storage, but it must begin invalid and receive a revision certificate only after initial state/BC coefficients are filled.

Explicit huge pages, alternate AoSoA layouts and broader SIMD specialization are **DEFER** items to be selected from full-grid profiles, not assumptions to add now.

## 4. Halo, physical boundary fill and FillPatch

### Public-source finding

AMReX distinguishes interior/periodic ghost exchange from physical-boundary fill. `MultiFab::FillBoundary` fills overlaps with valid data and optional periodic images; application boundary logic fills physical boundaries. See official [Basics: Boundary Conditions](https://amrex-codes.github.io/amrex/docs_html/Basics.html). Its communication metadata are cached by layout/distribution identity in pinned [`AMReX_FabArrayBase.H`](https://github.com/AMReX-Codes/amrex/blob/59d066aab774bc388cc6ed944f7beaf645607ed3/Src/Base/AMReX_FabArrayBase.H), and it exposes split begin/finish exchange in [`AMReX_FabArrayCommI.H`](https://github.com/AMReX-Codes/amrex/blob/59d066aab774bc388cc6ed944f7beaf645607ed3/Src/Base/AMReX_FabArrayCommI.H).

That AMReX source creates per-exchange state and send/receive request containers, so HUNDUN's fixed buffers and persistent MPI requests are not a source copy and are a stronger, project-specific hot-path optimization.

AMReX FillPatch exists primarily to combine same-level exchange, time interpolation, physical boundary fill and coarse/fine interpolation; see official [AmrCore FillPatch details](https://amrex-codes.github.io/amrex/docs_html/AmrCore.html) and [`FillPatcher`](https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1FillPatcher.html). HUNDUN v0.4 has one static Cartesian level, so adopting FillPatch, reflux or average-down would add the wrong lifecycle.

### Decision

**KEEP** explicit `begin/finish`, compute/communication overlap, revision publication only after `finish`, persistent addresses, and a separate compiled physical-BC stage.

**REVISE** "direct-neighbor-only" in Task 6. Define it as an exact stencil-overlap graph, not merely six face ranks. A width-four quadratic IBM donor or a multidimensional reconstruction can require edge/corner ghost values. The communication plan must do one of the following and test it:

- communicate directly with every rank whose owned box overlaps the grown local box, including edge/corner/periodic images; or
- use a fixed staged x/y/z exchange and prove later stages forward already-received data to fill all required corners.

Also reject a process grid whose minimum local extent cannot satisfy the registered staged reach, or compile a multi-hop schedule deliberately. Add periodic self-neighbor, same-peer-on-two-sides, unequal patch, edge/corner donor, ghost width 1--4, and simultaneous-field packing tests. Compile request/tag sets so only one instance of a persistent request is in flight; overlapping stages need distinct request objects or a proven serialized schedule.

## 5. Static IBM geometry versus AMReX EB2

### Public-source finding

AMReX EB2 builds a geometric database from an implicit function or STL and exposes long-lived geometry through an index space/level and an `EBFArrayBoxFactory`; see official [Embedded Boundary documentation](https://amrex-codes.github.io/amrex/docs_html/EB.html), pinned [`AMReX_EB2_Level.H`](https://github.com/AMReX-Codes/amrex/blob/59d066aab774bc388cc6ed944f7beaf645607ed3/Src/EB/AMReX_EB2_Level.H), and [`AMReX_EBFabFactory.H`](https://github.com/AMReX-Codes/amrex/blob/59d066aab774bc388cc6ed944f7beaf645607ed3/Src/EB/AMReX_EBFabFactory.H). The factory provides cut-cell quantities such as volume fraction, area fraction and centroids.

### Decision

**KEEP** the lifetime idea: STL parse/scan, region topology, donor stencils, surface quadrature, compact interface rows and wall geometry remain immutable for static geometry and rebuild only from explicit geometry/topology/partition/BC-plan changes.

**REVISE terminology and validation claims.** HUNDUN's planned positive-normal quadratic ghost-cell IBM is not AMReX's cut-cell finite-volume discretization. Do not claim AMReX EB conservation or stencil behavior. Borrow only the public ownership/lifetime pattern. Keep HUNDUN's own second-order reconstruction, penetration, traction, reaction and global conservation oracles.

Task 18's final halo reach must be derived after all IBM donor and compact-row registrations. A structured width-four halo is sufficient only if every stored donor lies in the certified grown box; otherwise an explicit irregular gather plan must register before freeze.

## 6. Operator, multigrid hierarchy and coarse numeric refresh

### Public-source finding

AMReX's `MLMG` owns a solve process around an `MLLinOp`; the operator type is chosen from the discretization, and coefficients/BCs are supplied before solve. EB uses a distinct `MLEBABecLap`. See official [MLMG and Linear Operator Classes](https://amrex-codes.github.io/amrex/docs_html/LinearSolvers.html).

AMReX-Hydro's pinned `MacProjector` owns its operator, `MLMG`, RHS, solution and flux arrays; it provides `updateBeta/updateCoeffs` to replace numeric coefficients without reconstructing the projector. See [`hydro_MacProjector.H`](https://github.com/AMReX-Fluids/amrex-hydro/blob/e49df248aabd2cc11865eb5be734a2f5f2f65ee5/Projections/hydro_MacProjector.H) and [`hydro_MacProjector.cpp`](https://github.com/AMReX-Fluids/amrex-hydro/blob/e49df248aabd2cc11865eb5be734a2f5f2f65ee5/Projections/hydro_MacProjector.cpp). IncFlo's MAC path calls those update methods when its density coefficient changes, confirming that long-lived object identity does not mean stale coefficients are silently reused.

HYPRE Struct separately creates/assembles grid, stencil, matrix and vectors, then performs solver setup and solve; see official [Structured-Grid System Interface](https://hypre.readthedocs.io/en/latest/ch-struct.html) and [solver setup/solve overview](https://hypre.readthedocs.io/en/latest/ch-solvers.html). PETSc likewise exposes explicit preconditioner reuse even when the operator changes, and documents that this suppresses the normal update; see [`PCSetReusePreconditioner`](https://petsc.org/release/manualpages/PC/PCSetReusePreconditioner/). These are lifecycle controls, not permission to use an uncertified exact operator.

### Decision

**REVISE** Tasks 11, 12 and 18 to use these four non-overlapping states:

1. `SymbolicPlan/CoarseningPlan`: topology, ownership, stencil pattern, boundary position/type pattern, level shapes, transfer sparsity and line/color schedule. It changes only for structural identities.
2. `ExactNumericState`: fine operator and every coarse numeric coefficient required by the selected Galerkin/geometric policy, plus numeric BC coefficients. Every coefficient revision refills/restricts what the exact `apply` and residual require before certification. No ratio threshold may skip this.
3. `PreconditionerSetupState`: smoother factors, coarse-solver setup and optional HYPRE setup. A pre-registered coefficient-change policy may reuse this across solves. Reuse is safe only because the outer Krylov method applies the current exact operator and checks a current FP64 true residual. Record reuse, forced refresh and convergence degradation.
4. `SolverWorkspace`: Krylov basis, MG level vectors, reduction storage and backend work arrays at fixed addresses; numeric changes never replace it unless capacity/backend changes.

The existing `coefficient_change_rebuild_ratio` should control item 3, not whether item 2 is current. Rename the counters accordingly (`exact_numeric_refill`, `coarse_numeric_refresh`, `preconditioner_setup`, `symbolic_rebuild`). A threshold crossing may rebuild preconditioner setup once; it must not imply rebuilding immutable coarsening topology unless a structural identity changed.

For HYPRE Struct, persist grid/stencil/matrix/vector/solver handles; update matrix values box-wise and assemble on coefficient change. Decide `Setup` reuse by the preconditioner policy and validate it by exact outer residual. Struct is scalar, cell-centered and fixed-stencil; it can represent the regular Cartesian pressure preconditioner but not HUNDUN's irregular compact IBM rows. Task 17's exact FGMRES operator must remain HUNDUN-owned, with HYPRE/Native MG only as preconditioner for that system.

## 7. Production compile/freeze order

### Plan defect

The Task 18 sequence freezes `FieldSchema -> ArenaLayout -> ... -> CommunicationPlan -> FrozenExecutionGraph`, while Task 10 says graph compilation determines workspace liveness, scratch alias intervals and halo nodes. Arena size and communication capacity therefore depend on graph analysis, but an executable graph cannot be finalized until arena views and communication handles exist.

Numeric state and hierarchy setup also depend on initialized field/BC coefficients. They can be allocated before the first attempt, but cannot truthfully carry a valid coefficient identity at case-plan freeze.

### Decision

**KEEP** one public `freeze_product` operation and post-freeze immutability, but **REVISE** it internally to a two-phase compiler:

```text
all capability registrations
-> canonical logical IR and dependency validation
-> graph/resource analysis (liveness, max workspace, ghost overlap,
   collective budget, service snapshot capacity)
-> freeze FieldSchema and compute ArenaLayout
-> allocate/first-touch arenas
-> instantiate geometry/BC/IBM/symbolic/coarsening plans
-> instantiate empty numeric/preconditioner/workspace capacities
-> instantiate merged communication metadata/buffers/requests
-> bind field views, operator views, halo tickets and graph nodes
-> validate fingerprints and seal
```

At driver initialization, after initial/boundary fields exist:

```text
fill exact numeric coefficients -> refresh coarse numeric coefficients
-> perform required preconditioner setup -> publish numeric certificates
-> enter first attempt
```

Add mutation tests for every phase boundary: a graph-derived workspace or ghost requirement introduced after analysis must reject; binding cannot change analyzed capacity; a never-filled numeric state cannot solve; replay must produce identical fingerprints and offsets on 1/2/4 ranks.

## 8. Collective count and failure consensus

### Public-source finding

AMReX profiling reports minimum/average/maximum time across processes and notes that communication timers can be dominated by wait time caused by load imbalance. Its optional synchronization diagnostic is explicitly not recommended for production because it perturbs execution. See official [AMReX Profiling Tools](https://amrex-codes.github.io/amrex/docs_html/AMReX_Profiling_Tools.html).

### Decision

**REVISE** Tasks 10, 11, 18, 20 and 21 so `ResourceCounters` and evidence include, per step and per solver:

- point-to-point message count and payload bytes;
- blocking and nonblocking collective count, payload and operation kind;
- Krylov global reductions per iteration and total;
- MG/coarse-solver collectives;
- closed-mass Newton and timestep-control reductions;
- failure-consensus reductions not already piggybacked on a required collective;
- time in halo pack/start, interior compute, wait/finish, unpack, pressure operator, preconditioner, reductions and I/O.

Compile a maximum collective budget in the frozen graph. Coalesce compatible scalar diagnostics/status with an already required stage reduction when ordering and error semantics permit. Do not add an unconditional barrier per stage or per timer. Use one synchronization at the external benchmark-region boundary; inside the run, report per-rank timing and reduce/report min/mean/max after the measured region.

Advanced pipelined/communication-avoiding PCG/FGMRES is **DEFER** until full-grid profiles show reductions are material. The v0.4 plan should nevertheless prevent an accidental extra all-reduce for every residual print, cache certificate or stage success flag.

MPI point-to-point failure recovery must also avoid one rank returning before its peers enter a matching communication. The graph should distinguish local preflight failures (consensus before starting the next collective) from failures after communication has started (complete/cancel according to the registered protocol, then consensus). This is a correctness condition and a performance reason to keep consensus points few and explicit.

## 9. Restart and output lifecycle

### Public-source finding

AMReX intentionally leaves checkpoint content to each application. Its documented pattern uses a problem-specific header, parallel field I/O, rank-zero header writing, and broadcast of header data on restart; see official [I/O: Checkpoint File](https://amrex-codes.github.io/amrex/docs_html/IO.html). This supports HUNDUN's typed, application-owned minimal restart rather than adopting an upstream file format.

AMReX async output copies the state and may require `MPI_THREAD_MULTIPLE` when fewer output files than ranks are used; it can also oversubscribe compute threads. The same official [I/O documentation](https://amrex-codes.github.io/amrex/docs_html/IO.html) describes those constraints.

POSIX `rename` is atomic with respect to directory entries, but durability requires synchronization. Linux documents that file `fsync` alone does not make the containing directory entry durable; the directory must also be synchronized. See [`rename`](https://pubs.opengroup.org/onlinepubs/9799919799/functions/rename.html), [`fsync`](https://man7.org/linux/man-pages/man2/fsync.2.html), and the POSIX [directory durability rationale](https://pubs.opengroup.org/onlinepubs/9799919799/xrat/V4_xbd_chap01.html).

### Decision

**KEEP** committed-state-only services, bounded `keep_last=1`, pending-generation validation, rank-change restart, VTI/VTR+Visit, and BE on the first resumed step.

**REVISE** Task 19 to require this crash-consistent default sequence on one filesystem:

```text
create pending generation -> write all rank data and header
-> close and fsync files -> fsync generation directories
-> collective validate hashes/schema
-> atomically rename/switch current -> fsync parent directory
-> remove previous generation -> fsync parent directory again
```

Rank zero should read/broadcast canonical metadata before ranks allocate/publish restarted state. No state is published until every rank validates its assigned data and collective status succeeds. Make snapshot lifetime explicit: synchronous output completes before the corresponding accepted layer may become writable. Async output is **DEFER** because it needs a pinned layer or a full snapshot copy and a separate MPI/thread policy.

## 10. Fair performance and release evidence

### Plan strengths to keep

**KEEP** full `480x480x48/64-rank` testing, fixed candidate identity, no `24^3` performance claim, two-step directional screening only, candidate freeze before long tests, alternating HUNDUN/COAST launch order, and separation of COAST short performance from HUNDUN-only periodic experimental validation.

AMReX's profiler's min/average/max-rank reporting reinforces that a parallel step is governed by the slowest rank, not rank zero or the mean; see [AMReX Profiling Tools](https://amrex-codes.github.io/amrex/docs_html/AMReX_Profiling_Tools.html).

### Required revision

**REVISE** Tasks 20--21:

1. Define one external timed hot region with output/checkpoint/statistics serialization disabled for both products. Synchronize once immediately before it; wall time is the maximum rank elapsed time or launcher-observed job elapsed time. Internal stage timers report min/mean/max without barriers inside stages.
2. Pre-register warmup and measured steps. The 2-step run remains directional. For each formal 20-step run, state exactly which steps are warmup and which are measured; do not let first coefficient fill/hierarchy setup contaminate the hot median unless both products' metric explicitly includes it.
3. Compare alternating runs as **paired ratios** `HUNDUN_hot / COAST_hot` from the same resource/time block. Report median paired ratio, P90 step time and a pre-registered uncertainty interval. Global medians from two unpaired samples can hide drift. Use at least five pairs, with a pre-registered maximum extension if the interval is inconclusive.
4. Record maximum per-rank RSS and node aggregate, not only rank-zero RSS; record frequency governor/turbo state when readable, SMT, NUMA binding, compiler/linker and optimization flags, MPI library, filesystem/output state, and pre/post process inventory.
5. Record exact operator refills, coarse numeric refreshes, preconditioner setups, hierarchy structural rebuilds, Krylov/MG iterations, collectives, halo messages/bytes and rejected attempts. A faster run obtained by using stale exact coefficients or looser true-residual/scientific tolerances is invalid.
6. Keep an external total-step timing as the performance authority. Fine-grained counters diagnose why; they must not replace the paired end-to-end result.

The release predicate may remain "no slower than COAST", but it should be evaluated on the paired metric with its pre-registered uncertainty rule rather than two independent medians.

## 11. Concrete edits by plan task

| Task | Decision | Required edit |
|---|---|---|
| Global constraints / Task 1 | **REVISE** | Replace the current hot schedule with the explicit thermophysical-predictor/two-pressure-corrector contract; expand pressure intermediate identities. |
| Task 4 | **REVISE** | Add NUMA first-touch, cache-line-isolated thread scratch, full capacity accounting and invalid-until-filled numeric storage. |
| Task 6 | **REVISE** | Define exact grown-box overlap or staged corner-fill semantics; register persistent request single-flight/tag rules; reject insufficient local extents. |
| Task 8 | **REVISE** | Place closed-domain `p_ref` Newton after thermophysical prediction and before pressure correction; define `pi` gauge and invalidations. |
| Task 10 | **REVISE** | Split logical graph/resource analysis from executable binding; add collective budgets. |
| Task 11 | **REVISE** | Count reductions; distinguish exact numeric state, preconditioner setup and workspace; keep current true residual. |
| Task 12 | **REVISE** | Always refresh exact/coarse numeric coefficients required by current operator; use threshold only for preconditioner setup; persist HYPRE handles. |
| Task 13 | **KEEP/REVISE wording** | Keep static caches; explicitly distinguish ghost-cell IBM from cut-cell EB2 and certify donor halo reach. |
| Task 14 | **REVISE** | Add a second-order thermophysical predictor contract using committed flux histories; add a coupled `p,h,Y,rho` manufactured oracle. |
| Task 15 | **REVISE before implementation** | Define full BDF density defect + `drho/dp` pressure increment, two nonlinear corrections, `p_abs` bounds and terminal EOS/continuity/mass/gauge gates. |
| Task 18 | **REVISE** | Implement registration -> analysis -> allocation/instantiation -> binding -> seal; fill numeric coefficients during driver initialization. |
| Task 19 | **REVISE** | Add directory fsync and explicit synchronous snapshot lifetime; root metadata broadcast. |
| Tasks 20--21 | **REVISE** | Add max-rank timing, exact warmup, paired ratios/uncertainty, collective counts and fair-build/resource receipts. |

## 12. Items explicitly not recommended for v0.4

- Do not add AMR, FillPatchTwoLevels, reflux, average-down or nonmatching patches.
- Do not replace the chosen ghost-cell IBM with AMReX cut-cell EB during this version.
- Do not add SIMPLE, PIMPLE, a configurable pressure-corrector loop or a third pressure solve to conceal an unspecified coupling defect.
- Do not copy OpenFOAM GPL code or translate its control flow. Treat it solely as a mathematical/lifecycle reference.
- Do not adopt AMReX classes or data layout wholesale. Its general-purpose per-call communication/projector allocations are not HUNDUN's CPU performance target.
- Do not enable asynchronous output until snapshot ownership, MPI thread level and measured benefit are designed.
- Do not add pipelined Krylov variants before reduction counts and full-grid profiles justify the numerical-complexity tradeoff.

## Source identity and license boundary

The plan's fixed source identities remain appropriate:

- OpenFOAM-dev `b9da51ab0673423aa2af6a45a72a3fbec9c66f9f`, GPL-3.0-or-later: reference public mathematics and lifecycle only.
- AMReX `59d066aab774bc388cc6ed944f7beaf645607ed3`, BSD-3-Clause: independently implement selected data/lifetime ideas; no AMR machinery in v0.4.
- IncFlo `7307d8725c2a538f09cafbeacbfeb63e0fb11d22`, BSD-3-Clause: projection sequencing and coefficient-refresh evidence only; IncFlo solves variable-density incompressible equations, not HUNDUN's local-pressure EOS system.
- AMReX-Hydro `e49df248aabd2cc11865eb5be734a2f5f2f65ee5`, BSD-3-Clause: long-lived projector/operator/MLMG and explicit coefficient refresh as lifecycle references.
- HYPRE official documentation, dual Apache-2.0/MIT in current releases: optional backend interface/lifecycle reference.
- PETSc official documentation, BSD-2-Clause: optional corroboration for explicit operator/preconditioner reuse controls; no PETSc dependency is proposed.

No upstream implementation text has been added to HUNDUN by this research note. All recommendations concern public equations, ownership boundaries, invalidation rules, measurement protocol and independently implemented tests.

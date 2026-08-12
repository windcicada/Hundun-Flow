<!-- SPDX-License-Identifier: Apache-2.0 -->

# HUNDUN-FLOW v0.4 public lifecycle survey

Date: 2026-08-13

Status: frozen research record for v0.4

This record freezes the public mathematical, layout, and resource-lifecycle ideas used to
design HUNDUN-FLOW v0.4. It is not a source-import record: Task 1 vendors, links, translates,
or makes a runtime dependency of none of the referenced sources. The sole later reuse
authorization is stated explicitly in sections 5 and 6. The machine-reviewable decision record
is
[`2026-08-13-hundun-v04-adoption-ledger.tsv`](2026-08-13-hundun-v04-adoption-ledger.tsv),
and the resulting execution order is
[`v0.4-target-hot-loop.md`](../architecture/v0.4-target-hot-loop.md).

## 1. Source identities and legal boundary

The four public repositories were inspected at the exact commits below. COAST was inspected
in place and read-only at the fixed path. The COAST working tree was already modified by its
owner; this survey did not change it and therefore treats the path, not its current Git tree,
as the fixed identity.

| Project | Fixed repository/path and revision | Files and symbols inspected | License boundary | Public idea or observed function | HUNDUN destination and decision |
| --- | --- | --- | --- | --- | --- |
| OpenFOAM-dev | [`OpenFOAM/OpenFOAM-dev`](https://github.com/OpenFOAM/OpenFOAM-dev/tree/b9da51ab0673423aa2af6a45a72a3fbec9c66f9f), `b9da51ab0673423aa2af6a45a72a3fbec9c66f9f` | `applications/modules/incompressibleFluid/correctPressure.C`: `rAU`, `rAtU`, `HbyA`, `phiHbyA`, pressure-equation flux, final velocity update; `COPYING` | GPL-3.0-or-later; public mathematics and lifecycle are reference-only. No source, expression, control flow, comments, naming arrangement, or translation may enter HUNDUN. | PISO intermediate dependencies and the pressure-solve relationship among pressure flux, face flux, pressure gradient, and velocity | `solver_piso`; clean-room equations, types, schedule, and tests |
| AMReX | [`AMReX-Codes/amrex`](https://github.com/AMReX-Codes/amrex/tree/59d066aab774bc388cc6ed944f7beaf645607ed3), `59d066aab774bc388cc6ed944f7beaf645607ed3` | `Src/Base/AMReX_FabArrayBase.H`: `BDKey`, `FB`, `CPC`, `getFB`, `getCPC`; `Src/Base/AMReX_FabArrayCommI.H`: `FillBoundary_nowait`, `FillBoundary_finish`; `Src/AmrCore/AMReX_FillPatcher.H`: `FillPatcher`; `Src/EB/AMReX_EB2_Level.H`: `EB2::Level`; `Src/EB/AMReX_EBFabFactory.H`: `EBFArrayBoxFactory`; `LICENSE` | BSD-3-Clause-style license text. HUNDUN nevertheless independently reimplements the ideas and exposes no AMReX type. | Box-local storage, cached halo metadata with begin/finish, fill-patch metadata, and EB resource ownership | `core_field_storage`, `parallel_communication`, `mesh_eb`; AMR fill-patch, reflux, and average-down are excluded from v0.4 |
| IncFlo | [`AMReX-Fluids/incflo`](https://github.com/AMReX-Fluids/incflo/tree/7307d8725c2a538f09cafbeacbfeb63e0fb11d22), `7307d8725c2a538f09cafbeacbfeb63e0fb11d22` | `src/incflo_regrid.cpp`: `MakeNewLevelFromCoarse`, `RemakeLevel`, `ClearLevel`; `src/convection/incflo_compute_MAC_projected_velocities.cpp`: `compute_MAC_projected_velocities`; `src/projection/incflo_apply_nodal_projection.cpp`: `ApplyNodalProjection`; `src/embedded_boundaries/eb_stl.cpp`; `LICENSE` | BSD-3-Clause. HUNDUN uses an independent implementation, not IncFlo classes or control flow. | Projection and EB resources are reconstructed at grid/ownership boundaries while current values and coefficients have narrower refresh boundaries. | `solver_projection`, `mesh_eb`; regrid and all AMR mechanisms are excluded |
| AMReX-Hydro | [`AMReX-Fluids/amrex-hydro`](https://github.com/AMReX-Fluids/amrex-hydro/tree/e49df248aabd2cc11865eb5be734a2f5f2f65ee5), `e49df248aabd2cc11865eb5be734a2f5f2f65ee5` | `Projections/hydro_MacProjector.{H,cpp}`: `MacProjector`, `updateBeta`, `updateCoeffs`, `project`, `m_linop`, `m_mlmg`; `Projections/hydro_NodalProjector.H`: `NodalProjector`, `define`, `m_linop`, `m_mlmg`; `LICENSE` | BSD-3-Clause-style license text. HUNDUN independently implements project-owned operator, hierarchy, and workspace abstractions. | Long-lived projector/operator/solver objects, explicit coefficient refresh, retained solve fields, and multigrid reuse | `solver_linear`; no AMReX-Hydro object or API is a HUNDUN authority |
| COAST | `/home/wyf/code_dev/Coast_software` (fixed read-only location; observed 2026-08-13, Git `HEAD` `46afa4137d3d4bab96db40846b0580567b82e557` plus pre-existing tracked changes) | Functionality paths are listed in section 5. | User-specified read-only functional reference. No license grant or redistribution right is presumed. Old COAST Fortran must not be copied or translated. The sole source-reuse exception is the user's authorization for Task 5 to port/adapt the `imb_mesh_y.cpp` scanning mathematics/method without provenance boilerplate; it does not authorize other COAST C++. Task 1 copies no source. | Replacement capability inventory, SIMPLE/ICCG hot-cycle observations, array/precomputation concepts, and I/O surface | HUNDUN-native `app_`, `mesh_`, `bc_`, `physics_`, `solver_`, and `io_` implementations; SIMPLE and COAST boundary limits are not v0.4 requirements |

The public-repository license files were read at their fixed commits, not inferred from the
current default branches. The COAST tree was inspected only for path/symbol/function evidence.
This document deliberately contains no source excerpt or line-by-line paraphrase.

## 2. PISO intermediate and authority contract

The v0.4 dependency chain is exactly:

```text
momentum numeric revision -> rAU
rAU + consistent diagonal revision -> rAtU
momentum numeric revision + current trial U -> HbyA
current HbyA + current trial U/phi + time/geometry/BC -> phiHbyA
pressure equation flux -> only final face-mass-flux writer
pressure gradient from the same solve -> final U update
```

The chain is a HUNDUN contract, not copied implementation structure. Its consequences are:

1. `rAU` is certified by the momentum diagonal, constraint, and boundary-coefficient
   revisions. `rAtU` additionally carries the consistent diagonal/time-correction revision.
2. `HbyA` is not merely a function of `rAU`: it is certified against the complete momentum
   numeric state and the current trial `U` plus applicable constraints.
3. `phiHbyA` is certified against the current `HbyA`, trial `U/phi`, time discretization,
   geometry, and boundary revisions.
4. Corrector 1 changes trial `U/phi`. Corrector 2 must therefore rebuild `HbyA/phiHbyA` or
   prove, using their full dependency tuples, that each remains valid. A shared label such as
   “momentum cache valid” is insufficient.
5. `rAU/rAtU` alone may be reused across the two correctors when every coefficient revision in
   their respective tuples is unchanged. A new RHS or changed trial velocity does not by
   itself invalidate these two coefficient-derived objects.
6. Corrector 1 may write a trial flux. Only the final corrector's pressure-equation flux path
   may publish `final face mass flux`; recomputing it from final cell velocity would create a
   second authority and is forbidden.
7. Final `U` uses the pressure gradient from that same final pressure solve. Energy, species,
   and passive scalars consume the published final face mass flux revision.

## 3. Geometry and boundary-plan lifetimes

All identities name immutable inputs. A revision change outside the listed identity cannot be
used as a reason to rebuild the resource.

| Resource and owner | Identity tuple | Lifetime | Only legal rebuild causes |
| --- | --- | --- | --- |
| `EBTopology`, `mesh_eb` | `(geometry_content_revision,geometry_transform_revision,cartesian_mesh_coordinates_revision,partition_revision,classification_policy_revision,compact_index_layout_revision)` | initialization through case shutdown for static geometry | geometry content, geometry transform, Cartesian mesh/coordinates, partition, classification policy, or compact-index layout revision changes; requested ghost extent, field values, time steps, and solver coefficients do not rebuild EBTopology |
| `BoundaryStencilPlan`, `mesh_eb`/`bc_plan` | `(EBTopology identity, BoundaryPlan revision, reconstruction-plan revision, donor-policy revision, required ghost extent, field-layout revision)` | Initialization through case shutdown | Any member of that tuple changes. A current field value, residual, retry, or time-step change cannot trigger donor search or stencil rebuild. |
| `SurfaceQuadraturePlan`, `mesh_eb` | `(EBTopology identity, surface-set revision, quadrature-family/order revision, compact-interface-index revision)` | Initialization through case shutdown | EB topology, selected surface, quadrature family/order, or compact interface mapping changes. Thermodynamic and flow-state changes do not rebuild it. |

For all three resources, borrowed views expire when their owning plan is replaced. No STL/BVH
query, donor search, quadratic-weight generation, or surface-quadrature construction is legal
inside the production hot loop. Requested ghost extent belongs to the factory/stencil or
communication capacity that consumes the topology; changing it does not change `EBTopology`.

## 4. Halo and linear-system lifetimes

| Resource and owner | Identity tuple | Lifetime | Only legal rebuild/refill/replacement causes |
| --- | --- | --- | --- |
| Persistent halo metadata (`HALO_METADATA`), `parallel_communication` | `(partition_revision,field_schema_revision,stage_ghost_set_revision,periodicity,peer_map,pack_span_layout)` | case execution after all stage ghost sets freeze | partition, field schema, registered stage ghost set, periodicity, peer map, or pack-span layout changes; buffer/request replacement and field numeric revisions do not rebuild metadata |
| Persistent halo buffers (`HALO_BUFFERS`), `parallel_communication` | `(send_capacity,receive_capacity,memory_kind,numa_placement,alignment)` | case execution at registered maximum capacity | replace only when required send/receive capacity exceeds allocation or memory kind, NUMA placement, or alignment plan changes; metadata/request replacement and new field values do not replace buffers |
| Persistent MPI requests (`HALO_REQUESTS`), `parallel_communication` | `(communicator_generation,peer_ranks,message_tags,message_counts,mpi_datatypes,registered_buffer_bindings)` | case execution; persistent requests survive time steps; one registered exchange instance is single-in-flight | recreate only when communicator, peer ranks, tags, counts, MPI datatypes, or registered buffer bindings change; metadata/buffer replacement alone does not recreate requests unless one of those request dependencies changes; field numeric revisions only invalidate ghost revisions, which publish after finish |
| `SymbolicPlan`, `solver_linear` | `(operator kind, scalar/face location, mesh topology, partition, EB interface pattern, boundary position/type pattern, stencil pattern, backend)` | Across all assemblies and solves with the same structural identity | Rebuild only when a structural identity member changes. Coefficients, RHS, tolerances, iteration count, time step, retry, and residual do not rebuild it. |
| `NumericState`, `solver_linear` | `(SymbolicPlan identity, diagonal revision, off-diagonal revision, time-coefficient revision, material/transport revision, numeric boundary-coefficient revision, constraint revision)` | Across solves while its full coefficient tuple is unchanged | Refill only when a coefficient tuple member changes. RHS or initial-guess changes do not refill. A structural change first replaces `SymbolicPlan`, which necessarily creates a new `NumericState`. |
| `HierarchyState`, `solver_linear` | `(SymbolicPlan identity, coarsening-plan revision, transfer/smoother revision, coefficient-policy epoch, backend)` | Across solves and time steps under the registered coefficient-change policy | Structural identity, coarsening, transfer/smoother, backend, or policy change rebuilds it. A coefficient change rebuilds it only when the pre-registered policy says the change crosses its threshold; otherwise coefficients are refreshed without changing hierarchy identity. RHS, tolerance, and iteration history never rebuild it. |
| `SolverWorkspace`, `solver_linear` | `(solver algorithm, backend, precision policy, execution-plan revision, maximum registered shape, maximum Krylov/subspace capacity, reduction layout)` | Allocated once at registered maximum capacity and reused by all compatible solves | Replace only for algorithm/backend/precision/execution-plan change or a required shape/subspace/reduction capacity beyond the registered maximum. A different operator value, RHS, tolerance, iteration count, or retry does not replace it. |

The four linear layers must remain separately measurable. “Rebuild solver” is not a valid
event name: counters must distinguish symbolic rebuild, numeric refill, hierarchy rebuild, and
workspace replacement. Halo dependencies are likewise non-transitive: buffer replacement never
rebuilds metadata, and it rebuilds a persistent request only if the request's registered buffer
binding or another request identity member actually changes.

## 5. COAST replacement capability trace

The fixed COAST reference is evidence of functions that HUNDUN must replace, not an algorithm
specification. The exact observations are:

| Required capability | Read-only evidence path and symbol | Observation and v0.4 disposition |
| --- | --- | --- |
| Local absolute-pressure EOS | `SRC.Coast/densty.F90::densty`, `SRC.Coast/temperature.F90::Temperature` | The compressible path evaluates thermodynamic state from local pressure, enthalpy, and composition. HUNDUN independently implements its approved `p_abs = p_ref + pi` EOS contract. |
| Pressure density derivative | `SRC.Coast/densty.F90::densty`, `SRC.Coast/compress.F90::compress` | A local pressure derivative is maintained and enters pressure-system assembly. HUNDUN owns and tests `(partial rho / partial p)_(h,Y)` independently. |
| Sound speed and Mach | `SRC.Coast/speedofsound.F90::speedofsound`, `SRC.Coast/averaged.F90::averaged`, `SRC.Coast/app/coast_legacy_driver.F90` Mach-field updates | COAST exposes mixture sound speed and Mach diagnostics. Two observed diagnostic forms have different meanings: boundary averaging uses `|U|/c`, while the main field path uses `|U| sqrt(drhodp)`. HUNDUN must select and document one thermodynamically consistent authority rather than silently equating them. |
| NSCBC | `SRC.Coast/bndry_NSCBC.F90::bndry_NSCBC`, `SRC.Coast/boundary_NSCBC.F90::boundary_NSCBC` | Subsonic characteristic-boundary capability is part of the replacement surface. HUNDUN boundary equations and validation are clean-room work. |
| SIMPLE-style transient iterations | `SRC.Coast/app/coast_legacy_driver.F90::coast_multibrick_run_timesteps`, `coast_multibrick_state_update`, `coast_multibrick_pressure_correction`; `SRC.Coast/input.F90::input` reads `niter` | The observed time-step loop alternates state update and pressure correction for `niter`. It is a functional/hot-cycle reference only. v0.4 requires transient PISO with exactly two pressure correctors and does not adopt SIMPLE. |
| Vreman wall function | `SRC.Coast/gamma_vreman.F90::gamma_vreman`, `SRC.Coast/wall.F90::WALL`, `SRC.Coast/boundary_wall.F90::boundary_wall` | The task shorthand resolves to two distinct observed capabilities: a Vreman SGS path and a wall-law path. HUNDUN must not conflate their revisions or ownership. |
| Flat `.d` inputs | `SRC.Coast/input.F90::input`; `EXEC/input.d`, `EXEC/probe.d`, `EXEC/vtk_output.d` | COAST consumes flat case files. HUNDUN uses `case.json` as sole configuration authority and flat referenced `.d` files only for bulk arrays/tables; all parsing remains cold-path. |
| STL scan | `SRC.Coast/imb_mesh_y.cpp` exported scan entry points and `SRC.Coast/imb_mesh_minimal.F90::imb_mesh_driver` | STL ingestion and scan capability is confirmed. Task 1 copies no source. The unique exception reserved for Task 5 is the user's authorization to port/adapt the `imb_mesh_y.cpp` scan mathematics/method without provenance boilerplate. That exception authorizes neither old COAST Fortran nor any other COAST C++. |
| Restart | `SRC.Coast/finish.F90::finish`, `SRC.Coast/app/coast_legacy_driver.F90` restart read/write schedule | Restart read/write and mesh-identity-aware runtime state are part of replacement scope. HUNDUN uses its own validated, transactional accepted-state format and publication. |
| Visit | `SRC.Coast/vtk.F90::{vtk,vtkblocks}` | VTK pieces and a VisIt index are emitted as scheduled snapshots. HUNDUN owns its own `io_visit` format, staging, and collective decision. |
| screen | `SRC.Coast/output.F90::output`, `SRC.Coast/coast_screen_summary.F90::coast_screen_write_step_summary` | Human-readable per-step diagnostics, extrema, flow, timing, and solver summaries are replacement requirements, not solver-state authorities. |
| ICCG lifecycle | `SRC.Coast/cgsol.F90::{cgsol,cgsol_multibrick}`, `SRC.Coast/module_coast_pressure_solver.F90::{coast_pressure_solve,coast_pressure_solve_legacy_iccg}`, `SRC.Coast/numerics/coast_numerics_wrappers.F90::solve_pressure_cg` | COAST assembles the current pressure system, selects/calls a legacy ICCG path, iterates, and records residual/iterations. This is lifecycle evidence only. HUNDUN uses the separate `SymbolicPlan`, `NumericState`, `HierarchyState`, and `SolverWorkspace` contract above. |

The fixed architecture also records that COAST has no periodic-boundary product capability;
the inspected legacy boundary tables mark a periodic condition code but do not provide a
v0.4-grade periodic product path. Consequently COAST is eligible only for nonperiodic short
performance comparison. This limitation, like its SIMPLE loop, is a reference-selection rule
and never a HUNDUN v0.4 algorithm requirement. HUNDUN still implements both periodic and
nonperiodic boundary plans, and long periodic cylinder statistics compare to experiments, not
to COAST.

## 6. Adoption boundary

- OpenFOAM contributes only public PISO mathematics and lifecycle questions under a strict
  GPL reference barrier.
- AMReX, IncFlo, and AMReX-Hydro contribute permissively licensed public ideas, but every
  HUNDUN type, equation assembly, cache key, schedule, and test remains an independent
  implementation. Their AMR algorithms are not part of v0.4.
- COAST contributes an observed replacement checklist and hot-cycle/resource questions. Old
  COAST Fortran is never copied or translated, and its SIMPLE/ICCG layout is not a target
  architecture. Task 1 itself copies no COAST source.
- The sole reuse exception is explicitly reserved for Task 5: the user authorizes porting or
  adapting the `imb_mesh_y.cpp` scanning mathematics/method without provenance boilerplate.
  The authorization does not extend to other COAST C++ or to any old COAST Fortran.
- Any implementation task that cannot state its HUNDUN owner, full identity tuple, lifetime,
  invalidation rule, and single writer has not satisfied this frozen survey.

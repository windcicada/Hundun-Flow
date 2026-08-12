# HUNDUN-FLOW Post-Task-11 Semantic-Port Architecture Design

**Date:** 2026-08-08  
**Status:** User-approved architecture; implementation begins only after Task
11 is accepted and the authoritative Stage 3 plan is amended.  
**Scope:** Stage 3 Tasks 12--21 and post-Stage-3 performance foundations.  
**Controlling principle:** Preserve HUNDUN-FLOW ownership, contracts and
architecture while independently reimplementing selected mature algorithms
and design semantics from public references.

## 1. Decision

HUNDUN-FLOW will not adopt AMReX, OpenFOAM, Basilisk or another solver as a
runtime dependency or backend. It will remain a project-owned C++ solver with
its existing topology, geometry, fields, execution, MPI, linear algebra,
transaction, checkpoint and diagnostic contracts.

The project will use a **semantic-port** process:

1. identify a mature public algorithm or architectural invariant;
2. freeze the mathematical behavior, inputs, outputs and mutations in a
   HUNDUN-owned design note;
3. write an independent RED oracle against the HUNDUN interface;
4. implement the minimum behavior using HUNDUN naming, data layouts and coding
   style without copying upstream source text, comments, control flow or ABI;
5. verify scientific behavior, deterministic MPI behavior and performance
   counters;
6. retain a provenance note identifying the public paper/document/repository
   used for comparison.

Direct source copying is outside this design, including code under permissive
licenses. A future direct BSD-licensed transplant would require a separate
user decision, retained upstream copyright/license notices and an isolated
legal review.

## 2. Goals

- deliver a usable constant-density static-IBM Stage 3 application soon after
  Task 11 rather than waiting until the former Task 19 position;
- reduce reimplementation risk by using mature operator, model and data-flow
  semantics;
- keep one numerical authority for immersed reconstruction, pressure,
  conservative flux, force and final-state reporting;
- separate regular/background kernels from sparse irregular/interface work;
- ensure all time-step geometry, stencils and communication plans are static;
- make attempt-local coefficients and rollback state explicit;
- establish exact counters before performance tuning;
- preserve Apache-2.0 project ownership and a clean independent source history.

## 3. Non-goals

This design does not authorize:

- vendoring or linking AMReX, OpenFOAM, Basilisk or IncFlo;
- replacing `MeshTopology`, `MeshGeometry`, `FlowState`, checked views,
  `ExecutionContext` or project-owned MPI/linear algebra;
- introducing `MultiFab`, `EB2`, `MLMG`, OpenFOAM dictionaries, object
  registries or run-time selection macros;
- rewriting the accepted Task 11 LFP-GCIBM operator into a true cut-cell
  method;
- AMR, moving bodies, FSI, rank-changing restart or production GPU support in
  Stage 3;
- adding damping, filtering, case tuning, relaxed scientific thresholds or
  additional PISO correctors;
- any 96-cubed test.

## 4. Public References and Legal Boundary

### 4.1 Primary architectural references

| Reference | Intended use | License/source rule |
|---|---|---|
| [AMReX](https://github.com/AMReX-Codes/amrex) and [EB documentation](https://amrex-codes.github.io/amrex/docs_html/EB.html) | EB support levels, sparse irregular data, regular/cut/covered dispatch, conservative redistribution, execution and performance concepts | BSD-3-Clause; semantic reimplementation only under this design |
| [IncFlo](https://github.com/AMReX-Fluids/incflo) | incompressible projection organization, EB/no-EB separation, convergence inventories | BSD-3-Clause; semantic reimplementation only |
| [OpenFOAM](https://github.com/OpenFOAM/OpenFOAM-dev) | WALE model boundary, explicit flow composition concepts, scheduled diagnostics | GPL-3.0; equations and public behavior only, no source copying |
| [OpenFOAM WALE](https://cpp.openfoam.org/v13/classFoam_1_1LESModels_1_1WALE.html) | independent WALE API and tensor oracle comparison | GPL-3.0 implementation; no source copying |
| [Basilisk `embed.h`](https://basilisk.fr/src/embed.h) | topology-aware embedded operators, boundary gradients, forces and flux consistency | GPL-3.0; mathematical comparison only |
| [ghost-cell-IBM](https://github.com/SpencerSchwart/ghost-cell-IBM) | hybrid ghost momentum/cut-pressure design comparison | no clear repository license found; no source copying |
| [NASSLARD2D](https://github.com/University-CFD-group/NASSLARD2D) | small laminar cylinder/projection validation ideas | GPL-3.0; no source copying |
| [3D-NS-FSI](https://github.com/Frederik-Kris/3D-NS-FSI) | image-point and moving-boundary test ideas | no clear repository license found; no source copying |

No exact public repository named `3D-NS-FSI-WLSQ` was located during the
2026-08-08 survey. The closest public repository above documents trilinear
image-point interpolation. It must not be represented as an authoritative
WLSQ source.

### 4.2 Clean-source rules

- Implementation workers receive the HUNDUN interface, equations, invariants,
  mutation list and test fixtures, not copied upstream functions.
- Upstream identifiers, comments, data layouts, macros, error messages and
  control flow must not be translated or mechanically refactored.
- GPL and unlicensed repositories are read-only scientific references.
- AMReX/IncFlo source is also read-only for normal semantic-port work even
  though BSD-3 permits redistribution.
- Every semantic port records the exact upstream URL, revision or document
  date, the behavior adopted, the behavior rejected and the independent tests.
- Public provenance acknowledges algorithms and papers without claiming
  upstream source reuse.
- No public `REUSE_MANIFEST.yaml` is added, and the existing `NOTICE` contract
  remains unchanged.

## 5. Architectural Model

```text
applications/
  Stage3Driver: construction, restart, retry, scheduling and output

flow/
  Stage3FlowComposition
  DensityAttemptAdapter
  WaleAttemptWorkspace
  TrialTransaction

finite_volume/
  BackgroundOperator
  ImmersedBoundaryAuthority
  optional ConservativeFluxCorrectionPlan

immersed/
  ImmersedSurface / SurfaceQuery / ImmersedDomain
  ActiveCellLayout / ActiveBoundaryLayout
  GhostStencilPlan / WallQuadraturePlan
  LocalFlowPatternTransform
  read-only geometry-support views

runtime/
  fields / checked views / MPI / halo / execution / checkpoint codec

diagnostics/
  read-only providers / exact counters / scheduled collection
```

Geometry decides **where** computation is legal. Operator plans decide **how**
it is discretized. Flow composition decides **when** each operation executes.
Diagnostics observe published summaries and counters but never become a
second numerical authority.

## 6. AMReX Semantics to Port

### 6.1 Geometry support levels

Expose logical support tiers through existing HUNDUN objects rather than
building a parallel EB database:

```text
basic:
  cell region, active ownership, interface connectivity

volume:
  basic plus active volume and cell centroid facts

full:
  volume plus face support, wall position/area/normal and reconstruction plans
```

The implementation may introduce small read-only view types when a real
consumer exists. It must not duplicate `ImmersedDomain`, `GhostStencilPlan`
or `WallQuadraturePlan` state.

### 6.2 Regular/irregular execution separation

Use existing active layouts to materialize deterministic owned-first work
partitions:

```text
regular active cells    -> branch-free background kernels
interface active cells  -> sparse immersed authority kernels
inactive cells          -> never read as physical state
```

The regular path must remain bitwise unchanged when IBM/WALE is absent. Static
partitioning occurs during construction, not during a time-step.

### 6.3 Static plans and attempt-local values

```text
static:
  geometry, IDs, ownership, donors, coefficients, halo plan, fingerprints

attempt-local:
  donor values, rho_wall, D_wall, u_lag, nu_t, mu_sgs, mu_eff, corrections
```

All attempt-local storage is discarded on rollback. No WALE coefficient,
temporary wall value or provisional force enters committed `FlowState`.

### 6.4 Conservative redistribution

AMReX-style redistribution is not part of the immediate MVP. A future
HUNDUN-owned `ConservativeFluxCorrectionPlan` is permitted only after a
mutation-sensitive RED proves a real cut-cell conservation or stability
defect.

When authorized it must:

- consume the same final face-flux provenance as transport and PISO;
- compute an explicit local conservation defect;
- redistribute the defect with deterministic nonnegative geometric weights;
- close the global correction to roundoff;
- preserve 1/2/4-rank decomposition invariance;
- expose exact correction, payload and reduction counters.

It must not silently filter the solution or become a general stabilizer.

## 7. OpenFOAM and Basilisk Semantics to Port

### 7.1 WALE

Implement the standard WALE mathematics independently:

- a backend-neutral tensor kernel;
- an explicit `WaleModel` interface;
- attempt-local `nu_t` and `mu_sgs`;
- one evaluation per trial after `rho_attempt` and `u_lag` are frozen;
- the same `mu_eff` in predictor, both correctors, final residual and wall
  force;
- no virtual calls or string lookup inside kernels;
- no wall function in Stage 3.

OpenFOAM is used to check model responsibility and published formula behavior,
not implementation structure.

### 7.2 Embedded pressure, force and topology

Basilisk supplies independent comparison points for:

- topological-connectivity guards before higher-order interpolation;
- embedded Dirichlet/Neumann gradient semantics;
- pressure and viscous surface-force decomposition;
- conservative face and wall-flux closure.

HUNDUN retains its full-strain low-Mach stress, unique pressure authority,
shared wall gradients and signed-force conventions accepted by Task 11.

## 8. Stage 3 Usable MVP

Immediately after Task 11, the first deliverable is a user-runnable vertical
slice:

- schema v3 constant-density flow;
- static STL, outside/inside LFP-GCIBM where already accepted;
- two-corrector PISO and shared final `FaceMassFlux`;
- pressure/viscous/total force and consistency report;
- adaptive retry and collective rollback;
- 1/2/4-rank execution;
- the same `hundun` executable via an early Stage 3 driver;
- Checkpoint v3 IBM-only continuation;
- minimal continuity, residual, wall-penetration, force and failure
  diagnostics;
- fast 12/24 screens and the retained 12/24/48 formal selectors;
- no 96-cubed execution.

The MVP excludes WALE, material density, ideal gas, combined IBM+WALE, AMR,
moving bodies, production GPU, rank-changing restart and multigrid.

## 9. Task 12--21 Amendment

The current plan places the first Stage 3 application driver too late. The
authoritative plan should be amended after Task 11 acceptance to use these
clusters:

1. **Task 19A -- constant IBM driver MVP**
   - construct only accepted constant-density IBM capabilities;
   - provide validate/run/output and no-restart execution first.
2. **Task 17A -- Checkpoint v3 codec and IBM-only continuation**
   - isolated byte protocol, transactional read and exact identity checks.
3. **Task 18A -- minimal diagnostics and counter substrate**
   - read-only module summaries and zero-cost disabled path.
4. **Task 12 -- WALE core**
   - independent tensor, exact-zero, y-cubed, identity and MPI tests.
5. **Task 13 -- body-fitted WALE vertical slice**
   - variable face viscosity and frozen attempt coefficients.
6. **Task 19B -- add `none/wale` to the driver**.
7. **Task 14 -- material-density IBM**.
8. **Task 15 -- ideal-gas IBM**.
9. **Task 16 -- combined IBM+WALE hard gate**.
10. **Task 17B -- WALE/density/combined continuation**.
11. **Task 18B -- complete Stage 3 diagnostics**.
12. **Task 19C -- complete all driver combinations**.
13. **Task 20 -- performance evidence and capability ledger**.
14. **Task 21 -- final acceptance**.

The scientific join remains:

```text
Task 11
  -> driver/restart MVP
  -> Task 12 -> Task 13 -----------+
  -> Task 14 -> Task 15 -----------+-> Task 16
  -> Task 17A / Task 18A ----------+-> complete driver -> Task 20 -> Task 21
```

Task 14 and Task 15 remain sequential because they share flow composition
files and closure state. Checkpoint protocol work, diagnostics adapters and
WALE tensor work are bounded independent modules, but implementation-worker
concurrency must still obey repository `AGENTS.md` and resource limits.

## 10. Test and Acceptance Strategy

### Development tiers

```text
fast:
  unit, header, policy, transaction, 12/24 task-focused numerical screens

screen:
  stable two-level direction checks and selected 1/2-rank tests

acceptance:
  only stable candidates; required 12/24/48 and 1/2/4 matrices
```

- Every tier uses the same product numerical path.
- No repair mechanically repeats Debug/Release/ASan/UBSan and every numerical
  selector.
- ASan/UBSan use task-focused small cases; Release owns long numerical
  acceptance.
- Independent low-memory unit/header/configuration tests may run in parallel.
- Large numerical jobs receive explicit CPU/MPI/memory groups.
- No two internally 96-thread numerical jobs run concurrently unless the host
  has at least twice that effective CPU capacity.
- Every detached long job records exact HEAD, dirty-diff hash, binary SHA-256,
  command, environment, log hash and exit status.
- 96-cubed is permanently excluded.

### Semantic-port RED requirements

Each port must name mutations that would distinguish an independent
implementation from a plausible wrong one. At minimum cover:

- wrong geometry support or connectivity;
- wrong conservation sign or missing redistribution closure;
- stale attempt values after retry;
- a second model evaluation inside a corrector;
- inactive-cell reads;
- decomposition-dependent ordering;
- diagnostic collection changing state or counters;
- fallback/limiter behavior not authorized by the design.

## 11. Performance Roadmap

Performance work is evidence-driven and ordered:

1. exact counters for queries, donors, halo bytes/messages, allocations,
   model evaluations, matvecs and reductions;
2. construction-time caching and zero time-step geometry work;
3. regular/interface work partitioning;
4. allocation-free attempt workspaces;
5. packed donor halo exchange with deterministic global ownership;
6. loop fusion only when counters/profiles prove memory traffic dominates;
7. after Stage 3, a project-owned structured geometric-multigrid pressure
   preconditioner inspired by AMReX concepts but independently implemented;
8. AMR and GPU remain separate later-stage designs.

Wall-clock thresholds are not portable scientific gates. Exact work counts,
payloads and allocation identities are gates; time, RSS and throughput are
recorded with compatibility metadata.

## 12. Designs Explicitly Rejected

- copying AMReX `MultiFab`, `EB2`, `MLMG` or ParmParse;
- adopting OpenFOAM run-time selection macros or global object registry;
- Basilisk-style global operator macro overrides;
- a parallel second IBM geometry database;
- widespread `my-*` forks of core operators;
- per-wall-point dynamic Eigen solves in production apply paths;
- per-step STL query, classification, QR or donor search;
- global value Allgather during every time-step;
- cut-cell redistribution without a conservation/stability RED;
- changing Task 11 product operators merely to resemble an upstream project.

## 13. Entry and Exit Conditions

### Entry

- Task 11 is formally accepted on an exact HEAD;
- its retained formal selectors, decomposition, rollback, collective failure,
  PISO and force evidence are sealed;
- the authoritative Stage 3 plan records this amendment and preserves the old
  plan/history;
- no Task 11 background job or stale worker remains.

### Exit

- the early constant-IBM driver/restart MVP is user runnable;
- Tasks 12--16 complete WALE and all density compositions;
- Checkpoint/diagnostics/driver cover every legal presence combination;
- exact counters show no per-step static geometry rebuild;
- capability ledger distinguishes accepted, deferred and out-of-scope rows;
- Stage 3 final acceptance passes without 96-cubed execution.

## 14. First Three Actions After Task 11

1. amend the authoritative Stage 3 plan and capability ledger using this
   architecture while preserving historical plans;
2. write the Task 19A constant-IBM driver MVP RED and construction-order
   contract;
3. freeze Task 17A Checkpoint v3 IBM-only protocol and Task 18A minimal
   diagnostics interfaces so they can be implemented as bounded modules.


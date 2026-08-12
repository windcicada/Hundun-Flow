# HUNDUN-FLOW v0.4 Cartesian Performance Architecture Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an independently implemented, CPU-first HUNDUN-FLOW v0.4 Cartesian CFD core and release it only after the Re=3900 cylinder numerical-correctness and COAST-performance joint gate returns `ACCEPT`.

**Architecture:** Keep v0.3 as an exact tracked-tree oracle while v0.4 is a separate flat-source product line. v0.4 freezes field, geometry, boundary, communication, operator, solver, and execution plans before the hot loop; it uses padded SoA fields, direction-separated conservative face fluxes, persistent workspaces, explicit revisions, native Cartesian multigrid, exactly two PISO correctors, and static EB/WALE plans.

**Tech Stack:** C++17, CMake 3.21+, MPI 3, POSIX/Linux, pthread/OpenMP-compatible persistent CPU teams, FP64 numerical authority, optional HYPRE Struct fallback, CTest, ASan, UBSan, perf-compatible counters.

## Global Constraints

- Performance is the first implementation constraint, but no optimization may reduce the two PISO correctors, weaken convergence/scientific thresholds, add case-specific damping/filtering, or change conservation, rollback, MPI, force-sign, and Restart contracts.
- v0.4 supports only `UniformCartesianPlan` and `StretchedCartesianPlan`; body-fitted, multiblock curvilinear, AMR, moving EB, and production GPU backends are out of scope.
- Linux CPU is the production target. Default placement is one MPI rank per NUMA domain with a persistent thread team; MPI-only remains supported.
- Product implementation may use public mathematics, data-layout, and lifecycle ideas from the fixed OpenFOAM, AMReX, IncFlo, and AMReX-Hydro revisions and read-only COAST location recorded in the approved spec. Do not copy or translate GPL or COAST source.
- Keep all v0.4 implementation files flat in `versions/v0.4/src` with stable responsibility prefixes. Public headers remain in `versions/v0.4/include/hundun`.
- Runtime field registration is allowed only before `FieldSchema::freeze()`. The frozen hot loop has no heap allocation, string lookup, geometry/donor search, per-cell virtual call, implicit MPI, undeclared write, or whole-field trial copy.
- FP64 is authoritative. Mixed precision is disabled in the initial implementation and may only be added by a later approved plan with an FP64 residual check and fallback.
- Preserve `/home/wyf/code_dev/.worktrees/hundun-flow-cylinder3900-stage4`, its dirty/untracked state, and all benchmark evidence. Never reset, clean, overwrite, or commit that worktree.
- Candidate tests run in this order: focused, `24^3/20-step`, full `480x480x48/64-rank` 2-step pairing, frozen-candidate full-grid 20-step pairing, then long statistics. Do not start long statistics before candidate freeze.
- There is one release node only: the Re=3900 numerical-correctness, robustness, COAST-performance, physical-accuracy, and provenance conjunction. Internal checkpoints never create alpha, beta, RC, or release tags.
- Use TDD for every product change. Each task begins with a failing contract or mutation-sensitive RED and ends with fresh focused verification and a DCO commit.

---

## File Map

The following map is normative. New v0.4 implementation files remain in one flat directory; a task may split a file only by adding another file with the same prefix.

| Public header | Flat implementation files | Responsibility |
| --- | --- | --- |
| `v04_types.hpp` | header-only POD contracts | indices, extents, spans, geometry/model tags |
| `v04_status.hpp` | `platform_status.cpp` | compact hot status and cold-path explanation |
| `v04_case.hpp` | `app_case.cpp`, `app_validate.cpp` | `CaseSpec` and `ValidatedModel` |
| `v04_storage.hpp` | `storage_schema.cpp`, `storage_arena.cpp`, `storage_state.cpp` | frozen fields, padded SoA, arena, revisions, time layers |
| `v04_mesh.hpp` | `mesh_cartesian.cpp`, `mesh_decomposition.cpp`, `mesh_tiles.cpp` | Cartesian geometry, one patch/rank, CPU tiles |
| `v04_platform.hpp` | `platform_cpu.cpp`, `platform_numa.cpp`, `platform_team.cpp` | ISA, NUMA placement, persistent thread team |
| `v04_comm.hpp` | `comm_plan.cpp`, `comm_persistent.cpp` | merged persistent halo and ghost certification |
| `v04_boundary.hpp` | `boundary_compile.cpp`, `boundary_apply.cpp` | `BoundarySpec` to immutable `BoundaryPlan` |
| `v04_discretization.hpp` | `disc_scheme.cpp`, `disc_gradient.cpp`, `disc_flux.cpp`, `disc_diffusion.cpp` | frozen schemes and conservative Cartesian kernels |
| `v04_execution.hpp` | `integration_graph.cpp`, `integration_transaction.cpp`, `integration_time.cpp` | frozen graph, stage consensus, attempt commit/rollback |
| `v04_linear.hpp` | `linear_lifecycle.cpp`, `linear_krylov.cpp`, `linear_cartesian_mg.cpp`, `linear_hypre_struct.cpp` | four-layer solver state, Krylov, native MG, fallback |
| `v04_flow.hpp` | `operator_momentum.cpp`, `operator_pressure.cpp`, `integration_piso.cpp` | equation authority, `rAU/rAtU/HbyA/phiHbyA`, two PISO corrections |
| `v04_eb.hpp` | `eb_geometry.cpp`, `eb_stencil.cpp`, `eb_pressure.cpp`, `eb_force.cpp` | static EB topology, interface correction, final-state forces |
| `v04_model.hpp` | `model_transport.cpp`, `model_wale.cpp`, `model_contribution.cpp` | `mu_eff`, WALE, density/scalar contributions |
| `v04_service.hpp` | `service_restart.cpp`, `service_diagnostics.cpp`, `service_performance.cpp` | committed snapshots, evidence, resource budgets |
| `v04_app.hpp` | `app_main.cpp`, `app_driver.cpp` | product executable and cold-path orchestration |

Tests are flat in `versions/v0.4/tests` and use the same prefixes. Case assets are under `versions/v0.4/cases`; no generated evidence is committed there.

---

### Task 1: Freeze the Public-Reference and COAST Lifecycle Survey

**Files:**
- Create: `docs/references/2026-08-12-hundun-flow-v0.4-hot-loop-lifecycle-survey.md`
- Create: `docs/references/2026-08-12-hundun-flow-v0.4-adoption-ledger.tsv`
- Create: `docs/architecture/v0.4-target-hot-loop.md`

**Interfaces:**
- Consumes: fixed upstream revisions from the approved design and read-only `/home/wyf/code_dev/Coast_software`.
- Produces: source/symbol/revision evidence, stage mapping, adopted/reworked/rejected decisions, and a target lifecycle that every later task cites.

- [ ] **Step 1: Record exact reference identities and licenses.**

  Use a read-only external cache, never a product source directory:

  ```bash
  REF=/home/wyf/code_dev/.reference-cache/hundun-v04-20260812
  mkdir -p "$REF"
  git clone --filter=blob:none https://github.com/OpenFOAM/OpenFOAM-dev.git "$REF/OpenFOAM-dev"
  git -C "$REF/OpenFOAM-dev" checkout --detach b9da51ab0673423aa2af6a45a72a3fbec9c66f9f
  git clone --filter=blob:none https://github.com/AMReX-Codes/amrex.git "$REF/amrex"
  git -C "$REF/amrex" checkout --detach 59d066aab774bc388cc6ed944f7beaf645607ed3
  git clone --filter=blob:none https://github.com/AMReX-Fluids/incflo.git "$REF/incflo"
  git -C "$REF/incflo" checkout --detach 7307d8725c2a538f09cafbeacbfeb63e0fb11d22
  git clone --filter=blob:none https://github.com/AMReX-Fluids/amrex-hydro.git "$REF/amrex-hydro"
  git -C "$REF/amrex-hydro" checkout --detach e49df248aabd2cc11865eb5be734a2f5f2f65ee5
  ```

  If a directory already exists, verify its origin and fetch only the named commit instead of recloning. Run `git rev-parse HEAD` in each checkout and `realpath /home/wyf/code_dev/Coast_software`; fail if a revision or origin differs. The ledger columns must be:

  ```text
  project\trevision\tlicense\tfile\tsymbol\tproblem\tdecision\thundun_owner
  ```

- [ ] **Step 2: Trace OpenFOAM PISO/PIMPLE assembly and reuse.**

  Record exact files and call order for momentum assembly, `rAU/rAtU`, `HbyA`, `phiHbyA`, pressure correction, non-orthogonal correction, and final `U/phi`. State explicitly which values survive a corrector and which are invalidated by coefficient or pressure changes.

- [ ] **Step 3: Trace AMReX/IncFlo EB, FillPatch, operator, and MLMG lifetimes.**

  Record exact files/symbols for geometry/topology caching, MultiFab/box layout, fill-patch ordering, coefficient updates, hierarchy reuse, and workspace ownership. Mark AMR-only mechanisms `rejected-v0.4`.

- [ ] **Step 4: Trace COAST hot-step ownership without copying source.**

  Use `rg` only for grid/IBM setup, SIMPLE/ICCG assembly, coefficient/preconditioner workspace, halo, force, output, and timing boundaries. The report paraphrases lifecycle and array-access facts; it contains no copied routine body, translated loop, or private data.

- [ ] **Step 5: Map all projects to the HUNDUN target lifecycle.**

  The architecture document must contain this exact stage skeleton and fill each row with authority, inputs, outputs, cache revision, communication, and workspace:

  ```text
  accept -> halo -> derived fields -> momentum numeric update -> predictor
         -> PISO-1 pressure/final-flux update -> PISO-2 pressure/final-flux update
         -> final gradient/force -> consensus -> commit -> services
  ```

- [ ] **Step 6: Verify provenance and commit.**

  Run:

  ```bash
  ! rg -n 'TB[D]|TO[D]O|FIXM[E]|XX[X]' docs/references/2026-08-12-hundun-flow-v0.4-* docs/architecture/v0.4-target-hot-loop.md
  git diff --check
  ```

  Expected: no placeholder; every adoption-ledger row has a revision, file, symbol, decision, and HUNDUN owner. Commit `docs: map v0.4 solver lifecycles` with DCO.

### Task 2: Establish Version Isolation and a Buildable v0.4 Skeleton

**Files:**
- Modify: `CMakeLists.txt`
- Modify: `CMakePresets.json`
- Modify: `README.md`
- Modify: `VERSION`
- Create: `versions/v0.3/**` from the exact tracked tree only
- Create: `versions/v0.3/FROZEN_SOURCE.json`
- Create: `versions/v0.4/CMakeLists.txt`
- Create: `versions/v0.4/include/hundun/v04_types.hpp`
- Create: `versions/v0.4/include/hundun/v04_status.hpp`
- Create: `versions/v0.4/src/platform_status.cpp`
- Create: `versions/v0.4/src/app_main.cpp`
- Create: `versions/v0.4/tests/CMakeLists.txt`
- Create: `versions/v0.4/tests/platform_status_test.cpp`

**Interfaces:**
- Consumes: governance-repository commit `4ae4832ad00b5e4d1129ee978e2e49cbb33bb7ef`, tree `0a59ffa9d71cee79ef561c49657c9e90c5948e3e`.
- Produces: root `HUNDUN_SOURCE_VERSION={v0.4,v0.3}`, default v0.4 executable, independently buildable frozen v0.3, POD base types, and `Status {StatusCode code; uint32_t detail;}`.

- [ ] **Step 1: Add a configure RED for version selection.**

  Add a CTest/CMake script that configures default, explicit v0.4, invalid `v9`, and explicit v0.3. Expected initial failure: `HUNDUN_SOURCE_VERSION` is unknown.

- [ ] **Step 2: Export the exact v0.3 tracked tree without reading dirty files.**

  Before export, require:

  ```bash
  test "$(git -C /home/wyf/code_dev/.worktrees/hundun-flow-cylinder3900-stage4 rev-parse HEAD)" = 4ae4832ad00b5e4d1129ee978e2e49cbb33bb7ef
  test "$(git -C /home/wyf/code_dev/.worktrees/hundun-flow-cylinder3900-stage4 rev-parse HEAD^{tree})" = 0a59ffa9d71cee79ef561c49657c9e90c5948e3e
  git -C /home/wyf/code_dev/.worktrees/hundun-flow-cylinder3900-stage4 archive 4ae4832ad00b5e4d1129ee978e2e49cbb33bb7ef | tar -x -C versions/v0.3
  ```

  `FROZEN_SOURCE.json` records source repository realpath, HEAD, tree, archive SHA-256, and import date. It must not list or hash worktree dirty files.

- [ ] **Step 3: Implement root dispatch and the minimal v0.4 status API.**

  Use these public contracts:

  ```cpp
  struct Int3 { std::int32_t x, y, z; };
  struct Real3 { double x, y, z; };
  template<class T> struct Span { const T* data; std::size_t size; };
  using NumaNode = std::int32_t;
  enum class CartesianKind : std::uint8_t { uniform, stretched };
  enum class DensityKind : std::uint8_t { constant, material, ideal_gas };
  enum class TurbulenceKind : std::uint8_t { laminar, wale };
  enum class StatusCode : std::uint16_t { ok, invalid_case, invalid_plan,
    allocation_failure, mpi_failure, numerical_failure, rejected_step };
  struct Status { StatusCode code{StatusCode::ok}; std::uint32_t detail{};
    constexpr explicit operator bool() const noexcept { return code == StatusCode::ok; } };
  ```

  Root configuration defaults to v0.4. Both version lines compile as C++17; v0.4 uses its own non-owning POD view instead of requiring `std::span`.

- [ ] **Step 4: Pass both build lanes.**

  Run:

  ```bash
  cmake -S . -B build/v04-debug -DHUNDUN_SOURCE_VERSION=v0.4 -DHUNDUN_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Debug
  cmake --build build/v04-debug -j2
  ctest --test-dir build/v04-debug --output-on-failure
  cmake -S . -B build/v03-release -DHUNDUN_SOURCE_VERSION=v0.3 -DHUNDUN_BUILD_TESTS=OFF -DCMAKE_BUILD_TYPE=Release
  cmake --build build/v03-release -j2
  ```

  Expected: v0.4 status test passes; v0.3 `hundun --version` reports its frozen version. Commit `build: isolate v0.3 and v0.4 source lines` with DCO.

### Task 3: Compile CaseSpec into ValidatedModel and Frozen FieldSchema

**Files:**
- Create: `versions/v0.4/include/hundun/v04_case.hpp`
- Create: `versions/v0.4/include/hundun/v04_storage.hpp`
- Create: `versions/v0.4/src/app_case.cpp`
- Create: `versions/v0.4/src/app_validate.cpp`
- Create: `versions/v0.4/src/storage_schema.cpp`
- Create: `versions/v0.4/tests/app_case_test.cpp`
- Create: `versions/v0.4/tests/storage_schema_test.cpp`

**Interfaces:**
- Produces: `Result<ValidatedModel> validate_case(const CaseSpec&)`; `Result<FieldId> FieldSchema::declare_field(FieldDescriptor)`; `Status FieldSchema::freeze()`; immutable descriptor lookup by `FieldId`.

- [ ] **Step 1: Write REDs for invalid physics and post-freeze mutation.**

  Test duplicate names, zero/negative cell counts, body-fitted selection, missing pressure boundary closure, a non-two PISO count for the Re=3900 profile, and `declare_field` after freeze. Expected: compile failure before the APIs exist.

- [ ] **Step 2: Define stable identifiers and validation results.**

  Use integer IDs and explicit cold-path results:

  ```cpp
  template<class T> struct Result { Status status; std::optional<T> value; };
  using FieldId = std::uint32_t;
  enum class ScalarKind : std::uint8_t { float64 };
  enum class Centering : std::uint8_t { cell, face_x, face_y, face_z };
  enum class FieldRole : std::uint8_t { state, authority, derived, scratch };
  struct CartesianSpec { CartesianKind kind; Int3 cells; Real3 origin, length; };
  struct PhysicsSpec { DensityKind density; TurbulenceKind turbulence; bool immersed; };
  struct FieldDescriptor { std::string name; ScalarKind scalar; Centering centering;
    std::uint8_t components; std::uint8_t ghost_width; FieldRole role; };
  struct ValidatedModel { CartesianSpec mesh; PhysicsSpec physics;
    std::vector<FieldDescriptor> fields; std::uint8_t piso_correctors; };
  ```

  Strings remain cold-path schema data; frozen plans store only IDs and POD descriptors.

- [ ] **Step 3: Implement validation and schema freeze.**

  Validation must reject every unsupported geometry at startup, register velocity/pressure/density and x/y/z face flux authorities once, and return a status rather than partially constructing a model.

- [ ] **Step 4: Verify and commit.**

  Run `ctest --test-dir build/v04-debug -R 'v04_(case|schema)' --output-on-failure`. Expected: all invalid cases return exact status codes and no schema mutation succeeds after freeze. Commit `feat(v0.4): freeze validated cases and field schema`.

### Task 4: Implement NUMA-Local Arena, Padded SoA Views, and Revisions

**Files:**
- Modify: `versions/v0.4/include/hundun/v04_storage.hpp`
- Create: `versions/v0.4/src/storage_arena.cpp`
- Create: `versions/v0.4/src/storage_state.cpp`
- Create: `versions/v0.4/tests/storage_arena_test.cpp`
- Create: `versions/v0.4/tests/storage_revision_test.cpp`
- Create: `versions/v0.4/tests/support_allocation_guard.cpp`

**Interfaces:**
- Consumes: frozen `FieldSchema`.
- Produces: `Arena::reserve(bytes, alignment, NumaNode)`; `FieldStorage::create(schema, PatchShape)`; `FieldView<T>` with x-contiguous indexing; `RevisionToken`; accepted/trial/history `LayerHandle` rotation.

- [ ] **Step 1: Write alignment, padding, authority, and allocation REDs.**

  Require 64-byte base alignment, x stride 1, padded y stride divisible by eight doubles, ghost/interior in one allocation, component-separated SoA, stale-view rejection in checked builds, and zero allocation events between `begin_step()` and `commit()/rollback()`.

- [ ] **Step 2: Implement one large arena and deterministic layout planning.**

  Calculate all offsets with checked `size_t` arithmetic before allocation. Each entry is POD:

  ```cpp
  struct PatchShape { Int3 interior; std::uint8_t ghost_width; };
  struct FieldLayout { std::size_t offset, x_stride, y_stride, z_stride;
    Int3 interior, allocated; std::uint8_t components, ghost_width; };
  struct RevisionToken { std::uint64_t value{}; };
  using LayerHandle = std::uint32_t;
  ```

- [ ] **Step 3: Implement handle rotation and local transaction logs.**

  `commit()` swaps layer handles and publishes pending revisions; `rollback()` discards pending revisions and restores only logged scalar/handle metadata. It never copies a complete field.

- [ ] **Step 4: Verify and commit.**

  Run Debug, ASan, and UBSan storage tests. Expected: no hot allocation, OOB, stale publish, or full-field copy counter. Commit `feat(v0.4): add padded SoA arena storage`.

### Task 5: Build CartesianGeometryPlan, MeshPatch, and CpuTile

**Files:**
- Create: `versions/v0.4/include/hundun/v04_mesh.hpp`
- Create: `versions/v0.4/src/mesh_cartesian.cpp`
- Create: `versions/v0.4/src/mesh_decomposition.cpp`
- Create: `versions/v0.4/src/mesh_tiles.cpp`
- Create: `versions/v0.4/tests/mesh_cartesian_test.cpp`
- Create: `versions/v0.4/tests/mesh_decomposition_mpi_test.cpp`

**Interfaces:**
- Produces: `CartesianGeometryPlan::uniform(...)`, `::stretched(...)`, `MeshPatch`, `CpuTileRange`, cell/face metrics, owned global box, and static nonuniform rank decomposition.

- [ ] **Step 1: Write geometry and decomposition REDs.**

  Cover exact uniform metrics, monotone stretched coordinates, face/cell volume identities, invalid coordinates, 1/2/4-rank no-gap/no-overlap ownership, periodic neighbor identity, and EB-weighted static split inputs.

- [ ] **Step 2: Implement the two geometry variants without per-cell polymorphism.**

  Use a tagged immutable plan selected once:

  ```cpp
  struct AxisMetrics { Span<double> cell, face, delta, inverse_delta; };
  struct CartesianGeometryPlan { CartesianKind kind; AxisMetrics x, y, z;
    RevisionToken geometry_revision; };
  ```

  Uniform kernels read scalar spacing; stretched kernels read one-dimensional metric arrays.

- [ ] **Step 3: Implement one patch per rank and variable CPU tiles.**

  Tiles partition only the owned patch; communication and field ownership remain patch-based. Default tile sizes are plan inputs, not case-hardcoded constants.

- [ ] **Step 4: Verify and commit.**

  Run unit tests and `mpiexec -n 1/2/4` decomposition tests. Commit `feat(v0.4): add Cartesian geometry and patch plans`.

### Task 6: Add CPU/NUMA/ISA Plan and Attempt Transactions

**Files:**
- Create: `versions/v0.4/include/hundun/v04_platform.hpp`
- Create: `versions/v0.4/include/hundun/v04_execution.hpp`
- Create: `versions/v0.4/src/platform_cpu.cpp`
- Create: `versions/v0.4/src/platform_numa.cpp`
- Create: `versions/v0.4/src/platform_team.cpp`
- Create: `versions/v0.4/src/integration_transaction.cpp`
- Create: `versions/v0.4/src/integration_time.cpp`
- Create: `versions/v0.4/tests/platform_cpu_test.cpp`
- Create: `versions/v0.4/tests/integration_transaction_mpi_test.cpp`

**Interfaces:**
- Produces: `CpuExecutionPlan`, startup ISA selection, persistent team, `TimeSchemePlan`, `AttemptTransaction`, and stage-local collective consensus.

- [ ] **Step 1: Write REDs for bounded specialization and collective rollback.**

  Require `scalar`, `avx2`, and `avx512` enum choices only; unsupported forced ISA returns `invalid_plan`. On 1/2/4 ranks, one-rank failure must make every rank reject the attempt with identical status and unchanged accepted revisions.

- [ ] **Step 2: Implement startup-only CPU discovery and placement validation.**

  Record logical CPUs, NUMA node, rank affinity, selected ISA, alignment, tile shape, and thread count in `CpuExecutionPlan`. No CPUID, affinity query, or team creation occurs after freeze.

- [ ] **Step 3: Implement persistent team and transaction state machine.**

  States are `idle -> active -> consensus -> committed|rolled_back`; illegal transitions return status. One collective is permitted at declared stage boundaries only.

- [ ] **Step 4: Verify and commit.**

  Run unit, 1/2/4-rank rollback, ASan, and forced-scalar tests. Commit `feat(v0.4): freeze CPU execution and transaction plans`.

### Task 7: Compile Merged Persistent Halo Communication

**Files:**
- Create: `versions/v0.4/include/hundun/v04_comm.hpp`
- Create: `versions/v0.4/src/comm_plan.cpp`
- Create: `versions/v0.4/src/comm_persistent.cpp`
- Create: `versions/v0.4/tests/comm_plan_test.cpp`
- Create: `versions/v0.4/tests/comm_persistent_mpi_test.cpp`

**Interfaces:**
- Consumes: `MeshPatch`, frozen field layouts, stage read sets, and field revisions.
- Produces: `CommunicationPlan::compile(...)`, `begin(stage, storage)`, `finish(stage, storage)`, merged peer buffers, persistent requests, and certified ghost revisions.

- [ ] **Step 1: Write payload, revision, and peer-only REDs.**

  On 1/2/4 ranks, require exact canonical ghost values, one merged message per peer/stage, payload equal to declared field regions, no all-gather, stale ghost rejection, and zero request/buffer allocation on repeated exchanges.

- [ ] **Step 2: Compile pack/unpack spans and persistent requests.**

  Store field ID, component, source span, destination span, peer, and expected revision in immutable arrays. Reject overlapping writes and locally owned ghost requests at compile time.

- [ ] **Step 3: Implement begin/interior/finish/boundary sequencing.**

  `begin` packs and starts requests; the caller executes interior tiles; `finish` waits, unpacks, and publishes ghost revision; then boundary tiles are eligible.

- [ ] **Step 4: Verify and commit.**

  Run 1/2/4-rank tests with allocation guard and message/byte counters. Commit `feat(v0.4): add merged persistent halo plans`.

### Task 8: Compile BoundarySpec and SchemePlan

**Files:**
- Create: `versions/v0.4/include/hundun/v04_boundary.hpp`
- Create: `versions/v0.4/include/hundun/v04_discretization.hpp`
- Create: `versions/v0.4/src/boundary_compile.cpp`
- Create: `versions/v0.4/src/boundary_apply.cpp`
- Create: `versions/v0.4/src/disc_scheme.cpp`
- Create: `versions/v0.4/tests/boundary_plan_test.cpp`
- Create: `versions/v0.4/tests/disc_scheme_test.cpp`

**Interfaces:**
- Produces: immutable `BoundaryPlan`, `SchemePlan`, compressed boundary spans, and static kernel tags.

- [ ] **Step 1: Write REDs for closure and hot dispatch.**

  Cover periodic pairing, velocity inlet/pressure outlet, no-slip wall, scalar zero-gradient, incompatible pressure closure, and a counter proving zero string comparisons during repeated application.

- [ ] **Step 2: Implement cold-path compilation.**

  Convert names and JSON values to field IDs, face spans, coefficient arrays, and enum kernel tags. Reject contradictory or uncovered physical faces before any field allocation.

- [ ] **Step 3: Freeze the initial scheme set.**

  Include central gradient, conservative linear face interpolation, bounded upwind/linear convection needed by existing validated profiles, orthogonal Cartesian diffusion, and BDF1/BDF2 time tags. No generic runtime expression tree is introduced.

- [ ] **Step 4: Verify and commit.**

  Run boundary/scheme tests and allocation guard. Commit `feat(v0.4): compile boundary and scheme plans`.

### Task 9: Implement Conservative Cartesian Kernels and Face-Flux Authority

**Files:**
- Modify: `versions/v0.4/include/hundun/v04_discretization.hpp`
- Create: `versions/v0.4/src/disc_gradient.cpp`
- Create: `versions/v0.4/src/disc_flux.cpp`
- Create: `versions/v0.4/src/disc_diffusion.cpp`
- Create: `versions/v0.4/tests/disc_manufactured_test.cpp`
- Create: `versions/v0.4/tests/disc_conservation_mpi_test.cpp`

**Interfaces:**
- Produces: direction-separated `FaceFluxView {x,y,z}`, gradient, divergence, convection, and diffusion kernels specialized by geometry/scheme/ISA tag.

- [ ] **Step 1: Write manufactured and conservation REDs.**

  Require exact constant/linear gradients, second-order smooth-grid trends, zero global divergence for periodic constant flux, pairwise internal-face cancellation, and 1/2/4-rank equality within FP64 reduction tolerance.

- [ ] **Step 2: Implement x-contiguous interior kernels.**

  Accept only POD views and preselected tags. Use `restrict`-equivalent nonaliasing contracts, contiguous x loops, and no bounds/variant checks inside the innermost loop.

- [ ] **Step 3: Fuse face reconstruction and mass-flux production.**

  Write each directional authority once and consume it for divergence and scalar transport. Add counters for bytes read/written and forbid a second writable face-flux field.

- [ ] **Step 4: Verify and commit.**

  Run manufactured, conservation, MPI, ASan/UBSan, and allocation tests. Commit `feat(v0.4): add conservative Cartesian kernels`.

### Task 10: Freeze the Execution Graph and Resource Contracts

**Files:**
- Modify: `versions/v0.4/include/hundun/v04_execution.hpp`
- Create: `versions/v0.4/src/integration_graph.cpp`
- Create: `versions/v0.4/tests/integration_graph_test.cpp`
- Create: `versions/v0.4/tests/integration_hot_path_test.cpp`

**Interfaces:**
- Consumes: field, boundary, communication, kernel, and transaction plans.
- Produces: `FrozenExecutionGraph::compile(StageSpec...)`, topologically ordered `StageNode` arrays, liveness-based workspace offsets, and per-stage byte/message/allocation budgets.

- [ ] **Step 1: Write graph-contract REDs.**

  Reject cycles, undeclared writes, two writers for one authority revision, read-before-produce, missing ghost dependency, overlapping scratch lifetimes, and implicit collective calls.

- [ ] **Step 2: Implement graph compilation and workspace liveness.**

  A node stores only IDs and function pointers selected at startup:

  ```cpp
  using StageId = std::uint32_t;
  struct KernelContext;
  struct RevisionRequirement { FieldId field; RevisionToken revision; };
  using KernelFn = Status (*)(KernelContext&) noexcept;
  struct WorkspaceSlice { std::size_t offset, bytes; };
  enum class CollectiveKind : std::uint8_t { none, stage_consensus };
  struct StageNode { StageId id; KernelFn kernel; Span<FieldId> reads, writes;
    Span<RevisionRequirement> ghosts; WorkspaceSlice scratch; CollectiveKind collective; };
  ```

- [ ] **Step 3: Enforce hot-path resource budgets.**

  With the allocation guard armed, execute a synthetic step 100 times. Expected: zero heap events, stable workspace address, stable message count, and byte counts no greater than the compiled contract.

- [ ] **Step 4: Verify and commit.**

  Run graph and hot-path tests. Commit `feat(v0.4): freeze execution graph and budgets`.

### Task 11: Implement Four-Layer Linear Lifecycle and Krylov Solvers

**Files:**
- Create: `versions/v0.4/include/hundun/v04_linear.hpp`
- Create: `versions/v0.4/src/linear_lifecycle.cpp`
- Create: `versions/v0.4/src/linear_krylov.cpp`
- Create: `versions/v0.4/tests/linear_lifecycle_test.cpp`
- Create: `versions/v0.4/tests/linear_krylov_mpi_test.cpp`

**Interfaces:**
- Produces: `SymbolicPlan`, `NumericState`, `HierarchyState`, `SolverWorkspace`, `CgSolver`, `FgmresSolver`, `SolveControl`, and `SolveReport`.

- [ ] **Step 1: Write lifecycle and nonsymmetric-solver REDs.**

  Require symbolic identity stability across numeric updates, no refill when coefficient revision is unchanged, persistent workspace across repeated solves, transactional growth failure, CG rejection of a declared nonsymmetric operator, and FGMRES convergence on a fixed nonnormal matrix at 1/2/4 ranks.

- [ ] **Step 2: Implement explicit identities and invalidation.**

  Use:

  ```cpp
  struct LinearIdentity { std::uint64_t topology, coefficients, hierarchy; };
  struct SolveControl { double atol, rtol; std::uint32_t max_iterations, restart; };
  struct SolveReport { Status status; std::uint32_t iterations;
    double initial_residual, final_residual; };
  ```

- [ ] **Step 3: Implement CG and restarted right-preconditioned FGMRES.**

  FGMRES uses two-pass modified Gram-Schmidt, Givens least squares, and an FP64 true-residual recomputation at restart and claimed convergence. All vectors come from `SolverWorkspace`.

- [ ] **Step 4: Verify and commit.**

  Run 1/2/4-rank solver tests, forced allocation failure, ASan, and UBSan. Commit `feat(v0.4): add persistent linear solver lifecycle`.

### Task 12: Implement NativeCartesianMG and HYPRE Struct Isolation

**Files:**
- Modify: `versions/v0.4/include/hundun/v04_linear.hpp`
- Create: `versions/v0.4/src/linear_cartesian_mg.cpp`
- Create: `versions/v0.4/src/linear_hypre_struct.cpp`
- Create: `versions/v0.4/tests/linear_cartesian_mg_test.cpp`
- Create: `versions/v0.4/tests/linear_cartesian_mg_mpi_test.cpp`

**Interfaces:**
- Consumes: Cartesian metrics, coefficient revisions, persistent communication.
- Produces: native full-coarsening MG, semi-coarsening/line relaxation for strong stretch, hierarchy reuse policy, and optional `HypreStructAdapter` behind one public preconditioner interface.

- [ ] **Step 1: Write convergence and reuse REDs.**

  Test periodic/Dirichlet Poisson manufactured solutions, uniform and 100:1 stretched axes, 1/2/4-rank residual reduction, unchanged hierarchy address for unchanged topology, numeric refresh without symbolic rebuild, and no HYPRE symbols when disabled.

- [ ] **Step 2: Implement full and semi-coarsening selection.**

  Select at plan freeze from metric aspect ratios. Full coarsening uses tensor restriction/prolongation; semi-coarsening preserves the strong direction and applies a line solve there.

- [ ] **Step 3: Implement quantified reuse.**

  Rebuild numeric smoothers only when coefficient revision changes; rebuild hierarchy only when topology/coarsening identity changes or the configured coefficient-change norm crosses the validated policy threshold. Record every rebuild counter.

- [ ] **Step 4: Verify and commit.**

  Run native MG tests with HYPRE off, then configure/build with HYPRE on when available. Commit `feat(v0.4): add native Cartesian multigrid`.

### Task 13: Build Momentum and Pressure Authorities with Exactly Two PISO Correctors

**Files:**
- Create: `versions/v0.4/include/hundun/v04_flow.hpp`
- Create: `versions/v0.4/src/operator_momentum.cpp`
- Create: `versions/v0.4/src/operator_pressure.cpp`
- Create: `versions/v0.4/src/integration_piso.cpp`
- Create: `versions/v0.4/tests/operator_momentum_test.cpp`
- Create: `versions/v0.4/tests/integration_piso_mpi_test.cpp`

**Interfaces:**
- Produces: `MomentumOperatorPlan`, `PressureOperatorPlan`, `MomentumIntermediates {rAU,rAtU,HbyA,phiHbyA}`, `PisoPlan {correctors=2}`, and final `U/face_flux` revisions.

- [ ] **Step 1: Write authority and mutation REDs.**

  Require exactly two pressure corrector calls, one momentum numeric assembly per step when coefficients are unchanged across correctors, correct invalidation of `rAU/rAtU/HbyA/phiHbyA`, one final face-flux writer, and failure when corrector 2 consumes corrector-1 stale intermediates.

- [ ] **Step 2: Implement momentum numeric update and intermediates.**

  Each intermediate stores its dependency tuple `(accepted_U, rho, mu_eff, dt, boundary, geometry)`. Reuse is permitted only on exact tuple equality. `rAtU` exists only for the selected consistent correction.

- [ ] **Step 3: Implement pressure correction and final publication.**

  Each corrector forms pressure RHS from `phiHbyA`, solves the declared pressure authority, corrects directional face flux, then updates velocity consistently. Only corrector 2 publishes final `U` and face-flux revisions.

- [ ] **Step 4: Verify and commit.**

  Run constant-field, Taylor-Green, conservation, rollback, and 1/2/4-rank two-corrector tests. Commit `feat(v0.4): add unified two-corrector PISO`.

### Task 14: Compile Static EB Geometry, Donor, and Surface Plans

**Files:**
- Create: `versions/v0.4/include/hundun/v04_eb.hpp`
- Create: `versions/v0.4/src/eb_geometry.cpp`
- Create: `versions/v0.4/src/eb_stencil.cpp`
- Create: `versions/v0.4/tests/eb_geometry_test.cpp`
- Create: `versions/v0.4/tests/eb_stencil_mpi_test.cpp`

**Interfaces:**
- Produces: `GeometryModel -> EbTopology -> BoundaryStencilPlan -> SurfaceQuadraturePlan`, compressed interface indices, donor weights, normals, areas, and geometry revision.

- [ ] **Step 1: Port scientific contracts as REDs, not implementation.**

  Recreate public mathematical fixtures for periodic cylinder ownership, surface windows, partition invariance, donor reproduction, normal orientation, and positive-normal donor feasibility. Include the existing ghost-donor failing configuration unchanged.

- [ ] **Step 2: Implement initialization-only geometry queries.**

  STL/BVH queries, cell classification, donor search, QR/weight construction, and quadrature happen before graph freeze. The hot plan stores compact indices and coefficients only.

- [ ] **Step 3: Repair positive-normal donor selection by root cause.**

  Preserve the required normal-side inequality and polynomial reproduction rank. Expand or reorder the deterministic candidate search only when geometry supports it; return `invalid_plan` when no scientifically valid donor set exists. Do not weaken the inequality or tolerance.

- [ ] **Step 4: Verify and commit.**

  Run geometry, mutation, periodic, 1/2/4-rank decomposition, and `test_immersed_wale_constant`-equivalent contracts. Commit `feat(v0.4): compile static EB stencil plans`.

### Task 15: Implement Exact IBM Pressure Correction and Final-State Force Authority

**Files:**
- Modify: `versions/v0.4/include/hundun/v04_eb.hpp`
- Create: `versions/v0.4/src/eb_pressure.cpp`
- Create: `versions/v0.4/src/eb_force.cpp`
- Create: `versions/v0.4/tests/eb_pressure_mpi_test.cpp`
- Create: `versions/v0.4/tests/eb_force_oracle_test.cpp`
- Create: `versions/v0.4/tests/eb_force_mutation_red_test.cpp`
- Create: `versions/v0.4/tests/integration_eb_transaction_mpi_test.cpp`

**Interfaces:**
- Consumes: native Cartesian pressure body, compact EB interface correction, final trial velocity, pressure, `mu_eff`, surface quadrature.
- Produces: exact nonsymmetric IBM pressure authority with FGMRES and compact preconditioner; operator force, budget reaction, pressure traction, viscous traction, and consistency report.

- [ ] **Step 1: Write the independent final-state force oracle RED.**

  The oracle reconstructs velocity gradients directly from the committed/final trial velocity and immutable stencil plan, then integrates traction. It must not call `eb_force.cpp` production helpers.

- [ ] **Step 2: Prove mutation sensitivity before implementation.**

  Compile a test-only mutation that feeds corrector scratch gradient into production force collection. Expected: oracle comparison fails while operator pressure and continuity remain unchanged. Restore the unmutated source before the green run.

- [ ] **Step 3: Implement regular-body plus compact-interface pressure apply.**

  Define homogeneous action as `F(p,0)-F(0,0)` and retain affine `F(0,g)` in the RHS. Use FGMRES outside; native Cartesian MG/compact correction is an inexact preconditioner only.

- [ ] **Step 4: Implement single force authority and transactional cache publication.**

  Final gradient depends on final trial `U` revision. Successful consensus publishes its cache for the next accepted-state consumer; failed attempts discard it. Four force fields share sign and unit metadata.

- [ ] **Step 5: Verify and commit.**

  Run oracle, mutation RED/green, force sign, pressure superposition, rollback, exactly-two-corrector, and 1/2/4-rank tests. Do not run 24³ yet. Commit `feat(v0.4): add exact IBM pressure and force authority`.

### Task 16: Add WALE and the Single mu_eff Authority

**Files:**
- Create: `versions/v0.4/include/hundun/v04_model.hpp`
- Create: `versions/v0.4/src/model_transport.cpp`
- Create: `versions/v0.4/src/model_wale.cpp`
- Create: `versions/v0.4/tests/model_wale_test.cpp`
- Create: `versions/v0.4/tests/model_wale_eb_mpi_test.cpp`

**Interfaces:**
- Produces: `DerivedFieldPlan` for shared velocity gradient, `WalePlan`, `TransportState`, and one `mu_eff` revision consumed by momentum diffusion and EB traction.

- [ ] **Step 1: Write WALE identity and reuse REDs.**

  Require zero eddy viscosity for uniform velocity and pure rotation identities, nonnegative finite values, correct Taylor-Green reference values, one gradient reconstruction per eligible revision, and identical `mu_eff` token at momentum and force consumers.

- [ ] **Step 2: Implement static model binding.**

  Select laminar/WALE once in `ValidatedModel`; store a direct kernel function in the execution graph. No model registry lookup or virtual call occurs per cell.

- [ ] **Step 3: Implement shared gradient cache.**

  Cache identity is `(U revision, geometry revision, boundary revision, EB stencil revision)`. A correct final-state cache may be published on commit and reused at next-step accepted-state stages.

- [ ] **Step 4: Verify and commit.**

  Run WALE unit, EB, 1/2/4-rank, rollback, allocation, and sanitizer tests. Commit `feat(v0.4): add WALE and unified transport authority`.

### Task 17: Add Conservative Density and Scalar Contributions

**Files:**
- Modify: `versions/v0.4/include/hundun/v04_model.hpp`
- Create: `versions/v0.4/src/model_contribution.cpp`
- Create: `versions/v0.4/tests/model_contribution_test.cpp`
- Create: `versions/v0.4/tests/integration_variable_density_mpi_test.cpp`

**Interfaces:**
- Produces: `ContributionPlan {equation, units, reads, writes, explicit_kernel, jacobian_kernel}`, conservative density/scalar transport, and static composition into existing equation authorities.

- [ ] **Step 1: Write units, conservation, and composition REDs.**

  Reject mismatched units, duplicate authorities, undeclared writes, and contribution cycles. Require periodic global mass/scalar conservation and 1/2/4-rank agreement.

- [ ] **Step 2: Implement cold-path contribution compilation.**

  Convert model contributions to ordered POD records and direct kernels. Explicit and Jacobian terms declare their stage and revisions; the equation remains the sole assembler.

- [ ] **Step 3: Implement constant and material-density profiles first.**

  Reuse final face-mass-flux authority for density and scalar transport. Do not duplicate convection or create a model-owned pressure/flux update path.

- [ ] **Step 4: Verify and commit.**

  Run conservation, composition, rollback, and MPI tests. Commit `feat(v0.4): add conservative model contributions`.

### Task 18: Add Restart, Diagnostics, Performance Evidence, and Input Migration

**Files:**
- Create: `versions/v0.4/include/hundun/v04_service.hpp`
- Create: `versions/v0.4/include/hundun/v04_app.hpp`
- Create: `versions/v0.4/src/service_restart.cpp`
- Create: `versions/v0.4/src/service_diagnostics.cpp`
- Create: `versions/v0.4/src/service_performance.cpp`
- Create: `versions/v0.4/src/app_driver.cpp`
- Modify: `versions/v0.4/src/app_main.cpp`
- Create: `versions/v0.4/tests/service_restart_mpi_test.cpp`
- Create: `versions/v0.4/tests/service_product_path_test.cpp`
- Create: `versions/v0.4/cases/cylinder_re3900/**`

**Interfaces:**
- Produces: committed-snapshot Restart, bitwise same-partition continuation, structured diagnostics, resource counters, v0.3 input migration, and tests-off product executable.

- [ ] **Step 1: Write Restart and product-path REDs.**

  Require schema/fingerprint validation before state mutation, bitwise continuous-vs-restart on the same partition, collective rejection of a corrupt rank payload, no observation of trial state, and identical numerical output from tests-on/tests-off binaries.

- [ ] **Step 2: Implement snapshot-only services.**

  Define `CommittedSnapshot` as immutable field views plus accepted layer and revision tokens. Services receive only this type and immutable plan identities. Restart reads into scratch storage, validates all ranks, then publishes through one transaction.

- [ ] **Step 3: Implement v0.3 case migration and Re=3900 case identity.**

  Migrate inputs and Restart metadata through explicit versioned readers; do not preserve v0.3 internal C++ APIs. Record canonical case/STL hashes and reject body-fitted/multiblock fields.

- [ ] **Step 4: Prove tests-off isolation.**

  Inspect link maps and symbols to ensure test oracles, full-domain gather helpers, and mutation seams are absent from the tests-off `hundun` binary.

- [ ] **Step 5: Verify and commit.**

  Run Restart 1/2/4-rank, corrupt-input, product-path, CLI validation, and offline configure/build tests. Commit `feat(v0.4): add committed-state services and migration`.

### Task 19: Freeze the Focused Gate and Performance Candidate Workflow

**Files:**
- Create: `versions/v0.4/tests/re3900_focused_manifest.cmake`
- Create: `versions/v0.4/cases/cylinder_re3900/case_24cube_20step.json`
- Create: `versions/v0.4/cases/cylinder_re3900/case_full_2step.json`
- Create: `versions/v0.4/docs/re3900-evidence-schema.md`
- Create: `cmake/HundunV04CandidateIdentity.cmake`
- Create: `docs/verification/v0.4-re3900-candidate-ledger.md`

**Interfaces:**
- Consumes: exact product candidate, v0.3 oracle, paired COAST commands/evidence, resource counters.
- Produces: immutable evidence manifests and a machine-checkable ordered gate state: `focused -> screen24 -> full2 -> frozen -> full20 -> statistics`.

- [ ] **Step 1: Write a gate-order mutation RED.**

  Attempt to register `full2`, `full20`, or `statistics` evidence without the preceding accepted state. Expected: manifest validator returns `REJECT` and leaves the ledger unchanged.

- [ ] **Step 2: Implement exact candidate identity.**

  Manifest fields include HEAD, tree, compiler, flags, tests-on/off binary hashes, case/STL hashes, MPI version, rank map, CPU/NUMA topology, environment allowlist, command, timestamps, exit status, and evidence-file SHA-256 values.

- [ ] **Step 3: Define focused and numerical equivalence checks.**

  Encode exact tolerances from the design: field normalized L2 `1e-10`; diagnostic relative `1e-10`/absolute `1e-12`; force absolute `1e-11` or relative `1e-9`; residual/conservation no more than `1.05x` reference; 20-step iterations no systematic increase beyond `2%`; bitwise Restart/protocol tests.

- [ ] **Step 4: Define performance comparison without running long jobs.**

  Record total/init/hot median/P90/RSS/communication/iterations. `HUNDUN median <= COAST median` is the release performance predicate; `<=1.10x` records only `NEAR`, never `ACCEPT`.

- [ ] **Step 5: Verify manifest tooling and commit.**

  Run only synthetic manifest tests and focused low-cost contracts. Do not run 24³ or full grid in this task. Commit `test(v0.4): freeze Re3900 joint-gate workflow`.

### Task 20: Execute the Re=3900 Numerical and Performance Joint Release Gate

**Files:**
- Modify: `docs/verification/v0.4-re3900-candidate-ledger.md`
- Create after execution: `docs/verification/v0.4-re3900-final-acceptance.md`
- External evidence only: `/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/**`

**Interfaces:**
- Consumes: a frozen exact v0.4 candidate and Task 19 gate tooling.
- Produces: one `ACCEPT` or `REJECT`; only `ACCEPT` authorizes a v0.4 release tag.

- [ ] **Step 1: Seal pre-run state and process inventory.**

  Record clean status, complete diff, HEAD/tree, binary/case hashes, existing task processes, evidence directory inventory, COAST pairing identity, and available memory. Do not stop unrelated processes.

- [ ] **Step 2: Run the focused gate.**

  Run Release focused tests, force oracle and mutation RED/green, WALE positive-normal donor contract, Restart/rollback, public headers, 1/2/4-rank MPI, ASan, and UBSan. Any failure returns `REJECT` before case runs.

- [ ] **Step 3: Run paired `24^3/20-step`.**

  Compare v0.4 with the frozen HUNDUN reference and paired COAST case. Check all numerical tolerances, two correctors per step, no instability, iteration trend, stage time, allocation count, and memory. Reject candidates without clear performance direction.

- [ ] **Step 4: Run paired full-grid 2-step.**

  Use `480x480x48`, 64 ranks, the frozen rank/NUMA mapping, `dt=0.006`, and equivalent physics. Record initialization and both steps, but do not infer final hot-step performance from step 2 alone.

- [ ] **Step 5: Freeze the candidate or return to implementation.**

  Freeze only if focused, 24³, and full2 pass numerical gates and show a credible path to `median <= COAST`. Once frozen, no source/config/build change is allowed without restarting at Step 1.

- [ ] **Step 6: Run paired full-grid 20-step.**

  Check numerical equivalence, total iterations, continuity/conservation, hot median/P90, peak RSS, communication, and `HUNDUN median <= COAST median`. Failure returns `REJECT`; it does not authorize threshold changes.

- [ ] **Step 7: Run fully developed statistics only for the frozen passing candidate.**

  Compare `Cd_mean`, `Cl_rms`, `St`, pressure distribution, and wake statistics to experiment and COAST. HUNDUN error must be no worse on all registered primary metrics and clearly better on at least one before claiming superior accuracy.

- [ ] **Step 8: Perform final provenance and complete-diff review.**

  Verify no copied upstream/COAST code, all DCO sign-offs, tests-off isolation, exact evidence hashes, no task MPI process, and final clean status. The report computes:

  ```text
  release = numerical_correctness_accept
         && robustness_accept
         && coast_performance_accept
         && physical_accuracy_accept
         && provenance_accept
  ```

- [ ] **Step 9: Record the single release decision.**

  If the conjunction is false, commit only a `REJECT` evidence receipt and return to the owning task. If true, commit `docs: accept v0.4 Re3900 joint release gate` with DCO; release/tag/push remains a separate user-authorized repository operation.

---

## Task Dependency and Review Order

Execute strictly in numerical dependency order:

```text
1 -> 2 -> 3 -> 4 -> 5 -> 6 -> 7 -> 8 -> 9 -> 10
                                      10 -> 11 -> 12 -> 13
5 + 7 + 9 + 13 -> 14 -> 15 -> 16 -> 17 -> 18 -> 19 -> 20
```

After each task, review requirements first, then code quality and performance risks, then inspect the complete task diff. A rejected review returns to the same task; it does not advance the dependency graph. Tasks 1-19 are internal checkpoints only. Task 20 is the sole release decision.

## Plan Self-Review

- **Spec coverage:** Tasks 2-10 cover version isolation, flat layout, schema, arena, Cartesian geometry, CPU/NUMA, communication, boundary/schemes, kernels, frozen graph, transactions, and resource budgets. Tasks 11-13 cover the four solver layers, native MG, unified equation authority, intermediate lifetimes, and exactly two PISO correctors. Tasks 14-17 cover static EB, the positive-normal donor failure, exact IBM pressure, independent force oracle/mutation RED, WALE, one `mu_eff`, and conservative multiphysics. Tasks 18-20 cover Restart/services, tests-off isolation, ordered evidence, COAST pairing, statistics, and the one release gate.
- **Scope:** GPU, AMR, body-fitted/multiblock, moving EB, fine-grained plugins, mixed precision, and case-specific tuning are explicitly excluded.
- **Type consistency:** `Status`, `RevisionToken`, `FieldId`, `FieldSchema`, `FieldStorage`, `CartesianGeometryPlan`, `MeshPatch`, `CpuExecutionPlan`, `CommunicationPlan`, `FrozenExecutionGraph`, the four linear lifecycle objects, `FaceFluxView`, and final `U/face_flux` revisions are defined once and consumed under the same names.
- **Lifecycle consistency:** geometry and symbolic plans invalidate only on topology revision; numeric coefficients invalidate on their declared tuple; hierarchy and workspace have separate identities; failed attempts never publish trial fields or caches.
- **Release consistency:** `1.10x COAST` is never an acceptance condition; only hot median no higher than COAST plus numerical, physical, robustness, and provenance acceptance can release v0.4.

# HUNDUN-FLOW v0.4 Cartesian Low-Mach Architecture Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an independently implemented, CPU-first HUNDUN-FLOW v0.4 single-phase low-Mach Cartesian CFD product and release it only after the Re=3900 cylinder numerical, robustness, full-grid short-performance, literature-accuracy, and provenance gates all accept.

**Architecture:** v0.4 compiles flat case input into immutable geometry, field, boundary, operator, communication, solver, and execution plans. Hot fields use padded x-fast SoA storage; one final face-mass-flux authority drives every conservative equation; static IBM data, halo metadata, linear hierarchies, and workspaces persist until their explicit revision identity changes. All capabilities register before one production freeze, after which exactly two PISO correctors advance a unified local-absolute-pressure EOS/enthalpy system without hot-path allocation or case-specific source branches.

**Tech Stack:** C++17, CMake 3.21+, MPI 3 C++ bindings, POSIX/Linux CPU and NUMA facilities, Threads, yyjson, optional HYPRE Struct adapter, CTest, AddressSanitizer, UndefinedBehaviorSanitizer, VTK XML output.

## Global Constraints

- Work only in `/home/wyf/code_dev/hundun-flow`; preserve `/home/wyf/code_dev/.worktrees/hundun-flow-cylinder3900-stage4`, its dirty/untracked files, and all `/home/wyf/code_dev/.benchmarks/cylinder-re3900*` evidence.
- Never reset, clean, overwrite, or commit changes in the existing cylinder worktree.
- Freeze v0.3 from commit `4ae4832ad00b5e4d1129ee978e2e49cbb33bb7ef`, tree `0a59ffa9d71cee79ef561c49657c9e90c5948e3e`, using tracked archive content only.
- Use accepted Stage 4 code identity `6407cd7c591ce088db7f1dd7e296d77acd18da1c`, tree `2791a1cee7ac8114f1696670d30c8951212d6024`, only as a scientific-contract reference.
- Do not copy or translate OpenFOAM GPL source or legacy COAST Fortran. Independently implement public mathematics, layout, and lifecycle ideas. The user-authorized `imb_mesh_y.cpp` STL scanning method may be ported or adapted without provenance boilerplate.
- v0.4 supports only uniform Cartesian and global tensor-product stretched Cartesian geometry. Reject body-fitted, multiblock curvilinear, AMR, and nonmatching refinement patches.
- v0.4 has one single-phase subsonic low-Mach compressible path. Do not add a constant-density fast path, shock path, reacting chemistry, PDF, spray, multiphase, or real-gas branch.
- EOS always consumes local absolute pressure `p_abs = p_ref + pi`; pressure correction includes `drho/dp`; sound speed, Mach, and NSCBC remain available even when acoustic CFL is not the default limit.
- v0.4 implements transient PISO only, with exactly two pressure correctors. Do not add SIMPLE, PIMPLE, strong PIMPLE, or a user-adjustable corrector count.
- Production turbulence choices are `vreman_wall_function` (default), `wale`, and `none` for verification. Do not add other production turbulence models in this version.
- Use static enthalpy as energy authority with full `Dp/Dt`, heat conduction, and viscous dissipation. Use N-1 inert species and optional passive scalars; do not persist `rhoU` or `rhoh` state arrays.
- `versions/v0.4/src` is flat and uses only `app_`, `core_`, `mesh_`, `parallel_`, `bc_`, `physics_`, `solver_`, and `io_` source prefixes.
- Keep tests in the source tree but no product `cases/` or `examples/` under `versions/v0.4`. Runtime examples are generated into a user-selected run directory.
- A production `FieldSchema`, arena layout, merged halo plan, linear lifecycles, and execution graph freeze exactly once after every capability registers and before the first product case.
- Hot-loop contracts: zero heap allocation, zero string dispatch, zero STL/BVH/donor search, zero implicit MPI, explicit revision-checked caches, one writer per authority, no full-field rollback copy.
- Preserve the standard IBM stencil contract: 14--32 fluid donors, all positive-normal, at least three normal bands, four tangential quadrants, reach at most four, deterministic QR, and no lower-order fallback.
- Do not accept the current final-gradient/surface-force candidate before an independent final-state force oracle and mutation-sensitive RED pass. Fix the existing WALE positive-normal donor failure by root cause, never by relaxing scientific constraints.
- Run validation in this order: focused contracts; full `480x480x48/64-rank` two-step HUNDUN/COAST short pairing; candidate freeze; at least five alternating full-grid 20-step pairings; HUNDUN-only literature statistics. Do not use `24^3` for performance and do not start long statistics before the short gate accepts.
- COAST is a short-performance baseline only because it has no periodic boundary support. HUNDUN literature statistics compare directly with experiments, not COAST long statistics.
- Every task ends with requirements review, code/performance review, complete-diff inspection, relevant tests, `git diff --check`, and a DCO-signed commit. A failed review returns to the same task.

## Stable Public Types

These names remain stable across tasks; later tasks add registrations and implementations without renaming them:

```cpp
namespace hundun::v04 {
struct Int3 { std::int32_t x{}, y{}, z{}; };
struct Real3 { double x{}, y{}, z{}; };
template<class T> struct Span { T* data{}; std::size_t size{}; };
using FieldId = std::uint16_t;
using StageId = std::uint16_t;
using RevisionToken = std::uint64_t;
using PlanFingerprint = std::uint64_t;
enum class StatusCode : std::uint16_t {
  ok, invalid_case, invalid_plan, allocation_failure, mpi_failure,
  numerical_failure, rejected_step, io_failure
};
struct Status {
  StatusCode code{StatusCode::ok};
  std::uint32_t detail{};
  constexpr explicit operator bool() const noexcept {
    return code == StatusCode::ok;
  }
};
}
```

## Target File Map

```text
versions/v0.4/
  CMakeLists.txt
  include/hundun/
    v04_types.hpp             stable POD identifiers and views
    v04_status.hpp            compact hot-path status
    v04_case.hpp              CaseSpec and ValidatedModel
    v04_field.hpp             registration, schema, arena views, revisions
    v04_mesh.hpp              Cartesian geometry, decomposition, STL scan
    v04_parallel.hpp          CPU/NUMA and persistent halo contracts
    v04_boundary.hpp          boundary, scheme, and time plans
    v04_physics.hpp           thermo, transport, contribution, turbulence
    v04_execution.hpp         stages, transactions, frozen graph, counters
    v04_linear.hpp            four-layer linear lifecycle and solvers
    v04_ibm.hpp               EB, reconstruction, pressure, and force
    v04_flow.hpp              equation plans and exactly-two-corrector PISO
    v04_io.hpp                Restart, Visit, screen, monitor, evidence
    v04_app.hpp               driver and CLI
  src/
    app_*.cpp core_*.cpp mesh_*.cpp parallel_*.cpp bc_*.cpp
    physics_*.cpp solver_*.cpp io_*.cpp
  tests/
    CMakeLists.txt
    unit/*.cpp mpi/*.cpp numerical/*.cpp integration/*.cpp
  docs/
    architecture.md numerics.md input-schema.md evidence-schema.md
```

The plan deliberately keeps public headers coarse and implementation files focused. Internal helpers remain in flat `src` as `<prefix>_<owner>_detail.hpp`; they are never installed.

---

### Task 1: Freeze the Public Lifecycle Survey and Adoption Ledger

**Files:**
- Create: `docs/references/2026-08-13-hundun-v04-public-lifecycle-survey.md`
- Create: `docs/references/2026-08-13-hundun-v04-adoption-ledger.tsv`
- Create: `docs/architecture/v0.4-target-hot-loop.md`
- Modify: `THIRD_PARTY.md`

**Interfaces:**
- Consumes: the fixed upstream revisions and COAST path in Global Constraints.
- Produces: one reviewable record mapping each referenced public idea to a HUNDUN owner, lifetime, invalidation rule, and independent-implementation decision.

- [ ] **Step 1: Record exact upstream identities and license boundaries.**

  Add one table row per fixed project with repository URL/path, revision, relevant files/symbols, license, allowed idea, prohibited copying, and HUNDUN destination. Record OpenFOAM as GPL-reference-only, AMReX/IncFlo/AMReX-Hydro as independently reimplemented BSD ideas, and COAST as user-specified read-only functional reference.

- [ ] **Step 2: Record the PISO intermediate lifecycle.**

  Document this exact contract:

  ```text
  momentum numeric revision -> rAU
  rAU + consistent diagonal revision -> rAtU
  momentum numeric revision + current trial U -> HbyA
  current HbyA + current trial U/phi + time/geometry/BC -> phiHbyA
  pressure equation flux -> only final face-mass-flux writer
  pressure gradient from the same solve -> final U update
  ```

  State explicitly that corrector 1 changes trial `U/phi`, so corrector 2 rebuilds or revision-recertifies `HbyA/phiHbyA`; `rAU/rAtU` alone may be reused when their coefficients are unchanged.

- [ ] **Step 3: Record geometry, halo, and linear lifetimes.**

  The survey must distinguish `EBTopology`, `BoundaryStencilPlan`, `SurfaceQuadraturePlan`, persistent halo metadata/buffers/requests, `SymbolicPlan`, `NumericState`, `HierarchyState`, and `SolverWorkspace`. Give each an identity tuple and its only legal rebuild causes.

- [ ] **Step 4: Trace COAST functionality without copying source.**

  Record the observed COAST capabilities required for replacement: local absolute pressure EOS, pressure density derivative, sound speed/Mach, NSCBC, SIMPLE-style transient iterations, Vreman wall function, flat `.d` inputs, STL scan, Restart, Visit, screen, and ICCG lifecycle. Mark SIMPLE and COAST nonperiodic limitations as references, not v0.4 algorithm requirements.

- [ ] **Step 5: Draw the target hot-loop text schedule.**

  Put this machine-reviewable stage order in `v0.4-target-hot-loop.md`:

  ```text
  begin attempt -> update thermo/transport -> momentum predictor
  -> PISO corrector 1 -> PISO corrector 2/final flux
  -> enthalpy -> inert species/passive scalars
  -> diagnostics snapshot decision -> collective commit/rollback
  ```

  Annotate halo begin/finish, cache publication, and collective consensus points.

- [ ] **Step 6: Verify the research artifact and commit.**

  Run:

  ```bash
  test "$(awk -F '\t' 'NR>1 {print $1}' docs/references/2026-08-13-hundun-v04-adoption-ledger.tsv | sort -u | wc -l)" -ge 5
  rg -n 'rAU|rAtU|HbyA|phiHbyA|EBTopology|SymbolicPlan|SolverWorkspace|NSCBC' docs/references/2026-08-13-hundun-v04-* docs/architecture/v0.4-target-hot-loop.md
  ! rg -n 'TB[D]|TO[D]O|FIXM[E]|XX[X]' docs/references/2026-08-13-hundun-v04-* docs/architecture/v0.4-target-hot-loop.md
  git diff --check
  ```

  Expected: all lifecycle terms appear, no placeholder pattern appears, and only documentation/`THIRD_PARTY.md` changed. Commit with `git commit -s -m "docs: freeze v0.4 lifecycle survey"`.

### Task 2: Isolate v0.3 and Establish a Buildable Flat v0.4 Skeleton

**Files:**
- Modify: `CMakeLists.txt`
- Modify: `CMakePresets.json`
- Modify: `README.md`
- Modify: `VERSION`
- Create: `cmake/HundunVersionDispatch.cmake`
- Create: `cmake/tests/version_dispatch_test.cmake`
- Create: `versions/v0.3/**` from tracked archive `4ae4832...` only
- Create: `versions/v0.3/FROZEN_SOURCE.json`
- Create: `versions/v0.4/CMakeLists.txt`
- Create: `versions/v0.4/include/hundun/v04_types.hpp`
- Create: `versions/v0.4/include/hundun/v04_status.hpp`
- Create: `versions/v0.4/src/core_status.cpp`
- Create: `versions/v0.4/src/app_main.cpp`
- Create: `versions/v0.4/tests/CMakeLists.txt`
- Create: `versions/v0.4/tests/unit/core_status_test.cpp`

**Interfaces:**
- Consumes: Task 1 survey; v0.3 source commit/tree fixed in Global Constraints.
- Produces: `HUNDUN_SOURCE_VERSION={v0.4,v0.3}` with v0.4 default, standalone version subtrees, target `hundun`, and the Stable Public Types.

- [ ] **Step 1: Write the version-dispatch configure RED.**

  In `cmake/tests/version_dispatch_test.cmake`, configure default, explicit `v0.4`, explicit `v0.3`, and invalid `v9`; assert the first three choose the intended subtree and `v9` fails with `unsupported HUNDUN_SOURCE_VERSION`.

  Run:

  ```bash
  cmake -P cmake/tests/version_dispatch_test.cmake
  ```

  Expected before implementation: failure because `HUNDUN_SOURCE_VERSION` is not recognized.

- [ ] **Step 2: Export v0.3 without reading dirty files.**

  Verify identities, create `versions/v0.3`, and use archive output only:

  ```bash
  test "$(git -C /home/wyf/code_dev/.worktrees/hundun-flow-cylinder3900-stage4 rev-parse HEAD)" = 4ae4832ad00b5e4d1129ee978e2e49cbb33bb7ef
  test "$(git -C /home/wyf/code_dev/.worktrees/hundun-flow-cylinder3900-stage4 rev-parse HEAD^{tree})" = 0a59ffa9d71cee79ef561c49657c9e90c5948e3e
  mkdir -p versions/v0.3
  git -C /home/wyf/code_dev/.worktrees/hundun-flow-cylinder3900-stage4 archive 4ae4832ad00b5e4d1129ee978e2e49cbb33bb7ef | tar -x -C versions/v0.3
  ```

  `FROZEN_SOURCE.json` contains source realpath, commit, tree, archive SHA-256, and import date; it contains no dirty-file inventory.

- [ ] **Step 3: Implement root version dispatch.**

  Set root `VERSION` to `0.4.0`, declare a cache string defaulting to `v0.4`, validate its two legal values, and `add_subdirectory(versions/${HUNDUN_SOURCE_VERSION})`. Keep each subtree's `project()` and dependency discovery local so v0.3 remains independently configurable.

- [ ] **Step 4: Add the minimal v0.4 public types and executable.**

  Implement the Stable Public Types exactly. `core_status.cpp` contains cold-path `status_message(Status)` mapping; `app_main.cpp --version` prints `HUNDUN-FLOW 0.4.0 source=v0.4`. All v0.4 sources obey the eight-prefix rule.

- [ ] **Step 5: Verify both source lines.**

  Run:

  ```bash
  cmake -P cmake/tests/version_dispatch_test.cmake
  cmake -S . -B build/v04-debug -DHUNDUN_SOURCE_VERSION=v0.4 -DHUNDUN_BUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Debug
  cmake --build build/v04-debug -j2
  ctest --test-dir build/v04-debug -R '^v04_core_status$' --output-on-failure
  cmake -S . -B build/v03-release -DHUNDUN_SOURCE_VERSION=v0.3 -DHUNDUN_BUILD_TESTS=OFF -DCMAKE_BUILD_TYPE=Release
  cmake --build build/v03-release -j2 --target hundun
  build/v04-debug/versions/v0.4/hundun --version
  ```

  Expected: dispatch tests and status test pass, both executables build, v0.4 reports `0.4.0`, and the cylinder worktree status is byte-for-byte unchanged from a pre-task receipt.

- [ ] **Step 6: Review and commit.**

  Confirm `git diff --check`, no v0.4 product source outside `versions/v0.4`, no archived `.git` directory, and tests-off v0.4 registers zero tests. Commit with `git commit -s -m "build: isolate v0.3 and v0.4 source lines"`.

### Task 3: Compile Flat Case Input into a Validated Capability Model

**Files:**
- Create: `versions/v0.4/include/hundun/v04_case.hpp`
- Create: `versions/v0.4/include/hundun/v04_field.hpp`
- Create: `versions/v0.4/src/app_case.cpp`
- Create: `versions/v0.4/src/app_case_detail.hpp`
- Create: `versions/v0.4/src/core_field_registry.cpp`
- Create: `versions/v0.4/tests/unit/app_case_test.cpp`
- Create: `versions/v0.4/tests/mpi/app_case_broadcast_test.cpp`
- Create: `versions/v0.4/tests/data/case_minimal_valid.json`
- Create: `versions/v0.4/docs/input-schema.md`
- Modify: `versions/v0.4/CMakeLists.txt`
- Modify: `versions/v0.4/tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `Status`, `Int3`, `Real3` from Task 2 and bundled yyjson.
- Produces: `CaseSpec`, `ValidatedModel`, `CaseCompiler::load_and_compile(MPI_Comm, const std::filesystem::path&)`, `FieldRegistry::declare_field`, and a reusable `FieldRegistry::freeze()` mechanism; it does not freeze the production schema.

```cpp
enum class GeometryKind : std::uint8_t { uniform, tensor_stretched };
enum class TurbulenceKind : std::uint8_t { none, wale, vreman_wall_function };
enum class TimeControlKind : std::uint8_t { fixed, adaptive_flow, adaptive_acoustic };
struct CaseSpec { std::filesystem::path root; /* parsed cold data */ };
struct ValidatedModel {
  GeometryKind geometry;
  TurbulenceKind turbulence;
  TimeControlKind time_control;
  std::vector<std::filesystem::path> data_files;
  std::optional<std::filesystem::path> stl_file;
  PlanFingerprint fingerprint;
};
class FieldRegistry {
 public:
  Status declare_field(std::string_view stable_name, std::uint8_t components,
                       std::uint8_t ghost_width, FieldId& out);
  Status freeze(FieldSchema& out);
};
```

- [ ] **Step 1: Write invalid-input REDs.**

  Test missing `case.json`, duplicate JSON keys, unknown units, paths outside the case root, `.d` files hidden in subdirectories, STL hidden in a subdirectory, body-fitted/AMR selection, constant-density path selection, SIMPLE/PIMPLE, non-two PISO count, reacting physics, unsupported turbulence, missing pressure closure, and duplicate field names.

- [ ] **Step 2: Run the REDs.**

  Run:

  ```bash
  cmake --build build/v04-debug -j2 --target v04_app_case_test
  ctest --test-dir build/v04-debug -R '^v04_app_case$' --output-on-failure
  ```

  Expected before implementation: compile failure because `CaseCompiler` and `FieldRegistry` do not exist.

- [ ] **Step 3: Implement rank-zero parsing and canonical compilation.**

  Rank 0 resolves `case.json`, direct-root `.d` and STL references, units, enums, pressure closure, and feature compatibility. Serialize the validated typed representation into a bounded byte buffer; broadcast length, bytes, status category, and lowest failing rank. Non-root ranks never parse JSON or touch case files.

- [ ] **Step 4: Implement registration without production freeze.**

  `declare_field` assigns deterministic stable IDs by registration order, rejects duplicates and post-freeze mutation, and stores component/ghost requirements. `freeze` publishes the immutable schema; product code does not call it before Task 18's single bundle freeze.

- [ ] **Step 5: Verify MPI identity and root-only I/O.**

  Register `v04_app_case_broadcast_{1,2,4}`. Interpose a test file-open counter and assert only rank 0 opens JSON/`.d`/STL metadata. Assert every rank receives identical `ValidatedModel::fingerprint` and identical failure category.

  Run:

  ```bash
  ctest --test-dir build/v04-debug -R '^v04_app_case(_broadcast_[124])?$' --output-on-failure
  ```

- [ ] **Step 6: Review and commit.**

  Verify the input schema explicitly states `case.json` authority and flat case-root files, no case-specific field names exist in product code, and production schema remains unfrozen. Commit with `git commit -s -m "feat(v0.4): compile flat case input"`.

### Task 4: Implement Padded SoA Storage, Revisions, and Attempt Transactions

**Files:**
- Modify: `versions/v0.4/include/hundun/v04_field.hpp`
- Create: `versions/v0.4/include/hundun/v04_execution.hpp`
- Create: `versions/v0.4/src/core_arena.cpp`
- Create: `versions/v0.4/src/core_arena_detail.hpp`
- Create: `versions/v0.4/src/core_revision.cpp`
- Create: `versions/v0.4/src/core_transaction.cpp`
- Create: `versions/v0.4/tests/unit/core_arena_test.cpp`
- Create: `versions/v0.4/tests/unit/core_transaction_test.cpp`
- Create: `versions/v0.4/tests/mpi/core_transaction_mpi_test.cpp`
- Modify: `versions/v0.4/CMakeLists.txt`
- Modify: `versions/v0.4/tests/CMakeLists.txt`

**Interfaces:**
- Consumes: synthetic `FieldSchema` from Task 3.
- Produces: `ArenaLayout`, `FieldStorage`, `FieldView<T>`, `StateLayers`, `RevisionSet`, `AttemptTransaction`, and allocation counters consumed by every later task.

```cpp
struct FieldView {
  double* base{};
  Int3 interior{};
  Int3 ghosts{};
  std::size_t stride_y{}, stride_z{}, component_stride{};
  RevisionToken revision{};
};
class AttemptTransaction {
 public:
  Status begin(const StateLayers& accepted, StateLayers& trial) noexcept;
  Status publish_pending_cache(FieldId, RevisionToken) noexcept;
  Status collective_finish(MPI_Comm, Status local) noexcept;
  bool committed() const noexcept;
};
```

- [ ] **Step 1: Write alignment, layout, and allocation REDs.**

  Assert base alignment is 64 bytes, x stride is one, padded y stride is divisible by eight doubles, components are separate SoA spans, interior and ghosts share one allocation, checked views reject stale revisions, and no allocation occurs while an allocation guard surrounds `begin -> work -> commit/rollback`.

- [ ] **Step 2: Write transaction REDs.**

  Test successful handle rotation, retry without accepted-state mutation, pending-cache discard on failure, no whole-field copy counter, and 1/2/4-rank consensus choosing the same status and lowest failing rank.

- [ ] **Step 3: Implement deterministic arena planning.**

  Compute every field offset from schema, padded extents, alignment, NUMA node, and lifetime class before allocating one aligned arena per NUMA owner. Build immutable views; do not expose owning vectors in operator interfaces.

- [ ] **Step 4: Implement layer rotation and compact rollback.**

  Keep accepted `n`, accepted `n-1`, trial, and scratch handles. `begin` binds trial destinations without copying accepted arrays; commit rotates handles and increments revisions; rollback discards trial and pending cache metadata. Small controller values use a bounded transaction log.

- [ ] **Step 5: Verify locally and under MPI.**

  Run:

  ```bash
  ctest --test-dir build/v04-debug -R '^v04_core_(arena|transaction|transaction_mpi_[124])$' --output-on-failure
  ```

  Expected: all tests pass; allocation count inside guarded attempts is zero and accepted CRCs remain unchanged after injected failures.

- [ ] **Step 6: Review and commit.**

  Inspect for `new`, `malloc`, `std::vector::resize`, or string construction reachable from guarded hot methods. Commit with `git commit -s -m "feat(v0.4): add revisioned SoA state arenas"`.

### Task 5: Build CartesianGeometryPlan and the Authorized STL Scan

**Files:**
- Create: `versions/v0.4/include/hundun/v04_mesh.hpp`
- Create: `versions/v0.4/src/mesh_cartesian.cpp`
- Create: `versions/v0.4/src/mesh_focus.cpp`
- Create: `versions/v0.4/src/mesh_decomposition.cpp`
- Create: `versions/v0.4/src/mesh_stl_scan.cpp`
- Create: `versions/v0.4/src/mesh_stl_scan_detail.hpp`
- Create: `versions/v0.4/tests/unit/mesh_cartesian_test.cpp`
- Create: `versions/v0.4/tests/unit/mesh_focus_test.cpp`
- Create: `versions/v0.4/tests/unit/mesh_stl_scan_test.cpp`
- Create: `versions/v0.4/tests/mpi/mesh_decomposition_mpi_test.cpp`
- Create: `versions/v0.4/tests/data/cube_binary.stl`
- Create: `versions/v0.4/tests/data/cylinder_ascii.stl`
- Modify: `versions/v0.4/CMakeLists.txt`
- Modify: `versions/v0.4/tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `ValidatedModel`, `ArenaLayout`, and stable POD types.
- Produces: `CartesianGeometryPlan`, `MeshPatch`, `CpuTile`, `TriangleSoA`, and `StlScanPlan`.

```cpp
struct CartesianGeometryPlan {
  GeometryKind kind;
  Span<const double> x_faces, y_faces, z_faces;
  Span<const double> dx, dy, dz;
  Int3 global_cells;
  RevisionToken topology_revision;
  PlanFingerprint fingerprint;
};
struct MeshPatch { Int3 begin, cells, process_grid, process_coord; };
struct StlScanPlan {
  Status classify(const CartesianGeometryPlan&, Span<std::uint8_t> region) const noexcept;
};
```

- [ ] **Step 1: Write geometry and focus-region REDs.**

  Test exact uniform faces, monotone stretched faces, focus target attainment, growth limits, max-cell/memory rejection, deterministic coordinates, `exact_cells`, and explicit rejection of body-fitted/AMR/nonmatching input.

- [ ] **Step 2: Write STL algorithm-consistency REDs.**

  Translate the user-authorized method's mathematical cases into tests: axis Möller--Trumbore hit, transverse bounding rejection, inclusive edge/vertex barycentrics, same-position/same-normal duplicate removal, opposite-normal retention, sorted parity, reversed triangle order, ASCII/binary identity, and scan-axis permutations.

- [ ] **Step 3: Run the REDs.**

  Run:

  ```bash
  cmake --build build/v04-debug -j2 --target v04_mesh_cartesian_test v04_mesh_stl_scan_test
  ctest --test-dir build/v04-debug -R '^v04_mesh_' --output-on-failure
  ```

  Expected before implementation: missing `CartesianGeometryPlan` and `StlScanPlan` compile errors.

- [ ] **Step 4: Implement uniform and tensor-stretched plans.**

  Generate global one-dimensional face arrays on rank 0, validate them, broadcast typed arrays, and compute metrics once. `focus_regions` merge into target spacing envelopes subject to hard min-spacing/growth/cell/memory limits. No per-cell geometry polymorphism is permitted.

- [ ] **Step 5: Implement batched STL scan.**

  Read and canonicalize triangles once into SoA, bin triangles by transverse extents, assign scan lines to fixed thread chunks, reuse bounded per-thread intersection scratch, sort intersections deterministically, apply duplicate rules and parity, and produce compact region flags. Do not use per-cell nearest-triangle comparisons.

- [ ] **Step 6: Implement deterministic decomposition and tiles.**

  Choose one contiguous `MeshPatch` per rank from a process grid; allow unequal patch sizes while preserving exact global coverage. Generate rank-local CPU tiles separately from ownership. Verify 1/2/4 ranks give the same global geometry and classification fingerprints.

- [ ] **Step 7: Verify and commit.**

  Run all `^v04_mesh_` tests and a Release scan timing that asserts a single STL load, zero hot classification allocations, and deterministic CRC. Commit with `git commit -s -m "feat(v0.4): add Cartesian geometry and STL scan plans"`.

### Task 6: Implement CPU/NUMA Placement and the Persistent Halo Engine

**Files:**
- Create: `versions/v0.4/include/hundun/v04_parallel.hpp`
- Create: `versions/v0.4/src/parallel_cpu.cpp`
- Create: `versions/v0.4/src/parallel_halo.cpp`
- Create: `versions/v0.4/src/parallel_halo_detail.hpp`
- Create: `versions/v0.4/tests/unit/parallel_cpu_test.cpp`
- Create: `versions/v0.4/tests/mpi/parallel_halo_mpi_test.cpp`
- Create: `versions/v0.4/tests/mpi/parallel_halo_failure_test.cpp`
- Modify: `versions/v0.4/CMakeLists.txt`
- Modify: `versions/v0.4/tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `MeshPatch`, `FieldView`, `RevisionToken`, and allocation guard.
- Produces: `CpuExecutionPlan`, `HaloEngine`, `HaloFieldSpec`, `HaloTicket`; Task 18 later compiles the production `CommunicationPlan`.

```cpp
struct HaloFieldSpec { FieldId field; std::uint8_t width; std::uint8_t components; };
class HaloEngine {
 public:
  Status reserve(MPI_Comm, const MeshPatch&, Span<const HaloFieldSpec>) noexcept;
  Status begin(StageId, Span<const FieldView>, HaloTicket&) noexcept;
  Status finish(HaloTicket&, Span<FieldView>) noexcept;
};
```

- [ ] **Step 1: Write placement and halo REDs.**

  Test default one-rank-per-NUMA recommendation, fixed-team creation, pure-MPI fallback, duplicate-core warning, direct-neighbor-only messages, periodic and nonperiodic peers, merged same-peer payloads, persistent buffer/request addresses, ghost revision publication only after `finish`, and zero allocation during repeated exchange.

- [ ] **Step 2: Inject MPI failures.**

  At 1/2/4 ranks inject pack, `MPI_Startall`, completion, and unpack failures. Assert all ranks return one category/lowest rank, no ghost certification is published, and accepted state is unchanged.

- [ ] **Step 3: Implement cold CPU discovery and finite specialization.**

  Query affinity and NUMA placement only during startup, validate requested rank/thread layout, and bind one persistent thread team. Store a finite ISA/tile selection; never dispatch on CPUID inside kernels.

- [ ] **Step 4: Implement reusable halo transport.**

  `reserve` calculates maximum merged peer payloads and constructs persistent requests/buffers. `begin` packs current revisions and starts requests; callers compute interior tiles; `finish` waits, unpacks, and publishes exact ghost revision certificates. The engine has no knowledge of equations or implicit fill-patch behavior.

- [ ] **Step 5: Verify and commit.**

  Run:

  ```bash
  ctest --test-dir build/v04-debug -R '^v04_parallel_(cpu|halo.*_[124])$' --output-on-failure
  ```

  Inspect message/byte counters against compiled bounds and commit with `git commit -s -m "feat(v0.4): add persistent CPU halo engine"`.

### Task 7: Compile Boundary, Scheme, and Time-Control Plans

**Files:**
- Create: `versions/v0.4/include/hundun/v04_boundary.hpp`
- Create: `versions/v0.4/src/bc_compile.cpp`
- Create: `versions/v0.4/src/bc_apply.cpp`
- Create: `versions/v0.4/src/bc_nscbc.cpp`
- Create: `versions/v0.4/src/bc_time.cpp`
- Create: `versions/v0.4/tests/unit/bc_compile_test.cpp`
- Create: `versions/v0.4/tests/unit/bc_nscbc_test.cpp`
- Create: `versions/v0.4/tests/unit/bc_time_test.cpp`
- Create: `versions/v0.4/tests/mpi/bc_decomposition_mpi_test.cpp`
- Modify: `versions/v0.4/CMakeLists.txt`
- Modify: `versions/v0.4/tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `ValidatedModel`, `CartesianGeometryPlan`, `FieldRegistry` registration API.
- Produces: `BoundaryPlan`, `SchemePlan`, `TimeSchemePlan`, `TimeControllerState`.

```cpp
enum class BoundaryKind : std::uint8_t {
  velocity_inlet, mass_flow_inlet, static_state_inlet, total_state_inlet,
  pressure_outlet, nscbc_inlet, nscbc_outlet, no_slip_wall, moving_wall,
  slip, symmetry, periodic, adiabatic_wall, isothermal_wall, heat_flux_wall
};
enum class ConvectionScheme : std::uint8_t { central2, limited_central2, tvd2 };
enum class TimeScheme : std::uint8_t { backward_euler, variable_bdf2 };
```

- [ ] **Step 1: Write closure and compatibility REDs.**

  Cover every boundary kind, scalar/species values or fluxes, backflow state, illegal periodic pairing, missing absolute-pressure closure, conflicting thermal authorities, unsupported shocks/supersonic settings, and unsupported scheme names. Assert lookup strings disappear after compilation.

- [ ] **Step 2: Write NSCBC analytic REDs.**

  Use uniform subsonic inflow/outflow states and small characteristic perturbations. Check incoming/outgoing characteristic selection, local sound speed use, pressure/temperature closure, sign conventions, and finite behavior as Mach approaches zero without dividing by velocity.

- [ ] **Step 3: Write time-control REDs.**

  Check fixed, `adaptive_flow`, and `adaptive_acoustic`; convective/viscous/thermal/species minima; growth and retry limits; BE startup/restart/retry; variable-step BDF2 coefficients; and default absence of an acoustic hard limit in `adaptive_flow`.

- [ ] **Step 4: Implement cold compilation and static dispatch.**

  Build compact patch-index spans and function-free enums/parameters for each stage. Compile only `central2`, `limited_central2`, `tvd2`, and second-order diffusion. Register all required fields/ghost widths but do not freeze the production schema.

- [ ] **Step 5: Verify decomposition and commit.**

  Run all `^v04_bc_` tests at 1/2/4 ranks, inspect the hot apply path for strings or virtual calls, and commit with `git commit -s -m "feat(v0.4): compile boundary and time plans"`.

### Task 8: Add Unified Thermodynamics, Transport, and Contribution Contracts

**Files:**
- Create: `versions/v0.4/include/hundun/v04_physics.hpp`
- Create: `versions/v0.4/src/physics_thermo.cpp`
- Create: `versions/v0.4/src/physics_transport.cpp`
- Create: `versions/v0.4/src/physics_contribution.cpp`
- Create: `versions/v0.4/src/physics_derived.cpp`
- Create: `versions/v0.4/tests/unit/physics_thermo_test.cpp`
- Create: `versions/v0.4/tests/unit/physics_transport_test.cpp`
- Create: `versions/v0.4/tests/unit/physics_contribution_test.cpp`
- Create: `versions/v0.4/tests/numerical/physics_closed_mass_test.cpp`
- Modify: `versions/v0.4/CMakeLists.txt`
- Modify: `versions/v0.4/tests/CMakeLists.txt`

**Interfaces:**
- Consumes: field registration, geometry metrics, boundary plan.
- Produces: `ThermodynamicsPlan`, `TransportPlan`, `DerivedFieldPlan`, `ContributionRegistry`, single base `mu_eff` authority.

```cpp
struct ThermoState { double rho, temperature, cp, drho_dp_hY, sound_speed, mach; };
class ThermodynamicsPlan {
 public:
  Status evaluate(double p_abs, double h, Span<const double> independent_Y,
                  Real3 velocity, ThermoState&) const noexcept;
};
struct ContributionSpec {
  FieldId conserved_quantity;
  StageId stage;
  Span<const FieldId> reads;
  bool supplies_implicit_diagonal;
};
```

- [ ] **Step 1: Write EOS and enthalpy inversion REDs.**

  Test temperature-dependent ideal-gas mixtures, N-1 mass fractions, sum/positivity rejection, `p_abs=p_ref+pi`, finite-difference verification of `drho_dp_hY`, sound speed/Mach, `h(T,Y)` inversion, and constant-property specialization matching the generic path.

- [ ] **Step 2: Write transport and derived-cache REDs.**

  Verify `mu(T,Y)` and conductivity, one velocity-gradient computation per matching revision, one `mu_eff` writer, precise invalidation by `U/geometry/boundary/turbulence` revisions, and rejection of stale cached gradients.

- [ ] **Step 3: Write contribution registration REDs.**

  Require units/conserved quantity, read set, stage, explicit source, optional implicit diagonal, and deterministic order. Reject two writers, undeclared reads, chemistry/reacting contributions, or registration after production freeze.

- [ ] **Step 4: Implement one thermodynamic path.**

  Implement mixture molecular weight, NASA-polynomial-style tabulated `cp/h` evaluation from case `.d` data, safeguarded temperature inversion, local absolute-pressure EOS derivative, transport, and derived-state revision publication. Do not add a `DensityKind::constant` branch or persistent conservative products.

- [ ] **Step 5: Implement closed-domain mass Newton.**

  Given committed `h,Y,pi` and target total mass, update `p_ref` using global mass residual and summed `drho/dp`, with bounded Newton steps and collective convergence/failure. Open domains leave `p_ref` under pressure-boundary authority.

- [ ] **Step 6: Verify and commit.**

  Run all `^v04_physics_(thermo|transport|contribution|closed_mass)` tests, scan public enums for a constant-density path, and commit with `git commit -s -m "feat(v0.4): add low-Mach thermo and transport authorities"`.

### Task 9: Implement Conservative Cartesian Kernels and One Face-Flux Authority

**Files:**
- Create: `versions/v0.4/src/solver_cartesian.cpp`
- Create: `versions/v0.4/src/solver_cartesian_detail.hpp`
- Create: `versions/v0.4/src/solver_face_flux.cpp`
- Create: `versions/v0.4/src/solver_diffusion.cpp`
- Create: `versions/v0.4/tests/numerical/solver_manufactured_cartesian_test.cpp`
- Create: `versions/v0.4/tests/numerical/solver_conservation_test.cpp`
- Create: `versions/v0.4/tests/mpi/solver_decomposition_mpi_test.cpp`
- Modify: `versions/v0.4/include/hundun/v04_execution.hpp`
- Modify: `versions/v0.4/CMakeLists.txt`
- Modify: `versions/v0.4/tests/CMakeLists.txt`

**Interfaces:**
- Consumes: `FieldView`, `CartesianGeometryPlan`, `BoundaryPlan`, `SchemePlan`, `ContributionRegistry`, halo engine.
- Produces: `FaceFluxView`, fused reconstruction/divergence kernels, gradient and diffusion kernels, explicit single-writer token for final mass flux.

```cpp
struct FaceFluxView { FieldView x, y, z; RevisionToken revision; };
struct KernelInvocation { Span<const FieldView> reads; Span<FieldView> writes; };
Status reconstruct_mass_flux(const SchemePlan&, const CartesianGeometryPlan&,
                             const KernelInvocation&, FaceFluxView&) noexcept;
```

- [ ] **Step 1: Write manufactured-order REDs.**

  On uniform and tensor-stretched grids, run 16/32/64 sequences for gradient, divergence, convection, and diffusion. Require both adjacent observed orders at least 1.8 above the roundoff floor.

- [ ] **Step 2: Write conservation and authority REDs.**

  Verify shared interior face values cancel bitwise between adjacent cells, periodic global divergence is roundoff zero, boundary flux equals volume change, density/momentum/enthalpy/scalars consume the same `FaceFluxView::revision`, and a second final-flux writer is rejected.

- [ ] **Step 3: Implement x-contiguous static kernels.**

  Use uniform/stretched plan specializations selected before the loop. Fuse face reconstruction, interpolation, provisional mass-flux formation, and divergence where lifetime analysis permits. Separate x/y/z flux arrays and avoid gather/scatter temporaries over full fields.

- [ ] **Step 4: Integrate explicit halo sequencing.**

  Operators request a halo ticket from the driver, compute interior tiles, finish exchange, then compute boundary tiles. Kernels themselves cannot call MPI. Publish output revisions only after both tile sets succeed.

- [ ] **Step 5: Verify performance contracts and commit.**

  Run numerical and 1/2/4-rank decomposition tests, allocation guards for 100 repeated invocations, and byte counters. Commit with `git commit -s -m "feat(v0.4): add conservative Cartesian kernels"`.

### Task 10: Implement the Stage Graph Compiler without Freezing a Product Case

**Files:**
- Modify: `versions/v0.4/include/hundun/v04_execution.hpp`
- Create: `versions/v0.4/src/core_execution_graph.cpp`
- Create: `versions/v0.4/src/core_resource_contract.cpp`
- Create: `versions/v0.4/tests/unit/core_execution_graph_test.cpp`
- Create: `versions/v0.4/tests/unit/core_resource_contract_test.cpp`
- Modify: `versions/v0.4/CMakeLists.txt`
- Modify: `versions/v0.4/tests/CMakeLists.txt`

**Interfaces:**
- Consumes: field/stage registrations and halo specifications from Tasks 3--9.
- Produces: `StageSpec`, `ExecutionGraphCompiler`, `FrozenExecutionGraph`, `ResourceCounters`; Task 18 performs the only production freeze.

```cpp
struct StageSpec {
  StageId id;
  Span<const FieldId> reads, writes, ghosts;
  Span<const std::uint8_t> ghost_widths;
  Span<const FieldId> invalidates;
  std::size_t workspace_bytes;
  bool collective_consensus;
};
class ExecutionGraphCompiler {
 public:
  Status register_stage(const StageSpec&);
  Status freeze_for_test(FrozenExecutionGraph&) const;
};
```

- [ ] **Step 1: Write graph REDs.**

  Reject read-before-produce, two writers, missing ghost production, undeclared invalidation, workspace overlap while live, a commit before collective consensus, a service reading trial state, and registration after freeze.

- [ ] **Step 2: Write resource-budget REDs.**

  Compile synthetic stages and assert exact max live workspace, allocation allowance zero, merged peer message/byte ceilings, allowed numeric refill/hierarchy rebuild counts, and counter overflow rejection.

- [ ] **Step 3: Implement deterministic topological compilation.**

  Preserve declared numerical order, add dependency edges, compute liveness/alias intervals, assign preallocated scratch offsets, attach explicit halo begin/finish nodes, and insert only declared collective consensus points. Fingerprint every compiled property.

- [ ] **Step 4: Prove mechanism-only scope.**

  Tests may call `freeze_for_test` on synthetic complete registrations. No product driver, cylinder case, production `FieldSchema`, or production `CommunicationPlan` may be created in this task.

- [ ] **Step 5: Verify and commit.**

  Run `^v04_core_(execution_graph|resource_contract)$`, inspect graph dumps for determinism under registration replay, and commit with `git commit -s -m "feat(v0.4): add execution graph compiler"`.

### Task 11: Implement the Four-Layer Linear Lifecycle and Krylov Solvers

**Files:**
- Create: `versions/v0.4/include/hundun/v04_linear.hpp`
- Create: `versions/v0.4/src/solver_linear_state.cpp`
- Create: `versions/v0.4/src/solver_linear_state_detail.hpp`
- Create: `versions/v0.4/src/solver_krylov.cpp`
- Create: `versions/v0.4/src/solver_reduction.cpp`
- Create: `versions/v0.4/tests/unit/solver_linear_lifecycle_test.cpp`
- Create: `versions/v0.4/tests/mpi/solver_krylov_mpi_test.cpp`
- Create: `versions/v0.4/tests/mpi/solver_reduction_mpi_test.cpp`
- Modify: `versions/v0.4/CMakeLists.txt`
- Modify: `versions/v0.4/tests/CMakeLists.txt`

**Interfaces:**
- Consumes: geometry/topology revision, boundary locations, explicit halo engine, persistent arena workspace, resource counters.
- Produces: `SymbolicPlan`, `NumericState`, `HierarchyState`, `SolverWorkspace`, `LinearOperator`, PCG, right-preconditioned restarted FGMRES, BiCGSTAB, FP64 true-residual checks.

```cpp
struct LinearIdentity {
  RevisionToken topology, boundary_layout, coefficients, hierarchy;
  PlanFingerprint fingerprint;
};
class LinearOperator {
 public:
  virtual Status apply(const FieldView& x, FieldView& y) const noexcept = 0;
  virtual ~LinearOperator() = default; // one call per field, never per cell
};
struct LinearSolveResult {
  Status status;
  std::uint32_t iterations;
  double initial_true_residual, final_true_residual;
};
```

- [ ] **Step 1: Write lifecycle REDs.**

  Construct synthetic topology/coefficient revisions and assert: topology change rebuilds all four layers; coefficient change refills numeric state but preserves symbolic state/workspace; hierarchy policy alone decides coarse rebuild; repeated solve reuses stable workspace addresses; stale numeric/hierarchy identities are rejected rather than silently used.

- [ ] **Step 2: Write solver REDs.**

  Use SPD Poisson, nonsymmetric advection-diffusion, a breakdown matrix, a zero RHS, and an injected MPI rank failure. Require PCG only for certified SPD operators, FGMRES/BiCGSTAB for nonsymmetric systems, deterministic failure classification, finite iterates, and explicit true-residual acceptance.

- [ ] **Step 3: Implement persistent vector/reduction workspace.**

  Register the maximum Krylov basis and reduction buffers before schema freeze. Reuse them for every solve. Implement reproducible fixed-tree reduction and performance `MPI_Allreduce` modes selected in the cold plan; both report an FP64 true residual.

- [ ] **Step 4: Implement PCG, FGMRES, and BiCGSTAB.**

  Solvers receive non-owning views, operator, preconditioner, tolerances, iteration cap, and workspace. They perform no allocation, no string formatting, and no implicit operator rebuild. FGMRES is restarted and right-preconditioned so IBM can vary its preconditioner while keeping an exact outer operator.

- [ ] **Step 5: Verify at 1/2/4 ranks and commit.**

  Run:

  ```bash
  ctest --test-dir build/v04-debug -R '^v04_solver_(linear_lifecycle|krylov_mpi_[124]|reduction_mpi_[124])$' --output-on-failure
  ```

  Repeat every solve 100 times under the allocation guard, inspect workspace addresses and iteration counters, then commit with `git commit -s -m "feat(v0.4): add persistent Krylov lifecycle"`.

### Task 12: Implement NativeCartesianMG and an Isolated HYPRE Struct Adapter

**Files:**
- Modify: `versions/v0.4/include/hundun/v04_linear.hpp`
- Create: `versions/v0.4/src/solver_mg.cpp`
- Create: `versions/v0.4/src/solver_mg_detail.hpp`
- Create: `versions/v0.4/src/solver_hypre.cpp`
- Create: `versions/v0.4/tests/numerical/solver_mg_convergence_test.cpp`
- Create: `versions/v0.4/tests/unit/solver_mg_reuse_test.cpp`
- Create: `versions/v0.4/tests/mpi/solver_mg_mpi_test.cpp`
- Create: `versions/v0.4/tests/integration/solver_hypre_isolation_test.cpp`
- Modify: `versions/v0.4/CMakeLists.txt`
- Modify: `versions/v0.4/tests/CMakeLists.txt`

**Interfaces:**
- Consumes: four-layer lifecycle, geometry metrics, boundary layout, halo engine.
- Produces: `NativeCartesianMgPlan`, `MgHierarchyPolicy`, `LineRelaxationPlan`, `HypreStructAdapter` behind `LinearPreconditioner`.

```cpp
enum class CoarseningKind : std::uint8_t { full_xyz, semi_xy, semi_xz, semi_yz };
struct MgHierarchyPolicy {
  CoarseningKind coarsening;
  double coefficient_change_rebuild_ratio;
  std::uint8_t pre_sweeps, post_sweeps;
};
```

- [ ] **Step 1: Write convergence REDs.**

  Test constant and variable-coefficient Poisson on uniform grids, stretched grids with one and two strong anisotropy directions, Dirichlet/Neumann/periodic combinations, and a null-space pressure case. Require finite monotone residual decrease and mesh-independent or pre-registered bounded iteration growth.

- [ ] **Step 2: Write reuse and semi-coarsening REDs.**

  Verify near-isotropic metrics choose full coarsening; strong stretched directions remain uncoarsened and use line relaxation; unchanged coefficients reuse hierarchy; below-threshold coefficient changes update coarse numeric data without symbolic rebuild; above-threshold changes increment exactly one hierarchy rebuild counter.

- [ ] **Step 3: Implement NativeCartesianMG.**

  Precompute level shapes, metric restriction, transfer spans, halo plans, smoother coloring/lines, null-space projection, and maximum workspace. Implement full weighting/restriction and second-order prolongation consistent with tensor metrics. All levels live in the arena and reuse persistent communication metadata.

- [ ] **Step 4: Implement HYPRE isolation.**

  Make HYPRE optional at configure time. The adapter owns all HYPRE objects and converts from HUNDUN numeric state only when its coefficient identity changes. No HYPRE type appears in installed HUNDUN headers, and disabling HYPRE leaves all native tests/builds intact.

- [ ] **Step 5: Verify and commit.**

  Run all `^v04_solver_mg` tests at 1/2/4 ranks plus `v04_solver_hypre_isolation` with HYPRE off and, when available, on. Record convergence factors and rebuild counters. Commit with `git commit -s -m "feat(v0.4): add native Cartesian multigrid"`.

### Task 13: Compile Static IBM Topology, Quadratic Stencils, and Surface Plans

**Files:**
- Create: `versions/v0.4/include/hundun/v04_ibm.hpp`
- Create: `versions/v0.4/src/mesh_ibm_topology.cpp`
- Create: `versions/v0.4/src/mesh_ibm_reconstruction.cpp`
- Create: `versions/v0.4/src/mesh_ibm_qr_detail.hpp`
- Create: `versions/v0.4/src/mesh_ibm_surface.cpp`
- Create: `versions/v0.4/tests/unit/mesh_ibm_quadratic_test.cpp`
- Create: `versions/v0.4/tests/mpi/mesh_ibm_stencil_mpi_test.cpp`
- Create: `versions/v0.4/tests/numerical/mesh_ibm_order_test.cpp`
- Create: `versions/v0.4/tests/data/ibm_plane.stl`
- Create: `versions/v0.4/tests/data/ibm_sphere.stl`
- Create: `versions/v0.4/tests/data/ibm_cylinder.stl`
- Modify: `versions/v0.4/CMakeLists.txt`
- Modify: `versions/v0.4/tests/CMakeLists.txt`

**Interfaces:**
- Consumes: Cartesian geometry and STL scan, halo reach, field registration, boundary plans.
- Produces: `EBTopology`, `QuadraticStencilPlan`, `BoundaryStencilPlan`, `SurfaceQuadraturePlan`; registers compact interface storage before Task 18.

```cpp
struct QuadraticStencilQuality {
  std::uint8_t donor_count, normal_band_count, quadrant_mask, reach;
  double condition_estimate, functional_l1;
  PlanFingerprint pivot_fingerprint;
};
struct EBTopology {
  Span<const std::uint8_t> region;
  Span<const std::uint32_t> interface_cells;
  RevisionToken geometry_revision;
  PlanFingerprint fingerprint;
};
```

- [ ] **Step 1: Port scientific contracts as independent REDs.**

  Read the accepted reference headers/tests at commit `6407cd7` for contracts only: `ib_quadratic_reconstruction`, `ib_ghost_stencil_plan`, `ib_surface`, and their unit/MPI tests. Re-express the ten-term cell-average quadratic basis, independent 12-tetrahedron moment oracle, deterministic pivot/tie behavior, constrained Dirichlet/Neumann forms, and rejection cases in new v0.4 tests; do not copy implementations.

- [ ] **Step 2: Write standard-stencil geometry REDs.**

  Require 14--32 unique fluid donors, reach at most four, strictly positive wall-normal coordinates, at least three normal bands, all four tangential quadrants, full rank, condition at most the registered limit, deterministic global-ID ordering, and no linear/constant fallback. Test duplicate donors, cavity-side donors, rank loss, poor condition, reversed partition, and the known positive-normal WALE donor geometry.

- [ ] **Step 3: Implement static topology and deterministic QR.**

  Build region/interface lists from scan output; compute wall points/normals; collect candidates by bounded logical neighborhoods; filter by fluid region and positive normal; score 14--32 donor subsets; solve weights with an independently written deterministic column-pivoted Householder QR; store compact immutable weights and fingerprints.

- [ ] **Step 4: Compile surface quadrature separately.**

  Build oriented surface elements, normals, areas, wall points, owner interface cells, and any explicitly declared shared-row donor union. Keep shared quadrature plans distinguishable from standard link stencils so they cannot silently bypass the positive-normal/32-donor contract.

- [ ] **Step 5: Add true h-refinement order tests.**

  Run 12/24/48 or 16/32/64 sequences for translated plane, sphere, and cylinder on uniform and stretched grids. Task 13 measures the quantities its static reconstruction API owns: wall value, wall-normal gradient, and near-wall value/penetration reconstruction. Task 17 adds the independent final-state traction and integrated-force oracle/order gates after the pressure/operator/final-state APIs exist; those gates may not be fabricated or silently treated as already passed here. Require decreasing errors and both adjacent observed orders at least 1.8 for metrics whose analytic truncation is above the floor. Run the static-plan decomposition selectors at 1/2/4 ranks.

- [ ] **Step 6: Verify and commit.**

  Run all `^v04_mesh_ibm_(quadratic|stencil_mpi_[124]|order|plan_[124])$` tests. Inspect every failure path to ensure it rejects instead of lowering order. Commit with `git commit -s -m "feat(v0.4): compile static quadratic IBM plans"`.

### Task 14: Implement Unified Momentum, Enthalpy, Species, and Scalar Equations

**Files:**
- Create: `versions/v0.4/include/hundun/v04_flow.hpp`
- Create: `versions/v0.4/src/solver_equations.cpp`
- Create: `versions/v0.4/src/solver_momentum.cpp`
- Create: `versions/v0.4/src/solver_enthalpy.cpp`
- Create: `versions/v0.4/src/solver_species.cpp`
- Create: `versions/v0.4/tests/numerical/solver_low_mach_mms_test.cpp`
- Create: `versions/v0.4/tests/numerical/solver_enthalpy_terms_test.cpp`
- Create: `versions/v0.4/tests/numerical/solver_species_conservation_test.cpp`
- Create: `versions/v0.4/tests/mpi/solver_equations_mpi_test.cpp`
- Modify: `versions/v0.4/CMakeLists.txt`
- Modify: `versions/v0.4/tests/CMakeLists.txt`

**Interfaces:**
- Consumes: conservative kernels, thermo/transport/contribution plans, boundary/time plans, numeric linear state.
- Produces: `MomentumEquationPlan`, `EnthalpyEquationPlan`, `SpeciesEquationPlan`, `ScalarEquationPlan`, common `EquationAssemblyContext`.

```cpp
struct EquationAssemblyContext {
  double dt, bdf_a0, bdf_a1, bdf_a2;
  RevisionToken time, geometry, boundary, thermo, transport, face_flux;
  Span<const FieldView> accepted, previous, trial;
};
```

- [ ] **Step 1: Write low-Mach manufactured REDs.**

  Use analytic variable `p_abs`, temperature, density, velocity, and enthalpy with nonzero divergence. Check mass, momentum, enthalpy, and pressure-work residuals on uniform/stretched 16/32/64 grids; require both observed orders at least 1.8.

- [ ] **Step 2: Isolate energy-term REDs.**

  Create cases where only one of `Dp/Dt`, conduction, viscous dissipation, unsteady enthalpy, or advective enthalpy is nonzero. Compare each assembled cell integral with an independent analytic oracle and verify signs/units.

- [ ] **Step 3: Write species/scalar conservation REDs.**

  Advect/diffuse N-1 inert species and a passive scalar using a supplied final face flux. Check sum closure, positivity rejection, composition-dependent EOS, global conservation including boundary flux, and no chemistry/source registration.

- [ ] **Step 4: Implement common assembly ownership.**

  Equation plans declare reads/writes/contributions once and use the same Cartesian kernels. Fused loops form `rhoU/rhoh` contributions transiently from primitive and derived views. No PISO-specific duplicate equation implementation and no persistent conservative-product fields are allowed.

- [ ] **Step 5: Implement open/closed pressure reference hooks.**

  Assemble local compressibility terms with `drho/dp`. Open cases consume absolute-pressure boundary closure; closed cases call the Task 8 mass-Newton service at the declared stage. Keep the pressure solve itself for Task 15.

- [ ] **Step 6: Verify and commit.**

  Run all `^v04_solver_(low_mach_mms|enthalpy_terms|species_conservation|equations_mpi_[124])$` tests and inspect field registration to confirm no constant-density or conservative-product field. Commit with `git commit -s -m "feat(v0.4): add unified low-Mach equations"`.

### Task 15: Implement Exactly Two PISO Correctors and Final Flux Publication

**Files:**
- Modify: `versions/v0.4/include/hundun/v04_flow.hpp`
- Create: `versions/v0.4/src/solver_piso.cpp`
- Create: `versions/v0.4/src/solver_piso_detail.hpp`
- Create: `versions/v0.4/src/solver_pressure.cpp`
- Create: `versions/v0.4/tests/unit/solver_piso_authority_test.cpp`
- Create: `versions/v0.4/tests/unit/solver_piso_mutation_test.cpp`
- Create: `versions/v0.4/tests/numerical/solver_piso_checkerboard_test.cpp`
- Create: `versions/v0.4/tests/numerical/solver_piso_temporal_order_test.cpp`
- Create: `versions/v0.4/tests/mpi/solver_piso_mpi_test.cpp`
- Modify: `versions/v0.4/CMakeLists.txt`
- Modify: `versions/v0.4/tests/CMakeLists.txt`

**Interfaces:**
- Consumes: equation plans, `NativeCartesianMG`, Krylov solvers, transaction state, final face-flux writer token.
- Produces: `PisoPlan`, `MomentumIntermediates`, `PressureCorrectionState`, `advance_piso_attempt`.

```cpp
struct MomentumIntermediates {
  FieldView rAU, rAtU, HbyA;
  FaceFluxView phiHbyA;
  RevisionToken rAU_identity, rAtU_identity, HbyA_identity, phiHbyA_identity;
};
Status advance_piso_attempt(const PisoPlan&, AttemptTransaction&,
                            const EquationAssemblyContext&,
                            LinearSolveResult pressure_results[2]) noexcept;
```

- [ ] **Step 1: Write authority and call-count REDs.**

  Require exactly two pressure solves/corrector calls, one final-flux publication, corrector 1 writes trial only, corrector 2 publishes final `U/face_flux`, enthalpy/species see only final flux, and any configured count other than two is rejected during case compilation.

- [ ] **Step 2: Write intermediate-revision mutation REDs.**

  Mutate momentum diagonal, consistent diagonal, trial `U`, trial `phi`, boundary, time, and geometry separately. Assert the exact invalidations: `rAU`; `rAtU`; `HbyA`; and `phiHbyA`. Specifically mutate corrector-1 trial `U/phi` and prove a stale corrector-2 `HbyA/phiHbyA` path fails RED.

- [ ] **Step 3: Write pressure-coupling REDs.**

  Check pressure matrix includes `drho/dp`, Rhie--Chow-like momentum weighting removes checkerboard modes, open/closed null-space handling, final continuity uses pressure-equation flux, and retry/failure publishes neither flux nor intermediate cache.

- [ ] **Step 4: Implement momentum/intermediate lifecycle.**

  Assemble or refill momentum numeric state; form `rAU`, optional `rAtU`, `HbyA`, and `phiHbyA` with separate dependency fingerprints. Reuse only certified intermediates. Corrector 2 necessarily rebuilds/revalidates trial-state-dependent quantities after corrector 1.

- [ ] **Step 5: Implement the two pressure corrections.**

  Assemble pressure RHS/operator from `phiHbyA`, time coefficients, EOS derivative, and BCs; solve; correct trial face flux directly from pressure-equation flux; update trial `U` from the same pressure-gradient path. Only the second correction acquires the final-writer token and publishes revisions.

- [ ] **Step 6: Verify spatial/temporal/MPI contracts.**

  Run checkerboard, variable-step BDF2 order, conservation, failure injection, and 1/2/4-rank tests. Require BE startup then observed BDF2 temporal order at least 1.8 above the floor and exact two-call diagnostics.

- [ ] **Step 7: Review and commit.**

  Inspect for a third corrector, a user corrector loop, final flux reconstructed from final `U`, or a common `accepted_U` dependency tuple. Commit with `git commit -s -m "feat(v0.4): add exactly-two-corrector PISO"`.

### Task 16: Add WALE and Default Vreman Wall-Function Turbulence

**Files:**
- Modify: `versions/v0.4/include/hundun/v04_physics.hpp`
- Create: `versions/v0.4/src/physics_wale.cpp`
- Create: `versions/v0.4/src/physics_vreman.cpp`
- Create: `versions/v0.4/src/physics_wall_function.cpp`
- Create: `versions/v0.4/tests/unit/physics_wale_test.cpp`
- Create: `versions/v0.4/tests/unit/physics_vreman_test.cpp`
- Create: `versions/v0.4/tests/unit/physics_wall_function_test.cpp`
- Create: `versions/v0.4/tests/mpi/physics_turbulence_mpi_test.cpp`
- Modify: `versions/v0.4/CMakeLists.txt`
- Modify: `versions/v0.4/tests/CMakeLists.txt`

**Interfaces:**
- Consumes: shared `DerivedFieldPlan`, velocity-gradient cache, wall distance, `EBTopology`, boundary plans, one `mu_eff` authority.
- Produces: `TurbulencePlan` for `none`, `wale`, and `vreman_wall_function`; no additional model choices.

```cpp
struct TurbulencePlan {
  TurbulenceKind kind;
  double coefficient, turbulent_prandtl, turbulent_schmidt;
  Status update(const FieldView& rho, const GradientView& grad_u,
                const WallDistanceView&, FieldView& mu_eff) const noexcept;
};
```

- [ ] **Step 1: Write WALE invariant REDs.**

  Check zero viscosity for rigid rotation and pure shear limits where required by the formula, nonnegative finite output, correct cubic wall scaling on analytic gradients, uniform/stretched metric identity, and one shared gradient-cache evaluation per `U` revision.

- [ ] **Step 2: Write Vreman and wall-function REDs.**

  Use analytic velocity-gradient tensors for zero/nonzero Vreman invariants; test nonnegative finite eddy viscosity, external and IBM wall faces, smooth/log-law transition, rough invalid input rejection, tangential relative velocity for moving walls, heat/scalar wall flux consistency, and single `mu_eff` publication.

- [ ] **Step 3: Reproduce the positive-normal donor failure as RED.**

  Construct the previously failing ghost geometry and require the standard IBM donor contract. The test fails if turbulence/wall reconstruction selects a nonpositive-normal donor; it may not skip the wall cell, relax coverage, or fall back to lower order.

- [ ] **Step 4: Implement static model binding.**

  Compile a single function-level model selection before the cell loop. `none` copies molecular viscosity into `mu_eff`; WALE and Vreman consume the same gradient view; wall function consumes shared wall distance and boundary/EB surface plans. Do not duplicate gradient or viscosity fields.

- [ ] **Step 5: Fix donor selection by topology/reach root cause.**

  If the RED fails, diagnose whether classification, periodic image, candidate reach, wall normal orientation, or halo availability is wrong. Correct that upstream plan and retain all 14--32/three-band/four-quadrant/positive-normal requirements.

- [ ] **Step 6: Verify and commit.**

  Run all turbulence tests at 1/2/4 ranks plus repeated-update allocation/cache counters. Assert case compilation accepts only the three model choices and defaults to Vreman wall function. Commit with `git commit -s -m "feat(v0.4): add WALE and Vreman wall functions"`.

### Task 17: Implement Exact IBM Pressure Correction and Final-State Force Authority

**Files:**
- Modify: `versions/v0.4/include/hundun/v04_ibm.hpp`
- Create: `versions/v0.4/src/solver_ibm_pressure.cpp`
- Create: `versions/v0.4/src/solver_ibm_force.cpp`
- Create: `versions/v0.4/src/solver_ibm_force_detail.hpp`
- Create: `versions/v0.4/tests/unit/solver_ibm_pressure_operator_test.cpp`
- Create: `versions/v0.4/tests/numerical/solver_ibm_force_oracle_test.cpp`
- Create: `versions/v0.4/tests/numerical/solver_ibm_force_mutation_test.cpp`
- Create: `versions/v0.4/tests/mpi/solver_ibm_force_mpi_test.cpp`
- Modify: `versions/v0.4/CMakeLists.txt`
- Modify: `versions/v0.4/tests/CMakeLists.txt`

**Interfaces:**
- Consumes: static IBM plans, exact final-state PISO contract, FGMRES, Native MG preconditioner, shared gradient/`mu_eff`, transaction cache publication.
- Produces: `IbmPressureOperator`, `IbmPressurePreconditioner`, `FinalSurfaceState`, `SurfaceForce`, `FinalForceCache`.

```cpp
struct SurfaceForce { Real3 pressure, viscous, total, moment; RevisionToken revision; };
struct FinalSurfaceState {
  RevisionToken final_u, final_pressure, final_gradient, mu_eff, geometry;
};
Status evaluate_surface_force(const SurfaceQuadraturePlan&,
                              const FinalSurfaceState&, SurfaceForce&) noexcept;
```

- [ ] **Step 1: Write an independent final-state force oracle first.**

  Build analytic constant/linear/quadratic pressure and velocity fields on plane, sphere, and cylinder surfaces. Integrate traction directly from analytic functions and triangle quadrature, outside all product reconstruction/operator code. Compare product pressure, viscous, total force, and moment with this oracle.

- [ ] **Step 2: Prove mutation sensitivity while still RED.**

  Add mutations that substitute corrector-1 pressure gradient, stale `U`, stale `mu_eff`, reversed normal, omitted area, wrong pressure sign, and provisional flux. Each mutation must make a specifically named assertion fail. Record the RED logs before implementing the candidate path.

- [ ] **Step 3: Implement regular-body plus compact-interface pressure apply.**

  The exact operator applies regular Cartesian coefficients plus immutable compact IBM rows. FGMRES always applies this exact operator; Native MG/compact approximation is preconditioner only. Store no duplicate full sparse matrix and rebuild compact rows only on geometry/topology revision.

- [ ] **Step 4: Implement one final-state traction path.**

  Reconstruct final pressure and velocity gradients from the final accepted candidate revision, combine molecular+turbulent stress through the one `mu_eff`, integrate oriented quadrature once, and publish `FinalForceCache` only in the successful transaction commit. Reaction budget and diagnostics read this same result.

- [ ] **Step 5: Verify order, consistency, failure, and MPI.**

  Run analytic oracle/mutation tests, h-refinement traction/force tests, operator-vs-reaction consistency, failed-attempt cache invisibility, and 1/2/4-rank decomposition. Require every registered mutation to fail on the mutated binary and pass on the production binary.

- [ ] **Step 6: Review and commit.**

  Inspect that no force path reads `HbyA`, corrector scratch, or provisional flux; no surface normal is recomputed in the hot loop; and no diagnostic has a second force implementation. Commit with `git commit -s -m "feat(v0.4): add exact IBM pressure and final force"`.

### Task 18: Freeze the Complete Production Schema, Halo, Solvers, and Graph Once

**Files:**
- Modify: `versions/v0.4/include/hundun/v04_field.hpp`
- Modify: `versions/v0.4/include/hundun/v04_parallel.hpp`
- Modify: `versions/v0.4/include/hundun/v04_execution.hpp`
- Create: `versions/v0.4/src/core_product_freeze.cpp`
- Create: `versions/v0.4/src/core_product_freeze_detail.hpp`
- Create: `versions/v0.4/tests/integration/core_product_freeze_test.cpp`
- Create: `versions/v0.4/tests/integration/core_hot_resource_test.cpp`
- Create: `versions/v0.4/tests/mpi/core_product_freeze_mpi_test.cpp`
- Modify: `versions/v0.4/CMakeLists.txt`
- Modify: `versions/v0.4/tests/CMakeLists.txt`

**Interfaces:**
- Consumes: every registration and plan from Tasks 3--17.
- Produces: `ProductPlanBundle` containing final `FieldSchema`, `ArenaLayout`, `CommunicationPlan`, boundary/EB/operator/solver plans, `FrozenExecutionGraph`, and exact resource contracts. This is the only production-freeze API.

```cpp
struct ProductPlanBundle {
  FieldSchema fields;
  ArenaLayout arena;
  CommunicationPlan communication;
  FrozenExecutionGraph graph;
  PlanFingerprint fingerprint;
};
Status freeze_product(const ValidatedModel&, ProductPlanBundle&) noexcept;
```

- [ ] **Step 1: Write incomplete-registration REDs.**

  Omit one field, ghost requirement, cache dependency, operator coefficient, solver workspace, service snapshot, or collective point at a time. Require `freeze_product` to reject before allocation. Attempt post-freeze registration and require immutable rejection.

- [ ] **Step 2: Write the exact freeze-order test.**

  Instrument and require this order:

  ```text
  complete capability registration -> FieldSchema -> ArenaLayout
  -> Boundary/EB plans -> Symbolic/Numeric/Hierarchy/Workspace plans
  -> merged CommunicationPlan -> FrozenExecutionGraph -> resource seal
  ```

  Assert every stage and field is registered before the first freeze event and every bundle fingerprint is identical at 1/2/4 ranks.

- [ ] **Step 3: Implement merged production plans.**

  Gather all field/stage specs, compute maximum ghost width/payload/workspace, freeze schema and arena, bind every view, instantiate static EB/boundary/linear plans, reserve peer buffers/requests, compile the complete stage graph, and seal counters. No later task may add a product field or hot stage without reopening Task 18 and invalidating downstream evidence.

- [ ] **Step 4: Prove hot-resource contracts.**

  Execute a synthetic complete step 100 times. Require zero heap events, stable field/workspace/buffer/request addresses, exact two correctors, one final-flux writer, no geometry query, no strings, message/byte counts within seal, and only policy-authorized numeric refill/hierarchy rebuild counts.

- [ ] **Step 5: Verify tests-off isolation and commit.**

  Build tests-on and tests-off. Inspect link maps/symbols so full-domain gather, mutation seams, and test oracles are absent from tests-off `hundun`. Run product-freeze tests at 1/2/4 ranks and commit with `git commit -s -m "feat(v0.4): freeze complete production plans"`.

### Task 19: Implement Product Driver, Minimal Restart, and Runtime Output

**Files:**
- Create: `versions/v0.4/include/hundun/v04_io.hpp`
- Create: `versions/v0.4/include/hundun/v04_app.hpp`
- Create: `versions/v0.4/src/app_driver.cpp`
- Modify: `versions/v0.4/src/app_main.cpp`
- Create: `versions/v0.4/src/io_restart.cpp`
- Create: `versions/v0.4/src/io_visit.cpp`
- Create: `versions/v0.4/src/io_screen.cpp`
- Create: `versions/v0.4/src/io_monitor.cpp`
- Create: `versions/v0.4/src/io_evidence.cpp`
- Create: `versions/v0.4/tests/integration/app_driver_test.cpp`
- Create: `versions/v0.4/tests/mpi/io_restart_mpi_test.cpp`
- Create: `versions/v0.4/tests/integration/io_product_path_test.cpp`
- Create: `versions/v0.4/tests/integration/app_init_case_test.cpp`
- Create: `versions/v0.4/docs/evidence-schema.md`
- Modify: `versions/v0.4/CMakeLists.txt`
- Modify: `versions/v0.4/tests/CMakeLists.txt`

**Interfaces:**
- Consumes: immutable `ProductPlanBundle`, `AttemptTransaction`, committed snapshots only.
- Produces: `hundun validate`, `hundun run`, `hundun init-case`, `RestartWriter/Reader`, VTI/VTR+Visit output, screen/monitor/evidence records.

- [ ] **Step 1: Write driver/CLI REDs.**

  Test flat case validation, dry plan print, run-directory separation, refusal to write into the source tree, `init-case --output <run-dir>` generation, and no bundled product case/example under `versions/v0.4`.

- [ ] **Step 2: Write minimal Restart REDs.**

  Save only current `U`, `p_ref/pi` or `p_abs`, `h`, N-1 species, final face mass flux, controller/stat state, schema/plan fingerprints, units, time, and integrity metadata. Assert derived `rho/T/mu`, `rhoU/rhoh`, and BDF2 previous history are absent. Corrupt/mismatched input must validate on all ranks before any state publication.

- [ ] **Step 3: Write overwrite/rank-change REDs.**

  With default `keep_last=1`, write pending, fsync/close, validate, atomically switch `current`, then delete the previous generation. Inject failure at every boundary and ensure at least one valid current restart remains. Restart 1->2, 2->4, and 4->1 ranks; first resumed step must use BE and subsequent steps variable BDF2.

- [ ] **Step 4: Implement committed-snapshot services.**

  Services receive immutable field views plus accepted revisions only. The driver executes the frozen graph, performs collective attempt finish, and only then schedules Restart/Visit/screen/monitor. VTI handles uniform grids; VTR handles stretched grids; `.visit` indexes rank/time files without full-domain gather.

- [ ] **Step 5: Implement structured evidence and tests-off equality.**

  Record stage wall time, init, hot step, RSS, halo messages/bytes, iterations, refills/rebuilds, allocation count, compiler/build/case/STL/binary fingerprints. Run the same minimal case with tests-on/off products and require numerically identical committed outputs while oracle symbols remain absent tests-off.

- [ ] **Step 6: Verify and commit.**

  Run CLI, driver, Restart 1/2/4-rank, corruption, rank-change, product-path, Visit schema, and `init-case` tests. Confirm Restart storage remains bounded over repeated writes. Commit with `git commit -s -m "feat(v0.4): add committed-state runtime services"`.

### Task 20: Freeze Focused Gates, Literature Data, and Candidate Evidence Workflow

**Files:**
- Create: `versions/v0.4/tests/focused_manifest.cmake`
- Create: `cmake/HundunV04CandidateIdentity.cmake`
- Create: `tools/v04_evidence_validate.py`
- Create: `tools/v04_literature_extract.py`
- Create: `docs/references/cylinder-re3900-parnaudeau.json`
- Create: `docs/references/cylinder-re3900-norberg.json`
- Create: `docs/references/backward-step-driver-seegmiller.json`
- Create: `docs/verification/v0.4-re3900-candidate-ledger.md`
- Create: `docs/verification/v0.4-literature-data-receipt.md`
- Modify: `versions/v0.4/docs/evidence-schema.md`
- Modify: `versions/v0.4/tests/CMakeLists.txt`

**Interfaces:**
- Consumes: complete tests-off/on product, all focused selectors, benchmark and literature identities.
- Produces: immutable gate state `focused -> full2 -> frozen -> full20 -> literature -> final`, candidate manifest validator, pre-registered experimental metrics/tolerances, and exact evidence receipts.

- [ ] **Step 1: Write gate-order mutation REDs.**

  Attempt to register `full2`, `frozen`, `full20`, or `literature` without the preceding accepted receipt; attempt to change source/case/STL/build after freeze; attempt to overwrite an evidence path. Require `REJECT` and unchanged ledger.

- [ ] **Step 2: Define the focused manifest.**

  Include EOS/`drho_dp`, enthalpy terms, closed mass, schemes, two-corrector authority/mutations, checkerboard/temporal order, IBM quadratic/order, force oracle/all mutations, positive-normal donor, WALE/Vreman, Restart/rollback/rank-change, public headers, 1/2/4 MPI, ASan, and UBSan. Do not include a small-grid performance test.

- [ ] **Step 3: Freeze literature data before simulations.**

  Encode Parnaudeau/Norberg citation/DOI, definitions, units, experimental values/uncertainties, and digitized mean/fluctuation profiles at `x/D=1.06,1.54,2.02`. Store source page/figure/table and extraction uncertainty per datum. Encode Driver--Seegmiller/NASA geometry and reattachment/profile data for later Task 22. Hash the extraction script and source attachments; no threshold may depend on HUNDUN output.

- [ ] **Step 4: Define exact candidate identity.**

  Manifest fields: HEAD/tree, compiler, flags, tests-on/off hashes, case/`.d`/STL hashes, MPI version, rank/process grid, CPU/NUMA/affinity, environment allowlist, command, start/end timestamps, exit, stdout/stderr and evidence SHA-256. Include process inventory before/after without terminating unrelated processes.

- [ ] **Step 5: Define performance statistics.**

  Full-grid 2-step is directional only. Formal full20 requires at least five alternating HUNDUN/COAST pairs on frozen resources. Compare hot-step median and report P90, init, RSS, communication, iterations/refills/rebuilds. Release predicate is `HUNDUN median <= COAST median`; `<=1.10x` records `NEAR`, never `ACCEPT`.

- [ ] **Step 6: Verify tooling without launching cases and commit.**

  Run synthetic manifests, corrupted hashes, reordered gates, changed-candidate mutations, literature-schema validation, and focused selector enumeration. Do not run full grid or long statistics in this task. Commit with `git commit -s -m "test(v0.4): freeze Re3900 evidence workflow"`.

### Task 21: Execute the Re=3900 Cylinder Joint Release Gate

**Files:**
- Runtime input only: a user-selected run directory containing flat `case.json`, `.d`, and STL files
- Modify: `docs/verification/v0.4-re3900-candidate-ledger.md`
- Create only after execution: `docs/verification/v0.4-re3900-final-acceptance.md`
- External evidence only: `/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/**`

**Interfaces:**
- Consumes: exact frozen candidate and Task 20 workflow.
- Produces: one `ACCEPT` or `REJECT`; only `ACCEPT` authorizes a v0.4 release operation.

- [ ] **Step 1: Seal pre-run state.**

  Record clean status, complete diff, HEAD/tree, tests-off/on hashes, case/`.d`/STL hashes, compiler/MPI/CPU/NUMA/affinity, available memory, existing HUNDUN/COAST/MPI processes, and complete evidence-directory inventory. Do not stop unrelated processes and never overwrite an old run.

- [ ] **Step 2: Run the focused gate.**

  Execute the frozen focused manifest in Release plus ASan/UBSan subsets. Any failure produces `REJECT` before a cylinder run. Require the force mutation RED/green evidence and positive-normal donor repair explicitly.

- [ ] **Step 3: Run paired full-grid two-step.**

  Use `coast_pairing_short`: `20D x 20D x piD`, cylinder at `x=5D`, `480x480x48`, 64 ranks, frozen process grid/affinity, `Re_D=3900`, `dt U/D=0.006`, and nonperiodic slip/symmetry transverse/span boundaries that both products can represent. Run HUNDUN and COAST with equivalent short-test physics and record init, both steps, RSS, communication, and iterations. This is a performance-direction screen, not a final hot-step conclusion.

- [ ] **Step 4: Freeze or return to the owning task.**

  Freeze only if focused correctness accepts, no NaN/deadlock/unhandled failure exists, full2 resource counters match contracts, and the measured direction makes `median <= COAST` credible. Record the exact immutable candidate. Any code/config/build change restarts at Step 1.

- [ ] **Step 5: Run at least five alternating full-grid 20-step pairs.**

  Alternate launch order to reduce drift. Use identical resources, measure post-startup hot steps, and compute median/P90 with the frozen tool. Require HUNDUN median no higher than COAST, no scientific-threshold relaxation, stable continuity/conservation, and no systematic iteration/refill/rebuild regression outside the pre-registered contract.

- [ ] **Step 6: Run HUNDUN-only literature statistics.**

  Use `literature_statistics`: Parnaudeau `20D x 20D x piD`, cylinder at `x=5D`, transverse/spanwise periodic boundaries, `480x480x48`, `dt U/D=0.006`, `Re_D=3900`, WALE. Develop at least `150D/U`; then collect approximately `2020D/U` (about 420 shedding periods) only from the frozen candidate. Compare `St`, `Lr/D`, `Cd_mean`, `Cl_rms`, centerline mean velocity, and mean/fluctuation profiles at the three registered x stations against pre-frozen experiment data and uncertainty rules.

- [ ] **Step 7: Audit provenance and publish one decision.**

  Verify exact hashes, complete diff, no copied GPL/COAST implementation, DCO trailers, tests-off isolation, no task-owned MPI process, and evidence immutability. Compute:

  ```text
  release = numerical_correctness_accept
         && robustness_accept
         && coast_short_performance_accept
         && literature_physical_accuracy_accept
         && provenance_accept
  ```

  Write one `ACCEPT` or `REJECT` report. A tag, push, or public release remains a separate user-authorized repository operation.

### Task 22: Validate the Driver--Seegmiller Backward-Facing Step without Source Changes

**Files:**
- Runtime input only: a user-selected backward-step case directory with `case.json`, direct-root `.d` tables, and `backward_step.stl`
- Create after execution: `docs/verification/v0.4-backward-step-validation.md`
- External evidence only: `/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-backward-step/**`

**Interfaces:**
- Consumes: Re3900 `ACCEPT` candidate, frozen Driver--Seegmiller data, tensor-stretched mesh, Vreman wall function, typed inlet-profile tables.
- Produces: evidence that a second literature case is configured without modifying product source. This task does not block the Re3900 release decision.

- [ ] **Step 1: Build the case entirely from runtime input.**

  Use the frozen NASA/Driver--Seegmiller geometry: step-height/tunnel-exit ratio 1:9 encoded by `backward_step.stl`, width identity when modeled in 3D, `U_ref` and Reynolds number from the frozen receipt, tensor-product stretching near walls/step, measured or documented inlet boundary-layer profile from direct-root `.d`, pressure outlet, no-slip walls, and default `vreman_wall_function`.

- [ ] **Step 2: Validate input and run focused case preflight.**

  Require `hundun validate` and plan fingerprints without source changes. Check wall-function y-plus plan, mass/pressure closure, inlet profile integration, mesh limits, and initial solver hierarchy before any long run.

- [ ] **Step 3: Run and compare literature observables.**

  Advance to statistically stable behavior using the pre-registered convergence/statistics policy. Compare reattachment length (reference about `6.26h` with frozen uncertainty), lower-wall `Cf/Cp`, mean velocity, and Reynolds-stress profiles at registered stations. Report numerical uncertainty and grid sensitivity using only global tensor-product refinement/stretching.

- [ ] **Step 4: Enforce the generality rule.**

  If the case requires a product-source edit, stop validation and classify the missing capability. Implement only a generic correction through the owning earlier task, rerun all affected focused tests and the complete Re3900 release gate, then restart this task. Never add `backward_step`, station, or geometry-name branches to product code.

- [ ] **Step 5: Record evidence.**

  Seal input/data/binary/resource hashes, statistics, plots, and comparison tables. State whether the case passed without source change. Commit only the validation report/receipts with `git commit -s -m "docs: validate v0.4 backward-facing step"` when complete.

## Dependency and Parallel Review Order

```text
1 -> 2 -> 3
3 -> {4, 5, 7, 8}
4 + 5 -> 6
4 + 5 + 7 -> {9, 10}
4 + 6 + 9 -> 11 -> 12
5 + 6 + 7 + 9 -> 13
8 + 9 + 11 + 12 -> 14 -> 15
8 + 9 + 13 + 14 + 15 -> 16
12 + 13 + 15 + 16 -> 17
10 + 15 + 16 + 17 -> 18 -> 19 -> 20 -> 21 -> 22
```

Tasks inside braces may run concurrently only when workers own disjoint files and agree on the Stable Public Types. The following edges are hard serial dependencies: version isolation before product code; storage/mesh before actual halo compile; boundary/scheme before operators; linear lifecycle before MG; MG before PISO; final-state PISO plus static IBM plus shared gradient/`mu_eff` before force; all physics stages before the single production freeze; and `18 -> 19 -> 20 -> 21 -> 22`.

## Plan Self-Review Checklist

- **Spec coverage:** Tasks 1--3 cover references, version isolation, flat source/input, and registration. Tasks 4--10 cover SoA arenas, transactions, Cartesian/STL, CPU/halo, BC/NSCBC/time, EOS/transport/contributions, conservative kernels, and graph compilation. Tasks 11--17 cover the four linear layers, Native MG/HYPRE, true IBM second-order evidence, unified enthalpy/species equations, exactly two PISO correctors, WALE/Vreman wall functions, and final-state force authority. Tasks 18--20 cover the once-only production freeze, minimal Restart/I/O, tests-off isolation, literature data, and ordered evidence. Tasks 21--22 implement only the approved cylinder and backward-step validation scope.
- **Scope exclusions:** No constant-density fast path, SIMPLE/PIMPLE, body-fitted/AMR, reactions, additional turbulence model, moving EB, product GPU, source-tree cases/examples, small-grid performance gate, or pre-freeze long statistics appears as an implementation deliverable.
- **Type consistency:** `Status`, `FieldId`, `StageId`, `RevisionToken`, `PlanFingerprint`, `FieldRegistry`, `FieldSchema`, `FieldView`, `CartesianGeometryPlan`, `MeshPatch`, `HaloEngine`, `BoundaryPlan`, `ThermodynamicsPlan`, `FaceFluxView`, the four linear layers, `EBTopology`, `AttemptTransaction`, `PisoPlan`, `ProductPlanBundle`, and committed snapshots are defined before consumption and retain one spelling.
- **Lifecycle consistency:** Task 3 tests freeze mechanics but cannot freeze production. Task 6 reserves engine resources but Task 18 compiles the complete merged communication plan. Task 10 tests graph compilation but Task 18 performs the only production freeze. `rAU/rAtU/HbyA/phiHbyA` have distinct dependencies. Corrector 2 is the sole final `U/face_flux` publisher.
- **Validation consistency:** COAST is used only for nonperiodic full-grid short performance. Periodic Parnaudeau long statistics are HUNDUN-only and compare directly to pre-frozen experiments. `24^3` is absent as a performance gate. Backward step follows Re3900 and cannot introduce case-specific source.
- **Placeholder scan:** Before committing this plan, run the command below and treat any match as a plan defect:

  ```bash
  ! rg -n '\x54\x42\x44|\x54\x4f\x44\x4f|\x46\x49\x58\x4d\x45|\x58\x58\x58|implement\x20later|fill\x20in\x20details|similar\x20to\x20Task' docs/superpowers/plans/2026-08-12-hundun-flow-v0.4-cartesian-performance-architecture.md
  git diff --check
  ```

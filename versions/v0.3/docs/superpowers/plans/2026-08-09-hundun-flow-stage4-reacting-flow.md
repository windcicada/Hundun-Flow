# HUNDUN-FLOW Stage 4 Reacting Flow Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:executing-plans` for the default serial path. Use `superpowers:subagent-driven-development` only after the user explicitly authorizes parallel execution at a stage gate. Execute one task at a time and return to the main agent for review.

**Goal:** 在已接受的 Stage 3 CPU-reference 求解器上交付无 Python 的 Cantera C++ 反应流、全组分/总热化学焓守恒输运、固定两次 PISO、Checkpoint v4、diagnostics 和可安装包。

**Architecture:** HUNDUN 保持网格、MPI、PISO、守恒输运、时间推进和事务权威；Cantera 只实现内部 thermo/transport/chemistry backend。热化学变量使用分区二阶 `C(dt/2)-T(dt)-C(dt/2)`，momentum 保留 BDF2，第二个 PISO 是最终约束同步点。

**Tech Stack:** C++17、CMake 3.21+、MPI-3、GCC 11/libstdc++、Cantera C++ 3.2.0、SUNDIALS、CTest、Apache-2.0/DCO。

## Global Constraints

- 实施 parent 必须是执行时正式接受的 Stage 3 product/code HEAD；`4F-0` 重新登记，不使用本计划编写基线代替。
- Linux x86_64、Ubuntu 22.04/glibc 2.35+、GCC 11、C++17、`_GLIBCXX_USE_CXX11_ABI=1` 是唯一发行 profile。
- 正常 configure/build/install/runtime/test 不得要求 Python、Conda 或在线 fetch；不得关闭宿主机网络。
- public header 不得出现 Cantera/SUNDIALS 类型；plugin ABI v1 不扩展。
- 所有 `rho Y_k` 和 `rho h_tc` 使用同一最终面质量通量；`h_tc` 含 formation enthalpy。
- chemistry 不向 `rho h_tc` 添加独立 heat-release source。
- 每个成功 step 恰好两次 PISO；不得加 damping、滤波、第三 corrector 或逐 case 调参。
- 同一 mechanism 的每 rank runtime、每 thread workspace 完全独立。
- NativeChemistryBackend 是 post-v1；当前不实现，但 internal service 和 value types
  不得含 Cantera 类型或假设，必须通过 analytic substitute-backend contract test。
- task gate 只运行直接 RED、unit/header/policy、必要 small MPI 和至多一个 12^3 smoke。
- worker 不访问 BOFFIN、COAST、私有研究数据，不修改 shared registry root，不提交或签署；主 agent完成完整 diff、DCO 和接受。
- `Files` 中出现 integration-owned CMake、root dispatch、registry、`VERSION` 或
  `AGENTS.md` 时，该行是主 agent的集成责任；worker 只返回所需 registration entries，
  不直接修改中央文件。
- Every task's Step 6 commit/receipt/DCO action is main-agent-only. A worker stops after
  Step 5 and returns its diff, commands, outputs and risks without staging or committing.
- 不 push、不发布。

---

## 1. Planned File Map

### Stable public value/config/report types

```text
include/hundun/chem_composition.hpp
include/hundun/chem_reports.hpp
include/hundun/flow_reacting_state.hpp
include/hundun/cfg_resolved_case_v4.hpp
include/hundun/cfg_resolved_case_v4_loader.hpp
include/hundun/flow_checkpoint_v4.hpp
include/hundun/diag_provider_registry.hpp
include/hundun/diag_reacting.hpp
```

### Internal services and composition

```text
src/chem_thermodynamics_service_detail.hpp
src/chem_backend_detail.hpp
src/chem_cantera_backend.cpp
src/chem_workspace_detail.hpp
src/flow_reacting_transaction_detail.hpp
src/flow_reacting_transaction.cpp
src/flow_reacting_transport_detail.hpp
src/flow_reacting_transport.cpp
src/flow_reacting_coupling_detail.hpp
src/flow_reacting_coupling.cpp
src/flow_reacting_boundary_detail.hpp
src/flow_reacting_boundary.cpp
src/cfg_resolved_case_v4_loader.cpp
src/cfg_resolved_case_v4_loader_detail.hpp
src/app_resolved_case_v4_broadcast.cpp
src/app_resolved_case_v4_broadcast_detail.hpp
src/app_reacting_flow_driver.cpp
src/app_reacting_flow_driver_detail.hpp
src/flow_checkpoint_v4.cpp
src/flow_checkpoint_v4_detail.hpp
src/diag_provider_registry.cpp
src/diag_reacting.cpp
```

### Packaging/provenance

```text
cmake/HundunCanteraPackage.cmake
cmake/HundunRelocatablePackage.cmake
third_party/cantera/UPSTREAM.json
third_party/cantera/PATCHES.md
LICENSES/cantera-BSD-3-Clause.txt
THIRD_PARTY_NOTICES
docs/numerics/stage4-capability-ledger.md
```

### Tests

```text
tests/unit/test_chem_composition.cpp
tests/unit/test_chem_services.cpp
tests/support/chem_analytic_backend.hpp
tests/support/chem_analytic_backend.cpp
tests/unit/test_reacting_transaction.cpp
tests/unit/test_resolved_case_v4.cpp
tests/unit/test_checkpoint_v4.cpp
tests/unit/test_diag_provider_registry.cpp
tests/unit/test_cantera_backend.cpp
tests/unit/test_reacting_transport.cpp
tests/unit/test_reacting_coupling.cpp
tests/mpi/test_reacting_rollback.cpp
tests/mpi/test_reacting_decomposition.cpp
tests/mpi/test_checkpoint_v4_mpi.cpp
tests/numerical/test_reacting_mms.cpp
tests/numerical/test_reacting_smoke.cpp
tests/cmake/stage4_package_contract.cmake
tests/cmake/stage4_source_policy.cmake
tests/acceptance/stage4_acceptance.sh
```

Stage 3 最终树若已经提供同职责 v3 registry 或 driver seam，主 agent在 `4F-0`
把上述“Create”改为“Modify”，但不得产生第二个 authority。

## 2. Frozen Internal Interfaces

后续任务必须使用以下名字，不自行发明近义接口：

```cpp
namespace hundun::chemistry {

struct SpeciesIdentity {
  std::string name;
  double molecular_weight_kg_per_kmol;
  std::vector<std::int32_t> element_counts;
};

struct CompositionIdentity {
  std::vector<std::string> element_names;
  std::vector<SpeciesIdentity> species;
  std::uint64_t fingerprint;
};

struct ThermochemicalPoint {
  double p0_pa;
  double h_tc_j_per_kg;
  std::vector<double> mass_fractions;
};

struct ThermodynamicProperties {
  double temperature_k;
  double density_kg_per_m3;
  double cp_j_per_kg_k;
  double mixture_molecular_weight_kg_per_kmol;
};

struct TransportProperties {
  double viscosity_pa_s;
  double conductivity_w_per_m_k;
  std::vector<double> mixture_diffusivity_m2_per_s;
};

struct ChemistryIntervalRequest {
  ThermochemicalPoint state;
  double start_time_s;
  double duration_s;
};

enum class ChemistryStatus : std::uint32_t {
  success,
  invalid_input,
  composition_mismatch,
  state_inversion_failure,
  integration_failure,
  non_finite_output,
  conservation_failure,
  workspace_failure
};

struct ChemistryIntervalReport {
  ThermochemicalPoint final_state;
  std::vector<double> integrated_rho_y_delta_kg_per_m3;
  ChemistryStatus status;
  double completed_duration_s;
  std::uint32_t internal_step_count;

  [[nodiscard]] bool succeeded() const noexcept {
    return status == ChemistryStatus::success;
  }
};

}  // namespace hundun::chemistry
```

For a non-success status, `final_state` is the original request state and every integrated
delta is canonical `+0.0`; `completed_duration_s` and `internal_step_count` are diagnostic
only. Partial backend state/deltas remain workspace-local and can never enter a transaction.

Internal service signatures in `src/*_detail.hpp`:

```cpp
namespace hundun::chemistry {

class ThermodynamicsService {
 public:
  virtual ~ThermodynamicsService() = default;
  virtual const CompositionIdentity& composition() const noexcept = 0;
  virtual ThermodynamicProperties evaluate(
      const ThermochemicalPoint&) const = 0;
};

class TransportPropertyService {
 public:
  virtual ~TransportPropertyService() = default;
  virtual TransportProperties evaluate(
      const ThermochemicalPoint&,
      const ThermodynamicProperties&) const = 0;
};

class ChemistryBackend {
 public:
  virtual ~ChemistryBackend() = default;
  virtual const CompositionIdentity& composition() const noexcept = 0;
  virtual ChemistryIntervalReport integrate(
      const ChemistryIntervalRequest&) = 0;
};

}  // namespace hundun::chemistry
```

### Frozen schema v4 additions

Schema v4 extends the accepted v3 root and preserves all unchanged v3 keys. New keys are exact:

| JSON path | Type / allowed value | Rule |
|---|---|---|
| `schema_version` | integer `4` | required |
| `simulation.type` | `reacting_flow` | required |
| `simulation.density_model` | `ideal_gas` | other v3 density models do not enable chemistry |
| `reacting.chemistry.backend` | `cantera` | only v1 production backend |
| `reacting.chemistry.mechanism.file` | path | required, case-root confined |
| `reacting.chemistry.mechanism.sha256` | 64 lowercase hex | required and verified before allocation |
| `reacting.chemistry.mechanism.phase` | nonempty string | required |
| `reacting.chemistry.relative_tolerance` | finite positive real | required, no hidden default |
| `reacting.chemistry.absolute_tolerance` | finite positive real | required, no hidden default |
| `reacting.chemistry.maximum_internal_steps` | positive integer | required |
| `reacting.thermodynamics.initial_p0_pa` | finite positive real | required |
| `reacting.thermodynamics.initial_temperature_k` | finite positive real | required |
| `reacting.thermodynamics.initial_mass_fractions` | ordered array of `{species,value}` | exactly all mechanism species once |
| `reacting.transport.model` | `mixture_averaged` | required |
| `reacting.pressure_constraint.mode` | `open_fixed_p0`, `closed`, `partially_closed` | required |
| `boundaries[*].reacting.species` | `non_catalytic_impermeable` or typed inlet/outlet state | catalytic rejected |
| `boundaries[*].reacting.thermal.mode` | `adiabatic` or `isothermal` | required on walls |
| `boundaries[*].reacting.thermal.temperature_k` | finite positive real | present iff isothermal |

All new objects reject unknown/duplicate members. `--print-resolved` emits mass fractions in
mechanism order and includes mechanism/composition fingerprints. Backend tolerances are case
identity and enter Checkpoint/Restart compatibility；the loader never tunes them by case type.

## 3. Task Sequence

### Common command protocol

Before 4P-2, backend-neutral tasks use:

```bash
cmake -S . -B build/stage4-debug \
  -DCMAKE_BUILD_TYPE=Debug \
  -DHUNDUN_BUILD_TESTS=ON
cmake --build build/stage4-debug -j32
ctest --test-dir build/stage4-debug -R "${HUNDUN_TASK_TEST_REGEX}" \
  --output-on-failure -j24
```

4P-2 freezes `HUNDUN_ENABLE_CANTERA` and `HUNDUN_CANTERA_PACKAGE_ROOT`. Tasks that consume
the verified artifact use:

```bash
cmake -S . -B build/stage4-debug \
  -DCMAKE_BUILD_TYPE=Debug \
  -DHUNDUN_BUILD_TESTS=ON \
  -DHUNDUN_ENABLE_CANTERA=ON \
  -DHUNDUN_CANTERA_PACKAGE_ROOT="${HUNDUN_STAGE4_CANTERA_ROOT}"
cmake --build build/stage4-debug -j32
ctest --test-dir build/stage4-debug -R "${HUNDUN_TASK_TEST_REGEX}" \
  --output-on-failure -j24
```

The main agent exports `HUNDUN_STAGE4_CANTERA_ROOT` from the accepted 4P-2 receipt and sets
`HUNDUN_TASK_TEST_REGEX` to the exact test basename(s) listed in the task. Focused Release
replaces the build directory/type with `build/stage4-release`/`Release`. MPI tests are
registered through CTest；a task invokes raw `mpiexec` only when its step says so. A worker
records the expanded command, environment, exit and log path.

### Task 4F-0: Accepted Stage 3 Intake and Compatibility Inventory

**Depends on:** User starts Stage 4 after Stage 3 acceptance.

**Files:**
- Create: `.superpowers/sdd/stage4-4F-0-baseline-receipt.md`
- Create: `docs/numerics/stage4-capability-ledger.md`
- Modify: `AGENTS.md`
- Test: `tests/cmake/stage4_source_policy.cmake`

**Interfaces:** Consumes the accepted Stage 3 HEAD, public headers, schema v1--v3,
Restart v1--v3, diagnostics registry and installed product. Produces the exact Stage 4
parent, affected-file inventory and a no-product-change compatibility gate.

- [ ] **Step 1: Freeze evidence.** Record `git rev-parse HEAD`, parent/tree, `git status
  --porcelain=v1`, public-header list, exported symbols, schema/Restart/diagnostic IDs and
  background MPI processes. The receipt must say `stage4_product_changes=none`.
- [ ] **Step 2: Add the failing source-policy fixture.** Make
  `stage4_source_policy.cmake` require the future Stage 4 allowlist and reject `Python.h`,
  `pybind`, `Cantera::` in public headers, BOFFIN/COAST paths and a second plugin ABI.
- [ ] **Step 3: Run RED.** Run
  `ctest --test-dir build/stage4-debug -R stage4_source_policy --output-on-failure` after
  configuring tests; expect failure because the Stage 4 ledger/allowlist is absent.
- [ ] **Step 4: Write only governance data.** Add exact accepted hashes, frozen API table,
  execution preconditions and the default serial rule to the ledger/AGENTS amendment.
- [ ] **Step 5: Verify.** Run `git diff --check` and the policy fixture; expect PASS. Do not
  run numerical tests.
- [ ] **Step 6: Main-agent commit.** Commit `docs: freeze Stage 4 accepted baseline` with
  the authorized DCO trailer.

### Task 4F-1: Composition and Total-Thermochemical-Enthalpy Identities

**Depends on:** 4F-0.

**Files:**
- Create: `include/hundun/chem_composition.hpp`
- Create: `include/hundun/chem_reports.hpp`
- Create: `tests/unit/test_chem_composition.cpp`
- Create: `tests/unit/test_chem_composition_header_contract.cpp`
- Modify: `src/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:** Produces the value types listed in section 2 and
`validate_composition_identity(const CompositionIdentity&)`.

- [ ] **Step 1: Write RED.** Test two species/two elements, all-species storage, exact
  species ordering, positive molecular weights, element matrix shape, fingerprint
  sensitivity, `h_tc` identity and stable `ChemistryStatus` values. Include a mutation
  fixture that drops the last species.
- [ ] **Step 2: Run RED.** Configure Debug and run
  `ctest --test-dir build/stage4-debug -R 'chem_composition' --output-on-failure`;
  expect compile failure because the header is absent.
- [ ] **Step 3: Implement minimal value types.** Use `std::vector`, checked finite/positive
  validation and a deterministic fingerprint over names, molecular weights and element
  counts. Do not add a dependent-species elimination path.
- [ ] **Step 4: Kill mutations.** Verify tests fail if species order is sorted, the last
  `rhoY` is inferred, formation enthalpy is excluded, or fingerprint omits element counts.
- [ ] **Step 5: Run gate.** Run focused unit/header tests and `git diff --check`.
- [ ] **Step 6: Main-agent commit.** Commit `feat: define reacting composition identities`.

### Task 4F-2: Internal Services and Analytic Backend

**Depends on:** 4F-1.

**Files:**
- Create: `src/chem_thermodynamics_service_detail.hpp`
- Create: `src/chem_backend_detail.hpp`
- Create: `tests/support/chem_analytic_backend.hpp`
- Create: `tests/support/chem_analytic_backend.cpp`
- Create: `tests/unit/test_chem_services.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:** Produces product-private `ThermodynamicsService`,
`TransportPropertyService` and `ChemistryBackend` contracts. Test support produces
`make_analytic_reacting_backend_for_tests(...)` using constant `cp`, two species and an
exactly integrable reversible source; it never enters a product target or installation.

- [ ] **Step 1: Write RED.** Add tests for `(p0,h,Y)->T`, identity preservation,
  scalar/batch-independent ordering, integrated chemistry delta, virtual destructor and
  backend-composition mismatch rejection. Compile two independent analytic backend classes
  against the same services to prove no Cantera-specific construction leaks into callers.
- [ ] **Step 2: Run RED.** Run `ctest -R chem_services`; expect missing internal service
  types.
- [ ] **Step 3: Implement minimal services.** Keep all interfaces private, return value
  reports, reject non-finite/non-positive states and never expose a backend-owned pointer.
- [ ] **Step 4: Mutation check.** Swap species order, change `h_tc` to `cp*T`, return an
  endpoint rate instead of an integral, expose a Cantera concrete type and remove the
  virtual destructor; each mutation must fail.
- [ ] **Step 5: Run focused Debug tests.** Use `ctest --test-dir build/stage4-debug -R
  'chem_services|chem_composition' --output-on-failure`.
- [ ] **Step 6: Main-agent commit.** Commit `feat: add internal reacting services`.

### Task 4F-3: Reacting Attempt and Source Transaction

**Depends on:** 4F-1, accepted Stage 3 transaction API.

**Files:**
- Create: `src/flow_reacting_transaction_detail.hpp`
- Create: `src/flow_reacting_transaction.cpp`
- Create: `tests/unit/test_reacting_transaction.cpp`
- Create: `tests/mpi/test_reacting_rollback.cpp`
- Modify: `src/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:** Produces `ReactingAttemptState`, `ReactingSourceDelta` and
`ReactingSourceTransaction::add_species/add_enthalpy/commit/rollback`. Deltas carry source
identity and units; commit is publish-once.

- [ ] **Step 1: Write RED.** Test all-species source accumulation, zero-sum total species
  chemistry mass, source provenance, duplicate commit rejection and bytewise rollback of
  state/history/`p0`.
- [ ] **Step 2: Run RED.** Expect compile failure for missing transaction.
- [ ] **Step 3: Implement minimal transaction.** Store trial deltas separately, expose
  read-only totals, validate before publication and reuse Stage 3 collective decision.
- [ ] **Step 4: Mutation check.** Directly modify committed arrays, omit one species,
  commit enthalpy chemistry heat or publish before collective agreement; RED must fail.
- [ ] **Step 5: Run unit and `mpiexec -n 2` rollback tests.** No 12^3 case.
- [ ] **Step 6: Main-agent commit.** Commit `feat: add reacting source transaction`.

### Task 4F-4: Schema v4 and Broadcast

**Depends on:** 4F-1, 4P-1 mechanism identity fields may initially use hashes supplied by
the test fixture.

**Files:**
- Create: `include/hundun/cfg_resolved_case_v4.hpp`
- Create: `include/hundun/cfg_resolved_case_v4_loader.hpp`
- Create: `src/cfg_resolved_case_v4_loader.cpp`
- Create: `src/cfg_resolved_case_v4_loader_detail.hpp`
- Create: `src/app_resolved_case_v4_broadcast.cpp`
- Create: `src/app_resolved_case_v4_broadcast_detail.hpp`
- Create: `tests/unit/test_resolved_case_v4.cpp`
- Create: `tests/mpi/test_resolved_case_v4_broadcast.cpp`
- Modify: `src/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:** Produces `ResolvedReactingCaseV4`,
`load_resolved_reacting_case_v4(path)` and
`broadcast_resolved_reacting_case_v4(MPI_Comm,int,const ResolvedReactingCaseV4*)`.

- [ ] **Step 1: Write RED.** Cover mechanism path/hash/phase, species identity, open vs
  closed `p0`, molecular transport, non-catalytic wall, chemistry tolerance and unknown-key
  rejection. Test path escape and rank-consistent error.
- [ ] **Step 2: Run RED.** Expect missing loader and broadcast.
- [ ] **Step 3: Implement rank-0 parse and typed broadcast.** Reuse yyjson and existing
  collective error pattern; do not serialize Cantera types or unordered maps.
- [ ] **Step 4: Mutation check.** Accept unknown keys, reorder species on non-root, ignore
  mechanism hash or allow both adiabatic/isothermal values; tests must fail.
- [ ] **Step 5: Run focused unit plus 1/2-rank tests and standalone header test.**
- [ ] **Step 6: Main-agent commit.** Commit `feat: add reacting schema v4`.

### Task 4F-5: Checkpoint v4 and Diagnostics Provider Contracts

**Depends on:** 4F-1, 4F-3, accepted Stage 3 v3 persistence/diagnostics design.

**Files:**
- Create: `include/hundun/flow_checkpoint_v4.hpp`
- Create: `include/hundun/diag_provider_registry.hpp`
- Create: `src/flow_checkpoint_v4_detail.hpp`
- Create: `src/diag_provider_registry.cpp`
- Create: `tests/unit/test_checkpoint_v4.cpp`
- Create: `tests/unit/test_diag_provider_registry.cpp`
- Modify: `src/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:** Produces versioned section IDs, `CheckpointSectionProvider`,
`DiagnosticProvider`, `CheckpointV4Manifest` and append-only registry APIs. Absence means
unregistered, never a fake-zero record.

- [ ] **Step 1: Write RED.** Test duplicate section/kind rejection, deterministic ordering,
  presence tags, unknown mandatory section rejection, optional-section skip, no mutation
  during diagnostics and no hidden collective in local provider callbacks.
- [ ] **Step 2: Run RED.** Expect missing registry contracts.
- [ ] **Step 3: Implement registries and manifest value types.** Keep v1--v3 readers
  untouched. Validate all providers before publish.
- [ ] **Step 4: Mutation check.** Register duplicate IDs, synthesize absent zero record,
  iterate unordered provider order or let diagnostics write state; RED must fail.
- [ ] **Step 5: Run focused unit/header/tests-off configure.**
- [ ] **Step 6: Main-agent commit.** Commit `feat: define reacting persistence registries`.

### Task 4P-1: Cantera Provenance Lock

**Depends on:** 4F-0. May proceed before product code.

**Files:**
- Create: `third_party/cantera/UPSTREAM.json`
- Create: `third_party/cantera/PATCHES.md`
- Create: `LICENSES/cantera-BSD-3-Clause.txt`
- Create or extend: `THIRD_PARTY_NOTICES`
- Create: `.superpowers/sdd/stage4-4P-1-cantera-provenance-receipt.md`
- Test: `tests/cmake/stage4_source_policy.cmake`

**Interfaces:** Produces exact source archive identity, dependency/license table, patch
ledger and allowed binary paths consumed by 4P-2/3/4.

- [ ] **Step 1: Verify official archive.** Record tag `v3.2.0`, commit
  `4a8358eb80cfeb50474386b5f9ec0b3a83519889` and SHA
  `a94682ef3cb60dc57c8d14fc4cccd94e8f6bb74cab9c3f5465ee90832859360b`.
- [ ] **Step 2: Freeze every transitive dependency.** Record version, URL, source SHA,
  binary SHA, license, ABI flags and whether bundled; do not write a candidate version as
  accepted before its archive is verified.
- [ ] **Step 3: Write policy RED.** Mutate one archive hash, omit one license and insert a
  Cantera file under `src`; the fixture must reject each mutation.
- [ ] **Step 4: Run provenance fixture and `git diff --check`.** No build required.
- [ ] **Step 5: Main-agent legal review.** Confirm no endorsement language and mechanisms
  are not covered by Cantera license.
- [ ] **Step 6: Main-agent commit.** Commit `docs: lock Cantera provenance`.

### Task 4P-2: Prebuilt Cantera Artifact Pipeline and Network-Independent Consumer

**Depends on:** 4P-1.

**Files:**
- Create: `cmake/build_bundled_cantera_linux_cpu.sh`
- Create: `cmake/HundunCanteraPackage.cmake`
- Create: `third_party/cantera/PREBUILT-LINUX-X86_64.json`
- Create: `tests/cmake/stage4_cantera_builder_contract.cmake`
- Modify: `CMakeLists.txt`, `cmake/HundunProvenanceGuard.cmake`

**Interfaces:** Produces CMake cache variables `HUNDUN_ENABLE_CANTERA` and
`HUNDUN_CANTERA_PACKAGE_ROOT`, a maintainer-only pinned Linux artifact builder, and normal-build
imported target `hundun_third_party_cantera` plus explicitly named `hundun_third_party_*`
transitive targets.
Normal configure accepts no URL/download command and never invokes the artifact builder.

- [ ] **Step 1: Write RED.** Configure a fixture with an empty local package root; require a
  precise missing-artifact error. Add mutations containing `FetchContent`, a network URL,
  or an invocation of the maintainer builder from normal CMake; policy must reject each.
- [ ] **Step 2: Run RED.** `cmake -P tests/cmake/stage4_cantera_builder_contract.cmake` must
  fail before the module exists.
- [ ] **Step 3: Implement the maintainer artifact producer.** It consumes only locally
  pinned source archives inside a frozen Ubuntu 22.04/GCC 11 builder root and emits shared
  libraries, data, licenses and a complete build-type/ABI/ISA/hash manifest. The official
  artifact is Release. If upstream Cantera requires
  Python/SCons, they are confined to this producer environment and never enter a normal
  HUNDUN configure/build/install/runtime/formal-acceptance path.
- [ ] **Step 4: Implement imported-target discovery.** Validate headers, shared libraries,
  data directory, transitive dependencies, ABI/ISA manifest and SHA files; never invoke the
  producer, SCons or Python from user CMake.
- [ ] **Step 5: Audit network independence.** Use producer/configure traces and fetch-string
  scans in isolated child processes, then run tests-off configure with the verified local
  package. Do not disconnect host networking.
- [ ] **Step 6: Main-agent package/provenance review and commit.** Confirm that fallback B/C
  was not silently selected; either accept A or stop for a design amendment. Commit
  `build: add prebuilt Cantera artifact pipeline`.

### Task 4P-3: C++ Link, Thread and MPI Spike

**Depends on:** 4P-2.

**Files:**
- Create: `tests/integration/cantera_cpp_link_spike.cpp`
- Create: `tests/mpi/test_cantera_workspace_isolation.cpp`
- Modify: `tests/CMakeLists.txt`

**Interfaces:** Proves one runtime per rank and one complete workspace per thread without
adding product adapter code.

- [ ] **Step 1: Write the spike.** Load one tiny mechanism whose independent source,
  redistribution permission and SHA-256 were recorded by 4P-1, or resolve an external
  hash-verified test asset. Construct two independent thread workspaces, evaluate
  thermo/transport/reaction and print only numeric hashes.
- [ ] **Step 2: Run Release link RED.** Expect missing imported target before 4P-2 artifact
  root is supplied; then configure with exact package root.
- [ ] **Step 3: Add thread/MPI assertions.** Within each rank, verify every active thread
  owns distinct mutable objects; never compare virtual pointer values between processes.
  Verify numeric results are deterministic for fixed lane order and 1/2-rank hashes match.
- [ ] **Step 4: Run `ctest -R cantera_workspace_isolation` and `ldd` in Release and ordinary
  GCC Debug.** Reject Python, Conda, `_GLIBCXX_DEBUG`, ABI=0, exceptions/RTTI mismatch,
  Clang/libc++ mixing, unexpected libstdc++ or absolute build-tree RPATH.
- [ ] **Step 5: Run focused ASan/UBSan single-rank spike.** No large chemistry case.
- [ ] **Step 6: Main-agent commit.** Commit `test: prove Cantera C++ workspace isolation`.

### Task 4P-4: Relocatable Linux Package Prototype

**Depends on:** 4P-2, 4P-3.

**Files:**
- Create: `cmake/HundunRelocatablePackage.cmake`
- Create: `tests/cmake/stage4_package_contract.cmake`
- Modify: `CMakeLists.txt`, `src/CMakeLists.txt`, `THIRD_PARTY_NOTICES`

**Interfaces:** Produces an install tree with relative RPATH, bundled shared libraries,
Cantera data and ABI manifest; no product release is published.

- [ ] **Step 1: Write RED.** Install to a temporary prefix, move the prefix, run a tiny C++
  backend executable and assert no original build path appears in `readelf -d`, `ldd` or
  tracked text.
- [ ] **Step 2: Run RED.** Expect relocation failure before package rules exist.
- [ ] **Step 3: Implement install/RPATH rules.** Use bundled shared libraries as the v1
  default and relative RPATH. Bundle only audited ordinary dependencies; never bundle
  glibc, a compiler runtime outside the frozen policy or unlicensed mechanisms. Static
  linking is rejected unless a later design amendment audits all licenses and ABI effects.
- [ ] **Step 4: Run clean Ubuntu-profile package smoke.** `hundun --version`, validate a
  case and load backend without Python.
- [ ] **Step 5: Record package file/binary hashes and size.** Size is informational. Test an
  atomic whole-bundle upgrade and reject replacing only one dependency `.so` under an old
  ABI manifest.
- [ ] **Step 6: Main-agent commit.** Commit `build: prototype relocatable reacting package`.

### Task 4C-1: Cantera Runtime and Per-Thread Workspace

**Depends on:** 4F-2, 4P-3.

**Files:**
- Create: `src/chem_workspace_detail.hpp`
- Create: `src/chem_cantera_backend.cpp`
- Create: `tests/unit/test_cantera_backend.cpp`
- Modify: `src/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:** Produces `CanteraBackendRuntime`, `CanteraWorkspacePool` and
`make_cantera_backend(const ResolvedReactingCaseV4&, CanteraWorkspacePool&)`.

- [ ] **Step 1: Write RED.** Assert exact mechanism/phase/composition identity, per-thread
  object independence, deterministic lane assignment and explicit missing-workspace error.
- [ ] **Step 2: Run RED.** Expect missing product backend.
- [ ] **Step 3: Implement runtime/workspace ownership.** Cache immutable identity only;
  each thread owns Solution/Thermo/Kinetics/Transport/Reactor/Integrator.
- [ ] **Step 4: Mutation check.** Share one Solution, use thread ID as species order or let a
  workspace outlive runtime; focused tests/sanitizers must fail.
- [ ] **Step 5: Run unit plus 1/2-rank focused tests.**
- [ ] **Step 6: Main-agent commit.** Commit `feat: add Cantera backend runtime`.

### Task 4C-2: Cantera Thermodynamics Adapter

**Depends on:** 4C-1.

**Files:**
- Modify: `src/chem_cantera_backend.cpp`
- Create: `tests/unit/test_cantera_thermodynamics.cpp`

**Interfaces:** Implements `ThermodynamicsService::evaluate` for `(p0,h_tc,Y)->T,rho,cp,W`.

- [ ] **Step 1: Write RED.** Use two synthetic states and compare adapter results to direct
  C++ Cantera calls made only inside the test. Cover non-normalized/negative Y rejection,
  out-of-range inversion and exact species order.
- [ ] **Step 2: Run RED.** Expect unimplemented evaluation.
- [ ] **Step 3: Implement state setting/inversion.** Preserve input `p0`/`h_tc`; never use
  mechanical pressure `pi` and never normalize a materially invalid composition silently.
- [ ] **Step 4: Mutation check.** Feed `p0+pi`, use sensible enthalpy, reorder species or
  accept negative Y; tests must fail.
- [ ] **Step 5: Run focused Debug/Release unit tests.**
- [ ] **Step 6: Main-agent commit.** Commit `feat: adapt Cantera thermodynamics`.

### Task 4C-3: Cantera Transport Adapter

**Depends on:** 4C-2.

**Files:**
- Modify: `src/chem_cantera_backend.cpp`
- Create: `tests/unit/test_cantera_transport.cpp`

**Interfaces:** Implements `TransportPropertyService::evaluate` returning molecular
mixture-averaged coefficients in SI units.

- [ ] **Step 1: Write RED.** Compare viscosity, conductivity and every species diffusivity
  against direct C++ reference at two states; verify size/order/units and finite positivity.
- [ ] **Step 2: Run RED.** Expect missing transport report.
- [ ] **Step 3: Implement minimal adapter.** Convert units once at the boundary and retain
  exact composition ordering.
- [ ] **Step 4: Mutation check.** Use thermal diffusivity as conductivity, convert cm2/s
  incorrectly or omit the last species; tests must fail.
- [ ] **Step 5: Run focused unit tests and thread-isolation test.**
- [ ] **Step 6: Main-agent commit.** Commit `feat: adapt Cantera gas transport`.

### Task 4C-4: Stiff Chemistry Interval Adapter

**Depends on:** 4C-2.

**Files:**
- Modify: `src/chem_cantera_backend.cpp`
- Create: `tests/unit/test_cantera_chemistry_interval.cpp`

**Interfaces:** Implements `ChemistryBackend::integrate`; the interval keeps Eulerian cell
mass and `h_tc` fixed and returns integrated `rhoY` deltas plus final `Y,T` identity.

- [ ] **Step 1: Write RED.** Reaction-off identity, one finite-rate interval, two half-step
  vs one full interval tolerance, element/mass/enthalpy budgets, absolute stage time and
  every backend-neutral solver-failure status. Require completed duration/internal-step
  diagnostics without exposing a Cantera exception or error type.
- [ ] **Step 2: Run RED.** Expect unimplemented interval.
- [ ] **Step 3: Implement adapter ODE contract.** Do not select a Cantera reactor merely by
  name; set/recover state so HUNDUN cell `rho` and `h_tc` invariants hold. Reinitialize on
  every restored attempt. Map all Cantera failures once to `ChemistryStatus`; retain raw
  backend text only in an internal log record.
- [ ] **Step 4: Mutation check.** Change density, add heat source, integrate `dt` instead of
  `dt/2`, reuse failed ReactorNet, expose partial failed deltas or return endpoint rates;
  RED must fail.
- [ ] **Step 5: Run focused unit and small sanitizer tests.**
- [ ] **Step 6: Main-agent commit.** Commit `feat: integrate stiff cell chemistry`.

### Task 4C-5: 0D/PSR Conformance and Backend Milestone

**Depends on:** 4C-2, 4C-3, 4C-4.

**Files:**
- Create: `tests/numerical/test_reacting_zero_dimensional.cpp`
- Create: `tests/numerical/test_reacting_psr.cpp`
- Create: `.superpowers/sdd/stage4-4C-backend-milestone.md`
- Modify: `tests/CMakeLists.txt`

**Interfaces:** Produces accepted scalar/batch/backend conformance evidence consumed by
Stage 5; no new product interface.

- [ ] **Step 1: Add reaction-only 0D oracle.** Compare HUNDUN backend interval against a
  direct C++ Cantera reference for species, temperature and integrated source.
- [ ] **Step 2: Add PSR fixture.** Use a small deterministic residence-time/source problem;
  compare steady residual and failure classification, not a long flame case.
- [ ] **Step 3: Add backend-spy test.** Verify thermo, transport, chemistry and PSR report
  the same mechanism/composition fingerprint.
- [ ] **Step 4: Run Debug and focused Release.** No MPI or long sweep.
- [ ] **Step 5: Main-agent complete-diff review of 4C.** Check Cantera types remain private,
  workspace ownership and licensing.
- [ ] **Step 6: Commit.** Commit `test: accept Cantera backend conformance`.

### Task 4R-0: Coupled C-T-C and Two-PISO Proof

**Depends on:** 4F-3, 4C-5, accepted Stage 3 flow hooks.

**Files:**
- Create: `src/flow_reacting_coupling_detail.hpp`
- Create: `tests/unit/test_reacting_coupling.cpp`
- Create: `tests/numerical/test_reacting_mms.cpp`
- Create: `docs/numerics/stage4-coupling-proof.md`

**Interfaces:** Freezes `ReactingCouplingSchedule` and
`ReactingStepReport{chemistry_call_count=2, pressure_corrector_count=2,...}`. No full
driver is implemented yet.

- [ ] **Step 1: Write non-commuting 2x2 ADR RED.** Require local slope 3 and global slope 2
  for `C/2-T-C/2`; kill `C-T`, `T-C`, wrong half duration and wrong stage time.
- [ ] **Step 2: Write operator spy.** Require two chemistry calls of `dt/2`, one scalar
  transport, exactly two PISO and PISO #2 consuming post-chemistry-2 state.
- [ ] **Step 3: Compare fixed-two-PISO candidates.** Prove either midpoint predictor flux
  is sufficient or a conservative predictor-to-final delta flux is required inside PISO #2.
- [ ] **Step 4: Implement only the selected schedule value object/hook contract.** Do not
  add a production step yet.
- [ ] **Step 5: Run analytic/MMS RED and document rejected candidate.**
- [ ] **Step 6: Main-agent commit.** Commit `docs: freeze reacting coupling schedule`.

### Task 4R-1: Reacting State and Initialization

**Depends on:** 4F-1, 4F-4, 4C-2, 4R-0.

**Files:**
- Create: `include/hundun/flow_reacting_state.hpp`
- Create: `src/flow_reacting_state.cpp`
- Create: `tests/unit/test_reacting_state.cpp`
- Modify: `src/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:** Produces `ReactingFlowState` with all `rhoY`, `rho_h_tc`, derived-cache
epoch and committed/trial/history lifecycle.

- [ ] **Step 1: Write RED.** Test all-species allocation, initialization from `p0,T,Y`,
  cache invalidation, non-positive state rejection and no reinterpretation of Stage 2 `h`.
- [ ] **Step 2: Run RED.** Expect missing state.
- [ ] **Step 3: Implement state through FieldRegistry.** Store conservative fields only;
  derived properties are attempt-local and invalidated on every writer epoch.
- [ ] **Step 4: Mutation check.** Infer last species, store temperature as authority, reuse
  stale cache or alias Stage 2 enthalpy descriptor; RED must fail.
- [ ] **Step 5: Run unit/header tests.**
- [ ] **Step 6: Main-agent commit.** Commit `feat: add reacting conservative state`.

### Task 4R-2: Conservative Species and Enthalpy Transport

**Depends on:** 4R-1, 4C-3.

**Files:**
- Create: `src/flow_reacting_transport_detail.hpp`
- Create: `src/flow_reacting_transport.cpp`
- Create: `tests/unit/test_reacting_transport.cpp`
- Create: `tests/mpi/test_reacting_decomposition.cpp`
- Modify: `src/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:** Produces `assemble_reacting_transport(...)` consuming the accepted
`FaceMassFlux`, molecular transport report and optional WALE diffusivity.

- [ ] **Step 1: Write RED.** Manufactured advection/diffusion for every species and
  `h_tc`, exact `sum(F_k)=0`, element/global conservation and final-flux identity.
- [ ] **Step 2: Run RED.** Expect missing operator.
- [ ] **Step 3: Implement MUSCL/central diffusion on existing product path.** Apply
  correction velocity once; do not maintain a reacting-only face flux.
- [ ] **Step 4: Mutation check.** Omit correction velocity, use predictor flux, transport
  only `Ns-1`, independently clip species or use stale density; RED must fail.
- [ ] **Step 5: Run focused unit and 1/2-rank manufactured tests.**
- [ ] **Step 6: Main-agent commit.** Commit `feat: transport reacting scalars conservatively`.

### Task 4R-3: Open-Domain C-T-C Reacting Step

**Depends on:** 4F-3, 4R-0, 4R-1, 4R-2, 4C-4.

**Files:**
- Create: `src/flow_reacting_coupling.cpp`
- Extend: `tests/unit/test_reacting_coupling.cpp`
- Extend: `tests/numerical/test_reacting_mms.cpp`

**Interfaces:** Produces `attempt_open_reacting_step(...)` and complete
`ReactingStepReport`; `p0` remains constant.

- [ ] **Step 1: Extend RED.** Cover chemistry/transport degenerate limits, full C-T-C,
  two-PISO spy, integrated source/divergence identity and failure injection at all stages.
- [ ] **Step 2: Run RED.** Expect missing production step.
- [ ] **Step 3: Implement selected schedule exactly.** Use one transaction; PISO #2 sees
  post-C2 thermo and final source; commit only after collective validation.
- [ ] **Step 4: Mutation check.** Write intermediate to BDF history, update `p0`, reuse C1
  on retry or add a third PISO; RED must fail.
- [ ] **Step 5: Run unit, 1/2-rank rollback and <=12^3 open smoke.**
- [ ] **Step 6: Main-agent commit.** Commit `feat: advance open reacting flow`.

### Task 4R-4: Body-Fitted Reacting Boundary Conditions

**Depends on:** 4R-2, 4R-3.

**Files:**
- Create: `src/flow_reacting_boundary_detail.hpp`
- Create: `src/flow_reacting_boundary.cpp`
- Create: `tests/unit/test_reacting_boundary.cpp`
- Modify: `src/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:** Produces non-catalytic impermeable species wall, adiabatic/isothermal
enthalpy wall, typed inlet state and pressure-outlet thermo policy.

- [ ] **Step 1: Write RED.** Wall species flux zero, adiabatic heat flux zero,
  isothermal `(p0,T,Y)->h_tc`, inlet redundant-state consistency and unsupported catalytic
  key rejection.
- [ ] **Step 2: Run RED.** Expect missing boundary implementation.
- [ ] **Step 3: Implement through existing BoundaryPatch authority.** Do not create a
  chemistry-specific patch registry.
- [ ] **Step 4: Mutation check.** Use `pi` for thermo pressure, allow species penetration,
  accept inconsistent T/h or synthesize catalytic rates; RED must fail.
- [ ] **Step 5: Run focused unit and one body-fitted 12^3 smoke.**
- [ ] **Step 6: Main-agent commit.** Commit `feat: add reacting boundary conditions`.

### Task 4R-5: IBM and WALE Reacting Coupling

**Depends on:** 4R-3, 4R-4, accepted Stage 3 IBM/WALE authority.

**Files:**
- Create: `src/flow_reacting_immersed_detail.hpp`
- Create: `src/flow_reacting_immersed.cpp`
- Create: `tests/mpi/test_reacting_immersed.cpp`
- Modify: `src/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:** Produces reacting scalar wall rows and effective diffusivity using the same
Stage 3 geometry, pressure, final flux and WALE coefficient identities.

- [ ] **Step 1: Write RED.** Non-catalytic IBM species condition, adiabatic/isothermal wall,
  shared pressure authority, WALE evaluated once per trial and molecular-only bitwise path.
- [ ] **Step 2: Run RED.** Expect missing reacting IBM adapter.
- [ ] **Step 3: Implement thin composition layer.** Reuse Stage 3 row reconstruction and
  wall ownership; do not copy IBM geometry or force code.
- [ ] **Step 4: Mutation check.** Recompute wall pressure, use different flux, evaluate WALE
  per species or let wall source bypass transaction; RED must fail.
- [ ] **Step 5: Run focused 1/2-rank and one <=12^3 IBM smoke.**
- [ ] **Step 6: Main-agent commit.** Commit `feat: couple reacting flow to IBM and WALE`.

### Task 4R-6: Closed and Partially-Closed `p0`

**Depends on:** 4R-3, 4R-4.

**Files:**
- Create: `src/flow_reacting_pressure_constraint_detail.hpp`
- Create: `src/flow_reacting_pressure_constraint.cpp`
- Create: `tests/mpi/test_reacting_closed_pressure.cpp`
- Modify: `src/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:** Produces global `p0` predictor/corrector and closed-domain divergence source
from one MPI reduction; `p0` is committed/rollback/checkpoint state.

- [ ] **Step 1: Write RED.** Uniform closed reactor against independent constant-volume
  thermo pressure, zero boundary net flux, partially closed integral identity, 1/2-rank
  collective equality and rollback.
- [ ] **Step 2: Run RED.** Expect missing constraint implementation.
- [ ] **Step 3: Implement global scalar update.** Use the same integrated chemistry and
  transport source as scalar update; never update `p0` cellwise.
- [ ] **Step 4: Mutation check.** Omit `dp0/dt`, use endpoint rate, update open-domain `p0`
  or commit before global agreement; RED must fail.
- [ ] **Step 5: Run focused unit/MPI and small closed-box test.**
- [ ] **Step 6: Main-agent commit.** Commit `feat: add closed reacting pressure constraint`.

### Task 4R-7: Reacting Driver Combinations

**Depends on:** 4F-4, 4R-3--4R-6.

**Files:**
- Create: `src/app_reacting_flow_driver.cpp`
- Create: `src/app_reacting_flow_driver_detail.hpp`
- Create: `tests/unit/test_reacting_driver.cpp`
- Modify: integration-owned root dispatch files and `src/CMakeLists.txt`

**Interfaces:** Produces `run_reacting_flow_case(const ResolvedReactingCaseV4&, MPI_Comm)`;
the root `hundun` executable dispatches schema v4 without changing v1--v3 behavior.

- [ ] **Step 1: Write RED.** Spy legal combinations: open/closed, body-fitted/IBM,
  none/WALE and molecular transport. Test illegal catalytic/missing mechanism/mismatched
  species combinations before allocation.
- [ ] **Step 2: Run RED.** Expect missing driver.
- [ ] **Step 3: Implement composition root.** Construct services/runtime/workspaces once,
  honor fixed operator order and emit structured reports.
- [ ] **Step 4: Mutation check.** Duplicate backend, bypass resolved case, change old
  dispatch or silently downgrade WALE/IBM; RED must fail.
- [ ] **Step 5: Run driver unit, `--validate`, `--print-resolved` and one-step 12^3 smoke.**
- [ ] **Step 6: Main-agent commit.** Commit `feat: add reacting flow driver`.

### Task 4A-1: Checkpoint v4 Reacting Sections

**Depends on:** 4F-5, 4R-1, 4R-6, 4R-7.

**Files:**
- Create: `src/flow_checkpoint_v4.cpp`
- Extend: `src/flow_checkpoint_v4_detail.hpp`
- Extend: `tests/unit/test_checkpoint_v4.cpp`
- Create: `tests/mpi/test_checkpoint_v4_mpi.cpp`

**Interfaces:** Registers sections for composition identity, all `rhoY`, `rho_h_tc`, `p0`,
BDF/time history and backend/mechanism identity. Workspaces/caches are excluded.

- [ ] **Step 1: Write RED.** continuous-vs-restart, presence combinations, CRC corruption,
  mechanism mismatch, failed-read state unchanged and same-partition 1/2-rank continuation.
- [ ] **Step 2: Run RED.** Expect missing codecs.
- [ ] **Step 3: Implement validate-then-publish codecs.** Use explicit lengths/units and
  stable section order; do not modify v2/v3 readers.
- [ ] **Step 4: Mutation check.** Omit last species, persist cache, ignore mechanism hash or
  partially publish failed read; RED must fail.
- [ ] **Step 5: Run focused unit/MPI and Restart smoke.**
- [ ] **Step 6: Main-agent commit.** Commit `feat: checkpoint reacting state v4`.

### Task 4A-2: Reacting Diagnostics and Counters

**Depends on:** 4F-5, 4R-7.

**Files:**
- Create: `include/hundun/diag_reacting.hpp`
- Create: `src/diag_reacting.cpp`
- Create: `tests/unit/test_reacting_diagnostics.cpp`
- Modify: integration-owned diagnostics registry

**Interfaces:** Registers providers for composition/mechanism identity, final residuals,
EOS drift, element/species/enthalpy budgets, chemistry calls/failures, workspace counts and
rollback classes.

- [ ] **Step 1: Write RED.** Stable IDs/units, absence means unregistered, local snapshot
  is read-only, exact chemistry-call and two-PISO counters, no hidden collective.
- [ ] **Step 2: Run RED.** Expect missing providers.
- [ ] **Step 3: Implement adapters over final authority reports.** Do not recompute thermo,
  source or residual in diagnostics.
- [ ] **Step 4: Mutation check.** Fake zero provider, endpoint source, mutable callback,
  unstable unordered key order or missing failure class; RED must fail.
- [ ] **Step 5: Run focused unit plus 1/2-rank diagnostic session.**
- [ ] **Step 6: Main-agent commit.** Commit `feat: add reacting diagnostics`.

### Task 4A-3: Low-Cost Stage 4 Acceptance and Optional Detached Evidence

**Depends on:** 4P-4, 4C-5, 4R-7, 4A-1, 4A-2.

**Files:**
- Create: `tests/acceptance/stage4_acceptance.sh`
- Create: `tests/support/run_stage4_detached.sh`
- Create: `.superpowers/sdd/stage4-4A-3-acceptance-manifest.md`
- Modify: `VERSION`
- Modify: `docs/numerics/stage4-capability-ledger.md`

**Interfaces:** Produces frozen Stage 4 code candidate `C4` at version `0.3.0` and its
complete low-cost evidence manifest. Detached 24^3 screens are optional diagnostics and
never block subsequent development; 48^3 waits for v1 final gate.

- [ ] **Step 1: Finalize and freeze candidate.** Set the single version authority to
  `0.3.0`, finish product/test source and selector registration, run a low-cost banner/
  configure preflight, then create a signed candidate commit `C4`. Record its HEAD/tree,
  dirty state and binaries; do not modify product or test source after this point.
- [ ] **Step 2: Freeze selectors.** Include 0D/PSR, C-T-C spy/MMS, 12^3 open/closed,
  IBM+WALE smoke, small 1/2/4-rank, rollback, Checkpoint, diagnostics, package/RPATH,
  full affected Debug and focused Release/sanitizers.
- [ ] **Step 3: Exclude forbidden work.** Assert no 48^3/96^3, Flame D, TPDF, spray or
  sanitizer-large MPI selector is registered.
- [ ] **Step 4: Run exact low-cost matrix on `C4`.** Capture commands, binaries, exit codes and log
  hashes. One M resource group at a time.
- [ ] **Step 5: Optionally launch a 24^3 detached screen.** Bind exact HEAD/binary/env; do
  not wait before governance review and do not treat it as acceptance evidence unless same
  tree/config/hash.
- [ ] **Step 6: Main-agent review and governance receipt.** Review the full Stage 4
  product diff/provenance/callers once. If a correctness fix changes `C4`, create a new
  candidate and rerun only invalidated clusters. Otherwise commit only the evidence
  manifest/ledger as `test: record Stage 4 acceptance evidence`, with
  `accepted_code_head=C4`.

### Task 4A-4: Stage 4 Exact-HEAD Seal

**Depends on:** 4A-3.

**Files:**
- Create: `.superpowers/sdd/stage4-final-acceptance-report.md`
- Modify: `docs/numerics/stage4-capability-ledger.md`
- Modify: `docs/numerics/stage4-6-capability-root.md`
- Modify: `AGENTS.md`

**Interfaces:** Produces governance seal `G4` for tested code candidate `C4` at `0.3.0`;
no default product-repository projection.

- [ ] **Step 1: Recover the tested candidate identity.** Record `C4` HEAD/parent/tree/diff
  SHA, worktree status, build/package/binary/log hashes, toolchain, MPI and background
  process state. Assert all commits after `C4` are governance-only.
- [ ] **Step 2: Verify DCO/provenance.** Check every Stage 4 commit and third-party notice;
  do not add or repair sign-off without authorization.
- [ ] **Step 3: Verify capability limits.** State no TPDF/TCR/spray, no Python, no generic
  Linux claim and no unverified mechanism redistribution.
- [ ] **Step 4: Record governance receipt version `0.3.0`.** Do not modify `VERSION`,
  product code or tests, and do not project to product repo unless the user separately
  authorizes an intermediate projection.
- [ ] **Step 5: Run final read-only verification.** `git diff --check`, acceptance manifest
  hash audit, `ldd`, RPATH, tests-off and worktree/background-process checks.
- [ ] **Step 6: Main-agent decision and governance commit.** Write `STAGE4_ACCEPT` or
  `REJECT`; commit `docs: seal Stage 4 reacting flow` only for an actually accepted `C4`.
  The report must state `accepted_code_head=C4`; the governance commit is not relabeled as
  the numerically tested code HEAD.

## 4. Critical Path and Stage Gate

```text
4F-0 -> 4F-1 -> 4F-2/3/4/5
4P-1 -> 4P-2 -> 4P-3 -> 4P-4
4F-2 + 4P-3 -> 4C-1 -> 4C-2 -> 4C-3/4 -> 4C-5
4F-3 + 4C-5 -> 4R-0 -> 4R-1 -> 4R-2 -> 4R-3
4R-3 -> 4R-4/5/6 -> 4R-7 -> 4A-1/2 -> 4A-3 -> 4A-4
```

默认串行执行。`4F` 或 `4C` 节点完成时，主 agent可以向用户建议是否提前启动
Stage 5/6 的纯模块分支；没有明确批准则继续串行完成 Stage 4。

## 5. Task Reference Matrix

URL、DOI、license class 和 rejected architecture 以
`docs/references/2026-08-09-hundun-flow-stage4-6-reference-catalog.md` 为准。表中“参考点”
只授权数学、接口或测试思想；除登记为 third-party 的完整 Cantera component 外，不授权
复制上游源码。

| Task | Reference | Exact point to reuse | Explicitly avoid |
|---|---|---|---|
| 4F-0 | accepted Stage 3 ledgers；MPI standard | frozen authority/ID inventory and collective evidence | inferring acceptance from a dirty tree |
| 4F-1 | Day--Bell；Cantera thermo docs | all-species identity and total thermochemical enthalpy | dependent-species storage or sensible-only enthalpy |
| 4F-2 | Cantera C++ API；PeleLMeX service layering；NASA/CHEMKIN/SUNDIALS future route | narrow backend-neutral thermo/transport/chemistry services | third-party types or Cantera assumptions in service contracts |
| 4F-3 | Day--Bell integrated source；Stage 3 transaction | integrated deltas and publish-once rollback | endpoint rates or direct field writes |
| 4F-4 | MPI typed broadcast；existing yyjson loader | rank-0 parse, deterministic typed broadcast | backend object serialization |
| 4F-5 | accepted Checkpoint v3/diagnostics registry | append-only provider IDs and validate-then-publish | changing v1--v3 readers or fake-zero absence |
| 4P-1 | Cantera v3.2.0 release/license | archive, revision, patch and transitive-license provenance | assuming Cantera BSD covers mechanisms |
| 4P-2 | Cantera build docs；CMake imported targets | maintainer-produced pinned artifact and local consumer | configure-time fetch or user-build Python |
| 4P-3 | Cantera C++/ReactorNet docs | per-rank runtime and per-thread mutable workspace tests | shared mutable Solution/Reactor objects |
| 4P-4 | CMake RPATH/install/CPack model | relocatable shared-library package | bundled glibc or unaudited static link |
| 4C-1 | Cantera C++ ownership model | immutable runtime identity plus thread-local full workspace | cross-thread mutable caches |
| 4C-2 | Cantera thermo state API | `(p0,h_tc,Y)->T,rho,cp,W` black-box adapter | mechanical pressure or silent invalid-Y normalization |
| 4C-3 | Cantera mixture-averaged transport docs | ordered SI viscosity/conductivity/diffusivity report | multicomponent/Soret scope creep |
| 4C-4 | Cantera ReactorNet；Day--Bell | local stiff interval and integrated composition delta | reactor class name as conservation proof |
| 4C-5 | Cantera 0D/PSR C++ examples | low-cost black-box conformance and failure mapping | long flame or Python oracle |
| 4R-0 | Strang；Day--Bell；Nonaka | symmetric C-T-C proof and final constraint timing | importing SDC/MAC control flow or third PISO |
| 4R-1 | PeleLMeX state model；Cantera thermo | conservative all-species/total-enthalpy authority | temperature authority or Stage 2 `h` reinterpretation |
| 4R-2 | Day--Bell；PeleLMeX transport model | shared final mass flux and correction-velocity species diffusion | reacting-only flux authority |
| 4R-3 | Day--Bell；Nonaka EOS consistency | integrated source, second-half state and final PISO visibility | endpoint heat release source or ad hoc damping |
| 4R-4 | Cantera thermo；PeleLMeX boundary concepts | typed non-catalytic species and thermal states | catalytic chemistry not in v1 |
| 4R-5 | accepted Stage 3 IBM/WALE；PeleLMeX scalar layering | reuse geometry/pressure/WALE/final-flux authorities | copied IBM or per-species WALE |
| 4R-6 | Day--Bell/Nonaka low-Mach pressure split | one global thermodynamic-pressure constraint | cellwise `p0` or mechanical-pressure mixing |
| 4R-7 | HUNDUN composition-root pattern；Cantera service boundary | construct backend once and dispatch one schema | a second executable or direct backend bypass |
| 4A-1 | accepted Checkpoint v3 protocol | stable sections, CRC and validate-then-publish | workspace/cache persistence |
| 4A-2 | accepted diagnostic-provider pattern；PeleLMeX integrated-source idea | reports over final authorities and exact counters | recomputing science in diagnostics |
| 4A-3 | CTest；cross-stage integration plan | hash-bound compact evidence and detached diagnostics | 48/96 or repeated long matrix |
| 4A-4 | cross-stage exact-HEAD protocol；DCO | tested code head separated from governance seal | version/product mutation after testing |

## 6. Task-Level Review and Receipt Template

每个 task 的主 agent receipt 必须写明：

```text
task id and candidate commit
parent/tree/diff SHA-256
files changed and callers
RED and killed mutations
commands/exits/log hashes
API/schema/Restart/diagnostic impact
allocation/MPI/rollback review
reference URL/revision and rejected architecture
license/provenance result
deferred work
ACCEPT or REJECT
```

任务通过后立即结束其 implementation agent。跨模块科学判断、完整 diff 和最终
Stage 4 接受始终由主 agent完成。

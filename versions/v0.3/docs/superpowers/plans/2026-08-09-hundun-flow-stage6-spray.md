# HUNDUN-FLOW Stage 6 Dilute Spray Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:executing-plans` for the default serial path. Use `superpowers:subagent-driven-development` only after the user explicitly authorizes parallel execution at a stage gate. Execute one task at a time and return to the main agent for review.

**Goal:** 交付单组分稀相 point-parcel 喷雾、Schiller--Naumann/Ranz--Marshall/Abramzon--Sirignano 交换、原子双向守恒、MPI migration、static IBM rebound、TAB、双 surrogate 接口、Checkpoint v6 和 diagnostics。

**Architecture:** Parcel 只采样 MeanState 和 Stage 4 thermo/transport services，每个 trial 计算一份积分交换量并通过 SprayCouplingTransaction 同时更新 parcel 与气相。Stage 5 stochastic fields 接收共同 vapor/enthalpy source，不创建 N 套 parcels。Pure-SoA container、staged migration 和 event-aware integrator都参与统一 rollback。

**Tech Stack:** C++17、MPI-3、Stage 4/5 services、counter RNG、CTest、Linux CPU-reference、Apache-2.0/DCO。

## Global Constraints

- 默认串行实施 parent 为 Stage 5 `STAGE5_ACCEPT` exact HEAD。只有用户明确授权并行时，pure parcel/property kernels 可从 Stage 4 service-contract freeze base 开发；正式接受前仍须更新到 Stage 5 accepted head。
- v1 只支持 dilute、spherical、internally uniform-temperature、single-component droplets。
- 产品蒸发路径为 Abramzon--Sirignano；`d^2` law 只作 oracle。
- parcel 只通过 ThermodynamicsService、TransportPropertyService 和 MeanState 查询气相；不得直接调用 Cantera。
- 一个 parcel 同时只有一个 owner rank；migration/injection/breakup 都属于 trial。
- 气相 mass/momentum/species/total-enthalpy gain 与 parcel loss 必须成对提交。
- 每个 gas step 固定一个 parcel predictor 和一个 corrector；内部允许 event-aware substeps；不得增加第三次 PISO。
- ESF common source 对全部 fields 使用同一物理 exchange；不得 per-field 重算 parcels。
- static IBM 只支持 rebound；deposition/film/splash 延期。
- `none` 和 TAB 是 v1 breakup modes；KH--RT、LISA、collision/coalescence 延期。
- n-dodecane/iso-octane 只证明接口通用性，不宣称真实 kerosene/gasoline science。
- COAST `EXEC/Fuels` 正式读取前必须用户确认路径；无再分发许可的机制不进入产品。
- task gate 不运行 48^3/96^3 或大型 sanitizer/MPI；最终唯一 48^3 属于跨阶段计划。
- worker 不决定能量/动量源代数、Stage 5组合时序或完整 rollback；主 agent拥有这些任务。
- `Files` 中出现 integration-owned CMake、root dispatch、registry、`VERSION` 或
  `AGENTS.md` 时，该行由主 agent集成；worker 只返回所需 registration entries，
  不直接修改中央文件。
- Every task's Step 6 commit/receipt/DCO action is main-agent-only. A worker stops after
  Step 5 and returns its diff, commands, outputs and risks without staging or committing.
- 不 push、不发布。

---

## 1. Planned File Map

### Public value/config/report types

```text
include/hundun/spray_reports.hpp
include/hundun/spray_liquid_material.hpp
include/hundun/cfg_resolved_case_v6.hpp
include/hundun/cfg_resolved_case_v6_loader.hpp
include/hundun/flow_checkpoint_v6.hpp
include/hundun/diag_spray.hpp
```

### Internal parcel and physics modules

```text
src/spray_parcel_state_detail.hpp
src/spray_parcel_container.cpp
src/spray_grid_coupling_detail.hpp
src/spray_grid_coupling.cpp
src/spray_trajectory_detail.hpp
src/spray_trajectory.cpp
src/spray_migration_detail.hpp
src/spray_migration.cpp
src/spray_injector_detail.hpp
src/spray_injector.cpp
src/spray_liquid_properties_detail.hpp
src/spray_liquid_properties.cpp
src/spray_gas_film_detail.hpp
src/spray_gas_film.cpp
src/spray_drag_detail.hpp
src/spray_drag.cpp
src/spray_heat_mass_detail.hpp
src/spray_heat_mass.cpp
src/spray_evaporation_detail.hpp
src/spray_evaporation.cpp
src/spray_exchange_integrator_detail.hpp
src/spray_exchange_integrator.cpp
src/spray_coupling_transaction_detail.hpp
src/spray_coupling_transaction.cpp
src/spray_esf_source_detail.hpp
src/spray_esf_source.cpp
src/spray_ibm_impact_detail.hpp
src/spray_ibm_impact.cpp
src/spray_tab_breakup_detail.hpp
src/spray_tab_breakup.cpp
src/flow_spray_coupling_detail.hpp
src/flow_spray_coupling.cpp
src/cfg_resolved_case_v6_loader.cpp
src/cfg_resolved_case_v6_loader_detail.hpp
src/app_resolved_case_v6_broadcast.cpp
src/app_spray_driver.cpp
src/app_spray_driver_detail.hpp
src/flow_checkpoint_v6.cpp
src/flow_checkpoint_v6_detail.hpp
src/diag_spray.cpp
```

### Provenance and tests

```text
.superpowers/oracles/coast-fuels-source-manifest.md
tests/fixtures/liquid/n_dodecane_properties.txt
tests/fixtures/liquid/iso_octane_properties.txt
tests/unit/test_spray_*.cpp
tests/mpi/test_spray_*.cpp
tests/numerical/test_single_droplet.cpp
tests/numerical/test_reacting_spray_smoke.cpp
tests/cmake/stage6_fuel_provenance_contract.cmake
tests/cmake/stage6_source_policy.cmake
tests/acceptance/stage6_acceptance.sh
docs/numerics/stage6-capability-ledger.md
```

Mechanism files only appear in tracked fixtures after 6F-4 confirms redistribution rights.
Otherwise tests resolve them from an external, hash-verified governance asset path.

## 2. Frozen Internal Interfaces

The value/report types below live in `namespace hundun::spray`.

```cpp
namespace hundun::spray {

struct ParcelId {
  std::uint64_t high;
  std::uint64_t low;
};

struct LiquidMaterialIdentity {
  std::string name;
  std::string gas_species_name;
  std::uint64_t fingerprint;
};

struct SprayParcelSnapshot {
  ParcelId id;
  std::array<double, 3> position_m;
  std::array<double, 3> velocity_m_per_s;
  double droplet_mass_kg;
  double droplet_diameter_m;
  double multiplicity;
  double temperature_k;
  std::uint64_t liquid_material_fingerprint;
  std::uint64_t owner_global_cell;
  double age_s;
};

struct LiquidProperties {
  double density_kg_per_m3;
  double cp_j_per_kg_k;
  double latent_heat_j_per_kg;
  double saturation_pressure_pa;
  double surface_tension_n_per_m;
  double viscosity_pa_s;
};

enum class LiquidPropertyStatus : std::uint32_t {
  success,
  unknown_material,
  temperature_out_of_range,
  correlation_domain_error,
  non_finite_output
};

struct LiquidPropertyReport {
  LiquidProperties properties;
  LiquidPropertyStatus status;

  [[nodiscard]] bool succeeded() const noexcept {
    return status == LiquidPropertyStatus::success;
  }
};

struct GasFilmSample {
  std::array<double, 3> velocity_m_per_s;
  double temperature_k;
  double density_kg_per_m3;
  double viscosity_pa_s;
  double conductivity_w_per_m_k;
  double vapor_diffusivity_m2_per_s;
  double vapor_mass_fraction;
  double p0_pa;
};

struct SprayExchangeDelta {
  double liquid_mass_delta_kg;
  std::array<double, 3> parcel_momentum_delta_kg_m_per_s;
  double parcel_energy_delta_j;
  double vapor_mass_to_gas_kg;
  std::array<double, 3> momentum_to_gas_kg_m_per_s;
  double h_tc_to_gas_j;
};

}  // namespace hundun::spray
```

Internal service signature in `src/spray_liquid_properties_detail.hpp`:

```cpp
namespace hundun::spray {

class LiquidPropertyService {
 public:
  virtual ~LiquidPropertyService() = default;
  virtual const LiquidMaterialIdentity& material(
      std::uint64_t fingerprint) const = 0;
  virtual LiquidPropertyReport evaluate(
      std::uint64_t fingerprint,
      double temperature_k) const = 0;
};

}  // namespace hundun::spray
```

For a non-success property report, every numeric member is canonical `+0.0` and must not
be consumed by drag/heat/mass kernels. The status is propagated through collective failure
and diagnostics without silent extrapolation.

`SprayExchangeDelta` signs are fixed by `6F-1`: parcel deltas describe new minus old;
gas fields receive the separately reported positive/negative transfers. Implementers must not
derive one from the other ad hoc in later modules.

### Frozen schema v6 additions

Schema v6 extends accepted v5 and adds one `spray` object:

| JSON path | Type / allowed value | Rule |
|---|---|---|
| `schema_version` | integer `6` | required |
| `spray.model` | `none` or `dilute_point_parcel` | required |
| `spray.coupling` | `one_way` or `two_way` | required when enabled |
| `spray.evaporation` | `none` or `abramzon_sirignano` | required |
| `spray.breakup` | `none` or `tab` | required |
| `spray.liquid.name` | nonempty string | label only, never dispatch key |
| `spray.liquid.fingerprint` | unsigned 64-bit integer | required |
| `spray.liquid.property_file` | confined path | required |
| `spray.liquid.property_sha256` | 64 lowercase hex | required |
| `spray.liquid.gas_species` | exact mechanism species name | required and resolved before allocation |
| `spray.particle_cfl` | finite positive real | required, no hidden default |
| `spray.maximum_substeps` | positive integer | required |
| `spray.ibm_interaction.mode` | `none` or `rebound` | deposition/film/splash rejected |
| `spray.ibm_interaction.normal_restitution` | finite in `[0,1]` | required for rebound |
| `spray.ibm_interaction.tangential_restitution` | finite in `[0,1]` | required for rebound |
| `spray.injectors` | ordered array | injector IDs unique |
| `spray.injectors[*].type` | `point` or `cone` | required |
| `spray.injectors[*].position_m` | three finite reals | required |
| `spray.injectors[*].axis` | normalized three-vector | required for cone |
| `spray.injectors[*].half_angle_rad` | finite in `[0,pi]` | required for cone |
| `spray.injectors[*].mass_flow_kg_per_s` | finite nonnegative real | required |
| `spray.injectors[*].droplet_diameter_m` | finite positive real | required |
| `spray.injectors[*].droplet_temperature_k` | finite positive real | required |
| `spray.injectors[*].parcel_multiplicity` | finite positive real | required |
| `spray.injectors[*].start_time_s` / `end_time_s` | finite ordered interval | required |

Unknown/duplicate keys and film/KH--RT/LISA/collision/multicomponent-liquid settings are
rejected. Property/mechanism hashes, injector order and all coupling flags enter case and
Checkpoint fingerprints. `name` is diagnostic text only；dispatch uses fingerprints and the
resolved gas-species identity.

## 3. Task Sequence

### Common command protocol

Pure parcel/property/correlation tasks may disable Cantera. Driver, evaporation-to-gas and
combined smokes use the accepted Stage 4 package root:

```bash
cmake -S . -B build/stage6-debug \
  -DCMAKE_BUILD_TYPE=Debug \
  -DHUNDUN_BUILD_TESTS=ON \
  -DHUNDUN_ENABLE_CANTERA=ON \
  -DHUNDUN_CANTERA_PACKAGE_ROOT="${HUNDUN_ACCEPTED_CANTERA_ROOT}"
cmake --build build/stage6-debug -j32
ctest --test-dir build/stage6-debug -R "${HUNDUN_TASK_TEST_REGEX}" \
  --output-on-failure -j24
```

The main agent exports `HUNDUN_ACCEPTED_CANTERA_ROOT` from G4 and sets
`HUNDUN_TASK_TEST_REGEX` to the task's exact registered test basename(s). For a pure task,
set `HUNDUN_ENABLE_CANTERA=OFF` and omit the package root. Focused Release uses
`build/stage6-release` and `Release`. COAST fuel assets are resolved only through the 6F-4
manifest；normal tests never guess or scan `EXEC/Fuels`.

### Task 6F-0: Stage 4/5 Intake and Compatibility Freeze

**Depends on:** Default serial path requires Stage 5 exact acceptance.

**Files:**
- Create: `.superpowers/sdd/stage6-6F-0-baseline-receipt.md`
- Create: `docs/numerics/stage6-capability-ledger.md`
- Modify: `AGENTS.md`
- Test: `tests/cmake/stage6_source_policy.cmake`

**Interfaces:** Produces exact parent/tree, Stage 4 service identities, Stage 5 common-source
hook, state/checkpoint/diagnostic IDs and allowed-file table.

- [ ] **Step 1: Record baseline.** Capture accepted Stage 5 HEAD/tree/worktree, service and
  registry signatures, PISO/IBM/WALE/ESF authorities and background processes.
- [ ] **Step 2: Add source-policy RED.** Reject direct Cantera calls, per-field parcels,
  OpenFOAM/COAST source, unlicensed mechanism files, 96^3 and a second mesh/flux authority.
- [ ] **Step 3: Run RED.** Expect failure before the Stage 6 ledger exists.
- [ ] **Step 4: Write governance only.** State no-product-change and default serial gate.
- [ ] **Step 5: Run policy and `git diff --check`.**
- [ ] **Step 6: Main-agent commit.** Commit `docs: freeze Stage 6 accepted baseline`.

### Task 6F-1: Parcel/Gas Conservation and Coupling Proof

**Depends on:** 6F-0. Main-agent task.

**Files:**
- Create: `docs/numerics/stage6-spray-conservation-proof.md`
- Create: `include/hundun/spray_reports.hpp`
- Create: `tests/unit/test_spray_exchange_contract.cpp`

**Interfaces:** Freezes `SprayExchangeDelta`, source signs, multiplicity scaling, vapor
momentum, drag, sensible/latent/thermochemical enthalpy decomposition and final PISO source.

- [ ] **Step 1: Derive balances for one physical droplet and one parcel.** Explicitly
  distinguish single-droplet state from multiplicity-scaled totals.
- [ ] **Step 2: Write algebra RED.** Closed gas+parcel mass/momentum/energy totals, pure
  drag, pure heating, pure evaporation and complete event-to-zero.
- [ ] **Step 3: Freeze vapor injection velocity/enthalpy convention and rejected double-count
  alternatives.**
- [ ] **Step 4: Implement report value types only.** No ODE or source deposition yet.
- [ ] **Step 5: Kill omitted vapor momentum, duplicated latent heat, wrong multiplicity and
  missing divergence-source mutations.**
- [ ] **Step 6: Main-agent commit.** Commit `docs: freeze spray conservation algebra`.

### Task 6F-2: Parcel State, Services and Transaction Contracts

**Depends on:** 6F-1.

**Files:**
- Create: `include/hundun/spray_liquid_material.hpp`
- Create: `src/spray_parcel_state_detail.hpp`
- Create: `src/spray_liquid_properties_detail.hpp`
- Create: `src/spray_coupling_transaction_detail.hpp`
- Create: `tests/unit/test_spray_state_contract.cpp`

**Interfaces:** Produces `ParcelId`, `SprayParcelSnapshot`, `LiquidMaterialIdentity`,
`LiquidProperties`, `LiquidPropertyStatus/Report`, `GasFilmSample`, trial/commit ownership
and service signatures.

- [ ] **Step 1: Write RED.** Positive mass/diameter/multiplicity/T, material fingerprint,
  stable owner, committed/trial semantics and distinct liquid/gas species identities.
- [ ] **Step 2: Run RED.** Expect missing contracts.
- [ ] **Step 3: Implement value/descriptors only.** Container owns arrays; reports own no
  mutable backend pointer.
- [ ] **Step 4: Mutation check.** Use species index as identity, conflate parcel mass with
  droplet mass, let owner be two ranks or commit deletion early; RED must fail.
- [ ] **Step 5: Run header/tests-off compile.**
- [ ] **Step 6: Commit.** Commit `feat: define spray state contracts`.

### Task 6F-3: Schema v6, Checkpoint and Diagnostic Identities

**Depends on:** 6F-2, Stage 5 registries.

**Files:**
- Create: `include/hundun/cfg_resolved_case_v6.hpp`
- Create: `include/hundun/cfg_resolved_case_v6_loader.hpp`
- Create: `include/hundun/flow_checkpoint_v6.hpp`
- Create: `include/hundun/diag_spray.hpp`
- Create: `tests/unit/test_stage6_identity_ledger.cpp`

**Interfaces:** Freezes one-way/two-way, evaporation, IBM rebound, TAB, injector, liquid
material, mechanism mapping, checkpoint sections and diagnostic kinds.

- [ ] **Step 1: Write RED.** Duplicate IDs, unsupported combinations, source-coupling flags,
  property validity, mechanism species mapping and absence semantics.
- [ ] **Step 2: Run RED.** Expect missing v6 identities.
- [ ] **Step 3: Implement append-only definitions.** Preserve v1--v5 IDs.
- [ ] **Step 4: Mutation check.** Enable film/KHRT, accept unknown liquid, reuse v5 section
  or represent absent spray with fake-zero provider; RED must fail.
- [ ] **Step 5: Run focused unit/header tests.**
- [ ] **Step 6: Commit.** Commit `feat: freeze Stage 6 identities`.

### Task 6F-4: COAST Fuels and Dual-Surrogate Provenance Gate

**Depends on:** 6F-0 and explicit user confirmation of COAST root/`EXEC/Fuels` realpath.
The main agent owns source confirmation and the redistribution decision; a worker may
receive only the frozen asset identities and policy fixture.

**Files:**
- Create: `.superpowers/oracles/coast-fuels-source-manifest.md`
- Create: `tests/cmake/stage6_fuel_provenance_contract.cmake`
- Create or externalize: `tests/fixtures/liquid/n_dodecane_properties.txt`
- Create or externalize: `tests/fixtures/liquid/iso_octane_properties.txt`
- Modify: `THIRD_PARTY_NOTICES` only for redistributable assets

**Interfaces:** Produces two liquid property identities and two independent reduced
mechanism identities/hashes/license decisions.

- [ ] **Step 1: Main agent asks the user to confirm exact COAST source and Fuels realpath.** Do not inspect
  a guessed/dirty/old path.
- [ ] **Step 2: Inventory only fuel mechanism/property candidates.** Record format,
  species/phase, size, revision and hash; do not inspect unrelated cases or research data.
- [ ] **Step 3: Select n-dodecane and iso-octane reduced candidates.** Prefer low-cost
  mechanisms actually used by COAST; record any species alias mapping explicitly.
- [ ] **Step 4: Decide redistribution per file.** Unclear license means external
  hash-verified acceptance asset, never tracked product data.
- [ ] **Step 5: Run provenance mutations.** Swapped mechanism, absent species, wrong hash,
  unlicensed tracked file and fuel-name hardcoding must fail.
- [ ] **Step 6: Commit governance/allowed fixtures only.** Commit
  `docs: freeze Stage 6 fuel provenance`.

### Task 6P-1: Pure-SoA Parcel Container and Stable IDs

**Depends on:** 6F-2.

**Files:**
- Create: `src/spray_parcel_container.cpp`
- Create: `tests/unit/test_spray_parcel_container.cpp`
- Modify: `src/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:** Produces `SprayParcelContainer` with component arrays, deterministic ID
order, committed/trial snapshots and staged create/remove.

- [ ] **Step 1: Write RED.** Add/read/update N parcels, ID uniqueness/order, multiplicity,
  staged removal, rollback and no allocation during a reserved attempt.
- [ ] **Step 2: Run RED.** Expect missing container.
- [ ] **Step 3: Implement pure SoA arrays and ID index.** No unordered iteration enters
  numerical order; reserve capacity before trial.
- [ ] **Step 4: Mutation check.** Duplicate ID, AoS-only hidden state, early erase,
  multiplicity default or nondeterministic iteration; RED must fail.
- [ ] **Step 5: Run unit/ASan focused tests.**
- [ ] **Step 6: Commit.** Commit `feat: add transactional parcel container`.

### Task 6P-2: Parcel/Grid Coupling Stencil

**Depends on:** 6P-1, accepted mesh/IBM topology.

**Files:**
- Create: `src/spray_grid_coupling_detail.hpp`
- Create: `src/spray_grid_coupling.cpp`
- Create: `tests/unit/test_spray_grid_coupling.cpp`
- Create: `tests/mpi/test_spray_grid_coupling_mpi.cpp`

**Interfaces:** Produces `ParcelGridCouplingStencil locate_and_build_stencil(position)` used
for both gas interpolation and conservative source deposition.

- [ ] **Step 1: Write RED.** Uniform/linear interpolation, weight sum one, face/edge/corner
  stable tie, global-cell identity and 1/2-rank partition boundary.
- [ ] **Step 2: Run RED.** Expect missing stencil.
- [ ] **Step 3: Implement one paired stencil authority.** Reuse MeshTopology/Geometry;
  forbid separate sampling/deposition stencils.
- [ ] **Step 4: Mutation check.** Weight loss, owner-by-rank tie, different source weights or
  nearest-cell-only hidden fallback; RED must fail.
- [ ] **Step 5: Run unit and 1/2-rank tests.**
- [ ] **Step 6: Commit.** Commit `feat: couple parcels to mesh stencils`.

### Task 6P-3: Event-Aware Parcel Trajectory Integrator

**Depends on:** 6P-1, 6P-2.

**Files:**
- Create: `src/spray_trajectory_detail.hpp`
- Create: `src/spray_trajectory.cpp`
- Create: `tests/unit/test_spray_trajectory.cpp`

**Interfaces:** Produces a generic predictor/corrector over a supplied acceleration callback,
particle-CFL substeps and ordered trajectory segments.

- [ ] **Step 1: Write RED.** Ballistic, constant acceleration, midpoint result, Galilean
  invariance, cell crossing segmentation, dt=0 and event time.
- [ ] **Step 2: Run RED.** Expect missing integrator.
- [ ] **Step 3: Implement one gas-step predictor/corrector with bounded substeps.**
- [ ] **Step 4: Mutation check.** Use old velocity for position, double-advance dt, drop a
  segment or continue after event; RED must fail.
- [ ] **Step 5: Run focused unit/UBSan.**
- [ ] **Step 6: Commit.** Commit `feat: integrate parcel trajectories`.

### Task 6P-4: Staged MPI Migration

**Depends on:** 6P-1--6P-3.

**Files:**
- Create: `src/spray_migration_detail.hpp`
- Create: `src/spray_migration.cpp`
- Create: `tests/mpi/test_spray_migration.cpp`

**Interfaces:** Produces trial `MigrationPlan`, outgoing/incoming buffers and publish/rollback
operations preserving one owner and stable ID.

- [ ] **Step 1: Write RED.** Same-rank move, neighbor move, multi-rank crossing, zero
  parcels, duplicate/drop detection, rollback and deterministic order.
- [ ] **Step 2: Run RED.** Expect missing migration.
- [ ] **Step 3: Implement typed all-rank ownership resolution and staged buffers.** Optimize
  local-neighbor path only after correctness; no rank-changing Restart.
- [ ] **Step 4: Mutation check.** Publish on receive, omit multiplicity/TAB state, duplicate
  owner or rank-local ID; RED must fail.
- [ ] **Step 5: Run 1/2/4-rank small tests.**
- [ ] **Step 6: Commit.** Commit `feat: migrate parcels collectively`.

### Task 6P-5: Deterministic Injector and RNG Domain

**Depends on:** 6P-1, Stage 5 counter RNG service.

**Files:**
- Create: `src/spray_injector_detail.hpp`
- Create: `src/spray_injector.cpp`
- Create: `tests/unit/test_spray_injector.cpp`
- Create: `tests/mpi/test_spray_injector_mpi.cpp`

**Interfaces:** Produces point/cone injector, residual mass state and IDs keyed by injector,
accepted step and ordinal; uses `purpose_tag=spray_injection`.

- [ ] **Step 1: Write RED.** Zero/constant flow, noninteger parcel residual, cone bounds,
  stable ID, retry no duplicate, Restart and rank-owner agreement.
- [ ] **Step 2: Run RED.** Expect missing injector.
- [ ] **Step 3: Implement deterministic mass-to-parcel conversion and direction sampling.**
- [ ] **Step 4: Mutation check.** Consume ESF stream, use attempt counter, discard residual,
  inject on every rank or hardcode fuel name; RED must fail.
- [ ] **Step 5: Run unit and 1/2-rank tests.**
- [ ] **Step 6: Commit.** Commit `feat: add deterministic spray injection`.

### Task 6L-1: Liquid Property Service and Correlation Laws

**Depends on:** 6F-2.

**Files:**
- Create: `src/spray_liquid_properties.cpp`
- Create: `tests/unit/test_spray_liquid_properties.cpp`
- Modify: `src/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:** Implements `LiquidPropertyService::evaluate(fingerprint,T)` for constant,
temperature-polynomial, Antoine and Clausius--Clapeyron laws with the frozen validity/status
report.

- [ ] **Step 1: Write RED.** SI units, constant/polynomial values, Antoine log10, fallback,
  validity boundaries, non-finite/negative rejection and fingerprint sensitivity.
- [ ] **Step 2: Run RED.** Expect missing service.
- [ ] **Step 3: Implement typed correlation variants.** Never branch on material name.
- [ ] **Step 4: Mutation check.** Celsius/Kelvin mix, natural log, silent extrapolation,
  negative latent heat or species-index identity; RED must fail.
- [ ] **Step 5: Run unit/UBSan.**
- [ ] **Step 6: Commit.** Commit `feat: add liquid property service`.

### Task 6L-2: n-Dodecane and Iso-Octane Property/Mechanism Mapping

**Depends on:** 6F-4, 6L-1.

**Files:**
- Consume when redistributable: `tests/fixtures/liquid/n_dodecane_properties.txt`
- Consume when redistributable: `tests/fixtures/liquid/iso_octane_properties.txt`
- Consume otherwise: `.superpowers/oracles/coast-fuels-source-manifest.md` and its
  hash-verified external asset root
- Create: `tests/unit/test_spray_dual_fuel_identity.cpp`

**Interfaces:** Produces two `LiquidMaterialIdentity`/property packs and runtime gas-species
mapping resolved by name/fingerprint, not compile-time index.

- [ ] **Step 1: Write RED.** Load both packs, verify materially different properties and
  independent mechanism/composition fingerprints.
- [ ] **Step 2: Add cross-wiring RED.** n-dodecane with iso-octane mechanism and absent gas
  species must reject before parcel allocation.
- [ ] **Step 3: Implement config-to-service mapping.** Keep data outside product when
  redistribution is not allowed.
- [ ] **Step 4: Mutation check.** Hardcode species index/name, share one property table or
  ignore mechanism hash; RED must fail.
- [ ] **Step 5: Run source-policy/provenance/unit tests.**
- [ ] **Step 6: Commit.** Commit `test: add dual-surrogate spray identities`.

### Task 6X-1: Gas and One-Third Film Sampling

**Depends on:** 6P-2, Stage 4 thermo/transport services.

**Files:**
- Create: `src/spray_gas_film_detail.hpp`
- Create: `src/spray_gas_film.cpp`
- Create: `tests/unit/test_spray_gas_film.cpp`

**Interfaces:** Produces `GasFilmSample sample_gas_film(MeanStateSnapshot, stencil,
droplet_surface_state)` with the frozen one-third rule.

- [ ] **Step 1: Write RED.** Uniform/linear gas fields, exact film weighting, `p0` not `pi`,
  vapor species mapping and non-finite property failure.
- [ ] **Step 2: Run RED.** Expect missing sampler.
- [ ] **Step 3: Implement through Stage 4 services and shared stencil only.**
- [ ] **Step 4: Mutation check.** Direct Cantera call, wrong film weights, mechanical
  pressure, per-field ESF sample or stale cache; RED must fail.
- [ ] **Step 5: Run focused unit tests.**
- [ ] **Step 6: Commit.** Commit `feat: sample spray gas film state`.

### Task 6X-2: Schiller--Naumann Drag Kernel

**Depends on:** 6F-1, 6X-1.

**Files:**
- Create: `src/spray_drag_detail.hpp`
- Create: `src/spray_drag.cpp`
- Create: `tests/unit/test_spray_drag.cpp`

**Interfaces:** Produces drag acceleration/force report from slip, gas properties, diameter,
liquid mass and multiplicity.

- [ ] **Step 1: Write RED.** Zero slip exact zero, Stokes relaxation, finite-Re reference,
  force direction, rotation/Galilean invariance and validity range.
- [ ] **Step 2: Run RED.** Expect missing kernel.
- [ ] **Step 3: Implement correlation with explicit low-Re limit.**
- [ ] **Step 4: Mutation check.** Reverse slip, omit diameter, use liquid density in Re,
  divide by zero or multiply multiplicity twice; RED must fail.
- [ ] **Step 5: Run unit/UBSan.**
- [ ] **Step 6: Commit.** Commit `feat: add spray drag kernel`.

### Task 6X-3: Ranz--Marshall Heat and Mass Transfer

**Depends on:** 6X-1.

**Files:**
- Create: `src/spray_heat_mass_detail.hpp`
- Create: `src/spray_heat_mass.cpp`
- Create: `tests/unit/test_spray_heat_mass.cpp`

**Interfaces:** Produces `Nu`, `Sh`, heat-transfer and mass-transfer coefficients in SI
units; no evaporation state update.

- [ ] **Step 1: Write RED.** Re=0 limits, independent Pr/Sc, positive finite coefficients,
  heat/mass analogy and dimensional checks.
- [ ] **Step 2: Run RED.** Expect missing correlations.
- [ ] **Step 3: Implement pure functions with explicit applicability status.**
- [ ] **Step 4: Mutation check.** Swap Pr/Sc, wrong exponent, radius/diameter confusion or
  missing conductivity/diffusivity; RED must fail.
- [ ] **Step 5: Run focused unit tests.**
- [ ] **Step 6: Commit.** Commit `feat: add droplet heat and mass transfer`.

### Task 6X-4: Abramzon--Sirignano Evaporation and Terminal Event

**Depends on:** 6L-1, 6X-1, 6X-3, 6F-1.

**Files:**
- Create: `src/spray_evaporation_detail.hpp`
- Create: `src/spray_evaporation.cpp`
- Create: `tests/unit/test_spray_evaporation.cpp`
- Create: `tests/numerical/test_single_droplet.cpp`

**Interfaces:** Produces one-substep liquid mass/temperature/exchange derivative and exact
event-to-zero report for a single-component droplet.

- [ ] **Step 1: Write RED.** Heating-only, no evaporation, simplified d2 limit, Stefan-flow
  correction, complete evaporation event and n-dodecane/iso-octane trend difference.
- [ ] **Step 2: Run RED.** Expect missing model.
- [ ] **Step 3: Implement Spalding/film/equilibrium equations and root-located terminal
  event.** Do not clip negative mass after integration.
- [ ] **Step 4: Mutation check.** Wrong `log1p`, latent sign, surface vapor fraction,
  negative mass clipping or fuel-specific constants; RED must fail.
- [ ] **Step 5: Run unit and short single-drop numerical test.** No Python comparison.
- [ ] **Step 6: Commit.** Commit `feat: evaporate single-component droplets`.

### Task 6X-5: Fixed Predictor/Corrector Exchange Integrator

**Depends on:** 6P-3, 6X-2, 6X-4.

**Files:**
- Create: `src/spray_exchange_integrator_detail.hpp`
- Create: `src/spray_exchange_integrator.cpp`
- Create: `tests/unit/test_spray_exchange_integrator.cpp`

**Interfaces:** Produces one gas-step `SprayExchangeDelta` from bounded internal substeps,
one predictor and one corrector; no gas source publication.

- [ ] **Step 1: Write RED.** Ballistic/drag/heating/evaporation limits, integrated delta,
  event time, substep refinement and retry-from-committed behavior.
- [ ] **Step 2: Run RED.** Expect missing integrator.
- [ ] **Step 3: Implement fixed two-pass exchange with accumulated actual increments.**
- [ ] **Step 4: Mutation check.** Double-count predictor, use endpoint derivative, continue
  after event, reuse failed predictor or add arbitrary outer iterations; RED must fail.
- [ ] **Step 5: Run unit/small sanitizer.**
- [ ] **Step 6: Commit.** Commit `feat: integrate parcel exchange`.

### Task 6X-6: Mean-State Spray Coupling Transaction

**Depends on:** 6P-2, 6X-5, Stage 4 source transaction/PISO hooks. Main-agent algebra review.

**Files:**
- Create: `src/spray_coupling_transaction.cpp`
- Create: `tests/unit/test_spray_coupling_transaction.cpp`
- Create: `tests/mpi/test_spray_coupling_mpi.cpp`

**Interfaces:** Produces staged parcel changes and conservative cell sources; PISO #2 sees
final mass/momentum/enthalpy/divergence source.

- [ ] **Step 1: Write RED.** One/multiple parcels, stencil deposition, closed global budgets,
  vapor species mapping, final-source visibility and rollback.
- [ ] **Step 2: Run RED.** Expect missing transaction.
- [ ] **Step 3: Implement pairwise accumulation/validation/publish.** Publish parcel removal
  and gas source only after collective agreement.
- [ ] **Step 4: Mutation check.** Weight loss, drag/vapor momentum double count, missing
  enthalpy, direct field write or PISO #2 stale source; RED must fail.
- [ ] **Step 5: Run unit and small 1/2-rank tests.**
- [ ] **Step 6: Main-agent commit.** Commit `feat: couple spray sources conservatively`.

### Task 6X-7: Stage 5 ESF Common-Source Adapter

**Depends on:** 6X-6, accepted Stage 5 EsfState/element-consistency contracts.

**Files:**
- Create: `src/spray_esf_source_detail.hpp`
- Create: `src/spray_esf_source.cpp`
- Create: `tests/unit/test_spray_esf_source.cpp`

**Interfaces:** Applies one physical vapor/enthalpy source to all N fields, then invokes the
Stage 5 consistency operator and returns pre/post budgets.

- [ ] **Step 1: Write RED.** N=2/4 common source, no field-specific parcel, ensemble mean
  matches MeanState, element/enthalpy budgets and rollback on infeasible projection.
- [ ] **Step 2: Run RED.** Expect missing adapter.
- [ ] **Step 3: Implement deterministic common-source update.** Reuse one species mapping and
  one physical exchange report.
- [ ] **Step 4: Mutation check.** Recompute per field, inject only field0, draw RNG or skip
  element consistency; RED must fail.
- [ ] **Step 5: Run unit plus small Stage 5 regression.**
- [ ] **Step 6: Main-agent commit.** Commit `feat: couple spray source to ESF fields`.

### Task 6B-1: Static IBM Impact, Rebound and Outlet Ledger

**Depends on:** 6P-3, accepted Stage 3 IBM geometry.

**Files:**
- Create: `src/spray_ibm_impact_detail.hpp`
- Create: `src/spray_ibm_impact.cpp`
- Create: `tests/unit/test_spray_ibm_impact.cpp`
- Create: `tests/mpi/test_spray_ibm_impact_mpi.cpp`

**Interfaces:** Produces earliest segment/surface event, restitution rebound and explicit
outlet removal report; no deposition/film state.

- [ ] **Step 1: Write RED.** Plane/sphere hit/miss, earliest of two hits, normal orientation,
  restitution limits, multiple-event bound, outlet ledger and partition boundary.
- [ ] **Step 2: Run RED.** Expect missing impact handler.
- [ ] **Step 3: Implement using Stage 3 surface query/normal authority.**
- [ ] **Step 4: Mutation check.** Endpoint-only collision, flipped normal, silent delete,
  infinite rebound loop or hidden deposition; RED must fail.
- [ ] **Step 5: Run unit and 1/2-rank tests.**
- [ ] **Step 6: Commit.** Commit `feat: rebound parcels from static IBM`.

### Task 6B-2: TAB Deformation Kernel

**Depends on:** 6F-1, 6F-2.

**Files:**
- Create: `src/spray_tab_breakup_detail.hpp`
- Create: `src/spray_tab_breakup.cpp`
- Create: `tests/unit/test_spray_tab.cpp`

**Interfaces:** Produces committed/trial deformation, deformation rate, threshold/status and
breakup request; no child creation.

- [ ] **Step 1: Write RED.** No forcing, damped oscillation, Weber forcing, threshold,
  dt refinement and rollback.
- [ ] **Step 2: Run RED.** Expect missing TAB kernel.
- [ ] **Step 3: Implement independently from the public TAB equations.** OpenFOAM code is
  not an implementation source.
- [ ] **Step 4: Mutation check.** Wrong restoring/damping sign, missing surface tension,
  endpoint threshold only or advance state on reject; RED must fail.
- [ ] **Step 5: Run focused unit/UBSan.**
- [ ] **Step 6: Commit.** Commit `feat: add TAB deformation model`.

### Task 6B-3: TAB Child Generation and Conservation

**Depends on:** 6B-2, 6P-1, Stage 5 counter RNG domain service.

**Files:**
- Modify: `src/spray_tab_breakup.cpp`
- Create: `tests/unit/test_spray_tab_children.cpp`

**Interfaces:** Produces staged child parcels with IDs derived from parent/accepted step/
ordinal and exact mass/momentum/temperature/multiplicity budgets.

- [ ] **Step 1: Write RED.** Deterministic child count/diameter, total mass/momentum/energy,
  stable IDs, distinct RNG domain and rollback restoring parent/no children.
- [ ] **Step 2: Run RED.** Expect missing child generation.
- [ ] **Step 3: Implement staged split and budget report.** Do not erase parent until commit.
- [ ] **Step 4: Mutation check.** Consume ESF/injector RNG, lose mass, duplicate ID, change
  liquid material or publish before collective agreement; RED must fail.
- [ ] **Step 5: Run unit/small ASan.**
- [ ] **Step 6: Commit.** Commit `feat: generate conservative TAB children`.

### Task 6I-1: Spray Coupling Coordinator

**Depends on:** 6P-4/5, 6X-6/7, 6B-1/3, Stage 4/5 step hooks. Main-agent task.

**Files:**
- Create: `src/flow_spray_coupling_detail.hpp`
- Create: `src/flow_spray_coupling.cpp`
- Create: `tests/unit/test_flow_spray_coupling.cpp`

**Interfaces:** Produces `attempt_spray_step(...)` and `SprayStepReport` with fixed injector,
sampling, predictor/corrector, source, chemistry/PISO, impact/breakup/migration order.

- [ ] **Step 1: Write operator-spy RED.** Require one physical exchange per parcel, second
  chemistry half sees vapor, PISO #2 sees final source and migration/breakup publish last.
- [ ] **Step 2: Run RED.** Expect missing coordinator.
- [ ] **Step 3: Implement one atomic trial using approved Stage 4/5 hooks.**
- [ ] **Step 4: Mutation check.** Per-field parcels, source after chemistry/PISO, predictor
  double count, third PISO or early migration publish; RED must fail.
- [ ] **Step 5: Run unit and <=12^3 smoke with the Stage 4
  `tests/support/chem_analytic_backend.*` fixture.** The fixture is not linked into `hundun`.
- [ ] **Step 6: Main-agent commit.** Commit `feat: coordinate reacting spray step`.

### Task 6I-2: Collective Failure and Complete Rollback

**Depends on:** 6I-1.

**Files:**
- Create: `tests/mpi/test_spray_rollback.cpp`
- Modify: `src/flow_spray_coupling.cpp`

**Interfaces:** Finalizes rollback for gas, ESF, parcels, owners, injector residual, TAB,
created/deleted IDs and RNG accepted clocks.

- [ ] **Step 1: Inject failure at injection, locate, drag, evaporation, grid-source deposition, ESF
  projection, IBM, TAB, migration, chemistry and PISO #2.**
- [ ] **Step 2: Require bytewise committed-state equality and all-rank decision.**
- [ ] **Step 3: Implement missing rollback/collective mappings only.** Retry recomputes the
  entire parcel trial.
- [ ] **Step 4: Mutation check.** Keep child, consume injector residual, publish migration,
  remove evaporated parcel or let one rank return early; RED must fail.
- [ ] **Step 5: Run small 1/2/4-rank failure tests.**
- [ ] **Step 6: Commit.** Commit `fix: make spray trials collectively atomic`.

### Task 6I-3: Schema v6, Broadcast and Driver

**Depends on:** 6F-3, 6I-1/2.

**Files:**
- Create: `src/cfg_resolved_case_v6_loader.cpp`
- Create: `src/cfg_resolved_case_v6_loader_detail.hpp`
- Create: `src/app_resolved_case_v6_broadcast.cpp`
- Create: `src/app_spray_driver.cpp`
- Create: `src/app_spray_driver_detail.hpp`
- Create: `tests/unit/test_resolved_case_v6.cpp`
- Create: `tests/mpi/test_resolved_case_v6_broadcast.cpp`
- Create: `tests/unit/test_spray_driver.cpp`

**Interfaces:** Produces `ResolvedSprayCaseV6` loader/broadcast and
`run_spray_case(...)` within the same `hundun` executable.

- [ ] **Step 1: Write RED.** one/two-way, evaporation, rebound, none/TAB, injector,
  property/mechanism identities, TPDF combinations and unsupported features.
- [ ] **Step 2: Run RED.** Expect missing implementation.
- [ ] **Step 3: Implement rank-0 parse/broadcast and composition root.** Validate all
  mappings before allocating parcels.
- [ ] **Step 4: Mutation check.** Hardcode fuel, silently enable film/KHRT, bypass Stage 5
  common source or alter old schemas; RED must fail.
- [ ] **Step 5: Run unit, 1/2-rank, validate/print-resolved and one-step smoke.**
- [ ] **Step 6: Commit.** Commit `feat: add spray schema and driver`.

### Task 6I-4: Checkpoint v6 Parcel Sections

**Depends on:** 6F-3, 6I-2/3.

**Files:**
- Create: `src/flow_checkpoint_v6.cpp`
- Create: `src/flow_checkpoint_v6_detail.hpp`
- Create: `tests/unit/test_checkpoint_v6.cpp`
- Create: `tests/mpi/test_checkpoint_v6_mpi.cpp`

**Interfaces:** Registers parcel arrays/IDs/owners, injector residual/counter, TAB history,
material/mechanism identity and Stage 5 state references. Trial buffers excluded.

- [ ] **Step 1: Write RED.** continuous-vs-restart, zero/many parcels, migrated parcels,
  injector/TAB continuation, corruption, identity mismatch and failed-read unchanged.
- [ ] **Step 2: Run RED.** Expect missing codecs.
- [ ] **Step 3: Implement deterministic validate-then-publish sections.** Same partition only.
- [ ] **Step 4: Mutation check.** Omit multiplicity/owner/TAB state, persist outgoing buffer,
  accept material mismatch or partially publish; RED must fail.
- [ ] **Step 5: Run unit and 1/2/4-rank small continuation.**
- [ ] **Step 6: Commit.** Commit `feat: checkpoint spray state v6`.

### Task 6I-5: Spray Diagnostics and Exact Counters

**Depends on:** 6F-3, 6I-3/4.

**Files:**
- Create: `src/diag_spray.cpp`
- Create: `tests/unit/test_spray_diagnostics.cpp`
- Modify: integration-owned diagnostic registry

**Interfaces:** Providers expose parcel counts/mass, injection, migration, drag, heating,
evaporation, complete-event, IBM impact, TAB, source budgets, failures and fuel identities.

- [ ] **Step 1: Write RED.** Stable IDs/units, exact counters, zero-parcel legitimate values,
  absent module unregistered, read-only callbacks and no hidden collective.
- [ ] **Step 2: Run RED.** Expect missing providers.
- [ ] **Step 3: Implement adapters over authority reports.** Do not recompute exchange.
- [ ] **Step 4: Mutation check.** Fake-zero absent module, omit outflow ledger, collapse
  failure classes, mutate container or unordered output; RED must fail.
- [ ] **Step 5: Run unit plus 1/2-rank diagnostic session.**
- [ ] **Step 6: Commit.** Commit `feat: diagnose spray coupling`.

### Task 6A-1: Parcel Mechanics and Migration Gate

**Depends on:** 6P-1--6P-5, 6B-1.

**Files:**
- Create: `tests/acceptance/stage6_parcel_gate.sh`
- Create: `.superpowers/sdd/stage6-6A-1-parcel-gate.md`

**Interfaces:** Produces mechanics/injection/migration/IBM evidence independent of
evaporation chemistry.

- [ ] **Step 1: Freeze selectors.** Container, stencil, ballistic/Stokes/Galilean,
  injector, 1/2/4-rank migration, rollback, Restart and rebound.
- [ ] **Step 2: Assert no evaporation/chemistry/48^3 selector is required.**
- [ ] **Step 3: Run Debug/focused Release/small sanitizer evidence.**
- [ ] **Step 4: Hash commands/binaries/logs and verify no background MPI remains.**
- [ ] **Step 5: Main-agent parcel ownership/full-diff review.**
- [ ] **Step 6: Commit receipt.** Commit `test: accept parcel mechanics`.

### Task 6A-2: Evaporation and Two-Way Conservation Gate

**Depends on:** 6L-1, 6X-1--6X-6, 6A-1.

**Files:**
- Create: `tests/acceptance/stage6_evaporation_gate.sh`
- Create: `.superpowers/sdd/stage6-6A-2-evaporation-gate.md`

**Interfaces:** Produces single-drop and gas/parcel conservation evidence.

- [ ] **Step 1: Run heating, d2 oracle, A--S short single-drop, event-to-zero and substep
  refinement.**
- [ ] **Step 2: Run one/multiple-parcel mass/momentum/h_tc closure and PISO #2 source spy.**
- [ ] **Step 3: Run small 1/2-rank two-way case and focused sanitizer.**
- [ ] **Step 4: Verify no clipping/negative state and exact source ledger.**
- [ ] **Step 5: Main-agent energy/source algebra review.**
- [ ] **Step 6: Commit receipt.** Commit `test: accept spray evaporation coupling`.

### Task 6A-3: Dual-Fuel, TAB and Stage 5 Combination Gate

**Depends on:** 6L-2, 6X-7, 6B-2/3, 6I-1--6I-5, 6A-2.

**Files:**
- Create: `tests/acceptance/stage6_acceptance.sh`
- Create: `.superpowers/sdd/stage6-6A-3-framework-gate.md`
- Create: `tests/numerical/test_reacting_spray_smoke.cpp`
- Modify: `VERSION`
- Modify: `docs/numerics/stage6-capability-ledger.md`

**Interfaces:** Produces frozen Stage 6 code candidate `C6` at `0.5.0-rc.1` and the
complete low-cost Stage 6 development evidence manifest.

- [ ] **Step 1: Finalize and freeze candidate.** Set `VERSION=0.5.0-rc.1`, finish all
  product/test source and selectors, run a low-cost preflight, create signed `C6`, and
  record HEAD/tree/binary/dirty hashes. Do not modify product or tests afterward.
- [ ] **Step 2: Run n-dodecane and iso-octane short single-drop/backend handoff and
  cross-mechanism rejection.**
- [ ] **Step 3: Run two <=8^3 reacting spray smokes, one per surrogate; run N=2 and N=4
  common-source tests without duplicating parcel calculations.**
- [ ] **Step 4: Run TAB child conservation, small 1/2/4-rank, Checkpoint, diagnostics,
  driver, headers/tests-off where affected.**
- [ ] **Step 5: Run complete affected Debug and focused Release/ASan/UBSan on `C6`, then
  perform the main-agent full Stage 6 diff/provenance/caller/science review.** A product
  fix creates a new `C6` and invalidates only consuming evidence.
- [ ] **Step 6: Commit governance receipt.** Commit only the evidence manifest/ledger as
  `test: record Stage 6 development evidence`, with `accepted_code_head=C6`.

### Task 6A-4: `0.5.0-rc.1` Development-Complete Seal

**Depends on:** 6A-3.

**Files:**
- Create: `.superpowers/sdd/stage6-development-complete-report.md`
- Modify: `docs/numerics/stage6-capability-ledger.md`
- Modify: `docs/numerics/stage4-6-capability-root.md`
- Modify: `AGENTS.md`

**Interfaces:** Produces governance seal `G6` declaring tested code candidate `C6` as
`V1_DEVELOPMENT_COMPLETE` at `0.5.0-rc.1`; it is not yet `V1_ACCEPT` and does not sync
the product repository.

- [ ] **Step 1: Recover `C6` and record its Stage 5 parent/tree/diff, binaries, mechanisms,
  liquid packs, logs and DCO.** Assert all later commits are governance-only.
- [ ] **Step 2: Verify no mechanism/license leakage, direct Cantera call, per-field parcel,
  background MPI or uncommitted product/test change.**
- [ ] **Step 3: State deferred capabilities and dual-surrogate claim limits.**
- [ ] **Step 4: Record governance version `0.5.0-rc.1`.** Do not modify `VERSION`, product
  code or tests, and do not tag/project `1.0.0`.
- [ ] **Step 5: Run final low-cost manifest/hash audit and `git diff --check`.**
- [ ] **Step 6: Main-agent decision/governance commit.** Write
  `V1_DEVELOPMENT_COMPLETE` or `REJECT`; commit `docs: freeze HUNDUN-FLOW v1 release
  candidate` only when `C6` passes. Record `accepted_code_head=C6`; do not call the
  governance commit the tested code HEAD.

## 4. Critical Path and Optional Parallel Edge

```text
6F-0 -> 6F-1 -> 6F-2/3
6F-2 -> 6P-1 -> 6P-2 -> 6P-3 -> 6P-4
6P-1 -> 6P-5
6F-2 -> 6L-1
6F-4 + 6L-1 -> 6L-2
6P-2 + Stage4 services -> 6X-1 -> 6X-2/3
6L-1 + 6X-3 -> 6X-4
6P-3 + 6X-2/4 -> 6X-5 -> 6X-6
Stage5 contracts + 6X-6 -> 6X-7
6P-3 -> 6B-1
6F-2 -> 6B-2 -> 6B-3
6P-4/5 + 6X-6/7 + 6B-1/3 -> 6I-1 -> 6I-2/3/4/5
6A-1 -> 6A-2 -> 6A-3 -> 6A-4
```

默认在 Stage 5 接受后串行执行。若用户在 Stage 4 或 Stage 5 节点批准并行，允许
提前执行 `6F-1/2` 的设计准备和 `6P-1/2/3`、`6L-1`、`6X-2/3`、`6B-2` 纯核；
不得提前接受 Stage 6、修改 shared integration files 或运行 6X-7/6I。

## 5. Task Reference Matrix

URL、DOI、license class 和 fuel-source 边界以
`docs/references/2026-08-09-hundun-flow-stage4-6-reference-catalog.md` 为准。Pele/AMReX
用于架构和测试思想；Code_Saturne/OpenFOAM 的 GPL 源码不得复制。COAST fuel files 只有
在用户确认 realpath 且许可审计通过后才能成为外部或 tracked test asset。

| Task | Reference | Exact point to reuse | Explicitly avoid |
|---|---|---|---|
| 6F-0 | accepted Stage 4/5 receipts | exact services, common-source hook and registry intake | starting from unaccepted or mixed heads |
| 6F-1 | Abramzon--Sirignano；Code_Saturne coupling categories | parcel/gas mass, momentum, energy and source-sign proof | latent/vapor-momentum double count |
| 6F-2 | AMReX ParticleContainer；PeleMP spray | SoA identity, material service and trial transaction boundaries | AMR/GPU runtime dependency |
| 6F-3 | accepted v5 registries | append-only v6 schema/checkpoint/diagnostic IDs | film/KH--RT or fake-zero presence |
| 6F-4 | confirmed COAST `EXEC/Fuels`；LLNL/Cantera mechanism pages | dual-surrogate identity, hash and redistribution decision | assuming Cantera license covers mechanism data |
| 6P-1 | AMReX ParticleContainer | pure SoA, stable IDs and staged create/remove | AMReX container/API or nondeterministic order |
| 6P-2 | AMReX particle interpolation；Code_Saturne coupling split | one paired sample/deposit stencil with partition ties | separate sampling and source weights |
| 6P-3 | Code_Saturne Lagrangian architecture；standard second-order ODE methods | segment/event-aware predictor-corrector | opaque legacy integrator control flow |
| 6P-4 | AMReX particle redistribution | one owner, staged migration and deterministic packing | AMR DistributionMap or early publish |
| 6P-5 | Salmon et al. Philox | injector-specific counter domain and retry identity | shared ESF stream or mutable cursor |
| 6L-1 | PeleMP liquid property service；standard Antoine/Clausius laws | typed material-independent SI correlations and validity | fuel-name branches or silent extrapolation |
| 6L-2 | LLNL n-dodecane/iso-octane pages；confirmed COAST fuel assets | two materially distinct identities and gas-species mapping | fixed species index or unlicensed redistribution |
| 6X-1 | PelePhysics/PeleMP spray equations | one-third film rule and service-mediated gas properties | direct Cantera or per-field gas samples |
| 6X-2 | Schiller--Naumann | Stokes limit, finite-Re sphere drag and validity status | copied solver source or sign ambiguity |
| 6X-3 | Ranz--Marshall | independent Nu/Sh limits and SI transfer coefficients | radius/diameter or Pr/Sc swaps |
| 6X-4 | Abramzon--Sirignano；PeleMP equations | single-component Stefan-flow evaporation and terminal event | multicomponent scope or negative-mass clipping |
| 6X-5 | second-order predictor/corrector；PeleMP substep concepts | integrated actual exchange and bounded event substeps | endpoint derivative or arbitrary outer iteration |
| 6X-6 | Code_Saturne two-way category；6F-1 proof | paired conservative gas/parcel publication before PISO #2 | direct field writes or double deposition |
| 6X-7 | Valiño ESF；Xu element consistency | identical physical source to all fields then one consistency solve | N parcel copies or field-0 injection |
| 6B-1 | accepted Stage 3 IBM geometry；Code_Saturne boundary separation | earliest segment hit, rebound and explicit outlet ledger | deposition/film/splash implementation |
| 6B-2 | O'Rourke--Amsden TAB；OpenFOAM model catalog | public deformation ODE, threshold and rollback | GPL source translation or KH--RT expansion |
| 6B-3 | TAB paper；counter RNG | deterministic child IDs and exact budgets | shared ESF/injector stream or early parent erase |
| 6I-1 | PeleMP module ordering；accepted Stage 4/5 hooks | one parcel exchange, chemistry visibility and publish-last composition | per-field parcels or third PISO |
| 6I-2 | MPI standard；HUNDUN transactions | all-rank failure and rollback of every parcel-related state | rank-local early return or partial migration |
| 6I-3 | yyjson/MPI broadcast；existing composition root | typed v6 schema and one executable dispatch | fuel hardcoding or old-schema mutation |
| 6I-4 | AMReX particle checkpoint concepts；accepted v5 protocol | deterministic same-partition parcel sections | AMReX format or rank-changing Restart |
| 6I-5 | Code_Saturne statistics separation；HUNDUN diagnostics | read-only authority reports and exact counters | recomputing exchange or fake-zero absence |
| 6A-1 | AMReX particle/Code_Saturne mechanics concepts | compact ownership, migration, Restart and rebound gate | chemistry/48^3 expansion |
| 6A-2 | Ranz--Marshall；Abramzon--Sirignano | single-drop limits and exact two-way budgets | Python oracle or empirical case tuning |
| 6A-3 | dual-surrogate sources；TAB；Stage 5 service contracts | no-hardcoding evidence and N=2/4 common-source integration | real kerosene/gasoline validation claim |
| 6A-4 | cross-stage exact-HEAD protocol；DCO | `C6`/`G6` separation and development-complete limits | calling rc evidence `V1_ACCEPT` |

## 6. Receipt Requirements

每个 task receipt 除通用 HEAD/tree/diff/commands/log/DCO 外，必须记录：

```text
parcel identity/owner/multiplicity impact
liquid material and gas species fingerprints
mass/momentum/h_tc budget before and after
sampling/deposition stencil identity
event/substep/predictor-corrector counts
injection/TAB RNG domain and accepted-step semantics
migration/rollback/Checkpoint impact
MeanState/ESF/PISO authority impact
fuel/mechanism license and hash identity
deferred physical models and capability limits
```

worker 只可实现冻结的 container、pure correlation、codec、fixture 或 adapter。主
agent保留 conservation derivation、exchange composition、Stage 5/PISO timing、完整
rollback、dual-fuel claims 和最终接受。

# HUNDUN-FLOW Stage 5 ESF/TPDF/TCR Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:executing-plans` for the default serial path. Use `superpowers:subagent-driven-development` only after the user explicitly authorizes parallel execution at a stage gate. Execute one task at a time and return to the main agent for review.

**Goal:** 在 accepted Stage 4 reacting-flow 上实现低雷诺一致 ESF/TPDF、exact IEM、reaction ensemble、PSR shadow、用户 TCR、Checkpoint v5、诊断和受控 COAST 差分 oracle。

**Architecture:** MeanState 保持守恒/PISO/面通量权威；每个 spatial rank 保存全部偶数 N 个 stochastic fields。随机场跨步连续，只有 Wiener increment 每个 accepted stochastic interval 更新。Stage 5 只通过 ChemistryBackend 使用 Cantera，COAST 只在用户确认后作为进程外私有 oracle。

**Tech Stack:** C++17、MPI-3、Stage 4 internal services、Philox counter RNG、CTest、可选 standalone Fortran oracle driver、Apache-2.0/DCO。

## Global Constraints

- 实施 parent 必须是 Stage 4 `STAGE4_ACCEPT` exact HEAD。
- 默认串行执行；Stage 4 节点未获用户并行授权时，不提前启动 Stage 5 产品实现。
- N 必须为偶数，minimum `N=2`，transient recommended `N=4`；N=2 标记 `minimal_pair_sampling`。
- 随机场状态跨步连续；只更新 `DeltaW`。禁止每步重置 field composition/enthalpy。
- 每 spatial rank 保存 owned/ghost 空间的全部 N fields；v1 不做 field-index MPI decomposition。
- ESF Wiener key 不含 cell、species、rank、PISO iteration 或 failed-attempt ordinal。
- molecular diffusion 不进入 Wiener coefficient；`mu_t=+0.0` 必须精确退化。
- TPDF 模式使用单一 composition diffusivity，配置要求 `Pr_t=Sc_t`。
- IEM 与 TCR 共享一个 mixer；`kappa=1` 是严格 IEM 极限。
- MeanState chemistry closure 使用 ensemble integrated increments；不得用 field 0 或 endpoint rate。
- COAST source、case、data、executable 不进入 Git、产品、安装包或 ABI；正式读取前必须请用户确认 realpath/version。
- 不运行 Vblowoff、Flame D、96^3 或大型 sanitizer/MPI。
- worker 不决定 ESF 方程、feedback timing、root branch、element correction 或完整组合；主 agent拥有这些科学任务。
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
include/hundun/comb_esf_reports.hpp
include/hundun/comb_tcr_reports.hpp
include/hundun/cfg_resolved_case_v5.hpp
include/hundun/cfg_resolved_case_v5_loader.hpp
include/hundun/flow_checkpoint_v5.hpp
include/hundun/diag_esf_tcr.hpp
```

### Internal state, kernels and composition

```text
src/comb_esf_state_detail.hpp
src/comb_esf_state.cpp
src/rt_counter_rng_philox_detail.hpp
src/rt_counter_rng_philox.cpp
src/comb_esf_wiener_detail.hpp
src/comb_esf_wiener.cpp
src/comb_esf_transport_detail.hpp
src/comb_esf_transport.cpp
src/comb_esf_iem_detail.hpp
src/comb_esf_iem.cpp
src/comb_esf_element_consistency_detail.hpp
src/comb_esf_element_consistency.cpp
src/comb_esf_reaction_ensemble_detail.hpp
src/comb_esf_reaction_ensemble.cpp
src/comb_esf_psr_shadow_detail.hpp
src/comb_esf_psr_shadow.cpp
src/comb_tcr_algebra_detail.hpp
src/comb_tcr_algebra.cpp
src/comb_tcr_root_state_detail.hpp
src/comb_tcr_root_state.cpp
src/comb_tcr_controller_detail.hpp
src/comb_tcr_controller.cpp
src/flow_tpdf_tcr_coupling_detail.hpp
src/flow_tpdf_tcr_coupling.cpp
src/flow_tpdf_tcr_immersed_detail.hpp
src/flow_tpdf_tcr_immersed.cpp
src/cfg_resolved_case_v5_loader.cpp
src/cfg_resolved_case_v5_loader_detail.hpp
src/app_resolved_case_v5_broadcast.cpp
src/app_tpdf_tcr_driver.cpp
src/app_tpdf_tcr_driver_detail.hpp
src/flow_checkpoint_v5.cpp
src/flow_checkpoint_v5_detail.hpp
src/diag_esf_tcr.cpp
```

### Private-oracle governance and HUNDUN-owned drivers

```text
.superpowers/oracles/coast-esf-source-manifest.md
.superpowers/oracles/coast-tcr-source-manifest.md
tests/oracles/coast_oracle_protocol.hpp
tests/oracles/coast_esf_driver.F90
tests/oracles/coast_tcr_driver.F90
tests/oracles/coast_oracle_runner.cpp
tests/cmake/coast_oracle_contract.cmake
tests/fixtures/esf_oracle_vectors.txt
tests/fixtures/tcr_oracle_vectors.txt
```

No COAST source file appears in this map. Formal runners copy allowlisted source to an
untracked generated directory outside the source tree.

## 2. Frozen Internal Interfaces

The report/value types below live in `namespace hundun::combustion`; internal ESF/TCR
state and kernels use the same full domain namespace rather than `esf`/`tcr` abbreviations.

```cpp
namespace hundun::combustion {

enum class StochasticPurpose : std::uint32_t {
  esf_transport = 1,
  oracle_injected_sign = 2
};

struct EsfCounterAddress {
  std::uint64_t seed;
  std::uint64_t accepted_step;
  std::uint32_t stochastic_stage;
  std::uint32_t field_pair;
  std::uint32_t spatial_direction;
  StochasticPurpose purpose;
};

struct WienerIncrement {
  std::uint32_t field_id;
  std::array<double, 3> delta_w;
};

enum class TcrStatus : std::uint32_t {
  inactive,
  direct_root,
  similar_state_fallback,
  fold_unresolved,
  weak_denominator,
  algebraic_invalid,
  statistical_unresolved,
  branch_crossing
};

enum class TcrMode : std::uint32_t {
  off,
  shadow,
  experimental,
  validated
};

struct TcrAlgebraInput {
  double eta;
  double rate_ratio_r;
  double previous_signed_root;
};

struct TcrAlgebraReport {
  double discriminant;
  double signed_root;
  double kappa;
  TcrStatus status;
};

struct TcrRootHistory {
  std::array<double, 9> values;
};

}  // namespace hundun::combustion
```

Internal state addressing is fixed as
`(global_cell_id, stochastic_field_id, component_id)`. Persistent arrays use
component-major SoA; chemistry packs are attempt-local scratch.

### Frozen schema v5 additions

Schema v5 extends accepted v4. It adds one `combustion` object with exact keys:

| JSON path | Type / allowed value | Rule |
|---|---|---|
| `schema_version` | integer `5` | required |
| `combustion.model` | `none` or `esf_tpdf` | required |
| `combustion.stochastic_fields.count` | even integer `>=2` | default `4`; explicit `2` for debug/steady |
| `combustion.stochastic_fields.seed` | unsigned 64-bit integer | required when ESF enabled |
| `combustion.stochastic_fields.rng` | `philox4x32_10_v1` | required, no alternate engine in v1 |
| `combustion.mixing.model` | `iem_tcr` | one mixer authority |
| `combustion.mixing.time_scale_s` | finite positive real | required, no hidden default |
| `combustion.tcr.mode` | `off`, `shadow`, `experimental`, `validated` | default `shadow` |
| `combustion.tcr.oracle_fingerprint` | unsigned 64-bit integer | required iff mode is `validated` |

When `combustion.model=esf_tpdf`, accepted v4
`les.wale.turbulent_prandtl` and `les.wale.turbulent_schmidt` must be bitwise equal after
canonical parsing；the loader does not add duplicate diffusivity keys. Unknown/duplicate keys,
odd N, N<2, validated-without-matching-oracle and a field-index decomposition key are rejected.
`--print-resolved` includes the sampling label and RNG/oracle identities.

## 3. Task Sequence

### Common command protocol

Pure contract/RNG/algebra tasks may configure with Cantera disabled. Any chemistry, driver or
product smoke uses the accepted Stage 4 package root:

```bash
cmake -S . -B build/stage5-debug \
  -DCMAKE_BUILD_TYPE=Debug \
  -DHUNDUN_BUILD_TESTS=ON \
  -DHUNDUN_ENABLE_CANTERA=ON \
  -DHUNDUN_CANTERA_PACKAGE_ROOT="${HUNDUN_ACCEPTED_CANTERA_ROOT}"
cmake --build build/stage5-debug -j32
ctest --test-dir build/stage5-debug -R "${HUNDUN_TASK_TEST_REGEX}" \
  --output-on-failure -j24
```

The main agent exports `HUNDUN_ACCEPTED_CANTERA_ROOT` from G4 and sets
`HUNDUN_TASK_TEST_REGEX` to the task's exact test basename(s). For a Cantera-free pure task,
set `HUNDUN_ENABLE_CANTERA=OFF` and omit the package root. Focused Release uses
`build/stage5-release` and `Release`. Live COAST oracle commands are never implicit in CTest；
they run only after 5E-R1/5T-0 and record the confirmed manifest and external build root.

### Task 5F-0: Stage 4 Intake and Compatibility Freeze

**Depends on:** Stage 4 exact acceptance.

**Files:**
- Create: `.superpowers/sdd/stage5-5F-0-baseline-receipt.md`
- Create: `docs/numerics/stage5-capability-ledger.md`
- Modify: `AGENTS.md`
- Test: `tests/cmake/stage5_source_policy.cmake`

**Interfaces:** Produces exact parent/tree, Stage 4 service signatures, registry IDs and
the Stage 5 allowed-file table.

- [ ] **Step 1: Record accepted evidence.** Capture Stage 4 HEAD/tree, worktree status,
  CompositionIdentity, ChemistryBackend, source transaction, Checkpoint v4 and diagnostic
  provider identities.
- [ ] **Step 2: Add policy RED.** Reject per-cell ESF RNG keys, `Cantera::` outside Stage 4
  adapter, COAST source in tracked files, Vblowoff selectors and a second flow authority.
- [ ] **Step 3: Run RED.** Expect failure before the Stage 5 ledger exists.
- [ ] **Step 4: Write governance only.** Add default serial gate and no-product-change
  receipt; do not alter Stage 4 types.
- [ ] **Step 5: Run policy and `git diff --check`.**
- [ ] **Step 6: Main-agent commit.** Commit `docs: freeze Stage 5 accepted baseline`.

### Task 5F-1: ESF/TCR Equations, Units and Coupling Contract

**Depends on:** 5F-0.

**Files:**
- Create: `docs/numerics/stage5-esf-tcr-equations.md`
- Create: `include/hundun/comb_esf_reports.hpp`
- Create: `include/hundun/comb_tcr_reports.hpp`
- Create: `tests/unit/test_esf_contracts.cpp`
- Modify: `src/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:** Produces report enums/types, exact low-Re ESF equation, unity-Lewis rule,
N policy, IEM exponent, TCR equations, publication identifiers and failure taxonomy.

- [ ] **Step 1: Write dimension/limit RED.** Verify `mu_t=0`, `dt=0`, constant field,
  `kappa=1`, N=2 pair and every TcrStatus serialization.
- [ ] **Step 2: Run RED.** Expect missing report types.
- [ ] **Step 3: Implement only value/report types.** No transport or TCR algorithm yet.
- [ ] **Step 4: Mutation check.** Put molecular diffusion in noise, remove IEM `1/2`, merge
  statuses or allow odd N; RED must fail.
- [ ] **Step 5: Main-agent equation review against cited papers and approved spec.** Record
  stable bibliography keys for product-source comments. A COAST path is not a scientific
  citation and must not replace the user-confirmed TCR publication.
- [ ] **Step 6: Commit.** Commit `docs: freeze ESF and TCR equations`.

### Task 5F-2: State, RNG, Retry and Service Contracts

**Depends on:** 5F-1.

**Files:**
- Create: `src/comb_esf_state_detail.hpp`
- Create: `src/rt_counter_rng_philox_detail.hpp`
- Create: `src/comb_esf_wiener_detail.hpp`
- Create: `tests/unit/test_esf_state_contract.cpp`

**Interfaces:** Produces `EsfStateLayout`, `EsfAttemptState`, `EsfCounterAddress`,
`WienerIncrement`, `PsrShadowReport` and `TcrRootHistory` ownership contracts.

- [ ] **Step 1: Write RED.** N=2/4 layouts, all fields per rank, committed/trial history,
  retry counter reuse, no mutable RNG cursor and nine-value root history.
- [ ] **Step 2: Run RED.** Expect missing internal contracts.
- [ ] **Step 3: Implement descriptors/value types only.** FieldRegistry owns arrays;
  PsrShadow is non-persistent; TcrRootHistory is persistent.
- [ ] **Step 4: Mutation check.** Reset fields per step, include rank/cell in key, persist
  workspace or omit one history value; RED must fail.
- [ ] **Step 5: Run unit/header/tests-off compile.**
- [ ] **Step 6: Commit.** Commit `feat: define ESF state contracts`.

### Task 5F-3: Schema v5, Persistence and Diagnostic Identities

**Depends on:** 5F-1, 5F-2, Stage 4 registries.

**Files:**
- Create: `include/hundun/cfg_resolved_case_v5.hpp`
- Create: `include/hundun/cfg_resolved_case_v5_loader.hpp`
- Create: `include/hundun/flow_checkpoint_v5.hpp`
- Create: `include/hundun/diag_esf_tcr.hpp`
- Create: `tests/unit/test_stage5_identity_ledger.cpp`

**Interfaces:** Freezes schema keys, checkpoint section IDs, diagnostic kinds and
capabilities without implementing codecs/providers.

- [ ] **Step 1: Write RED.** N even/minimum/default, `Pr_t=Sc_t`, seed/mode/oracle identity,
  duplicate ID rejection and absence semantics.
- [ ] **Step 2: Run RED.** Expect missing v5 identities.
- [ ] **Step 3: Implement append-only value definitions.** Do not change v1--v4 IDs.
- [ ] **Step 4: Mutation check.** Default N=8, allow N=1, fake-zero absent provider or reuse
  a v4 section ID; RED must fail.
- [ ] **Step 5: Run focused unit/header tests.**
- [ ] **Step 6: Commit.** Commit `feat: freeze Stage 5 identities`.

### Task 5E-R1: COAST ESF Semantic Audit

**Depends on:** 5F-0 and explicit user confirmation of the actual COAST source root.
The main agent owns the confirmation gate; a worker may receive only the frozen allowlist.

**Files:**
- Create: `.superpowers/oracles/coast-esf-source-manifest.md`
- Create: `docs/numerics/stage5-coast-esf-semantic-audit.md`
- No product source changes.

**Interfaces:** Produces a versioned behavior table for field persistence, stochastic
increment timing, antithetic pairing, N=2/4, IEM, boundaries, Restart and failure handling.

- [ ] **Step 1: Main agent stops at the user gate.** Ask for/confirm realpath and acceptable mature
  version(s); record commit/status before reading implementation files.
- [ ] **Step 2: Hash only the confirmed allowlist.** Do not inspect cases, research data,
  binaries or unrelated source.
- [ ] **Step 3: Extract semantic facts.** Record equations/units/call timing and classify
  each fact as paper-derived, implementation-observed or HUNDUN decision.
- [ ] **Step 4: Compare with 5F contracts.** Any conflict returns to main-agent spec
  amendment; worker may not choose.
- [ ] **Step 5: Run leak scan.** Ensure no source excerpt, private path beyond governance
  manifest or legacy name enters public docs/product.
- [ ] **Step 6: Commit governance only.** Commit `docs: audit COAST ESF semantics`.

### Task 5E-R2: COAST ESF Differential Oracle

**Depends on:** 5E-R1, 5E-1--5E-6 interfaces frozen.

**Files:**
- Create: `tests/oracles/coast_oracle_protocol.hpp`
- Create: `tests/oracles/coast_esf_driver.F90`
- Create: `tests/oracles/coast_oracle_runner.cpp`
- Create: `tests/fixtures/esf_oracle_vectors.txt`
- Create: `tests/cmake/coast_oracle_contract.cmake`

**Interfaces:** Produces a standalone process protocol accepting synthetic field arrays,
injected Wiener signs and one-step parameters; normal tests consume source-free fixtures.

- [ ] **Step 1: Write protocol RED.** Reject wrong version/component count/endianness and
  require exact status plus precision-aware numeric fields.
- [ ] **Step 2: Write HUNDUN-owned driver.** It may call only confirmed allowlisted pure
  COAST modules copied byte-for-byte into an untracked generated build directory.
- [ ] **Step 3: Build/run outside source tree.** Record source/driver/binary/output hashes;
  never link the oracle into `hundun`.
- [ ] **Step 4: Generate synthetic N=2/4 fixtures.** Cover persistence, one ESF increment,
  IEM and Restart continuation; no research case.
- [ ] **Step 5: Run source/package leak scans and source-free fixture tests.**
- [ ] **Step 6: Commit only driver/protocol/fixtures/manifest.** Commit
  `test: add COAST ESF differential oracle`.

### Task 5E-1: Philox Counter Kernel

**Depends on:** 5F-2.

**Files:**
- Create: `src/rt_counter_rng_philox.cpp`
- Create: `tests/unit/test_rng_philox.cpp`
- Modify: `src/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:** Produces `std::array<std::uint32_t,4> philox4x32_10(const
EsfCounterAddress&, std::uint32_t lane)`; no engine object or cursor.

- [ ] **Step 1: Write RED from published vectors.** Cover all-zero, nonzero key/counter,
  endianness and algorithm version.
- [ ] **Step 2: Run RED.** Expect missing function.
- [ ] **Step 3: Implement ten-round integer kernel from the public algorithm.** Use project
  names/style; do not copy Random123 source text.
- [ ] **Step 4: Mutation check.** Alter round count, multiply constants, word permutation or
  omit algorithm version; golden tests must fail.
- [ ] **Step 5: Run Debug/Release unit and UBSan.**
- [ ] **Step 6: Commit.** Commit `feat: add counter-based Philox kernel`.

### Task 5E-2: Balanced Wiener Scheduler

**Depends on:** 5E-1.

**Files:**
- Create: `src/comb_esf_wiener.cpp`
- Create: `tests/unit/test_esf_wiener.cpp`

**Interfaces:** Produces `make_balanced_wiener_increments(N,dt,address)` returning N
increments with `DeltaW(2m+1)=-DeltaW(2m)` and magnitude `sqrt(dt)`.

- [ ] **Step 1: Write RED.** N=2/4, three directions, zero dt, retry same sign/new magnitude,
  next accepted step new signs and decomposition-independent vectors.
- [ ] **Step 2: Run RED.** Expect missing scheduler.
- [ ] **Step 3: Implement domain-separated sign extraction.** Never accept cell/species/rank
  as an argument.
- [ ] **Step 4: Mutation check.** Odd N, unpaired signs, per-PISO redraw, failed-attempt
  counter or `dt` instead of `sqrt(dt)`; RED must fail.
- [ ] **Step 5: Run focused tests and fixture comparison.**
- [ ] **Step 6: Commit.** Commit `feat: schedule balanced ESF increments`.

### Task 5E-3: ESF SoA State, Halo and Transaction

**Depends on:** 5F-2, accepted Stage 4 FieldRegistry/Halo.

**Files:**
- Create: `src/comb_esf_state.cpp`
- Create: `tests/unit/test_esf_state.cpp`
- Create: `tests/mpi/test_esf_state_halo.cpp`
- Modify: `src/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:** Produces `EsfState` construction/access/commit/rollback and field-major Halo
pack order for all stochastic components.

- [ ] **Step 1: Write RED.** Component-major address, N=2/4, owned/ghost extents, stable
  global IDs, bytewise rollback and 1/2-rank Halo equality.
- [ ] **Step 2: Run RED.** Expect missing state implementation.
- [ ] **Step 3: Implement via FieldRegistry.** Keep no raw owning allocation in the model;
  make chemistry packs scratch-only.
- [ ] **Step 4: Mutation check.** Transpose field/component order, omit last field, create
  field-owner communicator or publish trial on Halo completion; RED must fail.
- [ ] **Step 5: Run unit and 1/2-rank focused tests.**
- [ ] **Step 6: Commit.** Commit `feat: add transactional ESF state`.

### Task 5E-4: Unity-Lewis Deterministic Field Transport

**Depends on:** 5E-3, Stage 4 reacting transport.

**Files:**
- Create: `src/comb_esf_transport_detail.hpp`
- Create: `src/comb_esf_transport.cpp`
- Create: `tests/unit/test_esf_transport.cpp`
- Create: `tests/mpi/test_esf_transport_mpi.cpp`

**Interfaces:** Produces deterministic `advance_esf_transport(...)` using MeanState density,
final FaceMassFlux and one composition diffusivity.

- [ ] **Step 1: Write RED.** Constant/linear manufactured fields, final-flux identity,
  unity-Lewis species/enthalpy diffusion and N=2/4 independent field transport.
- [ ] **Step 2: Run RED.** Expect missing operator.
- [ ] **Step 3: Reuse Stage 4 FV kernels.** Add only field batching/handles; do not maintain
  an alternative flux or limiter.
- [ ] **Step 4: Mutation check.** Use predictor flux, independent enthalpy diffusivity,
  field-specific density or per-field PISO; RED must fail.
- [ ] **Step 5: Run unit plus small 1/2-rank tests.**
- [ ] **Step 6: Commit.** Commit `feat: transport stochastic fields`.

### Task 5E-5: Low-Reynolds Stochastic Transport Increment

**Depends on:** 5E-2, 5E-4.

**Files:**
- Modify: `src/comb_esf_transport.cpp`
- Create: `tests/unit/test_esf_stochastic_increment.cpp`
- Extend: `tests/mpi/test_esf_transport_mpi.cpp`

**Interfaces:** Adds exactly one stochastic increment stage consuming WALE `mu_t`, `Sc_t`,
spatial gradients and balanced increments.

- [ ] **Step 1: Write RED.** `mu_t=+0.0` bitwise zero, constant-field zero, linear-gradient
  oracle, same sign across cells, field-pair behavior and 1/2-rank equality.
- [ ] **Step 2: Run RED.** Expect deterministic-only result.
- [ ] **Step 3: Implement approved coefficient.** Molecular diffusion remains in the
  deterministic operator only; consume one scheduler output per stage.
- [ ] **Step 4: Mutation check.** Add molecular coefficient, include cell key, redraw inside
  PISO, use different signs per species or apply twice; RED must fail.
- [ ] **Step 5: Run focused unit/MPI and COAST source-free fixture.**
- [ ] **Step 6: Commit.** Commit `feat: add low-Re ESF stochastic transport`.

### Task 5E-6: Exact IEM Operator

**Depends on:** 5E-3, 5F-1.

**Files:**
- Create: `src/comb_esf_iem_detail.hpp`
- Create: `src/comb_esf_iem.cpp`
- Create: `tests/unit/test_esf_iem.cpp`

**Interfaces:** Produces `apply_iem(EsfAttemptState&, mean, dt, tau, kappa)`; `kappa=1`
uses `exp(-dt/(2*tau))` exactly.

- [ ] **Step 1: Write RED.** N=2/4 mean invariance, variance decay, `dt=0`, equal fields,
  `kappa=1`, invalid tau and component permutation.
- [ ] **Step 2: Run RED.** Expect missing mixer.
- [ ] **Step 3: Implement one affine update over all fields/components.** Accumulate mean in
  deterministic order and do not normalize each field independently.
- [ ] **Step 4: Mutation check.** Explicit Euler, missing half factor, field-specific mean,
  independent clipping or separate IEM/TCR branches; RED must fail.
- [ ] **Step 5: Run focused tests and COAST fixture comparison.**
- [ ] **Step 6: Commit.** Commit `feat: add exact IEM mixing`.

### Task 5E-7: Element and Ensemble-Moment Consistency

**Depends on:** 5E-4--5E-6, 5C-2. Main-agent math gate precedes implementation.

**Files:**
- Create: `src/comb_esf_element_consistency_detail.hpp`
- Create: `src/comb_esf_element_consistency.cpp`
- Create: `tests/unit/test_esf_element_consistency.cpp`
- Create: `docs/numerics/stage5-element-consistency-proof.md`

**Interfaces:** Produces a coupled projection/report preserving simplex, element matrix,
ensemble mean and `h_tc`; returns infeasible instead of clipping.

- [ ] **Step 1: Main agent derives the constrained operator.** Record variables,
  constraints, nullspace, tolerance budget and infeasible conditions from the public paper;
  freeze mutations before worker coding.
- [ ] **Step 2: Write RED.** Feasible N=2/4 perturbations, element/mean/enthalpy budgets,
  roundoff-only canonicalization and deliberately infeasible states.
- [ ] **Step 3: Run RED.** Expect missing projection.
- [ ] **Step 4: Implement the frozen solve only.** No species clipping, ad hoc redistribution
  or fuel-specific correction.
- [ ] **Step 5: Kill clipping/renormalization/element-row omission mutations and run focused
  sanitizer tests.**
- [ ] **Step 6: Main-agent commit.** Commit `feat: conserve ESF elements and moments`.

### Task 5C-1: Stochastic-Field Chemistry Batch Adapter

**Depends on:** 5E-3, Stage 4 ChemistryBackend.

**Files:**
- Create: `src/comb_esf_reaction_ensemble_detail.hpp`
- Create: `src/comb_esf_reaction_ensemble.cpp`
- Create: `tests/unit/test_esf_chemistry_batch.cpp`

**Interfaces:** Produces `integrate_esf_chemistry_half(...)` with deterministic
cell/field lane order and one thread-owned backend workspace per active lane.

- [ ] **Step 1: Write RED.** Scalar-vs-batch equivalence, N=2/4, field permutation,
  mechanism identity, half duration/stage time and one failing lane report.
- [ ] **Step 2: Run RED.** Expect missing adapter.
- [ ] **Step 3: Implement batching over the Stage 4 interface only.** No Cantera type or
  second chemistry implementation enters Stage 5.
- [ ] **Step 4: Mutation check.** Share workspace, reorder species, use full dt, ignore a
  failing lane or call direct Cantera; RED must fail.
- [ ] **Step 5: Run focused unit and small ASan.**
- [ ] **Step 6: Commit.** Commit `feat: integrate stochastic-field chemistry`.

### Task 5C-2: Reaction Ensemble and Mean Source Closure

**Depends on:** 5C-1, Stage 4 source transaction.

**Files:**
- Modify: `src/comb_esf_reaction_ensemble.cpp`
- Create: `tests/unit/test_esf_reaction_ensemble.cpp`

**Interfaces:** Produces `ReactionEnsembleReport` and
`add_ensemble_chemistry_source(ReactingSourceTransaction&, report)` using equal `1/N`
weights and integrated deltas.

- [ ] **Step 1: Write RED.** Identical fields, distinct fields, N=2/4 weights, total species
  mass/element budgets, no chemistry enthalpy source and mean transaction provenance.
- [ ] **Step 2: Run RED.** Expect missing ensemble closure.
- [ ] **Step 3: Implement deterministic reduction and transaction insertion.** Never copy
  arithmetic field mean over MeanState.
- [ ] **Step 4: Mutation check.** Field0 authority, `1/(N-1)`, endpoint rates, extra heat
  release or direct MeanState write; RED must fail.
- [ ] **Step 5: Run focused tests and Stage 4 transaction regression.**
- [ ] **Step 6: Commit.** Commit `feat: close mean chemistry from ESF ensemble`.

### Task 5C-3: PSR Shadow Reactor

**Depends on:** 5C-1, Stage 4 PSR/backend conformance.

**Files:**
- Create: `src/comb_esf_psr_shadow_detail.hpp`
- Create: `src/comb_esf_psr_shadow.cpp`
- Create: `tests/unit/test_esf_psr_shadow.cpp`

**Interfaces:** Produces `PsrShadowReport evaluate_psr_shadow(const MeanStateSnapshot&,
ChemistryBackend&)`; no transport/Halo/PISO/checkpoint state.

- [ ] **Step 1: Write RED.** Deterministic reconstruction, same mechanism/tolerance as PDF,
  reaction-off, failure mapping and no state change.
- [ ] **Step 2: Run RED.** Expect missing evaluator.
- [ ] **Step 3: Implement as a non-owning service call.** Rebuild from committed mean each
  step and return a value report.
- [ ] **Step 4: Mutation check.** Treat field0 as PSR, use separate mechanism, persist hidden
  reactor state or let PSR enter transport; RED must fail.
- [ ] **Step 5: Run focused unit and Restart reconstruction test.**
- [ ] **Step 6: Commit.** Commit `feat: add PSR shadow reaction`.

### Task 5T-0: Exact COAST TCR Authority Confirmation Gate

**Depends on:** 5F-0 and user input. This is a main-agent task.

**Files:**
- Create: `.superpowers/oracles/coast-tcr-source-manifest.md`
- No product changes.

**Interfaces:** Produces confirmed realpath, commit, manifest, theory status, exact
allowlist and hashes. Product citations are already frozen as
`10.1016/j.proci.2026.106128` and `10.1016/j.cja.2026.104123`; a candidate COAST path is
not oracle authority until the user confirms it in this task.

- [ ] **Step 1: Main agent presents candidate source identity to the user.** Include realpath, commit,
  manifest status and hashes; do not read/copy implementation before confirmation.
- [ ] **Step 2: Record explicit source confirmation and resolved symlink paths.** Verify
  the manifest uses the two frozen TCR DOIs; do not infer or replace citations from a
  private source comment.
- [ ] **Step 3: Verify allowlist only.** Root transport, source control, operator budget,
  branch transport, finite field, safe kappa and Wiener counter pure modules; block driver,
  restart, cases, data and executable.
- [ ] **Step 4: Record copyright/access classification.** Private test oracle only.
- [ ] **Step 5: Run tracked-text/private-path scan.**
- [ ] **Step 6: Commit manifest.** Commit `docs: confirm COAST TCR oracle authority`.

### Task 5T-1: Standalone COAST TCR Oracle Protocol

**Depends on:** 5T-0.

**Files:**
- Extend: `tests/oracles/coast_oracle_protocol.hpp`
- Create: `tests/oracles/coast_tcr_driver.F90`
- Extend: `tests/oracles/coast_oracle_runner.cpp`
- Create: `tests/cmake/coast_tcr_oracle_contract.cmake`

**Interfaces:** Versioned binary/text protocol for synthetic TCR input/history/injected
signs and exact status/output; separate executable, no ABI link.

- [ ] **Step 1: Write malformed-protocol RED.** Wrong version, length, NaN, history count
  and status must be rejected.
- [ ] **Step 2: Implement HUNDUN-owned driver.** It calls only generated allowlisted COAST
  modules; source remains byte-for-byte and outside tracked tree.
- [ ] **Step 3: Build and hash oracle.** Record compiler/options/source/driver/binary.
- [ ] **Step 4: Run synthetic direct/fold/denominator/invalid/fallback vectors.**
- [ ] **Step 5: Run source/package leak scans.**
- [ ] **Step 6: Commit protocol/driver/contract only.** Commit
  `test: add standalone COAST TCR oracle`.

### Task 5T-2: Independent TCR Algebra and Status Kernel

**Depends on:** 5F-1. May run before 5T-1 because it derives from papers/spec, not COAST
source.

**Files:**
- Create: `src/comb_tcr_algebra_detail.hpp`
- Create: `src/comb_tcr_algebra.cpp`
- Create: `tests/unit/test_tcr_algebra.cpp`
- Modify: `src/CMakeLists.txt`, `tests/CMakeLists.txt`

**Interfaces:** Produces `evaluate_tcr_algebra(const TcrAlgebraInput&)` and exact status
classification.

- [ ] **Step 1: Write RED.** Direct plus/minus branch, inactive, `D<0`, near fold, weak
  denominator, non-finite and sign-history cases.
- [ ] **Step 2: Run RED.** Expect missing kernel.
- [ ] **Step 3: Implement approved equations with explicit preconditions/status.** Add a
  concise product-source comment citing `10.1016/j.proci.2026.106128` and
  `10.1016/j.cja.2026.104123`, and state that the implementation is independently derived.
  Never cite a private path and never clamp D or eta.
- [ ] **Step 4: Mutation check.** Flip sign, use current state for branch, clamp negative D,
  merge weak denominator with invalid or adjust eta; RED must fail.
- [ ] **Step 5: Run unit/UBSan and source-free oracle vectors when available.**
- [ ] **Step 6: Commit.** Commit `feat: add independent TCR algebra`.

### Task 5T-3: Root History, Source Control and Finite-Field Statistics

**Depends on:** 5T-2, 5C-2, 5C-3, 5F-2.

**Files:**
- Create: `src/comb_tcr_root_state_detail.hpp`
- Create: `src/comb_tcr_root_state.cpp`
- Create: `tests/unit/test_tcr_root_state.cpp`

**Interfaces:** Produces staged `TcrRootState`, nine-value history transition,
`TcrSourceControlReport` and N-aware statistics status.

- [ ] **Step 1: Write RED.** History transition, source control changes R only, N=2
  minimal-pair status, N=4 regular path, fold event certificate and retry rollback.
- [ ] **Step 2: Run RED.** Expect missing root state/controller.
- [ ] **Step 3: Implement versioned committed/trial root state.** Keep scratch guesses out of
  persistence.
- [ ] **Step 4: Mutation check.** Modify eta, infer sign from current D, consume only eight
  history values, accept uncertified crossing or advance history on reject; RED must fail.
- [ ] **Step 5: Run focused unit and transaction tests.**
- [ ] **Step 6: Commit.** Commit `feat: evolve TCR root state`.

### Task 5T-4: COAST Differential Replay and Feedback-Timing Proof

**Depends on:** 5T-1--5T-3, 5E-R2, 5E-7, 5C-2/3. Main-agent task.

**Files:**
- Create: `tests/fixtures/tcr_oracle_vectors.txt`
- Create: `tests/unit/test_tcr_coast_replay.cpp`
- Create: `docs/numerics/stage5-tcr-feedback-timing-proof.md`

**Interfaces:** Freezes one schedule: accepted-state lag or one fixed
predictor/evaluate/correct. Produces the exact comparison tolerance/status policy.

- [ ] **Step 1: Generate point and multi-step state replay with injected signs.** Use N=4;
  include every status and history field.
- [ ] **Step 2: Compare candidate timing schedules.** Require state/status/effect agreement;
  reject unbounded inner iteration.
- [ ] **Step 3: Freeze the selected schedule and rejected alternative in the proof.**
- [ ] **Step 4: Add mutations for lag index, stage ordering, sign redraw and fallback reuse.**
- [ ] **Step 5: Run source-free fixtures plus optional live oracle; no Vblowoff.**
- [ ] **Step 6: Main-agent commit.** Commit `docs: freeze TCR feedback timing`.

### Task 5T-5: TCR Mode Controller and Promotion Gate

**Depends on:** 5T-4.

**Files:**
- Create: `src/comb_tcr_controller_detail.hpp`
- Create: `src/comb_tcr_controller.cpp`
- Create: `tests/unit/test_tcr_controller.cpp`

**Interfaces:** Produces `apply_tcr_mode(TcrMode,...)`; default shadow, validated feedback
requires an accepted oracle identity.

- [ ] **Step 1: Write RED.** off/unity, shadow no feedback, experimental labeled feedback,
  validated gating, mode-stable RNG and status counters.
- [ ] **Step 2: Run RED.** Expect missing controller.
- [ ] **Step 3: Implement explicit mode switch over one algebra/mixer path.** No duplicated
  IEM implementation.
- [ ] **Step 4: Mutation check.** Default validated, shadow writes kappa, mode changes RNG
  stream or ignores oracle mismatch; RED must fail.
- [ ] **Step 5: Run focused unit and replay tests.**
- [ ] **Step 6: Commit.** Commit `feat: gate TCR feedback modes`.

### Task 5I-1: TPDF/TCR Step Coordinator

**Depends on:** 5E-5--5E-7, 5C-2/3, 5T-5, Stage 4 coupling hooks. Main-agent task.

**Files:**
- Create: `src/flow_tpdf_tcr_coupling_detail.hpp`
- Create: `src/flow_tpdf_tcr_coupling.cpp`
- Create: `tests/unit/test_tpdf_tcr_coupling.cpp`

**Interfaces:** Produces `attempt_tpdf_tcr_step(...)` and `TpdfTcrStepReport` with fixed
chemistry/transport/noise/mixing/PSR/TCR/PISO order.

- [ ] **Step 1: Write operator-spy RED.** Verify field C1, deterministic/stochastic
  transport once, IEM/TCR mixing, field C2, ensemble source, PSR/TCR candidate and PISO #2
  final visibility.
- [ ] **Step 2: Run RED.** Expect missing coordinator.
- [ ] **Step 3: Implement one atomic trial using approved feedback schedule.**
- [ ] **Step 4: Mutation check.** Reorder PSR/rates, apply noise twice, overwrite MeanState,
  commit root before PISO or add coupling iteration; RED must fail.
- [ ] **Step 5: Run unit and <=12^3 N=2 smoke with the Stage 4
  `tests/support/chem_analytic_backend.*` fixture.** The fixture is not a product backend.
- [ ] **Step 6: Main-agent commit.** Commit `feat: coordinate TPDF and TCR`.

### Task 5I-2: IBM/WALE and Reacting Boundary Integration

**Depends on:** 5I-1, Stage 4 IBM+WALE reacting path.

**Files:**
- Create: `src/flow_tpdf_tcr_immersed_detail.hpp`
- Create: `src/flow_tpdf_tcr_immersed.cpp`
- Create: `tests/mpi/test_tpdf_tcr_immersed.cpp`

**Interfaces:** Extends field transport to non-catalytic body-fitted/IBM walls with shared
WALE/pressure/final-flux authority.

- [ ] **Step 1: Write RED.** Impermeable stochastic species, adiabatic/isothermal h,
  `mu_t=0` wall limit, WALE once per trial and N=2/4 field identities.
- [ ] **Step 2: Run RED.** Expect missing integration.
- [ ] **Step 3: Implement thin adapters over Stage 4/5 operators.**
- [ ] **Step 4: Mutation check.** Per-field WALE, separate IBM rows, catalytic wall or
  stochastic wall RNG by cell; RED must fail.
- [ ] **Step 5: Run focused 1/2-rank and <=12^3 IBM smoke.**
- [ ] **Step 6: Commit.** Commit `feat: couple TPDF to IBM and WALE`.

### Task 5I-3: Atomic Rollback, Collective Failure and MPI

**Depends on:** 5I-1, 5I-2.

**Files:**
- Create: `tests/mpi/test_tpdf_tcr_rollback.cpp`
- Create: `tests/mpi/test_tpdf_tcr_decomposition.cpp`
- Modify: `src/flow_tpdf_tcr_coupling.cpp`

**Interfaces:** Finalizes collective failure mapping and rollback of MeanState, all fields,
root history and RNG clock.

- [ ] **Step 1: Inject failure at each stage.** Field C1, transport, stochastic increment,
  mixing, C2, PSR, TCR, element projection and PISO #2.
- [ ] **Step 2: Require bytewise state/history/counter restoration and all-rank agreement.**
- [ ] **Step 3: Implement missing rollback hooks only.** Retry reuses signs with new
  magnitude and reruns the complete attempt.
- [ ] **Step 4: Mutation check.** Advance accepted step, preserve trial root, redraw signs or
  let one rank return early; tests must fail/hang-detect.
- [ ] **Step 5: Run small 1/2/4-rank focused tests.** No 24^3.
- [ ] **Step 6: Commit.** Commit `fix: make TPDF trials collectively atomic`.

### Task 5I-4: Schema v5, Broadcast and Driver

**Depends on:** 5F-3, 5I-1--5I-3.

**Files:**
- Create: `src/cfg_resolved_case_v5_loader.cpp`
- Create: `src/cfg_resolved_case_v5_loader_detail.hpp`
- Create: `src/app_resolved_case_v5_broadcast.cpp`
- Create: `src/app_tpdf_tcr_driver.cpp`
- Create: `src/app_tpdf_tcr_driver_detail.hpp`
- Create: `tests/unit/test_resolved_case_v5.cpp`
- Create: `tests/mpi/test_resolved_case_v5_broadcast.cpp`
- Create: `tests/unit/test_tpdf_tcr_driver.cpp`

**Interfaces:** Produces `ResolvedTpdfTcrCaseV5` loader/broadcast and
`run_tpdf_tcr_case(...)` in the same `hundun` executable.

- [ ] **Step 1: Write RED.** N=2/4, seed, modes, tau/kappa parameters, `Pr_t=Sc_t`, oracle
  identity, illegal odd N/validated-without-oracle and old-schema compatibility.
- [ ] **Step 2: Run RED.** Expect missing implementation.
- [ ] **Step 3: Implement rank-0 parse/broadcast and composition root.** Reuse Stage 4
  mechanism/backend and existing root dispatch registration.
- [ ] **Step 4: Mutation check.** Default N=8, allow mismatched diffusivity, hidden field0 or
  direct Cantera construction; RED must fail.
- [ ] **Step 5: Run unit, 1/2-rank broadcast, validate/print-resolved and one-step smoke.**
- [ ] **Step 6: Commit.** Commit `feat: add TPDF-TCR schema and driver`.

### Task 5I-5: Checkpoint v5

**Depends on:** 5F-3, 5I-3/4.

**Files:**
- Create: `src/flow_checkpoint_v5.cpp`
- Create: `src/flow_checkpoint_v5_detail.hpp`
- Create: `tests/unit/test_checkpoint_v5.cpp`
- Create: `tests/mpi/test_checkpoint_v5_mpi.cpp`

**Interfaces:** Registers sections for N/layout, all fields, RNG algorithm/seed/accepted
clock, TCR mode/root history and Stage 4 mechanism identity. PSR/workspaces/caches excluded.

- [ ] **Step 1: Write RED.** N=2/4 continuous-vs-restart, corruption, mode/oracle mismatch,
  failed-read unchanged and same-partition 1/2/4 continuation.
- [ ] **Step 2: Run RED.** Expect missing codecs.
- [ ] **Step 3: Implement validate-then-publish section codecs.** Reconstruct PSR after read.
- [ ] **Step 4: Mutation check.** Omit field, save mutable RNG cursor, persist PSR workspace,
  accept N/mode mismatch or partial publish; RED must fail.
- [ ] **Step 5: Run focused unit/MPI and retry-after-Restart.**
- [ ] **Step 6: Commit.** Commit `feat: checkpoint TPDF-TCR state v5`.

### Task 5I-6: Stage 5 Diagnostics and Exact Counters

**Depends on:** 5F-3, 5I-4/5.

**Files:**
- Create: `src/diag_esf_tcr.cpp`
- Create: `tests/unit/test_esf_tcr_diagnostics.cpp`
- Modify: integration-owned diagnostic registry

**Interfaces:** Providers expose N/sampling label, RNG identity, mean/variance, pre/post
element budgets, chemistry/PSR rates, every TcrStatus, kappa/root history and rollback counts.

- [ ] **Step 1: Write RED.** Stable IDs/units, N=2 label, no fake zero, read-only callbacks,
  no hidden collective and exact stage/counter identities.
- [ ] **Step 2: Run RED.** Expect missing providers.
- [ ] **Step 3: Implement report adapters only.** Do not recompute means/rates/root.
- [ ] **Step 4: Mutation check.** Collapse statuses, omit sampling label, mutate state or use
  unordered output; RED must fail.
- [ ] **Step 5: Run unit plus 1/2-rank diagnostic session.**
- [ ] **Step 6: Commit.** Commit `feat: diagnose ESF and TCR`.

### Task 5A-1: ESF/IEM Milestone Gate

**Depends on:** 5E-1--5E-7, 5C-1/2, 5E-R2.

**Files:**
- Create: `tests/acceptance/stage5_esf_gate.sh`
- Create: `.superpowers/sdd/stage5-5A-1-esf-gate.md`

**Interfaces:** Produces N=2/4 ESF/IEM acceptance evidence without TCR feedback.

- [ ] **Step 1: Freeze selectors.** Philox/Wiener, laminar limit, transport, exact IEM,
  chemistry ensemble, element consistency, COAST source-free fixtures, Restart/rollback and
  small 1/2/4-rank.
- [ ] **Step 2: Assert exclusions.** No Vblowoff, 48^3, mechanism science sweep or long MPI.
- [ ] **Step 3: Run Debug/focused Release/small sanitizer evidence.**
- [ ] **Step 4: Verify N=2 is labeled minimal and N=4 is transient recommendation.**
- [ ] **Step 5: Main-agent full ESF diff/science review.**
- [ ] **Step 6: Commit receipt.** Commit `test: accept Stage 5 ESF core`.

### Task 5A-2: TCR Oracle and Feedback Milestone

**Depends on:** 5T-0--5T-5, 5C-3, 5A-1.

**Files:**
- Create: `tests/acceptance/stage5_tcr_gate.sh`
- Create: `.superpowers/sdd/stage5-5A-2-tcr-gate.md`

**Interfaces:** Produces N=4 COAST differential evidence and N=2 HUNDUN steady canary.

- [ ] **Step 1: Run source-free point/state/status/root-history fixtures.** Optionally rerun
  live oracle only when exact source/binary hashes match the manifest.
- [ ] **Step 2: Run off/shadow/experimental/validated mode tests.**
- [ ] **Step 3: Run N=2 steady/time-average canary and confirm limited capability wording.**
- [ ] **Step 4: Run feedback timing/operator spy and rollback.**
- [ ] **Step 5: Main-agent TCR math/provenance/full-diff review.** No Vblowoff.
- [ ] **Step 6: Commit receipt.** Commit `test: accept TCR oracle equivalence`.

### Task 5A-3: Integrated Stage 5 Framework Gate

**Depends on:** 5I-1--5I-6, 5A-1/2.

**Files:**
- Create: `tests/acceptance/stage5_acceptance.sh`
- Create: `.superpowers/sdd/stage5-5A-3-framework-gate.md`
- Modify: `VERSION`
- Modify: `docs/numerics/stage5-capability-ledger.md`

**Interfaces:** Produces frozen Stage 5 code candidate `C5` at version `0.4.0` and the
complete low-cost Stage 5 evidence manifest.

- [ ] **Step 1: Finalize and freeze candidate.** Set `VERSION=0.4.0`, finish all
  product/test source and selectors, run a low-cost preflight, create signed `C5`, and
  record HEAD/tree/binary/dirty hashes. Do not modify product or tests afterward.
- [ ] **Step 2: Run 0D, 1D and one <=12^3 IBM+WALE reacting smoke for N=2/4.**
- [ ] **Step 3: Run small 1/2/4-rank, Checkpoint, diagnostics, driver, headers and tests-off
  only where affected.**
- [ ] **Step 4: Run complete affected Debug and focused Release/ASan/UBSan on `C5`.**
- [ ] **Step 5: Hash every command/binary/log, verify no high-load process remains, and
  complete the main-agent requirements/quality/caller/provenance review.** A product fix
  creates a new `C5` and invalidates only consuming evidence.
- [ ] **Step 6: Commit governance receipt.** Commit only the evidence manifest/ledger as
  `test: record Stage 5 acceptance evidence`, with `accepted_code_head=C5`.

### Task 5A-4: Stage 5 Exact-HEAD Seal

**Depends on:** 5A-3.

**Files:**
- Create: `.superpowers/sdd/stage5-final-acceptance-report.md`
- Modify: `docs/numerics/stage5-capability-ledger.md`
- Modify: `docs/numerics/stage4-6-capability-root.md`
- Modify: `AGENTS.md`

**Interfaces:** Produces governance seal `G5` for tested code candidate `C5` at `0.4.0`;
no default product projection.

- [ ] **Step 1: Recover `C5` and record parent/tree/diff/binary/log/oracle hashes.** Assert
  all later commits are governance-only.
- [ ] **Step 2: Verify DCO, COAST source exclusion, mechanism/Cantera identity and worktree/
  background-process state.**
- [ ] **Step 3: State capability limits.** No Vblowoff, no field MPI decomposition, N=2
  limited sampling and no Stage 6 spray.
- [ ] **Step 4: Record governance receipt version `0.4.0`; do not modify `VERSION`, product
  code or tests, and do not sync product repo.**
- [ ] **Step 5: Run final read-only manifest and `git diff --check` audit.**
- [ ] **Step 6: Main-agent decision/governance commit.** Write `STAGE5_ACCEPT` or
  `REJECT`; commit `docs: seal Stage 5 ESF TPDF TCR` only when `C5` is accepted. Record
  `accepted_code_head=C5`; do not call the governance commit the tested code HEAD.

## 4. Critical Path and Serial Gate

```text
5F-0 -> 5F-1 -> 5F-2/3
5F-2 -> 5E-1 -> 5E-2
5F-2 -> 5E-3 -> 5E-4 -> 5E-5
5E-3 -> 5E-6
5E-3 -> 5C-1 -> 5C-2/3
5E-5/6 + 5C-2 -> 5E-7
5T-0 -> 5T-1
5F-1 -> 5T-2 -> 5T-3
5E-R1 -> 5E-R2
5T-1/3 + 5E-R2 + 5E-7 -> 5T-4 -> 5T-5
5E-7 + 5C-2/3 + 5T-5 -> 5I-1 -> 5I-2/3/4/5/6
5A-1 -> 5A-2 -> 5A-3 -> 5A-4
```

默认串行执行本图。Stage 5 接受后，主 agent在进入 Stage 6 前向用户报告是否值得
并行启动尚未完成的 Stage 6 pure-parcel 分支；没有明确指示则继续串行。

## 5. Task Reference Matrix

URL、DOI、license class 和 COAST C-level 边界以
`docs/references/2026-08-09-hundun-flow-stage4-6-reference-catalog.md` 为准。公开论文控制
方程；COAST 只提供用户确认后的进程外差分 oracle，不是产品实现模板。

| Task | Reference | Exact point to reuse | Explicitly avoid |
|---|---|---|---|
| 5F-0 | accepted Stage 4 receipt/services | exact backend, transaction and registry intake | starting from an unaccepted Stage 4 tree |
| 5F-1 | Valiño；Valiño--Mustata--Letaief；TCR two DOIs | low-Re ESF, IEM limit, TCR equations/status vocabulary | using COAST path as citation |
| 5F-2 | Random123 paper；Stage 4 transaction | stateless address, retry identity and staged state | mutable RNG cursor or per-cell key |
| 5F-3 | accepted v4 registries | append-only v5 schema/checkpoint/diagnostic IDs | changing v1--v4 identities |
| 5E-R1 | confirmed COAST ESF；Valiño papers | observe persistence/timing and classify provenance | copying control flow or research cases |
| 5E-R2 | confirmed COAST ESF C-oracle protocol | synthetic one-step differential records | linking private code into `hundun` |
| 5E-1 | Salmon et al. Philox | ten-round counter algorithm and published vectors | copying Random123 source text or stateful engine |
| 5E-2 | Valiño ESF；counter RNG | antithetic N=2/4 Wiener pairs and domain separation | cell/species/rank/PISO keys |
| 5E-3 | HUNDUN FieldRegistry/Halo/transaction | component-major SoA and all fields per spatial rank | field-index MPI decomposition in v1 |
| 5E-4 | Valiño ESF；accepted Stage 4 FV kernels | unity-Lewis deterministic transport on final flux | alternative limiter/flux/PISO |
| 5E-5 | Valiño--Mustata--Letaief | turbulent-only stochastic coefficient and laminar limit | molecular diffusion in noise |
| 5E-6 | IEM literature；TCR `kappa=1` contract | exact affine exponential mixing | explicit Euler or separate TCR mixer |
| 5E-7 | Xu et al. finite-field consistency | coupled simplex/element/mean constraints | fieldwise clipping and renormalization |
| 5C-1 | Stage 4 ChemistryBackend；Cantera black box | deterministic field batching through backend service | Cantera types in Stage 5 |
| 5C-2 | PDF ensemble definition；Day--Bell integrated source | equal-weight integrated mean chemistry closure | field 0 or endpoint-rate authority |
| 5C-3 | Stage 4 PSR conformance；Cantera backend | stateless committed-mean shadow reactor | persistent hidden PSR state |
| 5T-0 | TCR two DOIs；confirmed COAST | exact theory citation and private oracle identity | guessed source root or inferred citation |
| 5T-1 | confirmed COAST TCR C-oracle | versioned synthetic numeric protocol | ABI linking, case/data access |
| 5T-2 | TCR two DOIs | independent algebra, branches and status kernel | translating private source or clamping discriminant |
| 5T-3 | TCR two DOIs；finite-field statistics | nine-value history, source-control and N-aware status | advancing history on rejected trial |
| 5T-4 | TCR two DOIs；COAST differential replay | choose and freeze feedback timing from state/status evidence | unbounded inner iteration or Vblowoff |
| 5T-5 | TCR papers and accepted oracle identity | off/shadow/experimental/validated promotion gate | validated default or RNG stream changes by mode |
| 5I-1 | Strang/Day--Bell；TCR schedule | fixed chemistry/transport/noise/mixing/PISO order | hidden coupling iteration or third PISO |
| 5I-2 | accepted Stage 3 IBM/WALE and Stage 4 reacting walls | share geometry, pressure, WALE and final flux | per-field IBM/WALE authorities |
| 5I-3 | MPI standard；HUNDUN transaction model | collective decision and bytewise rollback | rank-local early return or retry redraw |
| 5I-4 | yyjson/MPI broadcast；Stage 4 composition root | typed v5 schema and one executable dispatch | direct Cantera construction or old-schema mutation |
| 5I-5 | accepted Checkpoint v4 protocol | all-field/RNG/root-history validate-then-publish | PSR/workspace persistence |
| 5I-6 | accepted diagnostics provider model | read-only final reports and exact counters | recomputed statistics/rates in diagnostics |
| 5A-1 | ESF papers；source-free COAST fixtures | compact N=2/4 laminar/IEM/consistency gate | long combustion case |
| 5A-2 | TCR papers；confirmed COAST oracle | point/state/status/root feedback equivalence | Vblowoff or mechanism sweep |
| 5A-3 | cross-stage compact gate | exact `C5`, small MPI, <=12^3 and focused sanitizers | candidate mutation after evidence |
| 5A-4 | cross-stage exact-HEAD protocol；DCO | `C5`/`G5` separation and capability limitations | calling governance commit the tested code head |

## 6. Receipt Requirements

每个 task receipt 除通用 HEAD/tree/diff/command/log/DCO 信息外，必须记录：

```text
N and sampling label
RNG algorithm/key identity
accepted-step/retry semantics
MeanState/field/PSR/TCR authority impact
element/moment budgets
TCR status and root-history impact
COAST source/oracle manifest identity（若消费）
source/package leak scan
deferred scientific claims
```

worker 只可实现冻结的纯核、codec、fixture 或 adapter。main agent保留 ESF 方程、
element correction、TCR root science、feedback timing、组合时序和最终接受。

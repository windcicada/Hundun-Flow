# HUNDUN-FLOW Stage 3 Parallel Completion v2 Implementation Plan

> **Status: `PROPOSED_DO_NOT_EXECUTE`.** 本文件作为不可变候选保持该标记；只有用户
> 批准后执行 S3-A0、把本文件及 design/reference 的 SHA-256 写入 tracked activation
> receipt 并原子切换 ledger/AGENTS，后续 task 才获得执行 authority。A0 不回写本文件，
> 因此 receipt 所绑定的 hash 不会自我失效。

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 从已接受的 constant IBM、Checkpoint v3 IBM-only、minimal diagnostics 和
body-fitted WALE 基线，完成三种密度的 IBM/WALE 组合、Restart、diagnostics、driver、
capability ledger、最终科学验收和 product 0.2.0 投影。

**Architecture:** `FixedStepImmersedFlow` 继续独占 immersed momentum/pressure/force；
一个构造期冻结的 density adapter 提供 attempt density/transport/closure，一个
attempt-local WALE authority 提供唯一 `mu_eff`。科学组合在主 agent lane 顺序合入；
Checkpoint constant-profile codec 和 WALE diagnostics provider 在独立 worker lane 完成，
涉及 density/closure/driver/counters 的其余集成仍由主 agent完成。最终长测只在软件
code-complete 并冻结 candidate 后运行。

**Tech Stack:** C++20、CMake 3.21+、Clang 15/libc++、OpenMPI、CPU-reference execution、
matrix-free FVM/PISO、schema v3、Checkpoint v3、CTest、Bash acceptance runners。

## Global Constraints

- 开发工作树固定为 `/home/wyf/code_dev/.worktrees/hundun-flow-stage3-framework`。
- 起始 accepted HEAD 固定为 `7fc8c5080528f6ea0dbc787c51ca40d9e0fa4553`；每个 packet
  在派发前写入实际 parent HEAD。
- governance 是唯一开发仓库；product 只在 Stage 3 最终接受后投影一次。
- 永久禁止 96-cubed；48-cubed 只在 code-complete frozen candidate 上运行。
- 不修改 Task 11 的科学阈值、force sign、authority、PISO corrector 数量或正式
  selector；发现影响时停止并交回主 agent。
- 不引入 AMReX、AMReX-Hydro、incflo、OpenFOAM、Basilisk、gslib、PETSc、Trilinos
  运行时依赖；只按 reference 文档做语义复用。
- 不复制、翻译或机械改写上游源码、注释、控制流、宏、ABI、错误文本或数据布局。
- 不访问私有 BOFFIN/COAST/COAST-2 或研究数据，不检查或干扰研究进程。
- worker 不提交、不添加 DCO、不改任务目标、不扩文件范围、不启动后台/长测试。
- 同时最多一个 implementation worker；只有独立 linked worktree/build tree 建立后，
  主 agent 才可与它并行处理科学 lane。
- 每个 task 都按 RED → minimum implementation → fast → main full-diff review → signed
  commit；不为放心追加未冻结测试。
- 任何 product/test source 改动后必须先更新 codegraphf index，再由主 agent审查
  callers/impact；精确文本搜索使用 `rg`。
- 每个“写 RED”步骤同时在 P0 分配的 task registration fragment 注册可执行程序和
  CTest name；P0 后任何剩余 task 均不得再修改 `tests/CMakeLists.txt`。unit
  timeout 60 秒，1/2-rank fast timeout 120 秒，12-cubed screen timeout 300 秒，labels
  至少含 `red;stage3;<task-id>;fast`。首次 RED 命令必须能实际启动并以预期 assertion/
  unsupported-path 原因失败，不能只因 target 未注册或编译语法错误失败。

---

## Phase A0 — Approval-bound activation transaction

### Task S3-A0: Activate the immutable v2 candidate

**Ownership:** main agent only；这是用户批准后的治理事务，不委派、不修改产品/测试。

**Consumes:** explicit user approval、clean candidate-doc commit and the exact bytes of the v2
design、reference and plan。

**Produces:** one tracked activation receipt and an atomic AGENTS/ledger authority switch；the
three hashed candidate documents remain byte-identical。

**References:** Git/DCO/exact-hash governance only；不采用外部数值算法。

**Files:**

- Modify: `AGENTS.md`
- Modify: `.superpowers/sdd/2026-08-09-hundun-flow-stage3-framework-completion/ledger.md`
- Create:
  `.superpowers/sdd/2026-08-09-hundun-flow-stage3-framework-completion/stage3-v2-activation.md`

**Non-TDD boundary:** A0 不实现行为，因此不伪造 product RED。它以 hash/status mutation
check 代替：任一文档 byte、批准状态、active profile 或 superseded row 被改错，contract
必须失败。

- [ ] **Step 1: freeze approval inputs**

要求 worktree clean、HEAD 等于用户审阅的 candidate commit。记录批准消息的原文、时间、
candidate commit、parent，并运行：

```bash
git status --short
git rev-parse HEAD HEAD^
sha256sum \
  docs/superpowers/specs/2026-08-09-hundun-flow-stage3-parallel-completion-v2-design.md \
  docs/references/2026-08-09-hundun-flow-stage3-public-algorithm-reference.md \
  docs/superpowers/plans/2026-08-09-hundun-flow-stage3-parallel-completion-v2.md
```

- [ ] **Step 2: write the activation transaction**

Receipt 保存上述三个完整 hash、批准原文和 `active_profile=stage3-parallel-framework-v2`。
`AGENTS.md` 把 v2 从 candidate 改为 active，但保留 Task 11/Stage 2 authority；ledger 把
legacy pending rows 标成 `superseded-by-stage3-parallel-framework-v2`，并把 v2 rows 从
proposed 改成 planned。三份 hashed candidate docs 不得修改。

- [ ] **Step 3: run mutation-sensitive governance checks**

重新计算三个 hash并逐字比对 receipt；用 `rg` 证明恰有一个 active profile、所有 legacy
pending rows 已 superseded、A0 receipt 未授权 Stage 4/luna/private access/push。临时把
receipt 中任一 hash 改一个 hex digit 的 in-memory/fixture check 必须返回 nonzero；真实
文件检查必须返回 zero。最后运行 `git diff --check`。

- [ ] **Step 4: main review and signed activation commit**

完整审查实际 diff 只能含上述三项。主 agent创建签署 commit
`docs: activate Stage 3 parallel completion v2`。只有该 commit clean 后才进入 P0。

---

## Worker packet 固定格式

每次派发必须把以下内容连同具体 task 一起发给 worker；worker 不需要读取整个项目历史：

```text
repository: the task-specific linked worktree from the task's Worker Packet block
build_tree: the task-specific build tree from the task's Worker Packet block
baseline_head: the exact 40-hex output of `git rev-parse HEAD` immediately before dispatch
task_id: exactly one S3 task ID copied from this plan
allowed_files: the narrower worker-only list copied from that task's Worker Packet block;
               never copy a main-only task's complete Files list
forbidden: commit, sign-off, threshold changes, extra corrector, long tests,
           private sources, unrelated cleanup, public API invention
required_reading_after_activation:
  AGENTS.md
  2026-08-09 stage3-parallel-completion-v2 design
  2026-08-09 public-algorithm-reference
  this task block
deliverable:
  changed-file list
  RED command and observed failure
  GREEN commands and exit status
  assumptions and unresolved blockers
  no commit
```

只有标有 **Worker Packet** 的 task 可以派发；main-only task 不得使用本模板。主 agent
在 packet 指定 worktree 派发前运行：

```bash
git status --short
git rev-parse HEAD
codegraphf sync .
```

若 worktree 不干净但修改正属于同一 active task，可继续；否则不得让 worker覆盖它。
worker 不提交。worker 返回后，主 agent 在该 isolated worktree 做完整 diff/调用方审查、
运行 GREEN 并创建签署 handoff commit；集成到 main 后再次运行受影响 GREEN。receipt
同时记录 `worker_baseline_head`、`handoff_commit` 和 `integration_parent`。

---

## Phase F0.5 — Parallel execution foundation

### Task S3-P0: Registration shards and isolated infrastructure lane

**Ownership:** main agent only。纯治理/测试构建图任务；不修改产品数值代码。

**Consumes:** signed S3-A0 activation commit whose parent chain contains accepted HEAD
`7fc8c5080528f6ea0dbc787c51ca40d9e0fa4553`。

**Produces:** 后续新测试的五个 registration fragments；可验证的 worker reading
exemption；独立 infra branch/worktree/build contract。

**References:** Git linked-worktree isolation 和现有 CTest helper；不采用外部数值算法。

**Files:**

- Modify: `tests/CMakeLists.txt`
- Create: `tests/cmake/stage3_science_registration.cmake`
- Create: `tests/cmake/stage3_checkpoint_registration.cmake`
- Create: `tests/cmake/stage3_diagnostics_registration.cmake`
- Create: `tests/cmake/stage3_framework_registration.cmake`
- Create: `tests/cmake/stage3_acceptance_registration.cmake`
- Create: `tests/cmake/stage3_registration_contract.cmake`
- Create: `tests/cmake/stage3_registration_contract_fixture.cmake`

**Registration ownership:**

| Fragment | Exclusive writers |
| --- | --- |
| `stage3_science_registration.cmake` | C1、D1、C2、D2、C3、S1 |
| `stage3_checkpoint_registration.cmake` | R1、R2 |
| `stage3_diagnostics_registration.cmake` | O1、O2 |
| `stage3_framework_registration.cmake` | A1、E1 |
| `stage3_acceptance_registration.cmake` | G1、DOC contract、V0/V1 inventory |

- [ ] **Step 1: write an executable registration RED**

先创建 validator/fixture，并只在 `tests/CMakeLists.txt` 注册：

```text
test_stage3_registration_contract
test_stage3_registration_contract_mutation
```

validator 要求五个 exact include 各出现一次、对应 fragment 存在且有 SPDX/include
guard；fixture 只复制 `tests/CMakeLists.txt` 和 fragments 到 build-tree sandbox，绝不改
source。此时尚未新增五个 include/fragments，运行：

```bash
cmake --preset debug
ctest --test-dir build/debug -R '^test_stage3_registration_contract$' \
  --output-on-failure
```

Expected: FAIL with `missing stage3 registration include`，而不是 CMake syntax/target failure。

- [ ] **Step 2: add fragments/includes and make the mutation contract GREEN**

`tests/CMakeLists.txt` 在现有 `if(HUNDUN_BUILD_TESTS)` 内、helper definitions 之后只
include 五个 fragment。fragment 初始只含 SPDX、include guard 和注释，不移动、不重命名
任何 accepted test。mutation fixture 在 build-tree copy 中删除一项 include，要求
validator nonzero 并把“检测到 mutation”作为 fixture PASS；不得修改/删除真实 source。
未 mutation 的真实 tree 必须通过。

- [ ] **Step 3: configure/build graph GREEN**

```bash
cmake --preset debug
cmake --build build/debug -j32 --target test_stage3_flow_header_contract
ctest --test-dir build/debug -R \
  '^(test_stage3_registration_contract(_mutation)?|source_layout_fixture|cmake_include_authority_fixture|test_stage3_flow_header_contract)$' \
  --output-on-failure
```

不运行 MPI 或数值 selector。

- [ ] **Step 4: sign P0 and freeze the post-C1 worktree command**

主 agent签署 P0。Do not create the lane from the pre-C1 HEAD。C1 Step 8 executes the following
only after its signed commit is accepted：

```bash
git worktree add -b coast/stage3-infrastructure-lane \
  /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure \
  "$(git rev-parse HEAD)"
cmake -S /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure \
  -B /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure/build/debug \
  -DCMAKE_BUILD_TYPE=Debug -DHUNDUN_BUILD_TESTS=ON
```

创建前只读确认目标路径不存在、branch 不存在、main worktree clean；不删除或复用未知
目录。worker 永远不进入 main worktree。

---

## Phase F1 — Scientific composition spine

### Task S3-C1: Constant-density IBM+WALE authority

**Ownership:** main agent only；不得把整个 task 或测试子集按 worker packet 派发。

**Consumes:** accepted Task 11 `FixedStepImmersedFlow`、Task 12 `WaleModel`、Task
13 body-fitted frozen `mu_eff`。

**Produces:** legal `LFP-GCIBM/wale/constant` product path；IBM-aware WALE gradient；
same-attempt wall/operator `mu_eff`；driver fast case。

**References:** WALE paper；OpenFOAM WALE responsibility；Basilisk
`dirichlet_gradient/embed_force`；AMReX-Hydro projection staging。

**Files:**

- Modify: `include/hundun/flow_immersed.hpp`
- Modify: `src/flow_immersed.cpp`
- Create: `src/flow_immersed_wale_detail.hpp`
- Modify: `src/flow_immersed_access_detail.hpp`
- Modify: `src/app_immersed_flow_driver.cpp`
- Modify: `tests/support/flow_immersed_test_access.hpp`
- Create: `tests/mpi/test_immersed_wale_constant.cpp`
- Create: `tests/acceptance/stage3_flow_models_fast.sh`
- Modify: `tests/cmake/stage3_science_registration.cmake`

**Forbidden files:** material/ideal-gas implementation、Checkpoint、diagnostics provider
files、Task 11 MMS selectors。

**Interfaces:** public `ImmersedFlowStepAttemptReport::wale` 已存在。新增 test-only
accessors：

```cpp
static std::uint64_t wale_evaluation_count(const FixedStepImmersedFlow&);
static les::WaleCoefficientIdentity
wale_coefficient_identity(const FixedStepImmersedFlow&);
static std::uint64_t wall_effective_viscosity_fingerprint(
    const FixedStepImmersedFlow&);
```

At C1 the executable accepts only no argument for the 8-cubed direct case and `fast` for the
12-cubed task case；only direct 1-rank and fast 1/2-rank CTest rows are registered。Unknown selector
or extra argv returns 2。S1 later adds and exclusively registers all formal 24/48 modes after the
combined product path is stable。

- [ ] **Step 0: add compile-preserving test-access seams**

Declare the three test-only accessors in `flow_immersed_test_access.hpp` and route them through
`flow_immersed_access_detail.hpp` under the existing test-access macro。Before C1 implementation
they return an explicit unavailable identity/count (or throw the stable test-only capability error)
without changing production state；all symbols link。No public non-test API is added。This ensures
the next RED fails because combined IBM+WALE is unsupported，not because an accessor is missing。

- [ ] **Step 1: 写 one-evaluation 和 rollback RED**

在 `test_immersed_wale_constant.cpp` 构造 8x8x8 periodic/background + closed sphere
case，传入 `domain` 和 `WaleModel`。测试代码必须包含：

```cpp
auto failed = flow.attempt(state, physics, stencil, {}, {});
HUNDUN_CHECK(!failed.wale.has_value());
HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
    before_committed, state.snapshot(flow::FlowLayer::committed)));

auto accepted = flow.attempt(state, physics, stencil, {}, {});
HUNDUN_CHECK(accepted.wale.has_value());
HUNDUN_CHECK(std::get<flow::StepAttemptReport>(accepted.base)
                 .pressure_corrector_count == 2U);
HUNDUN_CHECK(flow::test::ImmersedFlowTestAccess::wale_evaluation_count(flow) ==
             1U);
HUNDUN_CHECK(accepted.force.has_value());
```

同时加入 mutations：final velocity 重算 WALE、force 使用 molecular `mu`、第二
corrector 重新 evaluate、interface gradient 读取 solid slot、替换 `report.wale` 后旧
`ImmersedFlowDiagnosticSource` 仍被接受。

- [ ] **Step 2: 运行 RED**

```bash
cmake --build build/debug -j32 --target test_immersed_wale_constant
ctest --test-dir build/debug -R '^test_immersed_wale_constant_1_rank$' \
  --output-on-failure
```

Expected: FAIL，当前 facade collective 拒绝 `domain != nullptr && wale != nullptr`。

- [ ] **Step 3: 实现 `ImmersedWaleAttemptAuthority`**

在 detail header 中实现拥有 workspace 的 move-only class。BE 使用 `u^n`，BDF2 使用
`u^n + dt/dt_previous*(u^n-u^(n-1))`。regular cell 用 accepted background gradient，
interface cell 用现有 immersed reconstruction；不访问 inactive/solid physical value。

- [ ] **Step 4: 把同一 authority 接入所有消费者**

`mu_eff` 必须进入 momentum predictor、两次 corrector 使用的 diagonal/traction、final
momentum residual、wall force。pressure/force authority 仍由 Task 11 对象提供，不新增
第二份 wall pressure reconstruction。

同一 task 扩展 `diagnostic_report_seal`：constant base、四字段 force、WALE identity/
min/max/l2/zero/active count 全部进入 seal；任何字段 mutation 使 source authentication
失败。这里只修改 report authentication，不新增 provider。

- [ ] **Step 5: 开放 driver 组合**

删除仅针对 combined IBM+WALE 的拒绝，保持其他未实现密度组合显式拒绝。fast script
检查 stdout 中 `correctors=2`、四字段 force 和 WALE identity/min/max。

- [ ] **Step 6: 运行 task gate**

```bash
cmake --build build/debug -j32 --target \
  test_immersed_wale_constant hundun test_stage3_flow_header_contract
ctest --test-dir build/debug -j24 --output-on-failure -R \
  '^test_stage3_flow_header_contract$'
ctest --test-dir build/debug --output-on-failure -R \
  '^(test_immersed_wale_constant_[12]_rank|test_task19a_immersed_flow_dispatch_[12]_rank|test_stage3_flow_models_fast_[12]_rank)$'
```

MPI/selector tests 注册 `PROCESSORS`、`RESOURCE_LOCK hundun_stage3_mpi_m` 和 120/300 秒
timeout；命令不使用 `-j4` 并发多个 M job。只增加一个 12-cubed fast selector；不运行
24/48 wake。

- [ ] **Step 7: 主 agent 审查并提交**

检查 `WaleModel::evaluate` product caller 数量、wall viscosity fingerprint、Task 11
authority、rollback、allocation 和 1/2-rank ordering。主 agent执行：

```bash
git diff --check
git commit -s -m "feat: couple WALE to immersed flow"
```

- [ ] **Step 8: create the isolated infrastructure lane**

Run the exact P0 Step 4 worktree/configure commands from the accepted C1 HEAD。Record branch、
worktree gitdir pointer、CMake cache source root and clean status before dispatching R1。

---

### Task S3-D1: Material-density IBM vertical slice

**Ownership:** main agent only；该 vertical slice 不委派。

**Consumes:** S3-C1、Stage 2 material transport、final `FaceMassFlux`、accepted Ghost
reconstruction。

**Produces:** `LFP-GCIBM/none/material` driver；nonzero normal density reconstruction；
active mass conservation；density adapter seam。

**References:** incflo `incflo_update_density.cpp` 与 advance stage；AMReX flux accounting；
Basilisk homogeneous wall flux semantics。

**Files:**

- Modify: `include/hundun/flow_immersed.hpp`
- Modify: `src/flow_immersed.cpp`
- Create: `src/flow_immersed_density_detail.hpp`
- Modify: `src/flow_density_closure_detail.hpp`
- Modify: `src/flow_material_density_transport.cpp`
- Modify only for private access: `include/hundun/flow_material_density_transport.hpp`
- Modify only for authenticated report bridge: `include/hundun/flow_material_density_piso.hpp`
- Modify only if report sealing needs an out-of-line bridge: `src/flow_material_density_piso.cpp`
- Modify: `src/app_immersed_flow_driver.cpp`
- Create: `tests/mpi/test_immersed_material_density.cpp`
- Create: `tests/mpi/test_immersed_material_transaction.cpp`
- Modify: `tests/acceptance/stage3_flow_models_fast.sh`
- Modify: `tests/cmake/stage3_science_registration.cmake`

**Interfaces:** add exactly the `ImmersedFlowDensitySetup` from the v2 design and an overload:

```cpp
static FixedStepImmersedFlow create(
    const runtime::StructuredDecomposition&,
    const mesh::MeshTopology&, const mesh::MeshGeometry&,
    const boundary::BoundaryRegistry&,
    const immersed::ImmersedDomain*,
    const immersed::GhostStencilPlan*,
    const immersed::WallQuadraturePlan*,
    const immersed::LocalFlowPatternTransform*,
    const les::WaleModel*, ImmersedFlowDensitySetup,
    const runtime::MpiContext&, execution::ExecutionContext&,
    runtime::HaloExchange&, const linear::LinearSolver&,
    std::array<linear::Preconditioner*, 3>,
    const linear::LinearSolver&, linear::Preconditioner&);
```

旧 overload 构造 canonical constant setup，确保现有 callers source-compatible。

扩展已经是 `MaterialDensityStepAttemptReport` friend 的
`detail::DensityClosureBridge`。它从 adapter 的 existing `StepAttemptReport`、material
transport report、attempt identity 和 failure reason 构造并 seal 既有类型；不得新增
`ImmersedDensityReportBridge` 或 immersed-only material report。

```cpp
static MaterialDensityStepAttemptReport make_material_report(
    StepAttemptReport flow,
    std::optional<MaterialDensityTransportReport> material,
    MaterialTransportFailureReason failure,
    std::uint64_t material_field_count,
    runtime::FieldId shared_face_mass_flux_field,
    std::uint64_t attempt_identity);
```

该函数从 nested material report 派生 finalization identity、final-flux provenance 和
residual/conservation availability，验证 nested attempt identity 后调用原 `seal()`。

- [ ] **Step 0: 添加 compile-preserving disabled seam**

先声明 setup/new overload 和 bridge factory。旧 overload 仍走 accepted constant path；
new overload 对 material setup collective 返回稳定 `invalid_input/unsupported` report，
不进入 time step、不改 state。header contract 必须先编译通过，使下一步 RED 是运行时
unsupported-path failure，而不是不存在符号导致的编译失败。

- [ ] **Step 1: 写 density-gradient RED**

使用：

```text
rho(x,y,z) = 1 + 0.10*x + 0.05*y - 0.03*z
```

在 sphere interface 比较 analytic normal derivative。homogeneous-Neumann mutation、
clipping、epsilon、stencil fallback 都必须失败。

- [ ] **Step 2: 写 conservation/rollback RED**

```cpp
HUNDUN_CHECK(std::holds_alternative<flow::MaterialDensityStepAttemptReport>(
    report.base));
HUNDUN_CHECK(wall_mass_flux_bits == std::bit_cast<std::uint64_t>(0.0));
HUNDUN_CHECK(relative_active_mass_error <= 5.0e-12);
HUNDUN_CHECK(all_owned_active_density_positive_finite);
HUNDUN_CHECK(inactive_slots_are_positive_zero);
```

注入 negative `rho_wall`、NaN donor、rank-local failure；committed/history/final flux/
metadata 必须 bitwise neutral。另逐字段 mutation authenticated material base，并断言
`diagnostic_report_seal` 拒绝替换后的 report。

- [ ] **Step 3: 运行 RED**

```bash
cmake --build build/debug -j32 --target \
  test_immersed_material_density test_immersed_material_transaction
ctest --test-dir build/debug -R \
  '^(test_immersed_material_(density|transaction)_1_rank)$' \
  --output-on-failure
```

Expected: FAIL，当前 facade 只接受 constant density。

- [ ] **Step 4: 实现 private density adapter**

用 `std::variant` 实现 constant/material 两种 coarse-grained adapter。material adapter
先做 predictor transport，第一次 corrector 后做 provisional transport，第二次 corrector
后只用最终 corrected face flux完成 final transport；每个 corrector 前由同一 Ghost row
authority 生成 `rho_wall` 和 `D_wall`，验证正且有限。FlowState、transport、pressure
authority 全部 prepare 成功后才跨 collective publish；所有 publish 必须 `noexcept`、
allocation-free、MPI-free。adapter rollback 必须在 FlowState rollback 前后均 noexcept/
幂等。

- [ ] **Step 5: 接入 driver**

driver 根据 schema v3 构造 fields/spec/setup。material IBM fast case 8x8x8、一步；
stdout 不声称 Restart 或 WALE。

- [ ] **Step 6: task gate**

```bash
cmake --build build/debug -j32 --target \
  test_immersed_material_density test_immersed_material_transaction \
  test_material_density_transport test_material_density_piso \
  test_material_density_transport_header_contract \
  test_material_density_piso_header_contract \
  test_stage3_flow_header_contract hundun
ctest --test-dir build/debug -j24 --output-on-failure -R \
  '^(test_material_density_(piso|transport)_header_contract|test_stage3_flow_header_contract)$'
ctest --test-dir build/debug --output-on-failure -R \
  '^(test_immersed_material_(density|transaction)_[12]_rank|test_material_density_(transport|piso)_1_rank|test_stage3_flow_models_fast_[12]_rank)$'
```

不运行 1/2/4 full matrix、24/48 或 complete sanitizer。

- [ ] **Step 7: 主 agent impact review 和提交**

用 codegraphf 检查 `MaterialDensityTransport::advance`、
`FixedStepMaterialDensityFlow::attempt`、final flux callers。确认 Stage 2 body-fitted path
未改行为，提交：

```bash
git commit -s -m "feat: couple material density to immersed flow"
```

---

### Task S3-C2: Material-density body-fitted/IBM+WALE join

**Ownership:** main agent。

**Consumes:** S3-C1 WALE authority、S3-D1 material adapter。

**Produces:** material `rho_attempt` → one WALE evaluation → same `mu_eff`/final flux，覆盖
`none/wale/material` 和 `LFP-GCIBM/wale/material` 两种 schema v3 path。

**References:** incflo variable-density stage order；WALE paper；AMReX-Hydro projection。

**Files:**

- Modify: `src/flow_immersed.cpp`
- Modify: `src/flow_immersed_density_detail.hpp`
- Modify: `src/flow_immersed_wale_detail.hpp`
- Modify: `src/flow_density_closure_detail.hpp`
- Create: `src/flow_body_fitted_wale_detail.hpp`
- Modify: `src/flow_constant_density_piso.cpp`
- Modify: `include/hundun/flow_material_density_piso.hpp`
- Modify: `src/flow_material_density_piso.cpp`
- Modify: `src/app_immersed_flow_driver.cpp`
- Modify: `tests/mpi/test_wale_body_fitted.cpp`
- Create: `tests/mpi/test_material_wale_composition.cpp`
- Modify: `tests/acceptance/stage3_flow_models_fast.sh`
- Modify: `tests/cmake/stage3_science_registration.cmake`

**Public interface:** none。`ImmersedFlowStepAttemptReport` 既有 material base variant +
optional WALE + force 已足够。私有 bridge 签名冻结为：

```cpp
static MaterialDensityStepAttemptReport attempt_with_optional_wale(
    const FixedStepMaterialDensityFlow&, FlowState&, double molecular_mu,
    const MomentumTimeStencil&, const linear::SolveControl&,
    const linear::SolveControl&, const DensityClosureHooks*,
    const les::WaleModel*, les::WaleSummary*);
```

`model` 和 `summary` 必须同时为空或同时非空。existing public material `attempt` 和 ideal
bridge 都委托该函数的 null-WALE 分支；不得复制 body-fitted material PISO。

`FixedStepMaterialDensityFlow::attempt_common` 的 private declaration 同步追加最后两个
参数 `const les::WaleModel*` 和 `les::WaleSummary*`。现有 public `attempt` 传
`closure=null/model=null/summary=null`；现有 ideal bridge 传 closure 且 WALE 两项为空；
上面的新 bridge 是唯一可同时传 closure/WALE 的 caller。不得新增另一条 public
variable-density attempt。

`flow_body_fitted_wale_detail.hpp` 从 accepted constant implementation 提取
`BodyFittedWaleAttemptData` 和一个 internal template helper。helper 必须显式接收
`runtime::FieldView<const double> rho_attempt`，而不是自行选择 committed/trial density；
它继续从 committed/history velocity 构造 BE/BDF2 lagged gradient，并返回同一组 cell
coefficient、face `mu_sgs/mu_eff` 和 summary。constant caller 传 committed density；
material/ideal caller 只能在 predictor transport 和 optional predictor closure 完成、
trial density halo 已同步后传 trial density。helper 不拥有 transport、PISO 或 commit。

Selector contract at C2: accepted `test_wale_body_fitted` keeps its existing 12/24 short product
screens，and `test_material_wale_composition` accepts only `fast`。C2 registers no new formal row；
S1 later adds `formal` argv/CTest registrations and owns their resource metadata。

- [ ] **Step 1: 写 stale-density/second-evaluation RED**

分别用 `domain=nullptr` 和 closed sphere 构造 committed density 与 predicted
`rho_attempt` 不同的 case。两条路径都断言 WALE `mu_sgs` 使用 predicted density；
强制使用 committed density、corrector 内第二 evaluate、final velocity refresh 的
mutations 必须失败。

先把 Task 13 的 regular-cell lagged gradient/face interpolation 移到
`flow_body_fitted_wale_detail.hpp`，但保持调用点和计算顺序不变。已有
`test_wale_body_fitted` 对 refactor 前保存的 field/report fingerprints 做 bitwise RED；
material seam 仍返回 unsupported，因此本步骤可以编译并实际运行 RED。

- [ ] **Step 2: 写 retry identity RED**

第一次 attempt 在 momentum 后 recoverable failure；第二次使用减半 dt。断言 state
rollback bitwise、第二次 identity 按新 dt/rho_attempt 重算、只发布第二次 force/WALE。

- [ ] **Step 3: 运行 RED**

```bash
cmake --build build/debug -j32 --target \
  test_wale_body_fitted test_material_wale_composition
ctest --test-dir build/debug -R '^test_wale_body_fitted_1_rank$' \
  --output-on-failure
ctest --test-dir build/debug -R '^test_material_wale_composition_1_rank$' \
  --output-on-failure
```

Expected: the accepted body-fitted characterization row PASSes before and after the refactor；the
material composition row FAILs because D1 路径尚未允许 material+WALE。A compile/link failure is
not an acceptable RED。

- [ ] **Step 4: 最小 join**

body-fitted 固定调用顺序：既有 material predictor transport → `rho_attempt` view → shared
regular gradient/WALE → 既有 material momentum/PISO #1 → provisional transport → PISO #2
→ final transport/assessment/commit。IBM 固定调用顺序相同，但 gradient/viscosity 使用 C1
authority。两者不得复制 material transport、pressure solve、final-flux finalizer 或
WALE tensor kernel。

- [ ] **Step 5: fast gate 和提交**

```bash
cmake --build build/debug -j32 --target \
  test_wale_body_fitted test_material_wale_composition test_material_density_piso \
  test_material_density_piso_header_contract test_adaptive_time_control hundun
ctest --test-dir build/debug --output-on-failure -R \
  '^(test_material_density_piso_header_contract|test_adaptive_time_control)$'
ctest --test-dir build/debug --output-on-failure -R \
  '^(test_wale_body_fitted_[12]_rank|test_material_wale_composition_[12]_rank|test_material_density_piso_1_rank|test_stage3_flow_models_fast_[12]_rank)$'
git diff --check
git commit -s -m "feat: combine material density IBM and WALE"
```

主 agent必须用 codegraphf 证明 profile 5 调用既有 material PISO core，profile 6 调用
immersed PISO；同一 profile 内只有一个 pressure/final-flux authority。

---

### Task S3-D2: Ideal-gas IBM vertical slice

**Ownership:** main agent。

**Consumes:** S3-D1 density adapter、Stage 2 `IdealGasClosure`、active volumes。

**Produces:** `LFP-GCIBM/none/ideal_gas` driver；open fixed-p0 与 closed dynamic-p0；
closure rollback。

**References:** Bell--Marcus variable-density projection；incflo stage separation；AMReX
active geometry facts。外部 low-Mach codes 只作方程背景，不复用 implementation。

**Files:**

- Modify: `include/hundun/flow_immersed.hpp`
- Modify: `src/flow_immersed.cpp`
- Modify: `src/flow_immersed_density_detail.hpp`
- Modify: `src/flow_density_closure_detail.hpp`
- Modify only for private bridge: `include/hundun/flow_ideal_gas_closure.hpp`
- Modify only for bridge: `src/flow_ideal_gas_closure.cpp`
- Modify only for authenticated report bridge: `include/hundun/flow_ideal_gas_piso.hpp`
- Modify only if report sealing needs an out-of-line bridge: `src/flow_ideal_gas_piso.cpp`
- Modify: `src/app_immersed_flow_driver.cpp`
- Create: `tests/mpi/test_immersed_ideal_gas.cpp`
- Create: `tests/mpi/test_immersed_ideal_gas_transaction.cpp`
- Modify: `tests/acceptance/stage3_flow_models_fast.sh`
- Modify: `tests/cmake/stage3_science_registration.cmake`

**Private interface:** use existing `detail::DensityClosureAdapter::bind(...)` for the
predictor/provisional/final hooks. Extend that already-friend adapter with:

```cpp
static IdealGasStepAttemptReport make_ideal_report(
    MaterialDensityStepAttemptReport,
    std::optional<IdealGasClosureReport>,
    std::uint64_t attempt_identity);
```

Do not add a second closure bridge or a second ideal-gas report type.

- [ ] **Step 0: add compile-preserving ideal branch seam**

Add setup validation and the private factory signature, but let ideal setup return a stable
collective unsupported report before opening a trial. Header tests must compile; the following
RED must execute rather than fail because a symbol is missing.

- [ ] **Step 1: 写 active-volume closed-domain RED**

```text
T = 300 + 20*E(x,y,z)*F_b(x,y,z)^2
h = cp*T
p0 = M_target*R / sum_active(V_i/T_i)
```

断言 inactive volumes 不进入 sum，mass error `<=5e-12`，`h=cp*T` 和
`rho=p0/(R*T)` relative error `<=1e-12`。

- [ ] **Step 2: 写 open-domain 和 rollback RED**

open case 固定 positive p0；pressure outlet 只约束 `pi`。注入 non-positive T/rho、
closure mismatch、p0 update 后 rank-local failure；state/history/controller/final flux/p0
bitwise rollback。逐字段 mutation ideal/material nested report 和 closure report，断言
outer seal 与 `diagnostic_report_seal` 同时拒绝。

- [ ] **Step 3: 运行 RED**

```bash
cmake --build build/debug -j32 --target \
  test_immersed_ideal_gas test_immersed_ideal_gas_transaction
ctest --test-dir build/debug -R \
  '^(test_immersed_ideal_gas(_transaction)?_1_rank)$' \
  --output-on-failure
```

Expected: FAIL，density setup 的 ideal-gas branch 尚未实现。

- [ ] **Step 4: 扩展 adapter**

ideal adapter 复用 material conservative transport，并在 predictor transport 后、第一次
corrector 后的 provisional transport 后、第二次 corrector 后的 final transport 后调用
既有 closure。`p0`/transport/pressure authority/FlowState 全部先 prepare；一次 collective
pass 后才执行 no-throw、allocation-free、MPI-free publish。不增加 outer enthalpy
iteration 或第二 density field。`DensityClosureAdapter::make_ideal_report` 验证 material/
closure identity 和 final closure stage 后调用原 `seal()`。

- [ ] **Step 5: driver 和 task gate**

```bash
cmake --build build/debug -j32 --target \
  test_immersed_ideal_gas test_immersed_ideal_gas_transaction \
  test_ideal_gas_closure test_ideal_gas_piso \
  test_ideal_gas_header_contract test_stage3_flow_header_contract hundun
ctest --test-dir build/debug -j24 --output-on-failure -R \
  '^(test_ideal_gas_header_contract|test_stage3_flow_header_contract)$'
ctest --test-dir build/debug --output-on-failure -R \
  '^(test_immersed_ideal_gas(_transaction)?_[12]_rank|test_ideal_gas_(closure|piso)_1_rank|test_stage3_flow_models_fast_[12]_rank)$'
```

- [ ] **Step 6: 主审查和提交**

检查所有 p0/target mass/revision、open/closed authority、Checkpoint-facing state 和
Stage 2 ideal-gas callers。提交：

```bash
git commit -s -m "feat: couple ideal gas to immersed flow"
```

---

### Task S3-C3: Ideal-gas body-fitted/IBM+WALE and combined Gate 5

**Ownership:** main agent；这是不可委派科学 verdict。

**Consumes:** C1、D1、C2、D2。

**Produces:** all three density models 的 legal body-fitted WALE 与 IBM+WALE；关闭旧
Task 16/Gate 5。

**References:** WALE paper；incflo variable-density stage order；Basilisk force/flux
consistency；HUNDUN Task 11 authority。

**Files:**

- Modify: `src/flow_immersed.cpp`
- Modify: `src/flow_immersed_density_detail.hpp`
- Modify: `src/flow_immersed_wale_detail.hpp`
- Modify: `src/flow_density_closure_detail.hpp`
- Modify: `src/app_immersed_flow_driver.cpp`
- Create: `tests/mpi/test_ideal_gas_wale_composition.cpp`
- Create: `tests/mpi/test_immersed_combined_retry.cpp`
- Modify: `tests/acceptance/stage3_flow_models_fast.sh`
- Modify: `tests/cmake/stage3_science_registration.cmake`

Selector contract: `test_ideal_gas_wale_composition` and
`test_immersed_combined_retry` accept only `fast` in C3。S1 later adds and exclusively registers
formal selectors only for `test_ideal_gas_wale_composition`；combined-retry remains a V0 low-cost
row whose exact-C evidence is reused by V1，so S1 does not modify that source。

- [ ] **Step 1: 写 combined order RED**

对六个 WALE-enabled legal paths（profiles 2/3/5/6/8/9）记录 stage trace，必须精确
等于下列 conditional trace；`*-if-immersed` 只在 profiles 3/6/9 出现，
`*-if-variable` 只在 profiles 5/6/8/9 出现，`*-if-ideal` 只在 profiles 8/9
出现。条件不满足时 event 必须 absent，不允许用 no-op/fake module 记录占位：

```text
begin
density-predict-if-variable
closure-predict-if-ideal
lagged-gradient
wale-evaluate-once
momentum-predict
pressure-wall-authority-1-if-immersed
pressure-corrector-1
provisional-density-transport-if-variable
closure-provisional-if-ideal
pressure-wall-authority-2-if-immersed
pressure-corrector-2
final-density-transport-if-variable
closure-final-if-ideal
post-closure-assessment-if-ideal
final-residual
force-if-immersed
prepare-flow-state
prepare-density-transport-if-variable
prepare-closure-if-ideal
prepare-pressure-authority-if-immersed
collective-ready
publish-noexcept
```

任何 third corrector、final WALE refresh、跳过 provisional closure、把 provisional flux
当作 final flux 或在 collective-ready 后执行可能失败的动作都使 RED 失败。

- [ ] **Step 2: 写三密度 retry/collective RED**

每条路径注入一个 rank-local recoverable failure；断言 lowest failing rank、bitwise
rollback、retry identity、force/report isolation 和最终 commit。

- [ ] **Step 3: 运行 combined RED**

```bash
cmake --build build/debug -j32 --target \
  test_ideal_gas_wale_composition test_immersed_combined_retry
ctest --test-dir build/debug -R \
  '^(test_(ideal_gas_wale_composition|immersed_combined_retry)_1_rank)$' \
  --output-on-failure
```

Expected: FAIL，当前 ideal+WALE path 仍返回 stable unsupported，或 stage trace 缺少上述
provisional/prepare/publish events；target 必须编译并实际启动。

- [ ] **Step 4: 实现 ideal+WALE join 并删除合法组合拒绝**

不再修改 public API。`domain=nullptr` 调用 C2 的 existing material PISO private WALE
hook，并把 D2 `DensityClosureHooks` 传入同一 material core；`domain!=nullptr` 调用 D2
immersed adapter。两条路径共享 WALE identity contract，但不得共享/复制 pressure
solver。driver 只接受完整 setup，缺少 cp/R/p0 或 WALE control 时在 time step 前
collective 拒绝。

- [ ] **Step 5: combined fast gate**

```bash
cmake --build build/debug -j32 --target \
  test_ideal_gas_wale_composition test_immersed_combined_retry hundun
ctest --test-dir build/debug --output-on-failure -R \
  '^(test_(ideal_gas_wale_composition|immersed_combined_retry)_[12]_rank|test_stage3_flow_models_fast_[12]_rank)$'
```

fast case 最大 12-cubed、2 steps；不运行 finite-wake long statistics。

- [ ] **Step 6: 主 agent Gate 5 review 和提交**

从 C1 parent 到 C3 candidate 审查一次完整 combined diff。搜索多次 WALE evaluate、
stale density、non-final flux、wall molecular viscosity、persistent attempt report。
同时证明 profile 8 进入既有 material PISO + closure hooks，profile 9 进入 immersed
PISO；两者没有第二个 body-fitted variable-density pressure/final-flux authority。

```bash
git diff --check
git commit -s -m "feat: complete immersed density and WALE composition"
```

---

### Task S3-S1: Freeze final scientific selectors without running long rows

**Ownership:** main agent only。The TGV self-convergence oracle and final matrix semantics are scientific
judgments and are not delegated。

**Dependency:** C3 accepted。

**Consumes:** accepted WALE tensor/y-cubed tests；body-fitted/immersed product paths from C1--C3；
Task 11 authority and current scientific thresholds。

**Produces:** mutation-sensitive selector/output schema for WALE TGV convergence、body-fitted
channel、constant IBM+WALE and variable-density correctness；only a 12-cubed smoke runs now。

**References:** Nicoud--Ducros WALE equations；smooth Taylor--Green self-convergence；HUNDUN
Task 11 force/flux authority。No external result is used as pass evidence。

**Files:**

- Create: `tests/support/stage3_scientific_row.hpp`
- Create: `tests/support/stage3_scientific_row.cpp`
- Create: `tests/unit/test_stage3_scientific_row.cpp`
- Create: `tests/numerical/test_wale_taylor_green.cpp`
- Modify: `tests/mpi/test_wale_body_fitted.cpp`
- Modify: `tests/mpi/test_immersed_wale_constant.cpp`
- Modify: `tests/mpi/test_material_wale_composition.cpp`
- Modify: `tests/mpi/test_ideal_gas_wale_composition.cpp`
- Modify: `tests/cmake/stage3_science_registration.cmake`

**Frozen row schema:** each executable emits one canonical line containing
`row_id/cells/ranks/process_grid/steps/dt/final_time/velocity_l2/pressure_l2/nu_t_l2/continuity/
conservation/closure/force/WALE_identity/status`。Unavailable fields use an explicit availability
bit, never a numeric zero placeholder。The launcher rejects duplicate/missing rows and aggregates
only rows with identical physical problem/final time。

S1 is the exclusive owner of every formal CTest registration in
`stage3_science_registration.cmake`；C1--C3 retain only direct/fast rows。All 24-cubed formal rows
set `PROCESSORS` to ranks、`RESOURCE_LOCK hundun_stage3_mpi_m` and timeout 7200 seconds；all
48-cubed rows set `RESOURCE_LOCK hundun_stage3_mpi_h` and timeout 43200 seconds；12-cubed formal
variable-density rows use the M lock and timeout 1800 seconds。G1 validates these exact values but
does not define or rewrite them。

**TGV formal contract:** use smooth periodic Taylor--Green initial data、constant physical
viscosity and WALE controls。12/24/48 single-rank runs reach the same fixed physical time with dt
scaled by h；the aggregator restricts fine fields to the coarse cell-average layout and computes
two Richardson segments for velocity、gauge-normalized pressure and `nu_t` L2。Each error must be
finite、strictly positive、strictly decreasing and each observed segment must satisfy the existing
second-order bound `>=1.8`。A separate 24-cubed 1/2/4 decomposition row compares the same final
fields/fingerprints。This is self-convergence through the product path, not an alternate solver。

**Other formal rows:** body-fitted channel runs 48-cubed single-rank with the previously frozen
WALE baseline checks；constant IBM+WALE runs 48 single-rank and 24 on 1/2/4；material and ideal
IBM+WALE run 12/24 short 1/2/4 correctness with conservation/closure/retry。Warped、prism and
96-cubed selectors are never registered。

- [ ] **Step 1: write selector/schema/mutation RED**

Mutations cover mismatched final time、point-sample instead of cell-average restriction、average
slope hiding one failed segment、zero/epsilon error、wrong rank/process-grid、stale density、second
WALE evaluation and missing availability bit。All formal executables accept only their frozen argv
and return 2 for unknown selectors。

- [ ] **Step 2: register formal rows but execute only the smoke**

```bash
cmake --build build/debug -j32 --target \
  test_stage3_scientific_row_contract test_wale_taylor_green test_wale_body_fitted \
  test_immersed_wale_constant test_material_wale_composition \
  test_ideal_gas_wale_composition
ctest --test-dir build/debug -j24 --output-on-failure -R \
  '^test_stage3_scientific_row_contract$'
ctest --test-dir build/debug --output-on-failure -R \
  '^test_wale_taylor_green_12_smoke_1_rank$'
```

Expected initial RED: the schema/mutation contract fails on a deliberately removed pressure row or
the smoke reports unsupported；all targets compile and execute。Do not run 24/48。

- [ ] **Step 3: implement test-only row/aggregation support**

Use the same product executable/test facades and accepted fields；add the frozen `formal` argv and
CTest rows only here。Test support may restrict and
compare output but may not solve a PDE、recompute a replacement WALE field、filter data or change
thresholds。Formal rows use the M/H locks and timeouts frozen above；G1 only validates them。

- [ ] **Step 4: 12-cubed GREEN and review**

Run only the contract and 12-cubed smoke above。Main agent reviews equations、cell-average
restriction、normalization、selectors and confirms the `formal` label is not consumed by any
development task。Commit `test: freeze Stage 3 scientific selectors`。

---

## Phase F2 — Infrastructure lane

### Task S3-R1: Additive constant WALE Checkpoint v3 profiles

**Ownership:** bounded default worker eligible；主 agent负责 protocol/DCO verdict。

**Scheduling:** C1 accepted 后开始；可与 D1/C2 主 lane 并行。不得修改 driver。

**Consumes:** accepted 17A byte fixture、C1 constant IBM+WALE、Task 13 WALE-only。

**Produces:** profile enum values 2/3、const-correct write/read module views、codec/rollback
tests。

**References:** HUNDUN 17A protocol only；upstream checkpoint formats不是 byte authority。

**Files:**

- Modify: `include/hundun/flow_checkpoint_v3.hpp`
- Modify: `src/flow_checkpoint_v3.cpp`
- Modify: `src/flow_checkpoint_v3_detail.hpp`
- Modify: `tests/unit/test_checkpoint_v3_header_contract.cpp`
- Modify: `tests/unit/test_checkpoint_v3_codec.cpp`
- Create: `tests/mpi/test_checkpoint_v3_wale.cpp`
- Modify: `tests/cmake/stage3_checkpoint_registration.cmake`

**Worker Packet `S3-R1-codec`:**

```text
repository: /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure
build_tree: /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure/build/debug
allowed_files: exactly the seven Files entries above
forbidden_files: driver, flow_immersed*, density/closure, diagnostics, any main worktree file
resource_group: L for codec/header; one locked M job for each WALE continuation test
public_reference: none; accepted HUNDUN 17A byte fixture/protocol is the only format authority
integration_owner: main agent after full protocol review
```

**Interfaces:** preserve old profile-1 overload and add exactly:

```cpp
struct CheckpointV3WriteModules final {
  CheckpointV3Presence presence{CheckpointV3Presence::constant_static_ibm};
  const immersed::ImmersedSurface* surface{};
  const immersed::SurfaceQuery* query{};
  const immersed::ImmersedDomain* domain{};
  const immersed::GhostStencilPlan* ghost_plan{};
  const immersed::WallQuadraturePlan* wall_plan{};
  const immersed::LocalFlowPatternTransform* transform{};
  const les::WaleModel* wale{};
  const FixedStepImmersedFlow* flow{};
  const IdealGasClosure* ideal_gas{};
};

struct CheckpointV3ReadModules final {
  CheckpointV3Presence presence{CheckpointV3Presence::constant_static_ibm};
  const immersed::ImmersedSurface* surface{};
  const immersed::SurfaceQuery* query{};
  const immersed::ImmersedDomain* domain{};
  const immersed::GhostStencilPlan* ghost_plan{};
  const immersed::WallQuadraturePlan* wall_plan{};
  const immersed::LocalFlowPatternTransform* transform{};
  const les::WaleModel* wale{};
  FixedStepImmersedFlow* flow{};
  IdealGasClosure* ideal_gas{};
};

CheckpointV3Report write_checkpoint_v3(
    const runtime::MpiContext&, const runtime::StructuredDecomposition&,
    const mesh::MeshTopology&, const mesh::MeshGeometry&,
    const boundary::BoundaryRegistry&, const config::ImmersedFlowCaseConfig&,
    const CheckpointV3WriteModules&, const FlowState&,
    CheckpointV3ControlState, const std::filesystem::path&);

CheckpointV3ReadResult read_checkpoint_v3(
    const runtime::MpiContext&, const runtime::StructuredDecomposition&,
    const mesh::MeshTopology&, const mesh::MeshGeometry&,
    const boundary::BoundaryRegistry&, const config::ImmersedFlowCaseConfig&,
    const CheckpointV3ReadModules&, FlowState&, const std::filesystem::path&);
```

Pointers are borrowed for one call. Profile 2 requires only WALE among optional static modules and
forbids IBM/ideal pointers；profile 3 requires all IBM pointers + WALE + flow and forbids ideal；
old profile 1 rules and overload remain byte-identical。

- [ ] **Step 0: add compile-preserving views and reject 2/3**

Declare enum values, views and overloads; new overloads validate pointer presence and then return
stable `presence` failure for profiles 2/3 before I/O. Old overload delegates through profile 1.
This lets codec/header tests compile and makes RED a behavioral failure.

- [ ] **Step 1: 扩展 header/codec RED**

断言 enum 1/2/3 exact values，unknown profile、missing/extra module、nonzero absent section、
WALE transient identity persistence 都失败。旧 constant IBM golden bytes 必须相同。

- [ ] **Step 2: 运行 RED**

```bash
cmake --build \
  /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure/build/debug \
  -j32 --target test_checkpoint_v3_codec test_checkpoint_v3_header_contract hundun
ctest --test-dir \
  /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure/build/debug \
  -R '^(test_checkpoint_v3_(codec|header_contract))$' --output-on-failure
```

Expected: FAIL，新 profile 尚不存在。

- [ ] **Step 3: 最小 additive implementation**

保留 version/endian/CRC/section order；只在 manifest 的既有 presence byte 接受 2/3，
并增加 profile-governed identity sections。不要保存 `nu_t/mu_sgs/mu_eff`。

- [ ] **Step 4: focused gate**

```bash
cmake --build \
  /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure/build/debug \
  -j32 --target \
  test_checkpoint_v3_codec test_checkpoint_v3_wale \
  test_checkpoint_v3_header_contract hundun
ctest --test-dir \
  /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure/build/debug \
  --output-on-failure -R \
  '^(test_checkpoint_v3_(codec|header_contract|wale_[12]_rank))$'
```

WALE MPI tests 注册 resource lock 和 120 秒 timeout。worker 返回 diff，不提交。主 agent
在 infra worktree 复核旧 fixture、CRC、publish-last、failed-read neutrality 后签署
handoff commit `feat: extend Checkpoint v3 for constant WALE`；集成 main 后重跑同一 gate。

---

### Task S3-O1: WALE diagnostics and additive module kinds

**Ownership:** bounded default worker eligible；在 R1 接受后由同一 worker lane 执行。

**Scheduling:** 可与 C2/D2 main lane 并行；不得修改 `flow_immersed.cpp`。

**Consumes:** `les::WaleSummary`、minimal diagnostics schema v1。

**Produces:** stable enum values 18--22、direct WALE provider、absence inventory tests。

**References:** AMReX support-level inventory/counter separation；OpenFOAM model summary
responsibility。

**Files:**

- Modify: `include/hundun/diag_structured.hpp`
- Create: `include/hundun/diag_les_wale.hpp`
- Create: `src/diag_les_wale.cpp`
- Modify: `src/CMakeLists.txt`
- Create: `tests/unit/test_wale_diagnostics.cpp`
- Create: `tests/unit/test_wale_diagnostics_header_contract.cpp`
- Modify: `tests/cmake/stage3_diagnostics_registration.cmake`

**Worker Packet `S3-O1-wale-diagnostics`:**

```text
repository: /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure
build_tree: /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure/build/debug
baseline: signed R1 handoff commit in that worktree
allowed_files: exactly the seven Files entries above
forbidden_files: flow_immersed*, driver, Checkpoint, density/closure, main worktree
resource_group: L only
public_references:
  AMReX EBFabFactory at 5273b558c13011573e1b6bf71860db4fcc1f0cfb:
    https://github.com/AMReX-Codes/amrex/blob/5273b558c13011573e1b6bf71860db4fcc1f0cfb/Src/EB/AMReX_EBFabFactory.H
    adopt only presence/support-level inventory separation
  OpenFOAM WALE.H at 627e2f909cb6a08dbb3a74e9a34aa632a975650e:
    https://github.com/OpenFOAM/OpenFOAM-dev/blob/627e2f909cb6a08dbb3a74e9a34aa632a975650e/src/MomentumTransportModels/momentumTransportModels/LES/WALE/WALE.H
    adopt only model-summary responsibility; GPL source/control flow remains forbidden
integration_owner: main agent after enum/API/read-only review
```

**Interface:**

```cpp
DiagnosticDescriptor describe_diagnostics(const les::WaleSummary&) noexcept;
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const les::WaleSummary&);
void collect_diagnostics(const les::WaleSummary&,
                         const DiagnosticRequest&, DiagnosticSink&);
```

- [ ] **Step 0: add compile-preserving provider seam**

Append enum 18--22 and declare the three overloads. Initial `collect_diagnostics` returns a stable
`capability` failure without mutating the summary; this is a temporary seam only so RED compiles and
runs. It must be removed by Step 4.

- [ ] **Step 1: enum/provider RED**

static_assert 0--17 unchanged and 18--22 exact。WALE summary fields include identity、
owned-active count、min/max/l2、exact-zero count；nonfinite uses status encoding。

- [ ] **Step 2: absence/read-only RED**

no `WaleSummary` means caller registers no provider；不得构造 fake zero summary。重复收集
canonical identical，不修改 summary 或 business counters。

- [ ] **Step 3: 运行 RED**

```bash
cmake --build \
  /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure/build/debug \
  -j32 --target test_wale_diagnostics test_wale_diagnostics_header_contract \
  test_structured_diagnostics \
  test_structured_diagnostics_header_contract hundun
ctest --test-dir \
  /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure/build/debug \
  -R '^(test_wale_diagnostics.*|test_structured_diagnostics|test_structured_diagnostics_header_contract)$' \
  --output-on-failure
```

Expected: FAIL with the stable capability assertion；target 必须编译并实际启动。

- [ ] **Step 4: 实现 additive enum/provider**

追加 18--22，不重排 0--17。provider 只读取 `WaleSummary`，用现有
`DiagnosticValueStatus` 表示 nonfinite，不注册 fake absent record。

- [ ] **Step 5: 运行 GREEN**

```bash
cmake --build \
  /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure/build/debug \
  -j32 --target \
  test_wale_diagnostics test_wale_diagnostics_header_contract \
  test_structured_diagnostics \
  test_structured_diagnostics_header_contract hundun
ctest --test-dir \
  /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure/build/debug \
  -R '^(test_wale_diagnostics.*|test_structured_diagnostics|test_structured_diagnostics_header_contract)$' \
  --output-on-failure
cmake -S /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure \
  -B /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure/build/tests-off \
  -DCMAKE_BUILD_TYPE=Release -DHUNDUN_BUILD_TESTS=OFF
cmake --build \
  /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure/build/tests-off \
  -j32 --target hundun
```

worker 不提交。主 agent在 infra worktree 复核 enum ABI、tests-off linkage 后签署
handoff commit `feat: add WALE diagnostics`；集成 main 后重跑同一 gate。

---

### Task S3-R2: All density/IBM/WALE continuation profiles

**Ownership:** main agent only。codec、closure restore 和 driver 是一个 transaction，
不拆给 worker。

**Dependency:** S1 accepted，R1/O1 signed handoff commits 已集成到 main。此任务开始后
不再与 A1 并行修改 driver。

**Consumes:** accepted Checkpoint v3 profile 1 bytes、R1 profiles 2/3、C1--C3 九种
runtime composition、D2 ideal-gas closure transaction 和 O1 enum values 18--22。

**Produces:** profiles 4--9；continuous-vs-restart bitwise for all legal combinations。

**References:** accepted HUNDUN Checkpoint v3 byte protocol；incflo checkpoint/restart 只作
组合测试思想，不作为格式或 publication authority。

**Files:**

- Modify: `include/hundun/flow_checkpoint_v3.hpp`
- Modify: `src/flow_checkpoint_v3.cpp`
- Modify: `src/flow_checkpoint_v3_detail.hpp`
- Modify: `include/hundun/flow_ideal_gas_closure.hpp`
- Modify: `src/flow_ideal_gas_closure.cpp`
- Modify: `src/flow_density_closure_detail.hpp`
- Modify: `src/app_immersed_flow_driver.cpp`
- Modify: `tests/unit/test_checkpoint_v3_codec.cpp`
- Modify: `tests/unit/test_checkpoint_v3_header_contract.cpp`
- Create: `tests/mpi/test_checkpoint_v3_density_profiles.cpp`
- Create: `tests/acceptance/stage3_restart_fast.sh`
- Modify: `tests/cmake/stage3_checkpoint_registration.cmake`

**Frozen profile matrix:**

| Profile | IBM objects | WALE | Material fields | Ideal closure | Diagnostics kinds |
| --- | --- | --- | --- | --- | --- |
| 1 | required | forbidden | forbidden | forbidden | 18--21 |
| 2 | forbidden | required | forbidden | forbidden | 22 |
| 3 | required | required | forbidden | forbidden | 18--22 |
| 4 | required | forbidden | required | forbidden | 18--21 |
| 5 | forbidden | required | required | forbidden | 22 |
| 6 | required | required | required | forbidden | 18--22 |
| 7 | required | forbidden | required | required | 18--21 |
| 8 | forbidden | required | required | required | 22 |
| 9 | required | required | required | required | 18--22 |

Every forbidden module is null and has zero sections/bytes. Every required module has exactly one
matching object/identity section. Ideal open domain stores fixed p0/no target mass；closed domain
stores dynamic p0 + positive target mass；both store non-wrapping closure revision.

**Selector/final-row contract:** `test_checkpoint_v3_density_profiles fast` runs the 8-cubed
development cases；`test_checkpoint_v3_density_profiles formal 12` runs all profiles 1--9 at a
fixed 12-cubed global mesh for four accepted steps，with a second run restarting after step 2 and
continuing to step 4。Unknown selector、cell count or extra argv returns 2。Formal output contains
one canonical row per profile with `profile/cells/ranks/process_grid/continuous_steps/restart_step/
state_bits/metadata_bits/pressure_authority_bits/optional_closure_bits/status`；every available
bit identity must match continuous vs restarted execution，and unavailable modules carry an
availability bit rather than zero evidence。

R2 registers exact CTest names `checkpoint-continuation-n12-r1`、`-r2`、`-r4` with
`PROCESSORS` equal to ranks、`RESOURCE_LOCK hundun_stage3_mpi_m`、timeout 1800 seconds and the
`formal;scientific;restart;stage3` labels。They are listable but never executed in R2/V0；G1 maps
these exact tests to its same-named inventory rows and V1 is the sole execution owner。

Closure restore private contract：

```cpp
struct IdealGasClosureCheckpointPreparedRestore;

struct IdealGasClosureCheckpointAccess final {
  static IdealGasClosureState snapshot(const IdealGasClosure&);
  static IdealGasClosureCheckpointPreparedRestore prepare_restore(
      const IdealGasClosure&, const FlowState& restored_state,
      IdealGasClosureState);
  static void publish_restore(
      IdealGasClosure&,
      IdealGasClosureCheckpointPreparedRestore&&) noexcept;
};
```

`prepare_restore` validates mode/p0/target/revision/EOS against the prepared FlowState but does not
mutate the live closure；`publish_restore` is allocation-free and MPI-free。

- [ ] **Step 0: compile-preserving enum/restore seam**

Append enum values 4--9 and extend header snapshots, while codec returns stable `presence` failure
for 4--9. Add private `IdealGasClosureCheckpointAccess` declarations whose initial prepare returns
unsupported without mutation. This makes all new tests compile before RED.

- [ ] **Step 1: profiles 4--9 RED**

每个 profile 精确验证 required/forbidden modules。material persists FlowState density/h/
scalar；ideal gas additionally persists `thermodynamic_pressure_pa`、optional target mass、
closure revision。WALE transient fields不得持久化。

- [ ] **Step 2: corruption/rollback RED**

覆盖 wrong profile、missing closure、extra IBM section、CRC、partition、inactive negative
zero、failure after restore prepare。failed read 对 fields/history/controller/p0/immersed
pressure authority bitwise neutral。

- [ ] **Step 3: run codec/transaction RED**

```bash
cmake --build build/debug -j32 --target \
  test_checkpoint_v3_codec test_checkpoint_v3_header_contract \
  test_checkpoint_v3_density_profiles test_ideal_gas_header_contract
ctest --test-dir build/debug --output-on-failure -R \
  '^(test_checkpoint_v3_(codec|header_contract|density_profiles_1_rank)|test_ideal_gas_header_contract)$'
```

Expected: FAIL with profile-presence/closure-restore unsupported assertion；all targets compile and
execute.

- [ ] **Step 4: implement additive profiles and transactional restore**

旧 profile 1 bytes 不变。read prepares FlowState replacement、immersed pressure authority
和 optional ideal closure replacement before one collective ready boundary；publish functions
are `noexcept`/allocation-free/MPI-free. Driver constructs static plans、WALE、initial closure
和 flow before read；successful restore replaces closure state in place, then next attempt
recomputes WALE identity。

- [ ] **Step 5: fast continuation gate and frozen formal registration**

```bash
cmake --build build/debug -j32 --target \
  test_checkpoint_v3_density_profiles test_checkpoint_v3_header_contract \
  test_ideal_gas_closure test_ideal_gas_header_contract hundun
ctest --test-dir build/debug --output-on-failure -R \
  '^(test_checkpoint_v3_(density_profiles_[12]_rank|header_contract)|test_ideal_gas_(closure_1_rank|header_contract)|test_stage3_restart_fast_[12]_rank)$'
```

每个 case 8-cubed、continuous 2 steps vs restart after step 1。主 agent提交
`feat: complete Checkpoint v3 module profiles`。运行 `ctest -N -R
'^checkpoint-continuation-n12-r[124]$'` 确认三条 formal row 可列出，但不得运行它们。

---

### Task S3-O2: Complete Stage 3 providers and inventory

**Ownership:** main agent only。该任务需要修改 flow diagnostic snapshot，不委派。

**Dependency:** C3、O1、R2 accepted。

**Consumes:** profile truth table、authenticated flow/force/WALE reports、Checkpoint presence
and the immutable IBM surface/domain/Ghost/Wall/LFP plans。

**References:** AMReX `EBSupport` 的 presence-driven inventory、初始化/每步成本分离；
HUNDUN diagnostics schema v1 是唯一 record authority。

**Files:**

- Modify: `include/hundun/diag_immersed_module.hpp`
- Modify: `src/diag_immersed_module.cpp`
- Modify: `src/flow_immersed.cpp`
- Modify: `include/hundun/diag_checkpoint_v3.hpp`
- Modify: `src/diag_checkpoint_v3.cpp`
- Create: `include/hundun/diag_immersed_static.hpp`
- Create: `src/diag_immersed_static.cpp`
- Modify: `src/CMakeLists.txt`
- Create: `tests/unit/test_stage3_provider_inventory.cpp`
- Create: `tests/mpi/test_stage3_diagnostics_mpi.cpp`
- Create: `tests/unit/test_immersed_static_diagnostics_header_contract.cpp`
- Modify: `tests/cmake/stage3_diagnostics_registration.cmake`

**Produces:** profile-driven providers for kinds 18--22；authenticated static/per-attempt
snapshots。Exact counter instrumentation and artifacts belong to E1。

**Interface:** `diag_immersed_static.hpp` defines value-only
`ImmersedStaticDiagnosticSummary` and:

```cpp
ImmersedStaticDiagnosticSummary summarize_immersed_static(
    const immersed::ImmersedSurface&, const immersed::SurfaceQuery&,
    const immersed::ImmersedDomain&, const immersed::GhostStencilPlan&,
    const immersed::WallQuadraturePlan&,
    const immersed::LocalFlowPatternTransform&);

std::vector<DiagnosticModuleKind>
stage3_added_provider_inventory(flow::CheckpointV3Presence);
```

The summary owns scalars/IDs only and borrows nothing after construction。Provider overloads for
kinds 18--20 consume that summary；kind 21 consumes authenticated
`ImmersedFlowDiagnosticSource`；kind 22 remains O1 `WaleSummary`。No function accepts a mutable
numerical object。

- [ ] **Step 0: add compile-preserving provider seam**

Create the header/source and exact inventory signature。Initial summarizer/provider collection
returns stable `capability` failure after validating inputs；tests compile and execute before RED。

- [ ] **Step 1: inventory RED**

```text
body-fitted/wale -> les
IBM/none          -> immersed_surface, ghost_stencil, local_flow_pattern, wall_force
IBM/wale          -> all five
absent module     -> no descriptor, instance, counter, fake record
```

- [ ] **Step 2: freeze provider field contract and RED**

| Kind | Stable fields/units | Source authority |
| --- | --- | --- |
| 18 immersed_surface | vertices/triangles/components=count; bbox=m; area=m2; closed-volume=m3; area-vector closure=m2; orientation/fingerprint=count | accepted `ImmersedSurface` + immutable domain summary |
| 19 ghost_stencil | links/donors/QR-rank=count; condition=1; max-halo-reach=cell; wall-points/triangle-coverage=count; plan fingerprints=count | `GhostStencilPlan`/`WallQuadraturePlan` public immutable data |
| 20 local_flow_pattern | algorithm/row fingerprints=count; replacement groups/occurrences=count; coefficient norm=1; limiting-case status=count | published attempt snapshot from `FixedStepImmersedFlow`, never a recomputation |
| 21 wall_force | four pressure/viscous/total vectors=N; moment=N m; area closure=m2; point count=count; lowest rank=count | accepted report + same-attempt `WallForceSample` snapshot |
| 22 les | identity=count; nu_t min/max/l2=m2/s; exact-zero/owned-active=count | accepted `WaleSummary` from O1 |

The snapshot seal includes density variant、force、moment/area and WALE. A rejected attempt cannot
publish a static/per-attempt source. Field IDs and order are frozen in the test, not derived from
`std::map` iteration.

- [ ] **Step 3: disabled/read-only RED**

disabled/not-due collection adds zero allocation/collective/full-field copy/business counter。
重复 local collection canonical identical；collective only when request says collective。

- [ ] **Step 4: 运行 combined RED**

```bash
cmake --build build/debug -j32 --target \
  test_stage3_provider_inventory test_stage3_diagnostics_mpi
ctest --test-dir build/debug --output-on-failure -R \
  '^(test_stage3_(provider_inventory|diagnostics_mpi_1_rank))$'
```

Expected: FAIL，static providers/profile inventory 尚不完整；targets compile and execute。

- [ ] **Step 5: 实现 providers/snapshots**

新增 friend-free static summary adapters；扩展既有 `ImmersedFlowDiagnosticSource` 和
Checkpoint report collection，但不重算 pressure/force/WALE。`flow_immersed.cpp` only copies
already-published immutable values into the source；disabled path 在 snapshot/dispatch 前返回。

- [ ] **Step 6: focused gate**

```bash
cmake --build build/debug -j32 --target \
  test_stage3_provider_inventory test_stage3_diagnostics_mpi \
  test_immersed_static_diagnostics_header_contract \
  test_immersed_diagnostics_header_contract \
  test_checkpoint_v3_diagnostics_header_contract
ctest --test-dir build/debug -j24 --output-on-failure -R \
  '^(test_stage3_provider_inventory|test_immersed_static_diagnostics_header_contract|test_immersed_diagnostics_header_contract|test_checkpoint_v3_diagnostics_header_contract)$'
ctest --test-dir build/debug --output-on-failure -R \
  '^test_stage3_diagnostics_mpi_[12]_rank$'
cmake -S . -B build/stage3-tests-off -DCMAKE_BUILD_TYPE=Release \
  -DHUNDUN_BUILD_TESTS=OFF
cmake --build build/stage3-tests-off -j32 --target hundun
```

主 agent复核没有第二数值 authority 后提交 `feat: complete Stage 3 diagnostics`。

---

## Phase F3 — Framework closure

### Task S3-A1: Complete same-executable combination matrix

**Ownership:** main agent。

**Dependency:** C3、R2、O2 accepted。

**Consumes:** the nine-profile runtime implementations、Checkpoint presence matrix、provider
inventory and accepted schema-v3 validation/broadcast path。

**References:** incflo `incflo_advance.cpp` 的高层 construction/advance 分离；HUNDUN
现有 `app_flow_driver.cpp` 和 accepted 19A 是 CLI/stdout/transaction authority。

**Files:**

- Modify: `src/app_immersed_flow_driver.cpp`
- Modify only if broadcast fields changed: `src/app_case_config_broadcast.cpp`
- Modify: `tests/acceptance/stage3_flow_models_fast.sh`
- Modify: `tests/acceptance/stage3_restart_fast.sh`
- Create: `tests/unit/test_stage3_dispatch_inventory.cpp`
- Modify: `tests/cmake/stage3_framework_registration.cmake`

**Produces:** all nine Checkpoint profiles validate/run/restart/diagnostics through one `hundun`。

- [ ] **Step 1: freeze legal matrix**

Test inventory has exactly the following nine rows；no runtime-generated fallback row：

| ID | Density | IBM | WALE | p0 mode | Checkpoint | Added diagnostics |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | constant | yes | no | n/a | yes | 18--21 |
| 2 | constant | no | yes | n/a | yes | 22 |
| 3 | constant | yes | yes | n/a | yes | 18--22 |
| 4 | material | yes | no | n/a | yes | 18--21 |
| 5 | material | no | yes | n/a | yes | 22 |
| 6 | material | yes | yes | n/a | yes | 18--22 |
| 7 | ideal gas | yes | no | fixed when open, dynamic when closed | yes | 18--21 |
| 8 | ideal gas | no | yes | fixed when open, dynamic when closed | yes | 22 |
| 9 | ideal gas | yes | yes | fixed when open, dynamic when closed | yes | 18--22 |

Mutation: copy row 6 but clear its IBM object while retaining presence value 6；construction must
collectively reject before any output/step. Unknown ID 10 and any incomplete setup are also
pre-step errors。

- [ ] **Step 2: construction-order RED**

```text
config broadcast
-> mesh/boundary
-> optional IBM static plans
-> field registry/state
-> optional ideal-gas closure
-> optional WALE
-> FixedStepImmersedFlow
-> optional restore
-> retry/diagnostics/checkpoint
```

IBM absent 不构造 surface；WALE absent 不构造 model；output/checkpoint 不在 commit 前发布。

- [ ] **Step 3: 运行 dispatch RED**

```bash
cmake --build build/debug -j32 --target hundun test_stage3_dispatch_inventory
ctest --test-dir build/debug --output-on-failure -R \
  '^(test_stage3_dispatch_inventory|test_stage3_flow_models_fast_1_rank)$'
```

Expected: FAIL，至少一个尚未接入的 profile/restart/diagnostic inventory 不匹配；若所有
profiles 已由 earlier vertical slices 接入，则 inventory mutation 必须使 test RED。

- [ ] **Step 4: 实现最终 construction/dispatch wiring**

只补 inventory 证明缺失的组合和顺序。禁止重写已接受 vertical slice；validation 在任何
output/step 前完成，Checkpoint restore 在 static plans/model 构造后、retry loop 前完成。

- [ ] **Step 5: 8-cubed matrix GREEN**

```bash
cmake --build build/debug -j32 --target hundun test_stage3_dispatch_inventory
ctest --test-dir build/debug -j24 --output-on-failure -R \
  '^test_stage3_dispatch_inventory$'
ctest --test-dir build/debug --output-on-failure -R \
  '^test_stage3_(flow_models|restart)_fast_[12]_rank$'
```

- [ ] **Step 6: build/link gates**

```bash
cmake -S . -B build/stage3-tests-off -DCMAKE_BUILD_TYPE=Release \
  -DHUNDUN_BUILD_TESTS=OFF
cmake --build build/stage3-tests-off -j32 --target hundun
if nm -C build/stage3-tests-off/src/hundun |
    rg -q 'TestAccess|ENABLE_TEST_ACCESS'; then
  echo "tests-off hundun contains a test-only symbol" >&2
  exit 1
fi
ldd build/stage3-tests-off/src/hundun
```

`ldd` 只允许系统 C/C++、MPI、thread/dl；无 Python/vendor solver/GPU runtime。

- [ ] **Step 7: full impact review and commit**

codegraphf 检查 `app_main`、all drivers、broadcast、checkpoint、diagnostics callers。提交
`feat: complete Stage 3 application combinations`。

---

### Task S3-E1: Exact counters and performance artifact producer

**Ownership:** main agent only。该任务跨 geometry、flow、WALE、Checkpoint 和 artifact
schema，不委派。

**Dependency:** A1 accepted；所有九 profiles 的 construction/runtime owner 已稳定。

**Consumes:** O2 read-only provider snapshots；现有 performance artifact schema v1；
execution/MPI/Halo/linear-solver counters。

**Produces:** additive performance artifact schema v2；small 8-cubed RED/GREEN harness；
frozen 24-cubed 1/2/4 final selector（本 task 不运行 24-cubed）。

**References:** AMReX initialization-vs-step accounting only；HUNDUN existing
`diag_performance_artifact` is the serialization authority。

**Files:**

- Modify: `include/hundun/diag_performance_artifact.hpp`
- Modify: `src/diag_performance_artifact.cpp`
- Create: `include/hundun/diag_stage3_performance.hpp`
- Create: `src/diag_stage3_performance.cpp`
- Modify: `include/hundun/ib_surface_query.hpp`
- Modify: `src/ib_surface_query.cpp`
- Modify: `include/hundun/ib_domain.hpp`
- Modify: `src/ib_domain.cpp`
- Modify: `include/hundun/ib_ghost_stencil_plan.hpp`
- Modify: `src/ib_ghost_stencil_plan.cpp`
- Modify: `include/hundun/flow_immersed.hpp`
- Modify: `src/flow_immersed.cpp`
- Modify: `include/hundun/les_wale.hpp`
- Modify: `src/les_wale.cpp`
- Modify: `include/hundun/flow_checkpoint_v3.hpp`
- Modify: `src/flow_checkpoint_v3.cpp`
- Modify: `src/CMakeLists.txt`
- Modify: `tests/unit/test_performance_artifact.cpp`
- Modify: `tests/unit/test_task25_performance_header_contract.cpp`
- Modify: `tests/unit/test_immersed_domain_header_contract.cpp`
- Modify: `tests/unit/test_ghost_stencil_header_contract.cpp`
- Modify: `tests/unit/test_stage3_flow_header_contract.cpp`
- Modify: `tests/unit/test_wale_header_contract.cpp`
- Modify: `tests/unit/test_checkpoint_v3_header_contract.cpp`
- Create: `tests/unit/test_stage3_exact_counters.cpp`
- Create: `tests/unit/test_stage3_performance_header_contract.cpp`
- Create: `tests/support/stage3_case_generator.cpp`
- Create: `tests/support/stage3_performance_evidence.hpp`
- Create: `tests/support/stage3_performance_evidence.cpp`
- Create: `tests/mpi/test_stage3_performance.cpp`
- Modify: `tests/cmake/stage3_framework_registration.cmake`

**Artifact compatibility:** add `CounterMap algorithmic_work` to `ExactCounterMaps`。When
`Artifact::schema_version == 1`, canonical JSON remains byte-for-byte unchanged and the new map must
be empty；schema version 2 requires the Stage 3 map and records HEAD/tree/binary、compiler/MPI、
build type、rank/process grid、cpuset/thread budget、case/profile/geometry fingerprints、warmup/
measured steps and repetitions。

**Frozen counter contract:**

| Stable ID | Class | Increment/snapshot owner | Formula for a no-retry measured interval |
| --- | --- | --- | --- |
| `init.surface.triangles` | static snapshot | `ImmersedSurface` | global triangle count, once |
| `init.query.closest-calls` | monotonic work | `SurfaceQuery` | exact calls during construction |
| `init.query.segment-calls` | monotonic work | `SurfaceQuery` | exact segment queries during construction |
| `init.classification.cells` | monotonic work | `ImmersedDomain::create` | cells actually classified, global sum |
| `init.ghost.qr-plans` | monotonic work | `GhostStencilPlan::create` | deterministic QR plans attempted |
| `init.ghost.rejected-plans` | monotonic work | same | candidates rejected before accepted plan |
| `init.ghost.donor-references` | static snapshot | accepted Ghost/Wall plans | exact stored donor references |
| `init.wall.points` | static snapshot | `WallQuadraturePlan` | owned quadrature points, global sum |
| `step.ghost.constraints` | monotonic work | immersed reconstruction | actual constraint evaluations |
| `step.lfp.transforms` | monotonic work | immersed operator | actual row transforms/substitutions |
| `step.immersed.rows` | monotonic work | immersed operator | actual immersed rows assembled/applied |
| `step.pressure.wall-constraints` | monotonic work | two-corrector pressure authority | actual link constraints; accepted no-retry case equals two passes |
| `step.wall.quadrature-evaluations` | monotonic work | final force authority | actual points evaluated |
| `step.force.reductions` | monotonic work | wall-force authority | separately registered pressure/viscous/total reductions |
| `step.wale.gradient-cells` | monotonic work | body-fitted/immersed WALE authority | active cells whose 9-gradient is built |
| `step.wale.evaluations` | monotonic work | `WaleModel::evaluate` caller | exactly one per WALE-enabled attempt |
| `checkpoint.logical-io-bytes` | monotonic work | Checkpoint v3 codec | encoded/decoded logical bytes, not filesystem block size |

Backend `allocated_bytes`、Halo payload/messages、collectives/logical payload、matvec、
preconditioner and I/O maps retain their existing authorities。Every addition is checked for
overflow。Static snapshots never increase after construction。Monotonic work includes work done by
failed attempts and is deliberately not rolled back；it must never influence numerical decisions。
Artifacts bracket snapshots before/after the measured phase, so warmup and unrelated prior work are
excluded。The formal case forbids injected failure/retry；a separate unit mutation proves one failed
attempt adds actual work while FlowState remains bitwise rolled back。

- [ ] **Step 0: compile-preserving schema-v2 seam**

Add the new map and empty counter accessors returning zero snapshots. Preserve schema-v1 bytes.
`test_stage3_performance_header_contract.cpp` includes only
`hundun/diag_stage3_performance.hpp` before the standard standalone-header checks。Register the new
unit/header/MPI targets with L/M resource metadata；the 24-cubed selector is listable but not
invoked。

- [ ] **Step 1: write mutation-sensitive RED**

Cover every stable ID、initialization-vs-step separation、retry semantics、checked overflow、
schema-v1 byte preservation and schema-v2 required metadata。Mutations that double WALE evaluation,
count inactive cells, count provisional flux as final I/O, or rollback actual-work counters fail。

- [ ] **Step 2: run executable RED**

```bash
cmake --build build/debug -j32 --target \
  test_performance_artifact test_stage3_exact_counters test_stage3_performance \
  test_stage3_performance_header_contract test_task25_performance_header_contract \
  test_immersed_domain_header_contract test_ghost_stencil_header_contract \
  test_stage3_flow_header_contract test_wale_header_contract \
  test_checkpoint_v3_header_contract
ctest --test-dir build/debug -j24 --output-on-failure -R \
  '^(test_performance_artifact|test_stage3_exact_counters|test_stage3_performance_header_contract|test_task25_performance_header_contract|test_immersed_domain_header_contract|test_ghost_stencil_header_contract|test_stage3_flow_header_contract|test_wale_header_contract|test_checkpoint_v3_header_contract)$'
ctest --test-dir build/debug --output-on-failure -R \
  '^test_stage3_performance_1_rank_fast$'
```

Expected: unit assertion identifies the first zero/stub counter；MPI target compiles and executes an
8-cubed one-step case, then fails exact-counter comparison。

- [ ] **Step 3: implement minimum instrumentation and producer**

Increment only at the owner operations listed above。No diagnostic collection may itself increment
business/work counters。`stage3_performance_evidence` converts one immutable before/after snapshot
to schema-v2 artifact and uses existing canonical serializer；no second JSON encoder。

- [ ] **Step 4: small GREEN and freeze final selector**

```bash
cmake --build build/debug -j32 --target \
  test_performance_artifact test_stage3_exact_counters test_stage3_performance \
  test_stage3_performance_header_contract test_task25_performance_header_contract \
  test_immersed_domain_header_contract test_ghost_stencil_header_contract \
  test_stage3_flow_header_contract test_wale_header_contract \
  test_checkpoint_v3_header_contract
ctest --test-dir build/debug -j24 --output-on-failure -R \
  '^(test_performance_artifact|test_stage3_exact_counters|test_stage3_performance_header_contract|test_task25_performance_header_contract|test_immersed_domain_header_contract|test_ghost_stencil_header_contract|test_stage3_flow_header_contract|test_wale_header_contract|test_checkpoint_v3_header_contract)$'
ctest --test-dir build/debug --output-on-failure -R \
  '^test_stage3_performance_[12]_rank_fast$'
cmake -S . -B build/stage3-tests-off -DCMAKE_BUILD_TYPE=Release \
  -DHUNDUN_BUILD_TESTS=OFF
cmake --build build/stage3-tests-off -j32 --target hundun
```

`test_stage3_performance_<r>_rank_formal` is registered for ranks 1/2/4, fixed global 24-cubed,
2 warmup + 3 measured steps + 1 repetition, `RESOURCE_LOCK hundun_stage3_performance_m`, timeout
1800 seconds, and writes one canonical artifact per rank group under the runner-provided evidence
directory。Do not run formal mode in E1。

- [ ] **Step 5: main review and commit**

Review every increment caller with codegraphf，prove counters are observational and checked，rerun
Stage 1 performance schema-v1 snapshots。Commit `feat: add Stage 3 exact performance evidence`。

---

### Task S3-G1: Capability ledger, inventory and evidence manifest

**Ownership:** main agent only。ledger、selector inventory 和 projection contract 不委派。

**Dependency:** E1 accepted。

**Consumes:** all accepted task receipts、S1 frozen selector contracts、E1 counter/artifact
schema、existing product-projection manifests and the exact current CTest inventory。

**References:** AMReX/incflo regression metadata 只作可复现 artifact 字段对照；
HUNDUN exact-HEAD manifest、capability ledger 和 CTest inventory 是最终 authority。

**Files:**

- Create: `docs/numerics/stage3-capability-ledger.md`
- Create: `tests/acceptance/stage3_acceptance.sh`
- Create: `tests/acceptance/stage3_acceptance_inventory.tsv`
- Create: `tests/cmake/stage3_acceptance_contract.cmake`
- Create: `tests/cmake/stage3_capability_ledger_contract.cmake`
- Create: `tests/cmake/stage3_product_projection_contract.cmake`
- Create: `tests/cmake/stage3_product_projection.cmake`
- Create: `tests/support/stage3_evidence_manifest.hpp`
- Create: `tests/support/stage3_evidence_manifest.cpp`
- Create: `tests/unit/test_stage3_evidence_manifest.cpp`
- Modify: `tests/cmake/stage3_acceptance_registration.cmake`

**Produces:** complete inventory but does not execute long scientific group during development。

- [ ] **Step 1: 写 inventory/ledger/manifest RED**

先创建三个 CMake contract 和 manifest unit test。它们要求 ledger 每行有 disposition、
task、test、final owner；launcher 必须有六个固定 mode；manifest 缺少 exit/log/binary
identity 必须拒绝；projection contract 必须读取现有
`.superpowers/product-projection-manifest-2026-08-09.tsv` 并拒绝 tests、`.superpowers`、
private path/token 或未登记 product-only override。

```bash
cmake --build build/debug -j32 --target test_stage3_evidence_manifest
ctest --test-dir build/debug --output-on-failure -R \
  '^(test_stage3_(acceptance_contract|capability_ledger|evidence_manifest|product_projection_contract))$'
```

Expected: FAIL，ledger/launcher/manifest implementation 尚不存在；the projection contract also
fails on its frozen illegal-path/absent-final-projector mutation rather than CMake syntax。

- [ ] **Step 2: capability rows**

Every row is exactly `implemented-and-accepted`、`deferred` or `out-of-scope`，并链接 spec、
task commit、test name、final owner。redistribution/multigrid/AMR/GPU/moving bodies/rank-changing
Restart 必须 deferred/out-of-scope，不得遗漏。

- [ ] **Step 3: acceptance launcher modes**

```text
--list
--group <low-cost|scientific|performance|sanitizer|governance>
--candidate-head <40-hex>
--debug-root <absolute-path>
--release-root <absolute-path>
--asan-root <absolute-path>
--ubsan-root <absolute-path>
--tests-off-root <absolute-path>
```

`--list` 不要求 build roots，只输出 canonical inventory。任何 `--group` invocation 都必须
提供 candidate 和五个 absolute roots；launcher 根据 row 的 `build_role` 选择 binary，
不得搜索 `$PATH`、`build/debug` 或最新 mtime。`stage3_acceptance_inventory.tsv` 固定列：
`row_id/group/required/resource_class/ranks/build_role/producer_target/executable/argv/timeout_s/
artifact_subdir`。`build_role` 只能是 `debug`、`release`、`asan`、`ubsan`、`tests-off`
或 `source-only`；scientific/performance 固定使用 `release`。至少含以下 required rows：

| Group | Exact row IDs and command contract |
| --- | --- |
| low-cost | one exact row for every Stage 3 task-focused Debug unit/direct/fast CTest named in its accepted receipt; `task11-authority-current-tree`; `profiles-fast-r{1,2,4}` via `stage3_flow_models_fast.sh`; `restart-fast-r{1,2,4}`; `diagnostics-fast-r{1,2,4}`; Stage 1 low-cost whitelist; Stage 2 core whitelist |
| scientific | `wale-tgv-convergence-r1` via the S1 12/24/48 aggregator; `wale-tgv-n24-r{2,4}`; `wale-channel-n48-r1`; `constant-ibm-wale-n48-r1`; `constant-ibm-wale-n24-r{1,2,4}`; `material-ibm-wale-n{12,24}-r{1,2,4}`; `ideal-ibm-wale-n{12,24}-r{1,2,4}`; `checkpoint-continuation-n12-r{1,2,4}` |
| performance | `stage3-exact-counters-n24-r{1,2,4}` via E1 formal selector |
| sanitizer | focused ASan/UBSan unit + 8/12-cubed one-rank affected paths only |
| governance | tests-off/offline/headers/policy/provenance/DCO/ledger/projection/`nm`/`ldd` |

Sanitizer required row IDs are frozen，all one-rank and at most 8/12-cubed：

```text
asan-immersed-wale-constant-fast
asan-material-wale-composition-fast
asan-ideal-wale-composition-fast
asan-checkpoint-density-profiles-fast
asan-evidence-manifest-unit
ubsan-immersed-wale-constant-fast
ubsan-material-wale-composition-fast
ubsan-ideal-wale-composition-fast
ubsan-checkpoint-density-profiles-fast
ubsan-evidence-manifest-unit
```

48-cubed rows use resource H and timeout 43200 seconds；24-cubed scientific/performance rows use
locked M and timeout 7200/1800 seconds；small MPI uses timeout 600 seconds。All artifacts go under
`${HUNDUN_STAGE3_EVIDENCE_DIR}/<row_id>/`。Launcher validates source HEAD、clean status、cache、
tests enabled/disabled as required、canonical absolute build root、`CMakeCache.txt` SHA、binary
SHA、expected target/argv and free resource lock before execution。It records every row exit and
continues independent rows, but group exit is nonzero if any required row fails or lacks a terminal
manifest。A missing/mismatched root or a binary not built from `--candidate-head` fails before the
row starts。The G1 contract requires every non-`source-only` required row to name a nonempty
`producer_target` present in that build role's CMake target graph；this specifically covers every
V1 Release scientific/performance producer。`--list` may run at any time；V0 alone invokes
low-cost/sanitizer/governance on the frozen candidate，and V1 alone invokes
scientific/performance。

- [ ] **Step 4: manifest schema**

记录 HEAD/tree/diff、build role、canonical build root、CMake cache SHA、binary SHA/inode、
compiler/libc++/MPI、argv/env/cpuset/ranks、start/end、exit、duration、peak RSS、log SHA、
artifact SHA and required-row ID。missing exit status、duplicate row ID 或 root/cache/binary
identity mismatch 不可 accepted。

- [ ] **Step 5: low-cost contract GREEN**

```bash
cmake --build build/debug -j32 --target test_stage3_evidence_manifest
ctest --test-dir build/debug -j8 --output-on-failure -R \
  '^(test_stage3_(acceptance_contract|evidence_manifest|product_projection_contract)|test_stage3_capability_ledger)$'
bash tests/acceptance/stage3_acceptance.sh --list
```

禁止在本 task 调用 `--group scientific`、`--group performance` 或 sanitizer。主 agent提交
`test: define Stage 3 acceptance inventory`。

---

### Task S3-DOC: Public documentation finalization

**Ownership:** main agent only。

**Dependency:** G1 accepted；在 V0 candidate freeze 前完成。

**Consumes:** accepted capability ledger、profile truth table、public configuration/schema、
diagnostic/Restart contracts and task provenance receipts。

**Produces:** public 0.2.0 technical documentation whose claims are mechanically bounded by the
ledger；no acceptance evidence is created by prose。

**References:** 本计划的 public-algorithm-reference 和 task provenance receipts；只写
HUNDUN 已实现的公式/能力，不从上游项目复制说明文字。

**Files:**

- Modify: `README.md`
- Modify: `docs/index.md`
- Modify: `docs/architecture/overview.md`
- Modify: `docs/architecture/data-flow.md`
- Modify: `docs/api/configuration-schema.md`
- Modify: `docs/api/diagnostics.md`
- Modify: `docs/api/restart-schema.md`
- Modify: `docs/user-guide/quick-start.md`
- Modify: `docs/user-guide/restart.md`
- Modify: `docs/numerics/stage3-contracts.md`
- Modify: `docs/numerics/applicability-and-limitations.md`
- Modify: `docs/releases/current-capabilities.md`
- Modify: `docs/verification/accepted-capabilities.md`
- Modify: `docs/verification/conservation-summary.md`
- Modify: `docs/verification/convergence-summary.md`
- Modify: `docs/verification/decomposition-summary.md`
- Modify: `docs/verification/restart-and-rollback-summary.md`
- Create: `tests/cmake/stage3_documentation_contract.cmake`
- Modify: `tests/cmake/stage3_acceptance_registration.cmake`

- [ ] **Step 0: documentation claim contract RED**

Contract parses the capability ledger and requires every public Stage 3 claim to name one of the
nine profiles or an explicit deferred/out-of-scope row；it checks exact config keys、diagnostic
kinds 18--22、Checkpoint values 1--9、units、version wording and permanent 96-cubed exclusion。
Before the draft it must fail on at least one missing current-capability row, not on CMake syntax。

```bash
cmake --preset debug
ctest --test-dir build/debug -R '^test_stage3_documentation_contract$' \
  --output-on-failure
```

Expected: FAIL with the exact missing capability/profile claim ID。Only after this observed RED may
Step 1 edit public documentation。

- [ ] **Step 1: technical draft**

从 capability ledger 生成能力和限制；命令、JSON keys、单位、profile values 与程序一致。
文档把 0.2.0 能力声明明确绑定到“只有 V1 ACCEPT 才发布”的候选条件；已有证据可写成
已完成，final 48-cubed 在 V1 前仍标为 final-candidate gate，不伪造结果。V1 的实际
命令、计数和哈希只进入 governance acceptance report，不在冻结后回写产品文档。

- [ ] **Step 2: language skills in fixed order**

主 agent依次使用 `humanizer-zh`、`shuorenhua`。只处理公共中文技术说明；法律文本、
内部 plans、logs、receipts 不润色。

- [ ] **Step 3: factual/legal recheck**

逐项复核公式、sign、units、CLI、JSON、enum、version、SHA、Apache-2.0、NOTICE。运行：

```bash
git diff --check
ctest --test-dir build/debug -R \
  '^(test_stage3_documentation_contract|source_policy_.*|provenance_.*)$' \
  --output-on-failure
```

主 agent提交 `docs: finalize Stage 3 user documentation`。

---

### Task S3-V0: Code-complete low-cost preflight and freeze

**Ownership:** main agent only；不使用实现 worker。

**Dependency:** DOC accepted。

**Consumes:** code-complete product/test/public-doc tree、all task receipts and G1 acceptance
inventory。

**Produces:** one clean exact-HEAD candidate `C` plus reusable low-cost/sanitizer/governance
manifests；it does not produce the formal scientific verdict。

**References:** HUNDUN accepted task receipts、DCO/source-policy、exact-tree evidence rules；
不以外部 solver 的测试结果替代本项目 preflight。

**Files:** Read entire Stage 3 diff and all task receipts；Modify none。发现 defect 时 V0
立即 REJECT，另建一个有精确文件白名单的 repair packet；repair accepted 后从 V0 Step 1
重启。V0 freezes governance code/test/public-doc candidate `C`；V2 later may add only the listed
governance receipt/manifest files and may change `VERSION` only in product projection。

- [ ] **Step 1: full diff and history review**

Review `Task11 accepted base..candidate`，覆盖 math、public API/ABI、schema、Restart、
diagnostics、ownership、units/sign、MPI、rollback、allocation、tests、copyright、DCO。

- [ ] **Step 2: fresh low-cost build matrix**

```bash
cmake -S . -B build/stage3-final-debug -DCMAKE_BUILD_TYPE=Debug \
  -DHUNDUN_BUILD_TESTS=ON
cmake --build build/stage3-final-debug -j32

cmake -S . -B build/stage3-final-release -DCMAKE_BUILD_TYPE=Release \
  -DHUNDUN_BUILD_TESTS=ON
cmake --build build/stage3-final-release -j32 --target \
  hundun test_stage3_scientific_row_contract test_wale_taylor_green \
  test_wale_body_fitted test_immersed_wale_constant \
  test_material_wale_composition test_ideal_gas_wale_composition \
  test_checkpoint_v3_density_profiles test_stage3_performance

cmake -S . -B build/stage3-final-tests-off -DCMAKE_BUILD_TYPE=Release \
  -DHUNDUN_BUILD_TESTS=OFF
cmake --build build/stage3-final-tests-off -j32 --target hundun

cmake -S . -B build/stage3-final-asan -DCMAKE_BUILD_TYPE=Debug \
  -DHUNDUN_BUILD_TESTS=ON -DHUNDUN_ENABLE_ASAN=ON
cmake --build build/stage3-final-asan -j32 --target \
  test_immersed_wale_constant test_material_wale_composition \
  test_ideal_gas_wale_composition test_checkpoint_v3_density_profiles \
  test_stage3_evidence_manifest

cmake -S . -B build/stage3-final-ubsan -DCMAKE_BUILD_TYPE=Debug \
  -DHUNDUN_BUILD_TESTS=ON -DHUNDUN_ENABLE_UBSAN=ON
cmake --build build/stage3-final-ubsan -j32 --target \
  test_immersed_wale_constant test_material_wale_composition \
  test_ideal_gas_wale_composition test_checkpoint_v3_density_profiles \
  test_stage3_evidence_manifest
```

Sanitizer roots only build inventory-listed small targets；do not run 48-cubed or any
scientific/performance formal selector。Step 2 only configures/builds；all tests execute exactly
once through the inventory in Step 3，avoiding duplicate header/policy/unit runs inside V0。

- [ ] **Step 3: small MPI preflight**

Set `HUNDUN_STAGE3_EVIDENCE_DIR` to an ignored exact-C directory，then run exactly once with all
build roots explicit：

```bash
stage3_candidate=$(git rev-parse HEAD)
export HUNDUN_STAGE3_EVIDENCE_DIR="$PWD/build/stage3-evidence/$stage3_candidate"
stage3_roots=(
  --candidate-head "$stage3_candidate"
  --debug-root "$PWD/build/stage3-final-debug"
  --release-root "$PWD/build/stage3-final-release"
  --asan-root "$PWD/build/stage3-final-asan"
  --ubsan-root "$PWD/build/stage3-final-ubsan"
  --tests-off-root "$PWD/build/stage3-final-tests-off"
)
bash tests/acceptance/stage3_acceptance.sh --group low-cost "${stage3_roots[@]}"
bash tests/acceptance/stage3_acceptance.sh --group sanitizer "${stage3_roots[@]}"
bash tests/acceptance/stage3_acceptance.sh --group governance "${stage3_roots[@]}"
```

low-cost includes all profile 8/12-cubed 1/2 ranks、selected 4-rank decomposition、Restart、
rollback、collective failure、two correctors、authority and Stage 1/2 whitelists。V1 validates and
reuses these exact-C manifests；it does not rerun them。

- [ ] **Step 4: freeze candidate**

Before Step 2 record HEAD/parent/tree/diff SHA and require clean worktree；after Step 3 record every
binary/log SHA and confirm no HUNDUN background process。This identity is frozen candidate `C`。
If any defect is found, leave V0, repair in a new signed commit, invalidate rows by the design table,
and restart Step 1；otherwise tracked product/test/docs remain immutable through V1/V2。

---

## Phase F4 — Final frozen-candidate validation

### Task S3-V1: Single formal scientific matrix

**Ownership:** main agent；可委派默认 worker 只读监看日志，不可修改源码。

**Dependency:** V0 frozen candidate。

**Consumes:** exact candidate `C`、V0 terminal manifests、S1 scientific rows and E1 performance
rows。

**Produces:** terminal scientific/performance manifests and a PASS/REJECT verdict bound to `C`；
no tracked source/test/doc changes。

**References:** WALE paper tensor/y-cubed；Basilisk embedded gradient/force 独立 oracle；
incflo/AMReX-Hydro projection/decomposition tests。所有 pass threshold 仍来自 HUNDUN spec。

**Files:** no source/test/doc edits。输出只进入 ignored evidence directory 和最终 report。

- [ ] **Step 1: validate and reuse V0 groups**

```bash
bash tests/acceptance/stage3_acceptance.sh --list
```

Verify every required low-cost/sanitizer/governance row has a terminal V0 manifest bound to C、the
same test source、build configuration and binary SHA。Do not rerun an identical row merely because
V1 started。

- [ ] **Step 2: run the exact scientific group**

```bash
stage3_candidate=$(git rev-parse HEAD)
export HUNDUN_STAGE3_EVIDENCE_DIR="$PWD/build/stage3-evidence/$stage3_candidate"
stage3_roots=(
  --candidate-head "$stage3_candidate"
  --debug-root "$PWD/build/stage3-final-debug"
  --release-root "$PWD/build/stage3-final-release"
  --asan-root "$PWD/build/stage3-final-asan"
  --ubsan-root "$PWD/build/stage3-final-ubsan"
  --tests-off-root "$PWD/build/stage3-final-tests-off"
)
bash tests/acceptance/stage3_acceptance.sh --group scientific "${stage3_roots[@]}"
```

Scientific inventory：

```text
WALE TGV: 12/24/48 single-rank self-convergence; 24 on 1/2/4 ranks
WALE body-fitted channel: 48 single-rank baseline
constant IBM+WALE: 48 single rank; 24 on 1/2/4 ranks
material IBM+WALE: 12/24 short 1/2/4 ranks
ideal-gas IBM+WALE: 12/24 short 1/2/4 ranks
Checkpoint v3: 12 continuous-vs-restart 1/2/4 ranks
```

driver/diagnostics/profile small rows belong to V0 low-cost manifests and are only revalidated in
Step 1；V1 does not rerun them。

48-cubed 同时只能一个；96-cubed 永不运行。每个 detached unit 使用 manifest 记录 exact
identity。H job 期间允许只读 diff/manifest 审查，不允许改 tracked source/test/docs。

- [ ] **Step 3: run the exact performance group**

```bash
stage3_candidate=$(git rev-parse HEAD)
export HUNDUN_STAGE3_EVIDENCE_DIR="$PWD/build/stage3-evidence/$stage3_candidate"
stage3_roots=(
  --candidate-head "$stage3_candidate"
  --debug-root "$PWD/build/stage3-final-debug"
  --release-root "$PWD/build/stage3-final-release"
  --asan-root "$PWD/build/stage3-final-asan"
  --ubsan-root "$PWD/build/stage3-final-ubsan"
  --tests-off-root "$PWD/build/stage3-final-tests-off"
)
bash tests/acceptance/stage3_acceptance.sh --group performance "${stage3_roots[@]}"
```

Run fixed global 24-cubed on 1/2/4 ranks。Exact algorithmic/backend counters are hard gates；
wall time/RSS/bandwidth/throughput must be positive finite and are comparability metadata only。
Performance M rows do not overlap an H job。

- [ ] **Step 4: verdict**

任何 correctness/conservation/MPI/Restart/rollback failure 都 REJECT；不通过调阈值、
滤波、阻尼或 extra corrector 修复。timing/RSS 只记录，不设跨机器阈值。REJECT exits V1；
create one repair packet, then restart V0。Do not keep later rows and silently substitute a changed
binary。

---

### Task S3-V2: Exact-HEAD seal and product 0.2.0 projection

**Ownership:** main agent only。

**Dependency:** V1 all required rows PASS。

**Consumes:** candidate `C`、all V0/V1 manifests、projection whitelist/base manifest and clean
governance/product repositories。

**Produces:** signed product commit `P` at version 0.2.0、tracked governance seal commit `G` and
distinct recorded identities for C/P/G。

**References:** HUNDUN Checkpoint v3/exact-HEAD/DCO/product-projection contracts；外部项目
不提供 acceptance seal 或 release metadata。

**Files:**

- Create: `.superpowers/sdd/2026-08-09-hundun-flow-stage3-framework-completion/stage3-final-acceptance.md`
- Modify: `.superpowers/sdd/2026-08-09-hundun-flow-stage3-framework-completion/ledger.md`
- Create: `.superpowers/product-projection-manifest-stage3-final.tsv`
- Read: `.superpowers/product-projection-manifest-2026-08-09.tsv`
- Read: `.superpowers/product-stable-identity-allowlist-2026-08-09.tsv`
- Modify in product projection only: `VERSION`
- Project only paths emitted by `tests/cmake/stage3_product_projection.cmake`；never carry tests or
  `.superpowers` into product。

不再改 accepted governance code tree。

- [ ] **Step 1: write one exact-HEAD seal**

First prepare an ignored seal draft bound to
`accepted_governance_code_head=C`。Record parent/tree/diff、all commands/exits/log/artifact hashes、
binary/compiler/MPI、reused V0 evidence、DCO、worktree/process and capability dispositions。
Do not yet write `accepted_product_head`；it does not exist until Step 2 succeeds。

- [ ] **Step 2: project and accept product 0.2.0**

Require `/home/wyf/code_dev/hundun-flow` clean with zero remotes。Create a non-destructive product
linked worktree only if the target path does not exist：

```bash
git -C /home/wyf/code_dev/hundun-flow worktree add \
  -b coast/stage3-product-projection \
  /home/wyf/code_dev/hundun-flow-stage3-product.candidate HEAD
cmake \
  -DHUNDUN_GOVERNANCE_ROOT=/home/wyf/code_dev/.worktrees/hundun-flow-stage3-framework \
  -DHUNDUN_PRODUCT_ROOT=/home/wyf/code_dev/hundun-flow-stage3-product.candidate \
  -DHUNDUN_BASE_MANIFEST=.superpowers/product-projection-manifest-2026-08-09.tsv \
  -DHUNDUN_OUTPUT_MANIFEST=.superpowers/product-projection-manifest-stage3-final.tsv \
  -P tests/cmake/stage3_product_projection.cmake
```

The script resolves both manifest arguments against `HUNDUN_GOVERNANCE_ROOT`，not against the
product candidate or the shell working directory；the product candidate never receives a
`.superpowers` directory。

The script copies only accepted public/product paths from C、preserves the tests-off product preset、
sets product-only `VERSION=0.2.0` and fails rather than deleting an unexpected tracked path。Run
pre-commit tracked-text/private-token/license/NOTICE scan, then independent configure/build/install、
`hundun --version`、template `--validate`/`--print-resolved` and one-step 1-rank smoke in the
candidate。Main agent creates one authorized signed product commit `P`，runs post-commit tracked-text
and complete product-history scans，then fast-forwards the product root：

```bash
git -C /home/wyf/code_dev/hundun-flow merge --ff-only \
  coast/stage3-product-projection
```

No governance history or remote enters product。Do not push or publish。

- [ ] **Step 3: commit the governance seal**

Now write the tracked seal/ledger with distinct fields：

```text
accepted_governance_code_head=C
product_projection_source_tree=<tree of C>
accepted_product_head=P
scientific_matrix_head=C
```

Commit only the three listed governance receipt/manifest files as signed governance commit `G`。
`G` is not a tested product/code HEAD and must never replace C。Because a commit cannot contain its
own future hash, the final coordinator response records `governance_report_commit=G` after commit；
do not amend G merely to insert its own ID。

- [ ] **Step 4: final audit and stop**

Confirm product root at P、governance code C unchanged、governance report G only、candidate
worktree clean、no HUNDUN background/MPI jobs、no push/publication/private access。Only after the
product root fast-forward and manifest hashes match may the exact product candidate worktree be
removed with `git worktree remove`。Report Stage 3 accepted and stop；do not create Stage 4 work。

---

## Execution order and parallel schedule

```text
foundation: A0 -> P0 -> C1 -> create isolated infra worktree
main lane:          D1 -> C2 -> D2 -> C3 -> S1 -----------> integrate -> R2 -> O2 -> A1 -> E1 -> G1 -> DOC -> V0 -> V1 -> V2
infra lane:         R1 -> O1 -> signed handoff commits ----/
```

R1 必须等 C1 public/report interface 稳定，但可与 D1 并行；O1 可与 C2/D2 并行。
R1/O1 only run in the infra linked worktree and use separate registration/build trees。Main agent
signs their handoff commits and integrates them after S1；R2/O2/A1/E1 then run sequentially on main。
If no verified default worker is available，main agent executes R1/O1 in the infra worktree；it does
not substitute a differently configured worker。

Resource scheduler：L 组 build `-j32` / CTest `-j24`，最多两个不同 source/build tree；M 组
8/12/24-cubed 或 MPI 1/2/4，每作业总 CPU/thread budget `<=96`，最多两个且绑定不同
NUMA set；H 组 48-cubed 全机最多一个。任何预计超过 10 分钟的开发期 screen 不属于
task synchronous gate；不得在 V0 前另启 detached long screen，只能正常停止并移交 V1。
禁止在同一 build tree
并发 link；CTest M/H rows use `PROCESSORS` + resource locks；performance M 独占；禁止两个 H
job，禁止 H job 期间重建其 binary。

## Per-task receipt checklist

每个 task 提交前由主 agent记录：

```text
task id and accepted parent
allowed files and actual files
public reference revision / adopted behavior / rejected behavior
RED mutation and observed failure
GREEN commands, exits, binary/log SHA-256
codegraphf callers/impact
API/ABI/schema/Restart/diagnostic impact
MPI/rollback/allocation review
copyright independence statement
DCO and clean worktree
background process state
deferred tests and their final owner
worker baseline / handoff commit / integration parent, when delegated
```

receipt 缺少 RED、调用方或延期 owner 时不得接受 task。

## Spec-to-task coverage

| Controlling requirement | Implementation owner | Direct evidence | Final owner |
| --- | --- | --- | --- |
| approval-bound immutable authority activation | A0 | document hashes、approval receipt、single active profile、DCO | P0 and every later receipt |
| isolated parallel execution and registration ownership | P0 | registration mutation、worktree/build-root proof | every delegated receipt |
| Task 11 pressure/operator/final-flux/force authority | C1--C3 consume without replacement | authority fingerprints、two-corrector、force consistency fast | V1 |
| WALE IBM-aware gradient and frozen `mu_eff` | C1 | one-evaluation、solid-read mutation、wall viscosity fingerprint | V1 |
| material wall density/final-flux conservation | D1 | nonzero-gradient、mass、positive density、rollback | V1 |
| material body-fitted/IBM WALE | C2 | stale-density、two-domain composition、retry identity | V1 |
| ideal-gas active-volume p0 and closure | D2 | open/closed、mass/closure、p0 rollback | V1 |
| ideal body-fitted/IBM WALE and Gate 5 | C3 | stage trace、three-density retry/collective | V1 |
| frozen scientific rows without development long waits | S1 | selector mutations、12-cubed smoke、formal rows list-only | V1 |
| Checkpoint v3 profiles and transaction | R1/R2 | golden bytes、corruption、failed-read、continuous-vs-restart | V1 |
| provider inventory | O1/O2 | enum、absence、authenticated snapshots、read-only/disabled path | V1 |
| same-executable legal matrix | A1 | nine-profile dispatch/restart/diagnostics 1/2-rank fast | V1 |
| exact counters and performance artifacts | E1 | schema-v1 preservation、counter mutations、8-cubed harness | V1 24-cubed 1/2/4 |
| capability traceability and reproducible runner | G1 | ledger、exact inventory、manifest/projection contracts | V1/V2 |
| public user documentation | DOC | documentation/source-policy/provenance checks | V0/V2 |
| scientific convergence/decomposition/sanitizers | no development task waits for it | frozen-candidate matrix only | V1 |
| exact-HEAD acceptance and product 0.2.0 | V0--V2 | seal、DCO、clean tree、independent product build/install/smoke | V2 |

Every original Stage 3 sections 11 and 13--22 has an implementation and a final owner in this
table. No row is closed by a fast case alone.

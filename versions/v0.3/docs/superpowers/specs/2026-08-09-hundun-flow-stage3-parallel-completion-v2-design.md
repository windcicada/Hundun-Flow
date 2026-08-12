# HUNDUN-FLOW Stage 3 并行框架收口设计 v2

**状态：** `PROPOSED_FOR_USER_REVIEW`

**执行 profile：** `stage3-parallel-framework-v2`

**设计基线：** `7fc8c5080528f6ea0dbc787c51ca40d9e0fa4553`

**Task 11 科学基线：** `66080e324089599711fdb26082af9b330bfdb5ce`

**适用范围：** 当前基线之后的 Stage 3 未完成能力；不重解释已经接受的 Tasks 1--13、
17A、18A、19A 和 19B。

## 1. 修订决定

Stage 3 后续改为“一条科学主干、一条隔离的基础设施支线、一个最终长测阶段”。先用
`S3-P0` 一次性拆分新测试的 CMake 注册文件，并建立独立 linked worktree；此后科学
主干只有主 agent 顺序修改共享 flow composition，基础设施支线最多由一个边界明确的
默认 worker 顺序完成 Checkpoint 和 WALE diagnostics。worker 不提交，主 agent 在该
独立 worktree 完整审查并签署 handoff commit，再在科学主干的冻结 integration point
cherry-pick。driver 不再等到旧 Task 19 才统一接入，而是在每个密度 vertical slice
通过后立即增加对应运行组合。

本设计只重排未完成工作：

- 保留 Task 11 的 pressure/operator/final-flux/force 单一 authority；
- 保留成功 attempt 恰好两次 PISO corrector；
- 保留 committed/history/trial、collective failure 和 bitwise rollback；
- 保留 Checkpoint v3 已接受的 byte/protocol 兼容边界；
- 永久禁止 96-cubed；
- 开发 task 不等待 48-cubed 或其他大算例长测；
- 48-cubed 正式运行只属于完整软件冻结后的 final candidate；
- 开发期不启动预计超过 10 分钟的可选 screen；对应 row 移交 frozen-candidate V1，不能
  阻塞下一 task，也不能用旧树历史结果冒充有效 acceptance。

## 2. 当前真实能力

基线已经具备：

- static LFP-GCIBM、共享 pressure/flux/force authority 和 Task 11 科学接受；
- constant-density IBM driver、Checkpoint v3 IBM-only continuation 和最小诊断；
- backend-neutral WALE core；
- body-fitted constant-density WALE、一次求值、冻结 `mu_eff`、12/24 screen；
- repository split、flat public layout、tests-off product build。

未具备：

- IBM 内的 WALE gradient、wall `mu_eff` 和 constant IBM+WALE；
- material-density IBM 的非零法向密度重构和 final-flux transport；
- ideal-gas IBM 的 active-volume `p0` 与 h/T/rho closure；
- material/ideal IBM+WALE 的 `rho_attempt` 联合时序；
- WALE/density/combined Checkpoint v3 profiles；
- 完整 Stage 3 provider inventory 和 exact counters；
- 全合法组合 driver/restart/diagnostics；
- code-complete 后的正式科学矩阵、capability ledger、0.2.0 product projection。

## 3. 选择的架构

### 3.1 组成关系

```text
immutable schema/config
  -> static mesh/boundary/optional IBM plans
  -> FlowState committed/history/trial
  -> geometry route
       body-fitted: accepted constant/material PISO core
       immersed: FixedStepImmersedFlow + ImmersedDensityAttemptAdapter
  -> optional ImmersedWaleAttemptAuthority
  -> predictor -> PISO #1 -> provisional transport/closure
               -> PISO #2 -> final transport/closure
  -> final FaceMassFlux + force + diagnostic sources
  -> prepare commit -> collective publish | rollback
```

`FixedStepImmersedFlow` 仍是唯一 immersed momentum/pressure owner。密度 adapter
只能在冻结的阶段写 trial density/h/scalars、返回 `rho_attempt` 并参与
prepare/publish/rollback；它不能建立第二个 pressure solve、第二份 face flux 或第二份
force。WALE authority 只读 lagged velocity gradient 和 adapter 提供的
`rho_attempt`，每次 attempt 只求值一次。

`domain == nullptr` 不走一套新写的 variable-density PISO。它固定复用既有
`FixedStepConstantDensityFlow` 或 `FixedStepMaterialDensityFlow`；ideal gas 通过既有
`DensityClosureHooks` 接入 material core。为支持 profiles 5 和 8，Task C2/C3 只给
`FixedStepMaterialDensityFlow::attempt_common` 增加私有 WALE hook，复用 Task 13 已接受
的 body-fitted gradient/coefficient builder。旧 public `attempt(..., double mu, ...)`
继续委托到 null hook，Stage 2 行为必须 bitwise 不变。任何为 body-fitted material/
ideal 另写 pressure solve、face-flux finalizer 或 transport loop 的实现均不合格。

### 3.2 密度 setup

在 `include/hundun/flow_immersed.hpp` 增加一个粗粒度、构造期冻结的 setup：

```cpp
struct ImmersedFlowDensitySetup final {
  config::DensityModel model{config::DensityModel::constant};
  const runtime::FieldRegistry* registry{};
  FlowFieldIds fields{};
  std::optional<MaterialDensityTransportSpec> material_transport;
  IdealGasClosure* ideal_gas_closure{};
};
```

规则：

- constant：registry/closure 必须为空、material spec disengaged、`fields` 必须是 canonical
  default identity；旧 `create(...)` overload 委托到该 canonical setup；
- material：registry、fields、material spec 必须同时存在，closure 为空；
- ideal gas：material 三项存在，closure 非空且由 driver 拥有，facade 只借用；
- 所有 raw pointer 都是显式 borrowed collaborator；API 只能验证非空、field identity、
  layout 和 rank 一致性，不能声称运行时证明 C++ lifetime；调用方必须保证 registry、
  closure、mesh、solver、halo 和 execution context 全部比 flow 活得更久，并按 flow →
  closure → registry 的顺序析构；
- 构造期 collective 验证 model、field identity、layout 和 pointer-presence contract；
- setup 不进入 cell kernel，不引入 run-time selection registry 或字符串查找。

`FixedStepImmersedFlow::Impl` 的 body-fitted storage 使用：

```cpp
std::variant<std::monostate,
             FixedStepConstantDensityFlow,
             FixedStepMaterialDensityFlow> body_fitted;
```

constant 构造第一种 flow，material/ideal 构造同一个 material flow；ideal 只另外借用
closure 并在 attempt 时绑定既有 hooks。`domain != nullptr` 时 variant 保持 monostate，
只构造 immersed operator。构造期禁止同时拥有 body-fitted 和 immersed pressure core。

实现细节放在新的 `src/flow_immersed_density_detail.hpp`。私有 adapter 固定提供：

```cpp
begin_attempt(state, stencil)
stage_predictor_transport(state)
evaluate_predictor_closure_if_present(state)
rho_attempt_for_wale()
stage_after_corrector_one(state, provisional_flux)
evaluate_provisional_closure_if_present(state)
finalize_from_corrector_two_flux(state, final_flux)
evaluate_final_closure_if_present(state)
assess_final_and_post_closure(state, final_flux)
prepare_commit()                 // may fail before the collective boundary
publish_commit() noexcept        // allocation-free, MPI-free, cannot fail
rollback() noexcept
```

这些是 coarse-grained orchestration hooks，不是公开 ABI，也不允许在逐 cell 循环中
虚调用。实现可用 `std::variant` 静态分派。

adapter 不另造 density report schema，也不新增一组 friend authority。扩展既有
`detail::DensityClosureBridge` 生成 authenticated material report，扩展已经是
`IdealGasStepAttemptReport` friend 的 `detail::DensityClosureAdapter` 生成 ideal report。
driver、adaptive control 和 diagnostics 继续消费同一 `DensityStepAttemptReport`，不会
出现 immersed-only 的第二报告格式。

bridge 的冻结签名为：

```cpp
static MaterialDensityStepAttemptReport make_material_report(
    StepAttemptReport,
    std::optional<MaterialDensityTransportReport>,
    MaterialTransportFailureReason,
    std::uint64_t material_field_count,
    runtime::FieldId shared_face_mass_flux_field,
    std::uint64_t attempt_identity);

static IdealGasStepAttemptReport make_ideal_gas_report(
    MaterialDensityStepAttemptReport,
    std::optional<IdealGasClosureReport>,
    std::uint64_t attempt_identity);
```

bridge 必须从 nested report 读取 finalization identity、final-flux provenance、residual/
conservation availability 并调用原类型的 `seal()`；不能接受由调用方重复传入、可能与
nested report 不一致的数值。

成功 attempt 的 publish boundary 固定为：FlowState、density transport、closure、
pressure authority 和 diagnostic snapshot 全部先 prepare；一次 collective status
确认所有 rank 均已准备；随后仅执行 `noexcept`、allocation-free、MPI-free 的 swap/scalar
publication。任何 prepare 失败都在 publication 前回滚全部对象。不得在
`FlowState::publish_commit_attempt()` 之后调用可能失败、分配或 collective 的操作。

### 3.3 immersed WALE authority

`src/flow_immersed_wale_detail.hpp` 拥有 attempt-local workspace：

```text
u_lag
IBM-aware 9-component grad(u)
WaleAttemptCoefficients
cell mu_sgs
face mu_sgs
face mu_eff
wall-point mu_eff
WaleSummary
```

顺序固定：

```text
density prediction/closure
-> choose BE/BDF2 u_lag
-> background gradient on regular cells
-> accepted immersed reconstruction on interface cells
-> halo exchange and active-order validation
-> WaleModel::evaluate exactly once with rho_attempt
-> one deterministic face interpolation
-> one wall-quadrature interpolation
```

predictor、两次 corrector、最终 residual、transport diffusion 和 surface force 使用
同一 workspace。失败 attempt 不发布 summary，不把 transient coefficient 写入
FlowState 或 Checkpoint。

body-fitted 与 immersed 只共享数学输入合同，不共享几何 kernel：

- `src/flow_body_fitted_wale_detail.hpp` 从 Task 13 实现中提取已接受的 regular-cell
  lagged gradient、deterministic face interpolation 和 coefficient identity；constant 与
  material/ideal body-fitted core 都调用它；
- `src/flow_immersed_wale_detail.hpp` 只增加 interface reconstruction、active/inactive
  partition 和 wall-quadrature interpolation；
- 两者都调用同一个 `les::WaleModel::evaluate`，且都在 density predictor/closure 后、
  momentum predictor 前求值；
- refactor 前后 constant body-fitted WALE 的 report、fields 和 fingerprints 必须 bitwise
  相同，防止为 profiles 5/8 引入第二套 WALE 数学或改变已接受 Task 13。

### 3.4 Checkpoint v3 兼容策略

17A 已接受 `CheckpointV3Presence::constant_static_ibm = 1`，因此 v2 不把它重写成
新的 tag 布局。采用 additive profile enum，保留值 1 和原字段顺序：

```cpp
enum class CheckpointV3Presence : std::uint8_t {
  constant_static_ibm = 1,
  constant_body_fitted_wale = 2,
  constant_static_ibm_wale = 3,
  material_static_ibm = 4,
  material_body_fitted_wale = 5,
  material_static_ibm_wale = 6,
  ideal_gas_static_ibm = 7,
  ideal_gas_body_fitted_wale = 8,
  ideal_gas_static_ibm_wale = 9
};
```

每个 profile 决定哪些 identity sections 必须存在。没有相关模块时 section count 和
byte count 必须为零。WALE 保存配置和 transient-field schema，不保存 `nu_t`、
`mu_sgs` 或 attempt identity；Restart 后按 history 和 attempted dt 重算。ideal gas
额外保存 dynamic `p0`/target mass/revision。17A 的 constant IBM roundtrip 必须保持
bitwise 和 byte fixture 不变。

write/read 使用两个 const-correct 的 one-call borrowed view；不使用一个同时混合 const
和 mutable pointer 的万能 bundle：

```cpp
struct CheckpointV3WriteModules final {
  CheckpointV3Presence presence;
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
  CheckpointV3Presence presence;
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
```

两套 view 都只在一次调用期间借用对象。read 先把 FlowState、immersed pressure authority
和 ideal-gas closure state 全部准备到临时对象，collective 成功后才执行 no-throw
publish。旧 profile-1 overload 委托给新 overload，字节不变。

冻结 profile truth table：

| Value | Density | IBM section | WALE section | Material fields | Ideal closure | Domain mode | Added diagnostics |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | constant | required | forbidden | forbidden | forbidden | open/closed | 18--21 |
| 2 | constant | forbidden | required | forbidden | forbidden | open/closed | 22 |
| 3 | constant | required | required | forbidden | forbidden | open/closed | 18--22 |
| 4 | material | required | forbidden | required | forbidden | open/closed | 18--21 |
| 5 | material | forbidden | required | required | forbidden | open/closed | 22 |
| 6 | material | required | required | required | forbidden | open/closed | 18--22 |
| 7 | ideal gas | required | forbidden | required | required | open=fixed p0; closed=dynamic p0 | 18--21 |
| 8 | ideal gas | forbidden | required | required | required | open=fixed p0; closed=dynamic p0 | 22 |
| 9 | ideal gas | required | required | required | required | open=fixed p0; closed=dynamic p0 | 18--22 |

`required` 必须恰有一个匹配 section/object；`forbidden` 必须是 null pointer、零 section
count 和零 byte count。profile、schema config、module view 和 payload 任一不一致都在
I/O publication 前 collective 拒绝。

每个 profile 的唯一 owner 冻结如下；后续 task 只能实现或验证自己的列，不能把缺口
推给 A1 临时补写：

| Value | Runtime/science owner | Codec owner | Provider components consumed | Unique provider-integration owner | Driver/inventory owner | Final owner |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | accepted Task 11 + 19A constant IBM | accepted 17A | O2 kinds 18--21 | O2 | accepted 19A；A1 revalidates | V1；V2 seal |
| 2 | accepted Task 13 + 19B body-fitted constant WALE | R1 | O1 kind 22 | O2 | accepted 19B；A1 revalidates | V1；V2 seal |
| 3 | C1 constant IBM+WALE | R1 | O1 kind 22 + O2 kinds 18--21 | O2 | C1；A1 revalidates | V1；V2 seal |
| 4 | D1 material IBM | R2 | O2 kinds 18--21 | O2 | D1；A1 revalidates | V1；V2 seal |
| 5 | C2 body-fitted material WALE | R2 | O1 kind 22 | O2 | C2；A1 revalidates | V1；V2 seal |
| 6 | C2 material IBM+WALE | R2 | O1 kind 22 + O2 kinds 18--21 | O2 | C2；A1 revalidates | V1；V2 seal |
| 7 | D2 ideal-gas IBM | R2 | O2 kinds 18--21 | O2 | D2；A1 revalidates | V1；V2 seal |
| 8 | C3 body-fitted ideal-gas WALE | R2 | O1 kind 22 | O2 | C3；A1 revalidates | V1；V2 seal |
| 9 | C3 ideal-gas IBM+WALE | R2 | O1 kind 22 + O2 kinds 18--21 | O2 | C3；A1 revalidates | V1；V2 seal |

O1 is a bounded kind-22 component producer，not a profile-level integration owner；O2 consumes it
and owns the complete presence-driven provider set for every profile。A1 owns only complete-matrix
dispatch verification and missing wiring proven by its RED；it does not own a fallback numerical
implementation、codec or provider。

### 3.5 diagnostics

在 `DiagnosticModuleKind::performance = 17` 后追加：

```cpp
immersed_surface   = 18,
ghost_stencil      = 19,
local_flow_pattern = 20,
wall_force         = 21,
les                = 22
```

旧 0--17 和 schema v1 不重编号。provider inventory 由 profile 决定：不存在的模块
不注册 provider，也不输出伪造零记录。diagnostics 只读已发布 report/authority；disabled
或 not-due 路径不采样、不 collective、不复制全场。

`ImmersedFlowStepAttemptReport` 的 diagnostic authentication 必须覆盖三种
`DensityStepAttemptReport` variant、四字段 force 和 optional WALE summary。C1、D1、
D2 分别在引入相应 report 时扩展 seal 并加入 mutation RED；O2 不得在最后才补一个会使
早期合法 profile diagnostics 永远 stale 的旁路。

## 4. 新阶段与任务顺序

### Phase F0：已接受基础

Tasks 1--13、17A、18A、19A、19B。基线是 `7fc8c508...`，不重跑其长矩阵。

### Phase A0：批准绑定的激活事务

用户批准后先执行 `S3-A0`：把 immutable design/reference/plan 的 SHA-256 和批准原文写入
tracked activation receipt，并在同一签署治理提交中切换 AGENTS/ledger。A0 不修改三份
被哈希文档，不修改产品/测试，也不产生 worker reading exemption 之外的新权限。没有
A0 clean commit 时 P0 和后续任务均不得开始。

### Phase F0.5：并行执行基础

`S3-A0` 已激活 bounded reading exemption；`S3-P0` 只拆分后续新增测试的 CMake
registration，并冻结独立 worktree/build 规则；不移动既有测试、不改产品、不运行数值
矩阵。C1 accepted 后创建：

```text
main:  /home/wyf/code_dev/.worktrees/hundun-flow-stage3-framework
infra: /home/wyf/code_dev/.worktrees/hundun-flow-stage3-infrastructure
```

infra 从 C1 accepted HEAD 建立，使用 branch `coast/stage3-infrastructure-lane` 和自己的
`build/debug`。worker 只在 infra 工作树修改 Checkpoint/diagnostics 分片；main 工作树只
修改 science 分片。

### Phase F1：科学组合主干

1. `S3-C1`：constant IBM+WALE，冻结 immersed gradient 和 wall `mu_eff` authority。
2. `S3-D1`：material-density IBM + driver vertical slice。
3. `S3-C2`：material-density body-fitted/IBM+WALE；body-fitted 复用 material PISO
   private WALE hook，IBM 复用 density adapter。
4. `S3-D2`：ideal-gas IBM + driver vertical slice。
5. `S3-C3`：ideal-gas body-fitted/IBM+WALE；body-fitted 复用 material PISO + closure
   hooks，关闭旧 Gate 5。
6. `S3-S1`：冻结不在开发期运行的 TGV 12/24/48、自收敛、channel 和 combined final
   scientific selectors；只同步运行 12-cubed smoke。

C1--C3 与 D1--D2 这五个科学组合任务顺序执行，因为都修改
`flow_immersed.cpp`；S1 在它们之后冻结 test-only selector，不修改产品 flow
composition。旧 Task 14、15、16 的映射为：

| 新任务 | 旧边界 |
| --- | --- |
| S3-C1 | Task 16 constant subcluster，提前 |
| S3-D1 | Task 14 + 对应 19C driver 子集 |
| S3-C2 | Task 13 variable-density extension + Task 16 material subcluster |
| S3-D2 | Task 15 + 对应 19C driver 子集 |
| S3-C3 | Task 13 variable-density extension + Task 16 ideal-gas subcluster；Gate 5 verdict |

### Phase F2：基础设施支线

支线最多一个实现 worker，任务之间顺序执行；可与 F1 主 agent 工作并行：

1. C1 后在 infra worktree 执行 `S3-R1`：Checkpoint v3 additive constant WALE profiles
   和 codec RED；不改 driver。
2. `S3-O1`：WALE diagnostics、enum append 和 presence inventory；不改 flow orchestration。
3. S1 完成后，主 agent 把 R1/O1 的已签署 handoff commits 集成到 main。
4. `S3-R2`：material/ideal/combined continuation、closure restore 与 driver wiring；main-only。
5. `S3-O2`：完整 IBM/density/WALE providers；main-only。

### Phase F3：框架收口

1. `S3-A1`：全合法 driver/restart/diagnostics matrix，旧 Task 19C 收口。
2. `S3-E1`：exact-counter/performance artifact producer 和 24-cubed final selector；开发期
   只运行 8-cubed harness RED/GREEN。
3. `S3-G1`：capability ledger、acceptance inventory、evidence manifest 和 product
   projection contract；不运行 scientific/performance 长组。
4. `S3-DOC`：公共文档技术初稿、语言技能、主 agent 事实/法律复核。
5. `S3-V0`：code-complete low-cost preflight，修完后冻结 candidate。

### Phase F4：冻结候选验收与投影

1. `S3-V1`：唯一正式科学矩阵；软件已经 code-complete。
2. `S3-V2`：exact-HEAD seal、Task 21 verdict、product 0.2.0 projection。
3. 停在 Stage 3 边界，不进入 Stage 4。

## 5. 依赖与并行图

```text
accepted F0 -> A0 -> P0 -> C1
  |
  +-> main:  D1 -> C2 -> D2 -> C3 -> S1 ------------+
  |                                                  |
  +-> infra: R1 -> O1 -> signed handoff commits -----+-> integrate
                                                          |
                                                          v
                                                   R2 -> O2 -> A1 -> E1
                                                                         |
                                                                         v
                                                                  G1 -> DOC -> V0
                                                                                |
                                                                                v
                                                                           V1 -> V2
```

允许并行：

- 主 agent 在 main linked worktree 做 D1--C3 时，一个默认 worker 在 infra linked
  worktree 做 R1/O1；两者使用不同 source/build tree 和不同 registration fragment；
- L 类 build/header/unit 与一个 M 类 8/12-cubed MPI case 并行；
- code-complete 后，V1 的一个 H 类 48-cubed job 可与只读 diff/manifest 审查并行。

不允许并行：

- 两个任务同时修改 `flow_immersed.cpp`、`app_immersed_flow_driver.cpp`、同一 test
  registration fragment 或同一 build tree；
- worker 直接在 main worktree 工作，或 main 把 infra 未审查的 dirty diff 混入提交；
- R2 与 A1 同时修改 Checkpoint/driver wiring；
- 任意两个 H 类数值作业；
- final candidate H job 期间修改产品或测试源码。

## 6. 测试与长作业政策

### 6.1 task gate

每个开发 task 只同步等待：

- mutation-sensitive unit/direct RED；
- affected Debug unit/header/policy；
- 一个 8-cubed 或 12-cubed同产品路径 fast；
- collective 改动时的 1/2-rank；
- public header 改动时 standalone header；
- build graph 或 test seam 改动时 tests-off。

Release、ASan、UBSan 只对受影响的小型路径各取一个代表，不组成笛卡尔积。task
commit 不等待 24/48 refinement、工程统计或完整 1/2/4 matrix。

### 6.2 milestone diagnostics

12/24 两层 screen 只有在保守估计 10 分钟内完成时才可作为同步 task gate。预计更长、
首次运行越过 10 分钟或需要 48-cubed 时，开发期不继续等待、也不另启 detached 长作业；
立即正常停止该 screen，保留日志并把对应 row 移交 V1。这样 formal/long evidence 只有
一个冻结 binary owner，不会与正在变化的主干竞争资源。

### 6.3 final matrix

只有 V0 完成、公共文档和测试源码冻结后才开始：

- WALE TGV 12/24/48 单-rank self-convergence；24-cubed 1/2/4；body-fitted channel
  48-cubed 单-rank baseline；
- constant IBM+WALE 48-cubed 单 rank；24-cubed 1/2/4；
- material/ideal IBM+WALE 12/24 短程 1/2/4；
- Checkpoint v3 continuous-vs-restart 12-cubed 1/2/4；
- driver/diagnostics/presence 小型 1/2/4；
- exact-counter/performance artifact：24-cubed 1/2/4，同一 case 的 bytes/messages/
  queries/quadrature/evaluations/collectives/matvec/I/O 为硬门；wall time、RSS、bandwidth
  和 throughput 只要求 positive finite 并记录兼容性元数据，不设置跨机器阈值；
- complete affected Debug、focused Release、小型 focused ASan/UBSan；
- Stage 1 低成本完整回归、Stage 2 core whitelist；
- tests-off、offline、headers、policy、provenance、DCO、`nm`、`ldd`。

任意时刻只有一个 48-cubed/H 作业。96-cubed 不注册、不运行、不作为未来默认。
其中 Debug/Release/sanitizer/governance 小型证据由 V0 在 exact candidate `C` 上各运行
一次，V1 校验并复用，不因进入正式矩阵而重复执行；V1 只新增 scientific 和
performance groups。

### 6.4 resource groups and non-blocking scheduler

本机资源基线为 256 logical CPUs、128 physical cores、251 GiB memory：

| Group | Work | Limit | Concurrency |
| --- | --- | --- | --- |
| L | compile、header、policy、non-MPI unit | build `-j32`、CTest `-j24` | 最多两个不同 source/build tree |
| M | 8/12/24-cubed、MPI 1/2/4 | 单作业总 CPU/thread budget `<=96` | 最多两个，绑定不同 NUMA CPU set；performance 测量独占 M |
| H | 48-cubed、formal long selector | 单作业独占一个 NUMA resource group | 全机最多一个 |

同一 source tree 有 active edit 时不启动另一个构建；同一 build tree 永不并发两个
build/link。M/H CTest 注册必须设置 `PROCESSORS`、固定 resource lock 和 timeout；调度器
不能用一个宽泛 `ctest -j24` 同时放出多个 MPI job。task 的同步 gate 若预计超过 10 分钟，
必须停止并移到 V1，主开发继续。V1 的 detached units 使用独立、clean、exact-HEAD
worktree/build tree，记录 HEAD/tree/binary SHA/command/env/cpuset/log/exit/RSS；V1 期间
tracked code/test/docs 已冻结。

## 7. 证据失效规则

| 变化 | 失效证据 |
| --- | --- |
| 仅 governance 文档/索引 | 不失效数值/二进制证据 |
| public header 或 build graph | affected build/header/tests-off/linkage |
| flow orchestration | 消费该顺序的 fast、MPI、Restart、combined matrix |
| density adapter | 对应 density/closure/conservation/rollback |
| WALE gradient/viscosity | WALE、combined、wall force、Restart recomputation |
| Checkpoint codec | 对应 profiles 的 codec/corruption/continuation |
| diagnostics only | provider/canonical/counter；不重跑数值收敛 |
| selector/test source | 仅该 selector 消费的证据 |
| shared registration/support/launcher | 所有直接或传递消费该 harness 的 test rows |
| counter producer/formula | exact-counter unit、MPI performance artifact 和 capability row |
| artifact schema/parser | artifact contract 与 seal；数值结果在原 raw log/identity 可验证时保留 |
| compiler/MPI/build type | binary-bound evidence；source-only review 不失效 |
| cpuset/thread/rank/process grid | performance comparability；正确性证据按 selector contract 判定 |

Task 11 authority 未改时继续复用其科学接受；若 C1--C3 修改其 pressure、operator、
final-flux 或 force authority，必须先停止并由主 agent 判定是否触发 Task 11 证据失效。
V1 任一必需 row 失败时退出 V1，建立独立 repair packet；修复后从 V0 重新冻结，并只按
本表重跑失效 rows，不能在原 frozen candidate 上边改边继续长测。

## 8. 开源语义复用边界

每个任务使用 `docs/references/2026-08-09-hundun-flow-stage3-public-algorithm-reference.md`
中的固定 revision 和 reference points。采用顺序固定为：

```text
公开数学/行为
-> HUNDUN-owned equation and interface note
-> independent mutation RED
-> HUNDUN naming/layout/control flow implementation
-> provenance receipt
```

不复制、翻译或机械改写上游代码。AMReX/AMReX-Hydro/incflo/gslib 虽为宽松许可，
本阶段也只做语义参考；OpenFOAM/Basilisk 为 GPL，只能比较公式、责任边界和测试思想。
PETSc/Trilinos 只用于评估未来 solver 方向，不成为 Stage 3 依赖。

## 9. Worker 合同

worker 每次只收到一个冻结 packet，至少包含：

- exact baseline HEAD；
- allowed files；
- consumes/produces 的完整签名；
- mutation list 和预期 RED 失败原因；
- fast commands、资源类别和 timeout；
- forbidden files/behaviors；
- 公开 reference URL/revision 和“只复用什么”；
- 返回格式：diff 摘要、测试 exit、未解决问题；
- 明确禁止 commit、sign-off、扩 scope、改阈值和启动长测。

只有 R1、O1 和计划中另行写明的纯测试/机械子包可派发；main-only task 不得把整段
`Files` 清单复制给 worker。每个 worker packet 必须另列比 task 白名单更窄的
`allowed_files`，并指定独立 worktree、独立 build tree、worker baseline 和预期
integration parent。worker 不提交；主 agent 在 worker worktree 审查后创建签署 handoff
commit，集成到 main 后重新运行消费该 integration HEAD 的 GREEN。

本设计获用户批准前不产生 worker reading exemption。批准时 governance receipt 必须记录
design/plan/reference SHA-256，并原子更新 `AGENTS.md`：冻结的 bounded packet 只需完整读
`AGENTS.md` 的安全/版权边界、本设计、公开参考和自身 task block，不必重读已被 packet
摘要且不可修改的历史规格；科学主干、跨模块组合、完整 diff 和最终验收仍必须由主 agent
读取全部 authority。没有该 activation receipt 时不得派发 v2 worker。

主 agent 必须重新检查 worker diff、实际调用方和测试，不以 worker 报告代替验收。

## 10. 明确延期

不阻塞 Stage 3：strict ordinary-host zero allocation、replicated-vector 消除、peer-only
极端规模优化、geometric multigrid、AMR、moving/multiple bodies、FSI、production GPU、
rank-changing Restart、wall functions、chemistry/species/spray/particles。

不能延期：基本网格崩溃/NaN、质量或动量明显失控、MPI mismatch/deadlock、rollback、
Restart corruption、重复 numerical authority、错误 force sign、错误 corrector count、
用户无法通过同一 `hundun` 运行合法组合。

## 11. 完成条件

Stage 3 只有在 A0、P0、C1--C3、D1--D2、S1、R1--R2、O1--O2、A1、E1、G1、DOC、V0--V2
全部 accepted，24-cubed exact-counter artifact 和 final exact-HEAD matrix 通过，
governance/product 0.2.0 投影签署且无遗留后台/MPI 作业时完成。最终能力声明必须逐行
来自 capability ledger，不能把 deferred 项目写成已验证。

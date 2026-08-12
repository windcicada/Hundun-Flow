# HUNDUN-FLOW Stage 3 框架完成执行计划

> **后续修订候选：**
> `2026-08-09-hundun-flow-stage3-parallel-completion-v2.md` 重新拆分并排序
> Task 13+19B 之后的工作。用户批准前继续以本计划为执行 authority；批准后，v2
> 取代本文第 7--11 节中尚未执行的任务顺序和测试调度，历史正文保留不删除。

> **执行状态：** ACTIVE
> **Profile：** `stage3-compact-scientific-v1`
> **Branch：** `coast/stage3-framework-completion`
> **Worktree：** `/home/wyf/code_dev/.worktrees/hundun-flow-stage3-framework`

本计划接续并保留 2026-07-27 原计划、2026-08-05 科学优先修订和 2026-08-08 semantic-port 设计。旧文件不删除；发生冲突时，本计划只在仓库拓扑、Task 11 后顺序、每 task 测试成本、里程碑和最终矩阵 owner 上优先。

## 0. 冻结基线

- Task 11 product/code HEAD：`66080e324089599711fdb26082af9b330bfdb5ce`
- Task 11 result：`CORE_ACCEPT`
- Task 11 scientific convergence：`ACCEPTED_FOR_CURRENT_REQUIREMENTS`
- governance split seal：`ee4d2b18d0c68b3080edcc2e132045175961cfb8`
- initial product HEAD：`ae3d08bbb220d1d3b28ec070d1cba9c33fb85877`
- initial product version：`0.1.0`
- 96-cubed：永久禁止

迁移、目录切换、8 个 linked-worktree 指针修复和产品首次签署提交已经完成。旧 Stage 3 脏工作树只读保留。

## 1. 通用 task 协议

每个 task 使用一个冻结 evidence packet：目标、数学/行为合同、allowed files、mutation、fast 命令和失效规则。执行顺序为：

```text
主 agent 冻结 packet
-> mutation RED
-> bounded implementation（必要时使用用户当前指定的默认 worker）
-> affected build/unit/header/policy
-> 12^3 fast（需要数值路径时）
-> 1/2-rank（只在 collective 改动时）
-> 主 agent requirements + quality + callers + complete diff
-> 签署 task commit
-> receipt
```

task receipt 至少记录 parent/HEAD/tree、allowed files、测试命令和 exit、binary/log SHA、复用证据、DCO、worktree 和进程状态。

## 2. MVP-1：19A constant IBM driver

### 目标

让同一个 `hundun` 可执行程序真正运行 schema v3 constant-density static IBM，并复用 Task 11 的唯一 pressure/operator/final-flux/force authority。

### 新接口

```cpp
int run_immersed_flow_case(
    const CliOptions&,
    runtime::MpiContext&,
    const config::ImmersedFlowCaseConfig&,
    const std::filesystem::path& authoritative_case_root);
```

### 初始 allowed files

- 新增 `src/app_immersed_flow_driver.cpp`
- 新增 `src/app_immersed_flow_driver_detail.hpp`
- 修改 `src/app_main.cpp`
- 修改 `src/app_case_config_broadcast.cpp` 及必要 detail 头
- 修改 `src/CMakeLists.txt`
- 新增/修改直接 unit、MPI smoke、fixture 和 CMake registration

超出该清单必须由主 agent 证明构造依赖后追加，不允许顺带改变 `FixedStepImmersedFlow` 数值实现。

### 构造顺序

```text
load/broadcast immutable schema v3
-> decomposition/topology/geometry/boundaries
-> STL surface + SurfaceQuery
-> ImmersedDomain + active layouts
-> GhostStencilPlan + wall quadrature
-> fields/access plans/state
-> Task 11 FixedStepImmersedFlow
-> fixed/adaptive trial
-> exactly two PISO correctors
-> final flux/force report
-> collective accept or rollback
```

### RED/mutations

- CLI 把 schema v3 错误落入 schema v2 driver；
- geometry/plans 在每步重建；
- force 从 provisional pressure 或第二 authority 读取；
- third corrector；
- rank 局部失败未集体回退；
- output 在 step commit 前发布；
- IBM absent 或 WALE present 被 19A 非法接受。

### fast gate

- schema v3 validate/print-resolved；
- `8^3` 或 `12^3` constant IBM 一步 smoke；
- force 四字段、corrector count=2、finite residual/continuity；
- rollback/collective failure 小型 1/2-rank；
- tests-off build。

19A 不实现 Restart、完整 diagnostics、WALE、material 或 ideal gas。

## 3. MVP-2：17A Checkpoint v3 IBM-only

### 目标

为 19A constant IBM 冻结可事务恢复的 Checkpoint v3，而不是把 IBM 状态塞入 v2。

### 新接口和文件

- `include/hundun/flow_checkpoint_v3.hpp`
- `src/flow_checkpoint_v3.cpp`
- `src/flow_checkpoint_v3_detail.hpp`
- 直接 unit/MPI/codec fixture 和 `src/tests` CMake registration

### 协议

- canonical presence tag：IBM present、WALE absent、constant density；
- header/version/endian、rank/process-grid/owned-box identity；
- schema/config/geometry/static-plan fingerprint；
- committed flow fields/history/time controller/final flux；
- per-rank payload CRC-64；
- temporary publish，所有 rank 完成后最后写 `COMPLETED`；
- read 到临时对象，全部校验成功后一次提交；
- failed read 对调用方 state bitwise neutral；
- 相同分区 continuous-vs-restart bitwise continuation。

### mutations

错误 presence、漏 geometry fingerprint、CRC 绕过、提前 publish、partial read commit、rank identity 忽略、attempt-local force/gradient 被持久化。

### fast gate

小型 codec corruption、failed-read rollback、`8^3/12^3` continuous-vs-restart、1-rank；只有 collective 实现改变时补 2-rank。

## 4. MVP-3：18A minimal diagnostics/counters

### 目标

提供用户判断运行是否可信所需的最小 read-only 诊断，并在性能优化前建立 exact counter substrate。

### 文件

- `include/hundun/diag_immersed_module.hpp`
- `src/diag_immersed_module.cpp`
- `include/hundun/diag_checkpoint_v3.hpp`
- `src/diag_checkpoint_v3.cpp`
- 直接 provider/session/counter tests 与 CMake registration

### 字段

- continuity、final momentum residual；
- maximum/mean wall penetration；
- operator force、budget reaction、surface traction、consistency；
- PISO corrector count；
- failure category 和 lowest failing rank；
- construction：surface queries、classified cells、donors、static bytes；
- per attempt：model evaluations、matvec、reductions、halo bytes/messages、allocations。

disabled path 必须不采样、不分配、不改变产品计数；provider 只能读取已发布 authority。

## 5. MVP milestone

19A、17A、18A 分别通过 task gate 后冻结 MVP commit，运行：

- constant IBM `12^3`；
- 1/2/4-rank 小型分解；
- Restart bitwise；
- rollback 和 collective failure；
- exactly two correctors；
- pressure/operator/final-flux/force authority；
- tests-off、public headers、linkage。

通过后 0.1.x governance 能力可标记为 `MVP_ACCEPTED_INTERNAL`，但 product 不同步，直到 Task 21。

## 6. WALE cluster

### Task 12

新增：

- `include/hundun/les_wale.hpp`
- `src/les_wale.cpp`
- tensor oracle、exact `+0.0`、`y^3`、旋转、量纲、stale-density mutations。

保留类型：`WaleModel`、`WaleAttemptCoefficients`、`WaleCoefficientIdentity`、`WaleSummary`。kernel 不做字符串查询或虚调用，不读取 inactive cells。

### Task 13 + 19B

一次完成 body-fitted variable-viscosity WALE 和 `none/wale` driver。每个 trial 只求值一次，所有消费者共享同一 `mu_eff`；molecular-only 路径 bitwise 不变。

WALE milestone：tensor/`y^3`、12/24 TGV screen、1/2-rank、focused ASan/UBSan。

## 7. Density and combined cluster

### Task 14

material-density IBM；非零壁面法向密度梯度 RED，最终 `FaceMassFlux`、正密度、守恒和 collective rollback。

### Task 15

ideal-gas IBM；active-volume `p0`、open/closed domain、h/T 单一 authority、closure rollback。

### Task 16

三密度 IBM+WALE 共用固定 trial 顺序。WALE 单求值、最终通量、wall `mu_eff`、retry identity 和两次 PISO 是不可拆硬门。

14、15、16 顺序执行。Combined milestone 使用 12/24 constant wake、三密度小型 1/2/4、守恒/closure/retry。

## 8. Framework completion

- 17B：所有合法 IBM/WALE/density presence continuation；非法组合显式拒绝。
- 18B：完整 providers、session、schema 和 exact counters。
- 19C：全部合法 driver combinations、validate/run/restart/output。
- 20：capability ledger、exact performance evidence、最终 launcher；只对 measurement RED 证明的问题做窄性能修复。
- 公共文档：Task 20 后技术初稿，依次使用 `humanizer-zh`、`shuorenhua`，主 agent 再做事实和法律复核。
- 21：冻结 exact HEAD 后的唯一最终 acceptance owner。

## 9. Final compact matrix

- 复用 Task 11 科学结论；当前树重跑低成本 authority、rollback、decomposition、smoke。
- WALE：12/24/48 单 rank，24 的 1/2/4。
- constant IBM+WALE：48 单 rank，24 的 1/2/4。
- material/ideal IBM+WALE：12/24 短程 1/2/4 correctness、conservation、closure、retry。
- Checkpoint v3：12 continuous-vs-restart bitwise，1/2/4。
- driver/diagnostics/presence：小网格 1/2/4。
- exact counters：24；wall time/RSS 只记录。
- 完整受影响 Debug、focused Release、小型 focused ASan/UBSan。
- Stage 1 低成本完整回归、Stage 2 core whitelist。
- tests-off、offline、headers、policy、provenance、DCO、ledger、`nm`、`ldd`。

最终候选运行 H 组期间不得修改产品或测试源码。没有第二轮重复 full matrix。

## 10. 提交和产品同步

- worker 不提交、不添加 sign-off；主 agent 使用已授权身份签署 task commit。
- 不改写 Task 11 及更早接受历史。
- governance receipt 指向实际测试的 code HEAD；治理提交不能冒充数值产品 HEAD。
- Task 21 接受后，从 accepted governance code tree 重新生成 product `0.2.0` 投影，执行提交前后 tracked-text/full-history scan、独立 build/install/smoke，再创建单一签署 product commit。
- 不 push、不发布。

## 11. 立即执行顺序

1. 完成 post-switch receipt、authority docs 和 SDD ledger；
2. 主 agent 用 `codegraphf` 与 `rg` 冻结 19A 调用图、构造顺序、allowed files 和 mutations；
3. 写 19A mutation-sensitive RED，确认失败后才实现 driver。

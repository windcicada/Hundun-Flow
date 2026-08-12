# Task 19A：constant IBM driver 冻结证据包

状态：`ACCEPTED`

基线：`f037fd7`（分支 `coast/stage3-framework-completion`）

## 目标与边界

让同一个 `hundun` 可执行程序读取、广播并运行 schema v3 的静止
LFP-GCIBM constant-density case。运行必须复用 Task 11 已接受的
`FixedStepImmersedFlow` 数值路径，并保持恰好两次 PISO corrector、四字段
signed-force、collective failure 和 rollback 语义。

19A 不实现 Checkpoint v3、完整 diagnostics、WALE、material density、ideal
gas、移动壁面或新的数值算法。`immersed_boundary=none` 与 `les=wale` 在正常
运行模式下必须明确拒绝；它们分别由 19B/19C 接入。

## 已核对调用与影响闭包

- v3 权威配置：`ImmersedFlowCaseConfig` / `ResolvedCaseV3`，定义于
  `include/hundun/cfg_resolved_case_v3.hpp`；loader 和 canonical serializer
  已由 `cfg_resolved_case_v3_loader.cpp` 提供。
- v3 collective 广播：`broadcast_resolved_case_v3`，定义于
  `src/app_resolved_case_v3_broadcast.cpp`；既有 1/2-rank fault-sensitive 测试。
- CLI 当前只调用 v1/v2 loader，且只分发 `CaseConfig` 和 `FlowCaseConfig`。
- Task 11 产品内核：`FixedStepImmersedFlow::create/attempt`，codegraphf 影响
  闭包覆盖 `flow_immersed.cpp`、transaction、PISO、signed-force 和 laminar
  IBM tests。19A 不修改该闭包内的数值实现。
- 构造 oracle 只使用本仓库已接受的 Task 11 tests/support 和现有
  `app_flow_driver.cpp`；不访问或采用私有参考源码。

## 文件白名单

产品文件：

- 新增 `src/app_immersed_flow_driver_detail.hpp`
- 新增 `src/app_immersed_flow_driver.cpp`
- 修改 `src/app_main.cpp`
- 修改 `src/app_dispatch_order_detail.hpp`（RED 证明其硬编码
  `ResolvedCase`，需类型无关地保留既有 root/config/dispatch 顺序）
- 修改 `src/CMakeLists.txt`

测试与治理文件：

- 新增 `tests/unit/test_immersed_flow_dispatch_contract.cpp`
- 新增小型 CLI/driver fixture 或脚本（RED 已证明粗三角形 fixture 会产生
  退化 surface patch，因此使用仓库自有细分 cube 测试生成器）
- 修改 `tests/CMakeLists.txt`
- 本证据包与 task receipt

`src/app_resolved_case_v3_broadcast.cpp` 和 v3 loader 原则上不改；若 RED 暴露
现有 collective/config 缺陷，主 agent 必须先在本文件追加 blocker、影响
分析和窄化白名单。`flow_immersed.cpp`、FVM/IBM 数值实现、科学阈值和 PISO
corrector 数量禁止修改。

## 对外接口

```cpp
int run_immersed_flow_case(
    const CliOptions&,
    runtime::MpiContext&,
    const config::ImmersedFlowCaseConfig&,
    const std::filesystem::path& authoritative_case_root);
```

## 固定构造顺序

```text
load/broadcast immutable schema v3
-> decomposition / topology / geometry / body-fitted boundaries
-> STL surface and SurfaceQuery
-> ImmersedDomain and active layouts
-> GhostStencilPlan / WallQuadraturePlan / LocalFlowPatternTransform
-> field registry / access / FlowState
-> execution / halo / solvers
-> FixedStepImmersedFlow
-> fixed or adaptive trial/retry loop
-> exactly two correctors inside accepted Task 11 attempt
-> committed final flux and four-field force report
-> collective accept or bitwise rollback
```

静态 geometry、surface query、domain 和 plans 必须在时间循环外只构造一次。
输出只能消费 committed report/state；19A 的最小 smoke 不新增持久化格式。

## Mutation-sensitive RED

1. schema v3 被旧 loader 拒绝或误落入 schema v2 driver；
2. v3 validate/print-resolved 未使用 v3 canonical serializer；
3. 19A 接受 IBM absent、WALE present 或非 constant density；
4. geometry/domain/plans 被放进 step/retry 循环；
5. 绕过 `FixedStepImmersedFlow` 或增加第三次 pressure corrector；
6. force 来自 provisional pressure、第二 authority 或只报告单一 signed 值；
7. rank-local failure 未形成 collective failure，或失败后 state 未 bitwise rollback；
8. attempt commit 之前发布输出。

首个 RED 冻结为 CLI dispatch contract：v1/v2 分支身份不变；v3
validate/print-resolved 成功；正常 v3 仅把合法 19A 组合分发给
`run_immersed_flow_case`，非法 presence/density 在进入数值构造前拒绝。

## task fast gate

- `git diff --check`
- dispatch contract RED/GREEN
- v3 broadcast 1/2-rank 既有 focused 测试
- `8^3` 或 `12^3` 单步 constant IBM 产品路径 smoke
- report：`pressure_corrector_count == 2`、四个 force 字段存在、continuity
  和最终 residual 有限
- 小型 rollback/collective failure；只有 collective seam 变化时补 2-rank
- 公共/私有头与 tests-off 受影响构建

不运行 sanitizer、48³、96³、warped/prism 或完整回归。里程碑再执行
1/2/4-rank 和更宽矩阵。

## 回退边界与版权独立性

- 允许回退仅限新 driver、app 分发和直接测试注册；不得回退或重写 Task 11
  accepted 历史和旧证据。
- 不复制、翻译或机械改写 AMReX、IncFlo、OpenFOAM、Basilisk 或私有项目
  源码；19A 只复用 HUNDUN-FLOW 自有 public API 和已接受的产品构造语义。
- 不新增运行时依赖，不访问研究数据，不 push 或发布。

## 接受摘要

- schema-v3 dispatch RED：旧路径以 `/schema_version: unsupported schema
  version 3` 失败；切换到 v3 loader/broadcast/serializer 后 GREEN。
- normal-run RED：在未接入 driver 时明确失败；接入静态构造和 Task 11
  facade 后，12³ 1/2-rank 产品路径均 GREEN。
- 粗 cube fixture 曾以 `surface patch is empty` 失败；只修复测试 fixture，
  产品几何/重构实现未改。
- focused 矩阵 11/11、fresh dispatch 4/4、Release tests-off 均通过。
- requirements、quality、caller-impact 和完整 task diff 由主 agent 合并审查；
  未发现未解决的 correctness blocker。

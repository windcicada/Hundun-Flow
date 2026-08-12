# Task 18A：最小 IBM / Checkpoint v3 诊断冻结证据包

状态：`TASK_GATE_ACCEPTED`

基线：`eafad00f020a616e512e8a3b61fba0cff3580a77`

## 目标与权威边界

18A 只增加 read-only structured diagnostics，不改变 Task 11 数值路径、
Checkpoint v3 事务或 driver 控制流。

`ImmersedFlowDiagnosticSource` 只能由 `FixedStepImmersedFlow` 为最近一次完成
提交的 `FlowState` 创建。source 固定绑定该次提交的 state mutation identity、
step、time、report seal 与 pressure-authority publication generation。active attempt、
失败/rollback、旧 report、被后续提交替代的 source 或 moved-from flow 都必须拒绝。

## 18A 字段

- final continuity normalized L2；
- final momentum normalized L2 三分量；
- maximum / mean wall-normal penetration velocity；
- operator force、budget reaction、surface traction、consistency，各三分量；
- pressure corrector count；
- disposition、failure category、lowest failing rank；
- construction counts：本 rank classified cells、active cells、immersed links、
  pressure-authority donors、wall quadrature points；
- exact snapshots：FP64 reduction calls/scalars/bytes、halo exchanges/bytes/messages、
  execution allocation events/bytes/live/peak。

Checkpoint v3 adapter 报告 operation、phase、presence、disposition、failure、
lowest failing rank、CRC/fingerprint/partition/rollback status 与 manifest CRC。

## 明确延期

surface-query 调用次数、按算子拆分的 matvec、模型求值、每个 attempt 的 delta、
以及精确静态 byte attribution 尚无不改变热路径的公开权威，延期到 18B/Task 20。
18A 不以估算值、测试 seam 或读取 trial cache 冒充这些计数。

## 文件白名单

产品：

- 新增 `include/hundun/diag_immersed_module.hpp`
- 新增 `src/diag_immersed_module.cpp`
- 新增 `include/hundun/diag_checkpoint_v3.hpp`
- 新增 `src/diag_checkpoint_v3.cpp`
- 修改 `include/hundun/flow_immersed.hpp`
- 修改 `src/flow_immersed.cpp`
- 修改 `src/CMakeLists.txt`
- 修改 `src/diag_structured.cpp`（把 SI force unit `N` 加入冻结词表）

测试与治理：

- 新增 `tests/unit/test_immersed_diagnostics_header_contract.cpp`
- 新增 `tests/unit/test_checkpoint_v3_diagnostics_header_contract.cpp`
- 新增 `tests/mpi/test_immersed_diagnostics.cpp`
- 新增 `tests/mpi/test_checkpoint_v3_diagnostics.cpp`
- 修改 `tests/mpi/test_immersed_transaction.cpp`（复用既有 fixture）
- 修改 `tests/CMakeLists.txt`
- 本 packet、receipt、coordinator ledger

## mutation-sensitive RED

1. failed/rollback attempt 覆盖最近一次 published report；
2. source 接受 active attempt、旧 state identity、旧 publication generation；
3. wall penetration 从 predictor/trial 或几何常量读取，而非 committed final flux；
4. force 四字段交换符号或丢失 consistency；
5. pressure corrector count 不等于 report 的 exact 值；
6. collection 修改 source/product counter；
7. checkpoint v3 将 failed report 标成 ok、丢失最低失败 rank 或 rollback status；
8. collective collection 在某 rank source/request 无效时悬挂或部分 submit。

## fast gate

- 两个 standalone public-header contract；
- 1-rank IBM accepted + failed rollback + stale-source；
- 1/2-rank collective record agreement；
- checkpoint v3 completed/failed report adapter；
- `git diff --check`、focused Debug；
- 不运行 sanitizer、48^3 或完整矩阵。

## 版权与回退

只复用 HUNDUN-FLOW 自有 structured diagnostics、Task 11 report/authority 和
Checkpoint v3 report；无外部源码、无新依赖。允许回退范围仅限白名单文件，
且不得改变已经接受的 19A/17A 产品行为。

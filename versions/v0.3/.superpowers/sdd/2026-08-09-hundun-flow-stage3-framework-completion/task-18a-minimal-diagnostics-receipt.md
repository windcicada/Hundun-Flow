# Task 18A：最小 IBM / Checkpoint v3 诊断实施回执

基线：`eafad00f020a616e512e8a3b61fba0cff3580a77`

状态：`TASK_GATE_ACCEPTED`

## 已实现边界

- 新增 `ImmersedFlowDiagnosticSource` 和 local/collective structured adapter；
- source 只在调用 `diagnostic_source(state, report)` 时复制 report、读取 committed
  final flux 并采样 exact counters；未调用时不分配诊断对象、不采样 counter；
- flow 热路径只维护无分配的 attempt generation、state identity 和消费字段 seal；
- 任一新 attempt（包括 physics preflight failure）及 Checkpoint restore 都使旧 source
  失效；report mutation、错误 state、active attempt 和旧 source 均被拒绝；
- summary 覆盖 continuity、三分量 final momentum residual、maximum/mean wall
  penetration，以及 operator/budget/surface/consistency 的 pressure/viscous/total
  三分量力；
- counters 覆盖 classified/active cells、immersed links、donor references、wall
  quadrature points、FP64 reductions、halo exchanges/bytes/messages 与 execution
  allocation snapshots；
- collective adapter 先收敛 request/source readiness，验证 report/force 位级一致，
  再按 sum/max 规则聚合壁面统计和计数，record/sink 失败也 collective 收敛；
- 新增 Checkpoint v3 completed/failed structured adapter，保留 CRC、presence、
  fingerprint、partition、rollback、phase、failure 和最低失败 rank；
- structured diagnostics 的冻结 SI 单位词表增加 `N`，没有改变 record schema 版本。

## RED 与 mutation

- public headers 和 linker registration 缺失；
- corrector count mutation 仍被 diagnostic source 接受；
- 新 attempt / failed rollback 后旧 source 仍可读；
- physics preflight failure 不使旧 source 失效；
- 四字段 force 任一 authority 丢失；
- failed report 被标为 ok；
- collective record 使用各 rank 局部 wall/counter 值或跳过 report agreement；
- Checkpoint v3 completed write/read report 无法生成合法 structured record。

## focused 证据

- `test_checkpoint_v3_diagnostics_header_contract`：PASS；
- `test_checkpoint_v3_diagnostics`：PASS；
- `test_immersed_diagnostics_header_contract`：PASS；
- `test_structured_diagnostics`：PASS；
- `test_immersed_diagnostics_1_rank`：PASS；
- `test_immersed_diagnostics_2_rank`：PASS；
- 真实 `test_checkpoint_v3_1_rank` / `2_rank` adapter 接入：PASS；
- `test_checkpoint_v3_active_attempt_1_rank`：PASS；
- `test_checkpoint_v3_path_agreement_2_rank`：PASS；
- clang/libc++ Release tests-off build 与 `hundun --version`：PASS；
- `git diff --check`：PASS。

未运行 sanitizer、48³或完整矩阵；4-rank、MVP 12³ smoke、linkage 与完整
rollback/collective cluster 按计划进入 MVP gate。

## 延期与调用方

surface-query 调用次数、逐算子 matvec、模型求值、per-attempt delta 与精确静态
byte attribution 没有无侵入权威，延期到 18B/Task 20。19A 现有 console 输出已经
提供用户可见的 corrector/continuity/pressure/四字段 force；structured session
文件落盘与所有合法 driver combination 的 wiring 在 19C 完成。

## 审查与版权

主 agent 完成 requirements、完整 diff、caller、public API/schema、MPI、rollback、
allocation 和 tests-off 审查。因当前 Codex model catalog 明确要求重启前不得设置
子代理 model/reasoning override，本 task 未伪造 Luna review 证据。

实现只复用本项目 Apache-2.0 structured diagnostics、Task 11 report/authority 和
Checkpoint v3 report；未复制外部源码，未访问私有参考源码或研究数据，未新增
运行时依赖。

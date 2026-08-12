# Constant IBM MVP milestone 冻结证据包

状态：`MVP_ACCEPTED_INTERNAL`

候选基线：`24c85d2`

## 目标

证明 19A + 17A + 18A 在同一产品树上形成可运行、可 Restart、可诊断的
constant-density static LFP-GCIBM MVP。不得修改产品数值算法、阈值或 selector
语义。

## 最小矩阵

- 12³ `hundun` driver：1/2/4 ranks，均完成一步、两次 PISO、四字段 force、
  Checkpoint v3 write 与下一步 Restart；
- Checkpoint v3 continuous-vs-restart bitwise：1/2/4 ranks；
- decomposition fast：1x1x1、2x1x1、1x2x1、4x1x1、2x2x1；
- 4-rank immersed transaction：rollback、collective failure、authority rotation、
  corrector count 与 failed-state bitwise neutral；
- 当前树的 diagnostics 1/2-rank证据；4-rank driver/transaction 同时提供最终
  report/authority 的 rank coverage；
- public headers、tests-off、`hundun --version`、`ldd`；
- Task 11 signed-force mutation 与科学 acceptance 复用接受基线证据，因为
  `66080e3..24c85d2` 未修改 force/pressure/operator/final-flux 数值实现。

## 唯一测试图修改

- `tests/CMakeLists.txt`：Checkpoint v3 和 Task 19A driver 增加 4-rank注册；
- `tests/mpi/test_immersed_transaction.cpp`：Checkpoint v3 selector 接受 4 ranks；
- `tests/acceptance/task19a_immersed_flow_dispatch.sh`：检查最后一个 rank payload。

除上述测试/治理文件外不允许修改产品源码。高负载 48³、sanitizer、96³与完整
Task 11 科学矩阵不属于本 gate。

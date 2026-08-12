# Constant IBM MVP milestone 接受回执

产品基线：`24c85d2`

状态：`MVP_ACCEPTED_INTERNAL`

## 结论

19A、17A、18A 已在同一产品 commit 上形成可运行、可 Restart、可诊断的
constant-density static LFP-GCIBM MVP。该结论允许进入 Task 12，但不触发 product
仓库同步；product 仍等到 Task 21。

## 当前树证据

- `test_task19a_immersed_flow_dispatch_1_rank`：PASS；
- `test_task19a_immersed_flow_dispatch_2_rank`：PASS；
- `test_task19a_immersed_flow_dispatch_4_rank`：PASS，12³ `2x2x1`；
- `test_checkpoint_v3_4_rank`：PASS，continuous-vs-restart bitwise；
- `test_laminar_ibm_decomposition_4_rank_2x2x1_fast`：Clang/libc++ Release PASS；
- `test_immersed_transaction_4_rank`：Clang/libc++ Release PASS；
- Task 18A diagnostics 1/2-rank focused gate：PASS；
- Checkpoint v3 1/2-rank、rollback、path agreement：复用 17A/18A 当前产品
  authority 的已完成证据；
- Task 11 decomposition 五配置与 signed-force mutation：复用接受基线
  `66080e324089599711fdb26082af9b330bfdb5ce`，19A/17A/18A 未修改其数值
  authority、阈值或 selector；
- `test_checkpoint_v3_header_contract`、Task 18A 两个 diagnostics headers、
  `test_stage3_flow_header_contract`、`test_wall_force_header_contract`：PASS；
- clang/libc++ Release tests-off build、`hundun --version`、`ldd`：PASS；
- 无遗留 HUNDUN 或 MPI 进程。

## 调度修订记录

最初尝试在 Debug 下机械重跑五个 decomposition fast selector；首个 16³ 单 rank
运行约 4 分 40 秒仍未完成，诊断价值低。核对 PID、命令和工作目录后正常停止，
确认无残留进程。随后按紧凑计划改为：复用 Task 11 五配置证据，并在当前产品
commit 上运行 Release 4-rank `2x2x1` 单一确认。

4-rank driver 首次测试配置误用 `4x1x1`，导致每 rank x 方向 3 cells 小于冻结
halo width 4。未放宽产品约束；测试改用科学上合法且已经覆盖的 `2x2x1` 后通过。

## 审查与边界

主 agent 审查确认本 milestone 只修改测试注册、测试 process grid 和治理文档，
产品源码零修改。没有运行 48³、96³、sanitizer 或完整 Task 11 科学矩阵；没有
访问私有源码、研究数据、product 投影仓库，也没有 push 或发布。

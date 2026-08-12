# S3-E1 exact counters / performance artifact 验收回执

状态：`ACCEPTED`

accepted at：`2026-08-09T22:55:08+08:00`

accepted parent：`4e9a29e34a28eec1ecebfb15bc651a652ceb5b9a`

accepted commit：`d1f866ca38fcf0af93f1c61c1c3a575d2e78432f`

## Schema v2 与冻结 counter inventory

`diag_performance_artifact` 仍是唯一 canonical JSON authority。schema v2 只追加
`stage3_identity` 和 `exact_counters.algorithmic_work`；schema v1 要求这些字段为空且输出
顺序、字段名和编码路径不变。v2 强制记录 source commit/clean/tree/binary、compiler/MPI、
build type、rank/process grid、cpuset/thread budget、case/profile/geometry fingerprint、warmup/
measured/repetition，并精确接受冻结的 17 个 stable IDs，拒绝 missing/unknown IDs。

17 个 counter 分别由 SurfaceQuery、ImmersedDomain、GhostStencilPlan、WallQuadraturePlan、
FixedStepImmersedFlow、WaleModel 和 Checkpoint v3 codec 的实际 owner 累积。所有加法与 snapshot
delta 都检查 overflow/decrease；static construction snapshot 在 step interval 保持为零增量；
failed attempt 的实际 ghost/WALE work 保留，同时 `FlowState` step/time 回滚。diagnostic snapshot
本身不增加 business counters，也不参与数值分支。

## 8-cubed RED/GREEN 与 formal selector freeze

新增真实 constant + static IBM + WALE 8-cubed fixture。mutation 在 commit 前注入一次失败，证明
rollback state 与 retained work 同时成立；accepted interval 精确验证两次 PISO wall constraints、
wall quadrature、三次 force reductions、owned-active WALE gradients 和每 attempt 一次 WALE
evaluation。1/2-rank fast 均通过。

formal selector 已固定为 global 24-cubed、2 warmup + 3 measured steps + 1 repetition、ranks
1/2/4、timeout 1800、独占 `hundun_stage3_performance_m`。formal 路径禁止 failure injection，
以 `MPI_Wtime` 收集每 rank raw sample，只由 rank 0 在 runner 提供的
`HUNDUN_STAGE3_EVIDENCE_DIR` 原子写出一个 canonical schema-v2 artifact。E1 仅用 `ctest -N`
核对三项注册，没有执行 formal selector。

## Stage 1 schema-v1 performance reconciliation

复跑 Task 25 发现 Stage 3-capable flow scratch 已使旧 schema-v1 exact resource snapshot 增加
40 bytes/owned-cell transient peak、8 bytes/rank control block，以及每 measured step 32 bytes
FP64 logical payload。该差异来自已集成 Stage 3 flow path，不是 E1 serializer 字节变化。
独立 C++/MPI oracle 的公式和 shell snapshot 同步到 340 bytes/cell、24 bytes/rank、2360
logical bytes/measured-step，并修复 out-of-tree CMake 中 oracle/evidence helper 的 sibling
定位；1/2/4-rank 三项完整 Task 25 acceptance 全部通过。未修改 solver work、Halo、I/O、
matvec、preconditioner 或数值阈值。

## Verification

- focused build：全部 E1 unit/header/MPI targets 与 `hundun` PASS；
- unit/header/self gate：12/12 PASS；
- Stage 3 performance fast：1/2-rank 2/2 PASS，real 5.85 s；
- Task 25 schema-v1 performance：1/2/4-rank 3/3 PASS，real 154.95 s；
- formal inventory：1/2/4-rank 三项精确列出，未运行；
- Clang 15/libc++ Release tests-off `hundun` build 与 `--version`：PASS；
- tests-off `nm -C` 不含 `TestAccess`、`ENABLE_TEST_ACCESS` 或 test counter symbol：PASS；
- tests-off `ldd` 仅系统 C/C++、MPI、pthread/dl/hwloc/numa 等依赖：PASS；
- `git diff --check`：PASS。

内容 SHA-256：

- `diag_stage3_performance.hpp`：
  `058ad6b18fb353086beca882da175c8a804c9ae076b0258da6bbad83fcda4683`；
- `diag_stage3_performance.cpp`：
  `51f8a3a9b2065f29d64483e7e34ed1f8ed9412f15a7385824f5af30d6961e844`；
- `test_stage3_performance.cpp`：
  `96fb0902a910fa3d4138dec38c972782b16f86231a11edcb578ddadfb5ae3b9b`；
- `stage3_performance_evidence.cpp`：
  `2570b1de85fb1ade1e0f440a7977511205e359f8dc583df7c75f307aeb2b3d1b`；
- debug `test_stage3_performance`：
  `c21dea65966bbb68473d00ac40b780375109c05a96cf9fc4312d75ce443d3308`；
- tests-off `hundun`：
  `135dfd87ba8c4c0b658bdd71dc0987da084ab42c59c78e18af264a2509ece5b7`。

未运行 24/48/96-cubed 或任何正式矩阵；未访问私有源码、研究数据或研究进程，未
push/publish。Task 11 科学 authority、阈值、两次 PISO、force sign、rollback、Restart
与 MPI 一致性要求均未修改。

## 权威与 DCO

激活 design/reference/plan SHA-256 仍分别为
`4cc2052cf99874cefaffe1d11fc954195269fac77ea06fee06a92582ab5ea64e`、
`0349d4a6301b2213db98ba680557baab808bde62b40f09d99cf392af1e4058b8`、
`bed28804133174481d9ad41266a8d644a0ad056ff9be7ab48f44878b11c45f44`。

提交 subject：`feat: add Stage 3 exact performance evidence`

DCO：`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`

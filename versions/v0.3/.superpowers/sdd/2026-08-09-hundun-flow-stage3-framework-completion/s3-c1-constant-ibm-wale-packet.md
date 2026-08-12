# S3-C1 constant IBM+WALE 冻结执行包

状态：`TASK_GATE_ACCEPTED`

accepted parent：`5cc2f45fd3b6aa1adbc81eb8f5cd3e177fd68513`

activation commit：`2a6268d28666a65c176af0028e6c370ba12df85e`

## 目标

开放 `LFP-GCIBM/wale/constant` 产品路径。每次 attempt 只构造一个 move-only WALE
authority：BE 使用 `u^n`，BDF2 使用
`u^n + dt/dt_previous * (u^n - u^(n-1))`；regular gradient 使用 background
离散，interface gradient 使用既有 immersed reconstruction，inactive workspace 以 NaN
sentinel 证明不会读取 solid physical slot。

同一 authority 的 active-cell/face `mu_eff` 必须服务 momentum predictor、两次 PISO
corrector 的 diagonal/traction、final momentum residual 和 wall force。Task 11 的 pressure、
force sign、两次 PISO、rollback 与 MPI authority 不变。

## 允许文件

- `include/hundun/flow_immersed.hpp`
- `src/flow_immersed.cpp`
- `src/flow_immersed_wale_detail.hpp`
- `src/flow_immersed_access_detail.hpp`
- `src/app_immersed_flow_driver.cpp`
- `tests/support/flow_immersed_test_access.hpp`
- `tests/mpi/test_immersed_wale_constant.cpp`
- `tests/acceptance/stage3_flow_models_fast.sh`
- `tests/cmake/stage3_science_registration.cmake`
- 本 packet、对应 receipt 和 ledger status row

## 禁止边界

- 不修改 material/ideal-gas、Checkpoint、diagnostics provider 或 Task 11 MMS selector；
- 不运行 24/48/96³；正式 24/48³ 仍由 S3-V1 独占；
- IBM+WALE Checkpoint 读写仍由 S3-R1/R2 完成；C1 fast case 不触发 Checkpoint；
- 不访问私有源码/研究数据，不干扰研究进程，不 push/publish。

## RED / mutation contract

- accessor seams 先可编译并返回 unavailable；行为 RED 必须由 facade 拒绝
  `domain && wale` 产生；
- injected failure 后 history/committed/metadata bitwise rollback，且 WALE identity、count、
  wall fingerprint 不发布；
- repeated preparation 在 evaluate gate 被拒绝；
- wall force 若改用 molecular `mu`，实际 wall fingerprint 与 authority 不一致；
- interface inactive workspace 是 NaN，读取 solid slot 会失败；
- constant base、四组 force 或六字段 WALE summary 任一 mutation 都使 diagnostic source
  authentication 失败。

## 接受检查

- [x] executable behavioral RED 以 combined unsupported 失败；
- [x] BE/BDF2、rollback、恰好两次 PISO、单次 evaluate、force 和 WALE summary 通过；
- [x] repeated-evaluate、molecular-force、WALE-seal 三项临时 mutation 被执行测试杀死；
- [x] direct 8³、fast 12³ 1/2-rank、driver fast 1/2-rank 通过；
- [x] Task19A 1/2-rank 和 Stage 3 flow header contract 通过；
- [x] unknown/extra selector 返回 2；
- [x] 完整 diff、caller、allocation workspace、MPI ordering 和 Task 11 authority 由主 agent
  审查；
- [x] code graph sync、`git diff --check`、DCO、clean tree 与后台进程在提交前复核。

提交 subject：`feat: couple WALE to immersed flow`

签署 identity：`WANG YUDONG <wangyudong@buaa.edu.cn>`

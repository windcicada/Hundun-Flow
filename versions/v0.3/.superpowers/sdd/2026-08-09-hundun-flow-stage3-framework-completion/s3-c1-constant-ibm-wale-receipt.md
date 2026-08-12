# S3-C1 constant IBM+WALE 验收回执

状态：`ACCEPTED`

accepted parent：`5cc2f45fd3b6aa1adbc81eb8f5cd3e177fd68513`

accepted at：`2026-08-09T13:20:13+08:00`

## 实现边界

本 task 只修改 packet 白名单中的 9 个实现/测试文件，并增加本 packet、receipt 与 ledger
状态更新。未修改 material/ideal-gas、Checkpoint、diagnostics provider 或 Task 11 MMS
selector。

`ImmersedWaleAttemptAuthority` 持有唯一一次 `WaleAttemptCoefficients`、active-cell/face
`mu_eff`、summary 与 wall fingerprint。regular face 和 immersed wall 的 predictor/operator、
final residual、budget reaction 与 surface traction 消费同一冻结 authority。uniform molecular
路径仍使用构造期预分配 workspace，并保留 `mu == 0` 的无 viscosity exchange 快路。

driver 只开放不触发 Checkpoint 的 constant IBM+WALE 路径；本次运行会实际读写 WALE
Checkpoint 的组合仍显式拒绝，等待 S3-R1。material/ideal-gas 组合仍未开放。

## TDD 与 mutation 证据

行为 RED：

- seams/header 先成功编译；
- `test_immersed_wale_constant_1_rank` exit 8；
- 输出 `COMBINED_IBM_WALE_UNSUPPORTED reason=1`，失败点是 facade 对
  `domain != nullptr && wale != nullptr` 的现有拒绝，而非 compile/link/accessor。

临时 mutation（均由 `apply_patch` 施加并立即恢复）：

- 第二 corrector 前再次准备 WALE：exit 8，`correctors=1`，单次 evaluate gate 杀死；
- wall force 强制 molecular `mu`：exit 8，`correctors=2`，reason 9
  (`final_conservation_defect`)，wall fingerprint 杀死；
- diagnostic seal 忽略 WALE fields：exit 8，测试在 `rejected` 断言处杀死。

恢复后的 source 中 `wale->evaluate(` 恰好一个产品 caller，
`prepare_wale_authority(` 只有定义和单次 attempt 调用；没有 mutation 注入残留。

## GREEN / task gate

构建：

```text
cmake --build build/debug -j32 --target \
  test_immersed_wale_constant hundun test_stage3_flow_header_contract
```

exit 0。header gate `1/1 PASS`。

v2 原样 task-gate regex：`5/5 PASS`，0 failures，real 53.28 s：

- direct 8³ 1-rank；
- `test_stage3_flow_models_fast_[12]_rank`；
- `test_task19a_immersed_flow_dispatch_[12]_rank`。

注册的 executable 12³ fast 后缀行另行串行执行：`2/2 PASS`，0 failures，real
79.25 s；1-rank 51.28 s，2-rank 27.96 s。所有 MPI 行共享
`RESOURCE_LOCK hundun_stage3_mpi_m`，未并发多个 M job。

selector 合同：unknown selector exit 2；`fast` 加 extra argv exit 2。

提交前 fresh verification 重新 configure/build，并运行 layout、include authority、P0
registration contract/mutation、8³ direct 与 flow header：`6/6 PASS`，0 failures。

未运行 24/48/96³、Release、ASan、UBSan 或任何正式/长 selector。

## 证据身份

- implementation/test diff SHA-256（治理文件加入前）：
  `a562c6a62108f381738c01e49e94155899ea9281ed086220f290917fbfb32325`；
- Debug `test_immersed_wale_constant` SHA-256：
  `cfe3f28de69c7cd86e29bf283e3832ef5906b16cad28942e7f96e876eedf0570`；
- Debug `hundun` SHA-256：
  `39d94f531ccd11296a4ca6f4ab0bff9beb5bb71ba01059b359c812b1024d6af1`；
- Debug `test_stage3_flow_header_contract` SHA-256：
  `bdbfecbda0a821258d92e8e23d1f3c369d36ab0f81538ab995e9678be1a5db41`；
- final focused `LastTest.log` SHA-256：
  `00f5eebdb056dfa5206eb9870ab90b995f5f62217fe8c71cadec678710ae8118`；
- CMake `3.31.12`；codegraphf `0.9.8`，sync exit 0。

激活候选 design/reference/plan SHA-256 仍分别为
`4cc2052cf99874cefaffe1d11fc954195269fac77ea06fee06a92582ab5ea64e`、
`0349d4a6301b2213db98ba680557baab808bde62b40f09d99cf392af1e4058b8`、
`bed28804133174481d9ad41266a8d644a0ad056ff9be7ab48f44878b11c45f44`。

## 主 agent review

- Task 11 pressure/force authority、force sign、两次 PISO、rollback、Restart 与 MPI
  一致性要求未修改；
- inactive lagged workspace 的 NaN sentinel 在 BE/BDF2 direct case 中通过，证明 interface
  gradient 未读 solid physical slot；
- variable-viscosity diagonal 使用同一 face authority，wall diagonal/traction 使用同一
  active-cell authority；
- uniform molecular operator 的 staging buffers 在构造期分配，没有新增 steady-attempt
  临时 vector 分配；
- C1 没有委派 worker，完整 diff、科学判断、集成、测试、DCO 与提交均由主 agent 负责；
- 无私有路径、研究数据、研究进程、push 或 publish 操作。

提交 subject：`feat: couple WALE to immersed flow`

DCO：`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`

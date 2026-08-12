# Task 13 + 19B：body-fitted WALE flow 实施回执

基线：`43d8446a4de354c40749b94d1f2143a8d744a4e9`

状态：`TASK_GATE_ACCEPTED`

## 已实现边界

- `FixedStepConstantDensityFlow::attempt_with_wale` 在 attempt 开始后构造 BE/BDF2
  lagged velocity，只调用一次 `WaleModel::evaluate`；
- `mu_sgs` 从 cell active-order 映射到共享 face authority，形成冻结 `mu_eff`；
- predictor、accepted/history momentum residual、最终 momentum residual 和 physical
  boundary contribution 共用该 `mu_eff`；
- enthalpy/scalar 共用同一 `mu_sgs`，分别除以冻结的 `Pr_t`/`Sc_t`；
- molecular-only API 继续走原 scalar-viscosity 实现，constant face-viscosity RED 与
  旧路径 bitwise 等价；
- failed attempt 不发布 WALE summary，state/history/metadata bitwise rollback；retry
  使用相同输入得到相同 identity；
- body-fitted `none/wale` 已接入同一 `hundun`；输出包含 WALE identity、min/max 和
  exact-zero count；IBM+WALE 仍显式拒绝；
- driver 不为 WALE-only 写 Checkpoint v3，WALE presence continuation 明确延期至
  Task 17B。

## RED 暴露并修复的根因

最初 12/24 screen 的 molecular-only control 通过，而 WALE 路径出现 periodic
conservation failure。根因是一个 periodic face 的两个本地表示各自按不同 image
distance 插值 viscosity。修复后 periodic pair 采用对称 canonical 值；同 rank 中
较大 global-face ID bitwise 复制较小 ID 的 authority，跨 rank 则按同一 canonical
P/N 顺序求值。没有修改守恒阈值、PISO corrector 数量、滤波或阻尼。

## 最终 focused 证据

最终产品/测试源码树上的验证：

- Debug affected suite：11/11 PASS；包括 variable-viscosity FVM 1/2-rank、
  body-fitted WALE 1/2-rank、12/24 screen、driver 1/2-rank、PISO、WALE unit 和
  header contract；
- standalone/public header contracts：4/4 PASS；
- Release body-fitted WALE 1/2-rank + 12/24 screen：4/4 PASS；
- ASan focused `test_wale` + 1-rank body flow：2/2 PASS，
  `ASAN_OPTIONS=detect_leaks=0`；
- UBSan focused `test_wale` + 1-rank body flow：2/2 PASS；
- Release tests-off `hundun` build：PASS；
- `git diff --check`：PASS。

本节点实现与测试 diff SHA-256：
`499c2973fbff2f7f1385a3133fa02daaec01a94c417435f715b1cd721d1e1dd4`。

二进制 SHA-256：

- Debug `test_wale_body_fitted`：
  `3768fa378d557cfe78105d883ba35ff20d706be58d7567d5b7f5bd2c3637a286`；
- Release `test_wale_body_fitted`：
  `fa2793c0113f208b6833b3dc3b6f529c1ab1f16bee352fc7182a5f8c54c7f9e4`；
- ASan `test_wale_body_fitted`：
  `3f0888d379fd43ab35355413bce7c2b81a2ccac4cc28985ecac799e096f02a99`；
- UBSan `test_wale_body_fitted`：
  `e6c11074608655b21b748565b721c307c56c9c16413ae9a35d30beb3c77a4737`；
- tests-off `hundun`：
  `aefc6a33afc5ea50987b98b39e9be579681b5cbd31b33a07abeb109ef752c20c`；
- Debug `hundun`：
  `76b4f2a248ea56d4aea4a282ba435321d6fe995868c00eef68c232049a0de4c7`。

CTest 日志 SHA-256：

- Debug affected suite 完成时的 `LastTest.log`：
  `9f56278287a1fcd6d068b3330a713e6e7453018195d715fe8cc881c96f5ee6af`；
- 后续 Debug header-only run：
  `b753b223b44d0bfef6ae70932f7bd3fd340233f4621573d57fea254dcea41137`；
- Release：
  `97658532f89072058de2137ac4661d0a46144361c3328db309b1e07660680d2e`；
- ASan：
  `eaa61d567c3c785bd24a3a3c19583ce16fd7d7b438bc697060b735ec43358178`；
- UBSan：
  `d985eac287a76c23a126bc6c138e7f1a1f78c520001f6deb40f44dd46162d2c0`。

所有 focused build 使用 Clang 15、libc++、OpenMPI 和
`-DOMPI_SKIP_MPICXX`；Debug/Release/ASan/UBSan/tests-off 分属独立构建目录。
曾有一次同一 Debug archive 的并行 build 调度冲突，顺序重跑后通过，判定为 runner
竞争而非产品失败。没有 Task 13 后台进程遗留。

## 主 agent 审查

主 agent 完成数学/量纲、一次求值、多消费者 authority、PISO 顺序、public API/ABI、
ownership/lifetime、异常/collective rollback、allocation、periodic MPI、driver、测试
注册、完整 diff 和调用方审查。公开 API 增加 overload 和 report 字段，但未修改
schema、Restart 编码、diagnostics identity 或已有 field identity；项目仍处于 0.2.0
前的 Stage 3 开发期。

当前运行时提示 model catalog 变化，禁止重启前设置子代理 model/reasoning override；
因此本节点没有派生或伪造 Luna review 证据。审查由主 agent 完成。

## 延期与能力声明

- IBM+WALE：Task 16；
- WALE Checkpoint v3 presence：Task 17B；
- material/ideal-gas WALE composition：Tasks 14--16；
- 48-cubed、24-cubed 1/2/4-rank 正式科学矩阵：Task 21 冻结候选；
- 极端规模性能与非阻断整洁工作：Stage 3 末尾或后续维护。

本回执只接受 Task 13 + 19B 和 WALE milestone 的紧凑门，不把 12/24 screen 声称为
完整三层科学收敛证明。

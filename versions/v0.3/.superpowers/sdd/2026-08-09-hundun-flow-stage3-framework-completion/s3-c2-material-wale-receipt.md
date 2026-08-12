# S3-C2 material body-fitted/IBM+WALE 验收回执

状态：`ACCEPTED`

accepted parent：`3bdb3cf9ba3fe3c582c85af5ad40b9f759297641`

accepted at：`2026-08-09T15:58:21+08:00`

## 实现边界

本 task 实现 schema-v3 `none/wale/material` 与
`LFP-GCIBM/wale/material` 两条组合路径。没有新增 public variable-density
attempt；`ImmersedFlowStepAttemptReport` 继续复用 material base variant、optional WALE
与 optional force。冻结的 private
`detail::DensityClosureBridge::attempt_with_optional_wale` 是 body-fitted facade 进入既有
`FixedStepMaterialDensityFlow::attempt_common` 的唯一组合入口；public material attempt 和
ideal-gas bridge 均委托 null-WALE 分支。

Task 13 regular-cell lagged gradient、WALE evaluation、cell/face `mu_sgs` 与 `mu_eff`
已原样抽取到 `flow_body_fitted_wale_detail.hpp`。helper 显式接收 `rho_attempt`；constant
caller 传 committed density，material caller 只在 predictor transport、optional predictor
closure 与 trial-density halo 完成后传 trial density。helper 不拥有 transport、PISO、final
flux 或 commit。

body-fitted material 固定顺序保持 predictor transport → 一次 WALE → material momentum →
PISO #1 → provisional transport → PISO #2 → final corrected-flux transport/assessment →
commit。IBM material 在同一阶段显式把 trial density view 交给既有 C1 immersed WALE
authority；后续 corrector 不刷新系数。两条路径均只在成功 commit 后发布 WALE summary；
recoverable/non-retryable failure 保持 summary/force 未发布并 bitwise rollback。

shared helper 使用 attempt-local 一次性 evaluation token。material core 仅在 commit 后记录
真实 count，body-fitted facade 的既有 test access 读取该 private bridge 值；第二次调用在
运行期被拒绝。没有改动 PISO 阈值、exactly-two correctors、force sign、final-flux
finalizer、Restart 或 MPI 一致性 authority。

## TDD、RED 与 mutation 证据

初始可执行 RED：`test_wale_body_fitted_1_rank` 在 helper refactor 前后均 PASS；新增
`test_material_wale_composition_1_rank` 编译链接成功后 exit 8，并同时输出
`BODY_MATERIAL_WALE_UNSUPPORTED reason=1` 与
`IBM_MATERIAL_WALE_UNSUPPORTED reason=1`。失败来自 D1 组合能力未放行，不是
compile/link 或测试夹具错误。

composition test 对 `domain=nullptr` 与 closed sphere 各自构造 committed density/velocity
相同但 predictor face flux 不同的成对 case；正确实现的 WALE identity 必须随 predicted
`rho_attempt` 改变。另以相同 state、不同 molecular `mu` 构造最终 velocity 确实不同但
WALE identity 必须相同的成对 case，从而拒绝 corrector/final-velocity refresh。两条路径
都覆盖 first-attempt failure、history/committed/trial/metadata bitwise rollback、half-dt
retry identity 重算和只发布成功 attempt；IBM 另断言 wall `mu_eff` fingerprint、force 与
真实一次 evaluation count。

主 agent 用 `apply_patch` 施加并立即恢复以下可编译 mutation：

- body-fitted material 强制使用 committed density：CTest exit 8，2.73 s，死于
  collective-any identity-change 为零；
- IBM material 强制使用 committed density：CTest exit 8，19.40 s，死于同一独立
  identity authority；
- corrector/final transport 后第二次 evaluate：首次审查发现旧测试会漏过；加入 shared
  attempt token 与真实 count bridge 后重跑，CTest exit 8，37.88 s，body 路径以
  `invalid_input` 拒绝第二次调用；
- 使用独立 token 恶意绕过 count，并在 commit 前用 final trial velocity/density 刷新
  summary：CTest exit 8，2.91 s，精确死于 molecular-`mu` 变体的 WALE identity 不再
  相等。

最终源码无 mutation/debug 残留，恢复后的 composition 1/2-rank 与 characterization
重新通过。

## GREEN / task gate

冻结的 Debug targets 全部构建成功：

```text
test_wale_body_fitted
test_material_wale_composition
test_material_density_piso
test_material_density_piso_header_contract
test_adaptive_time_control
hundun
```

header/adaptive gate：`2/2 PASS`，0 failures，real 0.03 s。

v2 原样 focused gate：`7/7 PASS`，0 failures，real 77.59 s：

- material composition 1/2-rank：36.53 / 21.42 s；
- schema-v3 flow-model fast acceptance 1/2-rank：9.87 / 6.08 s；
- accepted body-fitted WALE characterization 1/2-rank：0.64 / 0.50 s；
- existing material PISO 1-rank：2.53 s。

fast acceptance 现在逐项执行 constant IBM+WALE、material body-fitted+WALE、material
IBM+WALE；每项要求 exactly two correctors 与 authenticated WALE identity，IBM 两项继续
要求完整 force，body-fitted 禁止伪造 force。material/WALE Restart 限制保持拒绝。

未运行 1/2/4 full matrix、正式 24/48³、96³、Release、ASan、UBSan 或 complete
sanitizer。

## 调用图与单一 authority

CodeGraphF `0.9.8` final sync exit 0、index up to date：

- `FixedStepMaterialDensityFlow::attempt_common` 唯一调用既有
  `PisoCoupler::correct_material_density`；源码精确匹配为 PISO #1/#2 两处；
- body-fitted facade 精确调用冻结的 `attempt_with_optional_wale`，该 bridge 唯一委托
  `attempt_common`；CodeGraphF 对这个新 private caller 漏报，已用 exact `rg` 交叉验证；
- IBM `attempt_immersed` 精确保留两处 `correct_active_pressure`、一处
  `prepare_wale_authority` 与一处 `finalize_from_corrector_two_flux`；
- material core 精确保留 predictor/provisional 两处 `stage_trial`、两处
  `correct_material_density` 和唯一 final `finalize_trial`/closure assessment 分支；同一
  profile 没有第二套 pressure 或 final-flux authority。

## 证据身份

- implementation/test staged diff SHA-256（治理文件加入前）：
  `882a73cff10c055c1cbf7fa240b068fe106a110d81ce46d23894b70f130464ef`；
- Debug `test_wale_body_fitted` SHA-256：
  `64f4d466a233b8a6f94d3196649092d2340bcec82796e95b382d29efd0a8da01`；
- Debug `test_material_wale_composition` SHA-256：
  `49bbdea8bc6263478384854bdd9b380d59ab0202e9d904eb83a1a00813c246cb`；
- Debug `test_material_density_piso` SHA-256：
  `85d060bca565ebf21123c29bc509ed14ee895d0f87cd1e6e5fc7723b01a35c8f`；
- Debug `test_material_density_piso_header_contract` SHA-256：
  `6af712d961fc7ec745a56d27ed1d50cad8898db600b0385ce7a8dfddf1929cd1`；
- Debug `test_adaptive_time_control` SHA-256：
  `04da33def978e0d3efb7b7ae23107dc37b698c8046f0d248fe6ba81ec18457e4`；
- Debug `hundun` SHA-256：
  `12eab988c81d2a435db1fb849965d06bf59ed37747b40cc841a6a058abb4e1a3`；
- final focused `LastTest.log` SHA-256：
  `0d23b8a00532b052a2b1ff799ea86113c70c98d3d57e680452bb7bab5a558c62`；
- CMake `3.31.12`；CodeGraphF `0.9.8`。

激活候选 design/reference/plan SHA-256 仍分别为
`4cc2052cf99874cefaffe1d11fc954195269fac77ea06fee06a92582ab5ea64e`、
`0349d4a6301b2213db98ba680557baab808bde62b40f09d99cf392af1e4058b8`、
`bed28804133174481d9ad41266a8d644a0ad056ff9be7ab48f44878b11c45f44`。

## 主 agent review

- public header contract 无新增接口；private bridge/model-summary 必须同时 null 或同时
  nonnull；public material 与 ideal-gas 路径的 null-WALE 行为重新通过；
- helper 的 constant caller 继续传 committed density，原 field/report fingerprints 位级
  PASS；material caller 明确位于 predictor trial-density halo 后、任何 momentum/PISO 前；
- IBM WALE authority 继续使用 C1 weighted-least-squares/active interpolation，只替换为
  显式 trial `rho_attempt` view；同一 authority 贯穿 momentum、wall viscous force 与
  final checks；
- exactly two PISO、force sign、rollback、Restart、MPI 失败 rank 与 Task 11 科学阈值均未
  修改；
- C2 未委派；R1 worker 结果仍隔离、未集成，O1 未启动；
- 无私有源码、研究数据、研究进程、push、publish 或 Stage 4–6 实现。

提交 subject：`feat: combine material density IBM and WALE`

DCO：`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`

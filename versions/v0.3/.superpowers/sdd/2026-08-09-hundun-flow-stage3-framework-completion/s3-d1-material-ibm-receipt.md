# S3-D1 material IBM vertical slice 验收回执

状态：`ACCEPTED`

accepted parent：`3735742b65d8372d54d6ac517eac932281f37fe2`

accepted at：`2026-08-09T15:08:50+08:00`

## 实现边界

本 task 实现 `LFP-GCIBM/none/material` vertical slice。public seam 只增加冻结的
`ImmersedFlowDensitySetup` 与 `FixedStepImmersedFlow::create` overload；旧 overload
继续构造 canonical constant setup。material 路径复用既有
`MaterialDensityTransport`、`MaterialDensityStepAttemptReport` 和
`detail::DensityClosureBridge`，没有新增 immersed-only material report、第二套 transport
或第三次 pressure corrector。

material adapter 的固定顺序为 predictor transport → PISO #1 → provisional transport →
PISO #2 → final corrected-flux transport。每次 corrector 使用同一 Ghost row 生成的
`rho_wall`、非零法向密度梯度和正 `D_wall` authority；final transport 只消费第二次
corrector 后的共享 `FaceMassFlux`。active density/transport/conservation 只遍历 fluid
cells，fluid--solid、solid--solid 与 inactive physical face 保持 canonical positive zero，
inactive physical boundary contribution 不进入 active conservation。

除 packet 文件外，修改了既有 private test-access 数据结构
`src/flow_immersed_access_detail.hpp`，只为从 facade 产品路径观测同一次 wall-density
authority，使 clipping/epsilon/fallback 与 homogeneous-Neumann mutation 可执行。没有增加
public API、diagnostics provider 或生产期 test seam。

## TDD、RED 与 mutation 证据

初始 behavior RED：注册与 header seam 可编译；material density/transaction 1-rank 两行
均 exit 8，并输出 `IMMERSED_MATERIAL_UNSUPPORTED`，失败点是 facade 尚未支持 material，
不是 compile/link failure。

主 agent 完整 diff review 另发现并关闭以下假阳性/漏洞：

- mutation 测试原先在 stale state 上调用 diagnostics，且 material base 被 seal 无条件拒绝；
  先移动测试并要求原始 report 可用，得到 exit 8、
  `immersed-flow diagnostic source is stale`，再扩展同一个
  `diagnostic_report_seal`。原 report GREEN，26 种逐字段 nested/outer mutation 全部被
  authenticated seal 与 diagnostic source 拒绝；
- 临时 `rho_wall += 1e-3` mutation exit 8，精确死于
  `maximum_product_wall_density_error <= 1e-12`；
- 临时 homogeneous-Neumann mutation（产品 authority 法向导数置零）exit 8，精确死于
  `maximum_product_wall_errors[1] <= 1e-10`；
- 两个 mutation 均由 `apply_patch` 施加并立即恢复，最终源码无 mutation/debug 残留。

transaction tests 覆盖 negative `rho_wall`、NaN `rho_wall`、NaN active donor、
rank-local corrector failure、late final-transport failure、before-commit failure。它们逐项
断言最低失败 rank 和 `non_positive_density`/`non_finite_state` 原因，并对
history/committed/trial/metadata 做 bitwise rollback 检查。无可靠 rank 的 MPI/未知异常
沿用 body-fitted authority：rollback 后抛出，不伪造 rank 或生成不可认证 report。

## GREEN / task gate

冻结的 8 个 Debug target 全部构建成功：

```text
test_immersed_material_density
test_immersed_material_transaction
test_material_density_transport
test_material_density_piso
test_material_density_transport_header_contract
test_material_density_piso_header_contract
test_stage3_flow_header_contract
hundun
```

fresh registration/layout/header gate：`9/9 PASS`，0 failures，real 2.20 s。它包括
source-layout、CMake include authority、P0 registration contract 及其 mutation fixtures，
以及三个 header contracts。

v2 原样 focused gate：`8/8 PASS`，0 failures，real 75.40 s：

- material IBM density 1/2-rank：8.69 / 5.48 s；
- material IBM transaction 1/2-rank：27.35 / 15.99 s；
- flow-model fast acceptance 1/2-rank：9.22 / 5.52 s；
- existing material transport/PISO 1-rank：0.64 / 2.48 s。

C1 constant preservation 另行串行验证：direct + 1/2-rank fast rows `3/3 PASS`，
0 failures，real 92.11 s。constant pressure coefficient path 不再分配临时 face-density
vector，也不增加 material-only conservation reduction。

未运行 1/2/4 full matrix、正式 24/48³、96³、Release、ASan、UBSan 或 complete
sanitizer。

## 证据身份

- implementation/test staged diff SHA-256（治理文件加入前）：
  `39e33f76cc59479e18ed36727cfcb0236d638266a89eea85957d8ba42531dd21`；
- Debug `test_immersed_material_density` SHA-256：
  `563c445c4eebde4f3e346d5697245a7e7002c0f0c8f8d06ea6c00391af64480e`；
- Debug `test_immersed_material_transaction` SHA-256：
  `801974ef7a439b49e1012368233fc2ad9669f60bc7bec6c8d8879506dfe50f28`；
- Debug `test_material_density_transport` SHA-256：
  `201c68359fdeadacb02bba5679a948244c4a53fc7a76b3a4885f701eb040007b`；
- Debug `test_material_density_piso` SHA-256：
  `45928f3ab00dc321af77f6b1848e618eb67bec515ec8dbd1f410117cf4fc69e1`；
- Debug `hundun` SHA-256：
  `1d9c9c81394ef3a2ec45df35cc3856f276c89158d7b1426bb9b5b8d9985e17c2`；
- Debug `test_stage3_flow_header_contract` SHA-256：
  `ffc548b06a56d1f071903a128463fc0c31acddac30a0f2e0a3df1af4a1f40c7a`；
- final focused `LastTest.log` SHA-256：
  `fbc61513e7a78461a8c9d98020fc088c35569144b871e2689553a6526ee1349b`；
- CMake `3.31.12`；CodeGraphF `0.9.8`，final sync exit 0、index up to date。

激活候选 design/reference/plan SHA-256 仍分别为
`4cc2052cf99874cefaffe1d11fc954195269fac77ea06fee06a92582ab5ea64e`、
`0349d4a6301b2213db98ba680557baab808bde62b40f09d99cf392af1e4058b8`、
`bed28804133174481d9ad41266a8d644a0ad056ff9be7ab48f44878b11c45f44`。

## 主 agent review

- exact `rg` caller audit 证明 body-fitted predictor/provisional/final transport callers 保留，
  private adapter bridge 是唯一新增 caller；CodeGraphF 对重载/大函数存在 caller
  under-report，因此未用空的 affected-tests 结果替代精确审阅和测试；
- pressure temporal RHS 使用 trial/committed/history 三层密度，pressure coefficient、
  exact Schur response、momentum diagonal 和 revision 消费同一 attempt authority；
- exactly two PISO、force sign、force-after-final-check、rollback、Restart 与 MPI 一致性
  authority 未修改；material+WALE 仍稳定拒绝，等待 S3-C2；material Restart 仍拒绝，
  等待 S3-R2；
- D1 未委派；R1 worker 结果仍隔离、未集成，O1 未启动；
- 无私有源码、研究数据、研究进程、push、publish 或 Stage 4–6 实现。

提交 subject：`feat: couple material density to immersed flow`

DCO：`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`

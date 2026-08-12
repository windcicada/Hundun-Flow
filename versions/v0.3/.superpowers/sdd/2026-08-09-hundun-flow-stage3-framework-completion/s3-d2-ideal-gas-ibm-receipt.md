# S3-D2 ideal-gas IBM vertical slice 验收回执

状态：`ACCEPTED`

accepted parent：`cd7f00c0d6749d40c200048f09144237de790244`

accepted at：`2026-08-09T17:27:57+08:00`

## 实现边界

本 task 实现 schema-v3 `LFP-GCIBM/none/ideal_gas` profile 7。driver 只放行完整
IBM、无 WALE、无 Restart 的 ideal-gas setup；body-fitted ideal 与 ideal+WALE 仍留给
S3-C3，density Restart 仍留给 S3-R2。

ideal adapter 只包装一个既有 `ImmersedMaterialDensityAdapter` 和既有
`DensityClosureHooks`，没有第二套 transport、第二 density field、第二 closure bridge 或
第三次 pressure corrector。固定阶段顺序为 predictor transport/closure → PISO #1 →
provisional transport/closure → PISO #2 → final corrected-flux transport/closure →
post-closure material assessment → final residual/force → 全部 prepare → 一次 collective →
noexcept publish。FlowState、material、p0/closure 与 pressure authority 均在 collective pass
后才发布。

closed p0 只使用 active fluid volume；inactive `rho/rho_h` 始终为 canonical positive zero。
open p0 与 configured value bitwise 固定，pressure outlet 只约束 `pi`。`cp`、`R` 与
configured p0 均通过 private closure adapter 与 attempt physics bitwise 匹配。public Stage 2
create/restore 仍以全部 owned cells 为 active，Checkpoint-facing closure state 未改变。

计划白名单之外的 `flow_material_density_{transport,piso}.cpp` 改动只为复用既有 friend
bridge 的 post-closure assessment 与 authenticated material report；三个 support header
改动只提供 test-only state/report mutation access。没有增加 public API 或生产期 test seam。

## TDD、RED 与 mutation 证据

初始可执行 behavior RED：两个新 target 均编译链接成功；1-rank direct/transaction 分别约
3.6 s 后 exit 8，输出 `IMMERSED_IDEAL_GAS_UNSUPPORTED`。失败来自 facade ideal branch
尚未实现，不是 missing symbol、CMake registration 或 fixture 错误。

主 agent 完整 diff review 另建立并关闭两项 RED：

- configured p0 one-ULP mismatch 原本被错误提交，1-rank exit 8、5.60 s；增加 private
  configured-p0 getter 与 bitwise setup 校验后，authenticated non-retryable rejection GREEN；
- `before_commit` collective failure 最初因 unavailable continuity/pressure 数值未归一而使
  material report bridge 自身无法认证；新增 prepare 后 failure case，按 authority availability
  归一失败报告，不改变成功路径或阈值，最终 1/2-rank bitwise rollback GREEN。

逐字段 report mutation 覆盖 194 种 material post-EOS evidence、28 个 closure report 字段，
以及 nested material、closure stage/candidate/seal、outer identity/presence/seal；原报告必须
authenticated，全部 mutation 同时被 outer seal 与 diagnostic source 拒绝。

主 agent 用 `apply_patch` 施加并立即恢复三项可编译 source mutation：

- 把 inactive volume 混入 closed-p0 分母：exit 8，86.49 s，死于
  `density_closure_failure`/未提交断言；
- 跳过 provisional closure：exit 8，5.28 s，死于 closure report unavailable；
- collective failure 前错误 publish closure：补齐 `before_commit` RED 后 exit 8，9.46 s，
  精确死于 closure state bitwise mismatch。

最终源码无 mutation/debug 残留。transaction 另覆盖 non-positive T/rho、closure layout
mismatch、p0 update 后 rank-local failure、prepare 后 rank-local failure，以及
history/committed/trial/metadata/final flux/controller/closure p0/revision bitwise rollback。

## GREEN / task gate

冻结的 Debug targets 全部无 warning 构建成功：

```text
test_immersed_ideal_gas
test_immersed_ideal_gas_transaction
test_ideal_gas_closure
test_ideal_gas_piso
test_ideal_gas_header_contract
test_stage3_flow_header_contract
hundun
```

header gate：`2/2 PASS`，0 failures，real 0.02 s。

v2 原样 focused gate：`8/8 PASS`，0 failures，real 44.94 s：

- ideal-gas IBM direct 1/2-rank：3.62 / 2.64 s；
- ideal-gas transaction 1/2-rank：9.75 / 6.56 s；
- schema-v3 flow-model fast acceptance 1/2-rank：12.13 / 7.60 s；
- existing ideal-gas closure/PISO 1-rank：1.02 / 1.62 s。

P0 registration/layout/mutation、Stage 3 source policy 与两项 header 的选定检查中 9 项
PASS。额外运行的历史 `test_task23_source_policy` 因寻找已经不存在的
`set_ideal_gas_restore` 源码标记而 FAIL；`git show HEAD` 证明 accepted parent 产品源码同样
没有该标记，来源是既有 policy 基线失配，不是 D2 diff。D2 不修改历史治理策略。

未运行 1/2/4 full matrix、正式 24/48³、96³、Release、ASan、UBSan 或 complete
sanitizer。环境没有 `clang-format`；`bash -n`、`git diff --check`、staged diff check 与编译器
warning policy 均 PASS。

## 调用图与单一 authority

CodeGraphF `0.9.8` final sync exit 0、index up to date。它能解析旧 closure callers，但对新
private bridge 与 `attempt_immersed` 大函数存在 caller under-report，因此用 exact `rg`
交叉验证：

- `make_ideal_report` 只有一个产品 caller；
- `attempt_immersed` 恰有两处 `correct_active_pressure`；
- predictor/provisional/final closure 各一处，final transport 与 post-closure assessment 各
  一处；
- `state.publish_commit_attempt()` 后只有一处 ideal adapter publish；同一 profile 没有第二套
  final-flux、pressure 或 closure authority。

## 证据身份

- implementation/test staged diff SHA-256（治理文件加入前）：
  `ee7fb42c077514abc6acf2dd028c87104b7c26dcadd377e1e26f069bac2687ff`；
- Debug `test_immersed_ideal_gas` SHA-256：
  `1c04a1a728a07748ea68e4ee27a5700c8f2439814ff701a7fe610d297a9d666a`；
- Debug `test_immersed_ideal_gas_transaction` SHA-256：
  `88f8b83a7620c7d46875fb6f36ecc1955af57ccd38118e82c265bfa69e6a1fae`；
- Debug `test_ideal_gas_closure` SHA-256：
  `75b80c3b9c0348c9c7ac18c85c3b80612652aadb0de7112cc704a0b6f26d31e7`；
- Debug `test_ideal_gas_piso` SHA-256：
  `704eb28175d196937462583ee349ad9be62962e64fe41c8c66b6e3429633de73`；
- Debug `hundun` SHA-256：
  `2957963198c7c0d8844e23f032c11da2e21e6cc037adb510847f3873bfac2e93`；
- Debug `test_stage3_flow_header_contract` SHA-256：
  `52486c9ac58b019e4f642f58c079e36dc5a0e6729962c69b133df4ab3b5ee2dd`；
- final focused `LastTest.log` SHA-256：
  `e3cfb313e5ba6323ee32ee1411dba0dca9c5e8c7c1ed9876255efa40f2f2e140`；
- CMake `3.31.12`；CodeGraphF `0.9.8`。

激活候选 design/reference/plan SHA-256 仍分别为
`4cc2052cf99874cefaffe1d11fc954195269fac77ea06fee06a92582ab5ea64e`、
`0349d4a6301b2213db98ba680557baab808bde62b40f09d99cf392af1e4058b8`、
`bed28804133174481d9ad41266a8d644a0ad056ff9be7ab48f44878b11c45f44`。

## 主 agent review

- active-volume mass、`h=cp*T` 与 EOS 分别满足 `5e-12`、`1e-12`、`1e-12`；inactive
  volume 的错误 p0 公式被独立排除；
- open p0 固定、pressure outlet 只约束 `pi`，closed p0/target mass/revision 与 rollback
  authority 完整；
- exactly two PISO、force sign、force-after-final-check、Task 11 阈值、Restart 与 MPI 最低
  失败 rank authority 未修改；
- D2 未委派；R1 worker 结果仍在 `coast/stage3-infrastructure-lane` 独立 dirty worktree，未
  集成，O1 未启动；
- 无私有源码、研究数据、研究进程、push、publish 或 Stage 4–6 实现。

提交 subject：`feat: couple ideal gas to immersed flow`

DCO：`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`

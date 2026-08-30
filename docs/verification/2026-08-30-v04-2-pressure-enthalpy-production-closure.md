# HUNDUN-FLOW V04-2 压力—焓产品闭环验收回执

日期：2026-08-30

裁决：`CONDITIONAL ACCEPT / V04-2 PRESSURE-ENTHALPY PRODUCT CLOSURE`

## 1. 裁决边界

本回执条件接受 pressure–enthalpy 主闭环、同目标时刻 refinement、失败回滚，以及
真实产品路径的 BE→BDF2 pressure–enthalpy 终端结果。以下三点不在本裁决中被提升为
“已完成”或“正式发布”：

- 两步 Re3900 的 momentum predictor limiter 均为
  `limited=true, theta=0, activations=1`，因此高阶 momentum 路径仍有数值质量问题；
- 长程稳定性、物理统计和正式发布门没有执行；
- 当前开发历史不 DCO-clean，融合前必须形成 signed squash 或少量逻辑 squash，并
  验证候选树与本回执验收的树等价。

## 2. 候选身份

- worktree：`/home/wyf/code_dev/hundun-flow-pressure-enthalpy-c1-production`
- branch：`codex/v04-pressure-enthalpy-c1-production`
- reopen snapshot：`2739f36cbb80886585eb5f7b32b625278b78f38b`
- pressure–enthalpy code commit：`d9006d261dd5247c7a2a8edc67a4fbb173544f15`
- tests-off `hundun` SHA-256：
  `c4e57d0fbf0f4076437781551fc5db4d3a48f6d8504bbbe4a96acedafbb6c60d`
- 最终 DCO-clean 候选 commit/tree：待封装和树等价校验后补录；本回执不虚构哈希。

对 `origin/main..candidate` 的 DCO 审计结果为 209 个提交：132 个没有 trailer，
75 个 `Signed-off-by` 与作者匹配，2 个只有 Codex trailer。因此当前分支历史不是
DCO-clean；不能因为最后一个提交带 sign-off 就声称整条分支合规。

## 3. 已证明的 pressure–enthalpy 主闭环

已实现并验证以下产品链路：

```text
accepted C2 state/flux
  -> exact nonlinear continuity/energy replay
  -> joint-L2 Armijo selection
  -> atomic p/h/rho/T/U/phi publication
  -> thermo + boundary/IBM + flux refresh
  -> fresh continuity-energy Jacobian/Schur solve at the same target
  -> independent EOS/C/E/mass/gauge terminal audit
  -> collective commit or exact rollback
```

IBM 产品路径使用 masked Cartesian spatial quasi-Newton `E_p/E_h`，不再使用导致能量
方向失真的双对角产品 fallback，也不再重复加入 `a0 V h rho_h`。证书仍明确声明
`full_nonlinear_jacobian=false` 和 `ibm_spatial_derivative=false`，没有把未实现导数冒充
完整 nonlinear Jacobian。

普通 PISO corrector 仍恰好为两个。若 C2 exact replay 未达到 continuity 和 energy
分量门，最多执行六次同目标 nonlinear refinement；每轮从刚接受的临时状态重建
thermo、边界/IBM、面通量、残差与 pressure–enthalpy 线性化，最终接受仍由 EOS、
continuity、energy、mass、gauge 五个分量门分别决定。

## 4. 本次封存契约修复

### 4.1 Continuity witness valid gate

CLI 失败输出现在先检查 `continuity_witness.valid`。无有效 witness 时不再把默认
`continuity=0`、局部索引零值和 `rank=-1` 伪装成真实连续性诊断；真实 CLI fixture
覆盖非连续性失败并验证 witness 字段被省略。详细 witness 仍只服务失败诊断，未进入
正常成功热路径。

### 4.2 联合 L2 globalization provenance

联合 L2 merit 使用新的 pressure–energy policy schema：

- schema tag：`0x7630347065676c32`（`v04pegl2`）；
- 新 provenance certificate：`0xe6d46d99f1321347`；
- 旧 L∞ 冻结 oracle：`0x9ef260037a8ea7d7`，保持字面不变。

mutation 测试证明旧 provenance 不能验证新联合 L2 policy；没有用当前候选代码重生成
旧 oracle 来证明自身正确。联合 L2 只负责选择耦合下降候选，continuity 和 energy 的
独立终端门没有被合并。

## 5. 真实 ProductDriver refinement 证据

现有小型确定性 fixture 直接运行真实 ProductDriver，没有手工构造 report。成功路径
实际触发两次 same-target refinement：

- `pressure_energy_refinement_solve_calls=2`，trajectory 长度为 4；
- ordinal 连续，target generation 不变，collective lineage 轮间互异且跨 rank 一致；
- 每轮 rank-local pressure state、numeric identity、flux provenance 和 linear identity
  均与已接受临时状态同步刷新；symbolic/hierarchy/workspace 的静态身份保持稳定；
- 每轮联合 L2 merit 下降，最后一次 refinement 首次同时通过 continuity 与 energy 门；
- 最终 EOS、continuity、energy、closed mass 和 gauge 五门分别通过，状态有限且正；
- 1-rank 与 2-rank 运行均通过必要的 collective target/lineage 语义检查。

成功 fixture 的终端分量为：EOS `0`、continuity 约 `2.2923e-14`、energy 约
`1.046043174e-10`、closed mass `0`，gauge 仅为浮点舍入量级。

失败 fixture 保留更严格的门，实际耗尽 capacity=6，记录六轮 refinement 和完整
trajectory 后拒绝；提交状态、历史和通量逐位精确回滚。既有 half retry、直接 half
及下一 BDF2 的逐位一致语义继续通过。由此分别证明了“真实 refinement 接受”和
“真实 refinement 耗尽后回滚”，不再用零 refinement 的 Re3900 startup 代替该证据。

## 6. 定向测试与 tests-off 构建

本轮受影响 focused 集合 9/9 通过：evidence V4 validator、真实 CLI witness、
globalization 单 rank/2-rank、PISO authority、product freeze、I/O product path，以及
ProductDriver refinement/retry 的 1-rank/2-rank。它覆盖联合 L2 schema/provenance
mutation、真实 refinement 成功与 capacity rollback 等本次改动。最后补入临时目录
删除哨兵和 trajectory→terminal C/E 位级绑定后，又复核 CLI 与 ProductDriver
1-rank/2-rank 共 3/3 通过。

Clang 测试使用了正确的 libc++ `LD_LIBRARY_PATH`，没有把运行库缺失误判为测试失败。
tests-off 增量构建成功，配置仍为 Release、Clang 15.0.6、libc++、MPI-3，
`HUNDUN_SOURCE_VERSION=v0.4`、`HUNDUN_BUILD_TESTS=OFF`，目标为
`build/c1-authority-release-tests-off/versions/v0.4/hundun`。

## 7. 64-rank Re3900 产品证据

### 7.1 已有独立 full 单步基线

独立审计此前已使用审计时的 tests-off 产品二进制、64 ranks 和冻结 core binding 运行
`case-full, dt=0.006`：

- 工件：
  `/home/wyf/code_dev/.benchmarks/v04-2-audit-full-one-step-2RDXYp/evidence.jsonl`；
- SHA-256：`1a18eb1b05fe544dabe56b5566d5e71262dfecc171bbd6061cb0b3e5dd0c65de`；
- `COMPLETED steps=1 time=0.006`，无 retry、无 refinement；
- continuity=`7.447131615171743e-7`，energy=`5.900106495213083e-7`，
  容差均为 `1e-6`；EOS、mass、gauge 均为 `0`；
- total linear iterations=`87`，max-rank step 约 `28.90 s`。

该证据确认正式 Re3900 对照目标为 `dt·U/D=0.006`；`dt=0.003` 是保守 half-step，
不是因为 `0.006` 过大而必须采用的时间步。

### 7.2 唯一两步 BE→BDF2 复核

代码和契约修复后只补了一次要求内的两步运行，没有重跑 `0.003/0.006` 单步、COAST
或完整 full/half 组：

- case：同一 `case-full, dt=0.006`；
- 工件：
  `/home/wyf/code_dev/.benchmarks/v04-2-conditional-two-step-Jr82mX/evidence.jsonl`；
- evidence SHA-256：
  `495b904f26f3dfb852e176dd84769746e59b5f8d1d310771202c7880d5a85e17`；
- schema：`HUNDUN_V04_EVIDENCE_V4`；
- 完成时间：step 1 为 `0.006`，step 2 为 `0.012`。

| step | 时间层 | fallback / retry | continuity | energy | EOS / mass / gauge | linear iterations | max-rank step |
|---:|---|---|---:|---:|---:|---:|---:|
| 1 | BE startup | false / false | 7.447131615171743e-7 | 5.900106495213083e-7 | 0 / 0 / 0 | 87 | 29.401 s |
| 2 | true BDF2 | false / false | 2.127108197142029e-7 | 1.636010211501308e-7 | 0 / 0 / 0 | 69 | 24.919 s |

两步的五个终端残差均低于 `1e-6`，`pressure_solve_calls=2`，
`pressure_energy_refinement_solve_calls=0`，termination 均为
`component_residuals_converged`。第二步 `requested_bdf_order=2`、`bdf_order=2`，且
`temporal_method_fallback=false`，所以真实 BDF2 时间层和 pressure–enthalpy 终端门
已经验证。该运行不负责证明 refinement；refinement 的产品证据由第 5 节的确定性
fixture 独立提供。

### 7.3 Momentum limiter 质量分类

两步的 `momentum_predictor_limiter` 均为：

```text
limited=true, theta=0, activations=1
```

接受规则没有人为要求 `theta>0`，所以这不否定已经通过的 pressure–enthalpy 五门和
BDF2 时间层。但 BDF2 步仍是全局 `theta=0`，必须分类为数值质量问题，不能静默把它
表述成高阶 momentum 产品路径已经完全验证。后续研究需要单独解释 limiter 为何在
这两个时间层退到全局低阶端点。

现有 evidence 不暴露实际 maximum convective CFL。本次没有仅为记录 CFL 扩张
evidence schema 或重构 fixed time control，因此回执不虚构该数值。

## 8. DCO 与可融合提交链

当前开发分支保留全部已有成果，但其 209-commit 历史不满足 DCO。正式交付前应：

1. 在明确的融合基线上形成一个 signed squash，或少量按逻辑拆分且每个均有有效
   `Signed-off-by` 的提交；
2. 比较封装前后 tree，证明与本回执通过测试和产品运行的代码树等价；
3. 对新候选范围重新审计全部 commit trailers；
4. 将最终候选 commit 和 tree 补录到本回执，再进入融合评审。

在这些动作完成前，本回执只给出 `CONDITIONAL ACCEPT`，不把某个末端签名提交当成
整条历史 DCO-clean 的替代证明。

## 9. 未执行项与范围声明

按用户约束，本次没有执行完整 CTest、长时间 Re3900、长统计、完整 full/half 时间
收敛组、COAST 重算、Stage 5/6 portable 或 Halo P0。既有且未受本修改影响的 MPI、
Restart、MMS、PISO temporal order 和 product pressure–energy temporal convergence
证据继续复用，但不在本回执中冒充新运行结果。

因此当前可以封存在 V04-2 task 节点的结论是：pressure–enthalpy 主闭环、真实
refinement/rollback 和 BE→BDF2 终端门已获证；momentum limiter 质量、长程物理验证、
正式发布门及 DCO-clean 融合链仍是明确保留条件。

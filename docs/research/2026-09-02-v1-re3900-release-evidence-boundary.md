<!-- SPDX-License-Identifier: Apache-2.0 -->

# HUNDUN-FLOW V1.0 Re=3900 文献与产品证据边界

日期：2026-09-02
状态：`RESEARCH DECISION / SUPERSEDED FOR V1.0 BY A USER-AUTHORIZED SCOPED POLICY`

> 2026-09-02 决策更新：用户随后明确授权“中短程测试 GREEN 后在主分支发布
> V1.0”。因此仓库在读取候选统计结果前新增了
> [`HUNDUN_V1_RE3900_MEDIUM_RELEASE_POLICY_V1`](../verification/v1.0-re3900-medium-release-policy.json)，
> 仅对 `v1.0.0` 建立新的、带 10D blockage/薄域限制的软件发布门。下文关于旧 v0.4
> Task 21 尚未完成的结论仍为真，但这些旧门不再阻断这一 scoped V1.0；相应能力也不得
> 被 V1.0 release 反向声称。

## 1. 裁决

当前 `20D x 10D x 3D`、128-rank 算例可以建立一个**可证伪的中短程
产品/物理 GREEN**：证明真实 ProductDriver 在这一域、网格和时间步上持续健康，形成完整、
不挑窗口的尾流统计，并与冻结 Parnaudeau 剖面和实验 Strouhal 作受限比较。

这个 GREEN 不能单独授权 `V1.0`。仓库现行权威发布谓词仍要求
`numerical correctness && robustness && COAST short performance && literature physical
accuracy && provenance`；规定顺序是
`focused -> full2 -> frozen -> full20 -> literature -> final`。V04-2 回执明确只接受数值产品
闭环，不接受 full20、长统计、完整 literature authority 或正式 release
（[`V04-2 closure plan`](../superpowers/plans/2026-08-30-v04-2-pressure-enthalpy-production-closure.md)，
“裁决边界”“V04-2 之外的正式发布门”；
[`V04-2 closure receipt`](../verification/2026-08-30-v04-2-pressure-enthalpy-production-closure.md)，
§1、§8）。把产品名从 v0.4 改为 V1.0 不会自动降低这些既有门槛。

如果项目决定让“中短程、10D 横向域”成为 V1 的最终物理门，必须先建立一份明确、版本化、
可机读且由人批准的新 V1 release policy，逐项说明它如何取代现行第 16 节/Task 21；不能在
看到 HUNDUN 结果后静默改写旧 predicate。本文没有作这种规格变更。

## 2. 一手来源与 authority 分层

### 2.1 控制规格

1. [`v0.4 Cartesian architecture`](../superpowers/specs/2026-08-12-hundun-flow-v0.4-cartesian-performance-architecture-design.md)
   §16 是现行 Re=3900 发布门：性能用 `20D x 20D x piD`、`480x480x48`、64 ranks；
   文献统计另用 HUNDUN periodic case，先发展 `150D/U`，再统计约 `2020D/U`
   （约 420 shedding periods）。
2. [`v0.4 implementation plan`](../superpowers/plans/2026-08-12-hundun-flow-v0.4-cartesian-performance-architecture.md)
   Task 20 Step 3、5 及 Task 21 Steps 1--7 冻结 gate 顺序、至少五组 full20 配对、
   `--require-complete` literature receipt 和最终 provenance 条件。
3. [`V04-2 formal-acceptance research`](2026-08-30-v04-2-coast-open-cfd-formal-acceptance.md)
   §8.2--§8.4 明确 full2、至少五组 paired full20、长程健康/统计以及缺失 force authority
   相互独立；两步成功不能替代它们。

### 2.2 冻结论文 authority

| 量 | 冻结身份 | 可作何种验收 | 不能怎样使用 |
|---|---|---|---|
| 15 组 Parnaudeau PIV 剖面 | `x/D=1.06,1.54,2.02`；每站 mean `u/Uc`、mean `v/Uc`、`u'u'/Uc^2`、`v'v'/Uc^2`、`u'v'/Uc^2` | scoped paper comparison；必须显示完整曲线、数字化坐标/数值误差和求解器采样区间 | 数字化 `Ey/Eq` 不是 PIV 实验置信区间，不能据此发明统一百分比门 |
| Parnaudeau 实验 `St` | `f_vs D/Uc = 0.208 +/- 0.002` | 当前唯一带论文报告不确定度的冻结 Re=3900 标量比较门 | 不能用后发表 CFD 表替代 |
| Parnaudeau `Lr/D=1.51`、`Umin/Uc=-0.34`、formation length `0.87` | 原论文直接表值，但没有各自报告的不确定度 | paper advisory；报告定义、点估计和 sampling interval | 不能自己补一个 pass band |
| Norberg `St` | `0.2108 @ Re=3704.6`、`0.2104 @ Re=4211.1` 为直接邻点；`0.208884 @ Re=3900` 为公式派生 | advisory corroboration | 非第二个 Re=3900 直接硬门 |
| Norberg `Cl_rms=0.083046` | sectional fit，公式派生 | advisory | 不等价于 periodic `3D` 或 `piD` 的积分升力 |
| Tadrist lift | Re=3000--7000 约 `0.02--0.08` 的正文范围 | target-neighbourhood advisory | 不是 Re=3900 点值，也不是实验不确定度 |

上述矩阵、公共 `y/D` 区间和每个 profile 的 `Ey/Eq` 见
[`thin-domain literature boundary`](2026-08-31-v04-re3900-thin-domain-literature-boundary.md)
“Accepted 15-profile matrix”“Meaning of the profile error bounds”。机器数据位于
[`cylinder-re3900-parnaudeau.json`](../references/cylinder-re3900-parnaudeau.json)
`experiment`、`observables`、`profile_authority.profiles`；Norberg 与 Tadrist 的限制分别位于
[`cylinder-re3900-norberg.json`](../references/cylinder-re3900-norberg.json)
`required_observables`/`pending_observables`/`excluded_as_release_oracles` 和
[`cylinder-re3900-tadrist.json`](../references/cylinder-re3900-tadrist.json)
`observables`/`excluded_as_release_oracles`。

冻结 Parnaudeau 原文为 14 页机构 PDF，SHA-256
`b8a775e5a5078e19fc47d9c5f47e95b81b4a68d4449e4444459088ed9befcdd4`。
PDF pp. 11--12、Figs. 11--15 是剖面来源；PDF p. 6、Table I 给出 LES 域、时间步、统计
时长和网格；PDF p. 7 说明先丢弃 `150D/U` 初始瞬态；PDF pp. 9--10 讨论统计收敛。

### 2.3 R4 receipt 仍不完整

[`v0.4-literature-data-receipt-r4-partial.json`](../verification/v0.4-literature-data-receipt-r4-partial.json)
字段 `complete=false`，`incomplete_references` 只有 Norberg；尚缺：

- Re=3900 的直接 total mean drag authority；
- 与目标 periodic span 可审计等价的 finite-span lift RMS authority。

本次只读重放结果是：三份 cylinder JSON `validate` 通过；R4 `receipt-validate` 通过；加入
`--require-complete` 后按预期以退出码 2 拒绝。当前文件身份是：

| 文件 | SHA-256 |
|---|---|
| R4 partial receipt | `1dfdeb16b92f4a29a5dc2a0c4daa146499e8decb1af5df197b652203e36b2e3f` |
| Parnaudeau JSON | `8a90682cf6d5d503538dabf47f8097a7e18f962eeeb5e53ce39858bcbb6361b9` |
| Norberg JSON | `4c8db4c241c72b383d68828d391e30b9f721128e4d1b0f83758b63df3ce6696f` |
| Tadrist JSON | `6788b727e72558dd87fd19b41e63e0030cc570ad772dc0a5e4cc5b2468e1705b` |

这允许 scoped comparison，不允许把 `literature`/`final` release gate 写成已通过。

## 3. Wang et al. (2024) 只作 peer reference

附件路径为
`/home/wyf/.codex/attachments/87217707-94dc-4958-abf8-330d4b71203c/An improved immersed boundary method with local flow pattern reconstruction and its validation _ Physics of Fluids _ AIP Publishing.pdf`。
该文件是 Wang Yudong et al., *Physics of Fluids* 36, 045145 (2024), DOI
`10.1063/5.0195598`；29 页，SHA-256
`f6d012a6bd62b2a755b4ecf2476d4783c03dbb5a51a36a60715aff53539d03e0`。完整审计另见
[`2026-09-01 Wang peer-reference audit`](2026-09-01-re3900-wang-peer-reference.md)。

原文对 cylinder case 明确给出的内容为：

| Wang 原文事实 | 页码 | 本项目用途 |
|---|---:|---|
| `Re=3900`、`D=0.02 m`、`Uc=2.89668 m/s` | PDF p. 11, §III A | 名义工况 peer identity |
| mesh spacing `0.2--0.8 mm`，即 `D/100--D/25`；`6.89e6` cells；64 blocks；圆柱附近坐标变换加密 | PDF pp. 11--12, Fig. 8 | 分辨率/结构化加密的 peer context，不是当前 face-coordinate certificate |
| 二阶中心空间离散、Crank--Nicolson、Smagorinsky | PDF pp. 10--11, §II C | 方法背景，不要求 HUNDUN 改算法 |
| mean profiles 位于 `x/D=1.06,1.54,2.02`；fluctuation profiles 同站位 | PDF pp. 13--15, Figs. 12--15 | 视觉 peer comparison；附件无可复用数值数组 |
| IB-method row：`Cd=0.9646`、`St=0.2344`、`umin/Uc=-0.2614`、`Lr/D=1.44` | PDF p. 15, Table I | 后发表数值计算结果；不是实验或 release oracle |

Wang 原文没有给出可审计的 cylinder 精确三维域、spanwise boundary、face coordinates、
时间步、发展时长、采样时长或统计不确定度。PDF p. 12 Fig. 9 只能视觉读出圆柱约位于
`x/D=5`、绘图区约到 `20D`；Fig. 8 也不能证明 span 或 periodicity。故它不能认证当前
`20D x 10D x 3D` case identity。

原文还存在必须保留的 `umin` 矛盾：PDF p. 12 正文给出最低平均速度 `-0.8867 m/s`，与
p. 11 的 `Uc=2.89668 m/s` 组合得到 `-0.306109062789`；p. 15 Table I 却给出
`-0.2614`。正文还把该最小值位置写成距圆柱中心 `0.0288 m = 1.44D`，数值恰等于 Table I
的 `Lr/D`，但 minimum-location 不是 mean-`u` zero crossing。两个 `umin` 必须分别标为
“正文自行归一化值”和“Table I 值”；不得平均或把差值伪装成 uncertainty，也不得从
`1.44` 的巧合倒推 `Lr` 定义。

## 4. 当前 10D x 3D case 的适用边界

Parnaudeau 实验的冻结身份是 blockage `4.3%`、effective span `20D`、PIV plane `z/D=0`
（Parnaudeau JSON `experiment`）。其 LES 是 `20D x 20D x piD`、圆柱距入口 `5D`、
spanwise periodic；LR LES 使用 `dt U/D=0.006`、`481x480x48`，统计 `2020D/U`
（约 420 cycles），HR LES 使用 `dt U/D=0.003`、`961x960x48`，统计 `250D/U`
（约 52 cycles）（Parnaudeau primary PDF p. 6, Table I）。

当前 smoke identity 则是 `20D x 10D x 3D`、约 10% blockage、tensor-stretched
`456x256x104`、128 ranks、`dt U/D=0.002`。网格/两步记录在
`/home/wyf/r39/report.md` 与 `/home/wyf/r39/manifest.json`：IBM strict-quadratic、MG、
两步 BE->BDF2、五个终端门、正性、双 CFL 和 AFC 均通过，但 `samples=0`、
`statistics_eligible=false`。因此现在只有 mesh/two-step product smoke，没有任何
时间平均或论文符合性结果。

域差异的后果必须显式保留：

- 横向 blockage 从实验 `4.3%`、标准 20D 计算域约 `5%` 变为约 `10%`；
- span 从 `piD` 变为 `3D`，且与实验 effective `20D` 均不等价；
- 当前 stretched smoke mesh 的 core/wake/outer 实际约为 `D/61`、`D/31`、`D/16`，
  不等于最初建议的 `D/80--D/100`、`D/40--D/50`、外部至少 `D/25`；
- 因而 `Cd`、`St`、`Umin`、`Lr` 或 profiles 的偏差不能直接归因为求解器。超出 scoped
  medium comparison 前，必须以只读旧 20D identity 做域敏感性对照；不得用 Wang
  Table I 抵消这一步。

## 5. 最小可证伪的中短程 GREEN

以下是**本任务的 scoped admission**，不是现行 Task 21 `literature`/`final` gate。
所有规则必须在读取该统计窗口结果前冻结。

### 5.1 窗口完整性

1. Fresh run 先发展至少 `20D/U`；若使用 exact candidate、mesh、BC、IBM、model 和
   accumulator provenance 均匹配的已发展 Restart，先重新稳定至少 `5D/U`。不满足 exact
   identity 时按 fresh 处理。
2. 随后连续采样至少 `50D/U`，并以实际观测的 shedding cycle 计数；至少包含 10 个完整
   周期。按冻结实验 `St=0.208`，`50D/U` 仅约 10.4 cycles，所以物理时间和实际周期数必须
   **同时**满足。
3. 统计按实际 `dt` 加权；失败 attempt、retry、temporal fallback 和 Restart recovery BE
   不计样。Restart 前后 accumulator、force/probe phase 和 committed time 必须连续。
4. 预先划分至少 5 个连续等长 blocks；报告每块、累计、前后半窗和 95% sampling interval，
   不得在看到结果后移动起点、删掉坏 block 或挑选有利窗口。

这只是最低 falsification window。Parnaudeau primary PDF pp. 9--10 明确警告：12 cycles
可以偶然产生“很好”的一、二阶剖面；近圆柱 mean statistics 通常需要超过 40 cycles；其
`Lr` 在 12/52/120 cycles 的粗略 sampling deviation 约为 `+/-0.5/+/-0.2/+/-0.1`，约
`+/-30%/+/-12%/+/-6%`，约 250 cycles 才看见 cumulative convergence。故 20+50 的
结果只能决定是否值得进入长程门，不能证明 converged literature statistics。

### 5.2 每步产品健康门

任一项失败即为 medium-short `RED`：

- Evidence V6 candidate/time/Restart identity 可验证，首条及连续时间锚点正确；
- 每个 accepted step 的 EOS、continuity、energy、closed mass、gauge 分量门分别通过，
  `p/T/rho` 有限且严格为正；
- provisional/committed CFL 都有效并低于配置上限；固定窗口不隐藏 retry 或 `dt` 改变；
- 正常统计步为真实 BDF2；Restart recovery BE 被明确标记并排除；
- AFC active 时 retained-correction L1 ratio 不能为零；`minimum_face_alpha=0` 和局部 limited
  face 本身不拒绝；无活动 correction 必须为 `not_applicable`；
- 没有质量/能量单向漂移、异常迭代突变、deadlock、非有限状态或未解释的 Restart 跳变。

这些门沿用
[`V04-2 formal-acceptance research`](2026-08-30-v04-2-coast-open-cfd-formal-acceptance.md)
§4.4、§5、§8.3 和最终 V04-2 Evidence V6 语义，不增设“所有面必须不受限”的伪门。

### 5.3 物理比较门

1. **Strouhal（唯一硬论文标量）**：用预注册 PSD/peak estimator 和 block bootstrap 给出
   HUNDUN `St` 的 95% sampling interval；该区间必须与 Parnaudeau 直接实验区间
   `[0.206,0.210]` 相交，且前/后半窗的主峰差必须能由各自 sampling interval 解释。
   不满足即 `RED`；相交只给 scoped GREEN，不代表 420-cycle 收敛。
2. **15 profiles（完整性与反例门）**：在精确站位、共同 `y/D` 区间和冻结归一化下输出全部
   15 组曲线；分别保存 paper extraction bound、solver interpolation error 和 block sampling
   interval。插值算子必须预注册；在每个 paper marker 处分别构造 paper extraction interval
   与 solver numerical-plus-sampling interval，并报告 overlap count、`L1/L2/Linf` residual 和
   对称/反对称 parity residual。任一站位/分量缺失、坐标/符号/归一化事后改变、选择性裁剪，
   某一完整 profile 在公共区间的 `overlap_count=0`，或冻结的 symmetry/parity 关系与 solver
   interval 明确冲突，均为 `RED`。由于 PIV 未报告逐点 experimental confidence interval，
   这一规则只排除明确反例；非零 overlap 也不能产生一个人工统一百分比“正式通过”阈值。
3. **统计稳定性**：St、profile norms、`Cd_mean`、`Cl_rms`、`Umin`、`Lr` 的 block sequence
   不得持续单向漂移。最小机器反例是 5 个 block point estimates 严格同向，且首末 block
   的 95% intervals 不相交；前后半窗 intervals 不相交也拒绝该 sampling window。此类结果
   记为 `SAMPLING INSUFFICIENT`，不是调整文献带或挑窗口。
4. **advisory / solver diagnostics**：Parnaudeau `Lr/Umin/formation length`、Norberg/Tadrist、
   Wang Table I、`Cd` 和 periodic-3D `Cl_rms` 全部单列。它们可以暴露反常，但不得冒充本
   medium gate 的论文硬阈值。

只有 5.1--5.3 全部成立，才可写：

```text
MEDIUM-SHORT 10D-BLOCKAGE Re3900 PRODUCT/PHYSICS GREEN
```

必须同句附带：`NOT CONVERGED LITERATURE STATISTICS / NOT V1 RELEASE ACCEPTANCE`。

## 6. V1 正式接受仍然缺少什么

| 阻断项 | 当前状态 | 关闭条件 |
|---|---|---|
| 正式 V1 policy identity | 仓库当前控制规格仍是 v0.4 §16/Task 21；没有一份已批准的新 policy 把 10D medium run 定义为 V1 final gate | 保持现行 predicate，或在看结果前批准、版本化并 mutation-test 一份明确替代它的 V1 policy |
| 新 exact machine candidate | 当前 10D smoke HEAD/mesh 是产品 smoke 身份，不是 Task 21 frozen candidate | clean exact HEAD/tree、tests-on/off、compiler/MPI、binary、case/STL/mesh、rank map、CPU、命令和 artifact hashes；候选改变后从 focused 重来 |
| robustness/focused | 10D 节点只跑了受影响 focused 和两步 smoke | 执行冻结 focused manifest，包括要求的 MPI、Restart/rollback、MMS/temporal、ASan/UBSan 和本候选直接受影响项 |
| COAST full2/full20 performance | 当前候选没有 resealed COAST equivalence，也没有至少五组正式 full20 | 按旧只读 COAST oracle 与现行 performance policy 完成等价性 reseal、full2、冻结及至少五组 paired full20；interval upper bound `<=1.0` |
| 域/网格 identity | 10D 横向、3D span、stretched smoke mesh 不等于冻结 `20D x 20D x piD` literature identity | 至少完成 10D/20D 域敏感性；正式门若不改规格，仍回到冻结 20D identity，并做必要空间分辨率敏感性 |
| 长程统计 | 当前 `samples=0`；即使 20+50 GREEN 也只有约 10 cycles | 现行门要求先发展 `150D/U`、后统计约 `2020D/U`/420 cycles，并通过长程健康、Restart accumulator 和 block convergence |
| complete literature authority | R4 `complete=false`，Norberg 两个 force rows pending | `receipt-validate --require-complete` 必须通过；Wang、sectional lift、pressure drag 或二手 CFD 表不能补位。若要删除这些量，必须先正式修改 release spec，而不是改 receipt 伪造 complete |
| final provenance/release | 尚无本候选的最终 ACCEPT receipt、V1 version/package/tag/release seal | 完整 diff、DCO、license/third-party、tests-off isolation、工件不可变、无遗留 MPI 进程、一个最终 ACCEPT；随后才可执行单独的 main/tag/push/GitHub release 操作 |

用户对“论文未提供的 total `Cd` 或严格等价 `Cl_rms` 不阻断”的要求适用于 scoped
thin/medium comparison，已经由
[`thin-domain literature boundary`](2026-08-31-v04-re3900-thin-domain-literature-boundary.md)
§“Comparison and verdict boundary”明确支持；它没有自动修改上表中的现行正式 release
predicate。两者必须保持可区分。

## 7. 禁止的结论

在上表全部关闭前，不得声称：

- `ACCEPT / HUNDUN-FLOW V1.0 RELEASE` 或 v0.4/V1 literature/final gate 已通过；
- 两步 smoke、20+50 medium window 或一个低 residual 已证明长程稳定/统计收敛；
- 10% blockage 的 10D 结果等价于 Parnaudeau 4.3% experiment、标准 20D case 或完整展向域；
- Wang `Cd/St/Umin/Lr` 是实验目标、可补齐 R4，或当前 case 已复现 Wang 的未披露域/采样；
- periodic `3D` integrated `Cl_rms` 等价于 Norberg sectional lift 或 Tadrist finite-span force；
- paper digitization error 是实验 uncertainty，或 profile 可以用事后统一百分比/选择性区间接受；
- HUNDUN 与 COAST 相互接近等于论文物理验证；COAST 始终只是性能/solver diagnostic，
  不是数值 oracle。

中短程 GREEN 的正确用途是：尽早证伪产品路径、域/网格和明显物理错误，并决定是否值得
投入正式域敏感性与长统计。它不是用较短运行重新命名正式发布门。

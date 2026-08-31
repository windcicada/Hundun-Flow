# HUNDUN-FLOW V04-2 压力—焓产品闭环计划

日期：2026-08-30；最终复核：2026-08-31

状态：`ACCEPT / V04-2 NUMERICAL PRODUCT CLOSURE`

## 裁决边界

本计划只关闭 V04-2 数值产品节点：同目标时刻 pressure–enthalpy 闭环、真实
ProductDriver refinement/rollback、MPI common-face AFC authority、provisional/committed 双
CFL，以及可证伪的 Evidence V6 provenance。

这不是 `HUNDUN-FLOW v0.4 Re=3900 JOINT RELEASE GATE`。本轮没有运行 full20、COAST
重算、长程稳定性/物理统计、完整文献 authority 或正式 release gate，也没有融合 Stage 5/6
portable 或 Halo P0。

## 闭合策略

IBM 不相容冷启动是原始压力失控的触发器；持续放大来自 pressure–enthalpy 分裂和 BDF2
历史链。修复后的产品路径在同一目标时间层同步更新 `p/h/rho/T/U/phi`，每轮 refinement
都从已接受的临时状态重建热力学、边界/IBM、质量通量、连续性—能量残差和
Jacobian/Schur 线性化。联合 L2 merit 只选择候选，EOS、continuity、energy、closed mass
和 gauge 仍分别终端验收。标准 BE/variable-BDF2 公式没有改变。

原有 pressure–enthalpy、refinement 和 capacity=6 回滚证据保留。本次重新打开任务所增加的
四个可证伪缺口按以下顺序关闭。

## 四个规格缺口与 RED→GREEN

| 缺口 | 旧实现的可证伪问题 | 修复与 GREEN 证据 |
|---|---|---|
| MPI partition-face AFC authority | 两侧只交换 cell ratio，之后各 rank 独立计算 shared-face alpha；没有唯一最终 authority，计数也可能重复 | `common_face_afc_v3_owner` 由负侧唯一 owner 发布三分量共享 alpha，正侧只消费 owner ghost；等量反号修正。partition-face mutation 故意污染 consumer alpha，1/2/4 ranks 仍得到相同 shared alpha、全局 RHS、守恒、active/limited/fraction |
| 不可变候选身份与首条时间 | V5 的 build/binary 是固定常量；validator 从第二条才验证 `delta(time)==dt`，restart 首条无可信锚点 | 新增 `HUNDUN_V04_EVIDENCE_V6`：绑定 exact HEAD/tree、canonical build manifest、运行中 executable SHA-256 和组合 identity；每行含 `previous_committed_time` 与 fresh/restart `run_start`。restart 必须提供实际 `manifest.bin`，并校验 FNV、SHA-256、step/time。首条错配、混合候选、三时间量与 manifest 一起平移均拒绝 |
| committed CFL 证书不完整 | 只有 out/abs max、limit、revision，不能证明 winner、rhoV、质量流量、dt 和真实 view | provisional/committed 共用单一 cell-CFL 内核；typed committed certificate 绑定 valid、density/final-flux revision、collective view identity、activity、dt、limit，以及 out/abs 各自 winner 的 global cell/rank、rhoV、outgoing/absolute mass flow。超限路径发布详细 witness，fixed fatal 与 adaptive retry/rollback 事务语义不变 |
| AFC 聚合量不足 | 只有 retained L1 和 limited faces；无法判断活动面分母、最小 alpha 和 N/A | V6 增加 `minimum_face_alpha`、`active_correction_faces`、`limited/active` fraction；只由唯一面 owner 聚合。无活动 correction 明确为 `not_applicable`，三个 ratio/min 字段为 null，不伪造 0 |

独立复审另外发现并关闭两条证据真实性反例：构建系统现在把与 executable identity 相同范围
的源码/index/ref 输入纳入自动重新 configure；外部 validator 对 restart/continuation 必须读取
真实 sidecar manifest，不能靠同步平移 JSON 中的三个时间量伪造合法 continuation。C++ 与
Python 对 restart SHA 的合法性判断也已统一。

## 冻结接受规则

- partition/periodic shared face 只有一个 alpha authority；两侧使用完全相同的 alpha，修正等量反号。
- active/limited 计数、retained ratio、minimum alpha 和 limited fraction 均按唯一物理面聚合。
- 成功步必须分别通过 EOS、continuity、energy、closed mass、gauge、有限性和严格正性门。
- 普通 pressure corrector 恰好为 2；same-target refinement 单独计数，容量为 6，失败逐位回滚。
- provisional 和 committed outward CFL 使用相同公式与配置 limit；absolute CFL 同时记录。
- fixed CFL 超限 fatal 且回滚，adaptive 超限 rejected-step 后缩步，不能隐藏 retry。
- V6 文件不得混 candidate identity；首行从 fresh 或经实际 restart manifest 认证的锚点验证
  `time - previous_committed_time == advective_cfl.dt`，后续逐行验证时间增量。
- 旧 V1–V5 schema/oracle 保持冻结，不以当前候选重生成旧 oracle 证明自身正确。

## 最小验证计划及完成状态

没有运行完整 CTest。最终代码候选执行并通过：

1. 14/14 focused：Evidence validator、solver-equations MPI 1/2/4、core-product-freeze MPI
   1/2/4、restart MPI-4、I/O product path、app-driver、restart rank-change、pressure–energy retry
   MPI 1/2/4；
2. `v04_product_pressure_energy_temporal_convergence`；
3. `v04_app_cli_continuity_witness`；
4. `python3 tools/v04_evidence_validate.py self-test`；
5. Release/Clang/libc++/MPI、tests-off 增量构建；
6. 最终 tests-off 二进制的 64-rank no-IBM fixed `dt=0.006` 两步 BE→BDF2；
7. 同一二进制的 64-rank IBM Re3900 fixed `dt=0.006` 两步 BE→BDF2。

独立 Spec 与 Standards 复审主动检查了 validator 错误接受、partition 两侧不同 alpha、restart
锚点平移、stale build identity 和 receipt/binary 不一致反例。发现的 no-IBM 时间步和旧封存文档
问题已修正；最终没有剩余 V04-2 阻断。

## 冻结候选

- 独立 worktree：`/home/wyf/code_dev/hundun-flow-v04-2-spec-closure`
- branch：`codex/v04-2-spec-closure-20260831`
- 运行代码提交：`2bc6fd5f89732dace0e7cd016ce4bd5dfdb64016`
- 运行代码 tree：`ff6eaf7668565df117e44b93c210519f8075e24a`
- build manifest SHA-256：`f367db108667550c718fc0f8bb59143be4fe6874c37ea95c31769387511c13ec`
- tests-off executable SHA-256：`13933c7d19c3659ef7e5aabab512f452e96bc6e4637e38f6b0a6bdb283831784`
- runtime candidate identity SHA-256：`2bebb5dee34375ffc55690230c7c02e4cc861d59032c96604f55fa6350fbc7c5`
- 封存 manifest：
  `/home/wyf/code_dev/.benchmarks/v04-2-final-re3900-VcRKDg3f/closure-manifest.json`
- 封存 manifest SHA-256：`718110b21db184439c9fe6a7ee38fd4c9fab58e7620deb8b0bfd466efadbd0e2`

文档封存提交位于运行代码提交之后；数值工件只绑定上述代码提交、tree、build manifest 和
executable bytes，不把后续纯文档提交冒充运行候选。

## V04-2 之外的正式发布门

后续正式 v0.4 release 仍须从 exact 候选重新执行规定顺序：new candidate、COAST scientific-
work equivalence reseal、`focused -> full2 -> frozen -> full20 -> literature -> final`、至少五组
paired full20、完整一手文献 authority、可 Restart 的统计 accumulator，以及长程正性、双 CFL、
AFC retained ratio、五门、retry/变步和守恒漂移监测。在这些门完成前不得创建或声称
`v0.4-re3900-final-acceptance.md`。

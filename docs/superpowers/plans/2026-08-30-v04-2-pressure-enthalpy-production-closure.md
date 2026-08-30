# HUNDUN-FLOW V04-2 压力—焓产品闭环计划

日期：2026-08-30

状态：`CONDITIONAL ACCEPT / V04-2 TASK NODE`

## 裁决含义

V0.4 的 pressure–enthalpy 主闭环、同目标时刻 refinement、失败回滚和 BE→BDF2
终端物理门已经得到产品路径证据，可以在 V04-2 节点作条件接受。这里的
`CONDITIONAL ACCEPT` 不是正式发布裁决；它保留的是数值质量和长期物理门，而非
提交历史不确定性：

- 两步 Re3900 的 pressure–enthalpy、BDF2 时间层和五个独立终端门已经通过；
- 两步 momentum predictor limiter 均为 `limited=true, theta=0, activations=1`，因此
  高阶 momentum 路径不能声明为已完全验证，这是明确保留的数值质量问题；
- 长程稳定性、物理统计和正式发布门没有执行；
- 原 209-commit 开发历史已保存在备份分支；可融合分支已从 `origin/main` 形成
  DCO-clean signed squash，并验证其 tree 与已验收树完全等价。

## 目标与边界

本阶段沿用既有工作树，不重做已经完成的压力—焓 RED、IBM 相容初始化、边界和 BDF
历史工作。目标仍是让真实 V0.4 产品路径在同一目标时间层上以一致的
`p/h/rho/T/U/phi` 状态通过 EOS、continuity、energy、mass 和 gauge 五个独立终端门，
并补齐两个封存契约以及真实 ProductDriver refinement 证据。

本阶段不融合 Stage 5/6 portable 或 Halo P0，不扩张时间控制和 evidence schema，
也不运行完整 CTest、长时间 Re3900、完整 full/half 收敛组或 COAST 重算。已有且
未受本次修改影响的测试证据继续复用。

## 根因裁决

IBM 不相容冷启动仍是压力失控的触发器，但不是持续放大的主因。有限差分 witness
把首个决定性不一致定位到 IBM 产品路径的能量—焓空间块：旧路径以双对角
quasi-Newton fallback 代替完整 Cartesian `E_p/E_h` 空间响应，并在 `E_h` 中再次
加入已经存在的 `a0 V h rho_h` 时间项。结果是 exact nonlinear replay 与 Schur
方向不一致，C2 只能选取下降候选，不能进入独立 continuity/energy 终端门。

本阶段没有修改标准 BDF2 系数，也没有把 IBM 导数冒充 full nonlinear Jacobian。
修复后的 IBM authority 明确标记为 masked Cartesian spatial quasi-Newton；缺失的
IBM donor-gradient 导数仍不作虚假声明，最终物理门继续作为接受 authority。

## 工作包与当前状态

| 工作包 | 状态 | 条件接受证据或剩余动作 |
|---|---|---|
| RED / Jacobian / Schur witness | 完成 | continuity 与 energy 分量独立比较，定位首个 `E_h` 空间不一致 |
| IBM pressure–enthalpy 空间块 | 完成 | activity-aware flux-only `E_p`、masked Cartesian `E_h`，移除时间项重复计数 |
| C2 same-target refinement | 完成 | 每轮从已接受临时状态刷新 thermo、边界/IBM、面通量、残差和线性化 |
| 原子状态同步 | 完成 | refinement 同步发布 `p/h/rho/T/U` 与物理质量通量，authority 单次消费 |
| globalization policy identity | 完成 | 联合 L2 使用新 schema/provenance；旧 L∞ 冻结 oracle 保持不变且 mutation 被拒绝 |
| CLI continuity witness | 完成 | 仅 `valid=true` 时输出 witness；其他失败不再伪装默认零值和 `rank=-1` 为真实诊断 |
| ProductDriver refinement | 完成 | 成功 fixture：calls=2、trajectory=4、五门通过；失败 fixture：capacity=6 后精确回滚；1/2-rank 通过 |
| rollback / retry | 完成 | 未接受 full step 精确回滚，half retry、直接 half 和下一 BDF2 语义保持一致 |
| Re3900 BE→BDF2 | 条件完成 | 唯一两步 `dt=0.006` 运行通过 pressure–enthalpy 五门且无 retry；两步无 refinement |
| momentum 高阶质量 | 待后续研究 | 两步均触发全局 `theta=0`；不得据此声明高阶 momentum 路径已验证 |
| DCO-clean 融合链 | 完成 | 从 `origin/main` 形成 signed squash；保留原历史备份并验证候选 tree 完全等价 |
| 长程与正式发布门 | 未执行 | 长程稳定性、物理统计和正式发布验收不属于本次短闭环 |

## 已冻结的接受规则

- globalization 可以接受联合 L2 下降，但不能替代 continuity 和 energy 分量终端门。
- 成功步必须分别满足 EOS、continuity、energy、closed mass 和 gauge 容差。
- `p_abs`、`rho` 和 `T` 必须在每个本地单元有限且严格为正。
- 普通 `pressure_solve_calls` 必须仍为 2；same-target refinement 使用独立计数，最大 6。
- refinement 报告必须形成从 1 开始的连续 ordinal，保持同一 target generation；每轮
  刷新 rank-local pressure/numeric/linear identity，并同步刷新跨 rank 一致的 collective
  state/flux provenance；collective lineage 跨 rank 一致且轮间互异，未使用 suffix 必须为空。
- 联合 L2 策略使用新的 policy schema/provenance identity；旧 L∞ 冻结 oracle 不得重生成。
- continuity witness 只在 `valid=true` 时可作为失败诊断输出，正常成功热路径不增加详细
  witness collective。
- V3 冻结 oracle 保持字面不变，V3 携带任一 V4 refinement 字段必须被拒绝。

## 验证分层

1. 主闭环证据：原 focused 9/9、PISO temporal order 和 product pressure–energy
   temporal convergence 已通过，未受本次契约修复影响的结果继续复用。
2. 契约与 refinement 证据：受影响 focused 集合 9/9 覆盖 L2 policy provenance
   mutation 和真实 ProductDriver 成功/回滚；最后的安全哨兵与 terminal 绑定补强后，
   CLI 及 ProductDriver 1/2-rank 又复核 3/3。
3. Re3900 基线证据：独立 `case-full, dt=0.006` 单步 startup 已通过，证明
   `dt=0.003` 只是保守 half-step，而非必要的时间步缩小。
4. Re3900 时间层证据：唯一一次 64-rank 两步 `dt=0.006` 运行完成真实
   BE→BDF2，无 temporal fallback、无 retry，第二步五个终端门通过。
5. 未证明项：实际 maximum convective CFL 未由现有 evidence 暴露；按范围约束不为此
   扩 schema。长程稳定性、物理统计、正式发布门及 momentum 高阶质量仍是后续门。

## 候选与 DCO 封装

- worktree：`/home/wyf/code_dev/hundun-flow-pressure-enthalpy-c1-production`
- branch：`codex/v04-pressure-enthalpy-c1-production`
- reopen snapshot：`2739f36cbb80886585eb5f7b32b625278b78f38b`
- pressure–enthalpy code commit：`d9006d261dd5247c7a2a8edc67a4fbb173544f15`
- 当前 tests-off `hundun` SHA-256：
  `c4e57d0fbf0f4076437781551fc5db4d3a48f6d8504bbbe4a96acedafbb6c60d`
- DCO-clean code commit：`59985dbf2d10dee167eddfdd38ae5c2752051605`
- 已验收 code tree：`7d172adbb1feb3ad4f692d303ef991d100e26252`
- 原始开发历史备份：`codex/v04-pressure-enthalpy-c1-production-pre-dco-20260830`

历史审计锚定 `origin/main..2739f36cbb80886585eb5f7b32b625278b78f38b`：209 个提交
中 132 个没有 `Signed-off-by` trailer，75 个 trailer 与作者匹配，2 个只有 Codex
trailer，因此没有机械修改旧作者历史。封装后的 code commit 与封装前候选
`abfdb7f2cf6f2514af721606ace625ee669c9c0b` 的 tree 均为上述 `7d172a...`；可融合范围只
保留由 WANG YUDONG 签署的 squash 和随后的签署回执提交。

具体命令、残差、运行代价、工件和条件接受边界见
`docs/verification/2026-08-30-v04-2-pressure-enthalpy-production-closure.md`。

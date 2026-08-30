# HUNDUN-FLOW V04-2 压力—焓产品闭环计划

日期：2026-08-30

状态：`ACCEPT / V04-2 NUMERICAL PRODUCT CLOSURE`

## 裁决边界

本计划已经完成 V04-2 数值产品节点。接受范围包括 pressure–enthalpy 同目标时刻闭环、
真实 ProductDriver refinement/rollback、common-face conservative AFC、provisional 与
committed 双 CFL authority，以及最终 tests-off Re3900 的 BE→BDF2 两步行为。

这不是 `HUNDUN-FLOW v0.4 Re=3900 JOINT RELEASE GATE`。正式发布仍需新的不可变候选、
COAST equivalence reseal、正式 full2、冻结、至少五组 paired full20、完整一手文献
authority、长程稳定性/物理统计和最终 provenance。V04-2 接受不能替代这些门。

## 已闭合的根因链

IBM 不相容冷启动是原始压力失控的触发器；持续放大来自 pressure–enthalpy 分裂和 BDF2
历史链。修复后的产品路径在同一目标时间层内同步更新 `p/h/rho/T/U/phi`，并用独立
EOS、continuity、energy、closed mass 和 gauge 门接受。C2 refinement 每轮从已接受的
临时状态重建 thermo、边界/IBM、质量通量、残差和 Jacobian/Schur 线性化，不复用旧方向。

标准 BE/variable-BDF2 公式没有被修改。所有方程装配、AssemblyEpoch、pressure-storage
直达入口和 PISO terminal 均验证 BDF 系数与证书 `dt` 一致；`dt` 进入 assembly certificate
和 momentum numeric identity。IBM Jacobian 仍诚实标记为 masked Cartesian spatial
quasi-Newton，没有冒充完整 nonlinear IBM donor-gradient Jacobian。

## 完成的工作包

| 工作包 | 结果 | 关键证据 |
|---|---|---|
| Pressure–enthalpy RED / Jacobian / Schur | 完成 | continuity、energy 有限差分与 Schur 预测独立比较；联合 L2 只选候选，分量门独立接受 |
| C2 same-target refinement | 完成 | 真实 ProductDriver 成功路径 calls=2；capacity=6 失败路径逐位回滚；ordinal/target/lineage/刷新 authority 通过 |
| CLI 终端诊断真实性 | 完成 | audit 前失败输出 `terminal_audit=unavailable`；五残差仅在 terminal audit present 时输出；witness 还要求 `valid=true` |
| Globalization provenance | 完成 | 联合 L2 使用 `v04pegl2`；旧 L∞ oracle 冻结，mutation 拒绝跨策略复用 |
| Momentum 正性限制 | 完成 | `common_face_afc_v2`，三分量共享 face alpha，owner/halo 成对守恒，smooth-extremum O(dx²) allowance，物理出流单边预算 |
| Momentum authority / MPI 原子性 | 完成 | velocity/geometry/boundary/BDF-dt、density exact view、三轴 flux exact view 绑定；本地错误在 FGMRES 和状态修改前 consensus |
| 双 CFL | 完成 | predictor provisional flux 与 terminal committed flux 分别计算，revision 不同、limit 相同；fixed/adaptive stage 60 策略和回滚有真实 driver fixture |
| Evidence V5 | 完成 | AFC L1 retained ratio、唯一受限面、双 CFL；文件内 schema/run identity/连续 step/static authority/time-dt/IBM activity 均 fail closed |
| 最终 no-IBM 对照 | 完成 | 64 ranks、两步、BE→BDF2、无 retry/fallback、五门与双 CFL 通过 |
| 最终 IBM Re3900 | 完成 | 64 ranks、`dt=0.006` 两步、无 retry/fallback；第二步 retained ratio `0.3649577`、双 CFL `<0.384` |
| DCO / 封存 | 完成 | 数值实现提交 signed；计划、回执和 ledger 使用独立 signed 文档提交 |

## 冻结的接受规则

- 联合 L2 merit 只能选择耦合下降候选；continuity 和 energy 必须分别过门。
- 成功步必须分别满足 EOS、continuity、energy、closed mass 和 gauge 容差。
- `p_abs`、`rho`、`T` 和全部提交状态必须有限且严格为正；失败 attempt 不得污染提交层。
- 普通 `pressure_solve_calls` 恰好为 2；same-target refinement 单独计数，容量为 6。
- momentum 限制量以 retained-correction L1 ratio 判断全域高阶保留；单个面 `alpha=0`
  不是全域坍缩，接受步的全局 retained ratio 必须严格大于 0。
- provisional 和 committed outward CFL 都必须不超过同一配置 limit；absolute CFL 作为诊断
  同时记录。fixed 超限 fatal，adaptive 超限 rejected-step 后事务缩步。
- CFL 证书必须绑定真实 density/flux rank-local view；相同 revision/storage/domain 的 foreign
  replica 不能冒充真实输入。
- V5 单文件不得混 schema 或候选；step 连续，`delta(time)` 与当前 CFL `dt` 一致，momentum
  plan、IBM activity 和 CFL limit 在运行内不变。
- 旧 V1–V4 工件和 oracle 保持冻结，不以当前实现重生成来证明自身正确。

## 最终最小验证集

最终候选没有运行完整 CTest、COAST 重算、长 Re3900 或 full/half 收敛组。执行并通过：

1. 12 项定向 focused 集合：CLI、V4/V5 validator、I/O、app-driver、pressure-storage、
   Schur/globalization、PISO authority、solver-equations 1/2-rank、refinement/retry 1/2-rank；
2. 当前 AFC v2/V5 契约直接修改的
   `v04_product_pressure_energy_temporal_convergence` 与
   `v04_restart_rank_change_continue`，定向复核 2/2；
3. `v04_core_product_freeze_mpi_test --mass-flow-only`，2 ranks；
4. Evidence V5 Python self-test；
5. tests-off Release/Clang/libc++/MPI 增量构建；
6. 64-rank、两步 no-IBM control；
7. 64-rank、两步 IBM Re3900 `case-full, dt=0.006`。

只有未与本轮 AFC/CFL/authority 契约重叠的既有 MMS、PISO temporal order 和
pressure–enthalpy/refinement 证据直接复用；被本轮改写的 temporal-convergence 与 Restart
检查使用当前版本重新运行，不以旧结果替代。

## 候选与最终工件

- worktree：`/home/wyf/code_dev/hundun-flow-pressure-enthalpy-c1-production`
- branch：`codex/v04-pressure-enthalpy-c1-production`
- 数值闭合提交：`287d6d7341f6fd7f61e553f7543f9e9615a90fab`
- 数值闭合 tree：`dfbedcf0646a173fbbeea9711d39eff60666e25a`
- tests-off `hundun` SHA-256：
  `babd0c407d5fd36d0663d4e94a462b2bdc54e09e96e33b20cb381a558dfe37f6`
- no-IBM V5：
  `/home/wyf/code_dev/.benchmarks/v04-2-afc2-noibm-accepted-ETZECl/run/evidence.jsonl`
- no-IBM evidence SHA-256：
  `91128e72607444270e93e0d98d8f2537d9f9b2305ef7f580ad248e6554c69891`
- IBM Re3900 V5：
  `/home/wyf/code_dev/.benchmarks/v04-2-afc2-re3900-accepted-SmS6Jy/run/evidence.jsonl`
- IBM evidence SHA-256：
  `be1f744b83a3fb33fd9446b67182c271ab9027228b6e45c05f3e7387f90ced58`

## V04-2 之外的后续门

V04-2 节点到此停止，不融合 Stage 5/6 portable 或 Halo P0。正式 v0.4 release 仍需：

1. 从最终 exact HEAD/tree/binary 建立新的 `HUNDUN_V04_CANDIDATE_V1`；
2. 对同一候选重新封印 COAST scientific-work equivalence；
3. 按 `focused -> full2 -> frozen -> full20 -> literature -> final` 顺序运行；
4. 至少五组 alternating paired full20 并满足冻结性能 policy；
5. 关闭 Re3900 total mean drag 和 periodic `pi D` `Cl_rms` 的一手 authority；
6. 修复/验证可 Restart 的统计 accumulator，并完成至少 `150D/U` 发展和约 `2020D/U`
   统计，持续监测双 CFL、AFC retained ratio、五门、正性、retry/变步和守恒漂移；
7. 完成最终 provenance、DCO、不可覆盖工件和发布回执审查。

在这些门完成前，不创建或声称 `v0.4-re3900-final-acceptance.md`。

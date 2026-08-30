# HUNDUN-FLOW V04-2 压力—焓产品闭环验收回执

日期：2026-08-30

裁决：`ACCEPT / V04-2 NUMERICAL PRODUCT CLOSURE`

## 1. 裁决范围

本回执正式接受 V04-2 数值产品节点：pressure–enthalpy 同目标时间层闭环、真实
refinement/rollback、common-face conservative AFC、双 revision CFL authority，以及最终
tests-off Re3900 的 BE→BDF2 两步路径。

本回执不是 HUNDUN-FLOW v0.4 正式发布接受。它不替代新的 machine candidate、COAST
equivalence/full2/full20、完整文献 authority、长程稳定性/统计和最终 release provenance。

## 2. 候选身份

- worktree：`/home/wyf/code_dev/hundun-flow-pressure-enthalpy-c1-production`
- branch：`codex/v04-pressure-enthalpy-c1-production`
- 数值闭合提交：`287d6d7341f6fd7f61e553f7543f9e9615a90fab`
- 数值闭合 tree：`dfbedcf0646a173fbbeea9711d39eff60666e25a`
- tests-off `hundun` SHA-256：
  `babd0c407d5fd36d0663d4e94a462b2bdc54e09e96e33b20cb381a558dfe37f6`
- 编译配置：Release、Clang 15.0.6、libc++、MPI-3、
  `HUNDUN_SOURCE_VERSION=v0.4`、`HUNDUN_BUILD_TESTS=OFF`
- rank layout：64 ranks，`ppr:32:socket:PE=1`，slot rank，core binding，
  `OMP/OPENBLAS/MKL_NUM_THREADS=1`

该数值提交由 `WANG YUDONG <wangyudong@buaa.edu.cn>` 作者签署。原始长开发历史仍保存在
pre-DCO 备份分支，没有机械补签；可融合链从既有 DCO-clean squash 继续叠加 signed commit。

## 3. Pressure–enthalpy 主闭环

已验证的真实产品链为：

```text
accepted same-target state
  -> exact nonlinear continuity/energy replay
  -> joint-L2 candidate selection
  -> atomic p/h/rho/T/U/phi publication
  -> thermo + boundary/IBM + flux refresh
  -> fresh continuity-energy Jacobian/Schur solve
  -> independent EOS/C/E/mass/gauge audit
  -> collective commit or exact rollback
```

联合 L2 merit 只选择下降候选，continuity 与 energy 仍分别验收。普通 PISO corrector 恰好
两个；最多六轮 refinement 使用独立计数。真实小型 ProductDriver fixture 已覆盖：

- 成功路径实际执行两轮 refinement，ordinal 连续、target 不变、lineage 跨 rank 一致且轮间
  不同，每轮 state/thermo/boundary/IBM/flux/Jacobian authority 刷新，最终五门通过；
- 更严格路径实际耗尽 capacity=6，拒绝后提交状态、历史和通量逐位回滚；
- 1/2-rank 成功、回滚、half retry、direct half 和下一 BDF2 语义通过。

最终 Re3900 两步本身不需要 refinement；它证明正常产品路径，不冒充 refinement fixture。

## 4. Momentum AFC 与 authority

旧 global-theta limiter 已替换为 `common_face_afc_v2`：

- 每个物理面由三个速度分量共享一个 alpha，面两侧使用同一修正，保持局部守恒和方向；
- smooth-extremum 使用二阶 allowance；物理出流只进入单边 budget；入口/回流保持 typed
  frozen boundary；IBM inactive faces 使用统一几何 activity；
- report 记录全局 retained-correction L1 ratio 和唯一受限面数。局部 `alpha=0` 不再被误报
  为整步全域高阶坍缩；只有全局 retained ratio 为零才代表修正全部丢失。

Authority 同时绑定 equation plan、velocity field、geometry、boundary revision、BDF-dt、
density 的完整 view identity，以及三轴 provisional flux 的 base/extents/stride/axis/storage/
revision-domain。相同 revision/storage/domain、不同 base/replica 的 foreign view mutation 被
拒绝。所有 rank-local preflight、boundary resolve 和 RHS 数值检查在 FGMRES/状态修改前
consensus；rank0-only authority 和 non-finite RHS mutation 均为零 solve、零状态修改。

## 5. CFL 与 V5 真实性

Momentum predictor 的 provisional flux 和 terminal pressure-corrected committed flux 分别
计算 outward/absolute convective CFL，绑定不同 flux revision 和同一配置 limit。真实
ProductDriver fixture 证明：

- fixed full step：provisional CFL 通过，五个 terminal 分量已完成，committed CFL 在 stage 60
  超限后返回 `invalid_case/10211`，提交状态、最终通量、压力参考和 closed mass 精确回滚；
- adaptive full step：第一次在同一 stage 60 返回 `rejected_step/10211`，缩到 half step 后
  接受；最终公开物理状态和 final mass flux 与 direct half control 逐位相同。

`HUNDUN_V04_EVIDENCE_V5` 记录 AFC scheme/retained ratio/limited faces、advective CFL 和
committed CFL。validator 还要求单文件 schema/run identity 不变、step 连续、time 严格递增、
`advective.dt == delta(committed time)`、momentum plan/IBM activity/CFL limit 静态，并要求
STL 与 activity collective 的零/非零语义一致。旧 V1–V4 oracle 保持冻结，V5 字段污染旧
schema 会被拒绝。

## 6. BDF–dt 与诊断真实性

BE 与 variable-BDF2 系数必须和 `dt` 一致。continuity、momentum、enthalpy、species、
AssemblyEpoch、direct pressure-storage 和 PISO terminal 均在写数值状态或签发 certificate 前
验证该关系。故不能再用一个时间步的 BDF 矩阵、另一个时间步的 CFL/证书通过验收。

CLI 只有在 terminal audit 确实完成时输出 EOS、continuity、energy、closed mass 和 gauge；
stage 44 等审计前失败输出 `terminal_audit=unavailable`，不输出默认零残差。
continuity witness 还要求 `valid=true`。真实 CLI fixture 保留失败阶段和压力求解路径断言。

## 7. 定向验证

最终 focused 集合为 12/12：

- `v04_evidence_v4_product_validator`
- `v04_app_cli_continuity_witness`
- `v04_solver_pressure_storage`
- `v04_solver_pressure_energy_schur`
- `v04_solver_pressure_energy_globalization_mpi_2`
- `v04_solver_piso_authority`
- `v04_solver_equations_mpi_1/2`
- `v04_io_product_path`
- `v04_app_driver`
- `v04_product_pressure_energy_retry_mpi_1/2`

另通过：

- `v04_core_product_freeze_mpi_test --mass-flow-only`，2 ranks；
- 当前 AFC v2/V5 契约直接修改的
  `v04_product_pressure_energy_temporal_convergence` 与
  `v04_restart_rank_change_continue`，定向复核 2/2；
- `python3 tools/v04_evidence_validate.py self-test`；
- tests-off 增量构建；
- 独立最终 P0/P1 只读审计，结论为无剩余确定 P0/P1。

定向集合第一次运行时，旧 PISO test fixture 用 `dt=1/a0` 构造 BDF2，按新契约被正确
拒绝；fixture 改为项目既有 `time_step_for_bdf()` 后最终 12/12 通过。生产校验没有放宽。
只有未与本轮 AFC/CFL/authority 契约重叠的既有 MMS、PISO temporal order 和
pressure–enthalpy/refinement 证据直接复用。被本轮改写的 temporal-convergence 与 Restart
检查使用当前版本重新运行；没有扩大为完整 CTest。

## 8. 最终 no-IBM control

工件：
`/home/wyf/code_dev/.benchmarks/v04-2-afc2-noibm-accepted-ETZECl/run/evidence.jsonl`

SHA-256：
`91128e72607444270e93e0d98d8f2537d9f9b2305ef7f580ad248e6554c69891`

64-rank 两步运行完成到 `t=0.012`，step 1 为 BE、step 2 为 BDF2，无 retry/fallback。
V5 validator 通过。

| step | C residual | E residual | retained L1 / limited faces | provisional / committed Co_out | max-rank step |
|---:|---:|---:|---:|---:|---:|
| 1 | `2.4932420565e-15` | `2.3063283973e-15` | `1 / 0` | `0.1440000000 / 0.1440000000` | `5.930 s` |
| 2 | `1.2330828250e-15` | `1.0130550697e-15` | `1 / 0` | `0.1440000000 / 0.1440000000` | `5.795 s` |

两步 EOS/mass/gauge 均为零，C/E 低于 `1e-6`；无 IBM activity，双 CFL 低于 `0.8`。

## 9. 最终 64-rank IBM Re3900

case：
`/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/c1-pressure-enthalpy-full-half-20260830-01/case-full`

配置：`480x480x48`，IBM cylinder，`Re_D=3900`，fixed `dt U/D=0.006`，两步。

工件：
`/home/wyf/code_dev/.benchmarks/v04-2-afc2-re3900-accepted-SmS6Jy/run/evidence.jsonl`

SHA-256：
`be1f744b83a3fb33fd9446b67182c271ab9027228b6e45c05f3e7387f90ced58`

V5 validator 通过，运行 `COMPLETED steps=2 time=0.012`。两步均无 retry、无 temporal
fallback，普通 pressure solves 为 2，refinement 为 0。ProductDriver positivity/finite gate
与五个 terminal 分量门均通过；V5 不序列化单元最小值，因此回执不虚构 p/T/rho minima。

| step | 时间层 | C residual | E residual | retained L1 | limited faces | Co_out provisional / committed | linear iterations | max-rank step |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | BE startup | `6.7110454556e-7` | `4.9656607156e-7` | `0.9684178508` | `70,436` | `0.3790477237 / 0.3797974726` | 87 | `30.053 s` |
| 2 | true BDF2 | `2.6262545883e-7` | `2.1919163778e-7` | `0.3649577176` | `25,865,378` | `0.3805492212 / 0.3830326558` | 72 | `26.010 s` |

两步 EOS、closed mass、gauge 均为零，所有 C/E 残差低于 `1e-6`，双 CFL 均低于 `0.8`。
step 2 虽广泛限制高阶面修正，但全局仍保留约 36.5% 的 L1 correction，原来的全域
`theta=0` 数值质量问题已经关闭。实际 maximum convective CFL 已被 V5 明确记录。

## 10. DCO 与正式发布剩余门

数值实现提交具有作者匹配的 `Signed-off-by`。完成本回执后再次审计 `origin/main..HEAD`
的作者/trailer，而不是以“最后一个提交有 sign-off”替代整条可融合链审查。

本轮按约束没有执行完整 CTest、COAST 重算、paired full20、长程 Re3900 或长统计。因此
以下仍阻止 `HUNDUN-FLOW v0.4 Re=3900 JOINT RELEASE GATE`：

1. 新 exact candidate 和 COAST equivalence reseal；
2. 正式 `focused -> full2 -> frozen -> full20 -> literature -> final` machine sequence；
3. 至少五组 alternating paired full20；
4. Re3900 total mean drag 与 periodic `pi D` `Cl_rms` 的完整一手 authority；
5. 至少 `150D/U` 发展和约 `2020D/U` 统计，以及可 Restart accumulator；
6. 长程双 CFL、AFC retained ratio、五门、正性、retry/变步和守恒漂移监测；
7. 最终 provenance、性能、物理统计和发布回执。

所以本节点的准确结论是 `ACCEPT / V04-2 NUMERICAL PRODUCT CLOSURE`，不是 v0.4 正式
发布接受。

# HUNDUN-FLOW V04-2 压力—焓产品闭环验收回执

日期：2026-08-30；最终复核：2026-08-31

裁决：`ACCEPT / V04-2 NUMERICAL PRODUCT CLOSURE`

## 1. 裁决范围

本回执接受 V04-2 数值产品节点：同目标时间层 pressure–enthalpy 闭环、真实
refinement/rollback、MPI 唯一面 AFC authority、完整的双 CFL certificate 和 Evidence V6
候选/时间真实性。

它不是 HUNDUN-FLOW v0.4 正式发布接受。本轮没有运行 full20、COAST 重算、长程
Re3900/物理统计、完整文献 authority 或正式 release gate。

## 2. 不可变运行候选

- worktree：`/home/wyf/code_dev/hundun-flow-v04-2-spec-closure`
- branch：`codex/v04-2-spec-closure-20260831`
- code commit：`2bc6fd5f89732dace0e7cd016ce4bd5dfdb64016`
- code tree：`ff6eaf7668565df117e44b93c210519f8075e24a`
- parent：`e483237b707ee46943f0276558e4a753d9ff0eab`
- build manifest SHA-256：`f367db108667550c718fc0f8bb59143be4fe6874c37ea95c31769387511c13ec`
- tests-off executable：
  `/home/wyf/code_dev/hundun-flow-v04-2-spec-closure/build/v04-2-spec-tests-off-final/versions/v0.4/hundun`
- executable SHA-256：`13933c7d19c3659ef7e5aabab512f452e96bc6e4637e38f6b0a6bdb283831784`
- runtime candidate identity SHA-256：`2bebb5dee34375ffc55690230c7c02e4cc861d59032c96604f55fa6350fbc7c5`
- build：Release、Clang/libc++、MPI，`HUNDUN_BUILD_TESTS=OFF`
- launcher：64 ranks，`ppr:32:socket:PE=1`、slot rank、core bind，三个线程变量均为 1

Evidence V6 每一行实际记录相同的 commit、tree、build manifest、executable 和组合 identity。
后续纯文档提交不改变上述运行候选，也不被表述为生成该二进制的 code commit。

## 3. Pressure–enthalpy 主闭环

真实产品路径为：

```text
accepted same-target state
  -> exact continuity/energy replay
  -> joint-L2 candidate selection
  -> atomic p/h/rho/T/U/phi publication
  -> thermo + boundary/IBM + mass-flux refresh
  -> fresh continuity-energy Jacobian/Schur solve
  -> independent EOS/C/E/mass/gauge and positivity audit
  -> collective commit or exact rollback
```

联合 L2 只选择下降候选，continuity 与 energy 仍分别过门。真实小型 ProductDriver fixture
实际执行两轮 refinement，验证 ordinal、target、lineage 和每轮 state/thermo/boundary/IBM/
flux/Jacobian 刷新；capacity=6 失败 fixture 验证提交状态、历史和通量逐位回滚。最终 Re3900
正常路径不需要 refinement，不能替代也没有被冒充为 refinement 证据。

CLI 仅在 terminal audit 已执行时输出五个真实残差；audit 前 stage-44 类失败输出
`terminal_audit=unavailable`。continuity witness 还要求 `valid=true`。

## 4. 四个重新打开的规格问题

### 4.1 MPI partition-face AFC

根因是旧代码交换 cell ratio 后仍由两侧分别计算 shared-face alpha。RED 将一个非平凡受限面
放在 MPI partition，并污染正侧 consumer 的本地 alpha。新 `common_face_afc_v3_owner` 由负侧
唯一 owner 发布最终三分量共享 alpha；两侧使用同值并向相邻控制体施加等量反号修正。

1/2/4-rank GREEN oracle 均得到：

- shared alpha `min == max`；
- retained L1 ratio `23/53`；
- active faces `18`，limited faces `1`，fraction `1/18`；
- 全局 RHS jump 精确为 `{0x1.ffffffffffff7p-1, -0x1.ffffffffffff8p-1}`；
- 全局守恒为零，且 active/limited/retained/min 只由 owner 计数。

### 4.2 Evidence V6 候选身份与首条时间

旧 V5 的 build/binary 常量不能区分不同候选，且首行没有可信前态锚点。V6 绑定 exact
HEAD/tree、canonical build manifest、运行中 `/proc/self/exe` SHA-256 和组合 identity；MPI
逐字节 consensus 后才写证据。构建依赖覆盖同一 executable-input pathspec，提交或相关源码
改变会自动重新 configure，不能保留 stale identity。

每行新增 `previous_committed_time` 和 `run_start`。fresh 由零步零时刻锚定；restart/
continuation 必须向 validator 提供实际 `manifest.bin`，校验长度、FNV、SHA-256、step/time。
以下 mutation 均 RED，修复后均 fail closed：

- 首行 `time=0.006, dt=0.003`；
- 同一文件混入不同 build/binary identity；
- 只平移首行时间；
- 同时平移 `run_start.previous_time`、`previous_committed_time` 和当前 time；
- restart 缺少、错误或畸形 sidecar；fresh 错传 sidecar。

合法 fresh 与合法 restart/rank-change continuation 均通过。旧 V1–V5 oracle 保持冻结。

### 4.3 committed CFL 完整证书

provisional 与 committed CFL 现在共用 `evaluate_cell_convective_cfl`，消除了公式重复。
committed certificate 绑定 valid、density/final-flux revision、density/flux collective view、IBM
activity、dt、limit，并分别记录 out/abs winner 的 global cell、rank、rhoV、outgoing mass
flow 和 absolute mass-flow sum。超限时 CLI 发布真实 winner；fixed fatal 和 adaptive retry/
缩步/回滚事务语义未改变。

### 4.4 AFC 聚合量

除 retained L1 和 limited faces 外，V6 记录 minimum face alpha、active correction faces 和
limited/active fraction，所有统计均按 MPI 唯一物理面 owner 聚合。无活动 correction 时明确
标记 `not_applicable`，retained/min/fraction 为 JSON null，而不是伪造为 0。

## 5. 定向测试

最终 code commit 上通过的 focused 集合为 14/14：

- `v04_evidence_v4_product_validator`；
- `v04_solver_equations_mpi_1/2/4`；
- `v04_core_product_freeze_mpi_1/2/4`；
- `v04_io_restart_mpi_4`；
- `v04_io_product_path`；
- `v04_app_driver`；
- `v04_restart_rank_change_continue`；
- `v04_product_pressure_energy_retry_mpi_1/2/4`。

另通过：

- `v04_product_pressure_energy_temporal_convergence`；
- `v04_app_cli_continuity_witness`；
- `python3 tools/v04_evidence_validate.py self-test`；
- Release/Clang/libc++/MPI tests-off configure/build。

没有运行完整 CTest。独立 Standards 最终复审为 0 findings；独立 Spec 复审识别出旧 no-IBM
工件不是 fixed `dt=0.006` 以及封存文档/命令仍指向旧候选。前者用最终二进制补跑一次，后者
由本回执和精确 artifact manifest 关闭；其余主动 validator/AFC 反例均正确拒绝或通过。

## 6. 64-rank no-IBM fixed full2

- root：`/home/wyf/code_dev/.benchmarks/v04-2-final-noibm-fixed-ZQnEQkc2`
- case SHA-256：`70a14e28529376c0b005957f0335fc85befd63d7d981cb9aa8a52e4d8bedc89b`
- thermophysics SHA-256：`85949346d7b3095d6ebab81a2e439a2bfe0f3ed78d21539b6c85bc71243dc910`
- evidence SHA-256：`860c550747e2f161a2313a1e764d1750cbc95c97be9422f086cfe9aeefc8f24f`
- log SHA-256：`0b3358e63c364a29a3dc3d336ebce5aab760c02b7f58a4fbddd5997ae14acd9d`
- restart manifest SHA-256：`f6d5bf06400b1fd01c7114481780c90a34bd92ac4f9f5d0b40311c35c17dd42c`

固定 `dt=0.006`，完成两步到 `t=0.012`，BE→true BDF2，无 retry/fallback。V6 runtime
validator 通过。第二步 EOS/mass/gauge 为零，continuity `1.2330828250412647e-15`，energy
`1.0130550696799174e-15`；provisional out/abs CFL 为
`0.1440000000000042/0.1440000000000049`，committed 为
`0.14400000000000424/0.14400000000000504`，均低于 `0.8`。无活动 correction，AFC 正确记录
`not_applicable`、active/limited 均为 0，而非伪造 retained ratio。

## 7. 64-rank IBM Re3900 fixed full2

- case：`/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/c1-pressure-enthalpy-full-half-20260830-01/case-full`
- root：`/home/wyf/code_dev/.benchmarks/v04-2-final-re3900-VcRKDg3f`
- case SHA-256：`ad2998ef0438d184c82d111697fdca477eb113a4f62efe60054070716ac86cf1`
- thermophysics SHA-256：`85949346d7b3095d6ebab81a2e439a2bfe0f3ed78d21539b6c85bc71243dc910`
- STL SHA-256：`bd264c586543de4ec330f53cb2d3d9dfba550823db69451ac3176c603d248f46`
- evidence SHA-256：`ec05c184735779d9497ce1d2d5075141fcb21e613f5a734e293959f50dba3f38`
- log SHA-256：`6d712ebf2173805b943c5907612fe7baa18601b61489987bca38877d7b91858c`
- restart manifest SHA-256：`dd52a252d2a51e781e2931aa8ad95561da5e51317f269b4612cb422450f75e05`

`480x480x48`、IBM cylinder、fixed `dt U/D=0.006`，两步到 `t=0.012`。V6 validator
通过；BE→true BDF2，无 retry/fallback，普通 pressure solves 每步为 2，refinement 为 0。
ProductDriver 的有限性/严格正性门和五个 terminal 分量门均通过；V6 没有序列化 p/T/rho
最小值，因此本回执不虚构这些数值。

| step | C residual | E residual | provisional out/abs CFL | committed out/abs CFL | retained L1 | min alpha | active / limited / fraction | linear iters | max-rank step |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | `6.7110454556e-7` | `4.9656607156e-7` | `0.3790477237 / 0.3790477237` | `0.3797974726 / 0.3797931894` | `0.9684178508` | `0` | `32,798,206 / 70,436 / 0.00214756` | 87 | `27.904 s` |
| 2 | `2.6262545874e-7` | `2.1919163712e-7` | `0.3805492212 / 0.3805406548` | `0.3830326558 / 0.3830353109` | `0.3649577176` | `0` | `32,834,706 / 25,865,378 / 0.78774508` | 72 | `24.622 s` |

两步 EOS/closed mass/gauge 为零，C/E 均低于 `1e-6`，双 CFL 均低于 `0.8`。step 2
存在大量局部 `alpha=0` 面，但全局仍保留 36.5% 的 L1 correction；这关闭了原全域
`theta=0` 坍缩问题，不等于证明长程 limiter 质量。

## 8. 封存 manifest、DCO 与最终裁决

- artifact manifest：
  `/home/wyf/code_dev/.benchmarks/v04-2-final-re3900-VcRKDg3f/closure-manifest.json`
- manifest SHA-256：`718110b21db184439c9fe6a7ee38fd4c9fab58e7620deb8b0bfd466efadbd0e2`

manifest 保存实际 binary/case/output 的绝对路径、完整 mpirun 命令和上述 artifact 哈希，不含
`<binary>`、`<case-full>` 或 `<root>` 占位符。代码提交含作者匹配的 `Signed-off-by`；封存
文档也使用 signed commit。原始含其他任务修改的工作树未 reset、stash 或覆盖。

最终裁决为 `ACCEPT / V04-2 NUMERICAL PRODUCT CLOSURE`。仍未完成的正式发布门包括新
machine candidate 与 COAST equivalence reseal、规定顺序的 full2/frozen/full20/literature/
final、至少五组 paired full20、长程稳定性/物理统计、完整一手文献 authority、可 Restart
统计 accumulator 和最终 release provenance。因此不得把本回执表述为 v0.4 正式发布接受。

# 严格 Re=3900 COAST 与普通可压缩 COAST 差异及 HUNDUN 优化含义

日期：2026-09-02
状态：只读研究记录；未修改 COAST、HUNDUN 或任何算例
用途：作为后续 HUNDUN-FLOW 性能优化计划的前置证据，不作为新的性能接受结论

## 1. 结论摘要

1. **[verified] 两者不是两套独立程序。** 普通可压缩 profiling 分支
   `84ea03dadaa67a446cbe60fdfe94c4cc471dc7bb` 直接建立在严格 COAST
   `3c22e0f029db1b2ca045ec9e212a95eacbcfe6a3` 上；普通运行通过
   `re3900_equivalence_enabled=.false.` 退出严格等价性路径。普通分支新增的主要是
   attribution 计时和固定 `dt=0.006` 的短测控制，不是另一套压力算法实现。
2. **[verified] 严格版的主导慢路径是压力 ICCG，不是终端审计。** 两步 attribution
   中，两个 SIMPLE pressure route 约 `33.23 s`，其中 ICCG inclusive 约
   `31.26 s`，占 COAST twin 总耗时约 `77%`；压力本地工作和 collective wait
   联合约 `32.50 s`。终端审计两次 twin 的两步总计仅 `0.146834 s` 和
   `0.130936 s`，平均约 `0.138885 s`，占 `40.646474 s` 的约 `0.342%`。
3. **[verified] 严格版同时锁死了三项互相放大的条件：** 压力后端
   `legacy_iccg`、禁止 fallback、所有主方程内层容差 `1e-6`/上限 500；压力真实残差
   每 128 次 ICCG 才复核一次。源码只允许在这些 128 次边界确认真实收敛，保存的
   full20 证据因此出现 `38 × 256 + 2 × 384` 次压力迭代。
4. **[verified] “严格版把 COAST 全部改成 FP64”不成立。** 两种模式共享 FP32
   主场、系数和工作数组；严格版增加的是 FP64 残差累加、FP64 `MPI_Allreduce`、
   压力增量补偿以及终端 EOS/continuity/gauge 审计。
5. **[verified] 固定两次外迭代不是两者差异。** 两边输入均为 `itnum=2`；严格版
   额外校验每步必须恰好完成两次 SIMPLE、每个主方程必须恰好求解两次，并在通过
   终端审计后才提交时间步。普通模式同样执行两次 SIMPLE，但严格终端函数立即返回。
6. **[verified] 普通短跑不是“已锁定 GMG”的稳态 auto 基线。** auto 需要 20 个
   pressure samples 才锁定；现有普通运行只有 3 步 × 2 次 SIMPLE，即 6 个压力样本，
   且没有 `pressure_solver_choice.dat`。它处于 legacy/GMG/Hypre 候选轮转阶段；Hypre
   不可用时回退 GMG。`0.453761 s` 应描述为“普通 auto 训练/回退阶段短测”，不能写成
   “稳定锁定 GMG 后的长期基线”。
7. **[inference] 对 HUNDUN 的直接含义是保留 FGMRES/F-cycle MG，借鉴普通 COAST
   的内外精度分层，而不是移植 legacy ICCG。** HUNDUN 应继续以 FP64 外层真实残差和
   终端物理门为权威，把早期内层容差放宽到 COAST 量级或采用 nonlinear forcing；
   现有配置 A/B 已证明 `1e-4` 能减少约 24% 单步耗时。
8. **[verified, supporting evidence] 已有同尺度但非同二进制的对照。** 严格 r10 与普通
   profile 均为 320×320×33、128 ranks；严格 BDF2 median 为 `1.192407 s`，普通为
   `0.453761 s`，前者约为后者 `2.6278×`。这支持“求解器/严格契约差异值得优先拆分”，
   但二进制和运行条件并未完全冻结一致，不能据此分配因果份额。
9. **[unknown] 目前还没有同一二进制、同一 320×320×33 case、逐项只切换
   backend/tolerance/strict-mode 的因子化 A/B。** 因而可以确认慢路径和工作膨胀机制，
   但不能把总差值精确拆成“ICCG、`1e-6`、FP64 audit、BDF/EOS seam”各占多少。

## 2. 身份冻结

### 2.1 严格 Re=3900 COAST

- 工作树：`/home/wyf/code_dev/Coast_software_worktrees/re3900-equivalence`
- 分支：`codex/re3900-equivalence-candidate`
- HEAD：`3c22e0f029db1b2ca045ec9e212a95eacbcfe6a3`
- tree：`ba449790e4918e7a3fd5c21e71c4e9f980a4691f`
- 冻结二进制 SHA-256：
  `c5e483a0442b430e58040a8e04f5adf42b70b3595319c643e9265c87341a9bc4`
- 身份证据：
  `/home/wyf/code_dev/.benchmarks/pa/runs/c0/candidate.identity:1-3`
- 运行命令：
  `/home/wyf/code_dev/.benchmarks/pa/runs/c0/command.txt:1`

严格 seam 提交链：

| 提交 | 时间 | 作用 |
|---|---|---|
| `2ae61043d6a9698eca605e47ab6fff92544f95d2` | 2026-08-21 16:19 +08:00 | 增加可审计 Re3900 equivalence mode |
| `725c90a99a13e8dbfbce6c4cec76059ca11ccf78` | 2026-08-21 19:05 +08:00 | 封存可靠残差、物理压力算子和瞬态 seam |
| `1075db10067bcbe2de0885f401d8bb50ab33f0d7` | 2026-08-21 19:11 +08:00 | 终端审计成为唯一时间步提交点 |
| `3c22e0f029db1b2ca045ec9e212a95eacbcfe6a3` | 2026-08-21 19:17 +08:00 | 修正 output-off finalization |

### 2.2 普通可压缩 COAST profiling

- 工作树：`/home/wyf/code_dev/cpt`
- 分支：`codex/coast-compressible-module-profile`
- profiling HEAD：`84ea03dadaa67a446cbe60fdfe94c4cc471dc7bb`
- profiling tree：`961a5eb3abf1660b649652fd915c40cc0676961b`
- 记录的 base source HEAD：仍为严格版 `3c22e0f...`
- 二进制 SHA-256：
  `cb6fa1297c78d618fb01b04f73421abec538877d44d626f4815a28b9f427d1f6`
- 身份证据：
  `/home/wyf/code_dev/.benchmarks/coast-hundun-compressible-modules-20260902/coast-candidate-seal.json:1-13`

`git merge-base codex/re3900-equivalence-candidate codex/coast-compressible-module-profile`
返回 `3c22e0f...`。因此本研究把差异分成两类：

- **运行时差异：** strict on/off、后端、fallback、容差和算例尺寸；
- **profiling overlay：** attribution 计时和短测时间步控制，不解释数值算法差异。

## 3. 配置差异

| 项目 | 严格 Re3900 | 普通可压缩 | 状态与含义 |
|---|---|---|---|
| `compressible` | `true` | `true` | **[verified] 相同**，不可用“普通版不可压缩”解释速度 |
| SIMPLE 次数 | 2 | 2 | **[verified] 相同** |
| u/v/w | `500, 1e-6` | `50, 1e-4` | **[verified]** 普通版更松、上限更低 |
| pressure | `500, 1e-6` | `500, 1e-4` | **[verified]** 上限相同，内层容差相差 100 倍 |
| enthalpy | `500, 1e-6` | `20, 1e-4` | **[verified]** 普通版显著更松 |
| pressure backend | `legacy_iccg` | `auto` | **[verified]** |
| fallback | `none` | `coast_gmg_pcg` | **[verified]** |
| strict seam | `.true.` | `.false.` | **[verified]** |
| fixed `dt` | `0.006` | `0.006` | **[verified] 相同** |
| strict terminal tolerance | `1e-6` | 配置仍存在，但 strict off | **[verified]** 普通模式不执行该审计 |
| 主场精度 | FP32 | FP32 | **[verified] 相同** |

一手配置：

- 严格 full20：
  `/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/candidates/h2c234b1-c3c22e0f-r1/coast_case/input.full20.d:48-53,95,102-111`
- 普通短测：
  `/home/wyf/code_dev/.benchmarks/coast-hundun-compressible-modules-20260902/runs/coast-ordinary-c2/case/input.d:30,48-53,95,102-111`

两个 input 的直接 diff 只在本研究关心的部分显示：步数/终止时间、五个方程的
cap/tolerance、pressure backend/fallback 和 strict enable 标志不同。

## 4. 压力求解路径

### 4.1 严格版：强制 legacy ICCG 和真实残差分段

**[verified]** 严格初始化会 fail-closed 检查：可压缩、无 restart/PDF/燃烧/喷雾、
固定材料和 Vreman 参数、恰好两次 SIMPLE、`legacy_iccg` 且无 fallback。任一条件不符
立即停止：

- `/home/wyf/code_dev/Coast_software_worktrees/re3900-equivalence/SRC.Coast/module_coast_re3900_equivalence.F90:195-230`

同一模块冻结：

- `re3900_linear_cap=500`；
- `re3900_internal_tolerance=1e-6`；
- FP64 `atol=1e-8`、`rtol=1e-6`；
- pressure true-residual interval = 128；
- BiCGStab true-residual interval = 8。

来源：`module_coast_re3900_equivalence.F90:119-128`。

压力 ICCG 的控制循环是：

1. 每次迭代做不完全 Cholesky 前后扫；
2. 一次 search-direction halo；
3. 至少三次显式 `MPI_Allreduce`（`sum1`、`sum2`、`err`）；
4. strict 模式下只有 `ncycl % 128 == 0` 时计算独立真实残差；
5. 若真实残差未通过，清空搜索方向并从可靠残差继续。

来源：
`/home/wyf/code_dev/Coast_software_worktrees/re3900-equivalence/SRC.Coast/cgsol.F90:187-351`。

这意味着严格版并非“每次迭代都能在真实残差达到阈值时退出”，而是将可能的接受点
量化到 128 的倍数。full20 pair-01 的 40 个压力求解中：

- 38 个为 256 iterations / 4 true-residual audits；
- 2 个为 384 iterations / 5 true-residual audits；
- 压力总迭代数 10,496，平均 262.4 次/solve。

原始证据首步与末步：

- `/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/f20_pibm_succ_p01_HC_fldtZG/C/re3900_evidence.jsonl:1-12`
- 同文件 `:210-220`

因此“128 审计节拍 + ICCG + 1e-6”是严格版工作量膨胀的**已验证机制**。但没有逐因子
A/B，尚不能把 256/384 次迭代分别归因到后端、容差、审计节拍或严格瞬态算子的
条件数。

### 4.2 普通版：auto 候选轮转和 GMG fallback

普通配置默认 auto，候选顺序为：

1. `legacy_iccg`；
2. `coast_gmg_pcg`；
3. Hypre Struct PFMG；
4. Hypre Struct SMG。

需要 20 个 pressure samples 才选择并锁定最快稳定后端：

- 候选及样本数：
  `/home/wyf/code_dev/cpt/SRC.Coast/module_coast_pressure_solver.F90:64-74`
- 轮转、记录和 20 次后锁定：同文件 `:183-250`
- 失败 fallback：同文件 `:87-167`

现有 ordinary c1/c2 每个只有 3 步、每步 2 次 SIMPLE，因此只有 6 个 pressure
samples，且运行目录没有 `monitor/pressure_solver_choice.dat`。screen 证据表明：

- 配置为 auto + GMG fallback：`.../coast-ordinary-c2/case/screen:113`；
- Hypre PFMG/SMG 返回 `status=-100` 并 fallback：同文件 `:163-164`；
- 每步最后记录的后端为 `coast_gmg_pcg`，迭代 7、4、0：
  `.../coast-ordinary-c2/case/monitor/pressure_solver.dat:1-6`。

**[verified limitation]** `pressure_solver.dat` 只保留每步最后一次 pressure solve，不能用
三行 GMG 证明 auto 已锁定 GMG。按源码确定性轮转，六次请求仍处在候选训练阶段。

### 4.3 普通 GMG 的实现特点

`coast_gmg_pcg` 是 PCG 外层加 block-local 两层结构 V-cycle：

- PCG：`module_coast_pressure_solver.F90:342-415`；
- operator halo：同文件 `:586-612`；
- FP64 dot accumulation + FP64 reduction：同文件 `:644-660`；
- 3 次 pre-Jacobi、3 次 post-Jacobi、局部粗网格：同文件 `:664-691`。

它仍会在每个 V-cycle 重复 allocation，并在每次 Jacobi operator 中做 halo，不能被当成
HUNDUN 的目标实现。但由于 ordinary 短测只需 0–7 次 PCG 迭代，它的工作量远低于
严格 ICCG 的 256/384 次。

## 5. 数据精度与真实残差

### 5.1 主场没有全局升为 FP64

**[verified]** 严格与普通模式共享相同数组模块：`field`、`coef`、`work`、`p`、`rho`、
`f` 等均声明为默认 `real`；driver 对主要场的显式 view 是 `real(kind=4)`：

- `/home/wyf/code_dev/Coast_software_worktrees/re3900-equivalence/SRC.Coast/module_arrays.F90:22-63`
- `/home/wyf/code_dev/Coast_software_worktrees/re3900-equivalence/SRC.Coast/app/coast_legacy_driver.F90:111-114`

`git diff 3c22e0f..84ea03d -- SRC.Coast/module_arrays.F90` 为空。因此普通 profile 没有另行
降低数据精度，严格版也没有全场 FP64 化。

### 5.2 严格版增加的 FP64 层

**[verified]** 严格版在每个线性 solve 前后独立形成 `b-Ax`，使用 double accumulation 和
`MPI_DOUBLE_PRECISION` reduction，并以 `max(1e-8, 1e-6*||b||)` 判定：

- solve begin/end：`module_coast_re3900_equivalence.F90:429-513`；
- pressure scaled reliable residual：同文件 `:515-604`；
- pressure compensated update：同文件 `:606-633`；
- 通用 reliable/true residual：同文件 `:661-785`。

full20 pair-01 每步 10 个线性求解；20 步共 200 solves、731 次 FP64 true-residual audits，
平均 36.55 次 audit/step。压力自身是 40 solves、162 audits。

**[inference]** 这些 audit 会增加全场访问、halo 和 reduction，但其独立耗时未被单独
计时，不能声称它占某个确定百分比。就压力路线而言，10,496 次 ICCG 迭代远多于
162 次真实残差复核，当前证据更支持“审计节拍导致迭代量化并暴露 ICCG 弱收敛”，
而不是“FP64 算术本身占了大多数时间”。

## 6. 外迭代、终端审计和提交语义

### 6.1 相同点

- 两边均 `itnum=2`；
- driver 均执行 `do iter=1,niter`；
- 每个 SIMPLE 包含 u/v/w/h/p 五个线性 solve。

严格 full20 证据每步固定 `simple_iterations=2`、`linear_solve_count=10`。

### 6.2 严格版额外语义

**[verified]** 严格 terminal audit：

- 校验五个方程的 solve inventory 都恰好为 2；
- 以 FP64 扫描流体单元的 EOS 和 BDF continuity；
- 检查压力出口 gauge、状态合法性和 outlet inventory；
- success path 做 6 次集合归约；
- 只有所有门通过后才更新历史场、`istep` 和 `tim`；
- 再写入 `re3900_evidence.jsonl`。

来源：
`/home/wyf/code_dev/Coast_software_worktrees/re3900-equivalence/SRC.Coast/module_coast_re3900_equivalence.F90:920-1214`。

普通模式调用同一函数，但在 `re3900_equivalence_enabled=.false.` 时于该函数首部返回。

严格 twin 的 `terminal_audit|terminal` 两次运行总计分别为：

- `0.146834355982719 s`：`/home/wyf/code_dev/.benchmarks/pa/summary.json:172-203`；
- `0.130936233981629 s`：同文件 `:375-406`。

**[verified]** 终端审计是重要的正确性差异，但它不是严格版十几秒慢问题的主因。

### 6.3 严格瞬态/EOS seam

严格模式还执行：

- BE→BDF2 的密度历史和 density defect；
- BDF2 storage rate/pressure reciprocal time；
- 固定材料/EOS/Vreman 参数；
- 每个压力 solve 的物理 operator snapshot/restore；
- FP32 主压力加 FP64 计算的增量补偿；
- accepted-state 历史场复制。

源码：

- `module_coast_re3900_equivalence.F90:330-408`；
- `app/coast_legacy_driver.F90:1905-1923,2288-2330`。

**[unknown]** 这些逻辑中的部分会增加数组流量，也可能改变压力算子条件数，但没有与
ordinary 同网格的逐项 A/B，不能给出独立性能份额。

## 7. MPI 与数组访问对比

| 路径 | 主要局部工作 | 主要同步 | 结论 |
|---|---|---|---|
| strict legacy ICCG | 每 iter 三角 preconditioner、operator、更新和残差全场 pass | 每 iter 至少 3 个显式 Allreduce + 1 个 halo；每 128 iter 再做 reliable audit | **[verified] 高迭代数放大本地带宽和同步** |
| ordinary GMG-PCG | 每 iter operator + block-local V-cycle（6 次 Jacobi operator） | dot/scaled-max reductions；operator/Jacobi halos | **[verified] 单 iter 不便宜，但只需少量迭代** |
| strict FP64 audit | 独立 `b-Ax` 全场 pass | halo + FP64 Allreduce | **[verified] 额外工作；单独耗时 unknown** |
| strict terminal | EOS/continuity/gauge 全场 pass | success path 6 reductions | **[verified] 约 0.34%，非主因** |

冻结 attribution 对严格 COAST 的总结：两个压力 route 约 `33.23 s`，ICCG inclusive
约 `31.26 s`；pressure local `20.81 s` + collective `11.69 s`，联合 `32.50 s`；
纯 halo 约 `1.12 s`、仅 `2.76%`：

- `/home/wyf/code_dev/.benchmarks/pa/report.md:111-115`

因此不能把 strict COAST 慢归咎于“halo payload 太大”；它与 HUNDUN 的相似之处恰恰是
压力路线的本地数组工作和细粒度同步共同主导。

## 8. 性能证据及不可直接相除的口径

### 8.1 严格 COAST

冻结两步 480×480×48、64 ranks、172,800 cells/rank：

- BE `21.927395 s`；
- BDF2 `17.431275 s`；
- 两步 `39.358670 s`。

来源：`/home/wyf/code_dev/.benchmarks/pa/report.md:37-52`。

更稳定的 5-pair full20 统计给出：

- strict COAST P90 `15.419434811 s/step`；
- 当时 HUNDUN `2c234b1` P90 `4.418378075 s/step`；
- 5 个配对的 HUNDUN/COAST median ratio `0.29128`。

来源：
`/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/formal_full20_h2c234b1_c3c22e0f_r1_03/PAIRED-05-STATS.json:1-19`。

因此“严格 COAST 与 HUNDUN 性能相近”只对某些后续 HUNDUN snapshot/两步口径成立；
它不是跨版本事实。后续 core-hotpath HUNDUN 两步为 `39.096166 s`，与冻结严格 COAST
`39.358670 s` 接近；早期 full20 HUNDUN 则明显更快。来源：
`/home/wyf/code_dev/.benchmarks/pa/report.md:3-12`。

### 8.2 普通可压缩 COAST

普通短测为 320×320×33、128 ranks、26,400 cells/rank，BDF2 median
`0.453761 s`；其配置为 compressible、auto、`1e-4`，计时窗口排除输出：

- `/home/wyf/code_dev/.benchmarks/coast-hundun-compressible-modules-20260902/report.md:1-12`
- 同文件 `:29-47` 为模块分解。

### 8.3 同尺度 supporting evidence

**[verified]** `/home/wyf/code_dev/.benchmarks/c3900s/r10` 是已接受的严格模式运行：

- 128 domains、320×320×33；
- `compressible=true`，五个 U/V/W/P/H 方程均为 cap 500、tolerance `1e-6`；
- `legacy_iccg`、无 fallback、`re3900_equivalence_enabled=.true.`；
- BDF2 steps 2–10 median `1.192407 s`；10 步、100 个 linear solves 和 10 个 terminal
  evidence rows 均通过保存的验证；
- 20 个 pressure solves 全部为 256 iterations；每步全部方程合计 34–37 次 FP64
  true-residual audits，除前 3 步外均为 36 次。

一手证据：

- `/home/wyf/code_dev/.benchmarks/c3900s/cause.md:14-28`；
- `/home/wyf/code_dev/.benchmarks/c3900s/r10/mesh_report.txt:5-7`；
- `/home/wyf/code_dev/.benchmarks/c3900s/r10/input.d:48-53,95,102-111`；
- `/home/wyf/code_dev/.benchmarks/c3900s/r10/re3900_evidence.jsonl`；
- `/home/wyf/code_dev/.benchmarks/c3900s/v10.json:1-15`。

普通 profile 同样为 320×320×33、128 ranks，BDF2 median `0.453761 s`；因此严格/普通
时间比约为 `1.192407 / 0.453761 = 2.6278×`。**[inference]** 网格和 rank 数相同使“尺度
差异”不再足以解释该比值，backend、tolerance 与 strict contract 是需要优先因子化的
候选。**[limitation]** 两次运行并非同一二进制，运行日期、计时工具和部分短测控制也不
完全相同，所以这只是 supporting evidence，不是可用于归因的受控 A/B。

### 8.4 比较限制

**[verified]** strict 与 ordinary 现有样本不同：

- 网格：480×480×48 vs 320×320×33；
- ranks：64 vs 128；
- cells/rank：172,800 vs 26,400，相差 6.545 倍；
- strict 是锁定 ICCG 的真实接受路径；ordinary 是未完成 auto 锁定的短测；
- 统计口径和运行日期不同。

**[inference only]** 用单个 BDF2 时间再除 cells/rank，strict 仍约慢 5.87 倍，但这不是
受控因果 A/B，不能用作接受指标。

## 9. 根因分级

### 已验证

1. strict pressure route 占总 twin 约 77%；
2. strict 强制 legacy ICCG/no fallback；
3. strict/ordinary 内层容差分别为 `1e-6`/`1e-4`；
4. strict 压力真实残差仅在 128 的倍数复核，实际压力迭代为 256/384；
5. strict ICCG 每迭代包含多次全场 pass、至少 3 个显式 reduction 和一次 halo；
6. 两边主场均为 FP32；
7. 两边均为两次 SIMPLE；
8. terminal audit 自身约 0.34%，不是主因；
9. ordinary short run 的 auto 尚未锁定。

### 有直接机制支持、但尚未独立计时

1. 每 solve 的 FP64 true-residual begin/end 和 periodic reliable correction；
2. operator snapshot/restore、pressure compensation 和 accepted history copy；
3. strict BDF/EOS storage 对压力算子条件数的影响；
4. strict 大网格和较少 ranks 对 ICCG/collective 扩展性的放大。

### 当前未知

1. legacy ICCG、`1e-6` 和 128-step audit cadence 各自的独立耗时份额；
2. ordinary auto 在完成 20 个相同工作负载样本后最终选择哪个 backend；
3. 同一 480×480×48 网格上 ordinary GMG/`1e-4` 的稳定 BDF2 时间；
4. 保留 strict physics/terminal contract，仅替换内层 backend/tolerance 后的正确性与速度；
5. FP32 preconditioner/FP64 outer residual 对 HUNDUN 的真实收益。

## 10. 对 HUNDUN-FLOW 的优化约束和机会

### 必须保留

1. HUNDUN 的 FGMRES + F-cycle native Cartesian MG 基线，不改成 legacy ICCG。
   当前 FGMRES 已把一个 Arnoldi column 的 basis dots 和 norm 合成一次 reduction：
   `/home/wyf/code_dev/hundun-flow-v04-ibm-adaptive-order/versions/v0.4/src/solver_krylov.cpp:2021-2091`。
2. 外层真实残差、终端 EOS/continuity/energy/gauge 和 rollback 仍使用 FP64 权威。
3. 两次 PISO 和最终候选态审计不变。
4. 当前 F-cycle MG 及非奇异 low-Mach storage 语义不变：
   `/home/wyf/code_dev/hundun-flow-v04-ibm-adaptive-order/versions/v0.4/src/solver_piso.cpp:6736-6767`。

### 优先优化

1. **P0：内外容差分层。** 把压力线性内层从固定 `1e-6` 改为经验证的 `1e-4`
   起步或 nonlinear-residual-aware forcing；最终 FP64 true residual 和终端门仍为
   `1e-6`。现有配置-only A/B 已节省 `24.0%`、压力迭代均值由 51 降到 36.5：
   `/home/wyf/code_dev/.benchmarks/coast-hundun-compressible-modules-20260902/report.md:41-47`。
2. **P0：保持 HUNDUN 的无量化停机，并 A/B 审计频率。** HUNDUN FGMRES 在递归残差
   达到 tolerance 时立即触发 FP64 true-residual audit，同时也在固定 interval 和 cycle
   end 审计；它没有 strict ICCG 那种“只能在 128 的倍数确认退出”的缺陷：
   `/home/wyf/code_dev/hundun-flow-v04-ibm-adaptive-order/versions/v0.4/src/solver_krylov.cpp:2298-2315`。
   可测试 interval 8、16 或自适应策略，但必须保留递归残差达标时的立即审计、cycle-end
   审计、FP64 终端接受条件，并以性能和数值 A/B 同时通过为准。
3. **P1：减少 pressure/MG/candidate 控制归约。** HUNDUN Arnoldi 已单归约，下一目标应是
   MG、candidate、AFC 和 fail-closed halo 的剩余控制同步，而不是重写 Krylov。
4. **P1：融合 pressure-energy/EOS/candidate/provenance 全场 pass。** strict COAST 说明
   低效压力路径会同时放大数组带宽和 collective；HUNDUN 应继续减少 Schur、candidate
   和审计重复遍历。
5. **P2：仅在预条件器内部试混合精度。** COAST 的快路径说明 FP32 工作集有带宽优势，
   但 HUNDUN 当前 FieldView 明确只接受 FP64：
   `/home/wyf/code_dev/hundun-flow-v04-ibm-adaptive-order/versions/v0.4/include/hundun/v04_execution.hpp:116-148`。
   因此只应在独立 MG workspace/coarse/halo 中试 FP32，外层 operator、解、真实残差和
   terminal audit 保持 FP64。

### 不应采取

1. 不把 HUNDUN 换成 strict COAST 的 legacy ICCG；
2. 不把普通 COAST 的 fixed-two SIMPLE 当作 HUNDUN 的最终停止准则；
3. 不为追求速度删除 FP64 true residual 或终端物理门；
4. 不把整个 HUNDUN 状态直接降为 FP32；
5. 不把 ordinary `0.453761 s` 当作同网格的目标值或稳定 locked-GMG 结果。

## 11. 后续因子化 A/B 与主计划关系

执行顺序以
`/home/wyf/code_dev/hundun-flow-v04-ibm-adaptive-order/docs/superpowers/plans/2026-09-02-re3900-performance-optimization.md`
的 R1 为唯一权威；本节只重述该矩阵，避免形成第二套并行计划。R1 先在同一 profiling
binary、同一 320×320×33 case、相同输出设置和 rank mapping 下运行：

| 组合 | Re3900 mode | 压力后端 | 线性容差 | 用途 |
|---|---|---|---|---|
| C0 | off | `auto` + GMG fallback | `1e-4` | 重现普通基线 |
| C1 | off | 强制 `coast_gmg_pcg`，无 fallback | `1e-4` | 去除 auto screening 噪声 |
| C2 | off | 强制 `legacy_iccg`，无 fallback | `1e-4` | 隔离后端差异 |
| C3 | off | 强制 `legacy_iccg`，无 fallback | `1e-6` | 隔离线性容差工作量 |
| C4 | on | 严格 `legacy_iccg`，无 fallback | `1e-6` | 测量严格物理/审计合并增量 |

每组至少记录 1 个 BE + 5 个 BDF2 步，并交替运行三轮。每个节点至少报告：

- backend 与真实迭代序列；
- true-residual audits、operator/preconditioner applies；
- Allreduce/Allgatherv/halo 次数与等待时间；
- pressure local、terminal audit 和整步时间；
- EOS/continuity/gauge；对 HUNDUN 还要报告 energy；
- 无输出计时窗口和不可变 binary/case hash。

只有当 C3 → C4 仍有显著差距时，才进入第二阶段 diagnostic-only patch，把 strict
cadence、high/low precision audit、terminal audit 等进一步拆开；补丁不得改写冻结
oracle。HUNDUN 的 interval 8/16 或自适应审计 A/B 属于主计划 H0，须在 R1 完成并通过
进入门后执行。

## 12. 可复现命令

```bash
# 身份与提交链
git -C /home/wyf/code_dev/Coast_software_worktrees/re3900-equivalence \
  status --short --branch
git -C /home/wyf/code_dev/Coast_software_worktrees/re3900-equivalence \
  log -4 --date=iso-strict --format='%H %cI %s'
git -C /home/wyf/code_dev/cpt merge-base \
  codex/re3900-equivalence-candidate codex/coast-compressible-module-profile

# 配置差异
diff -u \
  /home/wyf/code_dev/.benchmarks/pa/runs/c0/case/input.d \
  /home/wyf/code_dev/.benchmarks/coast-hundun-compressible-modules-20260902/runs/coast-ordinary-c2/case/input.d

# full20 压力迭代/审计频数
rg '"equation":"pressure"' \
  /home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/f20_pibm_succ_p01_HC_fldtZG/C/re3900_evidence.jsonl \
  | rg -o '"iterations":[0-9]+|"true_residual_audits":[0-9]+' \
  | sort | uniq -c

# 所有 solve 的 FP64 audit 总数
rg '"kind":"linear_solve"' \
  /home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/f20_pibm_succ_p01_HC_fldtZG/C/re3900_evidence.jsonl \
  | sed -E 's/^.*"true_residual_audits":([0-9]+).*$/\1/' \
  | awk '{sum += $1; n += 1} END {print n, sum, sum/n}'

# 严格契约静态测试（已 PASS）
cd /home/wyf/code_dev/Coast_software_worktrees/re3900-equivalence
PYTHONDONTWRITEBYTECODE=1 \
  python3 SRC.Coast/tests/check_re3900_equivalence_static.py
```

预期关键输出：pressure `38 × 256`、`2 × 384`；pressure audits `38 × 4`、
`2 × 5`；全部 `200 solves / 731 audits / 3.655 audits-per-solve`。

## 13. 本研究边界

- 未运行新的 MPI/长算例；
- 未修改 COAST、HUNDUN 或算例；
- 未把不同网格的时间比当作接受指标；
- 未声称 ordinary auto 已锁定 GMG；
- 未把 FP64 audit、strict transient seam 或 NUMA 的影响伪装成已独立量化；
- 后续实施应以
  `/home/wyf/code_dev/hundun-flow-v04-ibm-adaptive-order/docs/superpowers/plans/2026-09-02-re3900-performance-optimization.md`
  为执行顺序权威，并逐节点回填本文件列出的 unknown。

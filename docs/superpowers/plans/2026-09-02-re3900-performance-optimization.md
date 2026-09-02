# HUNDUN-FLOW Re3900 性能优化计划

日期：2026-09-02

状态：`COMPLETE / R0--R1 与 H0--H5 已闭环`

## 目标

在不降低最终数值验收精度的前提下，解释并缩小 HUNDUN-FLOW 与普通可压缩
COAST 的单步性能差距。先用严格 Re3900 COAST 的慢化机制建立因果对照，再逐项
优化 HUNDUN；每项只改变一个因素，未通过数值门和性能门不得进入下一项。

配套源码研究：
[严格 Re3900 COAST 与普通可压缩 COAST 差异](../../verification/2026-09-02-coast-strict-re3900-vs-ordinary-compressible.md)。

## 不可变边界

- COAST 与 HUNDUN 均保持 `compressible=true`。
- HUNDUN 保留两次 PISO、FGMRES/Native-MG、FP64 外层状态与真实残差。
- HUNDUN 的 EOS、continuity、energy、closed-mass、gauge、正性、CFL、commit/rollback
  终端验收不放宽；当前产品门仍为 `1e-6`。
- 不通过改变 20D x 10D x 3D 的物理边界、周期性、时间步或网格来获得提速。
- 不把 COAST `legacy_iccg`、固定两次 SIMPLE 或全局 FP32 状态直接迁入 HUNDUN。
- 所有新性能测试关闭计时窗口内输出，固定 MPI rank/线程/绑核和 executable identity。

## 当前证据

### 同尺度短算例

320 x 320 x 32/33、128 MPI ranks、`dt=0.006`：

| 程序/策略 | BDF2 中位耗时 | 说明 |
|---|---:|---|
| 普通可压缩 COAST | 0.453761 s | `auto` 探测阶段，混合 ICCG、GMG 和 Hypre->GMG fallback，线性容差 `1e-4` |
| 严格 Re3900 COAST | 1.192407 s | `legacy_iccg/no fallback`，线性容差 `1e-6`，FP64 真残差和终端审计 |
| HUNDUN 当前策略 | 3.049257 s | 压力相对容差 `1e-6` |
| HUNDUN 配置 A/B | 2.316371 s | 仅把压力相对容差改为 `1e-4`，终端门不变 |

严格 COAST 比普通 COAST 慢约 2.63 倍，但现有两组运行的源码身份和计时仪器不同，
只能证明“存在同量级慢化”，不能把全部差距归因于 ICCG 或审计。受控拆分由 R1 完成。

另一组冻结的两步 core-hotpath 证据为 HUNDUN 39.096166 s、严格 COAST
39.358670 s，差异约 0.67%，支持“严格 COAST 与慢 HUNDUN 可出现相近耗时”。但
full20 的 COAST P90 为 15.419435 s/step、当时 HUNDUN P90 为 4.418378 s/step，说明
这种相近只成立于该两步诊断，不能外推为长期等速。

证据：

- `/home/wyf/code_dev/.benchmarks/c3900s/cause.md`
- `/home/wyf/code_dev/.benchmarks/c3900s/r10/re3900_evidence.jsonl`
- `/home/wyf/code_dev/.benchmarks/coast-hundun-compressible-modules-20260902/report.md`
- `/home/wyf/code_dev/.benchmarks/coast-hundun-compressible-modules-20260902/summary.json`
- `/home/wyf/code_dev/.benchmarks/pa/final.md`
- `/home/wyf/code_dev/.benchmarks/pa/report.md`

### 20D x 10D x 3D 生产算例

- `456 x 256 x 104 / 128 ranks`，每核 94,848 cells；
- 49 个 BDF2 步中位 34.809398 s；
- Stage 50 中位 28.366345 s；
- pressure-energy refinement 最小/中位/最大为 3/6/7；
- 每增加一次 refinement 约增加 3.94 s，相关系数 0.9951；
- 代表步有 6,046 次阻塞集合通信，记录的 reduction 时间为 6.781 s。

证据：`/home/wyf/r39m/c12/x50/evidence.jsonl`。

主计划收口时，同一 step-50 restart 的 H0/H5 单步中位为
`41.473374 / 22.391036 s`，累计降低 `46.01%`；pressure-energy refinement 从 6 次
降为 3 次，blocking collectives 从 6,537 次降为 3,201 次。该对比用于表示本计划各节点
的累计工程收益；H5 自身的归因仍以 H3/H5 干净构建对照为准。

## 已验证的 COAST 差异及其含义

| 维度 | 普通可压缩 COAST | 严格 Re3900 COAST | 对 HUNDUN 的约束 |
|---|---|---|---|
| 压力后端 | `auto` 需 20 samples 才锁定；当前 6-solve 短跑仍混合探测，Hypre 不可用时回退 GMG | 强制 `legacy_iccg`，禁止 fallback | 保留 HUNDUN FGMRES/MG；只比较预条件器总工作量 |
| 压力迭代 | 每步末次 GMG 记录为 0--7 次，不能代表全部 auto 探测 solve | c3900s 的 20 个压力 solve 全部 256 次；full20 少量达到 384 次 | 优先消除过求解，而非更换为 ICCG |
| 线性控制 | U/V/W/H/P 多为 `1e-4` | 五类方程均强制 cap=500、内部 `1e-6` | 内层容差与终端产品门必须分离 |
| 数据精度 | 主场、系数和工作数组为默认 `REAL`/FP32 | 主存储仍为 FP32 | 严格模式并未全局升级为 FP64 |
| 真残差 | 求解器原生递归/缩放残差 | FP64 累加，解前/解后及周期中间审计 | HUNDUN 保留 FP64 外层真残差，可测试低精度预条件器 |
| 审计触发 | 无 Re3900 强制审计 | 压力每 128 次；BiCGStab 每 8 次 | 审计必须在“预测已收敛”时立即触发，不能只等固定间隔 |
| 产品验收 | 普通 legacy step 路径 | EOS/continuity/gauge/closed-mass 终端接受 | HUNDUN 终端门保留，审计 pass 可融合但不能删除 |
| 数值路径 | 普通可压缩 SIMPLE | BDF/密度缺陷/压力存储、补偿及边界等价钩子 | 2.63 倍差距不是纯压力后端 A/B |

严格 COAST 的压力循环只有在固定 128 次边界才执行真残差检查；现有证据中第一次审计
未通过后，只能继续到 256 次再退出。HUNDUN FGMRES 已在递归残差达到容差、固定间隔或
restart 末端三种条件下执行真残差，因此不存在完全相同的 128/256 量化停机问题；但当前
`true_residual_interval=4` 可能产生过密的完整 Schur operator/halo/归约审计。

严格 COAST 的 attribution 将约 77% 总时间归到 ICCG inclusive；终端审计两步均值约
0.1389 s，只占约 0.34%。因此保留终端精度门并优化求解工作量，比删除终端审计更合理。

## 推进节点

### R0：源码与历史证据研究

状态：`COMPLETE`

- [x] 固定严格 COAST 身份、输入和历史性能证据。
- [x] 固定普通可压缩 COAST 身份、输入和模块计时证据。
- [x] 核查 backend、fallback、线性容差、存储精度、真残差和终端审计。
- [x] 将 verified、inference、unknown 分开记录到配套研究文档。

### R1：同二进制 COAST 受控差异矩阵

状态：`COMPLETE`

在隔离工作树中，用同一 profiling binary、同一 320 x 320 x 33 case、相同输出设置和
相同 rank mapping 至少运行以下诊断组合：

| 组合 | Re3900 mode | 压力后端 | 线性容差 | 用途 |
|---|---|---|---|---|
| C0 | off | `auto` + GMG fallback | `1e-4` | 重现普通基线 |
| C1 | off | 强制 `coast_gmg_pcg`，无 fallback | `1e-4` | 去除 auto screening 噪声 |
| C2 | off | 强制 `legacy_iccg`，无 fallback | `1e-4` | 隔离后端差异 |
| C3 | off | 强制 `legacy_iccg`，无 fallback | `1e-6` | 隔离线性容差工作量 |
| C4 | on | 严格 `legacy_iccg`，无 fallback | `1e-6` | 测量严格物理/审计合并增量 |

每组至少记录 1 个 BE + 5 个 BDF2 步，并交替运行三轮。必须输出：

- max-rank step、pressure、momentum、terminal-audit 时间；
- 每类线性方程的迭代数、operator/preconditioner applies；
- FP64 真残差审计次数、全场遍历次数；
- `Allreduce/Allgatherv` 和 halo 次数；
- EOS/continuity/gauge/closed-mass 和压力真残差。

如果 C3 -> C4 仍有显著差距，再做严格 COAST 的诊断补丁：当递归残差首次达到门槛时
立即执行 FP64 真残差，固定间隔只承担 residual replacement。补丁只进入临时性能工作树，
不改变普通 COAST 主线。

**进入 H0 的门：**同二进制矩阵完整，能够把差距至少分成 backend、tolerance 和
strict-physics/audit 三类；否则不修改 HUNDUN 求解器。

完成证据（2026-09-02）：

- 同一 profiling binary、同一 320 x 320 x 33 case、128 ranks、固定绑核完成 C0--C4
  各三轮；每轮 1 BE + 5 BDF2，15/15 运行成功且 critical-rank self partition 闭合。
- BDF2 中位耗时 C0/C1/C2/C3/C4 分别为
  `0.398283 / 0.426525 / 0.415725 / 0.679622 / 1.203166 s`。
- C0 -> C1 强制 GMG 为 `+7.09%`，C1 -> C2 legacy ICCG 为 `-2.53%`；两者处于短跑
  auto 探测和后端噪声量级，不能解释主要慢化。
- C2 -> C3 将线性控制从 `1e-4` 收紧到 `1e-6` 后增加 `0.263898 s`（`+63.48%`）；
  C3 -> C4 严格物理/可靠残差增加 `0.523543 s`（`+77.03%`）。
- C4 terminal audit 中位仅 `0.012083 s`，不是 C3 -> C4 增量的主因；C4 压力
  operator/preconditioner 中位为 `520/512`，六步每 rank `MPI_Allreduce=17262`，而
  C3 为 `228/228` 和 `9381`。五组 `MPI_Allgatherv` 均为 0。
- 严格模式 18/18 步终端接受，180 次线性求解共 648 次 FP64 真残差审计，无 cap hit；
  EOS/continuity/gauge 最大残差为
  `1.102e-7 / 1.056e-7 / 2.139e-8`。

证据目录：`/home/wyf/code_dev/.benchmarks/coast-r1-controlled-matrix-20260902`；
`summary.json` SHA-256
`ceb8681f28265a88edcb1ad7c7c499918453c976e7bf53703bce06924c034902`，
`report.md` SHA-256
`e51c73efbdc5d188c7d4ab7931f56ead09f416212343aae27a9fa0e4236362d7`。

结论：R1 门已满足。backend 本身不是首要矛盾；线性过求解及严格路径增加的
operator/collective 工作量是可迁移到 HUNDUN 的主要优化信号。

### H0：HUNDUN 低风险线性控制

状态：`COMPLETE`

1. 以当前 `true_residual_interval=4` 为基线，A/B 测试 8 和 16。
   递归残差达到门槛时仍立即执行 FP64 真残差，restart 末端仍强制审计。
2. 压力内层相对容差以已实测的 `1e-4` 为首选，绝对容差及终端五门不变。
3. 增加按非线性残差收紧的 inexact forcing：早期 refinement 使用 `1e-4`，只有接近
   终端门或停滞时才收紧。

**验收：**短算例相对当前 `1e-6` 基线至少提速 20%；最终 FP64 真残差、五门、两次
PISO、迭代终止原因和 rollback 全部通过。任何 `interval` 造成 audit reject、restart 增加或
非线性 refinement 增加都不接受。

完成证据（2026-09-03）：

- 接受提交 `60cf89ec74ef2933d27067711587a9d315645bb3`，tree
  `2888797f341a7fd37459354588f00058b70c52cb`。保留 case `rtol=1e-6` 和全部
  `1e-6` 终端门，只在两次强制 PISO 方向及远离终端的 refinement 内使用受保护的
  `1e-4` inexact forcing；停滞、无效残差、接近终端或 case 更严格时恢复 case 控制。
- `true_residual_interval=8/16` 单变量收益不足；接受配置使用 interval 32，但递归残差
  达标和 restart 末端仍立即做 FP64 真残差，旧 checkpoint 的计划指纹不匹配时明确拒绝。
- 320 x 320 x 33、128 ranks 的 15 个 BDF2 步中位/P90 从
  `2.592414 / 2.807993 s` 降至 `2.033780 / 2.154805 s`，分别降低
  `21.55% / 23.26%`；线性迭代 `50 -> 36`、operator `65 -> 40`、blocking
  collectives `1372 -> 1196`，refinement 保持 0。
- 三轮 runtime validator 及 16 项原聚焦检查通过。产品时序收敛测试的 stage-44
  `2/5790` 已在精确的 H0 前提交 `faa8ad2` 上复现，登记为既有失败而非 H0 回归。

证据目录：`/home/wyf/code_dev/.benchmarks/hundun-h0-inexact-forcing-20260903`。

### H1：减少 pressure-energy refinement 数量

状态：`COMPLETE`

1. 用已保存的成功轨迹补丁确认当前 alpha、merit 和各终端残差收缩率。
2. 若仍为 `alpha=1` 且线性收敛，先测试 depth=2--3 的 safeguarded Aitken/Anderson。
3. 新方向必须走现有 exact candidate evaluator；不合格立即回退当前方向。

**目标：**20D x 10D x 3D 的 refinement 中位数由 6 降到不高于 3。按现有拟合，
降到 2 次的理论节省约 15.8 s/步；该数值是预测，不是验收承诺。

完成证据（2026-09-03）：

- 接受提交 `13ac181da171ed20ed544840ee423050bbbbdd3c`，tree
  `353647897ba51b783881d943cfad5bdf3b47c329`。以终端实际使用的
  `max(continuity, energy)` 为标量，使用去除上一轮 relaxation 影响的 safeguarded
  Aitken 系数；系数上限为 2，且 `alpha>1` 只生成候选，不获得发布权限。
- 所有外推态仍经过原 exact candidate evaluator、正性/有限性、严格 merit decrease、
  Armijo、谱系绑定和 commit/rollback；外推拒绝或数值失败时执行未改动的
  `alpha=1,1/2,...` ladder。
- 同一 step-50 restart 的诊断轨迹中，接受系数为 `1.13691 / 1.77344 / 2.0`；第 3 次
  refinement 后 continuity/energy 为 `9.940878e-7 / 7.947763e-7`，三次重复逐值一致。
- 干净构建的 20D x 10D x 3D 正式三轮中位/P90为 `24.262052 / 24.524881 s`；相对
  H0 同重启 `41.473374 s`，压力路径迭代 `240 -> 156`、operator `314 -> 205`、
  blocking collectives `6537 -> 4207`，refinement `6 -> 3`。Stage 40+50 中位为
  `22.610701 s`，reduction 中位为 `5.977023 s`。
- 干净构建的 globalization、PISO authority、1/2/4-rank 核心产品、故障注入回滚和
  pressure-energy retry 聚焦测试 18/18 通过。测试诊断另外修正了已存在的 controller
  ticket/已提交流场混淆、候选 alpha 元数据丢失，以及相消残差使用错误尺度的问题。

正式证据目录：
`/home/wyf/code_dev/.benchmarks/hundun-h1-pressure-energy-20260903/formal-runs`；三轮
`evidence.jsonl` SHA-256 分别为
`c8fbe69e3be7cdba620a54e82c16b39b8c5232e79df13474c6e601d1330abe47`、
`7d172ce606ec952852b870260bd44fb52339ff171165f0b4b0723ba56c4443cb`、
`439dd00a49e223728c81d1c395939ea864f67ba11ee1d0d458a5d183a9354f4d`。

### H2：减少压力路径集合通信

状态：`COMPLETE`

- 保留 MG/Schur 已有 prepared batch；
- 将 prepared epoch 扩展到候选态、AFC 和可共享失败 authority 的相邻阶段；
- 将状态发布搭载到不可避免的残差归约；
- 统计实际 PMPI `Allreduce/Allgatherv`，不只看 evidence 内 solver reduction counters。

**验收：**相同科学工作量下集合通信次数和 reduction 时间下降；failure injection、最低
失败 rank、persistent request 恢复以及 ghost revision 证书测试全部通过。

完成证据（2026-09-03）：

- 接受提交 `7bc64a60273a26d0cf883b7d24e605c3db3a6be1`，tree
  `6019b00a3defe40bb8df78812e04aa31d5dc2196`。PMPI 调用点审计确认热点不是 Schur
  apply：205 次 Schur apply 没有 `Allreduce`；主要冗余来自每次系数刷新时粗网格张量
  聚合的成功状态共识，以及 156 次 prepared-MG apply 的末端发布共识。
- 张量聚合保留每次 `MPI_Sendrecv` 前的集体前检，把其完成状态和本地映射状态延迟到
  下一次必需前检（最后一个轴延迟到函数末端前检），不跨越下一次 MPI；prepared-MG
  则把 Halo/workspace/counter 状态搭载到已有 FP64 投影归约，归约后的判断只依赖各 rank
  相同的全局投影。
- 同一诊断二进制和 step-50 restart 下，压力/候选作用域 PMPI `Allreduce` 从
  `3020 -> 2184`（`-27.68%`），`Allgatherv` 保持 `1096`；完整 evidence 的 blocking
  collectives 从 `4207 -> 3201`，恰好减少 `850` 次张量聚合共识和 `156` 次 prepared-MG
  发布共识。
- 干净 Release 构建的 20D x 10D x 3D 三轮中位/P90为
  `23.099072 / 23.168112 s`；相对 H1 中位 `24.262052 s` 再降低 `4.79%`。reduction
  中位为 `5.669419 s`，Stage 40+50 中位为 `21.455973 s`，分别较 H1 降低
  `5.15% / 5.11%`。
- 三轮均保持 162 次总线性迭代、3 次 refinement，以及完全相同的压力/精化
  operator/preconditioner/reduction 计数；terminal continuity/energy 仍逐值为
  `9.940878e-7 / 7.947763e-7`。干净构建的 MG 数值/复用、1/2/4-rank MG、奇数分区、
  更新契约、Krylov、Halo 故障与恢复、核心产品和 pressure-energy retry 共 27 项检查
  全部通过，并新增同一 persistent requests 在 deferred failure 后成功重试的断言。

正式证据目录：
`/home/wyf/code_dev/.benchmarks/hundun-h2-communication-20260903/formal-runs`；三轮
`evidence.jsonl` SHA-256 分别为
`471e718880441479c674dda3bb36605cdd3d7a513bfb8ed12a571f176955d8b1`、
`79090e7e93588bf3153c5ff307f93054e690edd72fbf174cc27c3de9f50a3edd`、
`a9e392a94d3b16441488290e1b07c2fda400d465ae5e951b00712ca4ede0ef61`；正式 runner
SHA-256 为 `43bad1cf70a3c98f15232d06bdd2a430a24078c21bce7252f0e4df1c48d6c848`。

### H3：减少全场数组 pass

状态：`COMPLETE`

- 在 EOS/transport/residual 生成 pass 中同步产生审计与 provenance；
- 合并相邻 zero、mask、EOS、transport、residual、hash pass；
- production 使用 revision/operation lineage，完整数值 hash 保留在诊断路径；
- dot product 使用多路 FP64 局部累加和确定性补偿合并；
- 内核缓存行首和 stride，不改变场布局。

**验收：**全场 pass/cell visits 减少，最终状态和终端审计在既定误差门内一致；不得只删除
certificate 或 hash 语义来获得提速。

完成证据（2026-09-03）：

- 接受提交 `6206bdc28c46b822ef6bd1d328ce87dc421b61cc`，tree
  `641fa4951d84a7a0cd366c089320e77add8e38ef`。冻结候选在生产核心中以完整
  pointer/layout/storage/revision-domain/revision lineage 绑定；测试核心仍逐值计算原 FP64
  位级 hash，因此无 revision 静默篡改的诊断覆盖没有删除。
- H2/H3 优化前后 gprof 均观察到每步每 rank 约 55 次冻结候选指纹调用；H2 每次扫描
  density、`HbyA`、`rAU`、压力/梯度（C1）及两组面通量，H3 生产路径的这些调用为
  O(1)，数组值访问从 55 次组合审计降为 0。优化后 profile 中该函数自耗时降为 0；
  完整数值 provenance、真实残差和终端审计仍保留在原路径。
- 干净 Release 构建的 pressure-energy Schur/globalization、PISO authority/mutation、
  1/2/4-rank PISO、核心产品及 pressure-energy retry 共 17 项检查全部通过。正式二进制
  SHA-256 为 `8b86c2359c51c65206436dba7fad93e5db25bae03babab5a9ff3d32e5f020286`，
  build manifest SHA-256 为
  `25e081b2155974cdd617bee8577fcef3532c77c77ae08231cad4b37e4d658949`。
- 清除一次遗留递归搜索造成的 I/O/load 污染后，采用 H2/H3、H3/H2、H2/H3 反序配对
  各三轮。H3/H2 单步中位为 `22.481344 / 23.077041 s`（`-2.58%`），三组配对
  相对变化中位为 `-4.46%`；Stage 40+50 中位为
  `20.843059 / 21.394500 s`（`-2.58%`），reduction 中位为
  `5.568769 / 5.779437 s`（`-3.65%`）。
- 六轮均为 162 次总线性迭代、205 次压力 operator、156 次 preconditioner、420 次
  solver reduction、3 次 refinement 和 3,201 次 blocking collective；H3 三轮终端
  continuity/energy 均逐值为 `9.940878e-7 / 7.947763e-7`。被遗留搜索污染的绝对
  时间三轮保留在证据目录，但明确排除在性能统计之外。

正式证据目录：
`/home/wyf/code_dev/.benchmarks/hundun-h3-array-passes-20260903/controlled-runs`；H3 三轮
`evidence.jsonl` SHA-256 分别为
`2121bc78c0af7e4fef877603781311e9956c278035e52b62d7c9062ac0aa9e47`、
`b77c61e7ca1f4c0bc2622005c72d1ab5d17e557fd21b0f3f137b898d19447f52`、
`41e4b52708055c30ef3f519601c01a5aa059d6ac7776108f0f7b56e401b986eb`。优化后 gprof
报告 SHA-256 为
`08b44a538260c758b954feca2400ed294d17d0de4b9c5d07684dde8e6d9a923b`。

### H4：预条件器工作量与分层精度

状态：`COMPLETE / 保留 FGMRES + F-cycle`

- 保留 FGMRES + F-cycle 作为可靠基线；历史证据中 F-cycle 快于 guarded V-cycle；
- A/B 测试 communication-light block-local/Schwarz 预条件器，但以总时间和总工作量判断；
- 仅在 MG 预条件器的工作数组、粗网格和 halo 中测试 FP32；
- 外层 operator、FGMRES 解、状态、FP64 真残差及终端审计保持 FP64；
- 不引入运行期 auto profiler，先离线选定明确后端并保留可靠 fallback。

完成证据（2026-09-03）：

- 以 H3 接受代码和同一 step-50 restart 做单变量诊断。最近的受控 F-cycle 样本为
  `22.048897 s`、162 次线性迭代、205/156 次 operator/preconditioner 和 3,201 次
  blocking collectives。
- 固定 V-cycle 为 `25.818467 s`（相对 F-cycle `+17.10%`），线性迭代升至 200，
  operator/preconditioner 为 `254/194`；首次系数刷新用 F、随后用 V 的混合方案为
  `26.566787 s`（`+20.49%`），迭代仍为 198。两者虽减少消息数，但增加的压力工作量
  抵消并超过通信收益，均拒绝。
- 四次红黑局部 sweep 的零重叠 rank-block Schwarz 将 structured messages 从
  `13,054,464` 降至 `2,170,368`，但线性迭代增至 1,046、operator/preconditioner
  增至 `1,323/1,040`、通信字节增至 `171,532,513,280`，单步为 `98.235333 s`
  （`+345.53%`），拒绝。
- 三个候选的终端 continuity/energy 仍在 `1e-6` 门内，因此拒绝依据是总工作量和总时间，
  不是精度失败。对应 evidence SHA-256 分别为
  `376015024ca3d9e0be2e2f4e5aa6746f8e4532666f5b8b8bf26eaab5409cecdb`、
  `71eece3bed99ae47fa7166b25df6e36f2bdd307264ed320d38579d697698a73c`、
  `ab3d093453813b80dd67d05fc13189ff069aa6f116926434f475b9e7fd96aec5`。
- 分层 FP32 审计确认现有 MG 数据契约从 `BasicFaceFieldView`、`MgWorkspaceRequirements`、
  `MgWorkspace`、层级/复制缓冲到 `HaloEngine` 均固定为 `double`/`MPI_DOUBLE`。把数值量化
  后仍存入 double 不能测试存储带宽或低精度 halo；真正的 FP32 候选需要跨层 typed
  workspace/halo 重构，超出本节点的低风险单变量改动范围，故不实施。外层状态、FGMRES
  解、Schur operator、FP64 真残差及终端审计维持原样。

诊断证据目录：
`/home/wyf/code_dev/.benchmarks/hundun-h4-preconditioner-20260903`。

### H5：动量 AFC 与硬件收尾

状态：`COMPLETE`

- 复用冻结的高阶面修正，减少 AFC 二次计算；
- 批量 limiter halo，保留唯一 common-face authority 和等量反号守恒；
- 最后测试 `-march/-mtune`、NUMA first-touch/绑核及节点内 halo。

Stage 30 明显小于 Stage 40+50；H0--H4 未完成前，不以 AFC 或编译器调优替代压力根因修复。

完成证据（2026-09-03）：

- 接受提交 `4476ac5a49da6037eb2ffab7285efbe12212091a`，tree
  `2e6db419a4ab4a51f09d77d6591ec404292a7394`。每个速度分量的高阶面重构由两次降为一次，
  三个分量分别保存在互不别名的 face workspace；第二次 AFC correction/metrics pass
  直接复用冻结结果。common-face alpha、owner authority、等量反号守恒、failure-before-
  mutation 和原始分量顺序不变。
- `phi_workspace` 从 4 个扩为 6 个 replica：0/1 仍服务时间/候选通量，2/3/4 保存三个
  分量的 AFC 面场，5 保存 common alpha。checkpoint 指纹仍只编码原有四个持久语义角色，
  两个纯临时缓存不改变 restart 兼容性；旧 step-50 restart 已在正式运行中通过。
- 只保留最后一个分量的首个候选在三组配对中 Stage 30 收益不稳定，已拒绝。完整三分量
  缓存的提交前诊断三组均使 Stage 30 降低 `4.18%--6.29%`，因此进入干净构建验收。
- 干净 Release 构建的关键回归共 21/21 通过，覆盖 pressure-energy Schur/globalization、
  PISO authority/mutation、solver equations、AFC 产品冻结、pressure-energy retry、应用驱动
  以及 MPI 1/2/4 ranks；`focused-21.log` SHA-256 为
  `b812daf1c4de9c1ae1cfe82a742c90fc07e97ef92a48a56c48bf2834f4e4418c`。
- H3/H5 三次正式单步中位如下；Stage 40+50 和 reduction 的小幅变化不归因于 AFC：

| 指标 | H3 | H5 | 变化 |
|---|---:|---:|---:|
| max-rank step | 22.481344 s | 22.391036 s | -0.40% |
| Stage 30 | 1.348177 s | 1.276689 s | -5.30% |
| Stage 40+50 | 20.843059 s | 20.809205 s | -0.16% |
| reduction | 5.568769 s | 5.489720 s | -1.42% |
| max-rank RSS | 409,706,496 B | 414,953,472 B | +5,246,976 B / +1.28% |

- 三轮 H5 均保持 162 次总线性迭代、205/156 次压力 operator/preconditioner、3 次
  refinement、3,201 次 blocking collectives，以及完全相同的 structured message/byte
  计数。AFC retained ratio 均为 `0.9976554563922444`，terminal continuity/energy 逐值为
  `9.940877864179103e-7 / 7.947762810986709e-7`。
- 另一次连续 step 51--53、128-rank 运行完整通过；三步的 refinement 均为 3，EOS、
  closed-mass、gauge 均为 0，continuity/energy 均低于 `1e-6`，CFL 约 0.339 且低于 0.8。
  该运行验证缓存 revision 随连续时间步推进，不纳入 H3/H5 性能中位数。
- 正式 runner SHA-256 为
  `fb9109bc0fd59f9cc2a4f833cddcbabe1cbc5455cb30d7cbdb8beefba0fe8d8f`，build manifest
  SHA-256 为
  `8bce3d923b1b0bbc6ac1af99bb02d5b64bb55dc3bf3286594ec9de2647899c76`。三轮正式
  `evidence.jsonl` SHA-256 分别为
  `fd1d17a21a4a9ea88c6fa0bb177c754ace4e5aace72a71d5e76a6f1cab149e73`、
  `291a601167697cec5d2e89035488aef44e9d3387dcc03c03873477e542296d9e`、
  `b86ee389aced682bfdc9bde8655a59ab844fdf3873344df2e27ac93f55f581cd`；连续三步证据
  SHA-256 为 `c5686b8121d4d279fd2ad586d2b5f3fee0550d3335968c857928a58fd75b4934`。
- `-march=znver3` 与额外 `-mtune=znver3` 的实际 runner `.text` 均为 4,593,509 bytes，
  SHA-256 同为 `df96cdbdee84744ded2f140e5d5e19b655f7721f215452e47950028596d4d777`，
  因此不增加无效编译开关。128-rank 命令已按两 socket 各 64 rank、每 rank 一个物理核
  绑核；分配发生在绑定进程内，Open MPI 提供 `vader/sm/self/tcp` 节点内传输。
- `codex/v04-intranode-halo-p0` 仍是 `P0_PLANNED_DEFERRED_NOT_IN_V1` 文档封存分支，
  没有可验收的节点内 halo 实现或测试，故不合并、不宣称能力，继续作为独立保留分支。

正式证据目录：`/home/wyf/code_dev/.benchmarks/hundun-h5-afc-20260903/formal-runs`。

## 每个节点的统一验收

1. 一次只提交一个可归因变化；保存 exact HEAD/tree、构建清单、binary SHA-256、case SHA-256。
2. 先跑 focused numerical/MPI tests，再跑 1 BE + 3--5 BDF2 短测；失败不进入长算例。
3. 对比 max-rank 中位时间、P90、迭代、refinement、operator/preconditioner applies、halo、
   collectives 和全场 pass；不能只报告单一 wall time。
4. 检查 EOS、continuity、energy、mass、gauge、CFL、AFC、正性、force/flux 及 rollback。
5. 性能变化小于噪声区间或依赖更弱终端门时，结论为 `REJECT`。
6. 接受后在本文件更新节点状态、commit/tree 和证据路径，再进入下一节点。

## 收口状态与后续

本计划已经完成；工作主线为 `codex/re3900-10d-smoke`，可靠后端保持
FGMRES + F-cycle，H5 三分量 AFC 缓存已进入主线。后续如需生产放大，应另立验证任务，
先运行更长的连续 BDF2 窗口并检查统计量、力/通量和内存高水位，不在本计划内继续叠加
求解器或精度变化。节点内 halo 只有在明确解除 P0 封存并重新完成 MPI/NUMA intake 后再评估。

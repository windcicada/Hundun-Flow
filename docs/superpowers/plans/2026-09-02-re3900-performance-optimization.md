# HUNDUN-FLOW Re3900 性能优化计划

日期：2026-09-02

状态：`ACTIVE / 当前节点 R1_COAST_CONTROLLED_DELTA`

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

状态：`PENDING / NEXT`

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

### H0：HUNDUN 低风险线性控制

状态：`PENDING`

1. 以当前 `true_residual_interval=4` 为基线，A/B 测试 8 和 16。
   递归残差达到门槛时仍立即执行 FP64 真残差，restart 末端仍强制审计。
2. 压力内层相对容差以已实测的 `1e-4` 为首选，绝对容差及终端五门不变。
3. 增加按非线性残差收紧的 inexact forcing：早期 refinement 使用 `1e-4`，只有接近
   终端门或停滞时才收紧。

**验收：**短算例相对当前 `1e-6` 基线至少提速 20%；最终 FP64 真残差、五门、两次
PISO、迭代终止原因和 rollback 全部通过。任何 `interval` 造成 audit reject、restart 增加或
非线性 refinement 增加都不接受。

### H1：减少 pressure-energy refinement 数量

状态：`PENDING`

1. 用已保存的成功轨迹补丁确认当前 alpha、merit 和各终端残差收缩率。
2. 若仍为 `alpha=1` 且线性收敛，先测试 depth=2--3 的 safeguarded Aitken/Anderson。
3. 新方向必须走现有 exact candidate evaluator；不合格立即回退当前方向。

**目标：**20D x 10D x 3D 的 refinement 中位数由 6 降到不高于 3。按现有拟合，
降到 2 次的理论节省约 15.8 s/步；该数值是预测，不是验收承诺。

### H2：减少压力路径集合通信

状态：`PENDING`

- 保留 MG/Schur 已有 prepared batch；
- 将 prepared epoch 扩展到候选态、AFC 和可共享失败 authority 的相邻阶段；
- 将状态发布搭载到不可避免的残差归约；
- 统计实际 PMPI `Allreduce/Allgatherv`，不只看 evidence 内 solver reduction counters。

**验收：**相同科学工作量下集合通信次数和 reduction 时间下降；failure injection、最低
失败 rank、persistent request 恢复以及 ghost revision 证书测试全部通过。

### H3：减少全场数组 pass

状态：`PENDING`

- 在 EOS/transport/residual 生成 pass 中同步产生审计与 provenance；
- 合并相邻 zero、mask、EOS、transport、residual、hash pass；
- production 使用 revision/operation lineage，完整数值 hash 保留在诊断路径；
- dot product 使用多路 FP64 局部累加和确定性补偿合并；
- 内核缓存行首和 stride，不改变场布局。

**验收：**全场 pass/cell visits 减少，最终状态和终端审计在既定误差门内一致；不得只删除
certificate 或 hash 语义来获得提速。

### H4：预条件器工作量与分层精度

状态：`PENDING`

- 保留 FGMRES + F-cycle 作为可靠基线；历史证据中 F-cycle 快于 guarded V-cycle；
- A/B 测试 communication-light block-local/Schwarz 预条件器，但以总时间和总工作量判断；
- 仅在 MG 预条件器的工作数组、粗网格和 halo 中测试 FP32；
- 外层 operator、FGMRES 解、状态、FP64 真残差及终端审计保持 FP64；
- 不引入运行期 auto profiler，先离线选定明确后端并保留可靠 fallback。

### H5：动量 AFC 与硬件收尾

状态：`PENDING`

- 复用冻结的高阶面修正，减少 AFC 二次计算；
- 批量 limiter halo，保留唯一 common-face authority 和等量反号守恒；
- 最后测试 `-march/-mtune`、NUMA first-touch/绑核及节点内 halo。

Stage 30 明显小于 Stage 40+50；H0--H4 未完成前，不以 AFC 或编译器调优替代压力根因修复。

## 每个节点的统一验收

1. 一次只提交一个可归因变化；保存 exact HEAD/tree、构建清单、binary SHA-256、case SHA-256。
2. 先跑 focused numerical/MPI tests，再跑 1 BE + 3--5 BDF2 短测；失败不进入长算例。
3. 对比 max-rank 中位时间、P90、迭代、refinement、operator/preconditioner applies、halo、
   collectives 和全场 pass；不能只报告单一 wall time。
4. 检查 EOS、continuity、energy、mass、gauge、CFL、AFC、正性、force/flux 及 rollback。
5. 性能变化小于噪声区间或依赖更弱终端门时，结论为 `REJECT`。
6. 接受后在本文件更新节点状态、commit/tree 和证据路径，再进入下一节点。

## 当前下一步

只执行 R1 的同二进制 COAST 受控差异矩阵。R1 完成并复核前，不开始 H0 代码修改。

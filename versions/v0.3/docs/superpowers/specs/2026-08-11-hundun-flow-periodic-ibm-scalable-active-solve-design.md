# HUNDUN-FLOW 周期 IBM 与可扩展活跃域求解设计

## 1. 状态与目标

本设计在已接受的 Stage 4 治理 seal
`033a685c90c1f9c674e93a4b82db10db4c381abe` 上继续工作。受测 Stage 4
产品代码仍为 `6407cd7c591ce088db7f1dd7e296d77acd18da1c`。目标是在不改变
Stage 4 反应流合同的前提下，融合 Re=3900 圆柱基准分支的周期 IBM 语义，并让
冻结的圆柱算例能够正常完成低成本验收。

正式物理配置保持 `Re=3900`、`D=U=rho=1`、`mu=1/3900`、
`dt=0.006`、两次 PISO corrector、WALE 和当前力权威。不得通过放宽残差阈值、
增加 corrector、滤波、阻尼或逐算例调参获得通过。96³永久不运行；正式统计长跑
只能在功能门通过后启动。

## 2. 已有证据与根因边界

原 Stage 3 产品树上的 64-rank 预检已经依次证明：

1. 封闭母 STL 可以在周期方向同时跨越一对周期面；分类使用完整母表面，coverage、
   wall quadrature、surface measure 和 force 只使用半开基本域内的表面片；
2. canonical cell identity、stencil-anchor-nearest image 和
   owned-box-nearest field image 必须分离；
3. LFP 周期邻点必须保持 canonical field identity，同时使用相邻周期镜像的几何；
4. 活跃压力和动量算子保留每 rank 全局 ID/hash 表会导致全尺寸构造 OOM；
5. 释放构造期索引后，48³首步能够越过周期几何、重构、collective 和内存阶段，
   但最终返回 `pressure_linear_solve`。

现有错误把 exact-predictor 内部动量 CG 与外层压力 BiCGStab 的失败折叠为同一原因。
失败耗时更接近首次 exact response 的三个内部动量求解，但在记录完整
`SolveReport` 之前不能把该推断当成结论。

## 3. 方案比较

### 3.1 采用：项目内 compact sparse Halo + 证据驱动求解修复

扩展现有 `GhostedVectorHalo`，使其除完整拓扑布局外还能接受任意“owned IDs 在前、
ghost IDs 在后”的 `VectorLayout`。构建期允许交换稀疏请求 ID；运行期只与实际 peer
交换边界值。活跃压力、动量和三分量动量对角线都复用该通道，不再为每次 matvec
执行全局值 `Allgatherv`。

先记录 exact predictor 内层和外层的独立报告，再只修被证明的求解层。优点是复用
Stage 2 已验收的 compact Halo 失败协议、chunking、生命周期和性能计数，不增加运行
时依赖，也不改变 exact Schur 数学映射。

### 3.2 不采用：引入 PETSc、HYPRE 或 Trilinos

成熟外部包能提供 AMG/ASM，但会立即改变 Stage 4 的无 Python/离线打包、ABI、RPATH、
许可证和验收范围。本基准修复不应顺带建立新的生产 solver backend。

### 3.3 不采用：绕过 exact predictor，直接求 compact Poisson

该方案开发成本最低，但会恢复 Task 11 已拒绝的近似压力权威，并可能破坏压力、最终
通量、operator force 和 surface consistency 的科学闭环。它不能作为 fallback。

## 4. 融合架构

新工作线从 Stage 4 seal 建立，保留原产品基准工作树原样。融合分为两层：

- 产品补丁：周期 surface window、periodic donor/LFP、活跃算子内存生命周期和失败
  诊断；
- 治理与基准资产：HUNDUN 圆柱 STL/case、C++ RED 和本设计记录。Python
  生成/runner/分析脚本及其测试只保存在仓库外 maintainer 根
  `/home/wyf/code_dev/.benchmarks/cylinder-re3900-stage4/maintainer-python`，不进入
  HUNDUN 正常 configure/build/install/test/runtime 的可见源码树。

旧 COAST case、patch、binary/source manifest 和双程序内部计划只保存在仓库外
`/home/wyf/code_dev/.benchmarks/cylinder-re3900-stage4/private-comparison-snapshot`。
它们不属于 HUNDUN 产品或治理历史，也不作为本修复的实现来源。

融合前后以 Stage 4 45-test 低成本矩阵和圆柱 focused RED 检查。Stage 4 的
Cantera、schema v4、Checkpoint v4 和 reacting coupling 文件如果没有调用方影响，
直接复用已接受证据，不运行化学长算例。

## 5. Sparse active Halo 合同

新增 `GhostedVectorHalo::create(..., VectorLayout)` 重载，合同如下：

- `VectorLayout` 的前 `owned_count` 个 ID 必须由本 rank 唯一拥有，其余 ID 必须是
  当前拓扑可见的 remote ghost；
- 每个 ghost 的 owner rank 由结构分解和 canonical cell coordinate 唯一确定；
- 构建期通过 request-ID 交换冻结 `peer -> send_indices/receive_indices`；运行期 payload
  与全局 active cell 数无关，仅与跨 rank 活跃邻接有关；
- 同一 global ID 不得重复、不得请求固体/缺失 ID、不得由错误 rank 发布；所有错误
  通过现有 collective failure 协议一致报告；
- full-topology 旧重载的 API、数值顺序、failure injection 和性能计数保持不变。

`ActivePressureOperator` 和 `ActiveMomentumOperator` 的 `VectorLayout` 已经包含 owned
与 ghost active IDs，因此它们各自持有一个 `GhostedVector` workspace 和 sparse Halo。
apply 顺序固定为：复制 owned input、`begin`、计算纯本地连接、`wait`、计算 remote
连接。压力 mobility 所需三分量动量对角线也通过同一 active Halo 顺序交换。

## 6. 求解失败可观测性与修复决策

每个 exact response 记录三分量内部动量 `SolveReport`；压力 correction 保留外层
BiCGStab 报告和独立重算残差。失败消息至少包含：

- `phase=affine_inner_momentum|homogeneous_inner_momentum|outer_pressure|independent_residual`；
- component；termination reason；iterations；initial/recursive/final residual；
- matvec、preconditioner 和 reduction 次数；最低失败 rank。

决策顺序固定：

1. 内层动量失败：先验证算子 SPD、对角正性和分解一致性。若仅 Jacobi 条件数不足，
   增加一个项目内、固定算法参数且 mutation-sensitive 的局部块预条件器；不改容差和
   最大迭代数来掩盖问题。
2. 外层压力失败但内层全部收敛：验证 exact operator 的线性/重复 apply、nullspace
   投影和 compact-Schur 预条件器质量。只有 RED 证明需要时，加入固定线性的 compact
   pressure V-cycle/Schwarz 近似；不切回 compact operator 作为产品方程。
3. 独立残差失败而 solver 报告收敛：修复 solver stopping/真实残差合同，不放宽科学
   阈值。

任何预条件器必须是通用网格/系数算法，不读取 case 名、Re、圆柱位置或网格层级。

## 7. 验收层级

### Fast RED

- sparse `VectorLayout` Halo：1/2/4 rank 值、请求顺序、错误 owner、缺失 ID、collective
  failure 和 payload 计数；
- periodic surface window、donor stencil/field image 和 LFP crossing；
- active pressure/momentum 与原全局交换 oracle 在小网格上逐位或 1 ULP 一致；
- nested solve failure provenance mutation；
- signed-force 四字段和两次 PISO 回归。

### Cylinder functional gate

- 48³、64 rank、固定 `dt=0.006` 完成至少一个完整时间步；
- 无崩溃、死锁、NaN、collective mismatch 或 retry；
- 恰好两次 PISO；压力、最终通量和 `operator_force/budget_reaction/
  surface_traction/consistency` 均存在且有限；
- 输出 `Cd`、`Cl`，但单步不解释为统计结论；
- active runtime exchange 峰值随本地 ghost 数增长，不随全局 active cell 数复制。

### Benchmark startup acceptance

- 冻结 480×480×48 正式 case 的 exact binary/input/command；
- 64 rank 完成初始化和最少一个完整步，记录 time/step、RSS、Cd/Cl 和求解计数；
- 只有该门通过后才能启动长时间 Cd/Cl/St 统计；长跑不阻塞代码开发。

## 8. 版权、并发和完成边界

实现只参考 HUNDUN 自有已验收 Halo/linear interfaces 和公开数学方法，不复制 COAST、
BOFFIN、GPL 或私有研究源码。现有 COAST snapshot 只保留为独立基准证据，不进入产品。
不停止或检查研究进程，不 push、不发布。

本任务完成时生成 exact-HEAD manifest，记录 Stage 4 parent、融合 diff、二进制和日志
SHA-256、MPI/toolchain、测试命令、后台进程状态及 ACCEPT/REJECT。只有冻结圆柱功能门
和正式 startup gate 均通过，才把目标标记为完成。

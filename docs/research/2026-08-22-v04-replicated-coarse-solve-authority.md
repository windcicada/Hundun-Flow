# HUNDUN-FLOW v0.4：小型最粗网格 rank-reduced / replicated coarse solve 的公开方法 authority

检索与访问日期：2026-08-22（UTC）。

本笔记是 `RESEARCH_ONLY`。它只回答公开的 multigrid 数学、布局和生命周期思想是否存在，
以及 HUNDUN 如何把一次 coarse solve 做成可审计的独立实现；不选择 HUNDUN 架构，不冻结
candidate，不给出 `ACCEPT/REJECT`，也没有运行 solver、full-grid、COAST 或测试。

检索范围限于 PETSc、hypre、Trilinos 的官方文档、MPI Forum 的 MPI-4.1 标准以及作者公开的
原始论文。没有读取、复制或引用 COAST/GPL 实现代码；下面的链接只用于核对公开接口语义和
方法描述，不授权移植任何第三方实现。

2026-08-22 由 main agent 重新打开并核对了 PETSc `PCMGSetLevels` /
`PCREDUNDANT` / `PCGAMGSetRankReductionFactors`、hypre 官方主干公开头文件、
May 等人的原始论文和 MPI-4.1 的 collective introduction/correctness 条文。
下文使用的 rank reduction、redundant coarse solve、collective 顺序和
type/count matching 主张与这些一手资料一致。

## 结论摘要

对 HUNDUN v0.4 native Cartesian MG 的小型最粗网格，以下两类想法都属于公开且已有先例的
multigrid 生命周期设计：

1. **rank-reduced / agglomerated coarse solve**：随着网格变粗，把 coarse unknowns 和
   coarse operator 重新布局到父 communicator 的较少 active ranks，在子 communicator 上
   solve，再把 correction 传回父布局。PETSc `PCMG` 允许每个 level 使用不同 communicator，
   `PCTELESCOPE` 明确运行在单个较小子 communicator 上，PETSc 的 GAMG 还公开了按 level 的
   rank-reduction schedule。见 [PCMGSetLevels](https://petsc.org/release//manualpages/PC/PCMGSetLevels/)、
   [PCTELESCOPE](https://petsc.org/release/manualpages/PC/PCTELESCOPE/) 和
   [PCGAMGSetRankReductionFactors](https://petsc.org/release/manualpages/PC/PCGAMGSetRankReductionFactors/)。

2. **replicated / redundant coarse solve**：把相同的 coarse problem 复制到多个 rank 或
   rank subgroup，各自进行 coarse solve，再在本地使用相同 correction；或者只把整个问题
   复制到若干 subgroup。PETSc `PCREDUNDANT` 是 subgroup redundant solve 的明确先例；
   hypre BoomerAMG 公开了在阈值以下让一个 rank 顺序 solve，或在 `redundant=1` 时让所有
   remaining active ranks 各自进行顺序 AMG 的选项。见 [PCREDUNDANT](https://petsc.org/release/manualpages/PC/PCREDUNDANT/)
   和 [hypre 官方 ParCSR 公开头文件](https://github.com/hypre-space/hypre/blob/master/src/parcsr_ls/HYPRE_parcsr_ls.h#L703-L721)。

这两个概念不能混为一谈：一个 active 子 communicator 上的单次 solve 是 rank reduction；
多个 subgroup 或多个 active rank 上的独立相同 solve 才是 redundancy。两者可以组成混合
方案，但本笔记不替 root 选择其中任何一种。

“公开方法思想”不等于“可以复制实现”。HUNDUN 必须独立实现自己的 Cartesian index map、
stencil/coarse operator、边界和 null-space 规则、MPI 调度、资源生命周期、错误传播和证据
计数；不能从 PETSc、hypre、Trilinos 或 COAST 旧实现复制代码。

## 1. 何时 coarse-grid communication/scalability 成为瓶颈

在顺序 multigrid 中，每级 unknown 数量随 coarsening 下降，因此计算量通常下降；在分布式
内存中，粗化后每个 rank 的计算量下降得更快，而通信 latency、global reduction、同步和
数据重分布的相对占比上升。May 等人的原始论文明确指出，越粗的 level，communication 相对
computation 越重要；一个很小的问题若仍分布在许多 MPI ranks 上，甚至 exact LU 的成本也
不再可忽略。[May et al., *Extreme-scale Multigrid Components within PETSc*, §1](https://arxiv.org/abs/1604.07163)

该论文列出的实际触发条件包括：通信无法与该 level 的计算重叠，或出现少于约一个 unknown
每 rank；此时可以截断层次、使用 inexact coarse solver，或把 unknowns 聚合到较少的
communicator。论文同时强调，理想策略依赖问题和机器，尤其是网络 latency、global reduction
成本和浮点吞吐的相对关系；不能把一个固定 rank factor 当作通用规律。[May et al., §1](https://arxiv.org/abs/1604.07163)

因此 HUNDUN 后续若考察 coarse rank reduction，触发条件应先作为可观测量冻结，而不是凭
“网格很小”判断。至少应记录：

* 每个 level 的全局/局部 active unknown 数和 unknown-per-active-rank；
* coarse operator apply、restriction、RHS gather/replication、coarse solve、correction
  scatter/replication、prolongation 的 max-rank wall time；
* collective 次数、参与 communicator 大小、发送/接收 payload bytes、等待/idle 时间和
  rank imbalance；
* setup/repartition 的 cold cost 及其对实际 solve 次数的摊销；
* 同一数学工作量下，rank-reduced、replicated 和当前全 communicator 路线的数值结果与
  iteration/convergence 记录。

PETSc 的原始 rank-agglomeration 论文给出了这一类设计的公开 performance 依据：coarse
communicator 可以逐级缩小，数据在 communicator 之间搬运，且在大规模算例中 repartitioning
可显著改善 time-to-solution；但它也指出 repartition setup/application cost 与 reduction
factor 有关，收益是问题和机器相关的。[May et al., §§1, 4](https://arxiv.org/abs/1604.07163)

## 2. 标准先例：rank-reduced 与 replicated/redundant

### 2.1 Rank-reduced / agglomerated

PETSc 的 `PCMGSetLevels()` 接口允许每个 multigrid level 提供独立 `MPI_Comm`，并明确允许
coarse level 使用比 fine level 更少的 processes；不参与 coarse solve 的 rank 以
`MPI_COMM_NULL` 表示。官方说明还要求 transfer 设计显式处理两个 communicator 之间的
布局变化：可以先在较大的 communicator 上做标准 restriction/interpolation，再用 MPI
搬运到较小/较大的布局。[PETSc `PCMGSetLevels`](https://petsc.org/release//manualpages/PC/PCMGSetLevels/)

`PCTELESCOPE` 是更直接的单子 communicator 先例。它把 coarse problem 的 RHS scatter 到
子 communicator 的 vector，在该子 communicator 上调用内部 KSP；未参与的父 communicator
ranks 在 solve 期间 idle。若 coarse `DM` 提供了子 communicator，PETSc 会检查它确实是父
communicator 的 sub-communicator，并使用它来定义 coarse solve 的参与者。[PETSc `PCTELESCOPE`](https://petsc.org/release/manualpages/PC/PCTELESCOPE/)
、[PCTelescopeSetUseCoarseDM](https://petsc.org/release/manualpages/PC/PCTelescopeSetUseCoarseDM/)

PETSc `PCGAMG` 公开了两种与此相同层次的控制：按每 rank equation 数量自动减少 coarse
processes 的 [PCGAMGSetProcEqLim](https://petsc.org/release//manualpages/PC/PCGAMGSetProcEqLim/)，
以及用户给定的逐 level reduction schedule [PCGAMGSetRankReductionFactors](https://petsc.org/release/manualpages/PC/PCGAMGSetRankReductionFactors/)。
`PCGAMGSetRepartition` 则明确说明，在减少 coarse MPI ranks 时重分配 degrees of freedom
可以改善各级负载，但会增加 setup cost。[PCGAMGSetRepartition](https://petsc.org/release//manualpages/PC/PCGAMGSetRepartition/)

Trilinos MueLu 官方用户指南也把 coarse-grid max size、coarse solver 和 coarse-level
repartition 作为参数化生命周期的一部分；其 `repartition: start level` 与
`repartition: min rows per proc` 用于从某个 coarse level 开始重新平衡。见
[MueLu User's Guide, §§5.4/Appendix B.2.4](https://trilinos.github.io/pdfs/mueluguide.pdf)。

### 2.2 Replicated / redundant

PETSc `PCREDUNDANT` 的公开定义是：把整个问题的 KSP solve 放到 processor subgroups 中，
例如 64 ranks 分成 4 个各含 16 ranks 的并行 solve。它证明“复制 coarse problem/solver 到
多个 subgroup”是标准 preconditioner 组织方式；但它的 subgroup 语义不自动等同于 HUNDUN
“每个 rank 都持有同一 coarse vector”，具体复制粒度仍需独立定义。[PETSc `PCREDUNDANT`](https://petsc.org/release/manualpages/PC/PCREDUNDANT/)

hypre 的 BoomerAMG API 给出更接近“小型最粗网格”的先例：
`HYPRE_BoomerAMGSetSeqThreshold` 规定阈值以下可让 process 0 做 sequential AMG；
`HYPRE_BoomerAMGSetRedundant` 与该阈值配合时，`redundant=1` 则让所有 remaining active
processes 各自做 sequential AMG。官方头文件将 `redundant` 默认值明确为关闭。
[hypre `HYPRE_BoomerAMGSetSeqThreshold` / `SetRedundant`](https://github.com/hypre-space/hypre/blob/master/src/parcsr_ls/HYPRE_parcsr_ls.h#L703-L721)

所以公开方法边界可以这样表述：

| 组织方式 | coarse operator/RHS 所在布局 | solve 执行者 | correction 返回方式 | HUNDUN 可借鉴的公开思想 |
|---|---|---|---|---|
| 全 communicator | 每个 active rank 持有其分块 | 全 communicator | 原有 prolongation | 基线 |
| rank-reduced | 聚合/重分区到 `C_c`，`|C_c| < |C_f|` | `C_c` 上一次 solve | `C_c -> C_f` scatter，再 prolongate | communicator hierarchy、布局变换、idle-rank 语义 |
| replicated subgroup | 每个 subgroup 有一份相同 coarse problem | 每个 subgroup 一次独立 solve | subgroup 内直接使用或选定 root 汇总 | redundant subgroup 语义 |
| replicated active ranks | 每个 active rank 有相同 coarse problem | 每个 rank 独立 solve | 各自直接 prolongate，或显式一致性检查 | sequential/redundant coarse solve 语义 |

表中的 HUNDUN 行只是术语和生命周期映射，不是架构选择。尤其是“相同 RHS”必须由确定的
gather/replicate 规则产生；不能以各 rank 的局部 RHS 未经核对地代替全局 coarse RHS。

## 3. HUNDUN 的最小安全方法边界

本节是结合公开先例和 HUNDUN 当前约束作出的实现安全要求；其中关于 HUNDUN interface、
allocation 和边界策略的内容是 HUNDUN 规范要求，不是 PETSc/hypre/Trilinos 的兼容性承诺。

### 3.1 不扩大原 MG interface

coarse solve 只应是现有 MG cycle 内部的 layout/lifecycle 变体，保留原有外部 MG interface
和 `C1 -> transport -> PISO1 -> C2 -> PISO2` seam。不要为了 expose communicator、root、
replication factor 或 staging buffers 而扩大 public API；这些应是内部 frozen policy/state。

对任意选定的 coarse layout，审计对象仍是同一个离散 correction contract：

\[
  r_c = R r_f,\qquad A_c e_c = r_c\quad\text{(精确或声明的 inexact solve)},\qquad
  x_f \leftarrow x_f + P e_c.
\]

若采用 Galerkin operator，公开 PETSc `PCMG` 文档给出的形式为
\(A_c = R A_f P\)；HUNDUN 的 stencil、系数、边界和数值顺序仍须独立实现并由自己的
oracle 验证。[PETSc `PCMG`](https://petsc.org/main/manualpages/PC/PCMG/)

### 3.2 冷/热生命周期

| 阶段 | 允许做的事 | 必须固定/可审计的状态 |
|---|---|---|
| cold setup | 决定是否 reduction/replication、active-rank schedule、父子 communicator、global-to-local map、coarse operator 结构、RHS/correction staging buffers、solver workspace、activity/null-space metadata | communicator 成员和 rank 顺序、global coarse index 顺序、root/replica 关系、counts/displacements、buffer capacity、operator structure hash |
| coefficient/BC refresh | 在系数、时间步或边界状态改变时刷新 coarse operator 的数值条目，并同步更新 boundary/activity/null-space metadata | 结构不变时不得隐式改变 map 或 communicator；刷新完成后才允许下一次 apply；若结构改变必须重新 cold setup |
| hot MG apply | 按固定顺序 local restriction → RHS gather/replicate → coarse solve → correction scatter/replicate → prolongation/add；只使用已有 workspace | 不创建/销毁 communicator，不构造 map，不做 heap allocation；所有 counts、rank 参与者和 collective ordinal 固定 |
| teardown/reconfigure | 释放 cold-owned communicator、operator、workspace 和 maps | 只在 cycle 外释放；不能把一次 apply 的临时 ownership 泄漏到下一次 apply |

`A_c` 的 coefficient refresh 与每次 apply 的 `r_c` 更新是两个生命周期：前者刷新 operator
数值和约束，后者每次由当前 fine residual 产生 RHS。复制方案不能以“上次 RHS”替代当前
RHS，也不能让 inactive rank 的旧 buffer 参与 correction。

### 3.3 RHS gather、solve、correction scatter/replication

建议把三段数据移动写成可编号的内部事件，而不是散落在 smoother/operator 代码中：

1. **restriction**：按已冻结的 global coarse ordering 计算 `r_c`；记录每个 owner 的
   `(global index, value)` 或等价确定性 packed layout。
2. **RHS gather/replicate**：rank-reduced 路线将分块 RHS gather 到 `C_c`；replicated 路线
   将同一 packed RHS 复制到每个 subgroup/rank。root、counts、displacements、datatype
   和 payload bytes 必须是证据的一部分。
3. **coarse solve**：在声明的 active participants 上执行相同的 coarse operator/solver
   policy；若是多个独立 replica，必须说明它们是数学上独立的相同 solve，而不是隐式平均。
4. **correction scatter/replicate**：把 `e_c` 按相同 global ordering 返回 fine owners；
   未持有 coarse unknown 的 fine rank 仍按父 transfer 规则得到零或相应 correction。
5. **prolongation/add**：保持原 MG 的 prolongation、fixed-DOF mask 和 additive 更新语义；
   不在 scatter 阶段偷偷改变数值 scaling。

### 3.4 Dirichlet、Neumann、periodic、null-space/activity

这些规则属于 HUNDUN 的 scientific-work contract，不能因换 communicator 而弱化：

* **Dirichlet/fixed DOF**：coarse restriction 不得把被消元/固定变量当作自由未知量；coarse
  correction 返回后 fixed DOF 必须保持既定值。若使用 row/column elimination，coarse
  operator、RHS 和 correction mask 必须遵循同一约定。
* **Neumann/pressure null-space**：保留兼容性条件和既定 null-space/projector 规则；若
  RHS 需要投影或 gauge fixing，必须在 replicated/rank-reduced 布局中执行同一数学操作，
  不能只在 root 上隐式 pin 一个值。
* **periodic**：global periodic neighbor mapping、wrap sign/weight 和 coarse global
  ordering 必须在父、子 communicator 间保持一致；local partition boundary 不能被误当作
  physical boundary。
* **activity/masks**：inactive cells、空 rank、切割边界和零贡献 DOF 必须在 operator、RHS、
  correction 三个阶段一致传播；复制的 operator/RHS 必须带有相同 active metadata。
* **odd/even 与退化层**：当 coarse extent 小于 active rank 数、某方向为 odd extent、或
  periodic wrap 退化时，必须走已验证的 generic mapping；不能依赖“每 rank 至少一个 cell”
  的隐式前提。

### 3.5 资源边界

`MPI_Comm_split`/`MPI_Comm_dup`、global maps、packing metadata、operator replicas、solver
factor/workspace 和 staging buffers 都属于 cold allocation。hot MG apply 只能复用这些
容量；任何 fallback 也必须提前分配其最大所需 workspace。若某一路线需要在 apply 中重新
分配或构造 child communicator，应视为违反 HUNDUN v0.4 性能政策，而不是“少量例外”。

## 4. MPI collective ordering、failure 与 accounting

### 4.1 MPI 标准约束

MPI-4.1 要求同一个 communicator 上的 collective initialization 按所有成员相同顺序
执行；collective 的参与者集合由 communicator 定义。[MPI-4.1, §2.4/§7.1](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node30.htm)
、[MPI-4.1, §7.2](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node114.htm)

MPI-4.1 的 correctness 章节进一步说明，collective 必须在所有成员中以相同顺序调用，
不同 overlapping communicators 之间还不能形成 cyclic dependency；blocking 与 nonblocking
collective 也不能互相匹配。所有 nonblocking collective 也必须按 communicator 保持顺序，
并在完成前保持输入/输出 buffer 的访问约束。[MPI-4.1, §7.14](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node172.htm)
、[MPI-4.1, §7.12](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node145.htm)

MPI 还要求 collective 接收方声明的数据量与发送方实际发送量匹配；`Gatherv`/`Allgatherv`
的 counts、displacements、datatype extent 和 packed global ordering 因此必须一起冻结，
不能只比较一个总 bytes 数。[MPI-4.1, §7.1](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node114.htm)

对 HUNDUN 的直接约束是：

* 父 communicator 上的 transfer/状态 collective 必须由全部规定参与者按固定 ordinal 调用；
  child coarse communicator 上的 solve/gather/scatter 只由该 child 的成员调用，parent 上的
  inactive rank 不能误调用 child collective；
* 每个 level、每个 MG apply 的 collective schedule 必须确定；若多个 communicator 重叠，
  要么固定无环顺序，要么采用已完成且有证据的 nonblocking 设计；
* `MPI_COMM_NULL` 只能表示该 rank 不属于某个 level communicator，不能作为可调用的
  communicator；
* 新增 collective 必须进入 evidence，至少包含 `ordinal/type/communicator/participants/`
  `root/counts/displacements/datatype/payload bytes`，并能与 old route 的 product/work
  counters 对账。

### 4.2 MPI failure 不是自动容错

MPI-4.1 规定，默认 predefined communicator 的 error handler 通常是
`MPI_ERRORS_ARE_FATAL`；`MPI_ERRORS_RETURN` 的作用是把错误码返回给用户，并不保证
后续 MPI 调用可以继续。MPI 标准本身也不提供 MPI process 意外永久停止后的恢复机制；某些
错误可能使后续 API 继续失败，用户必须决定终止还是使用明确设计的恢复路径。[MPI-4.1, §10.3](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node252.htm)
、[MPI-4.1, §3.8](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node47.htm)

因此 HUNDUN 的最小安全规则应是：

1. 每个 MPI return/error status 绑定到 `level/apply/collective ordinal`；
2. 任意 child gather/solve/scatter 出错时，不能将部分 `e_c` 写回 fine solution；父级应得到
   一致的 failure status，并按既定 policy 终止本次 apply 或整个 run；
3. 若启用 `MPI_ERRORS_RETURN` 仅用于诊断，也必须记录它不提供 process-failure recovery；
4. 不能以“某个 rank 先返回”作为成功判据。成功必须包括 active participants 的 collective
   completion、数值 status、correction checksum/residual 和资源审计。

### 4.3 建议的 HUNDUN evidence record

每一次 coarse transfer/solve 至少记录以下字段，供 root 自行纳入 frozen candidate 的证据矩阵：

```text
run_id, build/tree identity, level, mg_apply_id, phase_id
parent_comm_size, child_comm_size, active-rank list and global-to-local rank map
collective ordinal, collective type, communicator identity, root/replica group
datatype, counts, displacements, send bytes, receive bytes
coarse operator structure/value identity, RHS/correction checksum or exact oracle identity
max/min rank wall time, wait/idle time, iteration/status, returned error code
allocation count/bytes before and during apply, heap-allocation violation flag
```

这套字段是 HUNDUN 的审计建议，不是 MPI 标准定义；关键是把“新 collective 数量/类型/顺序/
参与者/bytes”与 numerical result、iteration、message/collective counters 和 hot-path heap
policy 放在同一份 evidence 中。只有这样，rank reduction 带来的通信变化才不会被误报成纯
solver speedup。

## 5. 建议的数值与资源验证矩阵

下表是后续实验建议，不替 root 决定路线、candidate freeze 或最终接受判定。

| 轴 | 最小覆盖 | 需要保存的证据 |
|---|---|---|
| communicator schedule | baseline 全 ranks；固定 factor 的 rank-reduced；single active rank；replicated subgroup；replicated active ranks；多级 reduction | 每 level 的 parent/child size、成员/rank map、collective schedule、idle/imbalance |
| Cartesian geometry | 1/2/4 ranks；even 与 odd global/local extents；coarse unknowns 少于 ranks；空或 inactive rank；各方向周期 wrap | global index map、owner map、restriction/prolongation oracle、periodic neighbor checksum |
| boundary/null-space | Dirichlet；pure Neumann/pressure null-space；mixed BC；periodic；activity mask/切割边界 | coarse compatibility/projector、fixed-DOF invariant、null-space residual、active mask identity |
| coefficients | 常系数；空间变系数；时间步/系数 refresh 后重复 apply；结构不变与结构改变两类 | `A_c` value/structure identity、refresh ordinal、冷重建与热复用的区别 |
| numerical equivalence | 与全 communicator baseline 比较 `r_c`、`e_c`、prolongated correction、post-cycle residual；相同 arithmetic order 时要求 bitwise，否则声明 tolerance | exact checksum 或误差范数、residual/convergence/iteration、nonfinite/status |
| solver policy | direct/iterative/inexact coarse policy 各自单独登记；不把不同 tolerance 或 iteration count 混成同一工作量 | coarse solver type、tolerance、iterations、convergence status、failure path |
| resources/performance | setup amortization；apply max-rank time；gather/solve/scatter 分项；collective bytes/count；hot-path allocation | wall-time breakdown、messages/collectives/bytes、heap/allocation counter、rank imbalance |
| failure/accounting | collective return error、coarse nonconvergence、replica mismatch、invalid map/shape；必要时只做受控 fault-injection | error code/stage、no-partial-correction invariant、全 rank status、可恢复性声明 |

最低限度应先证明：layout 变换不改变离散 coarse work；`r_c` 和 `e_c` 的 owner/global ordering
可复核；post-cycle residual 与 baseline 处于预先声明的等价范围；新增通信和内存资源完全
计数；coarse failure 不污染 fine state。性能是否足以进入 candidate、是否启动 formal
pairing 或长统计，仍由 root 按 governance 文档另行判断。

## 6. Primary-source receipt

| 来源 | 用于支持的关键主张 | 访问日期 |
|---|---|---|
| [May, Sanan, Rupp, Knepley, Smith, *Extreme-scale Multigrid Components within PETSc*](https://arxiv.org/abs/1604.07163) | 粗 level 上 communication/computation 比例上升；少于约一个 unknown/rank 或通信不能重叠时需要 truncation/agglomeration；rank agglomeration 与问题/机器相关 | 2026-08-22 |
| [PETSc `PCMGSetLevels`](https://petsc.org/release//manualpages/PC/PCMGSetLevels/) | 每 level 独立 communicator、`MPI_COMM_NULL` 非参与 rank、较大 communicator 上 transfer 后再搬运到较小布局 | 2026-08-22 |
| [PETSc `PCTELESCOPE`](https://petsc.org/release/manualpages/PC/PCTELESCOPE/) 与 [coarse-DM setup](https://petsc.org/release/manualpages/PC/PCTelescopeSetUseCoarseDM/) | 单个 reduced sub-communicator、RHS scatter、inactive parent ranks、coarse DM communicator hierarchy | 2026-08-22 |
| [PETSc `PCREDUNDANT`](https://petsc.org/release/manualpages/PC/PCREDUNDANT/) | processor subgroups 上的 redundant solves | 2026-08-22 |
| [PETSc `PCGAMGSetProcEqLim`](https://petsc.org/release//manualpages/PC/PCGAMGSetProcEqLim/)、[`PCGAMGSetRankReductionFactors`](https://petsc.org/release/manualpages/PC/PCGAMGSetRankReductionFactors/)、[`PCGAMGSetRepartition`](https://petsc.org/release//manualpages/PC/PCGAMGSetRepartition/) | coarse level 按 equation/rank 或 schedule 减少 ranks，并可重分配以改善负载 | 2026-08-22 |
| [hypre `HYPRE_BoomerAMGSetSeqThreshold` / `SetRedundant`](https://hypre.readthedocs.io/en/latest/api-sol-parcsr.html) | 小 coarse system 的 process-0 sequential solve，或 `redundant=1` 时所有 remaining active processes 的 sequential AMG | 2026-08-22 |
| [Trilinos MueLu User's Guide](https://trilinos.github.io/pdfs/mueluguide.pdf) | coarse max size、coarse solver 与从指定 level 开始的 repartition/min rows per proc | 2026-08-22 |
| [MPI-4.1 §2.4/§7.1](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node30.htm) 与 [§7.2](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node114.htm) | collective ordering、communicator participant scope、collective data matching | 2026-08-22 |
| [MPI-4.1 §7.14 correctness](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node172.htm) 与 [§7.12 nonblocking collectives](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node145.htm) | 同 communicator 同序、overlapping communicator deadlock/cyclic ordering、blocking/nonblocking matching 与 buffer lifetime | 2026-08-22 |
| [MPI-4.1 §10.3 error handling](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node252.htm) 与 [§3.8](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/node47.htm) | default fatal / `MPI_ERRORS_RETURN` 语义、后续调用可能继续失败、MPI 不提供 process-failure recovery | 2026-08-22 |

本 receipt 只证明公开方法先例和可审计边界存在；它不证明任何具体 HUNDUN 路线已经数值
等价、性能达标或可以进入 immutable candidate。

# HUNDUN-FLOW Stage 2 变密度流动设计规格

日期：2026-07-19

状态：本设计的方向与章节已经用户逐项批准；它只用于书面复核，尚不授权实现。
用户书面复核本规格后，必须另行编制并批准 Stage 2 实施计划，才能改动产品代码。

本规格承接已验收的 Stage 1，不修改 `stage1-runtime` 的接口、格式或科学声明。Stage 1 的现状与 Stage 2 强制规划顺序见 [Stage 1 runtime architecture](../../architecture/runtime.md)。

## 1. 已批准决定摘要

- Stage 2 采用低马赫、压力约束的变密度 Navier--Stokes 路径，不实现带声波、激波和 Riemann 求解器的全可压缩路径。
- Stage 2 的压力--速度耦合核心为独立实现的 PISO 与 Rhie--Chow；每个成功时间步固定执行两次 pressure corrector。
- SIMPLE、SIMPLEC 和稳态伪时间驱动不进入 Stage 2；未来可作为独立控制器增加，而不是改变本阶段方程算子。
- Stage 2 依次通过定密度、材料变密度、最小单组分焓--密度反馈三个数值门。
- 最小热膨胀闭合采用常 `cp`、常 `R`、`h -> T -> rho`，无物种和化学反应。
- Stage 2 同时求解质量、三分量动量、机械压力、焓和一个或多个通用标量。
- 质量、动量、焓与所有标量共享同一份经 PISO 修正的面质量通量。
- Stage 2 的开放边界只增加规定速度入口与机械压力出口；出口回流导致整步一致拒绝，不截断负流量，也不猜测回流热化学状态。
- 网格职责分离为 `MeshTopology`、`BoundaryPatch` 和 `MeshGeometry`。
- 字段安全分离为带 epoch/capability 的 checked view 和短生命周期 kernel view。
- 项目自有的后端中立线性代数与执行合同必须先于压力流实现；CPU reference 是本阶段数值基线，device 只保留接口语义，不实现生产 GPU 后端。

上述决定细化了总体设计中的 Stage 2 变密度流动范围，并落实了
[Stage 1 runtime architecture](../../architecture/runtime.md) 中记录的 Stage 2
强制前置项；总体阶段边界仍以[总体设计第 7、15、16 节](../../design/2026-07-16-hundun-flow-clean-cpp-solver-design.md)为准。

## 2. 范围与非目标

Stage 2 的交付范围是一个可独立验证的、瞬态、低马赫、变密度有限体积流动求解器。它在 Stage 1 的 C++17/MPI-3 运行时上建立新的组合路径，并保持 Stage 1 被动标量可执行程序作为独立回归路径。

本阶段包含：

- 质量守恒方程和三分量动量方程；
- 机械压力及压力修正方程；
- 显热焓和可配置数量的通用保守标量输运；
- 定密度、材料变密度和最小单组分理想气体热膨胀闭合；
- cell-centred 有限体积、PISO、Rhie--Chow 和自适应 BDF2；
- uniform adapter 与结构化曲线网格度量；
- Stage 2 已批准的 body-fitted 基础边界；
- 项目自有 CPU reference 线性求解与 matrix-free 压力算子。

本阶段明确不包含：

- 全可压缩、声学、激波、Riemann、WENO 或 DG；
- 物种方程、有限速率化学、反应热和多组分热力学；
- LES、IBM、TPDF-TCR、喷雾、颗粒和运动壁面；
- 多 STL、多部件、AMR、非结构网格或 rank-changing Restart；
- 生产 CUDA/HIP kernel、GPU runtime、PETSc、HYPRE、AMG 或 Kokkos；
- 复杂入口/出口选择器、时变运行日程和热化学回流状态。

Stage 2 不改变版权与独立性边界。BOFFIN 固定基线只可用于仓库外的私有独立性审计，COAST 与 COAST-2 只可作为黑盒科学需求参考，均不是源码上游。不得修改、清理、停止、打包或提交 `/home/wyf/code_dev/Coast_software` 及其研究算例、运行数据和进程。公开构建及运行仍不得依赖 Python，也不得发布或 push，除非用户以后另行明确授权。完整政策见 [source policy](../../development/source-policy.md) 和 [Stage 0/1 计划全局约束](../../plans/2026-07-16-hundun-flow-stage0-stage1.md)。

## 3. 控制方程与密度闭合

Stage 2 求解低马赫质量守恒：

```text
partial(rho)/partial(t) + div(rho*u) = 0
```

三分量动量使用同一守恒形式：

```text
partial(rho*u)/partial(t) + div(rho*u tensor u)
  = -grad(pi) + div(tau) + f
```

其中 `pi` 是维持连续性约束的机械压力，不是负责声学传播的全可压缩热力学压力。应力 `tau` 在 Stage 2 使用明确的层流黏性闭合，`f` 只接受本阶段已声明且可守恒审计的体源项。

显热焓采用保守输运：

```text
partial(rho*h)/partial(t) + div(mdot_f*h_f)
  = div(Gamma_h*grad(h)) + S_h
```

每个通用标量 `phi_s` 采用相同结构：

```text
partial(rho*phi_s)/partial(t) + div(mdot_f*phi_s,f)
  = div(Gamma_s*grad(phi_s)) + S_s
```

这里 `mdot_f` 必须是本时间步最终由 PISO 修正并通过连续性验收的唯一面质量通量。各方程不得维护彼此独立的面通量副本或用未修正预测通量提交新状态。

密度闭合按以下顺序引入并分别验收：

1. 定密度门：`rho = rho_ref`，先隔离验证压力--速度耦合和时间离散。
2. 材料变密度门：`rho` 由质量方程保守推进，焓暂不改变密度。
3. 最小热膨胀门：单组分、常 `cp`、常 `R`，由 `h` 得到 `T`，再由
   `rho = p0/(R*T)` 得到密度。

闭域热膨胀时，空间均匀的热力学压力 `p0` 必须作为全局状态更新，使闭域总质量与热力学闭合一致。`p0` 与机械压力 `pi` 职责分离，不能用压力修正场替代热力学状态。

开放域热膨胀的 `p0` 必须由 resolved case 配置提供，并且是正、有限、时间不变
的热力学参考量；机械压力出口只约束 `pi`，不能指定或更新 `p0`。最小热膨胀
入口必须在 `h` 与 `T` 中恰好选择一个权威热状态：选择 `h` 时按 `T = h/cp`
推导，选择 `T` 时按 `h = cp*T` 推导，两者都按 `rho = p0/(R*T)` 得到入口
密度。若配置还冗余给出其余 `h`、`T` 或 `rho`，这些值只能用于按同一单位和
配置容差严格交叉校验，不能按字段优先级覆盖推导值；任何不一致必须在运行前
拒绝。

本阶段不求解物种，不定义 `cp(T,Y)`、混合物气体常数、反应热或真实气体状态方程。多组分热化学密度反馈属于 Stage 4。非正或非有限的 `rho`、`T` 必须使全部 rank 一致拒绝本次试算，不得静默截断为下限后继续提交。

## 4. PISO、SIMPLE 与长期耦合策略

PISO 是本项目瞬态压力--速度耦合的内核，不是另一套控制方程。一次 corrector 根据临时面质量不平衡构造压力修正，随后一致修正机械压力、速度和面质量通量。Rhie--Chow 使用实际动量方程对角系数形成 collocated 面通量，抑制压力棋盘格。

Stage 2 的每个成功时间步固定执行两次 pressure corrector。固定次数形成可重复
的 reference 控制流，便于验证时间精度、守恒、MPI 分解不变性和 collective
预算；本阶段不在 corrector 外加入残差驱动的热化学非线性循环。

SIMPLE 及其相关稳态算法不作为 Stage 2 主路径。它们更适合稳态求解、伪时间
初始化或允许大量外迭代与欠松弛的工作流，今后可在相同算子和残差合同上作为
独立控制器实现，不得迫使基础算子依赖某一种耦合算法。

Stage 4 在焓、组分、热物性、密度与速度压力形成强非线性反馈时，可在 PISO
外增加残差驱动的热化学外层迭代。长期方向因此是“PISO 瞬态内核加可选外层
耦合”，而不是宣称 PISO 在所有问题上无条件优于 SIMPLE。

PISO 与 Rhie--Chow 的公开科学来源已经在
[scientific sources](../../references/scientific-sources.md) 中登记；旧程序的
控制流、消息和数组更新顺序不得作为实现来源。

## 5. 数值基线

- 空间离散为 collocated、cell-centred、二阶守恒有限体积。
- 时间推进首个有效步使用 backward Euler，具备历史状态后使用固定步长
  BDF2；固定步长路径通过后才启用 ratio-limited adaptive BDF2。
- 动量对流在共享面质量通量上使用动能友好的中心形式，并定义有界回退条件。
- 密度、焓和通用标量使用 MUSCL 重构与 MC limiter。
- 黏性和标量扩散使用二阶中心离散；曲线网格使用明确的非正交修正。
- 压力--速度耦合使用时间一致的 Rhie--Chow 与两次 PISO corrector。
- 压力算子优先采用 matrix-free 表达，CPU reference 路径使用 FP64。
- SPD 压力系统使用 CG；非对称或非正交修正路径使用 BiCGStab。
- Identity 与 Jacobi 是 Stage 2 基础预条件器，不引入外部线性代数依赖。
- 周期或全 Neumann 压力系统使用零均值投影显式处理 nullspace。

时间步采用试算事务：旧状态在整步验收前保持可恢复，候选值写入 stage fields。
线性求解不收敛、NaN/Inf、非正温度或密度、边界回流或最终残差不满足合同时，
所有 rank 对该步得到一致结论；可恢复失败缩小 `dt` 后从旧状态完整重试。

## 6. 总体依赖架构

依赖只允许向下：

```text
application / time driver
        -> flow solver / transport / boundary / density closure
        -> finite-volume operators / pressure coupling
        -> project-owned linear algebra contracts
        -> mesh / fields / execution / Halo / MPI
```

`application` 只负责类型化配置、对象组装、时间循环和输出，不包含数值 kernel。
`flow solver` 管理质量、动量、机械压力、焓和标量状态；`PisoCoupler` 只管理
压力--速度修正。`transport` 提供共享通量上的守恒输运，`boundary` 提供本阶段
边界物理，`DensityClosure` 提供密度与连续性所需的闭合量。

线性求解与执行合同必须使用项目自有类型，不得泄漏 CUDA、HIP、PETSc、
HYPRE、MPI 缓冲区或外部 allocator。Stage 2 的生产基线是 `cpu_reference`；
`cpu_optimized` 只能在计数证据定位热点后增加；device 类别仅冻结 execution
space、驻留、显式传输、event、生命周期与错误语义。CPU 与未来 device 后端
都只能替换粗粒度执行机制，不得改变项目自有的算子、残差、守恒、收敛、
分解不变性或失败合同。

## 7. Field epoch、capability 与 checked/kernel view

`FieldStorage` 拥有存活控制块和单调 generation。销毁、存储替换、网格重建、
重新分区或 Restart 读取事务进入时，旧 epoch 必须先失效；generation 回绕必须
明确拒绝。控制块的生命周期必须长于被检查的数据字节，确保陈旧 view 检查
本身没有未定义行为。

checked view 在暴露元素地址前检查 liveness、generation、数值类型、空间范围、
分量和访问能力。字段访问在组装期声明 actor、field 和 read/write 意图，冻结后
不可修改；reader 只获得 const view，同一阶段同一字段只允许一个已声明 writer。
这里的 capability 是组件间的正确性合同，不是安全权限系统。

kernel view 是独立且 trivially-copyable 的指针、extent 和 stride 借用类型，不含
epoch、边界、权限、共享所有权、虚调用或逐元素原子检查。它只能从已验证的
checked view 为一次词法受限的 kernel 调用生成，不得由 solver、Halo、I/O、
插件回调或跨步 lambda 长期保存。

kernel view 的安全边界由受限构造作用域、静态类型性质和 kernel review 建立，
而不是靠运行时 sanitizer 探测越界。测试必须静态验证其 trivially-copyable 且
不含 owner token，验证只有受限构造入口能够产生它，并比较 checked/kernel
路径的数值等价性；测试不得为验证失败行为而解引用陈旧或越界 kernel view。

checked view 与 kernel view 必须在所有 Stage 2 标量、分量及 ghost 布局上产生
相同数值。Stage 1 的现有 `FieldView` 仍作为冻结回归接口处理，迁移不得改变
Stage 1 被动标量数值结果或插件 ABI v1 的 metadata-only 边界。

## 8. MeshTopology、BoundaryPatch 与 MeshGeometry

`MeshTopology` 拥有逻辑 cell/face ID、owner/neighbour 连接、owned/ghost 分类、
分区无关全局 ID、周期配对和 patch membership。它组合或查询
`StructuredDecomposition`，但不拥有物理坐标、边界状态或 MPI 生命周期。

`BoundaryPatch` 只包含稳定 ID/name、拓扑 boundary-face range 或 compact set，
以及 pairing kind。内部面不属于物理 patch，rank 间分区面不得误分类为物理
边界；周期面必须具有唯一、互反的配对关系。

`MeshGeometry` 引用一份兼容 topology，拥有 cell/face centre、正 cell volume、
定向 face area vector、Jacobian 和曲线度量。共享面的两侧面积向量必须互反，
每个 cell 满足 closure，曲线映射必须通过 metric identity 与 free-stream
preservation 验收。

Stage 2 先提供 uniform adapter，精确保持 Stage 1 的 spacing、centre、volume、
local extent 与 owned box 合同；随后增加结构化曲线几何。网格层不包含入口
速度、出口压力、壁面通量、热物性、模型源项或时间日程。

## 9. Stage 2 边界与后续阶段分配

Stage 2 只实现下列 body-fitted 边界物理：

- `periodic`：由 topology 的互反配对和一致面通量实现；
- 静止无滑移壁面：不可渗透、壁面速度为零；
- `symmetry`：法向速度和相应法向通量满足对称约束；
- 规定速度入口：指定速度、密度或闭合所需热状态、焓和通用标量；
- 机械压力出口：只指定机械压力 `pi` 参考并使用一致的出流通量处理，不约束
  热力学参考量 `p0`。

规定速度入口是 Stage 2 唯一开放入口形式，不实现总质量流量控制。最小热膨胀
入口的 `h/T/rho/p0` 必须遵守第 3 节的唯一推导和冗余值交叉校验规则。机械压力
出口以第 8 步得到的最终面质量通量为准；第 9 步一旦在该通量上检测到回流，
必须报告 patch、时间步和法向质量通量证据，并由所有 rank 一致拒绝整个试算。
不得根据第 7 步的 provisional 通量作出通过或回流拒绝结论，不得把最终负通量
截为零，也不得在缺少热化学状态时临时把出口当入口。

复杂边界已经写入总体路线，不在本阶段遗漏或提前实现：Stage 3 交付静态 IBM
表面；Stage 4 交付反应气相多入口/多出口、质量流量入口、热化学回流、温度、
焓、组分及热壁面；Stage 5 交付 TPDF-TCR 边界状态；Stage 6 交付喷雾注入、
parcel 进出及静态壁碰撞；Stage 7 交付复杂选择器、优先级、冲突拒绝及时变
日程；Stage 8 冻结已经交付的 schema、文档与 SDK 示例。

运动或变形壁面、多部件运动几何、全可压缩/激波边界和 WENO/DG 专属边界
属于首版之后。各阶段必须随新增物理同步交付类型化 JSON schema、交叉校验和
unsupported-combination 拒绝测试。本节给出这些复杂边界的完整阶段分配。

## 10. 执行资源与数据驻留合同

`ExecutionContext` 是一次粗粒度数值操作所使用执行后端的唯一入口。它声明
backend identity、`ExecutionSpace::host/device`、操作的有序性、能力查询和统一
错误语义，但不暴露厂商 stream、queue、allocator 或设备指针。上下文之间不得
凭名称假定能够共享数据；跨上下文使用必须由明确兼容性合同允许。

`Buffer` 拥有一段线性存储及其 allocation identity，并且任一时刻只有一个
权威驻留位置。host 和未来 device 副本不能同时被视为最新状态。所有权转移、
重分配及销毁都必须等待引用该 allocation 的未完成 event，或把保持 allocation
存活的责任显式转移给 event；调用方不能在异步操作尚未完成时释放底层存储。

`VectorView` 是对 `Buffer` 子区间的短生命周期借用，至少携带：

- backend identity 与 execution space；
- scalar format、元素数和 stride；
- allocation identity、offset 和 epoch；
- const/read-write capability。

view 不拥有数据，也不能改变权威驻留位置。backend、scalar、范围、epoch 或
写能力不匹配必须在提交操作前失败，不能依靠 kernel 内部猜测或隐式修复。

`ExecutionEvent` 表示一次已提交操作的完成状态。`ready()` 只能无阻塞查询，
`wait()` 必须传播该操作的确定结果；重复等待具有相同语义。event 在完成前保持
其全部输入输出 allocation 存活。host reference 的同步操作返回已经完成的
inline event，从而让同一调用路径不需要以空句柄或特殊分支表示同步执行。

`transfer(source, destination, context)` 是改变权威驻留位置的唯一公共操作，
返回 `ExecutionEvent`。框架不提供隐式 host/device coherence、不在算子调用中
偷偷复制，也不允许 host 代码直接解引用 device view。需要回迁时必须显式
transfer 并等待其 event。

Halo 数据交换路径必须由运行时能力查询选择，而不能只凭编译宏推断：host 数据
使用 `host-direct`；未来 device 数据只有在运行时确认该 context、buffer 与 MPI
路径共同支持时才可使用 `device-direct`；否则必须选择显式
`device-host-staged` fallback，在通信前后分别提交 transfer 并按 event 依赖顺序
等待。若 direct 与 staged 路径都不满足能力、驻留或生命周期合同，则在开始 Halo
前明确拒绝，不得尝试隐式复制或未经确认的直接路径。

Stage 2 只交付生产 `cpu_reference` context。device context 只冻结上述类型、
能力拒绝、event 与生命周期语义，并使用不执行真实设备计算的 test double 验证
能力选择、transfer/event 顺序、staged fallback、无可用路径时的拒绝及错误传播。
本阶段不实现真实 device buffer、device-aware MPI 或生产 GPU 后端；device test
double 不能注册为可运行生产后端，也不能被性能或科学声明称为 GPU 支持。公共
合同和测试不得出现 CUDA、HIP、SYCL、PETSc、HYPRE 或厂商句柄。未来任何 CPU
或 device 生产后端都必须通过同一项目自有 FP64 residual、conservation、
convergence、decomposition-invariance 和 failure 合同。

## 11. 线性代数与 matrix-free 压力算子

线性代数接口采用项目自有的后端中立类型，并在整次算子或求解调用边界进行
运行时替换。逐 cell kernel 不使用虚调用。冻结的职责如下：

- `LinearOperator` 声明 domain/range layout、`ExecutionContext`、单调 revision，
  并提供 `apply(x, y)`；可选 `diagonal(d)` 以能力查询决定是否可用，不要求
  暴露 CSR、稀疏矩阵或外部库对象。
- `Preconditioner` 通过 `update(operator, revision)` 刷新状态，并通过
  `apply(r, z)` 产生预条件结果；revision 未匹配时不得使用旧状态。
- `LinearSolver` 只提供
  `solve(operator, preconditioner, b, x, control)`，不读取网格字段、Halo
  buffer 或厂商后端对象。
- `SolveControl` 至少包含 absolute tolerance、relative tolerance、maximum
  iterations 和独立残差重算周期；所有值在 collective 求解开始前验证。
- `SolveReport` 至少记录结束原因、迭代数、初始残差、递推残差、项目自有
  FP64 路径独立重算的最终残差，以及 matvec、preconditioner apply 和 global
  reduction 次数。

结束原因必须区分 converged、zero right-hand side、maximum iterations、
numerical breakdown、non-finite value、invalid control 和 collective failure。
递推残差只能用于迭代过程；成功声明必须由独立 FP64 最终残差满足同一
absolute/relative 合同支持。

Stage 2 的 reference 实现包括 Identity、Jacobi、CG 和 BiCGStab。CG 只用于
满足 SPD 合同的算子；非对称离散或显式非正交处理使用 BiCGStab。所有算法使用
专用连续本地 workspace，不让 solver 直接迭代 `FieldStorage`，也不在每次迭代
动态创建临时字段。

最小 `VectorOps` 包含 fill、copy、scale、axpy、线性组合、norm 和
`dot_batch`。`dot_batch` 将一次迭代中相互独立的内积合并到一个明确的
collective 批次，并在 `SolveReport` 中如实计数；Stage 2 不以隐藏同步或改变
数学停止条件来减少计数。

压力 `LinearOperator` 默认使用 matrix-free 表达。它拥有 cell-to-vector 映射、
边界贡献描述和可复用 ghost workspace，但不拥有压力或速度主字段。一次
`apply()` 的固定数据流是：

```text
Halo begin
-> 计算不依赖远端值的 interior cells
-> Halo wait
-> 计算 partition-boundary cells
```

operator revision 在系数、geometry、boundary 或分区相关 layout 改变时递增。
Jacobi 对角线和 Halo workspace 只能在 revision 与 layout 匹配时复用。

周期或全 Neumann 压力问题显式声明常数 nullspace。右端先投影到兼容子空间，
迭代向量及最终机械压力采用零均值规范化；不通过任意固定一个 cell 隐式改变
分区不变性。开放机械压力出口提供压力参考时，算子不得再施加重复零均值约束。

未来混合精度以 `ResidualCorrectionSolver` 包装器加入：低精度后端只计算修正量，
外层在项目自有 FP64 `VectorOps` 中重算残差并决定收敛。Stage 2 不实现该包装器
或低精度生产路径，但现有接口不得要求 operator、preconditioner 和 solver 使用
相同内部标量格式之外的厂商类型。

## 12. 时间步试算事务

一个 Stage 2 时间步严格执行以下九步，旧提交状态在第九步成功前保持不变：

1. 根据 CFL、扩散限制、上一步收敛证据和重试状态选择 `dt`。第一个有效步用
   backward Euler；有完整历史后使用 BDF2，固定步长验证通过后才允许受限步长比。
2. 从已经提交的历史状态和最终修正面质量通量构造二阶历史通量预测；不得使用
   未验收的旧 predictor flux 作为新步历史。
3. 在 stage fields 中预测 `rho`、`rho*h` 和每个 `rho*phi`，三者使用同一预测
   面质量通量及其相同时间层。
4. 执行当前 `DensityClosure`：材料变密度直接使用守恒预测；热膨胀路径执行
   `h -> T -> rho`，闭域还更新候选 `p0`，随后检查有限性、正性和总质量。
5. 由候选密度、黏性、边界和历史项组装三分量动量方程，求得动量 predictor
   与实际动量对角系数；任一分量失败都使整步失败。
6. 用实际动量对角系数进行时间一致的 Rhie--Chow 面插值，求解并应用第一次
   pressure corrector，一致更新候选机械压力、速度和面质量通量。
7. 使用第一次修正的面质量通量校正 `rho`、`rho*h` 和 `rho*phi`，刷新闭合、
   连续性右端和压力算子系数。该结果只是不可提交的 provisional iterate，既
   不能作为本步新状态验收，也不能从它再次推进一个时间步；不能只修正速度而
   提交旧输运状态。
8. 求解并应用第二次 pressure corrector，得到本步唯一候选最终机械压力、速度
   和面质量通量，并再次应用压力 nullspace 规范化。随后从本步起始提交状态及
   其历史层出发，以该最终通量重新完成 `rho`、`rho*h` 和全部 `rho*phi` 的
   保守 transport finalization，再刷新 `DensityClosure`、候选 `p0`、连续性
   右端和供验收使用的压力算子系数；这不是把第 7 步 provisional 状态再推进
   一个时间步。
9. 在 finalization 和 closure 刷新后，使用最终通量独立重算质量、三分量动量、
   焓、全部标量及各线性求解的离散残差，并检查所有物理、边界和性能合同。
   若刷新使连续性或压力残差超过合同，必须 collective 回滚并按确定规则缩小
   `dt`；只有全部 rank 一致成功后才原子提交 fields、历史、`dt`、最终面通量
   和 `p0`。

步骤 8 的 transport finalization 和 closure 刷新不会增加第三次 pressure
corrector，也不是残差驱动的热化学外循环。Stage 2 的固定 reference 路径仍然
恰好执行两次 pressure corrector；Stage 4 才可在步骤 3--8 外增加明确的热化学
非线性控制器，但不能改变这些算子的独立残差和守恒合同。

任一步的非有限值、非正 `rho/T`、线性求解失败、边界拒绝或 residual failure
都先通过 collective status 汇总，使所有 rank 得到相同失败类别与最低失败 rank
证据。任何 rank 都不得单独提交、单独缩小 `dt` 或提前进入下一次 collective。

失败时丢弃全部 stage fields、候选历史、候选 `p0` 和临时线性状态，恢复到该步
开始前的最后提交状态。可恢复失败按确定规则缩小 `dt` 后从步骤 1 完整重试；
达到配置的 minimum `dt` 或 maximum retries 时，以最后一次 collective 证据
明确终止。不可恢复的配置、layout 和数据损坏错误不参与重试。

机械压力出口在步骤 9 对最终面质量通量验收时出现回流属于边界拒绝：报告 patch
stable ID、step/time、最小法向质量通量及所在全局 face，并 collective 回滚
整个试算、缩小 `dt`。Stage 2 不以第 7 步 provisional 通量作出边界通过或拒绝
结论，不得截断最终负通量或构造未配置的回流状态。

## 13. Checkpoint v2 与诊断输出

Stage 2 新增独立、版本化的 Checkpoint v2；Stage 1 Restart v1 的字节、目录和
读取合同保持冻结。v2 至少持久化：

- `rho`、三分量动量、机械压力、焓和全部 persistent 通用标量；
- BDF2 的 `n`、`n-1` 状态、所需的历史 `dt`，以及时间积分器的 startup/history
  readiness 和当前 order 状态；
- 所有会影响 Restart 后下一次 `dt` 选择的 adaptive-controller 状态，包括当前
  limiter/hysteresis、上一步收敛或误差证据及重试相关记忆；控制器不得保留任何
  未持久化却会改变下一次选择的隐式状态；
- 最终修正面质量通量和动态闭域热膨胀状态 `p0`；固定的开放域 `p0` 不作为动态
  字段恢复，但至少必须进入 resolved-case fingerprint；
- topology、geometry、boundary 和 resolved-case fingerprint；
- 字段 schema、step、time、rank count、process grid 和每 rank owned boxes；
- 每 rank 文件名、逻辑 byte size、实际 byte size 和 CRC。

全局 manifest 给出 format/schema version、所有 rank 文件记录、共同 fingerprint
和恢复限制。各 rank 先写临时文件、完成 flush/close 与 CRC 校验，再通过一致
成功判定发布 rank 文件和 manifest；`COMPLETED` 标记必须最后写入。缺少标记、
缺 rank 文件、size/CRC/fingerprint 不一致或存在未声明 trailing bytes 时，整组
checkpoint 不可恢复。

读取采用事务语义：进入读取事务时立即使旧 checked view epoch 失效，所有文件、
schema、fingerprint 和数值范围全部验证后才替换字段；任一失败均保持原有字段
值和最后提交时间步不变。Checkpoint v2 首版只支持相同 rank 数、process grid
和 owned boxes 的恢复，并在不匹配时给出明确诊断。

Stage 1 primitive VTK 输出继续按原合同回归，不扩展其格式来承载曲线几何。
Stage 2 提供独立的结构化曲线几何诊断输出，包含足以检查 cell/face geometry、
守恒量和 partition identity 的版本元数据；其格式在实施计划中以无 Python 的
reader/writer 测试冻结。本阶段不加入 HDF5、MPI-IO 或 rank-changing Restart。

## 14. 测量、性能证据与优化准入

每次性能测量 artifact 至少记录：commit、工作树 clean/dirty 状态、
resolved-case fingerprint、编译器版本、完整 compiler/link flags、构建类型、MPI
实现、node 身份、CPU affinity、rank placement、rank/thread 数、process grid、
网格与每 rank owned cells、执行后端、数值配置、warmup steps、measured steps、
repetitions 和 correctness 结果。dirty 运行还必须保存可唯一识别该工作树差异的
摘要，不能与 clean commit 基线混同。未先通过相同残差、守恒、分解不变性及
输出合同的运行，不得作为性能基线。

wall-clock 使用 `MPI_Wtime` 包围定义清楚的 phase，排除或单独报告初始化与冷启动；
warmup 不计入测量。对 rank 数为 `p` 的每次重复 `r`，artifact 保存每个 rank 的
原始 measured-phase elapsed time 及换算后的 step-time sample `t[p,r,k]`；其中
phase elapsed time 除以 measured steps 得到该 rank 的 step time。该次重复值定义
为 maximum-rank step time `t[p,r] = max_k(t[p,r,k])`，报告值
`T_p = median_r(t[p,r])`。报告不得只保留聚合值，必须保留所有 raw per-rank
samples、重复编号和换算所需元数据，使 `t[p,r]` 与 `T_p` 可以重新计算。报告
同时保留迭代次数和工作量，避免把少做一次 corrector 或放宽容差误记为优化。

strong scaling 使用同一全局问题，定义 `S_p = T_1/T_p`、
`E_p = T_1/(p*T_p)`；weak scaling 保持每 rank 问题规模与数值配置一致，定义
`W_p = T_1/T_p`。两种 `T_1` 都必须来自与对应 `T_p` 满足下述兼容性条件的
单 rank 测量，并使用同一 warmup、measured steps、repetitions 与聚合规则。

Stage 2 必须测量并输出：

- strong scaling 与 weak scaling 的 case、rank 数和 max-rank time；
- allocated bytes per owned cell，并在可用平台记录 RSS bytes per owned cell；
- Halo payload bytes、pack/unpack bytes、message count、等待时间和有效带宽；
- collective 的 kind、count 与 logical payload bytes；
- Checkpoint/诊断 I/O 的 logical bytes、on-disk bytes、持续时间与吞吐；
- 压力和动量 `SolveReport` 的 iterations、matvec、preconditioner apply、
  reduction、递推残差与独立 FP64 最终残差。

确定性 exact counters 可作为普通 CI 硬门禁，包括 allocation bytes、Halo payload
与 message 数、collective 次数、matvec、preconditioner apply 和逻辑 I/O bytes。
这些门禁必须与具体 case、rank 数、process grid 和算法配置绑定。

wall-clock、RSS、带宽和实际文件系统吞吐只允许与兼容的同机基线比较。兼容性
要求至少覆盖硬件、MPI、编译器、构建选项、rank placement、case fingerprint、
数值容差和测量方法；任一关键元数据不匹配时结果标记为 incomparable，而不是
通过固定跨机器阈值判定回归。

`cpu_reference` 先建立正确性与工作量基线。只有计数或同机 phase 证据定位出
可重复热点后，才允许提出 `cpu_optimized` 实现任务；每项优化必须继续通过同一
FP64 残差、守恒、分解不变性、Checkpoint 和边界拒绝合同。Stage 2 不以未来
device 预留为由提前引入厂商依赖或无法测量收益的复杂缓存。

## 15. TDD 与验收矩阵

每个实现任务必须先加入能够因缺少目标行为而失败的最小测试，再实现使该测试
通过的代码。失败原因必须指向本任务合同，不能依赖尚未批准的后续模块。下列
矩阵是 Stage 2 实施计划必须逐项展开的最低验收集合：

- 字段安全：陈旧 epoch 在销毁、替换、重建、重分区和 Checkpoint 读取事务进入
  后均被拒绝；read/write capability、重复 writer、类型、范围和 generation
  回绕均明确失败；ASan/UBSan 只覆盖 stale checked view 的安全拒绝。kernel view
  通过受限构造作用域、trivially-copyable 且无 owner token 的静态性质、
  checked/kernel 数值等价和逐 kernel review 验收；测试不得解引用 stale 或
  out-of-bounds kernel view。
- 网格拓扑：cell/face owner-neighbour、owned/ghost、全局 ID 和分区面分类一致；
  `BoundaryPatch` stable ID、membership 与周期配对唯一且互反，分区面不能成为
  物理 patch。
- 网格度量：uniform adapter 保持 Stage 1 结果；曲线网格验证正 Jacobian、正
  volume、共享面面积向量互反、逐 cell closure、metric identity 和坏网格拒绝。
- 线性代数：CG 覆盖 SPD 制造系统，BiCGStab 覆盖非对称制造系统；两者覆盖零
  RHS、最大迭代、numerical breakdown、非有限输入、错误 control 和 collective
  失败。周期/全 Neumann 压力问题验证 RHS 投影、常数 nullspace 与零均值解；
  每次成功都以独立 FP64 残差复算，而非只相信递推残差。
- 定密度基础：先独立验证 matrix-free Poisson 算子、边界贡献、revision 和
  Jacobi 对角线，再验证固定两次 corrector 的 PISO。Taylor--Green 必须展示
  时间及空间阶，checkerboard 扰动必须受抑制；第二次 corrector 后的最终面质量
  通量必须用于 transport finalization，并在 closure/系数刷新后独立满足连续性、
  压力残差与分区不变性，测试还必须断言没有第三次 corrector。
- 材料变密度：密度波平移覆盖阶数、正性和总质量守恒；制造解与
  variable-density vortex 覆盖质量、三分量动量、机械压力和共享通量耦合；
  不允许焓在此门隐式改变密度。
- 最小热膨胀：分别检查常 `cp` 的 `h -> T`、常 `R` 的 `T -> rho`、非正状态
  拒绝和动态闭域 `p0` 更新；开放域验证配置 `p0` 正、有限且时间不变，机械压力
  出口只约束 `pi`，入口 `h/T/rho/p0` 的唯一推导与冗余值严格交叉校验会在运行
  前拒绝矛盾配置。闭域总质量、`h/T/rho/p0` 一致性、制造解及重试后的时间历史
  必须通过。
- 输运：焓与至少一个通用标量在最终 PISO 面质量通量上通过守恒、MUSCL/MC、
  扩散制造解、正性适用范围和 MPI 分解不变性；finalization 必须从本步起始提交
  状态及其历史层重做，不得提交或再次推进第 7 步 provisional iterate，也不得以
  predictor flux 提交。finalization/closure 刷新造成连续性或压力残差超限时，
  必须在 1/2/4 rank 上一致回滚并缩小 `dt`。
- 边界：periodic、静止无滑移壁面和 symmetry 分别验证速度与守恒通量；规定
  速度入口验证速度、密度或闭合热状态、焓和标量；机械压力出口验证只约束
  `pi` 与纯出流。第 9 步在最终面质量通量发现的任一出口回流必须在 1/2/4 rank
  上 collective 拒绝整步，证据包含 patch、时间步、最小法向质量通量和全局
  face，且提交状态保持不变。
- 曲线有限体积：验证 cell closure、free-stream preservation、非正交扩散与
  流动制造解，并比较不同合法 process grid 的守恒量和收敛结果。
- MPI：所有适用数值案例至少覆盖 1/2/4 rank 的 decomposition invariance；
  Halo 使用 `begin -> interior -> wait -> boundary`，并核对 payload、message 和
  collective 计数；device test double 覆盖 runtime capability 选择、direct/staged
  分支、transfer/event 顺序与无可用路径拒绝；rank-local 数值、边界、I/O 失败
  必须汇总成一致结果且无挂起。
- Checkpoint v2：连续运行与中断续算结果一致；BDF2 的 `n/n-1`、历史 `dt`、
  startup/order、全部 adaptive-controller 状态、最终面质量通量和动态闭域 `p0`
  连续，固定开放域 `p0` 由 fingerprint 校验；Restart 后下一次选定的 `dt` 和
  积分阶次必须与不中断运行相同。缺失 `COMPLETED`、截断、CRC/size/fingerprint
  错误及 trailing bytes 均被拒绝；失败读取不改变字段值或已提交 step/time，
  但进入事务时旧 view 已失效；不同 rank 数或分区必须明确拒绝。
- Stage 1 回归：高风险基础合同变更后运行相关 Stage 1 测试；Stage 2 最终出口
  必须重新通过已验收的完整 85 项 Stage 1 gate，不得更新基线来掩盖回归。

具体网格、初值、步长、误差范数和阈值不在本设计中临时猜定。它们必须在实施
计划中依据公开方程与 reference 计算冻结，并在对应 RED 测试写入前完成复核。

## 16. 八个实施硬门

Stage 2 只能按以下顺序推进；任一门的测试、requirements review、code-quality
review 或主 agent 独立复验未通过，均不得开始后一门的实现：

1. **规格与测量门**：冻结方程符号、残差与守恒定义、边界符号、测试参数、
   source references、resolved-case schema 增量和确定性性能计数。
2. **字段门**：实现 epoch/liveness、capability、checked/kernel view 分层及 stale
   checked view sanitizer 测试，同时保持 Stage 1 `FieldView` 回归。
3. **网格门**：实现 `MeshTopology`、`BoundaryPatch`、uniform adapter 与独立
   `MeshGeometry`，再通过曲线 topology/metrics 的纯几何验收。
4. **执行与线性门**：冻结 `ExecutionContext`、buffer/view/event/transfer 合同，
   实现 CPU reference `VectorOps`、Identity/Jacobi、CG/BiCGStab 和 matrix-free
   Poisson 单元测试。
5. **边界与有限体积门**：实现本阶段五类边界及共享面通量上的对流、黏性和
   标量扩散算子，尚不组装完整变密度时间步。
6. **定密度 PISO 门**：实现 Rhie--Chow、固定两次 PISO corrector、压力
   nullspace 和定密度瞬态验收；不以 SIMPLE 外循环替代。
7. **变密度门**：先验收材料变密度，再验收最小单组分焓--温度--密度与闭域
   `p0`；两者是连续但独立的数值门，不能合并后只报告最终案例。
8. **集成出口门**：固定步长 BDF2 通过后加入 ratio-limited adaptive BDF2，
   随后完成试算/重试、Checkpoint v2、曲线集成、MPI 与性能证据及全量回归。

`cpu_optimized` 不是 Stage 2 科学出口的阻断项；只有第 14 节要求的热点证据充分
时才可另立任务。device 仍只有不可注册为生产计算的合同和 test double，缺少
生产 device 后端不是失败，也不得因此宣称已经具备 GPU 计算能力。

## 17. 实施、复验与审查模式

未来 Stage 2 实施计划必须继续采用 subagent-driven development。每项功能按
以下不可跳过的顺序处理：

```text
RED test
-> fresh implementation worker
-> 主 agent 检查完整 diff 并独立重新运行测试
-> fresh requirements reviewer
-> 修复并重新进行 requirements review
-> fresh code-quality reviewer
-> 修复并重新进行 code-quality review
-> 主 agent 接受该任务并关闭全部相关 worker
```

同一时刻最多保留五个 worker，且只能有一个 implementation worker 活跃；review
worker 也必须是 fresh worker。worker 不得联系用户或向用户索取审批，完成任务
或不再需要后必须立即关闭。不得用实现者的完成陈述替代仓库证据。每次接受前
检查改动文件、提交边界、残留进程、构建依赖、诊断输出和未解释的测试跳过。
不得 publish 或 push。

正式实施计划必须在不扩大本规格范围的前提下，逐项冻结精确容差、误差范数、
网格与映射、rank/process-grid 集合、时间步和步长比限制、失败类别、maximum
retry、`dt` 缩减比例、输出格式及性能 artifact 字段。任何需要改变方程闭合、
增加边界类型或引入依赖的决定都必须先回到设计审批，不能埋入实现任务。

每个门的最终证据必须给出测试命令、退出码、计数或数值结果、accepted commit
以及两类 review 的结论。最终验收必须由主 agent 从已接受的 commit 重新构建，
不能复用 worker 未经核验的工作目录结果。

## 18. 独立性、公开来源与 Stage 2 出口

Stage 2 的方程和算法只从公开论文、教材、MPI 标准和本规格的测试合同推导。
私有黑盒能力矩阵只能列出需要覆盖的可观测科学能力、配置语义和验收证据；它
不能包含或引导读取旧源码实现。不得复制、翻译、机械改写、模仿或兼容旧程序的
源码、控制流、ABI、数组布局、注释、消息、输入格式、Decomp、Restart 或兼容层。

本阶段至少登记并在相应任务中引用以下公开来源：

- C. M. Rhie and W. L. Chow (1983), DOI `10.2514/3.8284`：collocated
  pressure--velocity coupling 与面通量设计的科学背景；
- R. I. Issa (1986), DOI `10.1016/0021-9991(86)90099-9`：PISO
  operator-splitting 的公开算法来源；
- R. J. LeVeque (2002), DOI `10.1017/CBO9780511791253`：守恒有限体积框架；
- B. van Leer (1977), DOI `10.1016/0021-9991(77)90095-X`，及 B. van Leer
  (1979), DOI `10.1016/0021-9991(79)90145-1`：limiter 与 MUSCL 重构；
- M. R. Hestenes and E. Stiefel (1952), *Methods of Conjugate Gradients for
  Solving Linear Systems*：CG；H. A. van der Vorst (1992), DOI
  `10.1137/0913035`：BiCGStab；
- E. Hairer and G. Wanner, *Solving Ordinary Differential Equations II*：
  BDF 方法与步长变化分析；
- MPI Forum, *MPI: A Message-Passing Interface Standard, Version 3.1*
  (2015)：communicator、nonblocking communication 与 collective 合同。

曲线度量、变密度低马赫闭合和时间一致 Rhie--Chow 的具体公开推导来源必须在
实施计划的规格冻结门补齐并经审查后，相关 RED 测试才能开始。来源登记只说明
科学依据，不授权照搬任何公开软件实现。

Stage 2 只有同时满足以下条件才可宣布完成：

- 第 16 节八个硬门全部由主 agent 接受，且第 15 节矩阵无跳过或未解释失败；
- 定密度、材料变密度、最小热膨胀三个门分别留有残差、守恒、阶数与 MPI 证据；
- 每个成功时间步确实固定两次 PISO corrector，以最终修正面质量通量从步初提交
  状态及历史层完成 transport finalization，并且不提交第 7 步 provisional 状态；
- 边界能力与拒绝行为符合第 9 节，出口回流在所有 rank 上事务回滚；
- Checkpoint v2 只按相同 rank 数、process grid 和 owned boxes 恢复，连续性、
  下一次 `dt`/积分阶次精确续算与损坏拒绝通过；Stage 1 Restart v1 与 primitive
  VTK 合同保持不变；
- `cpu_reference` 通过共同 FP64 残差、守恒和收敛合同，device 只保留诚实接口；
- 完整 85 项 Stage 1 gate、source-policy、离线构建、sanitizer 和 linkage 检查
  全部通过，公开构建及运行不依赖 Python；
- 性能 artifact 完整记录确定性计数和可比元数据，不以放宽数值合同换取结果；
- 所有实现与文档提交具有 DCO sign-off，工作树边界清楚，且未 publish/push。

Stage 2 的完成声明不得包含 LES、IBM、化学、物种、多组分热力学、TPDF-TCR、
喷雾、颗粒、生产 GPU、PETSc/HYPRE/AMG、WENO/DG、运动壁面、复杂热化学边界、
rank-changing Restart、全可压缩、声学或激波能力。这些仍按总体路线进入后续
阶段或首版之后；接口预留和 test double 不能被表述为已经实现相应功能。

不得访问、修改、清理、停止、打包或提交私有软件、研究算例、运行数据或研究
进程。私有独立性审计仍在公开仓库外进行，且不得向公共构建或运行引入 Python。

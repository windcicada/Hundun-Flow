# HUNDUN-FLOW Stage 3 静态 LFP-GCIBM 与 WALE LES 设计规格

日期：2026-07-27

状态：本设计的方向与章节已经用户逐项批准；本文件用于书面复核，尚不授权
修改产品代码。用户书面复核本规格后，必须另行编制并批准 Stage 3 实施计划。

本规格承接已经验收的 Stage 2 accepted HEAD
`7b9dbcfedfe93918bdf6cf635b71a563103c412b`。Stage 1 的 schema v1、Restart
v1、primitive VTK、banner 和 passive-scalar 行为保持冻结；Stage 2 的 schema
v2、变密度流动、PISO、Checkpoint v2、诊断与数值合同保持冻结。

## 1. 已批准决定摘要

- Stage 3 采用风险最低的科学推进顺序：
  `几何 -> 层流 IBM -> 独立 WALE -> IBM+WALE -> 集成出口`。
- IBM 是真正的 sharp-interface Ghost-Cell IBM，不采用扩散型体积惩罚，也不把
  求解后的字段覆盖称作 Ghost-Cell。
- 用户论文 Wang et al. (2024) 的局部流型重构、局部缩放/旋转和离散系数变换
  是方法的特色核心；HUNDUN-FLOW 在其上增加物理坐标二次约束、仿射 Ghost
  消元和 PISO 一致压力条件，形成 `LocalFlowPattern Ghost-Cell IBM`，简称
  `LFP-GCIBM`。
- 论文中的工程验证不能单独证明二阶。Stage 3 必须以制造解、系统网格加密和
  圆柱/球体壁面力收敛独立证明速度、压力、壁面穿透和积分受力均达到观测阶
  `>= 1.8`。
- IBM 首版只支持一个静态 STL、一个静态几何 part、不可渗透和无滑移壁面。
  焓与通用标量固定使用零法向扩散通量。
- 不使用 wall function 作为层流二阶壁面力闭合；壁面梯度与受力由同一个二次
  重构在真实三角面求积点计算。
- WALE 是 Stage 3 唯一生产 LES 模型。常系数 Smagorinsky、dynamic
  Smagorinsky 和其他亚格子模型均延后。
- IBM 与 WALE 产品路径必须兼容 Stage 2 的 `constant`、`material` 和
  `ideal_gas` 三种密度模型。
- 定密度案例承担完整且昂贵的二阶、壁面力和 LES 基准；材料变密度与
  ideal-gas 路径仍必须通过守恒、正性、失败回滚、Restart 和 MPI 集成门。
- 保留同一个 `hundun` 可执行程序。新增 schema v3 和独立 Stage 3 组合路径，
  不改变 schema v1/v2 的解析、canonical JSON 或运行行为。
- Stage 3 不实现化学、物种、多组分热力学、TPDF-TCR、喷雾、颗粒、生产 GPU、
  moving IBM、multi-part IBM、cut-cell、AMR、WENO/DG 或复杂热壁。

## 2. 论文评估与方法定位

本设计依据用户提供的论文：

> Y. Wang, F. Wang, J. Zhou, and J. Jin, “An improved immersed boundary
> method with local flow pattern reconstruction and its validation,”
> *Physics of Fluids* 36, 045145 (2024), DOI `10.1063/5.0195598`.

所复核 PDF 的 SHA-256 为：

```text
f6d012a6bd62b2a755b4ecf2476d4783c03dbb5a51a36a60715aff53539d03e0
```

论文方法通过几何扫描获得真实壁面位置和法向，在边界邻近单元建立缩放并旋转
到壁面对齐的虚拟局部网格，变换七点离散方程系数，再通过边界源项补充真实
壁面信息。该方法有三个适合 HUNDUN-FLOW 的特点：

1. 它直接作用于 cell-centred 离散方程，与当前有限体积布局一致；
2. 静态几何允许把扫描、局部坐标和系数映射预计算；
3. 它能够在不建立 body-fitted 网格的情况下保留真实壁面位置和法向。

论文采用二阶中心空间格式和 Crank--Nicolson 时间格式，并通过圆柱、GTMC 和
真实燃烧室展示工程有效性，但没有给出制造解、三层系统加密、速度/压力误差
范数、壁面穿透收敛或压力/黏性力收敛。因此，本规格不能把论文结果表述为 IBM
二阶证明。

原方法中的线性局部变化假设、局部缩放前后单元平均量不变假设、wall-function
剪切闭合，以及未单独冻结的压力修正边界，不足以直接承担本项目的严格二阶
合同。HUNDUN-FLOW 保留论文的局部流型和系数变换思想，但作以下明确增强：

- 线性场假设升级为物理坐标中的受约束二次重构；
- 边界源项表达升级为可审计的仿射 Ghost 约束；
- Ghost 约束在动量、压力和输运算子内部消元；
- 压力修正条件由零修正壁面法向质量通量离散推导；
- 层流壁面梯度和受力使用真实表面求积，不使用 wall function；
- 背景流体--流体面的共享通量仍保持成对互反，不被局部变换破坏。

这一路线保留论文方法的辨识度，同时使“二阶”成为必须由测试证明、也可能被
测试否决的科学合同。

## 3. 范围与非目标

Stage 3 交付以下能力：

- 一个静态、闭合、可定向的三角化 STL 表面；
- 分区无关的表面 fingerprint、空间查询、fluid/solid 分类和壁面截点；
- 每个流体--固体算子 link 的二次 `GhostStencilPlan`；
- 局部流型坐标、缩放/旋转和边界对齐系数变换；
- 速度 Dirichlet、焓/标量 homogeneous Neumann 和 PISO 压力 Ghost 条件；
- 真实表面上的压力力、黏性力和总力积分；
- WALE 亚格子黏度及焓/通用标量的梯度扩散闭合；
- constant/material/ideal-gas 三种密度模型下的 IBM 与 WALE 组合；
- Stage 3 独立 Checkpoint、结构化诊断、性能计数和 1/2/4-rank 验收；
- 一个不可注册为生产能力的静态 wall-kinematics/motion test double，用于证明
  后续接口能够表达非零壁速，但 Stage 3 产品配置只能得到零壁速。

Stage 3 明确不包含：

- immersed cut-cell volume/area fractions、小 cut-cell 稳定化或 conservative
  remap；
- 多 STL、多 part、相交/嵌套表面和 thin-shell 双侧流动；
- 移动、旋转、变形壁面，新覆盖/新暴露单元或 moving-boundary GCL；
- 壁面粗糙度、壁函数、等温壁、规定热通量、共轭传热或反应表面；
- 多入口/多出口扩展、热化学回流、物种、化学和反应热；
- dynamic/constant Smagorinsky、DES、RANS、TPDF-TCR；
- 生产 device kernel、GPU-aware MPI、混合精度生产求解器；
- 一般非结构网格、AMR、WENO、DG、PETSc/HYPRE/AMG；
- rank-changing Checkpoint 或将 Stage 3 文件伪装成 Restart v1/Checkpoint v2。

Stage 4 继续拥有可配置热壁、反应气相多入口/多出口、温度/焓/组分边界和
热化学回流。Stage 3 的零法向焓/标量通量不能被宣传为完整热壁能力。

## 4. 依赖架构与组件边界

依赖方向固定为：

```text
application / time driver
  -> IBM/LES flow composition / wall force / Stage 3 diagnostics
  -> immersed operator adapters / WALE / finite-volume / PISO
  -> GhostStencilPlan / local-flow-pattern transform / surface queries
  -> mesh / fields / execution / linear / Halo / MPI
```

新增职责分离为：

- `ImmersedSurface`：不可变三角几何、稳定 triangle ID、方向、面积与 fingerprint；
- `SurfaceQuery`：不可变空间索引、最近点、射线和线段相交；
- `ImmersedClassification`：全局 cell 分类、interface link 和 active-cell 映射；
- `ActiveBoundaryLayout`：从六个背景 patch 过滤仍与 active fluid 相邻的
  boundary faces，并决定入口、出口、周期配对和压力 reference 是否有效；
- `GhostStencilPlan`：donor、局部基、二次权重、秩、条件估计和依赖 halo；
- `WallQuadraturePlan`：真实 triangle quadrature points 的 fluid-side 二次
  reconstruction、稳定 triangle ownership 和只读受力求积权重；
- `LocalFlowPatternTransform`：论文来源的局部缩放/旋转与系数空间变换；
- `GhostConstraintTemplate`：从静态几何权重和当前壁面数据生成仿射约束；
- `ImmersedOperatorAdapter`：在粗粒度算子边界消元，不在逐 cell 路径虚调用；
- `WallForceIntegrator`：从同一二次重构计算真实表面 traction；
- `WaleModel`：只读速度/密度/几何并产生瞬时亚格子输运系数；
- `Stage3FlowDriver`：组合已有 Stage 2 事务、IBM、WALE、Checkpoint 和诊断。

这些组件不得共同继承一个侵入式调试基类。诊断通过只读 provider/adapter
接入。几何、stencil 和局部变换不拥有 `FlowState`；WALE 不拥有 MPI
communicator；应用层不包含逐 cell 数值 kernel。

Stage 2 的 `MeshTopology` 和 `MeshGeometry` 仍描述背景结构化网格。IBM 不把
三角表面伪装为 `BoundaryPatch`，也不把 fluid/solid marker 写入背景 topology。
二者通过稳定 global cell/face ID 和物理坐标关联。

## 5. Schema v3 与兼容入口

现有类型和函数保持不变：

```cpp
using ResolvedCase = std::variant<CaseConfig, FlowCaseConfig>;

ResolvedCase load_resolved_case(const std::filesystem::path&);
std::string to_resolved_json(const ResolvedCase&);
ResolvedCase broadcast_resolved_case(
    MPI_Comm, int root, const ResolvedCase* root_case);
```

Stage 3 使用加法式版本接口，避免改变 `ResolvedCase` 的 C++ 类型布局和 v1/v2
行为：

```cpp
struct ImmersedFlowCaseConfig;
using ResolvedCaseV3 =
    std::variant<CaseConfig, FlowCaseConfig, ImmersedFlowCaseConfig>;

ResolvedCaseV3 load_resolved_case_v3(const std::filesystem::path&);
std::string to_resolved_json_v3(const ResolvedCaseV3&);
ResolvedCaseV3 broadcast_resolved_case_v3(
    MPI_Comm, int root, const ResolvedCaseV3* root_case);
```

同一个 `hundun` dispatch shell 改用 v3-aware 入口。该入口对 schema v1/v2
必须委托冻结的原解析路径，并产生与 Stage 2 accepted HEAD 字节相同的
canonical resolved JSON。原 `load_resolved_case()` 仍拒绝 schema v3。

Schema v3 继续使用：

```text
simulation.type = "variable_density_flow"
density_model = constant | material | ideal_gas
mesh.mapping = uniform_box | analytic_warped_box
```

并增加：

```text
immersed_boundary.model = none | local_flow_pattern_ghost_cell
immersed_boundary.geometry.format = stl
immersed_boundary.geometry.file
immersed_boundary.geometry.length_scale_to_m
immersed_boundary.geometry.fluid_side = inside | outside
immersed_boundary.wall.velocity_m_per_s = [0,0,0]
immersed_boundary.wall.enthalpy = zero_normal_diffusive_flux
immersed_boundary.wall.scalars = zero_normal_diffusive_flux

les.model = none | wale
les.wale.coefficient
les.turbulent_prandtl
les.turbulent_schmidt
```

四种组合的对象形状固定为：

- `ibm=none` 时，`immersed_boundary` 对象只能包含 `model`；`geometry` 和
  `wall` 必须缺失；
- `ibm=local_flow_pattern_ghost_cell` 时，`geometry` 和 `wall` 全部上述成员
  都是 required，不能省略，也不能增加 motion、wall-function 或 thermal 键；
- `les=none` 时，`les` 对象只能包含 `model`；`wale`、
  `turbulent_prandtl` 和 `turbulent_schmidt` 必须缺失；
- `les=wale` 时，`wale.coefficient`、`turbulent_prandtl` 和
  `turbulent_schmidt` 全部 required；本阶段不使用隐式默认值；
- `none/none` 组合拒绝；`none/wale`、`LFP-GCIBM/none` 和
  `LFP-GCIBM/wale` 是三个合法组合；
- canonical JSON 对合法组合只输出该组合 required 的键，并使用上文固定顺序。

合法组合的下游对象存在性同样固定：

- `none/wale` 使用覆盖全部背景 cell 的 identity `ActiveCellLayout`，并直接
  使用 Stage 2 boundary 语义；不构造 surface、classification、interface-link、
  Ghost stencil、LFP、wall quadrature 或 wall-force provider，也不产生这些
  模块的计数；
- `LFP-GCIBM/none` 构造全部 IBM 对象，但不构造 WALE transient field、
  coefficient identity、provider 或计数；
- diagnostics 和 Checkpoint v3 对未启用模块写入下文定义的 canonical
  absence，而不是伪造空 surface、空 stencil 或零 fingerprint。

其余约束如下：

- schema v3 中 IBM 与 LES 至少启用一个，避免建立没有新增能力的 v2 别名；
- `local_flow_pattern_ghost_cell` 恰好引用一个 case-root 内 STL 文件；
- `length_scale_to_m` 必须正且有限；
- `fluid_side` 必须显式给出，不能靠三角方向猜测用户意图；
- wall velocity 三个分量只能在数学上等于零，canonical JSON 统一输出正零；
  非零、时变或 motion 键在 rank 0 预检阶段拒绝；
- h/标量壁面模式只能是固定的零法向扩散通量；
- WALE 系数和 turbulent Prandtl/Schmidt 必须正且有限；实施计划冻结其允许
  范围和验收基准值；
- schema v2 不接受这些新增键，v1/v2 未知键行为保持不变；
- 文件路径仍不得为绝对路径、包含 `..` 或逃逸 case root；
- 每个合法/非法组合、冗余子对象、缺失 required key 和 forbidden key 都必须
  有 exact JSON Pointer 测试及 typed MPI broadcast/canonical round-trip。

配置不暴露 polynomial degree、一阶 fallback、任意 regularization、wall
function 或“忽略坏 stencil”开关。二次重构是固定科学合同，不能被运行配置
降级。

## 6. 静态三角表面合同

### 6.1 输入与规范化

项目自有 C++17 reader 读取批准计划冻结的 STL 子集，不引入外部几何库或
Python。输入坐标先乘 `length_scale_to_m`，随后要求：

- 所有坐标、法向输入和计算量有限；
- 三角形面积严格为正且超过尺度相关退化阈值；
- 顶点焊接使用与几何尺度相关的确定性派生容差；它不是用户 schema 或
  canonical resolved JSON 键，而是 surface validation 的派生量，必须进入
  diagnostics、surface fingerprint 和 Checkpoint v3 compatibility；
- 每条无向边恰好连接两个三角形；
- 相邻三角方向一致；
- 表面闭合、可定向、具有非零封闭体积；
- 只有一个 connected component；
- 非相邻三角形不得自交；
- 表面不得与背景域的开放入口或压力出口相交；
- 表面与背景域边界之间必须满足实施计划冻结的最小可分辨距离。

STL 中存储的法向仅作为一致性诊断；产品法向由经过验证的顶点顺序重新计算。
闭合表面先建立一致的几何外法向，再根据 `fluid_side` 确定本规格统一使用的
“固体指向流体”法向：fluid 在表面外部时取几何外法向，fluid 在表面内部时取
其反向。稳定 triangle ID 来自规范化后仍保持输入顺序的三角记录。实现不得按
本地 rank 或空间索引遍历顺序重新编号。

### 6.2 SurfaceQuery

静态空间索引只在初始化构建。其 public query 至少支持：

```text
closest point and triangle ID
segment/surface intersections
ray parity classification
surface bounding box
bounded candidate-triangle query
```

空间索引内部节点顺序必须由稳定 triangle ID 和确定性分割规则决定。同一 STL、
映射和请求在不同 rank/process grid 上必须产生相同分类、截点、法向和
fingerprint。

射线落在边、顶点或与三角面近共面时不能任意计数。分类采用固定方向集合和
半开相交规则；多个确定方向不能得到一致结果时，几何预检 collective 拒绝，
不得按 rank、本地浮点偶然性或随机扰动选择结果。

### 6.3 分类与 interface link

每个背景 cell centre 分类为 `fluid` 或 `solid`。`ghost` 是 solid cell 针对一个
相邻 fluid operator row 的代数角色，不是第三种物理区域。每条由背景
owner/neighbour 连接产生的 fluid--solid link 必须：

- 在两 cell centre 的物理线段上与表面恰好相交一次；
- 记录稳定 `ImmersedLinkId`、fluid global cell ID、solid global cell ID、
  triangle ID、截点、固体指向流体的单位法向和无量纲截距；
- 对 MPI 分区两侧得到相同记录；
- 不进入六个 body-fitted `BoundaryPatch`；
- 在物理壁面质量通量预算中恰好出现一次。

一个 solid cell 可服务多个 interface link。Ghost 约束按 link 和 fluid row
存储，不能把多个不同截点强行压成一个全局 ghost 值。

fluid cells 构成 `ActiveCellLayout`。线性向量只包含 active fluid unknowns；
solid rows 不使用 identity 方程混入残差，也不进入收敛范数或守恒量。

cell-centre 分类本身不能证明表面已经被网格解析。初始化还必须通过
`SurfaceCoverage` 预检：

- 全局 fluid cell、solid cell 和 interface link 数均非零；
- 每个有面积 triangle 的全部受力 quadrature points 都关联至少一个合法
  interface link 和 active fluid row；
- 每个 quadrature point 在冻结的尺度相关半径内同时具有 fluid-side 与
  solid-side resolution witness，并且 witness 连线与该表面局部相交；
- surface 到最近 interface link 截点的最大距离不超过实施计划冻结的
  local-mesh-size bound；
- 所有 surface quadrature contributions 只关联一次，且任何参与算子的
  interface link 都具有已验证的 wall intercept；
- 小于可解析网格尺度、完全落入单个 cell、因相位位置没有 solid centre，或
  局部薄特征未产生 link 的几何必须 collective 拒绝，而不是被当作空 IBM。

测试必须平移同一 STL 穿过多个 mesh phase，并包含“闭合小物体完全落入一个
cell”的明确拒绝 oracle。

### 6.4 ActiveBoundaryLayout

`ActiveBoundaryLayout` 在 classification 后由 Stage 2 `BoundaryRegistry` 和
背景 patch 构造，只保留 owner 为 active fluid 的 background boundary faces。
Stage 2 registry 本身不修改；Stage 3 driver、PISO、回流检查和边界通量只查询
该 additive view。

- velocity inlet 和 pressure outlet 各自若配置，必须至少有一个 active face；
  零 active face 的开放 patch 在初始化 collective 拒绝；
- pressure outlet 只有在具有 active face 时才提供机械压力 reference；
- validation 必须先于压力算子构造；配置但零 active-face 的 pressure outlet
  没有可运行分支，必须先拒绝；
- 通过预检且根本未配置 pressure outlet 的合法闭域，压力算子使用 constant
  nullspace 和零均值规范化；
- periodic active faces 必须按原 topology 一一互反配对；一对 patch 同时零
  active faces 可以保持 inactive，但只有一侧 active 或 active pairing 不完整
  必须拒绝；
- no-slip/symmetry patch 的 inactive faces 不参与 flux、诊断、回流或计数；
- `fluid_side=inside` 的闭合 cavity 通常使六个背景 patch 全部 inactive，此时
  必须得到闭域 pressure-nullspace 语义，不能从背景 patch 猜测为开放域。

`ActiveBoundaryLayout` 的 ordered global face IDs、patch counts、pairing 和
pressure-reference decision 进入 diagnostics 与 Checkpoint v3 fingerprint。

## 7. GhostStencilPlan 与二次约束

### 7.1 物理坐标多项式

每个 interface link 在真实壁面截点 `x_b` 建立正交局部坐标
`(n,t1,t2)`。`n` 为固体指向流体的表面法向；`t1/t2` 由确定性轴选择产生，
不依赖三角遍历顺序。局部坐标按 stencil 尺度归一化。

重构使用完整三维二次基：

```text
1, n, t1, t2, n^2, n*t1, n*t2, t1^2, t1*t2, t2^2
```

donor 只能是物理 fluid cell。候选按物理距离、方向覆盖和 global cell ID
确定性选择。分区边界 donor 通过计划声明的 structured Halo 获取，不能因 donor
属于另一 rank 而改变 stencil。

Stage 2 的 `CellAverage` 数值不能被当作 cell-centre point sample。对 donor
cell `V_j` 和局部二次基 `b_m(x)`，设计矩阵行固定为：

```text
M_jm = (1/volume(V_j)) * integral_over_Vj(b_m(x) dV)
```

uniform cell 的这些 moments 必须在 FP64 roundoff 内精确；曲线背景使用与
`MeshGeometry` 表示同一物理 cell 的确定性体积分/多面体 moment，不得以逻辑
cell centre 代替。计划必须证明完整二次 cell-average polynomial 可重现 wall
value、normal/tangential gradient 和 ghost-centre value，覆盖 uniform 与
`analytic_warped_box`。

小型稠密系统使用项目自有、确定性的列主元 QR 或等价 rank-revealing FP64
方法预计算权重。必须记录：

- donor global IDs 和归一化坐标；
- 设计矩阵 rank；
- 条件估计；
- QR/pivot fingerprint；
- 所需最大 structured halo reach；
- Dirichlet、Neumann 和梯度求值权重；
- 计划 revision 与 geometry/classification fingerprint。

设计矩阵秩不足、条件估计超过冻结上限、方向覆盖不足、donor 越过壁面或 halo
reach 不可满足时，case 在开始时间推进前 collective 拒绝。禁止：

- 降为一次或常数重构；
- 用最近 cell 值替换；
- 未声明的 Tikhonov regularization；
- 忽略坏 link；
- 只在失败 rank 删除该 link；
- 让 product path 使用 test-only mutation seam。

### 7.2 仿射 Ghost 约束

对 interface link `P--G`，Ghost sample 的坐标固定为 solid neighbour 的实际
background cell centre `x_G`；它既不是 wall intercept，也不是另造的对称
image point。每个变量在该位置的 link-local 值写成：

```text
q_G = sum_j(w_j q_Dj) + c_D q_wall + c_N g_wall
```

其中 donor 权重来自静态 plan，`q_wall` 是 Dirichlet 数据，`g_wall` 是沿 `n`
的 Neumann 数据。静态无滑移速度使用 `q_wall=0`；焓和通用标量使用
`g_wall=0`。

材料密度不具有由不可渗透条件推出的 homogeneous Neumann 边界。需要
`rho_wall` 或近壁密度梯度时，只使用不带 wall derivative constraint 的
fluid-side 单边二次 cell-average 外推；唯一物理密度边界合同是质量壁面通量
为零。必须用合法的非零法向密度梯度证明该外推没有偷偷施加零梯度。

ideal-gas 路径由零法向焓通量、常 `cp/R` 和空间均匀 `p0` 得到一致的
`T/rho` 零法向延拓，不能独立覆盖由 closure 得到的密度。

仿射约束由算子在组装或 matrix-free apply 中消元。产品代码不得先求解含任意
solid 值的方程，再覆盖 ghost 或 fluid 值。

## 8. LocalFlowPatternTransform

`LocalFlowPatternTransform` 保存论文方法的特色结构：

1. 以 fluid cell、wall intercept 和局部法向构造虚拟 boundary-aligned cell；
2. 表达从背景局部轴到 `(n,t1,t2)` 的旋转；
3. 表达使虚拟控制面落在真实壁面的法向尺度变换；
4. 提供论文定义的完整局部 coefficient-row 变换，以及产品使用的守恒
   masked/marginal wall replacement；
5. 将 wall closure 写回 link-local 的算子贡献。

Stage 3 的变换不是任意 post-hoc source overwrite。它必须满足：

- 常数场保持；
- 线性场在平面壁、均匀正交网格上精确；
- 旋转前后量纲和符号一致；
- standalone `T_LFP(full coefficient row)` 在论文的平面/均匀/线性假设极限
  下逐系数退化到论文公式；
- 曲线背景网格使用物理几何向量，不能把逻辑坐标距离误当物理距离；
- fluid--fluid 共享 face flux 的两侧值保持互反；
- 所有由变换产生的非共享贡献都明确归属真实 immersed-wall flux；
- 与 `GhostStencilPlan` 的二次约束组合后，不能重新引入一阶 wall closure。

计划应提供两个相互独立的测试 oracle：

- standalone coefficient-map oracle，验证
  `T_LFP(full coefficient row)` 与论文 Eqs. 13--16；
- direct quadratic evaluator，在产品 masked/marginal replacement 的局部点
  用 cell-average 二次多项式直接求值。

产品 replacement 必须与第二个 evaluator 在 FP64 容差内一致。两个 oracle
只用于测试，不进入时间推进热点。产品为保持 active--active shared-face
守恒而采用 masked/marginal HUNDUN 扩展；除下文明确列出并由测试证明的特殊
极限外，不宣称产品最终 row 等于论文对完整 coefficient array 做一次
`T_LFP` 后的 row。

## 9. 算子数据流与守恒

### 9.1 唯一 immersed residual 代数

对每个 active row `P` 和每个离散方程，令 `L(P)` 为该 row 的全部 immersed
links。一个 row 可以有 1--6 个此类 link，不能假设只有一个 solid neighbour。
在任何消元前冻结一次不可变的完整背景行快照：

```text
R_bg_full(P)
  = R_shared(P) + R_immersed_bg(P)

R_immersed_bg(P)
  = sum_{g in G(P)} R_g,bg(P)
```

`R_shared(P)` 包含全部 active--active shared-face 项且不依赖任何 inactive
slot。`R_immersed_bg(P)` 是 ordinary background form 中所有指向 link-local
Ghost symbol 的 direct、orthogonal、non-orthogonal、deferred 和 gradient
依赖的并集；它不是从 canonical-zero inactive field slot 求值出来的数值。

`G(P)` 是在 plan 构造时确定的 replacement groups。每个 group `g` 保存：

- 非空稳定 link 集 `links(g) subseteq L(P)`；
- 一个 `mask_g`，选择 `R_immersed_bg(P)` 中属于该 group 的代数项
  occurrence；
- 所有关联几何、二次 reconstruction weights 和 joint-evaluator fingerprint。

全部 group masks 必须对 `R_immersed_bg(P)` 的代数项 occurrence 构成完整、
互斥的 partition。一个可分项属于恰好一个 singleton group；一个不可分的
cross-link reconstruction/deferred 项属于恰好一个 multi-link group，不能
复制给多个 groups。每个 `L(P)` 中的 link 必须被至少一个 group 引用，且其
唯一 link-local Ghost symbol 只有一份 affine constraint 定义。

产品 public operation 是一次 row-local 调用：

```text
W_P = wall_replacements(
  snapshot(R_bg_full(P)),
  {(group_g, mask_g, links(g), geometries(g)) | g in G(P)})
```

其数学定义为：

```text
W_P = sum_{g in G(P)} W_g

if links(g) = {l}:
  W_g = Eval(
    T_LFP,l(Coeff(R_bg_full(P)))
    - T_LFP,l(Coeff(R_bg_full(P) - R_g,bg(P))),
    q_local,l)

if |links(g)| > 1:
  W_g = JointEval_g(
    snapshot(R_bg_full(P)),
    mask_g,
    {(q_local,l, geometry_l) | l in links(g)})
```

这里的差分是 masked/marginal replacement 的数学定义，不是产品对 row 执行
两次 transform；实现利用 coefficient map 对 row coefficients 的线性性，在
一次 row-local 调用中直接形成同一结果。`JointEval_g` 是由 group 的二次
cell-average reconstruction 和各 link LFP geometry 在 plan 构造时生成的
pure linear evaluator；它同时读取所有关联 link symbols，只读同一个 immutable
snapshot，不允许按 link 顺序逐次修改 row。它对输入 link permutation
数学不变；实现按排序后的稳定 `ImmersedLinkId` tuple 固定 FP64 累加和
fingerprint。multi-link group 是 HUNDUN masked/marginal 扩展，不宣称等于
论文 whole-row transform。

最终 row 只有一个定义：

```text
Q_G(P) = {
  q_G,l -> affine_ghost_l(q_D, q_wall, g_wall)
  for every l in L(P)
}

R_final(P)
  = R_bg_full(P)
  - R_immersed_bg(P)
  + simultaneous_substitute(W_P, Q_G(P))
```

`R_shared(P)` 因而保持原来的 owner/neighbour shared-flux 贡献；LFP 引起的
其他 local donor coefficients 全部属于唯一 `W_P` wall reaction，而不是
第二份 shared-face flux。完整 `W_P` 先构造，再对 `Q_G(P)` 做恰好一次
row-wide simultaneous substitution；每个 link-local Ghost symbol 只有一个
权威 affine definition，即使它出现在多个 group terms 中也不能被顺序覆盖。
改变输入 group/link 遍历顺序不得改变代数、FP64 canonical order 或
fingerprint。

等价的产品残差链固定为：

```text
snapshot immutable R_bg_full at active fluid cell P
-> partition R_immersed_bg into complete, disjoint group masks G(P)
-> call wall_replacements(R_bg_full, group-set) exactly once for row P
-> remove R_immersed_bg exactly once from the actual row
-> apply Q_G(P) once as a row-wide simultaneous substitution in complete W_P
-> add W_P exactly once to row P
```

不存在第二个 wall forcing、第二次 Ghost substitution 或 solve 后覆盖。LFP
transform 是产品 `W_P` 的实际组成部分，不是只读诊断；每个仿射 Ghost
constraint 只在对应 link-local symbol 出现的位置代入。
`WallQuadraturePlan` 不向 momentum row 写任何 source，因此不能与该 residual
双计数。

standalone `T_LFP(full coefficient row)` 在 uniform、正交、平面壁和线性场
极限必须逐系数复现 Wang et al. 论文 Eqs. 13--16 的缩放/三轴旋转变换。
产品 `R_final` 则明确是为保留 active--active shared-face 互反性而定义的
masked/marginal HUNDUN 扩展；一般情况下不声称它等于
`T_LFP(R_bg_full)`。只有单 link、变换对 `R_shared` 为 identity、且测试逐系数
证明 `T_LFP(R_shared)=R_shared` 的特殊极限，才能声明产品最终 row 与论文完整
row 等价。

Eq. 17 的 solid-neighbour coefficient 消去结构和 Eq. 18 的
zero-normal-scalar closure 分别作为论文 coefficient-map 等价证据；HUNDUN
的 affine Ghost 消元和二次 cell-average reconstruction 是单独的增强证据。
Eq. 17 中依赖 wall function 的剪切项明确不采用，改由同一二次 polynomial
的 fluid-side velocity gradient 给出。测试必须分开记录论文等价部分、
masked/marginal 守恒扩展和这项有意差异，不能把三者合成一个模糊的“完全
复现论文”声明。

曲线/非正交网格上，ordinary `P--G` 的 orthogonal 与 explicit non-orthogonal
两部分必须一起移除；变换后的 wall value/gradient 由同一个物理坐标二次
polynomial 一次产生。不得保留旧 non-orthogonal correction 后再加入新的
wall-gradient correction。

独立测试必须对每个方程分别证明：

- 在普通 row 中只变动 `P--G` closure 时，fluid--fluid shared flux 完全不变；
- transformed-row 加一次 Ghost substitution 与 direct cell-average quadratic
  evaluator 的完整 row residual 一致；
- 删除 LFP transform 或重复 Ghost/wall term 的 mutation oracle 会失败；
- 同一 active row 含 2/3 个 solid neighbours 时，完整 group-mask partition、
  group/link 顺序置换不变性、漏 group、重复 group、cross-link term 复制和
  sequential-substitution mutation oracle；
- singleton `W_g` 与 marginal evaluator 一致，multi-link `JointEval_g` 与
  direct quadratic full-row evaluator 一致；
- constant、linear 和 quadratic polynomial reproduction；
- pressure、viscous、h/scalar 各自没有 double counting；
- uniform 与 warped、interior 与 partition-boundary multi-link row 使用同一
  代数。

### 9.2 Interface-band reconstruction 与执行顺序

`InterfaceBand` 包含所有拥有 immersed link 的 active rows，以及其
reconstruction/gradient stencil 触及的 active rows 和 active--active faces。
以下产品路径必须统一使用同一个 IBM-aware reconstruction provider：

- Rhie--Chow、pressure operator 和 pressure-correction gradient；
- momentum convection、viscous gradient 和 non-orthogonal/deferred term；
- density、h 和 scalar 的 MUSCL、diffusion 与 non-orthogonal term；
- WALE velocity gradient；
- wall value、wall gradient 和 surface traction。

该 provider 只能读取 active fluid `CellAverage` donors、预计算二次 moments/
weights 和 link-local affine Ghost symbols。它不得读取 owned 或 ghosted
inactive cell field slot，即使该 slot 满足 canonical positive-zero invariant。
canonical zero 是 state hygiene 和 Checkpoint 合同，不是数值边界条件。

因此 interface-band 的 `R_bg_full` 在 snapshot 前就必须由 IBM-aware provider
形成；ordinary `P--G` 依赖以 link-local symbolic Ghost 出现，而不能先调用
Stage 2 的 full-topology gradient 读取 solid zero、再试图只替换直接 face 项。
若未来加入的新 reconstruction/deferred 路径无法证明其依赖被统一 provider
覆盖，则该路径属于 layout/contract failure，不能静默进入产品。

active--active shared face 的 gradient、reconstruction 和 flux 仍只计算一次并
以互反符号加入 owner/neighbour；wall-derived correction 只进入第 9.1 节的
`W_P`。测试必须包括：

- uniform/warped 二次 cell-average polynomial 的 pressure、velocity、rho、
  h 和 scalar value/gradient reproduction；
- instrumented read-only test adapter 对 inactive-slot access 立即失败；
- 隔离 kernel oracle 改变非产品 fixture 的 inactive backing payload，而
  active row 和 active--active shared-face flux bitwise 不变；
- interior、partition-boundary 和 multi-link interface rows 的相同依赖图。

instrumented adapter 只存在于测试，不暴露产品 mutation seam；端到端产品案例
仍严格拒绝任何非 canonical-zero inactive state。

静态 plan 构造完成后，一次 IBM-aware matrix-free 操作固定为：

```text
Halo begin
-> 计算不依赖远端 donor 的 active interior rows
-> Halo wait
-> 更新需要远端 donor 的 link-local ghost constraints
-> 计算 partition-boundary rows
-> 计算 immersed-interface-band rows
```

不同 row 集合必须互斥且覆盖全部 active rows。一次 apply 中不查询 STL、
不遍历 BVH、不动态选择 donor、不分配 `FieldStorage`，也不改变 plan revision。

所有 fluid--fluid face 的对流、黏性和扩散通量继续以共享 face record 成对加入
owner/neighbour。fluid--solid 背景 face 的普通 neighbour flux 按第 9.1 节完整
移除，不能与 transformed closure 重复计算。质量壁面通量数学上为零；h/标量
扩散壁面通量数学上为零。动量的唯一 immersed contribution 是第 9.1 节最终
row residual。

Stage 3 的守恒量定义在 active fluid cells 使用的背景 cell volumes 上。该定义
是 Ghost-Cell 离散合同，不宣称实现 cut-cell 几何体积积分。制造解和 MPI
比较必须使用同一 active-cell measure；真实表面力单独在三角表面求积。

全局动量预算使用唯一 immersed row residual 的负和，称为
`operator_reaction_force`，并按 row 中的 pressure 与 viscous 项分别累计。
真实三角求积得到的 `surface_traction_force` 是独立物理精度量。二者不得彼此
替代：前者闭合离散动量，后者承担用户批准的 pressure/viscous/total force
二阶合同；两种表达的 pressure、viscous 和 total 差值也必须分别随网格至少
二阶收敛。

### 9.3 Full storage 与 inactive slots

Stage 3 通过 additive adapter 使用现有完整背景 `FlowState`/`FieldStorage`，
不改变 Stage 2 对象的布局、轮换或提交行为。active vector/残差使用
`ActiveCellLayout`，但每个 committed/history/trial layer 仍具有完整 cell/face
slots。

所有 owned inactive cell 的全部分量固定为 canonical FP64 positive zero；所有
solid--solid、fluid--solid 和 inactive body-boundary face 的 `FaceMassFlux` 及
其他 transient face slots 同样固定为 positive zero。正性、闭合、残差和物理
诊断只遍历 active IDs，不能把这些规范零误判为物理 `rho/T`。

初始化、trial 创建、retry rollback、commit、Halo workspace 准备和 Checkpoint
restore 都必须重新验证该 invariant。任何非零或负零 inactive owned slot 是
layout/state failure；产品不能用未初始化 solid 值作为 donor。这样 v3 可以
继续序列化完整背景布局并做 bitwise 比较，而 v1/v2 路径完全不变。

## 10. PISO 与压力 Ghost 条件

Stage 2 每个成功试算恰好两次 pressure corrector 的合同保持不变。LFP-GCIBM
不能增加第三次 corrector，也不能在 PISO 后覆盖速度来伪造无穿透。

对静止 wall link，压力修正边界由修正后的法向质量通量为零推导：

```text
mdot_wall^(k+1)
  = mdot_wall^* - D_wall * normal_gradient(pi')_wall
  = 0
```

`mdot_wall^*` 必须来自完整动量 predictor 和时间一致 Rhie--Chow 路径；
`D_wall` 是完整的正 mass-flux correction coefficient，包含同一 link 的
`rho_wall`、有效壁面 measure、真实截距以及由该步实际动量方程对角系数得到的
速度修正系数。该 effective measure 来自第 9.1 节唯一 transformed row，
不来自只读 `WallQuadraturePlan`。不能遗漏面积量纲，也不能用固定几何常数或
仅当前压力项替代。

material-density 的 fluid-side 外推 `rho_wall` 必须在形成 `D_wall` 前逐 link
验证为正且有限；constant/ideal-gas 的对应 wall density 也必须正且有限。
最终 `D_wall` 必须正且有限。失败按最低 rank collective 汇总为可恢复的
`non_positive_state` 或 `non_finite_state`，不得 clipping、改用零梯度、
替换 stencil 或把系数设为 epsilon 后继续。

几何部分和 Neumann gradient 权重由 `GhostStencilPlan` 预计算；每次 corrector
只更新依赖 `mdot_wall^*` 和 `D_wall` 的仿射常数。压力算子 revision
在这些动态系数、密度、动量对角、geometry 或 active layout 改变时按冻结规则
递增，Jacobi diagonal 不得跨 revision 误用。

应用压力修正后：

- fluid--solid background `FaceMassFlux` slot 保持 canonical positive zero，
  derived wall-link mass flux 按上述边界方程得到数学零，而不是对负值截断；
- fluid velocity 的 wall-normal 与 tangential no-slip 约束来自同一 Ghost
  reconstruction；
- pressure gauge/nullspace 仍按 Stage 2 周期/Neumann 规则处理；
- immersed closed wall 不提供额外机械压力 reference；
- 只有 `ActiveBoundaryLayout` 中至少含一个 active face 的 body-fitted pressure
  outlet 才提供 pressure reference 并禁用 nullspace；完全 inactive 的 outlet
  不能改变压力约束；
- checkerboard、最终 continuity 和两次 corrector 计数继续使用 Stage 2 阈值。

若最终独立重算的 wall penetration、continuity 或 pressure residual 超过合同，
该 trial 是可恢复数值失败，必须 collective 回滚并减小 `dt`，不能增加
corrector 或覆盖最终 flux。

## 11. 三种密度路径和输运边界

### 11.1 Constant

定密度路径是完整二阶 IBM 和 LES reference。它承担制造解、圆柱、球体、
wall-force、channel 和 IBM+WALE 的全部高成本科学矩阵。

### 11.2 Material density

材料密度仍由 Stage 2 质量方程保守推进：

```text
partial(rho)/partial(t) + div(mdot_f) = 0
```

所有 immersed-wall mass flux 为零；`rho*h` 和 `rho*phi` 使用同一最终 PISO
面质量通量。必须验证：

- active-fluid 总质量守恒；
- `rho` 正且有限；
- 合法非零法向 `rho` 梯度由 fluid-side 单边二次外推保持，而不是被改成
  homogeneous Neumann；
- 正 donor 但二次外推产生非正/非有限 `rho_wall` 的 trial collective 拒绝，
  且没有 clipping 或 stencil fallback；
- h/标量零法向扩散通量；
- final-flux provenance；
- 失败 trial 的 rho/h/scalar/history bitwise rollback；
- 1/2/4-rank decomposition invariance。

### 11.3 Ideal gas

Stage 2 最小闭合保持：

```text
h = cp*T
rho = p0/(R*T)
```

静态绝热 IBM 不改变该关系。闭域 `p0` 仍按 active-fluid volume 与目标总质量
更新；开放域 `p0` 仍为配置常量。IBM 分类、active volume 和 geometry
fingerprint 必须进入闭域质量目标及 Checkpoint compatibility。

非正 `T/rho`、非有限 closure 或 IBM boundary reconstruction 产生无效状态时，
沿用 Stage 2 collective failure 和 retry。Stage 3 不增加焓反馈外循环，也不
实现物种依赖 `cp/R`。

## 12. 壁面 traction 与受力

`WallQuadraturePlan` 使用与 `GhostStencilPlan` 相同的物理坐标二次基、donor
合法性、rank-revealing QR、条件上限和无 fallback 合同，但在真实 triangle
quadrature points 独立预计算只读求值/梯度权重。它没有写入 momentum row 的
API；任一有面积的 triangle 无合法 fluid-side quadrature stencil 时，初始化
collective 拒绝，不能跳过其受力。
在每个批准的三角求积点重构：

```text
pi
u
grad(u)
mu_eff
```

固体指向流体的单位法向记为 `n_s`。流体应力为：

```text
sigma = -pi*I + tau
```

流体作用于固体的力固定为：

```text
F_body = integral_surface(sigma * n_s) dA
F_pressure = integral_surface(-pi * n_s) dA
F_viscous = integral_surface(tau * n_s) dA
```

闭合低马赫 `p0` 为空间常数，其在完整闭合表面的解析净力为零；受力报告使用
机械压力 `pi`，并单独诊断三角面积向量 closure，不能把数值不闭合的 `p0`
贡献混入 drag。

压力、黏性和总力必须分别做 MPI reduction，并记录最低失败 rank。三者的
surface quadrature 与 triangle ownership 必须分区无关，每个 triangle 恰好由
一个 rank 积分。报告必须同时给出第 9.2 节的 `operator_reaction_force` 和本节
`surface_traction_force`；前者用于离散动量闭合，后者用于 pressure/viscous/
total force coefficient。二者的 pressure、viscous 和 total 差值是三个独立
consistency residual，必须有限并按第 18 节分别达到二阶，不能通过把
quadrature force 回填 solver 强制变成零。

时间步验收中的 `ForceAttemptReport` 是 attempt-local 派生证据，不属于
committed/history state，也不持久化“last accepted force report”。失败时丢弃
该 report；rollback 后若需要比较，只能从未改变的 committed state 重新收集，
并要求确定性相同。对 committed state 的只读 force/diagnostic 收集不得修改
flow state、solver cache 或业务性能计数器。

层流二阶 gate 禁止 wall function。未来高 Reynolds 壁面模型必须另立批准设计，
不能作为 Stage 3 精度失败的 fallback。

## 13. WALE LES

### 13.1 数学闭合

速度梯度记为：

```text
g_ij = partial(u_i)/partial(x_j)
S_ij = 0.5*(g_ij + g_ji)
G_ij = g_ik*g_kj
Sd_ij = 0.5*(G_ij + G_ji) - delta_ij*G_kk/3
```

WALE 运动黏度采用 Nicoud--Ducros 形式：

```text
nu_t = (Cw*Delta)^2
       * (Sd_ij*Sd_ij)^(3/2)
       / ((S_ij*S_ij)^(5/2) + (Sd_ij*Sd_ij)^(5/4))
```

当分子与分母的数学不变量均为零时，`nu_t` 精确为正零。实现不得用会在静止
流中产生非零黏度的固定 epsilon。任何尺度保护必须保持齐次性，并在计划中
冻结。

```text
Delta = cbrt(active background cell volume)
mu_sgs = rho*nu_t
mu_eff = mu + mu_sgs
```

亚格子各向同性应力吸收到机械压力 `pi`，不增加独立 SGS kinetic-energy
方程。动量方程只把 `mu_sgs` 加入与 Stage 2 一致的 deviatoric strain-rate
应力，不能把各向同性 SGS 分量重复加入机械压力和黏性项。焓和通用标量采用
梯度扩散：

```text
Gamma_h,sgs = mu_sgs / Pr_t
Gamma_phi,sgs = mu_sgs / Sc_t
```

WALE 只读本节冻结的 lagged velocity、指定 trial density 和 geometry，写入
声明的 transient `nu_t/mu_sgs` 字段。它不能直接修改速度、质量通量、压力或
committed state。速度梯度必须使用 IBM-aware 二次 reconstruction，不能跨过
wall 读取任意 solid 值。

### 13.2 变密度解释

材料变密度和 ideal-gas 路径使用 Favre-filtered 动量/标量语义：质量与所有
保守输运仍共享 Stage 2 最终 `FaceMassFlux`，`mu_sgs=rho*nu_t` 使用同一 trial
density。WALE 不提供新的密度闭合，不改变 `p0`，也不把体积平均速度替代为
另一份独立速度状态。

### 13.3 验证边界

WALE 单元证据至少覆盖：

- 零速度梯度得到 bitwise 正零；
- independent tensor evaluator 数值一致；
- 坐标正交旋转下标量 `nu_t` 不变；
- 网格尺度和速度尺度的量纲缩放；
- canonical near-wall expansion 的 `y^3` 渐近行为；
- 非有限梯度、密度或配置的 collective 分类；
- `mu_eff >= mu >= 0` 且有限；
- constant/material/ideal-gas 的同一路径一致性。

独立 LES 使用 periodic Taylor--Green/decaying turbulence 和 body-fitted
channel。IBM+LES 使用完全位于背景域内、由一个 closed connected STL 定义的
有限长圆柱或球体尾流。论文的 `Re=3900` 跨展向圆柱结果作为方法背景和未来
工程比较资料，不作为 Stage 3 强制出口；本阶段不为复刻该算例而引入穿越周期
边界的开放/周期商表面语义。任何以后与论文数据的正式比较都必须先冻结几何、
展向边界、统计窗口和可比性，不能把有限长物体的结果直接当作跨展向圆柱结果。

Stage 3 只声明“提供通过批准基准的 WALE baseline”，不声明 DNS 精度、通用壁模
能力或所有 Reynolds 数上的实验吻合。

### 13.4 WALE 时层与 lagging

Stage 3 不增加 LES 非线性外迭代。每次 attempt 的 WALE 系数只计算一次：

1. Stage 2 transport prediction 和本次 density closure 完成后，得到
   `rho_attempt`；
2. startup 使用 `u_lag=u^n`；具有 BDF2 history 时使用变量步长线性外推：

   ```text
   r = dt_attempt / dt_previous
   u_lag = u^n + r*(u^n - u^(n-1))
   ```

   这里使用与本次 attempt 相同且已经通过 `[0.5,2.0]` 限制的 `r`；
3. 在 momentum predictor 前，由 `rho_attempt` 和该 lagged velocity 计算
   `nu_t/mu_sgs/mu_eff`；
4. 同一 attempt 的 momentum predictor、两次 PISO corrector、transport
   finalization、最终 residual 和 wall traction 全部使用这份冻结系数；
5. 第二次 corrector 后不按 final velocity 重新计算 WALE，也不追加第三次
   momentum/pressure solve。

retry 使用相同 committed/history，但按新的 `dt` time-stencil 和新的
`rho_attempt` 重新计算一份 attempt-local WALE 状态。`WaleModel` 不拥有或提交
monotonic model revision；每份系数只携带由 WALE 配置、attempt step/dt/order、
committed/history state fingerprints 和 `rho_attempt` fingerprint 确定性派生的
`WaleCoefficientIdentity`。失败和成功 attempt 结束后都不把该 identity 变成
新的 committed state。

momentum/pressure operator 仍按 Stage 2 coefficient replacement 管理自己的
revision；它不能把 attempt-local WALE identity 冒充为 operator revision。
`nu_t` 和 `WaleCoefficientIdentity` 都不作为权威 Checkpoint 字段，Restart 后
由保存的 history、controller 和相同 trial prediction bitwise 重算。最终
residual 与
`surface_traction_force` 必须使用本 attempt 实际冻结的 `mu_eff`，不能用
final velocity 临时生成另一份黏度。

## 14. 时间事务、失败分类和 collective 一致性

Stage 2 的 committed/history/trial、adaptive BDF2、两次 PISO corrector、
transport finalization 和 retry 控制流保持不变。IBM/WALE 插入点固定为：

1. 初始化前完成 surface、classification、active layout、stencil 与 transform；
2. transport prediction/closure 后按第 13.4 节唯一 lagging 规则生成
   attempt-local WALE coefficients；
3. momentum predictor 使用 IBM-aware `mu_eff` 和 Ghost constraints；
4. 每次 pressure corrector 更新动态 pressure Ghost constants；
5. 第二次 corrector 后使用最终共享 flux 完成 transport finalization；
6. final closure 后独立重算 wall penetration、residual、conservation 和 force；
7. 全部合同通过后才提交 Stage 2 state；static plan 本身没有 trial mutation。

不可重试失败包括：

- STL 语法、路径、单位、拓扑、方向、自交或域相交错误；
- 分类歧义、零/多 wall intercept、surface coverage 或两侧 resolution witness
  不满足；
- stencil 秩、条件、donor 或 halo reach 错误；
- schema、active layout、fingerprint、capability 或 Checkpoint 完整性错误；
- MPI 操作错误和已经提交的输入损坏。

可恢复失败包括：

- trial 中非有限/非正 flow state；
- 派生 `rho_wall` 或 `D_wall` 非有限/非正；
- WALE 从有效配置和当前 trial 得到非有限系数；
- 线性不收敛或 numerical breakdown；
- 最终 continuity、momentum、h/scalar、wall penetration 或 pressure residual
  超限；
- Stage 2 已定义的最终出口回流。

所有失败必须给出统一分类、最低失败 rank 和稳定 invariant ID。任何 rank
不得单独改变 active layout、删除 stencil、提交状态或进入下一次 retry。

## 15. Checkpoint v3

Checkpoint v3 是独立格式；Restart v1 与 Checkpoint v2 的 API、字节和目录保持
冻结。v3 manifest 对 IBM 和 LES 分别保存 stable presence tag；每个 tag 只能
是 `absent` 或该模块冻结的 algorithm/version ID。`absent` 是 canonical
absence：对应 section count 和 byte count 均为零，且不得写入伪造的空对象
fingerprint。v3 在 v2 状态合同上增加：

- schema v3 canonical fingerprint；
- IBM presence tag；仅当 IBM enabled 时保存 STL 规范化 byte/content
  fingerprint、length scale、`fluid_side`、surface topology/orientation、
  classification、interface-link、active-layout、triangle ownership、
  SurfaceCoverage、ActiveBoundaryLayout、GhostStencilPlan、
  WallQuadraturePlan、LocalFlowPatternTransform 和 wall-force
  integrator/diagnostic fingerprints；
- IBM absent 时只保存 canonical absence tag；运行时使用 all-active identity
  layout 和原 Stage 2 boundary view，不能构造或比较 surface/stencil/LFP
  fingerprint；
- LES presence tag；仅当 WALE enabled 时保存 WALE 配置和 transient-field
  schema；
- LES absent 时只保存 canonical absence tag，不能构造或比较
  `WaleCoefficientIdentity` 或 WALE field fingerprint；
- Stage 2 全部 committed/history、dt/controller、final flux 和动态 `p0`。

BVH 节点、原始地址、平台相关缓存、派生 wall-force report 和 `nu_t` 瞬时值
不作为权威持久化状态。
恢复时先读取并验证文件，再从相同 STL、背景 mesh 和配置确定性重建
surface/classification/stencil/transform；只有 fingerprint 全部相等才可替换
flow state。重建失败不能用 checkpoint 中缓存权重绕过当前验证。

v3 延续：

- little-endian、binary64、CRC-64/ECMA-182；
- manifest、每-rank 文件、最终 `COMPLETED`；
- 相同 rank count、process grid、owned boxes 和 active layout；
- 完整背景 fields/face slots 按 Stage 2 field order 持久化；读取时所有 inactive
  owned cell/face doubles 必须为 canonical positive zero，否则语义拒绝；
- 进入读取事务即使旧 checked view 失效；
- 失败读取不改变字段值、committed step/time/history/controller；
- 成功续算的 fields、history、next `dt/order`、final flux、dynamic `p0` 和
  下一步按第 13.4 节重算的 IBM/WALE 结果与不中断运行 bitwise 相同。

Stage 3 不支持 geometry-changing、rank-changing 或 repartitioning restore。

## 16. 统一诊断与性能证据

`DiagnosticModuleKind` 以加法方式在现有 `performance=17` 之后追加稳定值：

```text
immersed_surface   = 18
ghost_stencil      = 19
local_flow_pattern = 20
wall_force         = 21
les                = 22
```

现有 0--17 的数值、`kDiagnosticRecordSchemaV1=1` 和旧 canonical records
保持不变。`ActiveBoundaryLayout` 继续使用既有 `boundary` kind，并以稳定
module/instance ID 区分。standalone header、underlying-value snapshot 和
Stage 2 canonical-record regression 必须证明追加枚举没有重编号旧值。

新增 provider/adapter 接入 `hundun::diagnostics`，至少覆盖：

- surface：triangle/vertex/component 数、bbox、面积、闭合、orientation 和
  fingerprint；
- query/classification：fluid/solid 数、interface links、歧义数和 active layout；
- stencil/quadrature：donor count min/max、rank、condition estimate、halo
  reach、revision、triangle coverage 和 fingerprint；
- LFP transform：scale/rotation invariant、coefficient norm 和 limiting-case
  状态；
- PISO IBM：predictor/final wall flux、wall penetration、dynamic pressure
  constraint 和 lowest failing rank；
- wall force：pressure/viscous/total force、moment、area closure 和单位；
- WALE：`nu_t/mu_sgs/mu_eff` min/max/norm、zero/active cell count、config 和
  attempt-local `WaleCoefficientIdentity`；
- driver/checkpoint：phase、step、retry、state fingerprint 和 v3 compatibility。

provider inventory 必须服从 presence tag：IBM absent 时不注册 18--21 的
provider、instance 或 counters；LES absent 时不注册 kind 22 的 provider、
instance 或 counters。缺失模块不是一条 status 为零的伪记录。driver 和
checkpoint provider 只报告稳定 presence/absence tag，并且同一合法组合在
1/2/4 ranks 上得到确定性相同的 module inventory。

诊断继续满足 Stage 2 已冻结行为：

- local scope 不隐式 collective；
- collective scope 必须显式且 rank 一致；
- 不修改 state、revision、generation、cache、allocation identity、controller
  或业务计数器；
- bounded sample 使用稳定 global cell/link/triangle ID；
- 非有限 FP64 以 value-status 和 bit representation 表达；
- 不输出原始地址、临时路径或完整场；
- 相同状态和请求产生确定 canonical records；
- 关闭诊断时，热点路径不增加隐式 allocation、collective 或全场复制。

exact counters 至少包括：

- surface triangle queries 和 candidate intersections；
- classification rays/segments；
- stencil donor、QR、rejected-plan 和 halo-reach counts；
- 每步 ghost-constraint evaluation 和 LFP transform counts；
- immersed interface rows、pressure wall constraints 和 wall quadrature counts；
- WALE gradient/cell counts；
- IBM 额外 Halo payload/message、linear matvec/reduction 和 I/O logical bytes。

初始化预计算与每步成本分开记录。Portable CI 对确定性 counters 硬判；wall
clock、RSS、带宽和吞吐只在兼容元数据下保存为 baseline，不设置跨机器阈值。

## 17. 后端中立与未来 GPU 兼容

Stage 3 只实现生产 `cpu_reference`。为了不封死后续 GPU：

- surface query 只在初始化使用；时间推进 kernel 不访问 STL parser 或 BVH；
- `GhostStencilPlan` 转换为连续、项目自有 buffer 和 trivially-copyable kernel
  view；
- donor IDs、offsets、weights 和 row ranges 使用显式 layout，不暴露
  CUDA/HIP/SYCL 类型；
- IBM 和 WALE 通过粗粒度 `ExecutionContext` 调用，不在逐 cell 路径虚调用；
- matrix-free apply 保留 `begin -> interior -> wait -> boundary/interface`；
- direct/staged Halo 仍由 runtime capability 决定；
- CPU 与未来 device 使用相同 FP64 residual、conservation、wall-force、
  convergence、rollback 和 decomposition contracts；
- test double 只能验证 capability/lifetime/rejection，不能注册为生产 GPU。

Stage 3 不实现 GPU-aware MPI、device STL scanning、混合精度或 vendor backend，
也不得据此宣称 GPU 计算已经可用。

## 18. 二阶精度与科学验收

### 18.1 误差语义

对 refinement levels `h, h/2, h/4`，其中 `h` 是实施计划冻结并由实际
`MeshGeometry` 重算的 characteristic spacing，连续两段观测阶定义为：

```text
p_1 = log(E_h/E_h2)/log(2)
p_2 = log(E_h2/E_h4)/log(2)
```

所有被列为二阶 hard gate 的量都要求 `p_1 >= 1.8` 且 `p_2 >= 1.8`。不能只报告
拟合直线平均斜率来掩盖其中一段失败。

除第 18.2 节 wall-penetration exact-enforcement 分支外，每个 `E` 必须有限、
严格为正、随两次 refinement 严格下降，并高于实施计划按 reference scale
冻结的 FP64 roundoff floor。reference solution、force coefficient
normalization 和每个误差的非零尺度必须在运行前确定，不能从候选结果反推。
不得给误差加 epsilon、截断为正数或在 `log()` 后忽略 NaN/Inf。若一个原本应
非零的 manufactured quantity 落入 roundoff floor，该案例无判定力，必须更换
制造解幅值而不是宣称无限阶。

数值离散误差使用批准容差；rollback、failed trial 和 Checkpoint 连续性使用
FP64 bitwise 比较；MPI decomposition 使用计划冻结的 max-field threshold。

### 18.2 必须通过的二阶量

在含斜平面局部面的闭合棱柱、圆柱和球体几何上分别验证：

- velocity volume-weighted L2 和 L-infinity；
- pressure gauge-normalized volume-weighted L2 和 L-infinity；
- 一层固定物理厚度近壁带内的 velocity/pressure L2；
- surface L2 和 L-infinity wall-normal penetration；
- integrated pressure force coefficient；
- integrated viscous force coefficient；
- integrated total force coefficient。
- `operator_reaction_force - surface_traction_force` 的 pressure、viscous 和
  total 三个归一化 consistency residual。

压力、黏性、总力和两种 force 表达的 consistency residual 必须分别达到
`>=1.8`，不能只要求它们相加后的总力。制造解必须使这些参考分量具有非零、
可归一化尺度，不能用接近零的 lift 掩盖误差。

若 wall-normal penetration 因仿射 Dirichlet 约束在三个网格上都不超过实施
计划冻结的 FP64 roundoff-scaled enforcement bound，则按“exact enforcement”
通过，这是强于二阶的结果；不得对 `0/0` 伪造收敛阶。只要任一层超过该 bound，
就必须按上述公式得到两段 `>=1.8` 的有限观测阶。

### 18.3 几何与制造解

二阶案例必须使用实际 STL reader 和产品分类/stencil/operator 路径。测试夹具由
纯 C++ 生成随背景网格同步加密的闭合斜棱柱、圆柱和球体 STL，使
faceting/chord error 至少为二阶并记录 geometry error；不能用 test-only
analytic geometry 绕过 STL path。

制造解体源只存在于测试，不进入 public schema。解析解满足静态 no-slip 和批准
的 h/scalar Neumann 条件，并提供解析 pressure、velocity gradient 和表面
traction。

至少覆盖：

- uniform background；
- `analytic_warped_box` 的批准受限映射；
- interface 穿过不同相位的 cell links，避免只测网格对齐壁面；
- 1/2/4 ranks 和至少两种合法 process grid；
- surface/mesh 相对位置平移后的稳定阶数；
- 二次 `CellAverage` polynomial 的 wall value/gradient/ghost-centre
  reproduction，不能用 point-sample oracle；
- `fluid_side=inside` 时全部 background patches inactive、pressure nullspace
  正确；配置零-active-face outlet 时 collective 拒绝；
- 小闭合物体落入单 cell、无 solid-centre 或 surface coverage 超限时
  collective 拒绝；
- stencil conditioning 接近但未超过阈值的合法案例；
- 秩不足、病态、分类歧义和多截点的 collective 拒绝案例。

### 18.4 物理基准

除制造解外至少包括：

- 一个 `fluid_side=inside` 的闭合 cavity 制造解，用于验证全部背景 patch
  inactive 时的 active-boundary/nullspace 语义；
- 低 Reynolds 圆柱流，检查 drag、lift、separation/recirculation 趋势；
- 球体低 Reynolds 流，检查压力/黏性 drag 分解；
- immersed Taylor--Green 或等价封闭瞬态；
- WALE body-fitted channel；
- 完全位于背景域内的有限长圆柱或球体 IBM+WALE 尾流。

工程案例允许使用公开接受区间，不承担形式二阶证明。完整分辨率、运行时间、
统计窗口和阈值由实施计划在 RED 测试前冻结。

### 18.5 密度、事务与 MPI

constant 路径执行完整二阶和 LES 矩阵。material/ideal-gas 至少分别覆盖：

- IBM wall mass flux 和 active-fluid mass conservation；
- `rho/T` 正性及 closure identity；
- h/标量零法向扩散；
- shared final-flux provenance；
- 所有 inactive cell/face slots 在 trial、rollback、commit 和 v3 restore 后
  都是 canonical positive zero；
- 一次 rank-local injected trial failure 后全体相同失败类别和最低失败 rank；
- fields、nested transported fields、history、controller 和 p0 的 bitwise
  rollback；失败后从原 committed state 重算的 attempt-local
  `operator_reaction_force/surface_traction_force` 必须与失败前确定性相同，
  被拒绝 trial 的派生 report 不属于 committed state；
- 1/2/4-rank decomposition invariance；
- Checkpoint v3 连续/续算 bitwise 一致。

## 19. 实施硬门

未来实施计划必须保持科学闭环，并按以下硬门推进：

1. **规格、schema 与几何门**：数值合同、schema v3、STL、surface validation、
   deterministic query、SurfaceCoverage、active cell/boundary layout；
2. **Ghost 计划门**：cell-average moments、二次 QR plan、donor/Halo、LFP
   limiting-case、仿射约束和明确拒绝，不组装完整流动；
3. **层流 IBM 门**：唯一 residual 链、算子消元、PISO pressure Ghost、共享
   通量、operator/surface force 和 constant-density 全部二阶矩阵；
4. **独立 WALE 门**：tensor formula、near-wall scaling、冻结 lagging、
   body-fitted LES、density-path 单元合同；
5. **IBM+WALE 与变密度门**：合法 closed finite-body wake、
   material/ideal-gas 守恒/正性/事务/MPI；
6. **集成出口门**：canonical inactive slots、Checkpoint v3、诊断、性能、
   曲线背景、工程案例和 Stage 1/2 全量回归。

不得为缩短 diff 机械拆开必须联合验证的几何、Ghost、压力、共享通量、守恒或
壁面力闭环。每个 task 可包含多个 acceptance cluster，但最终 verdict 只能对
完整 task 给出 accepted/rejected。

## 20. TDD、审查和 coordinator 协议

Stage 3 继续采用 subagent-driven development：

```text
RED test
-> fresh implementation worker
-> 主 agent 检查完整 task diff 并独立复验
-> fresh requirements reviewer
-> 必要 repair 和重新 requirements review
-> fresh code-quality reviewer
-> 必要 repair 和重新 code-quality review
-> 主 agent exact-HEAD 验收并关闭全部 worker
```

每个 task 实现前冻结 evidence matrix，建立双向 traceability：

```text
requirement
-> implementation
-> positive test
-> failure/rollback test
-> MPI/numerical acceptance
```

reviewer 必须审查从上一个 accepted task 到 candidate HEAD 的完整 task diff，
使用 `codegraphf` 查符号、调用方、影响范围和相关测试，使用 `rg` 搜索精确文本
与同类断言。连续两轮出现同类缺口时，主 agent 先做只读 closure sweep，再派
fresh repair worker 一次关闭完整问题类别。

Stage 3 继承 coordinator acceptance acceleration protocol v2：

- 中间循环运行 RED、最小 GREEN、finding closure 和 task-focused 矩阵；
- 每个最终 task candidate 由主 agent 在 exact HEAD 至少运行一次完整 Debug；
- Release/ASan/UBSan 使用 evidence matrix 冻结的聚焦矩阵；
- hard gate 与 Stage 3 最终出口执行计划规定的全配置和全量回归；
- 测试证据绑定 HEAD、preset、toolchain、MPI、rank、参数、产物和日志 SHA；
- reviewer 的直接证据、主 agent 完整 Debug 和最终出口命令不得跨角色省略。

worker/reviewer 交接包必须包含 accepted base、candidate HEAD、brief SHA-256、
required reading、允许文件、evidence matrix、不变量、findings、closure 条件、
测试矩阵、commit subject、exact DCO 和禁止范围。worker 不得联系用户、扩大
范围、访问私有源码目录、发布或 push。

统一测试 helper 必须区分 bitwise、数学精确和容差比较；nested collections
比较 outer size、inner size 和全部元素。helper 必须有 mutation-sensitive
oracle，证明 exact copy 通过、普通字段变化和嵌套字段变化均失败。

## 21. 独立性与公开来源

Stage 3 只从公开论文、教材、本规格和测试合同独立推导。用户论文是公开数学
来源，不授权复制任何旧程序实现。BOFFIN 固定基线仍只用于仓库外私有独立性
审计；COAST/COAST-2 仍只能作为黑盒能力和科学趋势参考。

至少登记：

- Wang et al. (2024), DOI `10.1063/5.0195598`：局部流型、局部坐标与离散系数
  重构的科学来源；
- Möller and Trumbore (1997), DOI `10.1080/10867651.1997.10487468`：射线--
  三角形相交；
- Tseng and Ferziger (2003), DOI `10.1016/j.jcp.2003.07.024`：非交错网格上的
  second-order Ghost-Cell 边界条件；
- Pan and Shen (2009), DOI `10.1002/fld.1942`：速度和压力 Ghost-Cell
  收敛验证；
- Seo and Mittal (2011), DOI `10.1016/j.jcp.2011.06.003`：sharp-interface IBM
  的质量守恒与压力振荡问题；
- Nicoud and Ducros (1999), DOI `10.1023/A:1009995426001`：WALE 模型；
- Stage 2 已登记的 Rhie--Chow、PISO、BDF2、有限体积和 MPI 来源继续适用。

实施计划必须在写对应 RED 测试前核对公式、符号、适用范围和可公开引用的
验证数据。不得读取或模仿旧源码的 IBM、LES、数组布局、控制流、消息、输入、
Decomp、Restart 或兼容层。

## 22. Stage 3 完成与能力声明

只有同时满足以下条件才能宣布 Stage 3 完成：

- 第 19 节六个硬门全部 accepted；
- LFP-GCIBM 在批准的闭合斜棱柱、圆柱和球体矩阵中，速度、压力及
  pressure/viscous/total force 两段观测阶全部 `>=1.8`，wall penetration
  达到 `>=1.8` 或满足第 18.2 节更强的 exact-enforcement 条件；
- 没有一阶 fallback、忽略 stencil、post-solve overwrite 或 wall-function
  精度替代；
- 每个 interface row 严格使用第 9.1 节唯一 LFP-transform/Ghost-substitution
  residual，WallQuadraturePlan 不写 solver row，double-count mutation oracle
  通过；
- SurfaceCoverage 能拒绝未被 cell-centre/link 解析的几何，
  ActiveBoundaryLayout 能正确过滤 inactive patch 和 pressure reference；
- pressure Ghost 由最终零 wall-normal corrected mass flux 推导，并保持恰好
  两次 PISO corrector；
- constant/material/ideal-gas 三条路径的守恒、正性、shared flux、rollback、
  Restart 和 MPI 合同均通过；
- WALE 的 tensor、near-wall、LES 和 IBM+LES 基准通过；
- Checkpoint v3 连续性、损坏拒绝和相同分区限制通过，Restart v1/Checkpoint
  v2 保持冻结；
- Stage 3 diagnostics provider coverage、determinism、no-state-mutation 和
  exact counters 通过；
- Debug、Release、ASan、UBSan、tests-off、header、policy、provenance、
  linkage、Stage 1 acceptance 和 Stage 2 acceptance 按计划全量通过；
- 所有 task commits 具有 DCO，全部 worker 关闭，无残留测试进程；
- 未访问或干预私有软件、研究算例、运行数据和研究进程；
- 未发布、未 push，公开构建与运行不依赖 Python。

Stage 3 能力声明只能说：

```text
单静态 STL 上通过二阶合同的 LFP-GCIBM，
以及通过批准基准的 CPU-reference WALE LES。
```

不得声称已经支持 cut-cell、移动/多部件 IBM、热壁、壁函数、化学、TPDF-TCR、
喷雾、生产 GPU、一般复杂燃烧室工作流或完全替代 COAST。

# HUNDUN-FLOW v0.4 Cartesian 低速可压缩性能架构设计

初始日期：2026-08-12

最终修订：2026-08-13

状态：最终批准

目标版本：HUNDUN-FLOW v0.4

唯一发布节点：Re=3900 三维圆柱绕流数值与性能门通过

## 1. 产品目标和不可变约束

v0.4 是 HUNDUN-FLOW 的新主线。目标不是为圆柱算例编写特化程序，而是建立可用于
外流、内流、射流和分离流等单相湍流问题的通用 Cartesian CFD 底层。圆柱和后台阶
只承担验证职责；更换几何、边界或物理参数必须只修改 case 输入，不得修改产品源码。

本版采用以下不可变约束：

1. 只支持单相、亚音速、无激波流动，但保留由流动、压力功和传热导致的密度与温度
   变化；不以“不可压缩”删去局部绝对压力对 EOS 的作用。
2. 只有一条低速可压缩产品路径，不保留常密度快速路径，不维护第二套方程或数组。
3. 只支持均匀 Cartesian 和全局张量积拉伸 Cartesian；不支持 body-fitted、多块曲线
   网格、AMR 或非匹配 refinement patches。
4. v0.4 只实现瞬态 PISO，每个时间步恰好两次 pressure corrector。稳态问题使用物理
   或伪瞬态推进；SIMPLE、PIMPLE 和 strong PIMPLE 不进入 v0.4。
5. Linux CPU 是生产后端；性能设计优先服务 MPI、NUMA、固定线程团队和 SIMD。
   GPU、移动 EB 和动态负载迁移不进入本版。
6. v0.3 单独封存为备份；v0.4 可破坏旧内部 API，但必须保留已验证科学合同，并对
   尚未真正证明的能力重新建立证据。
7. 性能优化不得减少 PISO corrector、放宽守恒/残差/IBM 几何约束、增加算例特定
   阻尼、滤波或硬编码分支。
8. 性能候选冻结前不运行长统计；小网格不得用作性能结论。

## 2. 参考边界与独立实现

第一工作单元固定核对下列项目的公开数学、数据布局和生命周期思想：

| 项目 | 固定身份或位置 | 参考范围 |
| --- | --- | --- |
| OpenFOAM-dev | `b9da51ab0673423aa2af6a45a72a3fbec9c66f9f` | PISO 时序、`rAU/rAtU`、`HbyA`、`phiHbyA`、最终 `U/phi` 权威关系 |
| AMReX | `59d066aab774bc388cc6ed944f7beaf645607ed3` | box-local 布局、halo/fill-patch 元数据、EB factory 和资源生命周期 |
| IncFlo | `7307d8725c2a538f09cafbeacbfeb63e0fb11d22` | projection、EB、operator/projector 重建边界 |
| AMReX-Hydro | `e49df248aabd2cc11865eb5be734a2f5f2f65ee5` | multigrid、operator 和 workspace 复用 |
| COAST | `/home/wyf/code_dev/Coast_software` | 功能范围、输入输出、SIMPLE/ICCG 热循环、数组和预生成数据思想 |
| STL 扫描参考 | 用户附件 `imb_mesh_y.cpp` | 用户授权直接复用或小改的扫描算法 |

OpenFOAM 为 GPL，只能参考公开数学和生命周期，不复制或翻译源码。AMReX/IncFlo
采用宽松许可证，但 HUNDUN 仍独立实现。除用户明确授权的 STL 扫描 C/C++ 方法外，
不得复制旧 COAST Fortran。用户明确要求该 STL 方法不记录原文件、作者或权利声明。

参考调研形成以下设计结论：

- PISO corrector 之间可复用的是经 revision 证明未变的动量系数及 `rAU/rAtU`，不是
  所有中间量；trial `U/phi` 改变后必须重建或重新认证 `HbyA/phiHbyA`。
- 最终 face mass flux 只能由最终压力方程 flux 通路发布，不能从最终速度另算一份。
- 静态 EB 几何、halo 元数据、线性 symbolic plan、MG hierarchy 和最大容量 workspace
  必须跨时间步持久化，并具有彼此独立的失效规则。
- AMR 的双层 fill-patch、reflux 和 average-down 不进入 v0.4。

## 3. 版本、源码与运行目录

### 3.1 版本隔离

仓库采用并列版本目录，根构建默认 v0.4，并允许显式构建 v0.3：

```text
hundun-flow/
  CMakeLists.txt
  VERSION
  versions/
    v0.3/
    v0.4/
      CMakeLists.txt
      include/hundun/
      src/
      tests/
      docs/
```

- `versions/v0.3` 只从圆柱工作树 commit
  `4ae4832ad00b5e4d1129ee978e2e49cbb33bb7ef` 的 tracked tree 导出，不读取或纳入
  该工作树当前 dirty/untracked 内容。
- Stage 4 已接受科学能力的独立参考身份为
  `6407cd7c591ce088db7f1dd7e296d77acd18da1c`，tree
  `2791a1cee7ac8114f1696670d30c8951212d6024`。
- 现有圆柱工作树和 `/home/wyf/code_dev/.benchmarks/cylinder-re3900*` 证据只读保留，
  不 reset、clean、覆盖或提交其中现有修改。

### 3.2 v0.4 扁平源码

`versions/v0.4/src` 是统一扁平产品源码目录，只采用八类稳定前缀：

| 前缀 | 所有权 |
| --- | --- |
| `app_` | CLI、case 编译、driver |
| `core_` | 类型、状态、arena、revision、事务、执行图 |
| `mesh_` | Cartesian 网格、STL 扫描、EB/IBM 静态计划 |
| `parallel_` | MPI、NUMA、线程、halo |
| `bc_` | 物理边界与冷路径编译 |
| `physics_` | EOS、输运、湍流、组分和源项 |
| `solver_` | 离散、方程、PISO、线性系统、MG |
| `io_` | Restart、Visit、screen、monitor、证据 |

公开头位于 `include/hundun`；非公开声明尽量与实现留在 `src`。源码树保留测试，但
不保存产品 `cases/` 或 `examples/`。安装后的运行工具可在用户指定运行目录创建
`examples`；案例与输出不污染源码。

## 4. Case 输入和冷路径编译

### 4.1 扁平 case 目录

一个 case 的输入直接位于其根目录：

```text
case-name/
  case.json
  mesh-focus.d
  inlet-profile.d
  cylinder.stl
```

`case.json` 是唯一主配置，负责启用功能、引用文件并规定单位。`*.d` 只保存大型数组、
坐标、剖面和表格，不得暗中覆盖 JSON。STL 和所有 `.d` 文件直接放在 case 根目录，
不设置 `geometry/`、`data/` 等必需子目录。

rank 0 完成读取、schema 校验、单位换算、文件哈希、引用闭合和物理兼容性检查，随后
广播紧凑类型化模型。其他 rank 不重复解析文本。无效输入必须在分配大字段或创建 MPI
persistent request 前集体失败。

### 4.2 初始化编译链

```text
CaseSpec
  -> ValidatedModel
  -> capability registration
  -> FieldSchema + CartesianGeometryPlan
  -> BoundaryPlan + EBTopology
  -> OperatorPlans + SolverPlans
  -> CpuExecutionPlan + CommunicationPlan
  -> FrozenExecutionGraph
```

字段、workspace、stage 和缓存需求可以在初始化 capability registration 阶段登记。
生产 `FieldSchema` 不能在 IBM、湍流、组分和 solver 尚未登记前提前冻结。所有能力登记
完成后，Task 18 按固定顺序一次性冻结布局、arena、halo、线性层和生产执行图。

热路径禁止：

- heap 分配、字符串查找、文本解析、STL/BVH 查询和 donor 搜索；
- 单元级虚函数、插件回调或运行时模型分派；
- 未声明字段写入、隐式 MPI、无 revision 缓存；
- accepted/trial 全字段回滚复制；
- 同一物理量存在两个可写 authority。

## 5. Cartesian 网格与 STL 扫描

### 5.1 CartesianGeometryPlan

v0.4 只有两种静态几何计划：

- `UniformCartesianPlan`：三个方向常量间距；
- `StretchedCartesianPlan`：全局一维 `x/y/z` face coordinates 和 metric arrays 的
  张量积。

JSON 可声明 `focus_regions`、目标尺寸和最大增长率，网格生成器据此生成全局一维坐标。
硬约束包括最大 cells、内存、最小间距和增长率。验证或配对可使用 `exact_cells`，但
不得生成局部非匹配 patch。一个 MPI rank 持有一个连续 `MeshPatch`；rank 内 `CpuTile`
只用于线程调度，不是所有权或通信单位。冷路径的 process-grid 与非均匀切分同时最小化
规则 halo surface、最大 rank cell work 和静态 IBM interface/donor work；仍允许用户显式
覆盖。分解策略、权重和结果进入 case fingerprint，不能为某个算例写入产品分支。

### 5.2 STL 扫描一致性

STL 在初始化期一次读取和规范化。扫描算法保持用户附件方法的数学一致性：

- 轴向 Möller–Trumbore 相交；
- 横向包围拒绝；
- inclusive barycentric 判定；
- 同位置同法向交点去重，反向法向交点保留；
- 排序后使用奇偶规则完成内外标记。

优化只改变数据组织和批处理：triangle bins、scan-line batch、SoA triangle data、一次性
scratch 和线程分块。不会改为逐单元最近三角形比较；无需逐单元与旧曲线坐标方法对比。

## 6. 字段布局、状态与 CPU 执行

### 6.1 权威布局

Eulerian 热字段采用 64-byte 对齐、x 方向连续的 padded SoA；interior 和 ghost 位于
同一分配，向量分量分别连续。三个方向的 face mass flux 分开存储。EB 紧凑接口数据和
小型固定宽度状态可使用 AoSoA。

不长期保存 `rhoU`、`rhoh` 等守恒乘积；在组装中由 `rho` 与 primitive authority 融合
产生。连续运行保留当前与上一时间层所需的 `U、p/h、Y`，派生的 `rho、T、mu、k` 由
EOS/transport revision 管理。不存在常密度专用数组或分支。

### 6.2 Arena、revision 和事务

最终网格、静态几何分类和分解确定后，初始化期按完整 logical graph 的最大需求建立
NUMA-local arenas，划分字段、solver workspace 和 step scratch。页面由其最终 worker/NUMA
owner 并行 first-touch；每线程计数和 scratch 按 cache line 隔离。
时间层通过 handle rotation 推进。每个 attempt 只修改 trial state、pending cache 和紧凑
事务日志；集体成功才发布 accepted revisions，失败直接丢弃 trial/pending 状态，不复制
恢复全场。

### 6.3 CPU 与 MPI

默认每个 NUMA domain 一个 MPI rank，并建立固定线程团队；纯 MPI 保持支持。MPI 只由
固定 communication thread 调用，v0.4 要求并验证 `MPI_THREAD_FUNNELED` 或更高等级，
不依赖 `MPI_THREAD_MULTIPLE`。启动时
选择有限的 ISA、对齐、tile 和 uniform/stretched kernel specialization。正式运行不在
热循环探测硬件或生成无界模板组合。

`hundun tune` 可离线测定 tile、线程、coarsening 和有限 kernel 变体，只允许调整硬件
执行参数；不得调整物理模型、离散 scheme、PISO 次数或收敛阈值。

## 7. Halo 和执行图

`CommunicationEngine` 先提供持久 buffer/request、pack/unpack span 和
`begin -> interior -> finish -> boundary` 调度。规则 Cartesian 算子默认只交换需要的
face slabs；静态 IBM 在 donor plan 完成后编译去重的远端 donor gather。只有确实需要
dense edge/corner ghost 的 stage 才编译完整 grown-box overlap 或固定 staged exchange，
不得假设六个面消息已经填满角点。完整 `CommunicationPlan` 只有在所有 stage 的
field/ghost sets 和 IBM donor reach 登记后才编译，按 peer 和 stage 合并消息，并用 ghost
revision 证明数据有效。每个 persistent request 实例只能 single-in-flight；重叠 stage
必须拥有不同 request 或由冻结图证明串行。

`FrozenExecutionGraph` 为每个 stage 固定：

- read/write sets 和单写者；
- required/published ghost revisions；
- cache dependency tuple 和 invalidation；
- workspace liveness 与别名边界；
- MPI begin/finish、collective epoch 和可合并的 failure-status reduction；
- 允许的 operator refill、hierarchy rebuild 和 cache publish。

资源合同记录每步 allocation 次数、halo bytes/messages、blocking/nonblocking collective
次数与时间、Krylov/MG reduction、内存峰值、operator refill、coarse numeric refresh、
preconditioner setup、linear iterations 和 stage wall time。产品热步必须为零 heap
allocation。失败状态尽量搭载在不可避免的 solver/closed-mass reduction 上；不得为每个
stage timer 或本地成功状态无条件增加 barrier/all-reduce。

## 8. 低速可压缩方程组

### 8.1 状态 authority

主状态为：

```text
U, pi, p_ref, h, Y[0..N-2], passive scalars, final face mass flux
p_abs = p_ref + pi
rho, T, cp, mu, lambda, sound_speed, Mach = derived state
```

EOS 始终使用局部绝对压力 `p_abs`。动量压力梯度使用扰动压力 `pi`，避免大背景压力
损害条件数。开放域由压力边界闭合绝对压力；封闭域使用全局质量约束和包含
`(partial rho / partial p)_(h,Y)` 的 `p_ref` Newton 更新。

### 8.2 热力学与输运

v0.4 支持理想气体惰性混合物，温度相关 `cp(T)、h(T)、mu(T)、lambda(T)`；可为常
物性输入选择计划期专用 kernel，但仍属于同一方程路径。组分采用 N−1 独立质量分数，
最后一个由和约束得到，局部组分参与 EOS 和物性。支持 passive scalars；不支持反应、
PDF、喷雾或真实气体。

### 8.3 守恒方程

- 连续性、动量、静焓、惰性组分和被动标量共享同一最终 face mass flux；
- 静焓是能量 authority，方程包含完整 `Dp/Dt`、热传导和黏性耗散；
- `rhoU/rhoh` 只在融合组装中形成，不成为第二状态；
- 扩散通量、边界通量和源项声明单位、守恒对象、Jacobian/显式部分和适用 stage；
- 所有方程默认 FP64，混合精度只允许在线性求解内部且必须有 FP64 true residual 与
  fallback。

### 8.4 时间步

默认 variable-step BDF2；首次启动、默认 Restart 后第一步和 retry 使用 backward Euler，
随后恢复 BDF2。时间控制模式包括：

- `fixed`；
- `adaptive_flow`（默认）：对流、黏性、热和标量限制；
- `adaptive_acoustic`：用户明确研究声学时启用声速限制。

默认亚音速非声学模式不采用显式声速 CFL 硬限制，但仍监控 sound speed、Mach 和
NSCBC 适用性。

## 9. 边界与离散

边界在冷路径编译成按 patch/stage 的 plan。v0.4 覆盖：

- velocity/mass-flow inlet；
- static/total pressure and temperature inlet；
- pressure outlet、backflow state；
- NSCBC 亚音速入口/出口；
- no-slip/moving/slip/symmetry/periodic；
- adiabatic、isothermal 和 specified heat flux walls；
- 组分与被动标量的 Dirichlet/flux 边界。

第一版 scheme 集只包括二阶中心、limited central、二阶 TVD 和二阶扩散。面重构、
质量通量和散度尽可能融合。所有 scheme 在初始化期静态绑定，内循环不做字符串分派。

浸没边界不是只作用于压力和表面力的附加诊断。初始化期必须把每条流固 link 编译成
唯一的界面面所有权和方程替换行：法向质量通量严格为零；动量按 no-slip/moving-wall
Dirichlet 二次重构替换常规扩散面项；静焓和标量按各自的 value/adiabatic/flux 条件替换；
固体控制体不参与流体守恒、闭域质量或终态连续性范数。替换逐 link 作用于流体行，禁止
通过写共享 solid ghost 的方式实现，因为同一 solid cell 可能对应多个不同壁面点。所需
远程 donor 只用预编译、去重的 compact gather 获取，不扩大常规 face halo。

## 10. 线性系统生命周期

每个线性系统严格拆为四个互不重叠的资源族：

1. `SymbolicPlan/CoarseningPlan`：拓扑、稀疏模式、边界位置、EB 接口、level shapes、
   transfer sparsity 和 line/color schedule；
2. `ExactNumericState`：当前 fine operator、全部必需 coarse numeric coefficients、numeric
   BC 和 diagonal/off-diagonal revision；
3. `PreconditionerSetupState`：smoother factor、coarse solver setup 和可选 HYPRE setup；
4. `SolverWorkspace`：Krylov vectors、MG level vectors、归约缓冲和 scratch。

四族分别拥有 identity 和失效条件。拓扑不变时禁止重建 symbolic/coarsening plan；每个
系数 revision 都必须刷新当前 exact fine/coarse numeric data 后才能认证；预注册的
coefficient-change policy 只能决定 preconditioner setup 是否复用，不能跳过 exact/coarse
numeric refresh；workspace 按最大 shape 持久化。计数器分别记录 symbolic rebuild、
exact numeric refill、coarse numeric refresh、preconditioner setup/reuse 和 workspace
replacement。

规则压力系统默认 `NativeCartesianMG`。近各向同性网格全方向粗化；强拉伸方向采用
semi-coarsening 与 line relaxation。Krylov 提供 PCG、FGMRES 和 BiCGSTAB。HYPRE
Struct 是隔离后备，不成为核心数据模型。IBM 非对称接口使用 exact operator + FGMRES，
紧凑/规则 MG 只能作为 preconditioner。

## 11. 两次 PISO 与中间量生命周期

v0.4 每个 accepted step 恰好执行两次 pressure corrector，不开放用户修改次数。为避免
corrector 2 后再改变 `h/Y` 使 EOS 密度与 continuity 失配，先用 committed face-flux
history、variable-step 系数和已声明 source/diffusion 形成二阶 thermophysical predictor；
startup/Restart/retry 使用 BE predictor。`h*/Y*` 在两次压力校正中冻结：

```text
predict h*/Y* from committed flux history
evaluate EOS/transport; closed-domain p_ref Newton and pi gauge
assemble/update momentum numeric state; predict trial U
corrector 1: build current intermediates -> pressure solve -> trial flux/U
corrector 2: rebuild/revalidate intermediates -> pressure solve -> final flux/U
terminal EOS/continuity/closed-mass/gauge audit
collective commit
```

中间量不能共享一个笼统 dependency tuple：

- `rAU` 依赖动量 diagonal、当前 density、implicit-source diagonal、约束和 boundary
  coefficient revisions；
- `rAtU` 依赖 `rAU` 与一致时间/对角修正 revision；
- `HbyA` 依赖 momentum numeric state、当前 trial `U` 和相关约束；
- pressure face coefficient 依赖 reciprocal diagonal、face-density interpolation、EOS 和
  `drho/dp` revisions；
- `phiHbyA` 依赖当前 `HbyA`、pressure face coefficient、trial/committed flux history、
  BDF、`p_ref`/gauge、geometry 和 numeric boundary revisions。

corrector 1 改写 trial `U/phi` 后，corrector 2 必须重新生成或严格重新认证
`HbyA/phiHbyA`；只有动量系数确实不变时才复用 `rAU/rAtU`。压力修正方程包含局部
`partial rho / partial p`。每次 correction 组装 full BDF density defect 与
`a0*V*(partial rho/partial p)_(h,Y)*delta_pi`，并在压力更新后重新评估 local EOS；不得把
`drho/dp` 当作脱离 BDF density residual 的独立附加项。最终压力方程 flux 是 final face
mass flux 的唯一 writer，最终 `U` 使用同一压力梯度通路更新。提交前分别报告并门控 EOS、
discrete continuity、闭域 total mass 和 pressure gauge/boundary closure residual。

## 12. 静态 IBM 与二阶能力

### 12.1 生命周期

```text
GeometryModel
  -> EBTopology
  -> BoundaryStencilPlan
  -> SurfaceQuadraturePlan
```

初始化期完成 cell classification、intersection、wall point/normal、donor search、quadratic
weights、surface quadrature 和 compact interface indices。热路径只读取计划和 field views。
v0.4 只支持静止 EB；geometry revision 改变才整体重建。

### 12.2 二次重构合同

保留 v0.3 的完整三维二次 cell-average basis、确定性 pivoted QR、condition/rank 门和
无低阶 fallback。标准 link stencil 要求 14--32 fluid donors、最大四层逻辑 halo、全部在
正法向、至少三层 normal bands、覆盖四个 tangential quadrants。共享 quadrature row 若
需要联合 donor set，必须作为不同 plan 类型明确声明，不能悄悄放宽标准 stencil 合同。

EB 与外部 Cartesian 边界相交时采用显式的 boundary-intersection quadratic plan，而非把
标准 link 合同全局放宽。它仍要求 14--32 个唯一 fluid donors、正法向、至少三层 normal
bands、完整十项二次秩与 condition 门；象限证书则必须覆盖该物理域在已声明
`symmetry`/`slip` 边界内几何可达的全部象限，并把 required/actual mask 纳入 plan
fingerprint。其他边界上的单侧 stencil 继续拒绝。surface quadrature 先按 Cartesian 开域
裁剪 STL，域外封口和与外边界共面的面积不属于 immersed traction authority。

v0.3 已证明合法 stencil 上的二次多项式复现，但现有公开 pressure 收敛证据不足以声称
任意几何全局二阶。v0.4 必须增加 uniform/stretched、平移几何和 1/2/4-rank 的真实
h-refinement sequence；二阶能力以预注册的 observed-order 门验收，不靠固定网格复现
代替。

### 12.3 压力与力 authority

IBM 压力采用规则 Cartesian 主体，并按 EB topology 精确移除不可穿透 fluid-solid face
link；solid row 为隔离 identity。该压力离散必须与 final face flux 的同一界面零通量定义
一致，不能在 pressure operator 中另引入 final flux 不会发布的二次 donor mass-flux。
二次重构继续用于速度、扩散和 surface traction；正值材料属性只在严格正 donor 包络内保留
二次值，overshoot 投影到包络，非法 donor 使 attempt 失败。pressure operator、momentum
reaction、final face flux 和 surface force 共享同一 EB topology 和 final-state revision。
surface-force 在进入力的 Allreduce 前先完成全 rank 状态共识，避免局部重构失败导致 collective
次序分叉。

当前 final-gradient/surface-force 候选不会因迁移自动接受，必须先满足：

1. 独立 final-state force oracle；
2. 可区分 final gradient 与 corrector scratch 的 mutation RED；
3. failed attempt 不发布 pending reconstruction/cache；
4. operator force、reaction budget、pressure contribution 和 surface traction 的符号、量纲
   与一致性检查。

`test_immersed_wale_constant` 的 ghost donor positive-normal 失败必须根因修复，不得放宽
donor 几何约束或科学阈值。

## 13. 湍流模型

case schema 只提供三种批准组合：

- `vreman_wall_function`：默认，作用于外部和 IBM 壁面；
- `wale`：保留壁面解析 LES，用于 Re=3900 文献长统计；
- `none`：层流/制造解验证。

内部仍拆为 `SubgridPlan {none,wale,vreman}` 与
`WallTreatmentPlan {resolved,wall_function}`，只允许上述组合。wall treatment 与 subgrid
model 拥有独立 fingerprint/revision，因为前者还控制 heat/species wall flux。WALE 与
Vreman 共享唯一 velocity-gradient cache、wall distance 和 `mu_eff` authority。
`mu_eff` revision 同时约束动量扩散、壁面处理、IBM traction 和诊断。不得为 force、IBM
或不同湍流模型各建一份梯度缓存。

## 14. Restart、输出和运行服务

所有服务只能读取 immutable committed snapshot，不得观察 trial state。输出位于独立
run directory：

- `Restart/`：默认 `keep_last=1`；先写 pending、校验并原子切换 current，再删除旧文件；
- `Visit/`：Cartesian 使用 VTI/VTR 和 `.visit` 索引；
- `screen`：人可读进度、残差、守恒、Mach、资源和错误摘要；
- `monitor/`：JSONL/CSV 力、流量、统计与性能数据。

默认 Restart 保存当前 `U、p_ref/pi（或 p_abs）、h、独立组分、final face mass flux` 及
继续控制所需的最小状态，不保存派生 `rho/T/mu`、长期 `rhoU/rhoh` 或上一 BDF2 历史层。
重启第一步用 BE，允许该步出现可控时间离散差异，之后恢复 BDF2。协议支持 rank-changing
重分解；读取必须先在 scratch 验证全部 rank，再一次性发布。默认同步 publication 在
文件与 generation directory `fsync`、原子 current 切换及 parent-directory `fsync` 完成后
才成功。Restart 后的 BE recovery step 标记在 evidence 中，不进入长期统计样本。

异常热路径只返回紧凑 status；冷路径解释文本。冻结图定义少量固定 collective epoch：
本地 preflight failure 在进入下一项 collective 前合并，已启动的通信先按注册协议完成或
取消，再统一决定；最终所有 rank 同时 commit、retry 或终止。正式 tests-off 二进制不得链接 full-domain gather、
mutation seam 或测试 oracle。

## 15. 实施结构与最终冻结

实施分为五个依赖阶段：

1. 基线与产品骨架：参考、版本隔离、CaseSpec；
2. 性能底座：arena、Cartesian mesh、CPU/halo、边界、热力学、守恒 kernel、graph compiler；
3. 求解与 IBM：线性四层、MG、EB plan、统一方程、两次 PISO、湍流和 final force；
4. 生产冻结与服务：先注册 driver/Restart/I/O 的 snapshot schema、容量和 stage，再由
   每个 `ValidatedModel` 的 `ProductCompiler` 执行 logical analysis、allocation/
   instantiation、binding 和一次 seal；Task 19 只实现 sealed service adapter 与格式；
5. 验收：focused、完整网格短测、候选冻结、20-step 性能、文献长统计、后台阶。

每个 compiled case 实例进入首次 attempt 前必须冻结并证明：单写者 authority、每个 cache 的精确 revision tuple、
accepted/trial 事务、EB plans、persistent halo、四层 solver identity、最大容量 workspace、
完整 stage liveness、两次 PISO 发布语义，以及 allocation/message/byte/collective/
refill/refresh/setup counters。freeze 可以分配 numeric/preconditioner capacity，但初始
field/BC coefficient 填充前必须保持 uncertified；driver initialization 完成 exact/coarse
numeric refresh 和必要 setup 后才可进入第一次 solve。

## 16. Re=3900 圆柱唯一发布门

### 16.1 两个 profile，不混淆用途

**`coast_pairing_short`** 只用于完整网格性能短测：

- `20D x 20D x piD`，圆柱中心距入口 `5D`；
- `480 x 480 x 48` cells，64 ranks，冻结绑核/NUMA；
- `Re_D=3900`，`dt U/D=0.006`；
- 由于 COAST 没有周期边界，HUNDUN 使用与 COAST 可对应的非周期 slip/symmetry 外边界；
- COAST 只提供性能基线，不是 v0.4 数值 oracle，也不参与长统计。

**`literature_statistics`** 只用于 HUNDUN 物理精度：

- 采用 Parnaudeau 等人的 `20D x 20D x piD` Cartesian/IBM 身份、圆柱距入口 `5D`、
  transverse/spanwise periodic 配置和 `dt U/D=0.006`；
- `Re_D=3900`、WALE、无算例专用源码；
- 至少发展 `150D/U`，随后按论文约 `2020D/U`（约 420 shedding periods）统计；
- 运行前冻结论文表格、数字化剖面、单位、误差和提取脚本哈希。

### 16.2 强制执行顺序

1. focused：EOS、焓、PISO、守恒、IBM 二阶、force oracle/mutation、WALE donor、Restart、
   1/2/4-rank、ASan 和 UBSan；
2. 完整 `480x480x48/64-rank` 2-step HUNDUN/COAST 短测，验证真实初始化、内存和性能方向；
3. 满足正确性且存在可信性能方向后冻结 exact candidate；
4. 冻结候选进行至少五组交替完整网格 20-step 短测；
5. 只有短测性能接受后才运行 HUNDUN 文献长统计；
6. 完成 provenance、完整 diff、进程和证据审查，输出一个 `ACCEPT` 或 `REJECT`。

取消 `24^3` 性能门。小网格只服务 focused 数值测试，不产生性能结论。

### 16.3 判据

```text
release = numerical_correctness_accept
       && robustness_accept
       && coast_short_performance_accept
       && literature_physical_accuracy_accept
       && provenance_accept
```

- 每步恰好两次 PISO，最终 continuity、守恒、true residual 和 transaction 合同通过；
- final force oracle 和 mutation RED 接受，positive-normal donor 失败已按根因修复；
- 正式性能只用 tests-off binary；记录 init、hot-step median/P90、RSS、communication、
  iterations 和 rebuild/refill counters；
- 正式 timed region 在两侧关闭 Restart/Visit/screen serialization，预注册 warmup 与测量
  steps，以 launcher elapsed 或 max-rank elapsed 为 authority；内部 timer 不引入 barrier；
- 至少五组交替运行按同一资源/时间块形成 paired ratio
  `HUNDUN_hot/COAST_hot`，报告 median paired ratio、P90 和预注册 uncertainty interval；
  接受时 paired ratio 不高于 `1.0` 且 uncertainty rule 通过；`1.10x` 只能记作内部
  `NEAR`，不构成发布；
- 文献主要指标为 `St、Lr/D、Cd_mean、Cl_rms`、centerline mean velocity，以及
  `x/D=1.06、1.54、2.02` 的 mean/fluctuation profiles；阈值和数字化误差必须在运行
  前登记，禁止查看结果后修改；
- COAST 与 HUNDUN 方程和周期能力不同，不要求全场 `1e-10` 等价，也不拿 COAST
  长统计决定 HUNDUN 精度；
- exact HEAD/tree、compiler/flags、binary、case、STL、rank map、CPU topology、命令和
  evidence hashes 全部封存；不得覆盖历史证据。

较快初始化、较低单次残差或 2-step 较快均不能替代 hot-step median 和文献统计。

## 17. 后台阶验证

圆柱发布门通过后，采用 Driver--Seegmiller/NASA 后台阶文献配置验证通用性：

- case 只通过根目录 JSON/`.d` 输入构建，不新增后台阶专用源码；
- 使用全局张量积拉伸 Cartesian 和默认 Vreman + wall function；
- 按文献设置 step height、expansion ratio、入口边界层和低速气体状态；
- 比较 reattachment length、`Cf/Cp`、mean velocity 和 Reynolds-stress profiles；
- 该验证不阻塞本次 Re=3900 发布节点；若后台阶暴露源码缺口，修复必须保持通用，且
  完整重跑圆柱发布门后才可继续声称 v0.4 已接受。

v0.4 正式验证只包括 Re=3900 圆柱和后台阶，不以射流或其他算例扩大本版门。

## 18. 不在 v0.4 范围内

- 常密度快速路径、SIMPLE、PIMPLE、strong PIMPLE；
- body-fitted、多块曲线网格、AMR、非匹配 refinement patches；
- 激波、超声速、声学研究默认模式、真实气体；
- 反应、化学、PDF、TPDF/TCR、喷雾和多相；
- 除 WALE、Vreman + wall function 之外的生产湍流模型；
- 移动 EB、动态负载迁移、生产 GPU 后端；
- 热循环细粒度插件 ABI；
- 源码树内产品 cases/examples；
- 为圆柱或后台阶编写特化求解流程；
- 候选冻结前的长统计和小网格性能结论。

## 19. 完成条件

设计实现完成必须同时满足：

1. v0.3 可独立构建且现有脏圆柱工作树与证据未被改动；
2. v0.4 默认构建采用扁平、独立的 Cartesian 低速可压缩架构；
3. 所有热路径符合 authority、revision、allocation、communication 和 solver lifecycle 合同；
4. 更换合法 case 不修改产品源码；
5. Re=3900 联合门按规定顺序返回 `ACCEPT`；
6. 圆柱通过后，后台阶可由纯输入构造并按文献验证。

除 Re=3900 联合门外，不设置 alpha、beta、RC 或其他公开发布节点。

# HUNDUN-FLOW v0.4 Cartesian 性能架构设计

日期：2026-08-12

状态：已批准设计，书面复核已确认

目标版本：HUNDUN-FLOW v0.4
唯一发布门：Re=3900 三维圆柱绕流数值正确性门与性能门同时通过

## 1. 目标与约束

v0.4 是性能优先、允许破坏旧内部 API 的新主线。目标是在保持可验证数值正确性
的前提下，建立适合湍流 CFD 的通用、鲁棒、高效率底层，而不是继续在 v0.3
热路径上叠加局部补丁。

本版采用以下强约束：

1. 只实现 CartesianGeometryPlan，包括均匀 Cartesian 网格和各坐标方向可拉伸
   的 Cartesian 网格；不支持 body-fitted、多块曲线网格或 AMR。
2. Linux CPU 是正式生产后端；优先支持 MPI、NUMA、固定线程团队和 SIMD。
   GPU 与 AMR 只保留粗粒度扩展边界，不进入 v0.4 验收范围。
3. v0.3 冻结为备份、输入迁移参考和独立数值 oracle；v0.4 是默认构建与运行主线。
4. 每个版本使用独立目录。v0.4 源码采用扁平目录和稳定前缀，减少跨目录跳转；
   扁平化不得牺牲模块所有权和依赖方向。
5. 只复用公开数学、数据布局和生命周期思想；不得复制 OpenFOAM GPL 源码，
   不得复制旧 COAST Fortran，也不得把上游默认参数当作 HUNDUN 调参依据。
6. 候选冻结前不启动长统计。性能优化不得减少 PISO corrector、放宽收敛阈值、
   增加经验阻尼或滤波，也不得逐算例调参。

## 2. 公开参考边界

设计参考固定到以下公开或用户授权只读基线：

| 项目 | 固定 revision / 位置 | 允许参考 | 许可证边界 |
| --- | --- | --- | --- |
| OpenFOAM-dev | `b9da51ab0673423aa2af6a45a72a3fbec9c66f9f` | PISO/PIMPLE 方程时序、`rAU/rAtU`、`HbyA`、`phiHbyA` 和最终 `U/phi` 权威关系 | GPL，只参考公开数学和生命周期，不复制实现 |
| AMReX | `59d066aab774bc388cc6ed944f7beaf645607ed3` | box/field 数据布局、halo/fill-patch、持久化通信与资源管理 | BSD-3-Clause，仍采用独立实现 |
| IncFlo | `7307d8725c2a538f09cafbeacbfeb63e0fb11d22` | projection、EB 与求解阶段组织 | BSD-3-Clause，仍采用独立实现 |
| AMReX-Hydro | `e49df248aabd2cc11865eb5be734a2f5f2f65ee5` | operator、multigrid hierarchy 和系数生命周期 | BSD-3-Clause，仍采用独立实现 |
| COAST | `/home/wyf/code_dev/Coast_software` | SIMPLE/ICCG 时序、数组布局、预生成数据和热循环复用 | 用户指定只读参考，不复制 Fortran |

HUNDUN-FLOW 保持独立方程语义和实现。COAST 的 SIMPLE 外迭代、松弛和压力修正
次数不能机械移植到 PISO；AMReX 的 AMR 层级语义和 OpenFOAM 的运行时对象模型
也不成为 v0.4 的依赖。

## 3. 产品目录与构建边界

仓库采用并列版本目录：

```text
hundun-flow/
  CMakeLists.txt
  README.md
  VERSION
  versions/
    v0.3/
      CMakeLists.txt
      include/hundun/
      src/
      tests/
      cases/
      docs/
    v0.4/
      CMakeLists.txt
      include/hundun/
      src/
      tests/
      cases/
      docs/
```

根构建默认选择 v0.4，并允许显式选择 v0.3。两个版本不共享产品源文件，避免
无意间改变备份 oracle。可共享的内容仅限仓库级许可证、构建辅助和不可执行的
证据说明。

`versions/v0.4/src` 是统一扁平源码目录。文件以职责前缀分组：

```text
platform_*    CPU、NUMA、MPI、ISA 选择
storage_*     FieldSchema、Arena、view、revision
mesh_*        CartesianGeometryPlan、MeshPatch、CpuTile
comm_*        halo、persistent request、通信计算调度
boundary_*    BoundarySpec、BoundaryPlan
eb_*          EBTopology、stencil、surface quadrature
disc_*        面通量、梯度、扩散、时间离散
operator_*    动量、压力、标量算子及系数更新
linear_*      Krylov、CartesianMG、HYPRE 适配
model_*       密度、输运、WALE、反应与源项贡献
integration_* 方程、耦合、时间步、事务和执行图
service_*     Restart、I/O、诊断、性能证据
app_*         配置解析、校验、启动和命令行
```

公开安装头保留在 `include/hundun`；非公开声明与实现尽量留在扁平 `src`。
源码位置不影响运行性能，性能由编译单元边界、链接时优化和数据访问决定。

## 4. 初始化流水线与冻结执行计划

运行前的唯一构造链是：

```text
CaseSpec
  -> ValidatedModel
  -> FieldSchema + CartesianGeometryPlan
  -> BoundaryPlan + EBTopology
  -> OperatorPlans + SolverPlans
  -> CpuExecutionPlan + CommunicationPlan
  -> FrozenExecutionGraph
```

`CaseSpec` 只表示用户输入；`ValidatedModel` 完成单位、物理组合、边界闭合、字段
需求和求解器兼容性检查。字段允许在初始化阶段动态注册，但 `FieldSchema` 一旦
冻结，热循环不得再增删字段、改布局或触发重新分配。

`FrozenExecutionGraph` 把每个 stage 的读集、写集、所需 ghost revision、缓存
失效条件、通信依赖、工作区和 collective failure 点编译成稳定计划。热路径禁止：

- heap 分配、字符串查找、几何搜索和 donor 搜索；
- 每单元虚函数或插件回调；
- 未声明写入、未带 revision 的缓存和隐式 MPI；
- trial/accepted 状态的全字段复制；
- 同一物理量存在两个可写 authority。

## 5. 网格、分解与数据布局

### 5.1 CartesianGeometryPlan

v0.4 只有两种初始化期静态方案：

- `UniformCartesianPlan`：常量间距，优先采用编译期/计划期简化；
- `StretchedCartesianPlan`：每个坐标方向使用一维坐标和 metric 数组。

几何、拓扑与热字段分离。一个 MPI rank 持有一个连续 `MeshPatch`，静态分解可
非均匀以平衡 EB 工作量。rank 内使用尺寸可调的 `CpuTile` 分配线程工作；tile
不是所有权或通信单位。

### 5.2 权威字段布局

Eulerian 热字段使用 padded SoA，x 方向连续，interior 与 ghost 位于同一分配。
向量分量各自连续。局部 stencil、EB 重构和小型逐单元状态可使用固定宽度 AoSoA。

速度、压力、密度、标量使用 collocated cell state；三个方向的 face flux 分开存储，
最终 conservative face flux 是质量守恒和后续输运的唯一通量 authority。派生量要么
即时融合计算，要么拥有明确 revision 和失效规则。

### 5.3 内存与 NUMA

每个 NUMA domain 默认一个 MPI rank，并建立持久线程团队；纯 MPI 是兼容模式。
大块 NUMA-local arena 在初始化期分配，按字段、solver workspace 和 step scratch
划分。时间层用 handle rotation，失败步用局部 transaction log 回退，避免全场复制。

## 6. 通信与 CPU 执行

`CommunicationPlan` 在初始化期合并同一邻居、同一阶段需要的字段，建立持久 MPI
request、pack/unpack 区间和 ghost revision 合同。执行顺序由图显式规定：启动 halo、
计算不依赖 ghost 的 tile、完成通信、计算边界 tile。不会在算子内部偷偷触发 halo。

`CpuExecutionPlan` 在启动时选择有限数量的 kernel specialization，包括 ISA、对齐、
tile 尺寸、均匀/拉伸网格和常见 stencil 组合。内层循环保持静态分派、可向量化和
确定的别名约束。不会生成无界模板组合，也不会在热循环做能力探测。

性能资源合同包括：每步字节搬运预算、halo 字节与消息数、融合边界、arena 峰值、
solver iteration 和 stage 时间。违反合同属于性能回归证据，而不只是诊断信息。

## 7. 离散、方程与压力速度耦合

### 7.1 方程 authority 与调度分离

方程对象定义离散数学、系数、边界贡献和读写集合；耦合调度只决定何时执行和
复用哪些已认证结果。PISO、SIMPLE 或其他算法若以后出现，必须调用同一方程
authority，不能复制一套离散实现。v0.4 的 Re=3900 正式路径固定为两次 PISO
corrector。

### 7.2 统一中间量语义

动量预测器和压力修正使用统一生命周期：

- 组装或更新动量算子的数值系数；
- 生成对角逆语义 `rAU`，需要修正对角时生成 `rAtU`；
- 从动量非对角与显式项生成 `HbyA`；
- 从 `HbyA` 和一致面插值生成 `phiHbyA`；
- 压力方程修正 face flux，并由最终压力通路更新 `U`；
- corrector 结束后只发布最终 `U`、最终 face flux 和相应 revision。

这些名称描述 HUNDUN 自己的数学合同，不复制 OpenFOAM 类型或源码。每个中间量
必须注明创建 stage、依赖 revision、可跨 corrector/时间步复用条件和失效原因。

### 7.3 通量与 scheme

所有守恒方程共享 direction-separated face flux authority。离散选项在初始化期编译成
`SchemePlan`，运行时不做字符串分派。第一版只实现 Re=3900 及后续通用湍流必需的
有限 scheme 组合，优先融合面重构、质量通量与散度，避免重复读取整个字段。

## 8. 线性算子与求解器生命周期

每个线性系统分为四层：

1. `SymbolicPlan`：拓扑、稀疏模式、边界位置和 EB 接口集合；
2. `NumericState`：随系数 revision 更新的数值系数；
3. `HierarchyState`：multigrid hierarchy、restriction、prolongation 和 smoother；
4. `SolverWorkspace`：Krylov 向量、归约缓冲和临时数组。

四层分别拥有 identity 与失效条件。拓扑未变时不得重建 symbolic plan；系数未变时
不得重填；层级可依据量化的 coefficient-change policy 复用；workspace 按最大需求
持久化。

规则 Cartesian 压力系统默认使用 `NativeCartesianMG`。均匀或近各向同性网格采用
全方向 coarsening；强拉伸方向使用 semi-coarsening 与 line relaxation。HYPRE Struct
作为独立适配后备，不成为核心数据模型。

动量系统默认使用精确张量动量算子配合标量 MG 预条件器；需要时保留 block 路径。
IBM 压力使用“规则 Cartesian 主体 + 紧凑接口修正”的 exact authority，外层使用
FGMRES；紧凑 MG 只能作为 preconditioner，不能替代非对称接口算子。

## 9. EB/IBM 生命周期

静止几何的初始化链为：

```text
GeometryModel
  -> EBTopology
  -> BoundaryStencilPlan
  -> SurfaceQuadraturePlan
```

初始化期完成分类、最近表面查询、donor 搜索、权重、法向、面积、quadrature 与接口
索引压缩。热路径只读取紧凑 plan 和字段 view，不做 STL/BVH 查询。几何 revision
变化才允许整体重建；v0.4 正式范围是静止 EB。

压力 operator、最终 face flux、动量 reaction 和 surface force 共享同一 EB topology
与最终状态 authority。当前 final-gradient/surface-force 候选不能因架构迁移自动获得
接受，必须先通过：

1. 独立 final-state force oracle；
2. 能区分 final gradient 与 corrector scratch 的 mutation-sensitive RED；
3. 失败 attempt 不发布 pending reconstruction 的 transaction 检查；
4. operator force、budget reaction、pressure contribution 和 surface traction 的符号、
   量纲与一致性检查。

既有 `test_immersed_wale_constant` ghost donor positive-normal 失败必须按根因修复，
不得放宽 donor 几何约束、科学阈值或 mutation 灵敏度。

## 10. 湍流、物性与多物理贡献

湍流模型在初始化期静态绑定。WALE 和以后模型通过共享派生字段计划读取速度梯度，
避免各模型重复重构。`mu_eff` 只有一个 authority，其 revision 同时约束动量扩散、
壁面处理、EB traction 和相关诊断。

密度、能量、组分、反应和源项通过 conservative contribution 接口加入既有方程：
贡献必须声明单位、守恒对象、Jacobian/显式项、读写集合和适用 stage。核心热循环
采用静态 composition；动态插件只允许在粗粒度模型构造或离线服务边界出现，不能
向单元循环注入虚调用。

数值权威默认 FP64。混合精度只允许用于有独立残差校验和 FP64 fallback 的 solver
内部，并且必须单独通过 Re=3900 数值门。归约提供 reproducible 模式和 performance
模式；Restart 与明确协议测试仍执行 bitwise 门。

## 11. 时间推进、错误处理与服务隔离

`TimeSchemePlan` 明确 accepted、trial 和历史时间层。每个 attempt 只修改 trial、局部
缓存和事务日志；所有 rank 在 stage 边界进行一次 collective consensus。成功后旋转
arena handle 并发布 revision，失败则丢弃 trial 与 pending cache，不复制恢复整场。

热内核返回紧凑 status code，不构造异常字符串。错误上下文在冷路径解释，并在约定
stage 做 collective 汇总，保证所有 rank 作出相同 commit/rollback 决定。

Restart、可视化输出、统计、调试 oracle 和性能证据属于 `service_*`，不进入算子
authority。异步或聚合 I/O 只能读取已提交快照；测试访问宏不得把全域 gather 或宽
记录注入正式产品二进制。

## 12. 开发顺序

以下是内部开发检查点，不是发布节点：

1. 冻结 v0.3：只从 `4ae4832ad00b5e4d1129ee978e2e49cbb33bb7ef` 的 tracked tree
   建立版本目录，不纳入该 worktree 的 dirty/untracked 内容；保留现有脏圆柱工作树
   和全部证据，不 reset、clean、覆盖或提交其中修改。
2. 建立 v0.4 基座：版本构建、CaseSpec、FieldSchema、arena、view、revision、
   CartesianGeometryPlan、MeshPatch、CpuTile。
3. 建立执行底层：NUMA/线程计划、persistent halo、FrozenExecutionGraph、状态码和
   transaction。
4. 建立 Cartesian 离散：边界计划、面通量、梯度、扩散、SchemePlan 和守恒标量。
5. 建立线性层：四层生命周期、NativeCartesianMG、Krylov 和 HYPRE Struct 后备。
6. 建立统一流动层：动量、`rAU/rAtU`、`HbyA/phiHbyA`、两次 PISO、最终 `U/phi`。
7. 建立 EB/IBM：静态 geometry plan、接口修正、force authority 和独立 oracle。
8. 接入 WALE、变密度与共享 `mu_eff`；随后迁移反应和更多标量能力。
9. 接入 Restart、I/O、诊断、配置迁移和性能证据服务。
10. 冻结性能候选，依次执行正式门；未冻结前不启动长统计。

每个检查点都必须保留可构建、可审查的纵向切片，但不得以“alpha”“beta”“RC”或
其他名称对外发布。

## 13. 唯一发布门：Re=3900 数值与性能联合门

v0.4 不设置分散的发布里程碑。只有下列联合判定全部满足，才允许发布：

```text
release = numerical_correctness_accept
       && robustness_accept
       && coast_performance_accept
       && physical_accuracy_accept
       && provenance_accept
```

### 13.1 候选身份与前置条件

- 冻结 HEAD、tree、编译器、构建参数、二进制、case、STL、rank 映射和证据清单；
- tests-on 与 tests-off 产品路径分离，正式性能只使用 tests-off；
- 每个 HUNDUN case 绑定同资源、等价设置的 COAST pairing；
- final-state force oracle 与 mutation RED 已先接受；已知 WALE donor 失败已按科学
  约束修复；
- 不存在本任务遗留的 MPI/测试进程，证据目录可追溯且不覆盖旧运行。

### 13.2 强制执行顺序

1. focused：单元、mutation/RED、Restart/rollback、`1/2/4-rank`、sanitizer 和公共
   合同；
2. `24^3/20-step`：稳定性、数值等价、迭代趋势和低成本性能 screen；
3. 完整 `480x480x48/64-rank` `2-step`：与 COAST 成对计时，验证真实数据布局、
   初始化和首个热步；
4. 候选满足前三项且性能方向冻结后，才允许完整网格 `20-step` pairing；
5. 只有冻结候选才运行充分发展的长时统计，并完成实验精度比较。

任一步失败都回到对应实现检查点；不得跳级，也不得用长测替代 focused 证据。

### 13.3 数值正确性

- 每步恰好两次 PISO corrector；IBM/WALE、边界、时间阶次、收敛阈值、守恒、
  rollback、collective failure 和 MPI 语义不变；
- velocity、pressure、face flux 相对已验证 HUNDUN reference 的归一化 L2 差
  `<= 1e-10`，近零量使用预先登记的物理绝对容差；
- 汇总关键有限值相对差 `<= 1e-10`、绝对差 `<= 1e-12`；
- Cd、Cl 与四字段力满足绝对差 `<= 1e-11` 或相对差 `<= 1e-9`，近零 Cl 以绝对
  容差为主；真实 correctness repair 由独立 oracle 建立新 reference，不能维持错误旧值；
- continuity、pressure residual 和守恒缺陷通过原科学阈值，且不系统性超过 reference
  的 `1.05` 倍；
- 20 步总线性迭代不得系统性增加超过 `2%`，不接受 NaN、崩溃、死锁或未处理异常；
- Restart、相同分区续算和明确协议测试保持 bitwise；其他数学等价性能优化允许受控
  舍入差异。

### 13.4 性能与物理精度

- 相同核数、等价网格、时间步和尽可能一致的物理设置；记录总时间、初始化时间、
  热步 median/P90、峰值 RSS、通信量和 solver iterations；
- HUNDUN 热步中位耗时必须不高于 COAST；`1.10x` 只作为内部接近线，不构成发布；
- HUNDUN continuity、守恒和稳定性不得劣于 COAST；
- `Cd_mean`、`Cl_rms`、`St`、压力分布或尾流统计相对实验的误差不得大于 COAST；
  至少一个主要指标明确改善后，才允许宣称精度优于 COAST；
- 较小线性残差、较快初始化或较快前两步均不能单独替代热步和物理统计结论。

### 13.5 发布判定

联合门只输出 `ACCEPT` 或 `REJECT`。内部检查点、短程 screen、单一性能提升和单一
物理指标均不得产生发布标签。`ACCEPT` 报告必须绑定 exact candidate、全部配对证据、
数值 oracle、性能统计、实验比较、DCO/许可证审查和最终干净工作树状态。

## 14. 不在 v0.4 范围内

- body-fitted、多块曲线网格和 AMR；
- 生产级 GPU 后端；
- 动态负载迁移和移动 EB；
- 热循环细粒度插件 ABI；
- 为单一 Re=3900 case 编写不可复用的特化求解流程；
- 在性能候选冻结前启动长时统计；
- 未通过联合门的公开版本发布。

## 15. 设计完成条件

本设计在以下事实同时成立时视为实现完成：v0.3 可独立构建并作为冻结 oracle；v0.4
默认构建使用扁平、独立的 Cartesian 性能架构；热循环符合资源与 authority 合同；
Re=3900 联合门按规定顺序返回 `ACCEPT`；除此之外没有其他发布节点。

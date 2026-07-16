# HUNDUN-FLOW 全新 C++ 燃烧 CFD 软件总体设计与开发交接

日期：2026-07-16

状态：总体路线、名称、首版范围和版权边界已确认；等待用户书面审阅后编制实施计划

文档性质：私有设计与开发交接文件，不作为 HUNDUN-FLOW 公共仓库的历史沿革说明

## 1. 项目定位

HUNDUN-FLOW 是一款从零建立、版权边界清晰、适合 GitHub 开源的 C++ 燃烧 CFD 软件。它面向燃烧模拟研究人员，重点支持通过稳定接口开发新的湍流、燃烧、蒸发、颗粒和方程源项模型。

项目不是对 BOFFIN、COAST 或 COAST-2 的重命名、翻译或渐进重构。现有 COAST 继续承担当前科研计算，HUNDUN-FLOW 在隔离的新仓库中独立开发。现有程序只提供科学需求、算例行为和整体趋势参照，不能成为源代码实现模板。

毕业阶段的目标功能与当前 COAST 大致一致，但不要求逐行、逐格式或逐浮点结果一致。首要功能范围包括：

- 亚音速可压缩、低马赫数变密度反应流；
- 结构化和曲线坐标网格；
- IBM 与 LES 基线；
- 有限速率化学反应；
- 与 TCR 紧密结合的 TPDF；
- IEM 作为 TCR 参数极限；
- 喷雾、蒸发、颗粒迁移和双向耦合；
- MPI 并行、Restart、诊断与可视化输出；
- 面向研究模型的插件和源项扩展接口。

首版不以超音速、激波捕捉、多部件运动 IBM 或 WENO/DG 实现为交付目标，但架构不得封死这些后续方向。

## 2. 名称与文化含义

项目正式名称为：

> **HUNDUN-FLOW**

中文名称为：

> **浑敦流体模拟框架**

正式仓库建议使用 `hundun-flow`，主可执行文件使用 `hundun`，C++ 根命名空间使用 `hundun`。

HUNDUN 可作技术释义：

> **H**igh-fidelity **U**nified **N**umerics for **D**ensity-varying **U**nsteady **N**avier-Stokes

项目英文副标题为：

> An extensible C++ framework for immersed reacting-flow simulation

《山海经·西山经》以“其状如黄囊，赤如丹火……浑敦无面目”描述帝江。项目介绍必须区分古籍原文和现代工程阐释：古籍提供“黄囊”和“丹火”的原始意象；HUNDUN-FLOW 将其工程化地理解为一个由边界包围的囊状密闭腔体，内部承载炽烈燃烧、湍流输运与相互作用的多物理过程。这一意象与燃烧室数值模拟相呼应。原文电子版本参见中国哲学书电子化计划：`https://ctext.org/text.pl?if=gb&node=83583&show=parallel`。

“浑敦”同时表示尚未被分解和建模的复杂整体。数值模拟的任务是从这种复杂性中建立守恒关系、可识别结构和可验证的物理认识。

公共 README 可使用以下中文短文：

> HUNDUN-FLOW（浑敦）之名源于《山海经》中的帝江意象。古籍以“状如黄囊，赤如丹火”描述其形态。本项目将这一意象工程化地理解为一个由边界包围的囊状密闭腔体，内部承载炽烈燃烧、湍流输运与多物理过程，恰与燃烧室数值模拟的研究对象相呼应。“浑敦”也象征复杂流动尚未被分解和建模的原始状态，而数值模拟的任务正是从这种复杂性中建立守恒、结构与可验证的物理认识。

公共 README 可使用以下英文短文：

> HUNDUN-FLOW takes its name from Hundun, an archaic image associated with Dijiang in the Classic of Mountains and Seas, described as pouch-like and red as cinnabar fire. The project reinterprets this image as a bounded chamber containing intense combustion, turbulent transport, and interacting physical processes. It also represents the initially undivided complexity of reacting flows, from which numerical simulation seeks to recover conservation, structure, and verifiable physical understanding.

单独的 `HUNDUN` 已被其他软件包使用，因此正式项目、论文标题和 GitHub 仓库统一使用 `HUNDUN-FLOW`；正文中完成首次定义后可简称 HUNDUN。

## 3. 版权与独立性边界

### 3.1 固定法律比较基线

内部版权审计只使用以下目录作为 BOFFIN 最终法律基线：

```text
/home/wyf/Compress_boffin/AECSC_WDQ/SRC-Spray
```

不再寻找或增加其他 BOFFIN 比较版本。该基线包含的内容可能比 W. P. Jones 教授持有的早期纯 Fortran BOFFIN 更广，因此采用它属于更严格的独立性检查。

### 3.2 现有程序的角色

- BOFFIN 基线：只用于私有版权相似性审计和禁止项识别；
- COAST/COAST-2：只用于私有科学需求、算例、输入输出和黑盒结果参照；
- HUNDUN-FLOW：新仓库、新架构、新数据模型、新方程组织和新实现。

不得从 BOFFIN、COAST 或 COAST-2 复制、翻译、逐句改写或机械重构以下内容：

- Fortran 求解器过程及其控制流；
- 离散方程的旧程序表达和数组更新顺序；
- 边界、MPI、Restart、TPDF、IBM、喷雾和化学耦合的 legacy 适配层；
- 全局变量布局、公共块、Fortran 下标语义和 ABI；
- 注释、错误信息、输出格式和内部命名，只要它们体现旧实现结构。

方程、算法和经验模型必须从论文、教材、标准或许可证兼容的公开资料重新推导，并在设计记录中注明来源。不能把观察旧代码当作算法来源。

### 3.3 C++ 资产的选择性复用

现有 COAST 中由用户独立编写的通用 C++ 资产可以逐文件审计，但默认动作仍是重写。只有同时满足以下条件才可进入 HUNDUN-FLOW：

1. 可证明由 WANG YUDONG 独立创作；
2. 不包含 BOFFIN/Fortran 数据布局、调用顺序或适配逻辑；
3. 对固定 BOFFIN 基线通过文本、标识符和结构相似性检查；
4. 接口符合 HUNDUN-FLOW 新架构，不为复用而扭曲架构；
5. 在私有审计报告中记录文件哈希、结论和审查日期。

候选范围仅限通用类型、FieldView、MPI 域分解思想、Halo schedule、CPU affinity、FPE、网格或 Restart 数据类型等。以下内容必须重写：`domxch`、Fortran indexing/ABI、`local_domain`、旧 runtime mesh、`coalesced_legacy_block`、旧 Decomp/Restart 适配器和旧 `.d` 输入路径。

公共仓库不发布 `REUSE_MANIFEST.yaml`，因为前置程序不会对外发布。私有证据建议保存为 `BOFFIN_INDEPENDENCE_AUDIT.md`，但不得把 COAST 描述为 HUNDUN-FLOW 的公共上游。

### 3.4 公共许可证

- 主许可证：Apache License 2.0；
- 贡献约束：Developer Certificate of Origin，使用 Signed-off-by；
- 第三方许可证：保存在 `LICENSES/` 和 `THIRD_PARTY.md`；
- 不采用 CLA 作为首版贡献门槛。

公共 `NOTICE` 仅包含：

```text
Copyright (c) 2026 WANG YUDONG
Contact: wangyudong@buaa.edu.cn
```

NOTICE 中不写项目、资助或前置软件说明。

## 4. 设计原则

1. **科学功能与代码实现分离**：可以要求物理能力等价，不能要求实现过程相似。
2. **守恒优先**：质量、动量、能量、组分和两相交换必须显式追踪守恒关系。
3. **模型可替换**：研究模型通过类型化服务和源项接口接入，不直接控制 MPI 或内部数组。
4. **依赖轻量**：并行、域分解、Halo、粒子迁移和基础线性求解器自主实现，仅依赖标准 MPI。
5. **首版可完成**：优先建立可靠二阶有限体积基线；WENO、DG 和运动 IBM 只预留接口。
6. **测试先于大算例**：每个阶段通过独立数值基准后才进入综合燃烧室。
7. **错误显式化**：不支持的模型组合、几何和配置必须明确失败，不静默降级。
8. **无 Python 运行依赖**：构建和运行主程序不得要求 Python、Conda 或在线下载。

## 5. 总体软件架构

建议目录结构：

```text
applications/
  combustor/
sdk/
  include/hundun/sdk/
models/
  turbulence/
  combustion/
  tpdf_tcr/
  particles/
  boundaries/
solver/
  flow/
  transport/
  pressure/
  time/
  chemistry/
mesh/
  structured/
  metrics/
  ibm/
runtime/
  mpi/
  decomposition/
  exchange/
  fields/
  particles/
  io/
tests/
examples/
docs/
LICENSES/
```

依赖方向固定为：

```text
applications -> models/sdk -> solver/mesh -> runtime
```

下层不得反向依赖应用或具体模型。核心模块职责如下：

- `applications`：JSON 配置、对象组装、主时间循环和运行生命周期；
- `sdk`：稳定模型接口、FieldView、服务对象、源项累加器和注册表；
- `models`：LES、燃烧、TPDF-TCR、喷雾、边界和 IBM 模型；
- `solver`：控制方程、离散算子、压力速度耦合、时间推进和化学后端；
- `mesh`：结构化块、曲线坐标度量、IBM 几何查询；
- `runtime`：MPI、域分解、Halo、字段所有权、粒子迁移、I/O 和诊断。

外部插件提供稳定 C ABI；仓库内 C++ SDK 在 C ABI 上封装类型安全接口。外部 `.so` 插件要求与主程序使用兼容的 C++ 编译器、MPI 和 ABI 版本。

## 6. 数据所有权与运行时

### 6.1 字段系统

字段由 `FieldRegistry` 声明并由运行时统一分配。模型只获得有生命周期限制的 `FieldView`，不得持有裸指针跨越网格重建、Restart 或重新分区。

每个字段必须声明：

- 名称、单位和数值类型；
- `FunctionSpace`；
- 分量数和守恒属性；
- Ghost 宽度和交换需求；
- Restart 和输出策略；
- 所有者与只读/可写权限。

首版 `FunctionSpace` 至少预留：

```text
CellAverage
FaceValue
VertexValue
ElementDof
QuadraturePoint
Particle
```

首版实际实现以 `CellAverage`、`FaceValue` 和 `Particle` 为主；其余用于高阶接口编译和 mock 测试。

### 6.2 MPI 与域分解

基础并行框架自主实现，依赖 MPI-3，不引入 PETSc、OpenFOAM runtime、Kokkos 或通用分布式网格框架。

运行时至少包含：

- 结构化域分解和全局/局部索引；
- 任意 Ghost 宽度的 `ExchangePlan`；
- 面、边、角邻居通信；
- 非阻塞 Halo 调度；
- 粒子跨 rank 迁移；
- MPI 错误传播和一致终止；
- 分区无关的全局 cell ID。

`ExchangePlan` 不得硬编码二阶格式的单层 Ghost。接口必须能表达未来 WENO 宽模板和 DG 面 trace DOF。

### 6.3 随机数

TPDF 随机数由以下键确定：

```text
global_cell_id + stochastic_field_id + time_step + user_seed
```

同一物理问题在不同 MPI 分解下应得到一致的随机序列。Restart 必须恢复足以继续同一序列的状态。

## 7. 控制方程与数值基线

### 7.1 物理范围

首版求解亚音速可压缩、低马赫数、变密度反应流。压力速度耦合负责声学刚性过滤和质量守恒，不要求超音速、激波或全可压缩 Riemann 求解器。

主要守恒量包括：

- 总质量；
- 三分量动量；
- 显热焓或与热化学后端一致的能量变量；
- 物种质量分数；
- 可选用户标量；
- 两相交换的质量、动量和能量源项。

### 7.2 首版离散

- 空间框架：现代二阶守恒、Cell-centered finite volume；
- 动量对流：具有动能友好性质的中心通量，在 IBM 邻域和坏网格处使用有界回退；
- 组分、焓和标量：MUSCL 重构，默认 MC limiter；
- TPDF 随机场：有界 MUSCL、正性保护和守恒归一化；
- 黏性与扩散：二阶中心格式和非正交修正；
- 压力速度耦合：Rhie-Chow 插值和变密度 PISO；
- 时间推进：BDF2 和自适应 CFL；
- 化学耦合：支持 Strang splitting；
- 线性求解：自主 CG/BiCGStab 和轻量预条件器。

旧 BOFFIN/COAST 离散格式不作为代码或控制流模板。数值公式应从公开文献重新推导，基于制造解、守恒和标准算例验收。

### 7.3 高阶扩展边界

物理方程由 `EquationSystem` 表示，离散由 `SpatialDiscretization` 表示。首版实现 `CellCenteredFVM`，预留：

```text
WenoFiniteVolume
DiscontinuousGalerkin
```

模型源项通过投影接口进入离散空间：FV 对应 cell source，DG 对应基函数系数。首版只要求 WENO/DG 接口可编译并通过 mock backend 测试，不实现完整 WENO/DG，也不承诺 DG 与 IBM/喷雾耦合。

## 8. 模型扩展接口

### 8.1 三层扩展能力

1. **简单源项 UDF**：显式 `Su`、隐式 `Sp`、物性、入口扰动和诊断；
2. **完整物理模型**：湍流、燃烧、化学、喷雾、蒸发、破碎、边界、IBM forcing；
3. **新方程或离散后端**：作为仓库内编译模块开发。

### 8.2 模型生命周期

统一生命周期为：

```text
configure
declare_fields
initialize
pre_step
add_sources
correct
post_step
checkpoint
restore
```

每个模型只能通过类型化服务访问网格、热化学、时间、统计和源项。模型不得：

- 自行分配或替换核心状态数组；
- 修改 MPI communicator 和 Halo schedule；
- 绕过 `SourceAccumulator` 直接改守恒方程；
- 在未声明的情况下写其他模型字段。

### 8.3 守恒事务

`SourceAccumulator` 为每个守恒方程累计显式和隐式源项，并记录来源。跨相交换使用 `CouplingTransaction`：气相增加的质量、动量和能量必须与颗粒相减少量成对提交。调试模式下可逐 cell 和全局审计不平衡量。

首批 SDK 示例应包括：

- 自定义标量源项；
- 自定义 LES 模型；
- 自定义液滴蒸发关联式。

## 9. TPDF 与 TCR 的统一结构

### 9.1 不拆分为松散插件

TCR 是毕业设计关键模型，并与随机场 PDF 反馈形成递归结构。首版采用紧密组合的：

```text
TpdfTcrCombustionModel
```

其内部拥有并协调：

```text
MeanStateTransport
StochasticFieldSet
StochasticTransport
ReactionEnsemble
PsrShadowReaction
TcrKappaSolver
MicromixingOperator
PdfStatistics
```

这些部件内部可单测和替换，但不能由应用层任意拼接成物理上不完整的调用顺序。

### 9.2 状态命名

禁止继续使用含义不清的 `field0` 表达。统一使用：

- `MeanState`：求解守恒输运方程的滤波平均状态；
- `StochasticFieldSet[N]`：N 个随机 PDF realization；
- `PsrShadowState`：用于 PSR/mean-field 参考反应率的非输运影子状态。

`PsrShadowState` 不参与对流、扩散、Halo、IBM 标量输运或压力求解，只用于与随机场相同化学后端的反应率诊断。

### 9.3 每步数据流

实现必须从公开 TPDF/TCR 方程重新推导。以下列表规定不可缺失的数据依赖，不预先指定唯一的算子分裂或子步先后；阶段 5 开始前必须另写方程与时间离散设计，并用 IEM 极限、守恒和时间步收敛选择最终顺序：

1. `MeanState` 进行守恒输运预测；
2. 随机场执行理论规定的随机输运增量、小尺度混合和有限速率化学；
3. `PsrShadowState` 由当前平均状态构造，使用同一化学后端计算 PSR 参考反应率；
4. `ReactionEnsemble` 汇总 PDF 反应率；
5. `TcrKappaSolver` 由 PDF/PSR reaction-rate ratio 求 `kappa_eff`；
6. `MicromixingOperator` 使用 `kappa_eff` 修正混合时间尺度；
7. 随机场统计量和平均态执行守恒一致性校正；
8. 可选耦合迭代在同一时间步内重复反应率、`kappa_eff` 和混合更新。

不得把随机场化学结果直接替代平均态守恒输运，也不得把 `PsrShadowState` 当成被输运的额外随机场。

### 9.4 IEM 与 TCR 模式

IEM 不是独立、重复实现的混合器，而是 TCR 参数极限：

```text
kappa_mode = unity_iem
kappa_eff  = 1
beta       = beta_iem
```

CD-PhySO-TCR 主路径为：

```text
kappa_mode = cdphyso
rate_ratio = abs(rdot_pdf) / max(abs(rdot_psr), eps_rate)
beta_tcr   = beta_iem * pow(kappa_eff, 1.0 / 3.0)
```

保留 `damkohler_debug` 作为消融模式，但不得作为主物理结论。`kappa_eff` 求解至少区分 inactive、direct root、eta projection 和 similar-state fallback，并输出状态统计。

### 9.5 TPDF 验收重点

- IEM 极限能严格恢复 `kappa_eff=1`；
- PDF/PSR 使用同一化学和热力学路径；
- 随机场数、均值和方差在 Restart 前后连续；
- 随机序列对 MPI 分解不敏感；
- `kappa_eff`、reaction ratio 和状态分类有独立参考测试；
- mean/PDF 质量分数保持正性、归一性和守恒；
- nf4 CD-PhySO-TCR 是正式综合验证配置。

## 10. 化学反应后端

提供稳定的 `ChemistryBackend` 接口，至少支持：

- 物种和热力学查询；
- 单 cell/批量反应积分；
- 反应率和热释放率；
- Jacobian 或近似 Jacobian 能力查询；
- 初猜、缓存和统计；
- 失败恢复和物理范围检查。

首版可使用 Cantera C++ 作为许可证兼容的后端。主程序运行不依赖 Python；发布构建可静态链接 Cantera C++，或通过包内相对 RPATH 提供动态库。Cantera 的 Python/SCons 构建需求不得传播为 HUNDUN-FLOW 的运行依赖。

为性能和长期独立性，接口应允许未来增加原生生成 kinetics 与 SUNDIALS/CVODE 后端。任何后端都必须通过同一 0D reactor、PSR、点火延迟、热释放和批量积分合同。

## 11. IBM 与网格

### 11.1 首版 IBM

IBM+LES 是首版基线。首版只实现和验证：

- 单个静态 STL；
- 固定结构化或曲线坐标背景网格；
- solid/fluid marker、法向、壁面距离或 SDF；
- 静态无滑移/热边界；
- 组分与 TPDF 的静态壁面处理；
- 液滴与静态壁面的碰撞接口。

首版明确不实现和不验证：

- 多 STL 和多部件组合；
- 运动或变形壁面；
- 新覆盖/新暴露单元；
- GCL；
- 活塞、阀门和曲轴运动；
- 动态壁面粒子碰撞。

不支持的 JSON 配置必须在预检阶段报错。

### 11.2 未来接口

内部 API 预留：

```text
GeometryScene
GeometryPart
RigidTransform
WallKinematics
MotionLaw
```

几何查询返回 `part_id` 和 `wall_velocity`；首版始终只有一个 part，且 wall velocity 为零。提供 mock `MotionLaw` 编译与注入测试，以证明接口可扩展，但不声称支持运动壁面。

## 12. 喷雾、蒸发和两相耦合

毕业交付必须包含可运行的 Eulerian-Lagrangian 喷雾能力：

- parcel 注入；
- 阻力和速度松弛；
- 液滴加热；
- d-squared 基线蒸发；
- 物种、质量、动量和能量双向耦合；
- MPI 粒子迁移；
- Restart；
- 静态 IBM 壁面碰撞；
- 破碎模型接口和至少一个经过验证的基线实现。

喷雾模型通过 `CouplingTransaction` 提交源项。任何粒子丢失、跨 rank 重复、气液质量不闭合或 Restart 数目跳变都必须触发诊断或失败。

多组分蒸发、稠密喷雾、液膜、运动壁面和高级破碎模型属于后续扩展，不应阻塞首个气相燃烧版本。

## 13. 配置、I/O 与可复现性

HUNDUN-FLOW 使用一个类型化 JSON case 配置，不支持旧 `.d` 文件自动回退。配置解析应在 rank 0 完成，执行 schema、范围、路径和交叉字段校验后广播类型化配置。

配置至少包含：

- 网格和 MPI 预期拓扑；
- 物理方程和时间推进；
- 边界与初始条件；
- LES、化学、TPDF-TCR 和喷雾模型；
- IBM 几何；
- Restart 和输出；
- 用户插件和参数。

未知键默认报错。错误信息包含 JSON Pointer、期望类型和实际值。配置和输出路径不得逃逸 case 根目录，除非维护者显式启用受信任的开发模式。

标准命令建议为：

```bash
mpirun -np N ./hundun case.json
./hundun case.json --validate
./hundun case.json --print-resolved
./hundun --version
```

Restart 必须保存：

- 全部守恒状态；
- TPDF 随机场和随机序列继续所需状态；
- 粒子；
- 时间、步数和自适应时间步状态；
- 模型 checkpoint；
- 网格、机制、配置和二进制版本元数据。

首版可使用自有二进制 Restart 和 VTK 输出。HDF5 是可选依赖，不得成为最小构建的强制条件。

## 14. 依赖与构建

### 14.1 强制依赖

- C++17 编译器；
- CMake；
- MPI-3；
- POSIX/Linux 和 C++ 标准库；
- 一个体积小、许可证兼容、随源码固定版本的 JSON parser。

### 14.2 可选依赖

- Cantera C++；
- SUNDIALS；
- HDF5；
- OpenMP。

### 14.3 不进入基础依赖

- OpenFOAM runtime；
- PETSc/HYPRE；
- Boost；
- Kokkos；
- Python runtime；
- Conda；
- 构建时在线拉取依赖；
- GPU 平台。

建议构建目标：

```text
hundun
hundun_mesh_tool
hundun_case_validate
hundun_plugin_sdk
unit_tests
mpi_tests
```

首发平台为 Linux CPU。旧 glibc 兼容通过目标环境构建实现，不捆绑 glibc。

## 15. 验证体系

### 15.1 单元测试

- 曲线坐标度量和几何一致性；
- limiter、重构、通量和扩散；
- Rhie-Chow 和压力校正；
- CG/BiCGStab；
- 热力学和反应后端；
- TCR 四类 `kappa_eff` 状态；
- 分区无关随机数；
- 粒子阻力、传热、蒸发和破碎；
- 喷雾交换守恒；
- IBM marker、法向和壁面距离。

### 15.2 数值验证

- 制造解二阶收敛；
- 全局质量、动量、能量和组分守恒；
- 正性和物种归一；
- 时间步收敛；
- MPI 分解不变性；
- Restart 连续性。

### 15.3 模块基准

流动与 LES：

- Taylor-Green vortex；
- channel flow；
- variable-density vortex；
- pressure checkerboard suppression。

IBM：

- IBM channel；
- cylinder；
- sphere。

化学：

- 0D reactor；
- PSR；
- ignition delay；
- chemistry backend comparison。

TPDF-TCR：

- IEM limit；
- nf4 stochastic ensemble；
- direct root、projection、fallback、inactive；
- PDF/PSR reaction-rate ratio；
- 随机数与 Restart 连续性。

喷雾：

- Stokes relaxation；
- d-squared evaporation；
- droplet heat transfer；
- MPI particle migration；
- two-way coupling conservation；
- static IBM collision。

### 15.4 综合验证

综合算例按以下顺序推进：

1. Flame D；
2. GTMC 双旋流；
3. Vblowoff/TCR 熄火；
4. 两相燃烧室。

HUNDUN-FLOW 与 COAST 的比较只用于黑盒科学合理性：比较守恒量、均值、趋势、火焰位置、热释放和统计范围，不要求 bitwise 或逐步一致。

### 15.5 工程质量

- GCC 和 Clang；
- Debug 和 Release；
- ASan/UBSan；
- MPI 正常和异常退出；
- `ldd` 确认无 Python runtime；
- 不出现 BOFFIN 符号、ABI 和 legacy 文件名；
- 小型测试数据进入 Git，大型 Restart 和输出外置。

## 16. 分阶段开发与出口条件

开发按阶段推进，不按周排期。每一阶段只有达到出口条件后才进入下一阶段。

### 阶段 0：版权和科学规格冻结

工作：

- 建立隔离的新 Git 仓库；
- 固定 BOFFIN 法律基线；
- 保存方程、算法和参考文献清单；
- 建立私有相似性审计；
- 固化许可证、DCO、NOTICE 和第三方规则。

出口：公共仓库不含前置代码；基础版权扫描和人工审查通过。

### 阶段 1：C++ 运行时

工作：

- MPI lifecycle；
- 域分解；
- FieldRegistry/FieldView；
- Halo exchange；
- 结构化网格和基础 I/O；
- JSON 配置；
- 插件 mock；
- 被动标量求解。

出口：多 rank 被动标量守恒、二阶收敛、分解不变和 Restart 连续。

### 阶段 2：变密度流动

工作：

- 质量、动量、压力、焓和标量；
- PISO/Rhie-Chow；
- 曲线坐标度量；
- 自适应时间步；
- 线性求解器。

出口：标准流动基准、守恒和无 checkerboard 通过。

### 阶段 3：静态 IBM 与 LES

工作：

- 单静态 STL；
- IBM marker、距离、法向和 forcing；
- 至少一个 LES 基线模型；
- IBM+LES 组合验证；
- 未来 motion API mock。

出口：IBM channel/cylinder/sphere 和 LES 基准通过；不支持配置能明确失败。

### 阶段 4：反应流与化学

工作：

- ChemistryBackend；
- 热力学、物种、焓和反应源项；
- 0D/PSR；
- Strang splitting；
- Flame D 气相验证。

出口：反应器基准和 Flame D 无 NaN、守恒、趋势合理。

### 阶段 5：TPDF 与绑定 TCR

工作：

- 随机场集合；
- 随机输运和微混合；
- IEM 极限；
- PDF/PSR 反应率；
- CD-PhySO `kappa_eff`；
- nf4 Restart 和诊断。

出口：TPDF-TCR 单元/数值合同通过，Vblowoff 熄火趋势具备可解释性。

### 阶段 6：喷雾与蒸发

工作：

- parcel、迁移、阻力、传热、蒸发；
- 双向守恒；
- Restart；
- 静态 IBM 碰撞；
- 基线破碎模型。

出口：单液滴、迁移、守恒和两相小算例通过。

### 阶段 7：综合功能和科学等价性

工作：

- Flame D、GTMC、Vblowoff 和两相燃烧室；
- 性能剖析；
- 黑盒科学参照；
- 文档和示例模型。

出口：四个综合算例可复现运行，主要物理趋势、守恒和稳定性达到研究使用要求。

### 阶段 8：GitHub 发布

工作：

- 公共文档；
- API/SDK 示例；
- CI；
- 版本和发布资产；
- 最终 BOFFIN 独立性审计。

出口：源码可复现构建，核心基准通过，四个综合算例有公开可运行配置，版权审计通过。

## 17. 首版非目标

- 超音速、激波和全可压缩燃烧；
- 完整 WENO 或 DG；
- 多 STL、多部件 IBM；
- 运动/变形壁面和活塞发动机；
- GPU；
- 液膜和稠密喷雾；
- OpenFOAM/PETSc 运行时集成；
- 对 BOFFIN/COAST 数值结果逐步复制；
- 将 COAST/COAST-2 转换成新仓库起点。

这些非目标不取消未来接口，但不允许以“预留未来”为理由扩大首版实现范围。

## 18. 交给下一位 Agent 的执行规则

下一位 agent 开始工作前必须：

1. 阅读本文件；
2. 确认新仓库与 Coast_software 物理隔离；
3. 将 `/home/wyf/Compress_boffin/AECSC_WDQ/SRC-Spray` 登记为只读审计基线；
4. 不打开旧 Fortran 过程来编写对应新过程；
5. 优先根据公开方程、测试和接口契约实现；
6. 对候选复用 C++ 文件逐个完成私有审计；
7. 每个阶段先写可失败的测试和出口验收；
8. 未经用户确认不得扩大首版范围。

代码查看优先使用本地 `codegraphf`。对现有 COAST 的查看只服务于资产归属审计、功能清单和黑盒输入输出识别，不用于复制控制流。

建议新 agent 首个实施计划只覆盖阶段 0 和阶段 1，不同时启动流动、IBM、化学、TPDF 和喷雾开发。阶段 1 完成前，禁止搬运综合燃烧算例或建立与 COAST 相似的主循环。

## 19. 最终毕业完成定义

同时满足以下条件才可称为 HUNDUN-FLOW 毕业版本：

- 源码版权和依赖许可证清晰，可在 GitHub 公开；
- 标准环境可从源码复现构建，无 Python 运行依赖；
- MPI、变密度流、IBM+LES、反应流、TPDF-TCR 和喷雾均有独立验证；
- Flame D、GTMC、Vblowoff 和两相燃烧室能够运行；
- 用户可通过 SDK 示例增加源项、LES 或蒸发模型；
- BOFFIN 独立性私有审计通过；
- 公共仓库不存在旧 Fortran ABI、legacy adapter 或来源不清代码；
- 文档诚实区分已实现功能、实验功能和仅预留接口。

HUNDUN-FLOW 的成功标准不是拥有最多模型，而是形成一个版权清晰、物理守恒、可验证、可扩展，并能让后续研究人员可靠加入新模型的燃烧模拟基础设施。

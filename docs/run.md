# 启动与运行控制

在算例目录放置 `case.json`、物性资产与 `run.json`，执行
`mpirun -np 4 /path/to/hundun`。`hundun init-case --output c` 生成示例。

```json
{"mode":"new","steps":100,"output":"../out","monitor_interval":1,
 "output_interval":20,"restart_interval":50,"diagnostics_interval":10}
```

`mode` 为 `new` 或 `restart`。恢复时添加 `"restart":"../out/Restart"`。
`steps` 表示本次新增步数；`end_time` 表示第一个达到该物理时间的接受步，
两者同时存在时取先达到的条件。时间步继续遵循算例时间控制。
`initial_state` 为 `[p,T,Ux,Uy,Uz,q...]`，用于指定新算均匀场。
所有相对运行路径以 `run.json` 所在目录为基准，输出目录与输入目录采用同级布局。

显式启动保持 `hundun run CASE --output OUT --steps N`。
`hundun run CASE --config FILE --steps N --output OUT` 先读运行文件，
再按具名键应用命令行覆盖。根进程解析后广播同一组最终参数。
启动信息记录最终步数、周期、方法指纹、进程数与输入／输出／恢复路径。

- `--monitor-interval N`：屏幕摘要和 `monitor.jsonl`，默认每接受步。
- `--output-interval N`：三维场输出。
- `--restart-interval N`：原生检查点。
- `--diagnostics-interval N`：详细守恒账本。
- `--until T`：显式命令的物理终止时间，配合 `--steps` 使用。

各周期独立配置，零值表示关闭该类定期输出。接受步摘要含 dt、
CFL 数值及定义、方程残差、迭代次数和最大进程推进壁钟时间。
`cn_phases.seconds` 保存 CN/BE 的准备、标量、物性、动量组装／求解、
压力组装／预处理／求解、校正及审核耗时。各阶段累计全部候选与重试，
包含阶段内部通信；每阶段取最大进程值，其所属进程可以随阶段变化。
整步时间使用 `seconds`，BE/PISO 的 CN 阶段项为零。

`physics_modules.seconds` 提供四类完整调用耗时：

| 项目 | 计时范围 |
|---|---|
| `reaction_sources` | 平均场有限速率／PaSR 的本步源准备，含物性查询与源写入 |
| `mean_reaction` | 输运后平均场化学区间推进，含状态检查、响应缓存及源更新 |
| `esf_reaction` | 随机场化学推进，含 PSR 参考、热状态重构与共同源更新 |
| `tcr_statistics` | 动态 TCR 终态统计调用，含调用内部的通信与历史准备 |

各项累计本步全部候选、重试及失败调用，并按最大进程壁钟输出。
计时位于 `advance` 内，部分调用同时属于 `cn_phases`；初始化独立于
本步计时，关闭的模型对应零值。ESF 的输运和隐式混合归入其实际执行
阶段，`esf_reaction` 聚焦化学与其直接状态处理。计时属于运行观测，
物理状态与 Restart 身份保持原定义。

`communication_observations.seconds` 分别提供 `structured_wait`
（已登记结构化 halo 的 MPI_Waitall）、`structured_control`
（这些 halo 的集体状态检查）、`linear_reductions`（公共线性归约对象
登记的通信）。这些观测覆盖各自具名操作，并与方程／物理模块耗时
重叠。IBM 供体、TCR 内部和直接 MPI 调用的通信仍由所属模块壁钟
承载。各项 scope 与 containment 字段记录统计口径；性能比较以
整步 `seconds` 为总耗时，逐模块观测用于定位成本。

在输出目录创建 `stop` 文件，本步接受后保存 Restart 并停止；
创建 `output` 文件，本步接受后额外保存场数据。请求由根进程接收，
全体进程执行，写入成功后登记 `control.jsonl` 回执并消费请求。
写入期间采用 `.pending` 标记保存请求身份，新到请求留给下次处理。
`hundun status OUT` 查看 `status.json` 的阶段、接受步和更新时间。

状态在每次步内尝试开始时更新。`phase=solving` 表示首次尝试，
`phase=retrying` 表示正在重试；`step/time` 对应已接受状态，
`target_step` 对应当前目标步。`attempt` 和 `coupling_sweep` 从 1
开始，`dt` 表示此次实际时间步；`retry_kind=time_step` 表示时间步
回退，`scalar_coupling` 表示同 dt 的标量耦合重算。
`previous_failure` 保存上次失败的 code、detail、stage、attempt 和 dt。
根进程通过临时文件和原子替换发布状态；写入错误在推进返回后由
应用统一报告，接受／回退始终由求解事务决定。状态发布耗时计入
实际推进壁钟。结束、保存停止和写盘阶段继续使用各自明确状态。

三维场默认采用 legacy 二进制 `.vtk`，在 `Visit/solution.visit` 中按真实
物理时间组织。向量 `Velocity` 与标量 `Temperature`、`Pressure`、
`PressureGauge`、`Density`、`Enthalpy` 使用 SI 单位。Pressure 为
参考压力加 pi，PressureGauge 为 pi。双状态 ESF 的 Density 对应
压力耦合 field0 密度；Temperature 与 Enthalpy 对应物理均值热状态。
双状态 ESF 的场输出同时提供以下诊断目录：

| legacy VTK 名称 | XML 名称 | 定义／单位 |
|---|---|---|
| DensityMeanEOS | rho_mean_eos | 物理均值组分、焓的 EOS 密度，kg/m³ |
| DensityStatistical | rho_pdf_mean | 随机场比容均值的倒数 `1/mean(1/rho_f)`，kg/m³ |
| DensityField0 | rho_field0 | field0 压力闭合密度，kg/m³，与 Density 同一状态 |
| TemperatureField0 | T_field0 | field0 正权 EOS 查询温度，K |
| EnthalpyField0 | h_field0 | 已接受的原始 field0 焓，J/kg |
| MassFractionsField0 | Y_field0 | 原始 field0 组分权重，多分量数组，kg/kg |

组分数组依照 `thermophysics.d` 的完整 species 顺序排列，包含余组分。
field0 正权 EOS 查询采用 `M=sum(max(Y0,0))`、`Y+=max(Y0,0)/M`
和 `h+=h0/M`，压力闭合密度为 `rho(p,h+,Y+)/M`；原始权重与焓
保持其独立状态。固定热力学压力配置使用 p0 查询全部 EOS，机械
压力继续由 Pressure／PressureGauge 表示。物理均值组分沿用各组分
具名标量及组成闭合定义。

上述诊断从接受态重建，按三维场周期输出；工作区和写盘容量进入
编译资源预算。固体格的诊断值采用零占位，按算例 IBM 几何选择流体区域。
原有变量目录继续通过 `--visit-format xml` 或运行文件
`"visit_format":"xml"` 选择。

采样点为实际单元中心，正方向分区交界增加相邻中心点，并以
`vtkGhostType` 标记重复点。格子内部、交界面、棱和角的值由输出专用
交换取得。时间索引在各分区数据完成后原子发布；每个输出目录保存
固定进程数的序列，切换进程数时选用新的输出目录。

`hundun restart-info OUT/Restart` 报告版本、步号、时间、dt、源进程数、
网格、方法／模型历史身份和文件校验结果。默认逐文件流式核对校验和，
工作缓冲区为 64 KiB；`--metadata-only` 读取并核对清单。
具体算例的模型兼容性由原生恢复入口结合 case.json 检查。

PDF 当前态迁移已接入正式程序：

```sh
mpirun -np 4 ./hundun import TRANSFER --format pdf-transfer-v1 --case CASE --output SEED
mpirun -np 4 ./hundun run CASE --restart SEED --output RUN --steps 1
```

TRANSFER 使用 `state.txt` 的 `HUNDUN_PDF_TRANSFER 1` 头，依次记录
`nx ny nz step time dt pressure_reference field_count species_count`，
再列出完整物种名称。物种顺序对应算例热物性资产；独立组分通过名称
映射，余组分使用算例定义。
`flow.f64` 每格保存 `u v w p_mechanical`，`rho_ref.f64` 保存参考
密度；`pdf0.f64` 等文件每格保存完整 Y 和 h；`fluid.u8` 使用
0／1 固体／流体标记。浮点载荷为小端 FP64，单元按 x 最快排列。

机械压力定义为 `pressure_reference+pi`。固定热力学压力配置下，
机械压力按源值保存，EOS 和 `rho*(h+K)-p_eos` 迁移账本共同使用
p0；耦合 EOS 配置使用正的绝对压力。固体格重建静止 295 K 空气
占位状态，转换器要求物种目录包含 O2 和 N2。

CASE 为配置完备的原生算例，TRANSFER 为已映射到目标网格的传输数据。
SEED 使用独立的新目录；路径准入保护源数据、算例与已有检查点。
新入口完整读回和原生初始化通过后，将报告写入 `SEED/import.json`。
独立工具 `v04_pdf_import CASE TRANSFER SEED` 共用同一实现，并延续
`TRANSFER/native.json` 的兼容报告路径。

转换报告记录压力分工、完整／独立物种顺序、
源步号及时间、物性重建前后库存。该入口从随机场物理均值重建 field0，
以 V1 当前态检查点登记方法历史恢复；原版完整统计历史、守恒网格
映射和源算例自动绑定按专项后续条目接续。首步推进还需
源场与目标边界、时间步及物理模型匹配。

携带被动标量时使用 `--format pdf-transfer-v2`：`state.txt` 的版本
设为 2，在 `field_count species_count` 后增加 `passive_count`，物种
名称之后列出被动标量名称。`passive.f64` 按单元保存这组标量，采用
小端 FP64、x 最快、分量交错的布局。源名称顺序可以与原生配置不同，
转换按名称映射并逐值读回。标量保持源值，支持有符号示踪量。
报告中的 `passive_scalar_order` 为原生保存顺序。

V2 扩展的是传输数据目录，输出继续采用 V1 当前态 Restart 和明确的
方法历史恢复。动态 TCR 所需的统计历史迁移单独接线；当前动
态模型配合这一当前态桥接器时，在恢复阶段返回 `invalid_plan/10217`。
完整动态模型的原生 Restart 继续使用已实现的模型历史格式。

Restart V6 保存当前流场、面通量和固定长度模型记录，时间历史明确
采用首步 BE 恢复。它用于同时满足“保留模型统计”与“重建时间历史”
的迁移状态；统计时钟、分支及累计量通过模型自身的类型检查。
V4／V5 继续承载完整时间历史，V1 保存普通当前态。V6 的固定记录
沿通用 Restart 分区转换，颗粒等变长记录采用其对应格式。

气相时间配置支持 `time.convective_cfl_definition`：`outgoing_sum` 为
单元向外面质量通量之和，`directional_max` 为六个面质量通量绝对值
的最大值；两者均乘以 `dt/(rho*V)`。原版气相对齐输入显式选用
`directional_max`，既有输入延续 `outgoing_sum`。目标及浮动范围由
`convective_cfl`、`convective_cfl_margin` 控制，推荐 0.30 和 0.05。

`hundun check`、接受步摘要和运行证据报告实际定义。监看同时记录
向外通量和、绝对通量和及方向最大值；准入和自适应 dt 使用所选指标。
固定 dt 超过所选上限时返回定位信息；自适应尝试超过上限时沿统一
回退／重试流程处理。Restart 的方法历史签名携带 CFL 定义，原生
连续恢复沿用来源定义。

气相参考对照可设置 `solver.reference_outer_iterations`。取值 1–64
表示每次候选尝试的固定外迭代次数；0 或省略表示按原方程残差停止。
该选项适用于 `cn_be`／`outer_corrected`，与完整压力求解配置放在同一
`solver` 对象内，可同时配置 `cold_stopping` 的参考尺度和阈值。
每次外迭代共享本步接受历史，末次执行完整方程及热力学审计，达到验收
阈值后提交本步；残差超限沿统一失败／回退流程处理。
`hundun check` 与运行证据分别记录请求次数和实际次数，Restart 的
算例与方法身份包含该配置。该模式用于同次数算法对照，日常运行继续
采用残差停止模式。


`cn_be` 的动量采用 CN；普通单流体焓采用守恒中点离散。配置组分
输运、反应、ESF 或喷雾时，焓与组分共同采用 BE，保持组成输运与
形成焓的一致时间层。`backward_euler` 的各方程采用 BE。
`hundun check`、接受步摘要及监看中的 `enthalpy_scheme` 报告实际
热方程策略，Restart 方法签名记录相应时间历史。
旧版普通焓检查点可通过显式 `--restart-method-recovery` 重建方法
历史；迁移运行另保存新输出目录及来源说明。新版本原生续算沿用
检查点所登记的方法身份。


压力参考后端可在 `solver.pressure_linear` 中设置 `"algorithm":"pcg"`
及 `"krylov_restart":0`，配合 `cn_be`／`outer_corrected`。程序使用
体积缩放、不完全 Cholesky 和 PCG，实际矩阵在每次求解前通过面系数
精确对称、弱对角占优及全域连通分量约束检查。该充分条件保留原矩阵
系数；一般压力系统使用 `fgmres` 或 `bicgstab`。

ICCG 的缩放行 L2 阈值由原始配置与最小单元体积换算，候选解同时
通过原单位 L2 和连续性审核。监看及证据中的 `pressure_iccg` 表示
实际后端组合，`pressure_original_l2` 与 `pressure_original_l2_limit`
保留原方程残差及门槛。`hundun check` 报告后端与预分配工作区字节数。
相同配置可使用原生 Restart 在不同进程数之间续算。

`cn_be`／`outer_corrected` 支持 `transported_scalars` 中的
`passive_scalar`。被动标量采用 CN 空间中点与完整 rho*q 守恒存储，
使用本步终态质量通量、公共边界矩阵消元和 IBM 固体行。中心格式直接
使用中点值；限幅格式从当前端点计算 VLS 面系数，再作用于中点。
校正矩阵包含空间项的 1/2 响应，每次候选共享已接受的标量历史。
有符号示踪量和多个被动标量沿用同一接口。`backward_euler` 入口继续
使用其 BE 标量调度。`hundun check`、运行证据和监看的 `passive_scheme`
记录当前 CN/BE 入口的被动标量时间身份。

CN 接线使用独立的 Restart 方法标记。此前 BE 被动标量检查点通过
`--restart-method-recovery` 显式重建方法历史后进入 CN 计算；新写出的
检查点可直接跨进程数续算。方法恢复的初始化身份与来源随运行记录保存。

运行证据的 `cold` 与监看的 `payload` 中，`passive_scalar_count`、
`passive_solve_calls`、`passive_iterations` 记录实际工作量；
`passive_residual` 为局部时间尺度
归一化的最大原方程残差，`passive_balance_defect` 为体积积分收支的
归一化缺陷。两者提交门槛均为 128 倍 FP64 epsilon。
`hundun check` 的 `passive_workspace_bytes` 报告各标量复用的预分配
校正工作区，输运耗时计入 CN 的 scalar 阶段。


### 近零物种／元素收支

`diagnostics.jsonl` 的 `composition_balance.roundoff_rule` 为
`fp64-local-storage-v1`。每项收支保留原 `defect` 与 `relative_defect`，
另给出同速率单位的 `storage_roundoff_bound` 和 `roundoff_applied`。
审核使用 `relative_defect < 1e-6`，或在原相对门槛之外满足
`abs(defect) <= storage_roundoff_bound`；后者将 `roundoff_applied` 设为 true。

2026-09-17 用户批准此近零收支定义。每个流体单元按实际 V、ρ、
随机场平均 Y 传播输入分辨率 `u(x)=epsilon_FP64*abs(x)+denorm_min`。
VρY 的包络用乘积展开计算，避免大数相减；前后两层包络相加后除以
本步 dt，再按 MPI 求和。随机场物理 Y 为非负值，其平均输入包络
同样为 `epsilon_FP64*mean(Y)+denorm_min`。总质量的组成 1 作为精确
常数；余组分采用总质量和独立组分包络之和；元素按原子数／分子量
加权传播。界限仅用于物种／元素的存储分辨率审核。

质量与总能量时间项在单元内形成增量后汇总，保留原始物理存量的
输出格式。小变化率使用 `mass_bdf_rate_kg_s` 和
`total_energy_bdf_rate_W`；直接相减已输出的绝对存量会再次引入
FP64 抵消误差。

### COAST 输入参数目录

`tools/input.py` 按已登记的 `input.F90` SHA-256 选择 dyn711、GTMC
或 624CF 的 `input.d` 行布局。用法如下：

```sh
python3 tools/input.py source/input.d --reader source/input.F90 --output check/input.json
```

输出含每项参数的原始行、行号、单位、读取值、尾随值及读取器默认值，
同时记录输入与读取器的路径和 SHA-256。数值表示源文件的十进制值；
源程序的 REAL 精度由配套构建清单另行确定。`steps` 保存源 `lstep`；
实际终止步数由原版启动／续算语义与源检查点共同确定。

扩展参数与 namelist 以原始记录保存，模型接线阶段按对应扩展读取器
解析。当前标量语法支持逗号、D 指数、逻辑值、引号字符串及行注释；
空值和重复语法先展开为明确的标量值，namelist 名称单独占一行。

目录标记为 `source_catalog_ready_native_case_binding_pending`。
压力模型和 CFL 映射单独列出；非对称 CFL 窗口附带策略选择状态。
完整原生算例还需绑定机理、几何、边界、同步 Restart／PDF 及历史，
并应用已批准的模型选择。例如 624CF 的源 SGS 名称作为来源保留，
Hundun 的 Vreman 选择在原生模型配置中登记。

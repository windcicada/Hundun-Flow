# 气相专项执行台账

入口为 [gas.md](gas.md)。当前实施承接 main eba8a1b 的参考冻结与运行接口节点。

| 项目 | 本次实现与证据 | 接续工作 |
|---|---|---|
| G0 | gas-ref.json：131 文件哈希、关键调用表；gas-cfl.json：完整 courant 例程的 259 输入对照 | 实际参考构建、算例后端与端到端计时窗口 |
| G1 | CFL 定义贯通配置、准入、自适应 dt、运行证据与 Restart 方法身份；固定外迭代参考模式贯通末次完整审计、方法身份和证据；定向检查见文末 | 普通焓 CN 已接线；完整 jstep 调度及非均匀场自适应重试 |
| G2 | Vreman、ICCG 原生接线及 CN/BE 被动标量 CN 守恒输运，证据见文末 | 矩阵／通量、壁函数、精度及大规模同场性能；近静态精度与效率；已批准舍入界及 300 K 接线见文末 |
| G3 | dyn711 原生 CN/BE 调度、单区间守恒重新混合、Bilger CN、终态局部统计及 4/8/16 场跨分区恢复；详见文末 | 热释放触发的重新混合轨迹、完整时间阶及实场动态模型组合 |
| G4 | 继承 JL4 原时间步组合检查 | 化学筛选／任务均衡、热源／边界及实场轨迹 |
| I1 | 独立监看、CN／化学／TCR 计时及通信包含关系、接受点 stop/output、步内重试状态；2→4 进程恢复通过 | 完整运行写盘成本与接口交付说明 |
| I2 | 当前目录启动、run.json、显式覆盖、步数／时间结束，严格模式检查通过 | 简洁物理输入模板及最终方法说明扩展 |
| I3 | legacy FP64 VTK、中心坐标、分区正向重叠、原子时间索引；3/8 进程 VTK 读取与流线、ESF 双状态字段目录通过 | VisIt 应用打开和迁移序列关联 |
| I4 | Restart 查看、流式校验和、损坏检测；正式 import 的 PDF 具名物种／固定 p0 迁移与 2→4 进程首步恢复；三个源版本的 input.d 参数目录及行级来源 | 源模型／边界自动绑定、扩展读取语义与统计历史导入、实场守恒映射及转换格式扩展 |

本台账逐项保留专项出口。各条证据按实际覆盖范围登记，完整专项验收
由 G0–G4 和 I1–I4 的后续条目共同闭合。

CFL 复现命令：

```sh
mkdir -p check/gc/f4 check/gc/f8
cp /home/wyf/code_dev/src.TCR.dyn711/courant.F90 check/gc/courant.F90
gfortran -O2 -Jcheck/gc/f4 tools/gas_cfl.f90 check/gc/courant.F90 -o check/gc/f4/probe
gfortran -O2 -fdefault-real-8 -Jcheck/gc/f8 tools/gas_cfl.f90 check/gc/courant.F90 -o check/gc/f8/probe
python3 tools/gas_cfl.py --float check/gc/f4/probe --double check/gc/f8/probe \
  --hundun b3/versions/v0.4/tests/v04_gas_cfl_probe --runner 'bash check/jam.sh' \
  --source check/gc/courant.F90 --output docs/gas-cfl.json
```

外围串行驱动把全局最大操作限定为单进程恒等操作，保持原例程的全部
局部循环和公式。默认 REAL 与 FP64 的比较分别记录，程序文件哈希
随结果保存。

接口定向检查：`check/gas-io.log`（5/5）、`check/gas-timing.log`
（CN/BE 监看）、`check/gas-fault.log`（输出路径与 MPI 故障收敛 2/2）。
`check/gas-vtk.log` 与独立 VTK 读取覆盖新版索引顺序保护；日志与
测试目录为本地生成资产，复现入口随源码保存。

G1 CFL：`check/gas-cfl-final.log` 通过原生定向检查；
`check/gas-cfl-check.log` 的输入、MPI 广播、时间控制与 PISO 四项回归
通过。外部证据校验器 self-test 返回成功。原版完整 courant 例程的
259 输入比较沿用 gas-cfl.json；新增接线保持局部公式与 FP64 运算顺序。

G1 固定外迭代：`check/gas-outer-final.log` 原生入口检查通过。覆盖 1／3
次完整外迭代、末次原方程残差审计、入口温差场单次残差超限、2→4
进程恢复、配置身份漂移、输入范围及外部证据篡改。`gas-outer-check.log`
的输入与 MPI 广播两项回归通过；证据校验器 self-test 返回成功。
正式构建完成后执行上述检查；原生残差停止模式与旧输入身份保持原行为。

G0/G1 时间算子基准：`gas-time.json` 登记完整 `cmod.F90`／`step.F90`
的 257 组输入、三条路径及源码／程序哈希。CN 与 Hundun 现用
`time_centre_cold_row` 逐系数、RHS 比较；BE 与压力行分别与解析式比较。
FP64 三条路径本组差异均为 0，原版默认 REAL 最大归一化差异为
1.65e-7。RHS 按输入项幅度归一化，避免抵消后的近零分母放大舍入。
此证据覆盖冻结行算子，普通焓完整接线及时间阶由 G1/G2 继续验收。

已定位的普通焓接线范围：原版 `boffin` 的普通焓调用
`condif → source → bndry2 → cmod → step(rhobar)`，其中 `source` 的
`dpdt` 保持全量；PDF 焓调用 `source_pdf → step(rho)`，两次 `jstep`
共享 `fold`。接线前 Hundun 高效能量组装与末次审计均使用 BE。
本次接线范围同时覆盖空间矩阵／RHS、密度时间层、压力功、IBM、边界
消元、总能量收支及 Restart 身份。既有 CN 动量行算子已具备原版对照。

```sh
mkdir -p check/gt/f4 check/gt/f8
cp /home/wyf/code_dev/src.TCR.dyn711/cmod.F90 check/gt/cmod.F90
cp /home/wyf/code_dev/src.TCR.dyn711/step.F90 check/gt/step.F90
gfortran -O2 -Jcheck/gt/f4 tools/gas_time.f90 check/gt/cmod.F90 check/gt/step.F90 -o check/gt/f4/probe
gfortran -O2 -fdefault-real-8 -Jcheck/gt/f8 tools/gas_time.f90 check/gt/cmod.F90 check/gt/step.F90 -o check/gt/f8/probe
python3 tools/gas_time.py --float check/gt/f4/probe --double check/gt/f8/probe \
  --hundun b3/versions/v0.4/tests/v04_gas_time_probe --runner 'bash check/jam.sh' \
  --cmod check/gt/cmod.F90 --step check/gt/step.F90 --output docs/gas-time.json
```

G1 普通焓 CN 公共算子：`EnthalpyMidpointView` 显式承载已接受／端点
热变量的平均值。完整组装、独立残差、热扩散面系数、焓校正矩阵和
全域能量账本共用该空间状态；`rho*h` 存储、压力功、动能项与注册
显式热源保持各自原定义。连续性约化使用
`R_h - h_mid*C`，其瞬态系数为 `rho_mid/dt`，端点响应为 1/2。
PDF／既有端点调用继续通过原接口默认值选择 BE。

`check/gas-mid-final.log` 定向检查通过：普通导热和 unity-Lewis 焓扩散、
中点对流、压力功、显式热源、变密度存储、半权扩散矩阵、连续性约化、
实体外边界的全域能量账本，以及完整／独立残差逐位一致。局部多项式
检查与实体边界的全域检查分别构造有效输入。原有 BE、解析阶数与
热边界检查随同通过。

本节点完成公共方程组件。生产调度仍需预分配中点工作区、刷新已接受
及试探状态的边界／halo、接入 IBM 外部热修正、冻结实际逐方程策略，
并登记配置与 Restart 方法身份；随后执行原生入口及时间阶检查。

正式程序重建成功；`check/gas-mid-app.log` 原生固定外迭代／残差拒绝／
2→4 进程恢复检查通过。该检查确认当前生产入口沿用端点热策略，
中点热算子的生产启用作为下一接线节点单独验收。


G1 普通焓 CN 生产接线：普通单流体的 `cn_be` 使用焓／温度中点，
组分输运、反应、ESF 与喷雾使用 BE 焓。原始 `input.F90:172` 对 PDF
路径设置 `icycl(nvh)=0`；单场组成／反应同属 PDF 输运分支，焓与组分
采用共同时间层以维持形成焓参考协变性。显式 `backward_euler` 保持 BE。

中点复用已预分配的候选工作区；接受场与端点的边界／halo 进入同一
空间状态。组装、独立末次审计、IBM 绝热热修正与全域能量账本共同
使用该状态。实际策略写入 check、监看、证据及 Restart 方法签名。
旧方法历史由显式恢复／迁移入口识别，新方法按自身身份原生续算。

原生周期热扩散检查 `v04_product_cn_heat_time` 采用三个 dt，弱温度
扰动的一阶傅里叶幅值分离线性时间误差与热膨胀二阶谐波。观测时间阶
为 2.00034、2.00010；离散 CN 放大因子与实算幅值差异为
3.05e-14、9.36e-15、1.44e-14 K。能量缺陷绝对值范围为
1.70e-11–5.25e-11 W。原方程与守恒阈值沿用原值，线性绝对容差设为
1e-17，以解析弱扰动的动量校正。数值见 `check/gas-heat-order.log`。

`check/gas-mid-wire.log` 的 CN 动量与原生入口两项检查通过；
`check/gas-mid-close.log` 中普通焓时间阶、ESF 热源 MPI2、温差入口及
2→4 进程恢复通过。通用编译计划测试识别出早期快照数量断言：四个
SGS 派生字段已加入当前输出，现按五个原生字段＋四个派生字段检查。

`check/gas-mid-ibm.log` 两项通过：编译计划／方法身份检查，以及
2 进程 IBM 温差质量入口检查。IBM 覆盖三个方向、内流／外流两类
流体侧，原流固零通量、固体速度、EOS、连续性和能量阈值保持。
外部证据校验器 self-test 返回成功；正式程序完成重建。


G2 ICCG 前置：PCG 增加与 FGMRES／BiCGStab 共用的原方程收敛审核。
检查先复现全局线性残差通过、低密度单元连续性仍为约 1e-5 的输入；
新增 PCG 用例在接线前失败（`check/gas-pcg-red.log`）。接线后继续
一次校正，连续性达到原阈值。原始线性相对／绝对阈值保持。

`check/gas-pcg-green.log` 四项通过：1／2／4 进程连续性审核和 2
进程 Krylov 合同检查。覆盖零 RHS、精确初值、迭代后审核拒绝、终止
拒绝、审核错误及失败时解数组保持。ICCG 的矩阵缩放、对称正定准入、
不完全 Cholesky 和生产入口仍属于接续工作。


G0/G2 ICCG 参考组件：`solver_iccg.hpp` 增加逐行体积缩放和块局部
不完全 Cholesky。分解使用三个低向系数的平方；回代显式取下三角的
转置，正枢轴形成 SPD 预条件器。矩阵 SPD 仍由独立准入负责。
原始 `cgsol.F90` 保持字节一致，外围程序提供单进程 MPI、周期 halo
及工作区，运行其全部分解、CG 与原残差停止循环。

[gas-iccg.json](gas-iccg.json) 记录 19 个确定性输入和程序／源码哈希。
覆盖三个方向周期、非周期、变体积、单格方向、弱对角占优连通域及
隔离固体单位行。二进制可精确表示的系数使矩阵体积缩放差异为零。
原版明确声明 REAL(kind=4) 工作指针；单精度解最大归一化误差为
1.85e-5，FP64 为 4.94e-13，不完全分解逆枢轴最大相对差异为 3.07e-7。
原版 `coef(bpc)` 保留原 RHS，体积缩放进入随后被更新的 CG 残差。
两种 RHS 存储约定分别核对，解另与同一制造解比较。

本组采用原版原行最大残差阈值 2e-6、Hundun 缩放行 L2 阈值
atol=1e-12／rtol=1e-13，服务于正确性对照；迭代数按各自停止标准
记录，效率比较接续使用共同标准。原生探针在全部矩阵上检查精确
对称、M 矩阵行及每个连通分量的约束后才声明 SPD。

`check/gas-iccg-unit.log` 通过预条件器双线性对称、正二次型、独立
稠密 L D L^T 作用、失败枢轴、失效后调用及体积缩放失效保持检查。
此节点交付参考组件；MPI 全局矩阵准入、原行残差反缩放和配置入口
继续按 G2 完成，当前生产后端选择沿用既有行为。

```sh
mkdir -p check/gi
cp /home/wyf/code_dev/src.TCR.dyn711/cgsol.F90 check/gi/cgsol.F90
mpifort -O2 -cpp -Jcheck/gi tools/gas_iccg.f90 check/gi/cgsol.F90 -o check/gi/ref
python3 tools/gas_iccg.py --reference check/gi/ref \
  --hundun b3/versions/v0.4/tests/v04_gas_iccg_probe --runner 'bash check/jam.sh' \
  --source check/gi/cgsol.F90 --output docs/gas-iccg.json
```


G2 分布式 ICCG 准入组件：对体积缩放后的行检查相邻面精确对称、
非负邻接系数、弱对角占优及每个连通分量的严格占优行。先合并本地
连通分量，再通过 halo 传播约束标志。矩阵系数保持原值；此项作为
充分条件使用，一般矩阵沿用适配的 Krylov 方法。

`check/gas-iccg-mpi.log` 的 1／2／4 进程检查通过。覆盖周期、单格
分区、跨分区压力约束传播、IBM 隔离后的多连通域、缺少压力参考的
分量、非对称面系数，以及单一进程出现非有限矩阵行时的集体返回。
全部工作数组由调用方预分配。生产后端配置和原方程单位下的残差
审核继续接线。


G2 ICCG 原生接线：`pressure_linear.algorithm=pcg` 选择 CN/BE 外迭代
中的体积缩放 IC/PCG。分布式 SPD 检查和正枢轴分解通过后进入求解；
原方程 L2 阈值由配置的 atol／rtol 与原 RHS 范数确定，缩放行使用
`min(V)*原 L2 门槛`，同时按各单元体积还原残差并执行原 L2／连续性
审核。缩放求解的残差和原单位残差分别报告。

新增行副本、连通域及分解工作数组在启动时预分配，计入每进程内存
预算。矩阵／系数保持原值，固定热力学压力与一般变密度压力的实际
矩阵分别接受检查。既有 FGMRES／MG 与 DILU 选择保持原配置。

`check/gas-iccg-wire.log` 中 1／2／4 进程原残差检查通过，新增缩放
反例分别隔离原 L2 与低密度连续性两种门槛，并覆盖局部无效体积和
1e±200 数量级范数。原生 CLI 的初次测试定位了输入广播分层：求解器
控制先解码、时间／耦合随后规范化；PCG 的调度准入现由完整生产计划
负责，保持分层解码合同。

`check/gas-iccg-pair.log` 原生入口检查通过：2 进程新算、4 进程恢复、
非对称矩阵拒绝、原残差证据篡改及真实场输出。ICCG 与 FGMRES 使用
相同网格／dt／起始场，首轮原矩阵和 RHS 哈希一致；整步场的最大
归一化差异：U=2.61e-15、pi=4.53e-15、T／rho／h=0，见
`check/gas-iccg-pair-values.log`。两种后端均通过各自真实残差和公共
方程／守恒审核。

`check/gas-iccg-reg.log` 的两项现有入口回归通过：CN/BE 冷态／平均
反应／PaSR 及跨进程数恢复，BE/PISO／SIMPLE，以及固定外迭代／温差
入口。外部证据校验器 self-test 返回成功。该节点完成 ICCG 原生正确性
接线；严格 SPD 充分条件及换算阈值会影响适用矩阵和求解成本，实际
冷态 10＋100 步性能窗口仍按 G2 原出口执行。

G3 原始统计定义：完整 `statistics.F90` 的实际调用条件为
`boffin` 的 `turbstat` 分支。κ 窗口先累计 8 次统计调用的带符号
`dt*1e5*rate`，第 9 次调用计算并清零；该次区间速率在窗口之外。
初次窗口的 κ 为 0.2。Cφ 的计数器从 0 开始，更新发生于第
1、2、9、16……次统计调用。模型时钟与流动外迭代分别管理。

新增 `models_tcr_dyn711` 内核保留该节奏、弱累计速率的单位比值、
带符号比值的正值规则、η 上限及化学／流动时间比超过 100 时的上支。
混合系数使用 `min(κ,1)`；小根采用有理化表达式，η=1 使用退化
线性方程的有限根，η=0 保留其有限上下根。Cφ 独立实现原始双重
空间滤波之间的原值乘积、OH→燃料→混合分数的优先级、范围判断
和邻居平滑。现有 `cdphyso_dynamic_v1` 的 624CF 定义保持原身份。

`tools/gas_tcr.f90` 连接逐字复制的完整 `statistics.F90`，包括
`Dynamic_Cphi` 与 `test_avg`。外围提供 2×2×2 均匀内点、坐标梯度
和单进程归约，当前对照范围为累计窗口、代数／弱反应分支、均匀场
时间尺度及 Cφ 时钟／均匀回退。40 组、每组 27 次调用，共 1,080
次调用的最大归一化差异为 FP64 5.55e-16、原生 REAL4 1.74e-7；
对应日志 `check/gas-tcr-compare.log`，哈希及覆盖范围见
[gas-tcr.json](gas-tcr.json)。空间变化的完整滤波和重新混合单列接续。

```sh
mkdir -p check/gd
cp /home/wyf/code_dev/src.TCR.dyn711/statistics.F90 check/gd/statistics.F90
gfortran -O2 -fcheck=all -Jcheck/gd tools/gas_tcr.f90 check/gd/statistics.F90 -o check/gd/ref4
gfortran -O2 -fcheck=all -fdefault-real-8 -fdefault-double-8 -Jcheck/gd \
  tools/gas_tcr.f90 check/gd/statistics.F90 -o check/gd/ref8
python3 tools/gas_tcr.py --float check/gd/ref4 --double check/gd/ref8 \
  --hundun b3/versions/v0.4/tests/v04_gas_tcr_probe --runner 'bash check/jam.sh' \
  --source check/gd/statistics.F90 --output docs/gas-tcr.json
```

`Dyn711History` 以独立版本身份保存每个单元的累计速率、κ 上下支、
有限退化根状态、Cφ 及统计时钟。所有候选计算从已接受窗口起步，
显式完成各单元／组分后才形成待提交记录，提交使用预分配缓冲交换。
原生逐单元记录格式支持通用恢复器按全局单元重分发。
`check/gas-tcr-history.log` 通过 20 步候选丢弃／新 dt 重算、8→9
窗口交界、恢复后字节一致、单单元记录重分配、固体保持、部分提交
检测及损坏时钟／分支拒绝。这是历史组件验证；原生运行入口及
跨进程数动态模型组合接续接线。

G3 非均匀滤波对照扩展：完整 `Dynamic_Cphi` 在 67 组密度／体积／
混合分数模板、536 个内点上执行原始两次滤波、系数选择与平滑。
原始 `test_D2` 的供体项为 `rho*V^(5/3)`，`test_RD2G2` 的供体项
为 `rho^2*V^(5/3)*|grad(q)|^2`；二者均除以共同质量权重和。
初次实现把后一项的密度也纳入 5/3 次幂，完整例程对照直接暴露了
差异，现已按源表达式的括号作用域修正。

`check/gas-mix-compare.log`：C++ 与 FP64 原例程的滤波比值及最终
Cφ 差异均为 0；REAL4 对应最大归一化差异分别为 6.62e-5、1.83e-6。
REAL4 的原值乘积存在相消，按精度单列记录。均匀梯度／变化密度／
变化体积与随机标量均纳入模板；当前外围 ghost 乘积固定初始化，
MPI／IBM 的生产空间交换继续由后续接线检查。原累计窗口与历史
定向检查仍通过，见 `check/gas-mix-unit.log`。

G3 公共滤波修正：继续以 624CF 冻结来源
`check/s6/statistics.F90` 的完整 `Dynamic_Cphi` 和其原辅助函数
独立对照。三条标量通道采用共同模板及单位分子量，以隔离权重、
有界供体乘积及二次滤波；67 组模板的 C++／FP64 差异为 0，
见 `check/gas-cf-compare.log` 与 gas-tcr.json 的 `cf_filter`。
这确认两类参考采用相同的 `rho*V^(5/3)` 与
`rho^2*V^(5/3)*|grad(q)|²` 权重。

现有 `DynamicTcrPlan` 对上述两项使用了 `(rho*V)^(5/3)` 与
`(rho²*V)^(5/3)`，密度随体积一同乘幂。当前生产路径已复用经
完整例程对照的公共 `dynamic_filter_moments`，保持 624CF 有界
乘积、组分分组、四步更新及共同事务。动态模型方法身份加入
`favre-volume-power-v2`；既有检查点继续保留其原来源身份，后续
迁移按方法恢复合同处理。

`check/gas-filter-green.log` 的四项检查通过。1／2／4 进程直接
检查生产 `DynamicTcrPlan`，覆盖全局前向边界、周期面、奇数
分区、单格分区和固体供体剔除。与整域串行结果的差异为 0；
密度整体乘 2 后的系数差异为 0，500 个系数位于裁剪区间内部。
JL4 真实机理原生五步计算及 1→2 进程恢复也通过，继续执行原
组分、元素及能量门槛。初次测试的分区夹具使用了超出周期供体
准入范围的小网格，现采用满足两格周期供体要求的 7×5×5 网格；
原生比较工具与程序按共同方法身份重新构建后完成恢复对照。

624CF 完整例程的额外构建入口（源文件保持原内容）：

```sh
python3 - <<'PY'
from pathlib import Path
source = Path('check/s6/statistics.F90').read_bytes()
Path('check/gd/cf.F90').write_bytes(source[source.index(b'      subroutine Dynamic_Cphi'):])
Path('check/gd/cf_stub.f90').write_text('subroutine statistics\n error stop 73\nend subroutine\n')
PY
gfortran -O2 -fcheck=all -fdefault-real-8 -fdefault-double-8 -Jcheck/gd \
  tools/gas_tcr.f90 check/gd/cf.F90 check/gd/cf_stub.f90 -o check/gd/cf8
```

`gas_tcr.py` 的原命令追加
`--cf-double check/gd/cf8 --cf-source check/s6/statistics.F90`。
该入口仅调用完整动态滤波例程，外围为串行归约和坐标梯度；字段
分子量转换及完整反应调度仍按 G3 生产组合出口接续。

G2 被动标量生产接线：原 CN/BE 输入在 ProductCompiler 返回 10201；
公共标量方程和恢复字段已具备。当前接入终态质量通量下的 BE 守恒
方程，已接受的 rho/q 历史保持固定；完整空间残差沿用配置的重构，
迎风校正矩阵配合 DILU/FGMRES。候选、固体行、边界和最终非对流率
进入共同事务。多个标量复用 72 字节/单元的矩阵与因子工作区，内存
预算、实际迭代、scalar 阶段计时及外部证据审核同步登记。

首次接线检查定位并修正两处合同：工作区预算采用已创建的标量角色
目录，方程计划在其后编译；普通被动标量采用变密度守恒行，ESF 继续
使用冻结密度的非守恒形式。标量校正复用的热方程临时数组在全域能量
审核完成后使用。最终原方程残差先除单元体积再以局部 rho/dt 标度
归一化，收支按体积积分，门槛均为 128 epsilon。

`v04_gas_passive_cli` 使用两个有符号仿射关联示踪量，检查入口输运、
温差驱动的变化密度、周期面、真实 VTK 值、2→4 进程恢复以及证据
篡改。仿射误差 2.22e-16，连续/恢复最大归一化差异 1.54e-11；
原方程残差 1.84e-15，收支缺陷 1.75e-16。IBM 用满足周期供体范围
的 12³ 网格，216 个固体单元的两个值逐位保持，流体侧扩散残差
5.46e-16。记录为 `check/gas-passive-values.log`、
`check/gas-passive-ibm.log`；共享 ESF 输运的 2 进程检查通过，
见 `check/gas-passive-test.log`。输入夹具采用原 STL 扫描所需预算、
周期范围和正向入口速度；这些准入条件与生产算法分别核查。

G3 混合分数来源复核：`fieldpdf.F90:421–479` 先由统计组分的比摩尔数
重建 Bilger 或单元素混合分数；`boffin.F90:396–429` 随后在
`icycl(nvf)>0` 时调用普通 `condif → cmod → step(rhobar)`。
统计位于本步压力/速度校正之后。rsfz-143 GTMC 的实际 `input.d:53`
为 `20 1e-3`，与已冻结 SHA-256
`b38285918f14e2ea7eea31f42ee9cd99903a3e60b7dd25473f3d28707444b849`
一致，只读副本在 `check/gas-gtmc-input.d`。因此原版动态滤波的混合
分数来源包含元素重建、同一本步历史的 CN 输运及终态取样。
本节点的通用 BE 被动标量按 G2 能力交付；G3 接续显式接入该 CN
时间层与元素定义，并固定实际参考程序的对应调度。

`check/gas-passive-final.log` 的原生标量组合与公共方程数值检查
2/2 通过；公共检查包含冻结密度/变密度行、物理边界、IBM、守恒
及校正回退。外部证据校验器 self-test 成功。正式构建记录为
`check/gas-passive-final-build.log`。本节点覆盖普通外边界及密闭
IBM 壁面；切面注入口的独立标量状态沿现有显式准入合同接续扩展。

监看收尾：`monitor.jsonl` 的 `payload` 显式输出五个 passive 指标，
接受步屏幕摘要输出实际被动标量迭代数。`check/gas-passive-monitor.log`
通过原生 2→4 恢复与 IBM 组合，逐接受步核对监看与证据的五项值完全
一致；监看条目直接复用本步报告。

G3 混合分数 CN 公共算子：`ScalarMidpointView` 为公共被动标量方程
提供已接受/端点平均值。`rho*q` 保持完整守恒存储；对流、扩散和
隐式源项使用中点，端点对空间矩阵及隐式项的响应系数为 1/2。
完整校正行同时缩放对流响应，IBM 固体行保持单位对角和零修正，
流固连接沿公共约束移除。中点字段身份进入方程证书，ESF 冻结密度
接口继续使用其独立时间合同。端点调用沿用既有 BE 运算顺序。

`check/gas-scalar-cn-final.log` 的公共方程与原生被动标量检查 2/2
通过，详值保存在 `check/gas-scalar-cn-values.log`。变密度二次函数
的独立解析存储/对流/扩散行差异 1.04e-17；完整残差差分与校正矩阵
的归一化差异 9.05e-15。IBM 及周期 IBM 的对应差异 2.85e-14。
时间细化采用正弦不变子空间，每一步通过两次独立公共残差组装求解
幅值，并核对全部单元原方程；与半离散精确指数解比较，时间阶为
2.00305、2.00076。该检查隔离时间误差，空间精度另按专项门槛登记。

本节点完成 CN 公共算子。生产被动标量当前继续采用前一节点的 BE
时间层；G3 接续工作为中点工作区与公共面 VLS 系数、逐次调用顺序、
终态审核及 Restart 方法身份的共同接线。元素重建及源检查点的
混合分数身份随实际 GTMC 配置一起处理。


G2/G3 被动标量 CN 生产接线：中点存储复用已完成热收支审核的候选
工作区，VLS 面系数由端点迭代计算，再用于中点输运。物理扩散与 VLS
取共同最大系数，校正行包含 1/2 响应；固体值、原方程审核与终态
非对流率继续进入共同事务。各标量依次复用既有面数组及预分配矩阵。
`passive_scheme=CN` 同步进入 check、证据和 monitor；BE 调度身份保持。

完整只读 `src.TCR.dyn711/vls.F90` 的 44 个模板涵盖双向通量、物理
扩散、常值／痕量、组分闭合、热变量及有符号单标量。原版 FP64 与
当前 upwind_constraint 公共面策略的差异为 0；默认 REAL 在近乎平坦
梯度处产生离散分支差异，逐例保存，其他模板差异处于 2e-6 内。
该入口覆盖均匀笛卡尔内面；IBM、MPI 与边界由原生组合检查承接。
数据见 [gas-vls.json](gas-vls.json)，复现：

```sh
python3 tools/vls.py --dyn711 --coast /home/wyf/code_dev/src.TCR.dyn711 \
  --native b3/versions/v0.4/tests/v04_solver_species_conservation_test \
  --runner 'bash /home/wyf/code_dev/flow/check/jam.sh' --output check/gv
```

`check/gas-passive-cn-test.log` 的原生 CN 入口通过变密度、有符号仿射
示踪量、周期边界、2→4 进程恢复和 IBM 组合。仿射误差 4.44e-16，
连续／恢复最大归一化差异 1.54e-11，原方程残差 6.79e-16，积分
收支缺陷 2.08e-16；216 个固体单元值逐位保持。详细输出保存于
`check/gas-passive-cn-values.log`。时间阶证据沿用前述公共算子细化，
完整原生瞬态细化继续按 G2 范围登记。

方法身份加入 `passive-cn-vls-v1`。此前原生 BE 标量检查点
`check/ps-i4m6spwj/r3/Restart` 的直接恢复返回 10213；显式方法
恢复完成 step 2 并写入新身份，随后 4 进程直接恢复完成 step 3。
三份记录为 `check/gas-passive-identity.log`、
`check/gas-passive-recovery.log`、`check/gas-passive-recovered-exact.log`。
来源检查点保持原样。该节点为原版混合分数 CN 输运提供公共生产能力；
元素重建、dyn711 统计调度和独立动态模型历史继续按 G3 出口接线。


G3 物种坐标对照：扩展完整 `statistics.F90`／`Dynamic_Cphi` 外围
驱动，使燃料、OH、产物和氧气使用实际非单位分子量。逐字原版与
C++ 原式诊断的 FP64 最大归一化差异分别为 dyn711 `3.68e-13`、
624CF `3.22e-15`。继承的统计时钟／根、混合分数滤波及单位分子量
对照继续通过。原式仅对 RD2G2 的组分梯度乘分子量平方，而 G2 保留
比摩尔数梯度；统一质量坐标后的系数变化见 gas-tcr.json 的
`species_coordinates`。原版程序源码保持原样，诊断只改变外围输入。

用户于 2026-09-17 明确采用一致的质量分数坐标。Hundun 的现有
624CF 滤波保持此定义，dyn711 接线同样采用该定义；原式作为具名
诊断对照，系数差异作为模型差异登记。运行时的质量／元素／能量
审核继续采用既定阈值。日志为 `check/gas-tcr-species.log`。

G3 历史与空间计划接线：`ProductTcrHistory` 现在承载独立 dyn711
记录身份，转发预提交、提交、回退、固定 V4 恢复及组合 V5 内层
恢复。部分统计窗口、各组分分支和 C_phi 时钟随完整单元记录迁移。
定向检查覆盖窗口恢复、提交前回退、异模型记录与损坏时钟拒绝。

`DynamicTcrPlan` 复用预分配坐标／梯度／乘积工作区和全局格号供体
交换，新增 dyn711 的原始乘积、OH→燃料→混合分数选择和中心权重
为 2 的邻域平滑。双方共享相同的密度／体积矩实现，624CF 保留其
缩放乘积及分组处理。dyn711 更新依赖已接受的调用时钟；交换前
共同核对跨进程时钟，时钟冲突返回统一错误。流体梯度在物理边界
取单侧差分；dyn711 后向滤波采用域内流体供体，边界和固体外的
供体权重为零，该处理作为笛卡尔／IBM 扩展明确登记。

`check/gas-dyn-final-test.log` 的五项检查通过：dyn711 历史接口，
1／2／4 进程滤波，以及正式程序 JL4／4 场 ESF 的五步新算／恢复。
滤波覆盖 7×5×5 网格、奇数分区、单格分区、周期面、固体掩码、
非均匀密度、整体密度缩放和重复候选。两种滤波的分区及缩放差异
均为 0；dyn711 首两次更新、第三次保持以及跨分区时钟冲突检查
通过。详细输出为 `check/gas-dyn-final-values.log`，正式构建为
`check/gas-dyn-final-build.log`。此前历史夹具显式设置完整 V4 历史
标识后通过；默认缺失历史标识按恢复合同正确拒绝。

该节点完成分布式动态滤波和通用历史接口。dyn711 的模型配置、
终态统计调度、化学／流动时间尺度、η 与混合分数 CN 来源、重新
混合及 8／16 场原生组合仍由 G3 后续节点交付。


G3 局部时间尺度与统计接线（2026-09-17）：用户明确批准
`tau_I=rho*k/epsilon`、`tau_K=sqrt(mu/epsilon)`、
`tau_flow=sqrt(tau_I*tau_K)`。k、epsilon 来自现有 Vreman，epsilon
使用 W/m³，mu 使用 Pa·s。零耗散取无限时间，正耗散且 k=0 取零
积分／流动时间；根选择继续采用 chemical_time > 100*tau_flow。
局部 SI 算子覆盖宽数值范围及量纲共同缩放。

只读完整 statistics.F90 的 flow 外围驱动复现：末格速度 3→12，
八格 tau_flow 从 0.0577350269 降至 0.0288675135；rho=4 时
原式为 0.0408248290。原版额外密度因子与末格索引的影响分别记录。
`tools/gas_tcr.py` 将这些反例及 67 组本地 SI 时间尺度输入纳入
同一可复现入口，局部公式最大归一化差异为 2.22e-16。
既有 1080 次完整原版统计、空间滤波和非单位分子量对照继续通过。
记录及程序哈希为 gas-tcr.json，日志为 check/gas-flow-reference.log。

DynamicTcrPlan 新增集体 stage_rates，读取终态物理质量分数、候选
密度／动力黏度／速度梯度以及各自单次区间 PDF／PSR 比摩尔率。
逐格查询 Vreman SGS 状态，逐组分归约有效化学时间及最大比摩尔数，
然后暂存 dyn711 率窗口。既有统计工作区复用一个分量保存流动时间；
空间滤波在其后复用同一内存并封装候选历史。输入检查、跨进程统计
时钟／组分数／dt／阈值一致性先于率写入，错误由全体共同返回。
原版化学弱反应规则保持：全域均弱时为 1e-12 s，1 s 为回退上限。

check/gas-flow-release-build.log 正式构建通过；
check/gas-flow-release-test.log 五项定向检查通过：局部时间／历史，
1／2／4 进程统计与空间滤波，以及正式 JL4／4 场 ESF 新算／恢复。
MPI 覆盖两个九次调用窗口、生成／消耗／全域弱反应、静止／简单
剪切／Vreman 状态、变密度、IBM 掩码、候选回退重算、单进程坏值、
跨进程 dt 和时钟冲突。分区与 rho/mu 共同缩放差异均为 0；完整
输出为 check/gas-flow-release-values.log。测试预期中的全域弱反应
时间已按原版初始化下限修正，新的独立上下分支检查通过。

本节点交付局部时间尺度及公共 MPI 统计计划。dyn711 原生模型配置、
PSR/PDF 一次区间与守恒重新混合、混合分数 CN 来源、终态调用位置
及 8／16 场组合保持下一接线节点。现有 cdphyso_dynamic_v1 的 JL4
回归用于守住已接入的生产行为，证据范围与新模型完整验收分别登记。


G3 原生 dyn711 生产接线（2026-09-17）：新增独立 `dyn711_v1`
配置、MPI 广播、方法指纹及运行摘要，进入正式 CN/BE 外迭代入口。
Bilger 混合分数使用命名的被动标量及公共 CN/VLS；η 使用已接受的
物理均值比摩尔状态。隐式混合使用已接受的逐物种 κ 和单格 Cφ，
焓采用燃料 κ。终态压力校正、物性和 Z 输运之后执行统计，候选
历史随全部参与者共同提交／回退。颗粒本地预提交状态先集体确认，
随后各进程统一进入统计通信。

PDF 与 PSR 各推进同一单次化学区间。PSR 以物理平均组成和焓作
PH 初态；源版本的 field0 焓保留为方法差异记录。重新混合将随机场
合并到 PSR 终态及共同平均焓，field0、物种源和模型率窗口承接同一
实际增量。原源版重新混合后的重复增量由该单区间定义替代。方法身份
同时包含终态统计、Bilger CN 和守恒重新混合；运行说明见 tcr.md。

八场 MPI 敏感性调查：同一第八步检查点，2／4 进程在全局单元
(6,3,0) 得到 k=6.03038e-16／6.61047e-16 m²/s，流动时间分别为
5.36751e-6／5.61975e-6 s。相同化学时间 5.5889988844e-4 s
跨过 100*tau_flow 阈值，κ 因此为 0.9999982573／0.3709445656。
独立 60 位九子式计算复现 SGS 数值，确认差异来自极弱横向速度梯度
及选根阈值。率窗口、η、Cφ、网格映射和化学时间在两种分区间一致。
诊断记录为 check/gas-dyn-replay2.log、gas-dyn-replay4.log、
gas-dyn-minors.json；原八步输入及检查点保留在 check/o8、o8s、
o8r、o8x，修复重放输出为 check/dk2、dk4。

用户批准方程分辨率一致的 k=0 局部极限。实际两组 nu_t 为
4.33226e-22／4.97218e-22 m²/s，FP64 中 mu+rho*nu_t 均等于 mu。
`dyn711_resolved_flow_times` 对此使用 k=0，其余状态使用原局部公式；
SGS 计算与输运系数保持原实现。新增规则登记到方法身份
`resolved-sgs-k-limit-v1`。原失败检查点先在同一旧方法身份下重放，
随后正式身份执行全新 4/8/16 场测试，清晰区分数值修复证据与方法迁移。
原检查点重放历史差异为 0，流场最大归一化差异 1.394e-14，
面质量通量最大差异 6.939e-18 kg/s；检查日志为 gas-resolved-compare.log。
局部回归覆盖原分支分歧、已分辨 SGS、零耗散及坏输入；1/2/4 进程
统计检查通过，见 check/gas-resolved-unit-test.log 前四项。

原生 4/8/16 场均使用 JL4、8×8×8 网格及 GTMC 原 dt：
2 进程推进八步，同一检查点分别以 2／4 进程推进第九步。
历史头和离散状态精确相同，率／κ／Cφ 的归一化差异为 0；物理场
最大归一化差异 1.394e-14，面通量最大差异 6.939e-18 kg/s。
归一化质量缺陷最大 1.724e-16、能量缺陷最大 2.350e-14，
元素缺陷最大 3.543e-10，保持原门槛。零热释放检查采用 1100 K
O2/N2、入口 O2=0.23／初态 0.232 及成对 ±1e-4 的 O2 波动，
一个原 dt 后各随机场完全相同；质量、能量和元素审核通过。
证据汇总为 gas-native.json；可复现入口为 CTest 的 v04_dyn711_*_cli。
该节点覆盖小规模原生生产接线和跨分区恢复；热释放分支长期轨迹、
真实 GTMC 动态续算、BE 调度扩展和完整专项出口继续逐项推进。

近静态夹具另暴露两项后续 G2 诊断：300 K 均匀 O2/N2 的净物种
收支约 2e-16～7e-16 kg/s，使用净项构造的归一化分母将相对缺陷
报告为 1；入口组成变化后，冷焓近参考值时终态残差在约 1e-16
附近停留并触发现有上限。原始失败日志分别为
check/gas-resolved-unit-test.log 和 check/gas-native-final-test.log。
1100 K 的本节点重新混合检查沿用全部原阈值，300 K 诊断继续单列。

收尾构建 check/gas-native-commit-build.log 通过。集体预提交分支明确后，
check/gas-native-commit-test.log 的八场及重新混合两项通过；4／16 场、
app_case、既有 JL4 使用此前适用于本节点的通过证据。各运行程序哈希
随 gas-native.json 保存，测试参数、阈值及原检查点定位一并登记。


G2 近静态收支与审核调度（2026-09-17，承接 3cc8f94）：当前程序
分别复现均匀 300 K 场的近零物种收支拒绝，以及入口组成变化场的
64 次外迭代上限。输入位于 check/c3、c4，初始复现日志为
check/gas-300-uniform-red.log、gas-300-gradient-red.log。

发现的第一项具体精度损失是物理账本先将各时间层全域存量转为
FP64，再相减乘以 1/dt。已有方程组装按局部增量计算，因而两者
精度不同。新增真实账本最小反例：rho=0.213389，h=300000 J/kg，
dt=1e-5 s，h 仅增加一个 FP64 间隔；原时间率为 1.45519e-6 W，
独立增量为 1.24209e-6 W，误差约 17.16%。
check/gas-balance-red-test.log 记录失败。现按单元先形成质量／总能量
增量，再进行全域归约；绝对存量继续按原格式报告。焓、密度、压力
和动能的单间隔回归及现有 CN 能量／边界算子检查通过。

第二项为临时求解目标兼作接受门槛。全域能量失败会加严求解目标，
后续在小于 FP64 量级的目标上跳过实际收支复核。现保留正常区间
已验证有效的加严精度；审核调度和局部筛选目标为
min(原配置门槛,max(临时求解目标,epsilon_FP64))。本地修正仍使用
临时求解目标；候选提交同时通过原方程、EOS、物种及实际全域能量
1e-6 门槛。直接取消正常区间加严会使现有跨分区率比较退化，
check/gas-balance-final-test.log 保存该中间诊断；最终实现保留此精度，
原 1/2/4 进程组成账本及 2→4 恢复重新通过原门槛。

用户批准近零物种／元素收支采用可追溯 FP64 存储包络。
每格的 V、rho、随机场物理平均 Y 使用
u(x)=epsilon_FP64*abs(x)+denorm_min；乘积包络按展开式计算，
前后两层相加除以实际 dt。余组分按总质量与独立组分包络相加，
元素按原子数／分子量传播。全部包络在 MPI 中归约后，与原相对
门槛共同审核。报告保持原相对缺陷，新增 storage_roundoff_bound、
roundoff_applied 和 fp64-local-storage-v1 规则身份，详见 run.md。

1/2/4 进程回归覆盖解析包络、零存量痕量源、单进程坏输入和
缺失源项：1e-3 kg/s 的缺失源及零存量组分的 1e-20 kg/s 源均被
拒绝。均匀 300 K 原生算例以原 dt 完成一步，原始 O2／N2 相对
缺陷仍为 1，绝对缺陷约 2.06e-16／6.82e-16 kg/s，应用包络并
明确报告；物理随机场重新混合后相同。带入口组成变化的 300 K
算例也通过原门槛，终态全域能量相对缺陷约 7.90e-8，目前需要
42 次外迭代，其效率仍作为后续近静态优化项。两例均使用正式
Cantera、原 GTMC dt 和只读参考模型资产。

收尾证据：check/gas-balance-ulp-build.log、gas-balance-gate-build.log；
check/gas-balance-final-test.log 的算子与 MPI 四项通过；最新原生
组合、300 K 两例、组成账本跨分区恢复、固定外迭代及拒绝检查见
check/gas-balance-close-test.log。定向汇总与实际程序身份见
[gas-balance.json](gas-balance.json)。该节点完成上述精度和近零审核
接线，gas.md 全部 G0–G4／I1–I4 的剩余项继续按表推进。

### I1：反应模块与通信观测计时（2026-09-17）

监看现已记录 reaction_sources、mean_reaction、esf_reaction 和
 tcr_statistics 四类完整调用耗时。计时在每次候选开始时重置，随现有
尝试报告累计，覆盖同一步的候选、失败调用和重试；关闭的模型为零。
ESF 化学位于 CN 外迭代之前，终态 TCR 位于审核阶段，新增模块项
与 CN 阶段有明确的包含关系。初始化与 Restart 物理身份沿用原定义。

通信观测复用现有资源计数，分别记录结构化 halo 等待、halo 控制
通信和公共线性归约对象的耗时。各项采用最大进程壁钟，具名范围
与包含关系写入 monitor.jsonl 和 run.md。直接 MPI、IBM 供体及
TCR 内部通信继续由所属模块壁钟承载，后续按 I1 扩展细分范围。

旧程序的缺字段检查得到预期失败。正式构建上的冷态外迭代检查、
8 场 dyn711 2→4 进程续算，以及真实 Cantera 平均场／PaSR
1→4 进程续算通过。覆盖关闭模型零值、有效模块正值、计时位于
整步壁钟内，以及原来的方程、守恒与历史审核。512 格 dyn711
第 9 步在当前两核绑核诊断中，2 进程整步约 0.618 s，ESF 反应
约 0.309 s，TCR 统计约 0.00069 s。这些数值用于计时接线检查，
实场性能验收继续按既定统一条件开展。

额外选取的 v04_species_terminal_refinement 在第二步的第 64 次
外迭代触发 17837。独立构建已提交基线 2362f67 后复现相同结果：
E=8.1439440162779271e-9、Y=8.1775844727605651e-14、反应
输入残差=2.8923774308455583e-8；三者与新增计时程序完全相同。
该既有收敛问题保留原测试与门槛，列为下一项独立诊断。
对照源码位于 check/tb，构建位于 check/tbb；失败日志与二进制
身份保存在 [gas-timing.json](gas-timing.json)。I1 的重试中状态及
进一步通信／写盘细分继续推进，完整专项目标保持进行中。

### G2：公共面限制器在浮点平台区的稳定处理（2026-09-17）

I1 记录的两步平均场失败已完成根因定位与修正。逐格探针显示第二步
第 10 次外迭代的组分变化约 6e-15、温度刷新约 6e-13 K；公共面
(4,7,4) 的附加扩散却由 0 跳到 7.672582437423995e-4，下一次又回到
0。面中间两点的组分差处于一个 ULP 内，上游差约 1.37e-13。
原限制器按中间梯度的瞬时符号形成比值，使公共系数在中心与上风
处理之间切换；有梯度的焓场继承该系数，进而改变能量与压力迭代。
逐项观测见 check/gas-mean-probe.log、gas-mean-face2.log。

修正沿用现有 128*epsilon_FP64 的模板分辨率定义：整体模板仍有
可分辨变化、面中间差处于该范围时，采用中间零梯度的上风限制。
完全平坦坐标与已分辨梯度继续采用各自既有分支。公共面系数仍统一
进入全部组分、焓、矩阵和残差。只读对照为 COAST vls.F90 的梯度
比值和上风系数；此处作为 Hundun 的 FP64 可分辨性修正记录。

最小实场模板在旧实现中随 1 ULP 改动给出 0／0.5／0.5；新增检查
要求三个值均为 0.5，并覆盖三轴、双向流动。完整物种守恒算子检查、
原两步平均场 1/2/4 进程检查及 8 场 dyn711 2→4 进程续算全部通过。
原失败的第二步在三个分区均为 16 次外迭代，物种残差分别为
1.3399174584485115e-15、1.334882960466225e-15 和
1.3308639002962352e-15；原方程、守恒和续算门槛继续生效。
证据与程序哈希见 [gas-plateau.json](gas-plateau.json)。本项关闭
上一条 I1 新发现的两步收敛问题；实场性能和其余专项出口继续推进。

### I1：步内重试的实时状态（2026-09-17）

ProductDriver::advance 增加按调用借用的 DriverAttemptObserver。
每次候选开始时发布不可变的 proposal、attempt、coupling_sweep 和
previous_failure，回调在本进程同步执行；根进程应用回调负责状态
文件原子替换。原来的尝试、标量重算及事务接受／回退继续决定
模型状态。应用在 advance 返回后集体汇总观测写入状态，并保留
已接受步的权威报告。

status.json 现能在长步内部显示 solving／retrying、目标步、此次 dt、
时间步回退／同 dt 重算类型，以及上次失败阶段和原因码。终态、
写盘和保存停止沿用现有状态流程。字段与计时归属见 run.md。

1/2/4 进程定向检查使用真实压力—能量 stage 54 拒绝，观测到 dt
由 0.009 缩至 0.0045，在第二次尝试开始前直接读取 retrying 文件。
该目标与随后一步的完整物理状态同直接半步 BE 对照逐位一致。
状态临时文件占用错误被记录，同时求解仍得到相同物理状态。
普通 CLI 的独立监看、临时输出、保存停止和 2→4 进程恢复检查
也通过。四项汇总为 check/gas-progress-green.log，实际状态记录及
程序身份见 [gas-progress.json](gas-progress.json)。此前 G2 的
公共面平台区修复继续作为基线。I1 其余计时范围与完整专项出口
继续按既定计划推进。


### I3：接受态 ESF 双状态输出目录（2026-09-17）

三维场增加 DensityMeanEOS、DensityStatistical、DensityField0、
TemperatureField0、EnthalpyField0 和 MassFractionsField0 六组数据。
物理均值 EOS 密度与统计密度分别采用 `rho(p,mean(h),mean(Y))` 和
`1/mean(1/rho_f)`；field0 按已接受的正权 EOS 定义查询，原始 h0／Y0
同时保留。完整名称、单位和组分顺序见 run.md。

重建仅发生于三维场输出，普通监看继续使用原生字段。工作区为每个
本地单元 `(ns+5)*8` 字节，加入现有编译内存和写盘容量预算。
SGS 输出复用原工作区；静态 IBM 固体格采用零占位。模型时间推进、
接受历史与方法身份沿用原定义。

`v04_esf_output_cli` 采用真实 JL4／4 场 ESF／动态 TCR 的 512 格用例：
2 进程两步新算，4 进程续算；独立 1／2／4 进程读取器直接从 Restart
的随机场和 field0 重建物性，对照输出值。最大归一化差异为
2.682e-15；两类物理密度的差异约 2.52e-8 kg/m³，确认其独立分工。
测试同时覆盖 Vreman 和层流输出分支。

开启／关闭场输出后的完整检查点逐字节一致；同源续算的 legacy／XML
输出对应原生检查点逐字节一致。二进制载荷按全局坐标逐值比较，
包括新增的七分量 field0 组分数组和分区重叠点。既有 dyn711 8 场
第 9 步检查点由当前程序直接读取，输出重建差异为 2.055e-16。
`v04_spray_fixed_tcr_cli` 完成固定 p0／IBM／煤油 ESF／THICK_EX／
SGS 破碎的输出与 1／2／4 进程恢复检查。

现有 I/O 综合检查定位到一条历史 JSON 字符串：073a5e9 已增加 CFL
定义及方向最大值，测试预期现同步包含这两个字段。该检查恢复通过，
生产证据格式保持原样。初次新增 CLI 检查的输出目录位于算例目录内，
触发已有路径保护；复现入口现采用同级短目录。

证据为 gas-state.json、check/gas-state-wire-test.log（1/1）与
check/gas-state-final-test.log（2/2）。本节点覆盖输出目录、状态定义、
写盘编码及恢复兼容性；VisIt 应用打开、实场迁移和完整气相专项按各自
条目接续。


### I4：PDF 迁移的压力分工与具名组分映射（2026-09-17）

原生固定 p0 已在 S29 接通，PDF 转换器仍按机械压力执行正压准入，
并把机械压力减入迁移能量账本。另一个实场适配条件来自煤油资产：
其物种顺序为 H2／H2O／CO／CO2／O2／N2／C12H23，转换器原先
假设 N2 位于末尾，按连续前缀读取独立组分。实际资产首先在 case
阶段触发 24110，源例程和资产保持原序。

`v04_pdf_import` 现通过物种名称构建独立组分映射，统一用于
热物性、输运物性和原生字段填充；PDF 保持完整源物种顺序。
机械压力按 `p_ref+pi` 保存，固定压力 EOS 使用 p0，参考／原生
能量采用同一 `rho*(mean_h+K)-p_eos` 定义。报告同时记录压力模型、
p0、机械压力基准、组分顺序、源时间及显式历史重建方式。

`v04_pdf_pressure_import_cli` 采用真实煤油热物性、2 场 ESF、
p0=790216.58 Pa 和 512 格。2／4 进程转换涵盖机械压力
−17942736 至 −17935736 Pa，原生压力和独立组分逐值读回一致。
正／负机械压力两组具有相同 h/Y，其固定 p0 质量与能量库存一致；
耦合 EOS 正压路径通过，负压路径按原准入规则返回 24106。

出口匹配的均匀场由 2 进程迁移，再由 4 进程从 V1 第 17 步推进至
第 18 步，t=0.00101 s，原方程与质量／能量审核通过。检查包含
迁移、完整随机场保存、field0 均值重建、原生读回及首步恢复。
记录为 gas-import.json 与 check/gas-import-close-test.log（1/1）。

额外首步探针保留了启动稳健性问题：将极值压力夹具直接连接到
100000 Pa 静压出口，并使用 dt=1e-5 s 时，外迭代速度和能量修正
持续增大，阶段 40 返回 10216，接受状态仍为第 17 步。
check/gas-import-resume.log 保存该窗口。它属于场／边界匹配及
启动收敛的后续分析；实际 624CF 的流体标记、源压力范围和守恒
映射场继续按实场入口审核。

本节点完成独立 PDF 桥接器的压力及物种适配。实际曲线网格到
笛卡尔／IBM 的守恒映射、原版 input.d 与统计历史迁移、统一
import 命令及实场轨迹保留各自验收出口。

### I4：版本化 input.d 参数目录（2026-09-17）

源读取器对照确认：原始 dyn711 采用固定记录顺序，GTMC 和 624CF
运行版具有可选 `restart_regrid_replace_mesh` 记录，并对统计起始步
采用三值读取／两值回退。`tools/input.py` 通过三个已冻结读取器的
SHA-256 选择布局，按 nf=6 解码公共记录，保存原始行、行号和单位。
读取器默认值另列原因；输入／读取器哈希进入输出来源清单。

GTMC 与 624CF 实际只读输入分别得到 43 项公共记录、58／26 项
扩展记录。读取结果分别为 4／2 场 ESF、压力参考 100000／790216.58 Pa，
CFL 窗口均为 0.20／0.30、目标 0.25。源 SGS 名称分别为 Vreman
和 Smagorinsky，624CF 的原生 Vreman 选择沿用已批准策略。

扩展参数与 namelist 保存原始记录，后续按实际模型读取器的搜索顺序、
覆盖规则和默认值接线。目录已给出压力与 CFL 的映射提示；模型、
资产、边界及历史绑定各自保留明确状态。数值按文本十进制记录，
源 REAL 舍入由完整构建精度清单覆盖。

`v04_input_catalog` 的 7 项检查覆盖当前布局、派生旧版布局、可选
记录、统计回退、行级来源、固定 p0／非对称 CFL、扩展保留和错误
准入。实际两例的读取器哈希准入及关键值另有实场检查。证据为
gas-input.json、check/gas-input-test.log（1/1）与 check/ig.json、
check/ic.json。旧版布局检查采用派生夹具；原始输入实场迁移继续
按完整算例入口核对。统一 import 命令及原生模型绑定仍在 I4 后续范围。

### I4：正式程序的 PDF 迁移入口（2026-09-17）

`hundun import TRANSFER --format pdf-transfer-v1 --case CASE --output SEED`
接通现有完整随机场迁移流程。实现归入 app_import.cpp，正式程序与
独立兼容工具共用同一套读取、物性重建、字段映射及 Restart 校验。
新入口绑定显式原生 CASE，格式标识区分已准备的传输数据与原版算例。

输出使用新的独立目录，集体准入检查源／目标路径和现有检查点。
正式入口保持传输文件原样，报告在完整读回及原生初始化成功后保存为
SEED/import.json。报告含质量、能量、压力分工及历史重建身份；能量
另记录绝对变化量与归一化库存，零归一化库存的相对值表示为 null。
空流体库存按 24113 返回，避免生成失去物理意义的归一化报告。
独立兼容工具继续使用 TRANSFER/native.json。

扩展 v04_pdf_pressure_import_cli 覆盖正式入口的 1／2／4 进程，
固定 p0 的带符号机械压力、具名煤油物种、耦合 EOS，以及 2→4
进程第 17→18 步恢复。源传输文件逐字节哈希保持一致；已有输出、
源路径、源内输出、算例路径及重复／未知选项共六类准入检查通过。
全固体夹具验证空库存准入。证据为 gas-import-cli.json 与
check/gas-import-cli-test.log。

本入口承担已映射 PDF 当前态到原生 Restart 的转换。源 input.d
扩展语义、模型／边界自动绑定、统计和颗粒历史、曲线网格守恒映射
及实场轨迹继续按专项各项出口推进。

### I4：具名被动标量迁移与动态历史依赖（2026-09-17）

实场入口核对确认，GTMC 的 dyn711 配置需要独立输运的混合分数 Z。
原 PDF 桥接器将标量数限定为 ns−1，并仅接收物种，因此动态配置
在 case 阶段就会终止。现增加 `pdf-transfer-v2`：头部明确被动
标量个数与名称，passive.f64 携带原始值；编译时分别建立物种与
被动标量的名称映射，字段保存使用各自原生目录。

定向检查使用真实煤油物性、2 场 ESF、固定 p0，以及 Z 和有符号
tracer 两个标量。源顺序 tracer／Z 与原生 Z／tracer 相反；2 进程
迁移后全部单元逐值读回一致，4 进程首步恢复达到第 18 步。
格式标识与头部版本匹配准入、非有限标量准入亦通过；既有 V1
压力和恢复检查继续通过。证据为 gas-import-z.json 与
check/gas-import-z-test.log（1/1）。

实际 dyn711 组合已推进到 restore 阶段，并稳定复现 10217：
ProductTcrHistory 要求完整模型历史及对应格式，当前态 V1 并未
携带这些记录。该窗口保存在 ip/dyn-history.log，测试明确记录为
后续模型历史迁移依赖；被动标量正向恢复检查使用 TCR off。
动态 TCR 迁移继续实施 typed 统计／分支历史与方法恢复的共同入口，
完整原生动态恢复证据保持其既有范围。本节点增加了实场所需数据
承载能力，实场动态续算的出口仍包含上述历史接口。

### I4：模型历史与时间历史分别恢复（2026-09-17）

10217 的接口分析确认，V4／V5 将模型记录与完整时间历史绑定，
而迁移可能拥有完整 TCR 统计、同时需要重建时间层。新增原生
Restart V6 保存当前字段、面通量和固定长度模型记录；其时间历史
标记为缺失，首步采用 BE 恢复。模型身份、逐单元统计时钟、累计
速率和分支仍由类型化恢复入口严格核对。

IO 层同步更新格式准入、布局估算、编码、流式读取、重分区和
检查点查看。V6 省略时间层及方法签名，保留完整模型记录。
V4／V5 继续表示完整时间历史；变长颗粒记录保持原格式。
原生恢复准入、运行证据及外部校验器同步识别两类历史状态。

定向 IO 检查覆盖 V6 的 1→4、4→1、4→4 进程转换、记录逐字节
一致、读取预算、身份／记录尺寸不匹配；同一检查继续覆盖原有
V1–V5。dyn711 类型化检查覆盖已累计 20 次统计的原样恢复、
时间恢复标记错误、统计时钟错误及失败后已接受状态保持。

生产闭环使用 JL4／4 场 dyn711 的第 8 步检查点，在 2 进程下
生成 V6 并核对全部模型记录，再由 4 进程完成第 9 步和原生写盘。
时间历史重建与模型统计保留分别记录，外部证据同时核对源清单。
这补齐了迁移承载格式与原生恢复能力；COAST 统计到原生类型记录
的转换、实际场绑定与完整实场轨迹继续按 I4／G4 出口推进。

证据：gas-history.json、check/gas-history-io.log（1/1）、
check/gas-history-native-test.log（2/2）及外部校验器自检。
新格式的伪造时间签名由外部校验器拒绝；现有模型记录与原检查点
均保持各自来源。首步物理推进通过后曾定位到内外证据校验器的
旧版本范围判断，最终闭环已同时覆盖这些报告路径。

### I4：源历史覆盖与显式模型初始化（2026-09-17）

逐项核对原始 dyn711 的 finish/start_pdf/config_boffin：PDF 文件
保存随机场焓与组分；普通统计文件保存 atime、fstat、ftau、fschem
和 phase_average。累计反应率 sum_w/sum_w0、dstep_ww0、init_ww0
及 dyn_C_count 属于另行分配的模型状态，配置阶段将计数置零。
运行版 finish 与 GTMC／624CF 冻结副本的 SHA-256 一致，另调用
TCR_ROOT_STATE_V2 写入器保存共享 s_kappa 根输运状态。该状态
与 dyn711 的逐物种累计统计分别对应不同模型；两份冻结输入的
root_transport 均为 off。来源与哈希登记在 gas-history-source.json。

正式 PDF 导入增加显式 `--model-history initialize`，支持
dyn711_v1 和 cdphyso_dynamic_v1。类型化构造器在源流动步号建立
目标初始统计，导入器采用目标编译器的模型身份并核对记录宽度，
V6 写入／读回逐字节检查这些记录。dyn711 使用零次统计、κ=0.2；
cdphyso_dynamic 使用 κ=1 和流动步号的更新相位；两者采用
配置的 2/c_z 作为初始混合系数。JSON 报告逐项记录这些定义。

省略显式选项时，动态模型继续要求完整历史来源；原生 V6 路径
保留已有统计。当前态导入、模型初始化与原生历史恢复分别保持
清楚身份。完整实场输入、守恒映射及实际续算按后续入口接续。

定向检查涵盖两种动态模型的全部 512 格初始记录、两种名称顺序
的被动标量、2 进程导入与 4 进程第 17→18 步推进。省略初始化
选项、格式不匹配及非有限标量保持严格准入。证据为
gas-history-init.json 与 check/gas-history-init-test.log（1/1）；
原生历史保留继续使用上一节点 gas-history.json 的独立证据。

### G4／I4：GTMC 实场 dyn711 输入与 IBM 标量入口（2026-09-17）

`tools/gtmc_seed.py` 准备第 32000 步的独立原生算例和 V2 传输数据。
逐文件核对既有 PDF 包哈希，并按 128 份冻结网格／Restart 的来源
清单读入原混合分数。变量 5 与写盘例程追加的 nvf 记录逐块一致，
每个目标单元恰由一个源块覆盖。新目录为 check/z0（算例）和
check/z1（传输）；原始算例、源检查点和 PDF 包保持原样。

几何沿用 153×325×151 的笛卡尔／IBM 资产，21 个质量流量入口、
静压出口、JL4、4 场 ESF、Vreman 和源初始 dt 保持原配置。
新增 dyn711 与具名 Z 输运；边界纯甲烷对应 Z=1，空气对应 Z=0。
局部 CFL 显式采用方向最大值、0.30±0.05。源 Z 的范围为
−7.6108472e−7 至 0.0122582391，完整保留这些原始数值。

128 进程 `hundun check` 返回 invalid_case/15906。逐条边界核对
排除了外边界温度／热类型／标量数值差异，命中的是内部入口
label=703741800、y_min、Z=1：core_patch_inlets_detail.hpp 当前
对 immersed patch 仅准入 species；IbmInterfaceInletState 和相应
字段枚举也只承载物种、焓、速度及动能。因而后续修复需要贯通
被动标量的内部入口状态、矩阵／通量与原方程审核，再恢复实场检查。
本节点保持该准入保护，实际推进接续上述接线完成后开展。

证据为 gas-gtmc-seed.json、check/gas-gtmc-prep.log 和
check/gas-gtmc-check.log。check/z1/prepare.json 保存全部 256 份
源网格／流场身份及传输来源；该数据准备与失败定位属于实场接线
证据，冷态计时及完整燃烧轨迹保持各自验收范围。

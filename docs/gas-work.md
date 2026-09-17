# 气相专项执行台账

入口为 [gas.md](gas.md)。当前实施承接 main eba8a1b 的参考冻结与运行接口节点。

| 项目 | 本次实现与证据 | 接续工作 |
|---|---|---|
| G0 | gas-ref.json：131 文件哈希、关键调用表；gas-cfl.json：完整 courant 例程的 259 输入对照 | 实际参考构建、算例后端与端到端计时窗口 |
| G1 | CFL 定义贯通配置、准入、自适应 dt、运行证据与 Restart 方法身份；固定外迭代参考模式贯通末次完整审计、方法身份和证据；定向检查见文末 | 普通焓 CN 已接线；完整 jstep 调度及非均匀场自适应重试 |
| G2 | Vreman 无散度反例定向回归 | ICCG 原生正确性接线见文末；矩阵／通量、壁函数、精度及大规模同场性能继续 |
| G3 | 继承 S28–S29；dyn711 独立代数／累计时钟、候选历史及完整统计例程的均匀场对照见文末 | 动态空间滤波、原生接线／重新混合，8/16 场动态模型组合（TCR off 复用 S24） |
| G4 | 继承 JL4 原时间步组合检查 | 化学筛选／任务均衡、热源／边界及实场轨迹 |
| I1 | 独立监看、CN 十阶段计时、接受点 stop/output、请求回执、状态查询；2→4 进程恢复通过 | 化学细分和通信包含关系，重试阶段状态通知 |
| I2 | 当前目录启动、run.json、显式覆盖、步数／时间结束，严格模式检查通过 | 简洁物理输入模板及最终方法说明扩展 |
| I3 | legacy FP64 VTK、中心坐标、分区正向重叠、原子时间索引；3/8 进程 VTK 读取与流线通过 | VisIt 应用打开、双状态附加字段目录和迁移序列关联 |
| I4 | Restart 查看、流式校验和、损坏检测通过 | 原版 input.d／历史检查点导入及来源说明 |

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

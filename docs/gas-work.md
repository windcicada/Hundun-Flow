# 气相专项执行台账

入口为 [gas.md](gas.md)。当前实施承接 main eba8a1b 的参考冻结与运行接口节点。

| 项目 | 本次实现与证据 | 接续工作 |
|---|---|---|
| G0 | gas-ref.json：131 文件哈希、关键调用表；gas-cfl.json：完整 courant 例程的 259 输入对照 | 实际参考构建、算例后端与端到端计时窗口 |
| G1 | CFL 定义贯通配置、准入、自适应 dt、运行证据与 Restart 方法身份；固定外迭代参考模式贯通末次完整审计、方法身份和证据；定向检查见文末 | 普通焓 CN 已接线；完整 jstep 调度及非均匀场自适应重试 |
| G2 | Vreman 无散度反例定向回归 | ICCG 分布式准入／生产接线、矩阵／通量、壁函数、精度及同场性能 |
| G3 | 继承 S28–S29 动态 TCR 与恢复证据 | 原始 dyn711 八步统计与分支／重新混合对齐，8/16 场动态模型组合（TCR off 复用 S24） |
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

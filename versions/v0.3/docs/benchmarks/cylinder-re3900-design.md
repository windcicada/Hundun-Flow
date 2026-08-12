# Re=3900 三维圆柱绕流基准设计

## 1. 用途和边界

该基准用于验证 HUNDUN-FLOW 在 Linux CPU 上的静止周期 IBM、WALE、
两次 PISO 修正、MPI 分解和四字段力报告能否在工程网格上共同工作。
当前验收门只要求冻结的全网格算例完成初始化和一个完整时间步。
它不证明时均阻力、脉动升力、Strouhal 数、尾流速度剖面或湍流统计已经
收敛。长时间统计只能在本功能门通过后单独启动。

## 2. 物理与数值配置

- 圆柱直径 `D=1 m`，来流 `U=1 m/s`，密度 `rho=1 kg/m3`；
- `Re_D=rho*U*D/mu=3900`，因此 `mu=1/3900 Pa s`；
- 计算域为 `[-5D,15D] x [-10D,10D] x [0,pi D]`；
- `x` 方向为速度入口/压力出口，`y` 和 `z` 方向周期；
- 圆柱为静止无滑移动边界，使用 LFP ghost-cell IBM；
- LES 使用 WALE，`Cw=0.5`；
- 固定时间步 `dt=0.006 s`，每步恰好两次 PISO corrector。

产品算例在
`benchmarks/cylinder_re3900/cases/hundun/case.json`。启动验收只把
`steps` 改为 1，把 diagnostics 改为每步记录并关闭网格输出；网格、
物理、边界、时间步、WALE 和 PISO 设置不变。

## 3. 网格尺度依据

冻结网格为 `480 x 480 x 48`：

- `Delta x = Delta y = D/24`，圆柱直径上有 24 个等距单元；
- `Delta z = pi*D/48 ≈ 0.06545D`，展向每个直径约 15.3 个单元；
- 展向长度 `pi D` 及三维 Re=3900 设置与公开圆柱尾流研究的常用问题
  定义对齐，但本项目不把“问题定义对齐”解释为统计已经验收；
- `48^3` 只用于快速暴露同一产品数值路径中的构造、压力求解和力诊断
  缺陷，不替代全网格启动门；
- `12^3` 会让当前圆柱 IBM 压力行缺少六个背景面，对该几何不是合法
  算例，不用于验收。

这一网格选择的目的是建立可运行、可追踪的工程基线，而不是借单一网格
宣称网格无关性。

## 4. 力系数和符号

驱动器输出的 `force_operator` 是物体受流体的力，
`force_budget_reaction` 与它反号。运算符力、表面牵引力和 consistency
必须同时有限，且保留 Task 11 已接受的单一符号权威。

对该配置，参考面积为 `D*Lz = pi D^2`，因此每步的瞬时系数为

```text
Cd = 2 * force_operator.x / (rho * U^2 * D * Lz)
Cl = 2 * force_operator.y / (rho * U^2 * D * Lz)
```

启动步的 `Cd/Cl` 只用来检查符号、量纲和输出链完整性，不与文献的
稳态或统计量直接对比。

## 5. 启动验收门

`480 x 480 x 48 / 64 ranks` 冻结算例必须：

1. 完成初始化和一个 committed step，无崩溃、死锁、NaN、collective
   mismatch 或未处理异常；
2. `attempts=1`，`correctors=2`，并通过产品自身的压力和独立残差合同；
3. continuity、pressure residual、WALE summary 和四字段力诊断全部有限；
4. 根据 `force_operator` 记录有限的瞬时 `Cd/Cl`；
5. 记录 exact HEAD/tree、binary/case/STL SHA-256、MPI/toolchain、命令、
   环境、退出状态、时间、峰值 RSS 和日志 SHA-256；
6. 验收后不留本基准的 systemd user unit 或 MPI 进程。

不允许为通过本门改变科学阈值、`dt`、物理配置、网格、PISO corrector
数量，或增加滤波、阻尼和按算例调参。

## 6. 公开科学参考

- Parnaudeau, Carlier, Heitz and Lamballais, *Experimental and numerical
  studies of the flow over a circular cylinder at Reynolds number 3900*,
  Physics of Fluids 20, 085101 (2008), DOI `10.1063/1.2957018`；
- Norberg, *An experimental investigation of the flow around a circular
  cylinder: influence of aspect ratio*, Journal of Fluid Mechanics 258
  (1994), DOI `10.1017/S0022112094003332`；
- Ma, Karamanos and Karniadakis, *Dynamics and low-dimensionality of a
  turbulent near wake*, Journal of Fluid Mechanics 410 (2000), DOI
  `10.1017/S0022112099007934`；
- 端木玉、万德成，《雷诺数为 3900 时三维圆柱绕流的大涡模拟》，《海洋工程》
  34(6), 2016，DOI `10.16483/j.issn.1005-9865.2016.06.002`。

这些文献用于定义问题、理解三维尾流和未来统计比较；本次一步启动
验收不把任何文献统计值设为通过阈值。

# v0.4 Re=3900 圆柱与后台阶一手验证资料核对

日期：2026-08-21

## 1. 使用范围与证据等级

本文只核对 v0.4 Task 20 所需的公开一手论文、作者公开稿和 NASA 原始数据。
为避免把图读数、经验拟合或后续论文的转引误当成直接实验值，以下采用三种标记：

- **已核验原始数值**：原论文表格、原始数据文件或官方验证档案明确给出的数值；
- **由原始公式计算**：公式来自原论文，但在目标 Reynolds 数处的值是本文代入计算所得；
- **待获取/数字化**：原文只有曲线，或者论文明确说明数值数据需向作者索取。

本文不是最终 acceptance receipt。任何用于发布阈值的数据仍应保存源文件、提取脚本、
页码/图号和 SHA-256，并在查看 HUNDUN 结果前冻结。

## 2. Parnaudeau 等人的 Re=3900 圆柱实验

### 2.1 主来源与定义

主来源是 P. Parnaudeau、J. Carlier、D. Heitz 和 E. Lamballais，
*Experimental and numerical studies of the flow over a circular cylinder at Reynolds number 3900*，
*Physics of Fluids* 20, 085101 (2008)，DOI
[`10.1063/1.2957018`](https://doi.org/10.1063/1.2957018)。可核对
[AIP 论文页](https://pubs.aip.org/pof/article/20/8/085101/256405/Experimental-and-numerical-studies-of-the-flow)
和[作者上传的全文](https://www.researchgate.net/publication/234130406_Experimental_and_numerical_studies_of_the_flow_over_a_circular_cylinder_at_Reynolds_number_3900)。

论文定义

\[
Re=U_cD/\nu,\qquad St=f_{vs}D/U_c,
\]

其中 `Uc` 是外部/来流速度，`D` 是圆柱直径，`fvs` 是涡脱落频率。坐标原点位于圆柱
中心，`x` 沿风洞轴，`z` 沿圆柱轴，`y` 按右手系定义。三站剖面横坐标为 `y/D`；
一阶量用 `Uc`、二阶量用 `Uc²` 归一化。

### 2.2 已核验实验配置

- 风洞试验段为 `28 cm × 28 cm` 方形截面、长 `100 cm`；来流湍流强度小于 `0.2%`，
  温度控制在 `±0.2 °C`。
- 圆柱直径 `D=12 mm`、总长 `280 mm`，两端板间距 `240 mm`，所以有效展弦比
  `L/D=20`；堵塞率 `4.3%`；圆柱中心距试验段入口 `3.5D`。
- PIV 在 `z=0` 平面进行，`Uc=4.8 m/s`，对应 `Re=3900`。每台相机使用 `5000`
  对图像；两幅图像的间隔为 `25 μs`。
- 统计剖面使用小视场 PIV case 2；论文给出的两个视场为 `3.6D × 2.9D` 和
  `1.9D × 1.6D`，Table I/正文中的 case 编号排版存在互换迹象，因此冻结原始数组时
  应以正文“statistical profiles come from case 2”及其 `0.025D` interrogation-window
  说明为准，不依据表格行序猜测。
- PIV 即时粒子位移误差小于 `0.3 pixel`，被判为错误并替换的速度向量少于 `0.1%`。
  这两个数不是三站统计剖面的逐点置信区间。
- HWA 的 `x/D=3,5,7,10` 统计与三站近尾流 PIV 不同：论文估计其一阶矩误差约
  `1%`、二阶矩误差约 `4%`。不能把这两个百分比直接赋给 `x/D=1.06,1.54,2.02`
  的 PIV 曲线。

### 2.3 可直接登记的标量

以下均来自原论文，不依赖 HUNDUN 输出：

| 量 | 已核验值 | 来源和边界 |
|---|---:|---|
| 实验 Strouhal 数 | `St = 0.208 ± 0.002` | HWA；误差主要来自约 `1%` 的 `Uc` 测量精度 |
| PIV 回流长度 | `Lr/D = 1.51` | Table II；`Lr` 定义为圆柱后缘到中心线平均流向速度变号位置的距离；表中没有为此 PIV 值给出独立 `±` 误差 |
| PIV 中心线最小速度 | `Umin/Uc = -0.34` | Table II |
| PIV 流向脉动形成长度 | `L_<u'u'>/D = 0.87` | Table II；定义为圆柱后缘到中心线 `⟨u'u'⟩` 最高峰位置的距离 |

论文对 **LES 统计收敛** 的专门估计是：12、52、120 个脱落周期时，`Lr` 的粗略误差
分别约 `±0.5`、`±0.2`、`±0.1`（约 `±30%`、`±12%`、`±6%`）；在 52 周期窗口中，
脉动速度方差最大值的相对不确定度约 `±10%`。这些是有限时长 LES 的收敛估计，
不是 PIV 仪器误差，也不应自动作为实验 acceptance band。

### 2.4 三个近尾流剖面站

原论文的三个站均为 `x/D=1.06, 1.54, 2.02`：

| 图 | 量 | 对称性/定性检查 |
|---|---|---|
| Fig. 11 | `⟨u⟩/Uc` | 近圆柱为 U 形，向下游演化成 V 形 |
| Fig. 12 | `⟨v⟩/Uc` | 关于 `y=0` 反对称 |
| Fig. 13 | `⟨u'u'⟩/Uc²` | 关于 `y=0` 对称；`x/D=1.06` 有两个剪切层峰 |
| Fig. 14 | `⟨v'v'⟩/Uc²` | 关于 `y=0` 对称，峰值位于中心线附近 |
| Fig. 15 | `⟨u'v'⟩/Uc²` | 关于 `y=0` 反对称 |

**待获取/数字化：** 论文没有为这些曲线提供数值表。脚注 31 明确说实验和数值统计可
通过联系作者获得。因此优先级应为：先向作者索取原数组；若无法取得，再从原 PDF 的
Fig. 11–15 数字化。数字化时必须分别保存 PIV 曲线，不得混入同图的 LES、Lourenco–Shih、
Kravchenko–Moin 或 Ma 等曲线，并登记像素标定、采点密度、对称性检查和数字化误差。
目前没有一手依据支持为每个 PIV 剖面点统一附加 `±10%`。

### 2.5 2026 年公开的补充 PIV 数据集

INRAE 在 Recherche Data Gouv 发布了 Georgeault 与 Heitz 的
*Non-time-resolved PIV dataset of flow over a circular cylinder at Reynolds
number 3900*，DOI `10.57745/DHJXM6`，版本 1（2026-02-19）。数据页说明其包含
2D2C PIV 瞬时速度场、均值与 Reynolds stress、偏度和峰度，并明确称其“补充和改进”
Parnaudeau 等（2008）的既有数据库。数据采用 Etalab Open License 2.0，可通过
Dataverse API 访问。

它是很有价值的一手补充数据，但不能仅凭描述就当作 2008 年 Fig. 11--15 的原始数组：
必须先核对实验批次、圆柱原点、`D/Uc`、空间窗口和统计样本，再决定能否在三站插值比较。
因此当前 receipt 只登记 DOI/API 身份，不把新数据场改名成 Parnaudeau 2008 曲线。

## 3. Norberg 圆柱力和 Strouhal 数据

### 3.1 直接实验来源

一手报告是 C. Norberg，*Effects of Reynolds Number and a Low-Intensity Freestream
Turbulence on the Flow Around a Circular Cylinder*，Chalmers Publikation 87/2 (1987)。
可核对[研究机构记录](https://portal.research.lu.se/en/publications/effect-of-reynolds-number-and-a-low-intensity-freestream-turbulen/)
和[作者报告全文](https://www.researchgate.net/publication/272090433_Effects_of_Reynolds_Number_and_Low-Intensity_Freestream_Turbulence_on_the_Flow_Around_a_Circular_Cylinder)。

报告附录把 `Tu=0.1%` 的离散 `St=f_sD/Uo` 和 `-3 dB` 相对带宽列成表。在 Re=3900
两侧，已核验的原始点为：

| `Re=UoD/ν` | `St` | `Δfs/fs` | 说明 |
|---:|---:|---:|---|
| `3704.6` | `0.2108` | `0.25%` | `D=3.99 mm`，lower-subcritical |
| `4211.1` | `0.2104` | `0.21%` | `D=3.99 mm`，lower-subcritical |

报告没有恰好 `Re=3900` 的离散点。若在这两点间做线性插值，结果必须标成“派生插值”而
不是实验记录。

同一报告 Table 3 的直接压力数据只在 `Re=3000` 和 `8000`（`D=6 mm`）列出：在
`Tu=0.1%` 下分别为 `St=0.213, C_Dp=0.98, -C_Pb=0.84` 和
`St=0.204, C_Dp=1.13, -C_Pb=1.05`。这里 `C_Dp` 是由表面压力积分得到的平均压力阻力，
不是包含摩擦的总平均阻力；不能把 `Re=3000` 的 `0.98` 直接登记为 Re=3900 总阻力。

### 3.2 Norberg 2003 经验式

补充来源是 C. Norberg，*Fluctuating lift on a circular cylinder: review and new
measurements*，*Journal of Fluids and Structures* 17 (2003) 57–96，DOI
[`10.1016/S0889-9746(02)00099-3`](https://doi.org/10.1016/S0889-9746(02)00099-3)，
可核对[作者单位记录](https://lup.lub.lu.se/search/publication/48a4572b-eeaf-4749-9823-a53e5ffe0813)
和[论文正文页](https://www.sciencedirect.com/science/article/pii/S0889974602000993)。

该文 Table 4 对实验汇编给出分段经验式。对 `1600 < Re < 1.5×10^5`：

\[
x=\log_{10}(Re/1600),\qquad
St\simeq0.1853+0.0261\exp(-0.9x^{2.3}).
\]

对 `1600 < Re < 5400` 的局部、单位展向长度极限下 sectional RMS lift：

\[
C'_{L,sectional}\simeq0.045+3.0x^{4.6}.
\]

代入 `Re=3900` 得到 **由原始公式计算** 的
`St≈0.208884`、`C'L,sectional≈0.083046`。这两项是跨多项实验拟合的中心估计，
不是 Re=3900 的直接测量；原式没有为这两个代入值提供置信区间。`C'L,sectional` 也不等于
有限周期展向域上整体力信号的 `Cl_rms`，除非明确处理展向相关长度和平均定义。

### 3.3 不能直接冻结的常见数值

后续 LES 文献经常把 `Cd≈0.98–0.99 ±0.05`、`Cl_rms≈0.10 ±0.05`、
`-Cpb≈0.88 ±0.05` 汇总为 “Norberg (1987), Re≈3900” 数据。这些数字没有在上述
Norberg 1987 可检索表格中以 `Re=3900` 直接列出；其中有些经 Kravchenko–Moin、
You–Moin 等论文再转引。若 Task 20 要使用它们，必须在 receipt 中明确标为“后续论文转引”，
或者回到 Norberg 原始曲线进行有误差记录的数字化。不能把它们与上面的直接离散点或
Norberg 2003 经验式混为同一证据等级。

## 4. Driver–Seegmiller/NASA 后台阶

### 4.1 主来源和标准零倾角配置

实验论文是 D. M. Driver 和 H. L. Seegmiller，*Features of a Reattaching Turbulent
Shear Layer in Divergent Channel Flow*，*AIAA Journal* 23(2), 163–171 (1985)，DOI
[`10.2514/3.8890`](https://doi.org/10.2514/3.8890)。NASA 的
[NPARC Backstep Study #1](https://www.grc.nasa.gov/www/wind/valid/backstep/backstep01/backstep01.html)
和 NASA/TMBWG 的
[2D Backward-Facing Step validation page](https://tmbwg.github.io/turbmodels/backstep_val.html)
给出了公开原始数据及后来修正的 CFD 解释。

v0.4 应使用顶壁角 `α=0°` 的平行壁配置；原实验还测量其他发散角，不能混入零倾角
validation 数据。

### 4.2 已核验几何与流动参数

- 台阶高 `H=12.7 mm (0.5 in)`。
- `H / tunnel-exit-height = 1/9`，`H / tunnel-width = 1/12`。因此按比例推得出口高度
  `9H=114.3 mm`、宽度 `12H=152.4 mm`；后两个是由已核验比例计算的派生值。
- TMBWG 修正说明：`Re_H≈36,000`，旧页面曾错误写成 `50,000`；台阶前动量厚度
  Reynolds 数为 `Re_θ=5000`，边界层厚度约 `1.5H`。
- 台阶前、约 `x/H=-4` 的中心通道参考速度 `Uref=44.2 m/s`；TMBWG 的可压缩 CFD
  复现以调背压获得约 `M=0.128`，原问题仍被归类为 essentially incompressible。
- NPARC 档案网格覆盖 `x/H=-105` 到 `+50`，这说明官方计算域选择，不等同于原风洞
  唯一强制域长。

### 4.3 可直接使用的数据和不确定度

NASA 保存了 Driver/Seegmiller 提供的
[原始数据说明与表格](https://www.grc.nasa.gov/WWW/wind/valid/backstep/backstep01/bstepdata.txt)。
其零倾角数据包括：

- `BACKSTP.DAT`：回流长度、`Cp`、`Cf` 和 Pitot 速度；
- `BACKSTP0.LDV`：`α=0°` 的 LDV 数据；
- `BACKSTP6.LDV`：`α=6°` 的 LDV 数据，仅用于发散壁研究，不能混入 v0.4 零倾角门。

**已核验原始数值：** 零倾角回流长度由 laser oil-flow interferometer 的壁面摩擦及
`Cf=0` 插值得到：

\[
x_r/H=6.26\pm0.10.
\]

原始 `Cf` 表逐点给出相对不确定度 `dCf/Cf`；`Cp` 表给出原值。TMBWG 还提供经过整理的
[Cp 数据](https://tmbwg.github.io/turbmodels/Backstep_validation/cp.expnew.dat)、
[Cf 数据](https://tmbwg.github.io/turbmodels/Backstep_validation/cf.exp.dat)和
[剖面数据](https://tmbwg.github.io/turbmodels/Backstep_validation/profiles.exp.dat)。
TMBWG 说明绘图时 `Cp` 被整体平移，使约 `x/H=40` 处为零；数据文件同时保留原始值和
平移值，比较时必须登记选用哪一列。

整理后的剖面文件包含 `x/H=-4,1,4,6,10` 五站；字段和缩放为：

- `y=y/H`, `u=U/Uref`, `v=V/Uref`；
- `uu=1000⟨u'u'⟩/Uref²`, `vv=1000⟨v'v'⟩/Uref²`,
  `uv=1000⟨u'v'⟩/Uref²`；
- `uuu, uvv, vuu, vvv` 等三阶矩乘以 `10000/Uref³`。

**不确定度边界：** TMBWG 的整理剖面文件没有给每个 LDV 剖面点提供误差条，不能自行
添加统一百分比。可直接冻结 `xr/H` 的 `±0.10` 和原始 `Cf` 表中的逐点相对误差；剖面
若要形成定量 acceptance band，需进一步从原论文/原始实验元数据取得测量不确定度，或把
数值/实验误差分开登记。

## 5. 对 Task 20 数据冻结的直接结论

1. 圆柱 `St` 首选 Parnaudeau HWA 的直接值 `0.208±0.002`；Norberg 2003 在 Re=3900
   的 `0.208884` 只能登记为独立经验式交叉检查。
2. 圆柱 `Cl_rms` 必须注明是 sectional 还是有限展向整体力。Norberg 经验式给出的
   `0.083046` 是 sectional 拟合值且无 Re=3900 专属置信区间，不能伪装成直接实验点。
3. Parnaudeau 三站的五类 PIV 剖面必须索取原数组或数字化 Fig. 11–15；当前公开论文没有
   可直接复制的数值表，也没有逐点 `±10%` 实验误差依据。
4. 圆柱 `Cd_mean` 的常见 Norberg 汇总值需要追溯/数字化后再冻结；现阶段可核验的一手
   Table 3 仅给出 Re=3000、8000 的压力阻力，不足以建立 Re=3900 总阻力 oracle。
5. 后台阶应直接使用 NASA 零倾角 ASCII 数据，无需从论文图中重画；核心标量是
   `xr/H=6.26±0.10`，剖面站固定为 `x/H=-4,1,4,6,10`，并保持 NASA 的归一化和
   `Cp` 平移约定可追踪。

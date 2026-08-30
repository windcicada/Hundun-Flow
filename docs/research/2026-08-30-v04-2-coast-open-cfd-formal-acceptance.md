<!-- SPDX-License-Identifier: Apache-2.0 -->

# V04-2 从动量全局 `theta=0` 到正式接受：COAST 与开源 CFD 一手源码研究及执行方案

研究日期：2026-08-30

研究基线：HUNDUN-FLOW `c7b13fd8eb70742713410922d73a418dee65de56`，分支
`codex/v04-pressure-enthalpy-c1-production`。研究时工作树中另有并行任务正在编写的
limiter/定向测试诊断修改；本文没有覆盖、修改或把这些未提交内容当成验收证据。

## 1. 结论

V04-2 目前仍只能是 `CONDITIONAL ACCEPT`。压力--焓闭环、同目标时刻 refinement、
回滚以及已有时间收敛证据可以保留；距离 V04-2 数值节点接受还缺两件实质工作：

1. 证明并消除 Re3900 BDF2 路径的**全域高阶动量修正坍缩**；
2. 用动量离散实际消费的质量通量计算**真实最大对流 CFL**，使固定 `dt U/D=0.006`
   也受到可审计的稳定性门约束。

这里不应使用“速度正性”或“动量正性”描述 `theta`。速度分量本来可以为负；当前
limiter 的目标是局部有界性或 invariant-domain 性。热力学正性仍由 `p/T/rho/h/Y`
的独立门负责。

推荐路线不是给全局 `theta` 加下限，也不是先减小 Re3900 时间步。应按以下顺序推进：

1. 用唯一、MPI 一致的 limiter winner 证书重放高阶端点、低阶端点、边界/IBM 约束和
   面修正分解，先找出第一个代数不相容项；
2. 若端点或边界/IBM 语义错误，只修这个错误并复测；
3. 若端点语义正确但单个局部极值仍把全域 `theta` 压到零，则把全域标量改为
   **面局部、守恒的高低阶通量修正 limiter**；
4. 用实际质量通量建立固定/自适应时间步共用的 CFL 审计；
5. 通过最小 focused、MPI、Restart、MMS/时间精度和短产品测试后冻结候选，再执行项目
   已注册的 full2、COAST full20 性能和 HUNDUN 长统计门。

正式接受**不要求每个局部面 limiter 都大于零**。圆柱 IBM 附近某些面取
`alpha_f=0` 可以是合法的局部保护。必须拒绝的是 BDF2 步中整个活动域的高阶修正
保留量为零，或 limiter 在光滑 MMS 上破坏既定收敛阶。

## 2. 当前实现给出的可证伪诊断

### 2.1 为什么一个局部单元会使全域 `theta=0`

当前 HEAD 在
[`versions/v0.4/src/solver_momentum.cpp`](../../versions/v0.4/src/solver_momentum.cpp)
的 `1097--1291` 行，对每个活动流体单元和速度分量构造

```text
M_i   = max(a_i, sum_j a_ij)
uH_i  = [b_i + (M_i-a_i) u_i] / M_i
uL_i  = [b_i + (M_i-a_i) u_i + delta_i] / M_i
```

然后以 `uL_i`、当前中心值和上风 donor 值形成局部区间，求每行允许的 `theta_i`，再用
MPI `checked_max(1-theta_i)` 得到一个全域 `theta`。最后所有活动行都按同一个值混合
高低阶 RHS。只要一个单元/分量满足“`uL_i` 已在区间端点，而高阶修正继续向区间外走”，
该行就给出 `theta_i=0`，整个域便退化到低阶端点。

因此，全局 `theta=0` 只证明有一个限制行，不证明所有高阶面修正都不合法。它是明确的
数值质量问题，因为它让一个局部 IBM/边界事件关闭全域高阶对流；它本身不是负密度或
压力--焓闭环再次失效的证据。

### 2.2 第一个必须验证的端点合同

高低阶对流差 `low_order_rhs_delta` 在同一文件 `961--980` 行形成；随后
`998--1007` 行才调用 `IbmEquationInterfacePlan::constrain_momentum`。IBM 实现在
[`versions/v0.4/src/solver_ibm_equations.cpp`](../../versions/v0.4/src/solver_ibm_equations.cpp)
的 `873--1123` 行修改界面流体行的 diagonal、residual 和 RHS，并把 solid 行改为单位行。

这不自动证明存在 bug：如果 IBM 变换对高阶和低阶方程完全相同，且与对流差可交换，
原 `delta_i` 仍然正确。必须用证书直接验证，而不是凭调用顺序推断。对 winner 行应逐项
检查：

```text
C_IBM(E_high + delta_conv) == C_IBM(E_high) + delta_conv
```

同时检查物理边界和活动面删除是否也满足相同关系。第一个不相等项才是根因。优先级为：

1. 高/低端点是否确实只相差已声明的一阶/高阶对流积分；
2. IBM 约束后的 diagonal/RHS 是否对两个端点一致；
3. inlet、outlet/backflow、MPI/periodic 和 IBM 面是否进入了正确的 donor/bound 集；
4. `delta_i` 是否能逐面重构，且单元和 MPI 面的符号、单位和和式完全一致；
5. 上述合同都成立后，才把问题分类为“合法的局部限制被错误地全域化”。

### 2.3 当前固定时间控制的证据缺口

[`versions/v0.4/src/bc_time.cpp`](../../versions/v0.4/src/bc_time.cpp) 的
`TimeSchemePlan::local_candidate` 在 `246--255` 行对 `fixed` 直接返回
`initial_dt`，不会读取任何局部稳定尺度。自适应模式虽在 `257--276` 行消费五类单位-CFL
尺度，但产品路径尚没有给 Re3900 固定步提供实际最大对流 CFL 证书。因此配置中的
`convective_cfl=0.8` 目前不是固定步的运行时门。

## 3. 一手源码对照

所有公开仓库都固定到完整 commit。这里只提取数学、数据所有权和验证原则；不复制其
代码、控制流、命名或架构。HUNDUN 的既有许可证边界见
[`THIRD_PARTY.md`](../../THIRD_PARTY.md)。

### 3.1 冻结 COAST：有 CFL 基础设施，但 Re3900 特例绕过了它

冻结性能基线为本地 COAST commit
`3c22e0f029db1b2ca045ec9e212a95eacbcfe6a3`。不可变内容应以如下命令读取，而不是把
当前可能已修改的 COAST 工作树当作证据：

```text
git -C /home/wyf/code_dev/Coast_software show \
  3c22e0f029db1b2ca045ec9e212a95eacbcfe6a3:<path>
```

一手位置如下：

- [`SRC.Coast/courant.F90`](/home/wyf/code_dev/Coast_software/SRC.Coast/courant.F90)，
  冻结版本 `35--69` 行：按三个方向用面质量通量、密度和单元 Jacobian 计算 Courant
  数，并以 `MPI_Allreduce(MPI_MAX)` 得到全局最大值；同时计算扩散数和 Peclet 数。
- [`SRC.Coast/app/coast_legacy_driver.F90`](/home/wyf/code_dev/Coast_software/SRC.Coast/app/coast_legacy_driver.F90)，
  冻结版本 `1790--1827` 行：普通路径调用 `courant`，CFL 越界时立即按目标 CFL 缩步，
  增长受 `time_step_growth_limit` 和 `time_step_max` 限制。
- 同一 driver 的 `1793--1805`、`1810--1824` 行：`coast_re3900_transient_enabled()`
  分支把 `cfl=0.0`，既不调用 `courant`，也不自适应时间步。
- [`SRC.Coast/coast_screen_summary.F90`](/home/wyf/code_dev/Coast_software/SRC.Coast/coast_screen_summary.F90)，
  冻结版本 `248--259`、`442--449` 行：保存并输出 CFL、三个方向 Courant/扩散数和调整前后
  `dt`。Re3900 特例输出的零值来自绕过路径，并非实际计算结果。

可直接移植的原则是“按活动控制体计算、MPI 取最大值、减步立即、增步受限、结果可见”。
不能照搬的部分包括 COAST 的方向面最大公式、旧 Fortran/SIMPLE 架构，以及 Re3900
特例的 `cfl=0`。特别是，现有冻结 COAST Re3900 证据**不能**证明 HUNDUN 的实际 CFL，
也不应为性能对照增加新的 COAST 科学工作量。COAST 仍只作为短性能基线，不是数值
oracle 或长统计 authority。

### 3.2 OpenFOAM-dev：低阶有界通量与面局部修正分离

固定源为 OpenFOAM-dev
[`09951d6e6b9e0fb45834507abf5bc5fea718330d`](https://github.com/OpenFOAM/OpenFOAM-dev/tree/09951d6e6b9e0fb45834507abf5bc5fea718330d)：

- [`MULESTemplates.C#L219-L385`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/09951d6e6b9e0fb45834507abf5bc5fea718330d/src/finiteVolume/fvMatrices/solvers/MULES/MULESTemplates.C#L219-L385)
  明确分开 upwind bounded flux `phiBD` 和 correction `phiCorr`，最终形成
  `phiBD + lambda*phiCorr`。
- [`MULESlimiter.C#L92-L165`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/09951d6e6b9e0fb45834507abf5bc5fea718330d/src/finiteVolume/fvMatrices/solvers/MULES/MULESlimiter.C#L92-L165)
  把正负修正通量按单元累积；
  [`#L305-L420`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/09951d6e6b9e0fb45834507abf5bc5fea718330d/src/finiteVolume/fvMatrices/solvers/MULES/MULESlimiter.C#L305-L420)
  由相邻单元预算选每个面的 limiter；
  [`#L517-L529`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/09951d6e6b9e0fb45834507abf5bc5fea718330d/src/finiteVolume/fvMatrices/solvers/MULES/MULESlimiter.C#L517-L529)
  在 coupled/MPI 面两侧取一致最小值。
- [`CourantNo.C#L39-L83`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/09951d6e6b9e0fb45834507abf5bc5fea718330d/src/functionObjects/field/CourantNo/CourantNo.C#L39-L83)
  从实际面通量计算 `0.5*dt*sum(|phi|)/V`，质量通量时再除密度。
- [`setDeltaT.H#L35-L55`](https://github.com/OpenFOAM/OpenFOAM-dev/blob/09951d6e6b9e0fb45834507abf5bc5fea718330d/src/finiteVolume/cfdTools/general/include/setDeltaT.H#L35-L55)
  对过大时间步立即减小，对增长采用 `1.2` 阻尼并受 `maxDeltaT` 限制。

可移植的是“有界低阶基线 + 面修正 + 相邻单元共同决定一个面系数 + coupled face
一致性”，不是 MULES 源码本身。MULES 面向有界标量，不是可直接复制的隐式向量动量
limiter；OpenFOAM 也是 GPL-3.0-or-later，仅可作公开数学参考。

### 3.3 SU2：局部 reconstruction limiter，不使用全域最小值控制所有点

固定源为 SU2
[`bc15466602a687d6fb796d5df7a12ce3fde0949a`](https://github.com/su2code/SU2/tree/bc15466602a687d6fb796d5df7a12ce3fde0949a)：

- [`computeLimiters_impl.hpp#L130-L243`](https://github.com/su2code/SU2/blob/bc15466602a687d6fb796d5df7a12ce3fde0949a/SU2_CFD/include/limiters/computeLimiters_impl.hpp#L130-L243)
  对每个点、每个变量从直接邻居构造 min/max 和 limiter；周期配对只同步配对点，MPI
  halo 从 owner 发布，不把所有点压成一个全域系数。
- [`CLimiterDetails.hpp#L140-L205`](https://github.com/su2code/SU2/blob/bc15466602a687d6fb796d5df7a12ce3fde0949a/SU2_CFD/include/limiters/CLimiterDetails.hpp#L140-L205)
  展示 Barth--Jespersen/Venkatakrishnan 一类局部限制函数以及光滑区弱化 limiting 的
  设计。
- [`CSolver.cpp#L1745-L1977`](https://github.com/su2code/SU2/blob/bc15466602a687d6fb796d5df7a12ce3fde0949a/SU2_CFD/src/solvers/CSolver.cpp#L1745-L1977)
  依据线性/非线性收敛和 under-relaxation 调整每点 local CFL，并独立全局归约
  min/max/average 供报告。

可移植的是 limiter 所有权局部化、owner/halo 一致性和 min/max/average 诊断。SU2 这里
限制的是非结构网格 MUSCL reconstruction，其 CFL 主要服务伪时间/非线性收敛；它不保证
HUNDUN 隐式动量 RHS 修正的面守恒，不能直接移植算法或阈值。

### 3.4 AMR-Wind：固定时间步仍计算实际 CFL，Restart/平均保留时间语义

固定源为 AMR-Wind 仓库当前 Kynema-SGF 代码
[`a7f6b6ddd4710bed5759a887d49c4fcd341b3783`](https://github.com/Exawind/amr-wind/tree/a7f6b6ddd4710bed5759a887d49c4fcd341b3783)：

- [`incflo_compute_dt.cpp#L33-L220`](https://github.com/Exawind/amr-wind/blob/a7f6b6ddd4710bed5759a887d49c4fcd341b3783/src/incflo_compute_dt.cpp#L33-L220)
  在活动 mask 上计算对流、扩散和源项谱率并全局取最大值；
  [`#L223-L340`](https://github.com/Exawind/amr-wind/blob/a7f6b6ddd4710bed5759a887d49c4fcd341b3783/src/incflo_compute_dt.cpp#L223-L340)
  为 prescribed/fixed `dt` 另设从 MAC 面速度计算 CFL 的路径。
- [`SimTime.cpp#L162-L280`](https://github.com/Exawind/amr-wind/blob/a7f6b6ddd4710bed5759a887d49c4fcd341b3783/src/core/SimTime.cpp#L162-L280)
  将谱率换成候选 `dt`、限制增长/min/max，并在固定模式仍保存实际 CFL；
  [`#L293-L323`](https://github.com/Exawind/amr-wind/blob/a7f6b6ddd4710bed5759a887d49c4fcd341b3783/src/core/SimTime.cpp#L293-L323)
  对违反 CFL 的 fixed `dt` 明确告警。
- [`IOManager.cpp#L330-L338`](https://github.com/Exawind/amr-wind/blob/a7f6b6ddd4710bed5759a887d49c4fcd341b3783/src/utilities/IOManager.cpp#L330-L338)
  和 [`io.cpp#L48-L68`](https://github.com/Exawind/amr-wind/blob/a7f6b6ddd4710bed5759a887d49c4fcd341b3783/src/utilities/io.cpp#L48-L68)
  在 checkpoint 中写回当前及前两层 `dt` 历史。
- [`TimeAveraging.cpp#L87-L126`](https://github.com/Exawind/amr-wind/blob/a7f6b6ddd4710bed5759a887d49c4fcd341b3783/src/utilities/averaging/TimeAveraging.cpp#L87-L126)
  以 `dt` 加权累积，并跳过 restart/非采样时间。
- [`verification.rst#L64-L105`](https://github.com/Exawind/amr-wind/blob/a7f6b6ddd4710bed5759a887d49c4fcd341b3783/docs/sphinx/developer/verification.rst#L64-L105)
  把 MMS 二阶收敛和两个 CFL 下的 Taylor vortex 作为离散验证。

可直接采用的原则是 fixed/adaptive 共用实际 CFL 观测、固体 mask、全局 winner、受限增长、
时间加权平均和 restart 不重复计样。不能照搬的是 `|U|/dx` 的不可压/AMR 谱率公式、
AMReX 数据结构和 EB/AMR 生命周期。HUNDUN 应从自己的质量通量与控制体质量定义 CFL。

### 3.5 PeleC：几何 mask、全局稳定步和运行时 effective-CFL 硬门

固定源为 PeleC
[`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`](https://github.com/AMReX-Combustion/PeleC/tree/bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8)：

- [`Timestep.H#L51-L87`](https://github.com/AMReX-Combustion/PeleC/blob/bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8/Source/Timestep.H#L51-L87)
  排除 covered EB cell，以 `dx/(c+|u|)` 形成局部 hydro 稳定步。
- [`PeleC.cpp#L793-L924`](https://github.com/AMReX-Combustion/PeleC/blob/bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8/Source/PeleC.cpp#L793-L924)
  汇总 hydro/扩散等局部最小步并做 MPI minimum；
  [`#L928-L997`](https://github.com/AMReX-Combustion/PeleC/blob/bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8/Source/PeleC.cpp#L928-L997)
  限制相邻步增长。
- [`Hydro.cpp#L291-L298`](https://github.com/AMReX-Combustion/PeleC/blob/bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8/Source/Hydro.cpp#L291-L298)
  对实际 advance 得到的 effective CFL 大于 1 明确告警，并可硬终止、要求从 checkpoint
  以较低 CFL 重启。
- [`_cpp_parameters#L263-L289`](https://github.com/AMReX-Combustion/PeleC/blob/bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8/Source/Params/_cpp_parameters#L263-L289)
  把 fixed/initial/min/max `dt`、CFL 和最大增长分别建模。

可移植的是活动几何 mask、运行前估计与运行后 effective-CFL 双重检查以及 fail-closed
行为。PeleC 的声学 CFL、守恒变量 Godunov 更新和 AMReX cut-cell redistribution 不适合
直接移植到 HUNDUN 的低马赫压力基方法与 ghost-cell IBM。

## 4. 推荐的 limiter 修复设计

### 4.1 先完成 failure-only winner 证书

正常成功热路径只保留低成本聚合量；当 `theta=0` 或候选被 limiter 质量门拒绝时，再形成
一个有效的、所有 rank 一致的 winner。推荐最小字段为：

- `valid`、attempt、BDF order、target time、global cell index、component、winning rank；
- `a_i`、`M_i`、`b_i`、`delta_i`、`u_center/uH/uL/lower/upper/theta_i`；
- 形成 bound 的每个 donor 值、面通量符号和 boundary/periodic/MPI/IBM 分类；
- mass-flux、geometry、IBM metric、boundary numeric 和 momentum assembly revision；
- 逐面 high-minus-low correction 的有向和，以及它与 `delta_i` 的差；
- IBM 前后高/低两行 replay 的 diagonal、RHS、residual 和 commutation 差。

winner 选择必须使用稳定的全序，例如 `(最大 depletion，最小 global linear index，最小
component，最小 rank)`。不能让多个并列 rank 各自向 `stderr` 打印局部坐标，也不能把
默认零值伪装成诊断。详细 donor/IBM 数据只在 failure path 收集；成功步不新增 collective。

该证书应先在已有小产品 fixture 的真实 `ProductDriver` 路径触发，不能手工拼 report。
RED 的终止条件是定位到第一个失败等式，而不是获得更多 line-search 候选。

### 4.2 端点语义有错时的最小修复

若 RED 证明错误来自 low-order endpoint，应让高、低两个临时方程通过完全相同的：

1. 物理 boundary numeric 更新；
2. IBM 界面 diagonal/RHS/residual 约束；
3. 活动面删除和 solid-row 处理；
4. MPI/periodic donor 发布。

然后以两者**约束后的 RHS 差**作为 limiter correction。可以保留等价的紧凑差分实现，
但必须有测试证明它与“双方程重放”逐位或在严格舍入界内一致。不要为修正一个端点合同
立即引入新的 limiter 家族。

### 4.3 端点正确但仍全域坍缩时，改成面局部守恒修正

如果 RED 证明当前端点和 IBM/边界合同正确，推荐采用 clean-room 的有限体积
flux-corrected transport 结构：

1. 以相同目标时刻状态组装单调低阶动量基线；
2. 对每个活动 fluid--fluid 面和分量形成有向的积分修正
   `A_f,k = Phi_high_f,k - Phi_low_f,k`；
3. 用低阶 predictor `uL_i,k`、相同 majorant `M_i,k` 和局部 admissible bounds 形成正负
   RHS 预算；
4. 对单元汇总 `A_f,k` 的正负贡献，得到不越过上下界的 `R_i,k+/- in [0,1]`；
5. 每个内部面从 owner 和 neighbour 的相应预算中取最小值。为保持向量修正方向，v0.4
   首选再对三个分量取一个共同 `alpha_f`；
6. 把 `alpha_f A_f` 以相反符号写入面两侧 RHS。MPI/periodic 面由唯一 owner 计算并发布
   同一个 `alpha_f`；IBM impermeable 面的对流修正为零；物理边界按 inflow/outflow
   语义一侧处理。

若把 `s_if A_f,k` 定义为面对单元 `i` 的有向 RHS 增量，则预算可写成：

```text
B_i,k+ = M_i,k (upper_i,k - uL_i,k)
B_i,k- = M_i,k (uL_i,k - lower_i,k)
S_i,k+ = sum_f max( s_if A_f,k, 0 )
S_i,k- = sum_f max(-s_if A_f,k, 0 )
R_i,k+ = min(1, B_i,k+ / (S_i,k+ + safe_zero))
R_i,k- = min(1, B_i,k- / (S_i,k- + safe_zero))
```

一个面的系数取决于该修正对 owner/neighbour 是正预算还是负预算。该结构使局部异常只
限制相邻面，并保持内部面修正的等量反号。全局 reduction 只计算诊断统计，不再参与
`alpha_f` 的数值选择。

局部 bounds 必须纳入真实 boundary/periodic/MPI donor；IBM 不应虚构 solid donor。
先使用尺度相关舍入容差。若光滑 MMS 的极值仍被降阶，再引入与空间截断误差同阶的
smooth-extremum allowance，并只用 MMS 冻结参数；不得按 Re3900 结果调参。

### 4.4 新 limiter 的验收量

原有单个 `theta` 应退役为方案 identity 与以下聚合量：

- `minimum_face_alpha`；
- `limited_faces / active_correction_faces`；
- `retained_correction_l1_ratio = sum|alpha_f A_f| / sum|A_f|`；
- limiter winner（仅失败诊断）；
- limiter 是否在光滑测试完全 inactive。

`minimum_face_alpha=0` 本身不拒绝。当 `sum|A_f|>0` 时，
`retained_correction_l1_ratio=0` 才表明整个活动修正域坍缩，任何 BDF2 产品步出现该状态
都拒绝；若高低阶通量本来完全相同，则该比值应标记为 `not_applicable`，不能伪装成零。
长程窗口还应拒绝持续的全域坍缩。无需先验发明一个“limited face 百分比”阈值；局部性
由合成测试和 MPI 分解不变性证明，数值耗散最终由 MMS、时间收敛和 Re3900 物理统计
共同约束。

## 5. 实际 CFL 与时间步方案

### 5.1 HUNDUN 的权威定义

对活动流体控制体 `i`，令 `dot_m_f` 为当前 face mass flux，`s_if` 把存储方向转换为
单元外法向，定义：

```text
Co_out,i = dt / (rho_i V_i) * sum_f max(s_if dot_m_f, 0)
Co_abs,i = dt / (2 rho_i V_i) * sum_f |dot_m_f|
```

对严格守恒、无体积质量源的单元，两者相等。`Co_out` 直接对应 upwind donor 的流出
驻留时间，建议作为动量对流稳定门；`Co_abs` 与 OpenFOAM 风格定义可比，作为交叉审计。
使用的是 HUNDUN 方程实际采用的 Cartesian 控制体体积，不引入 PeleC cut-cell 体积分数。
solid cell 排除，IBM impermeable 面质量通量必须为零，入口/出口和 periodic/MPI 面均
纳入。

应在两个 revision 上计算：

1. **advective CFL**：从本步动量对流组装实际消费的精确 mass-flux revision 计算，是
   稳定性 authority；
2. **committed CFL**：从压力--焓闭环接受后的 final mass flux 计算，证明下一历史层和
   终端物理状态一致。

每个最大值至少携带 `valid`、global cell、rank、`rho V`、outgoing/absolute flux sum、
`dt` 和 flux revision。全局 max 可以并入已有终端 audit reduction；详细 winner 只在
失败时发布，避免正常热路径增加 collective。

### 5.2 fixed 与 adaptive 的行为

固定 Re3900 正式算例继续使用 `dt U/D=0.006`，不得为了隐藏 limiter 或 CFL 问题改成
`0.003`。fixed 模式仍计算实际 CFL：

- 若 `max Co_out <= convective_cfl*(1+64 epsilon)`，保留配置时间步；
- 若超门，正式 benchmark 直接拒绝，不静默缩步，因为缩步会改变冻结的科学/性能工作；
- 非正式诊断运行可以给出稳定步建议，但其结果不能替代 `dt=0.006` 正式证据。

自适应模式把每个局部 `1/rate` 送入现有 `LocalTimeLimits`，MPI 取全局最小稳定步；减步
立即生效，增步受 `maximum_growth` 限制。非正状态或非线性终端门失败继续使用既有事务
回滚和 retry/BE 规则。不能同时用 CFL 预缩步和失败 retry 重复缩放而不记录 lineage。

对正式 Re3900 full/long 路径，任何 retry 或实际 `dt` 改变都使该段证据失效；应修复
算法后从已验证 snapshot 重跑，而不是把改变后的时间网格混入统计。

## 6. 最小实现与验证顺序

以下顺序把日常快速诊断与只对冻结候选执行的正式门分开。已有且不受 limiter/CFL 修改
影响的 pressure--enthalpy、refinement、capacity=6 rollback 和 provenance 证据直接复用。

| 阶段 | 实现或动作 | 最小证据 | 通过条件 |
|---|---|---|---|
| A | failure-only limiter winner 与高/低/IBM/边界 replay | 一个确定性小产品 fixture；一个 no-IBM control | 找到第一个不等式；所有 rank 报告同一 winner；默认值不输出 |
| B | 只修端点合同；若合同原本正确则实现面局部守恒 limiter | 两单元/两面合成测试，IBM/inlet/outlet fixture | 内部修正和为零；bounds 通过；局部坏面不降低无关面；`alpha=1` fast path |
| C | 实际 `Co_out/Co_abs` 计算与 fixed/adaptive 行为 | 均匀通量解析值、边界/IBM mask、MPI winner | 解析值、flux revision、rank winner 和配置门一致 |
| D | 真实 ProductDriver BE→BDF2 短测 | 现有小 fixture，短 no-IBM，短 IBM | 两步接受；第二步 BDF2；无全域修正坍缩；压力--焓五分量门和正性通过 |
| E | 受影响的精度与事务 | 光滑 MMS/temporal、失败回滚、Restart recovery | limiter 光滑区 inactive 且保持既定阶；失败不发布 trial；Restart BE 被标记并跳过统计 |
| F | 必要的 MPI 语义 | 日常只补一次 1/2-rank；最终沿用 4/64-rank 门 | owner/halo face alpha、CFL max、接受/回滚决定分解一致 |
| G | 候选短产品封存 | tests-off build；`dt=0.006` 64-rank full2 | BE 与 BDF2 均真实接受；无 retry；CFL 合格；终端门合格；记录成本 |
| H | 正式 release | 冻结 full20 COAST pairing、HUNDUN 长统计、provenance/DCO | 满足项目第 16 节全部门，不以短测替代 |

日常循环不运行完整 CTest、full20、COAST 或长统计。建议只构建受影响目标并运行：

1. limiter algebra/winner 定向测试；
2. IBM 边界产品 fixture；
3. CFL 解析 fixture；
4. 一个真实两步 ProductDriver reproduction。

只有修改 owner/halo/reduction 语义时才补一次 2-rank。最终候选才运行 tests-off、短
no-IBM、64-rank full2；Clang 测试使用正确 libc++ `LD_LIBRARY_PATH`。这不是减少正式门，
而是避免在诊断阶段重复昂贵、不能定位根因的运行。

## 7. Restart、MPI 和时间精度门

### 7.1 Restart

项目冻结设计在
[`v0.4 architecture design`](../superpowers/specs/2026-08-12-hundun-flow-v0.4-cartesian-performance-architecture-design.md)
第 8.4、14 节明确规定：Restart 不保存上一 BDF2 历史层，重启后第一步使用 BE，标记且
不计入长统计，然后恢复 BDF2。AMR-Wind 保存多层 `dt` 的做法说明另一种合法架构，但
不能静默改变 HUNDUN 已冻结的 Restart 身份。

本次必须验证：

- snapshot 只从 committed state 写出，包含当前 `U/p/h/Y/final mass flux` 和继续控制
  所需状态；
- limiter/CFL 失败后回滚不会污染 snapshot、历史或统计 accumulator；
- Restart recovery step 的 `origin/order` 明确为 BE，统计样本数不增加；下一步恢复 BDF2；
- same-rank 与既有 rank-changing 读取均形成同一物理接受决定；若 face-local limiter
  改变分区面语义，至少重跑一次相应 MPI Restart 测试；
- corrupt/mismatched provenance 继续 fail-closed。

### 7.2 MPI

每个 partition/periodic face 只能有一个 limiter authority。owner 计算或双方分别计算后
必须以确定性最小值合并，并向 halo 发布同一 `alpha_f`。两侧应用等量反号修正。全局
reduction 只用于 `minimum alpha`、计数、保留比、CFL max 和失败 winner，不得重新变成
控制全域解的 limiter。

1/2-rank 的通过条件不是日志相似，而是全局 index、受限面集合、守恒和终端接受决定
一致；浮点归约量按项目既有确定性/容差合同判断。仅当此次修改涉及 rank 语义时补一次
2-rank，4/64-rank 留给冻结候选。

### 7.3 时间精度

精度门分三层，不能只看 `bdf.order=2`：

1. 光滑 MMS 中 limiter inactive，空间/时间误差保持已注册二阶；
2. 产品 BE→BDF2 第二步实际进入 BDF2，且高阶动量修正保留比非零；
3. full/half 在同一物理时间比较 pressure、enthalpy、density、velocity 和 final flux，
   EOS、continuity、energy、mass、gauge 均独立通过。

局部 `alpha=0` 的非光滑/IBM 区域不要求点态二阶；但不能让它在光滑 MMS 中扩散到全域，
也不能用“BDF2 标志存在”掩盖整个对流修正已退化为低阶。

## 8. Re3900 正式接受门

### 8.1 V04-2 数值节点与 v0.4 release 必须分开

V04-2 数值节点可从 `CONDITIONAL ACCEPT` 提升为数值接受，当且仅当：

```text
v04_2_numerical_accept =
    pressure_enthalpy_closed_loop_accept
 && limiter_endpoint_certificate_accept
 && no_domain_wide_momentum_correction_collapse
 && actual_advective_cfl_accept
 && focused_mms_temporal_accept
 && rollback_restart_mpi_accept
 && dt_0p006_full2_BE_to_BDF2_accept
```

这仍不是 v0.4 正式发布。项目冻结的 release predicate 保持为：

```text
release = numerical_correctness_accept
       && robustness_accept
       && coast_short_performance_accept
       && literature_physical_accuracy_accept
       && provenance_accept
```

权威规格是
[`2026-08-12 v0.4 architecture design` 第 16 节](../superpowers/specs/2026-08-12-hundun-flow-v0.4-cartesian-performance-architecture-design.md)。

### 8.2 full2 与 COAST full20

冻结候选先执行 `480x480x48`、64 ranks、`20D x 20D x piD`、圆柱距入口 `5D`、
`Re_D=3900`、`dt U/D=0.006` 的两步 HUNDUN/COAST 短测。HUNDUN 必须在第二步真实走
BDF2，实际 CFL、limiter 保留量、终端五分量门、正性和无 retry 均通过。

然后冻结 exact HEAD/tree、编译器/flags、binary、case、STL、rank map 和绑核。正式性能
仍执行至少五组交替 full20 配对，使用 max-rank timed region、预注册 warmup 和 paired
uncertainty rule；`median(HUNDUN_hot/COAST_hot) <= 1.0` 才接受。加入 limiter/CFL 后应将
聚合归约并入已有阶段，不能用每步额外多次 collective 无意改变性能工作。COAST 的
Re3900 `CFL=0` 仅记录为其基线限制，不补算 COAST CFL，也不拿它作 HUNDUN 数值比较。

### 8.3 长程数值健康与统计稳态

`literature_statistics` 只运行 HUNDUN 的 transverse/spanwise periodic case，使用冻结
WALE、同一网格和 `dt=0.006`。按既有规格至少发展 `150 D/U`，随后统计约
`2020 D/U`（约 420 shedding periods）。每个接受步持续检查：

- `p/T/rho` 正且所有状态有限；
- EOS、continuity、energy、closed mass、gauge 分量门分别通过；
- 实际 advective/committed CFL 合格；
- 无 retry、无时间步改变、无全域 limiter 坍缩；
- 总质量、能量和入口/出口/IBM 通量无长期单调漂移；
- force、solver iteration、limiter limited fraction/retained ratio 无数值突变。

统计 accumulator 以实际 `dt` 加权，Restart BE recovery 和失败 attempt 不计样。不得只凭
瞬时 residual 宣称 statistically steady。无需新增一条长运行：对既定统计窗口做预注册的
连续 block 分析，报告 `Cd_mean`、`Cl_rms`、`St`、`Lr/D` 的 block 均值、积分自相关时间/
有效样本量和 95% sampling interval；窗口前后半段应在组合 sampling interval 内一致，
且不得存在显著线性漂移。判据和可能的延长规则必须在查看 HUNDUN 长统计结果前冻结。

主要物理输出保持：`St`、`Lr/D`、`Cd_mean`、`Cl_rms`、centerline mean velocity，以及
`x/D=1.06,1.54,2.02` 的 mean/fluctuation profiles。

### 8.4 当前独立的文献 authority 阻塞

Parnaudeau 的 Re3900 PIV/LES 一手论文为
[`doi:10.1063/1.2957018`](https://doi.org/10.1063/1.2957018)，项目已从
[`IRISA institution copy`](https://www.irisa.fr/fluminance/team/Carlier/publications/ParnaudeauCarlierHeitzLamballais_2008_POF.pdf)
冻结 15 组 profile。Tadrist 的直接全长 fluctuating-force 一手论文为
[`doi:10.1063/1.857804`](https://doi.org/10.1063/1.857804)，但当前只形成 target-neighborhood
advisory constraint，不能证明 periodic `pi D` span 等价或精确 Re3900 值。

当前
[`v0.4 literature data receipt`](../verification/v0.4-literature-data-receipt.md)
仍为 `complete=false`：直接测得、适用于 Re3900 的 total mean drag，以及能映射到 periodic
`pi D` 的 finite-span `Cl_rms` authority 尚未闭合。合法的下一步是继续取得并审计完整
Bishop--Hassan（[`doi:10.1098/rspa.1964.0004`](https://doi.org/10.1098/rspa.1964.0004)）
或其他满足 apparatus/span/normalization/uncertainty 条件的一手材料，然后更新冻结 machine
receipt。不能用二手 CFD 表、压力 drag、sectional lift、图上近似值或手工把
`complete=true` 代替。

这是一项与 limiter 修复独立的正式发布阻塞。可以并行完成代码、短产品和性能冻结；在
`receipt-validate --require-complete` 通过前，不得启动权威长统计或声称 v0.4 release
accepted。

## 9. 证据、provenance 与 DCO

limiter 从 global scalar 改为 face-local，或新增实际 CFL authority，都会改变数值/证据
语义。应复用现有通用字段和 reduction epoch；只有现有 schema 无法无歧义表达时才增加
最小字段，并升级 policy/schema/provenance identity。旧 oracle 保持冻结，不得由当前候选
重生成来证明自身正确。

最终 receipt 至少绑定：

- exact commit/tree 和 DCO-clean 可融合提交链；
- limiter scheme/version、bounds 规则、owner/halo 规则和源码/测试哈希；
- CFL 公式、gated flux revision、target、最大值/winner 和 fixed/adaptive policy identity；
- compiler、libc++ runtime、MPI、binary、case/STL、rank map、CPU topology 和命令；
- focused/full2/full20/long-stat artifact hash，及未重跑旧证据的可复用理由；
- 文献 machine receipt 和统计方法/阈值的冻结时间。

有 `Signed-off-by` 的最后一个提交不等于整条开发分支 DCO clean。数值候选稳定后再整理
提交链；诊断中不频繁改写历史。Stage 5/6 portable 与 Halo P0 不进入本节点。

## 10. 明确禁止的捷径

- 不设置 `theta_floor>0`，不把 `theta=0` 改成一个很小正数来绕过 bounds；
- 不逐单元 clip `U`，不使用不成对的 cell-local RHS theta 破坏面守恒；
- 不在计算实际 CFL 前降低正式 Re3900 `dt=0.006`；
- 不靠增加 corrector、line-search、阻尼层数或 retry 次数掩盖端点/IBM/边界不相容；
- 不要求所有局部 `alpha_f>0`，也不把单个局部零误判为全局降阶；
- 不把冻结 COAST Re3900 的 `cfl=0` 当作实际 CFL；
- 不把 COAST 当数值 oracle，不复制/翻译旧 COAST Fortran 或 OpenFOAM GPL 源码；
- 不以两步成功、一个低 residual 或最后一个有 sign-off 的提交代替 full20、长统计、完整
  文献 authority 和 provenance；
- 不在 `literature receipt complete=false` 时把长程物理统计写成正式接受。

按本方案完成后，`theta=0` 将从一个含义模糊的全局症状，转化为可定位的局部有界性
决策；固定 `dt` 将有真实 CFL 证书；V04-2 数值节点和 v0.4 正式 release 也将各自由明确、
可复现且不能互相替代的证据闭合。

## 11. 实施结果（2026-08-30）

推荐方案已经落地为 `common_face_afc_v2`。三个速度分量在同一物理面共享一个
face-local alpha，高阶反扩散修正按 owner/halo 规则成对保守；报告量改为全局保留修正
L1 比例和唯一受限面数，不再用全局最小 alpha 代表整步高阶保留量。CFL 同时在 momentum
predictor 的 provisional flux 和 pressure correction 后的 committed flux 上计算，二者绑定
不同 revision、同一配置上限和精确 rank-local view identity。

真实性补强还包括：BDF/variable-BDF2 系数必须与证书中的 `dt` 一致；rank-local
authority 或数值失败在任何 FGMRES collective 和状态修改前达成 MPI consensus；fixed
terminal-CFL 超限在 stage 60 fatal 并逐位回滚，adaptive 路径则缩到 half step，结果与
直接 half step 的公开物理状态及最终质量通量逐位一致。V5 validator 冻结单文件 schema、
run identity、momentum plan、IBM activity 和 CFL limit，并绑定连续 step 的 `delta(time)`
与 CFL `dt`。

最终 tests-off 二进制 SHA-256 为
`babd0c407d5fd36d0663d4e94a462b2bdc54e09e96e33b20cb381a558dfe37f6`。
64-rank `case-full, dt=0.006` 两步产品运行位于
`/home/wyf/code_dev/.benchmarks/v04-2-afc2-re3900-accepted-SmS6Jy/run/evidence.jsonl`，
SHA-256 为
`be1f744b83a3fb33fd9446b67182c271ab9027228b6e45c05f3e7387f90ced58`。
V5 validator 通过；step 1 为 BE、step 2 为真实 BDF2，均无 retry/fallback。第二步
continuity/energy 分别为 `2.626254588291675e-7` 和
`2.1919163778181067e-7`，provisional/committed outward CFL 分别为
`0.3805492212333011/0.3830326558478964`，均低于 `0.8`。第二步有
`25,865,378` 个唯一受限面，但 retained-correction L1 ratio 为
`0.36495771760827195 > 0`，直接排除了原来的全域 momentum 修正坍缩。

因此本研究支持的裁决是 `ACCEPT / V04-2 NUMERICAL PRODUCT CLOSURE`。它关闭本节点的
压力—焓、AFC、双 CFL 和短产品行为，不替代第 8 节列出的 full20、完整文献 authority、
长程稳定性/统计和最终 v0.4 joint release gate。

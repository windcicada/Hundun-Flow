# 方法优先的数值一致性审查

日期：2026-09-06。状态：已完成本报告范围内的方法审查；M1 完成有边界的实现与局部回归，整体数值正确性及 Re3900 长时间稳定性尚未验收。

## 范围和执行边界

按用户最新要求，本轮不以重现第 1297 步为前提，不先提交新算例。此前刚启动的 checkpoint 重放已人工停止，相关 MPI 进程已退出。中断记录见 [STOPPED.md](2026-09-06-1297-fix-evidence/STOPPED.md)；它不是新的数值失败或通过证据。原 checkpoint 和长测目录未修改。

实际工作区为 /home/wyf/code_dev/.worktrees/hundun-flow-strict-coast-parity。HEAD 为 b779bff4bbff7067e691952b25c90cbc86be7795；工作区包含此前 A–E 任务包的未提交修改，不能用 HEAD 单独标识当前程序。先完成静态方法审查，再对源代码和代数已经确认的 M1 编写局部解析检查、修复并运行相关小网格回归；没有把故障重放作为修复前提。未修改 Re3900 的物性、网格、边界、时间步、收敛阈值或迭代容量，未重新启动圆柱计算或长测，未推送 GitHub。

审查沿以下链条展开：守恒方程 → 时间离散 → 动量预测/面通量 → 压力—能量块 → 候选状态 → IBM/物理边界 → 最终验收/历史提交。以下 M1–M8 记录修改前的方法判断；仅文末“实际修复与验证”列出的项目在本轮运行过，不将静态阅读当成测试通过，也不声称覆盖仓库全部源文件。

证据分三级：

- **确认的实现事实/衔接缺口**：可以由调用链和代数直接确定。
- **明确采用的近似**：不能仅因不等于完整 Newton 就称为程序错误，需要说明省略项和适用范围。
- **待验证后果**：长期漂移、当前故障的归因、改动后的速度和稳定性，不能由静态检查单独确认。

## 先把实际求解的方程写清楚

令 Φ 为带方向的面质量流率，D 为面通量散度积分，B 为当前有效 BDF 算子。当前焓方程可概括为：

\[
R_C=V B(\rho)+D(\Phi),
\]

\[
R_E=V\{B(\rho h)-B(p)-\mathbf U\cdot\nabla p
-\nabla\cdot(\lambda\nabla T)-\tau:\nabla\mathbf U-S_h\}
+D(\Phi h_f).
\]

这里 h 是热力学焓，不是总焓；ρh−p 是内能密度。全域总能量账目还必须包含动能及相应边界通量，不能把 ρh−p 直接当总能量。

依据：[solver_enthalpy.cpp](../../versions/v0.4/src/solver_enthalpy.cpp)，evaluate_enthalpy_cell_terms / combine_enthalpy_cell_system，约 946–1054 行；[physics_thermo.cpp](../../versions/v0.4/src/physics_thermo.cpp)，complete_state_impl，约 810–857 行。

当前修正解的是下列近似块系统，而非完整的非线性方程 Jacobian：

\[
\begin{bmatrix}C_p&C_h\\ E_p&E_h\end{bmatrix}
\begin{bmatrix}\delta p\\\delta h\end{bmatrix}
=-\begin{bmatrix}R_C\\R_E\end{bmatrix}.
\]

实现从连续性行消去焓，使用

\[
S=E_p-E_h C_h^{-1}C_p,\quad
b_S=-R_E+E_h C_h^{-1}R_C,\quad
\delta h=-C_h^{-1}(R_C+C_p\delta p).
\]

这组消元符号与源码一致。关键不在于 Schur 公式有没有写反，而在于各块是否表达了候选状态实际发生的变化。源码也明确设置 full_nonlinear_jacobian=false。

## 主要发现

| 编号 | 性质 | 方法问题 | 当前 Re3900 的相关性 |
|---|---|---|---|
| M1 | 确认的衔接缺口；本轮已局部修复验证 | SIMPLE 新的动量预测没有进入 C2 内部面通量基线 | 直接相关；尚未证明是长测失败的唯一原因 |
| M2 | 确认的近似范围 | 候选通量更新密度，而 C_p/C_h 未包含全部通量密度响应 | 直接相关；收敛影响尚未定量 |
| M3 | 确认的近似范围 | IBM 能量残差和修正块使用不同的壁面空间响应 | 直接相关；不能把现有 spatial 模式称为完整 IBM Jacobian |
| M4 | 确认的验收缺口 | 最终没有重新检查动量方程，也没有累计的全域总能量收支 | 直接相关；目前不足以排除误差积累 |
| M5 | 确认的近似范围 | double-diagonal E_h 连局部热力学时间导数也采用了近似 | 直接相关；当前大小取决于实际焓基准和状态 |
| M6 | 确认的指标范围 | 能量门槛是时间项后向误差尺度，不是长期物理能量误差 | 直接相关；不应靠调整分母让故障通过 |
| M7 | 确认的数据流，后果待验证 | 标量按预测密度更新后，p/h 校正又改变密度，但未配套标量守恒校正 | 当前无输运标量，不是本算例直接原因 |
| M8 | 确认的闭合差异 | IBM 预测器存储速率用有界温度重构，目标能量残差用未截断重构 | 直接相关；可能改变预测质量和回退频率 |

### M1：SIMPLE 第二次动量预测与 C2 面通量脱节

调用链已确认：

1. [core_product_freeze.cpp](../../versions/v0.4/src/core_product_freeze.cpp) 约 13662–13672 行重新求解 SIMPLE 动量预测，更新 trial_velocity 和动量矩阵。
2. 同文件约 13681–13685 行把新速度与原 C1 修正后的 provisional_flux 一起交给 C2。
3. [solver_piso.cpp](../../versions/v0.4/src/solver_piso.cpp) 约 2901–2907 行以新速度作为 h_by_a；约 3175–3200 行却直接把输入的 C1 面通量复制为 C2 基线。
4. 约 3201–3238 行冻结偏移 B_f=Φ_C1−A_f I_f(ρ_0 U*_2)。候选通量再使用 A_f I_f(ρ_α U*_2)+B_f，加上压力增量响应（约 4522–4579 行）。

因此，固定 ρ、p、h，只改变第二次动量预测得到的 U*_2，在零压力/焓增量处仍有：

\[
\Phi_{C2}(0)=A_f I_f(\rho_0 U^*_2)
+[\Phi_{C1}-A_f I_f(\rho_0 U^*_2)]=\Phi_{C1}.
\]

这是代数上直接可见的抵消，不需要先跑到故障时间步。C2 连续性右端看不到这一轮新的内部速度预测；能量中的压力功、黏性耗散等却会看到新速度。物理边界 finalizer 能更新它负责的外边界，不能补回整个内部面的这项变化。

纯压力增量 refinement 延用前一次已修正通量有其合理性；**在重新求解动量之后仍采用同一条增量规则，则丢失了新动量预测的面响应**。这也是为何应区分 SIMPLE 新动量 sweep 和冻结动量下的 C2 refinement，而不是把二者都视为“corrector=2”。

外部交叉参考仅用于确认方法结构：OpenFOAM 的可压缩压力校正源码从当前动量方程重新构造 HbyA/phiHbyA，并包含瞬态修正，再组装压力方程。它不证明 HUNDUN 应照搬其全部离散。[OpenFOAM rhoPimpleFoam pEqn.H](https://api.openfoam.com/2312/compressible_2rhoPimpleFoam_2pEqn_8H_source.html)

修复边界：为新动量 sweep 重建与新动量矩阵、当前压力和 BDF 历史一致的面预测通量；压力—焓 refinement 则保留其已提交候选通量。不能只加一个 I(ρΔU)：rAU 改变时，Rhie–Chow 压力滤波和时间修正也必须更新。

进一步检查发现，旧 stationary 探测失败后的回退会再次调用普通 C2 refresh，但普通 C1→C2 授权本来就是一次性的；只有 typed refinement 可以合法再次进入。因此不能把这条旧回退误当成必须保留的可重入接口。本轮去掉了这次提前探测和重复 refresh，让 C2 一次组装、一次建立方向，再由现有精确候选循环判定是否已经收敛。这样也不会第二次把改写后的 face_aux 误当原时间历史。

**尚未证明**：这一缺口是第 1297 步超限的唯一原因，或者修好它就能消除所有长时间漂移。

### M2：密度相关的通量导数没有完整进入块系统

在 C2 的冻结动量候选附近，可把内部面通量写为：

\[
\Phi(\delta p,\delta h)=A_f I_f(\rho(p,h)U^*)+B_f
-d_f(\rho)\,[\![\delta p]\!].
\]

在零增量处，其一阶变化包括：

\[
\delta\Phi=A_f I_f[U^*(\rho_p\delta p+\rho_h\delta h)]
-d_f(\rho_0)[\![\delta p]\!].
\]

候选路径实际更新第一项中的密度；当前 C_p 主要包含 a0 V ρ_p 和压力梯度通量响应，C_h 则仅为 a0 V ρ_h。相应的通量密度响应没有完整进入连续性和能量块。

依据：[solver_pressure.cpp](../../versions/v0.4/src/solver_pressure.cpp) 约 209–268 行；[core_product_freeze.cpp](../../versions/v0.4/src/core_product_freeze.cpp) 约 9681–9726 行；[solver_piso.cpp](../../versions/v0.4/src/solver_piso.cpp) 约 4522–4579 行。

这是当前 quasi-Newton 的真实边界，不是 Schur 消元本身的符号错误。它意味着即使线性系统解得很准，也可能仍要多轮非线性修正。

修复/改进方向：先给出同一候选映射的完整块导数，再决定哪些项可作为受控近似。不能仅往 C_p 加项而忽略 C_h：完整的 C_h 将不再是纯单元对角块，现有按单元相除的消元前提需要重新处理。现有廉价 Schur 可继续作为预条件近似，不应未经分析就整体更换外层求解器。

### M3：IBM 空间响应仍不一致

目标残差确实包含 IBM 重构：

- 压力功由 correct_pressure_work 修正为壁面重构梯度对应的 U·∇p。
- 导热由 correct_zero_normal_diffusion 修正为重构温度 ghost 对应的热流。

而线性块：

- 一旦启用 IBM，E_p 的 pressure_work 绑定为空。这不是只省略壁面几个单元的导数，而是该模式没有装配整域的这组压力功线性响应。
- E_h 仅在活动笛卡尔面上装配温度增量导热，跳过 IBM 界面，没有对应的 donor 温度导数。

依据：[core_product_freeze.cpp](../../versions/v0.4/src/core_product_freeze.cpp) 约 9876–9906、11509–11550 行；[solver_ibm_equations.cpp](../../versions/v0.4/src/solver_ibm_equations.cpp) 1174–1225、1281–1330 行；[solver_pressure_energy.cpp](../../versions/v0.4/src/solver_pressure_energy.cpp) 约 3082–3136 行。

例如冻结 λ 和光滑 donor 分支时，壁面导热行包含 t_f(T_g−T_i)，故方向导数应包含 t_f(Σw_j δT_j−δT_i)，其中 δT_j=δh_j/cp_j。直接跳过该界面不等于这个重构算子的导数。

重要的排除项：当前 IBM 连续性压力算子通过移除不可穿透界面通量连接，与 zero_interface_flux 的质量通量约束相配。不能误报为“连续性也遗漏了同一套 donor 压力扩散”。见 [solver_ibm_pressure.cpp](../../versions/v0.4/src/solver_ibm_pressure.cpp) 约 243–246、347–383 行。

改进顺序：先补清楚冻结材料下的压力功和导热 donor 方向响应；再评估黏度、导热系数及湍流闭合导数是否需要进入同一轮。保持 NASA7 / COAST 原生输运，不使用常 cp/常 μ/常 λ 替代物理模型。

### M4：最终验收没有覆盖所有演化方程

最终路径重新组装了同一待提交通量上的能量残差，并检查连续性、EOS、压力边界、CFL 等；但没有在最终 ρ、U、p、μ 和最终 Φ 上重新组装动量残差。后面的速度 halo 和梯度更新服务于下一步的物性速率历史，不等于动量方程验收。

依据：[core_product_freeze.cpp](../../versions/v0.4/src/core_product_freeze.cpp) 14104–14260、14267–14540 行；[solver_piso.cpp](../../versions/v0.4/src/solver_piso.cpp) 的 audit_pending_final。

因此，当前记录可以证明某些离散方程的归一化残差达标，不能据此证明完整的 Navier–Stokes 离散系统在最终状态都已收敛。特别是 M1 中“新速度、旧内部通量”的组合，可能不被现有门槛区分。

补充方向：最终动量残差、面/单元动量插值关系、全域质量与总能量收支。全域能量应同时记录内能、动能、入口/出口输运、压力/黏性功、导热和 IBM 壁面贡献；把离散截断/数值耗散和未收敛误差分开，不能要求所有格式在任意状态下都逐位总能量守恒。

### M5–M6：时间导数近似和能量指标需要分开处理

精确的局部热力学时间导数为 E_h,time=a0 V(ρ+hρ_h)。spatial E_h 会从普通焓对角中扣除扩散代理、加入 a0 V hρ_h；double-diagonal 路径直接绑定普通焓对角 a0 Vρ+diffusion_diagonal，没有这项密度—焓时间响应。

依据：[solver_enthalpy.cpp](../../versions/v0.4/src/solver_enthalpy.cpp) 940 行；[solver_pressure_energy.cpp](../../versions/v0.4/src/solver_pressure_energy.cpp) 2420–2430 行；[core_product_freeze.cpp](../../versions/v0.4/src/core_product_freeze.cpp) 约 9699、9917–9923 行。

它被明确标成近似，并非错误地声明了完整 Newton。但“省空间 stencil”不必同时意味着省略一个廉价的局部热力学导数。后续应比较两种模式是否共享同一时间 Jacobian，而不是直接延长使用旧 diagonal 模式的轮数。

另一方面，能量归一化使用各 BDF 历史中 |ρh| 与 |p| 操作数的总尺度。它有助于识别大数相消下的后向误差，却不是累计物理能量漂移的上界。等步长 BDF2 的绝对压力部分约为 4pV/Δt；在本算例 p≈107136 Pa、Δt≈1.38×10⁻⁵ s 时，仅这部分乘 10⁻⁶ 就约为 3.1×10⁴ V W。这是阈值尺度的代数估算，不是测得的热源或能量泄漏。

依据：[core_product_freeze_detail.hpp](../../versions/v0.4/src/core_product_freeze_detail.hpp) 65–117 行。应保留后向误差指标，并另列有量纲残差、累计能量账目和时间/空间离散误差；不能更换分母来让当前失败被接受。

### M7：标量守恒不能只看取值范围和字段映射

预测器先得到 (ρq)*，随后以 ρ* 除成 q*。压力—焓候选复制并冻结组分质量分数，最终可把密度改为 ρ_new，而未同步演化标量守恒量。于是最终储存的量是 ρ_new q*，相对预测量多了 (ρ_new−ρ*)q*。

依据：[solver_thermophysical_predictor.cpp](../../versions/v0.4/src/solver_thermophysical_predictor.cpp) 约 3268–3279、3383–3412 行；[core_product_freeze.cpp](../../versions/v0.4/src/core_product_freeze.cpp) 10409–10418、13278–13296 行。

这与之前修复的 FieldId/role 混排不是同一个问题。即使所有字段身份正确、所有质量分数都处于 [0,1]，仍需说明最终密度、组分和质量流率采用怎样的保守时间分裂。当前 Re3900 没有这些输运标量，故不把它列作当前故障原因；扩展到变密度混合物之前必须补这项方法验收。

### M8：预测历史和最终能量残差的 IBM 温度重构不同

evaluate_thermophysical_rates 对温度扩散使用 correct_positive_bounded_zero_normal_diffusion；候选和最终能量残差使用 correct_zero_normal_diffusion。前者把重构温度限制在 donor 最小/最大值间，后者使用原重构值。

依据：[solver_thermophysical_rates.cpp](../../versions/v0.4/src/solver_thermophysical_rates.cpp) 323–330 行；[mesh_ibm_reconstruction.cpp](../../versions/v0.4/src/mesh_ibm_reconstruction.cpp) 1134–1172 行；目标残差调用位置见 M3。

由于能量最终会重新求解，这个差异本身尚不能证明最终能量方程被换掉。它首先意味着下一步的预测历史未必来自最终残差采用的同一空间算子；在重构越界或极值附近，可能影响预测质量、降阶和非线性成本。应明确区分“接受状态的物理速率”与“仅用于保正预测的限制修正”。

## 目前没有证据支持直接重写的部分

- BDF2 的变步长系数符号、ρh 与压力历史项，未发现基础公式错误。当前代码对预测器降阶会向动量及压力—能量传递同一个 effective_bdf；不是各方程随意使用不同 BDF。
- 持久化非对流焓速率已经不含 BDF(p)，避免在下一步再次外推压力时间导数；不能把旧版问题直接套回当前源文件。
- Schur RHS 和焓恢复的代数符号一致。
- 候选状态会重新求 EOS、温度和材料，并以实际候选通量验算能量；不是只看线性残差就提交。
- IBM 不可穿透通量确有显式约束。
- 现有测试已有 Cartesian E_p/E_h 方向导数、耦合时间阶数和重启覆盖。它们不能替代 SIMPLE 新动量 sweep、完整 IBM donor 响应和长期全域账目的验证，但也不能描述成“完全没有相关测试”。

## 求解策略和性能：放在方法一致性之后

当前 Aitken 初始步长由两个全局残差最大值组成的范数推算收缩率，并优先尝试外推候选；Armijo 检查的是合并范数下降，最终验收则要求各分量分别达标。因此合并范数下降、某一个分量上升并不矛盾，也不能单凭这一现象认定方程发散。应评估方向及各分量的实际响应，不先改成“一律 alpha=1”。

目前固定前半 refinement 使用 diagonal、后半使用 spatial。这个规则不判断所省略的项是否已经主导误差。先修 M1，并明确 M2/M3，再决定是否按方向预测误差或非线性停滞切换近似；不先增加 12 次容量。

FGMRES 的 true-residual 校验和恢复机制需保留。native MG 预条件的是压力代理，不是完整 Schur；加入空间 E_h 后，S 中会出现 E_h C_h⁻¹ C_p 的复合空间作用。其谱与 MG 代理是否匹配值得分析，但本轮没有把它的条件数或成本当成新测量结果。既有模块计时另见 [性能审查](2026-09-06-hundun-coast-module-review.md)，本轮不重复运行性能算例。

## 修复顺序与后续验证边界

| 顺序 | 修改目标 | 在方法明确后才安排的验证 |
|---|---|---|
| 1 | M1：区分新动量 sweep 与压力—焓 refinement，重建一致的面预测通量和时间/压力滤波项 | 公共耦合器入口：保持热力学状态不变、改变新动量预测，内部面通量及连续性 RHS 必须响应；重复进入 C2 不得重复加修正；PISO/refinement 行为单独检查 |
| 2 | M4：补同一最终状态/通量上的动量残差和累计质量/总能量观测 | 有独立收支的封闭/周期/开边界小问题；不以“跑完”代替守恒/误差检查 |
| 3 | M2/M3/M5：统一块导数的声明范围，补确认需要的热力学、密度通量及 IBM 方向响应 | 解析或独立差分 Jv：Cartesian、开边界、IBM 分开；冻结 limiter 分支和分支切换分开；不能仅比较两个共享同一错误实现的函数 |
| 4 | M7/M8：澄清标量保守分裂及预测速率闭合 | 非均匀标量、变密度、周期积分守恒；IBM 温度极值与有界分支 |
| 5 | 非线性策略和 MG/Krylov 性能优化 | 同物理推进时间、同模型/网格/边界的单轮比较，记录线性与非线性真实成本 |
| 6 | 回到 Re3900 | 先完成上述方法与局部验证，再安排有限时段的修复后运行；旧 checkpoint 只作带历史来源的继续调试起点，不自动升级为新的长期统计正确性证据 |

除下节明确记录的 M1 局部检查及相关回归外，这些仍是后续验证设计。M2–M8 尚未在本轮实现；不能用把第 1297 步“跑过去”代替方法验收，也不能用 M1 的局部通过宣布整个长测问题已经修复。

## 实际修复与验证：仅 M1

### 实现边界

- `solver_piso.cpp`：在现有 coupler 内区分 fresh SIMPLE C2 与压力增量 C2/refinement。前者用新 U、新 rAU、当前压力梯度和面压力差重建内部预测通量；后者继续使用前一次通量作为增量基线。
- C1 的 BDF 面历史偏移保存在原有预分配 face_aux 内，只允许同一时间、相同 BDF 系数的一次 fresh SIMPLE C2 消费；此后该存储转为候选通量偏移。没有新增全场数组，也没有增加公共入口。
- 固定物理边界仍由既有边界路径负责，IBM 仍应用原不可穿透面通量约束。没有改变 NASA7、输运模型或任何终止容差。
- `core_product_freeze.cpp`：移除前述 stationary 提前探测的重复 C2 入口。真正已收敛的基线仍可由精确候选路径接受；均匀场的第二次线性解仍是 zero_rhs、0 次迭代。其他接近收敛状态可能多做一次方向准备，尚未测量这项成本，不宣称本轮已提速。
- `v04_flow.hpp` 仅补充接口语义注释；测试沿用 `PressureVelocityCoupler` 和 `ProductDriver` 公共接口。`codebase-design` 的影响限于把这两种通量更新规则收进现有耦合器，没有向产品调用者增加模式开关。

### 独立解析断言

在 4×4×4 单位周期域上，使用常密度及已知动量矩阵来检查面离散恒等式。这是测试夹具，不是替换 Re3900 的物性模型。

1. **新动量脉冲**：保持热力学状态和 C1 通量不变，仅在 x=1 单元列令新 Ux=2。SIMPLE 相邻两个内部面的质量流率必须变为 1/16，且连续性 RHS 在相邻单元分别为 −1/16、+1/16；压力增量型 PISO 仍保留其零基线。生产修复前，此断言失败；修复后通过。
2. **同历史、新 rAU、新压力滤波**：BDF2 为 (1.5,−2,0.5)，已接受/前一层 Ux 分别为 3/0，C1 速度为零，得到 C1 面流率 1/4。C2 保持速度为零，把动量对角加倍，并设周期压力为 [0,1,0,−1]。指定面应得到 1/8+1/24−1/12=1/12。仅给旧通量加 I(ρΔU) 的补丁会仍得到 1/4，无法通过此断言。
3. **状态历史约束**：C1→C2 擅自更换 BDF 系数被拒绝；普通 C2 重复消费 C1 授权被拒绝。合法 typed IBM refinement 由产品级回归覆盖。

代码：[solver_piso_temporal_order_test.cpp](../../versions/v0.4/tests/numerical/solver_piso_temporal_order_test.cpp)，新增两个 `test_simple_*` 函数。

### 已执行结果

构建为 `build-review-fixes`，Clang 15.0.6 / libc++，Release `-O3 -DNDEBUG -march=znver3 -mno-fma -ffp-contract=off`，未启用 fast-math，HYPRE 关闭。仅重建相关测试目标；没有部署或替换原长测可执行文件。每种测试配置执行一轮，不做性能取中位数。

| 测试 | 最终结果 | 验证范围 |
|---|---|---|
| v04_solver_piso_temporal_order | PASS，0.33 s | 原有变步长时间阶数、BDF 面历史，以及上述新增解析断言 |
| v04_solver_piso_authority | PASS，0.87 s | 耦合器授权、状态及通量身份 |
| v04_solver_piso_checkerboard / mpi_2 / mpi_4 | 3/3 PASS，0.28 / 0.28 / 0.31 s | 原有压力棋盘模态检查 |
| v04_solver_piso_mpi_1 / 2 / 4 | 3/3 PASS，0.47 / 0.42 / 0.38 s | 原有跨分区 PISO 与通信/失败语义 |
| v04_core_product_freeze_mpi_1 / 2 / 4 | 3/3 PASS，1.25 / 0.86 / 3.69 s | 包含 SIMPLE 两次动量预测、IBM refinement、精确残差终止等产品回归 |

共 11 个不同测试配置的最终结果通过。这些耗时是小回归的执行时间，不能推算 Re3900 每步耗时或与 COAST 的速度比。

原始记录：

- [耦合器批次](2026-09-06-1297-fix-evidence/method-m1-coupler-ctest.log)：7 项通过，新增解析测试最初因本轮测试代码错误中止。
- [解析测试修正后](2026-09-06-1297-fix-evidence/method-m1-temporal-final-ctest.log)：1/1 通过。
- [产品批次](2026-09-06-1297-fix-evidence/method-m1-product-ctest.log)：3/3 通过。

验证过程问题未隐藏：新增测试曾写错证书成员名，编译报错后按现有 `dt` 字段修正；随后把夹具 ghost 宽度写死为 2，而 central2 夹具实际只分配 1，造成测试自身越界。改成读取字段的实际 ghost 宽度后通过。上述问题是新写测试的错误，不是原长测故障证据，也不计作数值方法的 RED。真正的 RED 是修复前 SIMPLE 新动量脉冲未进入面通量的解析断言。

未执行：Re3900 重放/中短测/长测、性能比较、M2–M8 的方向导数及累计守恒验证、全仓测试或 sanitizer 验收。当前结论是 M1 的通量衔接已修复并通过列出的局部检查，不是完整算法或长期统计已获验证。`git diff --check` 通过；此前未提交工作保留，未提交或推送本轮修改。

## 本次阅读的关键源文件身份

以下 SHA-256 记录于 2026-09-06 07:37 +08:00，是 M1 修复前的审查基线；“主要发现”的原始行号对应此阶段，不能把这些哈希当作修复后文件身份：

| 文件（versions/v0.4/src/ 下） | SHA-256 |
|---|---|
| core_product_freeze.cpp | 4280de935334c24d917d99a6342ba98b6641418b029502faf942e05f3725bcdd |
| solver_piso.cpp | 181d2d0e8355875f50938cb2300f09d15a9be5642e44607dc3927685289700dd |
| solver_pressure_energy.cpp | f3ed9c250c056ba23df891e3466fdb7abc6bef0630fdf5655e2d8d18e7883aa1 |
| solver_enthalpy.cpp | 87b6a1b319e71ed168a7b55b3c0553c0a1b0bbafb72642a10739972f389d7234 |
| solver_ibm_equations.cpp | 8cfa97e770f58c0e20714af8a0d8f349e2cb6f053fc975f0fcd8a9d3b9e4c153 |
| bc_time.cpp | 804484cbed8ef351dfb877d81ace8945ffda40fda14a8d82951548e15f1f3f6d |

M1 修复与最终解析测试之后的 SHA-256：

| 文件 | SHA-256 |
|---|---|
| versions/v0.4/src/solver_piso.cpp | c8aa72a52ba5fdfb4df23aea3f5b042cf62cdf04ece842bca38848c6eab9dc41 |
| versions/v0.4/src/core_product_freeze.cpp | 55acc77a629e102f4f9adc0660aa38b47cc7fd1d020b0fee7678000abb028a70 |
| versions/v0.4/include/hundun/v04_flow.hpp | 2d79beff8123a6ea06e183aa8079f4bd95ffab457e9bbaab45bc65b01cd3af3b |
| versions/v0.4/tests/support/piso_fixture.hpp | 75b4b5b18ca548b276f11aa714f37315df8e8c1f257be9fdf7cf5c7d8124b5fa |
| versions/v0.4/tests/numerical/solver_piso_temporal_order_test.cpp | b4ac4fd081a9a2352be22d72b67cec5ae185b738872e2e14916bd3e5f440487c |
| build-review-fixes/versions/v0.4/tests/v04_solver_piso_temporal_order_test | 3c91a66012d173f526a32d7829a448c86f3e01dcee9e3bd30b64a64ace1ebf80 |
| build-review-fixes/versions/v0.4/tests/v04_core_product_freeze_mpi_test | af43fe8331fd0704299aae5de6a0f80200217553ab88209727721daaa5321247 |

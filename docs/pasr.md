# 单场 PaSR 默认选择

反应配置中的 `model` 缺省为 `auto`，`ensemble.fields` 缺省为 1。
单场在输入编译阶段规范化为 `pasr_algebraic_v1`，不分配随机场或执行
随机数、IEM 和 TCR 统计。2/4/…/64 场仍规范化为 `esf_tpdf`。
`model: esf_tpdf` 配合 `fields: 1` 也遵循此规则。
显式的 `finite_rate_mean` 和原有 PaSR 配置保持各自方法。

```json
"reaction": {
  "model": "auto",
  "representation": "direct_cantera",
  "mechanism_file": "mechanism.yaml",
  "mechanism_sha256": "<机理文件的 SHA256>",
  "phase": "<机理中的 phase>",
  "chemistry_solver": {
    "relative_tolerance": 1e-8,
    "absolute_tolerance": 1e-14,
    "maximum_internal_steps": 5000
  },
  "ensemble": {"fields": 1},
  "mixing": {"c_z": 1.0, "turbulent_schmidt": 0.7}
}
```

`auto` 的 ensemble 和 mixing 均可省略，混合参数缺省为上述值。
ensemble 的 seed、initial_species_offsets、tcr 缺省为 0、空数组和 off。
单场拒绝非空偏移或启用 TCR，避免悄然忽略物理配置。
等价输入共享编译指纹与 Restart 方法身份；改变物理模型仍须重新初始化
其源项历史，不能把平均有限速率模型的 Restart 直接冒充 PaSR。

## 方程与参考

采用已存在的代数时间尺度 PaSR：

\[
\kappa=\frac{\tau_c}{\tau_c+\tau_m},\qquad
\overline{\dot\omega_s}=\kappa\dot\omega_s.
\]

κ 关系可见 Iavarone 等的
[An a priori assessment of the Partially Stirred Reactor (PaSR) model for MILD combustion](https://www.sciencedirect.com/science/article/pii/S1540748920303266)。
混合时间尺度是模型选择的一部分，不能由 κ 的公式唯一确定；
本实现保持现有 LES 标量耗散模型：

\[
\tau_m=\frac{C_z\Delta^2}{2(D+\nu_t/Sc_t)},\quad
D=\sum_sY_sD_s,\quad \Delta=V^{1/3},\quad
\tau_c=\frac{\rho\sum_{\dot\omega_s<0}Y_s}
{\sum_{\dot\omega_s<0}-\dot\omega_s}.
\]

D_s 取真实化学后端的混合平均扩散系数，ν_t 由选定的 SGS 模型给出。
该消耗质量时间尺度是 Hundun 现有 `pasr_algebraic_v1` 的明确选择；
不声称等同于论文中所有化学时间尺度估算或 COAST 的 progress-delta 定义。
所有组分使用同一 κ，保留机理源的质量与元素守恒。

零净反应率直接产生零源。纯组分的混合平均 D 可能严格为零；若 ν_t 也为零，
正 C_z 对应 τ_m=+∞、有限 τ_c 对应 κ=0。这是解析极限，不设置扩散率下限
或化学活跃阈值。近纯组分的非零次正规扩散率可能使 τ_m 超出 FP64；
此时用更宽的中间表达式计算无量纲比例，保留 FP64 可表示的微小 κ，
不会把微量扩散一律当作零。C_z=0 且 D+ν_t/Sc_t=0 的 0/0 仍拒绝。

## CN/BE 分步推进

无喷雾的 `cn_be` 与 `outer_corrected` 组合采用 COAST 的顺序：保存已接受状态，
进行两次组分/焓输运校正，再进行一次反应区积分，最后进行两次动量/压力
校正。两次输运校正共用旧时刻状态与同一 Δt，不是推进两个时间步。
单场复用公共输运工作区和确定性元组，不执行 Wiener 噪声、IEM 或 TCR。

κ 在输运后的 PH/Y 状态计算。Cantera 反应区按 COAST 的摩尔量积分方式，
在区间初始温度和密度下积分 Y/W，结束后按原焓闭合温度；随后组合

\[
Y_s^{n+1}=(1-\kappa)Y_s^{\mathrm{transported}}
             +\kappa Y_s^{\mathrm{reacted}}.
\]

全部组分共享 κ。原输入的质量分数相对、绝对容差保持不变；摩尔变量
绝对容差逐组分换算为 atol/W。零净消耗取 κ=0，不设置温度或燃料浓度
筛选阈值。外部化学提供者须提供有限区间积分及对应身份。
化学增量只计入一次，不再写入下一步的空间算子 EX2 历史。
失败尝试不提交状态；重试从相同旧状态重新输运和积分。

输运、动量和压力线性系统分别检查真实残差；静止初场采用首次压力校正
产生的动量作为本步固定参考尺度。仍计算并报告最终质量、总能量、
组分/元素收支。冻结载体密度的输运方程与后续密度更新产生的分裂项
单独列出，不能将这些已解释项误报为线性求解误差，也不能把阶段方程
通过解释为总体物理守恒验收通过。实际时间步精度和性能比较仍须完成。

新方法使用独立方法指纹与单场元组 Restart 布局。旧平均源布局不能冒充
新方法直接续算；必须明确迁移物理状态并重建历史。BE/PISO 保持原有
瞬时源定义；喷雾交换组合保持原有平均源耦合，其他模型保持各自方法。
新方法当前处于实场验收阶段，
不代表已达到 COAST 吞吐率或长测要求。

回归入口：`v04_reaction_selection` 覆盖默认/显式分派、非法输入、偶数多场
身份及跨输入拼写续算；`v04_pasr_pure` 覆盖真实 Cantera 纯组分初算与恢复；
`v04_models_combustion` 覆盖时间尺度及反应分数极限。

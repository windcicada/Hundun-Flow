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

CN/BE 使用包含压力与动能的公共总能量方程，并按实际组分通量和反应源
独立计算组分/元素账本。每次尝试从已接受状态重建黏度、SGS 与 PaSR 源，整个流动外迭代
复用同一源。失败尝试不提交状态；续算重新建立源历史。
该方法使用冻结的瞬时源，不是每次外迭代执行完整刚性反应器积分。
COAST `c4c1f788c0ca7359653a5acd621dd355bb5c7db3` 的 `fieldpdf.F90`
同样在流动迭代前处理化学，但其先积分再组合状态的完整 PDF 方法与这里的
代数 PaSR 有不同定义。大时间步的稳定性和时间细化仍须实际检验。

回归入口：`v04_reaction_selection` 覆盖默认/显式分派、非法输入、偶数多场
身份及跨输入拼写续算；`v04_pasr_pure` 覆盖真实 Cantera 纯组分初算与恢复；
`v04_models_combustion` 覆盖时间尺度及反应分数极限。

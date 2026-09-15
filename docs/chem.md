# 机理中的低浓度反应级数

对带负级数或分数级数的单向 Arrhenius 反应，可在机理 YAML 的相应
反应记录中配置浓度延拓。模型参数随机理 SHA-256 进入原生资产身份，
新算和 Restart 使用同一份机理文件。

```yaml
hundun-order-regularization:
  model: concentration_c1_v1
  negative-floor: 1e-12 mol/cm^3
  fractional-floor: 1e-14 mol/cm^3
```

上例对应 GTMC JL4 参考源码的默认阈值。浓度统一换算为 kmol/m³。
记反应级数为 a、实际浓度为 C、两个阈值为 e_n、e_f，则所选浓度因子为：

- a < 0：`max(C,e_n)^a`。
- 0 < a < 1 且 C < e_f：
  `C*e_f^(a-1)*(2-a+(a-1)*max(C,0)/e_f)`。
- 其余情况：原幂律因子。

分数级数的正浓度区间在 e_f 处连续匹配函数值和一阶导数；Newton
试探值使用其带符号延拓。各反应沿用原始 SI Arrhenius 常数、化学
计量系数和元素定义。整数幂及普通反应继续由 Cantera 的质量作用项
求值。该配置适用于单向、体相理想气体 Arrhenius 反应；输入检查
核对模型名、参数单位、阈值和反应类别。

Hundun 的浓度相关速率求值器接入每个工作区的 Kinetics，瞬时源项
查询和恒压反应器右端共享同一处理。当前反应器积分沿用稠密数值
Jacobian 的 Cantera/CVODES 后端。各工作区持有独立浓度缓存。

积分端点的组分界限处理采用整体 L1 修正预算
`min(8*absolute_tolerance,64*epsilon_FP64)`。越界幅度和总修正量
同时受此预算约束；超限状态返回显式失败，候选数组保持原值。
随后按处理后的组分增量继续审核质量与元素，焓坐标经公共物性接口
闭合。正常组分状态保持原值。

验证入口为 `v04_models_orders`、`v04_orders_input` 和
`v04_models_cantera_backend`。解析用例覆盖浓度单位、普通幂律区间、
低浓度区间、零浓度、积分解析解、两个半步及工作区复用。GTMC 的
JL4 逐状态对照和原生运行记录见 [算例接入记录](rc.md)。

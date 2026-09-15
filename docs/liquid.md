# 液相物性与绝对焓

液体资产采用 `HUNDUN_LIQUID_ASSET_V1` 和 SI 单位，绑定气相机理、
组分顺序、蒸气组分及焓参考。资产文件的完整字节内容生成材料身份，
该身份随颗粒进入运行状态与 Restart。

密度、比热、潜热、表面张力和黏度支持 `constant` 与 `cubic`。
`cubic` 采用 `theta = T - T_ref`，四个系数依次对应零至三次项。
煤油相关系提供以下固定系数选项：

| 资产字段 | 关系标识 | 输出单位 |
|---|---|---|
| `density` | `kerosene_density_v1` | kg/m³ |
| `cp` | `kerosene_cp_v1` | J/(kg·K) |
| `latent` | `kerosene_latent_v1` | J/kg |
| `saturation` | `kerosene_v1` | Pa |

对应的物性行写为：

```text
density kerosene_density_v1 298.15 0 0 0 0
cp kerosene_cp_v1 298.15 0 0 0 0
latent kerosene_latent_v1 298.15 0 0 0 0
```

`surface_tension` 和 `viscosity` 行随后给出所选关系与系数，
饱和蒸气压行采用 `saturation kerosene_v1`。固定关系的四个用户
系数保持为零；各标识对应表中指定的物理量。

煤油比热、潜热和饱和蒸气压采用亚临界温度区间
`43 K < T < 684.26 K`。算例在该区间内设置材料温度范围；材料
范围与气相机理的适用范围共同约束实际查询。

液相绝对焓由资产中的 `liquid_reference T_ref h_ref` 定义：

\[
h_l(T)=h_{ref}+\int_{T_{ref}}^T c_{p,l}(\theta)\,d\theta.
\]

煤油关系使用解析积分，并保留相邻浮点温度之间的小焓增量。
`h_ref` 按算例的气相形成焓与相变参考确定；颗粒库存、注入库存及
共同交换使用同一液相焓。通用气膜审核同时核对 `h_v = h_l + L`
及 `c_{p,v} = c_{p,l} + dL/dT`。煤油潜热导数使用其真实温度关系。

饱和蒸气压资产返回分段关系值；气膜模型按当地压力执行相平衡
处理。THICK_EX 的气膜查询与冻结物性核已有逐项参考证据，
颗粒推进、共同交换和 624CF 实场接线进度见 [开发记录](rc.md)。

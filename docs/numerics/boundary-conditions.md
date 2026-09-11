# 边界条件

六面名称固定为 `x_min`、`x_max`、`y_min`、`y_max`、`z_min`、`z_max`。流动类型、热边界、标量条件及倒流参数分别指定，合法枚举与组合见[当前 schema](../../versions/v0.4/docs/input-schema.md)。

周期边界必须配对。速度采用全局笛卡尔分量；边界方向、入口/出口权威和倒流处理不能按名称猜测。边界计划在输入类型检查之后另行验证跨面及参数一致性。

IBM 使用独立的静止封闭 STL 配置及流体侧、重构合同。不可穿透标量扩散的切面零通量与几何重构精度是两项要求，参见[离散方法](discretization.md)。不要沿用退休实现的 `velocity_m_per_s` 等 schema 字段。


冷态 `coast_cn_be` 支持质量流量入口与静压出口组合。`pressure_outlet`
给定绝对静压，`allow_backflow=true` 配合回流温度和组分。法向速度随
压力耦合求解，`backflow_velocity` 提供回流切向分量。出流 h/Y 外推，
回流 h/Y 取外界状态，边界 EOS 物性同步进入面通量和方程矩阵。
PISO 路径继续使用其既有速度和物性边界协议。
具体 GTMC 参数、迁移方式和验证记录见 [GTMC 开发记录](../gtmc.md)。

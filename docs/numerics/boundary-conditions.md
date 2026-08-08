# 边界条件

六个外边界 patch 固定命名为 `x_min`、`x_max`、`y_min`、`y_max`、`z_min`、`z_max`。可用类型为：

- `periodic`：必须成对出现；
- `no_slip_wall`：无滑移固壁；
- `symmetry`：法向速度和相应法向梯度受约束；
- `velocity_inlet`：速度为权威输入，热力学与标量字段按配置给定；
- `pressure_outlet`：压力扰动为权威输入，其余量采用出口规则。

schema 3 还可描述静止 STL 浸入壁面。当前只接受 `velocity_m_per_s=[0,0,0]`，焓和标量采用 `zero_normal_diffusive_flux`。`fluid_side` 明确 STL 的哪一侧是计算流体域。

边界方向以计算域外法向为准。配置速度使用全局笛卡尔分量，不随 patch 自动旋转。

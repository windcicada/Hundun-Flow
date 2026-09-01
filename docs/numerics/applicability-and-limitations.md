# 适用范围与限制

## 适合的用途

- Cartesian 六面体网格上的单相不可压/低马赫流动；
- 理想气体压力--焓耦合、温度相关状态与输运；
- 固定或自适应 BDF 时间推进；
- 静止、封闭、非相交 STL 的 immersed boundary；
- resolved-wall WALE 或 Vreman wall-function LES；
- MPI 并行、exact-history Restart、Visit 和 Evidence V6 审计。

IBM `strict_quadratic` 要求每个局部 stencil 能唯一确定三维二次多项式的 10 个系数。复杂曲率、窄缝、薄壁、尖角或特殊网格对齐应优先使用 `adaptive_order`：扩大同侧 donor 搜索后仍不能形成稳定二次基时，局部降为满秩线性；连线性也不能稳定构造时继续 fail-closed。

## V1.0 Re=3900 证据边界

V1.0 候选的物理门是 `20D x 10D x 3D`、约 10% blockage、spanwise periodic、tensor-stretched `456 x 256 x 104` 网格上的中短程测试。它要求至少 `20D/U` 发展、`50D/U` 连续采样和 10 个实际脱涡周期，并比较冻结 Parnaudeau Strouhal 与 15 组 PIV 剖面。

该结果不能被表述为：

- Parnaudeau 4.3% blockage / 20D effective span 的严格几何复现；
- 旧 v0.4 规格中 150+2020D/U、约 420 周期统计已经完成；
- total `Cd` 或 periodic-3D `Cl_rms` 已有一手实验硬门；
- 完整展向域、长程工程统计或普适 LES 精度已经验证。

## 不在当前范围内

- 可压缩激波、声学或高马赫数；
- 反应流、喷雾、颗粒、多相或辐射；
- 移动/变形 IBM、多个相交表面、非流形或未封闭 STL；
- AMR、嵌套网格、非结构网格和 GPU；
- 未经 case-specific preflight 的窄缝燃烧室或任意复杂几何；
- 不经网格、时间步、域和统计窗口敏感性检查的工程结论。

配置通过和单步残差合格都不等于物理设置合理。生产使用仍必须审计边界方向、IBM donor、CFL、正性、质量/能量漂移、limiter 活性和统计收敛。

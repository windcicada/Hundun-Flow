# 数据流

一次正常运行按以下顺序进行：

```text
case.json
  -> rank 0 解析并规范化
  -> MPI 广播并核对资源身份
  -> 建立拓扑、几何、边界和字段布局
  -> 按 profile 构造 IBM / density closure / WALE
  -> 可选读取 Checkpoint v2 或 v3
  -> 时间步试算
  -> 压力/速度修正与守恒检查
  -> 全 rank 接受或回退
  -> 写诊断和 Checkpoint
```

试算期间，新值与已接受值分开保存。任何 rank 报告不可恢复错误时，集体状态会收敛为失败；可重试错误则先回退，再按时间控制规则缩小步长。

浸入边界路径先从 STL 建立表面查询和 Ghost stencil，再由统一边界行参与动量、压力、最终通量和壁面力计算。诊断读取最终权威数据，不另行重算一套近似量。

WALE 每次试算只评估一次。浸入路径固定执行两次 PISO pressure corrector；任一阶段失败时，FlowState、密度闭合和可选模块一起回退，实际执行过的 performance work counter 不回退。

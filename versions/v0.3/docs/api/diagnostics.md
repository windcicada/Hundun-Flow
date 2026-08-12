# Diagnostics API

诊断文件采用 JSON Lines，每行是一条完整记录。公共读取约定是：

1. 按行解析 JSON；
2. 用记录中的描述符或字段名识别含义；
3. 保留未知字段，或安全忽略；
4. 不依赖对象成员顺序和记录顺序；
5. 用文件名中的 rank 与 step 做分组，不把不同步的文件误合并。

常见记录覆盖 MPI 身份、分区、字段布局、Halo、网格、边界、线性求解、最终通量、PISO、时间控制、Checkpoint 和流动驱动状态。并非每种配置都会产生全部记录。

Stage 3 追加 `DiagnosticModuleKind 18--22`：18 为 immersed surface，19 为 ghost stencil，20 为 local flow pattern，21 为 wall force，22 为 WALE。Checkpoint v3 沿用既有 checkpoint kind。记录单位沿用 descriptor：密度 `kg/m3`、动力黏度 `Pa*s`、力矩 `N*m`。failed/retried attempt 不发布 accepted-step module record。

诊断 JSON 适合监控和归档，不保证能恢复求解状态。Restart 必须使用完整检查点目录。

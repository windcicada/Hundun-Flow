# 架构概览

当前唯一实现以 `ApplicationService` 和 `ProductDriver` 组织算例编译与时间推进。普通 CLI 和专用 runner 共用产品内核，runner 另外负责试验统计和分模块观测。

模块划分见[实际源码表](modules.md)。实现约束包括：拓扑与几何分离；字段、工作区和通信请求有明确所有者；最终面通量及边界状态保持版本一致；试算不能直接覆盖已接受历史。

当前生产范围是 CPU/MPI 单相低马赫流动。燃烧、喷雾等规划不构成已集成能力；参见[能力与限制](../releases/current-capabilities.md)。

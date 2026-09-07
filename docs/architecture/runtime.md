# 并行运行时

当前实现由 `core_*`、`parallel_*` 和 `io_*` 管理存储、通信与持久化，公共接口见[模块结构](modules.md)。

字段和工作区在产品构建阶段规划；视图只借用存储，不拥有内存。halo/prepared epoch 与 IBM donor 交换绑定字段身份、分区和有效期，未完成请求不能跨越销毁或状态提交。

可能分配失败的本地阶段先转换为 Status，再由全部 rank 一致决定继续或退出。不能用单 rank 提前返回代替 collective 失败处理，也不能为减少通信次数删去生命周期和 epoch 关闭检查。

Restart 和 Visit 有独立发布及失败语义。求解提交后发生输出失败时，已接受状态不回退成未接受；来源 checkpoint 和历史附件保持只读。

预分配及 `owned_payload_bytes` 不等于完整进程 RSS 上限。[最新内存核查](../verification/2026-09-07-exclusive-module-acceptance.md)说明了实际覆盖的分配、阶段峰值和仍缺失的硬预算。

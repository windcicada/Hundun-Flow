# Restart 格式

当前产品使用 generation 目录和根目录的 `current` 指针发布 checkpoint，不是退休实现的 profile/presence 格式。

当前 V3 保存已接受状态、时间与速率历史、面通量、压力参考和闭域质量目标，并带方法历史签名。读取逐阶段验证身份、布局、文件长度、校验和与全 rank 状态；未验证的 owned image 不能发布为已接受状态。

V1 缺历史、V2 无方法签名、V3 精确历史及显式恢复的区别见[方法合同](../verification/v04-method-history-contract.md)。满足约束时可以跨 MPI 分区恢复；不能把“支持重新分区”理解为物理配置可以任意变化。

checkpoint 的持久化顺序、来源文件只读和一致失败合同不因仓库清理而改变。二进制布局不是外部工具可随意编辑的接口。

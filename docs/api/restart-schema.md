# Restart 格式

当前产品使用 generation 目录和根目录的 `current` 指针发布 checkpoint，不是退休实现的 profile/presence 格式。

当前 V3 保存已接受状态、时间与速率历史、面通量、压力参考和闭域质量目标，并带方法历史签名。读取逐阶段验证身份、布局、文件长度、校验和与全 rank 状态；未验证的 owned image 不能发布为已接受状态。

V1 缺历史、V2 无方法签名、V3 精确历史及显式恢复的区别见[方法合同](../verification/v04-method-history-contract.md)。满足约束时可以跨 MPI 分区恢复；不能把“支持重新分区”理解为物理配置可以任意变化。

checkpoint 的持久化顺序、来源文件只读和一致失败合同不因仓库清理而改变。二进制布局不是外部工具可随意编辑的接口。

## 读取容量与失败边界

`RestartReader::load` 的最后一个 C++ 参数为 `RestartReadLimits`：

| 容量 | 默认上限 | 约束范围 |
| --- | ---: | --- |
| `maximum_bulk_bytes` | 1 GiB / rank | 旧 out、新 owned image、解析块、原始读取块、清单和覆盖数组的同时存活上界 |
| `maximum_manifest_bytes` | 8 MiB | 原始 manifest 文件及其广播，读取/分配前检查 |
| `maximum_source_ranks` | 65536 | 来源 rank 目录，分配目录前检查 |

`current` 另外限制为 256 B。所有文件必须是普通文件；rank 文件在分配前要求
实际长度等于清单及字段布局推得的长度。完整性扫描、覆盖检查和方法历史检查仍保留。

容量可以 rank-local，不要求接收器指针一致；任一 rank 超限时全体失败，原 out 和
源文件不变。零值非法。嵌入调用方可显式提高有限容量；目前普通 CLI 和圆柱 runner
使用上述默认值，尚无命令行读取预算选项。不要以改小清单长度绕过限制。

`RestartReadReport` 的 `retained_image_bytes`、`new_image_bytes`、`peak_bulk_bytes`
分别记录旧 image 容量、新载荷和读取 bulk 上界。它们不包含小路径字符串、分配器及
MPI 内部资源，也不包括随后 ProductDriver 的物性恢复暂存、arena 外数组和 halo，
不是整个进程的 RSS 上限。普通日志关闭成功也不等于 checkpoint 的 fsync 持久化。

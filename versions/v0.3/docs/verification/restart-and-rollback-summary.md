# Restart 与回退摘要

Checkpoint v2 和 Checkpoint v3 已覆盖写出、读取、完整性校验、事务发布和失败回退。读取前会核对完成标记、manifest、文件大小、CRC、语义指纹、网格、边界、字段布局、分区和 optional module presence；检查失败不会部分覆盖正式状态。

时间推进路径还覆盖了局部异常收敛为 collective failure，以及失败试算回退到最近已接受状态。当前 Restart 只接受相同 rank 数和相同分区。

Checkpoint v3 continuation 已覆盖九个 profiles 的同分区恢复。rank-changing Restart 仍不支持。不要手工编辑检查点，也不要只复制其中一部分文件。

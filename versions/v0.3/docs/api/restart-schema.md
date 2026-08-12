# Restart 格式

Restart 是内部二进制事务格式，不是长期稳定的数据交换格式。产品中存在三条兼容路径：被动标量使用 Restart v1，schema 2 流动使用 Checkpoint v2，schema 3 profiles 使用 Checkpoint v3。

Checkpoint v2 目录包含：

- 一个描述全局状态和每个 rank 文件的 manifest；
- 每个 rank 的二进制 payload；
- 所有文件成功发布后写出的完成标记。

读取时检查 schema、文件大小、EOF、CRC、语义指纹、网格、边界、字段布局和分区。任何一项失败都会拒绝整次恢复，正式状态保持不变。

该格式当前要求相同 rank 数和相同分区。文件名、字节布局和内部枚举不是公共 API；外部工具不应修改 payload，也不应根据未文档化偏移读取数据。

Checkpoint v3 presence 1--9 精确对应九个 density/IBM/WALE profiles。manifest 记录 numerical config、geometry、optional module 和 flow state 的身份；恢复时必须先构造 presence 要求的对象。Checkpoint v3 仍要求相同 rank 数和 process grid，rank-changing Restart 不受支持。

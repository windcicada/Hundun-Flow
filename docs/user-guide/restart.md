# Restart

流动求解路径使用 Checkpoint v2 保存时间控制状态、场数据、网格与边界指纹和完整性信息。写出采用临时文件、发布和完成标记三个阶段；只有完成标记存在且清单、文件大小、CRC 和语义指纹都匹配时，目录才可读取。

配置项：

- `restart.read`：是否从已有目录读取；
- `restart.read_directory`：读取目录，`read=true` 时必填；
- `restart.write_directory`：新检查点的根目录；
- `restart.write_interval`：写出间隔。

Restart 要求当前配置、字段 schema、网格、边界和 MPI 分区与写出时兼容。当前不支持改变 rank 数后继续计算。读取失败时，程序不会把部分数据提交为新状态。

请把完成标记、manifest 和所有 rank 文件作为一个整体复制。不要编辑其中任何文件。

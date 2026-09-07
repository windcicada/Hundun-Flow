# 外部自动化操作指南

当前公共操作按[CLI](../api/cli.md)执行：

1. 记录源码、构建身份、输入哈希和 MPI 版本。
2. 按[当前 schema](../api/configuration-schema.md)准备 case root。
3. 使用 `hundun validate <case-dir> --dry-plan` 检查。
4. 确认资源允许后以明确 rank 数和独立输出目录运行。
5. 区分数值终止、输出失败、外部终止和未完成日志。
6. 只从校验通过的 checkpoint 续算；未知方法历史不默认兼容。
7. 保留原始证据；完整归因必须覆盖预期 ranks、步范围及来源身份。

不得猜测单位、静默修改物性/网格/dt/阈值、把配置校验当作科学验收，或停止其他作业。旧 CLI 与退休源码 schema 不属于当前接口。

# 模型与配置准备

外部工具修改 JSON 时应遵守以下顺序：

1. 选择与能力匹配的 `schema_version`；
2. 明确网格尺寸、原点、长度和 SI 单位；
3. 选择密度模型并补齐对应热力学字段；
4. 为六个外边界各写一条记录；
5. 明确 rank 数和进程网格；
6. 为 Restart、diagnostics 和性能输出使用互不冲突的目录；
7. 运行校验和规范化输出。

schema 3 的浸入边界只使用封闭、静止 STL，壁面速度必须严格为零。当前 `0.2.0 candidate` 已接入 schema 3 driver、Checkpoint v3、组合诊断和 `les.model=wale`。自动化工具仍不得因为 `--validate` 成功就直接启动计算；先把输入映射到 capability ledger 的 profile-1 至 profile-9，并核对 rank、process grid、STL 路径和适用范围。

工具必须把用户给定的物性、边界和误差标准视为受控输入。发现矛盾时应报告，不得靠经验值悄悄替换。

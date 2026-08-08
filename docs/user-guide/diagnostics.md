# 诊断输出

`diagnostics.directory` 指定输出目录，`diagnostics.write_interval` 指定时间步间隔，`diagnostics.write_mesh` 控制是否写出网格信息。

结构化诊断按 rank 和时间步写成 JSON Lines 文件，文件名形如：

```text
diagnostics.v1.rank-000000.step-00000000000000000001.jsonl
```

每行是一条独立 JSON 记录。记录包含描述符、请求范围和对应数据；字段集合由当前计算分支决定。处理程序应按字段名读取，忽略自己不认识的新增字段，不要依赖行顺序。

诊断输出用于观察状态，不是 Restart。要继续计算，必须使用完整的 Checkpoint 目录。

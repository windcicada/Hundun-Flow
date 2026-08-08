# 配置说明

计算输入是一个 JSON 文件。程序按 `schema_version` 选择解析规则，并拒绝未知 key、错误类型、重复成员和不满足约束的组合。

使用下面两条命令可以在不推进时间步的情况下检查输入：

```sh
hundun case.json --validate
hundun case.json --print-resolved
```

`--print-resolved` 输出经过规范化的配置，适合存档和比对。路径仍应按原始 `case.json` 的位置理解。

配置中的物理量采用 SI 单位，单位直接写在 key 中，例如 `length_m`、`initial_dt_s`、`dynamic_viscosity_pa_s`。不要依靠隐含单位换算。

完整字段表见[配置 schema](../api/configuration-schema.md)。

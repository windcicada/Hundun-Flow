# 配置说明

输入是 case root 中的 `case.json`。当前 `schema_version=1`；根对象包含网格、流动、求解器、物性、输运标量、六面边界、离散格式与时间控制。字段遵循关闭式校验，未知或重复 key、错误类型和不合法组合被拒绝。

```sh
mpirun -n 1 hundun validate /path/to/case --dry-plan
```

JSON 是配置权威；`.d` 只承载类型化表格，STL 承载几何。引用文件必须满足直接子文件、唯一普通文件及扩展名等约束，不支持路径逃逸。物理量采用 SI，不能假设旧式带单位后缀的 key 仍然可用。

完整字段与枚举见[输入 schema](../../versions/v0.4/docs/input-schema.md)，可运行模板见[最小算例](../../examples/minimal/README.md)。均匀初场由[CLI](../api/cli.md)显式指定，不再依靠最后遍历的边界。

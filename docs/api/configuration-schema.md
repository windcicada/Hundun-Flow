# 配置 schema

当前输入规范在[源码随附的 input schema](../../versions/v0.4/docs/input-schema.md)，以 `CaseCompiler`、后续计划编译和实际测试为最终实现依据。

当前 JSON 的 `schema_version` 为 1。产品版本、JSON schema、内部 typed-model wire 版本、Restart 版本和方法历史签名是不同概念，不能互换。退休实现的 schema 2/3、九个 profile 和 presence 编号不适用于本程序。

使用 `hundun validate <case-dir> --dry-plan` 校验，或从[当前测试模板](../../examples/minimal/README.md)开始。未知、重复、缺失字段及跨字段不合法组合不会被静默兼容。

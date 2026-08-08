# 已接受能力

本页记录当前公开投影的能力边界。对应的已接受源码版本为：

```text
source_commit=66080e324089599711fdb26082af9b330bfdb5ce
source_tree=ab071a61f00eba9ec973beb0fe600066a33ef74f
source_diff_sha256=e98f97ffae590778ac6ff37c4220fb843a4f168c7c8afb9c82ebbc85a25ca2dd
evidence_record_sha256=ba9c32eacf01380a9b656cf81a645231d44d5ff78a267f48761f440359a4c307
```

已接受范围：

- C++17、MPI 3、CPU 参考实现的独立构建和链接；
- 结构化均匀盒与解析扭曲盒的基础网格路径；
- 常密度、物质密度和理想气体密度闭合；
- 固定步长和带集体回退的自适应时间控制；
- 五类外边界条件；
- 单个封闭、静止 STL 的 Local Flow Pattern/Ghost-Cell 浸入边界；
- 两次 PISO pressure corrector；
- 压力、算子、最终通量和壁面力的统一权威链；
- signed-force 四字段语义；
- 1/2/4 rank 小规模分解一致性；
- Checkpoint v2、rollback、collective failure 和结构化 diagnostics。

这里的浸入边界条目是数值核心和公开 C++ 接口的接受结论，不表示 0.1.0 已经提供完整的 schema 3 命令行 driver。Checkpoint v2 和现有 diagnostics 适用于已集成的 schema 2 流动路径；IBM 专用 driver、Checkpoint v3 和组合诊断仍未发布。

未列出的能力不能从本页推定。具体限制见[适用范围与限制](../numerics/applicability-and-limitations.md)。

# 已接受能力

本页记录 `0.2.0 candidate` 的公开边界。精确 task commit、test 和 owner 以
[Stage 3 capability ledger](../numerics/stage3-capability-ledger.md)为准；冻结候选 HEAD 与
formal evidence 只写入治理验收报告，不在候选冻结后回写本页。

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
- Checkpoint v2/v3、rollback、collective failure 和结构化 diagnostics；
- profile-1 至 profile-9 的同一 executable dispatch；
- performance artifact schema v2 的 17 个 exact algorithmic work counters。

这些条目表示产品源码已接通 schema 3 driver。发布结论仍需 S3-V1 对同一个 clean HEAD 执行 formal scientific/performance rows。

未列出的能力不能从本页推定。具体限制见[适用范围与限制](../numerics/applicability-and-limitations.md)。

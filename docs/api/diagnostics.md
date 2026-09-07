# 诊断接口

普通应用的 Evidence 为版本化 JSON Lines，格式见[Evidence schema](../../versions/v0.4/docs/evidence-schema.md)。读取按字段名和 schema 分派，不依赖 JSON 成员顺序，不沿用退休实现的模块编号或文件命名规则。

专用 runner 的观测 CSV 由 `tools/v04_solver_observe.py` 校验。其 V3 完整性依赖冻结元数据中的预期 ranks、步范围与来源身份；重复行、非法计数/时间、缺 rank/loop 和截断尾步不能当作完整归因。

字段可用性、单位、归一化与采样窗口必须一起解释。诊断不替代 checkpoint，也不证明实验吻合，参见[诊断使用说明](../user-guide/diagnostics.md)。

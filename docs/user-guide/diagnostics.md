# 诊断输出

普通应用通过 `--output` 选择 run root，并以 `--output-interval`、`--restart-interval` 控制场输出和 checkpoint 周期，见[CLI](../api/cli.md)。Evidence 不会随这两项置零而关闭。

开发测试可显式启用 `--diagnostics-interval 1`，在不输出 Visit 的情况下逐提交步写入 `diagnostics.jsonl`。默认 0 不改变既有输出。每行具有 step/time/plan 和 payload；payload schema 为 `HUNDUN_V04_DEVELOPMENT_DIAGNOSTICS_V1`，包含 BDF 阶数、实际 dt、最终面通量 revision、质量/内能/动能、边界质量与能量通量、BDF 收支缺陷及累计缺陷，字段后缀给出 SI 单位。报告复用已提交步审计，不再扫描场或增加归约；无效或非有限报告拒绝写出。

边界通量向外为正，热量/黏性功输入向内为正。重启开始新的 `epoch_start_step`；累计量必须在同一 epoch 内解释。此开发台账明确保留 `statistics_eligible=false`，不改变 Evidence schema 或统计准入资格，也不能代替求解验收条件。中测的“500 步发展 + 2500 步观测”是后处理窗口，不是统计收敛声明。

专用圆柱 runner 另有 health、force、conservation、performance 和 solver-loop CSV，以及冻结运行元数据。完整性能归因必须按预期 rank 集合、完整步范围和来源身份校验，不得让两份同时缺行的日志互相证明完整。

求解残差、失败阶段、attempt/dt retry、组成 sweep、输出耗时和整步耗时须分开解释。末尾写入中的行或失败运行只能作为标记清楚的 partial 诊断，不能冒充完整统计。

Visit 用于可视化，Evidence/CSV 用于观测；恢复计算必须使用完整 Restart。历史记录不应因新程序或新脚本而覆写。格式见[Evidence schema](../../versions/v0.4/docs/evidence-schema.md)。

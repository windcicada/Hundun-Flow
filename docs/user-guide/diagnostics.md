# 诊断输出

普通应用通过 `--output` 选择 run root，并以 `--output-interval`、`--restart-interval` 控制场输出和 checkpoint 周期，见[CLI](../api/cli.md)。Evidence 不会随这两项置零而关闭。

专用圆柱 runner 另有 health、force、conservation、performance 和 solver-loop CSV，以及冻结运行元数据。完整性能归因必须按预期 rank 集合、完整步范围和来源身份校验，不得让两份同时缺行的日志互相证明完整。

求解残差、失败阶段、attempt/dt retry、组成 sweep、输出耗时和整步耗时须分开解释。末尾写入中的行或失败运行只能作为标记清楚的 partial 诊断，不能冒充完整统计。

Visit 用于可视化，Evidence/CSV 用于观测；恢复计算必须使用完整 Restart。历史记录不应因新程序或新脚本而覆写。格式见[Evidence schema](../../versions/v0.4/docs/evidence-schema.md)。

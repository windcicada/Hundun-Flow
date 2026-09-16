# ESF 随机场配置

原生 ESF/TPDF 使用成对 Wiener 随机增量。`reaction.ensemble.fields`
指定偶数场数；本轮生产入口验证覆盖 2、4、8、16 场。示例配置片段：

```json
"ensemble": {
  "fields": 16,
  "seed": 1234,
  "initial_species_offsets": [],
  "tcr": {"mode": "off"}
}
```

空偏移数组将各场初始化为给定组成；显式数组按“场、独立组分”排列，
长度为 `fields × 独立组分数`，每个组分的跨场偏移和为零。
随机数地址包含种子、已接受步、阶段、场对和空间方向，同一场对
沿用现有 Philox 地址规则。重算复用同一地址，原生 Restart 保留
各场、独立 field0、方法历史及公共输运缓存。

随机场统计均值承载物理组分与焓，field0 承担压力耦合的降噪状态。
空间输运、隐式 IEM、逐场化学和公共气液源项沿用当前调度与守恒定义。
Vreman 继续作为已选定的 LES 模型。

场数准入同时核对偶数性、原生 Restart 的 64 个主字段预算及
`mesh.limits.max_memory_bytes_per_rank`。主字段预算包含速度、压力、
焓、标量、随机场、field0 和输运缓存；实际可用场数随组分及资源配置
确定。使用 `hundun check <case>` 检查具体算例。

工作区按实际场数提前分配，包括 Wiener 增量、化学密度、压力与
统计缓存。新增缓存进入模型内存计量，整步及恢复的 halo 视图容量
跟随实际交换字段数。两场封存算例的复算保持全部比较值、颗粒记录
与面通量逐位一致。

验证入口为 `v04_kerosene_esf_8/16`、`v04_esf_mixture_mpi_1/2/4`
和 `v04_spray_thick_esf16_cli`。详细范围及证据见
[算例接入记录 S24](rc.md)。
本节点覆盖 TCR off 下的真实煤油化学、空间输运、静态 IBM、
THICK_EX、共同气液交换和原生恢复；实际 GTMC／624CF 的长期
验收沿各自实场计划开展。

# Restart 与方法恢复

```sh
hundun run /path/to/case --restart /path/to/old-run/Restart \
  --output /path/to/new-run --steps 10 --output-interval 0 --restart-interval 10
```

输出使用独立目录。来源根目录的 `current`、对应 generation 的 manifest、全部载荷及统计附件共同组成续算资产，按原始内容保存。

## 精确续算

V3 原生检查点携带方法历史签名；V4／V5 还承载固定／可变长度模型记录。默认入口核对物理 plan、schema、geometry、方法历史及文件完整性，恢复已接受场、前一时刻场、方程速率、面质量通量和质量目标。

跨分区读取遵循网格与计划约束；续算一致性以相同起点和有效时间历史进行比较。生产时间格式采用 `cn_be` 或 `backward_euler`。

## 化学精度细化

```sh
hundun run /path/to/fine --restart /path/to/old-run/Restart \
  --restart-refine-chemistry /path/to/source \
  --output /path/to/fine-run --steps 10 --output-interval 0 --restart-interval 10
```

`source` 为写出检查点的原算例。`fine` 保持全部物理输入及资产内容，调整 `reaction.chemistry_solver` 中的误差控制：相对、绝对容差保持或收紧，其中至少一项收紧；`maximum_internal_steps` 保持或增加。入口适用于真实 Cantera 反应后端，包括具有明确模型标识的煤油关系。

该入口使用 V3 及以上原生历史、严格存储布局和相同方法签名。当前场、历史场、速率及面通量直接恢复，后续步骤采用目标化学精度。来源算例身份写入 Evidence 的 `run_start.history.chemistry_source_case`，策略为 `refine_chemistry`。目标算例写出的新检查点继续使用普通精确续算入口。

## 显式方法恢复

```sh
hundun run /path/to/case --restart /path/to/old-run/Restart \
  --restart-method-recovery --output /path/to/new-run \
  --steps 10 --output-interval 0 --restart-interval 10
```

该策略按当前方法重建方程速率并保存来源闭域质量目标，执行对应的 BE 恢复步骤后继续配置的生产时间格式。V1 的缺省历史重建、V2 的方法签名检查与主动方法恢复分别记录。

圆柱 runner 使用 `--restart-root`。主动方法恢复建立新的统计 epoch；同方法精确续算可继承有效统计。物理统计窗口另行检查流场发展状态。

## 已登记存储迁移

`--restart-storage-compatibility mg-bundle-ghost-v1` 处理已登记的 MG 外层 ghost 布局，采用对应来源身份及存储布局检查，独立于主动方法恢复和化学精度细化策略。

签名组件、迁移表及链式回归见[方法历史合同](../verification/v04-method-history-contract.md)；文件提交规则见[Restart 格式](../api/restart-schema.md)。

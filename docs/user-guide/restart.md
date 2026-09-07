# Restart 与方法恢复

```sh
hundun run /path/to/case --restart /path/to/old-run/Restart \
  --output /path/to/new-run --steps 10 --output-interval 0 --restart-interval 10
```

输出必须使用独立目录。保留来源根目录的 `current`、所指 generation 的 manifest 和全部载荷、统计附件；不要手工改写文件或校验哈希。

## 精确续算

当前 V3 checkpoint 带方法历史签名；默认要求物理 plan、schema、geometry 和历史含义兼容，并验证文件完整性。有效状态、速率、质量目标和 BDF 历史按合同恢复。源格式版本不等于方法版本，也不等于产品版本。

支持满足网格/计划约束的跨分区读取，不承诺任意 rank 数都可用。重启一致性必须在相同有效时间历史下比较。

## 显式方法恢复

```sh
hundun run /path/to/case --restart /path/to/old-run/Restart \
  --restart-method-recovery --output /path/to/new-run \
  --steps 10 --output-interval 0 --restart-interval 10
```

该策略重新构造当前方法速率，保留合法来源的闭域质量目标，先执行必要的 BE 恢复，再恢复 BDF2。完整 V2 无方法签名时不能默认当作同方法历史；V1 确实缺历史的处理与主动恢复分开。

圆柱 runner 使用 `--restart-root`；主动方法恢复会建立新统计 epoch，不混入旧 accumulator 样本，一个 BE 步不表示流场已充分发展。同方法精确续算可以继承有效统计。

## 已登记存储迁移

`--restart-storage-compatibility mg-bundle-ghost-v1` 只处理已登记的 MG 外层 ghost 布局差异，不能与主动方法恢复一起使用。未知计划、物性、边界、几何或损坏载荷仍拒绝，不得靠修改 manifest 绕过。

签名组件、已知迁移表及链式回归见[方法历史合同](../verification/v04-method-history-contract.md)；文件提交规则见[Restart 格式](../api/restart-schema.md)。

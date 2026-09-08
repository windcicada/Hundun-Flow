# 命令行接口

```text
hundun --version
hundun validate <case-dir> [--dry-plan]
hundun init-case --output <case-dir>
hundun run <case-dir> --output <run-dir> --steps <N>
  [--output-interval <N>] [--restart-interval <N>]
  [--diagnostics-interval <N>]
  [--initial-state p,T,Ux,Uy,Uz[,q...]]
  [--restart <restart-dir>] [--restart-method-recovery]
  [--restart-storage-compatibility mg-bundle-ghost-v1]
```

`case-dir` 包含 `case.json` 及其直接引用的数据文件。`validate` 不推进时间步；`init-case` 创建模板。旧式 `hundun case.json`、`--validate` 和 `--print-resolved` 不属于当前 CLI。

均匀初场使用 Pa、K、m/s；q 按冻结标量目录排列。初场和 restart 互斥。方法恢复必须显式同时给出 restart，不能与存储迁移混用，见[Restart](../user-guide/restart.md)。

普通应用的 `--output-interval 0` 关闭 Visit/screen/monitor，`--restart-interval 0` 关闭 checkpoint；Evidence 仍启用。两项周期缺省均为 1，正式运行应明确指定。成功退出为 0，失败为非零；保留完整标准错误和结构化失败报告。

`--diagnostics-interval N` 独立输出已提交步的 `diagnostics.jsonl`，默认 0（关闭）；非零时按全局步号取样，并包含本次运行的最后一步。不需要启用 Visit，详见[诊断输出](../user-guide/diagnostics.md)。

圆柱专用 `v04_thin_domain_runner` 使用自己的 `--spec`、`--case-root`、`--run-root`、`--restart-root`、`--visit-interval` 参数，不应与普通应用参数混用。

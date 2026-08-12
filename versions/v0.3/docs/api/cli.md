# 命令行接口

```text
hundun <case.json> [--validate|--print-resolved] | hundun --version
```

| 形式 | 行为 |
| --- | --- |
| `hundun --version` | 输出程序版本后退出 |
| `hundun case.json` | 读取配置并运行 |
| `hundun case.json --validate` | 完整解析和校验配置，成功时输出 `VALID` |
| `hundun case.json --print-resolved` | 输出规范化 JSON，不推进时间步 |

成功返回 `0`，配置、I/O、MPI 或求解失败返回非零值。错误消息写到标准错误；正常信息写到标准输出。MPI 运行时，根 rank 负责用户可见的主要错误信息，失败状态会在 communicator 内统一。

CLI 不支持把配置从标准输入传入，也不接受未列出的额外参数。

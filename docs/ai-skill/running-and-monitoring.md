# 启动与监控

启动前记录工作目录、配置及程序 SHA-256、build manifest、rank 数、命令、环境和日志路径，确认没有超额占用计算资源。

```sh
mpirun -n 4 /absolute/path/to/hundun run /absolute/path/to/case \
  --output /absolute/path/to/new-run --steps 10 \
  --output-interval 10 --restart-interval 10
```

按 PID、工作目录和命令确认目标进程，再读取 health/Evidence、标准错误和 checkpoint 更新情况。尚未写完的尾行不作为完整接受记录，暂时无新行不自动判定为死锁。

停止、替换或重启冻结作业须有相应授权；不得按程序名批量杀进程，也不得让两个作业共写目录。专用圆柱 runner 的参数和普通应用不同，见[CLI](../api/cli.md)。

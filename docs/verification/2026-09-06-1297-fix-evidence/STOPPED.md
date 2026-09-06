# 重放人工中断记录

用户要求改为方法优先审查，不先重现故障。因此停止本目录 replay-baseline.sh 刚启动的重放。

- 操作对象：经命令行核实的 mpirun PID 161801；发送 SIGTERM。
- 运行目录：/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52/replay-baseline-1000to1296-20260906。
- 源 checkpoint：long-review-749to35000-20260905/Restart；未改动。
- 启动会话 64943 已结束，外层命令退出码 1；time 记录 491.28 秒。
- 停止后检查未发现该重放对应的 MPI runner 残留进程。
- 日志和可能存在的部分结果保留，不作为数值失败、性能基线、回归通过或完成证据。

退出码对应本次人工中断，不能解释为新出现的数值发散、MPI 程序缺陷或 OOM。

未再次启动重放；后续转入 [方法审查](../2026-09-06-method-first-numerical-audit.md)。

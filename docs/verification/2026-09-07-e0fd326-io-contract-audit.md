# e0fd326 之后的 I/O 与入口合同核查

基线：`e0fd326d3ed6cccd52f812c113712f34a08cec95`，实际工作区
`/home/wyf/code_dev/.worktrees/hundun-flow-strict-coast-parity`。
原工作区在本轮开始时干净。未改物性、网格、时间步、残差门槛、
refinement 容量或 checkpoint 发布/持久化协议。

## 1. 运行器日志结束状态：实际复现并修正

`tools/v04_thin_domain_runner.cpp::run()` 的原结束路径没有检查四个
root 流的 close 状态；performance 和逐 rank 的 solver 日志依赖析构关闭。
真实两 rank CLI 在持久化 step 1 之后，向 root `force.csv` 的 libc
`fclose` 注入 EIO：原程序两个 rank 均返回 0，仍打印 `COMPLETED`。
这是实际注入复现，不是一次真实磁盘故障记录。

最小修改：`complete_logs()` 显式处理所有实际打开的流，逐流 flush/close，
保存首个本地错误，再选出失败 rank 并广播错误。disabled 流不参与。
运行主体提前返回也先经过相同清理边界；原求解/输出失败退出码优先，
不因清理失败变成成功，也不改写已接受物理状态。只有全体关闭成功，
才完成日志只读化与目录同步并打印 `COMPLETED`。

“日志 flush/close 成功”只表示缓冲输出和关闭成功，不声称普通 CSV
已经 fsync 持久化。checkpoint 的文件同步、目录同步和 current 发布合同不变。

回归入口：`v04_runner_log_completion_mpi_2`。39 次独立 CLI 检查通过：
普通/观测开关、六类文件的 EIO/ENOSPC、末次显式刷新、关闭、
先 EIO 后 ENOSPC 的首错保留；含 rank 1 的 solver 文件。
31 个注入失败用例全 rank 返回 6，不打印完成；8 个成功/未启用用例正常完成。
每个用例的持久化 rank checkpoint 字节 SHA-256 与正常基线完全相同。

实际证据目录：
`/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52/module-io-audit-20260907/`。
`close-red-checkout/results.json`、`force-close-eio.log` 记录 RED；
`close-green/results.json` 和每次 CLI 日志记录 GREEN。

## 独立运行与构建观察

本轮开始前，冻结服务 `hundun-re3900-module-reviewed-20260907.service`
已于 15:42:58 CST 退出，状态 failed/exit-code 7；尝试第 7232 步时
返回 `status=5/804 stage=15`，已接受步为 7231，最新持久化 checkpoint 为 7000。
本轮未停止、替换或重启该程序。此退出未经过正常日志完成分支，
不能据此把它归因于上述 close 缺口；其数值根因尚未确定。

从基线 `git archive` 解压的无 `.git` 源码可独立 configure、构建
`hundun`、`v04_thin_domain_runner` 和 `v04_app_case_test`。
但真实 runner 初始执行返回 `candidate_identity_status=2/10505`：
当前运行证据身份要求有效 Git commit/tree，而归档构建记录为 unavailable。
本轮未伪造 Git 身份，也未把旧冻结程序的验收标签移给此构建。
这是与旧入口清理无关的归档运行限制，需要单独定义归档来源身份合同。

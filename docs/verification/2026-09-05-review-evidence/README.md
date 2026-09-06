# 本轮静态审阅实施证据

基线 `b779bff4bbff7067e691952b25c90cbc86be7795`；实现仍为本地未提交修改，
不把基线 SHA 当成本次二进制的完整身份。独立完整源码快照、工作区 patch、
自动构建 manifest、编译命令、输入副本及清单在：

`/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52/review-frozen-20260905/`

冻结二进制 SHA-256：
`3737f83420ce255c0268cf655855284a0690100a2af6456a85e075604952ac38`。
源码 archive：
`source-v04.tar.gz`，其 SHA 和每个源文件/输入 SHA 记录在 `inventory.json`。
这些文件用于重现本次本地候选，不用于改写历史二进制、receipt 或 checkpoint 的身份。

`performance-comparison.csv` 是当前同物理算例比较清单。PISO、低马赫 COAST、
普通可压缩 COAST 没有本轮同条件实测，不填入历史异条件数字，不声称超过 COAST。
单轮短窗只用于诊断与耗时估计，不代表湍流统计收敛或正式性能发布验收。

本轮完整说明见 `../2026-09-05-static-review-implementation.md`。原始测试启动输出
当前保存在 `build-review-fixes/`，通过/失败的独立日志均保留，最终摘要另外收集于此。

长测已于 2026-09-05 23:38:10 +08:00 启动，128 ranks，从 749 步推进至
35,000 步。`launch-long.sh` 是实际启动命令，`launch-status.json` 是启动后的
核验快照，不是完成凭据；持续运行数据位于独立的
`long-review-749to35000-20260905/` 目录。此轮没有提交或推送 GitHub。

# MPI 分解摘要

同一均匀球形产品路径在 `1x1x1`、`2x1x1` 和 `4x1x1` 进程网格上完成小规模分解检查。三种配置均正常退出，没有 collective mismatch、死锁或非有限诊断，严格分解 oracle 通过。

这项证据覆盖 1/2/4 rank 的基础一致性和共享浸入边界行的全局 donor 语义。它不等同于强扩展测试，也不说明更大 rank 数下的性能。

Stage 3 fast matrix 还覆盖九个 profiles、Checkpoint v3 continuation 和 diagnostics 的 1/2 rank，并选取 4 rank decomposition rows。正式 24-cubed 1/2/4 rank 是 S3-V1 gate，不是这里的小规模证据。

对应日志摘要 SHA-256：

```text
4a3669fd59c0ff0e964cb3e0ace9f423469318857b52455effe3215c8c193558
```

# 624CF 气液模型示例

煤油四步反应、2 场 ESF、Vreman、动态 TCR、IBM、THICK_EX 与 SGS 破碎，采用 624CF 原时间步和固定热力学压力。
该目录提供小型原生模型演示。实场网格和检查点沿各自的迁移流程组织。
正式构建方式见 [构建说明](../../docs/user-guide/build.md)。

在仓库根目录运行：

```sh
hundun check examples/cf
mpiexec -n 1 hundun run examples/cf --output cf1 --steps 3 \
  --initial-state 790216.58,900,0,0,0,0.001,0.05,0.005,0.01,0.2,0.02 --restart-interval 1 --diagnostics-interval 1
mpiexec -n 2 hundun run examples/cf --output cf2 --steps 1 \
  --restart cf1/Restart --restart-interval 1 --diagnostics-interval 1
```

[动态 TCR](../../docs/tcr.md) · [压力分工](../../docs/p0.md) · [运行记录](../../docs/rc.md)

# GTMC 气相模型示例

JL4、4 场 ESF、Vreman、动态 TCR，采用 GTMC 原时间步。
该目录提供小型原生模型演示。实场网格和检查点沿各自的迁移流程组织。
正式构建方式见 [构建说明](../../docs/user-guide/build.md)。

在仓库根目录运行：

```sh
hundun check examples/g
mpiexec -n 1 hundun run examples/g --output g1 --steps 3 \
  --initial-state 100000,1100,10,0,0,0.05,0.2,0.01,0.005,0.05,0.001 --restart-interval 1 --diagnostics-interval 1
mpiexec -n 2 hundun run examples/g --output g2 --steps 1 \
  --restart g1/Restart --restart-interval 1 --diagnostics-interval 1
```

[动态 TCR](../../docs/tcr.md) · [压力分工](../../docs/p0.md) · [运行记录](../../docs/rc.md)

实场运行可使用 `bash tools/hot.sh <程序> <算例> <Restart> <输出目录> <步数> 128`。
默认每 100 步输出可视化、每 500 步保存 Restart、每 10 步输出诊断。
`HF_PLOT_EVERY`、`HF_SAVE_EVERY`、`HF_DIAG_EVERY` 分别配置三种频率。

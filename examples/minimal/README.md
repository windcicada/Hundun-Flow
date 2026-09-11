# 最小运行模板

本目录与 `hundun init-case` 生成的输入一致：8×8×8 网格、CN/BE 时间推进、中央动量格式、单组分理想气体、速度入口与静压出口。入口速度为 1 m/s，温度为 300 K，出口压力为 101325 Pa。

在仓库根目录执行：

```sh
mpirun -n 1 build/versions/v0.4/hundun validate examples/minimal
mpirun -n 1 build/versions/v0.4/hundun run examples/minimal \
  --output /tmp/hf-run --steps 10 --output-interval 10 --restart-interval 10
```

每次运行使用独立输出目录。CN/BE 采用 CN 动量与 BE 质量、焓、组分输运；显式 `backward_euler` 用于兼容的 PISO/SIMPLE 配置。重启动保存速度端点、热力学状态、通量及算法历史签名。

[air.d](../air.d) 提供 O2/N2 空气的 NASA 热力学与 Perry 输运数据，O2 独立质量分数为 0.23291751145757963，N2 为补足组分。新算例见 [Re3900](../cyl/README.md)，开发记录见 [cn.md](../../docs/cn.md)。

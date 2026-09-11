# 快速开始

`examples/minimal` 提供与 `hundun init-case` 一致的 CN/BE 运行模板。入口速度为 1 m/s，温度为 300 K，出口静压为 101325 Pa。

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DHUNDUN_BUILD_TESTS=OFF
cmake --build build -j 2 --target hundun
build/versions/v0.4/hundun --version
mpirun -n 1 build/versions/v0.4/hundun validate examples/minimal
mpirun -n 1 build/versions/v0.4/hundun run examples/minimal \
  --output /tmp/hf-run --steps 10 --output-interval 10 --restart-interval 10
```

CLI 接收包含 `case.json` 的算例目录。数据文件遵循[输入规则](../../versions/v0.4/docs/input-schema.md)，输出放入独立目录。初始组分取各入口一致的指定值，也可通过 `--initial-state p,T,Ux,Uy,Uz[,q...]` 显式指定。

默认时间格式为 `cn_be`，耦合标记为 `CN_BE`；动量采用 CN 中央格式，质量、焓和组分采用 BE。计算检查包括物理状态、守恒、IBM 固体速度与实际时间步，统计计算还需确定网格和采样窗口。

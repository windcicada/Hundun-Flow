# 快速开始

最小模板位于 `examples/minimal`，与当前公共应用回归使用的输入一致，只用于检查安装和运行路径。

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DHUNDUN_BUILD_TESTS=OFF
cmake --build build/release -j 2 --target hundun
build/release/versions/v0.4/hundun --version
mpirun -n 1 build/release/versions/v0.4/hundun validate examples/minimal --dry-plan
mpirun -n 1 build/release/versions/v0.4/hundun run examples/minimal \
  --output run-minimal --steps 10 --output-interval 10 --restart-interval 10 \
  --initial-state 101325,300,0.1,0,0
```

CLI 接收包含 `case.json` 的目录，不是旧式的 `hundun case.json --validate`。所有数据文件按[当前输入规则](../../versions/v0.4/docs/input-schema.md)放在 case root。输出使用独立目录，不修改来源文件。

通过配置校验不等于通过精度或稳定性验证。正式计算前还须核对边界、网格、物性、dt、守恒与统计窗口；不要与现有长测并发占满资源。

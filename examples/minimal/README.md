# 最小运行模板

本目录的 JSON 和物性文件与当前回归夹具逐字节一致。网格为 8×8×8，单相低马赫理想气体、PISO，无 IBM、LES 或输运标量；它不是 Re3900 科学算例。

从仓库根目录运行：

```sh
mpirun -n 1 build/release/versions/v0.4/hundun validate examples/minimal --dry-plan
mpirun -n 1 build/release/versions/v0.4/hundun run examples/minimal \
  --output run-minimal --steps 10 --output-interval 10 --restart-interval 10 \
  --initial-state 101325,300,0.1,0,0
```

输出写入 `run-minimal`，不要让不同计算共用该目录。细节见[快速开始](../../docs/user-guide/quick-start.md)。

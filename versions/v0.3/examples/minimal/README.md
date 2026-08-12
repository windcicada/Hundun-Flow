# 最小运行模板

这个目录只展示输入结构、启动命令和输出位置，不是科学验证算例。

```sh
../../build/release/src/hundun case.json --validate
../../build/release/src/hundun case.json --print-resolved
mpiexec -n 1 ../../build/release/src/hundun case.json
```

运行会在本目录生成 `output` 和 `Restart`。再次运行前，请先把旧输出移到其他位置，或在副本中修改输出目录。

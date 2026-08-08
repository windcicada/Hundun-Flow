# 自动化安装与构建

自动化工具应先探测 `cmake`、C++ 编译器、`mpicxx` 和 `mpiexec`，再执行：

```sh
cmake -S . -B build/release \
  -DCMAKE_BUILD_TYPE=Release \
  -DHUNDUN_BUILD_TESTS=OFF
cmake --build build/release -j 2
build/release/src/hundun --version
```

每次构建记录完整命令、退出状态、生成二进制的 SHA-256、编译器版本和 MPI 版本。不要在构建失败后沿用旧二进制，也不要在不同 MPI 实现之间混用编译器包装器和运行器。

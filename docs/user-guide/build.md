# 构建

推荐使用独立构建目录：

```sh
cmake -S . -B build/release \
  -DCMAKE_BUILD_TYPE=Release \
  -DHUNDUN_BUILD_TESTS=OFF
cmake --build build/release -j 2
```

可用选项：

| 选项 | 默认值 | 作用 |
| --- | --- | --- |
| `HUNDUN_BUILD_TESTS` | `OFF` | 产品包不含测试源码，应保持关闭 |
| `HUNDUN_ENABLE_ASAN` | `OFF` | 为调试构建启用 AddressSanitizer |
| `HUNDUN_ENABLE_UBSAN` | `OFF` | 为调试构建启用 UndefinedBehaviorSanitizer |

构建完成后，执行：

```sh
build/release/src/hundun --version
```

不要在同一构建目录中来回切换编译器、MPI 或 sanitizer。需要更换配置时，新建构建目录更可靠。

# 构建

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DHUNDUN_BUILD_TESTS=OFF
cmake --build build/release -j 2 --target hundun
build/release/versions/v0.4/hundun --version
```

也可使用 `cmake --preset release` 和 `cmake --build --preset release -j 2`。仓库只构建当前实现；旧脚本的 `HUNDUN_SOURCE_VERSION=v0.4` 仍可用，其他值明确拒绝。

| 选项 | 默认 | 用途 |
| --- | --- | --- |
| `HUNDUN_BUILD_TESTS` | OFF | 启用随源码提供的回归；测试使用独立 test core |
| `HUNDUN_ENABLE_ASAN` | OFF | AddressSanitizer 检查构建 |
| `HUNDUN_ENABLE_UBSAN` | OFF | UndefinedBehaviorSanitizer 检查构建 |
| `HUNDUN_ENABLE_HYPRE` | OFF | 可选隔离的 HYPRE 适配器，需要可用依赖 |

更换编译器、MPI 或 sanitizer 时使用新的构建目录。性能测量使用不带 sanitizer 的构建，不使用 `-ffast-math` 绕过有限性和浮点语义。不要在已占满核数的长测旁启动编译或回归。

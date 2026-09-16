# 构建

先阅读[跨主机环境与运行](environment.md)。以下最小构建不包含真实
Cantera 化学后端，不能直接运行 `examples/g` 或 `examples/cf`：

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DHUNDUN_BUILD_TESTS=OFF
cmake --build build/release -j 2 --target hundun
build/release/versions/v0.4/hundun --version
```

反应计算使用 `HUNDUN_CANTERA_ROOT=/absolute/path/to/sdk cmake --preset release`
和 `cmake --build --preset release -j 2`。SDK 必须符合锁定的库和许可身份，
不是 Python `pip install cantera` 的安装目录。仓库只构建当前实现；
旧脚本的 `HUNDUN_SOURCE_VERSION=v0.4` 仍可用，其他值明确拒绝。

| 选项 | 默认 | 用途 |
| --- | --- | --- |
| `HUNDUN_BUILD_TESTS` | OFF | 启用随源码提供的回归；测试使用独立 test core |
| `HUNDUN_ENABLE_REACTING_CANTERA` | OFF | 启用真实化学；release preset 为 ON |
| `HUNDUN_CANTERA_PACKAGE_ROOT` | 空 | 已校验的 Cantera C++ SDK 根目录；release preset 从 `HUNDUN_CANTERA_ROOT` 读取 |
| `HUNDUN_ENABLE_ASAN` | OFF | AddressSanitizer 检查构建 |
| `HUNDUN_ENABLE_UBSAN` | OFF | UndefinedBehaviorSanitizer 检查构建 |
| `HUNDUN_ENABLE_HYPRE` | OFF | 可选隔离的 HYPRE 适配器，需要可用依赖 |

更换编译器、MPI 或 sanitizer 时使用新的构建目录。性能测量使用不带 sanitizer 的构建，不使用 `-ffast-math` 绕过有限性和浮点语义。不要在已占满核数的长测旁启动编译或回归。

## 源码归档与身份

Git checkout 和带有效 `.git_archival.txt` 的 GitHub/Git 源码归档均可构建运行。
归档通过 `export-subst` 保存原 commit/tree，实际源码另计算内容摘要。
不要自行填写来源 SHA，也不要把旧程序的身份或验收记录用于新程序。

2026-09-16 已从提交 `8235e59` 的 GitHub 归档独立构建，在另一台
glibc 2.31 主机上通过私有运行库完成 JL4 四场示例的新算和跨进程数恢复。
历史归档限制的修复见[独立运行包](../pkg.md)。

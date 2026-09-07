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

## 源码归档的当前限制

无 `.git` 的 `git archive` 源码可 configure、构建生产/测试目标并执行输入校验。
但普通 CLI 与圆柱 runner 的运行证据还要求有效的 Git commit/tree 身份；归档构建
会记录 unavailable，并在运行身份阶段拒绝推进。当前请从 Git checkout 构建正式运行程序。

归档支持需要单独定义随包来源身份与实际源码内容的验证合同，不能伪造 Git SHA、
取消身份检查或继承另一可执行程序的冻结验收标签。检查记录见
[干净构建与 I/O 核查](../verification/2026-09-07-e0fd326-io-contract-audit.md)。

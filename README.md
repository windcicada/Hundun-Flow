# HUNDUN-FLOW

HUNDUN-FLOW 是一个使用 C++17 实现的低马赫数有限体积求解程序。当前默认源码线是正在独立实现的 v0.4 Cartesian 架构；冻结的 v0.3 产品源码保存在 `versions/v0.3`，可通过同一根构建入口显式选择。

根构建默认选择 `HUNDUN_SOURCE_VERSION=v0.4`，只接受 `v0.4` 或 `v0.3`。两个版本目录各自拥有独立的 CMake project 和依赖发现逻辑，不会混合链接源码。

## 构建

v0.4 骨架需要 CMake 3.21 或更新版本和支持 C++17 的编译器：

```sh
cmake -S . -B build/v04-release \
  -DHUNDUN_SOURCE_VERSION=v0.4 -DCMAKE_BUILD_TYPE=Release
cmake --build build/v04-release -j 2
build/v04-release/versions/v0.4/hundun --version
```

冻结 v0.3 仍保留其原有 MPI 3 和其他构建要求：

```sh
cmake -S . -B build/v03-release \
  -DHUNDUN_SOURCE_VERSION=v0.3 -DCMAKE_BUILD_TYPE=Release
cmake --build build/v03-release -j 2 --target hundun
```

两个源码线也都可以直接将各自的 `versions/v0.x` 目录作为 CMake 源目录配置。构建过程不需要网络，也没有 Python 运行时依赖。

## 运行

v0.4 当前只提供稳定公开基础类型、状态消息和 `--version` 命令；数值求解与运行接口将在后续任务中实现。需要运行既有求解器时请选择冻结的 v0.3 源码线，并阅读其中的用户指南和适用范围说明。

## 文档入口

- [v0.4 Cartesian 架构实施计划](docs/superpowers/plans/2026-08-12-hundun-flow-v0.4-cartesian-performance-architecture.md)
- [v0.4 公开生命周期调研](docs/references/2026-08-13-hundun-v04-public-lifecycle-survey.md)
- [v0.4 目标热循环](docs/architecture/v0.4-target-hot-loop.md)
- [冻结 v0.3 用户指南](versions/v0.3/docs/user-guide/quick-start.md)
- [冻结 v0.3 配置格式](versions/v0.3/docs/api/configuration-schema.md)
- [冻结 v0.3 适用范围与限制](versions/v0.3/docs/numerics/applicability-and-limitations.md)

## 许可证

HUNDUN-FLOW 采用 Apache License 2.0。第三方组件及其许可证见 [THIRD_PARTY.md](THIRD_PARTY.md)。

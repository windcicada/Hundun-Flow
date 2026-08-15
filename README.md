# HUNDUN-FLOW

HUNDUN-FLOW 是一个使用 C++17 实现的低马赫数有限体积求解器，面向带浸入边界（IBM）的不可压 / 低马赫流动问题，支持 MPI 并行分解与 LES 湍流建模。

仓库当前包含两条独立源码线：

| 源码线 | 状态 | 说明 |
| --- | --- | --- |
| `versions/v0.3` | 冻结产品 | 可运行的既有求解器：静止周期 IBM、WALE、两次 PISO 修正、精确 Schur 压力通路、MPI 分解与四字段力报告 |
| `versions/v0.4` | 活跃开发 | 独立重写的 Cartesian 性能架构，当前为默认构建目标 |

两条源码线各自拥有独立的 CMake project 和依赖发现逻辑，通过根构建入口的 `HUNDUN_SOURCE_VERSION` 显式选择，不会混合链接源码。构建过程不需要网络，也没有 Python 运行时依赖。

## v0.4 Cartesian 架构进展

v0.4 按照《[Cartesian 性能架构实施计划](docs/superpowers/plans/2026-08-12-hundun-flow-v0.4-cartesian-performance-architecture.md)》分 22 个任务推进，当前已完成 13 个：

- 公开生命周期调研冻结与 v0.3/v0.4 构建隔离
- 扁平算例输入编译为经过验证的能力模型
- 带填充的 SoA 存储、修订版本与尝试事务
- Cartesian 几何计划与授权 STL 扫描
- CPU/NUMA 放置与持久 Halo 引擎
- 边界、离散格式与时间控制计划编译
- 统一热力学、输运与贡献契约
- 守恒 Cartesian 内核与单一面通量权威
- 阶段执行图编译器
- 四层线性生命周期与 Krylov 求解器（FGMRES 等）
- 原生 Cartesian 多重网格与隔离的 HYPRE Struct 适配器
- 静态 IBM 拓扑、二次插值模板与表面计划编译

后续任务（PISO 修正器、湍流模型、产品驱动、Re=3900 圆柱联合放行门等）见实施计划的 Task 15–22。

## 构建

要求 CMake 3.21 或更新版本，以及支持 C++17 的编译器。

构建默认的 v0.4 源码线：

```sh
cmake -S . -B build/v04-release \
  -DHUNDUN_SOURCE_VERSION=v0.4 -DCMAKE_BUILD_TYPE=Release
cmake --build build/v04-release -j 2
build/v04-release/versions/v0.4/hundun --version
```

构建冻结的 v0.3 产品线（另有 MPI 3 等要求）：

```sh
cmake -S . -B build/v03-release \
  -DHUNDUN_SOURCE_VERSION=v0.3 -DCMAKE_BUILD_TYPE=Release
cmake --build build/v03-release -j 2 --target hundun
```

两个源码线也都可以直接将各自的 `versions/v0.x` 目录作为 CMake 源目录独立配置。

## 运行

v0.4 当前提供稳定公开基础类型、状态消息和 `--version` 命令，数值求解与运行接口正在按任务计划实现。需要运行既有求解器时请选择冻结的 v0.3 源码线，并阅读其[用户指南](versions/v0.3/docs/user-guide/quick-start.md)和[适用范围与限制](versions/v0.3/docs/numerics/applicability-and-limitations.md)。

## 验证状态

能力声明以仓库内的验证文档为准，不作超出证据范围的外推：

- [收敛验证摘要](docs/verification/convergence-summary.md)
- [守恒验证摘要](docs/verification/conservation-summary.md)
- [已接受能力清单](docs/verification/accepted-capabilities.md)
- [Re=3900 圆柱绕流基准设计](versions/v0.3/docs/benchmarks/cylinder-re3900-design.md)

## 文档入口

- [v0.4 Cartesian 架构实施计划](docs/superpowers/plans/2026-08-12-hundun-flow-v0.4-cartesian-performance-architecture.md)
- [v0.4 目标热循环](docs/architecture/v0.4-target-hot-loop.md)
- [v0.4 公开生命周期调研](docs/references/2026-08-13-hundun-v04-public-lifecycle-survey.md)
- [冻结 v0.3 用户指南](versions/v0.3/docs/user-guide/quick-start.md)
- [冻结 v0.3 配置格式](versions/v0.3/docs/api/configuration-schema.md)

## 许可证

HUNDUN-FLOW 采用 Apache License 2.0，见 [LICENSE](LICENSE)。第三方组件及其许可证见 [THIRD_PARTY.md](THIRD_PARTY.md)。

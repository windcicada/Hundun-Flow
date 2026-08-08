# HUNDUN-FLOW

HUNDUN-FLOW 是一个使用 C++17 和 MPI 实现的低马赫数有限体积求解程序。目前 `hundun` 命令行程序可运行结构化网格上的常密度与变密度流动，支持固定或自适应时间推进、同分区 Restart 和结构化诊断输出。源码还提供已验收的静止 STL 浸入边界数值核心；0.1.0 尚未把它接入可运行、可 Restart 的 schema 3 命令行流程。

本仓库是可发布的产品源码。设计记录、测试程序和原始验收日志不随产品分发；公开文档只陈述已经具备且有证据支持的能力。

## 构建

需要 CMake 3.21 或更新版本、支持 C++17 的编译器，以及 MPI 3 实现。

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j 2
build/release/src/hundun --version
```

构建过程不需要网络，也没有 Python 运行时依赖。

## 运行

先检查配置，再启动计算：

```sh
build/release/src/hundun examples/minimal/case.json --validate
mpiexec -n 1 build/release/src/hundun examples/minimal/case.json
```

`examples/minimal` 只是输入和运行方式示例，不是科学验证算例。开始正式计算前，请阅读[适用范围与限制](docs/numerics/applicability-and-limitations.md)。

## 文档入口

- [安装、配置和运行](docs/user-guide/quick-start.md)
- [配置格式](docs/api/configuration-schema.md)
- [数值方法](docs/numerics/discretization.md)
- [当前能力与验证边界](docs/verification/accepted-capabilities.md)
- [面向外部自动化工具的操作说明](docs/ai-skill/index.md)
- [源码目录与命名规范](docs/development/naming-and-style.md)

## 许可证

HUNDUN-FLOW 采用 Apache License 2.0。第三方组件及其许可证见 [THIRD_PARTY.md](THIRD_PARTY.md)。

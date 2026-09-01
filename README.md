# HUNDUN-FLOW

HUNDUN-FLOW 是一个 C++17/MPI 低马赫数有限体积求解器，面向 Cartesian 网格上的单相流、LES 和静止浸入边界（IBM）。当前默认源码线位于 `versions/v0.4`；这里的 `v0.4` 是内部源码/API 命名空间，公开产品版本为 `1.0.0`。

当前分支是 V1.0 发布候选。只有预注册的 Re=3900 中短程产品/物理门、候选身份和发布 provenance 全部通过后，才会创建 `v1.0.0` tag；候选失败时不会发布。

## 主要能力

- tensor-stretched Cartesian 网格与 MPI Cartesian 分解；
- 单相低马赫理想气体热力学、温度/密度/输运同步更新；
- 两次 PISO pressure corrector，以及同目标时间层的 continuity--energy pressure--enthalpy nonlinear refinement；
- conservative final face mass-flux authority、EOS/continuity/energy/closed-mass/gauge 独立终端门；
- common-face owner AFC、provisional/committed CFL 证书和事务化 retry/rollback；
- 静止、封闭 STL IBM，支持严格三维二次和可审计的局部自适应降阶；
- WALE 与 Vreman wall-function 路径；
- exact-history MPI Restart、Visit 输出和 Evidence V6 provenance。

完整边界见[当前版本能力](docs/releases/current-capabilities.md)、[已接受能力](docs/verification/accepted-capabilities.md)和[适用范围与限制](docs/numerics/applicability-and-limitations.md)。

## 构建

要求 CMake 3.21、C++17 编译器和 MPI 3。

```sh
cmake -S . -B build/release \
  -DHUNDUN_SOURCE_VERSION=v0.4 \
  -DHUNDUN_BUILD_TESTS=OFF \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build/release -j 2 --target hundun
build/release/versions/v0.4/hundun --version
```

冻结的旧 `versions/v0.3` 产品线仍可通过 `-DHUNDUN_SOURCE_VERSION=v0.3` 独立构建；两条源码线不会混合链接。

## 运行

先验证算例，再运行：

```sh
mpirun -np 4 build/release/versions/v0.4/hundun validate case --dry-plan
mpirun -np 4 build/release/versions/v0.4/hundun run case \
  --output run --steps 10 --output-interval 0 --restart-interval 10
```

128-rank Re=3900 发布候选采用每 rank 一个物理核心，任何单个 rank 的 owned-cell 数不超过 100,000。精确门槛见 [`HUNDUN_V1_RE3900_MEDIUM_RELEASE_POLICY_V2`](docs/verification/v1.0-re3900-medium-release-policy.json)。该门只授权带明确 10% blockage/薄展向域限制的 V1.0 软件发布，不声称已完成旧 v0.4 规格中的 420-cycle 文献统计、完整展向域验证或所有 CFD 应用验证。

## 验证与文档

- [V04-2 数值产品闭环](docs/verification/2026-08-30-v04-2-pressure-enthalpy-production-closure.md)
- [V1.0 Re=3900 文献与产品证据边界](docs/research/2026-09-02-v1-re3900-release-evidence-boundary.md)
- [Re=3900 冻结文献边界](docs/research/2026-08-31-v04-re3900-thin-domain-literature-boundary.md)
- [文档入口](docs/index.md)

## 许可证

HUNDUN-FLOW 采用 Apache License 2.0，见 [LICENSE](LICENSE)。第三方组件及其许可证见 [THIRD_PARTY.md](THIRD_PARTY.md)。

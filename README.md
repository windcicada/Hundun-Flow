# HUNDUN-FLOW 1.0.0

HUNDUN-FLOW 是一个使用 C++17 和 MPI 编写的低马赫数有限体积求解器，面向 Cartesian 网格上的单相流、LES 和静止浸没边界（IBM）计算。它适合需要明确数值门限、可恢复计算和运行证据的 CPU 集群任务。

公开产品版本是 `1.0.0`。当前实现仍保留 `versions/v0.4` 目录和 `hundun::v04` 命名空间，便于已有算例和 API 继续使用；冻结的 `versions/v0.3` 源码线也可以独立构建。

## 为什么开发 HUNDUN-FLOW

这个项目有三个直接目的：

- 独立实现一套可检查的低马赫数压力、焓和密度耦合流程；
- 在 MPI 并行计算中保留终端残差、通量、Restart 和候选身份，失败时回滚，不把异常状态写成已接受结果；
- 在相同网格、物性和边界输入下对比 COAST 等既有程序的计算成本，同时保留温度相关热力学和输运模型，不用常 `cp`、常分子黏度或常导热系数换取速度。

## 已实现功能

- tensor-stretched Cartesian 网格、MPI Cartesian 分解，以及 COAST runtime axes 格式导入；
- 单相低马赫理想气体，支持温度和组分相关的 `cp`、密度、分子黏度与导热系数；
- 可选择 `PISO` 或 `SIMPLE` 耦合路径，两次 pressure corrector 和 same-target continuity--energy nonlinear refinement；
- BDF 时间推进、固定或自适应时间步、retry/rollback 和 provisional/committed CFL 检查；
- 入口、出口、周期、对称等外边界；
- 静止封闭 STL IBM，支持 `strict_quadratic` 和可审计的 `adaptive_order` 重构；
- WALE、Vreman 和 wall-function 路径；
- common-face owner AFC、最终面质量通量，以及 EOS、continuity、energy、closed-mass、gauge 五项终端检查；
- exact-history MPI Restart、Visit 输出；

[当前版本能力](docs/releases/current-capabilities.md)  [已接受能力](docs/verification/accepted-capabilities.md)。

## 计划功能

- 反应流、燃烧化学、喷雾、颗粒、多相和辐射；
- 移动或变形 IBM、相交表面、非流形和未封闭几何；
- AMR、嵌套网格、非结构网格和 GPU 后端；
- 可压缩激波、声学和高马赫数求解；
- 更广泛的网格/时间步独立性研究、长时间统计和复杂工程几何验证。

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
## 运行

先检查算例和 MPI 分解，再启动计算：

```sh
mpirun -np 4 build/release/versions/v0.4/hundun validate case --dry-plan
mpirun -np 4 build/release/versions/v0.4/hundun run case \
  --output run --steps 10 --output-interval 0 --restart-interval 10
```

`--output-interval 0` 关闭 Visit、screen 和 monitor 输出；Evidence 仍会写入运行目录。正式计算前应检查边界方向、IBM donor、CFL、正性、质量与能量残差、limiter 活性和统计窗口。

## 已完成测试

源码测试覆盖单元、数值、MPI、Restart、I/O、回滚和命令行路径。与 1.0.0 SIMPLE 优化直接相关的检查包括：

- SIMPLE diagonal Schur 选择策略的 RED→GREEN 测试；
- 1、2、4 ranks 的 SIMPLE+IBM 非零 refinement 测试；
- pressure--energy Schur、globalization、Krylov 和 pressure-energy retry 测试；
- 128 ranks 的严格算例静态检查、dry-plan 和一次正式性能运行；
- 每步无 retry，所有线性求解完成，EOS、continuity、energy、closed-mass、gauge 门限通过。

严格性能算例使用以下共同输入：

- `D=0.02 m`，`Uc=2.89668 m/s`，`Re=3900`；
- 计算域 `[-0.10,0.30] × [-0.10,0.10] × [0,0.06] m`，即 `20D × 10D × 3D`；
- `456 × 256 × 104 = 12,140,544` cells，HUNDUN 与 COAST 使用同一组坐标；
- 128 MPI ranks，每 rank 一个物理核心；
- COAST 原生空气模型：NASA7 `cp(T)`、Yoon--Thodos/Wilke `mu(T,Y)`、`k=cp*mu/0.70`；
- 每条路径只运行一次，其中 1 步启动，随后 5 步计时。

| 路径 | 五个热步中位数 | 结论 |
| --- | ---: | --- |
| HUNDUN PISO，同提交基线 | 27.624453 s | 基线 |
| HUNDUN SIMPLE，同提交基线 | 26.437491 s | 比 PISO 低 4.30% |
| HUNDUN SIMPLE，1.0.0 优化路径 | 17.540536 s | 比普通可压缩 COAST 低 17.68% |
| 普通可压缩 COAST | 21.306931 s | 对照 |

优化集中在 SIMPLE+IBM 的 pressure-energy refinement。旧路径在首次 C2 后切回空间 `E_p/E_h` Schur，而预条件器仍是连续性压力 MG；1.0.0 让后续 refinement 继续使用 double-diagonal Schur，省去了相应的空间 stencil 和焓 halo。Stage 50 中位数从 21.319866 s 降到 13.268423 s。refinement 次数和线性迭代没有减少，因此这项结果不能归因于放宽收敛条件或减少物理计算。

本次 HUNDUN 最大 continuity 和 energy 残差分别为 `7.451850e-7` 和 `9.369710e-7`，均低于 `1e-6`。性能证据对应优化提交 `fca607cf483e39996461619bb83295ee8cb05c98`，Evidence SHA-256 为 `ce9616316685ad3fc5f62b98be1d1d474677f0f8cea538a780007d24a607b274`。

## 文档

- [快速开始](docs/user-guide/quick-start.md)
- [配置说明](docs/user-guide/configuration.md)
- [控制方程](docs/numerics/governing-equations.md)
- [离散方法](docs/numerics/discretization.md)
- [Restart](docs/user-guide/restart.md)
- [诊断输出](docs/user-guide/diagnostics.md)
- [文档入口](docs/index.md)

## 许可证

HUNDUN-FLOW 采用 Apache License 2.0，见 [LICENSE](LICENSE)。第三方组件及许可证见 [THIRD_PARTY.md](THIRD_PARTY.md)。

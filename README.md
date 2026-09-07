# HUNDUN-FLOW

HUNDUN-FLOW 是面向流动与湍流燃烧研究的 C++17/MPI 数值模拟程序。开发目的，是把方程、离散、物性、并行通信和验证依据放在可检查的代码中，让研究者能够修改方法，并检验修改是否保持守恒、收敛与重启一致性。

当前产品版本为 **1.0.0**，已接入的生产路径是单相低马赫数流动。版本号不代表 Re3900 长期统计和实验对比已经完成。燃烧、TPDF/TCR 和喷雾仍是后续研发方向，不应当作当前可执行功能。

## 当前能力与限制

| 模块 | 已接入内容 | 仍需注意 |
| --- | --- | --- |
| 网格与 IBM | 均匀/张量拉伸笛卡尔网格、COAST 轴坐标文件、静止封闭 STL、局部二次/自适应降阶重构 | 不提供 AMR、移动 IBM 或任意曲面的统一二阶保证 |
| 流动与物性 | SIMPLE/PISO、压力—焓耦合、理想气体 EOS、温度相关热力学与输运物性 | 不面向激波、高马赫数或声学问题 |
| 时间推进 | 固定步长、变步长 BDF2、失败回退和必要的 BE 恢复 | 自动非对流时间尺度尚未完整接入 |
| 标量与 LES | 被动标量、参与 EOS 的组分输运、配对质量 remap、组成迭代、WALE/Vreman | 组分输运不等于已接入有限速率化学；曲面标量扩散精度仍有限制 |
| 并行与存储 | MPI 分区、持久 halo/donor 通信、共享归约与预分配工作区 | 完整进程峰值硬预算仍待完成 |
| 重启与输出 | 方法历史签名、精确续算、显式方法恢复、Visit、Evidence 和分模块观测 | 方法恢复不继承旧统计样本；观测溢出不能视为完整归因 |

## 源码组织

仓库只保留当前实现。下列 `versions/v0.4`、`v04_*` 和 `hundun::v04` 是当前源码、接口及证据格式沿用的兼容名称，**不是另一个旧版程序**。不为整理目录而改名，以免破坏构建、公共接口和证据追溯。退休实现可从 Git 历史取得。

| 位置 | 职责 |
| --- | --- |
| `versions/v0.4/include/hundun/` | 当前公共接口、字段视图、状态与计划 |
| `versions/v0.4/src/app_*`、`core_*` | 输入编译、应用驱动、资源分配、时间步状态的一致提交与回退 |
| `versions/v0.4/src/solver_*` | 有限体积算子、压力—能量耦合、标量、Krylov/MG |
| `versions/v0.4/src/mesh_*`、`bc_*`、`physics_*` | 网格、IBM、边界、热力学、输运与 LES |
| `versions/v0.4/src/parallel_*`、`io_*` | MPI 通信、Restart、Visit 与诊断 |
| `versions/v0.4/tests/` | 单元、MPI、制造解、故障注入与 CLI 回归 |
| `tools/` | 专用圆柱 runner、观测、证据校验与后处理 |
| `docs/` | 使用指南、数值说明及按提交记录的验证证据 |

## 构建与运行

需要 CMake 3.21、支持 C++17 的编译器、MPI 3 的 C 接口及 POSIX 线程库。源码随附 yyjson，默认构建不下载依赖。

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DHUNDUN_BUILD_TESTS=OFF
cmake --build build/release -j 2 --target hundun
build/release/versions/v0.4/hundun --version
mpirun -n 1 build/release/versions/v0.4/hundun validate examples/minimal --dry-plan
mpirun -n 1 build/release/versions/v0.4/hundun run examples/minimal \
  --output run-minimal --steps 10 --output-interval 10 --restart-interval 10 \
  --initial-state 101325,300,0.1,0,0
```

最小算例只验证输入和运行路径，不是圆柱算例或科学精度证明。输出目录应独立。需要回归时使用 `-DHUNDUN_BUILD_TESTS=ON`；运行测试前先确认没有占满机器的长测，不与它并发争抢资源。

## 已完成验证与结论

2026-09-07 的[模块验收报告](docs/verification/2026-09-07-exclusive-module-acceptance.md)记录了多标量容量、IBM 公共接口、初始化、方法恢复、统计 epoch、观测完整性和分配失败路径的检查。相关 1/2/4-rank MPI 与 ASan/UBSan 回归已完成，最终针对性验收为 27/27；这不是仓库全部测试或长期物理验收。

同一 Re3900 窗口、128 ranks、各一轮的[局部性能实验](docs/verification/experiments/2026-09-07-sparse-ibm-candidate.md)未显示总耗时收益，已撤回实验代码。当前没有据此宣称快于 COAST。

当前冻结圆柱长测采用 D=0.02 m、Uc=2.89668 m/s、Re=3900，计算域为 20D×10D×(π/2)D，网格 456×256×52，固定 dt=1.3808912271980336×10⁻⁵ s。变物性、守恒门槛及 checkpoint 合同保持不变。长测仍在进行，不能提前宣称统计收敛或实验吻合；它对应验收报告中的冻结程序，不会随 main 更新而自动替换。

后续工作包括曲面标量精度、完整内存预算、非对流时间尺度、长日志分段观测，以及在测量支持下逐项优化。燃烧与喷雾应在独立接口和回归闭合后再接入生产路径。

## 文档与许可证

- [快速开始](docs/user-guide/quick-start.md) · [命令行](docs/api/cli.md) · [输入配置](docs/api/configuration-schema.md)
- [重启与方法恢复](docs/user-guide/restart.md) · [当前能力](docs/releases/current-capabilities.md) · [文档入口](docs/index.md)

采用 Apache License 2.0，见 [LICENSE](LICENSE)。第三方来源见 [THIRD_PARTY.md](THIRD_PARTY.md)，贡献须遵守 [DCO](DCO.md)。

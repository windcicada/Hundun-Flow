# 当前模块结构

公共头文件位于 `versions/v0.4/include/hundun`，实现位于 `versions/v0.4/src`。不再保留根目录的重复实现。兼容命名空间为 `hundun::v04`。

| 前缀 / 位置 | 职责 |
| --- | --- |
| `app_*` | CLI、输入模型编译、应用运行和构建身份 |
| `core_*` | 字段目录、arena、资源计划、ProductDriver、提交与回退 |
| `mesh_*` | 笛卡尔网格、MPI 分解、STL、IBM donor 与重构 |
| `bc_*` | 外边界、热力学边界、时间控制 |
| `physics_*` | 热力学、输运、派生物性、LES 与贡献接口 |
| `solver_*` | 离散算子、压力—焓、动量、标量、Krylov/MG |
| `parallel_*` | CPU 资源、普通/prepared halo、IBM donor 通信 |
| `io_*` | 输出计划、Visit、Restart、监控和 Evidence |
| `tools/v04_thin_domain_runner.cpp` | 圆柱试验的运行、统计、性能记录 |
| `versions/v0.4/tests` | 独立 test core 和回归夹具 |

接口、计划、借用视图和容量共同构成调用合同。物性、通信和历史的权威来源不能由另一个模块静默覆盖。具体执行流见[数据流](data-flow.md)。

# 模块结构

公共头文件位于 `include/hundun`，实现位于 `src`。文件名前缀说明所属模块：

| 前缀 | 职责 |
| --- | --- |
| `cfg_` | JSON 配置、规范化和 schema |
| `rt_` | MPI、字段、Halo、分区、Restart 和 VTK I/O |
| `mesh_` | 结构化拓扑与几何 |
| `bc_` | 外边界条件 |
| `lin_` | 线性代数和迭代求解 |
| `fvm_` | 有限体积离散与 Poisson 算子 |
| `ib_` | 浸入表面、Ghost-Cell 重构和壁面力 |
| `flow_` | 时间推进、压力速度耦合和流动状态 |
| `diag_` | 结构化诊断 |
| `exec_` | CPU 执行与存储契约 |
| `sdk_` | 公共插件接口 |

`app_` 文件是可执行程序的私有编排层，不属于公共 C++ API。

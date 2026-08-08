# HUNDUN-FLOW 源码目录与命名规范

本文规定 HUNDUN-FLOW 产品源码的目录、文件、C++ 标识符和 CMake 目标命名。新增文件和重命名应遵守本文；已有公共接口若与本文不一致，以兼容性为先，不得仅为统一外观而改名。

## 1. 基本原则

命名要先把职责说清楚，再考虑缩短。文件名应直接表明所属领域和主要内容，不能靠目录层级补足含义。

产品源码采用扁平目录：

```text
include/
└── hundun/        公共头文件
src/               实现文件和私有头文件
```

`include/hundun/` 和 `src/` 下不得再按模块建立子目录。扁平结构依靠受控前缀区分领域，避免同名文件和含义不清的缩写。

## 2. 文件名前缀

产品文件必须使用下表中的前缀。

| 前缀 | 领域 | 示例 |
| --- | --- | --- |
| `app_` | 命令行程序、驱动和装配 | `app_main.cpp` |
| `cfg_` | 配置、解析和已解析配置 | `cfg_case_config.hpp` |
| `rt_` | 运行时、MPI、字段存储和重启基础设施 | `rt_mpi_context.hpp` |
| `exec_` | 执行空间和执行资源 | `exec_execution.hpp` |
| `mesh_` | 网格拓扑与几何 | `mesh_geometry.hpp` |
| `bc_` | 边界条件 | `bc_basic_boundary.hpp` |
| `ib_` | 浸入边界几何、重构和表面量 | `ib_surface.hpp` |
| `lin_` | 线性代数、求解器和预条件 | `lin_bicgstab.hpp` |
| `fvm_` | 有限体积离散和离散算子 | `fvm_cell_centered.hpp` |
| `flow_` | 流动状态、时间推进和压力速度耦合 | `flow_state.hpp` |
| `diag_` | 诊断、输出和性能记录 | `diag_session.hpp` |
| `sdk_` | 稳定的扩展接口与加载器 | `sdk_plugin_api.h` |

前缀只写一次。领域词已经在前缀中出现，主体名称就不再重复。例如写 `flow_state.hpp`，不写 `flow_flow_state.hpp`；写 `mesh_geometry.hpp`，不写 `mesh_mesh_geometry.hpp`。

文件名采用小写蛇形命名。可以使用行业通用缩写，不要自造缩写。名称只描述长期职责，不写临时序号、开发批次或完成状态。

以下前缀等到产品开始提供对应能力时再启用：

- `les_`：大涡模拟模型；
- `chem_`：化学动力学与物种反应；
- `comb_`：燃烧模型；
- `spray_`：喷雾模型；
- `part_`：离散颗粒模型。

预留前缀不能用来放占位文件。代码尚未形成独立领域时，也不要提前拆出新前缀。

## 3. 公共头文件

公共头文件只放在 `include/hundun/`，扩展名使用 `.hpp`；需要 C ABI 时可以使用 `.h`。包含路径统一写成：

```cpp
#include "hundun/flow_state.hpp"
#include "hundun/sdk_plugin_api.h"
```

公共头文件应满足以下要求：

- 单独包含时语义完整；
- 只暴露调用方需要的类型、函数和常量；
- 不包含 `src/` 中的头文件；
- 不暴露实现布局、故障注入接口或内部辅助类型；
- 不因实现重排而改变公共类型名、函数签名、命名空间或数据语义。

公共文件名应对应文件中的主要概念。一份头文件可以带上紧密相关的辅助类型，但不能成为跨领域的汇总入口。

## 4. 实现文件和私有头文件

实现文件放在 `src/`，扩展名使用 `.cpp`。通常，一个公共头文件对应一个同名实现文件：

```text
include/hundun/flow_state.hpp
src/flow_state.cpp
```

仅供实现使用的头文件也放在 `src/`，文件名以 `_detail.hpp` 结尾：

```text
src/flow_checkpoint_v2_detail.hpp
src/ib_surface_bvh_detail.hpp
```

`_detail.hpp` 不构成公共接口，不得被公共头文件包含。多个实现文件共享私有逻辑时，优先使用职责明确的私有头文件；只在单个实现文件中使用的辅助函数放入匿名命名空间。

`.cpp` 文件的包含顺序为：

1. 与本实现对应的公共头文件；
2. 本领域私有头文件；
3. 其他 HUNDUN-FLOW 公共头文件；
4. 第三方头文件；
5. C++ 标准库头文件。

各组之间留一个空行。不要依赖其他头文件的间接包含。

### 产品代码与测试支撑的边界

`tests/support/` 只放测试使用的适配器、故障注入、状态观察和测试夹具。这些文件不属于产品，也不会进入 tests-off 构建或产品仓库。

依赖方向只能从测试指向产品：

```text
tests/ -> include/hundun/
tests/ -> src/
src/   -X-> tests/
```

产品 `.cpp` 和 `src/*_detail.hpp` 不得包含、引用或通过宏拼出 `tests/` 下的路径。产品实现若需要与测试共享内部声明，应把产品侧的原始接口或中性桥接声明放在 `src/*_detail.hpp`，签名只使用产品类型、标准库类型或基础标量；`tests/support/` 再负责把这些接口转换成测试专用类型。测试专用枚举、typed snapshot、mutation 控制器和断言辅助类型不得进入产品头文件。

测试构建可以通过专用编译宏启用低成本观察或故障注入。宏关闭时，默认产品路径不得依赖测试状态，MPI collective 顺序、回滚、数值算法和公开接口也不能因此变化。CMake 中产品目标只允许公开使用 `include/`、私有使用 `src/`；仓库根目录和 `tests/` 不能成为产品目标的包含目录。

## 5. C++ 标识符

命名空间保留完整领域名称，例如：

```cpp
namespace hundun::finite_volume { }
namespace hundun::immersed { }
namespace hundun::diagnostics { }
```

文件名前缀是目录组织规则，不用于缩短命名空间。不得把 `finite_volume` 改成 `fvm`，也不得把 `immersed` 改成 `ib`。

C++ 标识符采用以下形式：

| 对象 | 形式 | 示例 |
| --- | --- | --- |
| 类型、类、枚举 | `PascalCase` | `FlowState` |
| 函数、方法 | `snake_case` | `advance_one_step` |
| 局部变量、参数 | `snake_case` | `cell_count` |
| 私有数据成员 | `snake_case_` | `local_size_` |
| 编译期常量 | `kPascalCase` | `kMaximumReach` |
| 枚举值 | `snake_case` | `collective_failure` |
| 宏 | `HUNDUN_UPPER_SNAKE_CASE` | `HUNDUN_ENABLE_ASAN` |

布尔名称应能直接读出真假含义，优先使用 `is_`、`has_`、`can_` 或明确的状态词。单位写进名称或类型，例如 `time_seconds`、`length_metres`；同一接口中不要依赖注释区分单位。

函数名通常由动作和对象组成。查询函数可以使用名词或 `find_`、`lookup_`；会改变状态的函数要写明动作。不要使用 `handle`、`process`、`do_work` 这类看不出结果的名称。

## 6. 接口与语义稳定性

重命名不能顺带改变以下内容：

- 公共类型和函数的含义；
- 配置键及其默认值；
- 重启数据的字段、版本和读取规则；
- 诊断字段、错误码、单位和符号约定；
- MPI collective 的参与顺序；
- 数值算子的权威数据来源。

确需改变公共名称或语义时，要先给出兼容策略，并把名称变化与行为变化分开处理。整理文件名时不要修改算法、阈值或默认运行参数。

## 7. CMake 目标

CMake 目标使用完整、稳定的领域名称，以 `hundun_` 开头：

```cmake
hundun_runtime
hundun_mesh
hundun_boundary
hundun_immersed
hundun_linear
hundun_fvm
hundun_flow
hundun_diagnostics
hundun_sdk
```

目标名面向构建依赖，文件名前缀面向源码导航，两者不要求逐字相同。可执行程序使用 `hundun`。

根目录 `CMakeLists.txt` 只负责项目声明、全局选项、依赖发现和子目录入口。产品目标及其源码清单放在 `src/CMakeLists.txt`。公共包含目录必须指向 `include/`；`src/` 只能作为目标的私有包含目录。

新增源码必须显式登记，不用递归通配自动收集产品文件。文件增删后，构建清单应能直接显示变化。

## 8. 新增文件的判断顺序

新增文件前按以下顺序判断：

1. 确认职责属于哪个已登记领域；
2. 选择对应前缀，并去掉主体中的重复领域词；
3. 判断调用方是否需要该声明，决定放入公共目录还是作为 `_detail.hpp`；
4. 检查名称是否表达物理量、算法对象或基础设施职责；
5. 将实现显式加入对应的 CMake 目标；
6. 检查公共接口、配置和持久化数据是否保持兼容。

若一个文件同时承担多个领域职责，应先拆清职责，再决定命名；不要用多个前缀拼接来回避设计问题。

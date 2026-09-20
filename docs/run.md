# 安装、运行与续算

HUNDUN-FLOW 使用 C++17、MPI 和 FP64。当前生产入口采用 CN/BE 外迭代压力耦合，
Vreman LES、有限速率化学、ESF/TPDF 和拉格朗日喷雾按算例显式配置。
动态混合提供 `cdphyso_dynamic_v1` 与 `dyn711_v1`，模型身份随 Restart 保存。

## 源码构建

准备 CMake 3.21+、Ninja、Clang 15、lld、libstdc++ ABI1、MPI 3+ 和
固定 Cantera 3.2 SDK。`VERSION` 提供产品版本，`versions/v0.4` 为实现目录。

```sh
HUNDUN_CANTERA_ROOT=/path/to/cantera cmake --preset release
cmake --build --preset release -j 2
export PATH="$PWD/b3/versions/v0.4:$PATH"
hundun --version
```

## 新算

```sh
hundun init-case --output case
mpiexec -n 4 hundun check case --dry-plan
mpiexec -n 4 hundun run case --output run --steps 10 \
  --output-interval 0 --restart-interval 10 --diagnostics-interval 1
```

`hundun check` 显示时间格式、压力耦合、化学后端、模型及资源计划。
局部 CFL 默认目标为 0.30，保持区间为 0.25–0.35；算例中的显式时间控制优先。

[气相示例](../examples/g/README.md) 使用 JL4、四场 ESF 和动态 TCR；
[气液示例](../examples/cf/README.md) 使用煤油四步反应、两场 ESF、
THICK_EX、SGS 破碎及共同气液交换。

## Restart

```sh
mpiexec -n 2 hundun run case --restart run/Restart --output next \
  --steps 10 --restart-interval 10 --diagnostics-interval 1
```

保持同一算例与物理资产。Restart 保存当前场、历史场、面通量、随机场、
模型时钟及颗粒状态。各分区写出并复读校验，随后发布完整清单和当前代指针。
进程数调整由原生场重分区入口处理。

一次性外部迁移通过格式名选择：

```sh
hundun import transfer --format pdf-transfer-v2 --case case --output seed \
  --model-history initialize
```

迁移报告记录场映射和模型历史初始化。`--model-history initialize` 从指定
流场建立本版本历史；原生 Restart 用于连续保留已接受的模型统计窗口。

## 输出与预算

`monitor.jsonl` 记录步长、方程迭代和模块时间；`diagnostics.jsonl` 记录
质量、能量和组分预算；`evidence.jsonl` 绑定程序、输入及接受状态身份。
预算使用物理随机场均值与 field0 耦合密度的明确分工。
原始预算缺陷与独立方程残差分别保存，支持查看有限迭代误差与离散账本差额。

`--output-interval`、`--restart-interval`、`--diagnostics-interval` 分别配置
云图、检查点和诊断频率。云图提供速度、压力、温度和 SGS 派生量。

# 配置 schema

程序按根对象中的 `schema_version` 选择格式。所有版本都拒绝未知 key 和重复成员。

## schema 1

用于保守被动标量路径。根对象包含 `case`、`resources`、`mesh`、`time`、`transport`、`initial_condition`、`restart` 和 `output`。该格式保留用于兼容已有输入；新流动计算应使用 schema 2 或 3。

## schema 2

用于低马赫数流动。主要对象如下：

| 对象 | 关键字段 |
| --- | --- |
| `case` | `name` |
| `simulation` | `type=variable_density_flow`、`density_model` |
| `resources` | 可选 `expected_ranks`、`process_grid` |
| `mesh` | `cells`、`origin_m`、`length_m`、`mapping`、可选 `warp_amplitude` |
| `time` | `mode`、`steps`、步长范围、CFL/扩散数目标和重试参数 |
| `physics` | 参考密度、动力黏度及所选密度模型需要的热力学量 |
| `scalars` | 标量名和 `diffusivity_m2_per_s` |
| `boundaries` | 六个 patch，每个 patch 恰好一条记录 |
| `restart` | 读取、写出目录和间隔 |
| `diagnostics` | 目录、间隔和网格开关 |
| `performance` | 是否启用及采样设置 |

`density_model` 可取 `constant`、`material`、`ideal_gas`。`mesh.mapping` 可取 `uniform_box` 或 `analytic_warped_box`。`time.mode` 可取 `fixed` 或 `adaptive`。

## schema 3（0.1.0 仅解析与校验）

在 schema 2 的共同流动字段上增加：

- `immersed_boundary`：`model` 为 `none` 或 `local_flow_pattern_ghost_cell`；启用后必须提供 STL `geometry` 和静止 `wall`；
- `les`：格式接受 `none` 或 `wale`，但当前公开验收能力仅覆盖 `none`。

当前静止壁面要求：

```json
{
  "velocity_m_per_s": [0.0, 0.0, 0.0],
  "enthalpy": "zero_normal_diffusive_flux",
  "scalars": "zero_normal_diffusive_flux"
}
```

`geometry.format` 必须是 `stl`，`length_scale_to_m` 必须为有限正数，`fluid_side` 为 `inside` 或 `outside`。

0.1.0 的配置库可以读取、规范化和广播 schema 3，但 `hundun` 命令行 driver 尚未装配该运行路径。使用 `--validate` 或 `--print-resolved` 成功，只证明输入符合 schema，不表示 IBM 算例能够由当前可执行程序推进或 Restart。

最可靠的字段检查方式是运行 `hundun case.json --print-resolved`，并把输出与输入一同存档。

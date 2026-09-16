# 液相物性与绝对焓

液体资产采用 `HUNDUN_LIQUID_ASSET_V1` 和 SI 单位，绑定气相机理、
组分顺序、蒸气组分及焓参考。资产文件的完整字节内容生成材料身份，
该身份随颗粒进入运行状态与 Restart。

密度、比热、潜热、表面张力和黏度支持 `constant` 与 `cubic`。
`cubic` 采用 `theta = T - T_ref`，四个系数依次对应零至三次项。
煤油相关系提供以下固定系数选项：

| 资产字段 | 关系标识 | 输出单位 |
|---|---|---|
| `density` | `kerosene_density_v1` | kg/m³ |
| `cp` | `kerosene_cp_v1` | J/(kg·K) |
| `latent` | `kerosene_latent_v1` | J/kg |
| `saturation` | `kerosene_v1` | Pa |

对应的物性行写为：

```text
density kerosene_density_v1 298.15 0 0 0 0
cp kerosene_cp_v1 298.15 0 0 0 0
latent kerosene_latent_v1 298.15 0 0 0 0
```

`surface_tension` 和 `viscosity` 行随后给出所选关系与系数，
饱和蒸气压行采用 `saturation kerosene_v1`。固定关系的四个用户
系数保持为零；各标识对应表中指定的物理量。

煤油比热、潜热和饱和蒸气压采用亚临界温度区间
`43 K < T < 684.26 K`。算例在该区间内设置材料温度范围；材料
范围与气相机理的适用范围共同约束实际查询。

液相绝对焓由资产中的 `liquid_reference T_ref h_ref` 定义：

\[
h_l(T)=h_{ref}+\int_{T_{ref}}^T c_{p,l}(\theta)\,d\theta.
\]

煤油关系使用解析积分，并保留相邻浮点温度之间的小焓增量。
`h_ref` 按算例的气相形成焓与相变参考确定；颗粒库存、注入库存及
共同交换使用同一液相焓。通用气膜审核同时核对 `h_v = h_l + L`
及 `c_{p,v} = c_{p,l} + dL/dT`。煤油潜热导数使用其真实温度关系。

饱和蒸气压资产返回分段关系值；气膜模型按当地压力执行相平衡
处理。THICK_EX 的气膜查询与冻结物性核已有逐项参考证据，
颗粒推进、共同交换和 624CF 实场接线进度见 [开发记录](rc.md)。

原生喷雾配置通过 `spray.evaporation` 选择蒸发模型：

```json
"evaporation": "thick_exchange"
```

`thick_exchange` 使用煤油密度、比热、潜热和饱和蒸气压关系，
气膜分子黏度来自流动求解器的输运物性计划。`abramzon_sirignano`
沿用已有输入默认值和 Restart 身份。

THICK_EX 的区间系数取子步起始状态，质量与温度遵循原有 BE
蒸发关系。速度采用现有 Schiller–Naumann 阻力的隐式更新，轨迹
采用端点平均速度。每个端点按新的液体密度由质量重算直径。
自适应误差控制使用质量、速度和温度；气液交换由相同端点的
质量、动量、绝对焓和动能库存确定，各段同时接受独立守恒审核。
A–S 的独立热通量积分沿用其交换增量误差控制。

原生 BE/PISO 与 CN/BE＋`outer_corrected` 入口支持真实煤油
物性、原生四步反应、喷射和气液共同交换。化学与颗粒组分源保留
各自注册身份，在目标组分方程中共同消费；CN/BE 的保守总能量
方程同时接收气液焓与动能交换。候选收支和最终收支使用同一步
交换账本，颗粒历史与流场完成共同提交。

原生验证入口为 `v04_spray_thick_cli` 和 `v04_spray_thick_cn_cli`，
覆盖 1→2→4 进程续算及同一检查点的 1/2 进程库存比较。
煤油机理资产通过 `hundun-chemistry` 指定
`kerosene_4_step_v1/frozen_material`，Cantera 提供热物性与积分器。
CN/BE＋ESF 的双随机场喷雾入口由 `v04_spray_thick_esf_cli` 验证。
每个随机场先接收同一保守交换，再完成输运、IEM 与化学；field0
承载降噪密度，随机场均值承载物理组分和焓。总能量方程同时消费
颗粒与统计输运贡献，两个账本分别记录。验证覆盖逐组分相间
库存、元素收支、完整场及颗粒记录的跨进程恢复。
`v04_spray_thick_ibm_cli` 覆盖静态立方体 IBM、双随机场、四步反应
和 THICK_EX 的共同推进，核对撞壁反射、流体侧颗粒位置及
1→2→4 进程续算。组分库存直接从实际随机场以扩展精度累计，
field0 继续提供压力耦合密度；守恒审核沿用原有阈值。
`v04_spray_thick_sgs_cli` 在该组合中启用 Vreman，逐分区读取
涡黏度、SGS 动能及两种耗散率，核对流体响应、固体零值和
同一步的跨进程输出。验证记录同时包含相间能量及组分／元素收支。
624CF 几何映射及实场验证继续按开发记录推进。

COAST 喷雾检查点的 SI 解码入口为 `tools/v04_spray_inventory.py`。
它按原生索引表区分当前槽和退役槽，分别保留有效液滴与待加入
子滴，将微米直径转换为米，并输出带来源身份和原始模型历史的
`parcels.jsonl`。完整记录及源文件哈希进入 `inventory.json`，
后续按目标 IBM 网格、物性和模型定义构建 Hundun Restart。

确定性／SGS 破碎核通过 `evaluate_sgs_breakup` 提供区间末候选。
它接收单位质量耗散率、两段暴露历史、已接受的 Poisson 倍率和
本次事件的随机数，返回破碎速率、候选历史及互补体积的子滴
直径比。调用方负责事件身份、子滴共同提交和气液账本。20 点
子滴分布与完整 COAST 例程的参考入口为
`tools/v04_sgs_breakup_reference.py`；生产接线及耗散率定义见开发记录。

SGS 破碎的双子滴候选采用互补体积和相同母滴权重，复用公共
质量、动量、平移动能和绝对焓审核。子滴身份支持复算及迁移
一致性；小子滴优先的质量分配保留 CDF 端点的可表示库存。
90 组分裂和 13 组异常请求通过，表面能单列于破碎账目。
生产历史与事务接线的阶段记录见 `docs/rc.md` S18。

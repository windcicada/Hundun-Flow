# 启动与运行控制

在算例目录放置 `case.json`、物性资产与 `run.json`，执行
`mpirun -np 4 /path/to/hundun`。`hundun init-case --output c` 生成示例。

```json
{"mode":"new","steps":100,"output":"../out","monitor_interval":1,
 "output_interval":20,"restart_interval":50,"diagnostics_interval":10}
```

`mode` 为 `new` 或 `restart`。恢复时添加 `"restart":"../out/Restart"`。
`steps` 表示本次新增步数；`end_time` 表示第一个达到该物理时间的接受步，
两者同时存在时取先达到的条件。时间步继续遵循算例时间控制。
`initial_state` 为 `[p,T,Ux,Uy,Uz,q...]`，用于指定新算均匀场。
所有相对运行路径以 `run.json` 所在目录为基准，输出目录与输入目录采用同级布局。

显式启动保持 `hundun run CASE --output OUT --steps N`。
`hundun run CASE --config FILE --steps N --output OUT` 先读运行文件，
再按具名键应用命令行覆盖。根进程解析后广播同一组最终参数。
启动信息记录最终步数、周期、方法指纹、进程数与输入／输出／恢复路径。

- `--monitor-interval N`：屏幕摘要和 `monitor.jsonl`，默认每接受步。
- `--output-interval N`：三维场输出。
- `--restart-interval N`：原生检查点。
- `--diagnostics-interval N`：详细守恒账本。
- `--until T`：显式命令的物理终止时间，配合 `--steps` 使用。

各周期独立配置，零值表示关闭该类定期输出。接受步摘要含 dt、
CFL 数值及定义、方程残差、迭代次数和最大进程推进壁钟时间。
`cn_phases.seconds` 保存 CN/BE 的准备、标量、物性、动量组装／求解、
压力组装／预处理／求解、校正及审核耗时。各阶段累计全部候选与重试，
包含阶段内部通信；每阶段取最大进程值，其所属进程可以随阶段变化。
整步时间使用 `seconds`，BE/PISO 的 CN 阶段项为零。

在输出目录创建 `stop` 文件，本步接受后保存 Restart 并停止；
创建 `output` 文件，本步接受后额外保存场数据。请求由根进程接收，
全体进程执行，写入成功后登记 `control.jsonl` 回执并消费请求。
写入期间采用 `.pending` 标记保存请求身份，新到请求留给下次处理。
`hundun status OUT` 查看 `status.json` 的阶段、接受步和更新时间。

三维场默认采用 legacy 二进制 `.vtk`，在 `Visit/solution.visit` 中按真实
物理时间组织。向量 `Velocity` 与标量 `Temperature`、`Pressure`、
`PressureGauge`、`Density`、`Enthalpy` 使用 SI 单位。Pressure 为
参考压力加 pi，PressureGauge 为 pi。双状态 ESF 的 Density 对应
压力耦合 field0 密度；Temperature 与 Enthalpy 对应物理均值热状态。
原有变量目录继续通过 `--visit-format xml` 或运行文件
`"visit_format":"xml"` 选择。

采样点为实际单元中心，正方向分区交界增加相邻中心点，并以
`vtkGhostType` 标记重复点。格子内部、交界面、棱和角的值由输出专用
交换取得。时间索引在各分区数据完成后原子发布；每个输出目录保存
固定进程数的序列，切换进程数时选用新的输出目录。

`hundun restart-info OUT/Restart` 报告版本、步号、时间、dt、源进程数、
网格、方法／模型历史身份和文件校验结果。默认逐文件流式核对校验和，
工作缓冲区为 64 KiB；`--metadata-only` 读取并核对清单。
具体算例的模型兼容性由原生恢复入口结合 case.json 检查。

气相时间配置支持 `time.convective_cfl_definition`：`outgoing_sum` 为
单元向外面质量通量之和，`directional_max` 为六个面质量通量绝对值
的最大值；两者均乘以 `dt/(rho*V)`。原版气相对齐输入显式选用
`directional_max`，既有输入延续 `outgoing_sum`。目标及浮动范围由
`convective_cfl`、`convective_cfl_margin` 控制，推荐 0.30 和 0.05。

`hundun check`、接受步摘要和运行证据报告实际定义。监看同时记录
向外通量和、绝对通量和及方向最大值；准入和自适应 dt 使用所选指标。
固定 dt 超过所选上限时返回定位信息；自适应尝试超过上限时沿统一
回退／重试流程处理。Restart 的方法历史签名携带 CFL 定义，原生
连续恢复沿用来源定义。

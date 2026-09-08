# Stage 4–6 原生耦合接线与圆柱续算验收

状态：P1–P8 原生接线、最终回归和干净程序验收通过，主线已合并。128 核短段对照通过，新程序圆柱长算已启动；启动前缀验收见下文。
本记录接续[源码归集记录](2026-09-08-integrated-stage46-re3900.md)，按用户
“完成耦合接线后才算合并完成，再续算圆柱”的新指示执行。

供体为 `codex/stage4-6-coast-replacement@f3ad77b0b0aa2105ff946a5a4680400ee41f1177`；
实际主线起点为 `de0e77ae715f43ba0e2cc0718f82029c80cb8c52`。
独立工作树为 `/home/wyf/code_dev/.worktrees/hundun-flow-integrated-stage4-6-re3900`。
其他开发工作树及供体未修改；本任务没有推送远程仓库。

## 原生消费合同

| 包 | 已接线的程序行为 | 验证入口 |
| --- | --- | --- |
| P1 | 解析/真实 Cantera 气相查询与推进；finite-rate/PaSR/ESF；完整 Ns 与原生 N−1、PH 状态、机理身份和化学源转换 | 1/2/4-rank 产品推进、真实 CLI、后端/组分/元素与失败检查 |
| P2 | 持续随机场、原生面通量/共同 Gamma 扩散/梯度、Philox/IEM、两个连续化学半区间、均值重整；物理边界与静态 IBM | 空间扰动、非零 WALE Dt、物理/IBM 边界和跨分区重启 |
| P3 | 实际 eta/R、shadow/experimental 反馈、显式折点证据消费、分支/折点/精确整数历史 | 原生 IEM 方差、六类错误证据撤回、V4/V5 持久状态 |
| P4–P6 | 实际液体资产与流体侧公共采样/沉积；A–S 自适应推进、几何/出口/壁面事件、TAB 及父删子建 | 两套合成资产、原生推进/反弹/TAB、库存与热/动量/动能账、热阶段分配检查 |
| P7 | 原生 owner/MPI 路由、跨单元迁移、注射器/粒子/TAB/ordinal/ID/RNG 时钟与恢复 | 1/2/4 ranks、全局 ID/容量/修订拒绝、真实 1↔4-rank 恢复 |
| P8 | 源先进入每个随机场与源后空间输运；Sm、组分、动量、总焓进入原生预测器、C1/C2 和最终压力/能量审计；所有参与者一次提交/回退 | 气液质量/总能量、当前步源只加入一次、PH/第二化学半区间/折点失败、实际 CLI 新算与续算 |

详细根因、RED/GREEN、增量提交和原始失败见[接线过程台账](../handoff/2026-09-08-product-coupling-progress.md)。

## 准入范围

ESF 与喷雾使用 backward Euler；静止物理壁面、出口、周期及静态 IBM 已接通。
移动壁面、AMR、稠密喷雾、液膜、多组分液滴等未在本供体合同内新增。
`validated` TCR 保持科学证据门；显式折点入口消费外部不可变证据，CLI 不自动识别折点。
当前资产/反应组合证明程序合同与数值收支，不替代真实燃料或湍流统计的独立科学验收。
内存检查覆盖原生 arena 和新增有界 owner/通信存储；不把它表述为全部后端/MPI 内部 RSS 硬上限。

## 最终验证与运行切换

最终 Clang [205/205](data/2026-09-09-native-stage46/clang-205.log)、GCC11/Cantera
[67/67](data/2026-09-09-native-stage46/gcc-cantera-67.log) 全部通过；两轮完整 JUnit
和构建目标/测试清单一并归档。IBM 增量首轮漏建后的补跑原始日志也保留。

干净源码提交 `83db41df2ec9c351c71e3aebe812509ddfd7262b`，tree
`5a1af1fca46e903402393bea63effce8ad086c3c`。Clang 和 GCC 构建清单的
core/target clean 均为 true；独立重算具名前缀 head/tree、核心内容和应用入口哈希一致。
GCC/Cantera 干净重建后，实际 IBM+ESF+喷雾 CLI 的 1→4 ranks 续算再次通过。
[冻结验收](data/2026-09-09-native-stage46/runtime/LOCAL_ACCEPTED.json)记录两个 ABI 的
程序哈希与 Cantera 环境；圆柱使用冻结的 Clang runner，SHA-256 为
`16dd8e972b1ac4459e2d517135993ed047251304ffb0acdedd2ad75ab361fcc3`。
程序始终绑定上述源码提交；后续验收文档提交不改变运行二进制。

主线 `codex/v04-restart-receipt-observability` 从 `de0e77a` 快进至 `83db41d`，
包含供体 `f3ad77b` 及全部原生接线提交，见[合并回执](data/2026-09-09-native-stage46/runtime/MERGE_VERIFIED.json)。
之后才停止旧暂停服务；130 个所属进程全部退出，两把 MPI 锁释放，原输出保留，
见[停止记录](data/2026-09-09-native-stage46/runtime/OLD_STOPPED.json)。
旧完整健康前缀为 10577；独立、完整的恢复权威为 10510，未落盘尾段由新程序重算。

圆柱保留原单相配置、网格、dt、容差和输出策略。封存的 10500、10510 各 133 个文件
再次通过 SHA-256；10500→10510 的 128-rank 短段仅运行一次，十步均为 BDF2、
一次接受、零重试，运行证据校验和完整求解器观测通过。
[短段回执](data/2026-09-09-native-stage46/runtime/PILOT_ACCEPTED.json)及
[原主线对照](data/2026-09-09-native-stage46/runtime/CONTROL_COMPARISON.json)确认：
260 个检查点、统计及 Visit 载荷逐字节一致，物理与稳定修订诊断 CSV 一致。
进程地址绑定的力证书哈希允许不同，原值完整保留；时间数据不作性能提升证据。

新长算于 2026-09-09 04:41:08 +08:00 启动，从 10510 续算 24490 步至 35000。
服务为 `hundun-re3900-coupled-long-20260909.service`；运行目录为
`/home/wyf/code_dev/.benchmarks/hundun-piso-simple-product-20260903/trial-D0p02-zpi2-52/long-coupled-35000-20260909`。
统计 epoch=7000、采样起点=17001 保留；方法签名兼容，不执行恢复降阶或统计重置。


[启动验收](data/2026-09-09-native-stage46/runtime/STARTUP_VERIFIED.json)通过：
健康前缀 10511–10522 共 12 步均为 BDF2、无重试，终端审计全部满足原容差；
连续性最大值 `2.8302423439211935e-7`、能量最大值 `9.296890369296903e-7`。
128 个 rank 全部属于新服务，实际 `/proc/PID/exe` 指向冻结 runner，无暂停或僵尸进程。
重启与冻结程序哈希复核未变；原始异步快照保留。观测流比健康 CSV 多完成一步，
验证了 13 步，唯一提示是仍在计算的下一步 10524 未完整，符合部分快照合同。

至此，本次原生耦合接线、主线合并和新版本圆柱续算切换完成。长算继续运行，
尚未到达 35000；不据此声称长期稳定性、网格独立性或湍流统计验收完成。

# 产品 TCR 历史与重启记录

当前产品入口是 `reaction.model=esf_tpdf` 下的 `ensemble.tcr`。本增量支持
周期、原生物理边界及静态 IBM、backward_euler 的 `shadow` 和 `experimental`；`validated`
仍因缺少科学验收证据而拒绝。产品接线验收见[原生组合验收](../verification/2026-09-09-native-stage46-coupling.md)。

## 统计与反馈的时间位置

每次尝试从已接受的 MeanState 与所有随机场取得 PH 状态；P8 组合先沉积本步气液交换，再从源后状态查询，以同一气相后端查询
净组分质量生成速率。进度速率为 `sum(progress_weights[s] * omega[s]) / rho`，
MeanState 的查询给出 PSR 参照速率，随机场速率等权平均。
反应物名称解析到已冻结的物种顺序，计算 `eta=sum(X_reactants)` 和
`R=mean(progress_rate_field)/progress_rate_PSR`。

本步对流、共同扩散、随机输运和 IEM 使用这个已接受状态的 TCR 控制量。
该显式时间位置写入 method history signature。TCR 只通过
`exp(-cbrt(kappa)*dt/(2*tau_mix))` 改变混合；两个连续化学半区间不乘 κ。
最后仍由 MeanState 的守恒方程与压力闭合决定气相状态，再重整随机场均值。

| 模式/事件 | 产品行为 |
| --- | --- |
| shadow，统计和根有效 | 更新观察到的分支历史，IEM 控制量保持 1 |
| shadow，弱速率或不可用根 | 保留上次有效历史，推进已接受步的修订号，IEM 控制量保持 1 |
| experimental，根有效 | 使用该分支 κ 改变 IEM，并暂存新历史 |
| experimental，弱速率、缺少分支、未解决折点 | 拒绝整个尝试，不改已接受状态 |
| 查询失败或非守恒/非有限化学速率 | 两种模式均拒绝整个尝试 |

分支首次选择遵循显式 `initialization_sign`；存在唯一容许根时允许 0。
后续保持已接受分支，不自动翻转。原生 `ProductCouplingBindings::tcr_fold` 可借用
`ProductTcrFoldProvider`，将显式 source-control continuation 证据传给 P3 FoldEvent。
提供者是对不可变证据的无分配、确定性查询，不维护独立的已接受时钟或分支。
其非零内容/算法身份加入产品与 Restart 合同，并在每次查询前后核对。

查询携带全局单元、已接受历史修订、本次输入修订及旧/新 eta/R；返回事件必须匹配
单元、两种修订和证据源身份，再由 P3 校验固定 eta、实际折点及显式离开分支。
迟到、错单元、非折点、改动身份或查询失败均拒绝整个尝试。只有气相共同接受后才
发布 fold count、last fold、分支和修订；V4/V5 恢复继续验证同一证据源身份。
CLI 没有内置自动折点识别器，未提供证据时沿原分支求解，未解决折点仍明确拒绝。
这个入口传递已有证据，不生成科学证据，不把实验反馈升级为 validated。

原生 1/2/4-rank 制造解验收使用 eta=1/4、R 从 1 到 4/3 再回到 1 的显式路径，
检验上分支 kappa=1 的 IEM 方差、六类错误证据的整步回退，以及真实磁盘恢复和
再次推进。该制造解是接口/事务证据，不是物理折点的自动识别验证。

## 提交与回退

TCR 历史及其编码字节在准备阶段预分配 accepted/trial 两套存储。每个单元计算完候选
根后先调用可移植 accept 校验，但只写 trial。所有单元化学完成后才封存候选；气相的
原生时间步接受后，以无分配的交换共同发布。压力、化学或最后重整失败都丢弃候选。
input revision 每个已接受步加一，重试不推进它，不用尝试编号替代历史身份。

恢复时先检查全部单元的代数、映射身份、初始分支、步号和编码，再参与 MPI 统一判定。
只有原生气相与通量恢复完成后才发布 TCR 历史。

## Restart V4 记录 ABI

每单元 120 字节，x 最快排序，所有整数及 IEEE FP64 位模式按小端编码。
记录身份绑定反应/ESF/TCR 配置与 `TCRHIS` 版本 1。整数从不经过浮点转换。

| 字节偏移 | 内容 |
| --- | --- |
| 0, 8, 16 | accepted step: u64；input revision: u64；algorithm version: u32 |
| 20 | initialized: u32，0 或 1 |
| 24, 32 | η、R：FP64 |
| 40, 48 | signed root、κ：FP64 |
| 56, 60 | 当前/初始分支 sign+1：u32，0/1/2 分别代表 -1/0/+1 |
| 64, 72 | 统计映射身份、fold count：u64 |
| 80, 88 | 最后折点的 η、R：FP64 |
| 96, 104, 112 | 最后折点 base revision 的 step/input/algorithm |
| 116 | 保留 u32，必须为 0 |

无折点时，最后折点输入为 0，base revision 为 `{0,0,1}`。读入时强制该规范编码。
通用 RestartReader 负责完整性、内存预算与分区搬运；ProductDriver 负责模型语义校验。
直接构造 RestartImage 也不能绕过产品校验。

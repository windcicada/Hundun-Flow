# 关键问题复核与分项处理

基线：`792f32efdfa71e2726fd0fc94513ee6f58ed0c03`。本轮只修改开发工作区；
`method-frozen-periodic-20260906` 的程序和输入保持独立、只读，不替换程序，
不修改物性、网格、dt、残差门槛、refinement 上限或 checkpoint 的持久化要求。
原定不打断长测，但实际发生主机重启，详见第 4 节。此后遵照用户新增的热负荷约束：
短测前先确认长测已经暂停，计算任务串行运行，总用核不超过 128；不自动恢复长测。

验证沿用用户确认的真实 runner CLI、ProductDriver、RestartReader/Writer 和公开数值入口。
每个缺陷先做存在性核对，再以单次配置回归区分实际复现与静态风险。
正确性修复、性能观测及未进入默认路径的算法实验分别提交。以下状态随工作更新。

| 项目 | 基线核查 | 实际证据 / 当前边界 |
|---|---|---|
| Visit 前路径构造 | 仍存在；Writer 内部已有保护，不重复改 Writer | rank 0、rank 1 分别注入真实分配失败均退出 134；栈落在 filesystem `operator/`，不是 Writer 内部 |
| 方法恢复与物理历史 | 仍通过改写 `backward_euler_recovery` 选择恢复；质量目标也随之重置 | 已复现并修复；非均匀闭域、签名兼容性、速率重建回归通过 |
| 方法恢复统计隔离 | 仍加载旧 accumulator，只排除第一个 BE 恢复步 | 已复现并修复；1/2/4 rank 非零样本 CLI 回归通过 |
| 新版性能与动量诊断 | multidot 已实现；最终方程/收支观测已存在 | 已补逐 loop/attempt 明细、归一化和区域分布；128-rank 单步观测通过，非正式速度验收 |
| 标量严格守恒 | 已记录有限 dt 的库存缺口；当前生产路径未配对校正最终标量通量 | 建立独立实验回归，不以缩放或归一化作为修复，不混入生产默认值 |

## 1. Visit：根因与验证

基线命令使用冻结 runner，`--visit-interval 1 --observe-performance`，
16³ 静止 IBM 小网格，2 ranks。测试通过现有弱 PMPI 观测入口在推进返回后
启动真实 `operator new` 分配失败扫描，不替换 Writer 或数值内核。
rank 0 和 rank 1 的第 2 个分配点（从 0 计）均触发未捕获 `bad_alloc` / SIGABRT，
退出 134；`addr2line` 确认为 `filesystem::operator/` → runner `run` → `main`。
前两个点是前置探针容器分配，已有本地阶段保护，可一致返回 7，不能把它们算成 Visit 缺陷。

修复将 Visit 目录在循环前的纯本地阶段构造，全 rank 汇总成功后复用；输出错误单独
记录已接受步与当前提交步/时间，并退出，不重试 advance、不回退已提交状态。
CLI 检查所有 rank 正常进入 MPI_Finalize、确实生成 Visit 索引、错误阶段以及提交步身份。
原 checkpoint 分配测试继续保留，不再把关闭 Visit 的测试当作此路径的证据。

早期测试脚本错误假定 VTK 扩展名为 vtu/vtr，而本均匀网格输出为 vti；已改为验证
公开 `.visit` 索引及块数。这两次夹具失败不计为程序 RED。

GREEN：修复版 2-rank CLI 对 rank 0 的 32 个、rank 1 的 60 个真实分配点逐一
注入一次，全部一致结束，无 SIGABRT/挂起。其中分别有 20/12 次 Writer 错误，
已接受步与错误后公开提交步均为 1，时间相同；所有成功输出与无故障参考逐字节相同。
非零 rank 的扫描窗口还包含其随后 checkpoint 阶段，不能把全部 60 点都称为 Visit 内分配。
日志和逐调用结果见 `2026-09-06-key-followup-evidence/visit-*`。

## 2. 方法恢复：源文件事实与消费者策略分开

实际 RED：非均匀周期闭域波动算例保存的目标质量为 1.1030949030378103，
旧恢复路径把它改为当前库存 1.1030948920068622。完整 V2 的旧 CLI 还可以不检查
方法身份便精确续算。恢复策略改为 `RestartHistoryPolicy`，不再改写输入 image。
完整源的质量目标始终保留；只有真正缺少历史的 V1 才使用冷启动库存定义。
主动恢复重建当前算子的速率，只在精确续算分支复制历史速率。

增加独立的语义 method/history signature，覆盖时间历史解释、热/动量速率、
周期通量、C1 联合目标与 SIMPLE 通量语义，不绑定 Git SHA。配置仍由原 plan/schema
验证。ProductDriver 写入带签名 V3；V1/V2 仍完整读取并校验全部源块与清单。
无签名的完整 V2 为“兼容性未知”，精确续算拒绝 10213，须显式方法恢复；V3
不同签名同样拒绝精确续算。V3 扩展了 common header，旧二进制不能读取 V3，
没有改动 fsync、原子发布、保留代数或 staging 预算合同。后续改变速率/离散语义时
必须更新对应签名成分；这不是自动识别任意源码变化的证明。

回归最初人为加入 1e-8 的相对质量差，仅放宽终端检查。这超过未改动的内部
predictor `ClosedMassPlan::certify_density_fields` 1e-12 容差，推进按设计返回 823。
该失败保留在 `method-history-diagnostic.log`，不能算重启算法失败。最终夹具撤回
容差改动，用 1e-13 的合法非零差异验证：目标位级保留，精确续算 BDF2，方法恢复
BE→BDF2。把源速率全部加/减 1e6 后，主动恢复得到的当前/前级速率与未污染源完全
相同，证明没有随后覆盖回旧速率。

新的 run_start.history 将源格式、源/目标签名和策略写入 evidence；验证器同时
核对被哈希绑定的源清单，不把 V2/V3 主动 BE 恢复谎称为源历史缺失。旧证据仍按
原格式验证。真实 CLI 覆盖损坏源块及篡改 history 元数据的拒绝。
暂不开放“旧 MG 存储迁移 + 主动方法恢复”的组合，保留既有约束。

## 3. 方法恢复：独立统计 epoch

实际 RED：源 checkpoint 有 2 个样本；旧方法恢复开始仍有 2 个，推进两步后
变为 3 个，旧方法样本与新方法混用。新路径先读取和验证旧附件，再清空运行内
累加器，原附件不修改。V2 accumulator 记录 epoch 起点、发展时长、最早采样步、
reset 原因、丢弃样本数及源清单 SHA256；继续支持旧 V1 附件。

默认新发展时长沿用 statistics spec 的 development_steps，但相对于恢复点计数。
`--restart-development-steps N` 仅在主动方法恢复时可显式覆盖；即便 N=0 也不采集
第一个 BE 恢复步。没有把一个 BE 步当成重新充分发展，窗口结束也不是湍流统计
已收敛的断言。绝对 collection_end_step 不变。

实际 CLI（1/2/4 ranks）各执行一次：源第 5 步有 2 样本；精确续算两步变为 4；
主动恢复两步仍为 0，新 epoch=5、采样最早 step=9；从恢复后第 7 步精确重启，
继续到第 10 步恰好 2 样本。显式发展时长 1 的覆盖也通过。源目录所有文件 SHA256
前后一致。新证据验证、三类篡改拒绝和源块损坏拒绝均通过。

本组实测：`restart-epoch-ctest-first.log`，10/10 PASS（77.70 s），包括
ProductDriver 方法恢复、1/2/4 rank 存储重启、4-rank V1/V2/V3 I/O 与 runner 自检。

## 4. 最终方法的性能与动量观测

### 运行边界与实际中断

主机在 14:30:31 重启。原冻结长测临时服务随之消失，health 最后一条完整记录为
1225，第 1226 条截断，最后完整 checkpoint 为 1000。没有向原长测发送停止信号，
也没有更换其程序/输入。第一次 128-rank 诊断与长测重叠且未完成，不作为性能证据。
用户随后明确热负荷限制：测试前暂停长测，计算总用核不超过 128。本轮后续测试
均在原长测停止时串行运行。不能把重启后 systemd 的 not-found/默认 success 当成
算例正常结束，也没有内核证据将该中断归因于数值发散。

有效诊断位于 `key-followup-observation-serial-20260906`：原冻结 Re3900 case/spec
只读复用，128 ranks，固定原 dt、物性和残差/refinement 设置，SIMPLE 冷启动 1 步。
Visit 关闭，runner 的最终 checkpoint 仍按原持久化合同写出。观察程序 SHA256：
`1473e092407392baa2cdaa3f7221b1c91c1efa7de204e12d9c149605b44a7e07`。
本次推进最大 rank 10.6833 s、完整步 20.1409 s（含终点 checkpoint 9.4458 s），
进程启动到结束 40.66 s。不能把这一步的 checkpoint 成本套到每个正常步。

先前原冻结程序第 966–1065 步的 100 步观测：平均完整步 9.1270 s，rank 平均
Krylov 3.6825 s、候选 2.0610 s；该版本已包含周期通量与 C1 修正，但当时有
低核数开发构建/回归并行，不能作为新热负荷约束下的独占性能验收。也不能用
冷启动单步与这段 BDF2 窗口计算加速比。

### 实测分层明细

以下时间为同一次单步诊断中各 rank 的平均值；只有一个 attempt。候选数为
baseline/extrapolation/ladder，线性收缩为 final_true/initial_true；完整原始分组
及连续性/能量收缩见 `re3900-loop-observation.json`。

| 修正 | 类型 | 迭代 | A/M 次数 | 候选 | 线性收缩 | solve (s) | A/M (s) | 候选 (s) |
|---|---|---:|---|---|---:|---:|---|---:|
| C1 | spatial | 41 | 58/41 | 1/0/1 | 7.29e-5 | 1.8267 | .9893/.4608 | .3227 |
| C2 | diagonal | 50 | 73/50 | 1/0/1 | 7.99e-6 | .9659 | .2583/.5908 | .2941 |
| C2 ref 1 | diagonal | 46 | 67/46 | 1/0/1 | 9.70e-7 | .8730 | .2366/.5060 | .3082 |
| C2 ref 2 | diagonal | 40 | 57/40 | 1/1/0 | 6.95e-6 | .8087 | .2027/.4344 | .2916 |
| C2 ref 3 | diagonal | 36 | 51/36 | 1/1/1 | 2.29e-5 | .6931 | .1819/.4001 | .4559 |
| C2 ref 4 | diagonal | 22 | 30/22 | 1/1/0 | 1.91e-4 | .4626 | .1083/.2395 | .3055 |

MG refill 合计约 .0530 s，其中复制约 .00776 s；不是当前 M apply 2.6316 s 的
主要成本。structured Halo Waitall .4798 s、控制归约 .3561 s 是部分通信成本，
不包括 IBM donor/MG 独立通信，不得声称是全部 MPI。额外 PMPI 记录显示
ReductionEngine packet Allreduce 2687 次/rank、约 1.1954 s/rank；普通 Halo
控制 Allreduce 1600 次/rank、约 .3544 s/rank；MG 粗层 Allgatherv 1410 次/rank、
约 .1348 s/rank。这些时间已嵌套在 A/M、候选等项内，不能再相加。

新增最终动量装配 .19944 s、终端指标 .00495 s、物理边界账本 .02482 s，分别计时。
新旧同一冷启动步 health 的所有共同数值字段一致（计时除外），conservation 的
共同字段全部字符串相同；力的数值一致，存储/状态证书 token 因进程绑定不同而变化。
观测验收 5/5 PASS，另 SIMPLE 时间收敛 1/1 PASS。

每个 advance 的详细记录使用固定 64-loop 存储，不随运行增长；超出时计数
`dropped_loops`，汇总器拒绝作完整归因，但不改数值成功/失败决定。该上限不是
refinement 上限。极端多次重试的完整逐 loop 外存记录仍是剩余工作。

### 归一化与下一项最小优化

动量归一化定义为 `abs(R_m)/(a0*rho*V*U_rms)`，`U_rms=sqrt(2*K/M)`，是
时间惯性力参考，不是线性求解残差，也不是完整矩阵后向误差。零速度参考时
normalization_valid=false，不伪造速度下限。本步 U_rms=2.91418 m/s，三分量
归一化最大值 .0648002/.0298098/.000753814；最差 gid=2870157/2869716/419169。
IBM 相邻区按 EB control-link 流体行定义，6500 单元；其他活跃单元 5987852。
X 分量邻近/内部 RMS=.0117354/.000279972。这里只描述局部分布，不新增经验门槛。

不再把 multidot、MG 已有 retained-defect 融合作为待办。依据新数据，下一项最小
候选限定为 MG W-cycle 同一 coarse RHS revision 再访问时复用已收集 RHS，避免重复
pack/Allgatherv。源码已有 continue_solution，但 RHS 仍重新收集；需先验证 RHS
revision、全 rank 缓存有效性和错误恢复分支，不能据 solution revision 单独跳过通信。
其收益上限只是上述 .1348 s 的一部分，不能声称能消除 M apply 的全部成本。
本轮仅选择并界定该候选，未把未经 A/B 与 MPI 故障验收的缓存接入生产。

## 5. 标量严格守恒：保留 RED 的独立项目

新增手动验收目标 `v04_product_scalar_contract_experiment`，默认不构建、不作为
绿色 CTest。实际运行 12 条单次轨迹：两类标量 × 均匀/非均匀 × 三个 dt。
正常推进、常量保持、有界性、组成/EOS 和约二阶自收敛通过；非均匀被动标量与
EOS 组分在粗 dt 的最大相对库存缺口分别 1.14038e-6、9.01814e-7，严格守恒失败。
没有改生产算法，没有做密度比缩放/全局归一化。

具体测试矩阵、真实失败证据、最终密度与标量面通量的配套校正设计及开放验收
见 [独立项目](2026-09-06-scalar-conservation-project.md)。该项仍未修复，不应因
Re3900 当前没有输运标量便宣称整个程序已经严格守恒。

## 6. 集成后的最后一轮验证

在长测停止、没有其他 MPI 作业的条件下，聚焦 CTest 以 `-j1` 串行执行，
13/13 PASS（50.65 s）。覆盖 runner 自检、2-rank 统计 epoch CLI、evidence workflow、
1/2/4-rank ProductDriver 存储重启与压力—能量重试、4-rank restart I/O、
PISO/SIMPLE 时间收敛及方法历史恢复。日志：
`2026-09-06-key-followup-evidence/final-focused-regressions.log`。

重试回归额外核对实际被拒绝的第一次 attempt 与随后接受的第二次 attempt：
1/2/4 ranks 均保留 14/2 条 loop 明细，分别绑定原 dt/减半 dt 和拒绝/接受状态，
`dropped_loops=0`。原有回退状态与直接小步推进的一致性检查继续通过。
对应原始测试输出摘录见 `final-retry-observation.log`；没有为观测而改变重试策略。

最终集成 runner 再对 Visit 开启时的 rank 0、rank 1 分别定点注入一次分配失败
（2 ranks、`--only-index 2`），两项通过；各自包含无故障参考运行。
见 `final-visit-rank-0.log`、`final-visit-rank-1.log`。这次定点复核不冒充重新完成
前述 32/60 个点的全扫描，也不把关闭 Visit 的 checkpoint 测试计入该证据。

独立标量实验的严格守恒失败不在 13 项绿色测试内，仍是明确未解决的验收项。
本轮没有完成全库测试、发布验收或新旧 COAST 的独占长窗口性能对比；MG RHS 缓存
只是下一项有测量依据的优化候选，尚未实施。最后确认原长测服务 inactive/not-found，
没有运行中的 mpirun，冻结程序与第 1000 步 checkpoint 保持原位。

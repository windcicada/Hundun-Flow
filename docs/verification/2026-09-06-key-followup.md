# 关键问题复核与分项处理

基线：`792f32efdfa71e2726fd0fc94513ee6f58ed0c03`。本轮只修改开发工作区；
`method-frozen-periodic-20260906` 及其 128-rank 长测保持独立、只读，不停止、
不替换程序，不修改物性、网格、dt、残差门槛、refinement 上限或 checkpoint 的持久化要求。

验证沿用用户确认的真实 runner CLI、ProductDriver、RestartReader/Writer 和公开数值入口。
每个缺陷先做存在性核对，再以单次配置回归区分实际复现与静态风险。
正确性修复、性能观测及未进入默认路径的算法实验分别提交。以下状态随工作更新。

| 项目 | 基线核查 | 实际证据 / 当前边界 |
|---|---|---|
| Visit 前路径构造 | 仍存在；Writer 内部已有保护，不重复改 Writer | rank 0、rank 1 分别注入真实分配失败均退出 134；栈落在 filesystem `operator/`，不是 Writer 内部 |
| 方法恢复与物理历史 | 仍通过改写 `backward_euler_recovery` 选择恢复；质量目标也随之重置 | 已复现并修复；非均匀闭域、签名兼容性、速率重建回归通过 |
| 方法恢复统计隔离 | 仍加载旧 accumulator，只排除第一个 BE 恢复步 | 已复现并修复；1/2/4 rank 非零样本 CLI 回归通过 |
| 新版性能与动量诊断 | multidot 已实现；最终方程/收支观测已存在 | 不复用旧错误周期/C1 算法的耗时占比；新冻结程序数据和新增观测待整理 |
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

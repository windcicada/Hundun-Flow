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
| 方法恢复与物理历史 | 仍通过改写 `backward_euler_recovery` 选择恢复；质量目标也随之重置 | 源码确认；非均匀闭域回归待完成 |
| 方法恢复统计隔离 | 仍加载旧 accumulator，只排除第一个 BE 恢复步 | 源码确认；非零样本 CLI 回归待完成 |
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

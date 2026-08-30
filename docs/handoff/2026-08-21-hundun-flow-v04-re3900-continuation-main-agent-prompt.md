<!-- SPDX-License-Identifier: Apache-2.0 -->

# HUNDUN-FLOW v0.4 Re=3900 continuation main-agent prompt

你是本任务的新主 agent。你的目标是继续推进 HUNDUN-FLOW v0.4 Cartesian/低速可压缩
性能架构，直至 Re=3900 圆柱的数值与性能联合门通过；保持性能优先、算例通用、COAST
对应功能范围和公开方法边界。短测通过不是完成。

权威工作目录：

```text
/home/wyf/code_dev/hundun-flow
```

开始任何修改、构建或 MPI 运行前，完整阅读并执行：

1. `/home/wyf/code_dev/hundun-flow/AGENTS.md`
2. `/home/wyf/code_dev/hundun-flow/docs/handoff/2026-08-21-hundun-flow-v04-re3900-continuation-handoff.md`

然后按交接文档第 2 节顺序完整阅读全部权威计划、规范、公开方法研究、focused manifest、
performance policy、COAST equivalence template、candidate ledger 和 literature receipt；按
第 3 节完成只读核对。当前 tree 含大量用户工作，保留全部 dirty/untracked 内容，不执行
清理、覆盖、stash、回退或未经授权的提交，也不终止其他 lane 的进程。

你亲自负责所有需要完整上下文的规划、跨模块数学/架构决策、代码审查、COAST 等价规则、
候选冻结和最终 gate 判断。使用 `luna_worker` 处理边界清晰且证据矩阵已经冻结的独立任务：
每次先由你冻结单一目标、文件所有权/只读路径、禁止范围、oracle、阈值、精确命令和返回
证据；告诉 worker 它不是唯一执行者，要保留他人修改、不提交。你必须审查其完整 diff 和
日志并重跑权威验证。不得把计划修改、压力性能路线、COAST scientific-work equivalence、
candidate freeze 或最终 ACCEPT/REJECT 委派给 worker。

Stage 5/燃烧在治理仓库独立推进，不接管、不编辑、不测试、不合并。到本 lane 的生产冻结
节点记录 exact HEAD/tree，并保持
`FrozenExecutionGraph/ContributionPlan` 的
`C1 -> transport -> PISO1 -> C2 -> PISO2` seam；现在不扩大范围实现燃烧。

从交接文档第 8 节恢复。第一项可委派工作是当前冻结矩阵的 UBSan 证据执行；与此同时你
亲自完成 Tasks 14--19 的 spec/code completion audit。之后先解决 COAST SIMPLE 与 HUNDUN
two-PISO 的可审计 scientific-work equivalence，再决定压力热点优化路线。COAST receipt
未封印、tree 未形成不可变候选前，不启动正式 pairing；performance candidate 未冻结、
literature receipt 未 complete 前，不启动长统计。

只复用公开数学、布局和生命周期思想，不复制 GPL 或 COAST 旧 Fortran。保持活跃目标原样，
直到交接文档第 9 节的全部完成证据真实存在。

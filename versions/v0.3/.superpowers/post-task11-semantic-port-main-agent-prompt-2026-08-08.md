# Main-Agent Handoff Prompt: Post-Task-11 Semantic Port

请在以下工作树接手 HUNDUN-FLOW：

```text
/home/wyf/code_dev/.worktrees/hundun-flow-stage3
```

首先完整阅读：

```text
docs/superpowers/specs/2026-08-08-hundun-flow-post-task11-semantic-port-architecture-design.md
```

这是用户已确认的 Task 11 之后架构方向。核心决定是：

- 保持 HUNDUN-FLOW 自有架构、命名、ABI、数据布局、MPI、execution、linear
  algebra、transaction、checkpoint 和 diagnostics；
- 不引入 AMReX/OpenFOAM/Basilisk/IncFlo 运行时依赖；
- 只做语义移植：数学与行为规格 -> mutation-sensitive RED -> HUNDUN 风格独立实现；
- 禁止复制、翻译、机械改写上游源码、注释、控制流、宏、错误消息和 ABI；
- AMReX/IncFlo 虽为 BSD-3，本阶段也只读参考；OpenFOAM/Basilisk/NASSLARD2D
  为 GPL，只能用作数学/行为参考；无明确许可证的 ghost-cell-IBM 与 3D-NS-FSI
  禁止复制；
- 不添加 public `REUSE_MANIFEST.yaml`，不修改冻结的 `NOTICE` 内容；
- 永远不运行 96-cubed。

执行边界：

1. 先只读核对 Task 11 是否已正式 accepted，不能因为单个 near-wall 修复通过就
   推断 Task 11 完成。
2. 在 Task 11 accepted 之前，不进入 Task 12 产品实现，也不修改 Task 11
   数值路径来迎合新架构。
3. Task 11 accepted 后，先修改权威 Stage 3 计划，保留历史计划和交接记录，说明
   任务重排理由；不要直接开始编码。
4. 建议把原 Task 19 拆为 19A/19B/19C，并把 19A constant-density IBM driver
   MVP 提前到 Task 11 后；把 Task 17 拆为 IBM-only codec/continuation 与后续
   combined continuation；把 Task 18 拆为最小 diagnostics/counter substrate 与
   完整 adapters。
5. 第一可用 MVP 必须包含 constant-density static LFP-GCIBM、两次 PISO、force、
   retry/rollback、1/2/4-rank、同一 `hundun` executable、Checkpoint v3 IBM-only
   continuation 和最小 diagnostics；不等待 WALE/三密度模型全部完成。
6. 随后执行 Task 12 WALE core、Task 13 body-fitted WALE、Task 14 material IBM、
   Task 15 ideal-gas IBM、Task 16 combined gate，再完成 checkpoint/diagnostics/
   driver/performance/final acceptance。

AMReX 语义移植优先级：

```text
P0  geometry facts 与 operator/flow composition 分层
P0  regular/interface/inactive 稀疏工作分区
P0  static plan 与 attempt-local values 分离
P0  exact counters、allocation-free apply、packed donor halo
P1  有 RED 后才实现 ConservativeFluxCorrectionPlan
P1  Stage 3 完成后规划 project-owned geometric multigrid
P2  AMR/GPU/moving-body/FSI 后续独立设计
```

执行方法：

- 主 agent 负责全局计划、数学判断、完整 diff 审查、版权独立性审查和最终验收；
- 优先使用已配置的 `ds_worker` 或等价配置处理边界明确的单模块任务；worker 不是
  单独的科学裁决者；
- 遵守 `AGENTS.md` 的 worker 数量、顺序和 requirement/code-quality review；
- 使用 codegraphf 检查符号、调用方和影响范围，使用 `rg` 精确搜索；
- 所有实现遵循 数学/行为规格 -> mutation-sensitive RED -> 最小实现 -> fast ->
  screen -> acceptance；
- 开发期使用 12/24 或等价低成本案例，稳定后才运行 12/24/48 acceptance；
- 长任务绑定 exact HEAD、dirty-diff hash、binary SHA-256、命令、环境、日志 hash
  和退出状态；
- 不同时运行多个会争用 CPU/内存/MPI 的大型矩阵；
- 不 push、不发布、不修改私有 COAST/BOFFIN 或研究数据/进程。

你的第一份回复应给用户：

1. Task 11 exact acceptance 状态与剩余 blocker；
2. 新设计与现有 Stage 3 Tasks 12--21 的差异表；
3. 建议的权威计划文件修改白名单；
4. 19A/17A/18A 的接口与依赖；
5. 哪些工作可由 bounded ds_worker 完成，哪些必须由主 agent 完成；
6. fast/screen/acceptance 调度与资源分组；
7. 版权独立性和上游 reference 记录方式；
8. 修改产品代码前的前三项工作。

不要假定本提示词覆盖工作区现状；必须先只读检查 HEAD、dirty/untracked、后台
测试、日志、当前交接和刚完成的 Task 11 证据。保护所有现有修改。


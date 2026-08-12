# HUNDUN-FLOW Stage 3 科学框架优先执行修订

> **2026-08-09 接续说明：** 本文件保留 Task 11 科学优先修订的历史和未被替代的
> 正确性边界。Task 11 接受后的仓库拆分、Tasks 12--21 顺序、紧凑测试分层及
> 最终验收 owner 由 `2026-08-09-hundun-flow-stage3-framework-completion.md`
> 接续控制。

状态：2026-08-05 经用户确认，立即控制 Stage 3 后续执行顺序、测试成本和
资源调度。它保留全部既有交接、拒绝、审查和证据记录，不追认 Task 11，
也不授权 Stage 4、发布、push 或私有参考源码访问。

## 1. 权威边界

本修订按以下优先级解释：

1. 2026-08-05 用户关于科学闭环、效率原则和计划重排的最新指令；
2. 冻结的 Stage 3 科学规格和 2026-08-05 force-consistency authority
   addendum；
3. 本执行修订；
4. 2026-07-27 Stage 3 实施计划中未被本修订明确替代的任务边界；
5. coordinator ledger、task brief、readiness 和历史交接记录。

本修订只替代下列执行条款：

- Task 11 后必须先执行 G1、G2、G3 才能进入 Task 12 的顺序；
- 每次 repair 机械重复完整 Debug、Release、ASan、UBSan 和 96³ 的做法；
- 所有数值网格只能用完整三层 acceptance 反馈开发方向的做法；
- 长时间测试只绑定 HEAD、却不绑定 dirty tree、二进制和 runner 退出状态的
  不完整证据身份。
- 旧计划将 Task 11 全部工作标为 main-agent-only 的 blanket rule：T11-S、跨模块
  判断、完整 diff 和最终 verdict 仍由主 agent 独占，但冻结后的 T11-M/T11-E
  局部 fixture、runner 和机械 repair 可交给一个 bounded worker；
- Task 11 中 strict ordinary-host zero-allocation、replicated-vector elimination
  和 peer-only/extreme-scale tuning 的原时点；它们仅按第 4.4/4.5 节条件移到
  Task 20，不删除。
- 旧 Task 20 仅允许 docs/tests/CMake 的 allowed-file 边界；它只按第 4.5 节
  measurement-RED-backed performance-repair subcluster 做窄扩展；
- 旧 Task 20 Step 5 与 Task 21 重复运行完整 Stage 3 acceptance inventory 的
  双重 owner；本修订由 Task 21 唯一拥有最终 exact-HEAD 全 inventory。

就后续顺序而言，它明确替代以下历史文件中的对应条款，而不改写其余证据：

- `.superpowers/sdd/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale/stage3-code-governance-amendment-readiness.md`
  的固定 `G1 -> G2 -> G3 -> Task 12` 顺序；
- `.superpowers/task-11-in-progress-handoff-2026-08-05.md` 的 post-Task-11
  future-order 指令，以及 Task 11 accepted 前不得在隔离 worktree 预备 Task 12
  的旧限制；
- `docs/superpowers/plans/2026-08-04-hundun-flow-task11-a22-a2-rhie-chow-interface-adjoint.md`
  和
  `docs/superpowers/plans/2026-08-04-hundun-flow-task11-a22-adjoint-pressure-coupling.md`
  中
  post-Task-11 先进入 G1/G2/G3 的条款；
- `.superpowers/sdd/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale/governance-g1-readiness.md`、
  同目录 `governance-g2-readiness.md` 和 `governance-g3-readiness.md` 中仅与
  旧 activation predecessor 有关的条款。

G1 readiness 必须在 Task 19 accepted HEAD 上重新生成路径/安装 inventory；旧的
provisional Task 11 identity 不能复用。G3 使用当时 accepted 的路径和名称，不再
等待 G2；未来执行 G2 时必须同步更新程序文档和 protected-name manifest。

它不降低或删除以下科学合同：

- Ghost-Cell/LFP 二阶重构；
- immersed pressure boundary；
- 守恒 Poisson、共享 face flux 和唯一 immersed residual；
- 成功 trial 恰好两次 PISO corrector；
- velocity、pressure、真实表面 pressure/viscous/total force；
- operator/surface pressure、viscous、total consistency；
- transaction rollback、collective failure 和 1/2/4-rank 分解一致性。

Stage 3 specification 第 18.2 节的所有正式二阶行继续要求两个相邻 refinement
segment 都满足原有 `>= 1.8` 门槛。coordinator ledger 中将 signed viscous net
vector 或 operator/surface consistency 降为普通 diagnostic 的记录，不再能
单独授权 Task 11 acceptance。若数学 RED 证明 force 名称、方向或比较运算与
离散残差约定矛盾，必须先建立单独的科学语义修订；不能通过测试侧换号、删除
行或改阈值解决。

2026-08-04 的 A22/A3/R1 规格、计划和证据保持原样。它们是当前 Task 11
候选链的一部分，不是 accepted Stage 3 capability，也不覆盖本节的 hard gate。

## 2. 当前接受边界

- Stage 0、Stage 1、Stage 2 已接受。
- Stage 3 Tasks 1--10 已接受；Task 10 接受提交是
  `0db56e463470dd1a605709ba05d8bd6a900f496b`。
- 当前工作分支 HEAD `28f5dd541a0e3ce9ecf852e53d83981add3a5be8`
  和其 tracked/untracked Task 11 修改均未接受。
- Task 11 阻塞 Gate 3、治理任务、Task 12 合入和后续科学验收。

任何 Task 11 诊断、RED、screen 或长测只说明其绑定候选的性质；它们不能把
dirty worktree、历史 candidate 或单段收敛结果升级为 acceptance。

## 3. 修订后的依赖图

```text
Gate 3
  Task 11-S science authority and second order
  Task 11-M MPI, transaction and collective failure
  Task 11-E evidence, selectors and detached runner
  Task 11-P bounded-memory/performance evidence
  -> one indivisible Task 11 verdict

Gate 4
  Task 12 standalone backend-neutral WALE
  -> Task 13 body-fitted WALE

Gate 5
  Task 14 material-density IBM
  -> Task 15 ideal-gas IBM
  -> Task 16 combined IBM+WALE

Gate 6 framework
  Task 17 Checkpoint v3
  -> Task 18 diagnostics and counters
  -> Task 19 same-executable driver

Post-framework closure
  G1 directory unification
  -> G3 program documentation system
  -> Task 20 performance evidence + capability ledger
  -> Task 21 final acceptance

Deferred maintenance
  G2 internal naming unification
```

Tasks 1--21 保留原编号和原科学边界。Task 11-S/M/E/P 是 acceptance clusters，
不是可独立 accepted 的新产品 tasks。任何一簇失败都使整个 Task 11 rejected。

G1 从原先的 Task 11 后前移门改到 Task 19 后：此时整体求解框架已存在，目录、
安装和导出路径可一次性稳定。G3 在 G1 后独立接受，生成程序文档、安装/运行说明
和非 Python consistency gate；Task 20 随后生成 capability ledger、性能证据并
验证 G3 文档 inventory。这样保留 G3 的独立治理 verdict，也避免与 Task 20 的
产品性能 repair 混成一个 diff。G2 只涉及内部命名，不影响数值正确性、用户接口
或运行能力，登记为 Stage 3 后维护债务，不再阻塞 Task 21。

## 4. Task 11 acceptance clusters

### 4.1 T11-S：科学权威

主 agent 独占数学判断、跨模块实现、完整 diff 审查和最终 verdict。顺序固定为：

```text
normal/stress/residual sign derivation
-> mutation-sensitive one-link RED
-> shared boundary-row authority RED
-> exact-state force/consistency decomposition RED
-> minimal product repair
-> fast
-> screen
-> frozen-candidate acceptance
```

禁止以滤波、额外阻尼、逐 case 调参、增加 PISO corrector、替代数值路径或放宽
阈值关闭 RED。若 RED 不敏感，先修测试；若数学未能唯一决定契约，停止产品修改
并提交显式科学语义修订。

### 4.2 T11-M：MPI 与 transaction

必须关闭：

- local count/extent failure 在 collective 前的统一 preflight；
- 命中产品路径的 donor-Halo underprovisioned fixture；
- 域内 near-periodic donor row；
- 成功/失败 trial 的唯一 bitwise equality authority；
- pressure、momentum、PISO 和 force 在 1/2/4 rank 的一致失败分类；
- rollback 后 committed/history/face flux/inactive canonical zero 不变。

该簇可由 bounded worker 实现已冻结的局部 fixture 或机械修复；主 agent 仍负责
collective 语义、跨模块审查和 exact-HEAD 验收。

### 4.3 T11-E：证据和 runner

必须关闭：

- expensive exact-A22 诊断默认关闭，只由显式 selector 开启；
- Task 11 相关 CTest 标签覆盖所有受影响测试；
- fast/screen/acceptance selectors 只改变规模和证据收集，不改变产品算法；
- detached runner 能可靠记录真实退出状态；
- 每项长证据绑定完整候选身份。

长测试 manifest 至少记录：

```text
accepted base and candidate HEAD
git status --porcelain=v1 and tracked patch SHA-256
untracked evidence manifest SHA-256
build preset/type, compiler, standard library and MPI
binary path, inode and SHA-256
argv, environment allow-list, rank count and resource class
systemd user unit/runner identity, PID tree and start time
log path, final exit status, duration, peak RSS and log SHA-256
```

只记录 HEAD 不能证明 dirty candidate；只记录 sidecar 不能替代 runner 的真实
退出状态。

### 4.4 T11-P：性能与分配

Task 11 必须证明：无泄漏、无 use-after-free、无无界 workspace 增长、diagnostic
关闭时没有与 selector 无关的大型 snapshot 保留、96³ 能在可用资源内完成。

以下项目在满足上述条件后移入 Task 20：普通 host allocation 严格零次、replicated
active vector 消除、peer-only scaling、极端规模性能微调。MPI collective safety、
整数溢出导致的 rank 分歧、Restart/rollback 和用户可运行性不属于可延期性能债。

### 4.5 Deferred-work ledger

| ID | 来源与风险 | 可延期理由 | 重新激活与关闭 |
|---|---|---|---|
| S3-D1 strict ordinary-host zero-allocation | Task 11 R2；风险是 jitter 和吞吐下降 | Task 11 先证明无泄漏、无 UAF、workspace 有界且 96³ 可运行；剩余分配不改变数值路径 | Task 20 performance matrix 开始时激活；由 Task 20 accepted 关闭 |
| S3-D2 replicated active vector elimination | Task 11 R2；风险是内存和强扩展效率 | acceptance rank/grid 内先证明资源有界且 decomposition 数值一致；不涉及 collective 顺序 | Task 20 memory/scaling gate 激活；由 Task 20 accepted 关闭 |
| S3-D3 peer-only/extreme-scale tuning | Task 11 R2；风险是大 rank 性能退化 | 1/2/4-rank correctness、collective safety 和用户运行路径不延期；仅优化通信集合/常数 | Task 20 scaling evidence 激活；由 Task 20 accepted 或 capability ledger 明确不支持该规模后关闭 |
| S3-D4 G2 internal naming | 原治理 G2；风险是内部命名漂移和以后文档同步成本 | public names、schema、CLI、diagnostic IDs 和 serialized fields bytewise 保护；不影响 science/runtime | Task 21 accepted 后的首个维护窗口激活为 `Stage3-M1/G2`，且必须在任何扩大 public surface 的重构前关闭；由独立 G2 acceptance 关闭 |

若 S3-D1/D2 在 Task 11 证明阶段暴露 leak、无界增长、96³ 不可运行或 rank-dependent
行为，立即撤销延期并返回 T11-P blocker。若 S3-D4 触及 public/frozen 名称，超出
maintenance 授权并停止等待用户决定。

旧 Task 20 allowed-file list 只覆盖 evidence/docs/tests。为关闭 S3-D1--D3，
本修订允许 Task 20 在 measurement RED 之后逐项开启 bounded performance-repair
subcluster。每个 subcluster 必须由主 agent 用 codegraphf 冻结最小产品文件清单、
单独 RED、资源上限和回归矩阵；不能用此授权做一般清理，也不能在一个 diff 中
混合 D1/D2/D3。该窄扩展明确替代旧 Task 20 allowed-file list，但不改变 Task 20
main-agent ownership 和完整验收责任。

## 5. 三层数值测试

### 5.1 Fast

目的：在一次实现迭代内证明 mutation、局部代数和方向。

- direct analytic/unit/row tests；
- 12³/24³，或证实能暴露同类错误的更小 case；
- Debug 用于合同、failure 和 bitwise rollback；
- Release 用于小规模数值方向；
- 只运行 changed-symbol 影响到的 1/2/4-rank 小矩阵。

Fast 必须走相同 production reconstruction、operator、PISO 和 force path。测试
可以注入解析 field 或选择网格，但不能维护替代 solver。

### 5.2 Screen

目的：用两层结果判断是否值得支付完整 acceptance 成本。

- 24³/48³ Release；
- 优先顺序：sphere uniform pressure/total consistency、warped prism viscous row、
  translated sphere multi-link/translation；
- 每个 formal row 必须有限、正、严格下降，并满足当前规格门槛；
- signed viscous force 只有在 candidate-independent non-degeneracy preflight
  冻结正式 fixture/分量后才能进入 screen；preflight 失败时停止并先修订 fixture，
  不能用候选误差选择分量；
- 任一关键行失败时停止 96³，回到 RED/根因，不做 case-specific tuning。

### 5.3 Acceptance

只对一个精确冻结候选运行：

- 九个正式 sequence 的 24³/48³/96³；
- 全部 formal rows 的两个相邻 order；
- 1/2/4-rank decomposition、closed transient 和 engineering path；
- complete Debug；
- 按 changed-symbol 影响范围冻结的 Release、ASan、UBSan、tests-off、header、
  policy、provenance、`nm`、`ldd` 矩阵；
- main-agent requirements review、code-quality review 和完整 Task 10..candidate
  diff 审查。

不得把 sanitizer 与每个大型数值 selector 做笛卡尔积；必须在 evidence matrix
中解释每个配置覆盖的缺陷类别。

### 5.4 Stage 3 final-inventory owner

- Task 20 运行 inventory list/contract、performance matrix、deferred-repair
  focused regressions、必要 screen 和该 task 的 complete Debug；不再次串行运行
  全部九序列和完整 Stage 3 acceptance shell。
- Task 21 是最终 exact-HEAD `stage3_acceptance.sh` 和完整 Release inventory 的
  唯一 owner，只运行一次。
- 只有 candidate source closure、binary SHA/inode、argv、rank、environment 和
  log identity 全部相同时才可复用 Task 20 证据；Task 20 后任何相关 tracked
  change 都使受影响证据失效。
- Task 11 的 Gate 3 formal matrix 仍在 Task 11 candidate 上独立执行；它不能
  替代最终集成后的 Task 21 evidence。

## 6. 资源调度

```text
L: compile, headers, policy, non-MPI unit
M: MPI 1/2/4, direct algebra, 12³/24³
H: 48³/96³, long engineering or full acceptance selector
```

- L 类可按实测内存并行，默认构建并发不超过机器安全余量。
- 同一 MPI 资源组内 M 类串行；系统资源充分时可与少量 L 类重叠。
- 任意时刻最多一个 H 类作业；禁止两个 96³/高内存矩阵争用。
- 所有 H 类运行时不得重建其 binary、改其 build tree、覆盖日志或改变环境。
- frozen/formal acceptance H 作业还要求其 source tree 完全不变；长作业期间只在
  独立 worktree 做只读审查、文档或不共享 binary/build tree 的轻量工作。
- 对已明确降级为 diagnostic、且已绑定 binary SHA/inode 的旧 H 作业，允许同一
  worktree 的计划/证据文档改动，但不得修改产品/测试 source 或 build tree，并
  必须记录 launch 后 dirty identity 已变化；这种作业永远不能升级为 acceptance。

Task 11 acceptance 候选冻结且 H 类 formal run 开始后，可在独立 worktree 预备
bounded-worker-eligible 的 Task 12。Task 12 不得使用 Task 11 dirty headers，
不得与 H 类作业争用 MPI/高内存资源，也不得在 Task 11 accepted 前合入 Stage 3
主线。

## 7. 后续任务执行策略

- Task 12/13：先形成独立 WALE 核心，再接 body-fitted flow；不等待 G2 命名。
- Task 14/15：数学推导和 fixture 可并行准备；涉及相同 flow/IBM 组合的产品修改
  顺序合入并由主 agent 完整审查。
- Task 16：保持 IBM、三种密度路径和 WALE 的不可分硬门。
- Task 17/18/19：先完成 transactional checkpoint，再稳定 diagnostics，最后建立
  同一 executable 的用户运行路径。
- G1：只做路径、CMake、install/export 和 policy 必需改动，不改变科学行为。
- G3：G1 accepted 后独立完成程序文档和非 Python consistency gate；它不等待
  G2，也不与性能产品 repair 混在一个 diff。
- Task 20：只在 G3 accepted 后开始，集中完成性能债、capability ledger、
  performance evidence，并把最终 G3 文档 inventory 纳入 contract test。
- Task 21：只在 Tasks 1--20、G1、G3 全部 accepted 且无大型后台作业后开始。

## 8. 计划修改和历史保留规则

- 旧规格、计划、handoff、rejection 和 review 文件不删除、不重命名、不回写成
  “accepted”。
- 旧文件可在顶部增加一条指向本修订的 sequencing notice；原正文保持可读。
- coordinator ledger 只追加新决定，不删除旧结论。
- 每个 deferred item 必须记录来源 task、风险、为何不影响科学正确性/实际使用、
  重新激活 gate 和最终关闭 task。
- 任何影响 convergence、conservation、MPI consistency、Restart reliability 或
  user executability 的 finding 不能标为 deferred。

## 9. 当前立即执行顺序

1. 保留正在运行的旧 sphere-uniform 24/48/96 诊断至自然结束，封存其真实
   `/usr/bin/time` 状态和日志 SHA；它不作为 acceptance。
2. 冻结 Task 11 数学/RED 设计，先使现有候选以正确理由变 RED。
3. 只对证明的根因做最小实现，依次通过 fast、代表性 screen 和完整 acceptance。

在第 2 步的 RED 被主 agent 审查并实际失败前，不再进行新的 Task 11 产品 repair。

## 10. 非承诺性 planning estimate

以下区间从 2026-08-05 takeover 起算，用于资源和依赖规划，不是以牺牲 hard gate
换取的 deadline：

- Task 11：通常 2--4 周；若 force semantic/shared-row RED 触发正式科学修订，
  4--6 周；
- 核心框架至 Task 19：8--14 周；
- 完整 Stage 3（含 G1、G3、Task20、Task21）：11--18 周。

每关闭一个 hard gate 后按实际长测吞吐、remaining RED 和机器可用性更新区间。
任何估计都不授权并发多个 H 类作业、删除科学行或跳过 exact-HEAD acceptance。

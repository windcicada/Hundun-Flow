# HUNDUN-FLOW Stage 4--6 主 agent 交接

日期：2026-08-11

状态：Stage 3 已接受；Stage 4-P0 已融合；accepted-state intake 已完成；Stage 4 产品实现尚未开始。

## 1. 当前唯一开发入口

```text
governance repository=/home/wyf/code_dev/hundun-flow-governance
worktree=/home/wyf/code_dev/.worktrees/hundun-flow-stage4-reacting-flow
branch=coast/stage4-reacting-flow
handoff_parent=0f39799692a5da47e9810a4e2262143fb2466996
handoff_parent_tree=1124e82f63ba61ed1041c9e77218fd08d695f6bb
handoff_commit=由提交后的外部交接结果登记，避免 tracked file 自引用
product_repository=/home/wyf/code_dev/hundun-flow
product_head=22ed17b438ffbb121ccda97898580183bd0803f8
```

`/home/wyf/code_dev/hundun-flow` 是接受后的产品投影，不是开发工作树。Stage 4--6 的源码、
测试、治理和验收只在 governance 仓库推进，直到计划规定的产品投影节点。

## 2. 接受基线与已完成融合

| Identity | Commit / tree | Meaning |
| --- | --- | --- |
| Stage 3 tested code `C` | `0cbd3d5bde4be63bc6346b4b32db771d87c59ea2` / `d50c1236f67bd2bdde58c94a125e530ae0f2ffea` | 冻结的产品、测试和数值权威 |
| Stage 3 governance `G` | `36bebc292e825fa15272481c6a00c2273fa61ce0` / `897560c30d7d7049a81605a257702b4091a13f25` | Stage 4 branch 的第一 parent |
| product `P` | `22ed17b438ffbb121ccda97898580183bd0803f8` / `7fb9ce848238eeab5dc1ad0908092d8d115851b4` | 已接受 `0.2.0` 产品投影 |
| P0 seal | `910fb1f7fc3df2e0c596d3682db06db442c03ccf` | `PREFLIGHT_PASS` 候选输入，不是 Stage 4 接受 |
| integration merge | `d45ef02706a17f12d38e050497f076cc5002fb51` / `50efe6741f4b6b6bbf113ce899edac754ede4d10` | parents=`G,P0`；签署融合提交 |
| intake commit / handoff parent | `0f39799692a5da47e9810a4e2262143fb2466996` / `1124e82f63ba61ed1041c9e77218fd08d695f6bb` | accepted-state intake；无产品变更 |

融合审计结论：23 个变更路径全部属于批准的治理、计划、参考和 P0 证据；Stage 3 产品投影
manifest 的 272 个路径相对 `C` 零差异，`include/src/tests/CMake/presets/cmake/VERSION` 也零
差异。P0 的非冲突 blob 与其 seal 完全一致；唯一冲突 `AGENTS.md` 已显式保留 Stage 3
身份、权威和 P0 限制。融合 diff SHA-256 为
`d450195159fb706235340e5888758d7b73c228ab1cca6ffd3ee579e722342f01`。

完整 intake 证据见
`.superpowers/sdd/stage4-4F-0-baseline-receipt.md`。该 receipt 记录了 Stage 3 receipt、产品
manifest、DCO、linked-worktree、80 个公共头、符号/target inventory、schema v1--v3、
Checkpoint v1--v3、diagnostics、transaction、final-flux、两次 PISO 和当前进程状态。

## 3. 权威文档顺序

开始工作前完整阅读当前 `AGENTS.md`，并按其 Required Reading 顺序读文档。Stage 4--6 的
直接控制层为：

1. `docs/superpowers/specs/2026-08-09-hundun-flow-stage4-6-linux-cpu-v1-architecture-design.md`
2. `docs/references/2026-08-09-hundun-flow-stage4-6-reference-catalog.md`
3. `docs/superpowers/plans/2026-08-09-hundun-flow-stage4-reacting-flow.md`
4. `docs/superpowers/plans/2026-08-09-hundun-flow-stage5-esf-tpdf-tcr.md`
5. `docs/superpowers/plans/2026-08-09-hundun-flow-stage6-spray.md`
6. `docs/superpowers/plans/2026-08-09-hundun-flow-stage4-6-v1-integration.md`
7. `.superpowers/stage4-p0/final-receipt.md`
8. `.superpowers/sdd/stage4-4F-0-baseline-receipt.md`

P0 计划和设计解释候选 artifact 的形成过程；它们不能覆盖正式 Stage 4 task contract。

## 4. 下一步：先收口 Task `4F-0`

accepted-state intake 已完成，但 Task `4F-0` 的产品测试部分尚未开始。下一位主 agent 从
Stage 4 plan 的 `4F-0` Step 2 继续：

1. 以 intake receipt 冻结的真实 authority 为输入，先写 `stage4_source_policy` RED；
2. source policy 必须拒绝 public header 中的 Cantera/SUNDIALS 类型、第二套 plugin ABI、
   私有路径和未登记 Stage 4 产品前缀；
3. 创建 `docs/numerics/stage4-capability-ledger.md`，只添加 Stage 4 计划/未实现行，不改写
   Stage 3 disposition；
4. 如测试注册确需改 `tests/CMakeLists.txt`，由主 agent拥有中央修改并审查；worker不得改；
5. 只运行 `git diff --check`、source-policy fixture 和必要的 tests-off configure；不运行数值
   测试；
6. 完整审查后创建签署提交 `docs: freeze Stage 4 accepted baseline` 的后续实现提交或使用不
   混淆既有 receipt 的明确 subject，并在 receipt 中指出 intake commit 已存在；不得 amend
   已有提交。

完成后按 Stage 4 计划执行 `4F-1`，不要跳到 Cantera adapter 或 reacting flow composition。

## 5. Stage 4 分配与推进

Stage 4 共 27 个任务，默认串行，分为五簇：

```text
4F-0..5  contracts: composition, services, transaction, schema v4, persistence/diagnostics
4P-1..4  package: provenance, artifact, C++/MPI boundary, relocation
4C-1..5  backend: runtime/workspace, thermo, transport, chemistry, 0D/PSR
4R-0..7  reacting: C-T-C/two-PISO proof, state, transport, boundary, IBM/WALE, p0, driver
4A-1..4  acceptance: Checkpoint, diagnostics, compact acceptance, exact-head seal
```

风险顺序不可反转：`4R-0` 必须在生产 reacting step 前冻结 C-T-C 与两次 PISO 的数学闭环；
`4P` 必须把 P0 候选升级为正式 HUNDUN 输入后，`4C` 才能依赖它。每个 task 只运行其计划
中的 mutation RED、focused unit/header/policy、必要 small MPI 和最多一个 12-cubed smoke。

Task 结束由主 agent完成 caller-impact、requirements、quality、完整 task diff、版权独立性和
DCO。worker 只接收一项边界清晰且证据矩阵冻结的工作，不修改中央 registry、CMake root、
科学组合、提交或 sign-off。当前用户规则仍是优先默认 worker，不使用 `luna_worker`，直到
用户再次改变规则。

## 6. Stage 5、Stage 6 和最终集成

默认接受历史严格串行：

```text
Stage 4 accepted (0.3.0)
-> 向用户报告 stage node，提出串行/有限并行建议并等待指示
-> Stage 5 accepted (0.4.0)
-> 再次报告并等待指示
-> Stage 6 development complete (0.5.0-rc.1)
-> 用户授权 final acceptance
-> frozen C / M1 / M2 / one H1 / package projection / V1 seal (1.0.0)
```

Stage 5 共 32 tasks：先冻结 ESF/TCR 方程、state/RNG/retry、schema/persistence，再做 COAST ESF
oracle、Philox/antithetic Wiener、ESF/IEM/consistency、field chemistry/PSR、TCR、组合、
Checkpoint/diagnostics 和四个接受 task。随机场最小 `N=2`，非定常推荐 `N=4`；每个 accepted
stochastic interval 生成新 Wiener increment，retry 不推进 accepted clock。

Stage 6 共 31 tasks：先冻结 parcel/gas 守恒、state/service/transaction、schema 和双代理燃料
来源，再做 pure-SoA parcel、stencil/trajectory/migration/injector、liquid properties、drag/
heat/mass/evaporation、two-way coupling、ESF common source、static IBM/TAB、组合和接受。
低成本通用性使用 n-dodecane 与 iso-octane surrogate，不宣称真实 kerosene/gasoline。

Stage 5/6 的 COAST 例外尚未开启。任何 source-level oracle 或 `EXEC/Fuels` 读取前，必须让
用户确认 exact current realpath/version；只能把 allowlisted pure module 复制到 Git 外的
generated temp tree，用独立进程和合成输入测试。COAST 源码、case、数据、ABI、消息或控制流
不能进入 Git、产品、安装包或 public tests。

## 7. P0 复用边界

可以按 hash 复用：Cantera/传递依赖来源与许可、Ubuntu 22.04/GCC 11/ABI=1 artifact candidate、
standalone C++ thread/MPI/relocation 证据、公开数学 vectors、双 surrogate 候选 identity。

必须由正式 Stage 4 重建：HUNDUN CMake integration、CPack/RPATH package、public-header isolation、
`ChemistryBackend`、workspace failure mapping、reacting transport、Checkpoint v4、diagnostics v4、
driver、rollback 和科学验收。P0 的 fuel chain 是 `PARTIAL`；机制是否可再分发必须独立决定。

## 8. 测试、资源和停止边界

- 永久禁止 96-cubed；v1 唯一 48-cubed 是最终 H1。
- 开发 task 不等待长测。长 selector 只在软件和 public docs 冻结后运行，或在不消费结果的
  独立任务旁 detached 运行。
- 不启动 Vblowoff、Flame D 或重复 COAST 已完成的科学矩阵。
- M1/M2 可在最终 frozen candidate 上按资源准入并行；H1 必须等待二者通过且独占 H 资源。
- 正常 configure/build/install/runtime/formal tests 不得要求 Python、Conda 或在线 fetch；
  “network-independent” 不得通过断开宿主网络实现。
- 不访问 BOFFIN、私有研究数据，不检查或干扰研究进程；只按 HUNDUN exact path 检查作业。
- 不 push、不发布；不改写已接受历史；不伪造 DCO。

## 9. 当前非阻断记录

- governance 主工作树 `/home/wyf/code_dev/hundun-flow-governance` 仍有用户原有 untracked
  `docs/plans/2026-07-21-hundun-flow-stage7-neural-warm-start.md`；不要触碰。
- intake 保留了两次配置/安装诊断失败，详情在 receipt；最终 clean product scratch build/
  install 已通过并报告 `0.2.0`。
- Stage 3 code tree 的 `VERSION=0.1.0`，accepted product 的 `VERSION=0.2.0`，是 projection
  manifest 明确记录的 override，不要误报为身份不一致修复。
- 当前没有活跃 HUNDUN build/MPI job；历史 Task 11/M2 systemd units 仅为 failed/inactive。

## 10. 首次报告要求

新主 agent 首次回复应报告：当前 HEAD/tree/parent、worktree clean 状态、intake receipt hash、
Stage 3 product blob 不变检查、P0 可复用/不可复用边界、Task `4F-0` 剩余文件白名单、RED 和
focused command。完成 `4F-0` 前不得进入 `4F-1`，Stage 4 接受前不得开始 Stage 5。

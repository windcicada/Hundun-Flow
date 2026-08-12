# HUNDUN-FLOW Stage 4--6 新主 agent 提示词

```text
接手 HUNDUN-FLOW Stage 4--6。开发工作树固定为：
/home/wyf/code_dev/.worktrees/hundun-flow-stage4-reacting-flow

首先完整阅读：
1. AGENTS.md
2. docs/handoff/2026-08-11-hundun-flow-stage4-6-main-agent-handoff.md
3. .superpowers/sdd/stage4-4F-0-baseline-receipt.md
4. AGENTS.md Required Reading 中的 Stage 4--6 规格、reference catalog、Stage 4/5/6 计划、
   integration plan 和 P0 final receipt。

先只读核对当前状态：
- branch=coast/stage4-reacting-flow
- expected handoff commit parent=0f39799692a5da47e9810a4e2262143fb2466996
- expected handoff parent tree=1124e82f63ba61ed1041c9e77218fd08d695f6bb
- current signed handoff commit HEAD/tree 以外部交接结果为准；不得从文档猜自引用值
- Stage 3 code C=0cbd3d5bde4be63bc6346b4b32db771d87c59ea2
- Stage 3 governance G=36bebc292e825fa15272481c6a00c2273fa61ce0
- product P=22ed17b438ffbb121ccda97898580183bd0803f8
- P0 seal=910fb1f7fc3df2e0c596d3682db06db442c03ccf
- integration merge=d45ef02706a17f12d38e050497f076cc5002fb51

融合和 accepted-state intake 已完成。不要重新融合 P0，不要另行融合 planning branch，不要
重复 Stage 3 数值验收。融合审计已证明 Stage 3 product manifest 的 272 个路径及
include/src/tests/CMake/presets/cmake/VERSION 相对 C 不变。P0 只提供候选 provenance、artifact、
standalone boundary、oracle vector 和 intake pattern；它不是 Stage 4 product/science evidence。

你的第一个实现节点是收口 Stage 4 Task 4F-0 的剩余产品测试部分。accepted-state receipt 已有，
从 Stage 4 plan 的 4F-0 Step 2 开始：先写 stage4_source_policy mutation RED，再创建 Stage 4
capability ledger，必要时由主 agent修改 tests/CMakeLists.txt 注册 focused test。文件白名单先
冻结并报告；不得修改数值算法、schema v4 产品 loader、Checkpoint v4、driver 或 Cantera adapter。
只运行 diff check、source-policy fixture 和必要 tests-off configure，不运行数值测试。主 agent
完成 caller-impact、requirements、quality、完整 diff、版权和 DCO 审查后创建签署 task commit。

随后严格执行批准的 Stage 4 27-task 计划：4F contracts -> 4P package -> 4C backend ->
4R reacting -> 4A acceptance。4R-0 的 C-T-C/two-PISO 数学证明必须先于生产 reacting step。
成功 step 始终恰好两次 PISO；不得添加 damping/filter/corrector、放宽阈值或逐 case 调参。
task 内只跑 mutation RED、focused unit/header/policy、必要 small MPI 和最多一个 12^3 smoke；
开发期间不等待长测。

Stage 4 接受并形成 0.3.0 exact-head seal 后停下，向用户报告并建议串行或有限并行；没有用户
新指示不得开始 Stage 5。Stage 5 接受后同样停下再决定 Stage 6。accepted history 固定为
Stage 4 -> Stage 5 -> Stage 6。最终软件与 public docs 冻结后才运行 M1/M2/唯一 H1；永久禁止
96^3，不运行 Vblowoff、Flame D 或重复 COAST 大矩阵。

默认只把边界清晰、证据矩阵冻结的独立任务交给子 agent；中央 CMake/registry、跨模块数学、
完整 diff、来源许可和最终验收由主 agent承担。当前用户规则是优先默认 worker，不使用
luna_worker，直到用户再次改变规则。worker 不提交、不添加 DCO、不扩大范围。

Linux v1 production profile 固定为 Ubuntu 22.04/glibc 2.35+、GCC 11/libstdc++、C++17、
_GLIBCXX_USE_CXX11_ABI=1、generic x86-64。普通 configure/build/install/test/runtime 不依赖
Python/Conda/online fetch，也不要通过断开宿主网络验证。Cantera 3.2 P0 artifact 必须在正式
4P tasks 中逐 hash 升级为 HUNDUN accepted input；public header 不得暴露 Cantera/SUNDIALS 类型。

任何 Stage 5/6 COAST oracle 或 EXEC/Fuels 读取前，必须请用户确认 exact current realpath/version。
COAST 只能进入 Git 外 generated temp oracle，通过独立进程和合成输入使用；不得把源码、case、
数据、ABI、消息或控制流带入 HUNDUN。不要访问 BOFFIN 或研究数据，不检查或干扰研究进程。

governance 主工作树已有用户 untracked Stage 7 文档，不要触碰。product repo 只是投影目标，
不要直接开发。不要 push、发布、改写接受历史或伪造 sign-off。

首次回复先给出只读核验结果、4F-0 剩余白名单、RED/mutation、focused commands、P0 reuse
边界和 blocker；然后连续完成 4F-0，停在该 task 节点报告，不要同一节点自动进入 4F-1。
```

# Stage 3 parallel-completion v2 方案审查记录

**状态：** READY_FOR_USER_REVIEW，尚未激活

**代码基线：** 7fc8c5080528f6ea0dbc787c51ca40d9e0fa4553

**审查范围：** 设计、公开参考、执行计划、AGENTS/ledger 候选映射；不含产品实现

**执行边界：** 未运行数值测试、MPI、sanitizer 或长测，未修改产品/测试源码

## 1. 候选文件

- docs/superpowers/specs/2026-08-09-hundun-flow-stage3-parallel-completion-v2-design.md
- docs/references/2026-08-09-hundun-flow-stage3-public-algorithm-reference.md
- docs/superpowers/plans/2026-08-09-hundun-flow-stage3-parallel-completion-v2.md

| Candidate | SHA-256 |
| --- | --- |
| design | 4cc2052cf99874cefaffe1d11fc954195269fac77ea06fee06a92582ab5ea64e |
| public reference | 0349d4a6301b2213db98ba680557baab808bde62b40f09d99cf392af1e4058b8 |
| implementation plan | bed28804133174481d9ad41266a8d644a0ad056ff9be7ab48f44878b11c45f44 |

三份文件保持 proposed/immutable candidate。只有用户批准后执行 S3-A0，把三份文件的
hash 和批准原文写入 tracked activation receipt，v2 才替代 compact v1 的未完成任务
顺序。A0 之前不得执行 P0、派发 v2 worker 或进入产品实现。

## 2. Pass 1：真实基线与首轮缺口

主 agent 用 rg、现有 CodeGraphF index 和源码只读检查核对：

- FixedStepMaterialDensityFlow::attempt_common 的现有两个 caller及 adaptive-control 影响；
- DensityClosureBridge、DensityClosureAdapter 和 ideal-gas 三阶段 closure；
- FixedStepImmersedFlow 的 53 个直接/传递影响符号；
- diagnostic_report_seal 当前只认证 constant base、没有覆盖 WALE；
- Checkpoint v3 旧 profile-1 overload、publish-last 和 driver caller；
- DiagnosticModuleKind 追加值会影响的 diagnostics/header consumers。

两个默认审查 agent 的首轮报告共暴露以下实质问题，均已反馈到文档：

| 首轮问题 | 修订结果 |
| --- | --- |
| C3 trace 漏 provisional closure | 冻结 predictor/provisional/final closure 和 conditional event trace |
| profiles 5/8 无 body-fitted variable-density WALE seam | 在既有 material PISO attempt_common 增 private optional-WALE hook；禁止第二 PISO |
| 多对象 publish 不原子 | 所有对象先 prepare，一次 collective-ready，之后只允许 noexcept/allocation-free/MPI-free publish |
| material/ideal/WALE report 未认证 | C1/D1/D2 分阶段扩展同一 report seal，逐字段 mutation |
| 单 worktree 伪并行 | P0 registration 分片；C1 后建立独立 infra linked worktree/build tree |
| worker packet 过宽 | 只有 R1/O1 有窄白名单 Worker Packet，worker 不提交 |
| reading exemption 与 AGENTS 冲突 | A0 receipt 激活后才允许 R1/O1 使用 bounded reading exemption |
| Checkpoint API/profile 真值表不完整 | 冻结 const-correct read/write views、overloads、values 1--9 和 required/forbidden matrix |
| 24-cubed exact-counter 证据缺失 | 新增 main-only E1；开发期只跑 8-cubed，V1 执行 24-cubed 1/2/4 |
| V0/V1/V2 身份混淆 | 分离 code candidate C、product commit P、governance report G |
| public references 不可复现 | commit-pinned links、许可证、Basilisk 页面 hash、论文 DOI |

## 3. Pass 2：架构、DAG 与 profile owner

复核结果：

- DAG 为 A0 -> P0 -> C1，随后 scientific main lane 与 R1/O1 infra lane 并行，S1 后
  汇合，再顺序执行 R2/O2/A1/E1/G1/DOC/V0/V1/V2；无环。
- C1/D1/C2/D2/C3 顺序修改共享 flow composition；不并发写 flow_immersed.cpp 或 driver。
- profiles 1--9 都有唯一 runtime/science owner、codec owner、profile-level provider
  integration owner、driver/inventory owner 和 final owner。
- O1 只生产 kind-22 provider component；O2 是全部 profile 的唯一 provider integration
  owner。
- body-fitted material/ideal 复用既有 material PISO，ideal 复用 closure hooks；IBM 路径
  保持 FixedStepImmersedFlow 是唯一 momentum/pressure/final-flux/force owner。
- borrowed registry/closure lifetime 是调用方合同；文档不再声称 C++ 可以运行时证明
  lifetime，并冻结 flow -> closure -> registry 的析构顺序。
- Checkpoint profile 1 bytes 不变，values 2--9 只 additive；WALE transient fields 不持久化。

## 4. Pass 3：可执行 packet、路径与构建

逐 task 检查了 19 个任务的 Ownership、Consumes、Produces、References、Files 和步骤。
结果：

- 所有当前不存在的 Modify 路径都由明确前置 Create 产生；没有 create/modify 顺序
  逆转或重复文件条目。
- P0 是唯一修改 tests/CMakeLists.txt 的后续任务；其 contract 对 build-tree fixture
  删除 include 的 mutation 必须敏感。
- R1/O1 worker packet 都指定 exact worktree、build tree、baseline、allowed/forbidden
  files、资源组、公开参考和 integration owner。
- public header 变更都有 standalone header target；新增 diag_stage3_performance.hpp
  有独立 header contract。
- RED 必须编译并实际启动；C1 增 compile-preserving test seam，C2 先跑 accepted
  characterization 再观察 material unsupported，DOC 在改正文前实际运行 claim RED。
- G1 inventory 固定 build_role、producer_target 和五个 absolute build roots；V0
  构建 Debug/Release/tests-off/ASan/UBSan，所有测试只由 inventory 执行一次。
- Release root 显式构建 TGV、scientific-row aggregator、全部 scientific producer 和
  performance producer。
- R2 冻结 12-cubed、profiles 1--9、1/2/4-rank continuous-vs-restart formal rows，但开发期
  只 list，不执行。

## 5. Pass 4：科学 selector、资源和证据失效

- S1 独占 scientific formal CTest registration；C1--C3 只注册 direct/fast。
- TGV 12/24/48 使用同一产品路径、相同物理终止时间、cell-average restriction 和两个
  独立 Richardson segment；不允许替代 solver、filter、epsilon 假误差或平均斜率掩盖。
- 24-cubed 使用 M lock，48-cubed 使用 H lock；任意时刻只有一个 H job；96-cubed 永久
  不注册。
- 开发期不启动预计超过 10 分钟的作业；越时 screen 正常停止并移交 V1。
- V0 在 exact candidate C 上运行 low-cost/sanitizer/governance 一次；V1 只验证复用并
  新增 scientific/performance，不重复小矩阵。
- 证据失效表覆盖 product/header/build graph、flow ordering、density、WALE、Checkpoint、
  diagnostics、shared harness、counter formula、artifact schema、compiler/MPI 和 cpuset。
- Task 11 pressure/operator/final-flux/force authority 未改时复用其接受；若被影响必须停止
  并由主 agent判定证据失效，不能静默继续。

## 6. Pass 5：版权、治理和投影

- 公开参考只提供数学行为、责任分层和测试思想；不复制源码、注释、控制流、ABI、宏、
  数据布局或错误文本。
- OpenFOAM/Basilisk 的 GPL 代码只作公式与架构对照；所有项目均不成为运行时依赖。
- 私有 BOFFIN/COAST/COAST-2 与研究数据不作为实现来源，未访问。
- product projection 由 manifest/whitelist contract 驱动，tests 和 .superpowers 不进入
  product；product 0.2.0 在 V1 通过后才产生签署 commit P。
- humanizer-zh、shuorenhua 只在 G1 完成后的 DOC task 使用，且主 agent随后重新核对
  公式、命令、key、单位、SHA 和法律文本。

## 7. 机械检查

- git diff --check：PASS
- Markdown fence parity：PASS
- 19 task required-field audit：PASS
- task step sequence audit：PASS
- create/modify/read path-order audit：PASS
- duplicate task-file audit：PASS
- post-P0 root tests/CMakeLists.txt ownership audit：PASS
- unfinished-marker scan：PASS（仅保留有明确定义的模板字段 task-id、row-id）
- product/test source diff：none
- long/MPI/sanitizer execution during this documentation task：none

## 8. 独立复审

- 架构复审最终结论：READY，无 Critical/Important。
- 可执行性复审最终结论：READY，无 Critical/Important/新增 Minor。
- 两次复审均由默认 agent 完成，没有 model 或 reasoning-effort override。

## 9. Verdict

READY_FOR_USER_REVIEW。

这份修订减少开发期矩阵和共享文件冲突，但不降低两次 PISO、单一 authority、守恒、
closure、force、rollback、Restart、MPI 或最终科学门槛。当前唯一外部前置条件是用户对
本 immutable candidate 的明确批准；批准后先执行 A0，不直接开始产品任务。

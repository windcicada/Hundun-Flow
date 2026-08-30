<!-- SPDX-License-Identifier: Apache-2.0 -->

# HUNDUN-FLOW v0.4 Cartesian/低速可压缩 Re=3900 后续工作交接

日期：2026-08-21

状态：`ACTIVE / NOT RELEASE-ACCEPTED`

本文件是下一主 agent 的当前状态单一事实源。详细数学、公开项目对照、验收规则和机器
证据格式仍以本文列出的权威文件为准；不要把本文件当作它们的副本。

## 1. 接管目标与完成定义

继续推进 HUNDUN-FLOW v0.4 Cartesian/低速可压缩性能架构，直至 Re=3900 圆柱的
数值与性能联合门通过。全过程保持以下优先级：

1. 性能优先，但数值精度、守恒、true residual 和事务语义不可降级；
2. 源程序保持算例通用，圆柱和后台阶只由运行目录中的 `case.json + *.d + STL`
   定义；更换射流、后台阶或其他 Cartesian 算例不应修改产品源码；
3. HUNDUN 的对应功能范围通常覆盖 COAST：单相亚音速低速可压缩、局部绝对压力
   EOS、压力修正中的 `drho/dp`、声速/Mach、NSCBC、Restart/Visit/screen；
4. 只复用公开数学、数据布局和生命周期思想，不复制 GPL 实现或 COAST 旧 Fortran；
5. v0.4 只支持均匀或全局张量积拉伸 Cartesian 网格、静止 STL IBM；不加入
   body-fitted 多块网格或非匹配 refinement patches；
6. 最终完成必须有不可变候选、匹配 COAST 的同机同资源性能证据、实验数值证据和
   provenance 审计，不能把“短测能跑通”改写成完成。

Stage 5/燃烧由治理仓库的另一条 lane 独立推进。本 lane 不实现或接管燃烧。到生产冻结
节点时，只记录可兼容接入的稳定 seam：
`FrozenExecutionGraph`、`ContributionPlan`，以及
`C1 -> transport -> PISO1 -> C2 -> PISO2`。

## 2. 权威阅读顺序

下一主 agent 在修改、构建或启动 MPI 前，完整阅读：

1. `AGENTS.md`；
2. 本文件；
3. `docs/superpowers/specs/2026-08-12-hundun-flow-v0.4-cartesian-performance-architecture-design.md`；
4. `docs/superpowers/plans/2026-08-12-hundun-flow-v0.4-cartesian-performance-architecture.md`；
5. `docs/architecture/v0.4-target-hot-loop.md`；
6. `docs/research/2026-08-20-v04-open-source-method-comparison.md`；
7. `docs/research/2026-08-21-v04-cylinder3900-backstep-primary-data.md`；
8. `versions/v0.4/tests/focused_manifest.cmake`；
9. `docs/verification/v0.4-re3900-performance-policy.json`；
10. `docs/verification/v0.4-re3900-coast-equivalence-template.json`；
11. `docs/verification/v0.4-re3900-candidate-ledger.md`；
12. `docs/verification/v0.4-literature-data-receipt.md` 和同名 JSON。

公开项目对照的引文和采用边界集中在第 6 项研究文档与
`docs/references/2026-08-13-hundun-v04-adoption-ledger.tsv`，不要在交接后凭记忆重新
解释 OpenFOAM、AMReX/AMReX-Hydro/IncFlo、HYPRE 或 PETSc 的实现。

## 3. 当前 Git 与工作区事实

当前 v0.4 权威工作目录是：

```text
/home/wyf/code_dev/hundun-flow
```

2026-08-21 完成交接文件后的最终核对：

```text
branch: main
HEAD:   8decc0d7e9f5815c091241cc2645c2a235e435d2
tree:   a024033c2b1b332ff413a3e4c9a037b9da2e361e
tracked modified/deleted entries: 45
untracked entries: 79（包含本交接文档和同目录提示词）
tracked diff stat: 45 files changed, 3568 insertions, 481 deletions
```

接管后的受控 checkpoint 取代上述初始 dirty-tree 身份作为当前推进依据：

```text
governance/docs branch: codex/re3900-v04-candidate
receipt-seal commit:    213c74c75744ace300edbe82b7253af08615e597
receipt-seal tree:      3e242e9f9deb6f5e1ea9ff69e06ce63050773f46
identity-fix commit:    7f0f0cc644b3155d86a8f1737cf714af1512b937
identity-fix tree:      3f3ac22b3bfdfa66f0b569f79102a021b41eca2f
HUNDUN source commit:   f0601781357d8e7c2098a359ca368dd00e7d9cf0
HUNDUN source tree:     8f6bf5d4297a234bc9eb56fec5605a6d9c6071ac
COAST source commit:    3c22e0f029db1b2ca045ec9e212a95eacbcfe6a3
COAST source tree:      ba449790e4918e7a3fd5c21e71c4e9f980a4691f
```

HUNDUN 与 COAST 的候选源码 worktree 均 clean。不可变 pre-full2 snapshot 是
`/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/candidates/hf060178-c3c22e0f-r2`；
其下无可写文件。旧 dirty-tree 事实和旧 snapshot 只保留为 provenance，不能晋级。
首个 `hf060178-c3c22e0f` snapshot 的 build identity 仍指向旧 snapshot 路径和旧 COAST
input/IBM hashes，已 fail-closed；`-r2` 中的 identity 由 clean HUNDUN source 重新生成，
其自身 SHA-256 为
`4ec0099b95aef1ab0ea117a7b2b8ec8e88c5409e38740ab529ba7ac43f8bae05`，
11 个绑定输入路径和哈希已逐项重验。

`.superpowers/sdd/` 下的入口指针受该目录 `.gitignore` 管理，不计入 status。

这些大量 dirty/untracked 内容是当前 v0.4 产品工作，不是垃圾。保持原样，不执行
`reset --hard`、`checkout --`、`clean`、stash、覆盖或自动提交。候选冻结工具要求 clean
tree，但如何形成 checkpoint/提交需要主 agent 完成整体验收后取得明确授权；不能为让
工具变绿而清理现有工作。

另有两个已登记 worktree 和一个可访问的旧 Stage 4 linked tree，含义不同：

- `/home/wyf/code_dev/.worktrees/hundun-flow-cylinder3900-stage4` 是旧 v0.3/Stage 4
  调试树，仍然 dirty；
- `/home/wyf/code_dev/.worktrees/hundun-flow-cylinder3900` 是旧 benchmark 分支；
- `/home/wyf/code_dev/.worktrees/hundun-flow-v04-cartesian-performance` 是较早的 v0.4
  分支快照。

不要把旧 Stage 4 树自动合并、清理或当成当前 v0.4 权威。Stage 5 在治理仓库的独立
worktree 中并行，亦不属于本交接范围。

### 接管后的第一轮只读核对

```bash
cd /home/wyf/code_dev/hundun-flow
git rev-parse --abbrev-ref HEAD
git rev-parse HEAD HEAD^{tree}
git status --short --branch
git diff --no-ext-diff
git ls-files --others --exclude-standard
git diff --check
git worktree list --porcelain
ps -eo pid,ppid,stat,etime,cmd | rg 'hundun|mpiexec|mpirun|orted' || true
find /home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900 \
  -maxdepth 3 -printf '%y %p -> %l\n' | sort
```

完整查看 diff 和重要 untracked 源文件后再判断状态；`git diff` 本身不会显示 untracked
文件内容。不得终止不属于本 lane 的进程。

## 4. 已实现并重新验证的架构

计划状态表把 Tasks 14--19 标为 `LOCALLY_ACCEPTED_UNCOMMITTED`；Task 20 为
`IN_PROGRESS`，Tasks 21--22 为 `PENDING`。各任务正文中的部分 checkbox 尚未同步，
不要仅凭 checkbox 推断实现缺失，应以代码、focused manifest 和验收结果完成逐项审计。

当前产品已具备：

- 扁平的 `versions/v0.4/src` 产品布局和单次 case compile/analyze/allocate/bind/seal；
- uniform/tensor-stretched Cartesian、STL 扫描、静态 IBM 二次重构与 compact donor
  exchange；
- primitive state + 派生 EOS/transport，不保存长期 `rhoU/rhoh`，无常密度快速路径；
- 局部绝对压力 EOS、`drho/dp` pressure storage、声速/Mach、开边界与闭域质量压力
  reference；
- thermophysical predictor 后冻结 `h*/Y*`，每 accepted step 恰好两次 PISO；
- 唯一 final face-mass-flux writer、terminal EOS/continuity/closed-mass/gauge 四门；
- WALE、默认 Vreman + wall function，统一 gradient/`mu_eff` authority；
- immutable linear plan、exact/coarse numeric refresh、setup reuse、固定 workspace；
- Restart 默认覆盖 current generation，Visit/screen/evidence 只读 committed state；
- accepted/trial 事务、rank-local failure consensus、零热路径 heap allocation 契约；
- C1/C2 model contribution seam，纯流 v0.4 绑定空 plan。

### 本轮为全网格正确性修复的关键点

1. Native MG 只对 fluid activity 求解：solid row 为 identity，fluid-solid face link 从
   preconditioner 中移除；activity map 是静态拓扑缓存。
2. coarse residual 对 finite-volume cell-integral residual 做守恒求和，不做 volume
   average。曾尝试的 normal-distance coarse-face scaling 使
   `v04_solver_mg_coarse_use` 退化，已经完整回退，不要恢复。
3. FGMRES 的 in-cycle FP64 true-residual audit 使用未占用 basis workspace 检查 trial
   solution，不再每个 audit interval 隐式重启 Arnoldi；restart 只发生在注册边界。
4. IBM pressure operator 与 final flux 采用同一不可穿透界面定义：精确移除
   fluid-solid regular link，solid row 隔离。旧 quadratic pressure-donor 路径与 final
   flux 不一致，已经移除；不要恢复 `ibm_pressure_donors`。
5. surface force 在 force Allreduce 前先进行全 rank 状态 consensus，避免局部重构失败
   造成 collective 次序分叉。
6. 正值材料属性使用 positive-bounded quadratic reconstruction：二次值在严格正 donor
   包络内原样保留，overshoot 投影到包络，非正/非有限 donor 使 attempt 失败。API 现在
   保证失败不改写调用者输出，并有 constant/overshoot/atomic rejection 回归。
7. momentum assembly 后先做 rank consensus，再进入 PISO collective；这修复了少量
   rank 局部 IBM 失败时与其他 rank 进入不同 collective 的死锁。
8. equation assemblers 已补齐完整 cell/face alias graph、checked storage intervals、
   dependency storage/revision domains、历史层实际 ghost reach、有限 RHS、tile coverage
   与只在 `AssemblyEpoch::finalize` 发布 certificate 的事务合同。

## 5. 当前验证证据

### 5.1 当前源码 focused

在最后一次 positive-property atomicity 修改之后执行：

```bash
cmake --build build/v04-task14-red-clang -j2
ctest --test-dir build/v04-task14-red-clang --output-on-failure -L v04_focused
```

结果：`107/107 passed`，real time `155.57 s`。覆盖 low-Mach MMS、enthalpy、species、
两次 PISO、IBM pressure/force oracle 与 mutation、MG/Krylov、Restart、100-step hot
resource、1/2/4-rank halo/equation/PISO/IBM 等。

当前源码的预注册 ASan 子集也已执行，`LastTest.log` 中 10 项全部 `Test Passed`：

```text
v04_public_headers
v04_solver_low_mach_mms
v04_solver_enthalpy_terms
v04_solver_species_conservation
v04_solver_piso_mutation
v04_solver_ibm_equation_interface
v04_solver_ibm_force_mutation
v04_core_product_freeze
v04_app_driver
v04_core_transaction
```

运行时需要：

```bash
LD_LIBRARY_PATH=/home/wyf/.local/opt/hundun-toolchain/clang/lib/x86_64-unknown-linux-gnu
ASAN_OPTIONS=detect_leaks=0:abort_on_error=1
```

日志：`build/v04-task14-asan-clang/Testing/Temporary/LastTest.log`。单独的
`v04_mesh_ibm_quadratic_test` 也在当前源码下通过 ASan。最后一次 atomicity 修改后尚未
重跑完整预注册 UBSan 子集，这是下一次候选 focused 前必须补的明确事项。

交接时 `git diff --check` 通过，未发现本 lane 的 HUNDUN/MPI 进程。

### 5.2 HUNDUN-only 全网格两步诊断

运行目录：

```text
/home/wyf/code_dev/.benchmarks/hundun-flow-v0.4-re3900/full2_clean_candidate_01
```

证据：`evidence.jsonl`，64 ranks，`480x480x48`，两步成功：

| step | mode | max-rank time | total linear iterations |
|---:|---|---:|---:|
| 1 | BE startup | 34.593771570 s | 247 |
| 2 | BDF2 | 23.918127001 s | 168 |

BDF2 step 的两个 pressure stages 分别为 `11.668151234 s` 和 `11.552785031 s`，即压力
路径占绝大多数 step 时间；两个 pressure solve 各约 84 iterations。step 2 记录
6,183,744 structured messages、38,701,168,128 structured bytes、3,087 blocking
collectives、736,935,765 reduction ns、零 heap allocations、零 setup 与两次 setup
reuse。

当次 identity：

```text
binary SHA-256: 63354dd6803afb310ba5f6172e4b085315202564e75831eeb7383a07b416bfff
case SHA-256:   7523c9a00e04190911b4d9cd106e4ab403c77db0e2a26ca79f55851fc4d9d3d7
.d SHA-256:     0cde4c77153be59455964ce778a01afaf92ecf95281fc90f67c95bd7f3510518
STL SHA-256:    bd264c586543de4ec330f53cb2d3d9dfba550823db69451ac3176c603d248f46
product:        18066577226237870531
CpuPlan:        8087156028938437934
```

这是修改前 dirty-tree 的 HUNDUN-only 诊断，不是可晋级 `full2` receipt。之后的
positive-property atomicity 修改已使 binary identity 变化；且当时没有匹配 COAST
运行。详细非晋级记录见 candidate ledger。

## 6. 尚未通过的硬门

### 6.1 不可变 pre-full2 候选已形成，性能 gate 尚未晋级

Tasks 14--19 的 spec/code completion audit、当前候选的 focused 107/107、预注册
ASan/UBSan 10/10、tests-on/tests-off 等价和 oracle isolation 已完成。HUNDUN
production source 固定在 `f0601781357d8e7c2098a359ca368dd00e7d9cf0`，tests-on 与
tests-off binary byte-identical，SHA-256 为
`dbe1d2bc13cca3337ec234825d0256fdd7e1d2e05faa68ad74e46ead865cacb4`。

这只形成不可变 pre-full2 候选。candidate ledger 的正式 `full2` 和后续 `frozen`
gate 尚未晋级；还缺 sealed receipt 下的 matched alternating full2。任何正式运行的
binary/input/case/CPU/MPI/affinity 身份偏移都会拒绝该样本。

### 6.2 COAST 等价基线已封印

`docs/verification/v0.4-re3900-coast-equivalence-template.json` 已在 immutable candidate
上 reseal，`--require-sealed` canonical digest 为
`23bf0beedd5dca508ebae79862376ab486bd44b95f9e4a495e8c262234d28750`。
scientific-work unit 保持“同一完整 accepted physical step”，不要求 SIMPLE/PISO
原生迭代数相等；HUNDUN 固定 two-PISO，COAST 固定 two-SIMPLE，所有原生工作均计时。

COAST 候选 binary SHA-256 为
`c5e483a0442b430e58040a8e04f5adf42b70b3595319c643e9265c87341a9bc4`。
压力 terminal/reliable residual 使用 ICCG scaling 前的 exact seven-coefficient physical
operator 和 high-plus-low pressure solution；BiCGStab/ICCG reliable replacement cadence
分别为 8/128。terminal audit 是唯一 accepted step/time commit authority，output-off
覆盖 per-step schedule 与 finalizer。

fresh candidate validation `diagnostic_coast_frozen_3c22e0f_output_off_02`
恰好执行 BE/BDF2 两个 accepted steps、每步 2 SIMPLE/10 solves；所有 true residual
和 terminal gates 通过，无 cap hit，Visit/Restart/Log/Results/monitor 文件数为零。
evidence SHA-256 为
`a6bbc871ae5f8c62792bbdab45081d3ed6e296f2bb560ae86e132055c8b6944e`。
这不是 matched full2；现在才允许启动正式 pairing。

### 6.3 性能门未通过

正式策略冻结在 `v0.4-re3900-performance-policy.json`：64 ranks、`480x480x48`、
full2 只作方向筛查；full20 用 steps 1--5 warmup、6--20 measured，至少五组交替 HC/CH
paired runs，接受要求 median ratio 的 deterministic bootstrap 95% upper bound `<=1.0`。

当前只有 HUNDUN 诊断，不能给出 ratio。压力路径是明确热点。下一主 agent 可以先做
架构级压力性能审查，但以下边界已冻结：

- Native MG 当前是 flexible preconditioner；其 residual-dependent 行为未被证明为
  symmetric/Galerkin，因此不可仅为进入 PCG 把 certificate 改成 `fixed_spd`；
- 若选择 PCG，必须先独立建立 symmetric/Galerkin MG 数学与 1/2/4-rank 证据；
- true residual、终态 continuity 和两次 corrector authority 不能为减少 collectives
  而放宽；
- 性能候选冻结前只做 focused 和完整网格短测，不启动长统计。

### 6.4 文献数值门 fail-closed

`v0.4-literature-data-receipt.json` 的 `complete=false` 是预期结果。已冻结的一手值包括
Parnaudeau `St=0.208±0.002`、`Lr/D=1.51`，Norberg 的 Re=3900 两侧直接 St 点和明确
标注的公式派生值，以及 NASA 后台阶 `xr/H=6.26±0.10` 与原始 ASCII 数据。

仍缺 Parnaudeau Fig. 11--15 三站原数组/受控数字化误差，以及 Re=3900 直接总阻力和
有限展向 `Cl_rms` authority。2026 INRAE 数据集可作为一手补充，但必须先核对实验批次、
坐标原点、归一化和 station mapping。以下命令目前必须失败：

```bash
python3 tools/v04_literature_extract.py receipt-validate \
  --require-complete docs/verification/v0.4-literature-data-receipt.json
```

不要用后续 CFD 论文的常见汇总值、LES 采样误差或查看 HUNDUN 输出后选择 tolerance 来
补齐缺口。literature receipt 完整前，不启动约 `2020 D/U` 的长统计。

## 7. luna_worker 委派合同

用户要求使用 `luna_worker` 完成边界清晰、证据矩阵已冻结的独立任务；所有需要完整
上下文的规划和审查由主 agent 完成。

主 agent 必须亲自承担：

- 总体计划/规范修改和跨模块依赖排序；
- PISO、压力算子、MG/Krylov、COAST scientific-work equivalence 的数学决策；
- 是否冻结候选、是否接受性能/文献门，以及所有 worker 输出的最终 review；
- 与 Stage 5 seam 的兼容性判断和 exact HEAD/tree 生产冻结记录。

只有当主 agent 先写出以下冻结 brief 时才调用 `luna_worker`：

1. 单一、可独立完成的目标；
2. 明确文件所有权或只读路径；
3. 明确禁止触碰的文件/lane；
4. 已冻结输入、oracle、accept/reject 阈值和精确命令；
5. 要返回的 diff、日志、hash 和未解决问题。

适合 `luna_worker` 的任务示例：

- 按固定矩阵执行 UBSan/1/2/4-rank 测试并整理只读证据；
- 对已冻结 COAST equivalence 字段做输入/能力 inventory，不替主 agent决定等价规则；
- 按固定 JSON schema 获取、校验和提取一手文献数据；
- 在唯一文件集内添加一个已给定 oracle 的回归测试；
- 在主 agent 已冻结公式、接口和验收命令后实现一个孤立组件。

每个 worker brief 都要说明：它不是代码库中的唯一执行者，要保留他人修改、适应当前
tree、不回退共享内容、不提交。主 agent 接收后必须查看完整 diff、核对越界修改、重跑
权威测试，不能把 worker 的“通过”直接当成发布证据。

以下任务不委派：改变计划、选择压力性能路线、修改 COAST 对应功能范围、修改 gate
policy、封印 candidate、解释 paired statistics、写最终 ACCEPT/REJECT。

## 8. 推荐恢复顺序

截至 identity-fix commit `7f0f0cc644b3155d86a8f1737cf714af1512b937`，原顺序中的
第 1--4 项和第 6 项的 pre-full2 identity 工作已经完成。第 5 项已有冻结 HUNDUN
候选，是否继续压力优化必须由正式 matched full2 的方向证据决定。当前唯一下一步是：
按 sealed receipt 和 frozen CPU plan 启动 fresh alternating matched full2；完成前不
修改 candidate identity。文献 receipt 仍不 complete，因此仍禁止 full20。

1. 完成第 2、3 节的只读核对，确认交接后没有新进程、新 diff 或外部证据变化。
2. 主 agent 对 Tasks 14--19 做一次 spec/code completion audit，重点查看最近的 MG、FGMRES、
   IBM pressure/force、positive-property 和 collective consensus 修改；不要仅重复测试。
3. 将当前预注册 UBSan 子集作为第一个冻结 brief 交给 `luna_worker` 执行，主 agent
   审查日志；随后复核 tests-off isolation/build identity。
4. 主 agent 解决 COAST SIMPLE 与 HUNDUN two-PISO 的 scientific-work equivalence 规则，
   再把只读 COAST capability/input inventory 作为独立 Luna 任务。
5. 在不改变科学工作量的前提下，围绕两个 pressure stages 的 iteration、collective、halo
   与 MG transfer/operator consistency 制定性能候选。每个候选按 focused -> 当前全网格
   2-step 的顺序验证；失败候选完整回退，不进入 ledger。
6. 代码和配置稳定后取得 clean candidate 授权，生成 tests-on/tests-off/inputs/CPU/MPI
   immutable identity，重新从 `focused` gate 开始。
7. 封印 COAST equivalence receipt 后执行交替 matched full2；方向门通过才冻结并进入至少
   五组 full20 pairing。
8. 可并行推进缺失的一手文献数据，但只有 complete receipt 才允许 HUNDUN-only 长统计。
9. Task 21 最终 ACCEPT 后，Task 22 才用纯运行输入建立 Driver--Seegmiller 后台阶；若需
   产品源码修改，回到对应早期任务实现通用能力并重跑完整 Re=3900 门。

## 9. 交接完成检查

本交接只完成“换主 agent”，没有完成活跃工程目标。下一 agent 接管后应持续推进，直到
以下证据同时存在：

- 当前不可变候选的 focused/ASan/UBSan 全绿；
- matched COAST/HUNDUN full2 与至少五组 full20 paired performance 接受；
- terminal EOS/continuity/closed-mass/gauge、力 oracle/mutation 和所有 provenance 门接受；
- complete literature authority 下的 Re=3900 物理统计接受；
- `docs/verification/v0.4-re3900-final-acceptance.md` 给出唯一 `ACCEPT`。

在此之前保持目标 active，不把 pending/NEAR/diagnostic 解释为发布通过。

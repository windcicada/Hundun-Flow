# HUNDUN-FLOW Agent Handoff

**Handoff date:** 2026-07-16
**Project directory:** `/home/wyf/code_dev/hundun-flow`
**Approved method:** subagent-driven development
**Approved implementation scope:** Stage 0 and Stage 1 only

## 1. Current State

HUNDUN-FLOW is a new, copyright-clean C++ combustion CFD project. The project
directory currently contains only this handoff, the approved design, the
approved Stage 0/1 plan, and project-level agent instructions.

No HUNDUN-FLOW source code has been written. No Git repository has been
initialized here. No worker has been dispatched. This state is intentional:
the dedicated HUNDUN agent must begin with Stage 0 of the approved plan.

The existing COAST code and ongoing research remain independent and untouched.
They may be used only for black-box scientific requirements and later
scientific comparisons, never as an implementation template.

## 2. Authoritative Documents

Read in this order:

1. [Project agent instructions](../../AGENTS.md)
2. [Approved clean C++ solver design](../design/2026-07-16-hundun-flow-clean-cpp-solver-design.md)
3. [Approved Stage 0/1 implementation plan](../plans/2026-07-16-hundun-flow-stage0-stage1.md)

Origin records in the private COAST workspace:

| Document | Origin commit | Handoff SHA-256 |
|---|---|---|
| Design | `50af30ae` | `710fe29fdcef3414743c307edee4336aa66b35c6f70f38146a25f54a4d078b1a` |
| Stage 0/1 plan | `33f486e4`, handoff adjustment `6d04fa4` | `cc0eac5dfbbdaf00fd4d6c329e1ae5d122dcec3b00e2e3edc9094be2703bc836` |

At execution start, compute both hashes again and record them in the first
coordinator report. The plan copy in this directory is authoritative because
it accounts for the documentation-only skeleton already existing.

## 3. Approved Product Direction

- Formal name: **HUNDUN-FLOW**
- Chinese name: **浑敦流体模拟框架**
- Repository name: `hundun-flow`
- Executable: `hundun`
- C++ namespace: `hundun`
- Subtitle: *An extensible C++ framework for immersed reacting-flow
  simulation.*
- Language: C++17
- Runtime foundation: MPI-3 plus a lightweight project-owned runtime
- Public license: Apache-2.0 with DCO
- Public runtime: no Python dependency

The name uses the ancient image of Hundun/Dijiang as an engineering
reinterpretation: a pouch-like bounded chamber containing intense combustion,
turbulent transport, and interacting physical processes. It is not presented
as a literal translation of the classical source.

## 4. Copyright Boundary

The final legal comparison baseline is fixed:

```text
/home/wyf/Compress_boffin/AECSC_WDQ/SRC-Spray
```

Rules:

1. Treat the baseline as read-only private audit evidence.
2. Do not use it as a source of algorithms expressed in code, naming, control
   flow, data layout, comments, diagnostics, or interfaces.
3. Treat COAST and COAST-2 as black-box scientific references only.
4. Stage 1 must reuse zero source files from BOFFIN, COAST, and COAST-2.
5. The private audit directory must remain outside the public repository:
   `/home/wyf/code_dev/hundun-flow-private-audit`.
6. Do not create a public `REUSE_MANIFEST.yaml`.
7. Preserve third-party licenses under `LICENSES/` and document them in
   `THIRD_PARTY.md`.
8. Do not claim legal certainty from an automated similarity scan. Record the
   scan as engineering evidence and retain a human review conclusion.

## 5. Stage 0/1 Deliverable

Stage 0 establishes:

- isolated Git history;
- Apache-2.0, DCO, exact NOTICE, third-party policy;
- fixed private baseline manifest;
- source provenance guard;
- private similarity-audit evidence.

Stage 1 establishes:

- strict typed JSON configuration;
- MPI lifecycle and collective failure propagation;
- three-dimensional structured decomposition;
- typed FieldRegistry, FieldStorage, and FieldView;
- arbitrary-width 26-neighbor halo exchange;
- uniform structured mesh;
- stable C plugin ABI;
- versioned per-rank Restart and primitive VTK output;
- conservative second-order MUSCL/SSPRK2 passive scalar;
- multi-rank conservation, convergence, decomposition invariance, and Restart
  continuity tests.

Stage 1 does not claim reacting flow, LES, IBM, TPDF-TCR, chemistry, particles,
spray, WENO, or DG.

## 6. Coordinator Protocol

The dedicated main agent must act as coordinator:

1. Read all three authoritative documents and inspect the documentation-only
   target directory.
2. Create an execution checklist from the 12 plan tasks.
3. Dispatch one fresh implementation worker for Task 1.
4. Review the worker diff and rerun its exact tests.
5. Dispatch a requirements reviewer focused only on plan compliance.
6. Dispatch a code-quality reviewer focused only on maintainability,
   regressions, and hidden dependencies.
7. Resolve findings, accept or reject the task, then close every worker.
8. Repeat sequentially for the next task.
9. Never keep more than five workers open; only one implementation worker may
   be active.
10. After Task 12, run every coordinator review gate and the private
    similarity audit before making a completion claim.

The coordinator must inspect repository state before and after every worker.
Worker reports are not acceptance evidence.

## 7. Required Worker Envelope

Every implementation/review worker prompt must state:

- do not contact the user or request user approval;
- do not touch COAST, COAST-2, research cases, jobs, or data;
- do not copy or translate legacy source;
- do not add deferred physics or unrelated refactors;
- do not publish or push;
- preserve unrelated changes;
- follow the exact task interfaces and TDD order;
- return one status: `DONE`, `DONE_WITH_CONCERNS`, `NEEDS_CONTEXT`, or
  `BLOCKED`;
- report changed files, commands run, exit codes, and unresolved concerns.

Questions are answered by the coordinator from the approved documents. Only a
dangerous or genuinely insufficient plan should be escalated to the user.

## 8. Execution Start Gate

Before dispatching Task 1:

```bash
cd /home/wyf/code_dev/hundun-flow
test ! -e .git
find . -mindepth 1 -type f -printf '%P\n' | sort
sha256sum \
  docs/design/2026-07-16-hundun-flow-clean-cpp-solver-design.md \
  docs/plans/2026-07-16-hundun-flow-stage0-stage1.md \
  docs/handoff/2026-07-16-hundun-flow-agent-handoff.md \
  AGENTS.md
ps -eo pid,ppid,stat,etime,pcpu,pmem,cmd |
  rg 'mpirun|(^|/)coast( |$)|aecsc|coastctl' || true
```

Expected files are exactly:

```text
AGENTS.md
docs/design/2026-07-16-hundun-flow-clean-cpp-solver-design.md
docs/handoff/2026-07-16-hundun-flow-agent-handoff.md
docs/plans/2026-07-16-hundun-flow-stage0-stage1.md
```

Do not stop or interact with any process found by the process check. The check
is informational and protects concurrent research.

## 9. Stop Conditions

Stop the current task and report to the coordinator when:

- the target already contains unapproved files or a Git history;
- the fixed baseline is missing or mutable;
- a task would require copying or translating legacy code;
- a dependency introduces a Python runtime or online build-time fetch;
- MPI/compiler constraints cannot satisfy the approved minimums;
- a worker proposes pressure coupling, IBM, chemistry, TPDF-TCR, spray, WENO,
  DG, moving walls, or multi-part geometry during Stage 1;
- tests cannot be made deterministic without changing an approved contract.

Do not delete or rewrite evidence to make a gate pass.

## 10. Final Acceptance Report

The coordinator report must include:

- accepted commit and `stage1-runtime` tag;
- all task and reviewer workers closed;
- exact build/test commands and pass counts;
- passive-scalar observed convergence orders;
- mass and decomposition error maxima;
- Restart continuation error;
- compiler/MPI versions;
- `ldd` result proving no Python or missing libraries;
- provenance-guard result;
- private similarity-audit result and human conclusion;
- confirmation that COAST and research data were untouched;
- explicit list of deferred stages.

Do not begin Stage 2 from this handoff. Stage 2 requires a new approved plan.

## 11. Bootstrap Prompt

The prompt below may be given directly to the dedicated HUNDUN agent:

```text
接手 /home/wyf/code_dev/hundun-flow 项目。当前目录只是文档交接骨架，不要把
COAST 当作源码上游。

首先依次阅读：
1. AGENTS.md
2. docs/handoff/2026-07-16-hundun-flow-agent-handoff.md
3. docs/design/2026-07-16-hundun-flow-clean-cpp-solver-design.md
4. docs/plans/2026-07-16-hundun-flow-stage0-stage1.md

执行方式固定为 subagent-driven development。主 agent 只负责拆任务、顺序
派 worker、审查、复验和最终验收；每个任务使用 fresh worker，最多同时保留
5 个 worker，且只能有 1 个 implementation worker 活跃。worker 不得联系我或
找我审批，完成或不用后立即关闭。每个实现任务后依次做 requirements review
和 code-quality review。

严格执行已批准的 Stage 0/1 计划，不进入 Stage 2。固定法律比较基线是
/home/wyf/Compress_boffin/AECSC_WDQ/SRC-Spray，只用于私有独立性审计。
Stage 1 不得复用、翻译、机械改写或模仿 BOFFIN、COAST、COAST-2 的源码、
控制流、ABI、数组布局、注释、消息、.d 输入、Decomp、Restart 或兼容层。
COAST/COAST-2 只能作为黑盒科学需求参考。

不要修改、清理、停止、打包或提交 /home/wyf/code_dev/Coast_software 及其
研究算例和运行数据；不要发布或 push。公开构建和运行不得依赖 Python。
本阶段只实现独立 C++17/MPI-3 运行时和被动标量验证，不实现变密度流、
LES、IBM、化学、TPDF-TCR、喷雾、WENO、DG、多部件或运动壁面。

开始前先执行 handoff 第 8 节的 start gate，核对目录只含四份批准文档，
记录 SHA-256 和现有研究进程，但不要干预进程。随后从计划 Task 1 开始。
每完成一个任务，向我报告 accepted commit、测试证据、review findings 和
已关闭 worker；不要只转述 worker 的 DONE。全部 Stage 0/1 gate 通过后停止，
提交最终验收报告，等待下一阶段指令。
```

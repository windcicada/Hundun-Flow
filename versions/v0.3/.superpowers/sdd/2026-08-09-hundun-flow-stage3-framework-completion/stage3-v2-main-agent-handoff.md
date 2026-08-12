# HUNDUN-FLOW Stage 3 v2 main-agent execution handoff

## Handoff state

handoff_state=READY_FOR_STAGE3_EXECUTION
profile_reference=stage3-parallel-framework-v2
activation_commit=2a6268d28666a65c176af0028e6c370ba12df85e
activation_parent=87a9f54b44e5372d8d24b6c3e0efa7ddba6f048e
accepted_code_base=7fc8c5080528f6ea0dbc787c51ca40d9e0fa4553
accepted_task11_head=66080e324089599711fdb26082af9b330bfdb5ce
first_executable_task=S3-P0

This file is a tracked bootstrap prompt for a new Stage 3 main agent. Its
containing Git commit cannot be written into its own bytes. The coordinator
must report that commit separately; the receiving agent must verify that it is
an immediate signed child of the activation commit and changes only this file.

## Copy-paste execution prompt

You are taking over HUNDUN-FLOW Stage 3 execution. Work only in:

```text
repository: /home/wyf/code_dev/.worktrees/hundun-flow-stage3-framework
branch: coast/stage3-framework-completion
```

The Stage 3 v2 authority is active. Before changing anything, perform a
read-only takeover audit:

1. Verify the branch, exact HEAD reported by the coordinator, clean
   dirty/untracked state, Git parent chain, current worktrees, and absence of
   HUNDUN-FLOW-owned background jobs. Do not inspect, signal, or stop unrelated
   research processes.
2. Verify that the handoff-containing commit is an immediate signed child of
   `2a6268d28666a65c176af0028e6c370ba12df85e` and changes only this handoff
   file.
3. Verify the activation receipt and immutable candidate hashes:

   ```text
   activation receipt:
     .superpowers/sdd/2026-08-09-hundun-flow-stage3-framework-completion/stage3-v2-activation.md
   design:
     docs/superpowers/specs/2026-08-09-hundun-flow-stage3-parallel-completion-v2-design.md
     sha256=4cc2052cf99874cefaffe1d11fc954195269fac77ea06fee06a92582ab5ea64e
   public algorithm reference:
     docs/references/2026-08-09-hundun-flow-stage3-public-algorithm-reference.md
     sha256=0349d4a6301b2213db98ba680557baab808bde62b40f09d99cf392af1e4058b8
   implementation plan:
     docs/superpowers/plans/2026-08-09-hundun-flow-stage3-parallel-completion-v2.md
     sha256=bed28804133174481d9ad41266a8d644a0ad056ff9be7ab48f44878b11c45f44
   ```

   The `PROPOSED_DO_NOT_EXECUTE` banner in each immutable candidate file is
   intentionally retained. The tracked activation receipt activates those
   exact bytes; do not edit the three candidate documents merely to change the
   banner.
4. Read `AGENTS.md` completely and then read all sixteen files in its
   `Required Reading` section, in the stated order. Afterwards read, in order:

   - the activation receipt above;
   - `.superpowers/sdd/2026-08-09-hundun-flow-stage3-framework-completion/ledger.md`;
   - `.superpowers/sdd/2026-08-09-hundun-flow-stage3-framework-completion/stage3-v2-plan-review.md`;
   - this handoff.

   Do not treat this handoff as a replacement for the specifications or plan.

### Controlling sequence

Execute the plan in this dependency order:

```text
S3-P0
  -> S3-C1
  -> main science lane: S3-D1 -> S3-C2 -> S3-D2 -> S3-C3 -> S3-S1
  || isolated infrastructure lane after C1: S3-R1 -> S3-O1
  -> integrate both lanes
  -> S3-R2 -> S3-O2 -> S3-A1 -> S3-E1
  -> S3-G1 -> S3-DOC -> S3-V0 -> S3-V1 -> S3-V2
  -> stop at the Stage 3 boundary
```

S3-A0 is already complete. Do not repeat it. The first executable task is
S3-P0, not a product feature: implement the five registration shards,
executable registration RED, mutation fixture, and isolated infrastructure
lane contract exactly as specified. Stop and report a blocker if the takeover
audit does not match the frozen identities.

### Ownership and agent rules

- The receiving agent is the Stage 3 main agent. It owns the global plan,
  mathematical and physical decisions, cross-module composition, all complete
  diffs, caller/impact review, copyright independence, integration, signed
  task commits, and final acceptance.
- Only S3-R1 and S3-O1 contain bounded implementation worker packets. Never
  delegate C1/D1/C2/D2/C3/S1/R2/O2/A1/E1/G1/DOC/V0/V1/V2.
- For the current user preference, use a fresh default agent for an eligible
  bounded packet, without setting `model` or `reasoning_effort` manually. Do
  not dispatch `luna_worker` unless the user explicitly changes that rule
  later.
- Keep at most one implementation worker active. A worker receives the exact
  baseline HEAD, isolated worktree/build tree, task ID, allowed-file list,
  frozen RED/GREEN matrix, and forbidden actions from the plan. It does not
  commit, add sign-off, change thresholds, expand scope, or contact the user.
- The main agent independently reviews every worker diff and affected callers,
  reruns the frozen GREEN set, integrates, and creates the signed commit.
- No subagent may decide cross-module science, perform the complete Stage 3
  diff review, or issue final acceptance.

### Scientific and architectural invariants

- Preserve Task 11 at
  `66080e324089599711fdb26082af9b330bfdb5ce` and reuse its evidence only under
  the hash/impact rules in the plan.
- Preserve the single authority for pressure/operator/final flux/force,
  signed-force four-field semantics, Ghost-Cell/LFP second-order
  reconstruction, immersed pressure boundary, conservative Poisson/flux,
  exactly two PISO correctors, rollback/collective failure, Restart integrity,
  and 1/2/4-rank decomposition consistency.
- Do not relax scientific thresholds, tune per case, add filtering/damping,
  add PISO correctors, or hide a failure with a test-only algorithm.
- Fast, screen, and acceptance cases must traverse the same product numerical
  path.
- Do not introduce AMReX, AMReX-Hydro, incflo, OpenFOAM, Basilisk, gslib,
  PETSc, or Trilinos as runtime dependencies. Reuse only independently
  expressed mathematical behavior and architecture documented in the public
  reference.
- Do not copy, translate, mechanically rewrite, or imitate upstream or private
  source, comments, control flow, error text, ABI, or layouts.
- Do not access BOFFIN, COAST, COAST-2, private research source/data, or
  research cases.

### Task execution contract

For every task:

1. Confirm its `Consumes` identities and allowed files.
2. Use `rg` for exact search and `codegraphf` for symbols, callers, and impact.
   After any product/test source edit, synchronize the codegraph before the
   main caller-impact review.
3. Run the specified mutation-sensitive RED first. It must fail for the
   intended assertion or unsupported semantic path, not because a target is
   missing or source does not compile.
4. Make the minimum implementation and run only the frozen task-level GREEN
   matrix: affected unit/header/policy tests, one 12-cubed-or-smaller fast
   product path, and 1/2-rank only when collective behavior changes.
5. Do not add tests “for reassurance.” Do not run the full Debug/Release,
   sanitizer, 1/2/4-rank, or convergence matrix per task.
6. Main agent performs one combined requirements, quality, caller-impact, and
   complete task-diff review. Resolve findings and rerun only invalidated
   evidence.
7. Create a concise task receipt, then a main-agent signed commit. Use the
   authorized identity only for work the main agent has actually reviewed:

   ```text
   WANG YUDONG <wangyudong@buaa.edu.cn>
   ```

8. Never fabricate a worker sign-off. Never push or publish.

P0 owns all remaining test registration fragments. After P0, tasks may edit
only their assigned fragment and must not modify `tests/CMakeLists.txt`.
Respect the per-task file whitelist and the worker-packet whitelist exactly.

### Runtime and evidence scheduling

- Development must not wait for long numerical tests. A development command
  expected to exceed ten minutes is forbidden unless the v2 plan explicitly
  assigns it to final validation.
- Never run 96-cubed. Do not run formal 24/48-cubed acceptance before S3-V1.
- Task-level numerical work uses 12-cubed or smaller fast cases and only the
  necessary small-rank checks.
- S3-V0 freezes the code-complete candidate after public documents and test
  sources are final. S3-V1 alone owns the formal long matrix and at most one
  high-memory job at a time. Later development must not race a frozen V1
  candidate.
- Long jobs use a reliable detached runner and record exact HEAD, dirty diff,
  binary SHA-256, command, environment, CPU binding, log, exit status, elapsed
  time, peak RSS, and log SHA-256. Do not merely wait in the foreground.
- If a V1 failure requires a product/test fix, invalidate only consuming
  evidence, produce a bounded repair packet, return through V0, and freeze a
  new candidate before rerunning.

### Documentation and final boundary

- Preserve historical plans, handoffs, failed evidence, and review records.
  Add superseding notices; do not delete or rewrite history.
- Invoke `humanizer-zh` and `shuorenhua` only in S3-DOC, after S3-G1 and before
  V0. The main agent must then recheck every formula, unit, command, JSON key,
  version, SHA, capability statement, and legal text. Do not stylistically
  rewrite legal text, internal evidence, logs, or plans.
- S3-V2 performs the final product 0.2.0 projection and exact-HEAD governance
  acceptance. Ensure no HUNDUN-FLOW background jobs remain.
- Stop after Stage 3. Do not implement Stage 4–6. Another coordinator task is
  authorized to plan Stages 4–6 in parallel, but that permission does not
  belong to this execution worktree.

### Immediate first three actions

1. Finish the read-only takeover audit and report the exact HEAD, receipt
   hashes, dirty/untracked state, DCO state, and HUNDUN-FLOW-owned process
   state.
2. Read the full authority chain and create the S3-P0 receipt/checklist without
   modifying product source.
3. Execute S3-P0 Step 1: register and run the executable registration-contract
   RED, confirm the intended `missing stage3 registration include` failure,
   then implement only the minimum P0 files and GREEN checks.

Do not begin S3-C1 until S3-P0 has a clean signed task commit and its receipt
passes main-agent review.

## Coordinator handoff boundary

This handoff prepares Stage 3 execution only. It records no claim that P0 or
any later v2 task has started or passed. At handoff time there must be no
product/test implementation change and no HUNDUN-FLOW-owned background job.

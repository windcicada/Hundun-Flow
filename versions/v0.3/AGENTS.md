# HUNDUN-FLOW Agent Instructions

## Current State

Stage 1 is accepted at the annotated `stage1-runtime` tag, Stage 2 is accepted
at its recorded coordinator-ledger head, and Stage 3 is accepted by the seal at
`36bebc292e825fa15272481c6a00c2273fa61ce0`. Do not reinitialize the repository
or rewrite accepted history. The Stage 4--6 Linux CPU v1 architecture and
detailed plans are the approved implementation authority. The independent P0
preflight is sealed at `910fb1f7fc3df2e0c596d3682db06db442c03ccf` with
`PREFLIGHT_PASS`; its evidence remains candidate input until formal Stage 4
tasks revalidate it. Stage 4 is accepted at tested code head
`6407cd7c591ce088db7f1dd7e296d77acd18da1c`, version `0.3.0`, by the governance
commit containing `.superpowers/sdd/stage4-final-acceptance-report.md`. Stage 5
has not started. Stop at the Stage 5 user decision gate until the user gives
explicit authorization. Determine progress from Git history and tracked
receipts, not from plan dates or an unaccepted worktree.

## Required Reading

Read these files in order before changing anything:

1. `docs/handoff/2026-07-16-hundun-flow-agent-handoff.md`
2. `docs/design/2026-07-16-hundun-flow-clean-cpp-solver-design.md`
3. `docs/development/naming-and-style.md`
4. `docs/development/source-policy.md`
5. `docs/plans/2026-07-16-hundun-flow-stage0-stage1.md`
6. `docs/superpowers/specs/2026-07-19-hundun-flow-stage2-variable-density-flow-design.md`
7. `docs/superpowers/plans/2026-07-19-hundun-flow-stage2-variable-density-flow.md`
8. `docs/superpowers/specs/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale-design.md`
9. `docs/superpowers/plans/2026-07-27-hundun-flow-stage3-static-lfp-gcibm-wale.md`
10. `docs/superpowers/specs/2026-08-05-hundun-flow-task11-force-consistency-authority-addendum.md`
11. `docs/superpowers/plans/2026-08-05-hundun-flow-stage3-science-first-execution-amendment.md`
12. `docs/superpowers/plans/2026-08-05-hundun-flow-task11-science-closure-red-design.md`
13. `docs/superpowers/specs/2026-08-08-hundun-flow-post-task11-semantic-port-architecture-design.md`
14. `docs/superpowers/specs/2026-08-09-hundun-flow-stage3-compact-scientific-design.md`
15. `docs/superpowers/plans/2026-08-09-hundun-flow-stage3-framework-completion.md`
16. `docs/superpowers/specs/2026-08-09-hundun-flow-stage3-parallel-completion-v2-design.md`
17. `docs/references/2026-08-09-hundun-flow-stage3-public-algorithm-reference.md`
18. `docs/superpowers/plans/2026-08-09-hundun-flow-stage3-parallel-completion-v2.md`
19. `docs/superpowers/specs/2026-08-09-hundun-flow-stage4-6-linux-cpu-v1-architecture-design.md`
20. `docs/references/2026-08-09-hundun-flow-stage4-6-reference-catalog.md`
21. `docs/superpowers/plans/2026-08-09-hundun-flow-stage4-reacting-flow.md`
22. `docs/superpowers/plans/2026-08-09-hundun-flow-stage5-esf-tpdf-tcr.md`
23. `docs/superpowers/plans/2026-08-09-hundun-flow-stage6-spray.md`
24. `docs/superpowers/plans/2026-08-09-hundun-flow-stage4-6-v1-integration.md`
25. `docs/superpowers/specs/2026-08-09-hundun-flow-pre-stage4-p0-preflight-design.md`
26. `docs/superpowers/plans/2026-08-09-hundun-flow-pre-stage4-p0-preflight.md`
27. `.superpowers/stage4-p0/final-receipt.md`
28. `.superpowers/sdd/stage4-4F-0-baseline-receipt.md` when present

For an active unaccepted task, also read its current coordinator handoff when
present in `.superpowers/`. A handoff is takeover evidence, not a portable
replacement for the tracked specification, addendum or plan.

The overall design controls architecture and the complete capability roadmap.
The Stage 3 specification controls the current scientific and numerical
scope. The 2026-08-05 force-consistency addendum restores the Task 11 hard
gates and requires an executable sign/orientation RED before changing their
formula. The Stage 3 plan controls its exact additive interfaces, Tasks 1--21,
tests, commits, and acceptance except where the 2026-08-05 science-first
execution amendment explicitly changes ordering, bounded-worker eligibility,
performance-debt timing, Task 20 product allowed-file subclusters, testing
cost, resource scheduling, final-evidence ownership or governance timing. The
2026-08-09 compact-scientific design and framework-completion plan supersede
the post-Task-11 task order, per-task test cadence, repository topology and
product projection timing while preserving every scientific, MPI,
transaction, Restart and user-executability requirement not explicitly
deferred there. The
parallel-completion v2 design and plan were activated by the tracked
`stage3-v2-activation.md` receipt and are the controlling authority for work
after accepted Task 13+19B. They supersede only the remaining-task order,
worker packets, development test cadence and final long-test scheduling. They
never reinterpret an accepted task or the Task 11 scientific authority. The
tracked Task 11 RED design and any current coordinator handoff control takeover
evidence until Task 11 is accepted. The Stage 2 documents
continue to control every frozen Stage 2 contract; the Stage 0/1 plan and
handoff continue to control frozen Stage 1 contracts and legal boundaries.
After Stage 3 acceptance and explicit user authorization, the 2026-08-09
Stage 4--6 specification controls v1 science, state, service, packaging,
copyright and capability boundaries. The three stage plans control their
task-local files and tests; the integration plan controls stage entry,
central-file ownership, versions, final evidence and product projection.
Before Stage 3 acceptance, the P0 design and plan control only their isolated
non-product preflight lane. Their artifact and oracle outputs remain candidates
until the accepted Stage 3 intake and formal Stage 4 tasks revalidate them.

## Execution Model

- During Stage 3, use the main-agent-led execution protocol frozen in the
  Stage 3 plan.
- The main agent owns task decomposition, full-context implementation,
  requirements review, code-quality review, complete-diff inspection,
  exact-HEAD verification and final acceptance.
- During Stage 3, use a fresh implementation worker only for a bounded,
  independent module or test fixture explicitly marked worker-eligible by the
  Stage 3 plan or the 2026-08-05 execution amendment. Task 11 delegation is
  limited to frozen
  T11-M/T11-E fixtures or mechanical repairs; do not delegate T11-S science,
  complex reviews, cross-module numerical reasoning or long-context plans.
- During the default Stage 4--6 serial path, keep only one implementation
  worker active. If the user explicitly approves a parallel stage edge, use
  separate worktrees and disjoint file ownership exactly as recorded in the
  integration plan.
- Until the user explicitly changes this rule, use the default worker type for
  bounded delegation. Do not dispatch `luna_worker` and do not set a model or
  reasoning-effort override manually.
- Close each worker immediately after acceptance or rejection.
- Workers never contact the user or request user approval.
- Worker blockers return to the coordinator with repository evidence.
- After implementation, the main agent runs a complete requirements review
  followed by a complete code-quality review. Resolve findings and repeat the
  affected review before accepting the task.
- Stage 4--6 execution is serial by default. At each accepted stage node the
  main agent reports the exact HEAD and evidence, recommends serial or limited
  parallel execution, and waits for the user's instruction. A parallel edge in
  a plan is not permission to start it.
- Stage 4--6 workers receive one frozen task and file allowlist. They do not
  change central registries, scientific coupling, task goals, commits or DCO.
  Cross-module science, complete diff, provenance and acceptance remain with
  the main agent.
- The approved pre-Stage-4 P0 lane may overlap Stage 3 only through its
  governance worktree and the exact external generated root
  `/home/wyf/code_dev/.hundun-flow-preflight/stage4-p0`. No other external
  write root is authorized. It uses at most one
  third-party build and small standalone 1/2-rank spikes. Workers do not commit;
  the main agent owns P0 provenance, ABI, mathematics, complete diff and seal.

The tracked parallel-completion v2 activation receipt records the exact
design/plan/reference hashes. A bounded R1/O1 worker packet may now use the
reading exemption defined by that plan: it must read this file completely,
the v2 design, the public-reference document and its exact task block, but
need not reread historical documents that the packet neither changes nor
reinterprets. This exemption does not apply to a main-only task, scientific
composition, complete-diff review, copyright independence or final
acceptance.

## Copyright Boundary

- Fixed legal comparison baseline:
  `/home/wyf/Compress_boffin/AECSC_WDQ/SRC-Spray`.
- BOFFIN is for private independence audit only, not implementation reference.
- COAST and COAST-2 remain private scientific references and are not product
  source ancestors. The only source-level exception is the Stage 5/6 private
  oracle protocol described below.
- Stages 1--3 reuse zero source files from BOFFIN, COAST, or COAST-2.
- Do not copy, translate, mechanically refactor, imitate, or adapt legacy
  source, control flow, ABI, array layout, comments, messages, `.d` inputs,
  Decomp, Restart, or compatibility layers.
- Do not add a public `REUSE_MANIFEST.yaml`.
- Public license is Apache-2.0 with DCO.
- Cantera and its transitive dependencies are governed third-party components,
  not HUNDUN original source. Preserve their upstream identity, licenses,
  revisions, archive/binary hashes and local patch diffs.
- A COAST ESF/TCR or fuel-source oracle may be prepared only after the user
  confirms the exact current realpath and version. Only allowlisted pure
  mathematical modules may be copied into an untracked generated directory
  outside the repository and compiled as a separate-process oracle using
  synthetic inputs. No COAST source, case, data, ABI, message or control flow
  may enter Git, the product, the installed package or public tests.
- `NOTICE` must contain only:

  ```text
  Copyright (c) 2026 WANG YUDONG
  Contact: wangyudong@buaa.edu.cn
  ```

## Workspace Safety

- Do not edit, stage, clean, stop, package, or otherwise disturb
  `/home/wyf/code_dev/Coast_software` or its research cases.
- Do not enter or inspect BOFFIN, COAST, COAST-2, the fixed private baseline,
  private research directories, or old implementation source while developing
  Stage 3. During Stage 5/6, the narrowly defined COAST oracle/fuel audit is
  allowed only after the task-specific confirmation gate; all other private
  paths remain outside the implementation evidence boundary.
- Do not delete or modify research data, and do not inspect, stop, signal, or
  otherwise interfere with research processes.
- Do not publish or push HUNDUN-FLOW until the coordinator reports acceptance
  and the user explicitly authorizes publication.
- Prefer `codegraphf` for later code navigation. Use `rg` for exact text and
  file searches when appropriate.
- Use `apply_patch` for manual edits.
- Preserve unrelated changes and never use destructive Git commands.

## Scope

Stage 3 product code is accepted at
`0cbd3d5bde4be63bc6346b4b32db771d87c59ea2`; its governance seal is
`36bebc292e825fa15272481c6a00c2273fa61ce0`, and its projected product is
`22ed17b438ffbb121ccda97898580183bd0803f8`. These identities and all accepted
Stage 3 numerical, transaction, Restart and MPI contracts are frozen.

The accepted Stage 4 scope is frozen at
`6407cd7c591ce088db7f1dd7e296d77acd18da1c`; its governance seal is the commit
containing the Stage 4 final acceptance report. The active scope is the Stage 5
user decision gate, not Stage 5 implementation. Do not begin Stage 5 without
explicit user authorization.
The 27/32/31 Stage 4--6 tasks and integration gates authorize Linux CPU v1,
bundled Cantera, reacting flow, ESF/TPDF/TCR and dilute single-component
point-parcel spray. Stage 4 starts with formal `4F-0`; P0 evidence may be reused
only after exact-hash intake and remains outside product history.

AMR, moving IBM, dense spray, liquid films, KH--RT, multicomponent droplets,
production GPU, rank-changing Restart and NativeChemistryBackend remain
post-v1 non-goals.

No Python dependency may enter the normal HUNDUN configure, build, install,
test, formal-acceptance or runtime path. The source policy's maintainer-only
Cantera artifact producer is not a normal HUNDUN build path and never enters an
installed package. A private similarity-audit script outside the public
repository may use Python.

# HUNDUN-FLOW Agent Instructions

## Current State

This repository was initialized from a documentation-only handoff skeleton.
Do not reinitialize it. Execute only the approved Stage 0/1 plan after explicit
user authorization, and determine progress from Git history and the coordinator
ledger.

## Required Reading

Read these files in order before changing anything:

1. `docs/handoff/2026-07-16-hundun-flow-agent-handoff.md`
2. `docs/design/2026-07-16-hundun-flow-clean-cpp-solver-design.md`
3. `docs/plans/2026-07-16-hundun-flow-stage0-stage1.md`

The design document controls architecture and scientific scope. The Stage 0/1
plan controls task order, exact interfaces, tests, commits, and acceptance.
The handoff document controls orchestration and legal boundaries.

## Execution Model

- Use subagent-driven development.
- The main agent is coordinator, reviewer, and final verifier only.
- Use a fresh implementation worker for each task, sequentially.
- Keep at most five workers open and only one implementation worker active.
- Close each worker immediately after acceptance or rejection.
- Workers never contact the user or request user approval.
- Worker blockers return to the coordinator with repository evidence.
- After implementation, run a requirement review followed by a code-quality
  review. Resolve findings before accepting the task.

## Copyright Boundary

- Fixed legal comparison baseline:
  `/home/wyf/Compress_boffin/AECSC_WDQ/SRC-Spray`.
- BOFFIN is for private independence audit only, not implementation reference.
- COAST and COAST-2 are private black-box scientific references only.
- Stage 1 reuses zero source files from BOFFIN, COAST, or COAST-2.
- Do not copy, translate, mechanically refactor, imitate, or adapt legacy
  source, control flow, ABI, array layout, comments, messages, `.d` inputs,
  Decomp, Restart, or compatibility layers.
- Do not add a public `REUSE_MANIFEST.yaml`.
- Public license is Apache-2.0 with DCO.
- `NOTICE` must contain only:

  ```text
  Copyright (c) 2026 WANG YUDONG
  Contact: wangyudong@buaa.edu.cn
  ```

## Workspace Safety

- Do not edit, stage, clean, stop, package, or otherwise disturb
  `/home/wyf/code_dev/Coast_software` or its research cases.
- Do not delete research data.
- Do not publish or push HUNDUN-FLOW until the coordinator reports acceptance
  and the user explicitly authorizes publication.
- Prefer `codegraphf` for later code navigation. Use `rg` for exact text and
  file searches when appropriate.
- Use `apply_patch` for manual edits.
- Preserve unrelated changes and never use destructive Git commands.

## Scope

The approved execution scope is Stage 0 and Stage 1 only. Do not implement
variable-density flow, pressure coupling, LES, IBM, chemistry, TPDF-TCR,
spray, WENO, DG, moving walls, or multi-part geometry in this run.

No Python dependency may enter the public build or runtime. A private
similarity-audit script outside the public repository may use Python.

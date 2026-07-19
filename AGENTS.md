# HUNDUN-FLOW Agent Instructions

## Current State

Stage 1 is accepted at the annotated `stage1-runtime` tag. Do not reinitialize
the repository or rewrite its accepted history. Stage 2 variable-density flow
is now authorized only through the approved specification and implementation
plan listed below. Determine progress from Git history and the coordinator
ledger, and do not infer authorization for Stage 3 or later work.

## Required Reading

Read these files in order before changing anything:

1. `docs/handoff/2026-07-16-hundun-flow-agent-handoff.md`
2. `docs/design/2026-07-16-hundun-flow-clean-cpp-solver-design.md`
3. `docs/plans/2026-07-16-hundun-flow-stage0-stage1.md`
4. `docs/superpowers/specs/2026-07-19-hundun-flow-stage2-variable-density-flow-design.md`
5. `docs/superpowers/plans/2026-07-19-hundun-flow-stage2-variable-density-flow.md`

The overall design controls architecture and the complete capability roadmap.
The Stage 2 specification controls the approved scientific and numerical
scope. The Stage 2 plan controls its exact interfaces, Tasks 1--26, hard-gate
order, tests, commits, and acceptance. The Stage 0/1 plan and handoff continue
to control frozen Stage 1 contracts, orchestration, and legal boundaries.

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
- Stage 1 and Stage 2 reuse zero source files from BOFFIN, COAST, or COAST-2.
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
- Do not enter or inspect BOFFIN, COAST, COAST-2, the fixed private baseline,
  private research directories, or old implementation source while developing
  Stage 2. They are outside the implementation evidence boundary.
- Do not delete or modify research data, and do not inspect, stop, signal, or
  otherwise interfere with research processes.
- Do not publish or push HUNDUN-FLOW until the coordinator reports acceptance
  and the user explicitly authorizes publication.
- Prefer `codegraphf` for later code navigation. Use `rg` for exact text and
  file searches when appropriate.
- Use `apply_patch` for manual edits.
- Preserve unrelated changes and never use destructive Git commands.

## Scope

The approved execution scope is exactly Stage 2 Tasks 1--26, in the eight hard
gates and order fixed by the Stage 2 plan. It includes low-Mach
variable-density flow, the constant/material/ideal-gas density gates, PISO and
Rhie--Chow, the five approved body-fitted boundary types, separated structured
topology/geometry, project-owned CPU-reference linear algebra and execution
contracts, adaptive BDF2 trial/retry, Checkpoint v2, curvilinear diagnostics,
MPI invariance, and performance evidence.

Do not enter Stage 3. Do not implement LES, IBM, chemistry, species,
multicomponent thermodynamics, TPDF-TCR, spray, particles, production GPU,
PETSc/HYPRE/AMG, WENO/DG, moving walls, complex thermochemical boundaries,
rank-changing Restart, fully compressible flow, acoustics, or shocks. Device
and GPU-aware paths remain non-production contract test doubles only.

No Python dependency may enter the public build or runtime. A private
similarity-audit script outside the public repository may use Python.

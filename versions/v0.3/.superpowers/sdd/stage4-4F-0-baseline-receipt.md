# HUNDUN-FLOW Stage 4 `4F-0` accepted-baseline intake

记录日期：2026-08-11

```text
schema=hundun.stage4.4f0.accepted_intake.v1
status=INTAKE_PASS
stage4_product_changes=none
stage4_product_implementation=not_started
stage4_scientific_acceptance=not_started
```

本 receipt 完成 Stage 4 开始前的 accepted-state 只读接收，并把 Stage 3、产品投影和
P0 候选证据绑定到精确 Git 对象。它不是 Stage 4 Task `4F-0` 的产品测试提交，也不接受
Stage 4 物理能力。下一位主 agent 必须从本 receipt 继续执行 Stage 4 计划中的 source-policy
RED、capability ledger 和后续 `4F-*` 任务；不得重新解释或重做已接受的 Stage 3 数值工作。

## 1. Accepted identities

```text
accepted_governance_head=36bebc292e825fa15272481c6a00c2273fa61ce0
accepted_governance_tree=897560c30d7d7049a81605a257702b4091a13f25
accepted_governance_parent=0cbd3d5bde4be63bc6346b4b32db771d87c59ea2
accepted_code_head=0cbd3d5bde4be63bc6346b4b32db771d87c59ea2
accepted_code_tree=d50c1236f67bd2bdde58c94a125e530ae0f2ffea
accepted_code_parent=fe9065f8559e1367e8e112505bdd565f108d217f
accepted_product_head=22ed17b438ffbb121ccda97898580183bd0803f8
accepted_product_tree=7fb9ce848238eeab5dc1ad0908092d8d115851b4
accepted_product_parent=ae3d08bbb220d1d3b28ec070d1cba9c33fb85877
accepted_receipt_sha256=eddfc417a91ab70fcd5a3da12b127a940ed24e8ca8416a3ad1820b51655617b7
accepted_receipt_commit=36bebc292e825fa15272481c6a00c2273fa61ce0
product_projection_manifest_sha256=224a3cdbb6fb104ad103256eb0de28732ff5e9dfbe47f8b120460bac2ea25f8c
accepted_version=code:0.1.0,product:0.2.0
accepted_banner=HUNDUN-FLOW 0.2.0
dco_result=PASS
```

| Role | Real path | Branch/topology | Status and remote |
| --- | --- | --- | --- |
| accepted governance seal | `/home/wyf/code_dev/.worktrees/hundun-flow-stage3-framework` | `coast/stage3-framework-completion` at `G` | clean；governance common repo has its existing `origin` |
| accepted tested code | `/home/wyf/code_dev/.worktrees/hundun-flow-stage3-accepted-code-intake` | detached exactly at `C` | clean；read-only intake worktree |
| accepted product | `/home/wyf/code_dev/hundun-flow` | `main` at `P` | clean；no product remote |
| Stage 4 integration | `/home/wyf/code_dev/.worktrees/hundun-flow-stage4-reacting-flow` | `coast/stage4-reacting-flow` | based on `G`; product implementation not started |

`C` retains the governance-tree version source `0.1.0`; the accepted product projection deliberately
overrides `VERSION` to `0.2.0`, as recorded by the projection manifest. Stage 4 must keep this distinction
explicit until its own version receipt; it must not call the governance code banner `0.2.0`.

All five identity commits (`C`, `G`, `P`, P0 seal and the integration merge) contain exactly the existing
authorized `Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>` trailer. The 21 commits after the P0
planning head through the P0 seal also pass the same DCO check. No trailer was synthesized for old history.

## 2. Linked-worktree topology

```text
worktree_git_pointer=gitdir: /home/wyf/code_dev/hundun-flow-governance/.git/worktrees/hundun-flow-stage3-accepted-code-intake
accepted_code_absolute_git_dir=/home/wyf/code_dev/hundun-flow-governance/.git/worktrees/hundun-flow-stage3-accepted-code-intake
accepted_code_common_git_dir=/home/wyf/code_dev/hundun-flow-governance/.git
accepted_code_toplevel=/home/wyf/code_dev/.worktrees/hundun-flow-stage3-accepted-code-intake
```

The pointer is a normal text `.git` file and agrees with Git's absolute paths.

## 3. P0 integration audit

```text
p0_seal=910fb1f7fc3df2e0c596d3682db06db442c03ccf
p0_evidence_candidate=401a5ff285aa1a6283d68a241155bbfab405b9e4
p0_final_receipt_sha256=e584a0f5609f8e7054cfe880f45919b28102ac056812e06a5909d06171a3b068
integration_head=d45ef02706a17f12d38e050497f076cc5002fb51
integration_parents=36bebc292e825fa15272481c6a00c2273fa61ce0,910fb1f7fc3df2e0c596d3682db06db442c03ccf
integration_tree=50efe6741f4b6b6bbf113ce899edac754ede4d10
integration_diff_sha256=d450195159fb706235340e5888758d7b73c228ab1cca6ffd3ee579e722342f01
integration_diff_boundary=23 governance/planning/evidence paths; 11336 insertions; 63 deletions
```

Audit result:

- merge used `--no-ff --no-commit` from accepted governance seal `G`;
- only `AGENTS.md` conflicted; resolution retains `C/G/P`, all Stage 3 authorities and every P0
  limitation while activating the user-approved Stage 4 intake;
- all other P0-touched blobs exactly equal the P0 seal;
- all 272 paths in the Stage 3 product projection manifest match `C` in the integration tree;
- `include/`, `src/`, `tests/`, central CMake, presets, `cmake/` and `VERSION` have zero integration
  changes relative to `C`;
- `git diff --check`, conflict-marker scan and P0 input-manifest validator pass;
- P0 remains governance/preflight history only and has not entered the product repository.

Reusable P0 evidence is limited to source/license hashes, the frozen Linux artifact candidate,
standalone C++ thread/MPI/relocation results, public mathematical vectors and intake command patterns.
It does not prove HUNDUN CMake integration, `ChemistryBackend`, reacting transport, Checkpoint v4,
diagnostics v4, Stage 4 science, real-fuel validity or COAST equivalence.

## 4. Build, install and API inventory

The accepted product `P` was freshly configured from a clean external build root with Clang 15/libc++,
`Release`, and `HUNDUN_BUILD_TESTS=OFF`. This reproduces the accepted product packaging path without
changing any repository file.

```text
accepted_product_scratch_build=/tmp/hundun-stage4-4f0-product-build.xuADvk
accepted_product_scratch_install=/tmp/hundun-stage4-4f0-product-install.z4VwYb
scratch_binary_sha256=2ce19de3f90cc5a67ded8c7f88dab02c8300d519a8d6b70a6e6b1c66b2e3f293
scratch_install_log_sha256=89351de3f125921eaff3f2a2b8a89c7ae9b8bcb9a9fd1b7dfef44e778468b4a7
scratch_install_tree_manifest_sha256=af1719ee741f8c2d272d1c903cfdece5ce506f14d467f75bfa276f2c0c9bb1b3
public_header_inventory_sha256=fc31763a30905ba2121a224a4b6a4f367b09ae46cfd701ed9721a94587a73192
installed_public_header_inventory_sha256=34a088c8f2a849c8634fe9373048f3733a7fc3f10f0969ebe4a94a9b7d39db93
exported_symbol_inventory_sha256=0f8c1a166e95107f31f797a516b75800110b9a233f306f12ee9d236bc0b817fe
cmake_target_inventory_sha256=7738ff549e7885e269f25597bca5153ee73291c8cc49254b4104c3cd91a53bb2
```

The source and installed inventories each contain 80 public headers. The installed binary reports
`HUNDUN-FLOW 0.2.0`. Its Stage 3 Clang/libc++ profile requires the frozen local libc++ directory in
the runtime environment and is not a relocatable GCC/libstdc++ distribution profile. This is an
intake limitation, not a Stage 3 regression; Stage 4 packaging must use its approved Ubuntu 22.04,
GCC 11/libstdc++, ABI=1 bundled profile.

Two retained setup failures are diagnostic evidence, not accepted results:

1. `/tmp/hundun-stage4-4f0-intake.OAzxBg` failed install because the old tests-off build had not built
   the SDK archive; building that target completed the installation graph without source changes.
2. `/tmp/hundun-stage4-4f0-product-build.JNuKGi` failed configure because GNU C was combined with a
   Clang-only `-stdlib=libc++` linker flag; the complete Clang 15 C/C++ profile passed.

## 5. Frozen authority inventory

```text
schema_v1_v2_v3_inventory=82bf5d2385de161678a06b5164404abbf0a36b11c24135895a4ac39ae2660881
checkpoint_v1_v2_v3_inventory=d595855174d3512a85fa8e9c1dc8126380a5a80871e3cdc53bcba732e110fbe3
diagnostics_inventory=91770008229d4f0091ac5889f9063a0e18be84af7fcaf020edf879570afa2f6a
capability_inventory=docs/numerics/stage3-capability-ledger.md@C
retry_rollback_authority=073c43dbde6e6ac47df889987014b4a87b915eedf972e664fd66c8b39d0499ad
collective_failure_authority=hundun::runtime::collective_status
final_flux_authority=863f11d62003e8ef03e5f19ad7b9836c4adb39fc3a67f6356092c3c7fbc173db
piso_corrector_authority=pressure_corrector_count==2U in accepted flow implementations and drivers
```

| Authority class | Accepted Stage 3 owner | Stage 4 action |
| --- | --- | --- |
| root executable/case dispatch | `src/app_main.cpp::run_case`, `dispatch_in_root_config_rank_order`, then `run_flow_case`/`run_immersed_flow_case` | extend the same root; no second executable dispatch |
| schema v1/v2/v3 | `load_resolved_case_v3`, `ResolvedCaseV3`, and existing v1/v2 loaders | add schema v4 through the existing loader/dispatch chain |
| Checkpoint v1/v2/v3 | runtime Restart v1, checkpoint-v2 protocol, `flow::write_checkpoint_v3`/`read_checkpoint_v3` and `CheckpointV3Presence` | add v4 sections through the existing persistence chain; preserve old IDs |
| diagnostics | `DiagnosticModuleKind`, structured record schema v1 and Stage 3 provider inventory (accepted kinds 18--22) | extend existing kind/provider authority; allocate new stable IDs once |
| capability ledger | `docs/numerics/stage3-capability-ledger.md` | create Stage 4 rows without rewriting Stage 3 dispositions |
| retry/rollback | `flow::FlowState::{begin,prepare_commit,publish_commit,rollback}_attempt` plus density-model hooks | chemistry trial state participates in the same transaction |
| collective failure | `runtime::collective_status` and accepted flow collective wrappers | map backend failures into this owner |
| final mass flux | accepted `FaceMassFlux`/`MaterialFaceMassFlux` bound to `state.fields().face_mass_flux` with `final_corrected` provenance | all reacting scalars consume the same final flux |
| PISO count | accepted reports and drivers require exactly `2U` | unchanged, including reacting paths |
| IBM wall/force | Stage 3 immersed domain/reconstruction/wall-force authority | reacting boundary extension only |
| WALE viscosity | accepted `WaleModel` and flow composition | consume through transport service; no duplicate LES evaluation |
| field identity/fingerprint | `FieldRegistry`, accepted state/config/checkpoint fingerprints | add composition/mechanism identity, preserving old identities |

The planned `chem_`, `comb_`, `spray_`, `hundun_chemistry`, `hundun_combustion` and
`hundun_spray` names have zero collisions in `C`. The only `rt_` textual match in the collision
scan is the existing source-policy prefix rule and is classified `same_authority_reuse`.

## 6. Process, limitations and next boundary

```text
background_hundun_jobs=none_active
deferred_scientific_validation=Stage 3 ledger dispositions remain unchanged
stage4_product_changes=none
```

No exact-path HUNDUN build, numerical or MPI process was active at intake completion. Four historical
Task 11/M2 systemd user units remain in failed/inactive state; none was started, stopped or changed.
No unrelated process was inspected or signalled.

The following remain outside v1 or outside the present gate: 96-cubed, AMR, moving/multiple IBM,
rank-changing Restart, production GPU, NativeChemistryBackend, dense spray, multicomponent droplets,
liquid films, KH--RT and unconfirmed private COAST oracle use. Formal COAST access still requires the
plan-specific user confirmation of exact realpath/version.

`INTAKE_PASS` authorizes the next main agent to continue the approved Stage 4 plan from the current
integration branch. It does not authorize a worker to change central ownership, skip mutation REDs,
start Stage 5/6, push, publish or claim reacting-flow capability.

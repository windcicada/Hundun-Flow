# Stage 4 capability ledger

This ledger starts the Stage 4 reacting-flow capability record without changing
any disposition in `docs/numerics/stage3-capability-ledger.md`. At Task `4F-0`,
no Stage 4 product or scientific capability is implemented or accepted.

```text
schema=hundun.stage4.capability_ledger.v1
stage=4
status=STAGE4_ACCEPT
accepted_code_head=6407cd7c591ce088db7f1dd7e296d77acd18da1c
accepted_code_tree=2791a1cee7ac8114f1696670d30c8951212d6024
version=0.3.0
governance_seal=commit-containing-stage4-final-acceptance-report
default_execution=serial
stage4_product_prefix_allowlist=app,bc,cfg,chem,diag,exec,flow,fvm,ib,les,lin,mesh,rt,sdk
plugin_abi_authority=HUNDUN_PLUGIN_METADATA_ABI_V1
```

The allowlist registers the new Stage 4 `chem_` product prefix and reuses the
existing authorities that Stage 4 is allowed to extend. It does not register the
future Stage 5 `comb_` or Stage 6 `spray_` prefixes. The existing `rt_` prefix
remains owned by the accepted runtime authority; this ledger does not create a
second runtime namespace.

## Accepted baseline identities

| Identity | Commit / tree or digest | Stage 4 meaning |
| --- | --- | --- |
| accepted Stage 3 code `C` | `0cbd3d5bde4be63bc6346b4b32db771d87c59ea2` / `d50c1236f67bd2bdde58c94a125e530ae0f2ffea` | frozen product, test, numerical and transaction authority |
| accepted Stage 3 governance `G` | `36bebc292e825fa15272481c6a00c2273fa61ce0` / `897560c30d7d7049a81605a257702b4091a13f25` | accepted governance parent of Stage 4 integration |
| accepted product `P` | `22ed17b438ffbb121ccda97898580183bd0803f8` / `7fb9ce848238eeab5dc1ad0908092d8d115851b4` | clean `0.2.0` product projection; not the development worktree |
| P0 seal | `910fb1f7fc3df2e0c596d3682db06db442c03ccf` | immutable candidate evidence only |
| P0 integration merge | `d45ef02706a17f12d38e050497f076cc5002fb51` / `50efe6741f4b6b6bbf113ce899edac754ede4d10` | signed governance/evidence merge; no product change |
| accepted-state intake | `0f39799692a5da47e9810a4e2262143fb2466996` / `1124e82f63ba61ed1041c9e77218fd08d695f6bb` | completed `4F-0` Step 1 evidence |
| execution handoff baseline | `2907092ff4435538c78fe0aa9fd2960191758261` / `39bf61808d64d9a1593baa220aecbc9d649046c2` | clean signed parent for the remaining `4F-0` work |
| accepted intake receipt | SHA-256 `76e3dc14937083ca7956663bff665891cae28ace5d18be646b7c75e5578e9627` | `INTAKE_PASS`; no Stage 4 product implementation |
| Stage 3 product projection manifest | SHA-256 `224a3cdbb6fb104ad103256eb0de28732ff5e9dfbe47f8b120460bac2ea25f8c` | 272 paths; zero mismatch against `C`, the integration tree and `P` |

The accepted code retains governance `VERSION=0.1.0`; the product projection
deliberately reports `0.2.0`. Stage 4 does not reinterpret that recorded
projection override.

## Frozen authority table

| Authority | Accepted owner | Stage 4 rule |
| --- | --- | --- |
| executable and case dispatch | `src/app_main.cpp::run_case` and the existing root dispatch chain | extend one root; never add a second executable dispatcher |
| schema v1--v3 | existing loaders through `load_resolved_case_v3` / `ResolvedCaseV3` | preserve old schema behavior; schema v4 is not implemented by `4F-0` |
| Checkpoint v1--v3 | existing Restart/checkpoint-v2/flow Checkpoint v3 chain | preserve readers and IDs; Checkpoint v4 is not implemented by `4F-0` |
| diagnostics | accepted provider/kind registry, including stable kinds 18--22 | append once through the existing registry; absence is not a fake-zero record |
| plugin ABI | metadata-only `HUNDUN_PLUGIN_METADATA_ABI_V1` | no model lifecycle ABI and no second plugin ABI |
| retry and rollback | `flow::FlowState` attempt lifecycle and density hooks | Stage 4 trial state must join the same transaction |
| collective failure | `hundun::runtime::collective_status` | backend failures must map into this authority |
| final face mass flux | accepted `FaceMassFlux` / `MaterialFaceMassFlux` with `final_corrected` provenance | every reacting scalar must consume the same final flux |
| pressure correction | exactly `pressure_corrector_count == 2U` | unchanged on every successful reacting path |
| IBM and WALE | accepted Stage 3 immersed geometry/wall-force and `WaleModel` | consume through thin extensions; do not duplicate geometry, pressure or LES evaluation |
| field identity | accepted `FieldRegistry` and state/config/checkpoint fingerprints | add composition/mechanism identity without changing old identities |

## Stage 4 capability rows

`planned-not-implemented` means the approved task exists but no product claim is
available. `deferred` and `out-of-scope` retain the architecture specification's
v1 limits.

| Capability | Disposition at `4F-0` | Planned task(s) | Required acceptance owner |
| --- | --- | --- | --- |
| all-species composition and total thermochemical enthalpy identities | implemented-and-verified-at-C4 | `4F-1` | Stage 4 exact-head seal |
| backend-neutral thermo, transport and chemistry services | implemented-and-verified-at-C4 | `4F-2` | Stage 4 exact-head seal |
| reacting source transaction and atomic rollback | implemented-and-verified-at-C4 | `4F-3` | Stage 4 exact-head seal |
| schema v4 and typed MPI broadcast | implemented-and-verified-at-C4 | `4F-4` | Stage 4 exact-head seal |
| Checkpoint v4 and diagnostics provider contracts | implemented-and-verified-at-C4 | `4F-5` | Stage 4 exact-head seal |
| audited no-Python Cantera package and relocation boundary | implemented-and-verified-at-C4 | `4P-1`--`4P-4` | Stage 4 exact-head seal |
| Cantera runtime, thermo, transport, stiff chemistry and 0D/PSR backend | implemented-and-verified-at-C4 | `4C-1`--`4C-5` | Stage 4 exact-head seal |
| fixed C-T-C and exactly-two-PISO coupling proof | implemented-and-verified-at-C4 | `4R-0` | Stage 4 exact-head seal |
| conservative reacting state, transport and open/closed flow | implemented-and-verified-at-C4 | `4R-1`--`4R-4`, `4R-6` | Stage 4 exact-head seal |
| reacting IBM and WALE composition | implemented-and-verified-at-C4 | `4R-5` | Stage 4 exact-head seal |
| schema-v4 reacting driver combinations | implemented-and-verified-at-C4 | `4R-7` | Stage 4 exact-head seal |
| Checkpoint v4 codecs and reacting diagnostics | implemented-and-verified-at-C4 | `4A-1`, `4A-2` | Stage 4 exact-head seal |
| Stage 4 compact acceptance at version `0.3.0` | accepted | `4A-3`, `4A-4` | tested `C4` plus governance seal |
| ESF/TPDF/TCR | deferred | Stage 5 after explicit user instruction | Stage 5 exact-head seal |
| dilute point-parcel spray | deferred | Stage 6 after explicit user instruction | Stage 6 development-complete seal |
| NativeChemistryBackend | out-of-scope | post-v1 | future independent plan |
| AMR, moving IBM, dense spray and rank-changing Restart | out-of-scope | post-v1 | future independent plan |
| 96-cubed execution | out-of-scope | permanent v1 exclusion | HUNDUN-FLOW governance |

## Execution preconditions and evidence boundary

- Continue serially as `4F contracts -> 4P package -> 4C backend -> 4R reacting
  -> 4A acceptance`. Completing `4F-0` does not authorize `4F-1` in the same
  task node.
- Stage 4 acceptance must stop for a user decision before Stage 5; Stage 5
  acceptance must stop again before Stage 6.
- The release profile is Ubuntu 22.04/glibc 2.35+, GCC 11/libstdc++, C++17,
  `_GLIBCXX_USE_CXX11_ABI=1`, generic x86-64. Normal configure, build,
  install, test and runtime paths require no Python, Conda or online fetch.
- Public headers expose only HUNDUN value/config/report types. Cantera,
  SUNDIALS, Python and pybind types remain outside the public API.
- P0 may be reused only for hash-bound source/license identities, the artifact
  candidate, standalone thread/MPI/relocation evidence, public mathematical
  vectors and intake command patterns. Formal Stage 4 tasks must establish all
  HUNDUN integration, backend, reacting-flow, persistence, diagnostics, driver
  and scientific evidence.
- BOFFIN and private research data remain outside the implementation boundary.
  COAST source access is not authorized in Stage 4. The later Stage 5/6 oracle
  exception requires explicit confirmation of the exact realpath/version and a
  generated, untracked, separate-process boundary.
- No Stage 4 capability becomes public or accepted until its mutation RED,
  focused tests, caller/requirements/quality/complete-diff review, provenance,
  DCO and exact-head receipt all pass.

# Stage 3 capability ledger

This ledger bounds every public Stage 3 claim. Disposition is exactly
`implemented-and-accepted`, `deferred`, or `out-of-scope`; the final owner is HUNDUN-FLOW unless
the row explicitly says otherwise.

| Capability | Disposition | Task / signed commit | Acceptance test | Final owner |
| --- | --- | --- | --- | --- |
| profile-1 constant body-fitted no-WALE | implemented-and-accepted | S3-A1 `3c1b8dc` | `profiles-fast-r{1,2,4}` | HUNDUN-FLOW |
| profile-2 material body-fitted no-WALE | implemented-and-accepted | S3-A1 `3c1b8dc` | `profiles-fast-r{1,2,4}` | HUNDUN-FLOW |
| profile-3 ideal-gas body-fitted no-WALE | implemented-and-accepted | S3-A1 `3c1b8dc` | `profiles-fast-r{1,2,4}` | HUNDUN-FLOW |
| profile-4 constant static-IBM no-WALE | implemented-and-accepted | S3-A1 `3c1b8dc` | `profiles-fast-r{1,2,4}` | HUNDUN-FLOW |
| profile-5 material static-IBM no-WALE | implemented-and-accepted | S3-A1 `3c1b8dc` | `profiles-fast-r{1,2,4}` | HUNDUN-FLOW |
| profile-6 ideal-gas static-IBM no-WALE | implemented-and-accepted | S3-A1 `3c1b8dc` | `profiles-fast-r{1,2,4}` | HUNDUN-FLOW |
| profile-7 constant static-IBM WALE | implemented-and-accepted | S3-C1/A1 `3c1b8dc` | `constant-ibm-wale-n24-r{1,2,4}` | HUNDUN-FLOW |
| profile-8 material static-IBM WALE | implemented-and-accepted | S3-C2/A1 `3c1b8dc` | `material-ibm-wale-n{12,24}-r{1,2,4}` | HUNDUN-FLOW |
| profile-9 ideal-gas static-IBM WALE | implemented-and-accepted | S3-C3/A1 `3c1b8dc` | `ideal-ibm-wale-n{12,24}-r{1,2,4}` | HUNDUN-FLOW |
| deterministic STL surface/query/classification | implemented-and-accepted | Tasks 3--5 / S3-O2 `0414194` | `diagnostics-fast-r{1,2}` | HUNDUN-FLOW |
| LFP ghost-cell reconstruction and static wall force | implemented-and-accepted | Tasks 6--11 / S3-O2 `0414194` | `task11-authority-current-tree` | HUNDUN-FLOW |
| body-fitted and immersed WALE | implemented-and-accepted | Tasks 12--13 / S3-O1 `463d048` | `wale-tgv-convergence-r1` | HUNDUN-FLOW |
| Checkpoint v3 profiles 1--9 | implemented-and-accepted | S3-R1/R2 `8f869aa` | `checkpoint-continuation-n12-r{1,2,4}` | HUNDUN-FLOW |
| diagnostics kinds 18--22 | implemented-and-accepted | S3-O1/O2 `0414194` | `diagnostics-fast-r{1,2,4}` | HUNDUN-FLOW |
| exact Stage 3 work counters / artifact schema v2 | implemented-and-accepted | S3-E1 `d1f866c` | `stage3-exact-counters-n24-r{1,2,4}` | HUNDUN-FLOW |
| redistribution | deferred | Stage 4 | none in Stage 3 | future coordinated task |
| multigrid | deferred | Stage 4 | none in Stage 3 | future coordinated task |
| AMR | out-of-scope | Stage 5+ | none in Stage 3 | future coordinated task |
| GPU backend | out-of-scope | Stage 5+ | none in Stage 3 | future coordinated task |
| moving immersed bodies | deferred | Stage 5+ | none in Stage 3 | future coordinated task |
| rank-changing Restart | out-of-scope | Stage 4+ | none in Stage 3 | future coordinated task |
| 96-cubed execution | out-of-scope | permanent Stage 3 exclusion | none | HUNDUN-FLOW governance |

The formal 24/48-cubed rows are candidate gates. A public 0.2.0 release claim is valid only after
S3-V1 returns ACCEPT for one frozen exact HEAD. No upstream or private implementation is an
acceptance authority.

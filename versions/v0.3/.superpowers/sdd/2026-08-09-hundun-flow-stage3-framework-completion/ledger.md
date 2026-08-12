# Stage 3 framework completion ledger

| Cluster | Status | Implementation owner | Requirements review | Quality/caller review | Commit/evidence |
| --- | --- | --- | --- | --- | --- |
| Repository split and flat layout | accepted | main + bounded Luna evidence | PASS | PASS | governance `ee4d2b1`; product `ae3d08b`; pre/post switch seals |
| 19A constant IBM driver | accepted | main agent | PASS | PASS | signed `17b8434`; 17A adds Restart wiring |
| 17A Checkpoint v3 IBM-only | accepted | main agent + Luna read-only review | PASS | PASS | this signed task commit; Task 17A receipt |
| 18A minimal diagnostics/counters | accepted | main agent | PASS | PASS | this signed task commit; Task 18A receipt |
| MVP milestone | accepted internal | main agent | PASS | PASS | MVP milestone receipt; test/governance commit |
| 12 WALE core | accepted | main agent | PASS | PASS | signed task commit; Task 12 receipt |
| 13 + 19B body-fitted WALE driver | accepted | main agent | PASS | PASS | this signed task commit; Task 13+19B receipt |
| 14 material IBM | superseded-by-stage3-parallel-framework-v2 | main agent | n/a | n/a | mapped to S3-D1 |
| 15 ideal-gas IBM | superseded-by-stage3-parallel-framework-v2 | main agent | n/a | n/a | mapped to S3-D2 |
| 16 combined IBM+WALE | superseded-by-stage3-parallel-framework-v2 | main agent | n/a | n/a | mapped to S3-C1/C2/C3/S1 |
| 17B/18B/19C framework completion | superseded-by-stage3-parallel-framework-v2 | packet-specific | n/a | n/a | mapped to S3-R1/R2/O1/O2/A1 |
| 20 ledger/counters/launcher | superseded-by-stage3-parallel-framework-v2 | main agent | n/a | n/a | mapped to S3-E1/G1 |
| Public documentation finalization | superseded-by-stage3-parallel-framework-v2 | main agent | n/a | n/a | mapped to S3-DOC |
| 21 exact-HEAD acceptance | superseded-by-stage3-parallel-framework-v2 | main agent | n/a | n/a | mapped to S3-V0/V1/V2 |

## Active parallel-completion v2 mapping

The user approved the immutable v2 candidate. The tracked
`stage3-v2-activation.md` receipt binds its exact design/reference/plan hashes.
These rows now control work after accepted Task 13+19B.

| v2 task | Legacy mapping | Status | Owner | Long-test owner |
| --- | --- | --- | --- | --- |
| S3-A0 immutable candidate activation | governance prerequisite after user approval | accepted | main agent | none |
| S3-P0 parallel execution foundation | governance/test registration prerequisite | accepted | main agent | none; this signed P0 commit and receipt |
| S3-C1 constant IBM+WALE | Task 16 constant subcluster | accepted | main agent | S3-V1; this signed C1 commit and receipt |
| S3-D1 material IBM vertical slice | Task 14 + driver subset | accepted | main agent | S3-V1; this signed D1 commit and receipt |
| S3-C2 material IBM+WALE | Task 16 material subcluster | accepted | main agent | S3-V1; this signed C2 commit and receipt |
| S3-D2 ideal-gas IBM vertical slice | Task 15 + driver subset | accepted | main agent | S3-V1; this signed D2 commit and receipt |
| S3-C3 ideal-gas IBM+WALE | Task 16 ideal subcluster / Gate 5 | accepted | main agent | S3-V1; this signed C3 commit and receipt |
| S3-S1 final scientific selector freeze | Task 16/20 final selector ownership | accepted | main agent | S3-V1; this signed S1 commit and receipt |
| S3-R1 Checkpoint constant profiles | Task 17B constant subset | accepted | bounded default worker in isolated lane; main reviewed/signed/integrated | S3-V1; handoff `c2a3b71`; integrated `efcdd25`; R1 receipt |
| S3-O1 WALE diagnostics | Task 18B WALE subset | accepted | bounded default worker in isolated lane; main reviewed/signed/integrated | S3-V1; handoff `8f927a0`; integrated `463d048`; O1 receipt |
| S3-R2 Checkpoint density profiles | Task 17B density/driver subset | accepted | main agent | S3-V1; signed `8f869aa`; R2 receipt |
| S3-O2 Stage 3 providers | Task 18B provider subset | accepted | main agent | S3-V1; signed `0414194`; O2 receipt |
| S3-A1 driver matrix | Task 19C | accepted | main agent | S3-V1; signed `3c1b8dc`; A1 receipt |
| S3-E1 exact counters/performance artifact | Task 18B counters + Task 20 performance | accepted | main agent | S3-V1 24-cubed; signed `d1f866c`; E1 receipt |
| S3-G1 capability/inventory/manifest | Task 20 ledger/launcher | accepted | main agent | signed `cdbc794`; G1 receipt |
| S3-DOC public documentation | documentation contract + public 0.2.0 candidate docs | accepted | main agent | signed `7abe828`; DOC receipt |
| S3-V0 | frozen-candidate preflight | accepted | main agent | exact code candidate `0cbd3d5`; low-cost/sanitizer/governance 31/31; prior rejected-candidate evidence retained only as history |
| S3-V1 WALE TGV repair | bounded acceptance repair | accepted | main agent | repair packet/receipt; signed repair commit; rerun only after new V0 freeze |
| S3-V1 Checkpoint inventory repair | bounded G1 identity repair | accepted | main agent | repair packet/receipt; signed repair commit; rerun only after new V0 freeze |
| S3-V1/V2 | Task 21 + product projection | accepted | main agent | exact-C scientific 23/23 + performance 3/3; total 57/57; product `22ed17b`; final exact-HEAD seal |

Only one implementation worker may be active. For the current task and until
the user changes the preference again, use only the default worker without a
manual model or reasoning override; do not dispatch `luna_worker`. No worker
may expand a frozen packet, sign commits or decide cross-module science.

# S3-V1 Checkpoint inventory repair packet

状态：`BOUNDED_REPAIR`

rejected candidate：`f69792e938eb33ae2a8f1798516aa2ae10bf20b2`

owner：main agent

## Failure and root cause

The single S3-V1 scientific launcher completed 22 of 23 required rows with
exit zero. The final `checkpoint-continuation-n12-r4` row exited 8 before the
test executable started because the inventory selected CTest
`test_checkpoint_v3_4_rank` and named a nonexistent producer target
`test_checkpoint_v3`. That CTest actually launches `test_immersed_transaction`,
which was not part of the frozen Release targeted build.

The same inventory cluster also bound the 1/2-rank rows to fast wrapper tests.
The authoritative Stage 3 registration already provides the intended formal
rows `checkpoint-continuation-n12-r{1,2,4}` through
`test_checkpoint_v3_density_profiles`, and that target is already part of the
frozen V0 Release build command. The defect is therefore limited to G1
inventory identity; product code, Checkpoint v3 semantics and CTest selectors
are unchanged.

## File boundary

Only these tracked files may change:

- `tests/acceptance/stage3_acceptance_inventory.tsv`;
- `tests/cmake/stage3_acceptance_contract.cmake`;
- this packet and its receipt;
- `.superpowers/sdd/2026-08-09-hundun-flow-stage3-framework-completion/ledger.md`.

Do not modify `tests/CMakeLists.txt`, Stage 3 registration fragments, product
code, Restart format or semantics, thresholds, PISO, force, rollback, or MPI
requirements.

## Required RED/GREEN

RED must compile and execute. The acceptance contract must reject the current
inventory because the three formal continuation rows are not bound to their
registered row IDs and producer, and because `restart-fast-r4` does not name
its actual producer.

GREEN requires:

- the acceptance contract and its existing mutations;
- inventory list-only review;
- Release targeted build of `test_checkpoint_v3_density_profiles`;
- the three formal continuation rows on 1/2/4 ranks;
- a restarted V0 and a new exact candidate before another V1 run.

All 22 successful scientific manifests from the rejected candidate remain
historical evidence only and are not silently reused for a changed candidate.

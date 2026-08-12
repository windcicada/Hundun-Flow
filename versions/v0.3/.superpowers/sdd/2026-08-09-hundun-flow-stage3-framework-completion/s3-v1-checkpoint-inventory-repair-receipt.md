# S3-V1 Checkpoint inventory repair receipt

状态：`ACCEPTED`

accepted parent：`f69792e938eb33ae2a8f1798516aa2ae10bf20b2`

## Result

The G1 inventory now binds all three formal Checkpoint continuation rows to
their registered CTest names and actual Release producer
`test_checkpoint_v3_density_profiles`. The low-cost legacy 4-rank row now names
its actual producer `test_immersed_transaction`, so evidence manifests cannot
silently substitute `CMakeCache.txt` for a missing producer binary.

No product source, test executable, CTest registration, Checkpoint format,
Restart semantic, threshold, PISO, force, rollback or MPI rule changed.

## TDD and focused evidence

Executable RED:

```text
test_stage3_acceptance_contract
  FAIL: Stage 3 restart-fast-r4 producer binding differs
```

GREEN on the final bytes:

```text
test_stage3_acceptance_contract          PASS
checkpoint-continuation-n12-r1           PASS 27.81 s
checkpoint-continuation-n12-r2           PASS 14.83 s
checkpoint-continuation-n12-r4           PASS  9.38 s
```

The rejected candidate completed 22 of 23 scientific rows with exit zero. Its
last row failed only because MPI could not execute the unbuilt legacy target;
performance was not run. Those manifests remain rejected history. The accepted
repair must restart V0 and produce a new clean exact candidate before V1.

No private source or research data was accessed; no research process was
inspected or signalled; no push, publication, 96-cubed run or Stage 4--6 action
occurred.

提交 subject：`fix: bind Stage 3 checkpoint acceptance rows`

DCO：`Signed-off-by: WANG YUDONG <wangyudong@buaa.edu.cn>`

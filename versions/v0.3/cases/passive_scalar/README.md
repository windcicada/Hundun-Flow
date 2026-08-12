# Stage 1 periodic passive scalar

Run the uninterrupted two-rank case from a configured build tree:

```sh
mpiexec -n 2 ./hundun /path/to/cases/passive_scalar/case.json
```

The run advances to absolute target step 20. It writes per-rank VTK files at
steps 10 and 20 below `output/` and complete Restart v1 checkpoints below
`Restart/step00000010/` and `Restart/step00000020/`.

`case_restart.json` reads the exact step-10 checkpoint, advances steps 11
through 20, writes VTK files below `output.resumed/`, and writes checkpoints
below `Restart.resumed/`. Restart v1 supports only the same rank count,
process grid, owned partition, and persistent-field schema. It does not
repartition data and makes no power-loss durability claim.

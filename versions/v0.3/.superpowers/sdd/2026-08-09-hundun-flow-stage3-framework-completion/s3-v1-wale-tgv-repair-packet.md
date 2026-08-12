# S3-V1 WALE Taylor--Green repair packet

状态：`BOUNDED_REPAIR`

rejected candidate：`4e0abc2cbb1e0731dc0234903d11f8707997c6b2`

owner：main agent

## Failure and scope

The first S3-V1 scientific run rejected the candidate before later rows were
allowed to continue:

- `wale-tgv-convergence-r1` failed the frozen first velocity Richardson
  segment (`1.7465267 < 1.8`);
- `wale-tgv-n24-r2` and `wale-tgv-n24-r4` failed the frozen `5e-12`
  decomposition comparison;
- `wale-channel-n48-r1` passed.

The launcher was stopped after rejection. No performance row ran, and no later
scientific result from the rejected candidate is reusable. The prior V0
manifests remain historical evidence only.

This repair may modify only:

- `src/flow_body_fitted_wale_detail.hpp`;
- `src/fvm_cell_centered.cpp`;
- `tests/mpi/test_cell_centered_fvm.cpp`;
- `tests/mpi/test_wale_body_fitted.cpp` (accepted deterministic snapshot only);
- this packet, its receipt, and `ledger.md`.

It must not change the TGV grids, final time, cell-average initial values or
restriction, either `>= 1.8` segment, the `5e-12` decomposition threshold,
WALE controls, solver, pressure-corrector count, force authority, rollback,
Restart, or MPI consistency requirements.

## Root causes and repair

The decomposition failure came from interpolating `mu_sgs` before its periodic
and remote ghost cells had a halo exchange. On a single-rank periodic boundary,
the ghost value remained zero; on multiple ranks it was populated remotely.
The repair releases the cell-viscosity write view, exchanges that scratch field,
then acquires the read view used by face interpolation.

The convergence failure was independent of solver tolerance, timestep scale,
molecular viscosity, velocity amplitude and WALE coefficient. The uniform-grid
momentum path reconstructed a face value by arithmetic averaging two quantities
whose declared function space is `cell_average`. That interpolation has a large
coarse-grid truncation term for the frozen 12/24/48 TGV sequence. For a uniform
spacing authority, the repair uses the symmetric four-cell finite-volume
formula

```text
q_f = (7(q_P + q_N) - (q_{P-1} + q_{N+1})) / 12.
```

It exactly reproduces cubic one-dimensional cell averages at the face and is
used only when `MeshGeometry::uniform_spacing_m()` is present and the existing
quality gate does not select MC limiting. Non-uniform and warped mappings, the
MC thresholds and physical-boundary handling retain their existing paths.

## Required RED/GREEN

RED must compile and execute:

- the unchanged formal TGV convergence row rejects the original interpolation;
- a uniform-grid quadratic cell-average face oracle rejects arithmetic
  averaging;
- the original 2/4-rank TGV comparison rejects the missing SGS halo exchange.

GREEN requires:

- the quadratic cell-average oracle on 1/2/4 ranks;
- the unchanged TGV convergence row and 24-cubed 1/2/4-rank rows;
- the existing Stage 2 Taylor--Green core row;
- affected source-policy/header and focused sanitizer rows in the restarted
  V0 inventory.

After this packet is accepted in one signed repair commit, restart S3-V0 Step 1
and freeze a new exact candidate before any S3-V1 rerun.

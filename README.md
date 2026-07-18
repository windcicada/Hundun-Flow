# HUNDUN-FLOW

HUNDUN-FLOW is an extensible C++ framework for immersed reacting-flow simulation.

The name draws on the archaic image of Hundun associated with Dijiang in the
Classic of Mountains and Seas: a pouch-like form, red as cinnabar fire. The
project reinterprets that image as a bounded chamber containing intense
combustion, turbulent transport, and interacting physical processes.

This engineering interpretation is not presented as a literal translation.
The primary text reference is:
https://ctext.org/text.pl?if=gb&node=83583&show=parallel

Stage 1 is complete at the repository gate. It provides an independently
implemented C++17/MPI-3 runtime and a conservative passive-scalar verification
kernel. Reacting flow, variable-density pressure coupling, LES, IBM,
TPDF-TCR, chemistry, particles, spray, WENO, DG, multi-part geometry, moving
walls, and accelerator implementations are roadmap items, not current
capabilities.

## Build and test

After installing a C++17 compiler, CMake 3.21 or newer, and an MPI-3
implementation, the vendored source tree builds without network access:

```sh
cmake --preset release
cmake --build --preset release -j 2
ctest --preset release --output-on-failure
bash tests/acceptance/stage1_acceptance.sh
```

The public configure, build, test, and run paths have no Python dependency.

## Run the Stage 1 case

```sh
mpiexec -n 2 build/release/hundun cases/passive_scalar/case.json
```

The canonical case writes per-rank legacy VTK scalar fields below
`cases/passive_scalar/output/` and collective Restart v1 checkpoints below
`cases/passive_scalar/Restart/step00000010/` and `step00000020/`. Each complete
checkpoint contains a manifest, two checksummed rank files, and a marker
written only after all files succeed. See
`cases/passive_scalar/README.md` for the exact continuation command and the
same-partition Restart limitation.

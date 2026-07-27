# HUNDUN-FLOW

HUNDUN-FLOW is a clean-room C++17/MPI-3 research CFD solver for
low-Mach-number variable-density flow. The current public baseline is
**0.1.0-stage2**.

中文简介：HUNDUN-FLOW 是面向低马赫数变密度流动的 C++17/MPI-3 开源 CFD
研究软件。Stage 2 已完成，包含常密度、材料密度和理想气体闭合路径，以及可复现的
MPI 验证与验收流程。

The project name draws on the image of Hundun and Dijiang in the *Classic of
Mountains and Seas*. It is an engineering name, not a literal translation or
a claim of historical interpretation.

## Project status

Stage 2 is the current implementation baseline. It delivers a CPU-reference
solver path and an evidence-led verification suite for structured,
low-Mach-number flow. The implementation and its accepted test matrix are
indexed in the [Stage 2 capability ledger](docs/numerics/stage2-capability-ledger.md).

The repository remains an early research release. It is intended for
reproducible numerical research and software development, not as a
drop-in production CFD package.

## What is available

- C++17, CMake 3.21+, and MPI-3 build with no Python runtime dependency.
- Structured-domain decomposition, MPI halo exchange, typed field access with
  epoch/capability checks, JSON case validation, and deterministic diagnostics.
- Project-owned CG and BiCGStab solvers, matrix-free Poisson operators, and
  cell-centred finite-volume transport.
- PISO pressure-velocity coupling with time-consistent Rhie--Chow treatment
  and two pressure correctors.
- Constant-density, material-density, and ideal-gas density-closure paths.
- Adaptive BDF2 retry/rollback control and same-layout transactional
  Checkpoint v2 restart.
- MPI verification at 1, 2, and 4 ranks, including conservation,
  decomposition invariance, rejected-state/rollback, diagnostic, and portable
  counter checks.
- The original periodic passive-scalar application, including VTK output and
  same-partition Restart v1, remains available as a frozen regression path.

## Requirements

- A C++17 compiler (GCC or Clang are covered by CI)
- CMake 3.21 or newer
- An MPI-3 implementation with `mpiexec`
- A POSIX shell for the acceptance scripts

## Quick start

Configure, build, and run the complete CTest suite:

```sh
cmake --preset release
cmake --build --preset release -j 2
ctest --preset release --output-on-failure
```

Run the public Stage 1 passive-scalar case:

```sh
mpiexec -n 2 build/release/hundun cases/passive_scalar/case.json
```

Run a Stage 2 case with one of the supported density models:

```sh
mpiexec -n 2 build/release/hundun /path/to/stage2-case.json
```

Stage 2 uses schema version 2 with
`simulation.type: "variable_density_flow"` and one of
`constant`, `material`, or `ideal_gas` density models. The repository's
acceptance scripts construct their canonical cases directly so the numerical
fixtures remain controlled by the test suite.

## Verification

Run the frozen Stage 1 application acceptance:

```sh
bash tests/acceptance/stage1_acceptance.sh
```

Run the selected Stage 2 acceptance matrix after a release build:

```sh
HUNDUN_STAGE2_BUILD_DIR=build/release \
  bash tests/acceptance/stage2_acceptance.sh
```

The Stage 2 gate verifies its exact registered test inventory before executing
the matrix serially. GitHub Actions builds GCC and Clang debug/release
configurations, runs CTest, runs the Stage 1 acceptance, and runs the Stage 2
acceptance on the release configuration.

## Repository layout

| Path | Purpose |
| --- | --- |
| `applications/` | Command-line application, Stage 1 regression driver, and Stage 2 flow driver |
| `config/` | Typed JSON schema validation and resolved-case serialization |
| `runtime/`, `execution/` | MPI lifetime, field storage, Halo exchange, and execution contracts |
| `mesh/`, `boundary/` | Structured topology/geometry and the approved boundary conditions |
| `linear/`, `finite_volume/` | Project-owned linear solvers and finite-volume/Poisson operators |
| `flow/` | PISO, density closures, time control, and Checkpoint v2 |
| `diagnostics/` | Deterministic diagnostic sessions, mesh records, and performance artifacts |
| `cases/` | Public passive-scalar regression inputs |
| `tests/` | Unit, MPI, contract, and acceptance tests |
| `sdk/` | Metadata-only plugin compatibility interface |

## Current limitations

The following are intentionally outside the Stage 2 scope:

- Reacting chemistry, species transport, LES, immersed boundaries, TPDF/TCR,
  spray, particles, moving walls, WENO, DG, AMR, and unstructured meshes.
- Production GPU execution, vendor solver backends, and a general model
  callback ABI.
- Rank-changing restarts; Checkpoint v2 requires identical rank count,
  process grid, and owned-box layout.
- Portable wall-clock, RSS, bandwidth, or throughput pass thresholds. The
  project gates exact portable counters and numerical correctness instead.

## License

HUNDUN-FLOW is distributed under the [Apache License 2.0](LICENSE). See
[THIRD_PARTY.md](THIRD_PARTY.md) and [LICENSES/](LICENSES/) for bundled
third-party notices.

## Contributing

Please keep changes independently implemented and preserve the repository's
source-policy and numerical-test contracts. Contributions are subject to the
[Developer Certificate of Origin](DCO.md).

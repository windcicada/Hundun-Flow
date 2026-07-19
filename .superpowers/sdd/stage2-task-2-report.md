# Stage 2 Task 2 Implementation Report

## Status and scope

Implemented Gate 1 Task 2 on accepted base
`030db24b0bd715a23744f0d70a5340bbaa07a97b` only. The change adds typed
schema-v2 configuration parsing, deterministic JSON I/O, typed MPI process
data transfer, and same-executable dispatch. Schema v2 normal execution is
explicitly deferred until Task 24. No mesh, field, execution, linear-system,
finite-volume, solver, or Checkpoint-v2 implementation was added.

Schema-v1 `CaseConfig`, `load_case_config()`, `validate_case_config()`,
`to_resolved_json(const CaseConfig&)`, `broadcast_case_config()`, Restart v1,
primitive VTK, version text, and the Stage 1 numerical body remain unchanged.
The numerical body was moved without kernel edits into a focused Stage 1
driver.

## Changed files

- `CMakeLists.txt`
- `CMakePresets.json`
- `applications/hundun/case_config_broadcast.cpp`
- `applications/hundun/case_config_broadcast.hpp`
- `applications/hundun/main.cpp`
- `applications/hundun/stage1_driver.cpp`
- `applications/hundun/stage1_driver.hpp`
- `config/include/hundun/config/resolved_case.hpp`
- `config/include/hundun/config/resolved_case_loader.hpp`
- `config/src/resolved_case_loader.cpp`
- `tests/acceptance/stage2_task2_dispatch.sh`
- `tests/mpi/test_resolved_case_broadcast.cpp`
- `tests/unit/test_resolved_case.cpp`
- `.superpowers/sdd/stage2-task-2-report.md`

## TDD evidence

RED was recorded before product implementation:

```text
/tmp/hundun-toolchain/cmake/bin/cmake --preset debug &&
/tmp/hundun-toolchain/cmake/bin/cmake --build --preset debug \
  --target test_resolved_case -j2
```

Configure succeeded, then compilation failed because the frozen API header
`hundun/config/resolved_case_loader.hpp` did not exist. Exit code: 2. This was
the expected RED for missing ResolvedCase/schema-v2 behavior, not a test
fixture error.

Focused GREEN:

```text
/tmp/hundun-toolchain/cmake/bin/cmake --build --preset debug \
  --target test_resolved_case -j2
/tmp/hundun-toolchain/cmake/bin/ctest --test-dir build/debug \
  -R '^test_resolved_case$' --output-on-failure
```

Result: build exit 0; 1/1 focused configuration test passed. The typed
broadcast and dispatch targets also compile successfully. MPI execution itself
is blocked by the local sandbox as described below.

## Configuration parsing and JSON I/O coverage

`test_resolved_case` covers:

- schema-v1 ResolvedCase load and byte-identical accepted resolved JSON;
- valid closed constant, material, and ideal-gas cases plus an open ideal-gas
  case;
- fixed/adaptive mode, uniform/analytic mapping, every density and boundary
  enum family, and both inlet thermal authorities;
- canonical root order, lexical scalar/inlet-scalar order, stable six-patch
  order, empty resolved resources, and byte-identical round trip;
- conditional warp, ideal-gas, Restart, inlet, and outlet members;
- numeric domains, fixed controller constants, cell/process-grid constraints,
  scalar names/duplicates/reserved names, diffusivity, intervals, and path
  containment;
- patch uniqueness, periodic pairing, open-boundary counts, scalar value sets,
  and constant/material/ideal-gas thermal-density cross-checks;
- unknown and duplicate key families, invalid types, malformed/non-finite JSON,
  and stable JSON-pointer diagnostics.

The implementation uses typed enums and optionals only for conditional JSON
members. Parsing canonicalizes general scalars, inlet scalar values, and the
six patches. JSON input remains strict and standards-conforming.

## Typed process data transfer and dispatch coverage

`test_resolved_case_broadcast` is registered for 1/2/4 ranks and covers
schema-v1/schema-v2 reconstructed JSON equality, invalid/null/intercommunicator
handles, inconsistent and out-of-range roots, root/non-root pointer discipline,
communicator/process-grid mismatch, and finalized MPI. The implementation
transfers a discriminant and every typed member, uses checked uint64 lengths
and chunked byte transfer for strings/paths, validates optionals/enums/counts,
and converges reconstruction/validation failures to the lowest relative rank.
It does not broadcast raw C++ object bytes.

`stage2_task2_dispatch.sh` covers schema-v2 `--validate`, canonical
`--print-resolved` round trip, no Stage 1 banner, and the stable normal-run
deferral:

```text
Stage 2 variable-density flow driver is not implemented before Task 24
```

## Debug, UBSan, and Stage 1 evidence

Debug:

- configure/build all targets: exit 0, no compiler warnings;
- `ctest -L unit -j1`: 10/10 passed;
- full serial `ctest`: 51/93 passed and 42 failed before test logic because
  OpenMPI could not create local communication sockets (`bind/socket` returned
  `EPERM`). Every failure was an MPI-initializing numerical, MPI, or dispatch
  test; non-MPI unit, acceptance-completeness, provenance, and source-policy
  tests passed.

UBSan:

- `ubsan` configure preset is present with Debug/test semantics and only
  `HUNDUN_ENABLE_UBSAN=ON` added;
- full UBSan build: exit 0;
- focused schema-v2 test: 1/1 passed without sanitizer output;
- UBSan `ctest -L unit -j1`: 10/10 passed without sanitizer output;
- focused/full MPI execution could not start because of the same sandbox
  socket restriction.

Stage 1:

- `build/debug/hundun --version` printed exactly
  `HUNDUN-FLOW 0.0.0-stage1`;
- the accepted Stage 1 script was run from a fresh temporary Clang 15/libc++
  build using the repository's required compiler/link flags; the clean build
  completed, then the first MPMD MPI gate could not produce its project
  diagnostic because MPI execution is unavailable in this sandbox;
- the Stage 1 numerical code moved to `stage1_driver.cpp` is otherwise
  byte-for-byte unchanged apart from receiving the already parsed typed v1
  configuration and referencing the same version constant.

The coordinator must rerun focused MPI, full Debug/UBSan, dispatch, and the
complete Stage 1 acceptance gate where local MPI communication is permitted.

## Source policy, provenance, and self-review

```text
/tmp/hundun-toolchain/cmake/bin/ctest --test-dir build/debug \
  -L provenance -j1 --output-on-failure
```

Result: 39/39 passed, including all 32 source-policy cases. Public build and
runtime paths add no Python, network fetch, dependency, device/vendor type, or
private-source reference.

Self-review checked the complete file boundary and confirmed:

- the old configuration loader/types and Stage 1 I/O formats are untouched;
- resolved JSON ordering is deterministic and round trips byte-identically;
- typed MPI transfer validates communicator/root/pointer state before payload
  collectives and converges rank-local preparation/allocation/reconstruction
  failures;
- v2 normal execution creates no deferred subsystem and prints no capability
  banner;
- the existing debug/release/asan preset meanings are unchanged;
- no Task 1 diagnostics, approved specification/plan, AGENTS, tag, remote, or
  private directory was modified.

`git diff --check` and staged source-policy/provenance checks are run again
immediately before the commit.

## DCO and concerns

The single commit uses subject
`feat: add stage two resolved case dispatch` and the configured DCO sign-off
`WANG YUDONG <wangyudong@buaa.edu.cn>`.

Concern: local MPI verification is incomplete solely because the managed
sandbox denies OpenMPI local socket setup. The coordinator must provide the
independent MPI and Stage 1 acceptance evidence before accepting this task.

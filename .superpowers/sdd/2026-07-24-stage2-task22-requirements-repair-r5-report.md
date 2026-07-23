# Stage 2 Task 22 Requirements Repair R5 Report

Date: 2026-07-24

Status: worker candidate complete; not an acceptance decision

## Scope and authority

- Accepted Task 21 base:
  `53362812364db043273c318fd8ad255eb74cba28`
- Starting Task 22 candidate:
  `be2de9ab2769e7498f1942419eb71c1c79d94b48`
- Starting tree:
  `cc3318911e5025edab8418d1e6f3d326d1e09f24`
- Frozen Task 22 brief SHA-256:
  `353db60c38ce3e6ea7d26f234c3fef402ceb76fc58e3dff730ea09043f00f0f1`
- R5 requirements review SHA-256:
  `4b2a3836b3de24e4e85f5e0c889fc981044a978f90cc903cb678074664a3b94e`
- R5 repair brief SHA-256:
  `723e0fde1dd2c3c9bcbe44051ec1fb250ab0af32f5b547a6bff665a402101d58`
- Accepted-base-to-candidate review package SHA-256:
  `58777e51fb2daa702acbfd615a8bb5d58e56944318dcd23be65cd02c948784cc`

Required reading was completed in the prescribed order: `AGENTS.md`, the
Stage 2 plan and Task 22, the frozen Task 22 brief, Requirements Reviews R4
and R5, the R4 repair brief and report, the R5 repair brief, and the complete
review package. No Task 23 or Stage 3 work was entered.

The protected untracked file
`docs/plans/2026-07-21-hundun-flow-stage7-neural-warm-start.md` was not opened,
modified, staged, or cleaned.

## Navigation and closure sweep

CodeGraphF and exact `rg` sweeps were used to navigate the three facade
implementations, `MaterialDiagonalOperator`, external Jacobi caches, test-only
seams, retry-BE restore sites, tests-off symbol contracts, and affected
Task 18--21 tests. The initial graph reported 232 files, 8,506 nodes, and
26,767 edges. An intermediate `codegraphf sync .` completed successfully
after the original eight-file repair. The final sync after the
coordinator-approved private linear evidence extension also passed and
reported ten changed graph files, nine modified and one removed, covering
733 nodes.

The same-class closure is:

```text
constant -> owned momentum operators + delegated external preconditioners
material -> owned momentum operators + delegated external preconditioners
ideal    -> delegated material facade operators/preconditioners
```

## RED evidence

### Cache authority RED

The direct cache-authority assertions were written first and the focused
target was built with:

```text
cmake --build build/debug --target test_adaptive_time_control_mpi -- -j2
```

The build failed because the old implementation had no
`ConstantDensityPisoTestAccess::facade_cache_snapshot`, and the material and
ideal summaries had no exact `operator_count` or `operators` inventory. Five
compile errors directly proved that the old `source_generation` and
aggregate-summary substitutes could not satisfy the authentic operator
authority oracle.

### Uninterrupted retry-BE RED

The uninterrupted-control requirement was written before its fixture and the
same focused build command failed on the intentionally missing
`require_retry_be_uninterrupted_control` assertion helper. This proved that
the prior fixture had no uninterrupted path to compare with the restored
retry-BE path.

During the minimum-GREEN refinement, the new behavioral oracles also caught
two invalid assumptions:

- the former constant early non-finite terminal did not reach momentum
  operator/cache work, so it was moved to the authentic final-continuity
  terminal;
- active and staging Jacobi buffers swap during a valid refresh, so persistent
  allocation authority is compared as the collision-free sorted two-buffer
  set, not by a false active-slot equality.

The former ideal early density-closure terminal likewise did not reach the
operators and was moved to the authentic final-conservation terminal.

## Implementation

The coordinator cache-contract addendum is implemented as written:

- physical/history/committed/trial/controller/closure state retains exact
  rollback assertions;
- persistent workspace and cache allocation identities and capacities remain
  exact;
- derived operator and preconditioner revisions are monotonic and are not
  rewound;
- every valid Jacobi cache revision is required to equal the exact owned
  momentum operator revision from which it was derived;
- clean successors prove changed diagonals refresh instead of reusing a stale
  inverse;
- ideal gas exposes an exact delegated material authority and explicitly marks
  the public-facing test snapshot as delegated;
- a fresh constant facade/preconditioner cold-path failure immediately after
  `begin_attempt()` proves that no operator revision is published and all
  three Jacobi caches remain invalid and unallocated, while the terminal
  fixtures prewarm authentic caches before testing full retry-limit
  chronology.

The new compile-gated, read-only facade snapshots expose an ordered,
collision-free vector of every previously inventoried persistent workspace
identity and capacity, a fixed array plus exact count for the three owned
momentum operator identities/revisions/diagonal values, and delegation state.
Allocating snapshot functions are not `noexcept`. No public header, public
`Preconditioner` transaction API, revision rewind, or numerical algorithm was
added.

The tests inspect the actual external `JacobiPreconditioner` objects using the
private compile-gated `PreconditionerTestAccess::jacobi_storage()` seam. The
two persistent Buffer identities and real Buffer byte sizes are observed
directly, and the active cached inverse is copied as exact FP64 values. Each
value is compared bitwise with independently evaluated `1 / operator
diagonal`; this proves the revision label and actual cached content are
coherent and that a stale inverse is not reused. Mutation checks prove that
changing a workspace identity/capacity, operator
identity/revision/diagonal value, Jacobi allocation identity/byte
size/revision, or one cached inverse value makes the oracle fail.

A completion code review initially found that the first version synthesized
Jacobi capacities, checked only revision labels rather than inverse contents,
and left `libhundun_linear.a` out of the registered tests-off archive scan.
The coordinator approved a minimal tracked-scope extension limited to
`linear/src/preconditioners_test_access.hpp` and
`linear/src/preconditioners.cpp`. Those files only extend the existing
`HUNDUN_LINEAR_ENABLE_TEST_ACCESS` private read-only snapshot; tests-off
production behavior and the public API are unchanged. The registered
contract now scans the linear archive and its test-seam names as well.

The new cold-path assertion initially expected a controller terminal failure,
but diagnostic instrumentation observed a committed result. Root-cause
tracing showed that the authentic controller entry intentionally resets all
product test faults before advance, clearing the injected constant stage. The
final cold-path oracle therefore calls the authentic constant facade
`attempt()` directly, where the after-`begin_attempt()` injection is active,
and verifies the recoverable result, unchanged operator authority, and three
invalid/unallocated Jacobi caches.

At the committed retry-BE boundary, the original controller, state, facade,
and preconditioners are retained for uninterrupted advance. A separate
`FlowState` is created from exact history, committed, trial, and metadata
values; a fresh controller is restored from the sealed exported state; and
fresh facade and Jacobi objects are used. Both authentic
controller-to-facade stencils are observed and compared bitwise for attempted
`dt`, previous `dt`, and all coefficients. Order, proposed-next-`dt`,
accepted order, history readiness, retry count, and post-step metadata are
also compared, while the independent mathematical stencil oracle remains in
place and does not call the product stencil factory.

## Changed tracked files

```text
CMakeLists.txt
flow/src/constant_density_piso.cpp
flow/src/constant_density_piso_test_access.hpp
flow/src/ideal_gas_closure_test_access.hpp
flow/src/ideal_gas_piso.cpp
flow/src/material_density_piso.cpp
flow/src/material_density_piso_test_access.hpp
linear/src/preconditioners.cpp
linear/src/preconditioners_test_access.hpp
tests/cmake/task22_test_access_contract.cmake
tests/mpi/test_adaptive_time_control.cpp
```

No other tracked file was changed. The two linear files are the
coordinator-approved allowlist extension; `CMakeLists.txt` was already in the
R5 brief's allowed subset.

## GREEN matrix

All commands used `/tmp/hundun-toolchain/cmake/bin/ctest` rather than the old
system CTest.

### Focused implementation and Debug

```text
cmake --build build/debug --target test_adaptive_time_control_mpi -- -j2
PASS

mpiexec -n 1 build/debug/test_adaptive_time_control_mpi acceptance
PASS
mpiexec -n 2 build/debug/test_adaptive_time_control_mpi acceptance
PASS
mpiexec -n 4 build/debug/test_adaptive_time_control_mpi acceptance
PASS

cmake --build build/debug -- -j2
PASS

/tmp/hundun-toolchain/cmake/bin/ctest --test-dir build/debug \
  -L task22 --output-on-failure
PASS 22/22
```

### Affected Task 18--21 Debug regressions

```text
/tmp/hundun-toolchain/cmake/bin/ctest --test-dir build/debug \
  --output-on-failure \
  -R '^(test_fixed_step_piso_[124]_rank|test_material_density_piso_[124]_rank|test_ideal_gas_piso_[124]_rank)$'
PASS 9/9
```

This covers constant, material-density, and ideal-gas fixed-flow tests at MPI
1/2/4 ranks.

### Release and sanitizers

```text
cmake --build build/release -- -j2
PASS
/tmp/hundun-toolchain/cmake/bin/ctest --test-dir build/release \
  -L task22 --output-on-failure
PASS 22/22

cmake --build build/asan -- -j2
PASS
ASAN_OPTIONS=detect_leaks=0 \
  /tmp/hundun-toolchain/cmake/bin/ctest --test-dir build/asan \
  -L task22 --output-on-failure
PASS 22/22

cmake --build build/ubsan -- -j2
PASS
/tmp/hundun-toolchain/cmake/bin/ctest --test-dir build/ubsan \
  -L task22 --output-on-failure
PASS 22/22
```

### Fresh tests-off Release boundary

A new independent build directory was used:

```text
/tmp/hundun-toolchain/cmake/bin/cmake -S . \
  -B /tmp/hundun-flow-task22-tests-off-r5-final \
  -DCMAKE_BUILD_TYPE=Release \
  -DHUNDUN_BUILD_TESTS=OFF \
  -DCMAKE_CXX_COMPILER=/tmp/hundun-toolchain/clang/bin/clang++ \
  '-DCMAKE_CXX_FLAGS=-stdlib=libc++ -DOMPI_SKIP_MPICXX'
PASS

/tmp/hundun-toolchain/cmake/bin/cmake \
  --build /tmp/hundun-flow-task22-tests-off-r5-final -- -j2
PASS full build
```

The Task22 standalone preprocessing/public-header syntax consumer was then
run manually against the fresh flow and diagnostics archives and passed.

`nm -C` was scanned across the fresh `libhundun_flow.a`,
`libhundun_linear.a`, and `libhundun_material_diagnostics.a`. It found none
of the test-only access classes, Jacobi storage seams, diagnostic mutation
seams, attempt-observer seams, `facade_cache_snapshot`,
`delegated_material_cache_snapshot`, or
`material_facade_cache_values_for_ideal`.

`ldd` with
`LD_LIBRARY_PATH=/tmp/hundun-toolchain/clang/lib/x86_64-unknown-linux-gnu`
resolved every dependency of the fresh tests-off `hundun` executable,
including the intended toolchain `libc++.so.1` and `libc++abi.so.1`.

Fresh artifact SHA-256 values:

```text
740240821bddd45ea4e16a29705d276193660107e5f2285dcafb36e6ec3cdea3  libhundun_flow.a
6ebc2230c99dbff75a5478fc44fc9d32531b89be7d82b426ef1ad6cda9e34ea3  libhundun_linear.a
b9fb93141b7eba1c360d08fedc3c532ef7d86f4d67d1c8c0a348cc97ae6ca63c  libhundun_material_diagnostics.a
cbc4f74dbe3e14e80cfdcf30cecd5418ff3a14aa9625e3d6acf2af0e03fdaa7c  hundun
```

### Global provenance

Global provenance remains a separate label from Task22:

```text
/tmp/hundun-toolchain/cmake/bin/ctest --test-dir build/debug \
  -L provenance --output-on-failure
PASS 39/39
```

`git diff --check` also passed before report generation.

## R4 report correction

This report explicitly supersedes the R4 report's repair-brief SHA typo:

```text
incorrect 3c5056ad...
correct   3c5056d28b1b9485ef0b664f85826e1e62f6c6a000e433f43c2a9af103664611
```

The correct value was independently reproduced from
`.superpowers/sdd/2026-07-24-stage2-task22-requirements-repair-r4-brief.md`.

## Remaining risk and ownership

No production numerical bug was exposed by the direct REDs; this repair is
tests-on evidence/introspection plus stronger authentic-path tests. The
independent completion review's three Important evidence findings and one
Minor full-state comparison suggestion were all addressed: real cached
content, real Buffer byte sizes, registered linear archive scanning, and
exact post-success `TimeControlState` equality are now asserted. A focused
strictly read-only re-review marked all four findings closed, confirmed the
cold-path oracle, found no new Critical or Important issue across the eleven
tracked files, and assessed the worker patch ready.

This worker matrix does not claim Task 22 acceptance. The main agent owns the
final exact-HEAD complete Debug acceptance and any subsequent independent
review. No push or publication was performed.

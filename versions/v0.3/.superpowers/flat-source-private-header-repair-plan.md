# HUNDUN-FLOW product-private header boundary repair plan

## Objective

Make the tests-free product projection independently buildable while preserving
the accepted product's public API and runtime behavior.

## Global constraints

- Public headers remain only under `include/hundun/`.
- Product-private headers remain flat under `src/` and use a registered prefix
  plus `_detail.hpp`.
- Test-only access, mutation, observation, and fault injection remain under
  `tests/support/`.
- No product source may include a path under `tests/`.
- Product targets may use `${PROJECT_SOURCE_DIR}/include` publicly and
  `${PROJECT_SOURCE_DIR}/src` privately; they must not require the repository
  root or `tests/` as an include directory.
- Do not change algorithms, constants, thresholds, PISO corrector count, public
  names, configuration keys, restart/diagnostic schemas, units, signs, MPI
  collective order, rollback behavior, or default runtime paths.
- Do not commit, sign off, push, publish, or switch repository directories.
- Verification is limited to layout/product-projection fixtures, builds,
  headers/linkage, and affected low-cost tests. No large numerical jobs.

## Task 1: Strengthen RED boundary fixtures

- Extend `tests/cmake/source_layout_fixture.cmake` so any `src/` include of
  `tests/` fails.
- Add a low-cost product-projection fixture that copies only the approved
  product roots to a temporary directory and verifies default configure plus a
  focused product build without `tests/`.
- Demonstrate the fixture fails on the current candidate for the known missing
  governance headers.

## Task 2: Reclassify the eight hard dependencies

Inspect and repair only these source/header pairs:

1. `src/app_case_config_broadcast.cpp`
2. `src/bc_basic_boundary.cpp`
3. `src/exec_execution.cpp`
4. `src/fvm_cell_centered.cpp`
5. `src/fvm_immersed_reconstruction.cpp`
6. `src/rt_field_storage.cpp`
7. `src/rt_halo_exchange.cpp`
8. `src/rt_mpi_context.cpp`

Move always-required implementation declarations into suitable existing or new
`src/*_detail.hpp` headers. Keep test-only APIs in the corresponding
`tests/support/` headers. Product sources must include only public or product-
private headers.

## Task 3: Make the source-boundary RED mutation-sensitive

- Reject any `tests/` path text in product C++ files, including a path hidden
  behind a macro and later consumed by `#include MACRO`.
- Add a self-contained low-cost mutation fixture proving that a macro-hidden
  `tests/support` include is rejected.
- Keep the real-tree fixture RED while known product/test dependencies remain;
  do not weaken or special-case the rule to obtain GREEN.

## Task 4: Reclassify diagnostics dependencies

- Execute this task as two sequential acceptance clusters: cluster A covers
  checkpoint-v2, structured, and material-density-PISO diagnostics; cluster B
  covers ideal-gas closure, material-density transport, and time-control
  diagnostics. Use a fresh implementation worker and independent review for
  each cluster.
- Remove all `src/diag_*.cpp -> tests/support` dependencies.
- Keep typed diagnostic test access under `tests/support`; use inline friend
  adapters or primitive/product-type raw hooks where implementation state is
  otherwise inaccessible.
- Verify clean tests-on consumers, clean tests-off diagnostic targets, and
  affected low-cost diagnostics tests only.

## Task 5: Reclassify flow dependencies

- Execute as three sequential acceptance clusters: A covers flow state and
  checkpoint-v2; B jointly covers adaptive time control plus the
  constant/material/ideal-gas PISO, density transport/closure, and momentum
  hooks it schedules; C covers the immersed flow authority adapters. Use a
  fresh implementation worker and independent review for each cluster.
- Remove all `src/flow_*.cpp -> tests/support` dependencies without changing
  time integration, pressure-velocity coupling, checkpoint, rollback,
  authority, diagnostics, or numerical behavior.
- Preserve the existing test-facing APIs and mutation-sensitive coverage in
  `tests/support`.
- Verify clean tests-on consumers, clean tests-off flow/application targets,
  and affected low-cost unit or small-rank tests only.

## Task 6: Reclassify remaining linear/FVM dependencies

- Remove the remaining `src/lin_*.cpp` and `src/fvm_*.cpp` dependencies on
  `tests/support`.
- Preserve preconditioner and immersed-operator test semantics, signs,
  ownership, MPI failure convergence, and rollback.
- Verify clean tests-on consumers, clean tests-off targets, and focused
  low-cost tests only.

## Task 7: Tighten CMake include authority

- Remove repository-root include directories from product targets.
- Keep `${PROJECT_SOURCE_DIR}/src` private and
  `${PROJECT_SOURCE_DIR}/include` public.
- Governance targets may add repository-root and test-support paths inside
  `tests/CMakeLists.txt` only.
- Make the product-projection fixture GREEN and rerun affected tests-off,
  header/linkage, rollback, collective, and test-access checks at low cost.

## Task 8: Audit and evidence refresh

- Verify no `src/` file includes `tests/`.
- Verify the tests-free projection configures and builds by default.
- Re-run normalized accepted-blob comparison and review every non-include
  product-source change.
- Regenerate the SHA-256 migration manifest with an explicit, truthful
  classification for the authorized private-boundary extraction.
- Update the audit report. Stop before commit or directory switching.

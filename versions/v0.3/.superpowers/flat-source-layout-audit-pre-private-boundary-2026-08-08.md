# HUNDUN-FLOW flat source-layout candidate audit

## Decision

Status: `PAUSED_FOR_APPROVAL`

The governance layout candidate is complete and its mixed-repository build is
verified. It is not yet safe to create the tests-free product projection. The
accepted product sources contain eight unconditional dependencies on governance
test-access or test-seam headers. Resolving that boundary requires source edits
outside the migration manifest's permitted change kinds, so the stop rule is
active.

No commit, DCO sign-off, repository-name switch, publication, or push was
performed.

## Frozen baseline

- acceptance profile: `task11-core-functional-v1`
- Task 11 result: `CORE_ACCEPT`
- accepted product HEAD: `66080e324089599711fdb26082af9b330bfdb5ce`
- accepted product tree: `ab071a61f00eba9ec973beb0fe600066a33ef74f`
- migration worktree: `/home/wyf/code_dev/.worktrees/hundun-flow-flat-layout`
- branch: `governance/flat-source-layout`
- current HEAD remains the accepted product HEAD
- governance layout commit: `PENDING_APPROVAL`
- product initial commit: `PENDING_APPROVAL`

## Completed candidate work

- Moved all product public headers into the flat `include/hundun/` directory.
- Moved product implementations and private headers into the flat `src/`
  directory.
- Applied the controlled file-prefix registry without changing C++ namespaces,
  public type names, function names, configuration keys, restart fields, or
  diagnostics schemas.
- Moved test-access and test-seam headers into `tests/support/` and marked them
  `not_exported`.
- Reduced the root CMake file to project-wide setup and subdirectory entry
  points.
- Moved product target registration into `src/CMakeLists.txt` and governance
  registration into `tests/CMakeLists.txt`.
- Set the product default to `HUNDUN_BUILD_TESTS=OFF`.
- Named the shortened unambiguous product targets `hundun_fvm` and
  `hundun_sdk`.
- Added the low-cost layout fixture
  `tests/cmake/source_layout_fixture.cmake`.
- Added the Chinese product manual
  `docs/development/naming-and-style.md`.

## Semantic diff audit

The 208 moved source/header/build-helper files were compared with their accepted
blobs:

- 28 files are byte-identical moves.
- 180 files differ only in preprocessor include lines.
- After removing include lines from both sides, all 208 pairs are byte-identical.
- Unexpected product content changes: 0.
- Old module-form includes in `include/`, `src/`, and `tests/`: 0.
- `git diff --check`: pass.

The mapping manifest contains 372 rows and uses only the approved change kinds:

- `cmake_registration_only`: 3
- `documentation_only`: 1
- `identical_move`: 20
- `include_path_only`: 162
- `not_exported`: 186

Manifest:
`.superpowers/flat-source-layout-manifest.tsv`

Manifest SHA-256:
`5399081939535217efe15a933a9a16fe03edb49c75472eadfa10309a21cc8a23`

Every old hash was checked against the accepted HEAD and every new hash was
checked against the current candidate file.

## Documentation gate

- technical draft SHA-256:
  `63d90e5f9e61761edabe5ae3aa82963f17afb8a6a24e94b3874fac28b222f0d1`
- after `humanizer-zh`:
  `52dd3b12280f8b4ed9cd7048c62e2606e415f9b9fe0def453deaa50c78dbc87e`
- after `shuorenhua`:
  `6f3813baf38e705d6a0b4dc6d05cec8855c02d2b724e538f0af63e0159c923f5`
- after main-agent technical correction:
  `c61d509e19061f458189b9fdfd3dbfdfcb35b48539812c51c38917630ff3cd63`
- prohibited history/process/private-source names in the product manual: 0

## Verification evidence

All builds used Clang 15.0.6, libc++, Open MPI 3.1, C++17, and at most two build
jobs.

1. Layout fixture: pass.
2. Debug tests-off configure and complete product build: pass.
3. Debug tests-on complete build of all registered targets: pass.
4. Focused Debug runtime set: 8/8 pass, including layout, public headers,
   plugin loading, one-rank transaction/rollback, signed-force authority, and
   collective status.
5. Non-large, non-acceptance unit set: 83 tests passed on the first run; eight
   governance path-authority failures were corrected and all eight passed on
   focused rerun. The additional fixture setup also passed.
6. Default Release configure with tests disabled and focused `hundun` plus
   `hundun_sdk` build: pass.
7. No 48-cube, 96-cube, warped, prism, convergence matrix, large MPI, or
   sanitizer numerical job was started.

The first focused CTest invocation could not load `libc++.so.1`; all seven
affected programs failed before `main()`. Repeating the same commands with the
project's documented Clang runtime library path produced the 8/8 pass above.

Accepted-tree warnings in `ib_ghost_stencil_plan.cpp` and
`fvm_immersed_operator.cpp` remain unchanged and were not repaired as part of
this layout-only candidate.

## Product-projection blocker

A clean temporary product projection was created without `tests/`. Its default
Release configure passed, but compilation stopped at missing governance headers:

- `src/exec_execution.cpp` requires
  `tests/support/exec_execution_test_access.hpp`.
- `src/rt_mpi_context.cpp` requires
  `tests/support/rt_mpi_context_test_seam.hpp`.

The diagnostic projection is located at:
`/tmp/hundun-flat-product-projection.fcyliJ`.

Across `src/`, 26 implementation files contain 36 governance-header include
lines. Most are already protected by existing test-access macros. The following
eight includes are unconditional and therefore form the hard product-export
boundary:

1. `src/app_case_config_broadcast.cpp`
2. `src/bc_basic_boundary.cpp`
3. `src/exec_execution.cpp`
4. `src/fvm_cell_centered.cpp`
5. `src/fvm_immersed_reconstruction.cpp`
6. `src/rt_field_storage.cpp`
7. `src/rt_halo_exchange.cpp`
8. `src/rt_mpi_context.cpp`

Shipping `tests/support/`, retaining tests in the product repository, or
generating hidden copies during CMake configuration would each violate the
approved product/governance boundary.

## Recommended authorized repair

Extract the small always-required internal declarations from those eight
governance headers into product-private, correctly prefixed `src/*_detail.hpp`
files. Keep fault injection, mutation, observation, and test-access APIs in
`tests/support/`, gated by the existing `HUNDUN_*_ENABLE_TEST_ACCESS` macros.

The repair must preserve every existing implementation body and default runtime
path. It must not alter numerical constants, algorithms, thresholds, PISO count,
units, signs, MPI collective ordering, rollback behavior, schemas, or public
interfaces.

Proposed file whitelist for the repair:

- the eight product `.cpp` files listed above;
- their eight corresponding headers in `tests/support/`;
- existing matching `src/*_detail.hpp` files where available;
- new, prefix-compliant `src/*_detail.hpp` files only where no suitable private
  header exists;
- `src/CMakeLists.txt` only if a private header must be listed explicitly;
- `tests/cmake/source_layout_fixture.cmake` and one product-projection build
  fixture.

This repair cannot truthfully be classified as only `include_path_only` without
first inspecting and splitting the declarations. Authorization is therefore
required before editing those product files.

## Approval boundary

If the repair above is approved, the next action is the minimum private-boundary
extraction followed by:

1. the failing tests-free product projection build turning green;
2. normalized product diff audit;
3. Debug tests-off and focused tests-on builds;
4. header/linkage and the eight affected low-cost governance checks;
5. manifest regeneration.

The governance migration commit and final directory-name switch remain separate
approval points.

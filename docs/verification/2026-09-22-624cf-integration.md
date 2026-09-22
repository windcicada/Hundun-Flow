# 624CF selective integration

## Scope and identity

Request: review `624CF` and integrate its valid changes, without reverting newer
mainline work. This is a selective port, not a merge of the entire donor history.

- Donor: `624CF@9adcdd804f20dd920d403cad4298224eb29ec6b1` (19 changed paths).
- Donor/main merge base: `8235e59d5a841cfa9c11921681167ea7cfc9689e`.
- Initial remote main: `98486e9c3e529bfc1fdb2bc39448215cb8149b76`.
- Updated integration base: `12e099fab3dfccebecc4342743e1c51ca6f06d4b`.
- Isolated worktree: `/home/wyf/code_dev/merge624`, branch
  `codex/merge-624cf-20260922`.
- Original `/home/wyf/code_dev/flow` worktree retained its modified `docs/alg.md`
  and untracked `docs/{agent,next,review}.md`; none belongs to this change.

The donor is owner-authored HUNDUN code with a DCO trailer. Incoming mainline
authors, committers and trailers were checked. New commits use
`WANG YUDONG <wangyudong@buaa.edu.cn>` with DCO. No COAST or OpenFOAM code was copied.

## Retained changes

| Area | Selected behavior and verification |
| --- | --- |
| Case assets | STL fingerprinting uses the existing mesh parser's 1 GiB budget; generic referenced data retains its 64 MiB budget. Sparse-file tests cover both limits. |
| ESF restart | Do not apply fresh-start stochastic offsets to temporary restored means. A two-species test restores pure B with configured offsets of ±0.1 and checks current/previous physical fields exactly. |
| Spray retry | Event-capacity exhaustion enters the existing bounded wave-retry path instead of being treated as permanently invalid input. Existing product spray tests cover 1/2/4 ranks. |
| Spray substeps | Grow an accepted substep by at most two, capped by the original bound, instead of restarting at the largest step after every acceptance. The analytic trajectory test also bounds attempt count. |
| Geometry tools | Add the identity-pinned STL repair utility and standalone 624CF geometry audit target. These neither import fields/parcels nor advance reacting flow. |

The all-solid MPI partition fix is **already supplied by mainline `3c82617`**,
including both projection routes and globally empty-domain rejection. Its donor
implementation and duplicate test were dropped in favor of that newer version.

The app-case test previously rejected `NOT_RUNTIME_AXES`, although the existing
parser deliberately accepts an alphabetic vendor prefix followed by
`_RUNTIME_AXES`. The negative fixture now uses `RUNTIME_AXES_INVALID`; positive
fixtures cover `COAST_RUNTIME_AXES` and `Vendor_RUNTIME_AXES`. No parser behavior
was changed for this correction.

## Standards

Source positions below refer to donor `9adcdd8`, unless otherwise noted.

1. **P1 — implicit TCR history initialization**, `tools/v04_pdf_import.cpp:739`:
   the old importer bypasses the explicit `--model-history initialize` policy
   documented in `docs/gas-work.md:866`. Excluded.
2. **P1 — dangling borrowed string**, importer lines 225–274: `file` is a
   `string_view` into a yyjson document freed before `root / file` is evaluated.
   Excluded with the old importer.
3. **P1 — unchecked parcel payload identity**, old importer/mapper: the mapper
   emits `mapped_sha256`, but the importer does not verify it. This does not
   satisfy the source/target identity contract in `docs/rc.md:3676`. Deferred.
4. **P2 — unbounded input staging**, old importer: row count is checked after
   accumulating all rows; the line limit is checked only after `getline` has
   allocated the line. Deferred with the old importer.
5. **P1 — incomplete method-recovery persistence**, candidate flux-history
   port: keeping distinct physical flux revisions makes immediate checkpoint
   publication structurally possible, but the checkpoint does not persist the
   pending first-step BE recovery flag. `TimeControllerState::restart` sets it;
   `restart_exact` clears it (`bc_time.cpp:360,388`). A write/read cycle can
   therefore select BDF2 using reconstructed history. The entire donor
   `complete_source_history` publication change is excluded. This finding is
   based on the controller/snapshot code path; a BDF2 round-trip advance was not
   run in this task.

Additional observations, not counted as hard violations:

- **Possible duplicated code:** the donor reinstates a monolithic importer,
  whereas mainline shares `app_import.cpp` between entry points
  (`docs/gas-work.md:779`). Preserve mainline's implementation boundary.
- **Dependency admission hold:** the new Cantera hash/manifest has no matching
  artifact or rebuild evidence available here. Do not add its allowlist entry.
  Verification uses the already admitted library SHA-256
  `093b62eadc4d44c3ef227c2d59554542820fdd8fde3497a0dcc46e3360040760`.

Standards total: five actionable findings, one design heuristic and one evidence
hold. Worst severity within this axis: P1.

## Spec

The applicable specifications are the tracked conservative continuation rules
in `docs/alg.md` and the migration identity/publication contract in `docs/rc.md`.

1. **P1 — nonconservative field transfer**, donor
   `tools/v04_aecsc_pdf_transfer.cpp:284`: nearest-value fill followed by a global
   density-only correction does not conserve momentum, species or enthalpy.
   Its own self-test produces x-momentum `88/3` instead of `26` (12.82% error),
   while merely printing the expected `26`. Exclude the transfer executable,
   its CMake wiring and its test; `docs/alg.md:27` requires conservative mapping.
2. **P1 — unsafe import path**, the dangling yyjson-backed view above makes the
   requested reliable continuation import unsafe. Exclude the old importer.
3. **P2 — incomplete conservation-receipt binding**, an altered finite parcel
   payload can be accepted without checking its recorded digest. Defer the
   spray import, remap changes and associated tests until the receipt is bound
   to the actual bytes consumed by the current shared importer.

Spec total: three findings; worst severity within this axis: P1.

## Verification

Build: Clang 15.0.6, Release, C++17, `-ffp-contract=off`, OpenMPI (MPI 3.1), and
the existing admitted Cantera 3.2.0 package. The existing Jammy rootfs/toolchain
was mounted read-only using a local `check/run` bubblewrap launcher; only this
isolated worktree was writable. No package installation or dependency admission
change was made.

Commands inside that environment:

```sh
cmake -S . -B check/b -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_MAKE_PROGRAM=/opt/ninja \
  -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
  -DCMAKE_C_FLAGS=-ffp-contract=off \
  '-DCMAKE_CXX_FLAGS=-stdlib=libstdc++ -D_GLIBCXX_USE_CXX11_ABI=1 -ffp-contract=off' \
  -DHUNDUN_BUILD_TESTS=ON -DHUNDUN_ENABLE_REACTING_CANTERA=ON \
  -DHUNDUN_CANTERA_PACKAGE_ROOT=/opt/ct
cmake --build check/b -j 4 --target v04_app_case_test \
  v04_solver_initialization_projection_mpi_test \
  v04_product_spray_advance_test v04_624cf_geometry_audit \
  v04_product_restart_storage_mpi_test v04_models_spray_events_test \
  v04_product_pressure_energy_temporal_convergence_test
ctest --test-dir check/b --output-on-failure --timeout 90 -j 2 \
  -R '^(v04_app_case|v04_models_spray_events|v04_product_restart_method_history|v04_solver_initialization_projection_mpi_[124]|v04_product_spray_advance_mpi_[124]|v04_product_restart_storage_mpi_[124])$'
git diff --check
```

Result: **12/12 focused tests passed** after updating to `12e099f` and excluding
the incomplete method-recovery port. Tests compare all restored physical ESF
fields exactly; the derived conductivity/cp transport cache is rebuilt and is
not asserted bit-identical (a one-ULP change was observed). Existing full-field
comparisons in the other restart cases remain unchanged.

An exploratory regression also demonstrated that the donor flux-history patch
alone mixes synthetic current revision 17 with restored previous revision 81,
failing with status `2/913`. A companion revision correction was evaluated, but
both it and that exploratory test were removed after the persistence problem
above was found. No incomplete fix is included.

Geometry commands (first uses host Python; audits use the build environment):

```sh
python3 tools/v04_repair_624cf_stl.py \
  --source /home/wyf/code_dev/flow/check/ji/case/Geometry/cf2222.STL \
  --output /home/wyf/code_dev/merge624/check/geometry/cf1.stl \
  --report /home/wyf/code_dev/merge624/check/geometry/repair.json
timeout 90s mpiexec -n 1 check/b/versions/v0.4/v04_624cf_geometry_audit check/geometry/cf1.stl
timeout 90s mpiexec -n 4 check/b/versions/v0.4/v04_624cf_geometry_audit check/geometry/cf1.stl
```

- Source SHA-256: `307879356b66179c1129dedd6f25213ecf240b5093d994cd0e9d361b505df65b`.
- Repaired SHA-256: `c1d3c7dfaa402931dae23e73677b7d4f42206a6977d02f6f308b841fc280b10f`,
  identical to the existing frozen repaired asset; neither original was changed.
- Repair: 280968 triangles, 421452 edges, six duplicate triangles removed, four
  seam pairs, relative enclosed-volume change `5.550986309287901e-7`.
- Both audits pass: 325846 fluid cells, 87767 immersed links, 330 inlet cells,
  1031 outlet cells. Surface volume `0.0020246215943216885 m^3`.
- Existing-output and incorrect-source-hash rejection checks pass without
  changing existing assets or producing a new output/report on those failures.

No full CTest, long CFD continuation, full spray/flame case, COAST comparison,
or performance run was performed. This result does not accept 624CF reacting
continuation, complete field/parcel migration, or full COAST replacement.

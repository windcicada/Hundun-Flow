# GTMC cold restart port — change and verification ledger

Status: active; user authorized non-equivalent VTK mean-field state transfer
on 2026-09-08 and explicitly scoped this work as Hundun development testing,
not a physical validation or exact continuation of COAST's PDF solution;
**the requested GTMC short and medium runs have not passed yet**. Unit tests and
geometry probes below are not acceptance runs.

Actual candidate `f10f15a077961a1be01027f613e317e7759eaf50` was built clean
with GCC 11.4 Release/tests OFF and run on 64 local ranks. The 100-step request
failed before the first new commit; a separate one-step replay reproduced the
same failure. The user subsequently authorized implementing the original
COAST zero-gradient/mass-closure outlet semantics for development testing.
That historical candidate was superseded: clean 0a9ed82 accepted one new GTMC
step before the G22 failure. Clean G23 `fcfc524` crossed that blocker and
committed steps 20001–20066; its requested short100 was deliberately stopped
after G25 established a stale-material defect, so it is NOT a short pass.
G26 below is the current repair candidate; short100 and medium3000 must use
a newly signed/clean build and must not concatenate candidates' step counts.
Detailed receipts are outside the source tree:
`/home/administrator/gtmc-hundun-port-20260908/candidate-f10f15a/RUN_STATUS.md`.

Baseline: `86542bb96678ec844ae5ac95d8f6391993da239e`, Hundun-Flow v0.4.
All code changes are local and unpushed. This ledger is updated with each
integration stage so mainline synchronization does not rely on chat history.

Current integration update: the temporary ProductCompiler rejection 15901 was
removed after source-aware frozen enthalpy, pressure-energy Schur and predictor
regressions passed. Earlier table entries describing that guard are historical
checkpoints, not the current executable behavior. Full native bridge and actual
64-rank GTMC advancement are the next gates. Before these runs the source will
be committed with DCO sign-off and rebuilt with honest clean identity; run
receipts will be stored outside the source tree until the candidate is stopped.

## Fixed input authority

The target is remote `/home/rsfz/data8t/gtmc/dev/d5_cold_20000`, step 20000,
time 0.015465333126485348 s, frozen 153 × 247 × 151 cells. The uint8 marker
contains 3,901,529 fluid and 1,804,912 solid cells; SHA-256:
`e5f7e9c64d88b2ac86ecaa255fb9a45bec775fca734f0fb07b5320adffb613b4`.
Twenty air patches and 250 internal fuel faces must retain their membership.
Whole Cartesian faces are authoritative; this is not a cut-cell-area model.

## Change inventory

### Outlet/patch development candidate follow-up

- G26 (development fix under verification; based on clean `fcfc524`):
  exchange generic-mixture thermal material on MPI/periodic neighbor slabs.
  G25 diagnostic root cause: candidate conductivity ghosts stayed at initial
  air values while received primitive ghosts and live conductivity changed.
  On rank 22 / global (77,82,59), replacing ONLY the x-minus conductivity in
  the production harmonic face formula predicts -1.02828700341196e-6 W;
  observed candidate-minus-terminal energy residual is -1.0282845193999855e-6
  W (2.484e-12 W remainder). The original first two attempts at step 20008
  exhaust 12 refinements with energy about 1.38e-10 despite continuity below
  16 eps. This is not a reason to relax either gate or raise iteration caps.
  Locations/reasons:
  `src/core_product_freeze.cpp:182` separates MPI exchange from COAST-only
  physical zero-gradient material handling, preserving generic physical EOS
  closure and G23 inlet-face material;
  `:8194`, `:9469`, `:11656`, `:14854` refresh live/candidate/final-rate k and
  k/cp for all transport kernels (COAST effective SGS conversion is unchanged);
  `:5570` explicitly retains the existing COAST cold-start behavior.
  `tests/mpi/product_pressure_energy_retry_mpi_test.cpp:1164` adds the real
  8^3 ProductDriver warm-restart constant/Sutherland comparison;
  `tests/CMakeLists.txt:1015` registers MPI 1/2/4 regression invocations.
  `src/core_product_freeze.cpp:68` bumps method-history with
  `generic-thermal-neighbor-material-v1`; the migration note in
  `docs/development/public-api-compatibility.md` explains why old V3 rates
  cannot be silently reused. This docs path is repository-relative.
  BEFORE fix, the new test fails Sutherland on all three decompositions
  (6/5792 globalization), while constant transport passes. AFTER the isolated
  halo repair, all MPI 1/2/4 cases pass in 1.35 s and the existing retry suite
  passes 3/3 in 7.57 s. Wider scalar/boundary/outlet checks pass 32/36; the
  four failures are the previously reproduced variable-thermo MPI 1/2/4 and
  stale-authority MPI 2 failures, not newly introduced regressions.
  The original GTMC native-20007 one-step replay then accepts step 20008 with
  BDF2, original dt=6.291342081924722e-9, no retry, C=1.942075401417315e-15,
  E=1.8219410364787955e-15, solver/launcher exit 0. Neither tolerances nor
  refinement capacity changed. Replay took 173.38 s max-rank step time in a
  distinct Release/no-LTO diagnostic build; it is not a performance comparison
  or short/medium acceptance. All 42 selected checks: 38 pass, 4 known fails.
  After the production method-history bump, the new seam and existing retry
  MPI 1/2/4 suites were rebuilt and repeated: 6/6 pass in 5.13 s.
  All source paths above are relative to `versions/v0.4/`. The independent
  worktree `thermal-halo-fix` leaves the active clean `fcfc524` short100 source
  and executable unchanged. External receipts:
  `/home/administrator/gtmc-hundun-port-20260908/pressure-refinement-diagnostic-v2/`.
  Diagnostic replay retained the old method-history identity only for
  controlled A/B; production now has the method-history bump and requires
  new native import before any acceptance run.

- G24/G25 (read-only diagnostic receipts, no numerical edits): G24 separates
  the independent total-energy balance defect (about -392 W) from the much
  smaller terminal equation defect; this audit remains unresolved and is not
  physical validation. G25 disproves Aitken/slow continuity as the primary
  retry blocker and temperature inversion as the energy-floor cause, then
  establishes the stale thermal ghost root cause described in G26.

- G23 (implementation under verification; user explicitly approved the
  boundary-contract change after G22): distinguish discrete Dirichlet p/h/Y
  mirrors from physical inlet-face thermodynamics. No composition clipping or
  EOS/tolerance relaxation. Locations and reasons:
  - `include/hundun/v04_boundary.hpp`, `src/bc_thermo.cpp`: append explicit
    `physical_inlet_face` closure kind, leaving the legacy entry default intact.
    Eligible locally owned fixed-composition inlets use the compiled Y target,
    resolved h mirror/owner face target, and extrapolated pressure. Validate
    physical owner and face by strict EOS, preserve p/h/Y mirrors, publish face
    material in the boundary slots and 2*T_face-T_owner for the T stencil.
    The shared inlet eligibility excludes walls, outlets and periodic/MPI faces.
    Corners extend one deterministic inlet trace with clamped tangential indices.
    Certificate kind, face-source/owner values and semantics bind the new
    dependencies; corrupt mirrors and same-revision owner changes fail closed.
  - Same BC files, `refresh_inlet_material`: generic later halo/zero-gradient
    and turbulence publication would otherwise overwrite face coefficients.
    Atomically reconstruct fixed-T/Y inlet material from owner pressure after
    those operations. Preserve the owner SGS dynamic-viscosity addition
    separately from face molecular viscosity. Optional output views restrict
    writes, avoid rewriting certified density, and retain physical EOS guards.
  - `include/hundun/v04_execution.hpp`, `include/hundun/v04_flow.hpp`,
    `src/solver_cartesian.cpp`, `src/solver_equations.cpp`: explicit opt-in
    equation/kernel material representation, bound by semantic fingerprints.
    Legacy standalone kernels default to exterior-cell coefficient samples.
  - `src/solver_diffusion.cpp`, `src/solver_equation_detail.hpp`: both the
    conservative residual and implicit transmissibility use the same inlet
    face coefficient over the mirror centre-to-centre distance. Primitive
    gradients and interior/other-face harmonic interpolation are unchanged.
    `src/solver_cartesian_detail.hpp`, `src/solver_momentum.cpp`,
    `src/core_conservation_detail.hpp`: inlet cross-stress and development
    viscous-work diagnostics consume the same face viscosity.
    `src/solver_face_flux.cpp`: inlet reconstruction uses rho_face*U_face*A,
    not interpolation of products involving mixed face/cell density samples.
  - `src/core_product_freeze.cpp`: activate both contracts together at fresh,
    live C1/C2 and candidate closures; propagate the kind in every consumer
    binding. Refresh inlet materials after cold restart, momentum publication,
    turbulence, candidate and final-rate derived-field boundary fills.
    Bump `method_history_signature` with physical-inlet-face-thermophysics-v1;
    do not reinterpret previous V3 histories under the changed boundary
    operator. Full acceptance needs a freshly imported V1 primitive checkpoint.
  - `src/solver_piso.cpp`, `src/solver_candidate_boundary.cpp`: reject a face
    material certificate consumed under an exterior-cell kernel (or vice versa).
  - `tests/unit/bc_thermo_test.cpp`: real scalar+enthalpy boundary fill followed
    by closure failed before the fix in 0.35 s; after the first repair serial
    and MPI-2 passed. Added owner/kind tamper, physical-pressure atomicity,
    separate SGS/molecular refresh, and actual conservative diffusion vs
    implicit-transmissibility oracles. Further integrated tests are in progress.
  Existing outlet MPI-1/2/4 plus BC tests passed after the initial integration
  (five tests, 1.94 s); later consumer/refresh revisions are not yet accepted.
  Actual GTMC short-100 and medium-3000 remain unpassed. Do not use a dirty
  build as acceptance or confuse these focused tests with the GTMC gates.
  G23 integrated follow-up: 25 wider regressions initially failed during
  initialization with 736. A launched GDB run of the real open scalar fixture
  showed all relevant ghost widths were 2 (falsifying the width hypothesis),
  but the new material refresh received pressure_reference=0 while its caller
  `rebuild_cold_velocity_dependents(cold_pressure_reference=101325)` had the
  correct state. The two new cold-path calls had accidentally used the not-yet
  initialized runtime member instead of the supplied cold pressure reference.
  Both now pass `cold_pressure_reference`; strict positive-pressure guards
  and view-width checks remain unchanged. Evidence: external
  `species-diagnostic-v1/refresh-width-gdb.log`. Wider tests must be rerun.
  G23 padding follow-up: after the cold-reference fix, 18/24 scalar tests
  passed; open EOS cases failed 737 at C1. GDB/targeted probes showed compiled
  stencil reach=1 but allocated/certified ghosts=2. At xmin, Y(owner0)=
  0.21916373024191385, Y(ghost1)=0.48083626975808613 correctly encode the
  0.35 face target, while unused ghost2=0.21913417161825449 was not filled by
  the one-layer Dirichlet operator. `bc_thermo.cpp::inlet_face_cell` now caps
  the trace source layer at the compiled reach and extends that valid trace
  into extra allocation padding. It does not alter stencil width or mirror
  guards. The compact-stencil real-fill unit regression failed first in
  0.35 s. All category probes/includes were removed after diagnosis.
  The three variable-thermophysics scalar failures (5792 in the smallest-dt
  passive case) also occur on an independently built unchanged 51e4a3a
  baseline; baseline open MPI-1/2/4 passes. This separate pre-existing
  globalization issue is not repaired or counted as passing here. Evidence:
  `baseline-51-comparison.log`, `padding-reach-gdb.log`, `padding-red-test.log`.
  After padding repair the open EOS case advances, but its old conservation
  oracle still used molecular viscosity at the interior owner composition
  instead of the newly authorized fixed inlet face composition/temperature.
  `tests/integration/product_scalar_contract_experiment.cpp::open_budget`
  now independently evaluates mu(T_inlet,Y_inlet) only on the fixed inlet;
  conditional pressure-outlet backflow retains the old owner oracle. The
  inventory tolerance remains 1e-12 and no computed solver flux is substituted
  into the analytic oracle. This follows the same residual/transmissibility
  face-coefficient regression above; rerun the open MPI-1/2/4 cases.
  Final focused verification: all 21 scalar uniform/near_pure/open/ibm/
  capacity/signed/restart cases on MPI-1/2/4 and both BC thermophysics tests
  passed (23/23, 48.35 s); external `g23-final-focused-tests.log` records it.
  Mass-balanced outlet MPI-1/2/4 also passed (3/3, 1.23 s), recorded in
  `g23-final-outlet-tests.log`: 26 focused cases passed in total.
  The old standalone stale-ghost-authority MPI-2 test also fails identically
  on the clean, independently built 51e4a3a baseline (collective stale rejection
  and rollback role-handle assertions, 0.38 s). See `baseline-authority-test.log`.
  Along with the three variable-thermophysics cases, this is a documented
  pre-existing regression, not a passing test or a speculative G23 repair.
  Temporary probes are absent and `git diff --check` passed. Next gate is
  clean Release rebuild of both CLI and importer, fresh V1 mean-state import,
  and real GTMC short-100; these results do not establish actual acceptance.

- G22 (diagnostic only): clean 0a9ed82 accepted the first real GTMC step,
  20001 at t=0.015465336638081572, dt=3.5115962230679085e-9. Its evidence row
  passed runtime validation against the exact initial manifest. Short-100
  then failed the second step at stage 40 / 802 (composition); predictor
  theta=0 and high margin=-1.4821969375237396e-323. A two-step replay preserved
  the native 20001 checkpoint (generation-20001-2080879002468335, manifest SHA
  c48880d4f7da078f2112197e3752584de53c80407ec06bdc30c902a4b7192b01).
  Exact restart from it reproduces stage 40/802 in 20.21 s. Temporary
  `src/solver_thermophysical_predictor.cpp` probe under
  HUNDUN_GTMC_SPECIES_PROBE / [DEBUG-gtmc-species] records any negative species
  certified by the predictor, including rhoY, endpoints/rates and selected
  theta. Purpose: distinguish density-product underflow in admissibility,
  boundary reconstruction, and inter-stage field mismatch. No physical
  acceptance guard is relaxed. Remove probe/includes after diagnosis.
  The first step's development energy ledger also reports -391.97986645 W
  balance defect (-1.37647501854e-6 J accumulated); this remains an audit item,
  not a physical validation claim. Short and medium tests have not passed.
  G22 follow-up: the 20.24 s predictor probe reproduced 802 with no negative
  certified predictor species, falsifying that proposed cause. Its code and
  includes were removed without changing predictor numerics. The probe now
  sat in `src/physics_thermo.cpp::ThermodynamicsPlan::composition`, printing only invalid
  tuples and a bounded diagnostic stack under the same environment flag,
  prefix [DEBUG-gtmc-thermo-composition]. Temporary cstdio/cstdlib/execinfo
  includes support it and must be removed after locating the actual caller.
  G22 root-cause confirmation: the 19.89 s invalid-composition/stack probe
  identifies `src/bc_thermo.cpp::evaluate_cell`, reached through
  `BoundaryThermophysicalFaceClosure::close` from ProductDriver's stage-40
  `refresh_coupled_state`. A follow-up ghost/owner probe reproduced it in
  19.99 s on rank 16 at local cell (23,-2,21), reflected owner (23,1,21),
  local shape (39,62,38). CH4 ghost=-0x0.000000000060cp-1022 and
  owner=+0x0.000000000060cp-1022, with exactly zero arithmetic face mean.
  All six retries preserve that antisymmetry. This is an admissible physical
  owner and imposed pure-air face, but a nonphysical Dirichlet stencil value.
  `src/bc_apply.cpp::apply_span<dirichlet>` correctly computes 2*face-owner;
  `bc_thermo.cpp::evaluate_cell` incorrectly assumes this stencil must also
  be a physical EOS composition for the intended inlet use case. Its current
  public contract explicitly rejects nonphysical ghost p/h/Y, so fixing the
  conflict requires distinguishing physical boundary-face thermodynamics
  from discrete ghost extension, including downstream/certificate semantics.
  A diagnostic harness reusing the real BC unit fixture and apply/close chain
  reproduces the same 802 in 0.34 s with ordinary, non-subnormal values:
  owner=0.8, prescribed face=0.35, ghost=-0.1. This rules out underflow as the
  general cause. Harness/logs/probe patch are outside the source in
  `/home/administrator/gtmc-hundun-port-20260908/species-diagnostic-v1/`.
  All predictor, thermodynamics-stack, and ghost probes/includes have now
  been removed from source; no production numerical fix was applied for G22.
  Do not clip ghost or interior Y, relax EOS bounds, or silently replace
  derived ghost values by face properties under the old certificate contract.
  The existing acceptance-build executable still contains the final probe
  and is DIAGNOSTIC ONLY until a fresh clean-source rebuild. No run is active.
  Proposed next step: explicit physical face-state authority, boundary
  consumers using that state, preserved stencil mirrors and strict physical
  EOS, with regression/certificate tests before further full GTMC runs.
  After probe cleanup, rebuilt `v04_bc_thermo_test`; existing serial/MPI-2
  contract tests pass (0.35/0.37 s). The separate diagnostic red loop repeated
  the identical 802 in 0.34 s. Existing tests therefore do not cover this
  valid-face/nonphysical-mirror integration conflict. Source diff is docs-only.

- G21: G20's 42.82 s probe reproduced norm=1 in trace-CH4 rows with q=0
  and next=0; extended-precision RHS is nonzero below binary64 range. Global
  absolute composition error was 2.77412637507e-17 and mass pairing
  1.68844377692e-15. This is an unrepresentable update, not failed physical
  mass closure. A real ScalarMassRemap regression with unit cell mass, zero
  inventory, fixed inlet fractions 0.2/0.3 and denorm_min correction flux
  reproduced 10215 / residual=1 / iteration=127 in 0.37 s.
  `src/solver_scalar_mass_remap_detail.hpp` now retains raw residual while
  separately computing residual beyond the unavoidable bound
  diagonal*denorm_min/2 in long double. The 32-epsilon relative gate applies
  to this excess; no mixture-mass tolerance, clipping, normalization or
  iteration-cap increase is used. Ordinary rows have negligible floor.
  Reduction capacity explicitly includes the fourth norm component.
  `include/hundun/v04_app.hpp::DriverScalarTransportReport` and
  `src/core_product_freeze.cpp::advance` retain both raw and convergence
  residuals. `tests/mpi/solver_scalar_mass_remap_mpi_test.cpp` plus
  `tests/CMakeLists.txt` exercise actual underflow closure and verify ordinary
  representable influx still requires an update to exact 0.0002/0.0003.
  MPI 1/2/4 pass (1.13 s total). The G20 probe and its includes are removed;
  logs remain in terminal-diagnostic-v1. Stored equation/history semantics
  are unchanged; this recognizes attainable floating-point convergence.
  Broader scalar conservation and full GTMC reruns remain pending.

- G20 (diagnostic only, no numerical repair): 9a418ce formal short-100 failed
  before any accepted step after 4:15.56. Final retry reports stage 53 / 5792,
  masking an earlier failure after the now-complete terminal audit. A separate
  fixed-dt=7.0231924461358171e-9 case/restart reproduced stage 65 / 10215 in
  41.66 s, one attempt. Terminal audit is available and passes: EOS=0,
  C=1.67053076098e-15, E=1.57039423988e-15, gauge/mass=0,
  committed absolute CFL=0.00472513880 < 0.3. The source
  `src/solver_scalar_mass_remap_detail.hpp` temporarily records the globally
  worst final-iteration cell/species, q/next/RHS/diagonal/masses/six correction
  fluxes under HUNDUN_GTMC_REMAP_PROBE / [DEBUG-gtmc-remap], with direct
  cstdio/cstdlib includes. Purpose: distinguish trace-species relative-error
  stagnation, fixed-source flux pairing, and outlet donor noncontraction.
  Remove after diagnosis. Formal v3 tolerances/input remain untouched;
  diagnostic-only artifacts live outside source in terminal-diagnostic-v1.

- G19: clean 25f4812 formal-v3 single-step reached both pressure-energy gates
  (C=1.67053076098e-15, E=2.92562756237e-11, attempt 7, dt=7.02319244614e-9),
  then failed stage 60 / invalid_plan 1503 after 3:40.56, zero accepted steps.
  A new 8x8x8 real ProductDriver regression in
  `versions/v0.4/tests/mpi/solver_mass_balanced_outlet_mpi_test.cpp` reproduces
  exactly stage 60 / 1503 in 0.40 s. Launching this test under GDB and breaking
  at `solver_piso.cpp:8904` observed closed=false, global_sum[4]=0 and
  global_max[3]=0: all preceding terminal bindings passed, but the mandatory
  open-boundary closure witness was absent. Root cause:
  `src/core_product_freeze.cpp::local_pressure_outlet_closure` only counted
  pressure_outlet. It now validates the compiled Neumann relation/value source
  and audits actual ghost-versus-owner absolute pressure for
  zero_gradient_mass_outlet. Ordinary pressure-outlet face-versus-prescribed
  pressure audit is unchanged. The mandatory nonempty witness and all terminal
  tolerances remain intact. No stored-history semantics or probe added.
  Real terminal-product plus existing outlet/patch regressions pass on MPI
  1/2/4 (0.39/0.41/0.42 s). Full-size rerun pending; this is not acceptance.

- G18: clean ae435c4 passed native import but actual formal-v3 one-step failed
  after 5:50.20, zero accepted steps (stage 54 / 10210, nine attempts). G17
  reduced energy to 7.3898770144e-12, below its 1e-10 gate; continuity stalled
  at 8.3106765680e-14 above the unchanged scalar-pairing 16-epsilon gate.
  The full-alpha trial improves C to 2.6501587219e-14 while E increases only
  to 7.3907623741e-12. A direct real-selector regression using this captured
  tuple failed in 0.36 s. Root cause: equal-weight merit ignores unequal
  component terminal targets, letting already-converged energy roundoff
  dominate continuity descent. `versions/v0.4/include/hundun/v04_flow.hpp`
  adds a default-one energy merit weight to samples/selection certificates;
  `src/solver_pressure_energy.cpp` validates matching positive finite weights,
  signs the distinct v04pewt1 policy, and revalidates weighted merit without
  changing raw residual evidence. Default-one legacy provenance is preserved.
  `src/core_product_freeze.cpp` sets weight=C_target/E_target, binds it on
  candidate replay, and uses the same scaled energy in Aitken/inexact forcing.
  Both independent terminal gates, physics and refinement capacity are intact.
  `tests/numerical/solver_pressure_energy_globalization_test.cpp` covers the
  captured tuple, raw evidence, weight tampering, invalid/foreign weights,
  extrapolation, and preserved legacy rejection/provenance. All four serial /
  MPI 1/2/4 selector tests pass (1.44 s). Product/full GTMC reruns pending.
  This changes selection policy, not stored rate/flux/history interpretation;
  method-history signature is deliberately unchanged per its contract.

- G17: the shortened G16 full-size probe reproduced the same failure in
  1:51.68 (one attempt instead of nine), with worst energy cell (85,80,90),
  on the internal methane source plane, residual -6.4836 W. Source inspection
  confirmed `src/solver_enthalpy.cpp::assemble_enthalpy_impl` applies
  `add_source_convection_correction`, but
  `assemble_target_coupled_enthalpy_residual` omitted it. The candidate replay
  therefore reconstructed h through a solid placeholder instead of using h_in,
  while the linear direction used the prescribed source. A real assembly
  regression in `tests/unit/solver_ibm_equation_interface_test.cpp` failed in
  0.48 s: full residual -1.254 W, candidate -0.00278667 W for the identical
  state and flux. Adding the same source correction to the candidate path
  made it pass (0.59 s). The regression also compares the complete residual
  arrays. The final no-probe batch passed all three tests: source interface,
  existing enthalpy terms, and app initialization/restart (8.49 s combined).
  The G16 probe and its direct includes were removed; its exact patch
  and rank logs are preserved outside source in `energy-diagnostic-v1`.
  `core_product_freeze.cpp::method_history_signature` bumps the source-state
  convection component v1 to v2 for patch cases. No tolerances were changed.
  Actual formal v3 advancement must now be rerun; this unit pass is not GTMC
  short/medium acceptance.

- G16 (active diagnosis, not a numerical repair): clean candidate 521232f
  imported/read back on 64 ranks, then its actual one-step run failed after
  7:35.88 with no accepted step. The earlier stage-44 outlet guard was passed;
  last failure is stage 54 / rejected-step 10210 after C2 and 12 pressure-energy
  refinements. Energy stalls near 2.04129538e-5 while continuity descends.
  `src/core_product_freeze.cpp` temporarily adds an environment-gated
  `[DEBUG-gtmc-energy]` probe at the already-scanned worst energy cell for
  alpha 0/1; it reports state, residual, correction and six-face flux with no
  extra collective or field mutation. Direct `<cstdio>/<cstdlib>` includes
  support the probe. It must be removed after diagnosis. A separate
  `case-energy-diagnostic-v1` clamps dt to the final failing dt and disallows
  smaller retries to shorten the reproducer; it is explicitly NOT acceptance.
  Formal v3 inputs/tolerances are unchanged. Current hypotheses are missing
  outlet-flux energy derivative, inlet energy residual/Jacobian mismatch, or
  mismatched live-versus-candidate energy state; none is yet confirmed.

- G15: the committed b6eb072 64-rank native import crashed before writing a
  restart. The same crash was minimized to the real 4-rank outlet/patch test
  (0.39 s). `addr2line` located `PressureEnergyCandidateBoundaryFinalizer::bind`;
  disassembly of the faulting load matched the parent mass-target lookup.
  `bc_compile.cpp` allocates parameter arrays only for locally owned faces,
  while the new global patch catalog is present on every rank. MPI 1/2 happened
  to leave every rank owning x_min; MPI 4 exposed non-owning ranks. The binder
  now checks parent targets only on owners, validates the local parameter
  index, and keeps collective descriptor equality/error consensus for all
  ranks. `tests/CMakeLists.txt` now includes MPI 4 for this regression.
  After the repair MPI 1/2/4 all passed (1.14 s combined).
  The b6eb072 import failure log is retained; its short/medium gates never ran.

- G12: COAST `s76/SRC.Coast/boundary2_dp.F90`, `bndry2dp.F90`,
  `Calcmassflowrate.F90` were read on the remote host. Marker -2 uses Neumann
  pressure correction and scales the existing signed outlet flux pool by
  `(flow + summ)/flow`, with `summ` equal to minus the global continuity defect
  including density storage. It is not the -60 characteristic/backflow model.
  `include/hundun/v04_case.hpp`, `src/app_case.cpp`, `src/bc_compile.cpp` add
  explicit `zero_gradient_mass_outlet` (enum appended, existing wire layout
  unchanged). Absolute pressure remains the EOS reference, not a Dirichlet
  correction. U/h/Y extrapolate. Old pressure-outlet/backflow behavior is kept.
  `src/core_product_freeze.cpp` enables the new finalizer path;
  `src/solver_piso.cpp` routes the new outlet as a mechanical transport face.
  `include/hundun/v04_flow.hpp` and `src/solver_piso.cpp` provide privately
  certified candidate BDF mass storage to `src/solver_candidate_boundary.cpp`.
  The latter rescales the positive net signed outlet pool, preserving local
  reverse flow, with no area fallback or clipping; unclosable pools reject
  atomically. A focused regression exposed that the existing `pressure_input`
  is deliberately corrector-2-only, so candidate history is separately bound
  at each successful pressure assembly. A second regression exposed the
  predictor's homogeneous-Neumann => fixed-flux shortcut; the new adjustable
  outlet now bypasses that shortcut without changing pressure response.
  Focused regression results are recorded below; actual GTMC gates remain open.
- G13: the two-patch real finalizer regression requested 0.125/0.250 kg/s but
  returned 0.1875/0.1875: the finalizer merged separately metered inlets on one
  Cartesian face. `v04_flow.hpp::PhysicalMassFlowPatch` introduces an immutable
  deep-copied membership/target binding. `solver_candidate_boundary.cpp`
  validates matching collective descriptors, local support ownership,
  non-overlap and parent total, fingerprints membership, and independently
  normalizes each patch using candidate EOS inlet density. Unlabelled faces
  stay zero. `core_product_freeze.cpp` passes the compiled label support and
  reserves enough reduction scalars for the global patch catalog. Its
  `method_history_signature` increments labelled-mass-inlets v1 to v2, so
  checkpoints from the old merged-patch finalizer are not treated as compatible.
  `tests/mpi/solver_mass_balanced_outlet_mpi_test.cpp`,
  `tests/support/candidate_boundary_fixture.hpp`, `tests/CMakeLists.txt` contain
  the real-chain MPI 1/2 regression. These tests are not GTMC acceptance.
- G14: rebuilding the legacy PISO MPI regression with GCC 11 failed because
  `tests/mpi/solver_piso_mpi_test.cpp::test_exact_eos_correction_collective_transaction`
  used a local constexpr variable as a lambda default argument. The test now
  passes that reference explicitly at each call. No production numerical code
  or test tolerance changed. The new parser/wire roundtrip in
  `tests/unit/app_case_test.cpp` passed on 2026-09-08.
- Focused verification: the new outlet/patch MPI 1/2 tests passed including
  nonzero BE density storage, Neumann pressure/pressure-correction ghosts,
  uniform signed reverse-flow scaling, and atomic zero/negative-pool rejection.
  Existing pressure-boundary tests (serial and MPI 1/2/4), IBM physical boundary
  flux MPI 1/2/4, app initialization/restart, and core patch-inlet tests passed.
  Legacy `v04_solver_piso_mpi_1` passed, but MPI 2/4 failed the existing negative-
  density / injected-Halo provenance assertions. A separate worktree at
  `18ce94c`, built with the same GCC 11 RelWithDebInfo configuration and only
  the identical G14 test compilation repair, reproduced those same failing
  groups. Both versions report candidate-boundary/frozen-candidate/exact-EOS
  groups passing. This is a verified pre-existing test-suite issue, not a
  passing regression suite; do not conceal it in mainline synchronization.
- Input v3 is a copy of the failed v2 with only `y_max.flow_kind` changed to
  `zero_gradient_mass_outlet`; all geometry, species, time controls and solver
  tolerances are unchanged. case.json SHA-256:
  `978089748ce71bafb436bd89ee49b21d597c81ce9fc422e1eb8befa552b51ebd`.
  Its updated `mean_field_case_receipt.json` links the unchanged v2 source
  fingerprint and the user's authorization. A new native import is required;
  the v2 checkpoint must not be relabelled as a v3 restart.


| ID | Cause / evidence | Files and symbols | Verification / synchronization notes |
|---|---|---|---|
| G01 | GCC 9.4 build failed because required standard declarations arrived only through transitive includes on the newer local compiler. | `versions/v0.4/include/hundun/v04_flow.hpp`; `src/bc_apply.cpp` (`<algorithm>`); `src/io_restart.cpp`; `src/solver_initialization_projection.cpp` (`<climits>`), all under `versions/v0.4`. | Remote GCC 9.4 Release build passed after direct includes. No numerical change. |
| G02 | Optimized COAST-axes parsing observed a rounded axis endpoint but unrounded declared endpoint despite nested casts. For 0.13710429672, the float32 authority is 0.13710430264472961. | `versions/v0.4/src/app_case.cpp`: `coast_binary32_value`, COAST axis parsing / domain normalization; `tests/unit/app_case_test.cpp`: non-binary endpoint regression. | Local Release GTMC parser and app-case unit test passed. The observation is verified; a standalone compiler-defect reproducer has not been established. Preserve explicit binary32 storage for COAST input only. |
| G03 | Whole-face boundary schema cannot preserve 20 separately metered air patches and the immersed methane inlet. Hashing a generic data file alone gives no numerical patch semantics. | `versions/v0.4/include/hundun/v04_case.hpp`: `PatchInletSpec`, `PatchInletsSpec`, `ValidatedModel::patch_inlets`; `src/app_case.cpp`: parser, secure label-file validation, immutable hash, extension wire 19; `tests/unit/app_case_test.cpp`: `test_patch_inlets_case_and_wire`. | Schema tests passed before adding G04. Product compiler currently deliberately rejects patch cases until source and boundary bindings are fully integrated. Existing case wire versions 9–18 remain unmodified. |
| G04 | STL is closed, not open: all 34,272 edges have incidence two. Global scan-hit deduplication merges two component crossings at y=-0.0135. Even component-local nonzero-winding union (bit-identical x/y/z scans) differs from frozen marker in 46,607 cells and changes 146 fuel backing solids. | `versions/v0.4/include/hundun/v04_case.hpp`: `ImmersedBoundarySpec::marker_file/marker_fingerprint`; `src/app_case.cpp`: explicit `imported_cartesian_marker` geometry mode, secure 0/1 file validation/hash/wire. Imported geometry compiler tracked separately. | In progress. Do not silently replace frozen marker with STL parity/union classification. Existing STL mode remains default; imported mode explicitly selects Cartesian cell-face wall geometry. STL is retained as hashed provenance, not reclassified geometry. |
| G05 | Existing IBM interface equation plan assumes every fluid–solid interface is impermeable, deleting prescribed fuel flux and using solid placeholder convection state. | `versions/v0.4/include/hundun/v04_ibm.hpp`, `src/solver_ibm_equations.cpp`, interface test; further equation consumers in progress. | Source-flux and full-state interface unit vertical slices passed; integration and acceptance pending. Detailed log: `docs/gtmc-port-ibm-inlet.md` when published. |
| G06 | One scale for an entire external face loses independent passage mass flows; accepting a label file without checking owner/backing geometry can inject into solids. | `versions/v0.4/src/core_patch_inlets_detail.hpp`: `read_case_binary`, `CompiledPatchInlets`, `compile_patch_inlets`; `src/core_product_freeze.cpp`: `resolve_static_boundary_values`, `ProductCompiler::compile`; `tests/unit/core_patch_inlets_test.cpp`; `tests/CMakeLists.txt`. | `ctest --test-dir build/gtmc-tests -R '^v04_core_patch_inlets$' --output-on-failure` passed (2026-09-08). Verifies real IBM link binding, one-face exact oriented source, outer support exclusion, changed-file rejection, missing negative backing rejection, and unsupported heterogeneous outer temperature rejection. Product integration remains guarded pending full equation wiring. Outer patches currently require shared parent-face thermal/scalar targets, explicitly rejected otherwise; each patch retains independent direction and mass scale. |
| G07 | PISO terminal CFL rejected any nonzero pressure-inactive face; pressure convergence audit returned zero there, dropping physical fuel inflow. A fixed mass source must be counted in continuity without opening a pressure correction link. | `versions/v0.4/src/solver_piso.cpp`: `PressureVelocityCoupler::audit_pending_final`, `audit_pressure_convergence`. | In progress; depends on exact source lookup API. Terminal permits only the exact prescribed value; convergence returns that fixed value for every pressure scale and rejects corrupted predictors. Pressure matrix / MG activity remains unchanged. Integration regression pending. |
| G08 | The physical boundary energy ledger reconstructed h and kinetic energy from ordinary Cartesian solid placeholders at the fuel face. | `versions/v0.4/src/core_conservation_detail.hpp`: `collect_boundary_balance`; `src/core_product_freeze.cpp`: audit call and `method_history_signature`. | Uses the same sparse source-face reconstruction correction as the equations for h and K, without adding mass twice. Source-patch cases get a distinct history signature; pre-existing cases retain the old signature. Integration regression pending. |
| G09 | D5 requests Vreman with wall functions off, but Hundun only exposed Vreman coupled to its equilibrium wall function. | `versions/v0.4/include/hundun/v04_case.hpp`: appended `TurbulenceKind::vreman`; `src/app_case.cpp`: `parse_turbulence` and wire enum validation; `src/physics_vreman.cpp`: `TurbulencePlan::compile`; `tests/unit/physics_vreman_test.cpp`, `app_case_test.cpp`. | Adds explicit `"model":"vreman"` using the existing Hundun Vreman kernel and resolved-wall treatment. Does not change the default or old enum values. Does not claim numerical identity to COAST's different Vreman discretization. Tests in progress. |
| G13 | The authorized VTK mean state must enter through a native checkpoint without relabeling COAST `dp/nvf/raw h`, and absent PDF/history must force one BE recovery step rather than be promoted to exact BDF history. | `tools/v04_coast_restart_import.cpp`; `versions/v0.4/CMakeLists.txt` (`v04_coast_restart_import`); `versions/v0.4/tests/integration/coast_restart_import_cli_test.py`; `versions/v0.4/tests/data/gtmc_hundun_thermophysics.d`; importer-side `gtmc-tools/prepare_gtmc_mean_case.py`, `gtmc_hundun_thermophysics.d`. | GCC 11 Release target and MPI 1/2 focused CTests passed. Full n64 native bridge/readback passed: V1 BE recovery, exact cell activity and field/phi payload; 250 fuel source faces sum to 0.0006966666667 kg/s. Restart `/tmp/hundun-gtmc.Itnz3O/gtmc-native-restart-n64-v1`; manifest SHA `08231ab...8465`, audit SHA `ade2b8...b68c`. Mass change is -1.78238605e-8 kg and total-energy change -24.2520370 J. The first driver is destroyed before readback compilation; max/rank conservative estimate is 224,701,632 B. This is bridge evidence, not the outstanding short/medium run acceptance. |

## G14 — independent committed development diagnostics

Ordinary application conservation observations were available in
`DriverStepReport` but not persisted unless a specialized runner was used;
the ordinary monitor was coupled to potentially large Visit output and did
not contain the boundary ledger. The default-off `--diagnostics-interval`
switch adds a small independent `diagnostics.jsonl`, with no new field scan,
reduction, numerical gate, or change to the V8 evidence schema.

Changed locations: `include/hundun/v04_app.hpp::ApplicationRunOptions` (appended
option); `src/app_main.cpp` (CLI parse/usage); `src/app_driver.cpp` (cold
collective control agreement, snapshot selection, committed finite-report
serialization through the existing MonitorWriter); all under `versions/v0.4`.
The public conservation comments now include prescribed IBM inlets. The new
schema is `HUNDUN_V04_DEVELOPMENT_DIAGNOSTICS_V1`, explicitly ineligible for
benchmark statistics. `docs/api/cli.md` and `docs/user-guide/diagnostics.md`
document defaults, units, sign and restart-epoch semantics.

Regression locations: `tests/integration/app_init_case_test.cpp` checks
default-off behavior and two committed diagnostic records with Visit disabled;
`tests/mpi/app_control_contract_mpi_test.cpp` adds divergent diagnostics
interval to cold collective rejection. The new test linked against a previous
archive observed the expected missing-ledger failure, but emitted an ODR
warning from the changed options type; this is **not accepted as a clean test
receipt**. `app_driver.cpp` compiled independently against current headers.
After the source-face implementation reached a compilable checkpoint, the
consistent GCC 11 build and `ctest -R '^v04_app_init_case$'` passed (8.77 s).
The added MPI control-agreement tests also passed on 2 and 4 ranks (0.76/0.74 s);
no GTMC acceptance step is represented by these application fixtures.

## COAST comparison evidence

The remote COAST implementation is consulted as the reference for field layout,
marker semantics, patch ownership, and geometry decisions. The immutable d5
marker is not assumed reproducible from the current mutable COAST source.
`imb_mesh_minimal.F90` computes the d5 default normal from signed six-neighbor
solid directions and wall distance from half the neighboring-center spacing;
the closest-STL metric is conditional, and `closest_triangle_metric` defaults
false. D5 disables wall functions and additional damping. These facts motivate
an explicit Cartesian-marker geometry mode, not a claim of smooth-wall numerical
identity between solvers. Source paths/line anchors and hashes are retained in
the geometry evidence receipts and will be consolidated before handoff.

## Required completion evidence

- Per-file patch/commit manifest, with reason, symbol/line anchors and test result.
- Native restart bridge: complete 128-block coverage, pressure semantics,
  species/enthalpy convention, and native readback proof.
- Exact frozen marker and inlet support; prescribed mass and energy ledgers.
- Local and remote compiler/test receipts; clean source/build runtime identity.
- Actual GTMC short 100-step gate and medium 3000-step gate with 500 development
  steps followed by reset and 2500 statistics steps. No unit/probe substitution.

### Development-run acceptance interpretation

The user explicitly confirmed software development testing. The ordinary
`hundun run` CLI intentionally writes `statistics_eligible=false`
(`src/app_driver.cpp`); it has no benchmark accumulator reset command. Preserve
that flag and do not impersonate Task-20 benchmark eligibility. For this GTMC
development medium test, "reset + 2500 statistics steps" means discard the first
500 committed development steps and start a fresh **postprocessed diagnostic
window** for the remaining 2500, not physical/statistical convergence proof.
Require actual contiguous committed steps, native restart readback, unchanged
candidate/input identities, the shipped runtime-evidence validator (with exact
run-start manifest), finite bounded thermodynamic/species state, and all enabled
terminal EOS/continuity/energy/CFL contracts. Record retries, limiter activity,
minimum dt, mass/energy ledgers and timing; never disable a gate to label a run
passed. The short/medium windows remain outstanding.

## Latest verified checkpoints

- Actual short-run failure and minimal replay (candidate f10f15a): both ended
  at step 20000, first numerical_failure/985 in predictor stage 10 at
  dt=4.4948431655269229e-7, last rejected_step/10464 in stage 44 on attempt 9
  at dt=1.7557981115339543e-9. No new step or Evidence row was committed.
  Wall times were 36.09/36.63 s, zero swap. Executable SHA-256
  `a48c2e5990b5898ead1c4f2fb5089704af49574fdfb1c39a7be37ccddd8bebd7`;
  target manifest SHA-256
  `d08b929a947e7319823791bde3dde1bd6460b0577a494ec33f8462ca6ab87cee`.
  `core_product_freeze.cpp::evaluate_pressure_energy_candidate` maps the
  finalizer status to 10210+254=10464. The only rejected_step returns in
  `solver_candidate_boundary.cpp::prepare_physical_boundary_flux` are the
  pressure-outlet backflow checks; with this case's allow_backflow=false,
  line 646 rejects negative outward mechanical flux. The frozen COAST input
  explicitly selects zeroGradientOutlet/marker -2 and documents its mass-
  closure path, not the pressure/characteristic -60 path. The previously
  documented allow_backflow=false observation was incomplete: ghost closure
  remains extrapolated, but final physical-flux acceptance rejects reversal.
  No rejection guard or numerical tolerance was relaxed.
- Source inspection also identifies an outstanding per-patch finalizer seam:
  `resolve_static_boundary_values` uses individual air-patch capacities, but
  `prepare_physical_boundary_flux` still recomputes one mass-flow capacity for
  the complete Cartesian face. Its finalization must gain patch membership
  and a focused regression before GTMC acceptance can be claimed. This is a
  source-inspection finding, not yet a separately executed failing fixture.
- Pre-candidate consistent-build regression batch passed 6/6:
  imported-marker MPI 1/2 ranks (0.67/0.60 s), app init/run/restart/diagnostics
  (8.23 s), resolved Vreman (0.34 s), compiled patch inlets (0.43 s), and case
  parser/wire (0.39 s). Source-interface, thermophysical IDP and pressure-energy
  Schur tests also passed after source certificate completion. The final
  acceptance executable will use a fresh GCC 11 Release build; existing test
  build directories are RelWithDebInfo, not performance evidence.
- Full native bridge (64 ranks, GCC 11 Release) passed after removal of 15901.
  Output `/tmp/hundun-gtmc.Itnz3O/gtmc-native-restart-n64-v1`, generation
  `generation-20000-2069446777267044`; manifest SHA-256
  `08231abeb01961c78d9dd216a633a480b8261b626f24c7d4a7ccafb2c2888465`.
  All 3,901,529 fluid / 1,804,912 solid cells were retained. Native V1 readback
  verified BE recovery and zero field/flux error. Exactly 250 previously zero
  internal faces received the prescribed fuel source, summing to
  0.0006966666667 kg/s. Native EOS reconstruction changed mass by
  -1.7823860524174596e-08 kg, rho*h inventory by -24.252034855706693 J,
  kinetic energy by -2.1625042897038037e-06 J, and total energy by
  -24.252037018210982 J. These are recorded transfer changes, not claimed
  conservation or physical-equivalence results. The receipt
  `mean-field-transfer.json` SHA-256 is
  `ade2b8c46f182fa1fc4caa222c3cf04b83532ac4d8d0914b77517015a3afb68c`.
  Conservative single-product storage estimate was 224,701,632 B/rank;
  elapsed 21.91 s, no swap. Synthetic importer CTest 1/2 ranks passed again.
  This verifies native initialization/readback, not yet multi-rank advancement.
- Root consistent-build `v04_solver_ibm_equation_interface` passed (0.58 s),
  including fixed h_in freeze, generic and compiled canonical +0 derivative,
  stale ordinary-certificate rejection, and the G10/G11 tests below.
- Input adapter tests `python3 -m unittest discover -s
  /tmp/hundun-gtmc.Itnz3O/gtmc-tools -p 'test_*.py'` passed 12/12 (18.085 s).
  A source-only adapter copy is retained in
  `/home/administrator/gtmc-hundun-port-20260908/adapters`; no credentials or
  raw production fields are embedded in the source patch.

- G12 (2026-09-08, case-input adapter): the COAST geometry label file includes
  3,105 air-labelled solid cells. It is retained byte-for-byte as provenance
  (`4af598b61caf328617e28990832a905d79b7a5dcc50b43879dc3fe4f5fb7b79c`).
  `gtmc-tools/prepare_gtmc_patch_base.py` creates typed `inlet_owner_labels.d`
  by intersecting only external air labels with the frozen fluid marker;
  fuel owner/backing labels must already match and are never silently filtered.
  Typed owner SHA-256 is
  `a9a566118864ea47a18bf7f98c2ddcdbedac80e3d23f01cdf93228fd94b1cc4a`.
  Input base is `/tmp/hundun-gtmc.Itnz3O/case-exact-marker-patches-v1`.
  The 8 inner and 12 outer name/label mappings and per-patch mdot were read
  from remote frozen `boundary_mesh_check.txt` / `boundary_conditions.d`.
  Air total is 0.019740000003995995 kg/s, fuel 0.0006966666667 kg/s.
  The old homogeneous geometry-only probe incorrectly combined these on y_min;
  it remains diagnostic-only. Marker and original labels were not modified.

- G10 (2026-09-08, inlet diffusion/gradient and limiter integration):
  `src/solver_ibm_equations.cpp::correct_velocity_gradient` treated full-state
  sources as zero-velocity walls. The new constant-prescribed-U regression in
  `tests/unit/solver_ibm_equation_interface_test.cpp` failed before the fix and
  passed afterward (`v04_solver_ibm_equation_interface_test`, exit 1 -> 0).
  `inlet_for_link` now binds fixed U in that gradient and `constrain_momentum`;
  source links do not apply the equilibrium wall law. Resolved Dirichlet
  viscous traction is retained (not skipped), and `inlet_viscous_work_input`
  reports the matching physical boundary work to
  `src/core_conservation_detail.hpp::collect_boundary_balance`.
  The explicit development inlet diffusion policy is zero thermal/species
  diffusive flux, with h/Y prescribed for advection; source thermal correction
  removes the ordinary solid-material coupling exactly. This is a declared
  inlet discretization, not a COAST-equivalence claim. Constant and linear
  velocity viscous-work oracles and both thermal closures with arbitrary hot
  solid placeholders passed in the source-bearing interface test. The compile
  API also now rejects inward phi with outward Cartesian normal U atomically;
  this regression failed before its direction check and passed afterward.
- G11 (2026-09-08, in progress): `src/solver_momentum.cpp` now subtracts the
  source-aware low-order upwind contribution when producing the high-low delta.
  The limiter's `MomentumPredictorLimiterWorkspace` in `include/hundun/v04_flow.hpp`
  carries an optional IBM authority; `src/core_product_freeze.cpp` binds it.
  Source flux is validated exactly and counted in CFL, while the pressure graph
  and antidiffusive source correction remain inactive (fixed high/low inlet
  states have zero difference). `EquationAssemblyCertificate::inlet_sources`
  and `src/solver_equations.cpp::AssemblyEpoch::record` bind/check that source
  identity, and the limiter rejects missing source authority. Source-bearing
  full momentum assembly and limiter regression passed: high-low delta is
  zero at the fixed inlet, and its exact absolute CFL is dt*phi/(2*rho*V).
  Existing `v04_solver_equations_mpi_1`, `_2`, `_4` also passed. Full-case
  integration is still pending; these are not GTMC acceptance steps.
- G11 follow-up: a full assembly regression with the prescribed source phi
  removed failed (assembly incorrectly succeeded). `assemble_momentum_impl`
  now validates source-bearing interface flux before writing any equation
  output; the same regression and the full source FCT/CFL test passed after
  this fix. Sealed-only legacy assembly entry behavior is unchanged.
- Full CH4/O2/N2 input is
  `/tmp/hundun-gtmc.Itnz3O/case-exact-marker-mean-patches-v2` (case SHA-256
  `b0cc3c30106fc0405b85f64352aec83f72d9be4fff82044a06a2c599c64cef3e`).
  v1 is retained as diagnostic input: duplicate reference to the typed labels
  in both `mesh.data_files` and `patch_inlets` was correctly rejected with
  case detail 9. v2 keeps exactly one typed reference plus original-label
  provenance and advances to the intentional ProductCompiler guard 15901.
  No parser/security validation was relaxed.
- Outlet preflight of the imported mean field: j=246 has 6,199 frozen-fluid
  cells, 171 with negative v; v ranges from -8.354767799377441 to
  19.48097038269043 m/s. The current pressure-outlet input keeps extrapolated
  U/scalars (`allow_backflow=false` means no prescribed reservoir-state switch
  in `bc_compile.cpp` / `resolve_static_boundary_values`, not a velocity clip).
  Do not silently replace these source velocities or invent a backflow inlet
  speed merely to suppress the observed recirculation.

- `v04_app_case` passed after G03/G04 marker and patch wire extensions (2026-09-08).
- `v04_core_patch_inlets` passed after G06 (2026-09-08).
- `v04_physics_vreman` passed including explicit resolved Vreman (2026-09-08).
- The marker-only full-size diagnostic command
  `mpiexec -n 16 build/gtmc-tests/versions/v0.4/hundun validate /tmp/hundun-gtmc.Itnz3O/case-exact-marker-probe --dry-plan`
  initially failed with `invalid plan detail=1304` (`kQuadraticCoverage`). After
  explicit single-normal-band closure and donor-selection repairs it passed
  with sealed=1 on 16 ranks: 218,175 links/groups, 176,204 quadratic, 41,971
  linear, 85,329 expanded-search groups, 41,236 rank fallback, 735 coverage
  fallback (including 285 explicit one-normal-band closures), zero condition
  or donor fallback. The intermediate 1305 failure was traced to truncating
  genuinely 3-D candidate sets (e.g. 217 candidates) to 32 coplanar donors;
  preserving tangential diversity fixed that without a 2-D closure. The
  probe intentionally lacks the final species/inlet physics and must never be
  used as a short/medium acceptance substitute.
- COAST `finish.F90` restart order is `f(1..6), p, f(nvf)`: `nvdp=4`,
  `nvf=5`, `nvh=6`. `coast_legacy_driver.F90` clears `dp` before its solve and
  `update.F90` adds it to physical `p`. Therefore record 4 is the last pressure
  correction, not physical pressure perturbation. Record 7 matches all
  5,706,441 VTK static-pressure values bitwise. The bridge must use
  `double(float32_absolute_pressure) - reference_pressure`, preserving record 4
  only as an audit field.
- The frozen run enables PDF but has zero `restart_pdf.*` files. Its record-5
  `nvf` working field is approximately zero everywhere and cannot reconstruct
  PDF mean CH4. Raw h + VTK mean Y has 138 cells without an inverse-temperature
  root in 200–6000 K; some valid inverses reach 1448.55 K despite a cold VTK
  temperature near 295 K. Replacing h by h(VTK T, VTK Y) changes the full-fluid
  `sum(rho*V*h)` by **-24.25214663167891 J** (absolute changes 24.417999023478465 J).
  This is a nonconservative state transfer, not an equivalent restart. The
  user subsequently **authorized the mean-field transfer** for Hundun software
  development testing. Native h/rho must be generated consistently by Hundun's
  thermodynamics and the resulting energy/mass changes recorded. Raw h remains
  preserved as audit data. Evidence and exact COAST source anchors:
  `/tmp/hundun-gtmc.Itnz3O/gtmc-tools/CHANGELOG.md`,
  `/tmp/hundun-gtmc.Itnz3O/d5-import-v1/thermo_validation.json`.

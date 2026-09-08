# GTMC IBM prescribed-inlet port note

## Scope and evidence

This slice addresses the GTMC methane annulus represented by 250 internal
fluid--solid Cartesian links.  The approved marker/patch probe places the
fluid owners at `j=80`, the masked neighbours at `j=79`, and the source on the
`y_negative` face with positive Cartesian mass flux.  Consequently this is
not a physical-domain face and cannot use the existing whole-face inlet BC.
It must remain disconnected from the pressure/MG graph while carrying a
fixed, non-zero continuity flux and fixed transported inlet state.

The pre-change equation interface had one authority only:
`IbmEquationInterfacePlan::zero_interface_flux` cleared every cut link to
zero.  The bulk convection kernels then reconstructed U/h/Y from the masked
solid placeholder.  Thus changing only continuity `phi` would conserve mass
but would inject arbitrary momentum, enthalpy, and composition.  A fixed
inlet state also has zero directional variation, so the pressure-energy Jv
must not differentiate the solid placeholder.

CodeGraphF was not available in this workspace.  Investigation used `rg` and
direct source inspection as the documented fallback.  No MPI launcher was
run.

## Implemented changes

### Stable source authority and state

- `versions/v0.4/include/hundun/v04_ibm.hpp`
  - `IbmInterfaceMassFluxSource` is the compatible flux-only binding.
  - `IbmInterfaceInletState` binds one stable `global_link` to exact `phi`, U,
    h, and all independent species.
  - `IbmInterfaceInletFieldKind` identifies velocity, enthalpy, independent
    species, and kinetic energy.  Kinetic energy is derived from the same U,
    rather than accepted as a second potentially inconsistent input.
  - `IbmEquationInterfacePlan::compile` retains every old overload and adds
    flux-only and full-state overloads.  The explicit independent-species
    count keeps a rank with zero local source links semantically bound.
  - `prescribed_face_flux` exposes only exact fixed-flux lookup.  A sealed cut
    face returns false and the output argument is left unchanged.
- `versions/v0.4/src/solver_ibm_equations.cpp`
  - `compile_sources` resolves stable global links, sorts them, rejects
    duplicates/unknown links/non-finite state, and verifies that `phi` points
    from solid to fluid for the link orientation.
  - `constrain_interface_flux` first performs the original blanket retirement
    of cut and solid-incident faces, then reapplies only prescribed source
    values.  `zero_interface_flux` remains as a compatibility spelling and
    delegates to it.
  - `validate_interface_flux` requires exact prescribed source values and zero
    on all remaining cut links.
  - The no-source fingerprint path is intentionally unchanged; source bytes
    enter the fingerprint only when a source binding exists.

### One reconstruction/correction seam

- `versions/v0.4/include/hundun/v04_execution.hpp` and
  `versions/v0.4/src/solver_diffusion.cpp`
  - `reconstruct_cartesian_convection_face` dispatches to the exact existing
    `reconstructed_face` templates for central2, limited-central2, and TVD2.
    No limiter or metric formula is copied into IBM code.
- `versions/v0.4/include/hundun/v04_ibm.hpp` and
  `versions/v0.4/src/solver_ibm_equations.cpp`
  - `add_source_convection_correction` and
    `add_source_first_order_upwind_correction` add the sparse difference

    `scale * s_face * phi * (q_inlet - q_cartesian) / V_owner`

    to the fluid owner only.  `s_face` is `+1` for an owner plus face and `-1`
    for an owner minus face.  Multiple prescribed faces on one owner are
    aggregated and preflighted before write.  The same stored state selects
    U, h, Y, or derived `0.5*|U|^2`.
  - `override_source_face_values` supplies the same values to face-based
    consumers and writes canonical positive zero for fixed-state directional
    variation.

### High-order equation assembly

- `versions/v0.4/src/solver_momentum.cpp`, symbol
  `assemble_momentum_impl`: applies the common correction to all three U
  components immediately after Cartesian convection.
- `versions/v0.4/src/solver_enthalpy.cpp`, symbol
  `assemble_enthalpy_impl`: applies it to h before convection scratch is copied
  into the fused enthalpy residual.
- `versions/v0.4/src/solver_species.cpp`, symbols `assemble_transport` and
  `assemble_species_impl`: passes the independent-species index to the common
  correction.  Passive scalars deliberately remain unchanged because the
  GTMC inlet contract provides no passive-scalar state.
- Enthalpy and species assembly certificates mix the IBM plan fingerprint when
  this state authority participates, matching the existing momentum
  certificate discipline.

### Low-order predictors, limiter, and source-side interface physics

- `versions/v0.4/src/solver_thermophysical_predictor.cpp`, symbols
  `predict_high_local`, `predict_low_from_accepted_rate`, and
  `predict_low_local_donor`: applies the same sparse h/Y correction to the
  high-order rates, accepted-rate BDF fallback, and limited local-donor
  fallback.  Predictor certificates include the source fingerprint, and both
  accepted and previous fluxes are checked against the same interface plan.
- `versions/v0.4/include/hundun/v04_execution.hpp` and
  `versions/v0.4/src/solver_face_flux.cpp`: the test-only
  `constrain_pending_face_flux_for_test` helper lets a fixture publish a real
  committed source-bearing flux without exposing pending writer internals in
  the production API.
- `versions/v0.4/src/solver_momentum.cpp`: first-order source corrections make
  fixed U identical in the high/low pair; the limiter counts prescribed flux
  in CFL but leaves its anti-diffusive and MG graph edge cut.  The assembly and
  limiter certificates reject a mismatched or omitted inlet authority.
- `versions/v0.4/src/solver_ibm_equations.cpp`, symbols
  `constrain_momentum`, `correct_velocity_gradient`,
  `correct_zero_normal_diffusion`, and `inlet_viscous_work_input`: source links
  use prescribed U rather than stationary-wall U, do not enable a wall
  function, contribute consistently signed viscous work, and use zero
  diffusive h/Y flux.  Non-source links retain the existing wall treatment.
- `compile_sources` additionally rejects a prescribed mass flux whose normal
  velocity points the opposite way.  This prevents one link from binding
  inward mass with outward momentum while preserving transactional compile.

### Frozen source state and pressure-energy Jv

The observed pre-change failure was reproduced, not inferred: an ordinary
`freeze_cartesian_target_convection_faces` succeeded, replacing its source
byte by `h_in` succeeded, and
`differentiate_frozen_cartesian_target_convection_faces` then returned
`invalid_plan` because its bitwise reconstruction check saw stale numeric
data.  Therefore an unsealed post-freeze overwrite was not a valid fix.

- `versions/v0.4/include/hundun/v04_execution.hpp` and
  `versions/v0.4/src/solver_diffusion.cpp`
  - `FrozenConvectionFixedFace` and the trailing
    `FrozenConvectionFaceField::{fixed_faces,fixed_face_authority}` describe a
    borrowed, sorted sparse fixed-face schedule.  Empty schedules retain the
    old revision/local-binding calculation exactly.
  - `seal_fixed_cartesian_target_convection_faces` rechecks every non-fixed
    face against the Cartesian reconstruction, invalidates the ordinary
    certificate before resealing, and mixes the rank-local schedule/authority
    into revision and local binding without changing the collective
    reconstruction token.
  - Generic central2, limited-central2, and TVD2 differentiation consume the
    same schedule and write canonical `+0` on a fixed face.  The compiled
    limited-central2 path records a dedicated fixed-branch sentinel and
    replays the same `+0`; fixed faces are differentiable constants, not
    semismooth/generalized faces.
- `versions/v0.4/include/hundun/v04_ibm.hpp` and
  `versions/v0.4/src/solver_ibm_equations.cpp`
  - `freeze_source_convection_faces` performs ordinary freeze, sparse inlet
    value replacement, reseal, and exact source validation as one member API.
  - `validate_frozen_source_face_values` requires the exact IBM plan-owned
    schedule and bitwise h/U/Y/K value, so a finite but forged fixed byte is
    rejected before an operator binds it.
- `versions/v0.4/src/solver_pressure_energy.cpp`
  - `PressureEnergyEnthalpyOperator::bind` requires that validation when a
    local source exists.  Generic and compiled E_h now obtain fixed `+0`
    directly from the common frozen-face machinery instead of independent
    solver-local overwrites.
  - `PressureEnergyPressureFluxOperator::apply_after_exchange` deliberately
    keeps source faces inactive: the absolute residual contains fixed
    `phi_source*h_in`, while the pressure Jv has
    `delta(phi_source)=0`.  The MG/pressure graph is never reopened.
  - `enthalpy_collective_fingerprint` no longer mixes the rank-local boolean
    `has_inlet_sources()`.  GTMC partitions may own fuel links on only some
    ranks; exact interface identity remains in the rank-local binding revision,
    while case-wide inlet semantics come from `equation_semantics`.
- `versions/v0.4/src/core_product_freeze.cpp`, symbol
  `assemble_pressure_energy_residual`: source-owning ranks call the IBM freeze
  member; other ranks retain the ordinary freeze.  The following consensus is
  common to both branches, and neither branch itself performs a collective.

### Focused tests

- `versions/v0.4/tests/unit/solver_ibm_equation_interface_test.cpp`, symbol
  `test_prescribed_interface_mass_flux`
  - exact source `phi` and canonical `+0` sealed links;
  - failure on missing source flux or non-zero sealed flux;
  - exact source-face lookup and sealed-face absence;
  - one U/h/Y state and derived kinetic energy;
  - fixed-state directional value `+0`;
  - exact correction of the fluid-owner divergence for all three production
    reconstruction schemes and first-order upwind;
  - the real legacy freeze/overwrite/stale-certificate failure;
  - exact frozen nonlinear `h_in`, canonical generic E_h `+0` for all three
    schemes, compiled limited-central2 `+0`, and rejection after h_in tamper;
  - inward `phi` with oppositely directed U is rejected atomically;
  - source-bearing momentum assembly/FCT/CFL and inlet
    velocity-gradient/viscous-work/zero-diffusion policies;
  - the pre-existing no-source IBM regression still runs afterward.
- `versions/v0.4/tests/numerical/solver_thermophysical_idp_test.cpp`, symbol
  `test_prescribed_ibm_source_predictor`: publishes an authoritative source
  flux and checks the analytic high-order rho/h/Y result, then forces the
  accepted-rate BDF low route and checks the same source state there.
- `versions/v0.4/tests/numerical/solver_pressure_energy_schur_test.cpp` remains
  the no-source regression for generic/compiled frozen convection and the
  pressure-energy/Schur contracts.  Source-specific frozen generic/compiled
  assertions live in the IBM interface fixture above.

The frozen-source test was first compiled before the new API and failed on the
missing members/fields.  The retained legacy probe then ran green while
demonstrating the actual stale-certificate rejection; the new source-aware
assertions passed after implementation.

## Verification performed

Build directory:
`/tmp/hundun-gtmc.Itnz3O/build-ibm-source-face-gcc11`, GCC/G++ 11,
RelWithDebInfo.  The project still enabled its configured IPO/LTO policy.

Commands and results:

```text
cmake --build . --target hundun_v04_test_core -j2
PASS (complete core static library)

cmake --build . --target v04_solver_ibm_equation_interface_test -j2
PASS (target built)

./versions/v0.4/tests/v04_solver_ibm_equation_interface_test
PASS; accepted-IBM-thermal-rate-mismatch=0,
IBM-donor-thermal-Jv-error=2.49793e-09,
IBM-solid-material-fluid-response=4.54747e-13

cmake --build . --target v04_solver_thermophysical_idp_test -j2
./versions/v0.4/tests/v04_solver_thermophysical_idp_test
PASS; source_bundle_oracle alpha=0.504881 expected=0.504881

cmake --build . --target v04_solver_pressure_energy_schur_test -j2
./versions/v0.4/tests/v04_solver_pressure_energy_schur_test
PASS; frozen-target p/h/mixed finite-difference and recovery diagnostics pass

git diff --check
PASS
```

All commands used direct executables with their `MPI_COMM_SELF`/single-process
paths.  No `mpiexec`/`mpirun` command was used.

## Mainline compatibility

- Existing no-argument IBM compilation and the public
  `zero_interface_flux` call remain source-compatible.
- A plan with no prescribed sources retains the old zero-cut-face behaviour
  and fingerprint, protecting restart/evidence identity for existing cases.
- A flux-only source plan can be used by continuity/PISO without inventing
  thermodynamic state.  Equation transport rejects a non-empty flux-only plan
  if a transported source value is requested, rather than reading a solid
  placeholder silently.
- Pressure/MG topology was not changed in this slice.  Source faces remain cut
  graph edges; `prescribed_face_flux` is a separate read-only authority for
  endpoint/CFL/ledger code.
- New frozen metadata is appended to existing aggregate structures, and an
  empty fixed schedule returns the pre-existing frozen revision/local-binding
  identities.  Existing no-source callers require no changes.
- No case JSON/schema, runner, restart, STL, or imported-marker format is
  defined here.  Those adapters bind into this core API rather than leaking
  case concepts into equation code.

## Remaining verification and out-of-scope work

- The source/no-source collective fingerprint rule is established directly by
  code structure, but this focused suite does not construct a real two-rank
  E_h Krylov bind with a source only on rank 0.  The imported-marker n16 bridge
  validates initialization, not this Krylov path.  The first multi-rank GTMC
  advance must assert equal collective E_h fingerprints and differing valid
  rank-local revisions before accepting a time step.
- The product bridge, case schema, imported marker, conservation ledger, and
  restart importer are owned by adjacent port slices and require their own
  signed-off tests.  This note records only the equation/source-face seam.
- Before a 128-rank GTMC attempt, run a serial/tiny non-production advance with
  exact mass, h, Y, kinetic-energy, and viscous-work ledgers, then a small
  multi-rank advance exercising source-owning and non-owning partitions.
- The COAST raw restart lacks compatible rho/T/Y/PDF state.  The user has
  authorized VTK mean-field transfer; conversion and provenance validation
  remain outside this IBM source-face slice.

// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_boundary.hpp"

#include <mpi.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>
#include <string_view>

namespace {

using namespace hundun::v04;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

BoundaryFaceSpec wall() {
  BoundaryFaceSpec value;
  value.flow_kind = BoundaryKind::no_slip_wall;
  value.thermal_kind = BoundaryKind::adiabatic_wall;
  value.mach_limit = 0.95;
  return value;
}

ValidatedModel open_model() {
  ValidatedModel value;
  value.fingerprint = 0x12345678U;
  value.pressure_reference = PressureReferenceKind::boundary_absolute;
  for (BoundaryFaceSpec& face : value.boundaries) {
    face = wall();
  }
  BoundaryFaceSpec& inlet = value.boundaries[0U];
  inlet.flow_kind = BoundaryKind::velocity_inlet;
  inlet.thermal_kind = BoundaryKind::none;
  inlet.velocity = Real3{2.0, 0.0, 0.0};
  inlet.temperature = 300.0;
  BoundaryFaceSpec& outlet = value.boundaries[1U];
  outlet.flow_kind = BoundaryKind::pressure_outlet;
  outlet.thermal_kind = BoundaryKind::none;
  outlet.pressure = 101325.0;
  value.schemes = SchemeSpec{};
  value.time = TimeControlSpec{};
  return value;
}

void add_scalar_catalog_and_closures(ValidatedModel& model,
                                     std::string_view stable_name,
                                     TransportedScalarRole role,
                                     CartesianFace selected_face,
                                     ScalarBoundaryKind selected_kind,
                                     double selected_value) {
  model.transported_scalars.push_back(
      TransportedScalarSpec{std::string(stable_name), role});
  for (std::size_t index = 0U; index < model.boundaries.size(); ++index) {
    BoundaryFaceSpec& face = model.boundaries[index];
    if (face.flow_kind == BoundaryKind::periodic) {
      continue;
    }
    const bool selected = index == static_cast<std::size_t>(selected_face);
    face.scalars.push_back(ScalarBoundarySpec{
        std::string(stable_name),
        selected ? selected_kind : ScalarBoundaryKind::zero_gradient,
        selected ? selected_value : 0.0,
        ScalarBoundaryKind::zero_gradient,
        0.0});
  }
}

CartesianMeshSpec mesh(Int3 cells) {
  CartesianMeshSpec value;
  value.kind = GeometryKind::uniform;
  value.lower = Real3{0.0, 0.0, 0.0};
  value.upper = Real3{6.0, 6.0, 6.0};
  value.has_exact_cells = true;
  value.exact_cells = cells;
  value.minimum_spacing = Real3{0.5, 0.5, 0.5};
  value.max_growth_ratio = 1.0;
  value.limits.max_global_cells =
      static_cast<std::uint64_t>(cells.x) *
      static_cast<std::uint64_t>(cells.y) *
      static_cast<std::uint64_t>(cells.z);
  value.limits.max_memory_bytes_per_rank = 1U << 20U;
  return value;
}

Status compile(const ValidatedModel& model,
               const CartesianGeometryPlan& geometry,
               const MeshPatch& patch, BoundaryPlan* retained = nullptr,
               BoundaryCompileDiagnostics* diagnostics = nullptr) {
  FieldRegistry registry;
  BoundaryPlan local_boundary;
  SchemePlan schemes;
  TimeSchemePlan time;
  BoundaryPlan& boundary = retained == nullptr ? local_boundary : *retained;
  return BoundaryCompiler::compile(MPI_COMM_SELF, model, geometry, patch,
                                   registry, boundary, schemes, time,
                                   diagnostics);
}

const BoundaryFacePlan* face(const BoundaryPlan& plan,
                             CartesianFace selected) {
  const BoundaryFacePlan* result = nullptr;
  return plan.face(selected, result) ? result : nullptr;
}

const BoundaryIndexSpan* span(const BoundaryPlan& plan,
                              CartesianFace selected,
                              BoundaryStage stage) {
  for (std::size_t index = 0U; index < plan.spans().size; ++index) {
    const BoundaryIndexSpan& candidate = plan.spans().data[index];
    if (candidate.face == selected && candidate.stage == stage) {
      return &candidate;
    }
  }
  return nullptr;
}

bool valid_flow_kind(BoundaryKind kind,
                     const CartesianGeometryPlan& geometry,
                     const MeshPatch& patch) {
  ValidatedModel candidate = open_model();
  BoundaryFaceSpec& face = candidate.boundaries[0U];
  face.flow_kind = kind;
  face.thermal_kind = BoundaryKind::none;
  face.velocity = Real3{2.0, 0.0, 0.0};
  face.direction = Real3{1.0, 0.0, 0.0};
  face.backflow_velocity = Real3{1.0, 0.0, 0.0};
  face.mass_flow_rate = 2.0;
  face.pressure = 101300.0;
  face.temperature = 300.0;
  face.total_pressure = 102000.0;
  face.total_temperature = 305.0;
  face.relaxation = 0.1;
  face.allow_backflow = kind == BoundaryKind::pressure_outlet ||
                        kind == BoundaryKind::nscbc_outlet;
  face.backflow_temperature = 300.0;
  if (kind == BoundaryKind::no_slip_wall) {
    face.velocity = Real3{};
  }
  return static_cast<bool>(compile(candidate, geometry, patch));
}

bool test_flow_and_thermal_matrix(const CartesianGeometryPlan& geometry,
                                  const MeshPatch& patch) {
  bool passed = true;
  const std::array<BoundaryKind, 11U> independent_flow_kinds{
      BoundaryKind::velocity_inlet, BoundaryKind::mass_flow_inlet,
      BoundaryKind::static_state_inlet, BoundaryKind::total_state_inlet,
      BoundaryKind::pressure_outlet, BoundaryKind::nscbc_inlet,
      BoundaryKind::nscbc_outlet, BoundaryKind::no_slip_wall,
      BoundaryKind::moving_wall, BoundaryKind::slip,
      BoundaryKind::symmetry};
  for (const BoundaryKind kind : independent_flow_kinds) {
    passed &= expect(valid_flow_kind(kind, geometry, patch),
                     "every supported non-periodic flow kind compiles");
  }

  ValidatedModel periodic = open_model();
  periodic.boundaries[2U].flow_kind = BoundaryKind::periodic;
  periodic.boundaries[2U].thermal_kind = BoundaryKind::none;
  periodic.boundaries[3U].flow_kind = BoundaryKind::periodic;
  periodic.boundaries[3U].thermal_kind = BoundaryKind::none;
  passed &= expect(static_cast<bool>(compile(periodic, geometry, patch)),
                   "a complete opposite periodic pair compiles");

  const std::array<BoundaryKind, 4U> invalid_flow_kinds{
      BoundaryKind::none, BoundaryKind::adiabatic_wall,
      BoundaryKind::isothermal_wall, BoundaryKind::heat_flux_wall};
  for (const BoundaryKind kind : invalid_flow_kinds) {
    ValidatedModel candidate = open_model();
    candidate.boundaries[2U].flow_kind = kind;
    candidate.boundaries[2U].thermal_kind = BoundaryKind::none;
    passed &= expect(!compile(candidate, geometry, patch),
                     "flow_kind accepts only flow boundary enumerators");
  }

  const std::array<BoundaryKind, 4U> valid_thermal_kinds{
      BoundaryKind::none, BoundaryKind::adiabatic_wall,
      BoundaryKind::isothermal_wall, BoundaryKind::heat_flux_wall};
  for (const BoundaryKind kind : valid_thermal_kinds) {
    ValidatedModel candidate = open_model();
    BoundaryFaceSpec& face = candidate.boundaries[2U];
    face.flow_kind = BoundaryKind::no_slip_wall;
    face.thermal_kind = kind;
    face.temperature = 310.0;
    face.heat_flux = 12.0;
    passed &= expect(static_cast<bool>(compile(candidate, geometry, patch)),
                     "every supported thermal kind compiles on a wall");
  }

  const std::array<BoundaryKind, 12U> invalid_thermal_kinds{
      BoundaryKind::velocity_inlet, BoundaryKind::mass_flow_inlet,
      BoundaryKind::static_state_inlet, BoundaryKind::total_state_inlet,
      BoundaryKind::pressure_outlet, BoundaryKind::nscbc_inlet,
      BoundaryKind::nscbc_outlet, BoundaryKind::no_slip_wall,
      BoundaryKind::moving_wall, BoundaryKind::slip,
      BoundaryKind::symmetry, BoundaryKind::periodic};
  for (const BoundaryKind kind : invalid_thermal_kinds) {
    ValidatedModel candidate = open_model();
    candidate.boundaries[2U].thermal_kind = kind;
    passed &= expect(!compile(candidate, geometry, patch),
                     "thermal_kind accepts only thermal enumerators");
  }

  ValidatedModel thermal_on_inlet = open_model();
  thermal_on_inlet.boundaries[0U].thermal_kind =
      BoundaryKind::isothermal_wall;
  thermal_on_inlet.boundaries[0U].temperature = 300.0;
  passed &= expect(!compile(thermal_on_inlet, geometry, patch),
                   "wall thermal authority is rejected on an inlet");

  ValidatedModel moving = open_model();
  moving.boundaries[2U].flow_kind = BoundaryKind::moving_wall;
  moving.boundaries[2U].velocity = Real3{0.5, 0.0, 0.0};
  passed &= expect(static_cast<bool>(compile(moving, geometry, patch)),
                   "moving wall accepts a nonzero velocity");
  ValidatedModel nonzero_no_slip = open_model();
  nonzero_no_slip.boundaries[2U].velocity = Real3{0.5, 0.0, 0.0};
  passed &= expect(!compile(nonzero_no_slip, geometry, patch),
                   "no-slip wall rejects a nonzero velocity");
  return passed;
}

bool test_closure_scheme_and_reach(const CartesianGeometryPlan& geometry,
                                   const MeshPatch& patch) {
  bool passed = true;
  ValidatedModel half_periodic = open_model();
  half_periodic.boundaries[2U].flow_kind = BoundaryKind::periodic;
  half_periodic.boundaries[2U].thermal_kind = BoundaryKind::none;
  passed &= expect(!compile(half_periodic, geometry, patch),
                   "a half periodic pair is rejected");

  ValidatedModel missing_pressure = open_model();
  missing_pressure.boundaries[1U] = wall();
  passed &= expect(!compile(missing_pressure, geometry, patch),
                   "open pressure reference requires an absolute authority");

  ValidatedModel closed_with_pressure = open_model();
  closed_with_pressure.pressure_reference = PressureReferenceKind::closed_mass;
  passed &= expect(!compile(closed_with_pressure, geometry, patch),
                   "closed mass reference rejects pressure authority");

  ValidatedModel closed;
  closed.fingerprint = 987U;
  closed.pressure_reference = PressureReferenceKind::closed_mass;
  for (BoundaryFaceSpec& face_spec : closed.boundaries) {
    face_spec = wall();
  }
  passed &= expect(static_cast<bool>(compile(closed, geometry, patch)),
                   "closed mass accepts an impermeable wall enclosure");
  closed.boundaries[0U].flow_kind = BoundaryKind::velocity_inlet;
  closed.boundaries[0U].thermal_kind = BoundaryKind::none;
  closed.boundaries[0U].velocity = Real3{1.0, 0.0, 0.0};
  passed &= expect(!compile(closed, geometry, patch),
                   "closed mass rejects a penetrating inlet without pressure authority");

  ValidatedModel supersonic = open_model();
  supersonic.boundaries[1U].mach_limit = 1.0;
  passed &= expect(!compile(supersonic, geometry, patch),
                   "supersonic boundary settings are rejected");

  ValidatedModel limiter = open_model();
  limiter.schemes.limiter = 1.000001;
  passed &= expect(!compile(limiter, geometry, patch),
                   "limiter coefficient above one is rejected");

  CartesianGeometryPlan narrow_geometry;
  MeshPatch narrow_patch;
  const Status narrow_status = CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, mesh(Int3{1, 6, 6}), GeometryBudget{0U, 1U},
      narrow_geometry, narrow_patch);
  passed &= expect(static_cast<bool>(narrow_status),
                   "narrow geometry fixture compiles");
  if (narrow_status) {
    passed &= expect(!compile(open_model(), narrow_geometry, narrow_patch),
                     "local owned extent covers compiled stencil reach");
  }
  return passed;
}

bool test_scalar_parameters(const CartesianGeometryPlan& geometry,
                            const MeshPatch& patch) {
  bool passed = true;
  BoundaryPlan no_scalar_plan;
  const Status no_scalar =
      compile(open_model(), geometry, patch, &no_scalar_plan);
  passed &= expect(static_cast<bool>(no_scalar) &&
                       face(no_scalar_plan, CartesianFace::x_min) != nullptr &&
                       face(no_scalar_plan, CartesianFace::x_min)->scalar_count ==
                           0U,
                   "flow and thermal parameters are excluded from scalar_count");

  const std::array<ScalarBoundaryKind, 4U> kinds{
      ScalarBoundaryKind::dirichlet, ScalarBoundaryKind::normal_flux,
      ScalarBoundaryKind::zero_gradient, ScalarBoundaryKind::convective};
  for (const ScalarBoundaryKind kind : kinds) {
    ValidatedModel candidate = open_model();
    add_scalar_catalog_and_closures(
        candidate, "mixture_fraction", TransportedScalarRole::passive_scalar,
        CartesianFace::x_min, kind, 0.25);
    BoundaryPlan plan;
    const Status status = compile(candidate, geometry, patch, &plan);
    passed &= expect(static_cast<bool>(status),
                     "every scalar value/flux kind compiles");
    if (status) {
      const BoundaryFacePlan* selected =
          face(plan, CartesianFace::x_min);
      passed &= expect(selected != nullptr && selected->scalar_count == 1U &&
                           selected->scalar_begin < plan.parameter_count(),
                       "scalar range contains exactly scalar parameters");
    }
  }
  return passed;
}

bool test_authoritative_patch_and_atomic_outputs(
    const CartesianGeometryPlan& geometry, const MeshPatch& patch) {
  ValidatedModel model = open_model();
  FieldRegistry registry;
  BoundaryPlan boundary;
  SchemePlan schemes;
  TimeSchemePlan time;
  bool passed = expect(static_cast<bool>(BoundaryCompiler::compile(
                           MPI_COMM_SELF, model, geometry, patch, registry,
                           boundary, schemes, time)),
                       "authoritative patch fixture compiles");
  if (!passed) {
    return false;
  }
  const auto unchanged = [&]() {
    return boundary.semantic_fingerprint() != 0U &&
           boundary.revision() == boundary.semantic_fingerprint() &&
           registry.fingerprint() != 0U && schemes.fingerprint() != 0U &&
           time.fingerprint() != 0U;
  };
  const PlanFingerprint old_boundary = boundary.semantic_fingerprint();
  const PlanFingerprint old_layout = boundary.local_layout_fingerprint();
  const PlanFingerprint old_registry = registry.fingerprint();
  const PlanFingerprint old_scheme = schemes.fingerprint();
  const PlanFingerprint old_time = time.fingerprint();
  const auto reject_unchanged = [&](MeshPatch candidate,
                                    std::string_view description) {
    const Status status = BoundaryCompiler::compile(
        MPI_COMM_SELF, model, geometry, candidate, registry, boundary,
        schemes, time);
    return expect(!status && unchanged() &&
                      boundary.semantic_fingerprint() == old_boundary &&
                      boundary.local_layout_fingerprint() == old_layout &&
                      registry.fingerprint() == old_registry &&
                      schemes.fingerprint() == old_scheme &&
                      time.fingerprint() == old_time,
                  description);
  };
  MeshPatch candidate = patch;
  candidate.begin.x += 1;
  passed &= reject_unchanged(candidate,
                             "shifted patch begin is rejected atomically");
  candidate = patch;
  candidate.begin.y = -1;
  passed &= reject_unchanged(candidate,
                             "negative patch begin is rejected atomically");
  candidate = patch;
  candidate.cells.z -= 1;
  passed &= reject_unchanged(candidate,
                             "patch hole is rejected atomically");
  candidate = patch;
  candidate.process_grid.z += 1;
  passed &= reject_unchanged(candidate,
                             "non-authoritative process grid is rejected atomically");
  candidate = patch;
  candidate.process_coord.x += 1;
  passed &= reject_unchanged(candidate,
                             "non-authoritative process coordinate is rejected atomically");
  return passed;
}

bool test_inward_directions(const CartesianGeometryPlan& geometry,
                            const MeshPatch& patch) {
  bool passed = true;
  const std::array<BoundaryKind, 3U> directed{
      BoundaryKind::mass_flow_inlet, BoundaryKind::static_state_inlet,
      BoundaryKind::total_state_inlet};
  for (const BoundaryKind kind : directed) {
    ValidatedModel zero = open_model();
    BoundaryFaceSpec& selected = zero.boundaries[0U];
    selected.flow_kind = kind;
    selected.velocity = {};
    selected.direction = {};
    selected.mass_flow_rate = 1.0;
    selected.pressure = 101325.0;
    selected.temperature = 300.0;
    selected.total_pressure = 102000.0;
    selected.total_temperature = 305.0;
    passed &= expect(!compile(zero, geometry, patch),
                     "directed inlet rejects a zero direction");
    selected.direction = Real3{-1.0, 0.0, 0.0};
    passed &= expect(!compile(zero, geometry, patch),
                     "directed inlet rejects an outward direction");
  }
  ValidatedModel nscbc = open_model();
  nscbc.boundaries[0U].flow_kind = BoundaryKind::nscbc_inlet;
  nscbc.boundaries[0U].velocity = Real3{-1.0, 0.0, 0.0};
  nscbc.boundaries[0U].pressure = 101325.0;
  nscbc.boundaries[0U].temperature = 300.0;
  passed &= expect(!compile(nscbc, geometry, patch),
                   "NSCBC inlet rejects an outward target velocity");

  ValidatedModel backflow = open_model();
  backflow.boundaries[1U].allow_backflow = true;
  backflow.boundaries[1U].backflow_temperature = 300.0;
  backflow.boundaries[1U].backflow_velocity = Real3{1.0, 0.0, 0.0};
  passed &= expect(!compile(backflow, geometry, patch),
                   "declared outlet backflow rejects an outward target velocity");
  backflow.boundaries[1U].backflow_velocity = Real3{-1.0, 0.0, 0.0};
  passed &= expect(static_cast<bool>(compile(backflow, geometry, patch)),
                   "declared outlet backflow accepts an inward target velocity");
  return passed;
}

bool test_resolver_descriptors(const CartesianGeometryPlan& geometry,
                               const MeshPatch& patch) {
  bool passed = true;
  const std::array<BoundaryKind, 4U> inlet_kinds{
      BoundaryKind::velocity_inlet, BoundaryKind::mass_flow_inlet,
      BoundaryKind::static_state_inlet, BoundaryKind::total_state_inlet};
  for (const BoundaryKind kind : inlet_kinds) {
    ValidatedModel candidate = open_model();
    BoundaryFaceSpec& inlet = candidate.boundaries[0U];
    inlet.flow_kind = kind;
    inlet.velocity = {};
    inlet.direction = Real3{1.0, 0.0, 0.0};
    inlet.mass_flow_rate = 1.0;
    inlet.pressure = 101325.0;
    inlet.temperature = 300.0;
    inlet.total_pressure = 102000.0;
    inlet.total_temperature = 305.0;
    BoundaryPlan plan;
    const Status status = compile(candidate, geometry, patch, &plan);
    const BoundaryIndexSpan* momentum =
        span(plan, CartesianFace::x_min, BoundaryStage::momentum);
    const BoundaryIndexSpan* enthalpy =
        span(plan, CartesianFace::x_min, BoundaryStage::enthalpy);
    const BoundaryIndexSpan* pressure =
        span(plan, CartesianFace::x_min, BoundaryStage::pressure);
    const BoundaryValueSource expected_momentum =
        kind == BoundaryKind::velocity_inlet
            ? BoundaryValueSource::compiled_vector
            : BoundaryValueSource::resolved_vector;
    const std::uint32_t expected_momentum_stride =
        kind == BoundaryKind::velocity_inlet ? 0U : 36U;
    passed &= expect(status && momentum != nullptr && enthalpy != nullptr &&
                         momentum->value_source == expected_momentum &&
                         enthalpy->value_source ==
                             BoundaryValueSource::resolved_scalar &&
                         momentum->resolved_stride ==
                             expected_momentum_stride &&
                         enthalpy->resolved_stride == 36U,
                     "all ordinary inlet forms compile a cellwise h resolver slice");
    if (kind == BoundaryKind::static_state_inlet ||
        kind == BoundaryKind::total_state_inlet) {
      passed &= expect(pressure != nullptr &&
                           pressure->value_source ==
                               BoundaryValueSource::resolved_scalar &&
                           pressure->resolved_stride == 36U,
                       "static/total inlet compiles a cellwise pi resolver slice");
    }
  }

  ValidatedModel backflow = open_model();
  backflow.boundaries[1U].allow_backflow = true;
  backflow.boundaries[1U].backflow_velocity = Real3{-1.0, 0.0, 0.0};
  backflow.boundaries[1U].backflow_temperature = 300.0;
  BoundaryPlan plan;
  const Status status = compile(backflow, geometry, patch, &plan);
  const BoundaryIndexSpan* momentum =
      span(plan, CartesianFace::x_max, BoundaryStage::momentum);
  const BoundaryIndexSpan* enthalpy =
      span(plan, CartesianFace::x_max, BoundaryStage::enthalpy);
  passed &= expect(status && momentum != nullptr && enthalpy != nullptr &&
                       momentum->value_source ==
                           BoundaryValueSource::resolved_vector &&
                       enthalpy->value_source ==
                           BoundaryValueSource::resolved_scalar,
                   "pressure outlet backflow compiles conditional U+h resolver slices");
  return passed;
}

bool test_inlet_temperature_and_species_composition(
    const CartesianGeometryPlan& geometry, const MeshPatch& patch) {
  bool passed = true;

  ValidatedModel cold_velocity = open_model();
  cold_velocity.boundaries[0U].temperature = 0.0;
  passed &= expect(!compile(cold_velocity, geometry, patch),
                   "velocity inlet rejects a missing absolute temperature");

  ValidatedModel composition = open_model();
  add_scalar_catalog_and_closures(
      composition, "O2", TransportedScalarRole::species,
      CartesianFace::x_min, ScalarBoundaryKind::dirichlet, 0.2);
  add_scalar_catalog_and_closures(
      composition, "N2", TransportedScalarRole::species,
      CartesianFace::x_min, ScalarBoundaryKind::dirichlet, 0.7);
  add_scalar_catalog_and_closures(
      composition, "tracer", TransportedScalarRole::passive_scalar,
      CartesianFace::x_min, ScalarBoundaryKind::dirichlet, -4.0);
  passed &= expect(static_cast<bool>(compile(composition, geometry, patch)),
                   "complete bounded inlet species and unrestricted passive scalar compile");

  composition.boundaries[0U].scalars[0U].value = -0.01;
  passed &= expect(!compile(composition, geometry, patch),
                   "inlet species rejects a negative Dirichlet fraction");
  composition.boundaries[0U].scalars[0U].value = 0.4;
  composition.boundaries[0U].scalars[1U].value = 0.7;
  passed &= expect(!compile(composition, geometry, patch),
                   "independent inlet species fractions sum to at most one");
  composition.boundaries[0U].scalars[0U].value = 0.2;
  composition.boundaries[0U].scalars[1U].value = 0.7;

  const std::array<BoundaryKind, 5U> inlet_kinds{
      BoundaryKind::velocity_inlet, BoundaryKind::mass_flow_inlet,
      BoundaryKind::static_state_inlet, BoundaryKind::total_state_inlet,
      BoundaryKind::nscbc_inlet};
  for (const BoundaryKind kind : inlet_kinds) {
    ValidatedModel incomplete = composition;
    BoundaryFaceSpec& inlet = incomplete.boundaries[0U];
    inlet.flow_kind = kind;
    inlet.velocity = Real3{2.0, 0.0, 0.0};
    inlet.direction = Real3{1.0, 0.0, 0.0};
    inlet.mass_flow_rate = 1.0;
    inlet.pressure = 101325.0;
    inlet.temperature = 300.0;
    inlet.total_pressure = 102000.0;
    inlet.total_temperature = 305.0;
    passed &= expect(static_cast<bool>(compile(incomplete, geometry, patch)),
                     "every inlet form accepts a complete species Dirichlet composition");
    inlet.scalars[1U].kind = ScalarBoundaryKind::zero_gradient;
    passed &= expect(!compile(incomplete, geometry, patch),
                     "every inlet form requires a complete species Dirichlet composition");
  }

  ValidatedModel wall_mixed = composition;
  wall_mixed.boundaries[2U].scalars[0U].kind =
      ScalarBoundaryKind::dirichlet;
  wall_mixed.boundaries[2U].scalars[0U].value = 0.8;
  wall_mixed.boundaries[2U].scalars[1U].kind =
      ScalarBoundaryKind::zero_gradient;
  wall_mixed.boundaries[2U].scalars[1U].value = 99.0;
  passed &= expect(static_cast<bool>(compile(wall_mixed, geometry, patch)),
                   "wall composition sums only species with Dirichlet authority");
  wall_mixed.boundaries[2U].scalars[0U].value = 1.01;
  passed &= expect(!compile(wall_mixed, geometry, patch),
                   "species Dirichlet authority is bounded on every face");
  return passed;
}

bool test_backflow_scalar_policy(const CartesianGeometryPlan& geometry,
                                 const MeshPatch& patch) {
  ValidatedModel candidate = open_model();
  add_scalar_catalog_and_closures(
      candidate, "O2", TransportedScalarRole::species,
      CartesianFace::x_min, ScalarBoundaryKind::dirichlet, 0.2);
  add_scalar_catalog_and_closures(
      candidate, "N2", TransportedScalarRole::species,
      CartesianFace::x_min, ScalarBoundaryKind::dirichlet, 0.7);
  add_scalar_catalog_and_closures(
      candidate, "mixture_fraction", TransportedScalarRole::passive_scalar,
      CartesianFace::x_min, ScalarBoundaryKind::dirichlet, 0.1);
  BoundaryFaceSpec& outlet = candidate.boundaries[1U];
  outlet.allow_backflow = true;
  outlet.backflow_velocity = Real3{-1.0, 0.0, 0.0};
  outlet.backflow_temperature = 300.0;

  bool passed = expect(!compile(candidate, geometry, patch),
                       "backflow outlet rejects missing scalar closures");
  outlet.scalars[0U].backflow_kind = ScalarBoundaryKind::dirichlet;
  outlet.scalars[0U].backflow_value = 0.2;
  outlet.scalars[1U].backflow_kind = ScalarBoundaryKind::dirichlet;
  outlet.scalars[1U].backflow_value = 0.7;
  outlet.scalars[2U].backflow_kind = ScalarBoundaryKind::dirichlet;
  outlet.scalars[2U].backflow_value = -3.0;
  BoundaryPlan plan;
  const Status valid = compile(candidate, geometry, patch, &plan);
  std::size_t resolved_scalars = 0U;
  for (std::size_t index = 0U; index < plan.spans().size; ++index) {
    const BoundaryIndexSpan& value = plan.spans().data[index];
    if (value.face == CartesianFace::x_max &&
        value.stage == BoundaryStage::scalar &&
        value.value_source == BoundaryValueSource::resolved_scalar) {
      ++resolved_scalars;
    }
  }
  passed &= expect(valid && resolved_scalars == 3U &&
                       plan.scalar_backflow_targets().size ==
                           plan.parameter_count(),
                   "backflow scalars compile cellwise conditional resolver slices");

  candidate.boundaries[1U].scalars[0U].backflow_value = 1.1;
  passed &= expect(!compile(candidate, geometry, patch),
                   "species backflow value is bounded by one");
  candidate.boundaries[1U].scalars[0U].backflow_value = 0.4;
  candidate.boundaries[1U].scalars[1U].backflow_value = 0.7;
  passed &= expect(!compile(candidate, geometry, patch),
                   "independent species backflow sum is bounded by one");
  return passed;
}

bool test_pressure_reference_and_diagnostics(
    const CartesianGeometryPlan& geometry, const MeshPatch& patch) {
  BoundaryPlan open;
  BoundaryCompileDiagnostics diagnostics{17};
  bool passed = expect(compile(open_model(), geometry, patch, &open,
                               &diagnostics) &&
                           open.pressure_reference() ==
                               PressureReferenceKind::boundary_absolute &&
                           diagnostics.lowest_failing_rank == -1,
                       "compiled plan exposes pressure reference and success diagnostics");

  ValidatedModel closed;
  closed.fingerprint = 987U;
  closed.pressure_reference = PressureReferenceKind::closed_mass;
  for (BoundaryFaceSpec& value : closed.boundaries) {
    value = wall();
  }
  BoundaryPlan closed_plan;
  passed &= expect(compile(closed, geometry, patch, &closed_plan) &&
                       closed_plan.pressure_reference() ==
                           PressureReferenceKind::closed_mass,
                   "closed-mass pressure policy is retained by the plan");

  ValidatedModel invalid = open_model();
  invalid.schemes.limiter = 2.0;
  diagnostics.lowest_failing_rank = -1;
  passed &= expect(!compile(invalid, geometry, patch, &open, &diagnostics) &&
                       diagnostics.lowest_failing_rank == 0,
                   "collective compile failure reports the lowest failing rank");
  return passed;
}

bool test_direct_scalar_catalog_limits(const CartesianGeometryPlan& geometry,
                                       const MeshPatch& patch) {
  bool passed = true;
  ValidatedModel invalid_name = open_model();
  add_scalar_catalog_and_closures(
      invalid_name, "bad/name", TransportedScalarRole::passive_scalar,
      CartesianFace::x_min, ScalarBoundaryKind::dirichlet, 0.0);
  passed &= expect(!compile(invalid_name, geometry, patch),
                   "direct model rejects an invalid scalar stable name");

  ValidatedModel too_many = open_model();
  for (std::size_t index = 0U; index < 65U; ++index) {
    add_scalar_catalog_and_closures(
        too_many, "s" + std::to_string(index),
        TransportedScalarRole::passive_scalar, CartesianFace::x_min,
        ScalarBoundaryKind::zero_gradient, 0.0);
  }
  passed &= expect(!compile(too_many, geometry, patch),
                   "direct model rejects more than 64 transported scalars");
  return passed;
}

bool test_safe_face_lookup(const CartesianGeometryPlan& geometry,
                           const MeshPatch& patch) {
  BoundaryPlan plan;
  bool passed = expect(static_cast<bool>(compile(open_model(), geometry, patch,
                                                 &plan)),
                       "safe face lookup fixture compiles");
  const BoundaryFacePlan* selected = nullptr;
  passed &= expect(static_cast<bool>(plan.face(CartesianFace::x_min,
                                               selected)) &&
                       selected != nullptr,
                   "valid face lookup returns a descriptor");
  selected = reinterpret_cast<const BoundaryFacePlan*>(1U);
  const Status invalid =
      plan.face(static_cast<CartesianFace>(255U), selected);
  passed &= expect(invalid.code == StatusCode::invalid_plan &&
                       selected == nullptr,
                   "invalid face enum returns invalid_plan and null output");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool passed = expect(
      static_cast<bool>(CartesianGeometryCompiler::compile(
          MPI_COMM_SELF, mesh(Int3{6, 6, 6}), GeometryBudget{0U, 1U},
          geometry, patch)),
      "boundary compiler geometry fixture compiles");
  if (passed) {
    passed &= test_flow_and_thermal_matrix(geometry, patch);
    passed &= test_closure_scheme_and_reach(geometry, patch);
    passed &= test_scalar_parameters(geometry, patch);
    passed &= test_authoritative_patch_and_atomic_outputs(geometry, patch);
    passed &= test_inward_directions(geometry, patch);
    passed &= test_resolver_descriptors(geometry, patch);
    passed &= test_inlet_temperature_and_species_composition(geometry, patch);
    passed &= test_backflow_scalar_policy(geometry, patch);
    passed &= test_pressure_reference_and_diagnostics(geometry, patch);
    passed &= test_direct_scalar_catalog_limits(geometry, patch);
    passed &= test_safe_face_lookup(geometry, patch);
  }
  if (passed) {
    std::cout << "v0.4 boundary compile tests passed\n";
  }
  MPI_Finalize();
  return passed ? 0 : 1;
}

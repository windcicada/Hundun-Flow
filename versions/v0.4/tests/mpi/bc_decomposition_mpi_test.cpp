// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_boundary.hpp"

#include <mpi.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <string_view>

namespace {

using namespace hundun::v04;

bool expect(bool condition, int rank, std::string_view description) {
  if (!condition) {
    std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  }
  return condition;
}

std::uint64_t mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= 1099511628211ULL;
  return hash;
}

std::uint64_t bits(double value) noexcept {
  std::uint64_t result = 0U;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

BoundaryFaceSpec wall() {
  BoundaryFaceSpec value;
  value.flow_kind = BoundaryKind::no_slip_wall;
  value.thermal_kind = BoundaryKind::adiabatic_wall;
  value.mach_limit = 0.95;
  return value;
}

ValidatedModel model(bool periodic_x) {
  ValidatedModel value;
  value.fingerprint = 1234567U;
  value.pressure_reference = PressureReferenceKind::boundary_absolute;
  for (BoundaryFaceSpec& face : value.boundaries) {
    face = wall();
  }
  value.boundaries[0].flow_kind = periodic_x ? BoundaryKind::periodic
                                              : BoundaryKind::velocity_inlet;
  value.boundaries[0].thermal_kind = BoundaryKind::none;
  value.boundaries[0].velocity = Real3{2.0, 0.0, 0.0};
  value.boundaries[0].temperature = 300.0;
  value.boundaries[1].flow_kind = periodic_x ? BoundaryKind::periodic
                                              : BoundaryKind::pressure_outlet;
  value.boundaries[1].thermal_kind = BoundaryKind::none;
  value.boundaries[1].pressure = 101325.0;
  if (periodic_x) {
    value.pressure_reference = PressureReferenceKind::closed_mass;
  }
  value.schemes = SchemeSpec{};
  value.time = TimeControlSpec{};
  return value;
}

ValidatedModel semantic_model() {
  ValidatedModel value = model(false);

  value.transported_scalars.push_back(TransportedScalarSpec{
      "mixture_fraction", TransportedScalarRole::passive_scalar});
  for (BoundaryFaceSpec& face : value.boundaries) {
    face.scalars.push_back(ScalarBoundarySpec{
        "mixture_fraction", ScalarBoundaryKind::zero_gradient, 0.0,
        ScalarBoundaryKind::zero_gradient, 0.0});
  }

  BoundaryFaceSpec& mass_inlet = value.boundaries[0U];
  mass_inlet.flow_kind = BoundaryKind::mass_flow_inlet;
  mass_inlet.velocity = Real3{};
  mass_inlet.direction = Real3{1.0, 0.0, 0.0};
  mass_inlet.mass_flow_rate = 2.0;
  mass_inlet.temperature = 300.0;
  mass_inlet.scalars[0U] = ScalarBoundarySpec{
      "mixture_fraction", ScalarBoundaryKind::dirichlet, 0.25,
      ScalarBoundaryKind::zero_gradient, 0.0};

  BoundaryFaceSpec& total_inlet = value.boundaries[2U];
  total_inlet.flow_kind = BoundaryKind::total_state_inlet;
  total_inlet.thermal_kind = BoundaryKind::none;
  total_inlet.direction = Real3{0.0, 1.0, 0.0};
  total_inlet.total_pressure = 102000.0;
  total_inlet.total_temperature = 305.0;

  BoundaryFaceSpec& nscbc_outlet = value.boundaries[3U];
  nscbc_outlet.flow_kind = BoundaryKind::nscbc_outlet;
  nscbc_outlet.thermal_kind = BoundaryKind::none;
  nscbc_outlet.pressure = 101325.0;
  nscbc_outlet.temperature = 300.0;
  nscbc_outlet.backflow_velocity = Real3{0.1, -0.2, 0.3};
  nscbc_outlet.backflow_temperature = 301.0;
  nscbc_outlet.relaxation = 0.1;
  nscbc_outlet.mach_limit = 0.90;
  nscbc_outlet.allow_backflow = true;
  nscbc_outlet.scalars[0U].backflow_kind = ScalarBoundaryKind::dirichlet;
  nscbc_outlet.scalars[0U].backflow_value = 0.1;

  BoundaryFaceSpec& moving = value.boundaries[4U];
  moving.flow_kind = BoundaryKind::moving_wall;
  moving.thermal_kind = BoundaryKind::heat_flux_wall;
  moving.velocity = Real3{0.5, 0.0, 0.0};
  moving.heat_flux = 12.0;
  return value;
}

CartesianMeshSpec mesh(int size, double upper_x = 4.0) {
  CartesianMeshSpec value;
  value.kind = GeometryKind::uniform;
  value.lower = Real3{0.0, 0.0, 0.0};
  value.upper = Real3{upper_x, 2.0, 2.0};
  value.has_exact_cells = true;
  value.exact_cells = Int3{4 * size, 4, 4};
  value.minimum_spacing = Real3{0.25, 0.25, 0.25};
  value.max_growth_ratio = 1.0;
  value.limits.max_global_cells =
      static_cast<std::uint64_t>(value.exact_cells.x) * 16U;
  value.limits.max_memory_bytes_per_rank = 1U << 20U;
  return value;
}

bool all_true(bool local) {
  const int input = local ? 1 : 0;
  int output = 0;
  return MPI_Allreduce(&input, &output, 1, MPI_INT, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS && output != 0;
}

std::uint64_t boundary_observable(const BoundaryPlan& boundary) noexcept {
  std::uint64_t hash = 1469598103934665603ULL;
  hash = mix(hash, boundary.semantic_fingerprint());
  hash = mix(hash, boundary.local_layout_fingerprint());
  hash = mix(hash, boundary.revision());
  hash = mix(hash, boundary.parameter_count());
  hash = mix(hash, boundary.required_ghost_width());
  hash = mix(hash, boundary.velocity_field());
  hash = mix(hash, boundary.pressure_field());
  hash = mix(hash, boundary.enthalpy_field());
  hash = mix(hash, static_cast<std::uint8_t>(boundary.pressure_reference()));
  const std::array<CartesianFace, 6U> faces{
      CartesianFace::x_min, CartesianFace::x_max, CartesianFace::y_min,
      CartesianFace::y_max, CartesianFace::z_min, CartesianFace::z_max};
  for (const CartesianFace selected : faces) {
    const BoundaryFacePlan* face = nullptr;
    if (!boundary.face(selected, face) || face == nullptr) {
      return 0U;
    }
    hash = mix(hash, static_cast<std::uint8_t>(face->flow_kind));
    hash = mix(hash, static_cast<std::uint8_t>(face->thermal_kind));
    hash = mix(hash, face->local_owner ? 1U : 0U);
    hash = mix(hash, face->periodic ? 1U : 0U);
    hash = mix(hash, face->flow_parameter);
    hash = mix(hash, face->thermal_parameter);
    hash = mix(hash, face->scalar_begin);
    hash = mix(hash, face->scalar_count);
  }
  for (std::size_t index = 0U; index < boundary.spans().size; ++index) {
    const BoundaryIndexSpan& span = boundary.spans().data[index];
    hash = mix(hash, static_cast<std::uint8_t>(span.stage));
    hash = mix(hash, static_cast<std::uint8_t>(span.relation));
    hash = mix(hash, static_cast<std::uint8_t>(span.face));
    hash = mix(hash, span.field);
    hash = mix(hash, span.component_begin);
    hash = mix(hash, span.component_count);
    hash = mix(hash, span.ghost_layers);
    hash = mix(hash, span.tangent_inner_count);
    hash = mix(hash, span.tangent_outer_count);
    hash = mix(hash, span.parameter);
    hash = mix(hash, span.resolved_begin);
    hash = mix(hash, static_cast<std::uint8_t>(span.value_source));
  }
  const std::array<Span<const double>, 22U> arrays{
      boundary.velocity_x(),          boundary.velocity_y(),
      boundary.velocity_z(),          boundary.direction_x(),
      boundary.direction_y(),         boundary.direction_z(),
      boundary.backflow_velocity_x(), boundary.backflow_velocity_y(),
      boundary.backflow_velocity_z(), boundary.scalar_targets(),
      boundary.scalar_backflow_targets(),
      boundary.pressure_targets(),    boundary.temperature_targets(),
      boundary.total_pressure_targets(),
      boundary.total_temperature_targets(),
      boundary.backflow_temperature_targets(),
      boundary.heat_flux_targets(), boundary.mass_flow_targets(),
      boundary.relaxation_rates(), boundary.mach_limits(),
      boundary.normal_distance_1(), boundary.normal_distance_2()};
  for (const Span<const double> values : arrays) {
    hash = mix(hash, values.size);
    for (std::size_t index = 0U; index < values.size; ++index) {
      hash = mix(hash, bits(values.data[index]));
    }
  }
  for (std::size_t index = 0U; index < boundary.allow_backflow().size;
       ++index) {
    hash = mix(hash, boundary.allow_backflow().data[index]);
    hash = mix(hash, static_cast<std::uint8_t>(
                         boundary.parameter_roles().data[index]));
    hash = mix(hash,
               static_cast<std::uint8_t>(boundary.scalar_kinds().data[index]));
  }
  const Int3 local = boundary.local_cells();
  hash = mix(hash, static_cast<std::uint32_t>(local.x));
  hash = mix(hash, static_cast<std::uint32_t>(local.y));
  hash = mix(hash, static_cast<std::uint32_t>(local.z));
  hash = mix(hash, boundary.resolved_scalar_count());
  hash = mix(hash, boundary.resolved_vector_count());
  hash = mix(hash, boundary.resolved_normal_gradient_count());
  return hash;
}

std::uint64_t scheme_observable(const SchemePlan& plan) noexcept {
  std::uint64_t hash = plan.fingerprint();
  hash = mix(hash, static_cast<std::uint8_t>(plan.momentum()));
  hash = mix(hash, static_cast<std::uint8_t>(plan.enthalpy()));
  hash = mix(hash, static_cast<std::uint8_t>(plan.species()));
  hash = mix(hash, static_cast<std::uint8_t>(plan.passive_scalar()));
  hash = mix(hash, static_cast<std::uint8_t>(plan.diffusion()));
  hash = mix(hash, bits(plan.limiter()));
  hash = mix(hash, plan.required_ghost_width());
  return hash;
}

std::uint64_t registry_observable(const FieldRegistry& registry,
                                  bool& ok) {
  FieldRegistry copy = registry;
  FieldSchema schema;
  const Status status = copy.freeze(schema);
  ok = static_cast<bool>(status);
  std::uint64_t hash = 1469598103934665603ULL;
  if (!status) {
    return hash;
  }
  hash = mix(hash, schema.size());
  for (const FieldDescriptor& field : schema) {
    hash = mix(hash, field.id);
    hash = mix(hash, field.components);
    hash = mix(hash, field.ghost_width);
    for (const unsigned char character : field.stable_name) {
      hash = mix(hash, character);
    }
  }
  return hash;
}

struct CompileState {
  FieldRegistry registry;
  BoundaryPlan boundary;
  SchemePlan schemes;
  TimeSchemePlan time;
};

Status compile(const ValidatedModel& input,
               const CartesianGeometryPlan& geometry,
               const MeshPatch& patch, CompileState& state,
               BoundaryCompileDiagnostics* diagnostics = nullptr) {
  return BoundaryCompiler::compile(MPI_COMM_WORLD, input, geometry, patch,
                                   state.registry, state.boundary,
                                   state.schemes, state.time, diagnostics);
}

bool same_collective_status(Status status, int rank) {
  const std::uint64_t packed =
      (static_cast<std::uint64_t>(status.code) << 32U) | status.detail;
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  const bool reduced =
      MPI_Allreduce(&packed, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                    MPI_COMM_WORLD) == MPI_SUCCESS &&
      MPI_Allreduce(&packed, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                    MPI_COMM_WORLD) == MPI_SUCCESS;
  return expect(reduced && !status && minimum == maximum, rank,
                "rank-divergent semantics are collectively rejected");
}

enum class Mutation : std::uint8_t {
  velocity,
  direction,
  backflow_velocity,
  mass_flow,
  pressure,
  temperature,
  total_pressure,
  total_temperature,
  heat_flux,
  relaxation,
  mach_limit,
  allow_backflow,
  scalar_kind,
  scalar_value,
  scalar_backflow_kind,
  scalar_backflow_value,
  scheme,
  time
};

void mutate(ValidatedModel& value, Mutation selected) {
  switch (selected) {
    case Mutation::velocity:
      value.boundaries[4U].velocity.x = 0.6;
      break;
    case Mutation::direction:
      value.boundaries[0U].direction = Real3{0.9, 0.1, 0.0};
      break;
    case Mutation::backflow_velocity:
      value.boundaries[3U].backflow_velocity.x = 0.15;
      break;
    case Mutation::mass_flow:
      value.boundaries[0U].mass_flow_rate = 2.1;
      break;
    case Mutation::pressure:
      value.boundaries[3U].pressure = 101300.0;
      break;
    case Mutation::temperature:
      value.boundaries[3U].temperature = 302.0;
      break;
    case Mutation::total_pressure:
      value.boundaries[2U].total_pressure = 102100.0;
      break;
    case Mutation::total_temperature:
      value.boundaries[2U].total_temperature = 306.0;
      break;
    case Mutation::heat_flux:
      value.boundaries[4U].heat_flux = 13.0;
      break;
    case Mutation::relaxation:
      value.boundaries[3U].relaxation = 0.2;
      break;
    case Mutation::mach_limit:
      value.boundaries[3U].mach_limit = 0.85;
      break;
    case Mutation::allow_backflow:
      value.boundaries[3U].allow_backflow = false;
      break;
    case Mutation::scalar_kind:
      value.boundaries[0U].scalars[0U].kind =
          ScalarBoundaryKind::zero_gradient;
      break;
    case Mutation::scalar_value:
      value.boundaries[0U].scalars[0U].value = 0.3;
      break;
    case Mutation::scalar_backflow_kind:
      value.boundaries[3U].scalars[0U].backflow_kind =
          ScalarBoundaryKind::zero_gradient;
      break;
    case Mutation::scalar_backflow_value:
      value.boundaries[3U].scalars[0U].backflow_value = 0.2;
      break;
    case Mutation::scheme:
      value.schemes.momentum = ConvectionScheme::central2;
      break;
    case Mutation::time:
      value.time.initial_dt = 2.0e-4;
      break;
  }
}

bool test_divergent_model(Mutation selected, std::string_view description,
                          const CartesianGeometryPlan& geometry,
                          const MeshPatch& patch, int rank, int size) {
  if (size == 1) {
    return true;
  }
  CompileState state;
  const ValidatedModel baseline = semantic_model();
  bool passed = expect(static_cast<bool>(compile(baseline, geometry, patch,
                                                  state)),
                       rank, "semantic baseline compiles");
  const std::uint64_t old_boundary = boundary_observable(state.boundary);
  const std::uint64_t old_schemes = scheme_observable(state.schemes);
  const std::uint64_t old_time = state.time.fingerprint();
  bool registry_ok = false;
  const std::uint64_t old_registry =
      registry_observable(state.registry, registry_ok);
  passed &= expect(registry_ok, rank, "registry baseline is observable");

  ValidatedModel divergent = baseline;
  if (rank == 1) {
    mutate(divergent, selected);
  }
  BoundaryCompileDiagnostics diagnostics;
  const Status status = compile(divergent, geometry, patch, state,
                                &diagnostics);
  passed &= same_collective_status(status, rank);
  passed &= expect(diagnostics.lowest_failing_rank == 1, rank,
                   "semantic divergence reports rank one as first failure");
  bool current_registry_ok = false;
  const std::uint64_t current_registry =
      registry_observable(state.registry, current_registry_ok);
  passed &= expect(boundary_observable(state.boundary) == old_boundary &&
                       scheme_observable(state.schemes) == old_schemes &&
                       state.time.fingerprint() == old_time &&
                       current_registry_ok && current_registry == old_registry,
                   rank, description);
  return passed;
}

bool test_registry_tail_divergence(const CartesianGeometryPlan& geometry,
                                   const MeshPatch& patch, int rank,
                                   int size) {
  if (size == 1) {
    return true;
  }
  CompileState state;
  const ValidatedModel baseline = semantic_model();
  bool passed = expect(static_cast<bool>(compile(baseline, geometry, patch,
                                                  state)),
                       rank, "registry-divergence baseline compiles");
  if (rank == 1) {
    FieldId unused = 0U;
    passed &= expect(static_cast<bool>(state.registry.declare_field(
                         "rank_one_tail", 1U, 2U, unused)),
                     rank, "rank-one registry tail is introduced");
  }
  const std::uint64_t old_boundary = boundary_observable(state.boundary);
  const std::uint64_t old_schemes = scheme_observable(state.schemes);
  const std::uint64_t old_time = state.time.fingerprint();
  bool registry_ok = false;
  const std::uint64_t old_registry =
      registry_observable(state.registry, registry_ok);
  passed &= expect(registry_ok, rank, "divergent registry is observable");

  const Status status = compile(baseline, geometry, patch, state);
  passed &= same_collective_status(status, rank);
  bool current_registry_ok = false;
  const std::uint64_t current_registry =
      registry_observable(state.registry, current_registry_ok);
  passed &= expect(boundary_observable(state.boundary) == old_boundary &&
                       scheme_observable(state.schemes) == old_schemes &&
                       state.time.fingerprint() == old_time &&
                       current_registry_ok && current_registry == old_registry,
                   rank,
                   "registry-tail rejection preserves all prior outputs");
  return passed;
}

bool test_geometry_identity(const CartesianGeometryPlan& baseline_geometry,
                            const MeshPatch& patch, int rank, int size) {
  if (size == 1) {
    return true;
  }
  CompileState state;
  const ValidatedModel baseline = semantic_model();
  bool passed = expect(
      static_cast<bool>(compile(baseline, baseline_geometry, patch, state)),
      rank, "geometry-divergence baseline compiles");
  const std::uint64_t old_boundary = boundary_observable(state.boundary);
  const std::uint64_t old_schemes = scheme_observable(state.schemes);
  const std::uint64_t old_time = state.time.fingerprint();
  bool registry_ok = false;
  const std::uint64_t old_registry =
      registry_observable(state.registry, registry_ok);
  passed &= expect(registry_ok, rank, "geometry baseline registry is observed");

  CartesianGeometryPlan divergent_geometry;
  MeshPatch unused_self_patch;
  const double upper_x = rank == 1 ? 4.25 : 4.0;
  const Status geometry_status = CartesianGeometryCompiler::compile(
      MPI_COMM_SELF, mesh(size, upper_x), GeometryBudget{0U, 1U},
      divergent_geometry, unused_self_patch);
  passed &= expect(static_cast<bool>(geometry_status), rank,
                   "rank-local divergent geometry compiles");
  const Status status = compile(baseline, divergent_geometry, patch, state);
  passed &= same_collective_status(status, rank);
  bool current_registry_ok = false;
  const std::uint64_t current_registry =
      registry_observable(state.registry, current_registry_ok);
  passed &= expect(boundary_observable(state.boundary) == old_boundary &&
                       scheme_observable(state.schemes) == old_schemes &&
                       state.time.fingerprint() == old_time &&
                       current_registry_ok && current_registry == old_registry,
                   rank,
                   "geometry-identity rejection preserves all prior outputs");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  CartesianGeometryPlan geometry;
  MeshPatch patch;
  bool passed = expect(
      static_cast<bool>(CartesianGeometryCompiler::compile(
          MPI_COMM_WORLD, mesh(size), GeometryBudget{0U, 1U}, geometry,
          patch)),
      rank, "geometry compiles for boundary decomposition");

  FieldRegistry registry;
  BoundaryPlan boundary;
  SchemePlan schemes;
  TimeSchemePlan time;
  ValidatedModel open = model(false);
  const Status status = BoundaryCompiler::compile(
      MPI_COMM_WORLD, open, geometry, patch, registry, boundary, schemes, time);
  passed &= expect(static_cast<bool>(status), rank,
                   "open boundary plan compiles collectively");
  const BoundaryFacePlan* x_min = nullptr;
  const BoundaryFacePlan* x_max = nullptr;
  passed &= expect(boundary.face(CartesianFace::x_min, x_min) &&
                       boundary.face(CartesianFace::x_max, x_max) &&
                       x_min != nullptr && x_max != nullptr &&
                       x_min->local_owner == (rank == 0) &&
                       x_max->local_owner == (rank + 1 == size),
                   rank, "only global edge ranks own x physical faces");
  std::uint64_t semantic = boundary.semantic_fingerprint();
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  MPI_Allreduce(&semantic, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(&semantic, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                MPI_COMM_WORLD);
  passed &= expect(semantic != 0U && minimum == maximum, rank,
                   "semantic fingerprint is decomposition independent");
  const std::uint64_t retained = boundary_observable(boundary);

  ValidatedModel periodic = model(true);
  FieldRegistry periodic_registry;
  BoundaryPlan periodic_boundary;
  SchemePlan periodic_schemes;
  TimeSchemePlan periodic_time;
  const Status periodic_status = BoundaryCompiler::compile(
      MPI_COMM_WORLD, periodic, geometry, patch, periodic_registry,
      periodic_boundary, periodic_schemes, periodic_time);
  bool has_periodic_physical_writer = false;
  for (std::size_t index = 0U; index < periodic_boundary.spans().size;
       ++index) {
    const CartesianFace face = periodic_boundary.spans().data[index].face;
    has_periodic_physical_writer |= face == CartesianFace::x_min ||
                                    face == CartesianFace::x_max;
  }
  passed &= expect(static_cast<bool>(periodic_status) &&
                       periodic_boundary.halo_topology().periodic_x &&
                       !has_periodic_physical_writer,
                   rank, "periodic pair becomes halo topology with no writer");

  ValidatedModel invalid = open;
  if (rank == 0) {
    invalid.schemes.limiter = std::numeric_limits<double>::quiet_NaN();
  }
  BoundaryCompileDiagnostics diagnostics;
  const Status rejected = BoundaryCompiler::compile(
      MPI_COMM_WORLD, invalid, geometry, patch, registry, boundary, schemes,
      time, &diagnostics);
  const std::uint64_t packed =
      (static_cast<std::uint64_t>(rejected.code) << 32U) | rejected.detail;
  MPI_Allreduce(&packed, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                MPI_COMM_WORLD);
  MPI_Allreduce(&packed, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                MPI_COMM_WORLD);
  passed &= expect(rejected.code == StatusCode::invalid_plan &&
                       minimum == maximum &&
                       diagnostics.lowest_failing_rank == 0 &&
                       boundary_observable(boundary) == retained,
                   rank, "collective invalid input preserves prior plan");

  const std::array<std::pair<Mutation, std::string_view>, 18U> mutations{{
      {Mutation::velocity, "velocity divergence preserves old outputs"},
      {Mutation::direction, "direction divergence preserves old outputs"},
      {Mutation::backflow_velocity,
       "backflow velocity divergence preserves old outputs"},
      {Mutation::mass_flow, "mass-flow divergence preserves old outputs"},
      {Mutation::pressure, "pressure divergence preserves old outputs"},
      {Mutation::temperature,
       "temperature divergence preserves old outputs"},
      {Mutation::total_pressure,
       "total-pressure divergence preserves old outputs"},
      {Mutation::total_temperature,
       "total-temperature divergence preserves old outputs"},
      {Mutation::heat_flux, "heat-flux divergence preserves old outputs"},
      {Mutation::relaxation,
       "relaxation divergence preserves old outputs"},
      {Mutation::mach_limit, "Mach-limit divergence preserves old outputs"},
      {Mutation::allow_backflow,
       "backflow switch divergence preserves old outputs"},
      {Mutation::scalar_kind,
       "scalar-kind divergence preserves old outputs"},
      {Mutation::scalar_value,
       "scalar-value divergence preserves old outputs"},
      {Mutation::scalar_backflow_kind,
       "scalar-backflow-kind divergence preserves old outputs"},
      {Mutation::scalar_backflow_value,
       "scalar-backflow-value divergence preserves old outputs"},
      {Mutation::scheme, "scheme divergence preserves old outputs"},
      {Mutation::time, "time-plan divergence preserves old outputs"}}};
  for (const auto& [selected, description] : mutations) {
    passed &= test_divergent_model(selected, description, geometry, patch,
                                   rank, size);
  }
  passed &= test_registry_tail_divergence(geometry, patch, rank, size);
  passed &= test_geometry_identity(geometry, patch, rank, size);

  passed = all_true(passed);
  if (rank == 0 && passed) {
    std::cout << "v0.4 boundary decomposition tests passed\n";
  }
  MPI_Finalize();
  return passed ? 0 : 1;
}

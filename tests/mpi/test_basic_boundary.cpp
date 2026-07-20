// SPDX-License-Identifier: Apache-2.0

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/finite_volume/poisson_boundary_adapter.hpp"

#include "boundary/src/basic_boundary_test_seam.hpp"
#include "hundun/config/resolved_case.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/field_access_plan.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/field_storage.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "runtime/src/mpi_context_test_seam.hpp"
#include "runtime/src/mpi_error.hpp"
#include "tests/support/allocation_attempt_guard.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using hundun::boundary::BoundaryKind;
using hundun::boundary::BoundaryRegistry;
using hundun::boundary::FinalFluxDecision;
using hundun::boundary::MassFluxRule;
using hundun::boundary::PressureRule;
using hundun::boundary::TransportRule;
using hundun::boundary::VelocityRule;
using hundun::config::BoundaryType;
using hundun::config::DensityModel;
using hundun::config::FlowBoundaryConfig;
using hundun::config::FlowCaseConfig;
using hundun::config::FlowScalarConfig;
using hundun::config::InletScalarValue;
using hundun::config::InletThermalAuthority;
using hundun::config::MeshMapping;
using hundun::config::PatchName;
using hundun::config::SimulationType;
using hundun::config::TimeMode;
using hundun::finite_volume::PressureConstraintMode;
using hundun::mesh::EntityOwnership;
using hundun::mesh::GlobalFaceId;
using hundun::mesh::MeshTopology;
using hundun::runtime::AccessMode;
using hundun::runtime::ActorId;
using hundun::runtime::DecompositionOptions;
using hundun::runtime::Error;
using hundun::runtime::FieldAccessPlan;
using hundun::runtime::FieldDescriptor;
using hundun::runtime::FieldId;
using hundun::runtime::FieldLayoutSet;
using hundun::runtime::FieldRegistry;
using hundun::runtime::FieldStorage;
using hundun::runtime::FunctionSpace;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::OutputPolicy;
using hundun::runtime::PhaseId;
using hundun::runtime::Real3;
using hundun::runtime::RestartPolicy;
using hundun::runtime::ScalarType;
using hundun::runtime::StructuredDecomposition;
namespace allocation_probe = hundun::test::allocation_probe;

constexpr std::array<PatchName, 6> kPatchNames{
    PatchName::x_min, PatchName::x_max, PatchName::y_min,
    PatchName::y_max, PatchName::z_min, PatchName::z_max};
constexpr std::array<const char *, 6> kPatchText{"x_min", "x_max", "y_min",
                                                 "y_max", "z_min", "z_max"};
constexpr PhaseId kPhase = 51U;
constexpr ActorId kActor = 73U;

FlowBoundaryConfig plain_boundary(std::size_t stable_id, BoundaryType type) {
  FlowBoundaryConfig boundary{};
  boundary.patch = kPatchNames.at(stable_id);
  boundary.type = type;
  return boundary;
}

FlowCaseConfig base_config(int ranks, DensityModel density_model,
                           bool open_domain, std::size_t scalar_count = 0U) {
  FlowCaseConfig config{};
  config.schema_version = 2;
  config.case_name = "task15_boundary";
  config.simulation_type = SimulationType::variable_density_flow;
  config.density_model = density_model;
  config.resources.expected_ranks = ranks;
  config.resources.process_grid = Int3{1, ranks, 1};
  config.mesh.cells = Int3{6, std::max(4, ranks * 2), 3};
  config.mesh.origin_m = Real3{0.0, 0.0, 0.0};
  config.mesh.length_m = Real3{1.0, 1.0, 1.0};
  config.mesh.mapping = MeshMapping::uniform_box;
  config.time.mode = TimeMode::fixed;
  config.time.steps = 2;
  config.time.initial_dt_s = 0.01;
  config.time.min_dt_s = 0.005;
  config.time.max_dt_s = 0.01;
  config.time.cfl_target = 0.5;
  config.time.diffusion_number_target = 0.25;
  config.time.growth_factor = 1.25;
  config.time.retry_factor = 0.5;
  config.time.max_retries = 8;
  config.physics.rho_ref_kg_per_m3 = 1.25;
  config.physics.dynamic_viscosity_pa_s = 1.0e-3;
  config.physics.inlet_consistency_rtol = 1.0e-12;

  const std::array<const char *, 3> names{"alpha", "beta", "zeta"};
  HUNDUN_CHECK(scalar_count <= names.size());
  for (std::size_t scalar = 0; scalar < scalar_count; ++scalar) {
    config.scalars.push_back(FlowScalarConfig{names[scalar], 0.0});
  }

  config.boundaries = {
      plain_boundary(0U, open_domain ? BoundaryType::velocity_inlet
                                     : BoundaryType::no_slip_wall),
      plain_boundary(1U, open_domain ? BoundaryType::pressure_outlet
                                     : BoundaryType::no_slip_wall),
      plain_boundary(2U, BoundaryType::symmetry),
      plain_boundary(3U, BoundaryType::symmetry),
      plain_boundary(4U, BoundaryType::periodic),
      plain_boundary(5U, BoundaryType::periodic)};

  if (open_domain) {
    auto &inlet = config.boundaries[0];
    inlet.velocity_m_per_s = Real3{2.0, -0.25, 0.5};
    inlet.thermal_authority = InletThermalAuthority::enthalpy;
    inlet.enthalpy_J_per_kg = 12.5;
    std::vector<InletScalarValue> values;
    for (std::size_t scalar = scalar_count; scalar > 0U; --scalar) {
      values.push_back(InletScalarValue{names[scalar - 1U], 0.125 * scalar});
    }
    inlet.scalar_values = std::move(values);
    if (density_model == DensityModel::material) {
      inlet.density_kg_per_m3 = 2.5;
    }
    config.boundaries[1].pressure_perturbation_pa = 17.0;
  }

  if (density_model == DensityModel::ideal_gas) {
    config.physics.cp_J_per_kg_K = 1005.0;
    config.physics.gas_constant_J_per_kg_K = 287.0;
    config.physics.thermodynamic_pressure_pa = 101325.0;
    if (open_domain) {
      auto &inlet = config.boundaries[0];
      inlet.thermal_authority = InletThermalAuthority::temperature;
      inlet.temperature_K = 300.0;
      inlet.enthalpy_J_per_kg = 301500.0;
      inlet.density_kg_per_m3 = 101325.0 / (287.0 * 300.0);
    }
  }

  config.restart.read = false;
  config.restart.write_directory = "checkpoints";
  config.restart.write_interval = 1;
  config.diagnostics.directory = "diagnostics";
  config.diagnostics.write_interval = 1;
  config.diagnostics.write_mesh = false;
  config.performance.enabled = false;
  config.performance.directory = "performance";
  config.performance.warmup_steps = 0;
  config.performance.measured_steps = 1;
  config.performance.repetitions = 1;
  return config;
}

std::array<bool, 3> periodic_axes() { return {false, false, true}; }

template <class Function> void expect_error(Function &&function) {
  bool threw = false;
  try {
    std::forward<Function>(function)();
  } catch (const Error &) {
    threw = true;
  }
  HUNDUN_CHECK(threw);
}

template <class Function>
void expect_collective_error(const MpiContext &mpi, Function &&function) {
  bool threw = false;
  try {
    std::forward<Function>(function)();
  } catch (const Error &) {
    threw = true;
  }
  int local = threw ? 1 : 0;
  int every = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local, &every, 1, MPI_INT, MPI_MIN, mpi.comm()) ==
               MPI_SUCCESS);
  HUNDUN_CHECK(every == 1);
}

template <class Function>
void expect_collective_error_message(const MpiContext &mpi, Function &&function,
                                     std::string_view expected_message) {
  bool matched = false;
  try {
    std::forward<Function>(function)();
  } catch (const Error &error) {
    matched = error.what() == expected_message;
  }
  int local = matched ? 1 : 0;
  int every = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local, &every, 1, MPI_INT, MPI_MIN, mpi.comm()) ==
               MPI_SUCCESS);
  HUNDUN_CHECK(every == 1);
}

template <class Function>
void expect_synchronous_mpi_operation_error(const MpiContext &mpi,
                                            Function &&function) {
  bool typed_error = false;
  try {
    std::forward<Function>(function)();
  } catch (const hundun::runtime::detail::MpiOperationError &) {
    typed_error = true;
  }
  int local = typed_error ? 1 : 0;
  int every = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local, &every, 1, MPI_INT, MPI_MIN, mpi.comm()) ==
               MPI_SUCCESS);
  HUNDUN_CHECK(every == 1);
}

std::uint64_t bits(double value) {
  std::uint64_t result{};
  static_assert(sizeof(result) == sizeof(value));
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

double range_safe_ratio_oracle(double numerator, double denominator_a,
                               double denominator_b) {
  return static_cast<double>(
      static_cast<long double>(numerator) /
      (static_cast<long double>(denominator_a) * denominator_b));
}

void check_scalar_copy(hundun::boundary::ScalarBoundaryValues values,
                       double interior) {
  HUNDUN_CHECK(bits(values.face) == bits(interior));
  HUNDUN_CHECK(bits(values.exterior) == bits(interior));
}

void check_evaluation_allocation_free(const BoundaryRegistry &registry,
                                      std::uint32_t patch_id,
                                      std::size_t scalar_count) {
  const Real3 interior_velocity{1.25, -2.5, 3.75};
  const Real3 area{2.0, -1.0, 3.0};
  allocation_probe::AllocationAttemptGuard guard;
  const auto velocity =
      registry.evaluate_velocity(patch_id, interior_velocity, area);
  const auto pressure = registry.evaluate_pressure(patch_id, 13.0);
  const auto density = registry.evaluate_density(patch_id, 1.75);
  const auto enthalpy = registry.evaluate_enthalpy(patch_id, 42.0);
  double scalar_sum = 0.0;
  for (std::size_t scalar = 0; scalar < scalar_count; ++scalar) {
    const auto value = registry.evaluate_scalar(
        patch_id, scalar, 0.25 + static_cast<double>(scalar));
    scalar_sum += value.face + value.exterior;
  }
  HUNDUN_CHECK(std::isfinite(velocity.face.x + velocity.exterior.y +
                             pressure.face + density.exterior + enthalpy.face +
                             scalar_sum));
  HUNDUN_CHECK(guard.attempts() == 0U);
}

FieldDescriptor face_descriptor(std::uint32_t components) {
  return FieldDescriptor{"final_face_mass_flux",
                         "kg/s",
                         "task15_test",
                         FunctionSpace::face_value,
                         ScalarType::float64,
                         components,
                         0,
                         true,
                         RestartPolicy::transient,
                         OutputPolicy::never};
}

void test_preflight_message_conversion_allocation_free() {
  const Error project_error(std::string(1024U, 'p'));
  const std::runtime_error standard_error(std::string(1024U, 's'));
  static_assert(noexcept(hundun::boundary::detail::fixed_preflight_message(
      std::declval<const Error &>())));
  static_assert(noexcept(hundun::boundary::detail::fixed_preflight_message(
      std::declval<const std::exception &>())));
  static_assert(
      noexcept(hundun::boundary::detail::fixed_unknown_preflight_message()));

  std::string_view project_message;
  std::string_view standard_message;
  std::string_view unknown_message;
  allocation_probe::AllocationAttemptGuard guard;
  project_message =
      hundun::boundary::detail::fixed_preflight_message(project_error);
  standard_message =
      hundun::boundary::detail::fixed_preflight_message(standard_error);
  unknown_message = hundun::boundary::detail::fixed_unknown_preflight_message();
  HUNDUN_CHECK(guard.attempts() == 0U);
  HUNDUN_CHECK(project_message ==
               "final pressure-outlet preflight rejected local data");
  HUNDUN_CHECK(standard_message ==
               "final pressure-outlet preflight failed locally");
  HUNDUN_CHECK(unknown_message ==
               "final pressure-outlet preflight failed unexpectedly");
}

void check_rules(const BoundaryRegistry &registry, std::uint32_t id,
                 BoundaryKind kind, VelocityRule velocity,
                 PressureRule pressure, TransportRule density,
                 TransportRule transport, MassFluxRule mass_flux) {
  const auto &descriptor = registry.patch(id);
  HUNDUN_CHECK(descriptor.stable_id() == id);
  HUNDUN_CHECK(descriptor.name() == kPatchText.at(id));
  HUNDUN_CHECK(descriptor.kind() == kind);
  HUNDUN_CHECK(descriptor.velocity_rule() == velocity);
  HUNDUN_CHECK(descriptor.pressure_rule() == pressure);
  HUNDUN_CHECK(descriptor.density_rule() == density);
  HUNDUN_CHECK(descriptor.enthalpy_rule() == transport);
  HUNDUN_CHECK(descriptor.scalar_rule() == transport);
  HUNDUN_CHECK(descriptor.mass_flux_rule() == mass_flux);
}

void test_closed_registry_and_zero_touch(const MpiContext &mpi, int size) {
  FlowCaseConfig config = base_config(size, DensityModel::constant, false, 2U);
  std::rotate(config.boundaries.begin(), config.boundaries.begin() + 2,
              config.boundaries.end());
  auto decomposition = StructuredDecomposition::create(
      mpi, config.mesh.cells, periodic_axes(),
      DecompositionOptions{config.resources.process_grid});
  MeshTopology topology(decomposition);
  BoundaryRegistry registry = BoundaryRegistry::create(config, topology);

  HUNDUN_CHECK(registry.scalar_count() == 2U);
  HUNDUN_CHECK(!registry.open_domain());
  HUNDUN_CHECK(!registry.velocity_inlet_patch_id().has_value());
  HUNDUN_CHECK(!registry.pressure_outlet_patch_id().has_value());
  HUNDUN_CHECK(registry.scalar_name(0U) == "alpha");
  HUNDUN_CHECK(registry.scalar_name(1U) == "beta");
  expect_error([&] { static_cast<void>(registry.scalar_name(2U)); });

  check_rules(registry, 0U, BoundaryKind::no_slip_wall,
              VelocityRule::prescribed_zero, PressureRule::zero_normal_gradient,
              TransportRule::copy_interior,
              TransportRule::zero_normal_diffusive_flux,
              MassFluxRule::identically_zero);
  HUNDUN_CHECK(!registry.patch(0U).inlet_state().has_value());
  HUNDUN_CHECK(!registry.patch(0U).pressure_value_pa().has_value());
  check_rules(registry, 2U, BoundaryKind::symmetry,
              VelocityRule::reflect_normal_copy_tangential,
              PressureRule::zero_normal_gradient, TransportRule::copy_interior,
              TransportRule::zero_normal_diffusive_flux,
              MassFluxRule::identically_zero);
  for (std::uint32_t id : {4U, 5U}) {
    check_rules(registry, id, BoundaryKind::periodic,
                VelocityRule::periodic_pair, PressureRule::periodic_pair,
                TransportRule::periodic_pair, TransportRule::periodic_pair,
                MassFluxRule::periodic_pair);
    HUNDUN_CHECK(registry.patch(id).paired_patch_id() ==
                 std::optional<std::uint32_t>{id == 4U ? 5U : 4U});
  }

  const auto poisson =
      hundun::finite_volume::make_poisson_boundary_spec(registry);
  HUNDUN_CHECK(poisson.mode == PressureConstraintMode::constant_nullspace);
  HUNDUN_CHECK(!poisson.pressure_reference_patch_id.has_value());

  const Real3 wall_interior{3.0, -4.0, 5.0};
  const auto wall =
      registry.evaluate_velocity(0U, wall_interior, Real3{-2.0, 1.0, 0.5});
  HUNDUN_CHECK(bits(wall.face.x) == bits(0.0));
  HUNDUN_CHECK(bits(wall.face.y) == bits(0.0));
  HUNDUN_CHECK(bits(wall.face.z) == bits(0.0));
  HUNDUN_CHECK_NEAR(wall.exterior.x, -wall_interior.x, 0.0);
  HUNDUN_CHECK_NEAR(wall.exterior.y, -wall_interior.y, 0.0);
  HUNDUN_CHECK_NEAR(wall.exterior.z, -wall_interior.z, 0.0);
  for (std::uint32_t patch_id : {0U, 2U}) {
    check_scalar_copy(registry.evaluate_pressure(patch_id, -17.0), -17.0);
    check_scalar_copy(registry.evaluate_density(patch_id, 1.75), 1.75);
    check_scalar_copy(registry.evaluate_enthalpy(patch_id, 91.0), 91.0);
    check_scalar_copy(registry.evaluate_scalar(patch_id, 0U, -0.25), -0.25);
    check_scalar_copy(registry.evaluate_scalar(patch_id, 1U, 0.75), 0.75);
  }
  check_evaluation_allocation_free(registry, 0U, registry.scalar_count());
  check_evaluation_allocation_free(registry, 2U, registry.scalar_count());

  FieldRegistry fields;
  const FieldId flux = fields.declare_field(face_descriptor(1U));
  fields.freeze();
  FieldAccessPlan access(fields);
  access.declare_access(kPhase, kActor, flux, AccessMode::read_write);
  access.freeze();
  FieldStorage storage(
      fields, FieldLayoutSet{Int3{1, 1, 1}, topology.local_face_count()});
  auto writer =
      storage.acquire_face_write<double>(access, kPhase, kActor, flux);
  writer(0U, 0) = -9.0;
  const FieldStorage &const_storage = storage;
  auto stale =
      const_storage.acquire_face_read<double>(access, kPhase, kActor, flux);
  storage.begin_rebuild();
  const auto counters_before = mpi.fp64_reduction_counters();
  const auto result =
      registry.assess_final_pressure_outlet_flux(topology, mpi, stale, 4U, 0.5);
  const auto counters_after = mpi.fp64_reduction_counters();
  HUNDUN_CHECK(result.decision == FinalFluxDecision::admissible);
  HUNDUN_CHECK(!result.evidence.has_value());
  HUNDUN_CHECK(counters_after.collective_calls ==
               counters_before.collective_calls);
}

void test_materialization_and_evaluation(const MpiContext &mpi, int size) {
  FlowCaseConfig ideal = base_config(size, DensityModel::ideal_gas, true, 2U);
  ideal.physics.cp_J_per_kg_K = 1.0;
  ideal.physics.gas_constant_J_per_kg_K =
      std::numeric_limits<double>::max() / 2.0;
  ideal.physics.thermodynamic_pressure_pa = std::numeric_limits<double>::max();
  auto &inlet_config = ideal.boundaries[0];
  inlet_config.temperature_K = 4.0;
  inlet_config.enthalpy_J_per_kg =
      std::nextafter(4.0, std::numeric_limits<double>::infinity());
  inlet_config.density_kg_per_m3 =
      std::nextafter(0.5, std::numeric_limits<double>::infinity());

  auto decomposition = StructuredDecomposition::create(
      mpi, ideal.mesh.cells, periodic_axes(),
      DecompositionOptions{ideal.resources.process_grid});
  MeshTopology topology(decomposition);
  BoundaryRegistry registry = BoundaryRegistry::create(ideal, topology);

  HUNDUN_CHECK(registry.open_domain());
  HUNDUN_CHECK(registry.velocity_inlet_patch_id() ==
               std::optional<std::uint32_t>{0U});
  HUNDUN_CHECK(registry.pressure_outlet_patch_id() ==
               std::optional<std::uint32_t>{1U});
  HUNDUN_CHECK(registry.patch(1U).pressure_value_pa() ==
               std::optional<double>{17.0});
  HUNDUN_CHECK(registry.scalar_count() == 2U);
  HUNDUN_CHECK(registry.scalar_name(0U) == "alpha");
  HUNDUN_CHECK(registry.scalar_name(1U) == "beta");

  check_rules(registry, 0U, BoundaryKind::velocity_inlet,
              VelocityRule::prescribed_inlet,
              PressureRule::zero_normal_gradient,
              TransportRule::prescribed_value, TransportRule::prescribed_value,
              MassFluxRule::prescribed_inlet_state);
  check_rules(registry, 1U, BoundaryKind::pressure_outlet,
              VelocityRule::pure_outflow, PressureRule::prescribed_value,
              TransportRule::pure_outflow, TransportRule::pure_outflow,
              MassFluxRule::outflow_only);

  const auto &inlet = registry.patch(0U).inlet_state();
  HUNDUN_CHECK(inlet.has_value());
  HUNDUN_CHECK(bits(inlet->temperature_K.value()) == bits(4.0));
  HUNDUN_CHECK(bits(inlet->enthalpy_J_per_kg) == bits(4.0));
  HUNDUN_CHECK(bits(inlet->density_kg_per_m3) == bits(0.5));
  HUNDUN_CHECK(inlet->scalar_values.size() == 2U);
  HUNDUN_CHECK_NEAR(inlet->scalar_values[0], 0.125, 0.0);
  HUNDUN_CHECK_NEAR(inlet->scalar_values[1], 0.25, 0.0);

  const Real3 interior{4.0, -2.0, 1.5};
  const Real3 area{2.0, -1.0, 3.0};
  const auto symmetry = registry.evaluate_velocity(2U, interior, area);
  const long double norm = std::sqrt(static_cast<long double>(area.x) * area.x +
                                     static_cast<long double>(area.y) * area.y +
                                     static_cast<long double>(area.z) * area.z);
  const std::array<long double, 3> normal{
      static_cast<long double>(area.x) / norm,
      static_cast<long double>(area.y) / norm,
      static_cast<long double>(area.z) / norm};
  const long double projection =
      static_cast<long double>(interior.x) * normal[0] +
      static_cast<long double>(interior.y) * normal[1] +
      static_cast<long double>(interior.z) * normal[2];
  const std::array<double, 3> expected_exterior{
      static_cast<double>(interior.x - 2.0L * projection * normal[0]),
      static_cast<double>(interior.y - 2.0L * projection * normal[1]),
      static_cast<double>(interior.z - 2.0L * projection * normal[2])};
  HUNDUN_CHECK_NEAR(symmetry.exterior.x, expected_exterior[0], 1.0e-14);
  HUNDUN_CHECK_NEAR(symmetry.exterior.y, expected_exterior[1], 1.0e-14);
  HUNDUN_CHECK_NEAR(symmetry.exterior.z, expected_exterior[2], 1.0e-14);
  const double face_normal = symmetry.face.x * static_cast<double>(normal[0]) +
                             symmetry.face.y * static_cast<double>(normal[1]) +
                             symmetry.face.z * static_cast<double>(normal[2]);
  HUNDUN_CHECK_NEAR(face_normal, 0.0, 2.0e-15);

  const auto wall = registry.evaluate_velocity(0U + 0U, interior, area);
  HUNDUN_CHECK_NEAR(wall.face.x, inlet->velocity_m_per_s.x, 0.0);
  HUNDUN_CHECK_NEAR(wall.exterior.x,
                    2.0 * inlet->velocity_m_per_s.x - interior.x, 0.0);
  const auto outlet_velocity = registry.evaluate_velocity(1U, interior, area);
  HUNDUN_CHECK_NEAR(outlet_velocity.face.x, interior.x, 0.0);
  HUNDUN_CHECK_NEAR(outlet_velocity.exterior.y, interior.y, 0.0);

  const auto inlet_pressure = registry.evaluate_pressure(0U, 33.0);
  HUNDUN_CHECK_NEAR(inlet_pressure.face, 33.0, 0.0);
  const auto outlet_pressure = registry.evaluate_pressure(1U, 33.0);
  HUNDUN_CHECK_NEAR(outlet_pressure.face, 17.0, 0.0);
  HUNDUN_CHECK_NEAR(outlet_pressure.exterior, 1.0, 0.0);
  const auto inlet_density = registry.evaluate_density(0U, 9.0);
  HUNDUN_CHECK_NEAR(inlet_density.face, 0.5, 0.0);
  HUNDUN_CHECK_NEAR(inlet_density.exterior, -8.0, 0.0);
  const auto inlet_h = registry.evaluate_enthalpy(0U, 2.0);
  HUNDUN_CHECK_NEAR(inlet_h.face, 4.0, 0.0);
  HUNDUN_CHECK_NEAR(inlet_h.exterior, 6.0, 0.0);
  HUNDUN_CHECK_NEAR(registry.evaluate_scalar(0U, 0U, 0.75).face, 0.125, 0.0);
  HUNDUN_CHECK_NEAR(registry.evaluate_scalar(0U, 1U, 0.75).face, 0.25, 0.0);
  HUNDUN_CHECK_NEAR(registry.evaluate_density(1U, 9.0).face, 9.0, 0.0);
  HUNDUN_CHECK_NEAR(registry.evaluate_enthalpy(1U, 2.0).exterior, 2.0, 0.0);
  HUNDUN_CHECK_NEAR(registry.evaluate_scalar(1U, 0U, 0.75).exterior, 0.75, 0.0);

  check_evaluation_allocation_free(registry, 0U, registry.scalar_count());
  check_evaluation_allocation_free(registry, 1U, registry.scalar_count());

  for (std::uint32_t periodic : {4U, 5U}) {
    expect_error([&] {
      static_cast<void>(registry.evaluate_velocity(periodic, interior, area));
    });
    expect_error(
        [&] { static_cast<void>(registry.evaluate_pressure(periodic, 1.0)); });
    expect_error(
        [&] { static_cast<void>(registry.evaluate_density(periodic, 1.0)); });
    expect_error(
        [&] { static_cast<void>(registry.evaluate_enthalpy(periodic, 1.0)); });
    expect_error([&] {
      static_cast<void>(registry.evaluate_scalar(periodic, 0U, 1.0));
    });
  }
  expect_error([&] {
    static_cast<void>(registry.evaluate_velocity(2U, interior, Real3{}));
  });
  expect_error([&] {
    static_cast<void>(registry.evaluate_velocity(
        2U, Real3{std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0}, area));
  });
  expect_error([&] {
    static_cast<void>(registry.evaluate_pressure(
        1U, std::numeric_limits<double>::infinity()));
  });
  expect_error(
      [&] { static_cast<void>(registry.evaluate_scalar(0U, 2U, 0.0)); });

  const auto poisson =
      hundun::finite_volume::make_poisson_boundary_spec(registry);
  HUNDUN_CHECK(poisson.mode ==
               PressureConstraintMode::pressure_reference_patch);
  HUNDUN_CHECK(poisson.pressure_reference_patch_id ==
               std::optional<std::uint32_t>{1U});

  FlowCaseConfig constant = base_config(size, DensityModel::constant, true, 0U);
  BoundaryRegistry constant_registry =
      BoundaryRegistry::create(constant, topology);
  const auto &constant_inlet = constant_registry.patch(0U).inlet_state();
  HUNDUN_CHECK(constant_inlet.has_value());
  HUNDUN_CHECK_NEAR(constant_inlet->density_kg_per_m3,
                    constant.physics.rho_ref_kg_per_m3, 0.0);
  HUNDUN_CHECK(!constant_inlet->temperature_K.has_value());
  HUNDUN_CHECK(constant_inlet->scalar_values.empty());

  FlowCaseConfig material = base_config(size, DensityModel::material, true, 1U);
  BoundaryRegistry material_registry =
      BoundaryRegistry::create(material, topology);
  HUNDUN_CHECK_NEAR(
      material_registry.patch(0U).inlet_state()->density_kg_per_m3, 2.5, 0.0);

  FlowCaseConfig enthalpy_authority = ideal;
  enthalpy_authority.physics.cp_J_per_kg_K = 4.0;
  enthalpy_authority.physics.gas_constant_J_per_kg_K =
      std::numeric_limits<double>::max() / 4.0;
  enthalpy_authority.physics.thermodynamic_pressure_pa =
      std::numeric_limits<double>::max();
  auto &enthalpy_inlet = enthalpy_authority.boundaries[0];
  enthalpy_inlet.thermal_authority = InletThermalAuthority::enthalpy;
  enthalpy_inlet.enthalpy_J_per_kg = 20.0;
  enthalpy_inlet.temperature_K =
      std::nextafter(5.0, std::numeric_limits<double>::infinity());
  const double enthalpy_density = range_safe_ratio_oracle(
      *enthalpy_authority.physics.thermodynamic_pressure_pa,
      *enthalpy_authority.physics.gas_constant_J_per_kg_K, 5.0);
  enthalpy_inlet.density_kg_per_m3 =
      std::nextafter(enthalpy_density, std::numeric_limits<double>::infinity());
  BoundaryRegistry enthalpy_registry =
      BoundaryRegistry::create(enthalpy_authority, topology);
  const auto &derived_enthalpy_state =
      enthalpy_registry.patch(0U).inlet_state().value();
  HUNDUN_CHECK(bits(derived_enthalpy_state.enthalpy_J_per_kg) == bits(20.0));
  HUNDUN_CHECK(bits(derived_enthalpy_state.temperature_K.value()) == bits(5.0));
  HUNDUN_CHECK(bits(derived_enthalpy_state.density_kg_per_m3) ==
               bits(enthalpy_density));

  auto temperature_bad_h = ideal;
  temperature_bad_h.boundaries[0].enthalpy_J_per_kg = 4.25;
  expect_error([&] {
    static_cast<void>(BoundaryRegistry::create(temperature_bad_h, topology));
  });
  auto temperature_bad_density = ideal;
  temperature_bad_density.boundaries[0].density_kg_per_m3 = 0.75;
  expect_error([&] {
    static_cast<void>(
        BoundaryRegistry::create(temperature_bad_density, topology));
  });
  auto enthalpy_bad_temperature = enthalpy_authority;
  enthalpy_bad_temperature.boundaries[0].temperature_K = 5.25;
  expect_error([&] {
    static_cast<void>(
        BoundaryRegistry::create(enthalpy_bad_temperature, topology));
  });
  auto enthalpy_bad_density = enthalpy_authority;
  enthalpy_bad_density.boundaries[0].density_kg_per_m3 =
      enthalpy_density * 1.25;
  expect_error([&] {
    static_cast<void>(BoundaryRegistry::create(enthalpy_bad_density, topology));
  });

  FlowCaseConfig huge = ideal;
  huge.boundaries[0].velocity_m_per_s =
      Real3{std::numeric_limits<double>::max(), 0.0, 0.0};
  BoundaryRegistry huge_registry = BoundaryRegistry::create(huge, topology);
  expect_error([&] {
    static_cast<void>(huge_registry.evaluate_velocity(
        0U, Real3{-std::numeric_limits<double>::max(), 0.0, 0.0}, area));
  });
}

void test_registry_rejections(const MpiContext &mpi, int size) {
  FlowCaseConfig valid = base_config(size, DensityModel::material, true, 2U);
  auto decomposition = StructuredDecomposition::create(
      mpi, valid.mesh.cells, periodic_axes(),
      DecompositionOptions{valid.resources.process_grid});
  MeshTopology topology(decomposition);

  auto duplicate = valid;
  duplicate.boundaries[1].patch = PatchName::x_min;
  expect_error([&] {
    static_cast<void>(BoundaryRegistry::create(duplicate, topology));
  });
  auto invalid_patch = valid;
  invalid_patch.boundaries[0].patch = static_cast<PatchName>(77);
  expect_error([&] {
    static_cast<void>(BoundaryRegistry::create(invalid_patch, topology));
  });
  auto invalid_kind = valid;
  invalid_kind.boundaries[3].type = static_cast<BoundaryType>(77);
  expect_error([&] {
    static_cast<void>(BoundaryRegistry::create(invalid_kind, topology));
  });
  auto periodic_mismatch = valid;
  periodic_mismatch.boundaries[5].type = BoundaryType::symmetry;
  expect_error([&] {
    static_cast<void>(BoundaryRegistry::create(periodic_mismatch, topology));
  });
  auto incomplete_inlet = valid;
  incomplete_inlet.boundaries[0].velocity_m_per_s.reset();
  expect_error([&] {
    static_cast<void>(BoundaryRegistry::create(incomplete_inlet, topology));
  });
  auto wrong_scalars = valid;
  wrong_scalars.boundaries[0].scalar_values->pop_back();
  expect_error([&] {
    static_cast<void>(BoundaryRegistry::create(wrong_scalars, topology));
  });
  auto invalid_open_pair = valid;
  invalid_open_pair.boundaries[1] =
      plain_boundary(1U, BoundaryType::no_slip_wall);
  expect_error([&] {
    static_cast<void>(BoundaryRegistry::create(invalid_open_pair, topology));
  });
  auto invalid_pressure = valid;
  invalid_pressure.boundaries[1].pressure_perturbation_pa =
      std::numeric_limits<double>::quiet_NaN();
  expect_error([&] {
    static_cast<void>(BoundaryRegistry::create(invalid_pressure, topology));
  });
  auto invalid_density = valid;
  invalid_density.boundaries[0].density_kg_per_m3 = 0.0;
  expect_error([&] {
    static_cast<void>(BoundaryRegistry::create(invalid_density, topology));
  });
  auto unsorted_scalars = valid;
  std::swap(unsorted_scalars.scalars[0], unsorted_scalars.scalars[1]);
  expect_error([&] {
    static_cast<void>(BoundaryRegistry::create(unsorted_scalars, topology));
  });

  FlowCaseConfig periodic_reference = valid;
  periodic_reference.boundaries[0] = plain_boundary(0U, BoundaryType::periodic);
  periodic_reference.boundaries[1] = plain_boundary(1U, BoundaryType::periodic);
  periodic_reference.boundaries[2] =
      plain_boundary(2U, BoundaryType::velocity_inlet);
  periodic_reference.boundaries[2].velocity_m_per_s = Real3{1.0, 0.0, 0.0};
  periodic_reference.boundaries[2].thermal_authority =
      InletThermalAuthority::enthalpy;
  periodic_reference.boundaries[2].enthalpy_J_per_kg = 1.0;
  periodic_reference.boundaries[2].density_kg_per_m3 = 1.0;
  periodic_reference.boundaries[2].scalar_values =
      valid.boundaries[0].scalar_values;
  periodic_reference.boundaries[3] =
      plain_boundary(3U, BoundaryType::pressure_outlet);
  periodic_reference.boundaries[3].pressure_perturbation_pa = 0.0;
  expect_error([&] {
    static_cast<void>(BoundaryRegistry::create(periodic_reference, topology));
  });
}

void check_counter_delta(hundun::runtime::Fp64ReductionCounters before,
                         hundun::runtime::Fp64ReductionCounters after,
                         std::uint64_t calls, std::uint64_t scalars,
                         std::uint64_t bytes_count) {
  HUNDUN_CHECK(after.collective_calls - before.collective_calls == calls);
  HUNDUN_CHECK(after.reduced_scalars - before.reduced_scalars == scalars);
  HUNDUN_CHECK(after.logical_payload_bytes - before.logical_payload_bytes ==
               bytes_count);
}

void test_final_flux_assessment(const MpiContext &mpi, int rank, int size) {
  FlowCaseConfig config = base_config(size, DensityModel::constant, true);
  auto decomposition = StructuredDecomposition::create(
      mpi, config.mesh.cells, periodic_axes(),
      DecompositionOptions{config.resources.process_grid});
  MeshTopology topology(decomposition);
  BoundaryRegistry registry = BoundaryRegistry::create(config, topology);

  FieldRegistry fields;
  const FieldId flux = fields.declare_field(face_descriptor(1U));
  fields.freeze();
  FieldAccessPlan access(fields);
  access.declare_access(kPhase, kActor, flux, AccessMode::read_write);
  access.freeze();
  FieldStorage storage(
      fields, FieldLayoutSet{Int3{1, 1, 1}, topology.local_face_count()});
  auto writer =
      storage.acquire_face_write<double>(access, kPhase, kActor, flux);
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    writer(face, 0) = 3.0 + static_cast<double>(face % 7U);
  }
  const auto &outlet_faces = topology.patch(1U).local_faces();
  HUNDUN_CHECK(!outlet_faces.empty());
  std::vector<hundun::mesh::LocalFaceId> owned_outlets;
  for (const auto face : outlet_faces) {
    if (topology.face_ownership(face) == EntityOwnership::owned) {
      owned_outlets.push_back(face);
    }
  }
  HUNDUN_CHECK(owned_outlets.size() >= 2U);
  writer(owned_outlets.front(), 0) = rank % 2 == 0 ? -0.0 : 0.0;

  const FieldStorage &const_storage = storage;
  auto reader =
      const_storage.acquire_face_read<double>(access, kPhase, kActor, flux);
  std::vector<std::uint64_t> before_bits(topology.local_face_count());
  for (std::size_t face = 0; face < before_bits.size(); ++face) {
    before_bits[face] = bits(reader(face, 0));
  }

  using TopologyObservation =
      hundun::boundary::detail::OutletTopologyObservationForTest;
  constexpr std::array<TopologyObservation, 6> topology_identity_changes{
      TopologyObservation::local_id_out_of_range,
      TopologyObservation::ownership_is_ghost,
      TopologyObservation::patch_id_is_not_outlet,
      TopologyObservation::outlet_patch_does_not_contain,
      TopologyObservation::global_id_differs,
      TopologyObservation::owned_cardinality_differs};
  for (const auto observation : topology_identity_changes) {
    if (rank == 0) {
      hundun::boundary::detail::set_next_outlet_topology_observation_for_test(
          observation);
    }
    const auto identity_counters_before = mpi.fp64_reduction_counters();
    expect_collective_error_message(
        mpi,
        [&] {
          static_cast<void>(registry.assess_final_pressure_outlet_flux(
              topology, mpi, reader, 6U, 1.0));
        },
        "final pressure-outlet preflight rejected local data");
    const auto identity_counters_after = mpi.fp64_reduction_counters();
    check_counter_delta(identity_counters_before, identity_counters_after, 0U,
                        0U, 0U);
    for (std::size_t face = 0; face < before_bits.size(); ++face) {
      HUNDUN_CHECK(bits(reader(face, 0)) == before_bits[face]);
    }
  }

  auto counters_before = mpi.fp64_reduction_counters();
  std::size_t admissible_assessment_allocations = 0U;
  auto result = [&] {
    allocation_probe::AllocationAttemptGuard guard;
    auto assessment = registry.assess_final_pressure_outlet_flux(
        topology, mpi, reader, 7U, 1.25);
    admissible_assessment_allocations = guard.attempts();
    return assessment;
  }();
  auto counters_after = mpi.fp64_reduction_counters();
  HUNDUN_CHECK(admissible_assessment_allocations == 0U);
  HUNDUN_CHECK(result.decision == FinalFluxDecision::admissible);
  HUNDUN_CHECK(!result.evidence.has_value());
  check_counter_delta(counters_before, counters_after, 1U, 1U, 8U);
  for (std::size_t face = 0; face < before_bits.size(); ++face) {
    HUNDUN_CHECK(bits(reader(face, 0)) == before_bits[face]);
  }

  writer(owned_outlets.front(), 0) = -2.0;
  for (std::size_t face = 0; face < before_bits.size(); ++face) {
    before_bits[face] = bits(reader(face, 0));
  }
  counters_before = mpi.fp64_reduction_counters();
  // Every rank arms this pre-call seam, which skips the MPI operation.  This
  // verifies typed operation-error propagation and sequence cutoff, not
  // recovery from a real partial-rank collective failure.
  hundun::runtime::detail::
      inject_synchronous_next_fp64_allreduce_pre_call_error_for_test(
          MPI_ERR_OTHER);
  expect_synchronous_mpi_operation_error(mpi, [&] {
    static_cast<void>(registry.assess_final_pressure_outlet_flux(
        topology, mpi, reader, 7U, 1.25));
  });
  counters_after = mpi.fp64_reduction_counters();
  check_counter_delta(counters_before, counters_after, 1U, 1U, 8U);
  for (std::size_t face = 0; face < before_bits.size(); ++face) {
    HUNDUN_CHECK(bits(reader(face, 0)) == before_bits[face]);
  }

  writer(owned_outlets.front(), 0) = -1.0;
  GlobalFaceId local_tied_id = std::numeric_limits<GlobalFaceId>::max();
  const bool contributes_most_severe =
      size == 1 || rank == size - 1 || (size >= 4 && rank == size - 2);
  if (contributes_most_severe) {
    writer(owned_outlets[0], 0) = -5.0;
    writer(owned_outlets[1], 0) = -5.0;
    local_tied_id = std::min(topology.global_face_id(owned_outlets[0]),
                             topology.global_face_id(owned_outlets[1]));
  }
  GlobalFaceId expected_id = local_tied_id;
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &expected_id, 1, MPI_UINT64_T,
                             MPI_MIN, mpi.comm()) == MPI_SUCCESS);
  before_bits.resize(topology.local_face_count());
  for (std::size_t face = 0; face < before_bits.size(); ++face) {
    before_bits[face] = bits(reader(face, 0));
  }
  counters_before = mpi.fp64_reduction_counters();
  result = registry.assess_final_pressure_outlet_flux(topology, mpi, reader, 8U,
                                                      1.5);
  counters_after = mpi.fp64_reduction_counters();
  HUNDUN_CHECK(result.decision == FinalFluxDecision::outlet_backflow);
  HUNDUN_CHECK(result.evidence.has_value());
  HUNDUN_CHECK(result.evidence->patch_id == 1U);
  HUNDUN_CHECK(result.evidence->step == 8U);
  HUNDUN_CHECK_NEAR(result.evidence->time_s, 1.5, 0.0);
  HUNDUN_CHECK(bits(result.evidence->minimum_outward_mass_flux_kg_per_s) ==
               bits(-5.0));
  HUNDUN_CHECK(result.evidence->global_face_id == expected_id);
  HUNDUN_CHECK(result.evidence->lowest_failing_rank == 0);
  check_counter_delta(counters_before, counters_after, 4U, 4U, 32U);
  for (std::size_t face = 0; face < before_bits.size(); ++face) {
    HUNDUN_CHECK(bits(reader(face, 0)) == before_bits[face]);
  }

  auto stale = reader;
  storage.begin_repartition();
  expect_collective_error(mpi, [&] {
    static_cast<void>(registry.assess_final_pressure_outlet_flux(
        topology, mpi, stale, 9U, 1.75));
  });

  FieldRegistry wrong_components_fields;
  const FieldId wrong_components =
      wrong_components_fields.declare_field(face_descriptor(2U));
  wrong_components_fields.freeze();
  FieldAccessPlan wrong_components_access(wrong_components_fields);
  wrong_components_access.declare_access(kPhase, kActor, wrong_components,
                                         AccessMode::read_write);
  wrong_components_access.freeze();
  FieldStorage wrong_components_storage(
      wrong_components_fields,
      FieldLayoutSet{Int3{1, 1, 1}, topology.local_face_count()});
  const FieldStorage &wrong_components_const = wrong_components_storage;
  auto wrong_components_view = wrong_components_const.acquire_face_read<double>(
      wrong_components_access, kPhase, kActor, wrong_components);
  expect_collective_error(mpi, [&] {
    static_cast<void>(registry.assess_final_pressure_outlet_flux(
        topology, mpi, wrong_components_view, 10U, 2.0));
  });

  FieldRegistry wrong_count_fields;
  const FieldId wrong_count =
      wrong_count_fields.declare_field(face_descriptor(1U));
  wrong_count_fields.freeze();
  FieldAccessPlan wrong_count_access(wrong_count_fields);
  wrong_count_access.declare_access(kPhase, kActor, wrong_count,
                                    AccessMode::read_write);
  wrong_count_access.freeze();
  FieldStorage wrong_count_storage(
      wrong_count_fields,
      FieldLayoutSet{Int3{1, 1, 1}, topology.local_face_count() - 1U});
  const FieldStorage &wrong_count_const = wrong_count_storage;
  auto wrong_count_view = wrong_count_const.acquire_face_read<double>(
      wrong_count_access, kPhase, kActor, wrong_count);
  expect_collective_error(mpi, [&] {
    static_cast<void>(registry.assess_final_pressure_outlet_flux(
        topology, mpi, wrong_count_view, 11U, 2.25));
  });

  FieldRegistry nonfinite_fields;
  const FieldId nonfinite = nonfinite_fields.declare_field(face_descriptor(1U));
  nonfinite_fields.freeze();
  FieldAccessPlan nonfinite_access(nonfinite_fields);
  nonfinite_access.declare_access(kPhase, kActor, nonfinite,
                                  AccessMode::read_write);
  nonfinite_access.freeze();
  FieldStorage nonfinite_storage(
      nonfinite_fields,
      FieldLayoutSet{Int3{1, 1, 1}, topology.local_face_count()});
  auto nonfinite_writer = nonfinite_storage.acquire_face_write<double>(
      nonfinite_access, kPhase, kActor, nonfinite);
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    nonfinite_writer(face, 0) = 1.0;
  }
  if (rank == 0) {
    nonfinite_writer(owned_outlets.front(), 0) =
        std::numeric_limits<double>::quiet_NaN();
  }
  const FieldStorage &nonfinite_const = nonfinite_storage;
  auto nonfinite_reader = nonfinite_const.acquire_face_read<double>(
      nonfinite_access, kPhase, kActor, nonfinite);
  std::vector<std::uint64_t> nonfinite_bits(topology.local_face_count());
  for (std::size_t face = 0; face < nonfinite_bits.size(); ++face) {
    nonfinite_bits[face] = bits(nonfinite_reader(face, 0));
  }
  expect_collective_error(mpi, [&] {
    static_cast<void>(registry.assess_final_pressure_outlet_flux(
        topology, mpi, nonfinite_reader, 12U, 2.5));
  });
  for (std::size_t face = 0; face < nonfinite_bits.size(); ++face) {
    HUNDUN_CHECK(bits(nonfinite_reader(face, 0)) == nonfinite_bits[face]);
  }
  nonfinite_writer(owned_outlets.front(), 0) =
      rank == 0 ? std::numeric_limits<double>::infinity() : 1.0;
  expect_collective_error(mpi, [&] {
    static_cast<void>(registry.assess_final_pressure_outlet_flux(
        topology, mpi, nonfinite_reader, 13U, 2.75));
  });
  expect_collective_error_message(
      mpi,
      [&] {
        static_cast<void>(registry.assess_final_pressure_outlet_flux(
            topology, mpi, nonfinite_reader, 14U, rank == 0 ? -1.0 : 3.0));
      },
      "final pressure-outlet preflight rejected local data");
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    MpiContext mpi = MpiContext::duplicate(MPI_COMM_WORLD);
    const int rank = mpi.rank();
    const int size = mpi.size();
    HUNDUN_CHECK(size == 1 || size == 2 || size == 4);
    test_preflight_message_conversion_allocation_free();
    test_closed_registry_and_zero_touch(mpi, size);
    test_materialization_and_evaluation(mpi, size);
    test_registry_rejections(mpi, size);
    test_final_flux_assessment(mpi, rank, size);
  });
}

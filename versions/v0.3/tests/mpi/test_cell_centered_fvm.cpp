// SPDX-License-Identifier: Apache-2.0

#include "tests/support/fvm_cell_centered_test_seam.hpp"
#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_access_plan.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_field_storage.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/allocation_attempt_guard.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using hundun::boundary::BoundaryRegistry;
using hundun::config::BoundaryType;
using hundun::config::DensityModel;
using hundun::config::FlowBoundaryConfig;
using hundun::config::FlowCaseConfig;
using hundun::config::FlowScalarConfig;
using hundun::config::InletScalarValue;
using hundun::config::InletThermalAuthority;
using hundun::config::PatchName;
using hundun::config::SimulationType;
using hundun::finite_volume::CellCenteredFvmOperators;
using hundun::finite_volume::declare_face_mass_flux;
using hundun::finite_volume::FaceMassFlux;
using hundun::finite_volume::FiniteVolumeQuantity;
using hundun::finite_volume::GradientScheme;
using hundun::finite_volume::monotonized_central;
using hundun::finite_volume::require_face_mass_flux_field;
using hundun::mesh::AnalyticWarpedBoxMapping;
using hundun::mesh::LocalCellId;
using hundun::mesh::MeshGeometry;
using hundun::mesh::MeshTopology;
using hundun::mesh::UniformBoxMapping;
using hundun::runtime::AccessMode;
using hundun::runtime::ActorId;
using hundun::runtime::DecompositionOptions;
using hundun::runtime::ExchangePlan;
using hundun::runtime::FieldAccessPlan;
using hundun::runtime::FieldDescriptor;
using hundun::runtime::FieldId;
using hundun::runtime::FieldLayoutSet;
using hundun::runtime::FieldRegistry;
using hundun::runtime::FieldStorage;
using hundun::runtime::FunctionSpace;
using hundun::runtime::HaloExchange;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::OutputPolicy;
using hundun::runtime::PhaseId;
using hundun::runtime::Real3;
using hundun::runtime::RestartPolicy;
using hundun::runtime::ScalarType;
using hundun::runtime::StructuredDecomposition;
namespace allocation_probe = hundun::test::allocation_probe;

Int3 process_grid_for(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw hundun::runtime::Error("unsupported Task 16 test rank count");
}

FlowCaseConfig periodic_case() {
  FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type = SimulationType::variable_density_flow;
  config.density_model = DensityModel::constant;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  config.scalars.push_back(FlowScalarConfig{"alpha", 0.125});
  constexpr std::array<PatchName, 6> names{PatchName::x_min, PatchName::x_max,
                                           PatchName::y_min, PatchName::y_max,
                                           PatchName::z_min, PatchName::z_max};
  for (std::size_t index = 0; index < names.size(); ++index) {
    FlowBoundaryConfig boundary{};
    boundary.patch = names[index];
    boundary.type = BoundaryType::periodic;
    config.boundaries[index] = std::move(boundary);
  }
  return config;
}

FlowCaseConfig mixed_case() {
  FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type = SimulationType::variable_density_flow;
  config.density_model = DensityModel::constant;
  config.physics.rho_ref_kg_per_m3 = 1.5;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  config.scalars.push_back(FlowScalarConfig{"alpha", 0.125});
  constexpr std::array<PatchName, 6> names{PatchName::x_min, PatchName::x_max,
                                           PatchName::y_min, PatchName::y_max,
                                           PatchName::z_min, PatchName::z_max};
  constexpr std::array<BoundaryType, 6> kinds{
      BoundaryType::velocity_inlet, BoundaryType::pressure_outlet,
      BoundaryType::no_slip_wall,   BoundaryType::symmetry,
      BoundaryType::periodic,       BoundaryType::periodic};
  for (std::size_t index = 0; index < names.size(); ++index) {
    config.boundaries[index].patch = names[index];
    config.boundaries[index].type = kinds[index];
  }
  auto &inlet = config.boundaries[0];
  inlet.velocity_m_per_s = Real3{1.0, 0.2, -0.1};
  inlet.thermal_authority = InletThermalAuthority::enthalpy;
  inlet.enthalpy_J_per_kg = 12.0;
  inlet.scalar_values =
      std::vector<InletScalarValue>{InletScalarValue{"alpha", 0.3}};
  config.boundaries[1].pressure_perturbation_pa = 17.0;
  return config;
}

FlowCaseConfig all_no_slip_case() {
  FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type = SimulationType::variable_density_flow;
  config.density_model = DensityModel::constant;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  constexpr std::array<PatchName, 6> names{
      PatchName::x_min, PatchName::x_max, PatchName::y_min,
      PatchName::y_max, PatchName::z_min, PatchName::z_max};
  for (std::size_t index = 0; index < names.size(); ++index) {
    config.boundaries[index].patch = names[index];
    config.boundaries[index].type = BoundaryType::no_slip_wall;
  }
  return config;
}

template <class Function> void expect_error(Function &&function) {
  bool threw = false;
  try {
    std::forward<Function>(function)();
  } catch (const hundun::runtime::Error &) {
    threw = true;
  }
  HUNDUN_CHECK(threw);
}

std::uint64_t bits(double value) {
  std::uint64_t result = 0U;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

FieldDescriptor cell_descriptor(std::string name, int components,
                                int ghost_width) {
  return FieldDescriptor{std::move(name),
                         "1",
                         "task16",
                         FunctionSpace::cell_average,
                         ScalarType::float64,
                         static_cast<std::uint32_t>(components),
                         ghost_width,
                         false,
                         RestartPolicy::transient,
                         OutputPolicy::never};
}

FieldDescriptor face_descriptor(std::string name, int components) {
  return FieldDescriptor{std::move(name),
                         "1",
                         "task16",
                         FunctionSpace::face_value,
                         ScalarType::float64,
                         static_cast<std::uint32_t>(components),
                         0,
                         false,
                         RestartPolicy::transient,
                         OutputPolicy::never};
}

struct CellYMoments final {
  double volume{};
  double y_integral{};
  double y_squared_integral{};
  double y_cubed_integral{};
};

Real3 subtract(Real3 left, Real3 right) noexcept {
  return {left.x - right.x, left.y - right.y, left.z - right.z};
}

double determinant(Real3 a, Real3 b, Real3 c) noexcept {
  return a.x * (b.y * c.z - b.z * c.y) -
         a.y * (b.x * c.z - b.z * c.x) +
         a.z * (b.x * c.y - b.y * c.x);
}

CellYMoments cell_y_moments(Int3 cell, const MeshGeometry &geometry) {
  constexpr std::array<Int3, 8> offsets{
      Int3{0, 0, 0}, Int3{1, 0, 0}, Int3{1, 1, 0}, Int3{0, 1, 0},
      Int3{0, 0, 1}, Int3{1, 0, 1}, Int3{1, 1, 1}, Int3{0, 1, 1}};
  std::array<Real3, 8> vertices{};
  Real3 reference{};
  for (std::size_t index = 0; index < vertices.size(); ++index) {
    vertices[index] = geometry.vertex_position_m(
        {cell.x + offsets[index].x, cell.y + offsets[index].y,
         cell.z + offsets[index].z});
    reference.x += 0.125 * vertices[index].x;
    reference.y += 0.125 * vertices[index].y;
    reference.z += 0.125 * vertices[index].z;
  }
  constexpr std::array<std::array<std::size_t, 3>, 12> triangles{{
      {{0, 4, 7}}, {{0, 7, 3}}, {{1, 2, 6}}, {{1, 6, 5}},
      {{0, 1, 5}}, {{0, 5, 4}}, {{3, 7, 6}}, {{3, 6, 2}},
      {{0, 3, 2}}, {{0, 2, 1}}, {{4, 5, 6}}, {{4, 6, 7}},
  }};
  CellYMoments result;
  for (const auto triangle : triangles) {
    const auto a = vertices[triangle[0]];
    const auto b = vertices[triangle[1]];
    const auto c = vertices[triangle[2]];
    const double volume =
        determinant(subtract(a, reference), subtract(b, reference),
                    subtract(c, reference)) /
        6.0;
    HUNDUN_CHECK(volume > 0.0 && std::isfinite(volume));
    const std::array<double, 4> y{reference.y, a.y, b.y, c.y};
    double sum = 0.0;
    double square_sum = 0.0;
    double cube_sum = 0.0;
    for (const double value : y) {
      sum += value;
      square_sum += value * value;
      cube_sum += value * value * value;
    }
    result.volume += volume;
    result.y_integral += volume * sum / 4.0;
    result.y_squared_integral +=
        volume * (sum * sum + square_sum) / 20.0;
    result.y_cubed_integral +=
        volume *
        (sum * sum * sum + 3.0 * sum * square_sum + 2.0 * cube_sum) /
        120.0;
  }
  HUNDUN_CHECK(result.volume > 0.0 && std::isfinite(result.volume));
  return result;
}

void check_limiter() {
  HUNDUN_CHECK_NEAR(monotonized_central(1.0, 1.0), 1.0, 0.0);
  HUNDUN_CHECK_NEAR(monotonized_central(-1.0, -1.0), -1.0, 0.0);
  HUNDUN_CHECK_NEAR(monotonized_central(1.0, -1.0), 0.0, 0.0);
  HUNDUN_CHECK_NEAR(monotonized_central(1.0, 3.0), 2.0, 0.0);
  HUNDUN_CHECK_NEAR(monotonized_central(3.0, 1.0), 2.0, 0.0);
  HUNDUN_CHECK(!std::signbit(monotonized_central(-0.0, -0.0)));
  HUNDUN_CHECK(
      monotonized_central(std::numeric_limits<double>::infinity(), 1.0) == 0.0);
  HUNDUN_CHECK(monotonized_central(std::numeric_limits<double>::quiet_NaN(),
                                   1.0) == 0.0);
  const double huge = std::numeric_limits<double>::max();
  HUNDUN_CHECK(std::isfinite(monotonized_central(huge, huge)));
  HUNDUN_CHECK(monotonized_central(huge, huge) == huge);
}

void run_gradient_case(const MpiContext &mpi, int ranks, bool warped) {
  constexpr Int3 extent{8, 8, 8};
  const std::array<bool, 3> periodic{true, true, true};
  auto decomposition = StructuredDecomposition::create(
      mpi, extent, periodic, DecompositionOptions{process_grid_for(ranks)});
  const auto fixture_local = decomposition.local_extent();
  HUNDUN_CHECK(fixture_local.x >= 2);
  HUNDUN_CHECK(fixture_local.y >= 2);
  HUNDUN_CHECK(fixture_local.z >= 2);
  MeshTopology topology(decomposition);
  MeshGeometry geometry =
      warped ? MeshGeometry(topology,
                            AnalyticWarpedBoxMapping(Real3{0.0, 0.0, 0.0},
                                                     Real3{1.0, 1.0, 1.0},
                                                     Real3{0.02, -0.015, 0.01}))
             : MeshGeometry(topology, UniformBoxMapping(Real3{0.0, 0.0, 0.0},
                                                        Real3{1.0, 1.0, 1.0}));
  BoundaryRegistry boundaries =
      BoundaryRegistry::create(periodic_case(), topology);
  CellCenteredFvmOperators operators =
      CellCenteredFvmOperators::create(topology, geometry);

  FieldRegistry registry;
  const FieldId q = registry.declare_field(cell_descriptor("q", 1, 2));
  const FieldId gradient =
      registry.declare_field(cell_descriptor("gradient", 3, 1));
  const FieldId velocity =
      registry.declare_field(cell_descriptor("gradient_velocity", 3, 2));
  const FieldId velocity_gradient = registry.declare_field(
      cell_descriptor("gradient_velocity_gradient", 9, 1));
  const FieldId mass_flux = declare_face_mass_flux(registry);
  const FieldId transported_face =
      registry.declare_field(face_descriptor("free_stream_face", 1));
  const FieldId momentum_face =
      registry.declare_field(face_descriptor("free_stream_velocity_face", 3));
  const FieldId gamma =
      registry.declare_field(face_descriptor("free_stream_gamma", 1));
  const FieldId scalar_residual =
      registry.declare_field(cell_descriptor("free_stream_scalar_r", 1, 0));
  const FieldId momentum_convective_residual = registry.declare_field(
      cell_descriptor("free_stream_momentum_convective_r", 3, 0));
  const FieldId momentum_viscous_residual = registry.declare_field(
      cell_descriptor("free_stream_momentum_viscous_r", 3, 0));
  const FieldId scalar_diffusive_residual = registry.declare_field(
      cell_descriptor("free_stream_scalar_diffusive_r", 1, 0));
  registry.freeze();
  FieldStorage storage(registry, FieldLayoutSet{decomposition.local_extent(),
                                                topology.local_face_count()});
  auto q_write = storage.view<double>(q);
  auto gradient_write = storage.view<double>(gradient);
  auto velocity_write = storage.view<double>(velocity);
  auto velocity_gradient_write = storage.view<double>(velocity_gradient);
  const auto owned = topology.owned_global_box();
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        const LocalCellId local =
            (static_cast<std::size_t>(k) *
                 static_cast<std::size_t>(decomposition.local_extent().y) +
             static_cast<std::size_t>(j)) *
                static_cast<std::size_t>(decomposition.local_extent().x) +
            static_cast<std::size_t>(i);
        const Real3 center = geometry.cell_center_m(local);
        q_write(i, j, k, 0) =
            2.0 + 1.25 * center.x - 0.75 * center.y + 0.5 * center.z;
        velocity_write(i, j, k, 0) =
            2.0 + 1.25 * center.x - 0.75 * center.y + 0.5 * center.z;
        velocity_write(i, j, k, 1) =
            -1.0 + 0.2 * center.x + 0.3 * center.y - 0.4 * center.z;
        velocity_write(i, j, k, 2) =
            0.5 - 0.6 * center.x + 0.7 * center.y + 0.8 * center.z;
        for (int component = 0; component < 3; ++component) {
          gradient_write(i, j, k, component) = -999.0;
        }
        for (int component = 0; component < 9; ++component) {
          velocity_gradient_write(i, j, k, component) = -999.0;
        }
      }
    }
  }
  auto halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 2));
  halo.exchange(storage, q);
  halo.exchange(storage, velocity);
  const FieldStorage &read_storage = storage;
  const auto q_read = read_storage.view<double>(q);
  const auto velocity_read = read_storage.view<double>(velocity);
  const GradientScheme linear_scheme =
      warped ? GradientScheme::weighted_least_squares
             : GradientScheme::green_gauss;
  operators.compute_gradient(linear_scheme, FiniteVolumeQuantity::density(),
                             boundaries, q_read, gradient_write);
  operators.compute_gradient(linear_scheme, FiniteVolumeQuantity::velocity(),
                             boundaries, velocity_read,
                             velocity_gradient_write);

  const double tolerance = 512.0 * std::numeric_limits<double>::epsilon() * 4.0;
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        const int gi = owned.begin.x + i;
        const int gj = owned.begin.y + j;
        const int gk = owned.begin.z + k;
        if (gi > 0 && gi + 1 < extent.x && gj > 0 && gj + 1 < extent.y &&
            gk > 0 && gk + 1 < extent.z) {
          HUNDUN_CHECK_NEAR(gradient_write(i, j, k, 0), 1.25, tolerance);
          HUNDUN_CHECK_NEAR(gradient_write(i, j, k, 1), -0.75, tolerance);
          HUNDUN_CHECK_NEAR(gradient_write(i, j, k, 2), 0.5, tolerance);
          constexpr std::array<double, 9> expected{1.25, -0.75, 0.5, 0.2, 0.3,
                                                   -0.4, -0.6,  0.7, 0.8};
          for (int component = 0; component < 9; ++component) {
            HUNDUN_CHECK_NEAR(velocity_gradient_write(i, j, k, component),
                              expected[static_cast<std::size_t>(component)],
                              tolerance);
          }
        }
      }
    }
  }

  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        q_write(i, j, k, 0) = 7.0;
        velocity_write(i, j, k, 0) = 4.0;
        velocity_write(i, j, k, 1) = -3.0;
        velocity_write(i, j, k, 2) = 2.0;
      }
    }
  }
  halo.exchange(storage, q);
  halo.exchange(storage, velocity);
  for (const GradientScheme scheme :
       {GradientScheme::green_gauss, GradientScheme::weighted_least_squares}) {
    operators.compute_gradient(scheme, FiniteVolumeQuantity::density(),
                               boundaries, q_read, gradient_write);
    operators.compute_gradient(scheme, FiniteVolumeQuantity::velocity(),
                               boundaries, velocity_read,
                               velocity_gradient_write);
    for (int k = 0; k < decomposition.local_extent().z; ++k) {
      for (int j = 0; j < decomposition.local_extent().y; ++j) {
        for (int i = 0; i < decomposition.local_extent().x; ++i) {
          HUNDUN_CHECK_NEAR(gradient_write(i, j, k, 0), 0.0, 0.0);
          HUNDUN_CHECK_NEAR(gradient_write(i, j, k, 1), 0.0, 0.0);
          HUNDUN_CHECK_NEAR(gradient_write(i, j, k, 2), 0.0, 0.0);
          for (int component = 0; component < 9; ++component) {
            HUNDUN_CHECK_NEAR(velocity_gradient_write(i, j, k, component), 0.0,
                              0.0);
          }
        }
      }
    }
  }
  if (!warped)
    return;

  constexpr Real3 free_stream_velocity{0.125, -0.25, 0.375};
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        q_write(i, j, k, 0) = 1.0;
        velocity_write(i, j, k, 0) = free_stream_velocity.x;
        velocity_write(i, j, k, 1) = free_stream_velocity.y;
        velocity_write(i, j, k, 2) = free_stream_velocity.z;
      }
    }
  }
  halo.exchange(storage, q);
  halo.exchange(storage, velocity);
  operators.compute_gradient(GradientScheme::weighted_least_squares,
                             FiniteVolumeQuantity::density(), boundaries,
                             q_read, gradient_write);
  operators.compute_gradient(GradientScheme::weighted_least_squares,
                             FiniteVolumeQuantity::velocity(), boundaries,
                             velocity_read, velocity_gradient_write);

  constexpr PhaseId write_phase = 701U;
  constexpr PhaseId read_phase = 702U;
  constexpr ActorId actor = 703U;
  FieldAccessPlan access(registry);
  access.declare_access(write_phase, actor, mass_flux, AccessMode::write);
  access.declare_access(read_phase, actor, mass_flux, AccessMode::read);
  access.declare_access(read_phase, actor, transported_face,
                        AccessMode::read_write);
  access.declare_access(read_phase, actor, momentum_face,
                        AccessMode::read_write);
  access.declare_access(read_phase, actor, gamma, AccessMode::read_write);
  access.freeze();
  operators.assemble_provisional_mass_flux(
      boundaries, q_read, velocity_read, registry, storage, access,
      write_phase, actor, mass_flux);
  const FieldStorage &free_stream_read_storage = storage;
  auto flux = FaceMassFlux::acquire(
      registry, free_stream_read_storage, access, read_phase, actor, mass_flux,
      topology);
  auto transported_face_write = storage.acquire_face_write<double>(
      access, read_phase, actor, transported_face);
  auto momentum_face_write = storage.acquire_face_write<double>(
      access, read_phase, actor, momentum_face);
  auto gamma_write =
      storage.acquire_face_write<double>(access, read_phase, actor, gamma);
  operators.reconstruct_momentum_faces(boundaries, flux, velocity_read,
                                       momentum_face_write);

  const auto flux_read = free_stream_read_storage.acquire_face_read<double>(
      access, read_phase, actor, mass_flux);
  const auto momentum_face_read =
      free_stream_read_storage.acquire_face_read<double>(
          access, read_phase, actor, momentum_face);
  const double eps = std::numeric_limits<double>::epsilon();
  const auto face_tolerance = [&](double scale) {
    return 512.0 * eps * std::max(1.0, scale);
  };
  for (hundun::mesh::LocalFaceId face = 0;
       face < topology.local_face_count(); ++face) {
    const auto area =
        geometry.face_area_vector_m2(face, hundun::mesh::FaceSide::owner);
    const double expected_flux =
        free_stream_velocity.x * area.x +
        free_stream_velocity.y * area.y +
        free_stream_velocity.z * area.z;
    HUNDUN_CHECK_NEAR(flux_read(face, 0), expected_flux,
                      face_tolerance(std::abs(expected_flux)));
    HUNDUN_CHECK_NEAR(momentum_face_read(face, 0), free_stream_velocity.x,
                      face_tolerance(std::abs(free_stream_velocity.x)));
    HUNDUN_CHECK_NEAR(momentum_face_read(face, 1), free_stream_velocity.y,
                      face_tolerance(std::abs(free_stream_velocity.y)));
    HUNDUN_CHECK_NEAR(momentum_face_read(face, 2), free_stream_velocity.z,
                      face_tolerance(std::abs(free_stream_velocity.z)));
  }

  const std::size_t cells = topology.owned_cell_count();
  const auto cell_index = [&](LocalCellId cell, std::size_t component,
                              std::size_t components) {
    return cell * components + component;
  };
  const auto independent_convective_oracle =
      [&](const std::vector<double> &constant_values) {
        const std::size_t components = constant_values.size();
        std::vector<double> residual(cells * components, 0.0);
        std::vector<double> magnitude(cells * components, 0.0);
        for (hundun::mesh::LocalFaceId face = 0;
             face < topology.local_face_count(); ++face) {
          const auto area = geometry.face_area_vector_m2(
              face, hundun::mesh::FaceSide::owner);
          const double analytic_flux =
              free_stream_velocity.x * area.x +
              free_stream_velocity.y * area.y +
              free_stream_velocity.z * area.z;
          const auto accumulate = [&](LocalCellId cell, double sign) {
            for (std::size_t component = 0; component < components;
                 ++component) {
              const auto index = cell_index(cell, component, components);
              const double contribution =
                  sign * analytic_flux * constant_values[component];
              residual[index] += contribution;
              magnitude[index] += std::abs(contribution);
            }
          };
          const auto owner = topology.owner(face);
          if (topology.cell_ownership(owner) ==
              hundun::mesh::EntityOwnership::owned)
            accumulate(owner, 1.0);
          const auto neighbour = topology.neighbour(face);
          if (!topology.periodic_pair(face).has_value() &&
              neighbour.has_value() &&
              topology.cell_ownership(*neighbour) ==
                  hundun::mesh::EntityOwnership::owned)
            accumulate(*neighbour, -1.0);
        }
        return std::pair{std::move(residual), std::move(magnitude)};
      };
  const auto verify_residual =
      [&](const hundun::runtime::FieldView<double> &product,
          const std::vector<double> &expected,
          const std::vector<double> &magnitude, std::size_t components) {
        std::vector<double> global_signed(components, 0.0);
        std::vector<double> global_magnitude(components, 0.0);
        for (LocalCellId cell = 0; cell < cells; ++cell) {
          const int i = static_cast<int>(
              cell % static_cast<std::size_t>(decomposition.local_extent().x));
          const std::size_t yz =
              cell /
              static_cast<std::size_t>(decomposition.local_extent().x);
          const int j = static_cast<int>(
              yz %
              static_cast<std::size_t>(decomposition.local_extent().y));
          const int k = static_cast<int>(
              yz /
              static_cast<std::size_t>(decomposition.local_extent().y));
          for (std::size_t component = 0; component < components;
               ++component) {
            const auto index = cell_index(cell, component, components);
            const double actual =
                product(i, j, k, static_cast<int>(component));
            const double tolerance =
                512.0 * eps * std::max(1.0, magnitude[index]);
            HUNDUN_CHECK_NEAR(actual, expected[index], tolerance);
            HUNDUN_CHECK(std::abs(actual) <= tolerance);
            global_signed[component] += actual;
            global_magnitude[component] += magnitude[index];
          }
        }
        HUNDUN_CHECK(MPI_Allreduce(
                         MPI_IN_PLACE, global_signed.data(),
                         static_cast<int>(components), MPI_DOUBLE, MPI_SUM,
                         mpi.comm()) == MPI_SUCCESS);
        HUNDUN_CHECK(MPI_Allreduce(
                         MPI_IN_PLACE, global_magnitude.data(),
                         static_cast<int>(components), MPI_DOUBLE, MPI_SUM,
                         mpi.comm()) == MPI_SUCCESS);
        for (std::size_t component = 0; component < components; ++component)
          HUNDUN_CHECK(
              std::abs(global_signed[component]) <=
              512.0 * eps * std::max(1.0, global_magnitude[component]));
      };
  const auto zero_residual =
      [&](const hundun::runtime::FieldView<double> &residual,
          std::size_t components) {
    for (int k = 0; k < decomposition.local_extent().z; ++k)
      for (int j = 0; j < decomposition.local_extent().y; ++j)
        for (int i = 0; i < decomposition.local_extent().x; ++i)
          for (std::size_t component = 0; component < components; ++component)
            residual(i, j, k, static_cast<int>(component)) = 0.0;
  };

  auto scalar_r = storage.view<double>(scalar_residual);
  const auto run_scalar_convective = [&](double constant_value,
                                         FiniteVolumeQuantity quantity) {
    for (int k = 0; k < decomposition.local_extent().z; ++k)
      for (int j = 0; j < decomposition.local_extent().y; ++j)
        for (int i = 0; i < decomposition.local_extent().x; ++i)
          q_write(i, j, k, 0) = constant_value;
    halo.exchange(storage, q);
    operators.reconstruct_transport_faces(quantity, boundaries, flux, q_read,
                                          transported_face_write);
    zero_residual(scalar_r, 1U);
    const auto transported_face_read =
        free_stream_read_storage.acquire_face_read<double>(
            access, read_phase, actor, transported_face);
    operators.accumulate_convective_residual(flux, transported_face_read,
                                             scalar_r);
    const auto [expected, magnitude] =
        independent_convective_oracle({constant_value});
    verify_residual(scalar_r, expected, magnitude, 1U);
  };
  run_scalar_convective(2.0, FiniteVolumeQuantity::enthalpy());
  run_scalar_convective(0.25, FiniteVolumeQuantity::scalar(0U));

  auto momentum_convective_r =
      storage.view<double>(momentum_convective_residual);
  zero_residual(momentum_convective_r, 3U);
  operators.accumulate_convective_residual(flux, momentum_face_read,
                                           momentum_convective_r);
  const auto [momentum_expected, momentum_magnitude] =
      independent_convective_oracle(
          {free_stream_velocity.x, free_stream_velocity.y,
           free_stream_velocity.z});
  verify_residual(momentum_convective_r, momentum_expected,
                  momentum_magnitude, 3U);

  auto momentum_viscous_r =
      storage.view<double>(momentum_viscous_residual);
  zero_residual(momentum_viscous_r, 3U);
  operators.accumulate_viscous_residual(
      boundaries, velocity_read,
      free_stream_read_storage.view<double>(velocity_gradient), 0.01,
      momentum_viscous_r);
  verify_residual(momentum_viscous_r, std::vector<double>(cells * 3U, 0.0),
                  std::vector<double>(cells * 3U, 0.0), 3U);

  for (hundun::mesh::LocalFaceId face = 0;
       face < topology.local_face_count(); ++face)
    gamma_write(face, 0) = 0.02;
  auto scalar_diffusive_r =
      storage.view<double>(scalar_diffusive_residual);
  zero_residual(scalar_diffusive_r, 1U);
  const auto gamma_read = free_stream_read_storage.acquire_face_read<double>(
      access, read_phase, actor, gamma);
  operators.compute_gradient(GradientScheme::weighted_least_squares,
                             FiniteVolumeQuantity::scalar(0U), boundaries,
                             q_read, gradient_write);
  operators.accumulate_scalar_diffusive_residual(
      FiniteVolumeQuantity::scalar(0U), boundaries, q_read,
      free_stream_read_storage.view<double>(gradient), gamma_read,
      scalar_diffusive_r);
  verify_residual(scalar_diffusive_r, std::vector<double>(cells, 0.0),
                  std::vector<double>(cells, 0.0), 1U);

  std::vector<std::array<double, 3>> closure(cells);
  std::vector<double> closure_scale(cells, 0.0);
  for (hundun::mesh::LocalFaceId face = 0;
       face < topology.local_face_count(); ++face) {
    const auto area =
        geometry.face_area_vector_m2(face, hundun::mesh::FaceSide::owner);
    const double norm = std::sqrt(area.x * area.x + area.y * area.y +
                                  area.z * area.z);
    const auto accumulate = [&](LocalCellId cell, double sign) {
      closure[cell][0] += sign * area.x;
      closure[cell][1] += sign * area.y;
      closure[cell][2] += sign * area.z;
      closure_scale[cell] += norm;
    };
    const auto owner = topology.owner(face);
    if (topology.cell_ownership(owner) ==
        hundun::mesh::EntityOwnership::owned)
      accumulate(owner, 1.0);
    const auto neighbour = topology.neighbour(face);
    if (!topology.periodic_pair(face).has_value() && neighbour.has_value() &&
        topology.cell_ownership(*neighbour) ==
            hundun::mesh::EntityOwnership::owned)
      accumulate(*neighbour, -1.0);
  }
  for (LocalCellId cell = 0; cell < cells; ++cell) {
    const double bound = 256.0 * eps * closure_scale[cell];
    for (const double component : closure[cell])
      HUNDUN_CHECK(std::abs(component) <= bound);
  }
  auto mutated_closure = closure;
  HUNDUN_CHECK(!mutated_closure.empty());
  mutated_closure[0][0] += 1.0;
  HUNDUN_CHECK(std::abs(mutated_closure[0][0]) >
               256.0 * eps * closure_scale[0]);
  if (mpi.rank() == 0)
    std::cout << "TASK25_WARPED_FREE_STREAM_FVM ranks=" << mpi.size()
              << " mutation_delta=1" << '\n';
}

void run_gradient_failure_contracts(const MpiContext &mpi, int ranks) {
  constexpr Int3 extent{6, 6, 6};
  auto decomposition = StructuredDecomposition::create(
      mpi, extent, std::array<bool, 3>{true, true, true},
      DecompositionOptions{process_grid_for(ranks)});
  MeshTopology topology(decomposition);
  MeshGeometry geometry(
      topology, UniformBoxMapping(Real3{0.0, 0.0, 0.0}, Real3{1.0, 1.0, 1.0}));
  BoundaryRegistry boundaries =
      BoundaryRegistry::create(periodic_case(), topology);
  CellCenteredFvmOperators operators =
      CellCenteredFvmOperators::create(topology, geometry);

  FieldRegistry registry;
  const auto q = registry.declare_field(cell_descriptor("failure_q", 1, 2));
  const auto q_narrow =
      registry.declare_field(cell_descriptor("failure_q_narrow", 1, 0));
  const auto q_vector =
      registry.declare_field(cell_descriptor("failure_q_vector", 3, 2));
  const auto gradient =
      registry.declare_field(cell_descriptor("failure_gradient", 3, 1));
  const auto gradient_vector =
      registry.declare_field(cell_descriptor("failure_gradient_vector", 9, 1));
  registry.freeze();
  FieldStorage storage(registry, FieldLayoutSet{decomposition.local_extent(),
                                                topology.local_face_count()});
  auto q_write = storage.view<double>(q);
  auto gradient_write = storage.view<double>(gradient);
  for (int k = -2; k < decomposition.local_extent().z + 2; ++k) {
    for (int j = -2; j < decomposition.local_extent().y + 2; ++j) {
      for (int i = -2; i < decomposition.local_extent().x + 2; ++i) {
        q_write(i, j, k, 0) = 2.0;
      }
    }
  }
  const auto reset_gradient = [&] {
    for (int k = 0; k < decomposition.local_extent().z; ++k) {
      for (int j = 0; j < decomposition.local_extent().y; ++j) {
        for (int i = 0; i < decomposition.local_extent().x; ++i) {
          for (int component = 0; component < 3; ++component) {
            gradient_write(i, j, k, component) = -17.0;
          }
        }
      }
    }
  };
  const auto require_gradient_unchanged = [&] {
    for (int k = 0; k < decomposition.local_extent().z; ++k) {
      for (int j = 0; j < decomposition.local_extent().y; ++j) {
        for (int i = 0; i < decomposition.local_extent().x; ++i) {
          for (int component = 0; component < 3; ++component) {
            HUNDUN_CHECK(bits(gradient_write(i, j, k, component)) ==
                         bits(-17.0));
          }
        }
      }
    }
  };
  const FieldStorage &read_storage = storage;
  const auto q_read = read_storage.view<double>(q);

  reset_gradient();
  expect_error([&] {
    operators.compute_gradient(static_cast<GradientScheme>(255U),
                               FiniteVolumeQuantity::density(), boundaries,
                               q_read, gradient_write);
  });
  require_gradient_unchanged();
  expect_error([&] {
    operators.compute_gradient(
        GradientScheme::green_gauss,
        FiniteVolumeQuantity{
            static_cast<hundun::finite_volume::FiniteVolumeQuantityKind>(255U),
            0U},
        boundaries, q_read, gradient_write);
  });
  require_gradient_unchanged();
  expect_error([&] {
    operators.compute_gradient(GradientScheme::green_gauss,
                               FiniteVolumeQuantity::scalar(1U), boundaries,
                               q_read, gradient_write);
  });
  require_gradient_unchanged();
  expect_error([&] {
    operators.compute_gradient(
        GradientScheme::green_gauss, FiniteVolumeQuantity::density(),
        boundaries, read_storage.view<double>(q_narrow), gradient_write);
  });
  require_gradient_unchanged();
  expect_error([&] {
    operators.compute_gradient(
        GradientScheme::green_gauss, FiniteVolumeQuantity::density(),
        boundaries, read_storage.view<double>(q_vector), gradient_write);
  });
  require_gradient_unchanged();
  expect_error([&] {
    operators.compute_gradient(GradientScheme::green_gauss,
                               FiniteVolumeQuantity::density(), boundaries,
                               q_read, storage.view<double>(gradient_vector));
  });
  require_gradient_unchanged();

  q_write(-1, 0, 0, 0) = std::numeric_limits<double>::quiet_NaN();
  expect_error([&] {
    operators.compute_gradient(GradientScheme::weighted_least_squares,
                               FiniteVolumeQuantity::density(), boundaries,
                               q_read, gradient_write);
  });
  require_gradient_unchanged();
  q_write(-1, 0, 0, 0) = 2.0;

  hundun::finite_volume::test::force_next_least_squares_singular();
  expect_error([&] {
    operators.compute_gradient(GradientScheme::weighted_least_squares,
                               FiniteVolumeQuantity::density(), boundaries,
                               q_read, gradient_write);
  });
  require_gradient_unchanged();

  FieldStorage wrong_extent_storage(
      registry, FieldLayoutSet{Int3{5, 6, 6}, topology.local_face_count()});
  const FieldStorage &wrong_extent_read = wrong_extent_storage;
  expect_error([&] {
    operators.compute_gradient(GradientScheme::green_gauss,
                               FiniteVolumeQuantity::density(), boundaries,
                               wrong_extent_read.view<double>(q),
                               wrong_extent_storage.view<double>(gradient));
  });

  FieldStorage stale_storage(registry,
                             FieldLayoutSet{decomposition.local_extent(),
                                            topology.local_face_count()});
  const FieldStorage &stale_read_storage = stale_storage;
  const auto stale_q = stale_read_storage.view<double>(q);
  auto stale_gradient = stale_storage.view<double>(gradient);
  stale_storage.begin_rebuild();
  expect_error([&] {
    operators.compute_gradient(GradientScheme::green_gauss,
                               FiniteVolumeQuantity::density(), boundaries,
                               stale_q, stale_gradient);
  });

  auto other_decomposition = StructuredDecomposition::create(
      mpi, Int3{7, 6, 6}, std::array<bool, 3>{true, true, true},
      DecompositionOptions{process_grid_for(ranks)});
  MeshTopology other_topology(other_decomposition);
  expect_error([&] {
    static_cast<void>(
        CellCenteredFvmOperators::create(other_topology, geometry));
  });
}

void run_shared_flux_case(const MpiContext &mpi, int ranks) {
  constexpr Int3 extent{12, 8, 6};
  const std::array<bool, 3> periodic{true, true, true};
  auto decomposition = StructuredDecomposition::create(
      mpi, extent, periodic, DecompositionOptions{process_grid_for(ranks)});
  MeshTopology topology(decomposition);
  MeshGeometry geometry(
      topology, UniformBoxMapping(Real3{0.0, 0.0, 0.0}, Real3{1.0, 1.0, 1.0}));
  BoundaryRegistry boundaries =
      BoundaryRegistry::create(periodic_case(), topology);
  CellCenteredFvmOperators operators =
      CellCenteredFvmOperators::create(topology, geometry);

  FieldRegistry registry;
  const FieldId density =
      registry.declare_field(cell_descriptor("density", 1, 2));
  const FieldId density_narrow =
      registry.declare_field(cell_descriptor("density_narrow", 1, 1));
  const FieldId velocity =
      registry.declare_field(cell_descriptor("velocity", 3, 2));
  const FieldId mass_flux = declare_face_mass_flux(registry);
  const FieldId density_face =
      registry.declare_field(face_descriptor("density_face", 1));
  const FieldId velocity_face =
      registry.declare_field(face_descriptor("velocity_face", 3));
  const FieldId gamma = registry.declare_field(face_descriptor("gamma", 1));
  const FieldId scalar_gradient =
      registry.declare_field(cell_descriptor("scalar_gradient", 3, 1));
  const FieldId velocity_gradient =
      registry.declare_field(cell_descriptor("velocity_gradient", 9, 1));
  const FieldId mass_residual =
      registry.declare_field(cell_descriptor("mass_residual", 1, 0));
  const FieldId scalar_residual =
      registry.declare_field(cell_descriptor("scalar_residual", 1, 0));
  const FieldId momentum_residual =
      registry.declare_field(cell_descriptor("momentum_residual", 3, 0));
  registry.freeze();
  require_face_mass_flux_field(registry, mass_flux);
  FieldStorage storage(registry, FieldLayoutSet{decomposition.local_extent(),
                                                topology.local_face_count()});
  auto rho = storage.view<double>(density);
  auto u = storage.view<double>(velocity);
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        const auto global = decomposition.global_cell(Int3{i, j, k});
        rho(i, j, k, 0) = 1.0 + 0.01 * global.x + 0.005 * global.y;
        const double gx = static_cast<double>(global.x);
        const double gy = static_cast<double>(global.y);
        const double gz = static_cast<double>(global.z);
        u(i, j, k, 0) = 0.5 + 0.001 * (gx * gx + 3.0 * gy + 5.0 * gz);
        u(i, j, k, 1) = 0.25 + 0.0015 * (2.0 * gx + gy * gy + 3.0 * gz);
        u(i, j, k, 2) = 0.125 + 0.002 * (gx + 2.0 * gy + gz * gz);
      }
    }
  }
  auto halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 2));
  halo.exchange(storage, density);
  halo.exchange(storage, velocity);

  constexpr PhaseId write_phase = 101U;
  constexpr PhaseId read_phase = 102U;
  constexpr ActorId actor = 201U;
  FieldAccessPlan access(registry);
  access.declare_access(write_phase, actor, mass_flux, AccessMode::write);
  access.declare_access(read_phase, actor, mass_flux, AccessMode::read);
  access.declare_access(read_phase, actor, density_face,
                        AccessMode::read_write);
  access.declare_access(read_phase, actor, velocity_face,
                        AccessMode::read_write);
  access.declare_access(read_phase, actor, gamma, AccessMode::read_write);
  access.freeze();
  const FieldStorage &read_storage = storage;
  const auto rho_read = read_storage.view<double>(density);
  const auto velocity_read = read_storage.view<double>(velocity);
  operators.assemble_provisional_mass_flux(boundaries, rho_read, velocity_read,
                                           registry, storage, access,
                                           write_phase, actor, mass_flux);
  FaceMassFlux flux = FaceMassFlux::acquire(
      registry, read_storage, access, read_phase, actor, mass_flux, topology);
  HUNDUN_CHECK(flux.field_id() == mass_flux);
  HUNDUN_CHECK(flux.face_count() == topology.local_face_count());
  auto rho_face = storage.acquire_face_write<double>(access, read_phase, actor,
                                                     density_face);
  auto u_face = storage.acquire_face_write<double>(access, read_phase, actor,
                                                   velocity_face);
  operators.reconstruct_transport_faces(FiniteVolumeQuantity::density(),
                                        boundaries, flux, rho_read, rho_face);
  std::optional<std::size_t> canonical_periodic_face;
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    const auto pair = topology.periodic_pair(face);
    if (pair.has_value() && topology.global_face_id(face) < *pair &&
        topology.find_local_face(*pair).has_value()) {
      canonical_periodic_face = face;
      break;
    }
  }
  HUNDUN_CHECK(canonical_periodic_face.has_value());
  hundun::finite_volume::test::override_next_face_metrics(
      topology.global_face_id(*canonical_periodic_face), 0.250001, 0.0);
  operators.reconstruct_momentum_faces(boundaries, flux, velocity_read, u_face);
  const auto flux_read = read_storage.acquire_face_read<double>(
      access, read_phase, actor, mass_flux);
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    HUNDUN_CHECK(std::isfinite(flux_read(face, 0)));
    HUNDUN_CHECK(std::isfinite(rho_face(face, 0)));
    HUNDUN_CHECK(rho_face(face, 0) > 0.0);
    HUNDUN_CHECK(std::isfinite(u_face(face, 0)));
    HUNDUN_CHECK(std::isfinite(u_face(face, 1)));
    HUNDUN_CHECK(std::isfinite(u_face(face, 2)));
    const auto pair = topology.periodic_pair(face);
    if (pair.has_value()) {
      const auto local_pair = topology.find_local_face(*pair);
      if (local_pair.has_value()) {
        HUNDUN_CHECK(bits(flux_read(face, 0)) ==
                     bits(-flux_read(*local_pair, 0)));
        HUNDUN_CHECK(bits(rho_face(face, 0)) == bits(rho_face(*local_pair, 0)));
        for (int component = 0; component < 3; ++component) {
          HUNDUN_CHECK(bits(u_face(face, component)) ==
                       bits(u_face(*local_pair, component)));
        }
      }
    }
  }

  struct FaceReplicaRecord final {
    std::uint64_t id;
    std::uint64_t flux;
    std::uint64_t transported;
    std::uint64_t velocity_x;
    std::uint64_t velocity_y;
    std::uint64_t velocity_z;
  };
  std::vector<FaceReplicaRecord> local_records;
  local_records.reserve(topology.local_face_count());
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    local_records.push_back(FaceReplicaRecord{
        topology.global_face_id(face), bits(flux_read(face, 0)),
        bits(rho_face(face, 0)), bits(u_face(face, 0)), bits(u_face(face, 1)),
        bits(u_face(face, 2))});
  }
  const int local_bytes =
      static_cast<int>(local_records.size() * sizeof(FaceReplicaRecord));
  std::vector<int> byte_counts(static_cast<std::size_t>(ranks));
  HUNDUN_CHECK(MPI_Allgather(&local_bytes, 1, MPI_INT, byte_counts.data(), 1,
                             MPI_INT, mpi.comm()) == MPI_SUCCESS);
  std::vector<int> byte_offsets(static_cast<std::size_t>(ranks), 0);
  int total_bytes = 0;
  for (int rank = 0; rank < ranks; ++rank) {
    byte_offsets[static_cast<std::size_t>(rank)] = total_bytes;
    total_bytes += byte_counts[static_cast<std::size_t>(rank)];
  }
  HUNDUN_CHECK(total_bytes % static_cast<int>(sizeof(FaceReplicaRecord)) == 0);
  std::vector<FaceReplicaRecord> gathered_records(
      static_cast<std::size_t>(total_bytes) / sizeof(FaceReplicaRecord));
  HUNDUN_CHECK(MPI_Allgatherv(local_records.data(), local_bytes, MPI_BYTE,
                              gathered_records.data(), byte_counts.data(),
                              byte_offsets.data(), MPI_BYTE,
                              mpi.comm()) == MPI_SUCCESS);
  std::sort(gathered_records.begin(), gathered_records.end(),
            [](const FaceReplicaRecord &lhs, const FaceReplicaRecord &rhs) {
              return lhs.id < rhs.id;
            });
  bool saw_partition_replica = false;
  for (std::size_t index = 1; index < gathered_records.size(); ++index) {
    if (gathered_records[index - 1U].id == gathered_records[index].id) {
      saw_partition_replica = true;
      HUNDUN_CHECK(gathered_records[index - 1U].flux ==
                   gathered_records[index].flux);
      HUNDUN_CHECK(gathered_records[index - 1U].transported ==
                   gathered_records[index].transported);
      HUNDUN_CHECK(gathered_records[index - 1U].velocity_x ==
                   gathered_records[index].velocity_x);
      HUNDUN_CHECK(gathered_records[index - 1U].velocity_y ==
                   gathered_records[index].velocity_y);
      HUNDUN_CHECK(gathered_records[index - 1U].velocity_z ==
                   gathered_records[index].velocity_z);
    }
  }
  HUNDUN_CHECK(saw_partition_replica == (ranks > 1));

  for (int k = -2; k < decomposition.local_extent().z + 2; ++k) {
    for (int j = -2; j < decomposition.local_extent().y + 2; ++j) {
      for (int i = -2; i < decomposition.local_extent().x + 2; ++i) {
        u(i, j, k, 0) = 0.0;
        u(i, j, k, 1) = 0.0;
        u(i, j, k, 2) = 0.0;
      }
    }
  }
  operators.assemble_provisional_mass_flux(boundaries, rho_read, velocity_read,
                                           registry, storage, access,
                                           write_phase, actor, mass_flux);
  operators.reconstruct_transport_faces(FiniteVolumeQuantity::density(),
                                        boundaries, flux, rho_read, rho_face);
  bool saw_local_zero_pair = false;
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    const auto pair = topology.periodic_pair(face);
    if (!pair.has_value() || topology.global_face_id(face) > *pair) {
      continue;
    }
    const auto local_pair = topology.find_local_face(*pair);
    if (!local_pair.has_value()) {
      continue;
    }
    saw_local_zero_pair = true;
    HUNDUN_CHECK(bits(flux_read(face, 0)) == bits(0.0));
    HUNDUN_CHECK(bits(flux_read(*local_pair, 0)) == bits(-0.0));
    HUNDUN_CHECK(bits(flux_read(face, 0)) == bits(-flux_read(*local_pair, 0)));
    HUNDUN_CHECK(bits(rho_face(face, 0)) == bits(rho_face(*local_pair, 0)));
  }
  HUNDUN_CHECK(saw_local_zero_pair);

  for (int k = -2; k < decomposition.local_extent().z + 2; ++k) {
    for (int j = -2; j < decomposition.local_extent().y + 2; ++j) {
      for (int i = -2; i < decomposition.local_extent().x + 2; ++i) {
        u(i, j, k, 0) = 0.75;
        u(i, j, k, 1) = -0.25;
        u(i, j, k, 2) = 0.125;
      }
    }
  }
  operators.assemble_provisional_mass_flux(boundaries, rho_read, velocity_read,
                                           registry, storage, access,
                                           write_phase, actor, mass_flux);

  expect_error([&] {
    operators.reconstruct_transport_faces(
        FiniteVolumeQuantity::density(), boundaries, flux,
        read_storage.view<double>(density_narrow), rho_face);
  });

  std::optional<std::size_t> threshold_face;
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    const auto neighbour = topology.neighbour(face);
    if (!topology.patch_id(face).has_value() && neighbour.has_value() &&
        topology.cell_ownership(topology.owner(face)) ==
            hundun::mesh::EntityOwnership::owned &&
        topology.cell_ownership(*neighbour) ==
            hundun::mesh::EntityOwnership::owned &&
        topology.logical_face(face).axis == hundun::mesh::FaceAxis::x) {
      const auto p = topology.global_cell(topology.owner(face));
      if (p.x >= 1 && p.x + 2 < extent.x) {
        threshold_face = face;
        break;
      }
    }
  }
  HUNDUN_CHECK(threshold_face.has_value());
  if (threshold_face.has_value()) {
    const auto selected_owner =
        topology.global_cell(topology.owner(*threshold_face));
    for (int k = 0; k < decomposition.local_extent().z; ++k) {
      for (int j = 0; j < decomposition.local_extent().y; ++j) {
        for (int i = 0; i < decomposition.local_extent().x; ++i) {
          const int gx = decomposition.global_cell(Int3{i, j, k}).x;
          double value = 2.0;
          if (gx == selected_owner.x - 1)
            value = 1.5;
          if (gx == selected_owner.x)
            value = 2.5;
          if (gx == selected_owner.x + 1)
            value = 1.5;
          if (gx == selected_owner.x + 2)
            value = 1.0;
          rho(i, j, k, 0) = value;
        }
      }
    }
    halo.exchange(storage, density);
    operators.reconstruct_transport_faces(FiniteVolumeQuantity::density(),
                                          boundaries, flux, rho_read, rho_face);
    HUNDUN_CHECK_NEAR(rho_face(*threshold_face, 0), 2.5, 0.0);
    for (int k = 0; k < decomposition.local_extent().z; ++k) {
      for (int j = 0; j < decomposition.local_extent().y; ++j) {
        for (int i = 0; i < decomposition.local_extent().x; ++i) {
          const auto global = decomposition.global_cell(Int3{i, j, k});
          rho(i, j, k, 0) = 1.0 + 0.1 * global.x;
        }
      }
    }
    halo.exchange(storage, density);
    operators.reconstruct_transport_faces(FiniteVolumeQuantity::density(),
                                          boundaries, flux, rho_read, rho_face);
    HUNDUN_CHECK_NEAR(rho_face(*threshold_face, 0),
                      1.0 + 0.1 * (selected_owner.x + 0.5), 1.0e-15);

    for (int k = 0; k < decomposition.local_extent().z; ++k) {
      for (int j = 0; j < decomposition.local_extent().y; ++j) {
        for (int i = 0; i < decomposition.local_extent().x; ++i) {
          const auto global = decomposition.global_cell(Int3{i, j, k});
          const double x = static_cast<double>(global.x);
          u(i, j, k, 0) = 1.0 + x * x + x + 1.0 / 3.0;
          u(i, j, k, 1) = 2.0;
          u(i, j, k, 2) = -3.0;
        }
      }
    }
    halo.exchange(storage, velocity);
    const auto face_id = topology.global_face_id(*threshold_face);
    hundun::finite_volume::test::override_next_face_metrics(face_id, 0.25,
                                                            70.0);
    operators.reconstruct_momentum_faces(boundaries, flux, velocity_read,
                                         u_face);
    const auto p = topology.global_cell(topology.owner(*threshold_face));
    const double x_p = static_cast<double>(p.x);
    const double q_p = 1.0 + x_p * x_p + x_p + 1.0 / 3.0;
    const double x_n = static_cast<double>(p.x + 1);
    const double q_n = 1.0 + x_n * x_n + x_n + 1.0 / 3.0;
    const double exact_face = 1.0 + x_n * x_n;
    HUNDUN_CHECK_NEAR(u_face(*threshold_face, 0), exact_face, 1.0e-14);
    hundun::finite_volume::test::override_next_face_metrics(face_id, 0.250001,
                                                            0.0);
    operators.reconstruct_momentum_faces(boundaries, flux, velocity_read,
                                         u_face);
    const double x_pm1 = static_cast<double>(p.x - 1);
    const double q_pm1 =
        1.0 + x_pm1 * x_pm1 + x_pm1 + 1.0 / 3.0;
    const double slope = monotonized_central(q_p - q_pm1, q_n - q_p);
    const double limited = q_p + 0.5 * slope;
    HUNDUN_CHECK_NEAR(u_face(*threshold_face, 0), limited, 0.0);
    HUNDUN_CHECK(limited != exact_face);
    hundun::finite_volume::test::override_next_face_metrics(face_id, 0.0,
                                                            70.000001);
    operators.reconstruct_momentum_faces(boundaries, flux, velocity_read,
                                         u_face);
    HUNDUN_CHECK_NEAR(u_face(*threshold_face, 0), limited, 0.0);

    for (int k = 0; k < decomposition.local_extent().z; ++k) {
      for (int j = 0; j < decomposition.local_extent().y; ++j) {
        for (int i = 0; i < decomposition.local_extent().x; ++i) {
          u(i, j, k, 0) = 0.75;
          u(i, j, k, 1) = -0.25;
          u(i, j, k, 2) = 0.125;
        }
      }
    }
    halo.exchange(storage, velocity);
    operators.reconstruct_momentum_faces(boundaries, flux, velocity_read,
                                         u_face);

    for (int k = 0; k < decomposition.local_extent().z; ++k) {
      for (int j = 0; j < decomposition.local_extent().y; ++j) {
        for (int i = 0; i < decomposition.local_extent().x; ++i) {
          const auto global = decomposition.global_cell(Int3{i, j, k});
          rho(i, j, k, 0) = 1.0 + 0.01 * global.x + 0.005 * global.y;
        }
      }
    }
    halo.exchange(storage, density);
    operators.reconstruct_transport_faces(FiniteVolumeQuantity::density(),
                                          boundaries, flux, rho_read, rho_face);
  }

  auto gamma_write =
      storage.acquire_face_write<double>(access, read_phase, actor, gamma);
  auto scalar_gradient_write = storage.view<double>(scalar_gradient);
  auto velocity_gradient_write = storage.view<double>(velocity_gradient);
  auto mass_r = storage.view<double>(mass_residual);
  auto scalar_r = storage.view<double>(scalar_residual);
  auto momentum_r = storage.view<double>(momentum_residual);
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    gamma_write(face, 0) = 0.0;
  }
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        mass_r(i, j, k, 0) = 3.0;
        scalar_r(i, j, k, 0) = -2.0;
        for (int component = 0; component < 3; ++component) {
          scalar_gradient_write(i, j, k, component) = 0.0;
          momentum_r(i, j, k, component) = 4.0;
        }
        for (int component = 0; component < 9; ++component) {
          velocity_gradient_write(i, j, k, component) = 0.0;
        }
      }
    }
  }
  std::vector<double> expected_mass(topology.owned_cell_count(), 3.0);
  std::vector<double> expected_scalar(topology.owned_cell_count(), -2.0);
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    const double m = flux_read(face, 0);
    const double transported = rho_face(face, 0);
    const auto owner = topology.owner(face);
    const auto neighbour = topology.neighbour(face);
    if (topology.cell_ownership(owner) ==
        hundun::mesh::EntityOwnership::owned) {
      expected_mass[owner] += m;
      expected_scalar[owner] += m * transported;
    }
    if (!topology.periodic_pair(face).has_value() && neighbour.has_value() &&
        topology.cell_ownership(*neighbour) ==
            hundun::mesh::EntityOwnership::owned) {
      expected_mass[*neighbour] -= m;
      expected_scalar[*neighbour] -= m * transported;
    }
  }
  operators.accumulate_mass_residual(flux, mass_r);
  const auto rho_face_read = read_storage.acquire_face_read<double>(
      access, read_phase, actor, density_face);
  operators.accumulate_convective_residual(flux, rho_face_read, scalar_r);
  const auto scalar_gradient_read = read_storage.view<double>(scalar_gradient);
  const auto velocity_gradient_read =
      read_storage.view<double>(velocity_gradient);
  const auto gamma_read =
      read_storage.acquire_face_read<double>(access, read_phase, actor, gamma);
  operators.accumulate_scalar_diffusive_residual(
      FiniteVolumeQuantity::density(), boundaries, rho_read,
      scalar_gradient_read, gamma_read, scalar_r);
  operators.accumulate_viscous_residual(
      boundaries, velocity_read, velocity_gradient_read, 1.25, momentum_r);
  const double rotation_tolerance =
      512.0 * std::numeric_limits<double>::epsilon();
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        const std::size_t cell =
            (static_cast<std::size_t>(k) *
                 static_cast<std::size_t>(decomposition.local_extent().y) +
             static_cast<std::size_t>(j)) *
                static_cast<std::size_t>(decomposition.local_extent().x) +
            static_cast<std::size_t>(i);
        HUNDUN_CHECK_NEAR(mass_r(i, j, k, 0), expected_mass[cell], 0.0);
        HUNDUN_CHECK_NEAR(scalar_r(i, j, k, 0), expected_scalar[cell], 0.0);
        HUNDUN_CHECK_NEAR(momentum_r(i, j, k, 0), 4.0, 0.0);
        HUNDUN_CHECK_NEAR(momentum_r(i, j, k, 1), 4.0, 0.0);
        HUNDUN_CHECK_NEAR(momentum_r(i, j, k, 2), 4.0, 0.0);
      }
    }
  }

  for (int k = -2; k < decomposition.local_extent().z + 2; ++k) {
    for (int j = -2; j < decomposition.local_extent().y + 2; ++j) {
      for (int i = -2; i < decomposition.local_extent().x + 2; ++i) {
        rho(i, j, k, 0) = 7.0;
      }
    }
  }
  for (int k = -1; k < decomposition.local_extent().z + 1; ++k) {
    for (int j = -1; j < decomposition.local_extent().y + 1; ++j) {
      for (int i = -1; i < decomposition.local_extent().x + 1; ++i) {
        for (int component = 0; component < 3; ++component) {
          scalar_gradient_write(i, j, k, component) = 0.0;
        }
      }
    }
  }
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    gamma_write(face, 0) = 1.0;
  }
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        scalar_r(i, j, k, 0) = 0.0;
      }
    }
  }
  operators.accumulate_scalar_diffusive_residual(
      FiniteVolumeQuantity::density(), boundaries, rho_read,
      scalar_gradient_read, gamma_read, scalar_r);
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        HUNDUN_CHECK_NEAR(scalar_r(i, j, k, 0), 0.0, 0.0);
      }
    }
  }

  const auto owned_box = topology.owned_global_box();
  for (int k = -2; k < decomposition.local_extent().z + 2; ++k) {
    for (int j = -2; j < decomposition.local_extent().y + 2; ++j) {
      for (int i = -2; i < decomposition.local_extent().x + 2; ++i) {
        const double x =
            (static_cast<double>(owned_box.begin.x + i) + 0.5) / extent.x;
        const double y =
            (static_cast<double>(owned_box.begin.y + j) + 0.5) / extent.y;
        u(i, j, k, 0) = -y;
        u(i, j, k, 1) = x;
        u(i, j, k, 2) = 0.0;
      }
    }
  }
  for (int k = -1; k < decomposition.local_extent().z + 1; ++k) {
    for (int j = -1; j < decomposition.local_extent().y + 1; ++j) {
      for (int i = -1; i < decomposition.local_extent().x + 1; ++i) {
        constexpr std::array<double, 9> rotation_gradient{
            0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        for (int component = 0; component < 9; ++component) {
          velocity_gradient_write(i, j, k, component) =
              rotation_gradient[static_cast<std::size_t>(component)];
        }
      }
    }
  }
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        for (int component = 0; component < 3; ++component) {
          momentum_r(i, j, k, component) = 0.0;
        }
      }
    }
  }
  operators.accumulate_viscous_residual(
      boundaries, velocity_read, velocity_gradient_read, 1.25, momentum_r);
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        HUNDUN_CHECK_NEAR(momentum_r(i, j, k, 0), 0.0, rotation_tolerance);
        HUNDUN_CHECK_NEAR(momentum_r(i, j, k, 1), 0.0, rotation_tolerance);
        HUNDUN_CHECK_NEAR(momentum_r(i, j, k, 2), 0.0, rotation_tolerance);
      }
    }
  }

  const auto counters_before = mpi.fp64_reduction_counters();
  {
    allocation_probe::AllocationAttemptGuard allocation_guard;
    operators.compute_gradient(GradientScheme::green_gauss,
                               FiniteVolumeQuantity::density(), boundaries,
                               rho_read, scalar_gradient_write);
    operators.reconstruct_transport_faces(FiniteVolumeQuantity::density(),
                                          boundaries, flux, rho_read, rho_face);
    operators.reconstruct_momentum_faces(boundaries, flux, velocity_read,
                                         u_face);
    operators.assemble_provisional_mass_flux(
        boundaries, rho_read, velocity_read, registry, storage, access,
        write_phase, actor, mass_flux);
    operators.accumulate_mass_residual(flux, mass_r);
    operators.accumulate_convective_residual(flux, rho_face_read, scalar_r);
    operators.accumulate_scalar_diffusive_residual(
        FiniteVolumeQuantity::density(), boundaries, rho_read,
        scalar_gradient_read, gamma_read, scalar_r);
    operators.accumulate_viscous_residual(
        boundaries, velocity_read, velocity_gradient_read, 1.25, momentum_r);
    HUNDUN_CHECK(allocation_guard.attempts() == 0U);
  }
  const auto counters_after = mpi.fp64_reduction_counters();
  HUNDUN_CHECK(counters_after.collective_calls ==
               counters_before.collective_calls);
  HUNDUN_CHECK(counters_after.reduced_scalars ==
               counters_before.reduced_scalars);
  HUNDUN_CHECK(counters_after.logical_payload_bytes ==
               counters_before.logical_payload_bytes);

  std::vector<std::uint64_t> scalar_snapshot(topology.owned_cell_count());
  std::vector<std::array<std::uint64_t, 3>> momentum_snapshot(
      topology.owned_cell_count());
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    scalar_snapshot[cell] = bits(scalar_r(i, j, k, 0));
    for (int component = 0; component < 3; ++component) {
      momentum_snapshot[cell][static_cast<std::size_t>(component)] =
          bits(momentum_r(i, j, k, component));
    }
  }
  gamma_write(0, 0) = -1.0;
  expect_error([&] {
    operators.accumulate_scalar_diffusive_residual(
        FiniteVolumeQuantity::density(), boundaries, rho_read,
        scalar_gradient_read, gamma_read, scalar_r);
  });
  gamma_write(0, 0) = 1.0;
  expect_error([&] {
    operators.accumulate_viscous_residual(
        boundaries, velocity_read, velocity_gradient_read, -1.0, momentum_r);
  });
  velocity_gradient_write(0, 0, 0, 0) =
      std::numeric_limits<double>::quiet_NaN();
  expect_error([&] {
    operators.accumulate_viscous_residual(
        boundaries, velocity_read, velocity_gradient_read, 1.25, momentum_r);
  });
  velocity_gradient_write(0, 0, 0, 0) = 0.0;
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    HUNDUN_CHECK(bits(scalar_r(i, j, k, 0)) == scalar_snapshot[cell]);
    for (int component = 0; component < 3; ++component) {
      HUNDUN_CHECK(
          bits(momentum_r(i, j, k, component)) ==
          momentum_snapshot[cell][static_cast<std::size_t>(component)]);
    }
  }

  const auto fill_move_sentinels = [&] {
    for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
      rho_face(face, 0) = 101.0;
      u_face(face, 0) = 103.0;
      u_face(face, 1) = 107.0;
      u_face(face, 2) = 109.0;
    }
    for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
      const int i = static_cast<int>(
          cell % static_cast<std::size_t>(decomposition.local_extent().x));
      const std::size_t yz =
          cell / static_cast<std::size_t>(decomposition.local_extent().x);
      const int j = static_cast<int>(
          yz % static_cast<std::size_t>(decomposition.local_extent().y));
      const int k = static_cast<int>(
          yz / static_cast<std::size_t>(decomposition.local_extent().y));
      scalar_gradient_write(i, j, k, 0) = 113.0;
      scalar_gradient_write(i, j, k, 1) = 127.0;
      scalar_gradient_write(i, j, k, 2) = 131.0;
      mass_r(i, j, k, 0) = 137.0;
      scalar_r(i, j, k, 0) = 139.0;
      momentum_r(i, j, k, 0) = 149.0;
      momentum_r(i, j, k, 1) = 151.0;
      momentum_r(i, j, k, 2) = 157.0;
    }
  };
  const auto require_move_sentinels = [&] {
    for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
      HUNDUN_CHECK(bits(rho_face(face, 0)) == bits(101.0));
      HUNDUN_CHECK(bits(u_face(face, 0)) == bits(103.0));
      HUNDUN_CHECK(bits(u_face(face, 1)) == bits(107.0));
      HUNDUN_CHECK(bits(u_face(face, 2)) == bits(109.0));
    }
    for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
      const int i = static_cast<int>(
          cell % static_cast<std::size_t>(decomposition.local_extent().x));
      const std::size_t yz =
          cell / static_cast<std::size_t>(decomposition.local_extent().x);
      const int j = static_cast<int>(
          yz % static_cast<std::size_t>(decomposition.local_extent().y));
      const int k = static_cast<int>(
          yz / static_cast<std::size_t>(decomposition.local_extent().y));
      HUNDUN_CHECK(bits(scalar_gradient_write(i, j, k, 0)) == bits(113.0));
      HUNDUN_CHECK(bits(scalar_gradient_write(i, j, k, 1)) == bits(127.0));
      HUNDUN_CHECK(bits(scalar_gradient_write(i, j, k, 2)) == bits(131.0));
      HUNDUN_CHECK(bits(mass_r(i, j, k, 0)) == bits(137.0));
      HUNDUN_CHECK(bits(scalar_r(i, j, k, 0)) == bits(139.0));
      HUNDUN_CHECK(bits(momentum_r(i, j, k, 0)) == bits(149.0));
      HUNDUN_CHECK(bits(momentum_r(i, j, k, 1)) == bits(151.0));
      HUNDUN_CHECK(bits(momentum_r(i, j, k, 2)) == bits(157.0));
    }
  };

  const auto moved_field = flux.field_id();
  const auto moved_face_count = flux.face_count();
  FaceMassFlux moved_flux(std::move(flux));
  HUNDUN_CHECK(flux.field_id() == moved_field);
  HUNDUN_CHECK(flux.face_count() == moved_face_count);
  HUNDUN_CHECK(moved_flux.field_id() == moved_field);
  HUNDUN_CHECK(moved_flux.face_count() == moved_face_count);
  fill_move_sentinels();
  expect_error([&] {
    operators.reconstruct_transport_faces(FiniteVolumeQuantity::density(),
                                          boundaries, flux, rho_read, rho_face);
  });
  expect_error([&] {
    operators.reconstruct_momentum_faces(boundaries, flux, velocity_read,
                                         u_face);
  });
  expect_error([&] { operators.accumulate_mass_residual(flux, mass_r); });
  expect_error([&] {
    operators.accumulate_convective_residual(flux, rho_face_read, scalar_r);
  });
  require_move_sentinels();
  operators.reconstruct_transport_faces(FiniteVolumeQuantity::density(),
                                        boundaries, moved_flux, rho_read,
                                        rho_face);
  HUNDUN_CHECK(std::isfinite(rho_face(0, 0)));

  std::vector<std::uint64_t> move_flux_snapshot(topology.local_face_count());
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    move_flux_snapshot[face] = bits(flux_read(face, 0));
  }
  fill_move_sentinels();
  CellCenteredFvmOperators moved_operators(std::move(operators));
  expect_error([&] {
    operators.compute_gradient(GradientScheme::green_gauss,
                               FiniteVolumeQuantity::density(), boundaries,
                               rho_read, scalar_gradient_write);
  });
  expect_error([&] {
    operators.reconstruct_transport_faces(FiniteVolumeQuantity::density(),
                                          boundaries, moved_flux, rho_read,
                                          rho_face);
  });
  expect_error([&] {
    operators.reconstruct_momentum_faces(boundaries, moved_flux, velocity_read,
                                         u_face);
  });
  expect_error([&] {
    operators.assemble_provisional_mass_flux(
        boundaries, rho_read, velocity_read, registry, storage, access,
        write_phase, actor, mass_flux);
  });
  expect_error([&] { operators.accumulate_mass_residual(moved_flux, mass_r); });
  expect_error([&] {
    operators.accumulate_convective_residual(moved_flux, rho_face_read,
                                             scalar_r);
  });
  expect_error([&] {
    operators.accumulate_scalar_diffusive_residual(
        FiniteVolumeQuantity::density(), boundaries, rho_read,
        scalar_gradient_read, gamma_read, scalar_r);
  });
  expect_error([&] {
    operators.accumulate_viscous_residual(
        boundaries, velocity_read, velocity_gradient_read, 1.25, momentum_r);
  });
  require_move_sentinels();
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    HUNDUN_CHECK(bits(flux_read(face, 0)) == move_flux_snapshot[face]);
  }
  moved_operators.reconstruct_transport_faces(FiniteVolumeQuantity::density(),
                                              boundaries, moved_flux, rho_read,
                                              rho_face);
  HUNDUN_CHECK(std::isfinite(rho_face(0, 0)));
}

double smooth_transport_reconstruction_error(const MpiContext &mpi, int ranks,
                                             int n) {
  const Int3 extent{n, n, n};
  auto decomposition = StructuredDecomposition::create(
      mpi, extent, std::array<bool, 3>{true, true, true},
      DecompositionOptions{process_grid_for(ranks)});
  MeshTopology topology(decomposition);
  MeshGeometry geometry(
      topology, UniformBoxMapping(Real3{0.0, 0.0, 0.0}, Real3{1.0, 1.0, 1.0}));
  BoundaryRegistry boundaries =
      BoundaryRegistry::create(periodic_case(), topology);
  CellCenteredFvmOperators operators =
      CellCenteredFvmOperators::create(topology, geometry);
  FieldRegistry registry;
  const auto density =
      registry.declare_field(cell_descriptor("smooth_density", 1, 2));
  const auto velocity =
      registry.declare_field(cell_descriptor("smooth_velocity", 3, 2));
  const auto flux_id = declare_face_mass_flux(registry);
  const auto face_density =
      registry.declare_field(face_descriptor("smooth_face_density", 1));
  registry.freeze();
  FieldStorage storage(registry, FieldLayoutSet{decomposition.local_extent(),
                                                topology.local_face_count()});
  auto rho = storage.view<double>(density);
  auto u = storage.view<double>(velocity);
  constexpr double two_pi = 6.283185307179586476925286766559;
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    const Real3 center = geometry.cell_center_m(cell);
    rho(i, j, k, 0) = 2.0 + 0.2 * std::sin(two_pi * center.x);
    u(i, j, k, 0) = 1.0;
    u(i, j, k, 1) = 0.0;
    u(i, j, k, 2) = 0.0;
  }
  auto halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 2));
  halo.exchange(storage, density);
  halo.exchange(storage, velocity);
  constexpr PhaseId write_phase = 211U;
  constexpr PhaseId read_phase = 212U;
  constexpr ActorId actor = 213U;
  FieldAccessPlan access(registry);
  access.declare_access(write_phase, actor, flux_id, AccessMode::write);
  access.declare_access(read_phase, actor, flux_id, AccessMode::read);
  access.declare_access(read_phase, actor, face_density,
                        AccessMode::read_write);
  access.freeze();
  const FieldStorage &read_storage = storage;
  const auto rho_read = read_storage.view<double>(density);
  operators.assemble_provisional_mass_flux(
      boundaries, rho_read, read_storage.view<double>(velocity), registry,
      storage, access, write_phase, actor, flux_id);
  const auto flux = FaceMassFlux::acquire(registry, read_storage, access,
                                          read_phase, actor, flux_id, topology);
  auto rho_face = storage.acquire_face_write<double>(access, read_phase, actor,
                                                     face_density);
  operators.reconstruct_transport_faces(FiniteVolumeQuantity::density(),
                                        boundaries, flux, rho_read, rho_face);
  double local_squared_error = 0.0;
  double local_count = 0.0;
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    if (topology.face_ownership(face) != hundun::mesh::EntityOwnership::owned ||
        topology.logical_face(face).axis != hundun::mesh::FaceAxis::x) {
      continue;
    }
    const double exact =
        2.0 + 0.2 * std::sin(two_pi * geometry.face_center_m(face).x);
    const double error = rho_face(face, 0) - exact;
    local_squared_error += error * error;
    local_count += 1.0;
  }
  double totals[2]{local_squared_error, local_count};
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, totals, 2, MPI_DOUBLE, MPI_SUM,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(totals[1] > 0.0);
  return std::sqrt(totals[0] / totals[1]);
}

void run_uniform_quadratic_diffusion_case(const MpiContext &mpi, int ranks) {
  constexpr Int3 extent{8, 8, 8};
  auto decomposition = StructuredDecomposition::create(
      mpi, extent, std::array<bool, 3>{true, true, true},
      DecompositionOptions{process_grid_for(ranks)});
  MeshTopology topology(decomposition);
  MeshGeometry geometry(
      topology, UniformBoxMapping(Real3{0.0, 0.0, 0.0}, Real3{1.0, 1.0, 1.0}));
  BoundaryRegistry boundaries =
      BoundaryRegistry::create(periodic_case(), topology);
  CellCenteredFvmOperators operators =
      CellCenteredFvmOperators::create(topology, geometry);
  FieldRegistry registry;
  const auto q = registry.declare_field(cell_descriptor("quadratic_q", 1, 2));
  const auto gradient =
      registry.declare_field(cell_descriptor("quadratic_gradient", 3, 1));
  const auto residual =
      registry.declare_field(cell_descriptor("quadratic_residual", 1, 0));
  const auto gamma =
      registry.declare_field(face_descriptor("quadratic_gamma", 1));
  registry.freeze();
  FieldStorage storage(registry, FieldLayoutSet{decomposition.local_extent(),
                                                topology.local_face_count()});
  auto q_write = storage.view<double>(q);
  auto gradient_write = storage.view<double>(gradient);
  auto residual_write = storage.view<double>(residual);
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    const Real3 center = geometry.cell_center_m(cell);
    q_write(i, j, k, 0) =
        center.x * center.x + center.y * center.y + center.z * center.z;
    gradient_write(i, j, k, 0) = 2.0 * center.x;
    gradient_write(i, j, k, 1) = 2.0 * center.y;
    gradient_write(i, j, k, 2) = 2.0 * center.z;
    residual_write(i, j, k, 0) = 0.0;
  }
  auto q_halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 2));
  q_halo.exchange(storage, q);
  auto gradient_halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 1));
  gradient_halo.exchange(storage, gradient);
  constexpr PhaseId phase = 311U;
  constexpr ActorId actor = 312U;
  FieldAccessPlan access(registry);
  access.declare_access(phase, actor, gamma, AccessMode::read_write);
  access.freeze();
  auto gamma_write =
      storage.acquire_face_write<double>(access, phase, actor, gamma);
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    gamma_write(face, 0) = 1.0;
  }
  const FieldStorage &read_storage = storage;
  operators.accumulate_scalar_diffusive_residual(
      FiniteVolumeQuantity::density(), boundaries, read_storage.view<double>(q),
      read_storage.view<double>(gradient),
      read_storage.acquire_face_read<double>(access, phase, actor, gamma),
      residual_write);
  const auto owned = topology.owned_global_box();
  bool checked = false;
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    const Int3 global{owned.begin.x + i, owned.begin.y + j, owned.begin.z + k};
    if (global.x == 0 || global.x + 1 == extent.x || global.y == 0 ||
        global.y + 1 == extent.y || global.z == 0 || global.z + 1 == extent.z) {
      continue;
    }
    checked = true;
    const double expected = -6.0 * geometry.cell_volume_m3(cell);
    const double tolerance = 512.0 * std::numeric_limits<double>::epsilon() *
                             std::max(1.0, std::abs(expected));
    HUNDUN_CHECK_NEAR(residual_write(i, j, k, 0), expected, tolerance);
  }
  int checked_int = checked ? 1 : 0;
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &checked_int, 1, MPI_INT, MPI_MAX,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(checked_int == 1);
}

struct WarpedDiffusionResult {
  double error{};
  double integral_defect{};
  std::vector<double> physical_laplacian;
  std::vector<hundun::mesh::GlobalCellId> global_ids;
};

WarpedDiffusionResult warped_diffusion_error(const MpiContext &mpi, int ranks,
                                             int n,
                                             bool reverse_nonorthogonal =
                                                 false) {
  const Int3 extent{n, n, n};
  const std::array<bool, 3> periodic{true, true, true};
  auto decomposition = StructuredDecomposition::create(
      mpi, extent, periodic, DecompositionOptions{process_grid_for(ranks)});
  MeshTopology topology(decomposition);
  MeshGeometry geometry(topology,
                        AnalyticWarpedBoxMapping(Real3{0.0, 0.0, 0.0},
                                                 Real3{1.0, 1.0, 1.0},
                                                 Real3{0.02, -0.015, 0.01}));
  BoundaryRegistry boundaries =
      BoundaryRegistry::create(periodic_case(), topology);
  CellCenteredFvmOperators operators =
      CellCenteredFvmOperators::create(topology, geometry);
  FieldRegistry registry;
  const FieldId q = registry.declare_field(cell_descriptor("q_diff", 1, 2));
  const FieldId gradient =
      registry.declare_field(cell_descriptor("gradient_diff", 3, 1));
  const FieldId residual =
      registry.declare_field(cell_descriptor("residual_diff", 1, 0));
  const FieldId gamma =
      registry.declare_field(face_descriptor("gamma_diff", 1));
  registry.freeze();
  FieldStorage storage(registry, FieldLayoutSet{decomposition.local_extent(),
                                                topology.local_face_count()});
  auto q_write = storage.view<double>(q);
  auto gradient_write = storage.view<double>(gradient);
  auto residual_write = storage.view<double>(residual);
  constexpr double two_pi = 6.283185307179586476925286766559;
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    const Real3 center = geometry.cell_center_m(cell);
    q_write(i, j, k, 0) = std::sin(two_pi * center.x);
    gradient_write(i, j, k, 0) = -999.0;
    gradient_write(i, j, k, 1) = -999.0;
    gradient_write(i, j, k, 2) = -999.0;
    residual_write(i, j, k, 0) = 0.0;
  }
  auto q_halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 2));
  q_halo.exchange(storage, q);
  const FieldStorage &gradient_source_storage = storage;
  operators.compute_gradient(
      GradientScheme::weighted_least_squares,
      FiniteVolumeQuantity::density(), boundaries,
      gradient_source_storage.view<double>(q), gradient_write);
  auto gradient_halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 1));
  gradient_halo.exchange(storage, gradient);
  constexpr PhaseId phase = 301U;
  constexpr ActorId actor = 302U;
  FieldAccessPlan access(registry);
  access.declare_access(phase, actor, gamma, AccessMode::read_write);
  access.freeze();
  auto gamma_write =
      storage.acquire_face_write<double>(access, phase, actor, gamma);
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    gamma_write(face, 0) = 1.0;
  }
  const FieldStorage &read_storage = storage;
  const auto q_read = read_storage.view<double>(q);
  const auto gradient_read = read_storage.view<double>(gradient);
  const auto gamma_read =
      read_storage.acquire_face_read<double>(access, phase, actor, gamma);
  struct NonorthogonalSignGuard final {
    explicit NonorthogonalSignGuard(bool reverse) : active(reverse) {
      if (active) {
        hundun::finite_volume::test::
            reverse_scalar_diffusion_nonorthogonal_contribution_for_test(true);
      }
    }
    ~NonorthogonalSignGuard() {
      if (active) {
        hundun::finite_volume::test::
            reverse_scalar_diffusion_nonorthogonal_contribution_for_test(
                false);
      }
    }
    bool active;
  } sign_guard(reverse_nonorthogonal);
  operators.accumulate_scalar_diffusive_residual(
      FiniteVolumeQuantity::density(), boundaries, q_read, gradient_read,
      gamma_read, residual_write);
  double local_error = 0.0;
  double local_volume = 0.0;
  double local_signed_residual = 0.0;
  double local_absolute_residual = 0.0;
  std::vector<double> physical_laplacian;
  std::vector<hundun::mesh::GlobalCellId> global_ids;
  physical_laplacian.reserve(topology.owned_cell_count());
  global_ids.reserve(topology.owned_cell_count());
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    const Real3 center = geometry.cell_center_m(cell);
    const double volume = geometry.cell_volume_m3(cell);
    const double analytic_negative_laplacian =
        two_pi * two_pi * std::sin(two_pi * center.x);
    const double analytic_physical_laplacian =
        -two_pi * two_pi * std::sin(two_pi * center.x);
    const double discrete_negative_laplacian =
        residual_write(i, j, k, 0) / volume;
    const double discrete_physical_laplacian =
        -residual_write(i, j, k, 0) / volume;
    const double error =
        discrete_negative_laplacian - analytic_negative_laplacian;
    const double physical_error =
        discrete_physical_laplacian - analytic_physical_laplacian;
    HUNDUN_CHECK_NEAR(
        error, -physical_error,
        64.0 * std::numeric_limits<double>::epsilon() *
            std::max({1.0, std::abs(error), std::abs(physical_error)}));
    local_error += error * error * volume;
    local_volume += volume;
    const double integrated = residual_write(i, j, k, 0);
    physical_laplacian.push_back(discrete_physical_laplacian);
    global_ids.push_back(topology.global_cell_id(cell));
    local_signed_residual += integrated;
    local_absolute_residual += std::abs(integrated);
  }
  double values[4]{local_error, local_volume, local_signed_residual,
                   local_absolute_residual};
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, values, 4, MPI_DOUBLE, MPI_SUM,
                             mpi.comm()) == MPI_SUCCESS);
  const double defect = std::abs(values[2]);
  HUNDUN_CHECK(
      defect <= 512.0 * std::numeric_limits<double>::epsilon() *
                    std::max(1.0, values[3]));
  return {std::sqrt(values[0] / values[1]), defect,
          std::move(physical_laplacian), std::move(global_ids)};
}

double reject_reversed_nonorthogonal_diffusion(
    const MpiContext &mpi, int ranks,
    const WarpedDiffusionResult &normal_error16) {
  const auto reversed16 = warped_diffusion_error(mpi, ranks, 16, true);
  const auto reversed32 = warped_diffusion_error(mpi, ranks, 32, true);
  const double reversed_order =
      std::log(reversed16.error / reversed32.error) / std::log(2.0);
  HUNDUN_CHECK(reversed_order < 1.8);
  HUNDUN_CHECK(reversed16.error > normal_error16.error);
  return reversed_order;
}

double no_slip_wall_traction_error(const MpiContext &mpi, int ranks, int cells,
                                   bool warped, double cubic_coefficient) {
  const Int3 extent{cells, cells, cells};
  auto decomposition = StructuredDecomposition::create(
      mpi, extent, std::array<bool, 3>{false, false, false},
      DecompositionOptions{process_grid_for(ranks)});
  MeshTopology topology(decomposition);
  MeshGeometry geometry =
      warped ? MeshGeometry(topology,
                            AnalyticWarpedBoxMapping(Real3{0.0, 0.0, 0.0},
                                                     Real3{1.0, 1.0, 1.0},
                                                     Real3{0.02, -0.015, 0.01}))
             : MeshGeometry(topology, UniformBoxMapping(Real3{0.0, 0.0, 0.0},
                                                        Real3{1.0, 1.0, 1.0}));
  BoundaryRegistry boundaries =
      BoundaryRegistry::create(all_no_slip_case(), topology);
  CellCenteredFvmOperators operators =
      CellCenteredFvmOperators::create(topology, geometry);

  FieldRegistry registry;
  const auto velocity =
      registry.declare_field(cell_descriptor("wall_order_velocity", 3, 2));
  const auto gradient = registry.declare_field(
      cell_descriptor("wall_order_velocity_gradient", 9, 1));
  const auto flux_id = declare_face_mass_flux(registry);
  const auto face_velocity =
      registry.declare_field(face_descriptor("wall_order_face_velocity", 3));
  registry.freeze();
  FieldStorage storage(registry, FieldLayoutSet{decomposition.local_extent(),
                                                topology.local_face_count()});
  auto velocity_write = storage.view<double>(velocity);
  auto gradient_write = storage.view<double>(gradient);
  constexpr double linear_coefficient = 0.75;
  constexpr double quadratic_coefficient = 1.25;
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    const auto moments = cell_y_moments(topology.global_cell(cell), geometry);
    const double y_average = moments.y_integral / moments.volume;
    const double y_squared_average =
        moments.y_squared_integral / moments.volume;
    const double y_cubed_average =
        moments.y_cubed_integral / moments.volume;
    velocity_write(i, j, k, 0) =
        linear_coefficient * y_average +
        quadratic_coefficient * y_squared_average +
        cubic_coefficient * y_cubed_average;
    velocity_write(i, j, k, 1) = 0.0;
    velocity_write(i, j, k, 2) = 0.0;
    for (int component = 0; component < 9; ++component)
      gradient_write(i, j, k, component) = 0.0;
    gradient_write(i, j, k, 1) =
        linear_coefficient + 2.0 * quadratic_coefficient * y_average +
        3.0 * cubic_coefficient * y_squared_average;
  }
  auto velocity_halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 2));
  velocity_halo.exchange(storage, velocity);
  auto gradient_halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 1));
  gradient_halo.exchange(storage, gradient);

  constexpr PhaseId phase = 351U;
  constexpr ActorId actor = 352U;
  FieldAccessPlan access(registry);
  access.declare_access(phase, actor, flux_id, AccessMode::read_write);
  access.declare_access(phase, actor, face_velocity, AccessMode::read_write);
  access.freeze();
  auto flux_write =
      storage.acquire_face_write<double>(access, phase, actor, flux_id);
  auto face_velocity_write = storage.acquire_face_write<double>(
      access, phase, actor, face_velocity);
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    flux_write(face, 0) = 0.0;
    for (int component = 0; component < 3; ++component)
      face_velocity_write(face, component) = 0.0;
  }

  const FieldStorage &read_storage = storage;
  const auto flux = FaceMassFlux::acquire(
      registry, read_storage, access, phase, actor, flux_id, topology);
  std::vector<hundun::finite_volume::PhysicalBoundaryMomentumContribution>
      contributions;
  operators.physical_boundary_momentum_contributions(
      boundaries, flux,
      read_storage.acquire_face_read<double>(access, phase, actor,
                                              face_velocity),
      read_storage.view<double>(velocity), read_storage.view<double>(gradient),
      1.0, contributions);
  double local_squared_error = 0.0;
  double local_count = 0.0;
  for (const auto &contribution : contributions) {
    const auto face = topology.find_local_face(contribution.global_face_id);
    HUNDUN_CHECK(face.has_value());
    if (topology.patch_id(*face) != std::optional<std::uint32_t>{2U} ||
        topology.face_ownership(*face) !=
            hundun::mesh::EntityOwnership::owned)
      continue;
    const double area = geometry.face_area_m2(*face);
    HUNDUN_CHECK(area > 0.0 && std::isfinite(area));
    const double recovered_derivative = contribution.viscous[0] / area;
    const double error = recovered_derivative - linear_coefficient;
    local_squared_error += error * error;
    local_count += 1.0;
  }
  double totals[2]{local_squared_error, local_count};
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, totals, 2, MPI_DOUBLE, MPI_SUM,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(totals[1] > 0.0);
  return std::sqrt(totals[0] / totals[1]);
}

void run_no_slip_quadratic_traction_order(const MpiContext &mpi, int ranks) {
  for (const bool warped : {false, true}) {
    const double exact_quadratic_error =
        no_slip_wall_traction_error(mpi, ranks, 8, warped, 0.0);
    HUNDUN_CHECK(exact_quadratic_error <=
                 32768.0 * std::numeric_limits<double>::epsilon());
    const double error8 =
        no_slip_wall_traction_error(mpi, ranks, 8, warped, 0.5);
    const double error16 =
        no_slip_wall_traction_error(mpi, ranks, 16, warped, 0.5);
    const double error32 =
        no_slip_wall_traction_error(mpi, ranks, 32, warped, 0.5);
    HUNDUN_CHECK(error8 > error16 && error16 > error32);
    HUNDUN_CHECK(std::log(error8 / error16) / std::log(2.0) >= 1.8);
    HUNDUN_CHECK(std::log(error16 / error32) / std::log(2.0) >= 1.8);
  }
}

void run_physical_boundary_case(const MpiContext &mpi, int ranks, bool warped) {
  constexpr Int3 extent{8, 8, 4};
  const std::array<bool, 3> periodic{false, false, true};
  auto decomposition = StructuredDecomposition::create(
      mpi, extent, periodic, DecompositionOptions{process_grid_for(ranks)});
  MeshTopology topology(decomposition);
  MeshGeometry geometry =
      warped ? MeshGeometry(topology,
                            AnalyticWarpedBoxMapping(Real3{0.0, 0.0, 0.0},
                                                     Real3{1.0, 1.0, 1.0},
                                                     Real3{0.02, -0.015, 0.01}))
             : MeshGeometry(topology, UniformBoxMapping(Real3{0.0, 0.0, 0.0},
                                                        Real3{1.0, 1.0, 1.0}));
  BoundaryRegistry boundaries =
      BoundaryRegistry::create(mixed_case(), topology);
  CellCenteredFvmOperators operators =
      CellCenteredFvmOperators::create(topology, geometry);
  const auto expected_descriptor =
      hundun::finite_volume::face_mass_flux_descriptor();
  HUNDUN_CHECK(expected_descriptor.name == "face_mass_flux");
  HUNDUN_CHECK(expected_descriptor.unit == "kg/s");
  HUNDUN_CHECK(expected_descriptor.owner == "flow");
  HUNDUN_CHECK(expected_descriptor.space == FunctionSpace::face_value);
  HUNDUN_CHECK(expected_descriptor.scalar_type == ScalarType::float64);
  HUNDUN_CHECK(expected_descriptor.components == 1U);
  HUNDUN_CHECK(expected_descriptor.ghost_width == 0);
  HUNDUN_CHECK(expected_descriptor.conservative);
  HUNDUN_CHECK(expected_descriptor.restart == RestartPolicy::persistent);
  HUNDUN_CHECK(expected_descriptor.output == OutputPolicy::never);

  FieldRegistry registry;
  const FieldId density =
      registry.declare_field(cell_descriptor("physical_density", 1, 2));
  const FieldId velocity =
      registry.declare_field(cell_descriptor("physical_velocity", 3, 2));
  const FieldId pressure =
      registry.declare_field(cell_descriptor("physical_pressure", 1, 2));
  const FieldId flux_id = declare_face_mass_flux(registry);
  const FieldId density_face =
      registry.declare_field(face_descriptor("physical_density_face", 1));
  const FieldId velocity_face =
      registry.declare_field(face_descriptor("physical_velocity_face", 3));
  const FieldId density_gradient = registry.declare_field(
      cell_descriptor("physical_density_gradient", 3, 1));
  const FieldId velocity_gradient = registry.declare_field(
      cell_descriptor("physical_velocity_gradient", 9, 1));
  const FieldId momentum_residual = registry.declare_field(
      cell_descriptor("physical_momentum_residual", 3, 0));
  const FieldId scalar_residual =
      registry.declare_field(cell_descriptor("physical_scalar_residual", 1, 0));
  const FieldId gamma =
      registry.declare_field(face_descriptor("physical_gamma", 1));
  registry.freeze();
  FieldStorage storage(registry, FieldLayoutSet{decomposition.local_extent(),
                                                topology.local_face_count()});
  auto rho = storage.view<double>(density);
  auto u = storage.view<double>(velocity);
  auto p = storage.view<double>(pressure);
  auto density_gradient_write = storage.view<double>(density_gradient);
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        const LocalCellId cell =
            (static_cast<std::size_t>(k) *
                 static_cast<std::size_t>(decomposition.local_extent().y) +
             static_cast<std::size_t>(j)) *
                static_cast<std::size_t>(decomposition.local_extent().x) +
            static_cast<std::size_t>(i);
        rho(i, j, k, 0) = 1.5 + geometry.cell_center_m(cell).x;
        p(i, j, k, 0) = 2.0;
        u(i, j, k, 0) = 1.0;
        u(i, j, k, 1) = 0.5;
        u(i, j, k, 2) = -0.25;
      }
    }
  }
  auto halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 2));
  halo.exchange(storage, density);
  halo.exchange(storage, velocity);
  constexpr PhaseId write_phase = 401U;
  constexpr PhaseId read_phase = 402U;
  constexpr ActorId actor = 403U;
  FieldAccessPlan access(registry);
  access.declare_access(write_phase, actor, flux_id, AccessMode::write);
  access.declare_access(read_phase, actor, flux_id, AccessMode::read);
  access.declare_access(read_phase, actor, density_face,
                        AccessMode::read_write);
  access.declare_access(read_phase, actor, velocity_face,
                        AccessMode::read_write);
  access.declare_access(read_phase, actor, gamma, AccessMode::read_write);
  access.freeze();
  const FieldStorage &read_storage = storage;
  const auto rho_read = read_storage.view<double>(density);
  const auto u_read = read_storage.view<double>(velocity);
  for (const GradientScheme scheme :
       {GradientScheme::green_gauss, GradientScheme::weighted_least_squares}) {
    operators.compute_gradient(scheme, FiniteVolumeQuantity::density(),
                               boundaries, rho_read, density_gradient_write);
    const auto owned = topology.owned_global_box();
    const double tolerance =
        512.0 * std::numeric_limits<double>::epsilon() * 2.5;
    for (int k = 0; k < decomposition.local_extent().z; ++k) {
      for (int j = 0; j < decomposition.local_extent().y; ++j) {
        for (int i = 0; i < decomposition.local_extent().x; ++i) {
          const int global_y = owned.begin.y + j;
          if (owned.begin.x + i == 0 && global_y > 0 &&
              global_y + 1 < extent.y &&
              (!warped || scheme == GradientScheme::weighted_least_squares)) {
            HUNDUN_CHECK_NEAR(density_gradient_write(i, j, k, 0), 1.0,
                              tolerance);
            HUNDUN_CHECK_NEAR(density_gradient_write(i, j, k, 1), 0.0,
                              tolerance);
            HUNDUN_CHECK_NEAR(density_gradient_write(i, j, k, 2), 0.0,
                              tolerance);
          }
        }
      }
    }
  }
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        rho(i, j, k, 0) = 2.0;
      }
    }
  }
  halo.exchange(storage, density);
  operators.assemble_provisional_mass_flux(boundaries, rho_read, u_read,
                                           registry, storage, access,
                                           write_phase, actor, flux_id);
  FaceMassFlux flux = FaceMassFlux::acquire(
      registry, read_storage, access, read_phase, actor, flux_id, topology);
  auto rho_face = storage.acquire_face_write<double>(access, read_phase, actor,
                                                     density_face);
  auto velocity_face_write = storage.acquire_face_write<double>(
      access, read_phase, actor, velocity_face);
  operators.reconstruct_transport_faces(FiniteVolumeQuantity::density(),
                                        boundaries, flux, rho_read, rho_face);
  operators.reconstruct_momentum_faces(boundaries, flux, u_read,
                                       velocity_face_write);
  const auto flux_read = read_storage.acquire_face_read<double>(
      access, read_phase, actor, flux_id);
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    const auto patch = topology.patch_id(face);
    if (!patch.has_value() ||
        topology.face_ownership(face) != hundun::mesh::EntityOwnership::owned) {
      continue;
    }
    const Real3 area =
        geometry.face_area_vector_m2(face, hundun::mesh::FaceSide::owner);
    if (*patch == 0U) {
      HUNDUN_CHECK_NEAR(flux_read(face, 0),
                        1.5 * (area.x + 0.2 * area.y - 0.1 * area.z), 0.0);
      HUNDUN_CHECK_NEAR(rho_face(face, 0), 1.5, 0.0);
      HUNDUN_CHECK_NEAR(velocity_face_write(face, 0), 1.0, 0.0);
      HUNDUN_CHECK_NEAR(velocity_face_write(face, 1), 0.2, 0.0);
      HUNDUN_CHECK_NEAR(velocity_face_write(face, 2), -0.1, 0.0);
    } else if (*patch == 1U) {
      HUNDUN_CHECK(flux_read(face, 0) > 0.0);
      HUNDUN_CHECK_NEAR(rho_face(face, 0), 2.0, 0.0);
      HUNDUN_CHECK_NEAR(velocity_face_write(face, 0), 1.0, 0.0);
      HUNDUN_CHECK_NEAR(velocity_face_write(face, 1), 0.5, 0.0);
      HUNDUN_CHECK_NEAR(velocity_face_write(face, 2), -0.25, 0.0);
    } else if (*patch == 2U) {
      HUNDUN_CHECK(!std::signbit(flux_read(face, 0)) &&
                   flux_read(face, 0) == 0.0);
      HUNDUN_CHECK_NEAR(velocity_face_write(face, 0), 0.0, 0.0);
      HUNDUN_CHECK_NEAR(velocity_face_write(face, 1), 0.0, 0.0);
      HUNDUN_CHECK_NEAR(velocity_face_write(face, 2), 0.0, 0.0);
    } else if (*patch == 3U) {
      HUNDUN_CHECK(!std::signbit(flux_read(face, 0)) &&
                   flux_read(face, 0) == 0.0);
      const double normal_dot = velocity_face_write(face, 0) * area.x +
                                velocity_face_write(face, 1) * area.y +
                                velocity_face_write(face, 2) * area.z;
      HUNDUN_CHECK_NEAR(normal_dot, 0.0, 0.0);
    }
  }

  auto gamma_write =
      storage.acquire_face_write<double>(access, read_phase, actor, gamma);
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    gamma_write(face, 0) = 1.0;
  }
  for (int k = -1; k < decomposition.local_extent().z + 1; ++k) {
    for (int j = -1; j < decomposition.local_extent().y + 1; ++j) {
      for (int i = -1; i < decomposition.local_extent().x + 1; ++i) {
        for (int component = 0; component < 3; ++component) {
          density_gradient_write(i, j, k, component) = 0.0;
        }
      }
    }
  }
  auto scalar_residual_write = storage.view<double>(scalar_residual);
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        scalar_residual_write(i, j, k, 0) = 0.0;
      }
    }
  }
  const auto gamma_read =
      read_storage.acquire_face_read<double>(access, read_phase, actor, gamma);
  operators.accumulate_scalar_diffusive_residual(
      FiniteVolumeQuantity::density(), boundaries, rho_read,
      read_storage.view<double>(density_gradient), gamma_read,
      scalar_residual_write);
  std::vector<double> expected_diffusion(topology.owned_cell_count(), 0.0);
  bool saw_prescribed_diffusion = false;
  for (const auto face : topology.patch(0U).local_faces()) {
    if (topology.face_ownership(face) != hundun::mesh::EntityOwnership::owned) {
      continue;
    }
    const LocalCellId owner = topology.owner(face);
    const Real3 area =
        geometry.face_area_vector_m2(face, hundun::mesh::FaceSide::owner);
    const Real3 center = geometry.cell_center_m(owner);
    const Real3 face_center = geometry.face_center_m(face);
    const Real3 d{2.0 * (face_center.x - center.x),
                  2.0 * (face_center.y - center.y),
                  2.0 * (face_center.z - center.z)};
    const double exterior = boundaries.evaluate_density(0U, 2.0).exterior;
    const double factor =
        (area.x * area.x + area.y * area.y + area.z * area.z) /
        (area.x * d.x + area.y * d.y + area.z * d.z);
    expected_diffusion[owner] -= (exterior - 2.0) * factor;
    saw_prescribed_diffusion = true;
  }
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    HUNDUN_CHECK_NEAR(scalar_residual_write(i, j, k, 0),
                      expected_diffusion[cell], 0.0);
  }
  int local_diffusion_saw = saw_prescribed_diffusion ? 1 : 0;
  int global_diffusion_saw = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_diffusion_saw, &global_diffusion_saw, 1,
                             MPI_INT, MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_diffusion_saw == 1);

  FieldAccessPlan wrong_access(registry);
  wrong_access.declare_access(read_phase, actor, flux_id, AccessMode::write);
  wrong_access.freeze();
  expect_error([&] {
    static_cast<void>(FaceMassFlux::acquire(registry, read_storage,
                                            wrong_access, read_phase, actor,
                                            flux_id, topology));
  });

  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        u(i, j, k, 0) = -1.0;
      }
    }
  }
  halo.exchange(storage, velocity);
  operators.assemble_provisional_mass_flux(boundaries, rho_read, u_read,
                                           registry, storage, access,
                                           write_phase, actor, flux_id);
  bool saw_outlet = false;
  for (const auto face : topology.patch(1U).local_faces()) {
    saw_outlet = true;
    HUNDUN_CHECK(flux_read(face, 0) < 0.0);
  }
  int local_saw = saw_outlet ? 1 : 0;
  int global_saw = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_saw, &global_saw, 1, MPI_INT, MPI_MAX,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_saw == 1);

  std::vector<std::uint64_t> flux_snapshot(topology.local_face_count());
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    flux_snapshot[face] = bits(flux_read(face, 0));
  }
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      for (int i = 0; i < decomposition.local_extent().x; ++i) {
        rho(i, j, k, 0) = -1.0;
      }
    }
  }
  halo.exchange(storage, density);
  expect_error([&] {
    operators.assemble_provisional_mass_flux(boundaries, rho_read, u_read,
                                             registry, storage, access,
                                             write_phase, actor, flux_id);
  });
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    HUNDUN_CHECK(bits(flux_read(face, 0)) == flux_snapshot[face]);
  }

  constexpr std::array<double, 9> affine_gradient{
      0.0, 2.0, 0.0, 0.0, 0.25, 0.0, 0.0, -0.3, 0.0};
  auto gradient_write = storage.view<double>(velocity_gradient);
  auto momentum_write = storage.view<double>(momentum_residual);
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    const Real3 x = geometry.cell_center_m(cell);
    u(i, j, k, 0) = affine_gradient[0] * x.x + affine_gradient[1] * x.y +
                    affine_gradient[2] * x.z;
    u(i, j, k, 1) = affine_gradient[3] * x.x + affine_gradient[4] * x.y +
                    affine_gradient[5] * x.z;
    u(i, j, k, 2) = affine_gradient[6] * x.x + affine_gradient[7] * x.y +
                    affine_gradient[8] * x.z;
    for (int component = 0; component < 9; ++component) {
      gradient_write(i, j, k, component) =
          affine_gradient[static_cast<std::size_t>(component)];
    }
    for (int component = 0; component < 3; ++component) {
      momentum_write(i, j, k, component) = 0.0;
    }
  }
  halo.exchange(storage, velocity);
  auto gradient_halo = HaloExchange::create(
      decomposition,
      ExchangePlan::create(decomposition, decomposition.local_extent(), 1));
  gradient_halo.exchange(storage, velocity_gradient);

  std::vector<std::array<double, 3>> expected(
      topology.owned_cell_count(), std::array<double, 3>{0.0, 0.0, 0.0});
  constexpr double mu = 1.75;
  const auto affine_velocity = [&](LocalCellId cell) {
    const Real3 x = geometry.cell_center_m(cell);
    return Real3{affine_gradient[0] * x.x + affine_gradient[1] * x.y +
                     affine_gradient[2] * x.z,
                 affine_gradient[3] * x.x + affine_gradient[4] * x.y +
                     affine_gradient[5] * x.z,
                 affine_gradient[6] * x.x + affine_gradient[7] * x.y +
                     affine_gradient[8] * x.z};
  };
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    const LocalCellId owner = topology.owner(face);
    const auto neighbour = topology.neighbour(face);
    const auto patch = topology.patch_id(face);
    const Real3 area =
        geometry.face_area_vector_m2(face, hundun::mesh::FaceSide::owner);
    Real3 d = geometry.face_displacement_m(face);
    const Real3 up = affine_velocity(owner);
    Real3 un{};
    if (patch.has_value() && !topology.periodic_pair(face).has_value()) {
      un = boundaries.evaluate_velocity(*patch, up, area).exterior;
      d = Real3{2.0 * d.x, 2.0 * d.y, 2.0 * d.z};
    } else {
      HUNDUN_CHECK(neighbour.has_value());
      un = affine_velocity(*neighbour);
    }
    std::array<double, 9> gface = affine_gradient;
    const double d2 = d.x * d.x + d.y * d.y + d.z * d.z;
    const std::array<double, 3> delta{un.x - up.x, un.y - up.y, un.z - up.z};
    const std::array<double, 3> dv{d.x, d.y, d.z};
    for (int component = 0; component < 3; ++component) {
      const std::size_t base = static_cast<std::size_t>(component) * 3U;
      const double projected =
          gface[base] * d.x + gface[base + 1U] * d.y + gface[base + 2U] * d.z;
      const double correction =
          (delta[static_cast<std::size_t>(component)] - projected) / d2;
      for (int direction = 0; direction < 3; ++direction) {
        gface[base + static_cast<std::size_t>(direction)] +=
            correction * dv[static_cast<std::size_t>(direction)];
      }
    }
    const double divergence = gface[0] + gface[4] + gface[8];
    const std::array<double, 3> sv{area.x, area.y, area.z};
    std::array<double, 3> traction{};
    for (int component = 0; component < 3; ++component) {
      for (int direction = 0; direction < 3; ++direction) {
        double stress = mu * (gface[static_cast<std::size_t>(component) * 3U +
                                    static_cast<std::size_t>(direction)] +
                              gface[static_cast<std::size_t>(direction) * 3U +
                                    static_cast<std::size_t>(component)]);
        if (component == direction) {
          stress -= mu * (2.0 / 3.0) * divergence;
        }
        traction[static_cast<std::size_t>(component)] +=
            stress * sv[static_cast<std::size_t>(direction)];
      }
    }
    if (topology.cell_ownership(owner) ==
        hundun::mesh::EntityOwnership::owned) {
      for (std::size_t component = 0; component < 3U; ++component) {
        expected[owner][component] -= traction[component];
      }
    }
    if (!patch.has_value() && neighbour.has_value() &&
        topology.cell_ownership(*neighbour) ==
            hundun::mesh::EntityOwnership::owned) {
      for (std::size_t component = 0; component < 3U; ++component) {
        expected[*neighbour][component] += traction[component];
      }
    }
  }
  const auto gradient_read = read_storage.view<double>(velocity_gradient);
  operators.accumulate_viscous_residual(boundaries, u_read, gradient_read, mu,
                                        momentum_write);
  std::vector<std::array<std::uint64_t, 3>> scalar_viscosity_bits(
      topology.owned_cell_count());
  bool saw_nonzero_traction = false;
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    for (int component = 0; component < 3; ++component) {
      const double oracle = expected[cell][static_cast<std::size_t>(component)];
      saw_nonzero_traction = saw_nonzero_traction || std::abs(oracle) > 1.0e-14;
      const double tolerance = 512.0 * std::numeric_limits<double>::epsilon() *
                               std::max(1.0, std::abs(oracle));
      HUNDUN_CHECK_NEAR(momentum_write(i, j, k, component), oracle, tolerance);
      scalar_viscosity_bits[cell][static_cast<std::size_t>(component)] =
          bits(momentum_write(i, j, k, component));
    }
  }
  local_saw = saw_nonzero_traction ? 1 : 0;
  global_saw = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_saw, &global_saw, 1, MPI_INT, MPI_MAX,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_saw == 1);

  for (std::size_t face = 0; face < topology.local_face_count(); ++face)
    gamma_write(face, 0) = mu;
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    for (int component = 0; component < 3; ++component)
      momentum_write(i, j, k, component) = 0.0;
  }
  operators.accumulate_viscous_residual(boundaries, u_read, gradient_read,
                                        gamma_read, momentum_write);
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    for (int component = 0; component < 3; ++component)
      HUNDUN_CHECK(bits(momentum_write(i, j, k, component)) ==
                   scalar_viscosity_bits[cell]
                                        [static_cast<std::size_t>(component)]);
  }
  for (std::size_t face = 0; face < topology.local_face_count(); ++face)
    gamma_write(face, 0) = 2.0 * mu;
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    for (int component = 0; component < 3; ++component)
      momentum_write(i, j, k, component) = 0.0;
  }
  operators.accumulate_viscous_residual(boundaries, u_read, gradient_read,
                                        gamma_read, momentum_write);
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    for (int component = 0; component < 3; ++component) {
      const double once = expected[cell][static_cast<std::size_t>(component)];
      HUNDUN_CHECK_NEAR(momentum_write(i, j, k, component), 2.0 * once,
                        1024.0 * std::numeric_limits<double>::epsilon() *
                            std::max(1.0, std::abs(once)));
    }
  }

  operators.reconstruct_momentum_faces(boundaries, flux, u_read,
                                       velocity_face_write);
  std::vector<hundun::finite_volume::PhysicalBoundaryMomentumContribution>
      momentum_contributions;
  std::vector<hundun::finite_volume::PhysicalBoundaryPressureContribution>
      pressure_contributions;
  std::vector<hundun::finite_volume::PhysicalBoundaryTransportContribution>
      transport_contributions;
  momentum_contributions.reserve(topology.local_face_count());
  pressure_contributions.reserve(topology.local_face_count());
  transport_contributions.reserve(topology.local_face_count());
  operators.physical_boundary_momentum_contributions(
      boundaries, flux,
      read_storage.acquire_face_read<double>(access, read_phase, actor,
                                              velocity_face),
      u_read, gradient_read, mu, momentum_contributions);
  operators.physical_boundary_pressure_contributions(
      boundaries, read_storage.view<double>(pressure), pressure_contributions);
  operators.physical_boundary_transport_contributions(
      FiniteVolumeQuantity::density(), boundaries, flux,
      read_storage.acquire_face_read<double>(access, read_phase, actor,
                                              density_face),
      rho_read, read_storage.view<double>(density_gradient), gamma_read,
      transport_contributions);
  HUNDUN_CHECK(momentum_contributions.size() ==
               pressure_contributions.size());
  HUNDUN_CHECK(momentum_contributions.size() ==
               transport_contributions.size());
  HUNDUN_CHECK(!momentum_contributions.empty());
  for (std::size_t index = 0; index < momentum_contributions.size(); ++index) {
    HUNDUN_CHECK(momentum_contributions[index].global_face_id ==
                 pressure_contributions[index].global_face_id);
    HUNDUN_CHECK(momentum_contributions[index].global_face_id ==
                 transport_contributions[index].global_face_id);
    for (std::size_t prior = 0; prior < index; ++prior) {
      HUNDUN_CHECK(momentum_contributions[index].global_face_id !=
                   momentum_contributions[prior].global_face_id);
    }
  }
  const auto momentum_capacity = momentum_contributions.capacity();
  operators.physical_boundary_momentum_contributions(
      boundaries, flux,
      read_storage.acquire_face_read<double>(access, read_phase, actor,
                                              velocity_face),
      u_read, gradient_read, mu, momentum_contributions);
  HUNDUN_CHECK(momentum_contributions.capacity() == momentum_capacity);
}

void run_extent_one_periodic_case(const MpiContext &mpi, int ranks) {
  constexpr Int3 extent{1, 8, 4};
  auto decomposition = StructuredDecomposition::create(
      mpi, extent, std::array<bool, 3>{true, true, true},
      DecompositionOptions{Int3{1, ranks, 1}});
  MeshTopology topology(decomposition);
  MeshGeometry geometry(
      topology, UniformBoxMapping(Real3{0.0, 0.0, 0.0}, Real3{1.0, 1.0, 1.0}));
  BoundaryRegistry boundaries =
      BoundaryRegistry::create(periodic_case(), topology);
  CellCenteredFvmOperators operators =
      CellCenteredFvmOperators::create(topology, geometry);
  FieldRegistry registry;
  const auto density = registry.declare_field(cell_descriptor("rho_one", 1, 2));
  const auto velocity = registry.declare_field(cell_descriptor("u_one", 3, 2));
  const auto flux_id = declare_face_mass_flux(registry);
  const auto residual = registry.declare_field(cell_descriptor("r_one", 1, 0));
  registry.freeze();
  FieldStorage storage(registry, FieldLayoutSet{decomposition.local_extent(),
                                                topology.local_face_count()});
  auto rho = storage.view<double>(density);
  auto u = storage.view<double>(velocity);
  auto r = storage.view<double>(residual);
  for (int k = -2; k < decomposition.local_extent().z + 2; ++k) {
    for (int j = -2; j < decomposition.local_extent().y + 2; ++j) {
      for (int i = -2; i < decomposition.local_extent().x + 2; ++i) {
        rho(i, j, k, 0) = 1.0;
        u(i, j, k, 0) = 1.0;
        u(i, j, k, 1) = 0.0;
        u(i, j, k, 2) = 0.0;
      }
    }
  }
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      r(0, j, k, 0) = 0.0;
    }
  }
  constexpr PhaseId wp = 501U;
  constexpr PhaseId rp = 502U;
  constexpr ActorId actor = 503U;
  FieldAccessPlan access(registry);
  access.declare_access(wp, actor, flux_id, AccessMode::write);
  access.declare_access(rp, actor, flux_id, AccessMode::read);
  access.freeze();
  const FieldStorage &read_storage = storage;
  operators.assemble_provisional_mass_flux(
      boundaries, read_storage.view<double>(density),
      read_storage.view<double>(velocity), registry, storage, access, wp, actor,
      flux_id);
  auto flux = FaceMassFlux::acquire(registry, read_storage, access, rp, actor,
                                    flux_id, topology);
  operators.accumulate_mass_residual(flux, r);
  for (int k = 0; k < decomposition.local_extent().z; ++k) {
    for (int j = 0; j < decomposition.local_extent().y; ++j) {
      HUNDUN_CHECK_NEAR(r(0, j, k, 0), 0.0, 0.0);
    }
  }
}

void run_topology_identity_mutation_contracts(const MpiContext &mpi,
                                              int ranks) {
  constexpr Int3 extent{6, 6, 6};
  auto decomposition = StructuredDecomposition::create(
      mpi, extent, std::array<bool, 3>{true, true, true},
      DecompositionOptions{process_grid_for(ranks)});
  MeshTopology topology(decomposition);
  MeshGeometry geometry(
      topology, UniformBoxMapping(Real3{0.0, 0.0, 0.0}, Real3{1.0, 1.0, 1.0}));
  BoundaryRegistry boundaries =
      BoundaryRegistry::create(periodic_case(), topology);
  CellCenteredFvmOperators operators =
      CellCenteredFvmOperators::create(topology, geometry);

  FieldRegistry registry;
  const auto flux_id = declare_face_mass_flux(registry);
  const auto density =
      registry.declare_field(cell_descriptor("identity_density", 1, 2));
  const auto velocity =
      registry.declare_field(cell_descriptor("identity_velocity", 3, 2));
  const auto transported =
      registry.declare_field(face_descriptor("identity_transported", 1));
  const auto transport_face =
      registry.declare_field(face_descriptor("identity_transport_face", 1));
  const auto momentum_face =
      registry.declare_field(face_descriptor("identity_momentum_face", 3));
  const auto mass_residual =
      registry.declare_field(cell_descriptor("identity_mass_residual", 1, 0));
  const auto convective_residual = registry.declare_field(
      cell_descriptor("identity_convective_residual", 1, 0));
  registry.freeze();

  constexpr PhaseId phase = 701U;
  constexpr ActorId actor = 702U;
  FieldAccessPlan access(registry);
  access.declare_access(phase, actor, flux_id, AccessMode::read_write);
  access.declare_access(phase, actor, transported, AccessMode::read_write);
  access.declare_access(phase, actor, transport_face, AccessMode::read_write);
  access.declare_access(phase, actor, momentum_face, AccessMode::read_write);
  access.freeze();
  FieldStorage storage(registry, FieldLayoutSet{decomposition.local_extent(),
                                                topology.local_face_count()});
  auto flux_write =
      storage.acquire_face_write<double>(access, phase, actor, flux_id);
  auto transported_write =
      storage.acquire_face_write<double>(access, phase, actor, transported);
  auto transport_face_write =
      storage.acquire_face_write<double>(access, phase, actor, transport_face);
  auto momentum_face_write =
      storage.acquire_face_write<double>(access, phase, actor, momentum_face);
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    flux_write(face, 0) = 0.0;
    transported_write(face, 0) = 1.0;
  }
  auto density_write = storage.view<double>(density);
  auto velocity_write = storage.view<double>(velocity);
  for (int k = -2; k < decomposition.local_extent().z + 2; ++k) {
    for (int j = -2; j < decomposition.local_extent().y + 2; ++j) {
      for (int i = -2; i < decomposition.local_extent().x + 2; ++i) {
        density_write(i, j, k, 0) = 1.0;
        velocity_write(i, j, k, 0) = 0.0;
        velocity_write(i, j, k, 1) = 0.0;
        velocity_write(i, j, k, 2) = 0.0;
      }
    }
  }
  auto mass_residual_write = storage.view<double>(mass_residual);
  auto convective_residual_write = storage.view<double>(convective_residual);
  const FieldStorage &read_storage = storage;
  const auto density_read = read_storage.view<double>(density);
  const auto velocity_read = read_storage.view<double>(velocity);
  const auto transported_read =
      read_storage.acquire_face_read<double>(access, phase, actor, transported);

  const auto verify_rejected_before_mutation = [&](const FaceMassFlux &flux) {
    for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
      transport_face_write(face, 0) = 19.0;
      momentum_face_write(face, 0) = 23.0;
      momentum_face_write(face, 1) = 29.0;
      momentum_face_write(face, 2) = 31.0;
    }
    expect_error([&] {
      operators.reconstruct_transport_faces(FiniteVolumeQuantity::density(),
                                            boundaries, flux, density_read,
                                            transport_face_write);
    });
    expect_error([&] {
      operators.reconstruct_momentum_faces(boundaries, flux, velocity_read,
                                           momentum_face_write);
    });
    for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
      HUNDUN_CHECK(bits(transport_face_write(face, 0)) == bits(19.0));
      HUNDUN_CHECK(bits(momentum_face_write(face, 0)) == bits(23.0));
      HUNDUN_CHECK(bits(momentum_face_write(face, 1)) == bits(29.0));
      HUNDUN_CHECK(bits(momentum_face_write(face, 2)) == bits(31.0));
    }

    for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
      const int i = static_cast<int>(
          cell % static_cast<std::size_t>(decomposition.local_extent().x));
      const std::size_t yz =
          cell / static_cast<std::size_t>(decomposition.local_extent().x);
      const int j = static_cast<int>(
          yz % static_cast<std::size_t>(decomposition.local_extent().y));
      const int k = static_cast<int>(
          yz / static_cast<std::size_t>(decomposition.local_extent().y));
      mass_residual_write(i, j, k, 0) = 37.0;
      convective_residual_write(i, j, k, 0) = 41.0;
    }
    expect_error(
        [&] { operators.accumulate_mass_residual(flux, mass_residual_write); });
    expect_error([&] {
      operators.accumulate_convective_residual(flux, transported_read,
                                               convective_residual_write);
    });
    for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
      const int i = static_cast<int>(
          cell % static_cast<std::size_t>(decomposition.local_extent().x));
      const std::size_t yz =
          cell / static_cast<std::size_t>(decomposition.local_extent().x);
      const int j = static_cast<int>(
          yz % static_cast<std::size_t>(decomposition.local_extent().y));
      const int k = static_cast<int>(
          yz / static_cast<std::size_t>(decomposition.local_extent().y));
      HUNDUN_CHECK(bits(mass_residual_write(i, j, k, 0)) == bits(37.0));
      HUNDUN_CHECK(bits(convective_residual_write(i, j, k, 0)) == bits(41.0));
    }
  };

  using Mutation =
      hundun::finite_volume::test::TopologySignatureMutationForTest;
  constexpr std::array<Mutation, 18> mutations{
      Mutation::cell_global_id,
      Mutation::cell_ownership,
      Mutation::face_ownership,
      Mutation::face_owner_local_id,
      Mutation::face_owner_global_id,
      Mutation::face_owner_ownership,
      Mutation::face_neighbour_presence,
      Mutation::face_neighbour_local_id,
      Mutation::face_neighbour_global_id,
      Mutation::face_neighbour_ownership,
      Mutation::logical_face,
      Mutation::face_patch_membership,
      Mutation::periodic_pair,
      Mutation::patch_stable_id,
      Mutation::patch_name,
      Mutation::patch_pairing_kind,
      Mutation::patch_paired_id,
      Mutation::patch_exact_membership};
  for (const Mutation mutation : mutations) {
    hundun::finite_volume::test::mutate_next_topology_signature(mutation);
    auto mutated_flux = FaceMassFlux::acquire(registry, read_storage, access,
                                              phase, actor, flux_id, topology);
    verify_rejected_before_mutation(mutated_flux);
  }

  auto same_shape_decomposition = StructuredDecomposition::create(
      mpi, extent, std::array<bool, 3>{false, true, true},
      DecompositionOptions{process_grid_for(ranks)});
  MeshTopology same_shape_topology(same_shape_decomposition);
  HUNDUN_CHECK(same_shape_topology.global_extent().x ==
               topology.global_extent().x);
  HUNDUN_CHECK(same_shape_topology.global_extent().y ==
               topology.global_extent().y);
  HUNDUN_CHECK(same_shape_topology.global_extent().z ==
               topology.global_extent().z);
  const auto same_shape_box = same_shape_topology.owned_global_box();
  const auto topology_box = topology.owned_global_box();
  HUNDUN_CHECK(same_shape_box.begin.x == topology_box.begin.x);
  HUNDUN_CHECK(same_shape_box.begin.y == topology_box.begin.y);
  HUNDUN_CHECK(same_shape_box.begin.z == topology_box.begin.z);
  HUNDUN_CHECK(same_shape_box.end.x == topology_box.end.x);
  HUNDUN_CHECK(same_shape_box.end.y == topology_box.end.y);
  HUNDUN_CHECK(same_shape_box.end.z == topology_box.end.z);
  HUNDUN_CHECK(same_shape_topology.local_face_count() ==
               topology.local_face_count());
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    HUNDUN_CHECK(same_shape_topology.global_face_id(face) ==
                 topology.global_face_id(face));
  }
  FieldStorage same_shape_storage(
      registry, FieldLayoutSet{same_shape_decomposition.local_extent(),
                               same_shape_topology.local_face_count()});
  auto same_shape_flux_write = same_shape_storage.acquire_face_write<double>(
      access, phase, actor, flux_id);
  for (std::size_t face = 0; face < same_shape_topology.local_face_count();
       ++face) {
    same_shape_flux_write(face, 0) = 0.0;
  }
  const FieldStorage &same_shape_read_storage = same_shape_storage;
  auto same_shape_flux =
      FaceMassFlux::acquire(registry, same_shape_read_storage, access, phase,
                            actor, flux_id, same_shape_topology);
  verify_rejected_before_mutation(same_shape_flux);
}

void run_failure_contracts(const MpiContext &mpi, int ranks) {
  constexpr Int3 extent{6, 6, 6};
  auto decomposition = StructuredDecomposition::create(
      mpi, extent, std::array<bool, 3>{true, true, true},
      DecompositionOptions{process_grid_for(ranks)});
  MeshTopology topology(decomposition);
  MeshGeometry geometry(
      topology, UniformBoxMapping(Real3{0.0, 0.0, 0.0}, Real3{1.0, 1.0, 1.0}));
  CellCenteredFvmOperators operators =
      CellCenteredFvmOperators::create(topology, geometry);

  FieldRegistry wrong_registry;
  auto wrong = hundun::finite_volume::face_mass_flux_descriptor();
  wrong.unit = "kg";
  const auto wrong_id = wrong_registry.declare_field(std::move(wrong));
  wrong_registry.freeze();
  expect_error([&] { require_face_mass_flux_field(wrong_registry, wrong_id); });
  FieldRegistry unfrozen;
  const auto unfrozen_id = declare_face_mass_flux(unfrozen);
  expect_error([&] { require_face_mass_flux_field(unfrozen, unfrozen_id); });
  expect_error([&] { static_cast<void>(declare_face_mass_flux(unfrozen)); });

  FieldRegistry registry;
  const auto flux_id = declare_face_mass_flux(registry);
  const auto transported =
      registry.declare_field(face_descriptor("overflow_transport", 1));
  const auto residual =
      registry.declare_field(cell_descriptor("overflow_residual", 1, 0));
  registry.freeze();
  expect_error([&] { require_face_mass_flux_field(registry, transported); });
  FieldStorage storage(registry, FieldLayoutSet{decomposition.local_extent(),
                                                topology.local_face_count()});
  constexpr PhaseId phase = 601U;
  constexpr ActorId actor = 602U;
  FieldAccessPlan access(registry);
  access.declare_access(phase, actor, flux_id, AccessMode::read_write);
  access.declare_access(phase, actor, transported, AccessMode::read_write);
  access.freeze();
  auto flux_write =
      storage.acquire_face_write<double>(access, phase, actor, flux_id);
  auto transported_write =
      storage.acquire_face_write<double>(access, phase, actor, transported);
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    flux_write(face, 0) = std::numeric_limits<double>::max();
    transported_write(face, 0) = std::numeric_limits<double>::max();
  }
  auto residual_write = storage.view<double>(residual);
  std::vector<std::uint64_t> original(topology.owned_cell_count());
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    residual_write(i, j, k, 0) = 1.0;
    original[cell] = bits(1.0);
  }
  const FieldStorage &read_storage = storage;
  using ConstructionFailure =
      hundun::finite_volume::test::FaceMassFluxConstructionFailureForTest;
  hundun::finite_volume::test::fail_next_face_mass_flux_construction(
      ConstructionFailure::length_error);
  expect_error([&] {
    static_cast<void>(FaceMassFlux::acquire(registry, read_storage, access,
                                            phase, actor, flux_id, topology));
  });
  hundun::finite_volume::test::fail_next_face_mass_flux_construction(
      ConstructionFailure::bad_alloc);
  expect_error([&] {
    static_cast<void>(FaceMassFlux::acquire(registry, read_storage, access,
                                            phase, actor, flux_id, topology));
  });
  auto flux = FaceMassFlux::acquire(registry, read_storage, access, phase,
                                    actor, flux_id, topology);
  FieldStorage wrong_face_layout(
      registry, FieldLayoutSet{decomposition.local_extent(),
                               topology.local_face_count() - 1U});
  const FieldStorage &wrong_face_layout_read = wrong_face_layout;
  expect_error([&] {
    static_cast<void>(FaceMassFlux::acquire(registry, wrong_face_layout_read,
                                            access, phase, actor, flux_id,
                                            topology));
  });
  expect_error(
      [&] { operators.accumulate_mass_residual(flux, residual_write); });
  const auto transported_read =
      read_storage.acquire_face_read<double>(access, phase, actor, transported);
  expect_error([&] {
    operators.accumulate_convective_residual(flux, transported_read,
                                             residual_write);
  });
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    HUNDUN_CHECK(bits(residual_write(i, j, k, 0)) == original[cell]);
  }

  FieldStorage stale_storage(registry,
                             FieldLayoutSet{decomposition.local_extent(),
                                            topology.local_face_count()});
  auto stale_flux_write =
      stale_storage.acquire_face_write<double>(access, phase, actor, flux_id);
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    stale_flux_write(face, 0) = 0.0;
  }
  const FieldStorage &stale_const = stale_storage;
  auto stale_flux = FaceMassFlux::acquire(registry, stale_const, access, phase,
                                          actor, flux_id, topology);
  stale_storage.begin_rebuild();
  expect_error(
      [&] { operators.accumulate_mass_residual(stale_flux, residual_write); });
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    HUNDUN_CHECK(bits(residual_write(i, j, k, 0)) == original[cell]);
  }

  auto other_decomposition = StructuredDecomposition::create(
      mpi, Int3{7, 6, 6}, std::array<bool, 3>{true, true, true},
      DecompositionOptions{process_grid_for(ranks)});
  MeshTopology other_topology(other_decomposition);
  FieldStorage other_storage(registry,
                             FieldLayoutSet{other_decomposition.local_extent(),
                                            other_topology.local_face_count()});
  auto other_flux_write =
      other_storage.acquire_face_write<double>(access, phase, actor, flux_id);
  for (std::size_t face = 0; face < other_topology.local_face_count(); ++face) {
    other_flux_write(face, 0) = 0.0;
  }
  const FieldStorage &other_read_storage = other_storage;
  auto other_flux =
      FaceMassFlux::acquire(registry, other_read_storage, access, phase, actor,
                            flux_id, other_topology);
  expect_error(
      [&] { operators.accumulate_mass_residual(other_flux, residual_write); });

  FaceMassFlux moved_flux(std::move(other_flux));
  HUNDUN_CHECK(moved_flux.field_id() == flux_id);
  HUNDUN_CHECK(moved_flux.face_count() == other_topology.local_face_count());

#ifdef HUNDUN_FINITE_VOLUME_ENABLE_TEST_ACCESS
  FieldStorage prepared_storage(
      registry, FieldLayoutSet{decomposition.local_extent(),
                               topology.local_face_count()});
  {
    auto prepared_flux_write = prepared_storage.acquire_face_write<double>(
        access, phase, actor, flux_id);
    for (std::size_t face = 0; face < topology.local_face_count(); ++face)
      prepared_flux_write(face, 0) = 0.0;
  }
  const FieldStorage &prepared_read_storage = prepared_storage;
  auto prepared =
      hundun::finite_volume::test::PreparedFaceMassFluxForTest::create(
          topology);
  const auto bind_prepared = [&] {
    return prepared.bind(registry, prepared_read_storage, access, phase, actor,
                         flux_id, topology);
  };
  {
    auto first = bind_prepared();
    expect_error([&] { static_cast<void>(bind_prepared()); });
    FaceMassFlux moved(std::move(first));
    expect_error(
        [&] { operators.accumulate_mass_residual(first, residual_write); });
    expect_error([&] { static_cast<void>(bind_prepared()); });
    HUNDUN_CHECK(moved.field_id() == flux_id);
    HUNDUN_CHECK(moved.face_count() == topology.local_face_count());
  }
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(decomposition.local_extent().x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(decomposition.local_extent().x);
    const int j = static_cast<int>(
        yz % static_cast<std::size_t>(decomposition.local_extent().y));
    const int k = static_cast<int>(
        yz / static_cast<std::size_t>(decomposition.local_extent().y));
    residual_write(i, j, k, 0) = 0.0;
  }
  {
    allocation_probe::AllocationAttemptGuard allocation_guard;
    auto second = bind_prepared();
    operators.accumulate_mass_residual(second, residual_write);
    HUNDUN_CHECK(allocation_guard.attempts() == 0U);
  }
  {
    auto third = bind_prepared();
    HUNDUN_CHECK(third.face_count() == topology.local_face_count());
  }
  expect_error([&] {
    static_cast<void>(prepared.bind(
        registry, wrong_face_layout_read, access, phase, actor, flux_id,
        topology));
  });
  expect_error([&] {
    static_cast<void>(prepared.bind(registry, prepared_read_storage, access,
                                    phase, actor, transported, topology));
  });
  expect_error([&] {
    static_cast<void>(prepared.bind(registry, other_read_storage, access,
                                    phase, actor, flux_id, other_topology));
  });
  FieldStorage prepared_stale_storage(
      registry, FieldLayoutSet{decomposition.local_extent(),
                               topology.local_face_count()});
  {
    auto stale_write = prepared_stale_storage.acquire_face_write<double>(
        access, phase, actor, flux_id);
    for (std::size_t face = 0; face < topology.local_face_count(); ++face)
      stale_write(face, 0) = 0.0;
  }
  const FieldStorage &prepared_stale_read = prepared_stale_storage;
  {
    auto stale_prepared = prepared.bind(
        registry, prepared_stale_read, access, phase, actor, flux_id,
        topology);
    prepared_stale_storage.begin_rebuild();
    expect_error([&] {
      operators.accumulate_mass_residual(stale_prepared, residual_write);
    });
  }
  auto moved_prepared = std::move(prepared);
  expect_error([&] { static_cast<void>(bind_prepared()); });
  {
    auto final = moved_prepared.bind(registry, prepared_read_storage, access,
                                     phase, actor, flux_id, topology);
    HUNDUN_CHECK(final.field_id() == flux_id);
  }
#endif
}

} // namespace

int main(int argc, char **argv) {
  MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    MpiContext mpi = MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
    const std::string_view mode =
        argc == 1 ? "full" : (argv[1] == nullptr ? "" : argv[1]);
    HUNDUN_CHECK(
        argc == 1 ||
        (argc == 2 &&
         (mode == "task25-warped-free-stream" ||
          mode == "task25-warped-diffusion-mutation" ||
          mode == "task11-no-slip-traction")));
    if (mode == "task25-warped-free-stream") {
      run_gradient_case(mpi, mpi.size(), true);
      return;
    }
    if (mode == "task25-warped-diffusion-mutation") {
      const auto normal16 = warped_diffusion_error(mpi, mpi.size(), 16);
      static_cast<void>(
          reject_reversed_nonorthogonal_diffusion(mpi, mpi.size(), normal16));
      return;
    }
    if (mode == "task11-no-slip-traction") {
      run_no_slip_quadratic_traction_order(mpi, mpi.size());
      return;
    }
    check_limiter();
    run_gradient_case(mpi, mpi.size(), false);
    run_gradient_case(mpi, mpi.size(), true);
    run_gradient_failure_contracts(mpi, mpi.size());
    run_failure_contracts(mpi, mpi.size());
    run_topology_identity_mutation_contracts(mpi, mpi.size());
    run_shared_flux_case(mpi, mpi.size());
    const double reconstruction_error16 =
        smooth_transport_reconstruction_error(mpi, mpi.size(), 16);
    const double reconstruction_error32 =
        smooth_transport_reconstruction_error(mpi, mpi.size(), 32);
    HUNDUN_CHECK(std::log(reconstruction_error16 / reconstruction_error32) /
                     std::log(2.0) >=
                 1.8);
    run_uniform_quadratic_diffusion_case(mpi, mpi.size());
    run_no_slip_quadratic_traction_order(mpi, mpi.size());
    run_physical_boundary_case(mpi, mpi.size(), false);
    run_physical_boundary_case(mpi, mpi.size(), true);
    run_extent_one_periodic_case(mpi, mpi.size());
    const auto error16 = warped_diffusion_error(mpi, mpi.size(), 16);
    const auto error32 = warped_diffusion_error(mpi, mpi.size(), 32);
    const auto error64 = warped_diffusion_error(mpi, mpi.size(), 64);
    const double order16_32 =
        std::log(error16.error / error32.error) / std::log(2.0);
    const double order32_64 =
        std::log(error32.error / error64.error) / std::log(2.0);
    HUNDUN_CHECK(order16_32 >= 1.8);
    HUNDUN_CHECK(order32_64 >= 1.8);
    std::array<double, 3> decomposition_differences{};
    if (mpi.size() > 1) {
      auto self = MpiContext::duplicate(MPI_COMM_SELF);
      const std::array<WarpedDiffusionResult, 3> references{
          warped_diffusion_error(self, 1, 16),
          warped_diffusion_error(self, 1, 32),
          warped_diffusion_error(self, 1, 64)};
      const std::array<const WarpedDiffusionResult *, 3> distributed{
          &error16, &error32, &error64};
      for (std::size_t level = 0; level < distributed.size(); ++level) {
        double local_max = 0.0;
        double local_reference_max = 0.0;
        for (std::size_t cell = 0;
             cell < distributed[level]->global_ids.size(); ++cell) {
          const auto global = static_cast<std::size_t>(
              distributed[level]->global_ids[cell]);
          HUNDUN_CHECK(global <
                       references[level].physical_laplacian.size());
          local_max = std::max(
              local_max,
              std::abs(distributed[level]->physical_laplacian[cell] -
                       references[level].physical_laplacian[global]));
          local_reference_max =
              std::max(local_reference_max,
                       std::abs(
                           references[level].physical_laplacian[global]));
        }
        double comparison[2]{local_max, local_reference_max};
        HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, comparison, 2, MPI_DOUBLE,
                                   MPI_MAX, mpi.comm()) == MPI_SUCCESS);
        decomposition_differences[level] = comparison[0];
        HUNDUN_CHECK(comparison[0] <=
                     5.0e-12 * std::max(1.0, comparison[1]));
      }
    }
    const double reversed_order =
        reject_reversed_nonorthogonal_diffusion(mpi, mpi.size(), error16);
    if (mpi.rank() == 0) {
      std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
                << "TASK25_WARPED_DIFFUSION ranks=" << mpi.size()
                << " error16=" << error16.error
                << " error32=" << error32.error
                << " error64=" << error64.error
                << " order16_32=" << order16_32
                << " order32_64=" << order32_64
                << " integral_defect16=" << error16.integral_defect
                << " integral_defect32=" << error32.integral_defect
                << " integral_defect64=" << error64.integral_defect
                << " decomposition16=" << decomposition_differences[0]
                << " decomposition32=" << decomposition_differences[1]
                << " decomposition64=" << decomposition_differences[2]
                << " reversed_nonorth_order=" << reversed_order << '\n';
    }
  });
}

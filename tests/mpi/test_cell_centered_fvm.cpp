// SPDX-License-Identifier: Apache-2.0

#include "finite_volume/src/cell_centered_fvm_test_seam.hpp"
#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/config/resolved_case.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/exchange_plan.hpp"
#include "hundun/runtime/field_access_plan.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/field_storage.hpp"
#include "hundun/runtime/halo_exchange.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
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
#include <string>
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
        u(i, j, k, 0) = 0.75;
        u(i, j, k, 1) = -0.25;
        u(i, j, k, 2) = 0.125;
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
  operators.reconstruct_momentum_faces(boundaries, flux, velocity_read, u_face);
  const auto flux_read = read_storage.acquire_face_read<double>(
      access, read_phase, actor, mass_flux);
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    HUNDUN_CHECK(std::isfinite(flux_read(face, 0)));
    HUNDUN_CHECK(std::isfinite(rho_face(face, 0)));
    HUNDUN_CHECK(rho_face(face, 0) > 0.0);
    HUNDUN_CHECK_NEAR(u_face(face, 0), 0.75, 0.0);
    HUNDUN_CHECK_NEAR(u_face(face, 1), -0.25, 0.0);
    HUNDUN_CHECK_NEAR(u_face(face, 2), 0.125, 0.0);
    const auto pair = topology.periodic_pair(face);
    if (pair.has_value()) {
      const auto local_pair = topology.find_local_face(*pair);
      if (local_pair.has_value()) {
        HUNDUN_CHECK(bits(flux_read(face, 0)) ==
                     bits(-flux_read(*local_pair, 0)));
        HUNDUN_CHECK(bits(rho_face(face, 0)) == bits(rho_face(*local_pair, 0)));
      }
    }
  }

  struct FaceReplicaRecord final {
    std::uint64_t id;
    std::uint64_t flux;
    std::uint64_t transported;
  };
  std::vector<FaceReplicaRecord> local_records;
  local_records.reserve(topology.local_face_count());
  for (std::size_t face = 0; face < topology.local_face_count(); ++face) {
    local_records.push_back(FaceReplicaRecord{topology.global_face_id(face),
                                              bits(flux_read(face, 0)),
                                              bits(rho_face(face, 0))});
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
          u(i, j, k, 0) = 1.0 + static_cast<double>(global.x * global.x);
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
    const double q_p = 1.0 + static_cast<double>(p.x * p.x);
    const double q_n = 1.0 + static_cast<double>((p.x + 1) * (p.x + 1));
    const double central = 0.5 * q_p + 0.5 * q_n;
    HUNDUN_CHECK_NEAR(u_face(*threshold_face, 0), central, 0.0);
    hundun::finite_volume::test::override_next_face_metrics(face_id, 0.250001,
                                                            0.0);
    operators.reconstruct_momentum_faces(boundaries, flux, velocity_read,
                                         u_face);
    const double q_pm1 = 1.0 + static_cast<double>((p.x - 1) * (p.x - 1));
    const double slope = monotonized_central(q_p - q_pm1, q_n - q_p);
    const double limited = q_p + 0.5 * slope;
    HUNDUN_CHECK_NEAR(u_face(*threshold_face, 0), limited, 0.0);
    HUNDUN_CHECK(limited != central);
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

double warped_diffusion_error(const MpiContext &mpi, int ranks, int n) {
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
    gradient_write(i, j, k, 0) = two_pi * std::cos(two_pi * center.x);
    gradient_write(i, j, k, 1) = 0.0;
    gradient_write(i, j, k, 2) = 0.0;
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
  operators.accumulate_scalar_diffusive_residual(
      FiniteVolumeQuantity::density(), boundaries, q_read, gradient_read,
      gamma_read, residual_write);
  double local_error = 0.0;
  double local_volume = 0.0;
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
    const double exact = two_pi * two_pi * std::sin(two_pi * center.x);
    const double error = residual_write(i, j, k, 0) / volume - exact;
    local_error += error * error * volume;
    local_volume += volume;
  }
  double values[2]{local_error, local_volume};
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, values, 2, MPI_DOUBLE, MPI_SUM,
                             mpi.comm()) == MPI_SUCCESS);
  return std::sqrt(values[0] / values[1]);
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

  constexpr std::array<double, 9> affine_gradient{1.0,  2.0, -0.5, -1.0, 0.25,
                                                  0.75, 0.4, -0.3, 2.0};
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
    u(i, j, k, 0) = 0.5 + affine_gradient[0] * x.x + affine_gradient[1] * x.y +
                    affine_gradient[2] * x.z;
    u(i, j, k, 1) = -0.75 + affine_gradient[3] * x.x +
                    affine_gradient[4] * x.y + affine_gradient[5] * x.z;
    u(i, j, k, 2) = 1.25 + affine_gradient[6] * x.x + affine_gradient[7] * x.y +
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
    return Real3{0.5 + affine_gradient[0] * x.x + affine_gradient[1] * x.y +
                     affine_gradient[2] * x.z,
                 -0.75 + affine_gradient[3] * x.x + affine_gradient[4] * x.y +
                     affine_gradient[5] * x.z,
                 1.25 + affine_gradient[6] * x.x + affine_gradient[7] * x.y +
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
    }
  }
  local_saw = saw_nonzero_traction ? 1 : 0;
  global_saw = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_saw, &global_saw, 1, MPI_INT, MPI_MAX,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_saw == 1);
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
}

} // namespace

int main(int argc, char **argv) {
  MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    MpiContext mpi = MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
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
    run_physical_boundary_case(mpi, mpi.size(), false);
    run_physical_boundary_case(mpi, mpi.size(), true);
    run_extent_one_periodic_case(mpi, mpi.size());
    const double error8 = warped_diffusion_error(mpi, mpi.size(), 8);
    const double error16 = warped_diffusion_error(mpi, mpi.size(), 16);
    const double error32 = warped_diffusion_error(mpi, mpi.size(), 32);
    HUNDUN_CHECK(std::log(error8 / error16) / std::log(2.0) >= 1.8);
    HUNDUN_CHECK(std::log(error16 / error32) / std::log(2.0) >= 1.8);
  });
}

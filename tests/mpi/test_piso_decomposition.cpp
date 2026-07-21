// SPDX-License-Identifier: Apache-2.0

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/config/resolved_case.hpp"
#include "hundun/execution/execution.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/flow/constant_density_piso.hpp"
#include "hundun/linear/conjugate_gradient.hpp"
#include "hundun/linear/preconditioners.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/exchange_plan.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/halo_exchange.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "tests/support/test_main.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kDt = 1.0e-4;
constexpr hundun::runtime::Int3 kExtent{8, 8, 4};

hundun::config::FlowCaseConfig periodic_case() {
  hundun::config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type =
      hundun::config::SimulationType::variable_density_flow;
  config.density_model = hundun::config::DensityModel::constant;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  constexpr std::array<hundun::config::PatchName, 6> names{
      hundun::config::PatchName::x_min, hundun::config::PatchName::x_max,
      hundun::config::PatchName::y_min, hundun::config::PatchName::y_max,
      hundun::config::PatchName::z_min, hundun::config::PatchName::z_max};
  for (std::size_t index = 0; index < names.size(); ++index) {
    config.boundaries[index].patch = names[index];
    config.boundaries[index].type = hundun::config::BoundaryType::periodic;
  }
  return config;
}

hundun::runtime::FieldDescriptor cell_field(const char *name,
                                            std::uint32_t components) {
  return {name,
          "1",
          "task18_decomposition",
          hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64,
          components,
          2,
          false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::runtime::FieldDescriptor face_field(const char *name,
                                            std::uint32_t components) {
  return {name,
          "1",
          "task18_decomposition",
          hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64,
          components,
          0,
          false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

std::array<double, 3> smooth_velocity(hundun::runtime::Real3 point) {
  return {std::sin(point.x) * std::cos(point.y),
          -std::cos(point.x) * std::sin(point.y), 0.0};
}

double smooth_pressure(hundun::runtime::Real3 point) {
  return 0.25 * (std::cos(2.0 * point.x) + std::cos(2.0 * point.y));
}

struct DecompositionResult final {
  std::vector<hundun::mesh::GlobalCellId> cell_ids;
  std::vector<double> velocity;
  std::vector<double> pressure;
  std::vector<hundun::mesh::GlobalFaceId> face_ids;
  std::vector<double> face_mass_flux;
};

DecompositionResult run_one_step(const hundun::runtime::MpiContext &mpi,
                                 hundun::runtime::Int3 process_grid) {
  auto decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, kExtent, {true, true, true},
      hundun::runtime::DecompositionOptions{process_grid});
  const hundun::mesh::MeshTopology topology(decomposition);
  const hundun::mesh::MeshGeometry geometry(
      topology, hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0},
                                                {2.0 * kPi, 2.0 * kPi, 1.0}));
  auto boundaries =
      hundun::boundary::BoundaryRegistry::create(periodic_case(), topology);

  hundun::runtime::FieldRegistry registry;
  hundun::flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell_field("rho", 1U));
  fields.velocity = registry.declare_field(cell_field("velocity", 3U));
  fields.mechanical_pressure = registry.declare_field(cell_field("pi", 1U));
  fields.face_velocity = registry.declare_field(face_field("u_face", 3U));
  fields.face_mass_flux =
      hundun::finite_volume::declare_face_mass_flux(registry);
  registry.freeze();

  const auto box = decomposition.owned_box();
  const hundun::runtime::Int3 local{box.end.x - box.begin.x,
                                    box.end.y - box.begin.y,
                                    box.end.z - box.begin.z};
  auto state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, kDt, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  hundun::flow::FlowLayerValues initial;
  initial.density.assign(topology.owned_cell_count(), 1.0);
  initial.velocity.resize(topology.owned_cell_count() * 3U);
  initial.mechanical_pressure.resize(topology.owned_cell_count());
  initial.face_velocity.resize(topology.local_face_count() * 3U);
  initial.face_mass_flux.resize(topology.local_face_count());
  for (hundun::mesh::LocalCellId cell = 0; cell < topology.owned_cell_count();
       ++cell) {
    const auto velocity = smooth_velocity(geometry.cell_center_m(cell));
    for (std::size_t component = 0; component < 3U; ++component) {
      initial.velocity[cell * 3U + component] = velocity[component];
    }
    initial.mechanical_pressure[cell] =
        smooth_pressure(geometry.cell_center_m(cell));
  }
  for (hundun::mesh::LocalFaceId face = 0; face < topology.local_face_count();
       ++face) {
    const auto velocity = smooth_velocity(geometry.face_center_m(face));
    for (std::size_t component = 0; component < 3U; ++component) {
      initial.face_velocity[face * 3U + component] = velocity[component];
    }
    const auto area =
        geometry.face_area_vector_m2(face, hundun::mesh::FaceSide::owner);
    initial.face_mass_flux[face] =
        velocity[0] * area.x + velocity[1] * area.y + velocity[2] * area.z;
  }
  state.seed_accepted_layers(initial, initial);

  auto halo = hundun::runtime::HaloExchange::create(
      decomposition,
      hundun::runtime::ExchangePlan::create(decomposition, local, 2));
  hundun::execution::CpuReferenceContext execution;
  hundun::linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  hundun::linear::ConjugateGradientSolver pressure_solver(execution, mpi);
  hundun::linear::JacobiPreconditioner mx(execution);
  hundun::linear::JacobiPreconditioner my(execution);
  hundun::linear::JacobiPreconditioner mz(execution);
  hundun::linear::JacobiPreconditioner pressure_preconditioner(execution);
  auto flow = hundun::flow::FixedStepConstantDensityFlow::create(
      decomposition, topology, geometry, boundaries, mpi, execution, halo,
      momentum_solver, {&mx, &my, &mz}, pressure_solver,
      pressure_preconditioner);
  const auto stencil = hundun::flow::make_momentum_time_stencil(
      hundun::flow::MomentumTimeOrder::backward_euler, kDt, 0.0);
  const auto report = flow.attempt(state, 1.0, 0.1, stencil, {}, {});
  HUNDUN_CHECK(report.disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(report.pressure_corrector_count == 2U);
  HUNDUN_CHECK(report.final_continuity_normalized_l2 <= 1.0e-10);

  const auto values = state.snapshot(hundun::flow::FlowLayer::committed);
  DecompositionResult result;
  result.velocity = values.velocity;
  result.pressure = values.mechanical_pressure;
  result.face_mass_flux = values.face_mass_flux;
  result.cell_ids.resize(topology.owned_cell_count());
  for (hundun::mesh::LocalCellId cell = 0; cell < topology.owned_cell_count();
       ++cell) {
    result.cell_ids[cell] = topology.global_cell_id(cell);
  }
  result.face_ids.resize(topology.local_face_count());
  for (hundun::mesh::LocalFaceId face = 0; face < topology.local_face_count();
       ++face) {
    result.face_ids[face] = topology.global_face_id(face);
  }
  return result;
}

DecompositionResult replicate_result(MPI_Comm communicator,
                                     const DecompositionResult &local) {
  int communicator_size = 0;
  HUNDUN_CHECK(MPI_Comm_size(communicator, &communicator_size) == MPI_SUCCESS);
  HUNDUN_CHECK(communicator_size > 0);
  HUNDUN_CHECK(local.cell_ids.size() <=
               static_cast<std::size_t>(std::numeric_limits<int>::max()));
  HUNDUN_CHECK(local.face_ids.size() <=
               static_cast<std::size_t>(std::numeric_limits<int>::max()));
  const int local_cell_count = static_cast<int>(local.cell_ids.size());
  const int local_face_count = static_cast<int>(local.face_ids.size());
  std::vector<int> cell_counts(static_cast<std::size_t>(communicator_size));
  std::vector<int> face_counts(static_cast<std::size_t>(communicator_size));
  HUNDUN_CHECK(MPI_Allgather(&local_cell_count, 1, MPI_INT, cell_counts.data(),
                             1, MPI_INT, communicator) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allgather(&local_face_count, 1, MPI_INT, face_counts.data(),
                             1, MPI_INT, communicator) == MPI_SUCCESS);
  std::vector<int> cell_displacements(cell_counts.size());
  std::vector<int> face_displacements(face_counts.size());
  int gathered_cell_count = 0;
  int gathered_face_count = 0;
  for (std::size_t rank = 0; rank < cell_counts.size(); ++rank) {
    cell_displacements[rank] = gathered_cell_count;
    face_displacements[rank] = gathered_face_count;
    gathered_cell_count += cell_counts[rank];
    gathered_face_count += face_counts[rank];
  }

  std::vector<hundun::mesh::GlobalCellId> gathered_cell_ids(
      static_cast<std::size_t>(gathered_cell_count));
  std::vector<hundun::mesh::GlobalFaceId> gathered_face_ids(
      static_cast<std::size_t>(gathered_face_count));
  HUNDUN_CHECK(MPI_Allgatherv(local.cell_ids.data(), local_cell_count,
                              MPI_UINT64_T, gathered_cell_ids.data(),
                              cell_counts.data(), cell_displacements.data(),
                              MPI_UINT64_T, communicator) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allgatherv(local.face_ids.data(), local_face_count,
                              MPI_UINT64_T, gathered_face_ids.data(),
                              face_counts.data(), face_displacements.data(),
                              MPI_UINT64_T, communicator) == MPI_SUCCESS);

  std::vector<double> local_cell_values(local.cell_ids.size() * 4U);
  for (std::size_t cell = 0; cell < local.cell_ids.size(); ++cell) {
    for (std::size_t component = 0; component < 3U; ++component) {
      local_cell_values[cell * 4U + component] =
          local.velocity[cell * 3U + component];
    }
    local_cell_values[cell * 4U + 3U] = local.pressure[cell];
  }
  std::vector<int> cell_value_counts(cell_counts.size());
  std::vector<int> cell_value_displacements(cell_counts.size());
  int gathered_cell_value_count = 0;
  for (std::size_t rank = 0; rank < cell_counts.size(); ++rank) {
    cell_value_displacements[rank] = gathered_cell_value_count;
    cell_value_counts[rank] = 4 * cell_counts[rank];
    gathered_cell_value_count += cell_value_counts[rank];
  }
  std::vector<double> gathered_cell_values(
      static_cast<std::size_t>(gathered_cell_value_count));
  HUNDUN_CHECK(MPI_Allgatherv(local_cell_values.data(), 4 * local_cell_count,
                              MPI_DOUBLE, gathered_cell_values.data(),
                              cell_value_counts.data(),
                              cell_value_displacements.data(), MPI_DOUBLE,
                              communicator) == MPI_SUCCESS);
  std::vector<double> gathered_face_values(
      static_cast<std::size_t>(gathered_face_count));
  HUNDUN_CHECK(MPI_Allgatherv(local.face_mass_flux.data(), local_face_count,
                              MPI_DOUBLE, gathered_face_values.data(),
                              face_counts.data(), face_displacements.data(),
                              MPI_DOUBLE, communicator) == MPI_SUCCESS);

  std::map<hundun::mesh::GlobalCellId, std::array<double, 4>> cells;
  for (std::size_t cell = 0; cell < gathered_cell_ids.size(); ++cell) {
    const std::array<double, 4> values{gathered_cell_values[cell * 4U],
                                       gathered_cell_values[cell * 4U + 1U],
                                       gathered_cell_values[cell * 4U + 2U],
                                       gathered_cell_values[cell * 4U + 3U]};
    HUNDUN_CHECK(cells.emplace(gathered_cell_ids[cell], values).second);
  }
  std::map<hundun::mesh::GlobalFaceId, double> faces;
  for (std::size_t face = 0; face < gathered_face_ids.size(); ++face) {
    const auto inserted =
        faces.emplace(gathered_face_ids[face], gathered_face_values[face]);
    if (!inserted.second) {
      HUNDUN_CHECK_NEAR(inserted.first->second, gathered_face_values[face],
                        5.0e-12);
    }
  }

  DecompositionResult result;
  result.cell_ids.reserve(cells.size());
  result.velocity.reserve(cells.size() * 3U);
  result.pressure.reserve(cells.size());
  for (const auto &[id, values] : cells) {
    result.cell_ids.push_back(id);
    result.velocity.insert(result.velocity.end(), values.begin(),
                           values.begin() + 3);
    result.pressure.push_back(values[3]);
  }
  result.face_ids.reserve(faces.size());
  result.face_mass_flux.reserve(faces.size());
  for (const auto &[id, value] : faces) {
    result.face_ids.push_back(id);
    result.face_mass_flux.push_back(value);
  }
  return result;
}

std::pair<double, double>
maximum_field_difference(const hundun::runtime::MpiContext &world,
                         const DecompositionResult &partitioned,
                         const DecompositionResult &reference) {
  std::unordered_map<hundun::mesh::GlobalCellId, std::size_t> cells;
  for (std::size_t index = 0; index < reference.cell_ids.size(); ++index) {
    HUNDUN_CHECK(cells.emplace(reference.cell_ids[index], index).second);
  }
  std::unordered_map<hundun::mesh::GlobalFaceId, std::size_t> faces;
  for (std::size_t index = 0; index < reference.face_ids.size(); ++index) {
    HUNDUN_CHECK(faces.emplace(reference.face_ids[index], index).second);
  }

  double comparison[2]{};
  for (std::size_t index = 0; index < partitioned.cell_ids.size(); ++index) {
    const auto found = cells.find(partitioned.cell_ids[index]);
    HUNDUN_CHECK(found != cells.end());
    const std::size_t reference_index = found->second;
    comparison[0] =
        std::max(comparison[0], std::abs(partitioned.pressure[index] -
                                         reference.pressure[reference_index]));
    comparison[1] =
        std::max({comparison[1], std::abs(partitioned.pressure[index]),
                  std::abs(reference.pressure[reference_index])});
    for (std::size_t component = 0; component < 3U; ++component) {
      const double left = partitioned.velocity[index * 3U + component];
      const double right = reference.velocity[reference_index * 3U + component];
      comparison[0] = std::max(comparison[0], std::abs(left - right));
      comparison[1] =
          std::max({comparison[1], std::abs(left), std::abs(right)});
    }
  }
  for (std::size_t index = 0; index < partitioned.face_ids.size(); ++index) {
    const auto found = faces.find(partitioned.face_ids[index]);
    HUNDUN_CHECK(found != faces.end());
    const double left = partitioned.face_mass_flux[index];
    const double right = reference.face_mass_flux[found->second];
    comparison[0] = std::max(comparison[0], std::abs(left - right));
    comparison[1] = std::max({comparison[1], std::abs(left), std::abs(right)});
  }
  world.allreduce_fp64_in_place(
      comparison, 2U, hundun::runtime::Fp64ReductionOperation::maximum);
  return {comparison[0], comparison[1]};
}

void run_decomposition_comparison(const hundun::runtime::MpiContext &world) {
  HUNDUN_CHECK(world.size() == 4);
  if (world.rank() == 0)
    std::cerr << "TASK18_DECOMP_STAGE world4\n" << std::flush;
  const auto four_rank = run_one_step(world, {2, 2, 1});

  if (world.rank() == 0)
    std::cerr << "TASK18_DECOMP_STAGE sub2\n" << std::flush;
  MPI_Comm raw_two_rank = MPI_COMM_NULL;
  HUNDUN_CHECK(MPI_Comm_split(MPI_COMM_WORLD, world.rank() / 2,
                              world.rank() % 2, &raw_two_rank) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Comm_set_errhandler(raw_two_rank, MPI_ERRORS_RETURN) ==
               MPI_SUCCESS);
  auto two_rank = hundun::runtime::MpiContext::duplicate(raw_two_rank);
  const auto two_rank_local = run_one_step(two_rank, {2, 1, 1});
  const auto two_rank_result = replicate_result(raw_two_rank, two_rank_local);
  HUNDUN_CHECK(MPI_Comm_free(&raw_two_rank) == MPI_SUCCESS);

  if (world.rank() == 0)
    std::cerr << "TASK18_DECOMP_STAGE self1\n" << std::flush;
  auto self = hundun::runtime::MpiContext::duplicate(MPI_COMM_SELF);
  const auto one_rank = run_one_step(self, {1, 1, 1});

  if (world.rank() == 0)
    std::cerr << "TASK18_DECOMP_STAGE compare\n" << std::flush;
  const auto four_to_two =
      maximum_field_difference(world, four_rank, two_rank_result);
  const auto two_to_one =
      maximum_field_difference(world, two_rank_result, one_rank);
  const auto four_to_one = maximum_field_difference(world, four_rank, one_rank);
  const auto accepted = [](std::pair<double, double> comparison) {
    return comparison.first <= 5.0e-12 * std::max(1.0, comparison.second);
  };
  HUNDUN_CHECK(accepted(four_to_two));
  HUNDUN_CHECK(accepted(two_to_one));
  HUNDUN_CHECK(accepted(four_to_one));
  if (world.rank() == 0) {
    std::cout << "TASK18_DECOMPOSITION max_4_to_2=" << four_to_two.first
              << " max_2_to_1=" << two_to_one.first
              << " max_4_to_1=" << four_to_one.first << '\n';
  }
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  auto world = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] { run_decomposition_comparison(world); });
}

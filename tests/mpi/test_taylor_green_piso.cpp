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

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kRho = 1.0;
constexpr double kNu = 0.1;

hundun::runtime::Int3 process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw hundun::runtime::Error("unsupported Taylor-Green rank count");
}

hundun::config::FlowCaseConfig periodic_case() {
  hundun::config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type =
      hundun::config::SimulationType::variable_density_flow;
  config.density_model = hundun::config::DensityModel::constant;
  config.physics.rho_ref_kg_per_m3 = kRho;
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
          "task18_taylor_green",
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
          "task18_taylor_green",
          hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64,
          components,
          0,
          false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

std::array<double, 3> exact_velocity(hundun::runtime::Real3 point,
                                     double time_s) {
  const double decay = std::exp(-2.0 * kNu * time_s);
  return {std::sin(point.x) * std::cos(point.y) * decay,
          -std::cos(point.x) * std::sin(point.y) * decay, 0.0};
}

double exact_pressure(hundun::runtime::Real3 point, double time_s) {
  return 0.25 * kRho * (std::cos(2.0 * point.x) + std::cos(2.0 * point.y)) *
         std::exp(-4.0 * kNu * time_s);
}

struct TrajectoryResult final {
  std::vector<double> velocity;
  double cell_volume_m3{};
  double velocity_l2_error{};
  double relative_mass_defect{};
  double max_continuity{};
  std::array<double, 3> max_momentum_conservation_defect{};
  std::uint64_t committed_steps{};
  hundun::flow::MomentumTimeOrder final_order{
      hundun::flow::MomentumTimeOrder::backward_euler};
  double final_dt_s{};
  double final_previous_dt_s{};
};

hundun::flow::FlowLayerValues
exact_layer(const hundun::mesh::MeshTopology &topology,
            const hundun::mesh::MeshGeometry &geometry, double time_s) {
  hundun::flow::FlowLayerValues values;
  const std::size_t cells = topology.owned_cell_count();
  values.density.assign(cells, kRho);
  values.velocity.resize(cells * 3U);
  values.mechanical_pressure.resize(cells);
  values.face_velocity.resize(topology.local_face_count() * 3U);
  values.face_mass_flux.resize(topology.local_face_count());
  for (hundun::mesh::LocalCellId cell = 0; cell < cells; ++cell) {
    const auto point = geometry.cell_center_m(cell);
    const auto velocity = exact_velocity(point, time_s);
    for (std::size_t component = 0; component < 3U; ++component) {
      values.velocity[cell * 3U + component] = velocity[component];
    }
    values.mechanical_pressure[cell] = exact_pressure(point, time_s);
  }
  for (hundun::mesh::LocalFaceId face = 0; face < topology.local_face_count();
       ++face) {
    const auto velocity = exact_velocity(geometry.face_center_m(face), time_s);
    for (std::size_t component = 0; component < 3U; ++component) {
      values.face_velocity[face * 3U + component] = velocity[component];
    }
    const auto area =
        geometry.face_area_vector_m2(face, hundun::mesh::FaceSide::owner);
    values.face_mass_flux[face] =
        kRho *
        (velocity[0] * area.x + velocity[1] * area.y + velocity[2] * area.z);
  }
  return values;
}

TrajectoryResult run_trajectory(const hundun::runtime::MpiContext &mpi,
                                int cells_xy, double dt_s,
                                double final_time_s) {
  const hundun::runtime::Int3 extent{cells_xy, cells_xy, 4};
  auto decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, extent, {true, true, true},
      hundun::runtime::DecompositionOptions{process_grid(mpi.size())});
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
      {0U, 0.0, dt_s, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  const auto initial = exact_layer(topology, geometry, 0.0);
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

  const int steps = static_cast<int>(std::llround(final_time_s / dt_s));
  HUNDUN_CHECK(steps > 0);
  HUNDUN_CHECK_NEAR(static_cast<double>(steps) * dt_s, final_time_s, 1.0e-15);
  double max_continuity = 0.0;
  std::array<double, 3> max_momentum_conservation_defect{};
  for (int step = 0; step < steps; ++step) {
    const auto order = step == 0
                           ? hundun::flow::MomentumTimeOrder::backward_euler
                           : hundun::flow::MomentumTimeOrder::bdf2;
    const auto stencil = hundun::flow::make_momentum_time_stencil(
        order, dt_s, step == 0 ? 0.0 : dt_s);
    const auto report = flow.attempt(state, kRho, kRho * kNu, stencil, {}, {});
    if (report.disposition !=
            hundun::flow::StepAttemptDisposition::committed &&
        mpi.rank() == 0) {
      std::cerr << "TASK18_TGV_GATE_FAILURE step=" << step
                << " reason=" << static_cast<int>(report.reason)
                << " momentum_residual="
                << report.final_momentum_normalized_l2[0] << ','
                << report.final_momentum_normalized_l2[1] << ','
                << report.final_momentum_normalized_l2[2]
                << " mass_defect="
                << report.final_mass_relative_conservation_defect
                << " momentum_defect="
                << report.final_momentum_relative_conservation_defect[0]
                << ','
                << report.final_momentum_relative_conservation_defect[1]
                << ','
                << report.final_momentum_relative_conservation_defect[2]
                << '\n';
    }
    HUNDUN_CHECK(report.disposition ==
                 hundun::flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(report.pressure_corrector_count == 2U);
    HUNDUN_CHECK(report.final_continuity_normalized_l2 <= 1.0e-10);
    max_continuity =
        std::max(max_continuity, report.final_continuity_normalized_l2);
    for (std::size_t component = 0; component < 3U; ++component) {
      HUNDUN_CHECK(
          report.final_momentum_relative_conservation_defect[component] <=
          5.0e-11);
      max_momentum_conservation_defect[component] = std::max(
          max_momentum_conservation_defect[component],
          report.final_momentum_relative_conservation_defect[component]);
    }
  }

  const auto result = state.snapshot(hundun::flow::FlowLayer::committed);
  double sums[4]{};
  for (hundun::mesh::LocalCellId cell = 0; cell < topology.owned_cell_count();
       ++cell) {
    const double volume = geometry.cell_volume_m3(cell);
    const auto exact =
        exact_velocity(geometry.cell_center_m(cell), final_time_s);
    for (std::size_t component = 0; component < 3U; ++component) {
      const double difference =
          result.velocity[cell * 3U + component] - exact[component];
      sums[0] += volume * difference * difference;
      sums[1] += volume * exact[component] * exact[component];
    }
    sums[2] += volume * result.density[cell];
    sums[3] += volume;
  }
  mpi.allreduce_fp64_in_place(sums, 4U,
                              hundun::runtime::Fp64ReductionOperation::sum);
  const auto metadata = state.metadata();
  return {result.velocity,
          geometry.cell_volume_m3(0U),
          std::sqrt(sums[0] / sums[1]),
          std::abs(sums[2] - kRho * sums[3]) / (kRho * sums[3]),
          max_continuity,
          max_momentum_conservation_defect,
          metadata.step,
          metadata.order,
          metadata.dt_s,
          metadata.previous_dt_s};
}

double difference_norm(const hundun::runtime::MpiContext &mpi,
                       const TrajectoryResult &left,
                       const TrajectoryResult &right) {
  HUNDUN_CHECK(left.velocity.size() == right.velocity.size());
  HUNDUN_CHECK_NEAR(left.cell_volume_m3, right.cell_volume_m3, 0.0);
  double sum = 0.0;
  for (std::size_t index = 0; index < left.velocity.size(); ++index) {
    const double difference = left.velocity[index] - right.velocity[index];
    sum += left.cell_volume_m3 * difference * difference;
  }
  mpi.allreduce_fp64_in_place(&sum, 1U,
                              hundun::runtime::Fp64ReductionOperation::sum);
  return std::sqrt(sum);
}

void run_taylor_green(const hundun::runtime::MpiContext &mpi) {
  constexpr std::array<int, 3> spatial_cells{16, 32, 64};
  constexpr double spatial_dt = 1.0e-4;
  std::array<TrajectoryResult, 3> spatial{
      run_trajectory(mpi, spatial_cells[0], spatial_dt, spatial_dt),
      run_trajectory(mpi, spatial_cells[1], spatial_dt, spatial_dt),
      run_trajectory(mpi, spatial_cells[2], spatial_dt, spatial_dt)};
  const double spatial_order_0 =
      std::log(spatial[0].velocity_l2_error / spatial[1].velocity_l2_error) /
      std::log(2.0);
  const double spatial_order_1 =
      std::log(spatial[1].velocity_l2_error / spatial[2].velocity_l2_error) /
      std::log(2.0);

  constexpr int temporal_cells = 64;
  constexpr std::array<double, 3> temporal_dt{0.02, 0.01, 0.005};
  constexpr double final_time = 0.04;
  std::array<TrajectoryResult, 3> temporal{
      run_trajectory(mpi, temporal_cells, temporal_dt[0], final_time),
      run_trajectory(mpi, temporal_cells, temporal_dt[1], final_time),
      run_trajectory(mpi, temporal_cells, temporal_dt[2], final_time)};
  const double coarse_medium = difference_norm(mpi, temporal[0], temporal[1]);
  const double medium_fine = difference_norm(mpi, temporal[1], temporal[2]);
  const double temporal_order =
      std::log(coarse_medium / medium_fine) / std::log(2.0);

  if (mpi.rank() == 0) {
    for (std::size_t level = 0; level < spatial.size(); ++level) {
      std::cout << "TASK18_TAYLOR_GREEN_SPATIAL cells=" << spatial_cells[level]
                << " dt=" << spatial_dt
                << " velocity_l2=" << spatial[level].velocity_l2_error
                << " mass_defect=" << spatial[level].relative_mass_defect
                << " continuity=" << spatial[level].max_continuity
                << " momentum_defect="
                << spatial[level].max_momentum_conservation_defect[0] << ','
                << spatial[level].max_momentum_conservation_defect[1] << ','
                << spatial[level].max_momentum_conservation_defect[2] << '\n';
    }
    std::cout << "TASK18_TAYLOR_GREEN_SPATIAL_ORDER coarse=" << spatial_order_0
              << " fine=" << spatial_order_1 << '\n';
    for (std::size_t level = 0; level < temporal.size(); ++level) {
      std::cout << "TASK18_TAYLOR_GREEN_TEMPORAL cells=" << temporal_cells
                << " dt=" << temporal_dt[level]
                << " velocity_l2=" << temporal[level].velocity_l2_error
                << " mass_defect=" << temporal[level].relative_mass_defect
                << " continuity=" << temporal[level].max_continuity
                << " momentum_defect="
                << temporal[level].max_momentum_conservation_defect[0] << ','
                << temporal[level].max_momentum_conservation_defect[1] << ','
                << temporal[level].max_momentum_conservation_defect[2]
                << " commits=" << temporal[level].committed_steps
                << " order=" << static_cast<int>(temporal[level].final_order)
                << " final_dt=" << temporal[level].final_dt_s
                << " previous_dt=" << temporal[level].final_previous_dt_s
                << '\n';
    }
    std::cout << "TASK18_TAYLOR_GREEN_TEMPORAL_ORDER coarse_medium="
              << coarse_medium << " medium_fine=" << medium_fine
              << " order=" << temporal_order << '\n';
  }

  HUNDUN_CHECK(spatial_order_0 >= 1.8);
  HUNDUN_CHECK(spatial_order_1 >= 1.8);
  HUNDUN_CHECK(temporal_order >= 1.8);
  for (const auto &item : spatial) {
    HUNDUN_CHECK(item.relative_mass_defect <= 5.0e-12);
    HUNDUN_CHECK(item.committed_steps == 1U);
    HUNDUN_CHECK(item.final_order ==
                 hundun::flow::MomentumTimeOrder::backward_euler);
  }
  for (std::size_t level = 0; level < temporal.size(); ++level) {
    const auto &item = temporal[level];
    HUNDUN_CHECK(item.relative_mass_defect <= 5.0e-12);
    HUNDUN_CHECK(item.committed_steps == (2U << level));
    HUNDUN_CHECK(item.final_order == hundun::flow::MomentumTimeOrder::bdf2);
    HUNDUN_CHECK(item.final_dt_s == temporal_dt[level]);
    HUNDUN_CHECK(item.final_previous_dt_s == temporal_dt[level]);
  }
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] { run_taylor_green(mpi); });
}

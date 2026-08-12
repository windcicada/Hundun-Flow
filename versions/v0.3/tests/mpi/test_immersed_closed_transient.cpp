// SPDX-License-Identifier: Apache-2.0

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/flow_immersed.hpp"
#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/ib_domain.hpp"
#include "hundun/ib_surface.hpp"
#include "hundun/ib_local_flow_pattern.hpp"
#include "hundun/ib_surface_query.hpp"
#include "hundun/lin_bicgstab.hpp"
#include "hundun/lin_conjugate_gradient.hpp"
#include "hundun/lin_preconditioners.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "src/ib_quadratic_reconstruction_detail.hpp"
#include "tests/support/stage3_mms.hpp"
#include "tests/support/stage3_test_contracts.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <thread>
#include <vector>

namespace {

using namespace hundun;
constexpr std::size_t kParallelWorkerBudget = 96U;

runtime::Int3 process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw runtime::Error("unsupported closed-transient rank count");
}

config::FlowCaseConfig closed_case(int cells, int ranks) {
  config::FlowCaseConfig result{};
  result.schema_version = 2;
  result.case_name = "stage3-task11-closed-transient";
  result.simulation_type = config::SimulationType::variable_density_flow;
  result.density_model = config::DensityModel::constant;
  result.resources.expected_ranks = ranks;
  result.resources.process_grid = process_grid(ranks);
  result.mesh.cells = {cells, cells, cells};
  result.mesh.origin_m = {0.0, 0.0, 0.0};
  result.mesh.length_m = {1.0, 1.0, 1.0};
  result.physics.rho_ref_kg_per_m3 = 1.0;
  result.physics.dynamic_viscosity_pa_s = 0.01;
  result.physics.inlet_consistency_rtol = 1.0e-12;
  constexpr std::array<config::PatchName, 6> names{
      config::PatchName::x_min, config::PatchName::x_max,
      config::PatchName::y_min, config::PatchName::y_max,
      config::PatchName::z_min, config::PatchName::z_max};
  for (std::size_t patch = 0U; patch < names.size(); ++patch) {
    result.boundaries[patch].patch = names[patch];
    result.boundaries[patch].type = config::BoundaryType::no_slip_wall;
  }
  return result;
}

runtime::FieldDescriptor cell_descriptor(const char *name,
                                         std::uint32_t components) {
  return {name,
          "1",
          "stage3_task11_closed_transient",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          4,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor face_descriptor(const char *name,
                                         std::uint32_t components) {
  return {name,
          "1",
          "stage3_task11_closed_transient",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::never};
}

void check_mpi(int code, const char *operation) {
  if (code != MPI_SUCCESS)
    throw runtime::Error(std::string(operation) + " failed");
}

void sum_in_place(const runtime::MpiContext &mpi, double *values,
                  std::size_t count) {
  HUNDUN_CHECK(count <=
               static_cast<std::size_t>(std::numeric_limits<int>::max()));
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, values, static_cast<int>(count),
                          MPI_DOUBLE, MPI_SUM, mpi.comm()),
            "MPI_Allreduce(closed-transient sum)");
}

void max_in_place(const runtime::MpiContext &mpi, double *values,
                  std::size_t count) {
  HUNDUN_CHECK(count <=
               static_cast<std::size_t>(std::numeric_limits<int>::max()));
  check_mpi(MPI_Allreduce(MPI_IN_PLACE, values, static_cast<int>(count),
                          MPI_DOUBLE, MPI_MAX, mpi.comm()),
            "MPI_Allreduce(closed-transient maximum)");
}

std::string collective_surface_path(
    const runtime::MpiContext &mpi,
    std::optional<test::Stage3TemporaryDirectory> &root_directory,
    const test::stage3::ManufacturedSurface &surface) {
  std::string path_text;
  if (mpi.rank() == 0) {
    root_directory.emplace("task11-closed-transient");
    const auto path = root_directory->path() / "sphere.stl";
    test::write_text(path, test::ascii_stl(surface.triangles, "sphere"));
    path_text = path.string();
  }
  std::uint64_t size = path_text.size();
  check_mpi(MPI_Bcast(&size, 1, MPI_UINT64_T, 0, mpi.comm()),
            "MPI_Bcast(closed-transient path size)");
  HUNDUN_CHECK(size <=
               static_cast<std::uint64_t>(std::numeric_limits<int>::max()));
  path_text.resize(static_cast<std::size_t>(size));
  check_mpi(MPI_Bcast(path_text.data(), static_cast<int>(path_text.size()),
                      MPI_BYTE, 0, mpi.comm()),
            "MPI_Bcast(closed-transient path)");
  return path_text;
}

bool finite(runtime::Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

double dot(runtime::Real3 left, runtime::Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

struct IntegralState final {
  double mass_kg{};
  double kinetic_energy_J{};
};

IntegralState integrals(const runtime::MpiContext &mpi,
                        const mesh::MeshTopology &topology,
                        const mesh::MeshGeometry &geometry,
                        const immersed::ImmersedDomain &domain,
                        const flow::FlowLayerValues &values) {
  std::array<double, 2> result{};
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    const double rho = values.density[cell];
    const double u = values.velocity[cell * 3U];
    const double v = values.velocity[cell * 3U + 1U];
    const double w = values.velocity[cell * 3U + 2U];
    const double volume = geometry.cell_volume_m3(cell);
    HUNDUN_CHECK(rho > 0.0 && std::isfinite(rho));
    HUNDUN_CHECK(std::isfinite(u) && std::isfinite(v) && std::isfinite(w));
    result[0] += rho * volume;
    result[1] += 0.5 * rho * (u * u + v * v + w * w) * volume;
  }
  sum_in_place(mpi, result.data(), result.size());
  return {result[0], result[1]};
}

double penetration_linf(const runtime::MpiContext &mpi,
                        const immersed::WallQuadraturePlan &wall_plan,
                        const runtime::FieldView<const double> &velocity) {
  double result = 0.0;
  for (const auto &point : wall_plan.local_points()) {
    runtime::Real3 wall_velocity{};
    wall_velocity.x = immersed::detail::value_with_origin_constraint(
        point.reconstruction, point.position_m, velocity, 0U, 0.0);
    wall_velocity.y = immersed::detail::value_with_origin_constraint(
        point.reconstruction, point.position_m, velocity, 1U, 0.0);
    wall_velocity.z = immersed::detail::value_with_origin_constraint(
        point.reconstruction, point.position_m, velocity, 2U, 0.0);
    result = std::max(
        result, std::abs(dot(wall_velocity, point.solid_to_fluid_normal)));
  }
  max_in_place(mpi, &result, 1U);
  return result;
}

void check_force(const flow::ForceAttemptReport &force) {
  for (const auto value :
       {force.operator_force.pressure_N, force.operator_force.viscous_N,
        force.operator_force.total_N, force.budget_reaction.pressure_N,
        force.budget_reaction.viscous_N, force.budget_reaction.total_N,
        force.surface_traction.pressure_N,
        force.surface_traction.viscous_N, force.surface_traction.total_N,
        force.consistency.pressure_N, force.consistency.viscous_N,
        force.consistency.total_N})
    HUNDUN_CHECK(finite(value));
}

void run(int argc, char **argv, const runtime::MpiContext &mpi) {
  HUNDUN_CHECK(argc == 2);
  const std::string selector = argv[1];
  HUNDUN_CHECK(selector == "smoke" || selector == "fast" ||
               selector == "acceptance");
  const int cells = selector == "smoke" ? 12 : selector == "fast" ? 48 : 96;
  const double final_time_s = selector == "smoke"  ? 0.0025
                              : selector == "fast" ? 0.02
                                                   : 0.05;
  const auto body = test::stage3::approved_body(test::stage3::BodyKind::sphere);
  const auto surface_fixture = test::stage3::make_manufactured_surface(
      body, 1.0 / static_cast<double>(cells));
  std::optional<test::Stage3TemporaryDirectory> root_directory;
  const auto surface_path =
      collective_surface_path(mpi, root_directory, surface_fixture);

  const auto decomposition = runtime::StructuredDecomposition::create(
      mpi, {cells, cells, cells}, {false, false, false},
      runtime::DecompositionOptions{process_grid(mpi.size())});
  const mesh::MeshTopology topology(decomposition);
  const mesh::MeshGeometry geometry(
      topology, mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  const auto boundaries = boundary::BoundaryRegistry::create(
      closed_case(cells, mpi.size()), topology);
  const auto surface = immersed::ImmersedSurface::load_collective(
      std::filesystem::path(surface_path), 1.0, mpi, 0);
  const auto query = immersed::SurfaceQuery::create(surface);
  const auto domain = immersed::ImmersedDomain::create(
      surface, query, config::ImmersedFluidSide::outside, topology, geometry,
      boundaries, mpi);
  const auto ghost_plan = immersed::GhostStencilPlan::create(
      surface, query, domain, topology, geometry, decomposition, mpi);
  const auto wall_plan = immersed::WallQuadraturePlan::create(
      surface, query, domain, topology, geometry, mpi);
  HUNDUN_CHECK(ghost_plan.maximum_halo_reach() <= 4U);

  double local_h_min = std::numeric_limits<double>::infinity();
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count(); ++cell)
    if (domain.region(cell) == immersed::CellRegion::fluid)
      local_h_min =
          std::min(local_h_min, std::cbrt(geometry.cell_volume_m3(cell)));
  double inverse_h_min = std::isfinite(local_h_min) ? 1.0 / local_h_min : 0.0;
  max_in_place(mpi, &inverse_h_min, 1U);
  HUNDUN_CHECK(inverse_h_min > 0.0 && std::isfinite(inverse_h_min));
  const double h_min = 1.0 / inverse_h_min;
  constexpr double rho = 1.0;
  constexpr double mu = 0.01;
  constexpr double speed_scale = 1.0;
  const double stable_dt =
      std::min(0.25 * h_min / speed_scale, 0.25 * rho * h_min * h_min / mu);
  const auto step_count =
      static_cast<std::uint64_t>(std::ceil(final_time_s / stable_dt));
  HUNDUN_CHECK(step_count > 0U);
  const double dt = final_time_s / static_cast<double>(step_count);
  HUNDUN_CHECK(speed_scale * dt / h_min <=
               0.25 * (1.0 + 64.0 * std::numeric_limits<double>::epsilon()));
  HUNDUN_CHECK(mu * dt / (rho * h_min * h_min) <=
               0.25 * (1.0 + 64.0 * std::numeric_limits<double>::epsilon()));

  runtime::FieldRegistry registry;
  flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell_descriptor("rho", 1U));
  fields.velocity = registry.declare_field(cell_descriptor("velocity", 3U));
  fields.mechanical_pressure =
      registry.declare_field(cell_descriptor("pi", 1U));
  fields.face_velocity =
      registry.declare_field(face_descriptor("face_velocity", 3U));
  fields.face_mass_flux = finite_volume::declare_face_mass_flux(registry);
  registry.freeze();

  flow::FlowLayerValues initial;
  initial.density.resize(topology.owned_cell_count(), 0.0);
  initial.velocity.resize(topology.owned_cell_count() * 3U, 0.0);
  initial.mechanical_pressure.resize(topology.owned_cell_count(), 0.0);
  initial.face_velocity.resize(topology.local_face_count() * 3U, 0.0);
  initial.face_mass_flux.resize(topology.local_face_count(), 0.0);
  std::vector<test::stage3::MmsCellAverage> initial_averages(
      topology.owned_cell_count());
  geometry.require_compatible(topology);
  const auto hardware_threads =
      std::max(1U, std::thread::hardware_concurrency());
  const std::size_t worker_count = std::min<std::size_t>(
      std::min<std::size_t>(
          hardware_threads,
          std::max<std::size_t>(1U, kParallelWorkerBudget /
                                        static_cast<std::size_t>(mpi.size()))),
      topology.owned_cell_count());
  std::vector<std::exception_ptr> failures(worker_count);
  std::vector<std::thread> workers;
  workers.reserve(worker_count);
  std::atomic<bool> cancel{false};
  try {
    for (std::size_t worker = 0U; worker < worker_count; ++worker) {
      workers.emplace_back([&, worker] {
        try {
          const auto begin =
              topology.owned_cell_count() * worker / worker_count;
          const auto end =
              topology.owned_cell_count() * (worker + 1U) / worker_count;
          for (mesh::LocalCellId cell = begin;
               cell < end && !cancel.load(std::memory_order_relaxed); ++cell) {
            if (domain.region(cell) == immersed::CellRegion::fluid)
              initial_averages[cell] =
                  test::stage3::detail::evaluate_cell_average_validated(
                      topology, geometry, cell, body, 0.0);
          }
        } catch (...) {
          failures[worker] = std::current_exception();
          cancel.store(true, std::memory_order_relaxed);
        }
      });
    }
  } catch (...) {
    cancel.store(true, std::memory_order_relaxed);
    for (auto &worker : workers)
      worker.join();
    throw;
  }
  for (auto &worker : workers)
    worker.join();
  for (const auto &failure : failures)
    if (failure)
      std::rethrow_exception(failure);
  for (mesh::LocalCellId cell = 0U; cell < topology.owned_cell_count();
       ++cell) {
    if (domain.region(cell) != immersed::CellRegion::fluid)
      continue;
    const auto &exact = initial_averages[cell];
    initial.density[cell] = rho;
    initial.velocity[cell * 3U] = exact.velocity_m_per_s.x;
    initial.velocity[cell * 3U + 1U] = exact.velocity_m_per_s.y;
    initial.velocity[cell * 3U + 2U] = exact.velocity_m_per_s.z;
  }
  for (mesh::LocalFaceId face = 0U; face < topology.local_face_count();
       ++face) {
    const bool owner_active =
        domain.region(topology.owner(face)) == immersed::CellRegion::fluid;
    const auto neighbour = topology.neighbour(face);
    const bool neighbour_active =
        neighbour.has_value() &&
        domain.region(*neighbour) == immersed::CellRegion::fluid;
    if (!owner_active || (neighbour.has_value() && !neighbour_active))
      continue;
    const auto exact =
        test::stage3::evaluate_mms(body, geometry.face_center_m(face), 0.0);
    initial.face_velocity[face * 3U] = exact.velocity_m_per_s.x;
    initial.face_velocity[face * 3U + 1U] = exact.velocity_m_per_s.y;
    initial.face_velocity[face * 3U + 2U] = exact.velocity_m_per_s.z;
    const auto area = geometry.face_area_vector_m2(face, mesh::FaceSide::owner);
    initial.face_mass_flux[face] = rho * dot(exact.velocity_m_per_s, area);
  }

  auto state = flow::FlowState::create(
      registry, {decomposition.local_extent(), topology.local_face_count()},
      fields, {0U, 0.0, dt, 0.0, flow::MomentumTimeOrder::backward_euler});
  state.seed_accepted_layers(initial, initial);
  auto halo = runtime::HaloExchange::create(
      decomposition, runtime::ExchangePlan::create(
                         decomposition, decomposition.local_extent(), 4));
  execution::CpuReferenceContext execution;
  linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  linear::BiCGStabSolver pressure_solver(execution, mpi);
  linear::JacobiPreconditioner mx(execution);
  linear::JacobiPreconditioner my(execution);
  linear::JacobiPreconditioner mz(execution);
  linear::JacobiPreconditioner pressure_pc(execution);
  immersed::LocalFlowPatternTransform transform;
  auto immersed_flow = flow::FixedStepImmersedFlow::create(
      decomposition, topology, geometry, boundaries, &domain, &ghost_plan,
      &wall_plan, &transform, nullptr, mpi, execution, halo, momentum_solver,
      {&mx, &my, &mz}, pressure_solver, pressure_pc);
  const flow::ImmersedFlowPhysics physics{config::DensityModel::constant,
                                    rho,
                                    mu,
                                    std::nullopt,
                                    std::nullopt,
                                    std::nullopt};
  linear::SolveControl control;
  control.max_iterations = 5000U;
  control.residual_recompute_interval = 20U;

  const auto initial_integrals =
      integrals(mpi, topology, geometry, domain, initial);
  HUNDUN_CHECK(initial_integrals.mass_kg > 0.0);
  HUNDUN_CHECK(initial_integrals.kinetic_energy_J > 0.0);
  double previous_energy = initial_integrals.kinetic_energy_J;
  for (std::uint64_t step = 0U; step < step_count; ++step) {
    const auto order = step == 0U ? flow::MomentumTimeOrder::backward_euler
                                  : flow::MomentumTimeOrder::bdf2;
    const auto stencil =
        flow::make_momentum_time_stencil(order, dt, step == 0U ? 0.0 : dt);
    const auto attempt =
        immersed_flow.attempt(state, physics, stencil, control, control);
    const auto &base = std::get<flow::StepAttemptReport>(attempt.base);
    if (mpi.rank() == 0 &&
        base.disposition != flow::StepAttemptDisposition::committed) {
      std::cerr << "closed_transient rejected step=" << step + 1U
                << " disposition=" << static_cast<int>(base.disposition)
                << " reason=" << static_cast<int>(base.reason)
                << " lowest_failing_rank=" << base.lowest_failing_rank
                << " correctors=" << base.pressure_corrector_count
                << " continuity=" << base.final_continuity_normalized_l2
                << " pressure_residual=" << base.final_pressure_residual_l2
                << " momentum_residuals="
                << base.final_momentum_normalized_l2[0] << ','
                << base.final_momentum_normalized_l2[1] << ','
                << base.final_momentum_normalized_l2[2] << '\n';
      for (std::size_t component = 0U;
           component < base.momentum.components.size(); ++component) {
        const auto &solve = base.momentum.components[component];
        std::cerr << "  momentum_solve[" << component
                  << "] reason=" << static_cast<int>(solve.reason)
                  << " iterations=" << solve.iterations
                  << " initial=" << solve.initial_residual
                  << " recursive=" << solve.recursive_residual
                  << " final=" << solve.final_residual
                  << " lowest_failing_rank=" << solve.lowest_failing_rank
                  << '\n';
      }
      for (std::size_t corrector = 0U; corrector < base.pressure.size();
           ++corrector) {
        const auto &solve = base.pressure[corrector];
        std::cerr << "  pressure_solve[" << corrector
                  << "] reason=" << static_cast<int>(solve.reason)
                  << " iterations=" << solve.iterations
                  << " initial=" << solve.initial_residual
                  << " recursive=" << solve.recursive_residual
                  << " final=" << solve.final_residual
                  << " lowest_failing_rank=" << solve.lowest_failing_rank
                  << '\n';
      }
    }
    HUNDUN_CHECK(base.disposition == flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(base.reason == flow::StepFailureReason::none);
    HUNDUN_CHECK(base.pressure_corrector_count == 2U);
    HUNDUN_CHECK(base.final_continuity_normalized_l2 <= 1.0e-10);
    HUNDUN_CHECK(base.final_pressure_residual_l2 <= 1.0e-10);
    for (const double residual : base.final_momentum_normalized_l2)
      HUNDUN_CHECK(residual <= 1.0e-9);
    HUNDUN_CHECK(attempt.force.has_value());
    check_force(*attempt.force);

    const auto committed = state.snapshot(flow::FlowLayer::committed);
    const auto current = integrals(mpi, topology, geometry, domain, committed);
    HUNDUN_CHECK(std::isfinite(current.kinetic_energy_J));
    const double energy_tolerance = 64.0 *
                                    std::numeric_limits<double>::epsilon() *
                                    std::max(1.0, previous_energy);
    HUNDUN_CHECK(current.kinetic_energy_J <=
                 previous_energy + energy_tolerance);
    const double relative_mass_error =
        std::abs(current.mass_kg - initial_integrals.mass_kg) /
        initial_integrals.mass_kg;
    HUNDUN_CHECK(relative_mass_error <= 5.0e-12);
    if (mpi.rank() == 0)
      std::cerr << "closed_transient step=" << step + 1U << '/' << step_count
                << " time=" << state.metadata().time_s
                << " kinetic_energy=" << current.kinetic_energy_J
                << " relative_mass_error=" << relative_mass_error
                << " continuity=" << base.final_continuity_normalized_l2
                << " pressure_residual=" << base.final_pressure_residual_l2
                << '\n';
    previous_energy = current.kinetic_energy_J;
  }

  const auto velocity =
      state.layer(flow::FlowLayer::committed).view<double>(fields.velocity);
  const double penetration = penetration_linf(mpi, wall_plan, velocity);
  HUNDUN_CHECK(penetration <= 8192.0 * std::numeric_limits<double>::epsilon());
  const auto metadata = state.metadata();
  HUNDUN_CHECK(metadata.step == step_count);
  HUNDUN_CHECK(std::abs(metadata.time_s - final_time_s) <=
               64.0 * std::numeric_limits<double>::epsilon() *
                   std::max(1.0, final_time_s));
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  const auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] { run(argc, argv, mpi); });
}

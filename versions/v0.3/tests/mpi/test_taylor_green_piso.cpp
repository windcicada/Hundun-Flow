// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "tests/support/flow_constant_density_piso_test_access.hpp"
#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/flow_constant_density_piso.hpp"
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
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
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
  hundun::flow::FlowLayerValues committed;
  std::vector<hundun::mesh::GlobalCellId> cell_ids;
  std::vector<hundun::mesh::GlobalFaceId> face_ids;
  std::vector<hundun::mesh::LocalFaceId> face_local_ids;
  double pressure_mean{};
  double cell_volume_m3{};
  double velocity_l2_error{};
  double relative_mass_defect{};
  double max_continuity{};
  double max_pressure_residual{};
  std::array<double, 3> max_momentum_residual{};
  std::array<double, 3> max_momentum_conservation_defect{};
  std::uint64_t max_solve_iterations{};
  std::uint64_t committed_steps{};
  hundun::flow::MomentumTimeOrder final_order{
      hundun::flow::MomentumTimeOrder::backward_euler};
  double final_dt_s{};
  double final_previous_dt_s{};
};

struct FlowStateSnapshot final {
  hundun::flow::FlowLayerValues history;
  hundun::flow::FlowLayerValues committed;
  hundun::flow::FlowLayerValues trial;
  hundun::flow::AcceptedStepMetadata metadata;
};

FlowStateSnapshot capture_state(const hundun::flow::FlowState &state) {
  return {state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata()};
}

bool state_bitwise_equal(const FlowStateSnapshot &left,
                         const FlowStateSnapshot &right) noexcept {
  return hundun::test::flow_layer_values_bitwise_equal(left.history,
                                                        right.history) &&
         hundun::test::flow_layer_values_bitwise_equal(left.committed,
                                                        right.committed) &&
         hundun::test::flow_layer_values_bitwise_equal(left.trial,
                                                        right.trial) &&
         hundun::test::accepted_step_metadata_bitwise_equal(left.metadata,
                                                             right.metadata);
}

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
                                double final_time_s,
                                bool analytic_warped = false) {
  const hundun::runtime::Int3 extent{cells_xy, cells_xy, 4};
  auto decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, extent, {true, true, true},
      hundun::runtime::DecompositionOptions{process_grid(mpi.size())});
  const auto fixture_local = decomposition.local_extent();
  HUNDUN_CHECK(fixture_local.x >= 2);
  HUNDUN_CHECK(fixture_local.y >= 2);
  HUNDUN_CHECK(fixture_local.z >= 2);
  const hundun::mesh::MeshTopology topology(decomposition);
  const hundun::mesh::MeshGeometry geometry =
      analytic_warped
          ? hundun::mesh::MeshGeometry(
                topology,
                hundun::mesh::AnalyticWarpedBoxMapping(
                    {0.0, 0.0, 0.0}, {2.0 * kPi, 2.0 * kPi, 1.0},
                    {0.02, -0.015, 0.01}))
          : hundun::mesh::MeshGeometry(
                topology,
                hundun::mesh::UniformBoxMapping(
                    {0.0, 0.0, 0.0}, {2.0 * kPi, 2.0 * kPi, 1.0}));
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
  hundun::linear::ConjugateGradientSolver momentum_cg(execution, mpi);
  hundun::linear::ConjugateGradientSolver pressure_cg(execution, mpi);
  hundun::linear::BiCGStabSolver momentum_bicgstab(execution, mpi);
  hundun::linear::BiCGStabSolver pressure_bicgstab(execution, mpi);
  const hundun::linear::LinearSolver &momentum_solver =
      analytic_warped
          ? static_cast<const hundun::linear::LinearSolver &>(
                momentum_bicgstab)
          : static_cast<const hundun::linear::LinearSolver &>(momentum_cg);
  const hundun::linear::LinearSolver &pressure_solver =
      analytic_warped
          ? static_cast<const hundun::linear::LinearSolver &>(
                pressure_bicgstab)
          : static_cast<const hundun::linear::LinearSolver &>(pressure_cg);
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
  const hundun::linear::SolveControl default_control{};
  const hundun::linear::SolveControl warped_control{
      1.0e-15, 1.0e-13, 500U, 20U};
  const auto &solve_control =
      analytic_warped ? warped_control : default_control;
  double max_continuity = 0.0;
  double max_pressure_residual = 0.0;
  std::array<double, 3> max_momentum_residual{};
  std::array<double, 3> max_momentum_conservation_defect{};
  std::uint64_t max_solve_iterations = 0U;
  static std::array<bool, 2> cancellation_contract_checked{};
  bool &mapping_cancellation_contract_checked =
      cancellation_contract_checked[analytic_warped ? 1U : 0U];
  for (int step = 0; step < steps; ++step) {
    const auto order = step == 0
                           ? hundun::flow::MomentumTimeOrder::backward_euler
                           : hundun::flow::MomentumTimeOrder::bdf2;
    const auto stencil = hundun::flow::make_momentum_time_stencil(
        order, dt_s, step == 0 ? 0.0 : dt_s);
    const bool check_cancellation_contract =
        !mapping_cancellation_contract_checked &&
        cells_xy == 64 && step == 0;
    double initial_zero_net_scale = 0.0;
    if (check_cancellation_contract) {
      using TestAccess =
          hundun::flow::test::ConstantDensityPisoTestAccess;
      HUNDUN_CHECK(
          hundun::test::flow_state_equality_oracle_is_mutation_sensitive());
      const auto rollback_baseline = capture_state(state);
      for (hundun::mesh::LocalCellId cell = 0;
           cell < topology.owned_cell_count(); ++cell) {
        initial_zero_net_scale +=
            std::abs(kRho * geometry.cell_volume_m3(cell) *
                     initial.velocity[cell * 3U]);
      }
      mpi.allreduce_fp64_in_place(
          &initial_zero_net_scale, 1U,
          hundun::runtime::Fp64ReductionOperation::sum);
      HUNDUN_CHECK(std::isfinite(initial_zero_net_scale));
      HUNDUN_CHECK(initial_zero_net_scale > 0.0);

      TestAccess::reset();
      constexpr double kInjectedRelativeDefect = 6.0e-11;
      TestAccess::set_final_mass_defect_perturbation(
          kInjectedRelativeDefect * 4.0 * kPi * kPi);
      const auto threshold_report =
          flow.attempt(state, kRho, kRho * kNu, stencil, solve_control,
                       solve_control);
      HUNDUN_CHECK(threshold_report.disposition ==
                   hundun::flow::StepAttemptDisposition::recoverable_failure);
      HUNDUN_CHECK(threshold_report.reason ==
                   hundun::flow::StepFailureReason::
                       final_conservation_defect);
      HUNDUN_CHECK(threshold_report.pressure_corrector_count == 2U);
      HUNDUN_CHECK(
          threshold_report.final_mass_relative_conservation_defect > 5.0e-11);
      HUNDUN_CHECK(
          threshold_report.final_mass_relative_conservation_defect < 1.0e-10);
      HUNDUN_CHECK(state.metadata().step == 0U);
      HUNDUN_CHECK(
          state_bitwise_equal(rollback_baseline, capture_state(state)));

      TestAccess::reset();
      constexpr double kZeroNetMutationRatio = 1.0e-10;
      const double injected_quantity =
          kZeroNetMutationRatio * initial_zero_net_scale;
      if (mpi.rank() == 0) {
        TestAccess::force_final_momentum_perturbation(
            0U, injected_quantity /
                    (kRho * geometry.cell_volume_m3(0U)));
      }
      for (std::size_t component = 0; component < 3U; ++component) {
        TestAccess::set_final_momentum_norm_squares(component, 0.0, 1.0);
      }
      const auto mutation_report =
          flow.attempt(state, kRho, kRho * kNu, stencil, solve_control,
                       solve_control);
      HUNDUN_CHECK(mutation_report.disposition ==
                   hundun::flow::StepAttemptDisposition::recoverable_failure);
      HUNDUN_CHECK(mutation_report.reason ==
                   hundun::flow::StepFailureReason::
                       final_conservation_defect);
      HUNDUN_CHECK(
          state_bitwise_equal(rollback_baseline, capture_state(state)));
      HUNDUN_CHECK(mutation_report.pressure_corrector_count == 2U);
      const auto mutation_diagnostic =
          TestAccess::last_momentum_conservation(0U);
      const double raw_over_cancellation_scale =
          std::abs(mutation_diagnostic.raw_defect) / initial_zero_net_scale;
      HUNDUN_CHECK(raw_over_cancellation_scale > 0.5e-10);
      HUNDUN_CHECK(raw_over_cancellation_scale < 1.5e-10);
      const double physical_scale =
          std::max({std::abs(mutation_diagnostic.quantity_n),
                    std::abs(mutation_diagnostic.quantity_np1),
                    stencil.dt_s *
                        mutation_diagnostic.absolute_boundary_flux /
                        stencil.alpha0,
                    std::abs(mutation_diagnostic.history_correction)});
      const double expected_denominator =
          physical_scale <=
                  64.0 * std::numeric_limits<double>::epsilon() *
                      initial_zero_net_scale
              ? initial_zero_net_scale
              : std::max(physical_scale,
                         std::numeric_limits<double>::min());
      const double expected_mutation_relative =
          std::abs(mutation_diagnostic.raw_defect) / expected_denominator;
      HUNDUN_CHECK_NEAR(
          mutation_report.final_momentum_relative_conservation_defect[0],
          expected_mutation_relative,
          512.0 * std::numeric_limits<double>::epsilon());
      HUNDUN_CHECK(
          mutation_report.final_momentum_relative_conservation_defect[0] >
          5.0e-11);
      HUNDUN_CHECK(state.metadata().step == 0U);
      if (mpi.rank() == 0) {
        std::cout
                  << (analytic_warped
                          ? "TASK25_CURVED_TGV_ZERO_NET_MUTATION raw="
                          : "TASK18_TGV_ZERO_NET_MUTATION raw=")
                  << mutation_diagnostic.raw_defect
                  << " cq=" << initial_zero_net_scale
                  << " denominator=" << expected_denominator
                  << " relative="
                  << mutation_report
                         .final_momentum_relative_conservation_defect[0]
                  << " reason=" << static_cast<int>(mutation_report.reason)
                  << " rank=" << mutation_report.lowest_failing_rank << '\n';
      }
      TestAccess::reset();
    }
    const auto report =
        flow.attempt(state, kRho, kRho * kNu, stencil, solve_control,
                     solve_control);
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
                << " qx="
                << hundun::flow::test::ConstantDensityPisoTestAccess::
                       last_momentum_conservation(0U)
                       .quantity_n
                << ','
                << hundun::flow::test::ConstantDensityPisoTestAccess::
                       last_momentum_conservation(0U)
                       .quantity_np1
                << " rawx="
                << hundun::flow::test::ConstantDensityPisoTestAccess::
                       last_momentum_conservation(0U)
                       .raw_defect
                << '\n';
    }
    HUNDUN_CHECK(report.disposition ==
                 hundun::flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(report.pressure_corrector_count == 2U);
    const auto verify_solve = [&](const hundun::linear::SolveReport &solve) {
      HUNDUN_CHECK(
          solve.reason == hundun::linear::SolveTerminationReason::converged ||
          solve.reason ==
              hundun::linear::SolveTerminationReason::zero_right_hand_side);
      HUNDUN_CHECK(solve.iterations <= solve_control.max_iterations);
      HUNDUN_CHECK(std::isfinite(solve.final_residual));
      max_solve_iterations =
          std::max(max_solve_iterations, solve.iterations);
    };
    for (const auto &solve : report.momentum.components)
      verify_solve(solve);
    for (const auto &solve : report.pressure)
      verify_solve(solve);
    HUNDUN_CHECK(report.final_continuity_normalized_l2 <= 1.0e-10);
    if (check_cancellation_contract) {
      const auto committed =
          state.snapshot(hundun::flow::FlowLayer::committed);
      double trial_zero_net_scale = 0.0;
      for (hundun::mesh::LocalCellId cell = 0;
           cell < topology.owned_cell_count(); ++cell) {
        trial_zero_net_scale +=
            std::abs(kRho * geometry.cell_volume_m3(cell) *
                     committed.velocity[cell * 3U]);
      }
      mpi.allreduce_fp64_in_place(
          &trial_zero_net_scale, 1U,
          hundun::runtime::Fp64ReductionOperation::sum);
      const double expected_scale =
          std::max(initial_zero_net_scale, trial_zero_net_scale);
      const auto natural_diagnostic =
          hundun::flow::test::ConstantDensityPisoTestAccess::
              last_momentum_conservation(0U);
      const double expected_relative =
          std::abs(natural_diagnostic.raw_defect) / expected_scale;
      HUNDUN_CHECK_NEAR(
          report.final_momentum_relative_conservation_defect[0],
          expected_relative,
          512.0 * std::numeric_limits<double>::epsilon());
      HUNDUN_CHECK(
          report.final_momentum_relative_conservation_defect[0] <= 5.0e-11);
      if (mpi.rank() == 0) {
        std::cout
                  << (analytic_warped
                          ? "TASK25_CURVED_TGV_ZERO_NET_NATURAL raw="
                          : "TASK18_TGV_ZERO_NET_NATURAL raw=")
                  << natural_diagnostic.raw_defect << " cq=" << expected_scale
                  << " denominator=" << expected_scale << " relative="
                  << report.final_momentum_relative_conservation_defect[0]
                  << '\n';
      }
      mapping_cancellation_contract_checked = true;
    }
    max_continuity =
        std::max(max_continuity, report.final_continuity_normalized_l2);
    max_pressure_residual =
        std::max(max_pressure_residual, report.final_pressure_residual_l2);
    for (std::size_t component = 0; component < 3U; ++component) {
      HUNDUN_CHECK(report.final_momentum_normalized_l2[component] <= 1.0e-9);
      max_momentum_residual[component] =
          std::max(max_momentum_residual[component],
                   report.final_momentum_normalized_l2[component]);
      HUNDUN_CHECK(
          report.final_momentum_relative_conservation_defect[component] <=
          5.0e-11);
      max_momentum_conservation_defect[component] = std::max(
          max_momentum_conservation_defect[component],
          report.final_momentum_relative_conservation_defect[component]);
    }
  }

  const auto result = state.snapshot(hundun::flow::FlowLayer::committed);
  double sums[5]{};
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
    sums[4] += volume * result.mechanical_pressure[cell];
  }
  mpi.allreduce_fp64_in_place(sums, 5U,
                              hundun::runtime::Fp64ReductionOperation::sum);
  const auto metadata = state.metadata();
  TrajectoryResult trajectory;
  trajectory.velocity = result.velocity;
  trajectory.committed = result;
  trajectory.cell_ids.reserve(topology.owned_cell_count());
  for (hundun::mesh::LocalCellId cell = 0;
       cell < topology.owned_cell_count(); ++cell)
    trajectory.cell_ids.push_back(topology.global_cell_id(cell));
  trajectory.face_ids.reserve(topology.owned_face_count());
  trajectory.face_local_ids.reserve(topology.owned_face_count());
  for (hundun::mesh::LocalFaceId face = 0;
       face < topology.local_face_count(); ++face) {
    if (topology.face_ownership(face) !=
        hundun::mesh::EntityOwnership::owned)
      continue;
    trajectory.face_ids.push_back(topology.global_face_id(face));
    trajectory.face_local_ids.push_back(face);
  }
  HUNDUN_CHECK(trajectory.face_ids.size() == topology.owned_face_count());
  HUNDUN_CHECK(trajectory.face_local_ids.size() ==
               topology.owned_face_count());
  HUNDUN_CHECK(std::set<hundun::mesh::GlobalFaceId>(
                   trajectory.face_ids.begin(), trajectory.face_ids.end())
                   .size() == trajectory.face_ids.size());
  trajectory.pressure_mean = sums[4] / sums[3];
  trajectory.cell_volume_m3 = geometry.cell_volume_m3(0U);
  trajectory.velocity_l2_error = std::sqrt(sums[0] / sums[1]);
  trajectory.relative_mass_defect =
      std::abs(sums[2] - kRho * sums[3]) / (kRho * sums[3]);
  trajectory.max_continuity = max_continuity;
  trajectory.max_pressure_residual = max_pressure_residual;
  trajectory.max_momentum_residual = max_momentum_residual;
  trajectory.max_momentum_conservation_defect =
      max_momentum_conservation_defect;
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &max_solve_iterations, 1,
                             MPI_UINT64_T, MPI_MAX, mpi.comm()) ==
               MPI_SUCCESS);
  trajectory.max_solve_iterations = max_solve_iterations;
  trajectory.committed_steps = metadata.step;
  trajectory.final_order = metadata.order;
  trajectory.final_dt_s = metadata.dt_s;
  trajectory.final_previous_dt_s = metadata.previous_dt_s;
  return trajectory;
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

std::array<double, 10> decomposition_field_differences(
    const hundun::runtime::MpiContext &mpi,
    const TrajectoryResult &distributed,
    const TrajectoryResult &single_rank) {
  constexpr std::size_t field_count = 10U;
  std::array<double, field_count * 2U> values{};
  const auto compare = [&](std::size_t field, double actual,
                           double reference) {
    values[field] = std::max(values[field], std::abs(actual - reference));
    values[field + field_count] =
        std::max(values[field + field_count], std::abs(reference));
  };
  std::vector<std::size_t> single_face_local;
  for (const auto id : single_rank.face_ids) {
    const auto needed = static_cast<std::size_t>(id) + 1U;
    if (single_face_local.size() < needed)
      single_face_local.resize(needed,
                               std::numeric_limits<std::size_t>::max());
  }
  for (std::size_t local = 0; local < single_rank.face_ids.size(); ++local)
  {
    const auto global =
        static_cast<std::size_t>(single_rank.face_ids[local]);
    HUNDUN_CHECK(single_face_local[global] ==
                 std::numeric_limits<std::size_t>::max());
    single_face_local[global] =
        static_cast<std::size_t>(single_rank.face_local_ids[local]);
  }
  HUNDUN_CHECK(single_face_local.size() <=
               static_cast<std::size_t>(std::numeric_limits<int>::max()));
  std::vector<int> distributed_ownership(single_face_local.size(), 0);
  for (const auto id : distributed.face_ids) {
    const auto global = static_cast<std::size_t>(id);
    HUNDUN_CHECK(global < distributed_ownership.size());
    HUNDUN_CHECK(distributed_ownership[global] == 0);
    distributed_ownership[global] = 1;
  }
  HUNDUN_CHECK(MPI_Allreduce(
                   MPI_IN_PLACE, distributed_ownership.data(),
                   static_cast<int>(distributed_ownership.size()), MPI_INT,
                   MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  for (const auto id : single_rank.face_ids)
    HUNDUN_CHECK(distributed_ownership[static_cast<std::size_t>(id)] == 1);
  HUNDUN_CHECK(distributed.committed.transported_cell_fields.size() ==
               single_rank.committed.transported_cell_fields.size());
  for (std::size_t nested = 0;
       nested < distributed.committed.transported_cell_fields.size();
       ++nested) {
    HUNDUN_CHECK(
        distributed.committed.transported_cell_fields[nested].size() ==
        distributed.committed.density.size());
    HUNDUN_CHECK(
        single_rank.committed.transported_cell_fields[nested].size() ==
        single_rank.committed.density.size());
  }
  for (std::size_t local = 0; local < distributed.cell_ids.size(); ++local) {
    const auto global =
        static_cast<std::size_t>(distributed.cell_ids[local]);
    compare(0U, distributed.committed.density[local],
            single_rank.committed.density[global]);
    for (std::size_t component = 0; component < 3U; ++component)
      compare(1U + component,
              distributed.committed.velocity[local * 3U + component],
              single_rank.committed.velocity[global * 3U + component]);
    compare(4U,
            distributed.committed.mechanical_pressure[local] -
                distributed.pressure_mean,
            single_rank.committed.mechanical_pressure[global] -
                single_rank.pressure_mean);
    for (std::size_t nested = 0;
         nested < distributed.committed.transported_cell_fields.size();
         ++nested)
      HUNDUN_CHECK_NEAR(
          distributed.committed.transported_cell_fields[nested][local],
          single_rank.committed.transported_cell_fields[nested][global],
          5.0e-12 *
              std::max(
                  1.0,
                  std::abs(single_rank.committed
                               .transported_cell_fields[nested][global])));
  }
  for (std::size_t local = 0; local < distributed.face_ids.size(); ++local) {
    const auto global =
        static_cast<std::size_t>(distributed.face_ids[local]);
    const auto distributed_local =
        static_cast<std::size_t>(distributed.face_local_ids[local]);
    HUNDUN_CHECK(distributed_local * 3U + 2U <
                 distributed.committed.face_velocity.size());
    HUNDUN_CHECK(distributed_local <
                 distributed.committed.face_mass_flux.size());
    HUNDUN_CHECK(global < single_face_local.size());
    const auto reference_local = single_face_local[global];
    HUNDUN_CHECK(reference_local !=
                 std::numeric_limits<std::size_t>::max());
    HUNDUN_CHECK(reference_local * 3U + 2U <
                 single_rank.committed.face_velocity.size());
    HUNDUN_CHECK(reference_local <
                 single_rank.committed.face_mass_flux.size());
    for (std::size_t component = 0; component < 3U; ++component)
      compare(5U + component,
              distributed.committed
                  .face_velocity[distributed_local * 3U + component],
              single_rank.committed
                  .face_velocity[reference_local * 3U + component]);
    compare(8U, distributed.committed.face_mass_flux[distributed_local],
            single_rank.committed.face_mass_flux[reference_local]);
  }
  values[9U] =
      std::abs(distributed.velocity_l2_error - single_rank.velocity_l2_error);
  values[9U + field_count] = std::abs(single_rank.velocity_l2_error);
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, values.data(),
                             static_cast<int>(values.size()), MPI_DOUBLE,
                             MPI_MAX, mpi.comm()) == MPI_SUCCESS);
  std::array<double, field_count> differences{};
  static constexpr std::array<const char *, field_count> names{
      "density",        "velocity_x",     "velocity_y",
      "velocity_z",     "pressure",       "face_velocity_x",
      "face_velocity_y", "face_velocity_z", "face_mass_flux",
      "result"};
  for (std::size_t field = 0; field < differences.size(); ++field) {
    differences[field] = values[field];
    if (mpi.rank() == 0) {
      std::cout << "TASK25_TAYLOR_GREEN_DECOMPOSITION_FIELD name="
                << names[field] << " difference=" << differences[field]
                << " reference_scale=" << values[field + field_count] << '\n';
    }
    HUNDUN_CHECK(differences[field] <=
                 5.0e-12 *
                     std::max(1.0, values[field + field_count]));
  }
  return differences;
}

void require_outside_prefix_face_mutation_rejected(
    const hundun::runtime::MpiContext &mpi,
    const TrajectoryResult &distributed,
    const TrajectoryResult &single_rank) {
  const auto missing = std::numeric_limits<std::uint64_t>::max();
  std::uint64_t local_global = missing;
  hundun::mesh::LocalFaceId local_face =
      std::numeric_limits<hundun::mesh::LocalFaceId>::max();
  for (std::size_t index = 0;
       index < distributed.face_ids.size(); ++index) {
    const auto candidate_local = distributed.face_local_ids[index];
    if (candidate_local < distributed.face_ids.size())
      continue;
    const auto candidate_global = static_cast<std::uint64_t>(
        distributed.face_ids[index]);
    if (candidate_global < local_global) {
      local_global = candidate_global;
      local_face = candidate_local;
    }
  }
  std::uint64_t selected_global = local_global;
  HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &selected_global, 1, MPI_UINT64_T,
                             MPI_MIN, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(selected_global != missing);

  auto mutated = distributed;
  if (local_global == selected_global) {
    HUNDUN_CHECK(local_face < mutated.committed.face_mass_flux.size());
    mutated.committed.face_mass_flux[local_face] += 1.0;
  }
  bool rejected = false;
  try {
    static_cast<void>(
        decomposition_field_differences(mpi, mutated, single_rank));
  } catch (const std::runtime_error &) {
    rejected = true;
  }
  const int local_rejected = rejected ? 1 : 0;
  int rejected_count = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_rejected, &rejected_count, 1, MPI_INT,
                             MPI_SUM, mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(rejected_count == mpi.size());
}

void run_taylor_green(const hundun::runtime::MpiContext &mpi) {
  constexpr std::array<int, 3> spatial_cells{16, 32, 64};
  const std::array<double, 3> warped_spatial_dt{
      1.0e-4 * std::pow(2.0 * kPi / spatial_cells[0], 2) / kNu,
      1.0e-4 * std::pow(2.0 * kPi / spatial_cells[1], 2) / kNu,
      1.0e-4 * std::pow(2.0 * kPi / spatial_cells[2], 2) / kNu};
  std::array<TrajectoryResult, 3> warped_spatial{
      run_trajectory(mpi, spatial_cells[0], warped_spatial_dt[0],
                     warped_spatial_dt[0], true),
      run_trajectory(mpi, spatial_cells[1], warped_spatial_dt[1],
                     warped_spatial_dt[1], true),
      run_trajectory(mpi, spatial_cells[2], warped_spatial_dt[2],
                     warped_spatial_dt[2], true)};
  const double warped_spatial_order_0 =
      std::log(warped_spatial[0].velocity_l2_error /
               warped_spatial[1].velocity_l2_error) /
      std::log(2.0);
  const double warped_spatial_order_1 =
      std::log(warped_spatial[1].velocity_l2_error /
               warped_spatial[2].velocity_l2_error) /
      std::log(2.0);

  std::array<std::array<double, 10>, 3> decomposition_differences{};
  if (mpi.size() > 1) {
    auto self =
        hundun::runtime::MpiContext::duplicate(MPI_COMM_SELF);
    for (std::size_t level = 0; level < warped_spatial.size(); ++level) {
      const auto single_rank =
          run_trajectory(self, spatial_cells[level],
                         warped_spatial_dt[level],
                         warped_spatial_dt[level], true);
      decomposition_differences[level] =
          decomposition_field_differences(
              mpi, warped_spatial[level], single_rank);
      if (level == 0U)
        require_outside_prefix_face_mutation_rejected(
            mpi, warped_spatial[level], single_rank);
    }
  }

  constexpr double uniform_spatial_dt = 1.0e-4;
  std::array<TrajectoryResult, 3> uniform_spatial{
      run_trajectory(mpi, spatial_cells[0], uniform_spatial_dt,
                     uniform_spatial_dt),
      run_trajectory(mpi, spatial_cells[1], uniform_spatial_dt,
                     uniform_spatial_dt),
      run_trajectory(mpi, spatial_cells[2], uniform_spatial_dt,
                     uniform_spatial_dt)};
  const double uniform_spatial_order_0 =
      std::log(uniform_spatial[0].velocity_l2_error /
               uniform_spatial[1].velocity_l2_error) /
      std::log(2.0);
  const double uniform_spatial_order_1 =
      std::log(uniform_spatial[1].velocity_l2_error /
               uniform_spatial[2].velocity_l2_error) /
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
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (std::size_t level = 0; level < uniform_spatial.size(); ++level) {
      std::cout << "TASK18_TAYLOR_GREEN_SPATIAL cells=" << spatial_cells[level]
                << " dt=" << uniform_spatial_dt
                << " velocity_l2="
                << uniform_spatial[level].velocity_l2_error
                << " mass_defect="
                << uniform_spatial[level].relative_mass_defect
                << " continuity=" << uniform_spatial[level].max_continuity
                << " pressure_residual="
                << uniform_spatial[level].max_pressure_residual
                << " max_solve_iterations="
                << uniform_spatial[level].max_solve_iterations
                << " momentum_residual="
                << uniform_spatial[level].max_momentum_residual[0] << ','
                << uniform_spatial[level].max_momentum_residual[1] << ','
                << uniform_spatial[level].max_momentum_residual[2]
                << " momentum_defect="
                << uniform_spatial[level]
                       .max_momentum_conservation_defect[0]
                << ','
                << uniform_spatial[level]
                       .max_momentum_conservation_defect[1]
                << ','
                << uniform_spatial[level]
                       .max_momentum_conservation_defect[2]
                << '\n';
    }
    std::cout << "TASK18_TAYLOR_GREEN_SPATIAL_ORDER coarse="
              << uniform_spatial_order_0
              << " fine=" << uniform_spatial_order_1 << '\n';
    for (std::size_t level = 0; level < warped_spatial.size(); ++level) {
      std::cout << "TASK25_CURVED_TAYLOR_GREEN_SPATIAL cells="
                << spatial_cells[level]
                << " dt=" << warped_spatial_dt[level]
                << " velocity_l2="
                << warped_spatial[level].velocity_l2_error
                << " mass_defect="
                << warped_spatial[level].relative_mass_defect
                << " continuity=" << warped_spatial[level].max_continuity
                << " pressure_residual="
                << warped_spatial[level].max_pressure_residual
                << " max_solve_iterations="
                << warped_spatial[level].max_solve_iterations
                << " momentum_residual="
                << warped_spatial[level].max_momentum_residual[0] << ','
                << warped_spatial[level].max_momentum_residual[1] << ','
                << warped_spatial[level].max_momentum_residual[2]
                << " momentum_defect="
                << warped_spatial[level]
                       .max_momentum_conservation_defect[0]
                << ','
                << warped_spatial[level]
                       .max_momentum_conservation_defect[1]
                << ','
                << warped_spatial[level]
                       .max_momentum_conservation_defect[2]
                << '\n';
    }
    std::cout << "TASK25_CURVED_TAYLOR_GREEN_SPATIAL_ORDER coarse="
              << warped_spatial_order_0
              << " fine=" << warped_spatial_order_1 << '\n';
    if (mpi.size() > 1) {
      for (std::size_t level = 0; level < warped_spatial.size(); ++level) {
        std::cout << "TASK25_TAYLOR_GREEN_DECOMPOSITION cells="
                  << spatial_cells[level] << " density="
                  << decomposition_differences[level][0] << " velocity="
                  << decomposition_differences[level][1] << ','
                  << decomposition_differences[level][2] << ','
                  << decomposition_differences[level][3] << " pressure="
                  << decomposition_differences[level][4]
                  << " face_velocity="
                  << decomposition_differences[level][5] << ','
                  << decomposition_differences[level][6] << ','
                  << decomposition_differences[level][7]
                  << " face_mass_flux="
                  << decomposition_differences[level][8] << " result="
                  << decomposition_differences[level][9] << '\n';
      }
    }
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

  HUNDUN_CHECK(uniform_spatial_order_0 >= 1.8);
  HUNDUN_CHECK(uniform_spatial_order_1 >= 1.8);
  HUNDUN_CHECK(warped_spatial_order_0 >= 1.8);
  HUNDUN_CHECK(warped_spatial_order_1 >= 1.8);
  HUNDUN_CHECK(temporal_order >= 1.8);
  for (const auto &item : uniform_spatial) {
    HUNDUN_CHECK(item.relative_mass_defect <= 5.0e-12);
    HUNDUN_CHECK(item.committed_steps == 1U);
    HUNDUN_CHECK(item.final_order ==
                 hundun::flow::MomentumTimeOrder::backward_euler);
  }
  for (const auto &item : warped_spatial) {
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

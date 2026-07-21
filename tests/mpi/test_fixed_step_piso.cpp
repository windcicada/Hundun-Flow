// SPDX-License-Identifier: Apache-2.0

#include "flow/src/constant_density_piso_test_access.hpp"
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
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

hundun::runtime::Int3 process_grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw hundun::runtime::Error("unsupported Task 18 rank count");
}

hundun::config::FlowCaseConfig periodic_case(double rho_ref) {
  hundun::config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type =
      hundun::config::SimulationType::variable_density_flow;
  config.density_model = hundun::config::DensityModel::constant;
  config.physics.rho_ref_kg_per_m3 = rho_ref;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  config.scalars.push_back(hundun::config::FlowScalarConfig{"alpha", 0.0});
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
          "task18",
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
          "task18",
          hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64,
          components,
          0,
          false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

class SecondComponentFailureSolver final : public hundun::linear::LinearSolver {
public:
  hundun::linear::SolveReport
  solve(const hundun::linear::LinearOperator &,
        hundun::linear::Preconditioner &,
        hundun::execution::VectorView<const double>,
        hundun::execution::VectorView<double>,
        const hundun::linear::SolveControl &) const override {
    hundun::linear::SolveReport report;
    if (calls_++ == 1U) {
      report.reason =
          hundun::linear::SolveTerminationReason::maximum_iterations;
      report.lowest_failing_rank = 0;
    } else {
      report.reason = hundun::linear::SolveTerminationReason::converged;
      report.final_residual = 0.0;
      report.lowest_failing_rank = -1;
    }
    return report;
  }

private:
  mutable std::size_t calls_{};
};

double checkerboard_amplitude(const hundun::runtime::MpiContext &mpi,
                              const hundun::mesh::MeshTopology &topology,
                              const hundun::flow::FlowLayerValues &values) {
  double parity[2]{};
  for (hundun::mesh::LocalCellId cell = 0; cell < topology.owned_cell_count();
       ++cell) {
    const auto global = topology.global_cell(cell);
    const double sign =
        ((global.x + global.y + global.z) % 2 == 0) ? 1.0 : -1.0;
    parity[0] += sign * values.mechanical_pressure[cell];
    parity[1] += 1.0;
  }
  mpi.allreduce_fp64_in_place(parity, 2U,
                              hundun::runtime::Fp64ReductionOperation::sum);
  return std::abs(parity[0] / parity[1]);
}

void run_zero_flow_transaction(const hundun::runtime::MpiContext &mpi) {
  constexpr double rho_ref = 1.25;
  constexpr hundun::runtime::Int3 extent{8, 6, 4};
  auto decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, extent, {true, true, true},
      hundun::runtime::DecompositionOptions{process_grid(mpi.size())});
  const hundun::mesh::MeshTopology topology(decomposition);
  const hundun::mesh::MeshGeometry geometry(
      topology,
      hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  auto boundaries = hundun::boundary::BoundaryRegistry::create(
      periodic_case(rho_ref), topology);

  hundun::runtime::FieldRegistry registry;
  hundun::flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell_field("rho", 1U));
  fields.velocity = registry.declare_field(cell_field("velocity", 3U));
  fields.mechanical_pressure = registry.declare_field(cell_field("pi", 1U));
  fields.face_velocity = registry.declare_field(face_field("u_face", 3U));
  fields.face_mass_flux =
      hundun::finite_volume::declare_face_mass_flux(registry);
  fields.transported_cell_fields.push_back(
      registry.declare_field(cell_field("alpha", 1U)));
  registry.freeze();

  const auto box = decomposition.owned_box();
  const hundun::runtime::Int3 local{box.end.x - box.begin.x,
                                    box.end.y - box.begin.y,
                                    box.end.z - box.begin.z};
  auto state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  const std::size_t cells = topology.owned_cell_count();
  hundun::flow::FlowLayerValues initial;
  initial.density.assign(cells, rho_ref);
  initial.velocity.assign(cells * 3U, 0.0);
  initial.mechanical_pressure.assign(cells, 0.0);
  initial.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  initial.face_mass_flux.assign(topology.local_face_count(), 0.0);
  initial.transported_cell_fields.resize(1U);
  initial.transported_cell_fields[0].resize(cells);
  for (hundun::mesh::LocalCellId cell = 0; cell < topology.owned_cell_count();
       ++cell) {
    const auto center = geometry.cell_center_m(cell);
    initial.transported_cell_fields[0][cell] =
        1.0 + 0.2 * std::sin(2.0 * 3.14159265358979323846 * center.x);
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
      pressure_preconditioner,
      {{fields.transported_cell_fields[0],
        hundun::finite_volume::FiniteVolumeQuantity::scalar(0U), 0.0}});
  const auto stencil = hundun::flow::make_momentum_time_stencil(
      hundun::flow::MomentumTimeOrder::backward_euler, 0.01, 0.0);
  const auto report = flow.attempt(state, rho_ref, 0.01, stencil, {}, {});
  if (report.disposition != hundun::flow::StepAttemptDisposition::committed) {
    std::cerr << "TASK18_ZERO_REPORT rank=" << mpi.rank()
              << " disposition=" << static_cast<int>(report.disposition)
              << " reason=" << static_cast<int>(report.reason)
              << " failing_rank=" << report.lowest_failing_rank
              << " correctors=" << report.pressure_corrector_count
              << " continuity=" << report.final_continuity_normalized_l2
              << " pressure=" << report.final_pressure_residual_l2 << '\n';
  }
  HUNDUN_CHECK(report.disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(report.pressure_corrector_count == 2U);
  HUNDUN_CHECK(report.final_continuity_normalized_l2 <= 1.0e-10);
  HUNDUN_CHECK(state.metadata().step == 1U);
  HUNDUN_CHECK(state.metadata().time_s == 0.01);
  HUNDUN_CHECK(state.snapshot(hundun::flow::FlowLayer::committed)
                   .transported_cell_fields[0] ==
               initial.transported_cell_fields[0]);

  using TestAccess = hundun::flow::test::ConstantDensityPisoTestAccess;
  TestAccess::reset();
  TestAccess::set_provisional_transport_sentinel(true);
  auto sentinel_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  sentinel_state.seed_accepted_layers(initial, initial);
  const auto sentinel_report =
      flow.attempt(sentinel_state, rho_ref, 0.0, stencil, {}, {});
  HUNDUN_CHECK(sentinel_report.disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(TestAccess::provisional_transport_calls() == 1U);
  HUNDUN_CHECK(TestAccess::final_transport_calls() == 1U);
  HUNDUN_CHECK(sentinel_state.snapshot(hundun::flow::FlowLayer::committed)
                   .transported_cell_fields[0] ==
               initial.transported_cell_fields[0]);

  TestAccess::reset();
  TestAccess::set_final_uniform_x_mass_flux(0.05);
  auto provenance_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  provenance_state.seed_accepted_layers(initial, initial);
  const auto provenance_report =
      flow.attempt(provenance_state, rho_ref, 0.0, stencil, {}, {});
  HUNDUN_CHECK(provenance_report.disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(provenance_report.final_continuity_normalized_l2 <= 1.0e-10);
  const auto advected =
      provenance_state.snapshot(hundun::flow::FlowLayer::committed)
          .transported_cell_fields[0];
  double final_flux_sensitivity = 0.0;
  for (std::size_t cell = 0; cell < advected.size(); ++cell) {
    final_flux_sensitivity = std::max(
        final_flux_sensitivity,
        std::abs(advected[cell] - initial.transported_cell_fields[0][cell]));
  }
  HUNDUN_CHECK(final_flux_sensitivity > 1.0e-8);

  TestAccess::reset();
  TestAccess::force_final_continuity_failure(true);
  auto rollback_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  rollback_state.seed_accepted_layers(initial, initial);
  const auto before_rollback =
      rollback_state.snapshot(hundun::flow::FlowLayer::committed);
  const auto rollback_report =
      flow.attempt(rollback_state, rho_ref, 0.0, stencil, {}, {});
  HUNDUN_CHECK(rollback_report.disposition ==
               hundun::flow::StepAttemptDisposition::recoverable_failure);
  HUNDUN_CHECK(rollback_report.reason ==
               hundun::flow::StepFailureReason::final_continuity_residual);
  HUNDUN_CHECK(rollback_report.pressure_corrector_count == 2U);
  HUNDUN_CHECK(rollback_report.suggested_dt_s == 0.005);
  HUNDUN_CHECK(rollback_state.metadata().step == 0U);
  const auto after_rollback =
      rollback_state.snapshot(hundun::flow::FlowLayer::committed);
  HUNDUN_CHECK(after_rollback.density == before_rollback.density);
  HUNDUN_CHECK(after_rollback.velocity == before_rollback.velocity);
  HUNDUN_CHECK(after_rollback.mechanical_pressure ==
               before_rollback.mechanical_pressure);
  HUNDUN_CHECK(after_rollback.face_velocity == before_rollback.face_velocity);
  HUNDUN_CHECK(after_rollback.face_mass_flux == before_rollback.face_mass_flux);
  HUNDUN_CHECK(after_rollback.transported_cell_fields ==
               before_rollback.transported_cell_fields);
  HUNDUN_CHECK(TestAccess::provisional_transport_calls() == 1U);
  HUNDUN_CHECK(TestAccess::final_transport_calls() == 1U);

  TestAccess::reset();
  TestAccess::force_final_pressure_failure(true);
  auto pressure_rollback_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  pressure_rollback_state.seed_accepted_layers(initial, initial);
  const auto pressure_before =
      pressure_rollback_state.snapshot(hundun::flow::FlowLayer::committed);
  const auto pressure_rollback_report =
      flow.attempt(pressure_rollback_state, rho_ref, 0.0, stencil, {}, {});
  HUNDUN_CHECK(pressure_rollback_report.disposition ==
               hundun::flow::StepAttemptDisposition::recoverable_failure);
  HUNDUN_CHECK(pressure_rollback_report.reason ==
               hundun::flow::StepFailureReason::final_pressure_residual);
  HUNDUN_CHECK(pressure_rollback_report.pressure_corrector_count == 2U);
  HUNDUN_CHECK(pressure_rollback_report.suggested_dt_s == 0.005);
  HUNDUN_CHECK(pressure_rollback_state.metadata().step == 0U);
  const auto pressure_after =
      pressure_rollback_state.snapshot(hundun::flow::FlowLayer::committed);
  HUNDUN_CHECK(pressure_after.density == pressure_before.density);
  HUNDUN_CHECK(pressure_after.velocity == pressure_before.velocity);
  HUNDUN_CHECK(pressure_after.mechanical_pressure ==
               pressure_before.mechanical_pressure);

  if (mpi.rank() == 0) {
    std::cout << "TASK18_TRANSACTION final_flux_sensitivity="
              << final_flux_sensitivity << " rollback_correctors="
              << rollback_report.pressure_corrector_count
              << " suggested_dt=" << rollback_report.suggested_dt_s
              << " provisional_calls="
              << TestAccess::provisional_transport_calls()
              << " final_calls=" << TestAccess::final_transport_calls() << '\n';
  }
  TestAccess::reset();

  SecondComponentFailureSolver component_failure_solver;
  hundun::linear::JacobiPreconditioner failure_mx(execution);
  hundun::linear::JacobiPreconditioner failure_my(execution);
  hundun::linear::JacobiPreconditioner failure_mz(execution);
  hundun::linear::JacobiPreconditioner failure_pressure(execution);
  auto component_failure_flow =
      hundun::flow::FixedStepConstantDensityFlow::create(
          decomposition, topology, geometry, boundaries, mpi, execution, halo,
          component_failure_solver, {&failure_mx, &failure_my, &failure_mz},
          pressure_solver, failure_pressure,
          {{fields.transported_cell_fields[0],
            hundun::finite_volume::FiniteVolumeQuantity::scalar(0U), 0.0}});
  auto component_failure_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  component_failure_state.seed_accepted_layers(initial, initial);
  const auto component_failure_report = component_failure_flow.attempt(
      component_failure_state, rho_ref, 0.0, stencil, {}, {});
  HUNDUN_CHECK(component_failure_report.disposition ==
               hundun::flow::StepAttemptDisposition::recoverable_failure);
  HUNDUN_CHECK(component_failure_report.reason ==
               hundun::flow::StepFailureReason::momentum_linear_solve);
  HUNDUN_CHECK(component_failure_report.lowest_failing_rank == 0);
  HUNDUN_CHECK(component_failure_report.pressure_corrector_count == 0U);
  HUNDUN_CHECK(component_failure_report.suggested_dt_s == 0.005);
  HUNDUN_CHECK(component_failure_state.metadata().step == 0U);
  HUNDUN_CHECK(
      component_failure_state.snapshot(hundun::flow::FlowLayer::committed)
          .velocity == initial.velocity);

  auto checker_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  auto checker = initial;
  for (hundun::mesh::LocalCellId cell = 0; cell < topology.owned_cell_count();
       ++cell) {
    const auto global = topology.global_cell(cell);
    checker.mechanical_pressure[cell] =
        ((global.x + global.y + global.z) % 2 == 0) ? 1.0 : -1.0;
  }
  checker_state.seed_accepted_layers(checker, checker);
  const auto checker_report =
      flow.attempt(checker_state, rho_ref, 0.0, stencil, {}, {});
  if (checker_report.disposition !=
      hundun::flow::StepAttemptDisposition::committed) {
    std::cerr << "TASK18_CHECKER_REPORT rank=" << mpi.rank()
              << " disposition=" << static_cast<int>(checker_report.disposition)
              << " reason=" << static_cast<int>(checker_report.reason)
              << " correctors=" << checker_report.pressure_corrector_count
              << " pressure0=" << checker_report.pressure[0].final_residual
              << " pressure1=" << checker_report.pressure[1].final_residual
              << " continuity=" << checker_report.final_continuity_normalized_l2
              << '\n';
  }
  HUNDUN_CHECK(checker_report.disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(checker_report.pressure_corrector_count == 2U);
  HUNDUN_CHECK(checker_report.final_continuity_normalized_l2 <= 1.0e-10);
  const auto checker_result =
      checker_state.snapshot(hundun::flow::FlowLayer::committed);
  double parity[2]{};
  for (hundun::mesh::LocalCellId cell = 0; cell < topology.owned_cell_count();
       ++cell) {
    const auto global = topology.global_cell(cell);
    const double sign =
        ((global.x + global.y + global.z) % 2 == 0) ? 1.0 : -1.0;
    parity[0] += sign * checker_result.mechanical_pressure[cell];
    parity[1] += 1.0;
  }
  mpi.allreduce_fp64_in_place(parity, 2U,
                              hundun::runtime::Fp64ReductionOperation::sum);
  const double parity_amplitude = std::abs(parity[0] / parity[1]);
  if (mpi.rank() == 0) {
    std::cout << "TASK18_CHECKERBOARD amplitude=" << parity_amplitude
              << " continuity=" << checker_report.final_continuity_normalized_l2
              << " correctors=" << checker_report.pressure_corrector_count
              << '\n';
  }
  HUNDUN_CHECK(parity_amplitude <= 1.0e-8);

  auto small_dt_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.001, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  small_dt_state.seed_accepted_layers(checker, checker);
  const auto small_dt_stencil = hundun::flow::make_momentum_time_stencil(
      hundun::flow::MomentumTimeOrder::backward_euler, 0.001, 0.0);
  const auto small_dt_report =
      flow.attempt(small_dt_state, rho_ref, 0.0, small_dt_stencil, {}, {});
  HUNDUN_CHECK(small_dt_report.disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(small_dt_report.pressure_corrector_count == 2U);
  const double small_dt_amplitude = checkerboard_amplitude(
      mpi, topology,
      small_dt_state.snapshot(hundun::flow::FlowLayer::committed));
  HUNDUN_CHECK(small_dt_amplitude <= 1.0e-8);
  if (mpi.rank() == 0) {
    std::cout << "TASK18_CHECKERBOARD_SMALL_DT amplitude=" << small_dt_amplitude
              << " continuity="
              << small_dt_report.final_continuity_normalized_l2 << '\n';
  }

  if (mpi.size() > 1) {
    auto mismatched_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
    mismatched_state.seed_accepted_layers(initial, initial);
    auto mismatched_stencil = stencil;
    if (mpi.rank() == 1)
      mismatched_stencil.alpha0 += 0.25;
    const auto mismatched_report = flow.attempt(mismatched_state, rho_ref, 0.0,
                                                mismatched_stencil, {}, {});
    HUNDUN_CHECK(mismatched_report.disposition ==
                 hundun::flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(mismatched_report.reason ==
                 hundun::flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(mismatched_report.lowest_failing_rank == 1);
    HUNDUN_CHECK(mismatched_report.suggested_dt_s == 0.0);
    HUNDUN_CHECK(mismatched_state.metadata().step == 0U);

    auto divergent_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
    divergent_state.seed_accepted_layers(initial, initial);
    const auto divergent_stencil =
        mpi.rank() == 1
            ? hundun::flow::make_momentum_time_stencil(
                  hundun::flow::MomentumTimeOrder::backward_euler, 0.02, 0.0)
            : stencil;
    const auto divergent_report =
        flow.attempt(divergent_state, rho_ref, 0.0, divergent_stencil, {}, {});
    HUNDUN_CHECK(divergent_report.disposition ==
                 hundun::flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(divergent_report.reason ==
                 hundun::flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(divergent_report.lowest_failing_rank == 1);
    HUNDUN_CHECK(divergent_state.metadata().step == 0U);

    auto divergent_transport_flow =
        hundun::flow::FixedStepConstantDensityFlow::create(
            decomposition, topology, geometry, boundaries, mpi, execution, halo,
            momentum_solver, {&mx, &my, &mz}, pressure_solver,
            pressure_preconditioner,
            {{fields.transported_cell_fields[0],
              hundun::finite_volume::FiniteVolumeQuantity::scalar(0U),
              mpi.rank() == 1 ? 0.125 : 0.0}});
    auto divergent_transport_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
    divergent_transport_state.seed_accepted_layers(initial, initial);
    const auto divergent_transport_report = divergent_transport_flow.attempt(
        divergent_transport_state, rho_ref, 0.0, stencil, {}, {});
    HUNDUN_CHECK(divergent_transport_report.disposition ==
                 hundun::flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(divergent_transport_report.reason ==
                 hundun::flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(divergent_transport_report.lowest_failing_rank == 1);
    HUNDUN_CHECK(divergent_transport_state.metadata().step == 0U);

    TestAccess::reset();
    TestAccess::force_local_derived_failure(mpi.rank() == 1);
    auto derived_failure_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
    derived_failure_state.seed_accepted_layers(initial, initial);
    const auto derived_failure_report =
        flow.attempt(derived_failure_state, rho_ref, 0.0, stencil, {}, {});
    HUNDUN_CHECK(derived_failure_report.disposition ==
                 hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(derived_failure_report.reason ==
                 hundun::flow::StepFailureReason::non_finite_trial);
    HUNDUN_CHECK(derived_failure_report.lowest_failing_rank == 1);
    HUNDUN_CHECK(derived_failure_report.pressure_corrector_count == 0U);
    HUNDUN_CHECK(derived_failure_report.suggested_dt_s == 0.005);
    HUNDUN_CHECK(derived_failure_state.metadata().step == 0U);

    TestAccess::reset();
    const auto derived_retry_stencil = hundun::flow::make_momentum_time_stencil(
        hundun::flow::MomentumTimeOrder::backward_euler,
        derived_failure_report.suggested_dt_s, 0.0);
    const auto derived_retry_report = flow.attempt(
        derived_failure_state, rho_ref, 0.0, derived_retry_stencil, {}, {});
    HUNDUN_CHECK(derived_retry_report.disposition ==
                 hundun::flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(derived_retry_report.reason ==
                 hundun::flow::StepFailureReason::none);
    HUNDUN_CHECK(derived_retry_report.pressure_corrector_count == 2U);
    HUNDUN_CHECK(derived_failure_state.metadata().step == 1U);
    HUNDUN_CHECK(derived_failure_state.metadata().time_s == 0.005);

    TestAccess::force_final_continuity_failure(mpi.rank() == 1);
    auto collective_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
    collective_state.seed_accepted_layers(initial, initial);
    const auto collective_report =
        flow.attempt(collective_state, rho_ref, 0.0, stencil, {}, {});
    HUNDUN_CHECK(collective_report.disposition ==
                 hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(collective_report.reason ==
                 hundun::flow::StepFailureReason::final_continuity_residual);
    HUNDUN_CHECK(collective_report.lowest_failing_rank == 1);
    HUNDUN_CHECK(collective_report.pressure_corrector_count == 2U);
    HUNDUN_CHECK(collective_state.metadata().step == 0U);
    TestAccess::reset();
  }

  auto bdf2_committed = initial;
  auto bdf2_history = initial;
  constexpr std::array<double, 3> velocity_n{0.5, 0.1, -0.2};
  constexpr std::array<double, 3> velocity_nm1{0.2, -0.3, 0.4};
  for (std::size_t cell = 0; cell < cells; ++cell) {
    for (std::size_t component = 0; component < 3U; ++component) {
      bdf2_committed.velocity[cell * 3U + component] = velocity_n[component];
      bdf2_history.velocity[cell * 3U + component] = velocity_nm1[component];
    }
  }
  for (hundun::mesh::LocalFaceId face = 0; face < topology.local_face_count();
       ++face) {
    for (std::size_t component = 0; component < 3U; ++component) {
      bdf2_committed.face_velocity[face * 3U + component] =
          velocity_n[component];
      bdf2_history.face_velocity[face * 3U + component] =
          velocity_nm1[component];
    }
    bdf2_committed.face_mass_flux[face] = 0.0;
    bdf2_history.face_mass_flux[face] = 0.0;
  }
  auto bdf2_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {1U, 0.015, 0.015, 0.01, hundun::flow::MomentumTimeOrder::bdf2});
  bdf2_state.seed_accepted_layers(bdf2_history, bdf2_committed);
  const auto bdf2_stencil = hundun::flow::make_momentum_time_stencil(
      hundun::flow::MomentumTimeOrder::bdf2, 0.02, 0.015);
  TestAccess::reset();
  const auto bdf2_report =
      flow.attempt(bdf2_state, rho_ref, 0.0, bdf2_stencil, {}, {});
  HUNDUN_CHECK(bdf2_report.disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  const double volume = geometry.cell_volume_m3(0U);
  for (std::size_t component = 0; component < 3U; ++component) {
    const double expected_diagonal =
        bdf2_stencil.alpha0 * rho_ref * volume / bdf2_stencil.dt_s;
    const double expected_rhs = -(rho_ref * volume / bdf2_stencil.dt_s) *
                                (bdf2_stencil.alpha1 * velocity_n[component] +
                                 bdf2_stencil.alpha2 * velocity_nm1[component]);
    if (mpi.rank() == 0) {
      std::cout << "TASK18_BDF2 component=" << component
                << " diagonal=" << TestAccess::last_momentum_diagonal(component)
                << " expected_diagonal=" << expected_diagonal
                << " rhs=" << TestAccess::last_momentum_rhs(component)
                << " expected_rhs=" << expected_rhs << '\n';
    }
    HUNDUN_CHECK_NEAR(TestAccess::last_momentum_diagonal(component),
                      expected_diagonal, 1.0e-13);
    HUNDUN_CHECK_NEAR(TestAccess::last_momentum_rhs(component), expected_rhs,
                      1.0e-13);
  }
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] { run_zero_flow_transaction(mpi); });
}

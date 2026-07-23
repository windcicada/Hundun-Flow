// SPDX-License-Identifier: Apache-2.0

#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/execution/execution.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/flow/adaptive_time_control.hpp"
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
#include "flow/src/adaptive_time_control_test_access.hpp"
#include "flow/src/constant_density_piso_test_access.hpp"
#include "flow/src/ideal_gas_closure_test_access.hpp"
#include "flow/src/material_density_piso_test_access.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/adaptive_time_control_test_support.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <string>

namespace {

hundun::runtime::Int3 grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw hundun::runtime::Error("unsupported Task22 rank count");
}

hundun::runtime::FieldDescriptor cell(const char *name,
                                      std::uint32_t components) {
  return {name, "1", "task22", hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64, components, 2, false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}
hundun::runtime::FieldDescriptor face(const char *name,
                                      std::uint32_t components) {
  return {name, "1", "task22", hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64, components, 0, false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::runtime::FieldDescriptor conservative_cell(
    const char *name, const char *unit) {
  return {name, unit, "task22",
          hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64, 1U, 2, true,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::config::FlowCaseConfig periodic_case() {
  hundun::config::FlowCaseConfig config{};
  config.schema_version = 2;
  config.simulation_type =
      hundun::config::SimulationType::variable_density_flow;
  config.density_model = hundun::config::DensityModel::constant;
  config.physics.rho_ref_kg_per_m3 = 1.0;
  config.physics.inlet_consistency_rtol = 1.0e-12;
  config.scalars.push_back({"alpha", 0.0});
  constexpr std::array<hundun::config::PatchName, 6> names{
      hundun::config::PatchName::x_min, hundun::config::PatchName::x_max,
      hundun::config::PatchName::y_min, hundun::config::PatchName::y_max,
      hundun::config::PatchName::z_min, hundun::config::PatchName::z_max};
  for (std::size_t i = 0; i < names.size(); ++i) {
    config.boundaries[i].patch = names[i];
    config.boundaries[i].type = hundun::config::BoundaryType::periodic;
  }
  return config;
}

class FirstCallFailureSolver final : public hundun::linear::LinearSolver {
public:
  explicit FirstCallFailureSolver(int lowest_failing_rank)
      : lowest_failing_rank_(lowest_failing_rank) {}

  hundun::linear::SolveReport
  solve(const hundun::linear::LinearOperator &,
        hundun::linear::Preconditioner &,
        hundun::execution::VectorView<const double>,
        hundun::execution::VectorView<double>,
        const hundun::linear::SolveControl &) const override {
    hundun::linear::SolveReport report;
    if (calls_++ == 0U) {
      report.reason =
          hundun::linear::SolveTerminationReason::maximum_iterations;
      report.lowest_failing_rank = lowest_failing_rank_;
      return report;
    }
    report.reason = hundun::linear::SolveTerminationReason::converged;
    report.final_residual = 0.0;
    report.lowest_failing_rank = -1;
    return report;
  }

private:
  int lowest_failing_rank_{};
  mutable std::size_t calls_{};
};

void run_fast(const hundun::runtime::MpiContext &mpi) {
  constexpr hundun::runtime::Int3 extent{8, 6, 4};
  auto decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, extent, {true, true, true},
      hundun::runtime::DecompositionOptions{grid(mpi.size())});
  const hundun::mesh::MeshTopology topology(decomposition);
  const hundun::mesh::MeshGeometry geometry(
      topology,
      hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  auto boundaries =
      hundun::boundary::BoundaryRegistry::create(periodic_case(), topology);

  hundun::runtime::FieldRegistry registry;
  hundun::flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell("rho", 1U));
  fields.velocity = registry.declare_field(cell("u", 3U));
  fields.mechanical_pressure = registry.declare_field(cell("pi", 1U));
  fields.face_velocity = registry.declare_field(face("uf", 3U));
  fields.face_mass_flux =
      hundun::finite_volume::declare_face_mass_flux(registry);
  fields.transported_cell_fields.push_back(
      registry.declare_field(cell("alpha", 1U)));
  registry.freeze();

  const auto box = decomposition.owned_box();
  const hundun::runtime::Int3 local{box.end.x - box.begin.x,
                                    box.end.y - box.begin.y,
                                    box.end.z - box.begin.z};
  auto state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0,
       hundun::flow::MomentumTimeOrder::backward_euler});
  const auto cells = topology.owned_cell_count();
  hundun::flow::FlowLayerValues initial;
  initial.density.assign(cells, 1.0);
  initial.velocity.assign(cells * 3U, 0.0);
  initial.mechanical_pressure.assign(cells, 0.0);
  initial.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  initial.face_mass_flux.assign(topology.local_face_count(), 0.0);
  initial.transported_cell_fields = {std::vector<double>(cells, 1.0)};
  state.seed_accepted_layers(initial, initial);
  auto other_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0,
       hundun::flow::MomentumTimeOrder::backward_euler});
  other_state.seed_accepted_layers(initial, initial);

  auto halo = hundun::runtime::HaloExchange::create(
      decomposition,
      hundun::runtime::ExchangePlan::create(decomposition, local, 2));
  hundun::execution::CpuReferenceContext execution;
  hundun::linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  hundun::linear::ConjugateGradientSolver pressure_solver(execution, mpi);
  hundun::linear::JacobiPreconditioner mx(execution), my(execution),
      mz(execution), pressure_preconditioner(execution);
  auto facade = hundun::flow::FixedStepConstantDensityFlow::create(
      decomposition, topology, geometry, boundaries, mpi, execution, halo,
      momentum_solver, {&mx, &my, &mz}, pressure_solver,
      pressure_preconditioner,
      {{fields.transported_cell_fields.front(),
        hundun::finite_volume::FiniteVolumeQuantity::scalar(0U), 0.0}});
  const hundun::config::FlowTimeConfig time{
      hundun::config::TimeMode::adaptive, 2, 0.01, 0.00125, 0.02,
      0.5, 0.25, 1.25, 0.5, 8};
  const int mismatch_rank = mpi.size() == 1 ? 0 : 1;
  const auto expect_factory_failure =
      [&](const hundun::config::FlowTimeConfig &candidate,
          hundun::config::DensityModel model, const char *prefix) {
        bool rejected = false;
        try {
          static_cast<void>(hundun::flow::Bdf2RetryController::create(
              candidate, model, topology, geometry, mpi, state));
        } catch (const hundun::runtime::Error &error) {
          rejected = true;
          HUNDUN_CHECK(std::string(error.what()).find(prefix) !=
                       std::string::npos);
          HUNDUN_CHECK(
              std::string(error.what()).find(std::to_string(mismatch_rank)) !=
              std::string::npos);
        }
        HUNDUN_CHECK(rejected);
        HUNDUN_CHECK(state.metadata().step == 0U);
      };
  {
    auto candidate = time;
    if (mpi.rank() == mismatch_rank)
      candidate.mode = static_cast<hundun::config::TimeMode>(255);
    expect_factory_failure(candidate, hundun::config::DensityModel::constant,
                           "time-control.create.invalid-config");
  }
  {
    auto candidate = time;
    if (mpi.rank() == mismatch_rank)
      candidate.steps = -1;
    expect_factory_failure(candidate, hundun::config::DensityModel::constant,
                           "time-control.create.invalid-config");
  }
  {
    auto candidate = time;
    if (mpi.rank() == mismatch_rank)
      candidate.initial_dt_s = std::numeric_limits<double>::infinity();
    expect_factory_failure(candidate, hundun::config::DensityModel::constant,
                           "time-control.create.invalid-config");
  }
  {
    auto candidate = time;
    if (mpi.rank() == mismatch_rank)
      candidate.cfl_target = 0.6;
    expect_factory_failure(candidate, hundun::config::DensityModel::constant,
                           "time-control.create.invalid-config");
  }
  {
    auto candidate = time;
    if (mpi.size() > 1) {
      if (mpi.rank() == mismatch_rank)
        ++candidate.steps;
      expect_factory_failure(candidate, hundun::config::DensityModel::constant,
                             "time-control.create.agreement");
    }
  }
  {
    auto model = hundun::config::DensityModel::constant;
    if (mpi.rank() == mismatch_rank)
      model = static_cast<hundun::config::DensityModel>(255);
    expect_factory_failure(time, model, "time-control.create.identity");
  }
  auto controller = hundun::flow::Bdf2RetryController::create(
      time, hundun::config::DensityModel::constant, topology, geometry, mpi,
      state);
  {
    auto wrong_model = hundun::flow::Bdf2RetryController::create(
        time, hundun::config::DensityModel::material, topology, geometry, mpi,
        state);
    const auto metadata_before = state.metadata();
    auto rejected = wrong_model.advance(state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(rejected.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
    HUNDUN_CHECK(rejected.attempt_count() == 0U);
    HUNDUN_CHECK(rejected.lowest_failing_rank() == 0);
    HUNDUN_CHECK(state.metadata().step == metadata_before.step);
  }
  {
    constexpr std::array<int, 4> control_fields{0, 1, 2, 3};
    if (mpi.size() > 1) {
      for (const int field : control_fields) {
        hundun::linear::SolveControl momentum{};
        if (mpi.rank() == mismatch_rank) {
          if (field == 0)
            momentum.atol *= 2.0;
          else if (field == 1)
            momentum.rtol *= 2.0;
          else if (field == 2)
            ++momentum.max_iterations;
          else
            ++momentum.residual_recompute_interval;
        }
        const auto rejected =
            controller.advance(state, facade, 1.0, 0.0, momentum, {});
        HUNDUN_CHECK(
            rejected.disposition() ==
            hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
        HUNDUN_CHECK(rejected.attempt_count() == 0U);
        HUNDUN_CHECK(rejected.lowest_failing_rank() == mismatch_rank);
        HUNDUN_CHECK(state.metadata().step == 0U);
      }
    }
  }
  {
    auto rejected =
        controller.advance(other_state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(rejected.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
    HUNDUN_CHECK(rejected.attempt_count() == 0U);
    HUNDUN_CHECK(rejected.lowest_failing_rank() == 0);
    HUNDUN_CHECK(other_state.metadata().step == 0U);
    HUNDUN_CHECK(state.metadata().step == 0U);
  }
  HUNDUN_CHECK(controller.state().accepted_step == 0U);
  {
    using TimeAccess = hundun::flow::test::AdaptiveTimeControlTestAccess;
    using Fault = hundun::flow::test::TimeControlPreflightFault;
    struct Case {
      Fault fault;
      std::uint8_t category;
    };
    constexpr std::array<Case, 3> cases{{
        {Fault::layout, 3U},
        {Fault::capability, 4U},
        {Fault::preparation, 7U},
    }};
    for (const auto &candidate : cases) {
      TimeAccess::reset_faults();
      TimeAccess::set_preflight_fault(candidate.fault, mismatch_rank);
      const auto metadata_before = state.metadata();
      auto rejected = controller.advance(state, facade, 1.0, 0.0, {}, {});
      HUNDUN_CHECK(
          rejected.disposition() ==
          hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
      HUNDUN_CHECK(rejected.attempt_count() == 0U);
      HUNDUN_CHECK(rejected.lowest_failing_rank() == mismatch_rank);
      HUNDUN_CHECK(TimeAccess::preflight_category(rejected) ==
                   candidate.category);
      HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
          metadata_before, state.metadata()));
      HUNDUN_CHECK(controller.state().accepted_step == 0U);
    }
    TimeAccess::reset_faults();
    const hundun::config::FlowTimeConfig fixed_time{
        hundun::config::TimeMode::fixed, 1, 0.01, 0.00001, 0.02,
        0.5, 0.25, 1.25, 0.5, 8};
    auto fixed_controller = hundun::flow::Bdf2RetryController::create(
        fixed_time, hundun::config::DensityModel::constant, topology,
        geometry, mpi, state);
    const auto reductions_before = mpi.fp64_reduction_counters();
    TimeAccess::set_preflight_fault(
        hundun::flow::test::TimeControlPreflightFault::preparation,
        mismatch_rank);
    const auto fixed_rejected =
        fixed_controller.advance(state, facade, 1.0, 0.0, {}, {});
    const auto reductions_after = mpi.fp64_reduction_counters();
    HUNDUN_CHECK(
        fixed_rejected.disposition() ==
        hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
    HUNDUN_CHECK(fixed_rejected.attempt_count() == 0U);
    HUNDUN_CHECK(TimeAccess::preflight_category(fixed_rejected) == 7U);
    HUNDUN_CHECK(fixed_rejected.lowest_failing_rank() == mismatch_rank);
    HUNDUN_CHECK(reductions_after.collective_calls ==
                 reductions_before.collective_calls);
    TimeAccess::reset_faults();
  }
  {
    using TimeAccess = hundun::flow::test::AdaptiveTimeControlTestAccess;
    auto malformed_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    malformed_state.seed_accepted_layers(initial, initial);
    auto malformed_controller = hundun::flow::Bdf2RetryController::create(
        time, hundun::config::DensityModel::constant, topology, geometry, mpi,
        malformed_state);
    const auto before =
        malformed_state.snapshot(hundun::flow::FlowLayer::committed);
    TimeAccess::set_recoverable_failures(1U, 0);
    TimeAccess::set_post_return_mutation(
        hundun::flow::test::TimeControlPostReturnMutation::attempted_dt,
        mismatch_rank);
    const auto malformed = malformed_controller.advance(
        malformed_state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(
        malformed.disposition() ==
        hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
    HUNDUN_CHECK(malformed.attempt_count() == 1U);
    HUNDUN_CHECK(TimeAccess::preflight_category(malformed) == 8U);
    HUNDUN_CHECK(malformed.lowest_failing_rank() == mismatch_rank);
    HUNDUN_CHECK(malformed_controller.state().accepted_step == 0U);
    HUNDUN_CHECK(malformed_state.metadata().step == 0U);
    HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
        before,
        malformed_state.snapshot(hundun::flow::FlowLayer::committed)));
    TimeAccess::reset_faults();
  }
  {
    using TimeAccess = hundun::flow::test::AdaptiveTimeControlTestAccess;
    constexpr std::array<hundun::flow::StepFailureReason, 6> reasons{
        hundun::flow::StepFailureReason::non_finite_trial,
        hundun::flow::StepFailureReason::final_momentum_residual,
        hundun::flow::StepFailureReason::final_transport_residual,
        hundun::flow::StepFailureReason::final_conservation_defect,
        hundun::flow::StepFailureReason::final_continuity_residual,
        hundun::flow::StepFailureReason::final_pressure_residual};
    for (const auto reason : reasons) {
      auto scheduled_state = hundun::flow::FlowState::create(
          registry, {local, topology.local_face_count()}, fields,
          {0U, 0.0, 0.01, 0.0,
           hundun::flow::MomentumTimeOrder::backward_euler});
      scheduled_state.seed_accepted_layers(initial, initial);
      auto scheduled_controller =
          hundun::flow::Bdf2RetryController::create(
              time, hundun::config::DensityModel::constant, topology,
              geometry, mpi, scheduled_state);
      TimeAccess::set_recoverable_failure_reason(reason);
      TimeAccess::set_recoverable_failures(1U, mismatch_rank);
      const auto scheduled = scheduled_controller.advance(
          scheduled_state, facade, 1.0, 0.0, {}, {});
      HUNDUN_CHECK(scheduled.disposition() ==
                   hundun::flow::TimeAdvanceDisposition::committed);
      HUNDUN_CHECK(scheduled.attempt_count() == 2U);
      HUNDUN_CHECK(scheduled.attempt(0).reason == reason);
      HUNDUN_CHECK(
          scheduled.attempt(0).disposition ==
          hundun::flow::StepAttemptDisposition::recoverable_failure);
      HUNDUN_CHECK(
          scheduled.attempt(1).disposition ==
          hundun::flow::StepAttemptDisposition::committed);
      HUNDUN_CHECK(scheduled_state.metadata().step == 1U);
      TimeAccess::reset_faults();
    }
  }
  {
    const auto run_linear_failure =
        [&](bool momentum_failure,
            hundun::flow::StepFailureReason expected_reason) {
          FirstCallFailureSolver one_shot(mismatch_rank);
          hundun::linear::JacobiPreconditioner fx(execution), fy(execution),
              fz(execution), fp(execution);
          const auto &selected_momentum =
              momentum_failure
                  ? static_cast<const hundun::linear::LinearSolver &>(one_shot)
                  : static_cast<const hundun::linear::LinearSolver &>(
                        momentum_solver);
          const auto &selected_pressure =
              momentum_failure
                  ? static_cast<const hundun::linear::LinearSolver &>(
                        pressure_solver)
                  : static_cast<const hundun::linear::LinearSolver &>(one_shot);
          auto failure_flow =
              hundun::flow::FixedStepConstantDensityFlow::create(
                  decomposition, topology, geometry, boundaries, mpi,
                  execution, halo, selected_momentum, {&fx, &fy, &fz},
                  selected_pressure, fp,
                  {{fields.transported_cell_fields.front(),
                    hundun::finite_volume::FiniteVolumeQuantity::scalar(0U),
                    0.0}});
          auto failure_state = hundun::flow::FlowState::create(
              registry, {local, topology.local_face_count()}, fields,
              {0U, 0.0, 0.01, 0.0,
               hundun::flow::MomentumTimeOrder::backward_euler});
          failure_state.seed_accepted_layers(initial, initial);
          auto failure_controller =
              hundun::flow::Bdf2RetryController::create(
                  time, hundun::config::DensityModel::constant, topology,
                  geometry, mpi, failure_state);
          const auto result = failure_controller.advance(
              failure_state, failure_flow, 1.0, 0.0, {}, {});
          HUNDUN_CHECK(result.disposition() ==
                       hundun::flow::TimeAdvanceDisposition::committed);
          HUNDUN_CHECK(result.attempt_count() == 2U);
          HUNDUN_CHECK(result.attempt(0).reason == expected_reason);
          HUNDUN_CHECK(
              result.attempt(0).disposition ==
              hundun::flow::StepAttemptDisposition::recoverable_failure);
          HUNDUN_CHECK(
              result.attempt(1).disposition ==
              hundun::flow::StepAttemptDisposition::committed);
        };
    run_linear_failure(
        true, hundun::flow::StepFailureReason::momentum_linear_solve);
    run_linear_failure(
        false, hundun::flow::StepFailureReason::pressure_linear_solve);
  }
  {
    auto signed_initial = initial;
    for (std::size_t face_id = 0;
         face_id < signed_initial.face_mass_flux.size(); ++face_id)
      signed_initial.face_mass_flux[face_id] =
          (face_id % 2U == 0U) ? 2.0 : -2.0;
    auto signed_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    signed_state.seed_accepted_layers(signed_initial, signed_initial);
    std::vector<double> flux_sum(topology.owned_cell_count());
    std::vector<double> geometry_sum(topology.owned_cell_count());
    for (hundun::mesh::LocalFaceId face_id = 0;
         face_id < topology.local_face_count(); ++face_id) {
      const auto area = geometry.face_area_m2(face_id);
      const auto d = geometry.face_displacement_m(face_id);
      const auto s = geometry.face_area_vector_m2(
          face_id, hundun::mesh::FaceSide::owner);
      const double factor =
          area * area / std::abs(s.x * d.x + s.y * d.y + s.z * d.z);
      const auto add = [&](std::size_t cell_id) {
        flux_sum[cell_id] +=
            std::abs(signed_initial.face_mass_flux[face_id]);
        geometry_sum[cell_id] += factor;
      };
      const auto owner = topology.owner(face_id);
      if (owner < topology.owned_cell_count())
        add(owner);
      const auto neighbour = topology.neighbour(face_id);
      if (neighbour && *neighbour < topology.owned_cell_count())
        add(*neighbour);
    }
    double expected[2]{};
    for (std::size_t cell_id = 0; cell_id < flux_sum.size(); ++cell_id) {
      const double denominator = 2.0 * geometry.cell_volume_m3(cell_id);
      expected[0] = std::max(expected[0], flux_sum[cell_id] / denominator);
      expected[1] =
          std::max(expected[1], 0.1 * geometry_sum[cell_id] / denominator);
    }
    HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, expected, 2, MPI_DOUBLE, MPI_MAX,
                               mpi.comm()) == MPI_SUCCESS);
    auto signed_controller = hundun::flow::Bdf2RetryController::create(
        time, hundun::config::DensityModel::constant, topology, geometry, mpi,
        signed_state);
    const auto signed_report =
        signed_controller.advance(signed_state, facade, 1.0, 0.1, {}, {});
    HUNDUN_CHECK(signed_report.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::committed);
    HUNDUN_CHECK(signed_report.convective_rate_per_s() == expected[0]);
    HUNDUN_CHECK(signed_report.diffusive_rate_per_s() == expected[1]);

    auto invalid_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    invalid_state.seed_accepted_layers(initial, initial);
    hundun::flow::test::MaterialDensityPisoTestAccess::
        set_accepted_face_mass_flux(
            invalid_state, 0U,
            std::numeric_limits<double>::quiet_NaN());
    auto invalid_controller = hundun::flow::Bdf2RetryController::create(
        time, hundun::config::DensityModel::constant, topology, geometry, mpi,
        invalid_state);
    const auto invalid_report =
        invalid_controller.advance(invalid_state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(invalid_report.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
    HUNDUN_CHECK(invalid_report.attempt_count() == 0U);
    HUNDUN_CHECK(
        hundun::flow::test::AdaptiveTimeControlTestAccess::preflight_category(
            invalid_report) == 5U);
  }
  {
    using TimeAccess = hundun::flow::test::AdaptiveTimeControlTestAccess;
    using Mutation = hundun::flow::test::TimeControlPostReturnMutation;
    constexpr std::array<Mutation, 5> mutations{
        Mutation::momentum_x_iterations, Mutation::momentum_y_iterations,
        Mutation::momentum_z_iterations, Mutation::pressure_one_iterations,
        Mutation::pressure_two_iterations};
    for (const auto mutation : mutations) {
      auto gate_state = hundun::flow::FlowState::create(
          registry, {local, topology.local_face_count()}, fields,
          {0U, 0.0, 0.01, 0.0,
           hundun::flow::MomentumTimeOrder::backward_euler});
      gate_state.seed_accepted_layers(initial, initial);
      auto gate_controller = hundun::flow::Bdf2RetryController::create(
          time, hundun::config::DensityModel::constant, topology, geometry,
          mpi, gate_state);
      TimeAccess::set_post_return_mutation(mutation, mpi.rank());
      const auto report =
          gate_controller.advance(gate_state, facade, 1.0, 0.0, {}, {});
      HUNDUN_CHECK(report.disposition() ==
                   hundun::flow::TimeAdvanceDisposition::committed);
      HUNDUN_CHECK(
          !report.attempt(0).all_linear_solves_within_half_limit);
      HUNDUN_CHECK(
          !gate_controller.state().last_all_linear_solves_within_half_limit);
      HUNDUN_CHECK(gate_controller.state().proposed_next_dt_s == 0.01);
      TimeAccess::reset_faults();
    }
  }
  {
    auto retry_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    retry_state.seed_accepted_layers(initial, initial);
    auto retry_controller = hundun::flow::Bdf2RetryController::create(
        time, hundun::config::DensityModel::constant, topology, geometry, mpi,
        retry_state);
    using TimeAccess = hundun::flow::test::AdaptiveTimeControlTestAccess;
    TimeAccess::reset_faults();
    TimeAccess::set_recoverable_failures(1U, 0);
    const auto retried =
        retry_controller.advance(retry_state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(retried.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::committed);
    HUNDUN_CHECK(retried.attempt_count() == 2U);
    HUNDUN_CHECK(retried.attempt(0).attempted_dt_s == 0.01);
    HUNDUN_CHECK(retried.attempt(1).attempted_dt_s == 0.005);
    HUNDUN_CHECK(retry_controller.state().last_retry_count == 1U);
    HUNDUN_CHECK(retry_state.metadata().step == 1U);
    TimeAccess::reset_faults();

    const hundun::config::FlowTimeConfig limit_time{
        hundun::config::TimeMode::fixed, 1, 0.01, 0.00001, 0.02,
        0.5, 0.25, 1.25, 0.5, 8};
    auto limit_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    limit_state.seed_accepted_layers(initial, initial);
    auto limit_controller = hundun::flow::Bdf2RetryController::create(
        limit_time, hundun::config::DensityModel::constant, topology,
        geometry, mpi, limit_state);
    const auto workspace_before =
        hundun::flow::test::ConstantDensityPisoTestAccess::
            mesh_workspace_snapshot(facade);
    const auto pressure_before =
        hundun::flow::test::ConstantDensityPisoTestAccess::
            pressure_operator_snapshot(facade);
    const auto reductions_before = mpi.fp64_reduction_counters();
    const auto exact_before = hundun::test::capture_adaptive_flow_state(
        limit_state, limit_controller.state(), {},
        workspace_before.data_identity, workspace_before.total_capacity,
        pressure_before.revision,
        {reductions_before.collective_calls,
         reductions_before.reduced_scalars,
         reductions_before.logical_payload_bytes, 0U});
    TimeAccess::set_recoverable_failures(9U, 0);
    const auto limited =
        limit_controller.advance(limit_state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(limited.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::retry_limit_reached);
    HUNDUN_CHECK(limited.attempt_count() == 9U);
    for (std::size_t index = 0; index < limited.attempt_count(); ++index)
      HUNDUN_CHECK(limited.attempt(index).disposition ==
                   hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(limit_state.metadata().step == 0U);
    const auto workspace_after =
        hundun::flow::test::ConstantDensityPisoTestAccess::
            mesh_workspace_snapshot(facade);
    const auto pressure_after =
        hundun::flow::test::ConstantDensityPisoTestAccess::
            pressure_operator_snapshot(facade);
    const auto reductions_after = mpi.fp64_reduction_counters();
    const auto exact_after = hundun::test::capture_adaptive_flow_state(
        limit_state, limit_controller.state(), {},
        workspace_after.data_identity, workspace_after.total_capacity,
        pressure_after.revision,
        {reductions_after.collective_calls,
         reductions_after.reduced_scalars,
         reductions_after.logical_payload_bytes, 0U});
    HUNDUN_CHECK(hundun::test::adaptive_flow_state_failed_attempt_preserved(
        exact_before, exact_after));
    TimeAccess::reset_faults();
    const auto fixed_success =
        limit_controller.advance(limit_state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(fixed_success.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::committed);
    HUNDUN_CHECK(fixed_success.attempt_count() == 1U);
    HUNDUN_CHECK(fixed_success.attempt(0).attempted_dt_s == 0.01);
    HUNDUN_CHECK(limit_controller.state().proposed_next_dt_s == 0.01);

    const hundun::config::FlowTimeConfig minimum_time{
        hundun::config::TimeMode::adaptive, 1, 0.01, 0.005, 0.02,
        0.5, 0.25, 1.25, 0.5, 8};
    auto minimum_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    minimum_state.seed_accepted_layers(initial, initial);
    auto minimum_controller = hundun::flow::Bdf2RetryController::create(
        minimum_time, hundun::config::DensityModel::constant, topology,
        geometry, mpi, minimum_state);
    TimeAccess::set_recoverable_failures(9U, 0);
    const auto minimum = minimum_controller.advance(
        minimum_state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(minimum.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::minimum_dt_failure);
    HUNDUN_CHECK(minimum.attempt_count() == 2U);
    HUNDUN_CHECK(minimum.attempt(1).attempted_dt_s == 0.005);
    HUNDUN_CHECK(minimum.limited_by_min_dt());
    TimeAccess::reset_faults();

    auto terminal_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    terminal_state.seed_accepted_layers(initial, initial);
    const auto maximum = std::numeric_limits<std::uint64_t>::max();
    hundun::flow::test::MaterialDensityPisoTestAccess::force_state_metadata(
        terminal_state,
        {maximum, 1.0, 0.01, 0.01,
         hundun::flow::MomentumTimeOrder::bdf2});
    hundun::flow::TimeControlState terminal_control;
    terminal_control.accepted_step = maximum;
    terminal_control.proposed_next_dt_s = 0.01;
    terminal_control.last_accepted_dt_s = 0.01;
    terminal_control.last_accepted_order =
        hundun::flow::MomentumTimeOrder::bdf2;
    terminal_control.history_ready = true;
    terminal_control.last_stability_metrics_available = true;
    terminal_control.revision = maximum;
    terminal_control.state_seal = TimeAccess::seal(
        time, hundun::config::DensityModel::constant, terminal_control);
    auto terminal_controller = hundun::flow::Bdf2RetryController::restore(
        time, hundun::config::DensityModel::constant, topology, geometry, mpi,
        terminal_state, terminal_control);
    const auto overflow = terminal_controller.advance(
        terminal_state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(overflow.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
    HUNDUN_CHECK(overflow.attempt_count() == 0U);
    HUNDUN_CHECK(overflow.lowest_failing_rank() == 0);
    HUNDUN_CHECK(TimeAccess::preflight_category(overflow) == 5U);
    HUNDUN_CHECK(terminal_controller.state().accepted_step == maximum);
  }
  const auto before = state.snapshot(hundun::flow::FlowLayer::committed);
  using TimeAccess = hundun::flow::test::AdaptiveTimeControlTestAccess;
  using PostMutation = hundun::flow::test::TimeControlPostReturnMutation;
  TimeAccess::reset_faults();
  TimeAccess::set_post_return_mutation(PostMutation::attempted_dt, mpi.rank());
  auto report = controller.advance(state, facade, 1.0, 0.0, {}, {});
  HUNDUN_CHECK(report.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  HUNDUN_CHECK(report.attempt_count() == 1U);
  HUNDUN_CHECK(report.attempt(0).order ==
               hundun::flow::MomentumTimeOrder::backward_euler);
  HUNDUN_CHECK(report.accepted_dt_s() == 0.01);
  HUNDUN_CHECK(report.stability_metrics_available());
  HUNDUN_CHECK(report.convective_rate_per_s() == 0.0);
  HUNDUN_CHECK(report.diffusive_rate_per_s() == 0.0);
  HUNDUN_CHECK(controller.state().accepted_step == 1U);
  HUNDUN_CHECK(state.metadata().step == 1U);
  HUNDUN_CHECK(TimeAccess::post_commit_observation_count() == 0U);
  TimeAccess::reset_faults();
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      before, state.snapshot(hundun::flow::FlowLayer::committed)));

  const auto restore_workspace =
      hundun::flow::test::ConstantDensityPisoTestAccess::
          mesh_workspace_snapshot(facade);
  const auto restore_pressure =
      hundun::flow::test::ConstantDensityPisoTestAccess::
          pressure_operator_snapshot(facade);
  const auto restore_before = hundun::test::capture_adaptive_flow_state(
      state, controller.state(), {}, restore_workspace.data_identity,
      restore_workspace.total_capacity, restore_pressure.revision);
  auto restored = hundun::flow::Bdf2RetryController::restore(
      time, hundun::config::DensityModel::constant, topology, geometry, mpi,
      state, controller.state());
  const auto restore_after = hundun::test::capture_adaptive_flow_state(
      state, restored.state(), {}, restore_workspace.data_identity,
      restore_workspace.total_capacity, restore_pressure.revision);
  HUNDUN_CHECK(hundun::test::adaptive_flow_state_bitwise_equal(
      restore_before, restore_after));
  HUNDUN_CHECK(restored.state().state_seal == controller.state().state_seal);
  auto second = restored.advance(state, facade, 1.0, 0.0, {}, {});
  HUNDUN_CHECK(second.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  HUNDUN_CHECK(second.attempt(0).order ==
               hundun::flow::MomentumTimeOrder::bdf2);
  HUNDUN_CHECK(state.metadata().step == 2U);
  TimeAccess::set_recoverable_failure_reason(
      hundun::flow::StepFailureReason::non_finite_trial);
  TimeAccess::set_recoverable_failures(2U, 0);
  const auto fallback =
      restored.advance(state, facade, 1.0, 0.0, {}, {});
  HUNDUN_CHECK(fallback.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  HUNDUN_CHECK(fallback.attempt_count() == 3U);
  HUNDUN_CHECK(fallback.attempt(0).order ==
               hundun::flow::MomentumTimeOrder::bdf2);
  HUNDUN_CHECK(fallback.attempt(1).order ==
               hundun::flow::MomentumTimeOrder::bdf2);
  HUNDUN_CHECK(fallback.attempt(2).order ==
               hundun::flow::MomentumTimeOrder::backward_euler);
  HUNDUN_CHECK(state.metadata().step == 3U);
  TimeAccess::reset_faults();
  const auto recovered =
      restored.advance(state, facade, 1.0, 0.0, {}, {});
  HUNDUN_CHECK(recovered.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  HUNDUN_CHECK(recovered.attempt(0).order ==
               hundun::flow::MomentumTimeOrder::bdf2);
  HUNDUN_CHECK(state.metadata().step == 4U);

  auto other_controller = hundun::flow::Bdf2RetryController::restore(
      time, hundun::config::DensityModel::constant, topology, geometry, mpi,
      state, restored.state());
  bool cross_controller_rejected = false;
  try {
    static_cast<void>(other_controller.diagnostic_source(state, second));
  } catch (const hundun::runtime::Error &) {
    cross_controller_rejected = true;
  }
  HUNDUN_CHECK(cross_controller_rejected);

  const auto out_of_band_stencil =
      hundun::flow::make_momentum_time_stencil(
          hundun::flow::MomentumTimeOrder::bdf2, state.metadata().dt_s,
          state.metadata().dt_s);
  const auto out_of_band =
      facade.attempt(state, 1.0, 0.0, out_of_band_stencil, {}, {});
  HUNDUN_CHECK(out_of_band.disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  const auto controller_before_rejection = restored.state();
  const auto rejected_after_out_of_band =
      restored.advance(state, facade, 1.0, 0.0, {}, {});
  HUNDUN_CHECK(rejected_after_out_of_band.disposition() ==
               hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
  HUNDUN_CHECK(rejected_after_out_of_band.attempt_count() == 0U);
  HUNDUN_CHECK(restored.state().state_seal ==
               controller_before_rejection.state_seal);
}

void run_acceptance(const hundun::runtime::MpiContext &mpi) {
  run_fast(mpi);
  constexpr hundun::runtime::Int3 extent{8, 4, 4};
  auto decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, extent, {true, true, true},
      hundun::runtime::DecompositionOptions{grid(mpi.size())});
  const hundun::mesh::MeshTopology topology(decomposition);
  const hundun::mesh::MeshGeometry geometry(
      topology,
      hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  auto config = periodic_case();
  config.scalars.clear();
  auto boundaries =
      hundun::boundary::BoundaryRegistry::create(config, topology);

  hundun::runtime::FieldRegistry registry;
  hundun::flow::FlowFieldIds fields;
  fields.density =
      registry.declare_field(conservative_cell("rho", "kg/m3"));
  fields.velocity = registry.declare_field(cell("u", 3U));
  fields.mechanical_pressure = registry.declare_field(cell("pi", 1U));
  fields.face_velocity = registry.declare_field(face("uf", 3U));
  fields.face_mass_flux =
      hundun::finite_volume::declare_face_mass_flux(registry);
  const auto rho_h =
      registry.declare_field(conservative_cell("rho_h", "J/m3"));
  fields.transported_cell_fields = {rho_h};
  registry.freeze();

  const auto box = decomposition.owned_box();
  const hundun::runtime::Int3 local{box.end.x - box.begin.x,
                                    box.end.y - box.begin.y,
                                    box.end.z - box.begin.z};
  const auto make_state = [&] {
    auto state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    hundun::flow::FlowLayerValues initial;
    initial.density.assign(topology.owned_cell_count(), 1.0);
    initial.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
    initial.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
    initial.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
    initial.face_mass_flux.assign(topology.local_face_count(), 0.0);
    initial.transported_cell_fields = {
        std::vector<double>(topology.owned_cell_count(), 300000.0)};
    state.seed_accepted_layers(initial, initial);
    return state;
  };

  auto halo = hundun::runtime::HaloExchange::create(
      decomposition,
      hundun::runtime::ExchangePlan::create(decomposition, local, 2));
  hundun::execution::CpuReferenceContext execution;
  hundun::linear::ConjugateGradientSolver momentum_solver(execution, mpi);
  hundun::linear::ConjugateGradientSolver pressure_solver(execution, mpi);
  hundun::linear::JacobiPreconditioner mx(execution), my(execution),
      mz(execution), pressure_preconditioner(execution);
  hundun::flow::MaterialDensityTransportSpec material_spec;
  material_spec.enthalpy_density = rho_h;
  material_spec.enthalpy_diffusivity_kg_per_m_s = 0.0;
  auto material_flow = hundun::flow::FixedStepMaterialDensityFlow::create(
      decomposition, topology, geometry, boundaries, mpi, execution, halo,
      momentum_solver, {&mx, &my, &mz}, pressure_solver,
      pressure_preconditioner, registry, fields, material_spec);
  const hundun::config::FlowTimeConfig time{
      hundun::config::TimeMode::adaptive, 2, 0.01, 0.00125, 0.02,
      0.5, 0.25, 1.25, 0.5, 8};
  using TimeAccess = hundun::flow::test::AdaptiveTimeControlTestAccess;
  using PostMutation = hundun::flow::test::TimeControlPostReturnMutation;
  const int failure_rank = mpi.size() == 1 ? 0 : 1;
  {
    auto state = make_state();
    const hundun::config::FlowTimeConfig terminal_time{
        hundun::config::TimeMode::fixed, 1, 0.01, 0.00001, 0.02,
        0.5, 0.25, 1.25, 0.5, 8};
    auto terminal_controller = hundun::flow::Bdf2RetryController::create(
        terminal_time, hundun::config::DensityModel::material, topology,
        geometry, mpi, state);
    const auto terminal_before = hundun::test::capture_adaptive_flow_state(
        state, terminal_controller.state());
    TimeAccess::set_recoverable_failures(9U, 0);
    const auto terminal = terminal_controller.advance(
        state, material_flow, 0.0, {}, {});
    HUNDUN_CHECK(terminal.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::retry_limit_reached);
    HUNDUN_CHECK(terminal.attempt_count() == 9U);
    const auto terminal_after = hundun::test::capture_adaptive_flow_state(
        state, terminal_controller.state());
    HUNDUN_CHECK(hundun::test::adaptive_flow_state_failed_attempt_preserved(
        terminal_before, terminal_after));
    TimeAccess::reset_faults();

    auto controller = hundun::flow::Bdf2RetryController::create(
        time, hundun::config::DensityModel::material, topology, geometry, mpi,
        state);
    TimeAccess::reset_faults();
    TimeAccess::set_post_return_mutation(PostMutation::attempted_dt,
                                         mpi.rank());
    auto report = controller.advance(state, material_flow, 0.0, {}, {});
    HUNDUN_CHECK(report.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::committed);
    HUNDUN_CHECK(controller.state().accepted_step == 1U);
    HUNDUN_CHECK(state.metadata().step == 1U);
    HUNDUN_CHECK(TimeAccess::post_commit_observation_count() == 0U);
    HUNDUN_CHECK(std::holds_alternative<
                 hundun::flow::MaterialDensityStepAttemptReport>(
        report.final_attempt()));
    const auto &material_final =
        std::get<hundun::flow::MaterialDensityStepAttemptReport>(
            report.final_attempt());
    HUNDUN_CHECK(material_final.flux_provenance() ==
                 hundun::flow::MaterialFluxProvenance::final_corrected);
    HUNDUN_CHECK(material_final.material_field_count() == 1U);
    TimeAccess::reset_faults();
  }
  {
    auto state = make_state();
    auto controller = hundun::flow::Bdf2RetryController::create(
        time, hundun::config::DensityModel::material, topology, geometry, mpi,
        state);
    TimeAccess::set_recoverable_failure_reason(
        hundun::flow::StepFailureReason::final_transport_residual);
    TimeAccess::set_recoverable_failures(1U, 0);
    const auto report = controller.advance(state, material_flow, 0.0, {}, {});
    HUNDUN_CHECK(report.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::committed);
    HUNDUN_CHECK(report.attempt_count() == 2U);
    HUNDUN_CHECK(
        report.attempt(0).reason ==
        hundun::flow::StepFailureReason::final_transport_residual);
    HUNDUN_CHECK(controller.state().last_retry_count == 1U);
    TimeAccess::reset_faults();
  }

  {
    auto terminal_state = make_state();
    auto terminal_closure = hundun::flow::IdealGasClosure::create(
        topology, geometry, boundaries, mpi, registry, fields, terminal_state,
        {rho_h, 1000.0, 287.0, 86100.0});
    auto terminal_flow = hundun::flow::FixedStepIdealGasFlow::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        momentum_solver, {&mx, &my, &mz}, pressure_solver,
        pressure_preconditioner, registry, fields, material_spec,
        std::move(terminal_closure));
    const hundun::config::FlowTimeConfig terminal_time{
        hundun::config::TimeMode::fixed, 1, 0.01, 0.00001, 0.02,
        0.5, 0.25, 1.25, 0.5, 8};
    auto terminal_controller = hundun::flow::Bdf2RetryController::create(
        terminal_time, hundun::config::DensityModel::ideal_gas, topology,
        geometry, mpi, terminal_state);
    const auto terminal_before = hundun::test::capture_adaptive_flow_state(
        terminal_state, terminal_controller.state(),
        terminal_flow.closure_state());
    TimeAccess::set_recoverable_failures(9U, 0);
    const auto terminal = terminal_controller.advance(
        terminal_state, terminal_flow, 0.0, {}, {});
    HUNDUN_CHECK(terminal.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::retry_limit_reached);
    HUNDUN_CHECK(terminal.attempt_count() == 9U);
    const auto terminal_after = hundun::test::capture_adaptive_flow_state(
        terminal_state, terminal_controller.state(),
        terminal_flow.closure_state());
    HUNDUN_CHECK(hundun::test::adaptive_flow_state_failed_attempt_preserved(
        terminal_before, terminal_after));
    TimeAccess::reset_faults();
  }

  auto ideal_state = make_state();
  auto closure = hundun::flow::IdealGasClosure::create(
      topology, geometry, boundaries, mpi, registry, fields, ideal_state,
      {rho_h, 1000.0, 287.0, 86100.0});
  auto ideal_flow = hundun::flow::FixedStepIdealGasFlow::create(
      decomposition, topology, geometry, boundaries, mpi, execution, halo,
      momentum_solver, {&mx, &my, &mz}, pressure_solver,
      pressure_preconditioner, registry, fields, material_spec,
      std::move(closure));
  {
    auto controller = hundun::flow::Bdf2RetryController::create(
        time, hundun::config::DensityModel::ideal_gas, topology, geometry, mpi,
        ideal_state);
    TimeAccess::set_post_return_mutation(PostMutation::reason, mpi.rank());
    auto report = controller.advance(ideal_state, ideal_flow, 0.0, {}, {});
    HUNDUN_CHECK(report.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::committed);
    HUNDUN_CHECK(controller.state().accepted_step == 1U);
    HUNDUN_CHECK(ideal_state.metadata().step == 1U);
    HUNDUN_CHECK(TimeAccess::post_commit_observation_count() == 0U);
    HUNDUN_CHECK(std::holds_alternative<
                 hundun::flow::IdealGasStepAttemptReport>(
        report.final_attempt()));
    const auto &ideal_final =
        std::get<hundun::flow::IdealGasStepAttemptReport>(
            report.final_attempt());
    HUNDUN_CHECK(ideal_final.flow().flux_provenance() ==
                 hundun::flow::MaterialFluxProvenance::final_corrected);
    HUNDUN_CHECK(ideal_final.closure_report().eos_max_relative_error() <=
                 1.0e-12);
    const auto closed_state = ideal_flow.closure_state();
    HUNDUN_CHECK(closed_state.mode ==
                 hundun::flow::IdealGasPressureMode::closed_dynamic);
    HUNDUN_CHECK(closed_state.target_mass_kg.has_value());
    TimeAccess::reset_faults();
    auto retry_controller = hundun::flow::Bdf2RetryController::restore(
        time, hundun::config::DensityModel::ideal_gas, topology, geometry, mpi,
        ideal_state, controller.state());
    TimeAccess::set_recoverable_failure_reason(
        hundun::flow::StepFailureReason::density_closure_failure);
    TimeAccess::set_recoverable_failures(1U, 0);
    const auto retried =
        retry_controller.advance(ideal_state, ideal_flow, 0.0, {}, {});
    HUNDUN_CHECK(retried.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::committed);
    HUNDUN_CHECK(retried.attempt_count() == 2U);
    HUNDUN_CHECK(
        retried.attempt(0).reason ==
        hundun::flow::StepFailureReason::density_closure_failure);
    HUNDUN_CHECK(retry_controller.state().last_retry_count == 1U);
    HUNDUN_CHECK(ideal_state.metadata().step == 2U);
    TimeAccess::reset_faults();
    for (const auto reason :
         {hundun::flow::StepFailureReason::transport_failure,
          hundun::flow::StepFailureReason::final_conservation_defect}) {
      TimeAccess::set_recoverable_failure_reason(reason);
      TimeAccess::set_recoverable_failures(1U, failure_rank);
      const auto scheduled = retry_controller.advance(
          ideal_state, ideal_flow, 0.0, {}, {});
      HUNDUN_CHECK(scheduled.disposition() ==
                   hundun::flow::TimeAdvanceDisposition::committed);
      HUNDUN_CHECK(scheduled.attempt_count() == 2U);
      HUNDUN_CHECK(scheduled.attempt(0).reason == reason);
      HUNDUN_CHECK(
          scheduled.attempt(1).disposition ==
          hundun::flow::StepAttemptDisposition::committed);
      TimeAccess::reset_faults();
    }
  }

  {
    constexpr hundun::runtime::Int3 open_extent{8, 4, 4};
    auto open_decomposition =
        hundun::runtime::StructuredDecomposition::create(
            mpi, open_extent, {false, false, false},
            hundun::runtime::DecompositionOptions{grid(mpi.size())});
    const hundun::mesh::MeshTopology open_topology(open_decomposition);
    const hundun::mesh::MeshGeometry open_geometry(
        open_topology,
        hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0},
                                        {1.0, 0.25, 0.25}));
    auto open_config = periodic_case();
    open_config.density_model = hundun::config::DensityModel::ideal_gas;
    open_config.scalars.clear();
    open_config.physics.cp_J_per_kg_K = 1000.0;
    open_config.physics.gas_constant_J_per_kg_K = 287.05;
    open_config.physics.thermodynamic_pressure_pa = 101325.0;
    open_config.boundaries[0].type =
        hundun::config::BoundaryType::velocity_inlet;
    open_config.boundaries[0].velocity_m_per_s = {1.0, 0.0, 0.0};
    open_config.boundaries[0].thermal_authority =
        hundun::config::InletThermalAuthority::temperature;
    open_config.boundaries[0].temperature_K = 300.0;
    open_config.boundaries[0].enthalpy_J_per_kg = 300000.0;
    const double open_density = 101325.0 / (287.05 * 300.0);
    open_config.boundaries[0].density_kg_per_m3 = open_density;
    open_config.boundaries[0].scalar_values =
        std::vector<hundun::config::InletScalarValue>{};
    open_config.boundaries[1].type =
        hundun::config::BoundaryType::pressure_outlet;
    open_config.boundaries[1].pressure_perturbation_pa = 0.0;
    for (std::size_t patch = 2U; patch < 6U; ++patch)
      open_config.boundaries[patch].type =
          hundun::config::BoundaryType::symmetry;
    auto open_boundaries = hundun::boundary::BoundaryRegistry::create(
        open_config, open_topology);

    hundun::runtime::FieldRegistry open_registry;
    hundun::flow::FlowFieldIds open_fields;
    open_fields.density =
        open_registry.declare_field(conservative_cell("rho", "kg/m3"));
    open_fields.velocity = open_registry.declare_field(cell("u", 3U));
    open_fields.mechanical_pressure =
        open_registry.declare_field(cell("pi", 1U));
    open_fields.face_velocity =
        open_registry.declare_field(face("uf", 3U));
    open_fields.face_mass_flux =
        hundun::finite_volume::declare_face_mass_flux(open_registry);
    const auto open_rho_h =
        open_registry.declare_field(conservative_cell("rho_h", "J/m3"));
    open_fields.transported_cell_fields = {open_rho_h};
    open_registry.freeze();
    const auto open_box = open_decomposition.owned_box();
    const hundun::runtime::Int3 open_local{
        open_box.end.x - open_box.begin.x,
        open_box.end.y - open_box.begin.y,
        open_box.end.z - open_box.begin.z};
    auto open_state = hundun::flow::FlowState::create(
        open_registry, {open_local, open_topology.local_face_count()},
        open_fields,
        {0U, 0.0, 0.001, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    hundun::flow::FlowLayerValues open_initial;
    open_initial.density.assign(open_topology.owned_cell_count(),
                                open_density);
    open_initial.velocity.assign(open_topology.owned_cell_count() * 3U, 0.0);
    for (std::size_t cell_id = 0;
         cell_id < open_topology.owned_cell_count(); ++cell_id)
      open_initial.velocity[cell_id * 3U] = 1.0;
    open_initial.mechanical_pressure.assign(
        open_topology.owned_cell_count(), 0.0);
    open_initial.face_velocity.assign(
        open_topology.local_face_count() * 3U, 0.0);
    open_initial.face_mass_flux.assign(open_topology.local_face_count(), 0.0);
    for (hundun::mesh::LocalFaceId face_id = 0;
         face_id < open_topology.local_face_count(); ++face_id) {
      const auto area = open_geometry.face_area_vector_m2(
          face_id, hundun::mesh::FaceSide::owner);
      open_initial.face_velocity[face_id * 3U] = 1.0;
      open_initial.face_mass_flux[face_id] = open_density * area.x;
    }
    open_initial.transported_cell_fields = {std::vector<double>(
        open_topology.owned_cell_count(), open_density * 300000.0)};
    open_state.seed_accepted_layers(open_initial, open_initial);
    auto open_halo = hundun::runtime::HaloExchange::create(
        open_decomposition,
        hundun::runtime::ExchangePlan::create(open_decomposition, open_local,
                                              2));
    hundun::flow::MaterialDensityTransportSpec open_material_spec;
    open_material_spec.enthalpy_density = open_rho_h;
    auto open_closure = hundun::flow::IdealGasClosure::create(
        open_topology, open_geometry, open_boundaries, mpi, open_registry,
        open_fields, open_state,
        {open_rho_h, 1000.0, 287.05, 101325.0});
    auto open_flow = hundun::flow::FixedStepIdealGasFlow::create(
        open_decomposition, open_topology, open_geometry, open_boundaries, mpi,
        execution, open_halo, momentum_solver, {&mx, &my, &mz},
        pressure_solver, pressure_preconditioner, open_registry, open_fields,
        open_material_spec, std::move(open_closure));
    const hundun::config::FlowTimeConfig open_time{
        hundun::config::TimeMode::adaptive, 1, 0.001, 0.000125, 0.002,
        0.5, 0.25, 1.25, 0.5, 8};
    auto open_controller = hundun::flow::Bdf2RetryController::create(
        open_time, hundun::config::DensityModel::ideal_gas, open_topology,
        open_geometry, mpi, open_state);
    TimeAccess::set_recoverable_failure_reason(
        hundun::flow::StepFailureReason::boundary_backflow);
    TimeAccess::set_recoverable_failures(1U, 0);
    const auto open_report =
        open_controller.advance(open_state, open_flow, 0.0, {}, {});
    HUNDUN_CHECK(open_report.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::committed);
    HUNDUN_CHECK(open_report.attempt_count() == 2U);
    HUNDUN_CHECK(
        open_report.attempt(0).reason ==
        hundun::flow::StepFailureReason::boundary_backflow);
    HUNDUN_CHECK(open_report.attempt(0).attempted_dt_s == 0.001);
    HUNDUN_CHECK(open_report.attempt(1).attempted_dt_s == 0.0005);
    HUNDUN_CHECK(open_report.attempt(0).order ==
                 hundun::flow::MomentumTimeOrder::backward_euler);
    HUNDUN_CHECK(open_report.attempt(1).order ==
                 hundun::flow::MomentumTimeOrder::backward_euler);
    const auto &open_final =
        std::get<hundun::flow::IdealGasStepAttemptReport>(
            open_report.final_attempt());
    HUNDUN_CHECK(open_final.flow().flux_provenance() ==
                 hundun::flow::MaterialFluxProvenance::final_corrected);
    HUNDUN_CHECK(open_final.closure_report().eos_max_relative_error() <=
                 1.0e-12);
    const auto open_closure_state = open_flow.closure_state();
    HUNDUN_CHECK(open_closure_state.mode ==
                 hundun::flow::IdealGasPressureMode::open_fixed);
    HUNDUN_CHECK(open_closure_state.thermodynamic_pressure_pa == 101325.0);
    HUNDUN_CHECK(!open_closure_state.target_mass_kg.has_value());
    TimeAccess::reset_faults();
  }
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  if (argc != 2)
    return 2;
  const std::string mode(argv[1]);
  if (mode == "fast")
    return hundun::test::run([&] { run_fast(mpi); });
  if (mode == "acceptance")
    return hundun::test::run([&] { run_acceptance(mpi); });
  return 2;
}

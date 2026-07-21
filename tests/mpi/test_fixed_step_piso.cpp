// SPDX-License-Identifier: Apache-2.0

#include "flow/src/constant_density_piso_test_access.hpp"
#include "hundun/boundary/basic_boundary.hpp"
#include "hundun/config/resolved_case.hpp"
#include "hundun/execution/execution.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/finite_volume/matrix_free_poisson.hpp"
#include "hundun/flow/constant_density_piso.hpp"
#include "hundun/linear/conjugate_gradient.hpp"
#include "hundun/linear/preconditioners.hpp"
#include "hundun/mesh/mesh_geometry.hpp"
#include "hundun/mesh/mesh_topology.hpp"
#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/exchange_plan.hpp"
#include "hundun/runtime/field_access_plan.hpp"
#include "hundun/runtime/field_registry.hpp"
#include "hundun/runtime/halo_exchange.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "linear/src/preconditioners_test_access.hpp"
#include "runtime/src/mpi_error.hpp"
#include "tests/support/test_main.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
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

hundun::config::FlowCaseConfig open_case(double rho_ref,
                                         double outlet_pressure,
                                         double inlet_velocity = 0.0) {
  auto config = periodic_case(rho_ref);
  config.scalars = {hundun::config::FlowScalarConfig{"alpha", 0.0}};
  config.boundaries[0].type = hundun::config::BoundaryType::velocity_inlet;
  config.boundaries[0].velocity_m_per_s =
      hundun::runtime::Real3{inlet_velocity, 0.0, 0.0};
  config.boundaries[0].thermal_authority =
      hundun::config::InletThermalAuthority::enthalpy;
  config.boundaries[0].enthalpy_J_per_kg = 1.0;
  config.boundaries[0].scalar_values =
      std::vector<hundun::config::InletScalarValue>{{"alpha", 2.0}};
  config.boundaries[1].type = hundun::config::BoundaryType::pressure_outlet;
  config.boundaries[1].pressure_perturbation_pa = outlet_pressure;
  for (std::size_t patch = 2U; patch < config.boundaries.size(); ++patch) {
    config.boundaries[patch].type = hundun::config::BoundaryType::symmetry;
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

class DishonestConvergedSolver final : public hundun::linear::LinearSolver {
public:
  hundun::linear::SolveReport
  solve(const hundun::linear::LinearOperator &,
        hundun::linear::Preconditioner &,
        hundun::execution::VectorView<const double>,
        hundun::execution::VectorView<double>,
        const hundun::linear::SolveControl &) const override {
    hundun::linear::SolveReport report;
    report.reason = hundun::linear::SolveTerminationReason::converged;
    report.final_residual = 0.0;
    report.lowest_failing_rank = -1;
    return report;
  }
};

class SelectedRankPressureFailureSolver final
    : public hundun::linear::LinearSolver {
public:
  hundun::linear::SolveReport
  solve(const hundun::linear::LinearOperator &,
        hundun::linear::Preconditioner &,
        hundun::execution::VectorView<const double>,
        hundun::execution::VectorView<double>,
        const hundun::linear::SolveControl &) const override {
    hundun::linear::SolveReport report;
    report.reason =
        hundun::linear::SolveTerminationReason::non_finite_value;
    report.final_residual = std::numeric_limits<double>::infinity();
    report.lowest_failing_rank = 1;
    return report;
  }
};

class MpiOperationThrowingSolver final : public hundun::linear::LinearSolver {
public:
  hundun::linear::SolveReport
  solve(const hundun::linear::LinearOperator &,
        hundun::linear::Preconditioner &,
        hundun::execution::VectorView<const double>,
        hundun::execution::VectorView<double>,
        const hundun::linear::SolveControl &) const override {
    throw hundun::runtime::detail::MpiOperationError(
        "injected Task 18 operation failure");
  }
};

class RuntimeErrorThrowingSolver final
    : public hundun::linear::LinearSolver {
public:
  hundun::linear::SolveReport
  solve(const hundun::linear::LinearOperator &,
        hundun::linear::Preconditioner &,
        hundun::execution::VectorView<const double>,
        hundun::execution::VectorView<double>,
        const hundun::linear::SolveControl &) const override {
    throw hundun::runtime::Error(
        "injected Task 18 unclassified runtime failure");
  }
};

class NonRuntimeThrowingSolver final
    : public hundun::linear::LinearSolver {
public:
  hundun::linear::SolveReport
  solve(const hundun::linear::LinearOperator &,
        hundun::linear::Preconditioner &,
        hundun::execution::VectorView<const double>,
        hundun::execution::VectorView<double>,
        const hundun::linear::SolveControl &) const override {
    throw std::logic_error("injected Task 18 non-runtime failure");
  }
};

std::uint64_t fp64_bits(double value) noexcept {
  std::uint64_t bits{};
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

bool fp64_vector_bitwise_equal(const std::vector<double> &left,
                               const std::vector<double> &right) noexcept {
  if (left.size() != right.size()) {
    return false;
  }
  for (std::size_t index = 0; index < left.size(); ++index) {
    if (fp64_bits(left[index]) != fp64_bits(right[index])) {
      return false;
    }
  }
  return true;
}

bool flow_layer_values_bitwise_equal(
    const hundun::flow::FlowLayerValues &left,
    const hundun::flow::FlowLayerValues &right) noexcept {
  if (!fp64_vector_bitwise_equal(left.density, right.density) ||
      !fp64_vector_bitwise_equal(left.velocity, right.velocity) ||
      !fp64_vector_bitwise_equal(left.mechanical_pressure,
                                 right.mechanical_pressure) ||
      !fp64_vector_bitwise_equal(left.face_velocity, right.face_velocity) ||
      !fp64_vector_bitwise_equal(left.face_mass_flux,
                                 right.face_mass_flux) ||
      left.transported_cell_fields.size() !=
          right.transported_cell_fields.size()) {
    return false;
  }
  for (std::size_t field = 0; field < left.transported_cell_fields.size();
       ++field) {
    if (!fp64_vector_bitwise_equal(left.transported_cell_fields[field],
                                   right.transported_cell_fields[field])) {
      return false;
    }
  }
  return true;
}

void check_layer_equal(const hundun::flow::FlowLayerValues &left,
                       const hundun::flow::FlowLayerValues &right) {
  HUNDUN_CHECK(flow_layer_values_bitwise_equal(left, right));
}

void check_flow_layer_bitwise_oracle() {
  hundun::flow::FlowLayerValues baseline;
  baseline.density = {0.0, 1.0};
  baseline.velocity = {2.0, 3.0, 4.0};
  baseline.mechanical_pressure = {5.0};
  baseline.face_velocity = {6.0, 7.0, 8.0};
  baseline.face_mass_flux = {9.0};
  baseline.transported_cell_fields = {{0.0, 10.0}, {11.0}};

  const auto exact_copy = baseline;
  HUNDUN_CHECK(flow_layer_values_bitwise_equal(baseline, exact_copy));

  auto ordinary_field_mutation = baseline;
  ordinary_field_mutation.density.front() = -0.0;
  HUNDUN_CHECK(
      !flow_layer_values_bitwise_equal(baseline, ordinary_field_mutation));

  auto transported_field_mutation = baseline;
  transported_field_mutation.transported_cell_fields.front().front() = -0.0;
  HUNDUN_CHECK(
      !flow_layer_values_bitwise_equal(baseline, transported_field_mutation));
}

void check_metadata_equal(const hundun::flow::AcceptedStepMetadata &left,
                          const hundun::flow::AcceptedStepMetadata &right) {
  HUNDUN_CHECK(left.step == right.step);
  HUNDUN_CHECK(fp64_bits(left.time_s) == fp64_bits(right.time_s));
  HUNDUN_CHECK(fp64_bits(left.dt_s) == fp64_bits(right.dt_s));
  HUNDUN_CHECK(fp64_bits(left.previous_dt_s) ==
               fp64_bits(right.previous_dt_s));
  HUNDUN_CHECK(left.order == right.order);
}

void check_solve_report_equal(const hundun::linear::SolveReport &left,
                              const hundun::linear::SolveReport &right) {
  HUNDUN_CHECK(left.reason == right.reason);
  HUNDUN_CHECK(left.iterations == right.iterations);
  HUNDUN_CHECK(fp64_bits(left.initial_residual) ==
               fp64_bits(right.initial_residual));
  HUNDUN_CHECK(fp64_bits(left.recursive_residual) ==
               fp64_bits(right.recursive_residual));
  HUNDUN_CHECK(fp64_bits(left.final_residual) ==
               fp64_bits(right.final_residual));
  HUNDUN_CHECK(left.matvec_count == right.matvec_count);
  HUNDUN_CHECK(left.preconditioner_apply_count ==
               right.preconditioner_apply_count);
  HUNDUN_CHECK(left.global_reduction_count == right.global_reduction_count);
  HUNDUN_CHECK(left.lowest_failing_rank == right.lowest_failing_rank);
}

void check_attempt_report_equal(const hundun::flow::StepAttemptReport &left,
                                const hundun::flow::StepAttemptReport &right) {
  HUNDUN_CHECK(left.disposition == right.disposition);
  HUNDUN_CHECK(left.reason == right.reason);
  HUNDUN_CHECK(left.lowest_failing_rank == right.lowest_failing_rank);
  HUNDUN_CHECK(left.pressure_corrector_count ==
               right.pressure_corrector_count);
  HUNDUN_CHECK(fp64_bits(left.attempted_dt_s) ==
               fp64_bits(right.attempted_dt_s));
  HUNDUN_CHECK(fp64_bits(left.suggested_dt_s) ==
               fp64_bits(right.suggested_dt_s));
  for (std::size_t index = 0; index < left.momentum.components.size();
       ++index) {
    check_solve_report_equal(left.momentum.components[index],
                             right.momentum.components[index]);
  }
  for (std::size_t index = 0; index < left.pressure.size(); ++index) {
    check_solve_report_equal(left.pressure[index], right.pressure[index]);
  }
  HUNDUN_CHECK(fp64_bits(left.final_continuity_normalized_l2) ==
               fp64_bits(right.final_continuity_normalized_l2));
  HUNDUN_CHECK(fp64_bits(left.final_pressure_residual_l2) ==
               fp64_bits(right.final_pressure_residual_l2));
  HUNDUN_CHECK(left.final_momentum_normalized_l2 ==
               right.final_momentum_normalized_l2);
  HUNDUN_CHECK(left.final_transport_normalized_l2 ==
               right.final_transport_normalized_l2);
  HUNDUN_CHECK(fp64_bits(left.final_mass_relative_conservation_defect) ==
               fp64_bits(right.final_mass_relative_conservation_defect));
  HUNDUN_CHECK(left.final_momentum_relative_conservation_defect ==
               right.final_momentum_relative_conservation_defect);
  HUNDUN_CHECK(left.final_transport_relative_conservation_defect ==
               right.final_transport_relative_conservation_defect);
  HUNDUN_CHECK(left.final_backflow_evidence.has_value() ==
               right.final_backflow_evidence.has_value());
}

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
  const auto rank_mismatch_extra_field =
      registry.declare_field(cell_field("rank_mismatch_extra", 1U));
  const auto direct_momentum_diagonal =
      registry.declare_field(cell_field("direct_momentum_diagonal", 3U));
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
  using TestAccess = hundun::flow::test::ConstantDensityPisoTestAccess;

  auto direct_coupler = hundun::flow::PisoCoupler::create(
      decomposition, topology, geometry, boundaries, mpi, execution, halo,
      pressure_solver, pressure_preconditioner);
  const auto direct_report = direct_coupler.correct(
      state, rho_ref,
      state.layer(hundun::flow::FlowLayer::committed)
          .view<double>(fields.velocity),
      {});
  HUNDUN_CHECK(
      direct_report.disposition ==
      hundun::flow::PressureCorrectionDisposition::non_retryable_failure);
  HUNDUN_CHECK(direct_report.reason ==
               hundun::flow::StepFailureReason::invalid_input);
  HUNDUN_CHECK(direct_report.lowest_failing_rank == 0);
  HUNDUN_CHECK(!direct_report.accepted);

  hundun::linear::JacobiPreconditioner direct_diagonal_preconditioner(
      execution);
  auto direct_diagonal_coupler = hundun::flow::PisoCoupler::create(
      decomposition, topology, geometry, boundaries, mpi, execution, halo,
      pressure_solver, direct_diagonal_preconditioner);
  auto direct_diagonal_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0,
       hundun::flow::MomentumTimeOrder::backward_euler});
  direct_diagonal_state.seed_accepted_layers(initial, initial);
  direct_diagonal_state.begin_attempt();

  const hundun::runtime::Int3 bad_diagonal_global{
      extent.x / 2, extent.y / 2, extent.z / 2};
  const auto bad_diagonal_global_id =
      topology.global_cell_id(bad_diagonal_global);
  const auto locally_observed_bad_diagonal =
      topology.find_local_cell(bad_diagonal_global_id);
  const auto observer_status = hundun::runtime::collective_status(
      mpi, !locally_observed_bad_diagonal.has_value(),
      "rank observes selected momentum diagonal");
  HUNDUN_CHECK(!observer_status.ok);

  const auto set_direct_diagonal = [&](double selected_value) {
    auto diagonal = direct_diagonal_state.trial_layer().view<double>(
        direct_momentum_diagonal);
    for (int k = -2; k < local.z + 2; ++k) {
      for (int j = -2; j < local.y + 2; ++j) {
        for (int i = -2; i < local.x + 2; ++i) {
          for (int component = 0; component < 3; ++component) {
            diagonal(i, j, k, component) = 1.0;
          }
        }
      }
    }
    if (locally_observed_bad_diagonal.has_value() &&
        topology.cell_ownership(*locally_observed_bad_diagonal) ==
            hundun::mesh::EntityOwnership::owned) {
      const auto direct_box = decomposition.owned_box();
      const hundun::runtime::Int3 index{
          bad_diagonal_global.x - direct_box.begin.x,
          bad_diagonal_global.y - direct_box.begin.y,
          bad_diagonal_global.z - direct_box.begin.z};
      for (int component = 0; component < 3; ++component) {
        diagonal(index.x, index.y, index.z, component) = selected_value;
      }
    }
  };

  const auto check_direct_diagonal_rollback =
      [&](const hundun::flow::FlowLayerValues &history_before,
          const hundun::flow::FlowLayerValues &committed_before,
          const hundun::flow::FlowLayerValues &trial_before,
          const hundun::flow::AcceptedStepMetadata &metadata_before) {
        check_layer_equal(
            direct_diagonal_state.snapshot(hundun::flow::FlowLayer::history),
            history_before);
        check_layer_equal(direct_diagonal_state.snapshot(
                              hundun::flow::FlowLayer::committed),
                          committed_before);
        check_layer_equal(
            direct_diagonal_state.snapshot(hundun::flow::FlowLayer::trial),
            trial_before);
        check_metadata_equal(direct_diagonal_state.metadata(), metadata_before);
      };

  const auto require_direct_diagonal_rejection = [&](double selected_value) {
    set_direct_diagonal(selected_value);
    halo.exchange(direct_diagonal_state.trial_layer(),
                  direct_momentum_diagonal);
    const auto history_before =
        direct_diagonal_state.snapshot(hundun::flow::FlowLayer::history);
    const auto committed_before =
        direct_diagonal_state.snapshot(hundun::flow::FlowLayer::committed);
    const auto trial_before =
        direct_diagonal_state.snapshot(hundun::flow::FlowLayer::trial);
    const auto metadata_before = direct_diagonal_state.metadata();
    const auto rejection = direct_diagonal_coupler.correct(
        direct_diagonal_state, rho_ref,
        direct_diagonal_state.layer(hundun::flow::FlowLayer::trial)
            .view<double>(direct_momentum_diagonal),
        {});
    HUNDUN_CHECK(
        rejection.disposition ==
        hundun::flow::PressureCorrectionDisposition::non_retryable_failure);
    HUNDUN_CHECK(rejection.reason ==
                 hundun::flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(rejection.lowest_failing_rank ==
                 observer_status.failing_rank);
    HUNDUN_CHECK(!rejection.accepted);
    check_direct_diagonal_rollback(history_before, committed_before,
                                   trial_before, metadata_before);
  };

  require_direct_diagonal_rejection(-0.25);
  require_direct_diagonal_rejection(
      std::numeric_limits<double>::quiet_NaN());

  set_direct_diagonal(1.0);
  halo.exchange(direct_diagonal_state.trial_layer(), direct_momentum_diagonal);
  const auto direct_diagonal_retry = direct_diagonal_coupler.correct(
      direct_diagonal_state, rho_ref,
      direct_diagonal_state.layer(hundun::flow::FlowLayer::trial)
          .view<double>(direct_momentum_diagonal),
      {});
  HUNDUN_CHECK(
      direct_diagonal_retry.disposition ==
      hundun::flow::PressureCorrectionDisposition::accepted);
  HUNDUN_CHECK(direct_diagonal_retry.reason ==
               hundun::flow::StepFailureReason::none);
  HUNDUN_CHECK(direct_diagonal_retry.lowest_failing_rank == -1);
  HUNDUN_CHECK(direct_diagonal_retry.accepted);
  direct_diagonal_state.rollback_attempt();

  auto invalid_quantity_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  invalid_quantity_state.seed_accepted_layers(initial, initial);
  hundun::runtime::FieldAccessPlan preflight_liveness(registry);
  preflight_liveness.declare_access(17U, 17U, fields.velocity,
                                    hundun::runtime::AccessMode::read);
  preflight_liveness.freeze();
  const auto preflight_trial =
      invalid_quantity_state.trial_layer().acquire_read<double>(
          preflight_liveness, 17U, 17U, fields.velocity);
  auto invalid_quantity_flow =
      hundun::flow::FixedStepConstantDensityFlow::create(
          decomposition, topology, geometry, boundaries, mpi, execution, halo,
          momentum_solver, {&mx, &my, &mz}, pressure_solver,
          pressure_preconditioner,
          {{fields.transported_cell_fields[0],
            hundun::finite_volume::FiniteVolumeQuantity::density(), 0.0}});
  const auto invalid_quantity_report = invalid_quantity_flow.attempt(
      invalid_quantity_state, rho_ref, 0.0, stencil, {}, {});
  HUNDUN_CHECK(invalid_quantity_report.disposition ==
               hundun::flow::StepAttemptDisposition::non_retryable_failure);
  HUNDUN_CHECK(invalid_quantity_report.reason ==
               hundun::flow::StepFailureReason::invalid_input);
  HUNDUN_CHECK(invalid_quantity_report.suggested_dt_s == 0.0);
  HUNDUN_CHECK(invalid_quantity_state.metadata().step == 0U);
  HUNDUN_CHECK_NEAR(preflight_trial(0, 0, 0, 0), 0.0, 0.0);

  const auto require_invalid_transport = [&](auto specs) {
    auto invalid_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    invalid_state.seed_accepted_layers(initial, initial);
    auto invalid_flow = hundun::flow::FixedStepConstantDensityFlow::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        momentum_solver, {&mx, &my, &mz}, pressure_solver,
        pressure_preconditioner, std::move(specs));
    const auto invalid_report =
        invalid_flow.attempt(invalid_state, rho_ref, 0.0, stencil, {}, {});
    HUNDUN_CHECK(invalid_report.disposition ==
                 hundun::flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(invalid_report.reason ==
                 hundun::flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(invalid_report.suggested_dt_s == 0.0);
    HUNDUN_CHECK(invalid_state.metadata().step == 0U);
  };
  require_invalid_transport(std::vector<hundun::flow::ConstantDensityTransportSpec>{
      {fields.transported_cell_fields[0],
       hundun::finite_volume::FiniteVolumeQuantity::velocity(), 0.0}});
  require_invalid_transport(std::vector<hundun::flow::ConstantDensityTransportSpec>{
      {fields.transported_cell_fields[0],
       hundun::finite_volume::FiniteVolumeQuantity::scalar(1U), 0.0}});
  require_invalid_transport(std::vector<hundun::flow::ConstantDensityTransportSpec>{
      {fields.transported_cell_fields[0],
       hundun::finite_volume::FiniteVolumeQuantity::scalar(0U), 0.0},
      {fields.transported_cell_fields[0],
       hundun::finite_volume::FiniteVolumeQuantity::scalar(0U), 0.0}});
  require_invalid_transport(std::vector<hundun::flow::ConstantDensityTransportSpec>{
      {fields.velocity, hundun::finite_volume::FiniteVolumeQuantity::scalar(0U),
       0.0}});
  require_invalid_transport(std::vector<hundun::flow::ConstantDensityTransportSpec>{
      {fields.transported_cell_fields[0],
       hundun::finite_volume::FiniteVolumeQuantity::scalar(0U), -1.0}});
  require_invalid_transport(std::vector<hundun::flow::ConstantDensityTransportSpec>{
      {fields.transported_cell_fields[0],
       hundun::finite_volume::FiniteVolumeQuantity::scalar(0U),
       std::numeric_limits<double>::quiet_NaN()}});
  auto malformed_enthalpy =
      hundun::finite_volume::FiniteVolumeQuantity::enthalpy();
  malformed_enthalpy.scalar_index = 1U;
  require_invalid_transport(std::vector<hundun::flow::ConstantDensityTransportSpec>{
      {fields.transported_cell_fields[0], malformed_enthalpy, 0.0}});

  const auto require_control_preflight =
      [&](hundun::linear::SolveControl momentum_control,
          hundun::linear::SolveControl pressure_control,
          int expected_failing_rank) {
        auto invalid_state = hundun::flow::FlowState::create(
            registry, {local, topology.local_face_count()}, fields,
            {0U, 0.0, 0.01, 0.0,
             hundun::flow::MomentumTimeOrder::backward_euler});
        invalid_state.seed_accepted_layers(initial, initial);
        const auto committed_before =
            invalid_state.snapshot(hundun::flow::FlowLayer::committed);
        const auto trial_before = invalid_state.trial_layer().acquire_read<double>(
            preflight_liveness, 17U, 17U, fields.velocity);
        const auto invalid_report = flow.attempt(
            invalid_state, rho_ref, 0.0, stencil, momentum_control,
            pressure_control);
        HUNDUN_CHECK(invalid_report.disposition ==
                     hundun::flow::StepAttemptDisposition::
                         non_retryable_failure);
        HUNDUN_CHECK(invalid_report.reason ==
                     hundun::flow::StepFailureReason::invalid_input);
        HUNDUN_CHECK(invalid_report.lowest_failing_rank ==
                     expected_failing_rank);
        HUNDUN_CHECK(invalid_report.suggested_dt_s == 0.0);
        HUNDUN_CHECK(invalid_report.pressure_corrector_count == 0U);
        HUNDUN_CHECK(invalid_state.metadata().step == 0U);
        check_layer_equal(
            invalid_state.snapshot(hundun::flow::FlowLayer::committed),
            committed_before);
        HUNDUN_CHECK_NEAR(trial_before(0, 0, 0, 0), 0.0, 0.0);
      };

  const int local_invalid_rank = mpi.size() == 1 ? 0 : 1;
  hundun::linear::SolveControl invalid_momentum_control;
  if (mpi.rank() == local_invalid_rank) {
    invalid_momentum_control.atol = -1.0;
  }
  require_control_preflight(invalid_momentum_control, {}, local_invalid_rank);
  hundun::linear::SolveControl invalid_pressure_control;
  if (mpi.rank() == local_invalid_rank) {
    invalid_pressure_control.residual_recompute_interval = 0U;
  }
  require_control_preflight({}, invalid_pressure_control, local_invalid_rank);
  if (mpi.size() > 1) {
    hundun::linear::SolveControl mismatched_momentum_control;
    if (mpi.rank() == 1) {
      mismatched_momentum_control.max_iterations =
          (UINT64_C(1) << 60U) + 500U;
    }
    require_control_preflight(mismatched_momentum_control, {}, 1);
    hundun::linear::SolveControl mismatched_pressure_control;
    if (mpi.rank() == 1) {
      mismatched_pressure_control.atol = -0.0;
    }
    require_control_preflight({}, mismatched_pressure_control, 1);
  }

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
  if (mpi.rank() == 0) {
    std::cout << "TASK18_FINAL_GATES momentum_residual="
              << report.final_momentum_normalized_l2[0] << ','
              << report.final_momentum_normalized_l2[1] << ','
              << report.final_momentum_normalized_l2[2]
              << " transport_residual="
              << report.final_transport_normalized_l2[0]
              << " mass_defect="
              << report.final_mass_relative_conservation_defect
              << " momentum_defect="
              << report.final_momentum_relative_conservation_defect[0] << ','
              << report.final_momentum_relative_conservation_defect[1] << ','
              << report.final_momentum_relative_conservation_defect[2]
              << " transport_defect="
              << report.final_transport_relative_conservation_defect[0]
              << '\n';
  }
  HUNDUN_CHECK(state.metadata().step == 1U);
  HUNDUN_CHECK(state.metadata().time_s == 0.01);
  HUNDUN_CHECK(fp64_vector_bitwise_equal(
      state.snapshot(hundun::flow::FlowLayer::committed)
          .transported_cell_fields[0],
      initial.transported_cell_fields[0]));

  const auto make_zero_state = [&] {
    auto candidate = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    candidate.seed_accepted_layers(initial, initial);
    return candidate;
  };

  if (mpi.size() > 1) {
    const auto warm_operator = TestAccess::pressure_operator_snapshot(flow);
    const auto warm_jacobi =
        hundun::linear::test::PreconditionerTestAccess::jacobi_storage(
            pressure_preconditioner);
    HUNDUN_CHECK(warm_operator.present);
    HUNDUN_CHECK(warm_operator.identity != 0U);
    HUNDUN_CHECK(warm_jacobi.cache_valid);
    HUNDUN_CHECK(warm_jacobi.revision == warm_operator.revision);

    auto refresh_state = make_zero_state();
    const auto refresh_history_before =
        refresh_state.snapshot(hundun::flow::FlowLayer::history);
    const auto refresh_committed_before =
        refresh_state.snapshot(hundun::flow::FlowLayer::committed);
    const auto refresh_metadata_before = refresh_state.metadata();
    const auto changed_stencil = hundun::flow::make_momentum_time_stencil(
        hundun::flow::MomentumTimeOrder::backward_euler, 0.005, 0.0);
    TestAccess::set_pressure_operator_refresh_failure_rank(1);
    const auto refresh_failure = flow.attempt(
        refresh_state, rho_ref, 0.0, changed_stencil, {}, {});
    HUNDUN_CHECK(
        refresh_failure.disposition ==
        hundun::flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(refresh_failure.reason ==
                 hundun::flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(refresh_failure.lowest_failing_rank == 1);
    HUNDUN_CHECK(refresh_failure.pressure_corrector_count == 0U);
    HUNDUN_CHECK(refresh_failure.suggested_dt_s == 0.0);
    check_metadata_equal(refresh_state.metadata(), refresh_metadata_before);
    check_layer_equal(
        refresh_state.snapshot(hundun::flow::FlowLayer::history),
        refresh_history_before);
    check_layer_equal(
        refresh_state.snapshot(hundun::flow::FlowLayer::committed),
        refresh_committed_before);
    const auto after_failure_operator =
        TestAccess::pressure_operator_snapshot(flow);
    const auto after_failure_jacobi =
        hundun::linear::test::PreconditionerTestAccess::jacobi_storage(
            pressure_preconditioner);
    HUNDUN_CHECK(after_failure_operator.present);
    HUNDUN_CHECK(after_failure_operator.identity == warm_operator.identity);
    HUNDUN_CHECK(after_failure_operator.revision == warm_operator.revision);
    HUNDUN_CHECK(after_failure_jacobi.cache_valid == warm_jacobi.cache_valid);
    HUNDUN_CHECK(after_failure_jacobi.revision == warm_jacobi.revision);
    HUNDUN_CHECK(after_failure_jacobi.allocation_identities ==
                 warm_jacobi.allocation_identities);

    TestAccess::reset();
    const auto retry_report = flow.attempt(
        refresh_state, rho_ref, 0.0, changed_stencil, {}, {});
    HUNDUN_CHECK(retry_report.disposition ==
                 hundun::flow::StepAttemptDisposition::committed);
    const auto after_retry_operator =
        TestAccess::pressure_operator_snapshot(flow);
    const auto after_retry_jacobi =
        hundun::linear::test::PreconditionerTestAccess::jacobi_storage(
            pressure_preconditioner);
    HUNDUN_CHECK(after_retry_operator.present);
    HUNDUN_CHECK(after_retry_operator.identity == warm_operator.identity);
    HUNDUN_CHECK(after_retry_operator.revision == warm_operator.revision + 1U);
    HUNDUN_CHECK(after_retry_jacobi.cache_valid);
    HUNDUN_CHECK(after_retry_jacobi.revision == after_retry_operator.revision);
    auto warm_allocations = warm_jacobi.allocation_identities;
    auto retry_allocations = after_retry_jacobi.allocation_identities;
    std::sort(warm_allocations.begin(), warm_allocations.end());
    std::sort(retry_allocations.begin(), retry_allocations.end());
    HUNDUN_CHECK(retry_allocations == warm_allocations);

    hundun::linear::JacobiPreconditioner reference_x(execution);
    hundun::linear::JacobiPreconditioner reference_y(execution);
    hundun::linear::JacobiPreconditioner reference_z(execution);
    hundun::linear::JacobiPreconditioner reference_pressure(execution);
    auto reference_flow = hundun::flow::FixedStepConstantDensityFlow::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        momentum_solver, {&reference_x, &reference_y, &reference_z},
        pressure_solver, reference_pressure,
        {{fields.transported_cell_fields[0],
          hundun::finite_volume::FiniteVolumeQuantity::scalar(0U), 0.0}});
    auto reference_state = make_zero_state();
    const auto reference_report = reference_flow.attempt(
        reference_state, rho_ref, 0.0, changed_stencil, {}, {});
    check_attempt_report_equal(retry_report, reference_report);
    check_metadata_equal(refresh_state.metadata(), reference_state.metadata());
    check_layer_equal(
        refresh_state.snapshot(hundun::flow::FlowLayer::history),
        reference_state.snapshot(hundun::flow::FlowLayer::history));
    check_layer_equal(
        refresh_state.snapshot(hundun::flow::FlowLayer::committed),
        reference_state.snapshot(hundun::flow::FlowLayer::committed));

    SelectedRankPressureFailureSolver direct_rank_failure_solver;
    hundun::linear::JacobiPreconditioner direct_rank_failure_preconditioner(
        execution);
    auto rank_failure_coupler = hundun::flow::PisoCoupler::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        direct_rank_failure_solver, direct_rank_failure_preconditioner);
    auto direct_rank_failure_state = make_zero_state();
    auto positive_diagonal_values = initial;
    std::fill(positive_diagonal_values.velocity.begin(),
              positive_diagonal_values.velocity.end(), 1.0);
    direct_rank_failure_state.seed_accepted_layers(positive_diagonal_values,
                                                   positive_diagonal_values);
    direct_rank_failure_state.begin_attempt();
    halo.exchange(direct_rank_failure_state.trial_layer(), fields.velocity);
    const auto direct_rank_failure_report = rank_failure_coupler.correct(
        direct_rank_failure_state, rho_ref,
        direct_rank_failure_state.layer(hundun::flow::FlowLayer::trial)
            .view<double>(fields.velocity),
        {});
    HUNDUN_CHECK(
        direct_rank_failure_report.disposition ==
        hundun::flow::PressureCorrectionDisposition::recoverable_failure);
    HUNDUN_CHECK(direct_rank_failure_report.reason ==
                 hundun::flow::StepFailureReason::pressure_linear_solve);
    HUNDUN_CHECK(
        direct_rank_failure_report.solve.lowest_failing_rank == 1);
    HUNDUN_CHECK(direct_rank_failure_report.lowest_failing_rank == 1);
    HUNDUN_CHECK(!direct_rank_failure_report.accepted);
    direct_rank_failure_state.rollback_attempt();

    SelectedRankPressureFailureSolver fixed_rank_failure_solver;
    hundun::linear::JacobiPreconditioner fixed_rank_failure_preconditioner(
        execution);
    auto fixed_rank_failure_flow =
        hundun::flow::FixedStepConstantDensityFlow::create(
            decomposition, topology, geometry, boundaries, mpi, execution,
            halo, momentum_solver, {&mx, &my, &mz},
            fixed_rank_failure_solver, fixed_rank_failure_preconditioner,
            {{fields.transported_cell_fields[0],
              hundun::finite_volume::FiniteVolumeQuantity::scalar(0U),
              0.0}});
    auto fixed_rank_failure_state = make_zero_state();
    const auto fixed_rank_failure_history =
        fixed_rank_failure_state.snapshot(hundun::flow::FlowLayer::history);
    const auto fixed_rank_failure_committed =
        fixed_rank_failure_state.snapshot(hundun::flow::FlowLayer::committed);
    const auto fixed_rank_failure_metadata = fixed_rank_failure_state.metadata();
    const auto fixed_rank_failure_report = fixed_rank_failure_flow.attempt(
        fixed_rank_failure_state, rho_ref, 0.0, stencil, {}, {});
    HUNDUN_CHECK(
        fixed_rank_failure_report.disposition ==
        hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(fixed_rank_failure_report.reason ==
                 hundun::flow::StepFailureReason::pressure_linear_solve);
    HUNDUN_CHECK(
        fixed_rank_failure_report.pressure[0].lowest_failing_rank == 1);
    HUNDUN_CHECK(fixed_rank_failure_report.lowest_failing_rank == 1);
    HUNDUN_CHECK(fixed_rank_failure_report.pressure_corrector_count == 0U);
    HUNDUN_CHECK(fixed_rank_failure_report.suggested_dt_s == 0.005);
    check_metadata_equal(fixed_rank_failure_state.metadata(),
                         fixed_rank_failure_metadata);
    check_layer_equal(
        fixed_rank_failure_state.snapshot(hundun::flow::FlowLayer::history),
        fixed_rank_failure_history);
    check_layer_equal(
        fixed_rank_failure_state.snapshot(hundun::flow::FlowLayer::committed),
        fixed_rank_failure_committed);

    MpiOperationThrowingSolver operation_failure_solver;
    hundun::linear::JacobiPreconditioner operation_failure_x(execution);
    hundun::linear::JacobiPreconditioner operation_failure_y(execution);
    hundun::linear::JacobiPreconditioner operation_failure_z(execution);
    hundun::linear::JacobiPreconditioner operation_failure_pressure(execution);
    auto operation_failure_flow =
        hundun::flow::FixedStepConstantDensityFlow::create(
            decomposition, topology, geometry, boundaries, mpi, execution,
            halo, operation_failure_solver,
            {&operation_failure_x, &operation_failure_y, &operation_failure_z},
            pressure_solver, operation_failure_pressure,
            {{fields.transported_cell_fields[0],
              hundun::finite_volume::FiniteVolumeQuantity::scalar(0U),
              0.0}});
    auto operation_failure_state = make_zero_state();
    const auto operation_failure_history =
        operation_failure_state.snapshot(hundun::flow::FlowLayer::history);
    const auto operation_failure_committed =
        operation_failure_state.snapshot(hundun::flow::FlowLayer::committed);
    const auto operation_failure_metadata = operation_failure_state.metadata();
    const auto operation_failure_report = operation_failure_flow.attempt(
        operation_failure_state, rho_ref, 0.0, stencil, {}, {});
    HUNDUN_CHECK(
        operation_failure_report.disposition ==
        hundun::flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(operation_failure_report.reason ==
                 hundun::flow::StepFailureReason::collective_operation);
    HUNDUN_CHECK(operation_failure_report.lowest_failing_rank == -1);
    HUNDUN_CHECK(operation_failure_report.pressure_corrector_count == 0U);
    HUNDUN_CHECK(operation_failure_report.suggested_dt_s == 0.0);
    check_metadata_equal(operation_failure_state.metadata(),
                         operation_failure_metadata);
    check_layer_equal(
        operation_failure_state.snapshot(hundun::flow::FlowLayer::history),
        operation_failure_history);
    check_layer_equal(
        operation_failure_state.snapshot(hundun::flow::FlowLayer::committed),
        operation_failure_committed);

    const auto require_unclassified_failure =
        [&](const hundun::linear::LinearSolver &throwing_solver) {
          hundun::linear::JacobiPreconditioner failure_x(execution);
          hundun::linear::JacobiPreconditioner failure_y(execution);
          hundun::linear::JacobiPreconditioner failure_z(execution);
          hundun::linear::JacobiPreconditioner failure_pressure(execution);
          auto failure_flow =
              hundun::flow::FixedStepConstantDensityFlow::create(
                  decomposition, topology, geometry, boundaries, mpi,
                  execution, halo, throwing_solver,
                  {&failure_x, &failure_y, &failure_z}, pressure_solver,
                  failure_pressure,
                  {{fields.transported_cell_fields[0],
                    hundun::finite_volume::FiniteVolumeQuantity::scalar(0U),
                    0.0}});
          auto failure_state = make_zero_state();
          const auto history_before =
              failure_state.snapshot(hundun::flow::FlowLayer::history);
          const auto committed_before =
              failure_state.snapshot(hundun::flow::FlowLayer::committed);
          const auto metadata_before = failure_state.metadata();
          const auto failure_report = failure_flow.attempt(
              failure_state, rho_ref, 0.0, stencil, {}, {});
          HUNDUN_CHECK(
              failure_report.disposition ==
              hundun::flow::StepAttemptDisposition::non_retryable_failure);
          HUNDUN_CHECK(failure_report.reason ==
                       hundun::flow::StepFailureReason::invalid_input);
          HUNDUN_CHECK(failure_report.lowest_failing_rank == -1);
          HUNDUN_CHECK(failure_report.pressure_corrector_count == 0U);
          HUNDUN_CHECK(failure_report.suggested_dt_s == 0.0);
          check_metadata_equal(failure_state.metadata(), metadata_before);
          check_layer_equal(
              failure_state.snapshot(hundun::flow::FlowLayer::history),
              history_before);
          check_layer_equal(
              failure_state.snapshot(hundun::flow::FlowLayer::committed),
              committed_before);
        };
    RuntimeErrorThrowingSolver runtime_failure_solver;
    require_unclassified_failure(runtime_failure_solver);
    NonRuntimeThrowingSolver nonruntime_failure_solver;
    require_unclassified_failure(nonruntime_failure_solver);
  }

  const auto warm_workspace = TestAccess::mesh_workspace_snapshot(flow);
  auto workspace_success_state = make_zero_state();
  const auto workspace_success_report = flow.attempt(
      workspace_success_state, rho_ref, 0.0, stencil, {}, {});
  HUNDUN_CHECK(workspace_success_report.disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  const auto after_success_workspace =
      TestAccess::mesh_workspace_snapshot(flow);
  HUNDUN_CHECK(after_success_workspace.vector_count ==
               warm_workspace.vector_count);
  HUNDUN_CHECK(after_success_workspace.total_capacity ==
               warm_workspace.total_capacity);
  HUNDUN_CHECK(after_success_workspace.data_identity ==
               warm_workspace.data_identity);

  TestAccess::reset();
  TestAccess::set_attempt_failure_stage(
      hundun::flow::test::AttemptFailureStage::after_corrector_2);
  auto workspace_retry_state = make_zero_state();
  const auto workspace_retry_report = flow.attempt(
      workspace_retry_state, rho_ref, 0.0, stencil, {}, {});
  HUNDUN_CHECK(workspace_retry_report.disposition ==
               hundun::flow::StepAttemptDisposition::recoverable_failure);
  const auto after_retry_workspace =
      TestAccess::mesh_workspace_snapshot(flow);
  HUNDUN_CHECK(after_retry_workspace.vector_count ==
               warm_workspace.vector_count);
  HUNDUN_CHECK(after_retry_workspace.total_capacity ==
               warm_workspace.total_capacity);
  HUNDUN_CHECK(after_retry_workspace.data_identity ==
               warm_workspace.data_identity);
  TestAccess::reset();

  if (mpi.size() > 1) {
    TestAccess::set_pressure_operator_construction_failure_rank(1);
    auto construction_failure_state = make_zero_state();
    const auto construction_failure_before = construction_failure_state.snapshot(
        hundun::flow::FlowLayer::committed);
    auto uncached_flow = hundun::flow::FixedStepConstantDensityFlow::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        momentum_solver, {&mx, &my, &mz}, pressure_solver,
        pressure_preconditioner,
        {{fields.transported_cell_fields[0],
          hundun::finite_volume::FiniteVolumeQuantity::scalar(0U), 0.0}});
    const auto construction_failure_report = uncached_flow.attempt(
        construction_failure_state, rho_ref, 0.0, stencil, {}, {});
    HUNDUN_CHECK(
        construction_failure_report.disposition ==
        hundun::flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(construction_failure_report.reason ==
                 hundun::flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(construction_failure_report.lowest_failing_rank == 1);
    HUNDUN_CHECK(construction_failure_report.pressure_corrector_count == 0U);
    HUNDUN_CHECK(construction_failure_report.suggested_dt_s == 0.0);
    HUNDUN_CHECK(construction_failure_state.metadata().step == 0U);
    check_layer_equal(
        construction_failure_state.snapshot(
            hundun::flow::FlowLayer::committed),
        construction_failure_before);
    TestAccess::reset();

    TestAccess::reset();
    if (mpi.rank() == 0) {
      TestAccess::set_final_momentum_norm_squares(
          0U, 1.21e-18, 1.0);
    } else if (mpi.rank() == 1) {
      TestAccess::set_final_momentum_norm_squares(0U, 0.0, 10000.0);
    } else {
      TestAccess::set_final_momentum_norm_squares(0U, 0.0, 0.0);
    }
    auto asymmetric_momentum_state = make_zero_state();
    const auto asymmetric_momentum_report = flow.attempt(
        asymmetric_momentum_state, rho_ref, 0.0, stencil, {}, {});
    HUNDUN_CHECK(asymmetric_momentum_report.disposition ==
                 hundun::flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(asymmetric_momentum_report.final_momentum_normalized_l2[0] <
                 1.0e-9);

    TestAccess::reset();
    if (mpi.rank() == 0) {
      TestAccess::set_final_transport_norm_squares(
          0U, 1.21e-18, 1.0);
    } else if (mpi.rank() == 1) {
      TestAccess::set_final_transport_norm_squares(0U, 0.0, 10000.0);
    } else {
      TestAccess::set_final_transport_norm_squares(0U, 0.0, 0.0);
    }
    auto asymmetric_transport_state = make_zero_state();
    const auto asymmetric_transport_report = flow.attempt(
        asymmetric_transport_state, rho_ref, 0.0, stencil, {}, {});
    HUNDUN_CHECK(asymmetric_transport_report.disposition ==
                 hundun::flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(asymmetric_transport_report.final_transport_normalized_l2[0] <
                 1.0e-9);
  }

  TestAccess::reset();
  if (mpi.rank() == 0) {
    TestAccess::set_final_momentum_norm_squares(0U, 4.0e-18, 1.0);
  } else {
    TestAccess::set_final_momentum_norm_squares(0U, 0.0, 0.0);
  }
  auto global_momentum_failure_state = make_zero_state();
  const auto global_momentum_failure_report = flow.attempt(
      global_momentum_failure_state, rho_ref, 0.0, stencil, {}, {});
  HUNDUN_CHECK(global_momentum_failure_report.disposition ==
               hundun::flow::StepAttemptDisposition::recoverable_failure);
  HUNDUN_CHECK(global_momentum_failure_report.reason ==
               hundun::flow::StepFailureReason::final_momentum_residual);
  HUNDUN_CHECK(global_momentum_failure_report.lowest_failing_rank == 0);

  TestAccess::reset();
  if (mpi.rank() == 0) {
    TestAccess::set_final_transport_norm_squares(0U, 4.0e-18, 1.0);
  } else {
    TestAccess::set_final_transport_norm_squares(0U, 0.0, 0.0);
  }
  auto global_transport_failure_state = make_zero_state();
  const auto global_transport_failure_report = flow.attempt(
      global_transport_failure_state, rho_ref, 0.0, stencil, {}, {});
  HUNDUN_CHECK(global_transport_failure_report.disposition ==
               hundun::flow::StepAttemptDisposition::recoverable_failure);
  HUNDUN_CHECK(global_transport_failure_report.reason ==
               hundun::flow::StepFailureReason::final_transport_residual);
  HUNDUN_CHECK(global_transport_failure_report.lowest_failing_rank == 0);

  if (mpi.size() > 1) {
    TestAccess::reset();
    if (mpi.rank() == 1) {
      TestAccess::set_final_momentum_norm_squares(
          0U, std::numeric_limits<double>::infinity(), 1.0);
    }
    auto non_finite_momentum_state = make_zero_state();
    const auto non_finite_momentum_report = flow.attempt(
        non_finite_momentum_state, rho_ref, 0.0, stencil, {}, {});
    HUNDUN_CHECK(non_finite_momentum_report.disposition ==
                 hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(non_finite_momentum_report.reason ==
                 hundun::flow::StepFailureReason::final_momentum_residual);
    HUNDUN_CHECK(non_finite_momentum_report.lowest_failing_rank == 1);

    TestAccess::reset();
    if (mpi.rank() == 1) {
      TestAccess::set_final_transport_norm_squares(
          0U, std::numeric_limits<double>::infinity(), 1.0);
    }
    auto non_finite_transport_state = make_zero_state();
    const auto non_finite_transport_report = flow.attempt(
        non_finite_transport_state, rho_ref, 0.0, stencil, {}, {});
    HUNDUN_CHECK(non_finite_transport_report.disposition ==
                 hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(non_finite_transport_report.reason ==
                 hundun::flow::StepFailureReason::final_transport_residual);
    HUNDUN_CHECK(non_finite_transport_report.lowest_failing_rank == 1);

    TestAccess::reset();
    if (mpi.rank() == 1) {
      TestAccess::force_momentum_conservation_aggregate_overflow(0U, true);
    }
    auto momentum_aggregate_state = make_zero_state();
    const auto momentum_aggregate_before =
        momentum_aggregate_state.snapshot(
            hundun::flow::FlowLayer::committed);
    const auto momentum_aggregate_report = flow.attempt(
        momentum_aggregate_state, rho_ref, 0.0, stencil, {}, {});
    HUNDUN_CHECK(momentum_aggregate_report.disposition ==
                 hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(momentum_aggregate_report.reason ==
                 hundun::flow::StepFailureReason::final_momentum_residual);
    HUNDUN_CHECK(momentum_aggregate_report.lowest_failing_rank == 1);
    HUNDUN_CHECK(momentum_aggregate_report.pressure_corrector_count == 2U);
    HUNDUN_CHECK(momentum_aggregate_report.suggested_dt_s == 0.005);
    check_layer_equal(
        momentum_aggregate_state.snapshot(
            hundun::flow::FlowLayer::committed),
        momentum_aggregate_before);
    HUNDUN_CHECK(momentum_aggregate_state.metadata().step == 0U);

    TestAccess::reset();
    if (mpi.rank() == 1) {
      TestAccess::force_transport_conservation_aggregate_overflow(0U, true);
    }
    auto transport_aggregate_state = make_zero_state();
    const auto transport_aggregate_before =
        transport_aggregate_state.snapshot(
            hundun::flow::FlowLayer::committed);
    const auto transport_aggregate_report = flow.attempt(
        transport_aggregate_state, rho_ref, 0.0, stencil, {}, {});
    HUNDUN_CHECK(transport_aggregate_report.disposition ==
                 hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(transport_aggregate_report.reason ==
                 hundun::flow::StepFailureReason::final_transport_residual);
    HUNDUN_CHECK(transport_aggregate_report.lowest_failing_rank == 1);
    HUNDUN_CHECK(transport_aggregate_report.pressure_corrector_count == 2U);
    HUNDUN_CHECK(transport_aggregate_report.suggested_dt_s == 0.005);
    check_layer_equal(
        transport_aggregate_state.snapshot(
            hundun::flow::FlowLayer::committed),
        transport_aggregate_before);
    HUNDUN_CHECK(transport_aggregate_state.metadata().step == 0U);

    TestAccess::reset();
    std::array<double, 8> globally_overflowing_parts{};
    globally_overflowing_parts[6] =
        0.75 * std::numeric_limits<double>::max();
    globally_overflowing_parts[7] =
        0.75 * std::numeric_limits<double>::max();
    TestAccess::set_momentum_conservation_parts(
        0U, globally_overflowing_parts);
    auto global_aggregate_state = make_zero_state();
    const auto global_aggregate_before =
        global_aggregate_state.snapshot(
            hundun::flow::FlowLayer::committed);
    const auto global_aggregate_report = flow.attempt(
        global_aggregate_state, rho_ref, 0.0, stencil, {}, {});
    HUNDUN_CHECK(global_aggregate_report.disposition ==
                 hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(global_aggregate_report.reason ==
                 hundun::flow::StepFailureReason::final_momentum_residual);
    HUNDUN_CHECK(global_aggregate_report.lowest_failing_rank == 0);
    HUNDUN_CHECK(global_aggregate_report.pressure_corrector_count == 2U);
    HUNDUN_CHECK(global_aggregate_report.suggested_dt_s == 0.005);
    check_layer_equal(
        global_aggregate_state.snapshot(
            hundun::flow::FlowLayer::committed),
        global_aggregate_before);
    HUNDUN_CHECK(global_aggregate_state.metadata().step == 0U);
    if (mpi.rank() == 0) {
      std::cout << "TASK18_R3_LOCAL_AGGREGATE momentum_rank="
                << momentum_aggregate_report.lowest_failing_rank
                << " transport_rank="
                << transport_aggregate_report.lowest_failing_rank
                << " global_rank="
                << global_aggregate_report.lowest_failing_rank
                << " correctors="
                << transport_aggregate_report.pressure_corrector_count
                << " suggested_dt="
                << transport_aggregate_report.suggested_dt_s << '\n';
    }
  }

  const auto run_cancellation_threshold_oracle = [&](double epsilon_multiple) {
    TestAccess::reset();
    const double epsilon = std::numeric_limits<double>::epsilon();
    std::array<double, 8> local_parts{};
    local_parts[4] = 1.0;
    local_parts[5] = 1.0;
    local_parts[6] = epsilon / stencil.dt_s;
    local_parts[7] = epsilon_multiple * epsilon / stencil.dt_s;
    TestAccess::set_momentum_conservation_parts(0U, local_parts);
    auto threshold_state = make_zero_state();
    const auto threshold_before = threshold_state.snapshot(
        hundun::flow::FlowLayer::committed);
    const auto threshold_report = flow.attempt(
        threshold_state, rho_ref, 0.0, stencil, {}, {});
    const auto diagnostic = TestAccess::last_momentum_conservation(0U);
    const double global_cq = static_cast<double>(mpi.size());
    const double expected_raw = epsilon * global_cq;
    const double expected_absolute_flux =
        epsilon_multiple * epsilon * global_cq / stencil.dt_s;
    const double expected_relative =
        epsilon_multiple < 64.0 ? epsilon : 1.0 / epsilon_multiple;
    const double tolerance =
        32.0 * epsilon * std::max(1.0, expected_relative);
    HUNDUN_CHECK_NEAR(diagnostic.quantity_n, 0.0, 0.0);
    HUNDUN_CHECK_NEAR(diagnostic.quantity_np1, 0.0, 0.0);
    HUNDUN_CHECK_NEAR(diagnostic.absolute_boundary_flux,
                      expected_absolute_flux,
                      32.0 * epsilon * expected_absolute_flux);
    HUNDUN_CHECK_NEAR(diagnostic.raw_defect, expected_raw,
                      32.0 * epsilon * expected_raw);
    HUNDUN_CHECK(diagnostic.raw_defect != 0.0);
    HUNDUN_CHECK_NEAR(
        threshold_report.final_momentum_relative_conservation_defect[0],
        expected_relative, tolerance);
    if (mpi.rank() == 0) {
      std::cout << "TASK18_R3_CQ_THRESHOLD multiple=" << epsilon_multiple
                << " raw=" << diagnostic.raw_defect
                << " absolute_flux=" << diagnostic.absolute_boundary_flux
                << " relative="
                << threshold_report
                       .final_momentum_relative_conservation_defect[0]
                << " disposition="
                << static_cast<int>(threshold_report.disposition) << '\n';
    }
    if (epsilon_multiple < 64.0) {
      HUNDUN_CHECK(threshold_report.disposition ==
                   hundun::flow::StepAttemptDisposition::committed);
    } else {
      HUNDUN_CHECK(threshold_report.disposition ==
                   hundun::flow::StepAttemptDisposition::recoverable_failure);
      HUNDUN_CHECK(threshold_report.reason ==
                   hundun::flow::StepFailureReason::final_conservation_defect);
      HUNDUN_CHECK(threshold_report.lowest_failing_rank == 0);
      HUNDUN_CHECK(threshold_report.pressure_corrector_count == 2U);
      HUNDUN_CHECK(threshold_report.suggested_dt_s == 0.005);
      check_layer_equal(
          threshold_state.snapshot(hundun::flow::FlowLayer::committed),
          threshold_before);
      HUNDUN_CHECK(threshold_state.metadata().step == 0U);
    }
  };
  run_cancellation_threshold_oracle(63.0);
  run_cancellation_threshold_oracle(65.0);

  TestAccess::reset();
  TestAccess::set_final_mass_defect_perturbation(7.5e-11);
  auto fine_conservation_state = make_zero_state();
  const auto fine_conservation_report = flow.attempt(
      fine_conservation_state, rho_ref, 0.0, stencil, {}, {});
  HUNDUN_CHECK(fine_conservation_report.disposition ==
               hundun::flow::StepAttemptDisposition::recoverable_failure);
  HUNDUN_CHECK(fine_conservation_report.reason ==
               hundun::flow::StepFailureReason::final_conservation_defect);
  HUNDUN_CHECK(fine_conservation_report.final_mass_relative_conservation_defect >
               5.0e-11);
  HUNDUN_CHECK(fine_conservation_report.final_mass_relative_conservation_defect <
               1.0e-10);
  TestAccess::reset();

  const auto first_committed =
      state.snapshot(hundun::flow::FlowLayer::committed);
  std::size_t exact_periodic_pairs = 0U;
  for (hundun::mesh::LocalFaceId face = 0;
       face < topology.local_face_count(); ++face) {
    const auto pair_global = topology.periodic_pair(face);
    if (!pair_global.has_value() ||
        topology.global_face_id(face) >= *pair_global) {
      continue;
    }
    const auto pair = topology.find_local_face(*pair_global);
    if (!pair.has_value())
      continue;
    HUNDUN_CHECK(first_committed.face_mass_flux[face] ==
                 -first_committed.face_mass_flux[*pair]);
    for (std::size_t component = 0; component < 3U; ++component) {
      HUNDUN_CHECK(first_committed.face_velocity[face * 3U + component] ==
                   first_committed.face_velocity[*pair * 3U + component]);
    }
    ++exact_periodic_pairs;
  }
  HUNDUN_CHECK(exact_periodic_pairs > 0U);

  auto offset_pressure = initial;
  std::fill(offset_pressure.mechanical_pressure.begin(),
            offset_pressure.mechanical_pressure.end(), 3.25);
  auto normalized_pressure_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0,
       hundun::flow::MomentumTimeOrder::backward_euler});
  normalized_pressure_state.seed_accepted_layers(offset_pressure,
                                                 offset_pressure);
  const auto normalized_pressure_report = flow.attempt(
      normalized_pressure_state, rho_ref, 0.0, stencil, {}, {});
  HUNDUN_CHECK(normalized_pressure_report.disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  const auto normalized_pressure =
      normalized_pressure_state.snapshot(hundun::flow::FlowLayer::committed);
  double pressure_moment[2]{};
  for (hundun::mesh::LocalCellId cell = 0;
       cell < topology.owned_cell_count(); ++cell) {
    const double volume = geometry.cell_volume_m3(cell);
    pressure_moment[0] += volume * normalized_pressure.mechanical_pressure[cell];
    pressure_moment[1] += volume;
  }
  mpi.allreduce_fp64_in_place(
      pressure_moment, 2U, hundun::runtime::Fp64ReductionOperation::sum);
  HUNDUN_CHECK(std::abs(pressure_moment[0] / pressure_moment[1]) <= 1.0e-12);

  auto dishonest_initial = initial;
  for (hundun::mesh::LocalCellId cell = 0;
       cell < topology.owned_cell_count(); ++cell) {
    const double velocity =
        0.2 * std::sin(2.0 * 3.14159265358979323846 *
                       geometry.cell_center_m(cell).x);
    dishonest_initial.velocity[cell * 3U] = velocity;
  }
  for (hundun::mesh::LocalFaceId face = 0;
       face < topology.local_face_count(); ++face) {
    const double velocity =
        0.2 * std::sin(2.0 * 3.14159265358979323846 *
                       geometry.face_center_m(face).x);
    dishonest_initial.face_velocity[face * 3U] = velocity;
    dishonest_initial.face_mass_flux[face] =
        rho_ref * velocity *
        geometry.face_area_vector_m2(face, hundun::mesh::FaceSide::owner).x;
  }
  DishonestConvergedSolver dishonest_pressure_solver;
  hundun::linear::JacobiPreconditioner dishonest_pressure_preconditioner(
      execution);
  auto dishonest_flow = hundun::flow::FixedStepConstantDensityFlow::create(
      decomposition, topology, geometry, boundaries, mpi, execution, halo,
      momentum_solver, {&mx, &my, &mz}, dishonest_pressure_solver,
      dishonest_pressure_preconditioner,
      {{fields.transported_cell_fields[0],
        hundun::finite_volume::FiniteVolumeQuantity::scalar(0U), 0.0}});
  auto dishonest_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0,
       hundun::flow::MomentumTimeOrder::backward_euler});
  dishonest_state.seed_accepted_layers(dishonest_initial, dishonest_initial);
  const auto dishonest_before =
      dishonest_state.snapshot(hundun::flow::FlowLayer::committed);
  const auto dishonest_report =
      dishonest_flow.attempt(dishonest_state, rho_ref, 0.0, stencil, {}, {});
  HUNDUN_CHECK(dishonest_report.disposition ==
               hundun::flow::StepAttemptDisposition::recoverable_failure);
  HUNDUN_CHECK(dishonest_report.reason ==
               hundun::flow::StepFailureReason::pressure_linear_solve);
  HUNDUN_CHECK(dishonest_report.pressure_corrector_count == 0U);
  HUNDUN_CHECK(dishonest_report.suggested_dt_s == 0.005);
  check_layer_equal(
      dishonest_state.snapshot(hundun::flow::FlowLayer::committed),
      dishonest_before);

  using FailureStage = hundun::flow::test::AttemptFailureStage;
  constexpr std::array<std::pair<FailureStage, std::uint32_t>, 8>
      transaction_stages{{
          {FailureStage::after_begin, 0U},
          {FailureStage::after_momentum, 0U},
          {FailureStage::after_face_predictor, 0U},
          {FailureStage::after_corrector_1, 1U},
          {FailureStage::after_provisional_transport, 1U},
          {FailureStage::after_corrector_2, 2U},
          {FailureStage::after_final_transport, 2U},
          {FailureStage::before_commit, 2U},
      }};
  for (const auto &[stage, expected_correctors] : transaction_stages) {
    TestAccess::reset();
    TestAccess::set_attempt_failure_stage(stage);
    auto staged_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    staged_state.seed_accepted_layers(initial, initial);
    const auto history_before =
        staged_state.snapshot(hundun::flow::FlowLayer::history);
    const auto committed_before =
        staged_state.snapshot(hundun::flow::FlowLayer::committed);
    const auto metadata_before = staged_state.metadata();
    const auto staged_report =
        flow.attempt(staged_state, rho_ref, 0.0, stencil, {}, {});
    HUNDUN_CHECK(staged_report.disposition ==
                 hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(staged_report.reason ==
                 hundun::flow::StepFailureReason::non_finite_trial);
    HUNDUN_CHECK(staged_report.lowest_failing_rank == 0);
    HUNDUN_CHECK(staged_report.pressure_corrector_count ==
                 expected_correctors);
    HUNDUN_CHECK(staged_report.suggested_dt_s == 0.005);
    check_layer_equal(
        staged_state.snapshot(hundun::flow::FlowLayer::history),
        history_before);
    check_layer_equal(
        staged_state.snapshot(hundun::flow::FlowLayer::committed),
        committed_before);
    check_metadata_equal(staged_state.metadata(), metadata_before);
  }
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
  HUNDUN_CHECK(fp64_vector_bitwise_equal(
      sentinel_state.snapshot(hundun::flow::FlowLayer::committed)
          .transported_cell_fields[0],
      initial.transported_cell_fields[0]));

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
  check_layer_equal(after_rollback, before_rollback);
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
  check_layer_equal(pressure_after, pressure_before);

  const int final_gate_failure_rank = mpi.size() == 1 ? 0 : 1;
  TestAccess::reset();
  if (mpi.rank() == final_gate_failure_rank) {
    TestAccess::force_final_momentum_perturbation(1U, 1.0e-4);
  }
  auto momentum_gate_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  momentum_gate_state.seed_accepted_layers(initial, initial);
  const auto momentum_gate_before =
      momentum_gate_state.snapshot(hundun::flow::FlowLayer::committed);
  const auto momentum_gate_report =
      flow.attempt(momentum_gate_state, rho_ref, 0.0, stencil, {}, {});
  HUNDUN_CHECK(momentum_gate_report.disposition ==
               hundun::flow::StepAttemptDisposition::recoverable_failure);
  HUNDUN_CHECK(momentum_gate_report.reason ==
               hundun::flow::StepFailureReason::final_momentum_residual);
  HUNDUN_CHECK(momentum_gate_report.lowest_failing_rank == 0);
  HUNDUN_CHECK(momentum_gate_report.pressure_corrector_count == 2U);
  HUNDUN_CHECK(momentum_gate_report.suggested_dt_s == 0.005);
  const auto momentum_gate_after =
      momentum_gate_state.snapshot(hundun::flow::FlowLayer::committed);
  check_layer_equal(momentum_gate_after, momentum_gate_before);
  HUNDUN_CHECK(momentum_gate_state.metadata().step == 0U);

  TestAccess::reset();
  if (mpi.rank() == final_gate_failure_rank) {
    TestAccess::force_final_transport_perturbation(0U, 1.0e-4);
  }
  auto transport_gate_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  transport_gate_state.seed_accepted_layers(initial, initial);
  const auto transport_gate_before =
      transport_gate_state.snapshot(hundun::flow::FlowLayer::committed);
  const auto transport_gate_report =
      flow.attempt(transport_gate_state, rho_ref, 0.0, stencil, {}, {});
  HUNDUN_CHECK(transport_gate_report.disposition ==
               hundun::flow::StepAttemptDisposition::recoverable_failure);
  HUNDUN_CHECK(transport_gate_report.reason ==
               hundun::flow::StepFailureReason::final_transport_residual);
  HUNDUN_CHECK(transport_gate_report.lowest_failing_rank == 0);
  HUNDUN_CHECK(transport_gate_report.pressure_corrector_count == 2U);
  HUNDUN_CHECK(transport_gate_report.suggested_dt_s == 0.005);
  const auto transport_gate_after =
      transport_gate_state.snapshot(hundun::flow::FlowLayer::committed);
  check_layer_equal(transport_gate_after, transport_gate_before);
  HUNDUN_CHECK(transport_gate_state.metadata().step == 0U);

  TestAccess::reset();
  TestAccess::force_final_conservation_failure(
      mpi.rank() == final_gate_failure_rank);
  auto conservation_gate_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0, hundun::flow::MomentumTimeOrder::backward_euler});
  conservation_gate_state.seed_accepted_layers(initial, initial);
  const auto conservation_gate_before =
      conservation_gate_state.snapshot(hundun::flow::FlowLayer::committed);
  const auto conservation_gate_report =
      flow.attempt(conservation_gate_state, rho_ref, 0.0, stencil, {}, {});
  HUNDUN_CHECK(conservation_gate_report.disposition ==
               hundun::flow::StepAttemptDisposition::recoverable_failure);
  HUNDUN_CHECK(conservation_gate_report.reason ==
               hundun::flow::StepFailureReason::final_conservation_defect);
  HUNDUN_CHECK(conservation_gate_report.lowest_failing_rank ==
               final_gate_failure_rank);
  HUNDUN_CHECK(conservation_gate_report.pressure_corrector_count == 2U);
  HUNDUN_CHECK(conservation_gate_report.suggested_dt_s == 0.005);
  const auto conservation_gate_after =
      conservation_gate_state.snapshot(hundun::flow::FlowLayer::committed);
  check_layer_equal(conservation_gate_after, conservation_gate_before);
  HUNDUN_CHECK(conservation_gate_state.metadata().step == 0U);
  TestAccess::reset();

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
  check_layer_equal(
      component_failure_state.snapshot(hundun::flow::FlowLayer::committed),
      initial);

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

    const auto require_rank_mismatched_transport =
        [&](hundun::flow::ConstantDensityTransportSpec local_spec) {
          auto mismatched_flow =
              hundun::flow::FixedStepConstantDensityFlow::create(
                  decomposition, topology, geometry, boundaries, mpi,
                  execution, halo, momentum_solver, {&mx, &my, &mz},
                  pressure_solver, pressure_preconditioner, {local_spec});
          auto candidate = hundun::flow::FlowState::create(
              registry, {local, topology.local_face_count()}, fields,
              {0U, 0.0, 0.01, 0.0,
               hundun::flow::MomentumTimeOrder::backward_euler});
          candidate.seed_accepted_layers(initial, initial);
          const auto mismatch_report =
              mismatched_flow.attempt(candidate, rho_ref, 0.0, stencil, {}, {});
          HUNDUN_CHECK(mismatch_report.disposition ==
                       hundun::flow::StepAttemptDisposition::
                           non_retryable_failure);
          HUNDUN_CHECK(mismatch_report.reason ==
                       hundun::flow::StepFailureReason::invalid_input);
          HUNDUN_CHECK(mismatch_report.lowest_failing_rank == 1);
          HUNDUN_CHECK(mismatch_report.suggested_dt_s == 0.0);
          HUNDUN_CHECK(candidate.metadata().step == 0U);
        };
    require_rank_mismatched_transport(
        {fields.transported_cell_fields[0],
         mpi.rank() == 0
             ? hundun::finite_volume::FiniteVolumeQuantity::enthalpy()
             : hundun::finite_volume::FiniteVolumeQuantity::scalar(0U),
         0.0});
    require_rank_mismatched_transport(
        {fields.transported_cell_fields[0],
         hundun::finite_volume::FiniteVolumeQuantity::scalar(
             mpi.rank() == 0 ? 0U : 1U),
         0.0});
    require_rank_mismatched_transport(
        {mpi.rank() == 0 ? fields.transported_cell_fields[0]
                         : rank_mismatch_extra_field,
         hundun::finite_volume::FiniteVolumeQuantity::scalar(0U), 0.0});

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
  HUNDUN_CHECK(cells > 1U);
  for (std::size_t component = 0; component < 3U; ++component) {
    bdf2_history.velocity[3U + component] +=
        static_cast<double>(cells) *
        (velocity_n[component] - velocity_nm1[component]);
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

  auto mutation_committed = initial;
  auto mutation_history = initial;
  for (hundun::mesh::LocalCellId cell = 0; cell < topology.owned_cell_count();
       ++cell) {
    const auto center = geometry.cell_center_m(cell);
    mutation_committed.velocity[cell * 3U] =
        0.25 + 0.05 * std::sin(2.0 * 3.14159265358979323846 * center.x);
    mutation_committed.velocity[cell * 3U + 1U] =
        -0.1 + 0.04 * std::cos(2.0 * 3.14159265358979323846 * center.x);
    mutation_committed.velocity[cell * 3U + 2U] =
        0.02 + 0.03 *
                   std::sin(2.0 * 3.14159265358979323846 * center.y);
    mutation_history.velocity[cell * 3U] =
        0.25 + 0.02 * std::cos(2.0 * 3.14159265358979323846 * center.x);
    mutation_history.velocity[cell * 3U + 1U] =
        -0.1 + 0.03 * std::sin(2.0 * 3.14159265358979323846 * center.x);
    mutation_history.velocity[cell * 3U + 2U] =
        0.02 + 0.01 *
                    std::cos(2.0 * 3.14159265358979323846 * center.y);
    mutation_committed.mechanical_pressure[cell] =
        0.07 * std::sin(2.0 * 3.14159265358979323846 * center.x);
    mutation_history.mechanical_pressure[cell] =
        0.03 * std::cos(2.0 * 3.14159265358979323846 * center.x);
    mutation_committed.transported_cell_fields[0][cell] =
        1.0 + 0.2 * std::sin(2.0 * 3.14159265358979323846 * center.x);
    mutation_history.transported_cell_fields[0][cell] =
        1.0 + 0.15 * std::cos(2.0 * 3.14159265358979323846 * center.x);
  }
  for (hundun::mesh::LocalFaceId face = 0; face < topology.local_face_count();
       ++face) {
    const auto area = geometry.face_area_vector_m2(
        face, hundun::mesh::FaceSide::owner);
    mutation_committed.face_mass_flux[face] = 0.08 * rho_ref * area.x;
    mutation_history.face_mass_flux[face] = 0.03 * rho_ref * area.x;
  }
  const auto mutation_stencil = hundun::flow::make_momentum_time_stencil(
      hundun::flow::MomentumTimeOrder::bdf2, 0.01, 0.01);
  const auto make_mutation_state = [&] {
    auto candidate = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {1U, 0.01, 0.01, 0.01, hundun::flow::MomentumTimeOrder::bdf2});
    candidate.seed_accepted_layers(mutation_history, mutation_committed);
    return candidate;
  };
  TestAccess::reset();
  auto mutation_baseline_state = make_mutation_state();
  const auto mutation_baseline_report = flow.attempt(
      mutation_baseline_state, rho_ref, 0.04, mutation_stencil, {}, {});
  if (mutation_baseline_report.disposition !=
          hundun::flow::StepAttemptDisposition::committed &&
      mpi.rank() == 0) {
    std::cerr << "TASK18_MUTATION_BASELINE reason="
              << static_cast<int>(mutation_baseline_report.reason)
              << " continuity="
              << mutation_baseline_report.final_continuity_normalized_l2
              << " momentum="
              << mutation_baseline_report.final_momentum_normalized_l2[0]
              << ','
              << mutation_baseline_report.final_momentum_normalized_l2[1]
              << ','
              << mutation_baseline_report.final_momentum_normalized_l2[2]
              << " transport="
              << (mutation_baseline_report.final_transport_normalized_l2.empty()
                      ? -1.0
                      : mutation_baseline_report
                            .final_transport_normalized_l2[0])
              << " conservation="
              << mutation_baseline_report.final_mass_relative_conservation_defect
              << ','
              << mutation_baseline_report
                     .final_momentum_relative_conservation_defect[0]
              << ','
              << mutation_baseline_report
                     .final_momentum_relative_conservation_defect[1]
              << ','
              << mutation_baseline_report
                     .final_momentum_relative_conservation_defect[2]
              << '\n';
  }
  HUNDUN_CHECK(mutation_baseline_report.disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  if (mpi.rank() == 0) {
    std::cout << "TASK18_MUTATION_BASELINE_FINAL_GATES momentum_residual="
              << mutation_baseline_report.final_momentum_normalized_l2[0]
              << ','
              << mutation_baseline_report.final_momentum_normalized_l2[1]
              << ','
              << mutation_baseline_report.final_momentum_normalized_l2[2]
              << " transport_residual="
              << mutation_baseline_report.final_transport_normalized_l2[0]
              << " mass_defect="
              << mutation_baseline_report.final_mass_relative_conservation_defect
              << " momentum_defect="
              << mutation_baseline_report
                     .final_momentum_relative_conservation_defect[0]
              << ','
              << mutation_baseline_report
                     .final_momentum_relative_conservation_defect[1]
              << ','
              << mutation_baseline_report
                     .final_momentum_relative_conservation_defect[2]
              << " transport_defect="
              << mutation_baseline_report
                     .final_transport_relative_conservation_defect[0]
              << '\n';
  }

  constexpr std::array momentum_mutations{
      hundun::flow::test::MomentumAssemblyMutation::omit_convection,
      hundun::flow::test::MomentumAssemblyMutation::omit_viscosity,
      hundun::flow::test::MomentumAssemblyMutation::omit_pressure,
      hundun::flow::test::MomentumAssemblyMutation::omit_alpha1,
      hundun::flow::test::MomentumAssemblyMutation::omit_alpha2,
      hundun::flow::test::MomentumAssemblyMutation::replace_history_flux};
  for (const auto mutation : momentum_mutations) {
    TestAccess::reset();
    TestAccess::set_momentum_assembly_mutation(mutation);
    auto mutation_state = make_mutation_state();
    const auto mutation_report = flow.attempt(
        mutation_state, rho_ref, 0.04, mutation_stencil, {}, {});
    HUNDUN_CHECK(mutation_report.disposition ==
                 hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(mutation_report.reason ==
                 hundun::flow::StepFailureReason::final_momentum_residual);
    HUNDUN_CHECK(mutation_report.pressure_corrector_count == 2U);
    HUNDUN_CHECK(mutation_state.metadata().step == 1U);
  }
  constexpr std::array transport_mutations{
      hundun::flow::test::TransportAssemblyMutation::omit_history_spatial,
      hundun::flow::test::TransportAssemblyMutation::use_provisional_flux};
  for (const auto mutation : transport_mutations) {
    TestAccess::reset();
    TestAccess::set_final_uniform_x_mass_flux(0.04 * rho_ref);
    TestAccess::set_transport_assembly_mutation(mutation);
    auto mutation_state = make_mutation_state();
    const auto mutation_report = flow.attempt(
        mutation_state, rho_ref, 0.04, mutation_stencil, {}, {});
    HUNDUN_CHECK(mutation_report.disposition ==
                 hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(mutation_report.reason ==
                 hundun::flow::StepFailureReason::final_transport_residual);
    HUNDUN_CHECK(mutation_report.pressure_corrector_count == 2U);
    HUNDUN_CHECK(mutation_state.metadata().step == 1U);
  }
  TestAccess::reset();

  constexpr double transport_velocity = 0.1;
  constexpr double transport_final_time = 0.2;
  constexpr std::array<double, 3> transport_dt{0.025, 0.0125, 0.00625};
  std::array<double, 3> transport_error{};
  TestAccess::reset();
  TestAccess::set_final_uniform_x_mass_flux(rho_ref * transport_velocity);
  for (std::size_t refinement = 0; refinement < transport_dt.size();
       ++refinement) {
    auto transported = initial;
    for (hundun::mesh::LocalCellId cell = 0;
         cell < topology.owned_cell_count(); ++cell) {
      const auto global = topology.global_cell(cell);
      const double mode_sign = (global.x % 2 == 0) ? 1.0 : -1.0;
      transported.transported_cell_fields[0][cell] = 1.0 + 0.2 * mode_sign;
    }
    auto transported_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, transport_dt[refinement], 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    transported_state.seed_accepted_layers(transported, transported);
    const auto steps = static_cast<std::uint32_t>(
        std::llround(transport_final_time / transport_dt[refinement]));
    HUNDUN_CHECK_NEAR(static_cast<double>(steps) * transport_dt[refinement],
                      transport_final_time, 1.0e-15);
    for (std::uint32_t step = 0; step < steps; ++step) {
      const auto order =
          step == 0U ? hundun::flow::MomentumTimeOrder::backward_euler
                     : hundun::flow::MomentumTimeOrder::bdf2;
      const auto transport_stencil = hundun::flow::make_momentum_time_stencil(
          order, transport_dt[refinement],
          step == 0U ? 0.0 : transport_dt[refinement]);
      const auto transport_report = flow.attempt(
          transported_state, rho_ref, 0.0, transport_stencil, {}, {});
      HUNDUN_CHECK(transport_report.disposition ==
                   hundun::flow::StepAttemptDisposition::committed);
      HUNDUN_CHECK(transport_report.pressure_corrector_count == 2U);
    }
    const auto final_transport =
        transported_state.snapshot(hundun::flow::FlowLayer::committed)
            .transported_cell_fields[0];
    const double dx = 1.0 / static_cast<double>(extent.x);
    const double decay =
        std::exp(-2.0 * transport_velocity * transport_final_time / dx);
    double error_square = 0.0;
    double volume_sum = 0.0;
    for (hundun::mesh::LocalCellId cell = 0;
         cell < topology.owned_cell_count(); ++cell) {
      const auto global = topology.global_cell(cell);
      const double mode_sign = (global.x % 2 == 0) ? 1.0 : -1.0;
      const double exact = 1.0 + 0.2 * mode_sign * decay;
      const double error = final_transport[cell] - exact;
      const double cell_volume = geometry.cell_volume_m3(cell);
      error_square += cell_volume * error * error;
      volume_sum += cell_volume;
    }
    double sums[2]{error_square, volume_sum};
    mpi.allreduce_fp64_in_place(
        sums, 2U, hundun::runtime::Fp64ReductionOperation::sum);
    transport_error[refinement] = std::sqrt(sums[0] / sums[1]);
    HUNDUN_CHECK(std::isfinite(transport_error[refinement]));
    HUNDUN_CHECK(transport_error[refinement] >
                 std::numeric_limits<double>::min());
  }
  const double transport_order_0 =
      std::log(transport_error[0] / transport_error[1]) / std::log(2.0);
  const double transport_order_1 =
      std::log(transport_error[1] / transport_error[2]) / std::log(2.0);
  if (mpi.rank() == 0) {
    std::cout << "TASK18_TRANSPORT_BDF2 errors=" << transport_error[0] << ','
              << transport_error[1] << ',' << transport_error[2]
              << " orders=" << transport_order_0 << ',' << transport_order_1
              << '\n';
  }
  HUNDUN_CHECK(transport_order_0 >= 1.8);
  HUNDUN_CHECK(transport_order_1 >= 1.8);
  TestAccess::reset();
}

void run_open_boundary_composition(const hundun::runtime::MpiContext &mpi) {
  constexpr double rho_ref = 1.25;
  constexpr double pressure_reference = -2.5;
  constexpr double transport_diffusivity = 0.01;
  constexpr hundun::runtime::Int3 extent{8, 6, 4};
  auto decomposition = hundun::runtime::StructuredDecomposition::create(
      mpi, extent, {false, false, false},
      hundun::runtime::DecompositionOptions{process_grid(mpi.size())});
  const hundun::mesh::MeshTopology topology(decomposition);
  const hundun::mesh::MeshGeometry geometry(
      topology,
      hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
  auto boundaries = hundun::boundary::BoundaryRegistry::create(
      open_case(rho_ref, pressure_reference), topology);

  hundun::runtime::FieldRegistry registry;
  hundun::flow::FlowFieldIds fields;
  fields.density = registry.declare_field(cell_field("open_rho", 1U));
  fields.velocity = registry.declare_field(cell_field("open_velocity", 3U));
  fields.mechanical_pressure =
      registry.declare_field(cell_field("open_pi", 1U));
  fields.face_velocity =
      registry.declare_field(face_field("open_u_face", 3U));
  fields.face_mass_flux =
      hundun::finite_volume::declare_face_mass_flux(registry);
  fields.transported_cell_fields.push_back(
      registry.declare_field(cell_field("open_alpha", 1U)));
  registry.freeze();

  const auto box = decomposition.owned_box();
  const hundun::runtime::Int3 local{box.end.x - box.begin.x,
                                    box.end.y - box.begin.y,
                                    box.end.z - box.begin.z};
  const std::size_t cells = topology.owned_cell_count();
  hundun::flow::FlowLayerValues initial;
  initial.density.assign(cells, rho_ref);
  initial.velocity.assign(cells * 3U, 0.0);
  initial.mechanical_pressure.assign(cells, pressure_reference);
  initial.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
  initial.face_mass_flux.assign(topology.local_face_count(), 0.0);
  initial.transported_cell_fields = {std::vector<double>(cells, 1.0)};

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
  using TestAccess = hundun::flow::test::ConstantDensityPisoTestAccess;

  TestAccess::reset();
  auto reference_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0,
       hundun::flow::MomentumTimeOrder::backward_euler});
  reference_state.seed_accepted_layers(initial, initial);
  const auto reference_report =
      flow.attempt(reference_state, rho_ref, 0.0, stencil, {}, {});
  HUNDUN_CHECK(reference_report.disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  HUNDUN_CHECK(TestAccess::last_pressure_constraint_mode() ==
               static_cast<int>(hundun::finite_volume::PressureConstraintMode::
                                    pressure_reference_patch));
  HUNDUN_CHECK(TestAccess::last_pressure_operator_mode() ==
               static_cast<int>(hundun::finite_volume::PressureConstraintMode::
                                    pressure_reference_patch));
  const auto reference_result =
      reference_state.snapshot(hundun::flow::FlowLayer::committed);
  for (double pressure : reference_result.mechanical_pressure) {
    HUNDUN_CHECK_NEAR(pressure, pressure_reference, 0.0);
  }

  auto advective_boundaries = hundun::boundary::BoundaryRegistry::create(
      open_case(rho_ref, pressure_reference, 0.1), topology);
  hundun::linear::JacobiPreconditioner advective_mx(execution);
  hundun::linear::JacobiPreconditioner advective_my(execution);
  hundun::linear::JacobiPreconditioner advective_mz(execution);
  hundun::linear::JacobiPreconditioner advective_pressure(execution);
  auto advective_flow = hundun::flow::FixedStepConstantDensityFlow::create(
      decomposition, topology, geometry, advective_boundaries, mpi, execution,
      halo, momentum_solver, {&advective_mx, &advective_my, &advective_mz},
      pressure_solver, advective_pressure,
      {{fields.transported_cell_fields[0],
        hundun::finite_volume::FiniteVolumeQuantity::scalar(0U),
        transport_diffusivity}});
  auto advective_initial = initial;
  for (hundun::mesh::LocalCellId cell = 0;
       cell < topology.owned_cell_count(); ++cell) {
    const double x = geometry.cell_center_m(cell).x;
    advective_initial.velocity[cell * 3U] = 0.1 + 0.04 * x;
    advective_initial.transported_cell_fields[0][cell] = 3.0 + 0.3 * x;
  }
  for (hundun::mesh::LocalFaceId face = 0;
       face < topology.local_face_count(); ++face) {
    const auto area = geometry.face_area_vector_m2(
        face, hundun::mesh::FaceSide::owner);
    advective_initial.face_mass_flux[face] = rho_ref * 0.1 * area.x;
    advective_initial.face_velocity[face * 3U] = 0.1;
  }
  auto advective_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0,
       hundun::flow::MomentumTimeOrder::backward_euler});
  advective_state.seed_accepted_layers(advective_initial, advective_initial);
  TestAccess::reset();
  TestAccess::set_final_uniform_x_mass_flux(0.1 * rho_ref);
  TestAccess::set_final_uniform_x_mass_flux_override(true);
  TestAccess::set_final_momentum_norm_squares(1U, 0.0, 1.0);
  TestAccess::set_final_momentum_norm_squares(2U, 0.0, 1.0);
  const auto advective_be_report = advective_flow.attempt(
      advective_state, rho_ref, 0.0, stencil, {}, {});
  if (advective_be_report.disposition !=
          hundun::flow::StepAttemptDisposition::committed &&
      mpi.rank() == 0) {
    std::cerr << "TASK18_OPEN_BE reason="
              << static_cast<int>(advective_be_report.reason)
              << " momentum="
              << advective_be_report.final_momentum_normalized_l2[0]
              << ',' << advective_be_report.final_momentum_normalized_l2[1]
              << ',' << advective_be_report.final_momentum_normalized_l2[2]
              << " transport="
              << (advective_be_report.final_transport_normalized_l2.empty()
                      ? -1.0
                      : advective_be_report.final_transport_normalized_l2[0])
              << " conservation="
              << advective_be_report
                     .final_momentum_relative_conservation_defect[0]
              << ','
              << (advective_be_report
                          .final_transport_relative_conservation_defect.empty()
                      ? -1.0
                      : advective_be_report
                            .final_transport_relative_conservation_defect[0])
              << '\n';
  }
  HUNDUN_CHECK(advective_be_report.disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  const auto advective_after_be =
      advective_state.snapshot(hundun::flow::FlowLayer::committed);
  const auto be_momentum_diagnostic =
      TestAccess::last_momentum_conservation(0U);
  const auto be_transport_diagnostic =
      TestAccess::last_transport_conservation(0U);

  const auto global_quantities =
      [&](const hundun::flow::FlowLayerValues &previous,
          const hundun::flow::FlowLayerValues &current,
          const hundun::flow::FlowLayerValues &next) {
        std::array<double, 12> values{};
        for (hundun::mesh::LocalCellId cell = 0;
             cell < topology.owned_cell_count(); ++cell) {
          const double mass = rho_ref * geometry.cell_volume_m3(cell);
          for (std::size_t component = 0; component < 3U; ++component) {
            values[component * 3U] +=
                mass * previous.velocity[cell * 3U + component];
            values[component * 3U + 1U] +=
                mass * current.velocity[cell * 3U + component];
            values[component * 3U + 2U] +=
                mass * next.velocity[cell * 3U + component];
          }
          values[9] += mass * previous.transported_cell_fields[0][cell];
          values[10] += mass * current.transported_cell_fields[0][cell];
          values[11] += mass * next.transported_cell_fields[0][cell];
        }
        mpi.allreduce_fp64_in_place(
            values.data(), values.size(),
            hundun::runtime::Fp64ReductionOperation::sum);
        return values;
      };
  const auto global_boundary_contributions =
      [&](const hundun::flow::FlowLayerValues &previous,
          const hundun::flow::FlowLayerValues &current,
          const hundun::flow::FlowLayerValues &pressure_state, bool bdf2) {
        std::array<double, 14> values{};
        for (hundun::mesh::LocalFaceId face = 0;
             face < topology.local_face_count(); ++face) {
          const auto patch = topology.patch_id(face);
          if (!patch.has_value() ||
              advective_boundaries.patch(*patch).kind() ==
                  hundun::boundary::BoundaryKind::periodic) {
            continue;
          }
          const auto owner = topology.owner(face);
          HUNDUN_CHECK(owner < topology.owned_cell_count());
          const auto area = geometry.face_area_vector_m2(
              face, hundun::mesh::FaceSide::owner);
          const double mass_flux = rho_ref * 0.1 * area.x;
          values[0] += mass_flux;
          values[1] += std::abs(mass_flux);
          const hundun::runtime::Real3 velocity_n{
              current.velocity[owner * 3U],
              current.velocity[owner * 3U + 1U],
              current.velocity[owner * 3U + 2U]};
          const hundun::runtime::Real3 velocity_nm1{
              previous.velocity[owner * 3U],
              previous.velocity[owner * 3U + 1U],
              previous.velocity[owner * 3U + 2U]};
          const auto face_velocity_n = advective_boundaries
                                           .evaluate_velocity(*patch,
                                                              velocity_n, area)
                                           .face;
          const auto face_velocity_nm1 =
              advective_boundaries
                  .evaluate_velocity(*patch, velocity_nm1, area)
                  .face;
          const double pressure_face =
              advective_boundaries
                  .evaluate_pressure(
                      *patch, pressure_state.mechanical_pressure[owner])
                  .face;
          const std::array<double, 3> area_components{area.x, area.y, area.z};
          const std::array<double, 3> velocity_components_n{
              face_velocity_n.x, face_velocity_n.y, face_velocity_n.z};
          const std::array<double, 3> velocity_components_nm1{
              face_velocity_nm1.x, face_velocity_nm1.y,
              face_velocity_nm1.z};
          for (std::size_t component = 0; component < 3U; ++component) {
            const double convection_n =
                mass_flux * velocity_components_n[component];
            const double convection_nm1 =
                mass_flux * velocity_components_nm1[component];
            const double convection =
                bdf2 ? 2.0 * convection_n - convection_nm1 : convection_n;
            const double pressure =
                pressure_face * area_components[component];
            const double effective = convection + pressure;
            values[2U + component] += effective;
            values[5U + component] += std::abs(effective);
            values[10U + component] +=
                std::abs(convection) + std::abs(pressure);
          }
          const double scalar_n =
              advective_boundaries
                  .evaluate_scalar(
                      *patch, 0U,
                      current.transported_cell_fields[0][owner])
                  .face;
          const double scalar_nm1 =
              advective_boundaries
                  .evaluate_scalar(
                      *patch, 0U,
                      previous.transported_cell_fields[0][owner])
                  .face;
          const auto scalar_diffusion = [&](double owner_value) {
            const auto evaluated = advective_boundaries.evaluate_scalar(
                *patch, 0U, owner_value);
            if (advective_boundaries.patch(*patch).scalar_rule() !=
                hundun::boundary::TransportRule::prescribed_value) {
              return 0.0;
            }
            const auto owner_center = geometry.cell_center_m(owner);
            const auto face_center = geometry.face_center_m(face);
            const hundun::runtime::Real3 d{
                2.0 * (face_center.x - owner_center.x),
                2.0 * (face_center.y - owner_center.y),
                2.0 * (face_center.z - owner_center.z)};
            const double area_square =
                area.x * area.x + area.y * area.y + area.z * area.z;
            const double area_dot_d =
                area.x * d.x + area.y * d.y + area.z * d.z;
            HUNDUN_CHECK(area_dot_d > 0.0);
            const double physical_flux =
                transport_diffusivity * (evaluated.exterior - owner_value) *
                area_square / area_dot_d;
            return -physical_flux;
          };
          const double convection_n = mass_flux * scalar_n;
          const double convection_nm1 = mass_flux * scalar_nm1;
          const double diffusion_n = scalar_diffusion(
              current.transported_cell_fields[0][owner]);
          const double diffusion_nm1 = scalar_diffusion(
              previous.transported_cell_fields[0][owner]);
          const double transport_convection =
              bdf2 ? 2.0 * convection_n - convection_nm1 : convection_n;
          const double transport_diffusion =
              bdf2 ? 2.0 * diffusion_n - diffusion_nm1 : diffusion_n;
          const double transport =
              transport_convection + transport_diffusion;
          values[8] += transport;
          values[9] += std::abs(transport);
          values[13] +=
              std::abs(transport_convection) + std::abs(transport_diffusion);
        }
        mpi.allreduce_fp64_in_place(
            values.data(), values.size(),
            hundun::runtime::Fp64ReductionOperation::sum);
        return values;
      };
  const auto check_diagnostic =
      [](const hundun::flow::test::ConservationDiagnostic &diagnostic,
         double quantity_nm1, double quantity_n, double quantity_np1,
         double signed_boundary, double absolute_boundary,
         const hundun::flow::MomentumTimeStencil &time_stencil) {
        const double history =
            time_stencil.alpha2 * (quantity_nm1 - quantity_n) /
            time_stencil.alpha0;
        const double boundary =
            time_stencil.dt_s * signed_boundary / time_stencil.alpha0;
        const double absolute_scale =
            time_stencil.dt_s * absolute_boundary / time_stencil.alpha0;
        const double raw = quantity_np1 - quantity_n + boundary + history;
        const double denominator = std::max(
            {std::abs(quantity_n), std::abs(quantity_np1), absolute_scale,
             std::abs(history), std::numeric_limits<double>::min()});
        const double tolerance =
            512.0 * std::numeric_limits<double>::epsilon() *
            std::max({1.0, std::abs(quantity_nm1), std::abs(quantity_n),
                      std::abs(quantity_np1), absolute_scale});
        HUNDUN_CHECK_NEAR(diagnostic.quantity_nm1, quantity_nm1, tolerance);
        HUNDUN_CHECK_NEAR(diagnostic.quantity_n, quantity_n, tolerance);
        HUNDUN_CHECK_NEAR(diagnostic.quantity_np1, quantity_np1, tolerance);
        HUNDUN_CHECK_NEAR(diagnostic.signed_boundary_flux, signed_boundary,
                          tolerance);
        HUNDUN_CHECK_NEAR(diagnostic.absolute_boundary_flux,
                          absolute_boundary, tolerance);
        HUNDUN_CHECK_NEAR(diagnostic.boundary_integral, boundary, tolerance);
        HUNDUN_CHECK_NEAR(diagnostic.history_correction, history, tolerance);
        HUNDUN_CHECK_NEAR(diagnostic.raw_defect, raw, tolerance);
        HUNDUN_CHECK_NEAR(diagnostic.relative_defect,
                          std::abs(raw) / denominator, tolerance);
      };
  const auto be_quantities = global_quantities(
      advective_initial, advective_initial, advective_after_be);
  const auto be_boundaries = global_boundary_contributions(
      advective_initial, advective_initial, advective_after_be, false);
  HUNDUN_CHECK(be_boundaries[5] < be_boundaries[10]);
  HUNDUN_CHECK(be_boundaries[9] < be_boundaries[13]);
  if (mpi.rank() == 0) {
    std::cout << "TASK18_R3_FACE_COMPOSITION momentum_qnm1="
              << be_momentum_diagnostic.quantity_nm1
              << " expected=" << be_quantities[0]
              << " momentum_abs="
              << be_momentum_diagnostic.absolute_boundary_flux
              << " expected_abs=" << be_boundaries[5]
              << " termwise_abs=" << be_boundaries[10]
              << " transport_abs="
              << be_transport_diagnostic.absolute_boundary_flux
              << " expected_abs=" << be_boundaries[9]
              << " termwise_abs=" << be_boundaries[13] << '\n';
  }
  check_diagnostic(be_momentum_diagnostic, be_quantities[0], be_quantities[1],
                   be_quantities[2], be_boundaries[2], be_boundaries[5],
                   stencil);
  check_diagnostic(be_transport_diagnostic, be_quantities[9],
                   be_quantities[10], be_quantities[11], be_boundaries[8],
                   be_boundaries[9], stencil);
  HUNDUN_CHECK(std::abs(be_momentum_diagnostic.quantity_n) > 0.0);
  HUNDUN_CHECK(std::abs(be_momentum_diagnostic.signed_boundary_flux) > 0.0);
  HUNDUN_CHECK(std::abs(be_transport_diagnostic.quantity_n) > 0.0);
  HUNDUN_CHECK(std::abs(be_transport_diagnostic.signed_boundary_flux) > 0.0);
  const auto advective_bdf2_stencil =
      hundun::flow::make_momentum_time_stencil(
          hundun::flow::MomentumTimeOrder::bdf2, 0.01, 0.01);
  const auto advective_bdf2_report = advective_flow.attempt(
      advective_state, rho_ref, 0.0, advective_bdf2_stencil, {}, {});
  if (advective_bdf2_report.disposition !=
          hundun::flow::StepAttemptDisposition::committed &&
      mpi.rank() == 0) {
    std::cerr << "TASK18_OPEN_BDF2 reason="
              << static_cast<int>(advective_bdf2_report.reason)
              << " momentum="
              << advective_bdf2_report.final_momentum_normalized_l2[0]
              << ',' << advective_bdf2_report.final_momentum_normalized_l2[1]
              << ',' << advective_bdf2_report.final_momentum_normalized_l2[2]
              << " transport="
              << (advective_bdf2_report.final_transport_normalized_l2.empty()
                      ? -1.0
                      : advective_bdf2_report
                            .final_transport_normalized_l2[0])
              << " conservation="
              << advective_bdf2_report
                     .final_momentum_relative_conservation_defect[0]
              << ','
              << (advective_bdf2_report
                          .final_transport_relative_conservation_defect.empty()
                      ? -1.0
                      : advective_bdf2_report
                            .final_transport_relative_conservation_defect[0])
              << '\n';
  }
  HUNDUN_CHECK(advective_bdf2_report.disposition ==
               hundun::flow::StepAttemptDisposition::committed);
  const auto advective_after_bdf2 =
      advective_state.snapshot(hundun::flow::FlowLayer::committed);
  const auto bdf2_momentum_diagnostic =
      TestAccess::last_momentum_conservation(0U);
  const auto bdf2_transport_diagnostic =
      TestAccess::last_transport_conservation(0U);
  const auto bdf2_quantities = global_quantities(
      advective_initial, advective_after_be, advective_after_bdf2);
  const auto bdf2_boundaries = global_boundary_contributions(
      advective_initial, advective_after_be, advective_after_bdf2, true);
  check_diagnostic(bdf2_momentum_diagnostic, bdf2_quantities[0],
                   bdf2_quantities[1], bdf2_quantities[2],
                   bdf2_boundaries[2], bdf2_boundaries[5],
                   advective_bdf2_stencil);
  check_diagnostic(bdf2_transport_diagnostic, bdf2_quantities[9],
                   bdf2_quantities[10], bdf2_quantities[11],
                   bdf2_boundaries[8], bdf2_boundaries[9],
                   advective_bdf2_stencil);
  HUNDUN_CHECK(std::abs(bdf2_momentum_diagnostic.history_correction) > 0.0);
  HUNDUN_CHECK(std::abs(bdf2_transport_diagnostic.history_correction) > 0.0);
  if (mpi.rank() == 0) {
    std::cout << "TASK18_OPEN_CONSERVATION be_momentum_flux="
              << be_momentum_diagnostic.signed_boundary_flux
              << " be_transport_flux="
              << be_transport_diagnostic.signed_boundary_flux
              << " bdf2_momentum_flux="
              << bdf2_momentum_diagnostic.signed_boundary_flux
              << " bdf2_transport_flux="
              << bdf2_transport_diagnostic.signed_boundary_flux << '\n';
  }
  TestAccess::reset();

  TestAccess::set_final_uniform_x_mass_flux(-0.1);
  TestAccess::set_final_uniform_x_mass_flux_override(true);
  auto backflow_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      {0U, 0.0, 0.01, 0.0,
       hundun::flow::MomentumTimeOrder::backward_euler});
  backflow_state.seed_accepted_layers(initial, initial);
  const auto backflow_before =
      backflow_state.snapshot(hundun::flow::FlowLayer::committed);
  const auto backflow_report =
      flow.attempt(backflow_state, rho_ref, 0.0, stencil, {}, {});
  HUNDUN_CHECK(backflow_report.disposition ==
               hundun::flow::StepAttemptDisposition::recoverable_failure);
  HUNDUN_CHECK(backflow_report.reason ==
               hundun::flow::StepFailureReason::boundary_backflow);
  HUNDUN_CHECK(backflow_report.pressure_corrector_count == 2U);
  HUNDUN_CHECK(backflow_report.suggested_dt_s == 0.005);
  HUNDUN_CHECK(TestAccess::final_transport_calls() == 1U);
  HUNDUN_CHECK(backflow_report.final_backflow_evidence.has_value());
  HUNDUN_CHECK(backflow_report.final_backflow_evidence->patch_id == 1U);
  HUNDUN_CHECK(backflow_report.final_backflow_evidence->step == 1U);
  HUNDUN_CHECK(backflow_report.final_backflow_evidence->time_s == 0.01);
  bool has_local_outlet = false;
  double expected_minimum = 0.0;
  hundun::mesh::GlobalFaceId local_minimum_face =
      std::numeric_limits<hundun::mesh::GlobalFaceId>::max();
  for (hundun::mesh::LocalFaceId face = 0;
       face < topology.local_face_count(); ++face) {
    if (topology.patch_id(face) != std::optional<std::uint32_t>{1U})
      continue;
    has_local_outlet = true;
    const double injected =
        -0.1 * geometry.face_area_vector_m2(
                   face, hundun::mesh::FaceSide::owner)
                   .x;
    const auto global_face = topology.global_face_id(face);
    if (injected < expected_minimum ||
        (injected == expected_minimum && injected < 0.0 &&
         global_face < local_minimum_face)) {
      expected_minimum = injected;
      local_minimum_face = global_face;
    }
  }
  const double local_expected_minimum = expected_minimum;
  double evidence_expectations[]{-local_expected_minimum,
                                 has_local_outlet
                                     ? -static_cast<double>(mpi.rank())
                                     : -static_cast<double>(mpi.size())};
  mpi.allreduce_fp64_in_place(
      evidence_expectations, 2U,
      hundun::runtime::Fp64ReductionOperation::maximum);
  expected_minimum = -evidence_expectations[0];
  double negative_lowest_face =
      local_expected_minimum == expected_minimum && expected_minimum < 0.0 &&
              local_minimum_face !=
                  std::numeric_limits<hundun::mesh::GlobalFaceId>::max()
          ? -static_cast<double>(local_minimum_face)
          : -std::numeric_limits<double>::max();
  mpi.allreduce_fp64_in_place(
      &negative_lowest_face, 1U,
      hundun::runtime::Fp64ReductionOperation::maximum);
  const auto expected_face = static_cast<hundun::mesh::GlobalFaceId>(
      -negative_lowest_face);
  const int expected_rank = static_cast<int>(-evidence_expectations[1]);
  double boundary_sums[2]{};
  for (hundun::mesh::LocalFaceId face = 0;
       face < topology.local_face_count(); ++face) {
    const auto patch = topology.patch_id(face);
    if (!patch.has_value() ||
        boundaries.patch(*patch).kind() ==
            hundun::boundary::BoundaryKind::periodic) {
      continue;
    }
    const double injected =
        -0.1 * geometry.face_area_vector_m2(
                   face, hundun::mesh::FaceSide::owner)
                   .x;
    boundary_sums[0] += injected;
    boundary_sums[1] += std::abs(injected);
  }
  mpi.allreduce_fp64_in_place(
      boundary_sums, 2U, hundun::runtime::Fp64ReductionOperation::sum);
  const auto mass_diagnostic = TestAccess::last_mass_conservation();
  const double boundary_tolerance =
      64.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, boundary_sums[1]);
  HUNDUN_CHECK_NEAR(mass_diagnostic.signed_boundary_flux, boundary_sums[0],
                    boundary_tolerance);
  HUNDUN_CHECK_NEAR(mass_diagnostic.absolute_boundary_flux, boundary_sums[1],
                    boundary_tolerance);
  HUNDUN_CHECK_NEAR(mass_diagnostic.boundary_integral,
                    stencil.dt_s * boundary_sums[0], boundary_tolerance);
  HUNDUN_CHECK_NEAR(mass_diagnostic.history_correction, 0.0, 0.0);
  HUNDUN_CHECK(mass_diagnostic.absolute_boundary_flux > 0.0);
  if (mpi.rank() == 0) {
    std::cout << "TASK18_BACKFLOW min="
              << backflow_report.final_backflow_evidence
                     ->minimum_outward_mass_flux_kg_per_s
              << " expected_min=" << expected_minimum << " face="
              << backflow_report.final_backflow_evidence->global_face_id
              << " expected_face=" << expected_face << " rank="
              << backflow_report.final_backflow_evidence->lowest_failing_rank
              << " expected_rank=" << expected_rank << '\n';
  }
  HUNDUN_CHECK_NEAR(
      backflow_report.final_backflow_evidence
          ->minimum_outward_mass_flux_kg_per_s,
      expected_minimum, 0.0);
  HUNDUN_CHECK(backflow_report.final_backflow_evidence->global_face_id ==
               expected_face);
  HUNDUN_CHECK(backflow_report.final_backflow_evidence->lowest_failing_rank ==
               expected_rank);
  HUNDUN_CHECK(backflow_report.lowest_failing_rank ==
               backflow_report.final_backflow_evidence->lowest_failing_rank);
  check_layer_equal(
      backflow_state.snapshot(hundun::flow::FlowLayer::committed),
      backflow_before);
  HUNDUN_CHECK(backflow_state.metadata().step == 0U);
  TestAccess::reset();
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] {
    check_flow_layer_bitwise_oracle();
    run_zero_flow_transaction(mpi);
    run_open_boundary_composition(mpi);
  });
}

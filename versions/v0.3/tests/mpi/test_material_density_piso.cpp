// SPDX-License-Identifier: Apache-2.0

#include "tests/support/flow_material_density_piso_test_access.hpp"
#include "tests/support/flow_material_density_transport_test_access.hpp"
#include "hundun/bc_basic_boundary.hpp"
#include "hundun/cfg_resolved_case.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/flow_material_density_piso.hpp"
#include "hundun/lin_bicgstab.hpp"
#include "hundun/lin_conjugate_gradient.hpp"
#include "hundun/lin_preconditioners.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_mpi_operation_error.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/task25_counter_checks.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <new>
#include <utility>

namespace task20_allocation_oracle {

std::atomic<const hundun::flow::FlowState *> observed_state{};
std::atomic<std::uint64_t> active_attempt_allocations{};

void record_active_attempt_allocation() noexcept {
  const auto *state = observed_state.load(std::memory_order_relaxed);
  if (state != nullptr &&
      hundun::flow::test::MaterialDensityPisoTestAccess::state_attempt_active(
          *state) &&
      hundun::flow::test::MaterialDensityPisoTestAccess::
          face_flux_path_observation_active())
    active_attempt_allocations.fetch_add(1U, std::memory_order_relaxed);
}

void *allocate(std::size_t size) {
  record_active_attempt_allocation();
  if (void *result = std::malloc(size == 0U ? 1U : size))
    return result;
  throw std::bad_alloc();
}

void *allocate_aligned(std::size_t size, std::size_t alignment) {
  record_active_attempt_allocation();
  void *result = nullptr;
  if (posix_memalign(&result, alignment, size == 0U ? 1U : size) == 0)
    return result;
  throw std::bad_alloc();
}

class Observation final {
public:
  explicit Observation(const hundun::flow::FlowState &state) {
    active_attempt_allocations.store(0U, std::memory_order_relaxed);
    observed_state.store(&state, std::memory_order_relaxed);
  }
  ~Observation() { observed_state.store(nullptr, std::memory_order_relaxed); }
  Observation(const Observation &) = delete;
  Observation &operator=(const Observation &) = delete;

  std::uint64_t count() const noexcept {
    return active_attempt_allocations.load(std::memory_order_relaxed);
  }
};

} // namespace task20_allocation_oracle

void *operator new(std::size_t size) {
  return task20_allocation_oracle::allocate(size);
}
void *operator new[](std::size_t size) {
  return task20_allocation_oracle::allocate(size);
}
void operator delete(void *pointer) noexcept { std::free(pointer); }
void operator delete[](void *pointer) noexcept { std::free(pointer); }
void operator delete(void *pointer, std::size_t) noexcept {
  std::free(pointer);
}
void operator delete[](void *pointer, std::size_t) noexcept {
  std::free(pointer);
}
void *operator new(std::size_t size, std::align_val_t alignment) {
  return task20_allocation_oracle::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}
void *operator new[](std::size_t size, std::align_val_t alignment) {
  return task20_allocation_oracle::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}
void operator delete(void *pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void *pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete(void *pointer, std::size_t, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void *pointer, std::size_t, std::align_val_t) noexcept {
  std::free(pointer);
}

namespace {

hundun::runtime::Int3 grid(int ranks) {
  if (ranks == 1)
    return {1, 1, 1};
  if (ranks == 2)
    return {2, 1, 1};
  if (ranks == 4)
    return {2, 2, 1};
  throw hundun::runtime::Error("unsupported Task 20 rank count");
}

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
  for (std::size_t patch = 0; patch < names.size(); ++patch) {
    config.boundaries[patch].patch = names[patch];
    config.boundaries[patch].type = hundun::config::BoundaryType::periodic;
  }
  return config;
}

hundun::config::FlowCaseConfig wall_case() {
  auto config = periodic_case();
  config.boundaries[0].type =
      hundun::config::BoundaryType::no_slip_wall;
  config.boundaries[1].type =
      hundun::config::BoundaryType::no_slip_wall;
  return config;
}

hundun::runtime::FieldDescriptor cell(const char *name, const char *unit,
                                      std::uint32_t components,
                                      bool conservative) {
  return {name,
          unit,
          "task20",
          hundun::runtime::FunctionSpace::cell_average,
          hundun::runtime::ScalarType::float64,
          components,
          2,
          conservative,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

hundun::runtime::FieldDescriptor face(const char *name, const char *unit,
                                      std::uint32_t components) {
  return {name,
          unit,
          "task20",
          hundun::runtime::FunctionSpace::face_value,
          hundun::runtime::ScalarType::float64,
          components,
          0,
          false,
          hundun::runtime::RestartPolicy::persistent,
          hundun::runtime::OutputPolicy::never};
}

struct ExactState final {
  hundun::flow::FlowLayerValues history;
  hundun::flow::FlowLayerValues committed;
  hundun::flow::FlowLayerValues trial;
  hundun::flow::AcceptedStepMetadata metadata;
};

std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

class TerminalMomentumSolver final : public hundun::linear::LinearSolver {
public:
  enum class Mode {
    rankless_result,
    thrown_operation,
    reliable_result,
    maximum_iterations,
    numerical_breakdown,
    non_finite_value
  };
  explicit TerminalMomentumSolver(Mode mode, std::size_t fail_on_call = 1U)
      : mode_(mode), fail_on_call_(fail_on_call) {}

  hundun::linear::SolveReport solve(
      const hundun::linear::LinearOperator &,
      hundun::linear::Preconditioner &,
      hundun::execution::VectorView<const double>,
      hundun::execution::VectorView<double> x,
      const hundun::linear::SolveControl &) const override {
    ++calls_;
    if (calls_ < fail_on_call_) {
      for (std::size_t index = 0; index < x.size(); ++index)
        x[index] = 0.0;
      hundun::linear::SolveReport success;
      success.reason =
          hundun::linear::SolveTerminationReason::zero_right_hand_side;
      success.initial_residual = 0.0;
      success.recursive_residual = 0.0;
      success.final_residual = 0.0;
      return success;
    }
    if (mode_ == Mode::thrown_operation)
      throw hundun::runtime::MpiOperationError(
          "injected material momentum operation failure");
    hundun::linear::SolveReport report;
    if (mode_ == Mode::maximum_iterations)
      report.reason =
          hundun::linear::SolveTerminationReason::maximum_iterations;
    else if (mode_ == Mode::numerical_breakdown)
      report.reason =
          hundun::linear::SolveTerminationReason::numerical_breakdown;
    else if (mode_ == Mode::non_finite_value)
      report.reason =
          hundun::linear::SolveTerminationReason::non_finite_value;
    else {
      report.reason =
          hundun::linear::SolveTerminationReason::collective_failure;
      report.lowest_failing_rank =
          mode_ == Mode::reliable_result ? 0 : -1;
    }
    return report;
  }

  std::size_t calls() const noexcept { return calls_; }

private:
  Mode mode_;
  std::size_t fail_on_call_;
  mutable std::size_t calls_{};
};

ExactState capture(const hundun::flow::FlowState &state) {
  return {state.snapshot(hundun::flow::FlowLayer::history),
          state.snapshot(hundun::flow::FlowLayer::committed),
          state.snapshot(hundun::flow::FlowLayer::trial), state.metadata()};
}

void check_equal(const ExactState &expected,
                 const hundun::flow::FlowState &state) {
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      expected.history, state.snapshot(hundun::flow::FlowLayer::history)));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      expected.committed, state.snapshot(hundun::flow::FlowLayer::committed)));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      expected.trial, state.snapshot(hundun::flow::FlowLayer::trial)));
  HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
      expected.metadata, state.metadata()));
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
    {
      const hundun::flow::test::MaterialMomentumConservationInput input{
          -2.0, 1.0, 4.0, -0.25, 0.5, 0.75,
          2.0, 1.0, 4.0, 0.25, 0.5, 0.75,
          0.2, 1.5, 0.5, true};
      const double history = (input.alpha2 / input.alpha0) *
                             (input.momentum_n_minus_1 - input.momentum_n);
      const double boundary = (input.dt_s / input.alpha0) *
                              (2.0 * input.boundary_n -
                               input.boundary_n_minus_1);
      const double source = (input.dt_s / input.alpha0) *
                            input.source_n_plus_1;
      const double raw = input.momentum_n_plus_1 - input.momentum_n +
                         history + boundary - source;
      const double boundary_abs =
          2.0 * input.boundary_abs_n + input.boundary_abs_n_minus_1;
      double denominator = std::max(
          {std::abs(input.momentum_n), std::abs(input.momentum_n_plus_1),
           std::abs(history),
           (input.dt_s / std::abs(input.alpha0)) * boundary_abs,
           (input.dt_s / std::abs(input.alpha0)) *
               input.source_abs_n_plus_1,
           std::numeric_limits<double>::min()});
      const double cancellation =
          std::max({input.momentum_abs_n_minus_1, input.momentum_abs_n,
                    input.momentum_abs_n_plus_1});
      if (denominator <= 64.0 * std::numeric_limits<double>::epsilon() *
                             cancellation)
        denominator = std::max(denominator, cancellation);
      const double expected = std::abs(raw) / denominator;
      const double actual = hundun::flow::test::MaterialDensityPisoTestAccess::
          momentum_conservation_defect(input);
      HUNDUN_CHECK(actual == expected);
      auto wrong_source = input;
      wrong_source.source_n_plus_1 = -wrong_source.source_n_plus_1;
      HUNDUN_CHECK(hundun::flow::test::MaterialDensityPisoTestAccess::
                       momentum_conservation_defect(wrong_source) != actual);
      auto cancelling = input;
      cancelling.momentum_n_minus_1 = 0.0;
      cancelling.momentum_n = 0.0;
      cancelling.momentum_n_plus_1 = 0.0;
      cancelling.boundary_n_minus_1 = 0.0;
      cancelling.boundary_n = 0.0;
      cancelling.source_n_plus_1 = 0.0;
      cancelling.momentum_abs_n_minus_1 = 1.0e12;
      cancelling.momentum_abs_n = 1.0e12;
      cancelling.momentum_abs_n_plus_1 = 1.0e12;
      cancelling.source_abs_n_plus_1 = 1.0e12;
      HUNDUN_CHECK(hundun::flow::test::MaterialDensityPisoTestAccess::
                       momentum_conservation_defect(cancelling) == 0.0);
    }
    {
      using hundun::flow::StepFailureReason;
      using Access = hundun::flow::test::MaterialDensityPisoTestAccess;
      const auto select = [&](bool low_rank_uses_recoverable_reason) {
        const bool low = mpi.rank() == 0;
        const bool high = mpi.rank() == mpi.size() - 1 && mpi.size() > 1;
        hundun::flow::test::MaterialPhaseFailureForTest local{};
        if (low) {
          local = {true,
                   static_cast<std::uint8_t>(
                       low_rank_uses_recoverable_reason
                           ? StepFailureReason::non_finite_trial
                           : StepFailureReason::invalid_input),
                   low_rank_uses_recoverable_reason};
        } else if (high) {
          local = {true,
                   static_cast<std::uint8_t>(
                       low_rank_uses_recoverable_reason
                           ? StepFailureReason::invalid_input
                           : StepFailureReason::non_finite_trial),
                   !low_rank_uses_recoverable_reason};
        }
        return Access::select_phase_failure(mpi, local);
      };
      const auto first = select(true);
      HUNDUN_CHECK(first.lowest_failing_rank == 0);
      HUNDUN_CHECK(first.reason == static_cast<std::uint8_t>(
                                       StepFailureReason::non_finite_trial));
      HUNDUN_CHECK(first.recoverable);
      const auto second_selection = select(false);
      HUNDUN_CHECK(second_selection.lowest_failing_rank == 0);
      HUNDUN_CHECK(second_selection.reason == static_cast<std::uint8_t>(
                                                  StepFailureReason::invalid_input));
      HUNDUN_CHECK(!second_selection.recoverable);

      bool exact_terminal = false;
      try {
        Access::require_reliable_collective_result(
            static_cast<std::uint8_t>(StepFailureReason::collective_operation),
            -1);
      } catch (const hundun::runtime::MpiOperationError &error) {
        exact_terminal =
            std::string_view(error.what()) ==
            "material-density flow collective failure has no reliable rank";
      }
      HUNDUN_CHECK(exact_terminal);
      Access::require_reliable_collective_result(
          static_cast<std::uint8_t>(StepFailureReason::collective_operation),
          0);
      Access::require_reliable_collective_result(
          static_cast<std::uint8_t>(StepFailureReason::invalid_input), -1);
    }
    constexpr hundun::runtime::Int3 extent{8, 4, 4};
    auto decomposition = hundun::runtime::StructuredDecomposition::create(
        mpi, extent, {true, true, true},
        hundun::runtime::DecompositionOptions{grid(mpi.size())});
    hundun::mesh::MeshTopology topology(decomposition);
    hundun::mesh::MeshGeometry geometry(
        topology,
        hundun::mesh::UniformBoxMapping({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
    auto boundaries = hundun::boundary::BoundaryRegistry::create(
        periodic_case(), topology);

    hundun::runtime::FieldRegistry registry;
    hundun::flow::FlowFieldIds fields;
    fields.density = registry.declare_field(cell("rho", "kg/m3", 1U, true));
    fields.velocity =
        registry.declare_field(cell("velocity", "m/s", 3U, false));
    fields.mechanical_pressure =
        registry.declare_field(cell("pi", "Pa", 1U, false));
    fields.face_velocity =
        registry.declare_field(face("face_velocity", "m/s", 3U));
    fields.face_mass_flux =
        hundun::finite_volume::declare_face_mass_flux(registry);
    const auto rho_h =
        registry.declare_field(cell("rho_h", "J/m3", 1U, true));
    const auto alternate_density =
        registry.declare_field(cell("rho_alternate", "kg/m3", 1U, true));
    fields.transported_cell_fields = {rho_h};
    registry.freeze();

    const auto box = decomposition.owned_box();
    const hundun::runtime::Int3 local{box.end.x - box.begin.x,
                                      box.end.y - box.begin.y,
                                      box.end.z - box.begin.z};
    const auto make_state = [&](double density) {
      auto state = hundun::flow::FlowState::create(
          registry, {local, topology.local_face_count()}, fields,
          {0U, 0.0, 0.01, 0.0,
           hundun::flow::MomentumTimeOrder::backward_euler});
      hundun::flow::FlowLayerValues initial;
      initial.density.assign(topology.owned_cell_count(), density);
      initial.velocity.assign(topology.owned_cell_count() * 3U, 0.0);
      initial.mechanical_pressure.assign(topology.owned_cell_count(), 0.0);
      initial.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
      initial.face_mass_flux.assign(topology.local_face_count(), 0.0);
      initial.transported_cell_fields = {
          std::vector<double>(topology.owned_cell_count(), 3.0)};
      state.seed_accepted_layers(initial, initial);
      return state;
    };

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
    hundun::flow::MaterialDensityTransportSpec specification;
    specification.enthalpy_density = rho_h;
    specification.enthalpy_diffusivity_kg_per_m_s = 0.0;
    {
      hundun::mesh::MeshGeometry warped_geometry(
          topology, hundun::mesh::AnalyticWarpedBoxMapping(
                        {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0},
                        {0.02, -0.015, 0.01}));
      hundun::linear::BiCGStabSolver warped_pressure_solver(execution, mpi);
      hundun::linear::JacobiPreconditioner warped_pressure_preconditioner(
          execution);
      auto warped = hundun::flow::FixedStepMaterialDensityFlow::create(
          decomposition, topology, warped_geometry, boundaries, mpi,
          execution, halo, momentum_solver, {&mx, &my, &mz},
          warped_pressure_solver, warped_pressure_preconditioner, registry,
          fields, specification);
      auto warped_state = make_state(1.0);
      const auto warped_report = warped.attempt(
          warped_state, 0.0,
          hundun::flow::make_momentum_time_stencil(
              hundun::flow::MomentumTimeOrder::backward_euler, 0.01, 0.0),
          {}, {});
      HUNDUN_CHECK(warped_report.flow().disposition ==
                   hundun::flow::StepAttemptDisposition::committed);
      HUNDUN_CHECK(warped_report.flow().pressure_corrector_count == 2U);
      HUNDUN_CHECK(warped_state.metadata().step == 1U);
    }
    {
      auto wall_decomposition =
          hundun::runtime::StructuredDecomposition::create(
              mpi, extent, {false, true, true},
              hundun::runtime::DecompositionOptions{grid(mpi.size())});
      hundun::mesh::MeshTopology wall_topology(wall_decomposition);
      hundun::mesh::MeshGeometry wall_geometry(
          wall_topology, hundun::mesh::UniformBoxMapping(
                             {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}));
      auto wall_boundaries = hundun::boundary::BoundaryRegistry::create(
          wall_case(), wall_topology);
      const auto wall_box = wall_decomposition.owned_box();
      const hundun::runtime::Int3 wall_local{
          wall_box.end.x - wall_box.begin.x,
          wall_box.end.y - wall_box.begin.y,
          wall_box.end.z - wall_box.begin.z};
      auto wall_halo = hundun::runtime::HaloExchange::create(
          wall_decomposition, hundun::runtime::ExchangePlan::create(
                                  wall_decomposition, wall_local, 2));
      bool rejected = false;
      try {
        static_cast<void>(hundun::flow::FixedStepMaterialDensityFlow::create(
            wall_decomposition, wall_topology, wall_geometry, wall_boundaries,
            mpi, execution, wall_halo, momentum_solver, {&mx, &my, &mz},
            pressure_solver, pressure_preconditioner, registry, fields,
            specification));
      } catch (const hundun::runtime::Error &error) {
        rejected = std::string_view(error.what()) ==
                   "material flow requires fully periodic boundaries";
      }
      HUNDUN_CHECK(rejected);
    }
    auto flow = hundun::flow::FixedStepMaterialDensityFlow::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        momentum_solver, {&mx, &my, &mz}, pressure_solver,
        pressure_preconditioner, registry, fields, specification);
    const auto initial_pressure_halo =
        flow.pressure_halo_performance_counters();
    HUNDUN_CHECK(initial_pressure_halo.completed_exchanges == 0U);
    HUNDUN_CHECK(initial_pressure_halo.begin_calls == 0U);
    HUNDUN_CHECK(initial_pressure_halo.wait_calls == 0U);
    {
      auto movable_flow =
          hundun::flow::FixedStepMaterialDensityFlow::create(
              decomposition, topology, geometry, boundaries, mpi, execution,
              halo, momentum_solver, {&mx, &my, &mz}, pressure_solver,
              pressure_preconditioner, registry, fields, specification);
      hundun::flow::FixedStepMaterialDensityFlow moved_flow(
          std::move(movable_flow));
      HUNDUN_CHECK(moved_flow.pressure_halo_performance_counters()
                       .completed_exchanges == 0U);
      bool moved_from_rejected = false;
      try {
        static_cast<void>(
            movable_flow.pressure_halo_performance_counters());
      } catch (const hundun::runtime::Error &) {
        moved_from_rejected = true;
      }
      HUNDUN_CHECK(moved_from_rejected);
    }
    const auto be = hundun::flow::make_momentum_time_stencil(
        hundun::flow::MomentumTimeOrder::backward_euler, 0.01, 0.0);

    auto state = make_state(1.0);
    const auto report = [&] {
      task20_allocation_oracle::Observation allocation_observation(state);
      auto result = flow.attempt(state, 0.0, be, {}, {});
      HUNDUN_CHECK(allocation_observation.count() == 0U);
      return result;
    }();
    HUNDUN_CHECK(report.flow().disposition ==
                 hundun::flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(report.flow().reason == hundun::flow::StepFailureReason::none);
    HUNDUN_CHECK(report.flow().pressure_corrector_count == 2U);
    HUNDUN_CHECK(report.material_report_available());
    HUNDUN_CHECK(report.material_report().disposition() ==
                 hundun::flow::MaterialTransportDisposition::finalized);
    HUNDUN_CHECK(report.flux_provenance() ==
                 hundun::flow::MaterialFluxProvenance::final_corrected);
    HUNDUN_CHECK(report.material_field_count() == 1U);
    HUNDUN_CHECK(report.final_continuity_residual_available());
    HUNDUN_CHECK(report.final_pressure_residual_available());
    HUNDUN_CHECK(state.metadata().step == 1U);
    const auto pressure_halo =
        flow.pressure_halo_performance_counters();
    const auto pressure_halo_parts =
        hundun::flow::test::MaterialDensityPisoTestAccess::
            material_pressure_halo_counters(flow);
    HUNDUN_CHECK(pressure_halo_parts.ordinary_available);
    HUNDUN_CHECK(pressure_halo_parts.material_final_available);
    const auto expected_pressure_halo =
        hundun::flow::test::MaterialDensityPisoTestAccess::
            combine_pressure_halo_counters(
                pressure_halo_parts.ordinary,
                pressure_halo_parts.material_final);
    HUNDUN_CHECK(hundun::test::task25_counters_equal(
        pressure_halo, expected_pressure_halo));
    HUNDUN_CHECK(pressure_halo.completed_exchanges > 0U);
    {
      auto overflow_left = pressure_halo_parts.ordinary;
      auto overflow_right = pressure_halo_parts.material_final;
      overflow_left.completed_exchanges =
          std::numeric_limits<std::uint64_t>::max();
      overflow_right.completed_exchanges = 1U;
      bool overflow_rejected = false;
      try {
        static_cast<void>(
            hundun::flow::test::MaterialDensityPisoTestAccess::
                combine_pressure_halo_counters(overflow_left,
                                               overflow_right));
      } catch (const hundun::runtime::Error &) {
        overflow_rejected = true;
      }
      HUNDUN_CHECK(overflow_rejected);
    }
    const auto pressure_evidence =
        hundun::flow::test::MaterialDensityPisoTestAccess::
            material_pressure_evidence(flow);
    HUNDUN_CHECK(pressure_evidence.rhs_raw.size() ==
                 topology.owned_cell_count());
    HUNDUN_CHECK(pressure_evidence.rhs_solve.size() ==
                 topology.owned_cell_count());
    HUNDUN_CHECK(pressure_evidence.correction.size() ==
                 topology.owned_cell_count());
    HUNDUN_CHECK(pressure_evidence.final_face_density.size() ==
                 topology.local_face_count());
    HUNDUN_CHECK(pressure_evidence.corrector_ordinal == 2U);
    HUNDUN_CHECK(!pressure_evidence.token_available);
    HUNDUN_CHECK(pressure_evidence.final_operator_available);
    HUNDUN_CHECK(pressure_evidence.final_operator_revision >= 2U);
    HUNDUN_CHECK(std::all_of(
        pressure_evidence.final_face_density.begin(),
        pressure_evidence.final_face_density.end(),
        [](double value) { return value > 0.0 && std::isfinite(value); }));
    double projected_integral{};
    for (std::size_t cell_index = 0;
         cell_index < topology.owned_cell_count(); ++cell_index)
      projected_integral += geometry.cell_volume_m3(cell_index) *
                            pressure_evidence.rhs_solve[cell_index];
    mpi.allreduce_fp64_in_place(
        &projected_integral, 1U,
        hundun::runtime::Fp64ReductionOperation::sum);
    HUNDUN_CHECK(std::abs(projected_integral) <=
                 64.0 * std::numeric_limits<double>::epsilon());
    HUNDUN_CHECK(
        hundun::flow::test::MaterialDensityPisoTestAccess::
            report_authenticated(report));
    for (const auto corruption : {
             hundun::flow::test::MaterialReportCorruptionForTest::
                 success_corrector_count,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 success_provenance,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 success_shared_field,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 parent_transport_size,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 material_transport_size,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 unavailable_numeric_value,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 material_reason_mapping,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 reliable_collective_rank,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 material_count_zero,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 material_count_plus_two,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 material_count_five,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 parent_transport_residual_value,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 nested_transport_residual_value,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 parent_transport_conservation_value,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 nested_transport_conservation_value,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 nested_density_residual_availability,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 nested_transport_residual_availability,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 nested_mass_conservation_availability,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 nested_transport_conservation_availability,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 nested_minimum_density_availability,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 nested_attempt_identity,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 nested_finalization_identity,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 nested_shared_field,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 nested_provenance,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 nested_residual_outer_size,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 nested_conservation_outer_size}) {
      auto mutated = report;
      hundun::flow::test::MaterialDensityPisoTestAccess::corrupt_report(
          mutated, corruption);
      HUNDUN_CHECK(!hundun::flow::test::MaterialDensityPisoTestAccess::
                        report_authenticated(mutated));
    }

    const auto exercise_terminal_solver = [&](TerminalMomentumSolver::Mode mode,
                                              std::string_view message) {
      TerminalMomentumSolver terminal_solver(mode);
      hundun::linear::JacobiPreconditioner tx(execution), ty(execution),
          tz(execution), tp(execution);
      auto terminal_flow = hundun::flow::FixedStepMaterialDensityFlow::create(
          decomposition, topology, geometry, boundaries, mpi, execution, halo,
          terminal_solver, {&tx, &ty, &tz}, pressure_solver, tp, registry,
          fields, specification);
      auto terminal_state = make_state(1.0);
      const auto before = capture(terminal_state);
      bool exact = false;
      try {
        static_cast<void>(
            terminal_flow.attempt(terminal_state, 0.0, be, {}, {}));
      } catch (const hundun::runtime::MpiOperationError &error) {
        exact = std::string_view(error.what()) == message;
      }
      HUNDUN_CHECK(exact);
      HUNDUN_CHECK(terminal_solver.calls() == 1U);
      check_equal(before, terminal_state);
    };
    exercise_terminal_solver(
        TerminalMomentumSolver::Mode::rankless_result,
        "material-density flow collective failure has no reliable rank");
    exercise_terminal_solver(
        TerminalMomentumSolver::Mode::thrown_operation,
        "injected material momentum operation failure");
    {
      TerminalMomentumSolver terminal_solver(
          TerminalMomentumSolver::Mode::reliable_result);
      hundun::linear::JacobiPreconditioner tx(execution), ty(execution),
          tz(execution), tp(execution);
      auto terminal_flow = hundun::flow::FixedStepMaterialDensityFlow::create(
          decomposition, topology, geometry, boundaries, mpi, execution, halo,
          terminal_solver, {&tx, &ty, &tz}, pressure_solver, tp, registry,
          fields, specification);
      auto terminal_state = make_state(1.0);
      const auto before = capture(terminal_state);
      const auto terminal_report =
          terminal_flow.attempt(terminal_state, 0.0, be, {}, {});
      HUNDUN_CHECK(terminal_report.flow().disposition ==
                   hundun::flow::StepAttemptDisposition::
                       non_retryable_failure);
      HUNDUN_CHECK(terminal_report.flow().reason ==
                   hundun::flow::StepFailureReason::collective_operation);
      HUNDUN_CHECK(terminal_report.flow().lowest_failing_rank == 0);
      HUNDUN_CHECK(
          hundun::flow::test::MaterialDensityPisoTestAccess::
              report_authenticated(terminal_report));
      HUNDUN_CHECK(terminal_solver.calls() == 1U);
      check_equal(before, terminal_state);
    }
    for (const auto &[mode, fail_on_call] :
         std::array{std::pair{TerminalMomentumSolver::Mode::maximum_iterations,
                             1U},
                    std::pair{TerminalMomentumSolver::Mode::numerical_breakdown,
                              2U},
                    std::pair{TerminalMomentumSolver::Mode::non_finite_value,
                              3U}}) {
      TerminalMomentumSolver terminal_solver(mode, fail_on_call);
      hundun::linear::JacobiPreconditioner tx(execution), ty(execution),
          tz(execution), tp(execution);
      auto terminal_flow = hundun::flow::FixedStepMaterialDensityFlow::create(
          decomposition, topology, geometry, boundaries, mpi, execution, halo,
          terminal_solver, {&tx, &ty, &tz}, pressure_solver, tp, registry,
          fields, specification);
      auto terminal_state = make_state(1.0);
      const auto before = capture(terminal_state);
      const auto terminal_report =
          terminal_flow.attempt(terminal_state, 0.0, be, {}, {});
      HUNDUN_CHECK(terminal_report.flow().disposition ==
                   hundun::flow::StepAttemptDisposition::recoverable_failure);
      HUNDUN_CHECK(terminal_report.flow().reason ==
                   hundun::flow::StepFailureReason::momentum_linear_solve);
      HUNDUN_CHECK(terminal_report.flow().lowest_failing_rank == 0);
      HUNDUN_CHECK(
          hundun::flow::test::MaterialDensityPisoTestAccess::
              report_authenticated(terminal_report));
      HUNDUN_CHECK(terminal_solver.calls() == fail_on_call);
      check_equal(before, terminal_state);
    }
    const auto exercise_terminal_pressure = [&](TerminalMomentumSolver::Mode mode,
                                                std::size_t fail_on_call,
                                                std::string_view message) {
      TerminalMomentumSolver terminal_solver(mode, fail_on_call);
      hundun::linear::JacobiPreconditioner tx(execution), ty(execution),
          tz(execution), tp(execution);
      auto terminal_flow = hundun::flow::FixedStepMaterialDensityFlow::create(
          decomposition, topology, geometry, boundaries, mpi, execution, halo,
          momentum_solver, {&tx, &ty, &tz}, terminal_solver, tp, registry,
          fields, specification);
      auto terminal_state = make_state(1.0);
      const auto before = capture(terminal_state);
      bool exact = false;
      try {
        static_cast<void>(
            terminal_flow.attempt(terminal_state, 0.0, be, {}, {}));
      } catch (const hundun::runtime::MpiOperationError &error) {
        exact = std::string_view(error.what()) == message;
      }
      HUNDUN_CHECK(exact);
      HUNDUN_CHECK(terminal_solver.calls() == fail_on_call);
      check_equal(before, terminal_state);
    };
    exercise_terminal_pressure(
        TerminalMomentumSolver::Mode::rankless_result, 1U,
        "material-density flow collective failure has no reliable rank");
    exercise_terminal_pressure(
        TerminalMomentumSolver::Mode::rankless_result, 2U,
        "material-density flow collective failure has no reliable rank");
    exercise_terminal_pressure(
        TerminalMomentumSolver::Mode::thrown_operation, 1U,
        "injected material momentum operation failure");
    for (const std::size_t fail_on_call : {1U, 2U}) {
      TerminalMomentumSolver terminal_solver(
          TerminalMomentumSolver::Mode::reliable_result, fail_on_call);
      hundun::linear::JacobiPreconditioner tx(execution), ty(execution),
          tz(execution), tp(execution);
      auto terminal_flow = hundun::flow::FixedStepMaterialDensityFlow::create(
          decomposition, topology, geometry, boundaries, mpi, execution, halo,
          momentum_solver, {&tx, &ty, &tz}, terminal_solver, tp, registry,
          fields, specification);
      auto terminal_state = make_state(1.0);
      const auto before = capture(terminal_state);
      const auto terminal_report =
          terminal_flow.attempt(terminal_state, 0.0, be, {}, {});
      HUNDUN_CHECK(terminal_report.flow().disposition ==
                   hundun::flow::StepAttemptDisposition::
                       non_retryable_failure);
      HUNDUN_CHECK(terminal_report.flow().reason ==
                   hundun::flow::StepFailureReason::collective_operation);
      HUNDUN_CHECK(terminal_report.flow().lowest_failing_rank == 0);
      HUNDUN_CHECK(
          hundun::flow::test::MaterialDensityPisoTestAccess::
              report_authenticated(terminal_report));
      HUNDUN_CHECK(terminal_solver.calls() == fail_on_call);
      check_equal(before, terminal_state);
    }
    for (const auto &[mode, fail_on_call] :
         std::array{std::pair{TerminalMomentumSolver::Mode::maximum_iterations,
                             1U},
                    std::pair{TerminalMomentumSolver::Mode::numerical_breakdown,
                              2U}}) {
      TerminalMomentumSolver terminal_solver(mode, fail_on_call);
      hundun::linear::JacobiPreconditioner tx(execution), ty(execution),
          tz(execution), tp(execution);
      auto terminal_flow = hundun::flow::FixedStepMaterialDensityFlow::create(
          decomposition, topology, geometry, boundaries, mpi, execution, halo,
          momentum_solver, {&tx, &ty, &tz}, terminal_solver, tp, registry,
          fields, specification);
      auto terminal_state = make_state(1.0);
      const auto before = capture(terminal_state);
      const auto terminal_report =
          terminal_flow.attempt(terminal_state, 0.0, be, {}, {});
      HUNDUN_CHECK(terminal_report.flow().disposition ==
                   hundun::flow::StepAttemptDisposition::recoverable_failure);
      HUNDUN_CHECK(terminal_report.flow().reason ==
                   hundun::flow::StepFailureReason::pressure_linear_solve);
      HUNDUN_CHECK(terminal_report.flow().lowest_failing_rank == 0);
      HUNDUN_CHECK(
          hundun::flow::test::MaterialDensityPisoTestAccess::
              report_authenticated(terminal_report));
      HUNDUN_CHECK(terminal_solver.calls() == fail_on_call);
      check_equal(before, terminal_state);
    }

    using TerminalMode =
        hundun::flow::test::MaterialTerminalModeForTest;
    using TerminalPoint =
        hundun::flow::test::MaterialTerminalPointForTest;
    using TestAccess = hundun::flow::test::MaterialDensityPisoTestAccess;
    const auto exercise_terminal_point = [&](TerminalPoint point,
                                             TerminalMode mode) {
      TestAccess::reset_terminal_fault();
      TestAccess::set_terminal_fault(point, mode);
      hundun::linear::JacobiPreconditioner tx(execution), ty(execution),
          tz(execution), tp(execution);
      auto terminal_flow = hundun::flow::FixedStepMaterialDensityFlow::create(
          decomposition, topology, geometry, boundaries, mpi, execution, halo,
          momentum_solver, {&tx, &ty, &tz}, pressure_solver, tp, registry,
          fields, specification);
      auto terminal_state = make_state(1.0);
      const auto before = capture(terminal_state);
      bool exact = false;
      try {
        static_cast<void>(
            terminal_flow.attempt(terminal_state, 0.0, be, {}, {}));
      } catch (const hundun::runtime::MpiOperationError &error) {
        exact = std::string_view(error.what()) ==
                (mode == TerminalMode::returned_rankless
                     ? "material-density flow collective failure has no "
                       "reliable rank"
                     : "injected material terminal operation failure");
      }
      HUNDUN_CHECK(exact);
      HUNDUN_CHECK(!TestAccess::has_diagnostic_report(terminal_flow));
      check_equal(before, terminal_state);
      HUNDUN_CHECK(TestAccess::terminal_point_calls(point) == 1U);
      for (std::size_t later = static_cast<std::size_t>(point) + 1U;
           later < static_cast<std::size_t>(TerminalPoint::count); ++later)
        HUNDUN_CHECK(TestAccess::terminal_point_calls(
                         static_cast<TerminalPoint>(later)) == 0U);
      TestAccess::reset_terminal_fault();
    };
    constexpr std::array terminal_points{
        TerminalPoint::predictor_stage,
        TerminalPoint::momentum_x,
        TerminalPoint::momentum_y,
        TerminalPoint::momentum_z,
        TerminalPoint::pressure_corrector_one,
        TerminalPoint::provisional_stage,
        TerminalPoint::pressure_corrector_two,
        TerminalPoint::public_finalizer,
        TerminalPoint::final_continuity_reduction,
        TerminalPoint::final_continuity_status,
        TerminalPoint::final_pressure_entry,
        TerminalPoint::final_pressure_gamma_sum,
        TerminalPoint::final_pressure_gamma_count,
        TerminalPoint::final_pressure_residual_reduction,
        TerminalPoint::final_momentum_residual_reduction,
        TerminalPoint::final_momentum_conservation_reduction,
        TerminalPoint::final_momentum_status,
        TerminalPoint::final_conservation_status};
    for (const auto point : terminal_points) {
      exercise_terminal_point(point, TerminalMode::returned_rankless);
      exercise_terminal_point(point, TerminalMode::thrown_operation);
    }
    for (const auto point : {TerminalPoint::predictor_stage,
                             TerminalPoint::public_finalizer,
                             TerminalPoint::final_pressure_entry}) {
      TestAccess::reset_terminal_fault();
      TestAccess::set_terminal_fault(point, TerminalMode::returned_reliable);
      hundun::linear::JacobiPreconditioner tx(execution), ty(execution),
          tz(execution), tp(execution);
      auto terminal_flow = hundun::flow::FixedStepMaterialDensityFlow::create(
          decomposition, topology, geometry, boundaries, mpi, execution, halo,
          momentum_solver, {&tx, &ty, &tz}, pressure_solver, tp, registry,
          fields, specification);
      auto terminal_state = make_state(1.0);
      const auto before = capture(terminal_state);
      const auto terminal_report =
          terminal_flow.attempt(terminal_state, 0.0, be, {}, {});
      HUNDUN_CHECK(terminal_report.flow().disposition ==
                   hundun::flow::StepAttemptDisposition::
                       non_retryable_failure);
      HUNDUN_CHECK(terminal_report.flow().reason ==
                   hundun::flow::StepFailureReason::collective_operation);
      HUNDUN_CHECK(terminal_report.flow().lowest_failing_rank == 0);
      HUNDUN_CHECK(TestAccess::report_authenticated(terminal_report));
      check_equal(before, terminal_state);
      TestAccess::reset_terminal_fault();
    }

    {
      using TransportAccess =
          hundun::flow::test::MaterialDensityTransportTestAccess;
      enum class LateGate { transport_residual, mass, transport_conservation };
      for (const auto gate : {LateGate::transport_residual, LateGate::mass,
                              LateGate::transport_conservation}) {
        hundun::linear::JacobiPreconditioner lx(execution), ly(execution),
            lz(execution), lp(execution);
        auto late_flow = hundun::flow::FixedStepMaterialDensityFlow::create(
            decomposition, topology, geometry, boundaries, mpi, execution,
            halo, momentum_solver, {&lx, &ly, &lz}, pressure_solver, lp,
            registry, fields, specification);
        auto late_state = make_state(1.0);
        const auto before = capture(late_state);
        const int failing_rank = mpi.size() - 1;
        const double injected =
            gate == LateGate::transport_residual ? 2.0e-9 : 1.0e-9;
        TransportAccess::reset();
        if (gate == LateGate::transport_residual)
          TransportAccess::set_transport_residual(0U, injected, failing_rank);
        else if (gate == LateGate::mass)
          TransportAccess::set_mass_conservation_defect(injected,
                                                        failing_rank);
        else
          TransportAccess::set_transport_conservation_defect(
              0U, injected, failing_rank);
        std::optional<hundun::flow::MaterialDensityStepAttemptReport>
            late_report;
        try {
          late_report.emplace(
              late_flow.attempt(late_state, 0.0, be, {}, {}));
        } catch (...) {
          TransportAccess::reset();
          throw;
        }
        TransportAccess::reset();
        const auto expected_reason =
            gate == LateGate::transport_residual
                ? hundun::flow::StepFailureReason::final_transport_residual
                : hundun::flow::StepFailureReason::final_conservation_defect;
        HUNDUN_CHECK(late_report->flow().disposition ==
                     hundun::flow::StepAttemptDisposition::
                         recoverable_failure);
        HUNDUN_CHECK(late_report->flow().reason == expected_reason);
        HUNDUN_CHECK(late_report->flow().lowest_failing_rank == failing_rank);
        HUNDUN_CHECK(late_report->flow().pressure_corrector_count == 2U);
        HUNDUN_CHECK(late_report->material_report_available());
        const auto &nested = late_report->material_report();
        HUNDUN_CHECK(nested.disposition() ==
                     hundun::flow::MaterialTransportDisposition::
                         recoverable_failure);
        HUNDUN_CHECK(
            nested.reason() ==
            (gate == LateGate::transport_residual
                 ? hundun::flow::MaterialTransportFailureReason::
                       final_transport_residual
                 : hundun::flow::MaterialTransportFailureReason::
                       final_conservation_defect));
        HUNDUN_CHECK(nested.lowest_failing_rank() == failing_rank);
        HUNDUN_CHECK(nested.flux_provenance() ==
                     hundun::flow::MaterialFluxProvenance::final_corrected);
        HUNDUN_CHECK(late_report->flux_provenance() ==
                     hundun::flow::MaterialFluxProvenance::final_corrected);
        HUNDUN_CHECK(
            late_report->flow().final_transport_normalized_l2.size() ==
            nested.transport_normalized_l2().size());
        HUNDUN_CHECK(late_report->flow()
                             .final_transport_relative_conservation_defect
                             .size() ==
                     nested.transport_relative_conservation_defect().size());
        for (std::size_t field = 0;
             field < nested.transport_normalized_l2().size(); ++field) {
          HUNDUN_CHECK(
              nested.transport_residual_availability()[field] == 1U);
          HUNDUN_CHECK(
              nested.transport_conservation_availability()[field] == 1U);
          HUNDUN_CHECK(bits(late_report->flow()
                                .final_transport_normalized_l2[field]) ==
                       bits(nested.transport_normalized_l2()[field]));
          HUNDUN_CHECK(
              bits(late_report->flow()
                       .final_transport_relative_conservation_defect[field]) ==
              bits(nested.transport_relative_conservation_defect()[field]));
        }
        HUNDUN_CHECK(late_report->mass_conservation_available() ==
                     nested.mass_conservation_available());
        HUNDUN_CHECK(late_report->mass_conservation_available());
        HUNDUN_CHECK(
            bits(late_report->flow()
                     .final_mass_relative_conservation_defect) ==
            bits(nested.mass_relative_conservation_defect()));
        if (mpi.rank() == failing_rank) {
          if (gate == LateGate::transport_residual)
            HUNDUN_CHECK(bits(nested.transport_normalized_l2()[0]) ==
                         bits(injected));
          else if (gate == LateGate::mass)
            HUNDUN_CHECK(bits(nested.mass_relative_conservation_defect()) ==
                         bits(injected));
          else
            HUNDUN_CHECK(
                bits(nested.transport_relative_conservation_defect()[0]) ==
                bits(injected));
        }
        HUNDUN_CHECK(
            hundun::flow::test::MaterialDensityPisoTestAccess::
                report_authenticated(*late_report));
        check_equal(before, late_state);
        auto late_source =
            late_flow.diagnostic_source(late_state, *late_report);
        HUNDUN_CHECK(late_source.report().attempt_identity() ==
                     late_report->attempt_identity());
        const auto invalidation =
            late_flow.attempt(late_state, -1.0, be, {}, {});
        HUNDUN_CHECK(invalidation.flow().reason ==
                     hundun::flow::StepFailureReason::invalid_input);
        bool stale = false;
        try {
          static_cast<void>(late_source.report());
        } catch (const hundun::runtime::Error &error) {
          stale = std::string_view(error.what()) ==
                  "material flow diagnostic source is stale";
        }
        HUNDUN_CHECK(stale);
        check_equal(before, late_state);
      }
    }

    auto source = flow.diagnostic_source(state, report);
    const auto expect_stale_query = [](auto &&query) {
      bool stale_query = false;
      try {
        query();
      } catch (const hundun::runtime::Error &error) {
        stale_query = std::string_view(error.what()) ==
                      "material flow diagnostic source is stale";
      }
      HUNDUN_CHECK(stale_query);
    };
    const auto expect_all_queries_stale = [&](auto &candidate) {
      expect_stale_query(
          [&] { static_cast<void>(candidate.fingerprint_field_count()); });
      expect_stale_query(
          [&] { static_cast<void>(candidate.fingerprint_field_id(0U)); });
      expect_stale_query(
          [&] { static_cast<void>(candidate.field_unit(0U)); });
      expect_stale_query(
          [&] { static_cast<void>(candidate.field_entity(0U)); });
      expect_stale_query(
          [&] { static_cast<void>(candidate.field_component_count(0U)); });
      expect_stale_query(
          [&] { static_cast<void>(candidate.field_item_count(0U)); });
      expect_stale_query(
          [&] { static_cast<void>(candidate.field_global_id(0U, 0U)); });
      expect_stale_query([&] {
        static_cast<void>(candidate.field_value(0U, 0U, 0U));
      });
      expect_stale_query([&] {
        static_cast<void>(candidate.owned_cell_layout_fingerprint());
      });
      expect_stale_query([&] {
        static_cast<void>(candidate.global_cell_layout_fingerprint());
      });
      expect_stale_query([&] {
        static_cast<void>(candidate.owned_face_layout_fingerprint());
      });
      expect_stale_query([&] {
        static_cast<void>(candidate.global_face_layout_fingerprint());
      });
      expect_stale_query(
          [&] { static_cast<void>(candidate.relative_rank()); });
      expect_stale_query(
          [&] { static_cast<void>(candidate.committed_step()); });
      expect_stale_query(
          [&] { static_cast<void>(candidate.committed_time_s()); });
      expect_stale_query(
          [&] { static_cast<void>(candidate.global_cell_extent()); });
      expect_stale_query(
          [&] { static_cast<void>(candidate.owned_global_box()); });
      expect_stale_query(
          [&] { static_cast<void>(candidate.owned_cell_count()); });
      expect_stale_query([&] {
        static_cast<void>(candidate.canonical_owned_face_count());
      });
      expect_stale_query(
          [&] { static_cast<void>(candidate.cell_volume_m3(0U)); });
      expect_stale_query([&] { static_cast<void>(candidate.report()); });
    };
    HUNDUN_CHECK(source.fingerprint_field_count() == 5U);
    HUNDUN_CHECK(source.fingerprint_field_id(0) == "face_mass_flux");
    HUNDUN_CHECK(source.fingerprint_field_id(4) == "velocity");
    HUNDUN_CHECK(source.committed_step() == 1U);
    HUNDUN_CHECK(source.field_item_count(3) == topology.owned_cell_count());
    HUNDUN_CHECK(source.field_value(3, 0U, 0U) == 1.0);

    auto moved_flow = std::move(flow);
    expect_all_queries_stale(source);
    auto fresh_source = moved_flow.diagnostic_source(state, report);
    HUNDUN_CHECK(fresh_source.committed_step() == 1U);
    bool moved_from_rejected = false;
    try {
      static_cast<void>(flow.diagnostic_source(state, report));
    } catch (const hundun::runtime::Error &error) {
      moved_from_rejected = std::string_view(error.what()) ==
                            "material flow object has been moved from";
    }
    HUNDUN_CHECK(moved_from_rejected);

    const auto expect_invalid = [&](hundun::flow::FlowState &candidate,
                                    const auto &candidate_stencil,
                                    const auto &momentum_control,
                                    const auto &pressure_control) {
      const auto before = capture(candidate);
      const auto invalid_report = moved_flow.attempt(
          candidate, 0.0, candidate_stencil, momentum_control,
          pressure_control);
      HUNDUN_CHECK(invalid_report.flow().disposition ==
                   hundun::flow::StepAttemptDisposition::
                       non_retryable_failure);
      HUNDUN_CHECK(invalid_report.flow().reason ==
                   hundun::flow::StepFailureReason::invalid_input);
      HUNDUN_CHECK(TestAccess::report_authenticated(invalid_report));
      check_equal(before, candidate);
    };
    {
      auto wrong_fields = fields;
      wrong_fields.density = alternate_density;
      auto wrong_state = hundun::flow::FlowState::create(
          registry, {local, topology.local_face_count()}, wrong_fields,
          {0U, 0.0, 0.01, 0.0,
           hundun::flow::MomentumTimeOrder::backward_euler});
      auto initial = make_state(1.0).snapshot(
          hundun::flow::FlowLayer::committed);
      wrong_state.seed_accepted_layers(initial, initial);
      expect_invalid(wrong_state, be, hundun::linear::SolveControl{},
                     hundun::linear::SolveControl{});
    }
    {
      hundun::runtime::FieldRegistry other_registry;
      hundun::flow::FlowFieldIds other_fields;
      other_fields.density =
          other_registry.declare_field(cell("rho", "kg/m3", 1U, true));
      other_fields.velocity = other_registry.declare_field(
          cell("velocity", "m/s", 3U, false));
      other_fields.mechanical_pressure = other_registry.declare_field(
          cell("pi", "Pa", 1U, false));
      other_fields.face_velocity =
          other_registry.declare_field(face("face_velocity", "m/s", 3U));
      other_fields.face_mass_flux =
          hundun::finite_volume::declare_face_mass_flux(other_registry);
      const auto other_rho_h = other_registry.declare_field(
          cell("rho_h", "J/m3", 1U, true));
      other_fields.transported_cell_fields = {other_rho_h};
      other_registry.freeze();
      auto other_state = hundun::flow::FlowState::create(
          other_registry, {local, topology.local_face_count()}, other_fields,
          {0U, 0.0, 0.01, 0.0,
           hundun::flow::MomentumTimeOrder::backward_euler});
      auto initial = make_state(1.0).snapshot(
          hundun::flow::FlowLayer::committed);
      other_state.seed_accepted_layers(initial, initial);
      expect_invalid(other_state, be, hundun::linear::SolveControl{},
                     hundun::linear::SolveControl{});
    }
    {
      const hundun::runtime::Int3 wrong_local{local.x + 1, local.y, local.z};
      auto wrong_layout_state = hundun::flow::FlowState::create(
          registry, {wrong_local, topology.local_face_count()}, fields,
          {0U, 0.0, 0.01, 0.0,
           hundun::flow::MomentumTimeOrder::backward_euler});
      const auto cells = static_cast<std::size_t>(wrong_local.x) *
                         static_cast<std::size_t>(wrong_local.y) *
                         static_cast<std::size_t>(wrong_local.z);
      hundun::flow::FlowLayerValues initial;
      initial.density.assign(cells, 1.0);
      initial.velocity.assign(cells * 3U, 0.0);
      initial.mechanical_pressure.assign(cells, 0.0);
      initial.face_velocity.assign(topology.local_face_count() * 3U, 0.0);
      initial.face_mass_flux.assign(topology.local_face_count(), 0.0);
      initial.transported_cell_fields = {std::vector<double>(cells, 3.0)};
      wrong_layout_state.seed_accepted_layers(initial, initial);
      expect_invalid(wrong_layout_state, be, hundun::linear::SolveControl{},
                     hundun::linear::SolveControl{});
    }
    {
      auto wrong_stencil = be;
      wrong_stencil.alpha0 += 1.0;
      auto candidate = make_state(1.0);
      expect_invalid(candidate, wrong_stencil,
                     hundun::linear::SolveControl{},
                     hundun::linear::SolveControl{});
      auto invalid_control = hundun::linear::SolveControl{};
      invalid_control.residual_recompute_interval = 0U;
      candidate = make_state(1.0);
      expect_invalid(candidate, be, invalid_control,
                     hundun::linear::SolveControl{});
      candidate = make_state(1.0);
      expect_invalid(candidate, be, hundun::linear::SolveControl{},
                     invalid_control);
    }

    auto invalid_state = make_state(1.0);
    const auto invalid_before = capture(invalid_state);
    const auto invalid =
        moved_flow.attempt(invalid_state, -1.0, be, {}, {});
    HUNDUN_CHECK(invalid.flow().disposition ==
                 hundun::flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(invalid.flow().reason ==
                 hundun::flow::StepFailureReason::invalid_input);
    check_equal(invalid_before, invalid_state);
    expect_all_queries_stale(fresh_source);

    auto negative_state = make_state(-1.0);
    const auto negative_before = capture(negative_state);
    const auto negative =
        moved_flow.attempt(negative_state, 0.0, be, {}, {});
    HUNDUN_CHECK(negative.flow().disposition ==
                 hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(negative.flow().reason ==
                 hundun::flow::StepFailureReason::transport_failure);
    HUNDUN_CHECK(negative.material_failure_reason() ==
                 hundun::flow::MaterialTransportFailureReason::
                     non_positive_density);
    HUNDUN_CHECK(negative.flow().suggested_dt_s == 0.005);
    HUNDUN_CHECK(
        hundun::flow::test::MaterialDensityPisoTestAccess::
            report_authenticated(negative));
    for (const auto corruption : {
             hundun::flow::test::MaterialReportCorruptionForTest::
                 material_reason_mapping,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 unavailable_numeric_value,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 material_count_zero,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 material_count_plus_two,
             hundun::flow::test::MaterialReportCorruptionForTest::
                 material_count_five}) {
      auto mutated = negative;
      hundun::flow::test::MaterialDensityPisoTestAccess::corrupt_report(
          mutated, corruption);
      HUNDUN_CHECK(!hundun::flow::test::MaterialDensityPisoTestAccess::
                        report_authenticated(mutated));
    }
    check_equal(negative_before, negative_state);

    auto allocation_state = make_state(1.0);
    const auto allocation_before = capture(allocation_state);
    hundun::flow::test::MaterialDensityPisoTestAccess::
        set_preflight_allocation_failure_rank(mpi.size() - 1);
    const auto identity_before_preflight_failure = negative.attempt_identity();
    std::optional<hundun::flow::MaterialDensityStepAttemptReport>
        allocation_report;
    try {
      allocation_report.emplace(
          moved_flow.attempt(allocation_state, 0.0, be, {}, {}));
    } catch (...) {
      hundun::flow::test::MaterialDensityPisoTestAccess::
          reset_preflight_allocation_failure();
      throw;
    }
    hundun::flow::test::MaterialDensityPisoTestAccess::
        reset_preflight_allocation_failure();
    HUNDUN_CHECK(allocation_report->flow().disposition ==
                 hundun::flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(allocation_report->flow().reason ==
                 hundun::flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(allocation_report->flow().lowest_failing_rank ==
                 mpi.size() - 1);
    HUNDUN_CHECK(allocation_report->flow().pressure_corrector_count == 0U);
    HUNDUN_CHECK(allocation_report->material_field_count() == 1U);
    HUNDUN_CHECK(!allocation_report->material_report_available());
    HUNDUN_CHECK(!allocation_report->final_continuity_residual_available());
    HUNDUN_CHECK(!allocation_report->final_pressure_residual_available());
    HUNDUN_CHECK(allocation_report->flow().attempted_dt_s == 0.01);
    HUNDUN_CHECK(allocation_report->flow().suggested_dt_s == 0.01);
    HUNDUN_CHECK(
        allocation_report->flow().final_transport_normalized_l2.size() == 1U);
    HUNDUN_CHECK(allocation_report->flow()
                     .final_transport_relative_conservation_defect.size() ==
                 1U);
    HUNDUN_CHECK(allocation_report->flow().final_transport_normalized_l2[0] ==
                 0.0);
    HUNDUN_CHECK(allocation_report->flow()
                     .final_transport_relative_conservation_defect[0] == 0.0);
    HUNDUN_CHECK(allocation_report->attempt_identity() ==
                 identity_before_preflight_failure + 1U);
    HUNDUN_CHECK(
        hundun::flow::test::MaterialDensityPisoTestAccess::report_authenticated(
            *allocation_report));
    check_equal(allocation_before, allocation_state);
    auto allocation_source =
        moved_flow.diagnostic_source(allocation_state, *allocation_report);
    HUNDUN_CHECK(allocation_source.report().attempt_identity() ==
                 allocation_report->attempt_identity());

    const auto bdf2 = hundun::flow::make_momentum_time_stencil(
        hundun::flow::MomentumTimeOrder::bdf2, 0.01, 0.01);
    const auto second = moved_flow.attempt(state, 0.0, bdf2, {}, {});
    HUNDUN_CHECK(second.flow().disposition ==
                 hundun::flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(state.metadata().step == 2U);
    expect_all_queries_stale(allocation_source);

    {
      hundun::linear::JacobiPreconditioner ax(execution), ay(execution),
          az(execution), ap(execution), bx(execution), by(execution),
          bz(execution), bp(execution);
      auto first_flow = hundun::flow::FixedStepMaterialDensityFlow::create(
          decomposition, topology, geometry, boundaries, mpi, execution, halo,
          momentum_solver, {&ax, &ay, &az}, pressure_solver, ap, registry,
          fields, specification);
      auto lifetime_state = make_state(1.0);
      const auto first_report =
          first_flow.attempt(lifetime_state, 0.0, be, {}, {});
      auto first_source =
          first_flow.diagnostic_source(lifetime_state, first_report);
      auto second_flow = hundun::flow::FixedStepMaterialDensityFlow::create(
          decomposition, topology, geometry, boundaries, mpi, execution, halo,
          momentum_solver, {&bx, &by, &bz}, pressure_solver, bp, registry,
          fields, specification);
      const auto second_report =
          second_flow.attempt(lifetime_state, 0.0, be, {}, {});
      expect_all_queries_stale(first_source);
      auto second_source =
          second_flow.diagnostic_source(lifetime_state, second_report);
      auto moved_state = std::move(lifetime_state);
      expect_all_queries_stale(second_source);
      HUNDUN_CHECK(moved_state.metadata().step == 2U);
    }
    {
      hundun::linear::JacobiPreconditioner sx(execution), sy(execution),
          sz(execution), sp(execution);
      auto reseed_flow = hundun::flow::FixedStepMaterialDensityFlow::create(
          decomposition, topology, geometry, boundaries, mpi, execution, halo,
          momentum_solver, {&sx, &sy, &sz}, pressure_solver, sp, registry,
          fields, specification);
      auto reseed_state = make_state(1.0);
      const auto reseed_report =
          reseed_flow.attempt(reseed_state, 0.0, be, {}, {});
      auto reseed_source =
          reseed_flow.diagnostic_source(reseed_state, reseed_report);
      const auto accepted =
          reseed_state.snapshot(hundun::flow::FlowLayer::committed);
      reseed_state.seed_accepted_layers(accepted, accepted);
      expect_all_queries_stale(reseed_source);
    }
    {
      auto epoch_state = make_state(1.0);
      const auto accepted =
          epoch_state.snapshot(hundun::flow::FlowLayer::committed);
      TestAccess::force_state_diagnostic_identity(
          epoch_state, std::numeric_limits<std::uint64_t>::max() - 1U);
      epoch_state.seed_accepted_layers(accepted, accepted);
      HUNDUN_CHECK(TestAccess::state_diagnostic_identity(epoch_state) ==
                   std::numeric_limits<std::uint64_t>::max());
      const auto epoch_before = capture(epoch_state);
      for (int repetition = 0; repetition < 2; ++repetition) {
        bool seed_exhausted = false;
        try {
          epoch_state.seed_accepted_layers(accepted, accepted);
        } catch (const hundun::runtime::Error &error) {
          seed_exhausted = std::string_view(error.what()) ==
                           "FlowState diagnostic mutation identity would wrap";
        }
        HUNDUN_CHECK(seed_exhausted);
        check_equal(epoch_before, epoch_state);
        HUNDUN_CHECK(!TestAccess::state_attempt_active(epoch_state));
      }
      for (int repetition = 0; repetition < 2; ++repetition) {
        bool begin_exhausted = false;
        try {
          epoch_state.begin_attempt();
        } catch (const hundun::runtime::Error &error) {
          begin_exhausted = std::string_view(error.what()) ==
                            "FlowState diagnostic mutation identity would "
                            "wrap";
        }
        HUNDUN_CHECK(begin_exhausted);
        check_equal(epoch_before, epoch_state);
        HUNDUN_CHECK(!TestAccess::state_attempt_active(epoch_state));
      }
    }
    {
      auto begin_state = make_state(1.0);
      TestAccess::force_state_diagnostic_identity(
          begin_state, std::numeric_limits<std::uint64_t>::max() - 1U);
      begin_state.begin_attempt();
      HUNDUN_CHECK(TestAccess::state_diagnostic_identity(begin_state) ==
                   std::numeric_limits<std::uint64_t>::max());
      begin_state.rollback_attempt();
    }

    hundun::flow::test::MaterialDensityPisoTestAccess::
        force_flow_attempt_identity(moved_flow,
                                    std::numeric_limits<std::uint64_t>::max() -
                                        1U);
    const auto maximum_identity =
        moved_flow.attempt(state, 0.0, bdf2, {}, {});
    HUNDUN_CHECK(maximum_identity.attempt_identity() ==
                 std::numeric_limits<std::uint64_t>::max());
    HUNDUN_CHECK(maximum_identity.flow().disposition ==
                 hundun::flow::StepAttemptDisposition::committed);
    HUNDUN_CHECK(
        hundun::flow::test::MaterialDensityPisoTestAccess::
            report_authenticated(maximum_identity));
    auto maximum_source =
        moved_flow.diagnostic_source(state, maximum_identity);
    const auto exhaustion_before = capture(state);
    for (int repetition = 0; repetition < 2; ++repetition) {
      bool exhausted = false;
      try {
        static_cast<void>(moved_flow.attempt(state, 0.0, bdf2, {}, {}));
      } catch (const hundun::runtime::Error &error) {
        exhausted = std::string_view(error.what()) ==
                    "material flow attempt identity would wrap";
      }
      HUNDUN_CHECK(exhausted);
      check_equal(exhaustion_before, state);
    }
    bool maximum_source_stale = false;
    try {
      static_cast<void>(maximum_source.report());
    } catch (const hundun::runtime::Error &error) {
      maximum_source_stale = std::string_view(error.what()) ==
                             "material flow diagnostic source is stale";
    }
    HUNDUN_CHECK(maximum_source_stale);
    check_equal(exhaustion_before, state);
  });
}

// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/diag_time_control.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/fvm_cell_centered.hpp"
#include "hundun/flow_adaptive_time_control.hpp"
#include "hundun/lin_conjugate_gradient.hpp"
#include "hundun/lin_preconditioners.hpp"
#include "hundun/mesh_geometry.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_registry.hpp"
#include "hundun/rt_halo_exchange.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "hundun/rt_mpi_operation_error.hpp"
#include "hundun/rt_structured_decomposition.hpp"
#include "tests/support/flow_adaptive_time_control_test_access.hpp"
#include "tests/support/flow_constant_density_piso_test_access.hpp"
#include "tests/support/flow_ideal_gas_closure_test_access.hpp"
#include "tests/support/flow_material_density_piso_test_access.hpp"
#include "tests/support/lin_preconditioners_test_access.hpp"
#include "tests/support/flow_state_equality.hpp"
#include "tests/support/adaptive_time_control_test_support.hpp"
#include "tests/support/allocation_attempt_guard.hpp"
#include "tests/support/test_main.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <string>
#include <string_view>
#include <utility>

namespace {

class OneRecordSink final : public hundun::diagnostics::DiagnosticSink {
public:
  void submit(const hundun::diagnostics::DiagnosticRecord &value) override {
    record = value;
    ++calls;
  }

  hundun::diagnostics::DiagnosticRecord record;
  int calls{};
};

void require_zero_trusted_tail_events() {
  const auto observed =
      hundun::flow::test::AdaptiveTimeControlTestAccess::
          trusted_tail_observation();
  HUNDUN_CHECK(observed.allocation_attempts == 0U);
  HUNDUN_CHECK(observed.controller_collectives == 0U);
  HUNDUN_CHECK(observed.field_state_traversals == 0U);
  HUNDUN_CHECK(observed.callbacks_or_sinks == 0U);
}

std::uint64_t fp64_bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

hundun::flow::MomentumTimeStencil independent_stencil(
    hundun::flow::MomentumTimeOrder order, double dt,
    double previous_dt) {
  if (order == hundun::flow::MomentumTimeOrder::backward_euler)
    return {order, dt, 0.0, 1.0, -1.0, 0.0};
  const double ratio = dt / previous_dt;
  return {order, dt, previous_dt, (1.0 + 2.0 * ratio) / (1.0 + ratio),
          -(1.0 + ratio), ratio * ratio / (1.0 + ratio)};
}

void require_stencil_bitwise_equal(
    const hundun::flow::MomentumTimeStencil &actual,
    const hundun::flow::MomentumTimeStencil &expected) {
  HUNDUN_CHECK(actual.order == expected.order);
  HUNDUN_CHECK(fp64_bits(actual.dt_s) == fp64_bits(expected.dt_s));
  HUNDUN_CHECK(fp64_bits(actual.previous_dt_s) ==
               fp64_bits(expected.previous_dt_s));
  HUNDUN_CHECK(fp64_bits(actual.alpha0) == fp64_bits(expected.alpha0));
  HUNDUN_CHECK(fp64_bits(actual.alpha1) == fp64_bits(expected.alpha1));
  HUNDUN_CHECK(fp64_bits(actual.alpha2) == fp64_bits(expected.alpha2));
}

template <class Snapshot>
bool facade_cache_exact(const Snapshot &left, const Snapshot &right) noexcept {
  if (left.operator_count != right.operator_count ||
      left.delegated != right.delegated ||
      left.workspaces.size() != right.workspaces.size())
    return false;
  for (std::size_t index = 0; index < left.workspaces.size(); ++index) {
    if (left.workspaces[index].identity != right.workspaces[index].identity ||
        left.workspaces[index].capacity != right.workspaces[index].capacity)
      return false;
  }
  for (std::size_t index = 0; index < left.operator_count; ++index) {
    if (left.operators[index].identity != right.operators[index].identity ||
        left.operators[index].revision != right.operators[index].revision ||
        left.operators[index].diagonal.size() !=
            right.operators[index].diagonal.size())
      return false;
    for (std::size_t value = 0;
         value < left.operators[index].diagonal.size(); ++value) {
      if (fp64_bits(left.operators[index].diagonal[value]) !=
          fp64_bits(right.operators[index].diagonal[value]))
        return false;
    }
  }
  return true;
}

template <class Snapshot>
bool facade_cache_allocation_exact(const Snapshot &left,
                                   const Snapshot &right) noexcept {
  if (left.operator_count != right.operator_count ||
      left.delegated != right.delegated ||
      left.workspaces.size() != right.workspaces.size())
    return false;
  for (std::size_t index = 0; index < left.workspaces.size(); ++index) {
    if (left.workspaces[index].identity != right.workspaces[index].identity ||
        left.workspaces[index].capacity != right.workspaces[index].capacity)
      return false;
  }
  for (std::size_t index = 0; index < left.operator_count; ++index) {
    if (left.operators[index].identity != right.operators[index].identity)
      return false;
  }
  return true;
}

template <class Left, class Right>
bool facade_cache_delegation_exact(const Left &delegated,
                                   const Right &owned) noexcept {
  if (!delegated.delegated || owned.delegated ||
      delegated.operator_count != owned.operator_count ||
      delegated.workspaces.size() != owned.workspaces.size())
    return false;
  for (std::size_t index = 0; index < delegated.workspaces.size(); ++index) {
    if (delegated.workspaces[index].identity !=
            owned.workspaces[index].identity ||
        delegated.workspaces[index].capacity !=
            owned.workspaces[index].capacity)
      return false;
  }
  for (std::size_t index = 0; index < delegated.operator_count; ++index) {
    if (delegated.operators[index].identity !=
            owned.operators[index].identity ||
        delegated.operators[index].revision !=
            owned.operators[index].revision ||
        delegated.operators[index].diagonal.size() !=
            owned.operators[index].diagonal.size())
      return false;
    for (std::size_t value = 0;
         value < delegated.operators[index].diagonal.size(); ++value) {
      if (fp64_bits(delegated.operators[index].diagonal[value]) !=
          fp64_bits(owned.operators[index].diagonal[value]))
        return false;
    }
  }
  return true;
}

template <class Snapshot>
void require_facade_cache_inventory(const Snapshot &snapshot,
                                    bool delegated) {
  HUNDUN_CHECK(snapshot.delegated == delegated);
  HUNDUN_CHECK(!snapshot.workspaces.empty());
  HUNDUN_CHECK(snapshot.operator_count == 3U);
  for (std::size_t index = 0; index < snapshot.operator_count; ++index) {
    HUNDUN_CHECK(snapshot.operators[index].identity != 0U);
    HUNDUN_CHECK(snapshot.operators[index].revision != 0U);
  }
}

template <class Snapshot>
void require_facade_cache_mutation_sensitive(const Snapshot &snapshot) {
  HUNDUN_CHECK(facade_cache_exact(snapshot, snapshot));
  {
    auto changed = snapshot;
    ++changed.workspaces.front().identity;
    HUNDUN_CHECK(!facade_cache_exact(snapshot, changed));
  }
  {
    auto changed = snapshot;
    ++changed.workspaces.front().capacity;
    HUNDUN_CHECK(!facade_cache_exact(snapshot, changed));
  }
  {
    auto changed = snapshot;
    ++changed.operators.front().revision;
    HUNDUN_CHECK(!facade_cache_exact(snapshot, changed));
  }
  {
    auto changed = snapshot;
    ++changed.operators.front().identity;
    HUNDUN_CHECK(!facade_cache_exact(snapshot, changed));
  }
  {
    auto changed = snapshot;
    HUNDUN_CHECK(!changed.operators.front().diagonal.empty());
    changed.operators.front().diagonal.front() =
        std::nextafter(changed.operators.front().diagonal.front(),
                       std::numeric_limits<double>::infinity());
    HUNDUN_CHECK(!facade_cache_exact(snapshot, changed));
  }
}

struct JacobiCacheAuthority final {
  std::array<hundun::linear::test::JacobiStorageSnapshot, 3> storage{};
};

JacobiCacheAuthority capture_jacobi_cache(
    std::array<const hundun::linear::JacobiPreconditioner *, 3>
        preconditioners) {
  JacobiCacheAuthority result;
  for (std::size_t component = 0; component < preconditioners.size();
       ++component) {
    result.storage[component] =
        hundun::linear::test::PreconditionerTestAccess::jacobi_storage(
            *preconditioners[component]);
  }
  return result;
}

bool jacobi_allocation_exact(const JacobiCacheAuthority &left,
                             const JacobiCacheAuthority &right) noexcept {
  for (std::size_t component = 0; component < left.storage.size();
       ++component) {
    std::array<std::pair<hundun::execution::AllocationIdentity, std::size_t>,
               2>
        left_buffers{}, right_buffers{};
    for (std::size_t slot = 0; slot < left_buffers.size(); ++slot) {
      left_buffers[slot] = {
          left.storage[component].allocation_identities[slot],
          left.storage[component].byte_sizes[slot]};
      right_buffers[slot] = {
          right.storage[component].allocation_identities[slot],
          right.storage[component].byte_sizes[slot]};
    }
    std::sort(left_buffers.begin(), left_buffers.end());
    std::sort(right_buffers.begin(), right_buffers.end());
    if (left_buffers != right_buffers)
      return false;
  }
  return true;
}

bool jacobi_authority_exact(const JacobiCacheAuthority &left,
                            const JacobiCacheAuthority &right) noexcept {
  if (!jacobi_allocation_exact(left, right))
    return false;
  for (std::size_t component = 0; component < left.storage.size();
       ++component) {
    if (left.storage[component].revision !=
            right.storage[component].revision ||
        left.storage[component].cache_valid !=
            right.storage[component].cache_valid ||
        left.storage[component].cached_inverse.size() !=
            right.storage[component].cached_inverse.size())
      return false;
    for (std::size_t index = 0;
         index < left.storage[component].cached_inverse.size(); ++index) {
      if (fp64_bits(left.storage[component].cached_inverse[index]) !=
          fp64_bits(right.storage[component].cached_inverse[index]))
        return false;
    }
  }
  return true;
}

template <class Snapshot>
void require_jacobi_operator_coherence(
    const JacobiCacheAuthority &jacobi, const Snapshot &facade) {
  HUNDUN_CHECK(facade.operator_count == jacobi.storage.size());
  for (std::size_t component = 0; component < facade.operator_count;
       ++component) {
    HUNDUN_CHECK(jacobi.storage[component].cache_valid);
    HUNDUN_CHECK(jacobi.storage[component].revision ==
                 facade.operators[component].revision);
    HUNDUN_CHECK(jacobi.storage[component].allocation_identities[0] != 0U);
    HUNDUN_CHECK(jacobi.storage[component].allocation_identities[1] != 0U);
    HUNDUN_CHECK(jacobi.storage[component].byte_sizes[0] != 0U);
    HUNDUN_CHECK(jacobi.storage[component].byte_sizes[1] != 0U);
    HUNDUN_CHECK(jacobi.storage[component].cached_inverse.size() ==
                 facade.operators[component].diagonal.size());
    for (std::size_t index = 0;
         index < jacobi.storage[component].cached_inverse.size(); ++index) {
      const double expected =
          1.0 / facade.operators[component].diagonal[index];
      HUNDUN_CHECK(
          fp64_bits(jacobi.storage[component].cached_inverse[index]) ==
          fp64_bits(expected));
    }
  }
}

void require_jacobi_mutation_sensitive(const JacobiCacheAuthority &authority) {
  HUNDUN_CHECK(jacobi_authority_exact(authority, authority));
  {
    auto changed = authority;
    changed.storage.front().cache_valid =
        !changed.storage.front().cache_valid;
    HUNDUN_CHECK(!jacobi_authority_exact(authority, changed));
  }
  {
    auto changed = authority;
    ++changed.storage.front().allocation_identities.front();
    HUNDUN_CHECK(!jacobi_authority_exact(authority, changed));
  }
  {
    auto changed = authority;
    ++changed.storage.front().byte_sizes.front();
    HUNDUN_CHECK(!jacobi_authority_exact(authority, changed));
  }
  {
    auto changed = authority;
    ++changed.storage.front().revision;
    HUNDUN_CHECK(!jacobi_authority_exact(authority, changed));
  }
  {
    auto changed = authority;
    if (changed.storage.front().cached_inverse.empty()) {
      changed.storage.front().cached_inverse.push_back(1.0);
    } else {
      changed.storage.front().cached_inverse.front() =
          std::nextafter(changed.storage.front().cached_inverse.front(),
                         std::numeric_limits<double>::infinity());
    }
    HUNDUN_CHECK(!jacobi_authority_exact(authority, changed));
  }
}

void require_cold_jacobi_unpublished(
    const JacobiCacheAuthority &authority) {
  for (const auto &storage : authority.storage) {
    HUNDUN_CHECK(!storage.cache_valid);
    HUNDUN_CHECK(storage.revision == 0U);
    HUNDUN_CHECK(storage.allocation_identities[0] == 0U);
    HUNDUN_CHECK(storage.allocation_identities[1] == 0U);
    HUNDUN_CHECK(storage.byte_sizes[0] == 0U);
    HUNDUN_CHECK(storage.byte_sizes[1] == 0U);
    HUNDUN_CHECK(storage.cached_inverse.empty());
  }
}

class JacobiColdFaultReset final {
 public:
  JacobiColdFaultReset() noexcept {
    hundun::linear::test::PreconditionerTestAccess::
        reset_cold_update_fault();
  }

  ~JacobiColdFaultReset() noexcept {
    hundun::linear::test::PreconditionerTestAccess::
        reset_cold_update_fault();
  }

  JacobiColdFaultReset(const JacobiColdFaultReset &) = delete;
  JacobiColdFaultReset &operator=(const JacobiColdFaultReset &) = delete;
};

struct AttemptBoundaryEvidence final {
  hundun::flow::Bdf2RetryController *controller{};
  std::array<hundun::flow::MomentumTimeStencil, 9> stencils{};
  std::array<hundun::test::AdaptiveFlowStateSnapshot, 9> states{};
  std::size_t count{};

  static void observe(void *context,
                      const hundun::flow::MomentumTimeStencil &stencil,
                      const hundun::flow::FlowState &state) {
    auto &self = *static_cast<AttemptBoundaryEvidence *>(context);
    HUNDUN_CHECK(self.count < self.stencils.size());
    self.stencils[self.count] = stencil;
    self.states[self.count] = hundun::test::capture_adaptive_flow_state(
        state, self.controller->state());
    ++self.count;
  }
};

void require_authentic_diagnostic(
    const hundun::runtime::MpiContext &mpi,
    const hundun::flow::Bdf2RetryController &controller,
    const hundun::flow::FlowState &state,
    const hundun::flow::TimeAdvanceReport &report,
    hundun::diagnostics::DiagnosticStatus expected_status,
    hundun::diagnostics::DiagnosticFailureClass expected_classification,
    std::string_view expected_code, int expected_collective_rank) {
  auto source = controller.diagnostic_source(state, report);
  for (const auto scope :
       {hundun::diagnostics::DiagnosticScope::local,
        hundun::diagnostics::DiagnosticScope::collective}) {
    OneRecordSink sink;
    const hundun::diagnostics::DiagnosticRequest request{
        hundun::diagnostics::DiagnosticLevel::summary, scope,
        {mpi.rank(), state.metadata().step, state.metadata().time_s,
         "time-control.advance-result"},
        {}, 0U};
    if (scope == hundun::diagnostics::DiagnosticScope::local)
      hundun::diagnostics::collect_diagnostics(source, request, sink);
    else
      hundun::diagnostics::collect_diagnostics(source, mpi, request, sink);
    HUNDUN_CHECK(sink.calls == 1);
    HUNDUN_CHECK(sink.record.status == expected_status);
    HUNDUN_CHECK(sink.record.failure.classification ==
                 expected_classification);
    HUNDUN_CHECK(sink.record.failure.code == expected_code);
    const int expected_rank =
        scope == hundun::diagnostics::DiagnosticScope::collective
            ? expected_collective_rank
            : -1;
    if (sink.record.failure.lowest_failing_rank != expected_rank)
      throw std::runtime_error(
          std::string("authentic diagnostic rank mismatch for ") +
          std::string(expected_code) + ": observed " +
          std::to_string(sink.record.failure.lowest_failing_rank) +
          ", expected " + std::to_string(expected_rank));
  }
}

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

class AlwaysFailureSolver final : public hundun::linear::LinearSolver {
public:
  explicit AlwaysFailureSolver(int lowest_failing_rank)
      : lowest_failing_rank_(lowest_failing_rank) {}

  hundun::linear::SolveReport
  solve(const hundun::linear::LinearOperator &,
        hundun::linear::Preconditioner &,
        hundun::execution::VectorView<const double>,
        hundun::execution::VectorView<double>,
        const hundun::linear::SolveControl &) const override {
    hundun::linear::SolveReport report;
    report.reason =
        hundun::linear::SolveTerminationReason::maximum_iterations;
    report.lowest_failing_rank = lowest_failing_rank_;
    return report;
  }

private:
  int lowest_failing_rank_{};
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
  const int mismatch_rank = mpi.size() == 1 ? 0 : 1;
  const hundun::config::FlowTimeConfig time{
      hundun::config::TimeMode::adaptive, 2, 0.01, 0.00125, 0.02,
      0.5, 0.25, 1.25, 0.5, 8};
  const auto make_layout_state =
      [&](hundun::runtime::FieldLayoutSet layout) {
        auto candidate = hundun::flow::FlowState::create(
            registry, layout, fields,
            {0U, 0.0, 0.01, 0.0,
             hundun::flow::MomentumTimeOrder::backward_euler});
        const auto candidate_cells =
            static_cast<std::size_t>(layout.cell_interior_extent.x) *
            static_cast<std::size_t>(layout.cell_interior_extent.y) *
            static_cast<std::size_t>(layout.cell_interior_extent.z);
        hundun::flow::FlowLayerValues values;
        values.density.assign(candidate_cells, 1.0);
        values.velocity.assign(candidate_cells * 3U, 0.0);
        values.mechanical_pressure.assign(candidate_cells, 0.0);
        values.face_velocity.assign(layout.face_count * 3U, 0.0);
        values.face_mass_flux.assign(layout.face_count, 0.0);
        values.transported_cell_fields = {
            std::vector<double>(candidate_cells, 1.0)};
        candidate.seed_accepted_layers(values, values);
        return candidate;
      };
  {
    const auto correct_layout =
        hundun::runtime::FieldLayoutSet{local, topology.local_face_count()};
    auto create_positive = hundun::flow::Bdf2RetryController::create(
        time, hundun::config::DensityModel::constant, topology, geometry, mpi,
        state);
    auto restore_positive = hundun::flow::Bdf2RetryController::restore(
        time, hundun::config::DensityModel::constant, topology, geometry, mpi,
        state, create_positive.state());
    HUNDUN_CHECK(restore_positive.state().state_seal ==
                 create_positive.state().state_seal);

    const auto require_layout_rejection =
        [&](bool restore, bool cell_extent_mismatch) {
          auto layout = correct_layout;
          if (mpi.rank() == mismatch_rank) {
            if (cell_extent_mismatch)
              ++layout.cell_interior_extent.x;
            else
              ++layout.face_count;
          }
          auto candidate = make_layout_state(layout);
          const auto before = hundun::test::capture_adaptive_flow_state(
              candidate, {});
          bool rejected = false;
          try {
            if (restore) {
              static_cast<void>(hundun::flow::Bdf2RetryController::restore(
                  time, hundun::config::DensityModel::constant, topology,
                  geometry, mpi, candidate, create_positive.state()));
            } else {
              static_cast<void>(hundun::flow::Bdf2RetryController::create(
                  time, hundun::config::DensityModel::constant, topology,
                  geometry, mpi, candidate));
            }
          } catch (const hundun::runtime::Error &error) {
            rejected = true;
            const std::string message(error.what());
            HUNDUN_CHECK(
                message.find(restore ? "time-control.restore.layout"
                                     : "time-control.create.layout") !=
                std::string::npos);
            HUNDUN_CHECK(
                message.find("lowest failing rank " +
                             std::to_string(mismatch_rank)) !=
                std::string::npos);
          }
          HUNDUN_CHECK(rejected);
          const auto after = hundun::test::capture_adaptive_flow_state(
              candidate, {});
          HUNDUN_CHECK(hundun::test::adaptive_flow_state_bitwise_equal(
              before, after));
        };
    require_layout_rejection(true, true);
    require_layout_rejection(true, false);
    require_layout_rejection(false, true);
    require_layout_rejection(false, false);
  }
  {
    auto cold_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    cold_state.seed_accepted_layers(initial, initial);
    hundun::linear::JacobiPreconditioner cold_mx(execution),
        cold_my(execution), cold_mz(execution),
        cold_pressure(execution);
    auto cold_facade = hundun::flow::FixedStepConstantDensityFlow::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        momentum_solver, {&cold_mx, &cold_my, &cold_mz}, pressure_solver,
        cold_pressure,
        {{fields.transported_cell_fields.front(),
          hundun::finite_volume::FiniteVolumeQuantity::scalar(0U), 0.0}});
    using ConstantAccess =
        hundun::flow::test::ConstantDensityPisoTestAccess;
    ConstantAccess::reset();
    const auto cold_cache_before =
        ConstantAccess::facade_cache_snapshot(cold_facade);
    const auto cold_jacobi_before =
        capture_jacobi_cache({&cold_mx, &cold_my, &cold_mz});
    require_cold_jacobi_unpublished(cold_jacobi_before);
    require_jacobi_mutation_sensitive(cold_jacobi_before);
    JacobiColdFaultReset cold_fault_reset;
    hundun::linear::test::PreconditionerTestAccess::
        arm_fail_next_cold_update_before_publication();
    const auto cold_report = cold_facade.attempt(
        cold_state, 1.0, 0.0,
        independent_stencil(
            hundun::flow::MomentumTimeOrder::backward_euler, 0.01, 0.0),
        {}, {});
    hundun::linear::test::PreconditionerTestAccess::
        reset_cold_update_fault();
    HUNDUN_CHECK(cold_report.disposition ==
                 hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(cold_report.reason ==
                 hundun::flow::StepFailureReason::momentum_linear_solve);
    HUNDUN_CHECK(cold_report.lowest_failing_rank == 0);
    const auto cold_cache = ConstantAccess::facade_cache_snapshot(cold_facade);
    const auto cold_jacobi =
        capture_jacobi_cache({&cold_mx, &cold_my, &cold_mz});
    HUNDUN_CHECK(cold_cache.operator_count == 3U);
    HUNDUN_CHECK(cold_cache.operators[0].revision >
                 cold_cache_before.operators[0].revision);
    require_cold_jacobi_unpublished(cold_jacobi);
    require_jacobi_mutation_sensitive(cold_jacobi);
    const auto cold_retry = cold_facade.attempt(
        cold_state, 1.0, 0.0,
        independent_stencil(
            hundun::flow::MomentumTimeOrder::backward_euler, 0.01, 0.0),
        {}, {});
    HUNDUN_CHECK(cold_retry.disposition ==
                 hundun::flow::StepAttemptDisposition::committed);
    const auto cold_cache_published =
        ConstantAccess::facade_cache_snapshot(cold_facade);
    const auto cold_jacobi_published =
        capture_jacobi_cache({&cold_mx, &cold_my, &cold_mz});
    for (std::size_t component = 0; component < 3U; ++component) {
      HUNDUN_CHECK(cold_cache_published.operators[component].revision >
                   cold_cache.operators[component].revision);
    }
    require_jacobi_operator_coherence(cold_jacobi_published,
                                      cold_cache_published);
    require_jacobi_mutation_sensitive(cold_jacobi_published);
  }
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
  {
    bool rejected = false;
    try {
      hundun::test::allocation_probe::FailNextAllocationGuard fail(
          mpi.rank() == mismatch_rank);
      static_cast<void>(hundun::flow::Bdf2RetryController::create(
          time, hundun::config::DensityModel::constant, topology, geometry,
          mpi, state));
    } catch (const hundun::runtime::Error &error) {
      rejected = true;
      HUNDUN_CHECK(std::string(error.what()).find(
                       "time-control.create.state") != std::string::npos);
      HUNDUN_CHECK(std::string(error.what()).find(
                       std::to_string(mismatch_rank)) != std::string::npos);
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(state.metadata().step == 0U);
  }
  for (std::size_t ordinal = 1U; ordinal <= 7U; ++ordinal) {
    using TimeAccess = hundun::flow::test::AdaptiveTimeControlTestAccess;
    const int origin_rank =
        mpi.size() == 4 && ordinal % 2U != 0U ? 2 : mismatch_rank;
    TimeAccess::reset_faults();
    TimeAccess::set_raw_fault(ordinal, origin_rank);
    bool typed = false;
    try {
      static_cast<void>(hundun::flow::Bdf2RetryController::create(
          time, hundun::config::DensityModel::constant, topology, geometry,
          mpi, state));
    } catch (const hundun::runtime::MpiOperationError &) {
      typed = true;
    }
    HUNDUN_CHECK(typed);
    const auto observed = TimeAccess::raw_fault_observation();
    HUNDUN_CHECK(observed.raw_operations == ordinal);
    HUNDUN_CHECK(observed.fault_ordinal == ordinal);
    HUNDUN_CHECK(observed.requested_rank == origin_rank);
    HUNDUN_CHECK(observed.local_origins ==
                 (mpi.rank() == origin_rank ? 1U : 0U));
    std::uint64_t origin_count = observed.local_origins;
    HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, &origin_count, 1, MPI_UINT64_T,
                               MPI_SUM, mpi.comm()) == MPI_SUCCESS);
    HUNDUN_CHECK(origin_count == 1U);
    HUNDUN_CHECK(state.metadata().step == 0U);
  }
  hundun::flow::test::AdaptiveTimeControlTestAccess::reset_faults();
  auto controller = hundun::flow::Bdf2RetryController::create(
      time, hundun::config::DensityModel::constant, topology, geometry, mpi,
      state);
  {
    const auto rejected =
        controller.advance(state, facade, 1.0, -1.0, {}, {});
    HUNDUN_CHECK(
        hundun::flow::test::AdaptiveTimeControlTestAccess::
            preflight_category(rejected) == 1U);
    require_authentic_diagnostic(
        mpi, controller, state, rejected,
        hundun::diagnostics::DiagnosticStatus::failed,
        hundun::diagnostics::DiagnosticFailureClass::invalid_input,
        "time-control.preflight.config", 0);
  }
  {
    using TimeAccess = hundun::flow::test::AdaptiveTimeControlTestAccess;
    const auto before = state.metadata();
    hundun::test::allocation_probe::FailNextAllocationGuard fail(
        mpi.rank() == mismatch_rank);
    const auto rejected =
        controller.advance(state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(
        rejected.disposition() ==
        hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
    HUNDUN_CHECK(rejected.attempt_count() == 0U);
    HUNDUN_CHECK(TimeAccess::preflight_category(rejected) == 7U);
    HUNDUN_CHECK(rejected.lowest_failing_rank() == mismatch_rank);
    HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
        before, state.metadata()));
  }
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
    require_authentic_diagnostic(
        mpi, wrong_model, state, rejected,
        hundun::diagnostics::DiagnosticStatus::failed,
        hundun::diagnostics::DiagnosticFailureClass::invalid_input,
        "time-control.preflight.identity", 0);
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
      hundun::diagnostics::DiagnosticFailureClass classification;
      const char *code;
    };
    constexpr std::array<Case, 3> cases{{
        {Fault::layout, 3U,
         hundun::diagnostics::DiagnosticFailureClass::layout,
         "time-control.preflight.layout"},
        {Fault::capability, 4U,
         hundun::diagnostics::DiagnosticFailureClass::capability,
         "time-control.preflight.capability"},
        {Fault::preparation, 7U,
         hundun::diagnostics::DiagnosticFailureClass::invalid_input,
         "time-control.preflight.preparation"},
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
      require_authentic_diagnostic(
          mpi, controller, state, rejected,
          hundun::diagnostics::DiagnosticStatus::failed,
          candidate.classification, candidate.code, mismatch_rank);
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
    auto source =
        malformed_controller.diagnostic_source(malformed_state, malformed);
    OneRecordSink sink;
    hundun::diagnostics::collect_diagnostics(
        source,
        {hundun::diagnostics::DiagnosticLevel::summary,
         hundun::diagnostics::DiagnosticScope::local,
         {mpi.rank(), malformed_state.metadata().step,
          malformed_state.metadata().time_s, "time-control.advance-result"},
         {}, 0U},
        sink);
    HUNDUN_CHECK(sink.calls == 1);
    HUNDUN_CHECK(sink.record.failure.code ==
                 "time-control.preflight.report");
    require_authentic_diagnostic(
        mpi, malformed_controller, malformed_state, malformed,
        hundun::diagnostics::DiagnosticStatus::failed,
        hundun::diagnostics::DiagnosticFailureClass::invalid_input,
        "time-control.preflight.report", mismatch_rank);
    TimeAccess::reset_faults();
  }
  {
    using TimeAccess = hundun::flow::test::AdaptiveTimeControlTestAccess;
    struct ReasonCase final {
      hundun::flow::StepFailureReason reason;
      hundun::diagnostics::DiagnosticFailureClass classification;
      const char *code;
    };
    constexpr std::array<ReasonCase, 6> reasons{{
        {hundun::flow::StepFailureReason::non_finite_trial,
         hundun::diagnostics::DiagnosticFailureClass::non_finite_state,
         "time-control.non-finite-trial"},
        {hundun::flow::StepFailureReason::final_momentum_residual,
         hundun::diagnostics::DiagnosticFailureClass::non_convergence,
         "time-control.final-momentum-residual"},
        {hundun::flow::StepFailureReason::final_transport_residual,
         hundun::diagnostics::DiagnosticFailureClass::non_convergence,
         "time-control.final-transport-residual"},
        {hundun::flow::StepFailureReason::final_conservation_defect,
         hundun::diagnostics::DiagnosticFailureClass::conservation,
         "time-control.final-conservation-defect"},
        {hundun::flow::StepFailureReason::final_continuity_residual,
         hundun::diagnostics::DiagnosticFailureClass::non_convergence,
         "time-control.final-continuity-residual"},
        {hundun::flow::StepFailureReason::final_pressure_residual,
         hundun::diagnostics::DiagnosticFailureClass::non_convergence,
         "time-control.final-pressure-residual"},
    }};
    for (const auto &reason_case : reasons) {
      auto scheduled_state = hundun::flow::FlowState::create(
          registry, {local, topology.local_face_count()}, fields,
          {0U, 0.0, 0.01, 0.0,
           hundun::flow::MomentumTimeOrder::backward_euler});
      scheduled_state.seed_accepted_layers(initial, initial);
      auto scheduled_controller =
          hundun::flow::Bdf2RetryController::create(
              time, hundun::config::DensityModel::constant, topology,
              geometry, mpi, scheduled_state);
      TimeAccess::set_recoverable_failure_reason(reason_case.reason);
      TimeAccess::set_recoverable_failures(1U, mismatch_rank);
      const auto scheduled = scheduled_controller.advance(
          scheduled_state, facade, 1.0, 0.0, {}, {});
      HUNDUN_CHECK(scheduled.disposition() ==
                   hundun::flow::TimeAdvanceDisposition::committed);
      HUNDUN_CHECK(scheduled.attempt_count() == 2U);
      HUNDUN_CHECK(scheduled.attempt(0).reason == reason_case.reason);
      HUNDUN_CHECK(
          scheduled.attempt(0).disposition ==
          hundun::flow::StepAttemptDisposition::recoverable_failure);
      HUNDUN_CHECK(
          scheduled.attempt(1).disposition ==
          hundun::flow::StepAttemptDisposition::committed);
      HUNDUN_CHECK(scheduled_state.metadata().step == 1U);
      TimeAccess::reset_faults();

      auto terminal_state = hundun::flow::FlowState::create(
          registry, {local, topology.local_face_count()}, fields,
          {0U, 0.0, 0.01, 0.0,
           hundun::flow::MomentumTimeOrder::backward_euler});
      terminal_state.seed_accepted_layers(initial, initial);
      auto terminal_controller =
          hundun::flow::Bdf2RetryController::create(
              time, hundun::config::DensityModel::constant, topology,
              geometry, mpi, terminal_state);
      TimeAccess::set_recoverable_failure_reason(reason_case.reason);
      TimeAccess::set_recoverable_failures(9U, mismatch_rank);
      const auto terminal = terminal_controller.advance(
          terminal_state, facade, 1.0, 0.0, {}, {});
      HUNDUN_CHECK(
          terminal.disposition() ==
          hundun::flow::TimeAdvanceDisposition::minimum_dt_failure);
      HUNDUN_CHECK(terminal.reason() == reason_case.reason);
      const bool globally_observed =
          reason_case.reason ==
              hundun::flow::StepFailureReason::non_finite_trial ||
          reason_case.reason ==
              hundun::flow::StepFailureReason::final_momentum_residual ||
          reason_case.reason ==
              hundun::flow::StepFailureReason::final_transport_residual;
      require_authentic_diagnostic(
          mpi, terminal_controller, terminal_state, terminal,
          hundun::diagnostics::DiagnosticStatus::failed,
          reason_case.classification, reason_case.code,
          globally_observed ? 0 : mismatch_rank);
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
    const auto run_terminal_linear_failure =
        [&](bool momentum_failure,
            hundun::flow::StepFailureReason expected_reason,
            const char *expected_code) {
          AlwaysFailureSolver always_fails(mismatch_rank);
          hundun::linear::JacobiPreconditioner fx(execution), fy(execution),
              fz(execution), fp(execution);
          const auto &selected_momentum =
              momentum_failure
                  ? static_cast<const hundun::linear::LinearSolver &>(
                        always_fails)
                  : static_cast<const hundun::linear::LinearSolver &>(
                        momentum_solver);
          const auto &selected_pressure =
              momentum_failure
                  ? static_cast<const hundun::linear::LinearSolver &>(
                        pressure_solver)
                  : static_cast<const hundun::linear::LinearSolver &>(
                        always_fails);
          auto terminal_flow =
              hundun::flow::FixedStepConstantDensityFlow::create(
                  decomposition, topology, geometry, boundaries, mpi,
                  execution, halo, selected_momentum, {&fx, &fy, &fz},
                  selected_pressure, fp,
                  {{fields.transported_cell_fields.front(),
                    hundun::finite_volume::FiniteVolumeQuantity::scalar(0U),
                    0.0}});
          auto terminal_state = hundun::flow::FlowState::create(
              registry, {local, topology.local_face_count()}, fields,
              {0U, 0.0, 0.01, 0.0,
               hundun::flow::MomentumTimeOrder::backward_euler});
          terminal_state.seed_accepted_layers(initial, initial);
          auto terminal_controller =
              hundun::flow::Bdf2RetryController::create(
                  time, hundun::config::DensityModel::constant, topology,
                  geometry, mpi, terminal_state);
          const auto result = terminal_controller.advance(
              terminal_state, terminal_flow, 1.0, 0.0, {}, {});
          HUNDUN_CHECK(
              result.disposition() ==
              hundun::flow::TimeAdvanceDisposition::minimum_dt_failure);
          HUNDUN_CHECK(result.reason() == expected_reason);
          require_authentic_diagnostic(
              mpi, terminal_controller, terminal_state, result,
              hundun::diagnostics::DiagnosticStatus::failed,
              hundun::diagnostics::DiagnosticFailureClass::non_convergence,
              expected_code, mismatch_rank);
        };
    run_terminal_linear_failure(
        true, hundun::flow::StepFailureReason::momentum_linear_solve,
        "time-control.momentum-linear-solve");
    run_terminal_linear_failure(
        false, hundun::flow::StepFailureReason::pressure_linear_solve,
        "time-control.pressure-linear-solve");
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

    const auto run_directional_stability =
        [&](const hundun::mesh::MeshGeometry &case_geometry,
            hundun::flow::FixedStepConstantDensityFlow *case_facade,
            std::size_t enabled_directions) {
          auto directional_initial = initial;
          for (std::size_t face_id = 0;
               face_id < directional_initial.face_mass_flux.size();
               ++face_id) {
            const auto area_vector = case_geometry.face_area_vector_m2(
                face_id, hundun::mesh::FaceSide::owner);
            const std::array<double, 3> magnitudes{
                std::abs(area_vector.x), std::abs(area_vector.y),
                std::abs(area_vector.z)};
            const auto dominant = static_cast<std::size_t>(
                std::max_element(magnitudes.begin(), magnitudes.end()) -
                magnitudes.begin());
            directional_initial.face_mass_flux[face_id] =
                dominant < enabled_directions
                    ? ((face_id % 2U == 0U) ? 1.0e-12 : -1.0e-12) *
                          static_cast<double>(dominant + 1U)
                    : 0.0;
          }
          auto directional_state = hundun::flow::FlowState::create(
              registry, {local, topology.local_face_count()}, fields,
              {0U, 0.0, 0.01, 0.0,
               hundun::flow::MomentumTimeOrder::backward_euler});
          directional_state.seed_accepted_layers(directional_initial,
                                                 directional_initial);
          std::vector<double> directional_flux(
              topology.owned_cell_count());
          std::vector<double> directional_geometry(
              topology.owned_cell_count());
          for (hundun::mesh::LocalFaceId face_id = 0;
               face_id < topology.local_face_count(); ++face_id) {
            const auto area = case_geometry.face_area_m2(face_id);
            const auto displacement =
                case_geometry.face_displacement_m(face_id);
            const auto area_vector = case_geometry.face_area_vector_m2(
                face_id, hundun::mesh::FaceSide::owner);
            const auto factor =
                area * area /
                std::abs(area_vector.x * displacement.x +
                         area_vector.y * displacement.y +
                         area_vector.z * displacement.z);
            const auto add = [&](std::size_t cell_id) {
              directional_flux[cell_id] +=
                  std::abs(directional_initial.face_mass_flux[face_id]);
              directional_geometry[cell_id] += factor;
            };
            const auto owner = topology.owner(face_id);
            if (owner < topology.owned_cell_count())
              add(owner);
            const auto neighbour = topology.neighbour(face_id);
            if (neighbour && *neighbour < topology.owned_cell_count())
              add(*neighbour);
          }
          double directional_expected[2]{};
          for (std::size_t cell_id = 0;
               cell_id < topology.owned_cell_count(); ++cell_id) {
            const auto denominator =
                4.0 * case_geometry.cell_volume_m3(cell_id);
            directional_expected[0] =
                std::max(directional_expected[0],
                         directional_flux[cell_id] / denominator);
            directional_expected[1] =
                std::max(directional_expected[1],
                         0.1 * directional_geometry[cell_id] / denominator);
          }
          if (case_facade != nullptr)
            HUNDUN_CHECK(MPI_Allreduce(
                             MPI_IN_PLACE, directional_expected, 2,
                             MPI_DOUBLE, MPI_MAX, mpi.comm()) == MPI_SUCCESS);
          auto directional_controller =
              hundun::flow::Bdf2RetryController::create(
                  time, hundun::config::DensityModel::constant, topology,
                  case_geometry, mpi, directional_state);
          if (case_facade != nullptr) {
            const auto directional_report = directional_controller.advance(
                directional_state, *case_facade, 2.0, 0.1, {}, {});
            HUNDUN_CHECK(directional_report.stability_metrics_available());
            HUNDUN_CHECK(directional_report.convective_rate_per_s() ==
                         directional_expected[0]);
            HUNDUN_CHECK(directional_report.diffusive_rate_per_s() ==
                         directional_expected[1]);
          } else {
            const auto rates = hundun::flow::test::
                AdaptiveTimeControlTestAccess::stability_rates(
                    directional_controller, directional_state, 2.0, false,
                    0.1);
            HUNDUN_CHECK(rates[0] == directional_expected[0]);
            HUNDUN_CHECK(rates[1] == directional_expected[1]);
          }
        };
    run_directional_stability(geometry, &facade, 1U);
    run_directional_stability(geometry, &facade, 2U);
    run_directional_stability(geometry, &facade, 3U);

    const hundun::mesh::MeshGeometry warped_geometry(
        topology, hundun::mesh::AnalyticWarpedBoxMapping(
                      {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0},
                      {0.01, -0.008, 0.006}));
    run_directional_stability(warped_geometry, nullptr, 3U);

    for (const double invalid_flux :
         {std::numeric_limits<double>::quiet_NaN(),
          std::numeric_limits<double>::infinity(),
          -std::numeric_limits<double>::infinity()}) {
      auto invalid_state = hundun::flow::FlowState::create(
          registry, {local, topology.local_face_count()}, fields,
          {0U, 0.0, 0.01, 0.0,
           hundun::flow::MomentumTimeOrder::backward_euler});
      invalid_state.seed_accepted_layers(initial, initial);
      hundun::flow::test::MaterialDensityPisoTestAccess::
          set_accepted_face_mass_flux(invalid_state, 0U, invalid_flux);
      auto invalid_controller = hundun::flow::Bdf2RetryController::create(
          time, hundun::config::DensityModel::constant, topology, geometry,
          mpi, invalid_state);
      const auto invalid_report = invalid_controller.advance(
          invalid_state, facade, 1.0, 0.0, {}, {});
      HUNDUN_CHECK(
          invalid_report.disposition() ==
          hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
      HUNDUN_CHECK(invalid_report.attempt_count() == 0U);
      HUNDUN_CHECK(
          hundun::flow::test::AdaptiveTimeControlTestAccess::
              preflight_category(invalid_report) == 5U);
      require_authentic_diagnostic(
          mpi, invalid_controller, invalid_state, invalid_report,
          hundun::diagnostics::DiagnosticStatus::failed,
          hundun::diagnostics::DiagnosticFailureClass::invalid_input,
          "time-control.preflight.state", 0);
    }
  }
  {
    using TimeAccess = hundun::flow::test::AdaptiveTimeControlTestAccess;
    using Mutation = hundun::flow::test::TimeControlPostReturnMutation;
    constexpr std::array<Mutation, 5> mutations{
        Mutation::momentum_x_iterations, Mutation::momentum_y_iterations,
        Mutation::momentum_z_iterations, Mutation::pressure_one_iterations,
        Mutation::pressure_two_iterations};
    struct GateCase final {
      std::uint64_t limit;
      std::uint64_t iterations;
      bool expected;
    };
    constexpr std::array<GateCase, 8> gate_cases{{
        {1U, 0U, true},
        {1U, 1U, false},
        {2U, 1U, true},
        {2U, 2U, false},
        {3U, 1U, true},
        {3U, 2U, false},
        {std::numeric_limits<std::uint64_t>::max(),
         std::numeric_limits<std::uint64_t>::max() / 2U, true},
        {std::numeric_limits<std::uint64_t>::max(),
         std::numeric_limits<std::uint64_t>::max() / 2U + 1U, false},
    }};
    for (const auto mutation : mutations) {
      for (const auto gate : gate_cases) {
        auto gate_state = hundun::flow::FlowState::create(
            registry, {local, topology.local_face_count()}, fields,
            {0U, 0.0, 0.01, 0.0,
             hundun::flow::MomentumTimeOrder::backward_euler});
        gate_state.seed_accepted_layers(initial, initial);
        auto gate_controller = hundun::flow::Bdf2RetryController::create(
            time, hundun::config::DensityModel::constant, topology, geometry,
            mpi, gate_state);
        hundun::linear::SolveControl control;
        control.max_iterations = gate.limit;
        TimeAccess::set_post_return_mutation(mutation, mpi.rank());
        TimeAccess::set_post_return_iteration_value(gate.iterations);
        const auto report = gate_controller.advance(
            gate_state, facade, 1.0, 0.0, control, control);
        HUNDUN_CHECK(report.disposition() ==
                     hundun::flow::TimeAdvanceDisposition::committed);
        HUNDUN_CHECK(
            report.attempt(0).all_linear_solves_within_half_limit ==
            gate.expected);
        HUNDUN_CHECK(
            gate_controller.state()
                .last_all_linear_solves_within_half_limit == gate.expected);
        HUNDUN_CHECK(gate_controller.state().proposed_next_dt_s ==
                     (gate.expected ? 0.0125 : 0.01));
        TimeAccess::reset_faults();
      }
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
    hundun::runtime::FieldAccessPlan retained_view_plan(registry);
    retained_view_plan.declare_access(
        2201U, 1U, fields.velocity,
        hundun::runtime::AccessMode::read_write);
    retained_view_plan.freeze();
    auto retained_trial_view =
        limit_state.trial_layer().acquire_write<double>(
            retained_view_plan, 2201U, 1U, fields.velocity);
    const auto constant_cache_before =
        hundun::flow::test::ConstantDensityPisoTestAccess::
            facade_cache_snapshot(facade);
    const auto constant_jacobi_before =
        capture_jacobi_cache({&mx, &my, &mz});
    require_facade_cache_inventory(constant_cache_before, false);
    require_facade_cache_mutation_sensitive(constant_cache_before);
    require_jacobi_operator_coherence(constant_jacobi_before,
                                      constant_cache_before);
    require_jacobi_mutation_sensitive(constant_jacobi_before);
    const auto reductions_before = mpi.fp64_reduction_counters();
    const auto exact_before = hundun::test::capture_adaptive_flow_state(
        limit_state, limit_controller.state(), {}, 0U, 0U, 0U,
        {reductions_before.collective_calls,
         reductions_before.reduced_scalars,
         reductions_before.logical_payload_bytes, 0U});
    TimeAccess::set_recoverable_failure_reason(
        hundun::flow::StepFailureReason::final_continuity_residual);
    TimeAccess::set_recoverable_failures(9U, 0);
    const auto limited =
        limit_controller.advance(limit_state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(limited.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::retry_limit_reached);
    HUNDUN_CHECK(limited.attempt_count() == 9U);
    for (std::size_t index = 0; index < limited.attempt_count(); ++index)
      HUNDUN_CHECK(limited.attempt(index).disposition ==
                   hundun::flow::StepAttemptDisposition::recoverable_failure);
    HUNDUN_CHECK(limited.reason() ==
                 hundun::flow::StepFailureReason::final_continuity_residual);
    HUNDUN_CHECK(limit_state.metadata().step == 0U);
    require_authentic_diagnostic(
        mpi, limit_controller, limit_state, limited,
        hundun::diagnostics::DiagnosticStatus::failed,
        hundun::diagnostics::DiagnosticFailureClass::non_convergence,
        "time-control.final-continuity-residual", 0);
    const auto constant_cache_after =
        hundun::flow::test::ConstantDensityPisoTestAccess::
            facade_cache_snapshot(facade);
    const auto constant_jacobi_after =
        capture_jacobi_cache({&mx, &my, &mz});
    HUNDUN_CHECK(facade_cache_allocation_exact(constant_cache_before,
                                               constant_cache_after));
    HUNDUN_CHECK(jacobi_allocation_exact(constant_jacobi_before,
                                         constant_jacobi_after));
    for (std::size_t component = 0; component < 3U; ++component) {
      HUNDUN_CHECK(constant_cache_after.operators[component].revision >
                   constant_cache_before.operators[component].revision);
    }
    require_jacobi_operator_coherence(constant_jacobi_after,
                                      constant_cache_after);
    const auto reductions_after = mpi.fp64_reduction_counters();
    const auto exact_after = hundun::test::capture_adaptive_flow_state(
        limit_state, limit_controller.state(), {}, 0U, 0U, 0U,
        {reductions_after.collective_calls,
         reductions_after.reduced_scalars,
         reductions_after.logical_payload_bytes, 0U});
    HUNDUN_CHECK(hundun::test::adaptive_flow_state_failed_attempt_preserved(
        exact_before, exact_after, limited.attempt_count(),
        {1332U, 2052U, 16416U, 0U}));
    bool retained_view_stale = false;
    try {
      static_cast<void>(retained_trial_view(0, 0, 0, 0));
    } catch (const hundun::runtime::Error &) {
      retained_view_stale = true;
    }
    HUNDUN_CHECK(retained_view_stale);
    TimeAccess::reset_faults();
    const auto fixed_success =
        limit_controller.advance(limit_state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(fixed_success.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::committed);
    HUNDUN_CHECK(fixed_success.attempt_count() == 1U);
    HUNDUN_CHECK(fixed_success.attempt(0).attempted_dt_s == 0.01);
    HUNDUN_CHECK(limit_controller.state().proposed_next_dt_s == 0.01);
    const auto constant_cache_success =
        hundun::flow::test::ConstantDensityPisoTestAccess::
            facade_cache_snapshot(facade);
    const auto constant_jacobi_success =
        capture_jacobi_cache({&mx, &my, &mz});
    HUNDUN_CHECK(facade_cache_allocation_exact(constant_cache_after,
                                               constant_cache_success));
    HUNDUN_CHECK(jacobi_allocation_exact(constant_jacobi_after,
                                         constant_jacobi_success));
    for (std::size_t component = 0; component < 3U; ++component) {
      HUNDUN_CHECK(constant_cache_success.operators[component].revision >
                   constant_cache_after.operators[component].revision);
    }
    require_jacobi_operator_coherence(constant_jacobi_success,
                                      constant_cache_success);
    const std::array<std::uint64_t, 2> no_stability_reduction{0U, 0U};
    HUNDUN_CHECK(TimeAccess::stability_reduction_observation() ==
                 no_stability_reduction);
    retained_view_stale = false;
    try {
      static_cast<void>(retained_trial_view(0, 0, 0, 0));
    } catch (const hundun::runtime::Error &) {
      retained_view_stale = true;
    }
    HUNDUN_CHECK(retained_view_stale);

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
    require_authentic_diagnostic(
        mpi, minimum_controller, minimum_state, minimum,
        hundun::diagnostics::DiagnosticStatus::failed,
        hundun::diagnostics::DiagnosticFailureClass::non_finite_state,
        "time-control.non-finite-trial", 0);
    TimeAccess::reset_faults();

    const hundun::config::FlowTimeConfig simultaneous_time{
        hundun::config::TimeMode::fixed, 1, 0.01,
        std::ldexp(0.01, -8), 0.02,
        0.5, 0.25, 1.25, 0.5, 8};
    auto simultaneous_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    simultaneous_state.seed_accepted_layers(initial, initial);
    auto simultaneous_controller =
        hundun::flow::Bdf2RetryController::create(
            simultaneous_time, hundun::config::DensityModel::constant,
            topology, geometry, mpi, simultaneous_state);
    TimeAccess::set_recoverable_failures(9U, 0);
    const auto simultaneous = simultaneous_controller.advance(
        simultaneous_state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(
        simultaneous.disposition() ==
        hundun::flow::TimeAdvanceDisposition::minimum_dt_failure);
    HUNDUN_CHECK(simultaneous.attempt_count() == 9U);
    for (std::size_t index = 0U; index < 9U; ++index) {
      HUNDUN_CHECK(simultaneous.attempt(index).attempted_dt_s ==
                   std::ldexp(0.01, -static_cast<int>(index)));
      HUNDUN_CHECK(
          simultaneous.attempt(index).order ==
          hundun::flow::MomentumTimeOrder::backward_euler);
      HUNDUN_CHECK(
          !simultaneous.attempt(index)
               .all_linear_solves_within_half_limit);
    }
    HUNDUN_CHECK(simultaneous.limited_by_min_dt());
    TimeAccess::reset_faults();

    const auto maximum = std::numeric_limits<std::uint64_t>::max();
    for (const auto accepted_step :
         {std::uint64_t{1} << 53U, (std::uint64_t{1} << 53U) + 1U,
          maximum}) {
      auto terminal_state = hundun::flow::FlowState::create(
          registry, {local, topology.local_face_count()}, fields,
          {0U, 0.0, 0.01, 0.0,
           hundun::flow::MomentumTimeOrder::backward_euler});
      terminal_state.seed_accepted_layers(initial, initial);
      hundun::flow::test::MaterialDensityPisoTestAccess::force_state_metadata(
          terminal_state,
          {accepted_step, 1.0, 0.01, 0.01,
           hundun::flow::MomentumTimeOrder::bdf2});
      hundun::flow::TimeControlState terminal_control;
      terminal_control.accepted_step = accepted_step;
      terminal_control.proposed_next_dt_s = 0.01;
      terminal_control.last_accepted_dt_s = 0.01;
      terminal_control.last_accepted_order =
          hundun::flow::MomentumTimeOrder::bdf2;
      terminal_control.history_ready = true;
      terminal_control.last_stability_metrics_available = true;
      terminal_control.revision = accepted_step;
      terminal_control.state_seal = TimeAccess::seal(
          time, hundun::config::DensityModel::constant, terminal_control);
      auto different_horizon = time;
      different_horizon.steps += 17;
      auto terminal_controller = hundun::flow::Bdf2RetryController::restore(
          different_horizon, hundun::config::DensityModel::constant, topology,
          geometry, mpi, terminal_state, terminal_control);
      const auto terminal_report = terminal_controller.advance(
          terminal_state, facade, 1.0,
          accepted_step == maximum ? 0.0 : -1.0, {}, {});
      HUNDUN_CHECK(
          terminal_report.disposition() ==
          hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
      HUNDUN_CHECK(terminal_report.attempt_count() == 0U);
      HUNDUN_CHECK(terminal_report.lowest_failing_rank() == 0);
      HUNDUN_CHECK(TimeAccess::preflight_category(terminal_report) ==
                   (accepted_step == maximum ? 5U : 1U));
      HUNDUN_CHECK(terminal_controller.state().accepted_step ==
                   accepted_step);

      auto source = terminal_controller.diagnostic_source(terminal_state,
                                                          terminal_report);
      OneRecordSink sink;
      hundun::diagnostics::collect_diagnostics(
          source, mpi,
          {hundun::diagnostics::DiagnosticLevel::bounded_state_sample,
           hundun::diagnostics::DiagnosticScope::collective,
           {mpi.rank(), accepted_step, 1.0,
            "time-control.advance-result"},
           {}, 256U},
          sink);
      HUNDUN_CHECK(sink.calls == 1);
      HUNDUN_CHECK(sink.record.samples.size() == 9U);
      const auto limb = [&](std::size_t index) {
        double value{};
        const auto value_bits = sink.record.samples[index].value.bits;
        std::memcpy(&value, &value_bits, sizeof(value));
        return static_cast<std::uint64_t>(value);
      };
      const auto accepted_reconstructed =
          limb(0U) | (limb(1U) << 32U);
      const auto revision_reconstructed =
          limb(7U) | (limb(8U) << 32U);
      HUNDUN_CHECK(accepted_reconstructed == accepted_step);
      HUNDUN_CHECK(revision_reconstructed == accepted_step);
      if (accepted_step == maximum)
        HUNDUN_CHECK(sink.record.failure.code ==
                     "time-control.preflight.state");
    }
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
  require_zero_trusted_tail_events();
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

  const std::array<std::function<void(hundun::flow::TimeControlState &)>, 13>
      unsealed_state_mutations{{
          [](auto &value) { ++value.schema_version; },
          [](auto &value) { ++value.accepted_step; },
          [](auto &value) {
            value.proposed_next_dt_s =
                std::nextafter(value.proposed_next_dt_s,
                               std::numeric_limits<double>::infinity());
          },
          [](auto &value) {
            value.last_accepted_dt_s =
                std::nextafter(value.last_accepted_dt_s,
                               std::numeric_limits<double>::infinity());
          },
          [](auto &value) {
            value.last_accepted_order =
                value.last_accepted_order ==
                        hundun::flow::MomentumTimeOrder::bdf2
                    ? hundun::flow::MomentumTimeOrder::backward_euler
                    : hundun::flow::MomentumTimeOrder::bdf2;
          },
          [](auto &value) { value.history_ready = !value.history_ready; },
          [](auto &value) {
            value.last_all_linear_solves_within_half_limit =
                !value.last_all_linear_solves_within_half_limit;
          },
          [](auto &value) {
            value.last_convective_rate_per_s =
                std::nextafter(value.last_convective_rate_per_s, 1.0);
          },
          [](auto &value) {
            value.last_diffusive_rate_per_s =
                std::nextafter(value.last_diffusive_rate_per_s, 1.0);
          },
          [](auto &value) {
            value.last_stability_metrics_available =
                !value.last_stability_metrics_available;
          },
          [](auto &value) { ++value.last_retry_count; },
          [](auto &value) { ++value.revision; },
          [](auto &value) { ++value.state_seal; },
      }};
  for (const auto &mutate : unsealed_state_mutations) {
    auto candidate = restored.state();
    mutate(candidate);
    bool rejected = false;
    try {
      static_cast<void>(hundun::flow::Bdf2RetryController::restore(
          time, hundun::config::DensityModel::constant, topology, geometry,
          mpi, state, candidate));
    } catch (const hundun::runtime::Error &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  }

  const std::array<
      std::function<void(hundun::config::FlowTimeConfig &)>, 9>
      sealed_config_mutations{{
          [](auto &value) {
            value.mode = value.mode == hundun::config::TimeMode::adaptive
                             ? hundun::config::TimeMode::fixed
                             : hundun::config::TimeMode::adaptive;
          },
          [](auto &value) {
            value.initial_dt_s =
                std::nextafter(value.initial_dt_s, value.max_dt_s);
          },
          [](auto &value) {
            value.min_dt_s = std::nextafter(value.min_dt_s, 0.0);
          },
          [](auto &value) {
            value.max_dt_s =
                std::nextafter(value.max_dt_s,
                               std::numeric_limits<double>::infinity());
          },
          [](auto &value) {
            value.cfl_target =
                std::nextafter(value.cfl_target, 1.0);
          },
          [](auto &value) {
            value.diffusion_number_target =
                std::nextafter(value.diffusion_number_target, 1.0);
          },
          [](auto &value) {
            value.growth_factor =
                std::nextafter(value.growth_factor, 2.0);
          },
          [](auto &value) {
            value.retry_factor =
                std::nextafter(value.retry_factor, 0.0);
          },
          [](auto &value) { --value.max_retries; },
      }};
  for (const auto &mutate : sealed_config_mutations) {
    auto candidate = time;
    mutate(candidate);
    bool rejected = false;
    try {
      static_cast<void>(hundun::flow::Bdf2RetryController::restore(
          candidate, hundun::config::DensityModel::constant, topology,
          geometry, mpi, state, restored.state()));
    } catch (const hundun::runtime::Error &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  }
  {
    auto different_horizon = time;
    ++different_horizon.steps;
    auto horizon_restore = hundun::flow::Bdf2RetryController::restore(
        different_horizon, hundun::config::DensityModel::constant, topology,
        geometry, mpi, state, restored.state());
    HUNDUN_CHECK(horizon_restore.state().state_seal ==
                 restored.state().state_seal);
  }
  if (mpi.size() > 1) {
    auto rank_distinct = restored.state();
    if (mpi.rank() == 1)
      rank_distinct.last_retry_count =
          rank_distinct.last_retry_count == 0U ? 1U : 0U;
    rank_distinct.state_seal = TimeAccess::seal(
        time, hundun::config::DensityModel::constant, rank_distinct);
    bool rejected = false;
    try {
      static_cast<void>(hundun::flow::Bdf2RetryController::restore(
          time, hundun::config::DensityModel::constant, topology, geometry,
          mpi, state, rank_distinct));
    } catch (const hundun::runtime::Error &error) {
      rejected = true;
      HUNDUN_CHECK(std::string_view(error.what()).find(
                       "time-control.restore.agreement") !=
                   std::string_view::npos);
      HUNDUN_CHECK(std::string_view(error.what()).find("rank 1") !=
                   std::string_view::npos);
    }
    HUNDUN_CHECK(rejected);
  }
  const auto second_step_start =
      hundun::test::capture_adaptive_flow_state(state, restored.state());
  AttemptBoundaryEvidence second_evidence{&restored};
  TimeAccess::set_attempt_observer(&AttemptBoundaryEvidence::observe,
                                   &second_evidence);
  auto second = restored.advance(state, facade, 1.0, 0.0, {}, {});
  HUNDUN_CHECK(second.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  HUNDUN_CHECK(second.attempt(0).order ==
               hundun::flow::MomentumTimeOrder::bdf2);
  HUNDUN_CHECK(second_evidence.count == 1U);
  HUNDUN_CHECK(
      fp64_bits(second.attempt(0).attempted_dt_s /
                second_step_start.metadata.dt_s) ==
      fp64_bits(1.25));
  require_stencil_bitwise_equal(
      second_evidence.stencils[0],
      independent_stencil(hundun::flow::MomentumTimeOrder::bdf2,
                          second.attempt(0).attempted_dt_s,
                          second_step_start.metadata.dt_s));
  require_authentic_diagnostic(
      mpi, restored, state, second,
      hundun::diagnostics::DiagnosticStatus::ok,
      hundun::diagnostics::DiagnosticFailureClass::none, "none", -1);
  TimeAccess::set_attempt_observer(nullptr, nullptr);
  HUNDUN_CHECK(state.metadata().step == 2U);
  const auto fallback_step_start =
      hundun::test::capture_adaptive_flow_state(state, restored.state());
  AttemptBoundaryEvidence fallback_evidence{&restored};
  TimeAccess::set_attempt_observer(&AttemptBoundaryEvidence::observe,
                                   &fallback_evidence);
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
  HUNDUN_CHECK(fallback_evidence.count == fallback.attempt_count());
  for (std::size_t index = 0; index < fallback_evidence.count; ++index) {
    require_stencil_bitwise_equal(
        fallback_evidence.stencils[index],
        independent_stencil(
            fallback.attempt(index).order,
            fallback.attempt(index).attempted_dt_s,
            fallback_step_start.metadata.dt_s));
    auto normalized = fallback_evidence.states[index];
    normalized.trial_generation = fallback_step_start.trial_generation;
    normalized.attempt_identity = fallback_step_start.attempt_identity;
    normalized.diagnostic_identity = fallback_step_start.diagnostic_identity;
    HUNDUN_CHECK(hundun::test::adaptive_flow_state_bitwise_equal(
        normalized, fallback_step_start));
    HUNDUN_CHECK(
        fallback_evidence.states[index].trial_generation ==
        fallback_step_start.trial_generation + index);
    HUNDUN_CHECK(fallback_evidence.states[index].attempt_identity ==
                 fallback_step_start.attempt_identity + index);
    HUNDUN_CHECK(fallback_evidence.states[index].diagnostic_identity ==
                 fallback_step_start.diagnostic_identity + index);
  }
  HUNDUN_CHECK(state.metadata().step == 3U);
  TimeAccess::set_attempt_observer(nullptr, nullptr);
  TimeAccess::reset_faults();
  const auto retry_be_export = restored.state();
  const auto retry_be_before_restore =
      hundun::test::capture_adaptive_flow_state(state, retry_be_export);
  const auto retry_be_history =
      state.snapshot(hundun::flow::FlowLayer::history);
  const auto retry_be_committed =
      state.snapshot(hundun::flow::FlowLayer::committed);
  const auto retry_be_trial = state.snapshot(hundun::flow::FlowLayer::trial);
  const auto retry_be_metadata = state.metadata();
  auto retry_be_restored_state = hundun::flow::FlowState::create(
      registry, {local, topology.local_face_count()}, fields,
      retry_be_metadata);
  retry_be_restored_state.seed_accepted_layers(retry_be_history,
                                               retry_be_committed);
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      retry_be_restored_state.snapshot(hundun::flow::FlowLayer::history),
      retry_be_history));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      retry_be_restored_state.snapshot(hundun::flow::FlowLayer::committed),
      retry_be_committed));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      retry_be_restored_state.snapshot(hundun::flow::FlowLayer::trial),
      retry_be_trial));
  HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
      retry_be_restored_state.metadata(), retry_be_metadata));
  hundun::linear::JacobiPreconditioner restored_mx(execution),
      restored_my(execution), restored_mz(execution),
      restored_pressure(execution);
  auto retry_be_restored_facade =
      hundun::flow::FixedStepConstantDensityFlow::create(
          decomposition, topology, geometry, boundaries, mpi, execution, halo,
          momentum_solver, {&restored_mx, &restored_my, &restored_mz},
          pressure_solver, restored_pressure,
          {{fields.transported_cell_fields.front(),
            hundun::finite_volume::FiniteVolumeQuantity::scalar(0U), 0.0}});
  auto retry_be_restored = hundun::flow::Bdf2RetryController::restore(
      time, hundun::config::DensityModel::constant, topology, geometry, mpi,
      retry_be_restored_state, retry_be_export);
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      retry_be_restored_state.snapshot(hundun::flow::FlowLayer::history),
      retry_be_before_restore.history));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      retry_be_restored_state.snapshot(hundun::flow::FlowLayer::committed),
      retry_be_before_restore.committed));
  HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
      retry_be_restored_state.snapshot(hundun::flow::FlowLayer::trial),
      retry_be_before_restore.trial));
  HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
      retry_be_restored_state.metadata(), retry_be_before_restore.metadata));
  HUNDUN_CHECK(hundun::test::time_control_state_bitwise_equal(
      retry_be_export, retry_be_restored.state()));
  AttemptBoundaryEvidence uninterrupted_evidence{&restored};
  TimeAccess::set_attempt_observer(&AttemptBoundaryEvidence::observe,
                                   &uninterrupted_evidence);
  const auto uninterrupted =
      restored.advance(state, facade, 1.0, 0.0, {}, {});
  HUNDUN_CHECK(uninterrupted.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  HUNDUN_CHECK(uninterrupted_evidence.count == 1U);
  AttemptBoundaryEvidence recovered_evidence{&retry_be_restored};
  TimeAccess::set_attempt_observer(&AttemptBoundaryEvidence::observe,
                                   &recovered_evidence);
  const auto recovered = retry_be_restored.advance(
      retry_be_restored_state, retry_be_restored_facade, 1.0, 0.0, {}, {});
  HUNDUN_CHECK(recovered.disposition() ==
               hundun::flow::TimeAdvanceDisposition::committed);
  HUNDUN_CHECK(recovered.attempt(0).order ==
               hundun::flow::MomentumTimeOrder::bdf2);
  HUNDUN_CHECK(recovered_evidence.count == 1U);
  require_stencil_bitwise_equal(
      recovered_evidence.stencils[0],
      independent_stencil(hundun::flow::MomentumTimeOrder::bdf2,
                          recovered.attempt(0).attempted_dt_s,
                          retry_be_before_restore.metadata.dt_s));
  HUNDUN_CHECK(fp64_bits(uninterrupted.attempt(0).attempted_dt_s) ==
               fp64_bits(recovered.attempt(0).attempted_dt_s));
  HUNDUN_CHECK(uninterrupted.attempt(0).order ==
               recovered.attempt(0).order);
  require_stencil_bitwise_equal(uninterrupted_evidence.stencils[0],
                                recovered_evidence.stencils[0]);
  require_stencil_bitwise_equal(
      uninterrupted_evidence.stencils[0],
      independent_stencil(hundun::flow::MomentumTimeOrder::bdf2,
                          uninterrupted.attempt(0).attempted_dt_s,
                          retry_be_before_restore.metadata.dt_s));
  HUNDUN_CHECK(fp64_bits(uninterrupted.proposed_next_dt_s()) ==
               fp64_bits(recovered.proposed_next_dt_s()));
  HUNDUN_CHECK(fp64_bits(restored.state().proposed_next_dt_s) ==
               fp64_bits(retry_be_restored.state().proposed_next_dt_s));
  HUNDUN_CHECK(restored.state().last_accepted_order ==
               retry_be_restored.state().last_accepted_order);
  HUNDUN_CHECK(restored.state().history_ready ==
               retry_be_restored.state().history_ready);
  HUNDUN_CHECK(restored.state().last_retry_count ==
               retry_be_restored.state().last_retry_count);
  HUNDUN_CHECK(hundun::test::time_control_state_bitwise_equal(
      restored.state(), retry_be_restored.state()));
  HUNDUN_CHECK(hundun::test::accepted_step_metadata_bitwise_equal(
      state.metadata(), retry_be_restored_state.metadata()));
  HUNDUN_CHECK(state.metadata().step == 4U);
  HUNDUN_CHECK(retry_be_restored_state.metadata().step == 4U);
  TimeAccess::set_attempt_observer(nullptr, nullptr);
  TimeAccess::set_attempt_observer(&AttemptBoundaryEvidence::observe,
                                   &uninterrupted_evidence);
  TimeAccess::exercise_trusted_tail_attempt_observer(
      state, uninterrupted_evidence.stencils[0]);
  HUNDUN_CHECK(
      TimeAccess::trusted_tail_observation().callbacks_or_sinks == 1U);
  TimeAccess::set_attempt_observer(nullptr, nullptr);

  {
    const hundun::config::FlowTimeConfig ratio_time{
        hundun::config::TimeMode::fixed, 5, 0.01, 0.00001, 0.02,
        0.5, 0.25, 1.25, 0.5, 8};
    auto ratio_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    ratio_state.seed_accepted_layers(initial, initial);
    auto ratio_controller = hundun::flow::Bdf2RetryController::create(
        ratio_time, hundun::config::DensityModel::constant, topology,
        geometry, mpi, ratio_state);
    TimeAccess::reset_faults();
    TimeAccess::set_recoverable_failures(1U, 0);
    const auto reduced_start = ratio_controller.advance(
        ratio_state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(reduced_start.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::committed);
    HUNDUN_CHECK(fp64_bits(ratio_state.metadata().dt_s) ==
                 fp64_bits(0.005));
    TimeAccess::reset_faults();

    AttemptBoundaryEvidence ratio_two{&ratio_controller};
    TimeAccess::set_attempt_observer(&AttemptBoundaryEvidence::observe,
                                     &ratio_two);
    const auto doubled =
        ratio_controller.advance(ratio_state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(doubled.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::committed);
    HUNDUN_CHECK(ratio_two.count == 1U);
    HUNDUN_CHECK(fp64_bits(doubled.attempt(0).attempted_dt_s / 0.005) ==
                 fp64_bits(2.0));
    require_stencil_bitwise_equal(
        ratio_two.stencils[0],
        independent_stencil(hundun::flow::MomentumTimeOrder::bdf2, 0.01,
                            0.005));

    AttemptBoundaryEvidence ratio_one{&ratio_controller};
    TimeAccess::set_attempt_observer(&AttemptBoundaryEvidence::observe,
                                     &ratio_one);
    const auto equal =
        ratio_controller.advance(ratio_state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(equal.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::committed);
    HUNDUN_CHECK(ratio_one.count == 1U);
    HUNDUN_CHECK(fp64_bits(equal.attempt(0).attempted_dt_s /
                          ratio_state.metadata().previous_dt_s) ==
                 fp64_bits(1.0));
    require_stencil_bitwise_equal(
        ratio_one.stencils[0],
        independent_stencil(hundun::flow::MomentumTimeOrder::bdf2, 0.01,
                            0.01));

    AttemptBoundaryEvidence ratio_half{&ratio_controller};
    TimeAccess::set_attempt_observer(&AttemptBoundaryEvidence::observe,
                                     &ratio_half);
    TimeAccess::set_recoverable_failures(1U, 0);
    const auto halved =
        ratio_controller.advance(ratio_state, facade, 1.0, 0.0, {}, {});
    HUNDUN_CHECK(halved.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::committed);
    HUNDUN_CHECK(ratio_half.count == 2U);
    HUNDUN_CHECK(fp64_bits(halved.attempt(1).attempted_dt_s /
                          equal.accepted_dt_s()) ==
                 fp64_bits(0.5));
    require_stencil_bitwise_equal(
        ratio_half.stencils[1],
        independent_stencil(hundun::flow::MomentumTimeOrder::bdf2, 0.005,
                            0.01));
    TimeAccess::set_attempt_observer(nullptr, nullptr);
    TimeAccess::reset_faults();
  }

  auto other_controller = hundun::flow::Bdf2RetryController::restore(
      time, hundun::config::DensityModel::constant, topology, geometry, mpi,
      state, retry_be_restored.state());
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
  const auto controller_before_rejection = retry_be_restored.state();
  const auto rejected_after_out_of_band =
      retry_be_restored.advance(state, facade, 1.0, 0.0, {}, {});
  HUNDUN_CHECK(rejected_after_out_of_band.disposition() ==
               hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
  HUNDUN_CHECK(rejected_after_out_of_band.attempt_count() == 0U);
  HUNDUN_CHECK(retry_be_restored.state().state_seal ==
               controller_before_rejection.state_seal);
  {
    auto bound_to_other = hundun::flow::Bdf2RetryController::create(
        time, hundun::config::DensityModel::constant, topology, geometry, mpi,
        other_state);
    auto moved_to = std::move(other_state);
    TimeAccess::reset_faults();
    bool rejected = false;
    try {
      static_cast<void>(
          bound_to_other.advance(moved_to, facade, 1.0, 0.0, {}, {}));
    } catch (const hundun::runtime::Error &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(TimeAccess::raw_operation_count() == 0U);

    TimeAccess::reset_faults();
    rejected = false;
    try {
      static_cast<void>(
          controller.advance(other_state, facade, 1.0, 0.0, {}, {}));
    } catch (const hundun::runtime::Error &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(TimeAccess::raw_operation_count() == 0U);
  }

  {
    TimeAccess::reset_faults();
    TimeAccess::set_active(controller, true);
    const auto state_before_overlap =
        state.snapshot(hundun::flow::FlowLayer::committed);
    const auto control_before_overlap = controller.state();
    bool rejected = false;
    try {
      static_cast<void>(
          controller.advance(state, facade, 1.0, 0.0, {}, {}));
    } catch (const hundun::runtime::Error &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(TimeAccess::raw_operation_count() == 0U);
    HUNDUN_CHECK(hundun::test::flow_layer_values_bitwise_equal(
        state_before_overlap,
        state.snapshot(hundun::flow::FlowLayer::committed)));
    HUNDUN_CHECK(hundun::test::time_control_state_bitwise_equal(
        control_before_overlap, controller.state()));
    TimeAccess::set_active(controller, false);
  }

  {
    auto live_controller = std::move(controller);
    TimeAccess::reset_faults();
    bool rejected = false;
    try {
      static_cast<void>(
          controller.advance(state, facade, 1.0, 0.0, {}, {}));
    } catch (const hundun::runtime::Error &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(TimeAccess::raw_operation_count() == 0U);

    auto live_facade = std::move(facade);
    TimeAccess::reset_faults();
    rejected = false;
    try {
      static_cast<void>(
          live_controller.advance(state, facade, 1.0, 0.0, {}, {}));
    } catch (const hundun::runtime::Error &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(TimeAccess::raw_operation_count() == 0U);
  }
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
  config.scalars.push_back({"beta", 0.0});
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
  const auto rho_alpha =
      registry.declare_field(conservative_cell("rho_alpha", "kg/m3"));
  const auto rho_beta =
      registry.declare_field(conservative_cell("rho_beta", "kg/m3"));
  fields.transported_cell_fields = {rho_h, rho_alpha, rho_beta};
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
        std::vector<double>(topology.owned_cell_count(), 300000.0),
        std::vector<double>(topology.owned_cell_count(), 0.5),
        std::vector<double>(topology.owned_cell_count(), 0.25)};
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
  material_spec.scalar_densities = {rho_alpha, rho_beta};
  material_spec.scalar_diffusivities_kg_per_m_s = {0.2, 0.4};
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
    using MaterialAccess =
        hundun::flow::test::MaterialDensityPisoTestAccess;
    auto cold_state = make_state();
    hundun::linear::JacobiPreconditioner cold_mx(execution),
        cold_my(execution), cold_mz(execution), cold_pressure(execution);
    auto cold_flow = hundun::flow::FixedStepMaterialDensityFlow::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        momentum_solver, {&cold_mx, &cold_my, &cold_mz}, pressure_solver,
        cold_pressure, registry, fields, material_spec);
    const auto cold_cache_before =
        MaterialAccess::facade_cache_snapshot(cold_flow);
    const auto cold_jacobi_before =
        capture_jacobi_cache({&cold_mx, &cold_my, &cold_mz});
    require_cold_jacobi_unpublished(cold_jacobi_before);
    require_jacobi_mutation_sensitive(cold_jacobi_before);
    JacobiColdFaultReset cold_fault_reset;
    hundun::linear::test::PreconditionerTestAccess::
        arm_fail_next_cold_update_before_publication();
    const auto cold_report = cold_flow.attempt(
        cold_state, 0.0,
        independent_stencil(
            hundun::flow::MomentumTimeOrder::backward_euler, 0.01, 0.0),
        {}, {});
    hundun::linear::test::PreconditionerTestAccess::
        reset_cold_update_fault();
    HUNDUN_CHECK(
        cold_report.flow().disposition ==
        hundun::flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(cold_report.flow().reason ==
                 hundun::flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(cold_report.flow().lowest_failing_rank == 0);
    const auto cold_cache = MaterialAccess::facade_cache_snapshot(cold_flow);
    const auto cold_jacobi =
        capture_jacobi_cache({&cold_mx, &cold_my, &cold_mz});
    HUNDUN_CHECK(cold_cache.operators[0].revision >
                 cold_cache_before.operators[0].revision);
    require_cold_jacobi_unpublished(cold_jacobi);
    require_jacobi_mutation_sensitive(cold_jacobi);

    const auto clean_report = cold_flow.attempt(
        cold_state, 0.0,
        independent_stencil(
            hundun::flow::MomentumTimeOrder::backward_euler, 0.01, 0.0),
        {}, {});
    HUNDUN_CHECK(clean_report.flow().disposition ==
                 hundun::flow::StepAttemptDisposition::committed);
    const auto published_cache =
        MaterialAccess::facade_cache_snapshot(cold_flow);
    const auto published_jacobi =
        capture_jacobi_cache({&cold_mx, &cold_my, &cold_mz});
    for (std::size_t component = 0; component < 3U; ++component) {
      HUNDUN_CHECK(published_cache.operators[component].revision >
                   cold_cache.operators[component].revision);
    }
    require_jacobi_operator_coherence(published_jacobi, published_cache);
    require_jacobi_mutation_sensitive(published_jacobi);
  }
  {
    using IdealAccess =
        hundun::flow::test::IdealGasClosureTestAccess;
    auto cold_state = make_state();
    hundun::linear::JacobiPreconditioner cold_mx(execution),
        cold_my(execution), cold_mz(execution), cold_pressure(execution);
    auto cold_closure = hundun::flow::IdealGasClosure::create(
        topology, geometry, boundaries, mpi, registry, fields, cold_state,
        {rho_h, 1000.0, 287.0, 86100.0});
    auto cold_flow = hundun::flow::FixedStepIdealGasFlow::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        momentum_solver, {&cold_mx, &cold_my, &cold_mz}, pressure_solver,
        cold_pressure, registry, fields, material_spec,
        std::move(cold_closure));
    const auto cold_cache_before =
        IdealAccess::facade_cache_snapshot(cold_flow);
    const auto cold_owned_before =
        IdealAccess::delegated_material_cache_snapshot(cold_flow);
    HUNDUN_CHECK(
        facade_cache_delegation_exact(cold_cache_before, cold_owned_before));
    const auto cold_jacobi_before =
        capture_jacobi_cache({&cold_mx, &cold_my, &cold_mz});
    require_cold_jacobi_unpublished(cold_jacobi_before);
    require_jacobi_mutation_sensitive(cold_jacobi_before);
    JacobiColdFaultReset cold_fault_reset;
    hundun::linear::test::PreconditionerTestAccess::
        arm_fail_next_cold_update_before_publication();
    const auto cold_report = cold_flow.attempt(
        cold_state, 0.0,
        independent_stencil(
            hundun::flow::MomentumTimeOrder::backward_euler, 0.01, 0.0),
        {}, {});
    hundun::linear::test::PreconditionerTestAccess::
        reset_cold_update_fault();
    const auto &cold_flow_report = cold_report.flow().flow();
    HUNDUN_CHECK(
        cold_flow_report.disposition ==
        hundun::flow::StepAttemptDisposition::non_retryable_failure);
    HUNDUN_CHECK(cold_flow_report.reason ==
                 hundun::flow::StepFailureReason::invalid_input);
    HUNDUN_CHECK(cold_flow_report.lowest_failing_rank == 0);
    const auto cold_cache = IdealAccess::facade_cache_snapshot(cold_flow);
    const auto cold_owned =
        IdealAccess::delegated_material_cache_snapshot(cold_flow);
    HUNDUN_CHECK(facade_cache_delegation_exact(cold_cache, cold_owned));
    HUNDUN_CHECK(cold_cache.operators[0].revision >
                 cold_cache_before.operators[0].revision);
    const auto cold_jacobi =
        capture_jacobi_cache({&cold_mx, &cold_my, &cold_mz});
    require_cold_jacobi_unpublished(cold_jacobi);
    require_jacobi_mutation_sensitive(cold_jacobi);

    const auto clean_report = cold_flow.attempt(
        cold_state, 0.0,
        independent_stencil(
            hundun::flow::MomentumTimeOrder::backward_euler, 0.01, 0.0),
        {}, {});
    HUNDUN_CHECK(clean_report.flow().flow().disposition ==
                 hundun::flow::StepAttemptDisposition::committed);
    const auto published_cache =
        IdealAccess::facade_cache_snapshot(cold_flow);
    const auto published_owned =
        IdealAccess::delegated_material_cache_snapshot(cold_flow);
    HUNDUN_CHECK(
        facade_cache_delegation_exact(published_cache, published_owned));
    const auto published_jacobi =
        capture_jacobi_cache({&cold_mx, &cold_my, &cold_mz});
    for (std::size_t component = 0; component < 3U; ++component) {
      HUNDUN_CHECK(published_cache.operators[component].revision >
                   cold_cache.operators[component].revision);
    }
    require_jacobi_operator_coherence(published_jacobi, published_cache);
    require_jacobi_mutation_sensitive(published_jacobi);
  }
  if (mpi.size() > 1) {
    using AuthorityMutation =
        hundun::flow::test::MaterialTransportAuthorityMutation;
    for (const auto mutation :
         {AuthorityMutation::omitted, AuthorityMutation::reordered,
          AuthorityMutation::same_maximum_different_sequence}) {
      hundun::linear::JacobiPreconditioner mismatch_mx(execution),
          mismatch_my(execution), mismatch_mz(execution),
          mismatch_pressure(execution);
      auto mismatch_flow =
          hundun::flow::FixedStepMaterialDensityFlow::create(
              decomposition, topology, geometry, boundaries, mpi, execution,
              halo, momentum_solver,
              {&mismatch_mx, &mismatch_my, &mismatch_mz}, pressure_solver,
              mismatch_pressure, registry, fields, material_spec);
      if (mpi.rank() == 1)
        hundun::flow::test::MaterialDensityPisoTestAccess::
            mutate_transport_authority(mismatch_flow, mutation);
      auto mismatch_state = make_state();
      auto mismatch_controller =
          hundun::flow::Bdf2RetryController::create(
              time, hundun::config::DensityModel::material, topology,
              geometry, mpi, mismatch_state);
      const auto mismatch_report = mismatch_controller.advance(
          mismatch_state, mismatch_flow, 0.0, {}, {});
      HUNDUN_CHECK(
          mismatch_report.disposition() ==
          hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
      HUNDUN_CHECK(mismatch_report.attempt_count() == 0U);
      HUNDUN_CHECK(mismatch_report.lowest_failing_rank() == 1);
      HUNDUN_CHECK(TimeAccess::preflight_category(mismatch_report) == 6U);
      require_authentic_diagnostic(
          mpi, mismatch_controller, mismatch_state, mismatch_report,
          hundun::diagnostics::DiagnosticStatus::failed,
          hundun::diagnostics::DiagnosticFailureClass::invalid_input,
          "time-control.preflight.transport-authority", 1);
    }
  }
  {
    TimeAccess::reset_faults();
    hundun::flow::FlowLayerValues variable_initial;
    const auto cell_count = topology.owned_cell_count();
    variable_initial.density.resize(cell_count);
    for (std::size_t cell_id = 0; cell_id < cell_count; ++cell_id)
      variable_initial.density[cell_id] =
          1.0 + 0.05 * static_cast<double>((cell_id % 7U) + 1U);
    variable_initial.velocity.assign(cell_count * 3U, 0.0);
    variable_initial.mechanical_pressure.assign(cell_count, 0.0);
    variable_initial.face_velocity.assign(topology.local_face_count() * 3U,
                                          0.0);
    variable_initial.face_mass_flux.resize(topology.local_face_count());
    for (std::size_t face_id = 0;
         face_id < variable_initial.face_mass_flux.size(); ++face_id)
      variable_initial.face_mass_flux[face_id] =
          face_id % 2U == 0U ? 1.0e-12 : -1.0e-12;
    variable_initial.transported_cell_fields = {
        std::vector<double>(cell_count, 300000.0),
        std::vector<double>(cell_count, 0.5),
        std::vector<double>(cell_count, 0.25)};
    auto variable_state = hundun::flow::FlowState::create(
        registry, {local, topology.local_face_count()}, fields,
        {0U, 0.0, 0.01, 0.0,
         hundun::flow::MomentumTimeOrder::backward_euler});
    variable_state.seed_accepted_layers(variable_initial, variable_initial);
    std::vector<double> flux_sum(cell_count);
    std::vector<double> geometry_sum(cell_count);
    for (hundun::mesh::LocalFaceId face_id = 0;
         face_id < topology.local_face_count(); ++face_id) {
      const auto area = geometry.face_area_m2(face_id);
      const auto displacement = geometry.face_displacement_m(face_id);
      const auto area_vector = geometry.face_area_vector_m2(
          face_id, hundun::mesh::FaceSide::owner);
      const auto factor =
          area * area /
          std::abs(area_vector.x * displacement.x +
                   area_vector.y * displacement.y +
                   area_vector.z * displacement.z);
      const auto add = [&](std::size_t cell_id) {
        flux_sum[cell_id] +=
            std::abs(variable_initial.face_mass_flux[face_id]);
        geometry_sum[cell_id] += factor;
      };
      const auto owner = topology.owner(face_id);
      if (owner < cell_count)
        add(owner);
      const auto neighbour = topology.neighbour(face_id);
      if (neighbour && *neighbour < cell_count)
        add(*neighbour);
    }
    double expected[2]{};
    for (std::size_t cell_id = 0; cell_id < cell_count; ++cell_id) {
      const auto denominator =
          2.0 * variable_initial.density[cell_id] *
          geometry.cell_volume_m3(cell_id);
      expected[0] =
          std::max(expected[0], flux_sum[cell_id] / denominator);
      expected[1] =
          std::max(expected[1], 0.4 * geometry_sum[cell_id] / denominator);
    }
    const auto reductions_before = mpi.fp64_reduction_counters();
    HUNDUN_CHECK(MPI_Allreduce(MPI_IN_PLACE, expected, 2, MPI_DOUBLE, MPI_MAX,
                               mpi.comm()) == MPI_SUCCESS);
    auto variable_controller = hundun::flow::Bdf2RetryController::create(
        time, hundun::config::DensityModel::material, topology, geometry, mpi,
        variable_state);
    const auto variable_report = variable_controller.advance(
        variable_state, material_flow, 0.0, {}, {});
    const auto reductions_after = mpi.fp64_reduction_counters();
    HUNDUN_CHECK(variable_report.stability_metrics_available());
    HUNDUN_CHECK(variable_report.convective_rate_per_s() == expected[0]);
    HUNDUN_CHECK(variable_report.diffusive_rate_per_s() == expected[1]);
    HUNDUN_CHECK(reductions_after.collective_calls >
                 reductions_before.collective_calls);
    HUNDUN_CHECK(reductions_after.reduced_scalars >=
                 reductions_before.reduced_scalars + 2U);
    const auto stability_reduction =
        TimeAccess::stability_reduction_observation();
    HUNDUN_CHECK(stability_reduction[0] == 1U);
    HUNDUN_CHECK(stability_reduction[1] == 2U);
  }
  {
    auto state = make_state();
    const hundun::config::FlowTimeConfig terminal_time{
        hundun::config::TimeMode::fixed, 1, 0.01, 0.00001, 0.02,
        0.5, 0.25, 1.25, 0.5, 8};
    auto terminal_controller = hundun::flow::Bdf2RetryController::create(
        terminal_time, hundun::config::DensityModel::material, topology,
        geometry, mpi, state);
    const auto material_cache_before =
        hundun::flow::test::MaterialDensityPisoTestAccess::
            facade_cache_snapshot(material_flow);
    const auto material_jacobi_before =
        capture_jacobi_cache({&mx, &my, &mz});
    require_facade_cache_inventory(material_cache_before, false);
    require_facade_cache_mutation_sensitive(material_cache_before);
    require_jacobi_operator_coherence(material_jacobi_before,
                                      material_cache_before);
    require_jacobi_mutation_sensitive(material_jacobi_before);
    const auto terminal_reductions_before = mpi.fp64_reduction_counters();
    const auto terminal_before = hundun::test::capture_adaptive_flow_state(
        state, terminal_controller.state(), {}, 0U, 0U, 0U,
        {terminal_reductions_before.collective_calls,
         terminal_reductions_before.reduced_scalars,
         terminal_reductions_before.logical_payload_bytes, 0U});
    TimeAccess::set_recoverable_failures(9U, 0);
    const auto terminal = terminal_controller.advance(
        state, material_flow, 0.0, {}, {});
    HUNDUN_CHECK(terminal.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::retry_limit_reached);
    HUNDUN_CHECK(terminal.attempt_count() == 9U);
    HUNDUN_CHECK(terminal.reason() ==
                 hundun::flow::StepFailureReason::final_continuity_residual);
    require_authentic_diagnostic(
        mpi, terminal_controller, state, terminal,
        hundun::diagnostics::DiagnosticStatus::failed,
        hundun::diagnostics::DiagnosticFailureClass::non_convergence,
        "time-control.final-continuity-residual", 0);
    const auto terminal_reductions_after = mpi.fp64_reduction_counters();
    const auto material_cache_after =
        hundun::flow::test::MaterialDensityPisoTestAccess::
            facade_cache_snapshot(material_flow);
    const auto material_jacobi_after =
        capture_jacobi_cache({&mx, &my, &mz});
    HUNDUN_CHECK(facade_cache_allocation_exact(material_cache_before,
                                               material_cache_after));
    HUNDUN_CHECK(
        jacobi_allocation_exact(material_jacobi_before, material_jacobi_after));
    for (std::size_t component = 0; component < 3U; ++component) {
      HUNDUN_CHECK(material_cache_after.operators[component].revision >
                   material_cache_before.operators[component].revision);
    }
    require_jacobi_operator_coherence(material_jacobi_after,
                                      material_cache_after);
    const auto terminal_after = hundun::test::capture_adaptive_flow_state(
        state, terminal_controller.state(), {}, 0U, 0U, 0U,
        {terminal_reductions_after.collective_calls,
         terminal_reductions_after.reduced_scalars,
         terminal_reductions_after.logical_payload_bytes, 0U});
    HUNDUN_CHECK(hundun::test::adaptive_flow_state_failed_attempt_preserved(
        terminal_before, terminal_after, terminal.attempt_count(),
        {1440U, 81432U, 651456U, 0U}));
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
    require_zero_trusted_tail_events();
    HUNDUN_CHECK(std::holds_alternative<
                 hundun::flow::MaterialDensityStepAttemptReport>(
        report.final_attempt()));
    const auto &material_final =
        std::get<hundun::flow::MaterialDensityStepAttemptReport>(
            report.final_attempt());
    HUNDUN_CHECK(material_final.flux_provenance() ==
                 hundun::flow::MaterialFluxProvenance::final_corrected);
    HUNDUN_CHECK(material_final.material_field_count() == 3U);
    const auto material_cache_success =
        hundun::flow::test::MaterialDensityPisoTestAccess::
            facade_cache_snapshot(material_flow);
    const auto material_jacobi_success =
        capture_jacobi_cache({&mx, &my, &mz});
    HUNDUN_CHECK(facade_cache_allocation_exact(material_cache_after,
                                               material_cache_success));
    HUNDUN_CHECK(
        jacobi_allocation_exact(material_jacobi_after, material_jacobi_success));
    for (std::size_t component = 0; component < 3U; ++component) {
      HUNDUN_CHECK(material_cache_success.operators[component].revision >
                   material_cache_after.operators[component].revision);
    }
    require_jacobi_operator_coherence(material_jacobi_success,
                                      material_cache_success);
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
    using MaterialAccess =
        hundun::flow::test::MaterialDensityPisoTestAccess;
    hundun::linear::JacobiPreconditioner collective_mx(execution),
        collective_my(execution), collective_mz(execution),
        collective_pressure(execution);
    auto collective_material_flow =
        hundun::flow::FixedStepMaterialDensityFlow::create(
            decomposition, topology, geometry, boundaries, mpi, execution,
            halo, momentum_solver,
            {&collective_mx, &collective_my, &collective_mz},
            pressure_solver, collective_pressure, registry, fields,
            material_spec);
    auto state = make_state();
    auto controller = hundun::flow::Bdf2RetryController::create(
        time, hundun::config::DensityModel::material, topology, geometry, mpi,
        state);
    MaterialAccess::reset_terminal_fault();
    MaterialAccess::set_terminal_fault(
        hundun::flow::test::MaterialTerminalPointForTest::
            final_pressure_entry,
        hundun::flow::test::MaterialTerminalModeForTest::returned_reliable);
    const auto report = controller.advance(
        state, collective_material_flow, 0.0, {}, {});
    MaterialAccess::reset_terminal_fault();
    HUNDUN_CHECK(
        report.disposition() ==
        hundun::flow::TimeAdvanceDisposition::non_retryable_failure);
    HUNDUN_CHECK(report.reason() ==
                 hundun::flow::StepFailureReason::collective_operation);
    require_authentic_diagnostic(
        mpi, controller, state, report,
        hundun::diagnostics::DiagnosticStatus::failed,
        hundun::diagnostics::DiagnosticFailureClass::collective_operation,
        "time-control.collective-operation", 0);
  }

  {
    auto terminal_state = make_state();
    hundun::linear::JacobiPreconditioner terminal_mx(execution),
        terminal_my(execution), terminal_mz(execution),
        terminal_pressure(execution);
    auto terminal_closure = hundun::flow::IdealGasClosure::create(
        topology, geometry, boundaries, mpi, registry, fields, terminal_state,
        {rho_h, 1000.0, 287.0, 86100.0});
    auto terminal_flow = hundun::flow::FixedStepIdealGasFlow::create(
        decomposition, topology, geometry, boundaries, mpi, execution, halo,
        momentum_solver, {&terminal_mx, &terminal_my, &terminal_mz},
        pressure_solver, terminal_pressure, registry, fields, material_spec,
        std::move(terminal_closure));
    const hundun::config::FlowTimeConfig terminal_time{
        hundun::config::TimeMode::fixed, 1, 0.01, 0.00001, 0.02,
        0.5, 0.25, 1.25, 0.5, 8};
    auto prewarm_controller = hundun::flow::Bdf2RetryController::create(
        terminal_time, hundun::config::DensityModel::ideal_gas, topology,
        geometry, mpi, terminal_state);
    const auto prewarm = prewarm_controller.advance(
        terminal_state, terminal_flow, 0.0, {}, {});
    HUNDUN_CHECK(prewarm.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::committed);
    const auto prewarm_export = prewarm_controller.state();
    auto terminal_controller = hundun::flow::Bdf2RetryController::restore(
        terminal_time, hundun::config::DensityModel::ideal_gas, topology,
        geometry, mpi, terminal_state, prewarm_export);
    const auto ideal_cache_before =
        hundun::flow::test::IdealGasClosureTestAccess::
            facade_cache_snapshot(terminal_flow);
    const auto ideal_owned_before =
        hundun::flow::test::IdealGasClosureTestAccess::
            delegated_material_cache_snapshot(terminal_flow);
    const auto ideal_jacobi_before =
        capture_jacobi_cache({&terminal_mx, &terminal_my, &terminal_mz});
    require_facade_cache_inventory(ideal_cache_before, true);
    require_facade_cache_mutation_sensitive(ideal_cache_before);
    HUNDUN_CHECK(
        facade_cache_delegation_exact(ideal_cache_before, ideal_owned_before));
    require_jacobi_operator_coherence(ideal_jacobi_before, ideal_cache_before);
    require_jacobi_mutation_sensitive(ideal_jacobi_before);
    const auto terminal_reductions_before = mpi.fp64_reduction_counters();
    const auto terminal_before = hundun::test::capture_adaptive_flow_state(
        terminal_state, terminal_controller.state(),
        terminal_flow.closure_state(), 0U, 0U, 0U,
        {terminal_reductions_before.collective_calls,
         terminal_reductions_before.reduced_scalars,
         terminal_reductions_before.logical_payload_bytes, 0U});
    TimeAccess::set_recoverable_failure_reason(
        hundun::flow::StepFailureReason::final_conservation_defect);
    TimeAccess::set_recoverable_failures(9U, 0);
    const auto terminal = terminal_controller.advance(
        terminal_state, terminal_flow, 0.0, {}, {});
    HUNDUN_CHECK(terminal.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::retry_limit_reached);
    HUNDUN_CHECK(terminal.attempt_count() == 9U);
    HUNDUN_CHECK(terminal.reason() ==
                 hundun::flow::StepFailureReason::final_conservation_defect);
    require_authentic_diagnostic(
        mpi, terminal_controller, terminal_state, terminal,
        hundun::diagnostics::DiagnosticStatus::failed,
        hundun::diagnostics::DiagnosticFailureClass::conservation,
        "time-control.final-conservation-defect", 0);
    const auto terminal_reductions_after = mpi.fp64_reduction_counters();
    const auto ideal_cache_after =
        hundun::flow::test::IdealGasClosureTestAccess::
            facade_cache_snapshot(terminal_flow);
    const auto ideal_owned_after =
        hundun::flow::test::IdealGasClosureTestAccess::
            delegated_material_cache_snapshot(terminal_flow);
    const auto ideal_jacobi_after =
        capture_jacobi_cache({&terminal_mx, &terminal_my, &terminal_mz});
    HUNDUN_CHECK(
        facade_cache_allocation_exact(ideal_cache_before, ideal_cache_after));
    HUNDUN_CHECK(
        facade_cache_delegation_exact(ideal_cache_after, ideal_owned_after));
    HUNDUN_CHECK(
        jacobi_allocation_exact(ideal_jacobi_before, ideal_jacobi_after));
    for (std::size_t component = 0; component < 3U; ++component) {
      HUNDUN_CHECK(ideal_cache_after.operators[component].revision >
                   ideal_cache_before.operators[component].revision);
    }
    require_jacobi_operator_coherence(ideal_jacobi_after, ideal_cache_after);
    const auto terminal_after = hundun::test::capture_adaptive_flow_state(
        terminal_state, terminal_controller.state(),
        terminal_flow.closure_state(), 0U, 0U, 0U,
        {terminal_reductions_after.collective_calls,
         terminal_reductions_after.reduced_scalars,
         terminal_reductions_after.logical_payload_bytes, 0U});
    HUNDUN_CHECK(hundun::test::adaptive_flow_state_failed_attempt_preserved(
        terminal_before, terminal_after, terminal.attempt_count(),
        {1602U, 82350U, 658800U, 0U}));
    TimeAccess::reset_faults();

    auto transport_controller = hundun::flow::Bdf2RetryController::restore(
        terminal_time, hundun::config::DensityModel::ideal_gas, topology,
        geometry, mpi, terminal_state, prewarm_export);
    TimeAccess::set_recoverable_failure_reason(
        hundun::flow::StepFailureReason::transport_failure);
    TimeAccess::set_recoverable_failures(9U, failure_rank);
    const auto transport_terminal = transport_controller.advance(
        terminal_state, terminal_flow, 0.0, {}, {});
    HUNDUN_CHECK(
        transport_terminal.disposition() ==
        hundun::flow::TimeAdvanceDisposition::retry_limit_reached);
    HUNDUN_CHECK(transport_terminal.reason() ==
                 hundun::flow::StepFailureReason::transport_failure);
    require_authentic_diagnostic(
        mpi, transport_controller, terminal_state, transport_terminal,
        hundun::diagnostics::DiagnosticStatus::failed,
        hundun::diagnostics::DiagnosticFailureClass::numerical_breakdown,
        "time-control.transport-failure", failure_rank);
    TimeAccess::reset_faults();
    const auto ideal_success = transport_controller.advance(
        terminal_state, terminal_flow, 0.0, {}, {});
    HUNDUN_CHECK(ideal_success.disposition() ==
                 hundun::flow::TimeAdvanceDisposition::committed);
    const auto ideal_cache_success =
        hundun::flow::test::IdealGasClosureTestAccess::
            facade_cache_snapshot(terminal_flow);
    const auto ideal_owned_success =
        hundun::flow::test::IdealGasClosureTestAccess::
            delegated_material_cache_snapshot(terminal_flow);
    const auto ideal_jacobi_success =
        capture_jacobi_cache({&terminal_mx, &terminal_my, &terminal_mz});
    HUNDUN_CHECK(
        facade_cache_allocation_exact(ideal_cache_after, ideal_cache_success));
    HUNDUN_CHECK(facade_cache_delegation_exact(ideal_cache_success,
                                               ideal_owned_success));
    HUNDUN_CHECK(
        jacobi_allocation_exact(ideal_jacobi_after, ideal_jacobi_success));
    for (std::size_t component = 0; component < 3U; ++component) {
      HUNDUN_CHECK(ideal_cache_success.operators[component].revision >
                   ideal_cache_after.operators[component].revision);
    }
    require_jacobi_operator_coherence(ideal_jacobi_success,
                                      ideal_cache_success);
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
    require_zero_trusted_tail_events();
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

    auto terminal_open_controller =
        hundun::flow::Bdf2RetryController::restore(
            open_time, hundun::config::DensityModel::ideal_gas, open_topology,
            open_geometry, mpi, open_state, open_controller.state());
    TimeAccess::set_recoverable_failure_reason(
        hundun::flow::StepFailureReason::boundary_backflow);
    TimeAccess::set_recoverable_failures(9U, 0);
    const auto terminal_open = terminal_open_controller.advance(
        open_state, open_flow, 0.0, {}, {});
    HUNDUN_CHECK(
        terminal_open.disposition() ==
        hundun::flow::TimeAdvanceDisposition::minimum_dt_failure);
    HUNDUN_CHECK(terminal_open.reason() ==
                 hundun::flow::StepFailureReason::boundary_backflow);
    const int outlet_coordinates[3]{grid(mpi.size()).x - 1, 0, 0};
    int outlet_rank = -1;
    hundun::runtime::check_mpi_result(
        MPI_Cart_rank(open_decomposition.comm(), outlet_coordinates,
                      &outlet_rank),
        "MPI_Cart_rank(Task22 outlet authority)");
    require_authentic_diagnostic(
        mpi, terminal_open_controller, open_state, terminal_open,
        hundun::diagnostics::DiagnosticStatus::failed,
        hundun::diagnostics::DiagnosticFailureClass::boundary,
        "time-control.boundary-backflow", outlet_rank);
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

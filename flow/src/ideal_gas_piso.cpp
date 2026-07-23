// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow/ideal_gas_piso.hpp"

#include "density_closure_detail.hpp"
#include "fixed_step_flow_detail.hpp"
#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#include "ideal_gas_closure_test_access.hpp"
#endif

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <utility>

namespace hundun::flow {
namespace {

constexpr std::uint64_t kIdealStepSeed = 0x696465616c737465ULL;
constexpr runtime::PhaseId kStatePhase = 1800U;
constexpr runtime::ActorId kStateActor = 1800U;

void mix(std::uint64_t &hash, std::uint64_t value) noexcept {
  hash ^= value + 0x9e3779b97f4a7c15ULL + (hash << 6U) + (hash >> 2U);
}

void require_index(std::size_t index, std::size_t size, const char *message) {
  if (index >= size)
    throw runtime::Error(message);
}

} // namespace

detail::DensityClosureHooks
detail::DensityClosureAdapter::bind(IdealGasClosure &closure,
                                    runtime::FieldId enthalpy_density,
                                    double enthalpy_rate_J_per_kg_s) {
  DensityClosureHooks hooks;
  hooks.object = &closure;
  hooks.enthalpy_density = enthalpy_density;
  hooks.enthalpy_rate_J_per_kg_s = enthalpy_rate_J_per_kg_s;
  hooks.begin = &DensityClosureAdapter::begin;
  hooks.evaluate = &DensityClosureAdapter::evaluate;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  hooks.prepare_attempt = &DensityClosureAdapter::prepare_attempt;
  hooks.begin_allocation_observation =
      &DensityClosureAdapter::begin_allocation_observation;
  hooks.end_allocation_observation =
      &DensityClosureAdapter::end_allocation_observation;
  hooks.outer_failure = &DensityClosureAdapter::outer_failure;
  hooks.after_halo = &DensityClosureAdapter::after_halo;
  hooks.before_post_assessment =
      &DensityClosureAdapter::before_post_assessment;
#endif
  hooks.before_outlet = &DensityClosureAdapter::before_outlet;
  hooks.before_prepare = &DensityClosureAdapter::before_prepare;
  hooks.prepare = &DensityClosureAdapter::prepare;
  hooks.publish = &DensityClosureAdapter::publish;
  hooks.rollback = &DensityClosureAdapter::rollback;
  return hooks;
}

void detail::DensityClosureAdapter::begin(void *object, const FlowState &state,
                                          std::uint64_t identity) {
  static_cast<IdealGasClosure *>(object)->begin_attempt(state, identity);
}

detail::DensityClosureEvaluation
detail::DensityClosureAdapter::evaluate(void *object, FlowState &state,
                                        DensityClosureStage stage) {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  auto &closure = *static_cast<IdealGasClosure *>(object);
  test::IdealGasClosureTestAccess::begin_allocation_observation(closure);
  struct ObservationGuard final {
    IdealGasClosure &closure;
    ~ObservationGuard() {
      test::IdealGasClosureTestAccess::end_allocation_observation(closure);
    }
  } guard{closure};
#endif
  const auto concrete_stage = stage == DensityClosureStage::predictor
                                  ? IdealGasClosureStage::predictor
                              : stage == DensityClosureStage::provisional
                                  ? IdealGasClosureStage::provisional
                                  : IdealGasClosureStage::final;
  const auto &report =
      static_cast<IdealGasClosure *>(object)->evaluate(state, concrete_stage);
  return {report.disposition() == IdealGasClosureDisposition::closed,
          report.disposition() ==
              IdealGasClosureDisposition::recoverable_failure,
          report.lowest_failing_rank()};
}

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
void detail::DensityClosureAdapter::prepare_attempt(void *object) {
  auto &closure = *static_cast<IdealGasClosure *>(object);
  test::IdealGasClosureTestAccess::consume_attempt_preparation_fault(closure);
}
void detail::DensityClosureAdapter::begin_allocation_observation(
    void *object) noexcept {
  test::IdealGasClosureTestAccess::begin_allocation_observation(
      *static_cast<IdealGasClosure *>(object));
}
void detail::DensityClosureAdapter::end_allocation_observation(
    void *object) noexcept {
  test::IdealGasClosureTestAccess::end_allocation_observation(
      *static_cast<IdealGasClosure *>(object));
}
int detail::DensityClosureAdapter::outer_failure(
    void *object, DensityClosureOuterPoint point) {
  return static_cast<IdealGasClosure *>(object)
      ->outer_failure_for_test(static_cast<std::uint8_t>(point));
}
void detail::DensityClosureAdapter::after_halo(
    void *object, DensityClosureStage stage, runtime::FieldId density,
    runtime::FieldId enthalpy_density) {
  static_cast<IdealGasClosure *>(object)->record_halo_for_test(
      static_cast<std::uint8_t>(
          stage == DensityClosureStage::predictor
              ? IdealGasClosureStage::predictor
          : stage == DensityClosureStage::provisional
              ? IdealGasClosureStage::provisional
              : IdealGasClosureStage::final),
      density, enthalpy_density);
}
void detail::DensityClosureAdapter::before_post_assessment(void *object,
                                                           FlowState &state) {
  static_cast<IdealGasClosure *>(object)->before_post_assessment_for_test(
      state);
}
#endif

void detail::DensityClosureAdapter::before_outlet(void *object,
                                                  FlowState &state) {
  static_cast<IdealGasClosure *>(object)->before_outlet(state);
}
int detail::DensityClosureAdapter::before_prepare(
    void *object, FlowState &state, AcceptedStepMetadata accepted) {
  return static_cast<IdealGasClosure *>(object)->before_prepare(state,
                                                                accepted);
}

int detail::DensityClosureAdapter::prepare(void *object) {
  return static_cast<IdealGasClosure *>(object)->prepare_commit();
}
void detail::DensityClosureAdapter::publish(void *object) noexcept {
  static_cast<IdealGasClosure *>(object)->publish_commit();
}
void detail::DensityClosureAdapter::rollback(void *object) noexcept {
  static_cast<IdealGasClosure *>(object)->rollback();
}
double detail::DensityClosureAdapter::gas_constant_J_per_kg_K(
    const IdealGasClosure &closure) noexcept {
  return closure.gas_constant_J_per_kg_K();
}
double detail::DensityClosureAdapter::cp_J_per_kg_K(
    const IdealGasClosure &closure) noexcept {
  return closure.cp_J_per_kg_K();
}

IdealGasStepAttemptReport::IdealGasStepAttemptReport()
    : flow_(detail::DensityClosureBridge::make_report()) {}

const MaterialDensityStepAttemptReport &
IdealGasStepAttemptReport::flow() const noexcept {
  return flow_;
}
bool IdealGasStepAttemptReport::closure_report_available() const noexcept {
  return closure_report_.has_value();
}
const IdealGasClosureReport &IdealGasStepAttemptReport::closure_report() const {
  if (!closure_report_)
    throw runtime::Error("ideal-gas closure report is unavailable");
  return *closure_report_;
}
std::uint64_t IdealGasStepAttemptReport::attempt_identity() const noexcept {
  return attempt_identity_;
}
std::uint64_t IdealGasStepAttemptReport::compute_seal() const noexcept {
  std::uint64_t hash = kIdealStepSeed;
  mix(hash, attempt_identity_);
  mix(hash, detail::DensityClosureBridge::report_seal(flow_));
  mix(hash, closure_report_.has_value() ? 1U : 0U);
  if (closure_report_)
    mix(hash, closure_report_->compute_seal());
  return hash;
}
bool IdealGasStepAttemptReport::semantic_valid() const noexcept {
  if (!detail::DensityClosureBridge::report_authenticated(flow_) ||
      attempt_identity_ == 0U ||
      attempt_identity_ != flow_.attempt_identity() ||
      (closure_report_ &&
       (!closure_report_->authenticated() ||
        closure_report_->attempt_identity() != attempt_identity_)))
    return false;
  const bool closure_failure =
      flow_.flow().reason == StepFailureReason::density_closure_failure;
  if (closure_failure)
    return closure_report_.has_value() &&
           closure_report_->disposition() != IdealGasClosureDisposition::closed;
  if (flow_.flow().disposition == StepAttemptDisposition::committed)
    return closure_report_.has_value() &&
           closure_report_->disposition() ==
               IdealGasClosureDisposition::closed &&
           closure_report_->stage() == IdealGasClosureStage::final &&
           closure_report_->final_metrics_available();
  return !closure_report_ ||
         closure_report_->disposition() == IdealGasClosureDisposition::closed;
}
void IdealGasStepAttemptReport::seal() noexcept {
  seal_ = detail::DensityClosureBridge::report_authenticated(flow_) &&
                  (!closure_report_ || closure_report_->authenticated()) &&
                  attempt_identity_ == flow_.attempt_identity()
              ? compute_seal()
              : 0U;
}
bool IdealGasStepAttemptReport::authenticated() const noexcept {
  return seal_ != 0U && seal_ == compute_seal() &&
         detail::DensityClosureBridge::report_authenticated(flow_) &&
         (!closure_report_ || closure_report_->authenticated()) &&
         attempt_identity_ == flow_.attempt_identity() && semantic_valid();
}

struct FixedStepIdealGasFlow::Impl final {
  Impl(const runtime::StructuredDecomposition &decomposition,
       const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
       const boundary::BoundaryRegistry &boundaries,
       const runtime::MpiContext &mpi, const runtime::FieldRegistry &registry,
       FlowFieldIds supplied_fields, FixedStepMaterialDensityFlow material,
       IdealGasClosure supplied_closure)
      : decomposition(&decomposition), topology(&topology), geometry(&geometry),
        boundaries(&boundaries), mpi(&mpi), registry(&registry),
        fields(std::move(supplied_fields)), material(std::move(material)),
        closure(std::move(supplied_closure)) {
    const auto box = topology.owned_global_box();
    const auto global = topology.global_extent();
    owned_cell_fingerprint =
        "cell.f64.owned." + std::to_string(box.begin.x) + "." +
        std::to_string(box.begin.y) + "." + std::to_string(box.begin.z) + "." +
        std::to_string(box.end.x) + "." + std::to_string(box.end.y) + "." +
        std::to_string(box.end.z);
    global_cell_fingerprint = "cell.f64.global." + std::to_string(global.x) +
                              "." + std::to_string(global.y) + "." +
                              std::to_string(global.z);
  }

  const runtime::StructuredDecomposition *decomposition;
  const mesh::MeshTopology *topology;
  const mesh::MeshGeometry *geometry;
  const boundary::BoundaryRegistry *boundaries;
  const runtime::MpiContext *mpi;
  const runtime::FieldRegistry *registry;
  FlowFieldIds fields;
  FixedStepMaterialDensityFlow material;
  IdealGasClosure closure;
  double enthalpy_rate_J_per_kg_s{};
  std::uint64_t source_generation{1U};
  const FlowState *last_state{};
  std::uint64_t last_state_identity{};
  std::uint64_t last_report_seal{};
  std::string owned_cell_fingerprint;
  std::string global_cell_fingerprint;
};

struct IdealGasClosureDiagnosticSource::Impl final {
  const FixedStepIdealGasFlow::Impl *flow{};
  const FlowState *state{};
  IdealGasStepAttemptReport report;
  std::uint64_t source_generation{};
  std::uint64_t state_identity{};
  std::uint64_t report_seal{};
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  mutable detail::DensityClosureDiagnosticTestState diagnostic_test_state;
#endif
};

FixedStepIdealGasFlow::FixedStepIdealGasFlow(
    std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
FixedStepIdealGasFlow::~FixedStepIdealGasFlow() noexcept = default;
FixedStepIdealGasFlow::FixedStepIdealGasFlow(
    FixedStepIdealGasFlow &&other) noexcept
    : impl_(std::move(other.impl_)) {
  if (impl_) {
    if (impl_->source_generation == std::numeric_limits<std::uint64_t>::max())
      impl_->source_generation = 0U;
    else
      ++impl_->source_generation;
  }
}

FixedStepIdealGasFlow FixedStepIdealGasFlow::create(
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const runtime::MpiContext &mpi,
    execution::ExecutionContext &execution_context,
    runtime::HaloExchange &cell_halo,
    const linear::LinearSolver &momentum_solver,
    std::array<linear::Preconditioner *, 3> momentum_preconditioners,
    const linear::LinearSolver &pressure_solver,
    linear::Preconditioner &pressure_preconditioner,
    const runtime::FieldRegistry &registry, FlowFieldIds fields,
    MaterialDensityTransportSpec specification, IdealGasClosure &&closure) {
  const auto match = runtime::collective_status(
      mpi,
      closure.matches(topology, geometry, boundaries, mpi, registry, fields),
      "ideal-gas closure collaborators do not match flow");
  if (!match.ok)
    throw runtime::Error("ideal-gas closure collaborators do not match flow");
  std::unique_ptr<Impl> prepared;
  bool construction_ok = true;
  try {
    auto material = detail::DensityClosureBridge::create_open_capable(
        decomposition, topology, geometry, boundaries, mpi, execution_context,
        cell_halo, momentum_solver, momentum_preconditioners, pressure_solver,
        pressure_preconditioner, registry, fields, specification);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (test::IdealGasClosureTestAccess::consume_facade_create_fault(
            closure, mpi.rank())) {
      throw std::bad_alloc();
    }
#endif
    prepared = std::make_unique<Impl>(
        decomposition, topology, geometry, boundaries, mpi, registry,
        std::move(fields), std::move(material), std::move(closure));
  } catch (const runtime::MpiOperationError &) {
    throw;
  } catch (...) {
    construction_ok = false;
  }
  const auto construction = runtime::collective_status(
      mpi, construction_ok, "ideal-gas flow facade construction failed");
  if (!construction.ok)
    throw runtime::Error("ideal-gas flow facade construction failed on rank " +
                         std::to_string(construction.failing_rank));
  if (!prepared)
    throw runtime::Error("ideal-gas flow facade construction produced no value");
  return FixedStepIdealGasFlow(std::move(prepared));
}

IdealGasStepAttemptReport FixedStepIdealGasFlow::attempt(
    FlowState &state, double mu, const MomentumTimeStencil &stencil,
    const linear::SolveControl &momentum_control,
    const linear::SolveControl &pressure_control) const {
  if (!impl_)
    throw runtime::Error("ideal-gas flow object has been moved from");
  impl_->last_state = nullptr;
  impl_->last_state_identity = 0U;
  impl_->last_report_seal = 0U;
  if (impl_->source_generation == 0U ||
      impl_->source_generation == std::numeric_limits<std::uint64_t>::max())
    throw runtime::Error("ideal-gas diagnostic source generation would wrap");
  ++impl_->source_generation;
  IdealGasStepAttemptReport result;
  const auto closure_hooks = detail::DensityClosureAdapter::bind(
      impl_->closure, impl_->fields.transported_cell_fields.front(),
      impl_->enthalpy_rate_J_per_kg_s);
  result.flow_ = detail::DensityClosureBridge::attempt(
      impl_->material, state, mu, stencil, momentum_control, pressure_control,
      closure_hooks);
  result.attempt_identity_ = result.flow_.attempt_identity();
  try {
    const auto &latest = impl_->closure.latest_report();
    if (latest.attempt_identity() == result.attempt_identity_)
      result.closure_report_ = latest;
    else
      result.closure_report_.reset();
  } catch (const runtime::Error &) {
    result.closure_report_.reset();
  }
  result.seal();
  impl_->last_state = &state;
  impl_->last_state_identity = state.diagnostic_mutation_identity();
  impl_->last_report_seal = result.seal_;
  return result;
}

IdealGasClosureState FixedStepIdealGasFlow::closure_state() const {
  if (!impl_)
    throw runtime::Error("ideal-gas flow object has been moved from");
  return impl_->closure.state();
}

MaterialDensityFlowDiagnosticSource
FixedStepIdealGasFlow::flow_diagnostic_source(
    const FlowState &state, const IdealGasStepAttemptReport &report) const {
  if (!impl_ || !report.authenticated() ||
      report.flow().flow().reason == StepFailureReason::density_closure_failure)
    throw runtime::Error("ideal-gas flow diagnostic source is unavailable");
  return impl_->material.diagnostic_source(state, report.flow());
}

IdealGasClosureDiagnosticSource
FixedStepIdealGasFlow::closure_diagnostic_source(
    const FlowState &state, const IdealGasStepAttemptReport &report) const {
  if (!impl_)
    throw runtime::Error("ideal-gas flow object has been moved from");
  const auto state_identity = diagnostic_state_identity(state);
  if (!report.authenticated() || impl_->last_state != &state ||
      !diagnostic_state_matches(state, state_identity) ||
      state_identity == 0U ||
      impl_->last_state_identity != state_identity ||
      impl_->last_report_seal != report.seal_)
    throw runtime::Error("ideal-gas closure diagnostic source is stale");
  auto source = std::make_unique<IdealGasClosureDiagnosticSource::Impl>();
  source->flow = impl_.get();
  source->state = &state;
  source->report = report;
  source->source_generation = impl_->source_generation;
  source->state_identity = state_identity;
  source->report_seal = report.seal_;
  return IdealGasClosureDiagnosticSource(std::move(source));
}

IdealGasClosureDiagnosticSource::IdealGasClosureDiagnosticSource(
    std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
IdealGasClosureDiagnosticSource::~IdealGasClosureDiagnosticSource() noexcept =
    default;
IdealGasClosureDiagnosticSource::IdealGasClosureDiagnosticSource(
    IdealGasClosureDiagnosticSource &&other) noexcept
    : impl_(std::move(other.impl_)) {}

void IdealGasClosureDiagnosticSource::validate() const {
  if (!impl_ || !impl_->flow || !impl_->state ||
      impl_->source_generation != impl_->flow->source_generation ||
      !FixedStepIdealGasFlow::diagnostic_state_matches(*impl_->state,
                                                       impl_->state_identity) ||
      impl_->report_seal != impl_->report.seal_ ||
      !impl_->report.authenticated())
    throw runtime::Error("ideal-gas closure diagnostic source is stale");
}

std::uint64_t FixedStepIdealGasFlow::diagnostic_state_identity(
    const FlowState &state) noexcept {
  return state.diagnostic_mutation_identity();
}

bool FixedStepIdealGasFlow::diagnostic_state_matches(
    const FlowState &state, std::uint64_t identity) noexcept {
  return identity != 0U &&
         state.diagnostic_mutation_identity() == identity &&
         !state.attempt_active();
}

std::size_t IdealGasClosureDiagnosticSource::fingerprint_field_count() const {
  validate();
  return 3U;
}
std::string_view
IdealGasClosureDiagnosticSource::fingerprint_field_id(std::size_t field) const {
  validate();
  constexpr std::string_view ids[]{"p0", "rho", "rho_h"};
  require_index(field, 3U, "ideal-gas fingerprint field is out of bounds");
  return ids[field];
}
std::string_view IdealGasClosureDiagnosticSource::fingerprint_field_unit(
    std::size_t field) const {
  validate();
  constexpr std::string_view units[]{"Pa", "kg/m3", "J/m3"};
  require_index(field, 3U, "ideal-gas fingerprint field is out of bounds");
  return units[field];
}
IdealGasDiagnosticEntity
IdealGasClosureDiagnosticSource::fingerprint_field_entity(
    std::size_t field) const {
  validate();
  require_index(field, 3U, "ideal-gas fingerprint field is out of bounds");
  return field == 0U ? IdealGasDiagnosticEntity::global_scalar
                     : IdealGasDiagnosticEntity::cell;
}
std::size_t IdealGasClosureDiagnosticSource::fingerprint_field_item_count(
    std::size_t field) const {
  validate();
  require_index(field, 3U, "ideal-gas fingerprint field is out of bounds");
  return field == 0U ? (impl_->flow->mpi->rank() == 0 ? 1U : 0U)
                     : impl_->flow->topology->owned_cell_count();
}
std::uint64_t IdealGasClosureDiagnosticSource::fingerprint_field_global_id(
    std::size_t field, std::size_t item) const {
  const auto count = fingerprint_field_item_count(field);
  require_index(item, count, "ideal-gas fingerprint item is out of bounds");
  return field == 0U ? 0U : impl_->flow->topology->global_cell_id(item);
}

double
IdealGasClosureDiagnosticSource::committed_cell_value(runtime::FieldId field,
                                                      std::size_t item) const {
  const auto values =
      detail::DensityClosureDiagnosticAccess::acquire_committed(*this);
  const auto extent = values.density.interior_extent();
  const std::size_t count = static_cast<std::size_t>(extent.x) *
                            static_cast<std::size_t>(extent.y) *
                            static_cast<std::size_t>(extent.z);
  require_index(item, count, "ideal-gas diagnostic cell is out of bounds");
  const auto plane =
      static_cast<std::size_t>(extent.x) * static_cast<std::size_t>(extent.y);
  const int k = static_cast<int>(item / plane);
  const auto within = item % plane;
  const int j = static_cast<int>(within / static_cast<std::size_t>(extent.x));
  const int i = static_cast<int>(within % static_cast<std::size_t>(extent.x));
  if (field == impl_->flow->fields.density)
    return values.density(i, j, k, 0);
  if (field == impl_->flow->fields.transported_cell_fields.front())
    return values.enthalpy_density(i, j, k, 0);
  throw runtime::Error("ideal-gas diagnostic field is invalid");
}

double IdealGasClosureDiagnosticSource::fingerprint_field_value(
    std::size_t field, std::size_t item) const {
  validate();
  if (field == 0U) {
    require_index(item, fingerprint_field_item_count(field),
                  "ideal-gas p0 item is out of bounds");
    return impl_->flow->closure.state().thermodynamic_pressure_pa;
  }
  if (field == 1U)
    return committed_cell_value(impl_->flow->fields.density, item);
  if (field == 2U)
    return committed_cell_value(
        impl_->flow->fields.transported_cell_fields.front(), item);
  throw runtime::Error("ideal-gas fingerprint field is out of bounds");
}

std::size_t IdealGasClosureDiagnosticSource::sample_field_count() const {
  validate();
  return 5U;
}
std::string_view
IdealGasClosureDiagnosticSource::sample_field_id(std::size_t field) const {
  validate();
  constexpr std::string_view ids[]{"enthalpy", "p0", "rho", "rho_h",
                                   "temperature"};
  require_index(field, 5U, "ideal-gas sample field is out of bounds");
  return ids[field];
}
std::string_view
IdealGasClosureDiagnosticSource::sample_field_unit(std::size_t field) const {
  validate();
  constexpr std::string_view units[]{"J/kg", "Pa", "kg/m3", "J/m3", "K"};
  require_index(field, 5U, "ideal-gas sample field is out of bounds");
  return units[field];
}
IdealGasDiagnosticEntity
IdealGasClosureDiagnosticSource::sample_field_entity(std::size_t field) const {
  validate();
  require_index(field, 5U, "ideal-gas sample field is out of bounds");
  return field == 1U ? IdealGasDiagnosticEntity::global_scalar
                     : IdealGasDiagnosticEntity::cell;
}
std::size_t IdealGasClosureDiagnosticSource::sample_field_item_count(
    std::size_t field) const {
  validate();
  require_index(field, 5U, "ideal-gas sample field is out of bounds");
  return field == 1U ? (impl_->flow->mpi->rank() == 0 ? 1U : 0U)
                     : impl_->flow->topology->owned_cell_count();
}
std::uint64_t IdealGasClosureDiagnosticSource::sample_field_global_id(
    std::size_t field, std::size_t item) const {
  const auto count = sample_field_item_count(field);
  require_index(item, count, "ideal-gas sample item is out of bounds");
  return field == 1U ? 0U : impl_->flow->topology->global_cell_id(item);
}
double
IdealGasClosureDiagnosticSource::sample_field_value(std::size_t field,
                                                    std::size_t item) const {
  validate();
  if (field == 1U) {
    require_index(item, sample_field_item_count(field),
                  "ideal-gas p0 sample item is out of bounds");
    return impl_->flow->closure.state().thermodynamic_pressure_pa;
  }
  const double rho = committed_cell_value(impl_->flow->fields.density, item);
  const double rho_h = committed_cell_value(
      impl_->flow->fields.transported_cell_fields.front(), item);
  if (field == 0U)
    return rho_h / rho;
  if (field == 2U)
    return rho;
  if (field == 3U)
    return rho_h;
  if (field == 4U)
    return (rho_h / rho) / impl_->flow->closure.cp_J_per_kg_K();
  throw runtime::Error("ideal-gas sample field is out of bounds");
}

std::string_view
IdealGasClosureDiagnosticSource::owned_cell_layout_fingerprint() const {
  validate();
  return impl_->flow->owned_cell_fingerprint;
}
std::string_view
IdealGasClosureDiagnosticSource::global_cell_layout_fingerprint() const {
  validate();
  return impl_->flow->global_cell_fingerprint;
}
int IdealGasClosureDiagnosticSource::relative_rank() const {
  validate();
  return impl_->flow->mpi->rank();
}
std::uint64_t IdealGasClosureDiagnosticSource::committed_step() const {
  validate();
  return impl_->state->metadata().step;
}
double IdealGasClosureDiagnosticSource::committed_time_s() const {
  validate();
  return impl_->state->metadata().time_s;
}
std::size_t IdealGasClosureDiagnosticSource::owned_cell_count() const {
  validate();
  return impl_->flow->topology->owned_cell_count();
}
double IdealGasClosureDiagnosticSource::cell_volume_m3(std::size_t cell) const {
  validate();
  require_index(cell, owned_cell_count(),
                "ideal-gas diagnostic volume is out of bounds");
  return impl_->flow->geometry->cell_volume_m3(cell);
}
IdealGasClosureState IdealGasClosureDiagnosticSource::closure_state() const {
  validate();
  return impl_->flow->closure.state();
}
const IdealGasStepAttemptReport &
IdealGasClosureDiagnosticSource::report() const {
  validate();
  return impl_->report;
}

detail::DensityClosureReadSession
detail::DensityClosureDiagnosticAccess::acquire_committed(
    const IdealGasClosureDiagnosticSource &source) {
  source.validate();
  const auto &state = *source.impl_->state;
  const auto &access = FlowStateSolverAccess::access(state);
  const auto &storage = state.layer(FlowLayer::committed);
  return {storage.acquire_read<double>(access, kStatePhase, kStateActor,
                                       source.impl_->flow->fields.density),
          storage.acquire_read<double>(
              access, kStatePhase, kStateActor,
              source.impl_->flow->fields.transported_cell_fields.front()),
          source.impl_->flow->topology->owned_global_box(),
          source.impl_->flow->topology->global_extent()};
}

double detail::DensityClosureDiagnosticAccess::gas_constant_J_per_kg_K(
    const IdealGasClosureDiagnosticSource &source) {
  source.validate();
  return DensityClosureAdapter::gas_constant_J_per_kg_K(
      source.impl_->flow->closure);
}

double detail::DensityClosureDiagnosticAccess::cp_J_per_kg_K(
    const IdealGasClosureDiagnosticSource &source) {
  source.validate();
  return DensityClosureAdapter::cp_J_per_kg_K(source.impl_->flow->closure);
}

const runtime::MpiContext &detail::DensityClosureDiagnosticAccess::mpi(
    const IdealGasClosureDiagnosticSource &source) {
  source.validate();
  return *source.impl_->flow->mpi;
}

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
detail::DensityClosureDiagnosticTestState &
detail::DensityClosureDiagnosticAccess::test_state(
    const IdealGasClosureDiagnosticSource &source) {
  source.validate();
  return source.impl_->diagnostic_test_state;
}
#endif

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
void test::IdealGasClosureTestAccess::set_uniform_enthalpy_rate(
    FixedStepIdealGasFlow &flow, double rate) {
  if (!flow.impl_ || !std::isfinite(rate))
    throw runtime::Error("ideal-gas test enthalpy rate is invalid");
  flow.impl_->enthalpy_rate_J_per_kg_s = rate;
}
bool test::IdealGasClosureTestAccess::report_authenticated(
    const IdealGasClosureReport &report) noexcept {
  return report.authenticated();
}
bool test::IdealGasClosureTestAccess::report_authenticated(
    const IdealGasStepAttemptReport &report) noexcept {
  return report.authenticated();
}
bool test::IdealGasClosureTestAccess::post_eos_evidence_authenticated(
    const MaterialDensityStepAttemptReport &report) noexcept {
  return ::hundun::flow::detail::DensityClosureBridge::
      post_eos_evidence_authenticated(report);
}
void test::IdealGasClosureTestAccess::exhaust_source_generation(
    FixedStepIdealGasFlow &flow) {
  if (!flow.impl_)
    throw runtime::Error("ideal-gas test flow has been moved from");
  flow.impl_->source_generation = std::numeric_limits<std::uint64_t>::max();
}
void test::IdealGasClosureTestAccess::force_finalization_identity_wrap(
    FixedStepIdealGasFlow &flow) {
  if (!flow.impl_)
    throw runtime::Error("ideal-gas test flow has been moved from");
  ::hundun::flow::detail::DensityClosureBridge::
      force_finalization_identity_wrap(flow.impl_->material);
}
std::vector<test::IdealGasHaloTraceEntry>
test::IdealGasClosureTestAccess::halo_trace(
    const FixedStepIdealGasFlow &flow) {
  if (!flow.impl_)
    throw runtime::Error("ideal-gas test flow has been moved from");
  std::vector<IdealGasHaloTraceEntry> result;
  for (const auto &entry : flow.impl_->closure.halo_trace_for_test())
    result.push_back(
        {static_cast<IdealGasClosureStage>(entry[0]),
         static_cast<runtime::FieldId>(entry[1]),
         static_cast<runtime::FieldId>(entry[2])});
  return result;
}
std::uint64_t test::IdealGasClosureTestAccess::source_generation(
    const IdealGasClosureDiagnosticSource &source) {
  source.validate();
  return source.impl_->source_generation;
}
void test::IdealGasClosureTestAccess::set_post_store_corruption(
    FixedStepIdealGasFlow &flow, int rank, bool enthalpy_density) {
  if (!flow.impl_)
    throw runtime::Error("ideal-gas test flow has been moved from");
  flow.impl_->closure.set_post_store_corruption_for_test(rank,
                                                        enthalpy_density);
}
void test::IdealGasClosureTestAccess::set_candidate_precedence_fault(
    FixedStepIdealGasFlow &flow, int rank) {
  if (!flow.impl_)
    throw runtime::Error("ideal-gas test flow has been moved from");
  flow.impl_->closure.set_candidate_precedence_fault_for_test(rank);
}
void test::IdealGasClosureTestAccess::set_stage_failure(
    FixedStepIdealGasFlow &flow, IdealGasClosureStage stage,
    IdealGasClosureFailureReason reason, int rank) {
  if (!flow.impl_)
    throw runtime::Error("ideal-gas test flow has been moved from");
  flow.impl_->closure.set_stage_failure_for_test(stage, reason, rank);
}
void test::IdealGasClosureTestAccess::set_outer_failure(
    FixedStepIdealGasFlow &flow, IdealGasOuterFailurePoint point, int rank) {
  if (!flow.impl_)
    throw runtime::Error("ideal-gas test flow has been moved from");
  flow.impl_->closure.set_outer_failure_for_test(
      static_cast<std::uint8_t>(point), rank);
}
void test::IdealGasClosureTestAccess::set_prepare_fault(
    FixedStepIdealGasFlow &flow, IdealGasPrepareFault fault, int rank) {
  if (!flow.impl_)
    throw runtime::Error("ideal-gas test flow has been moved from");
  flow.impl_->closure.set_prepare_fault_for_test(
      fault == IdealGasPrepareFault::state_prepare, rank);
}
void test::IdealGasClosureTestAccess::set_post_store_mpi_fault(
    FixedStepIdealGasFlow &flow, int rank) {
  if (!flow.impl_)
    throw runtime::Error("ideal-gas test flow has been moved from");
  flow.impl_->closure.set_post_store_mpi_fault_for_test(rank);
}
void test::IdealGasClosureTestAccess::set_post_assessment_fault(
    FixedStepIdealGasFlow &flow, IdealGasPostAssessmentFault fault, int rank) {
  if (!flow.impl_)
    throw runtime::Error("ideal-gas test flow has been moved from");
  if (fault == IdealGasPostAssessmentFault::non_finite_state ||
      fault == IdealGasPostAssessmentFault::non_positive_density) {
    flow.impl_->closure.set_post_assessment_fault_for_test(
        static_cast<std::uint8_t>(fault), rank);
    return;
  }
  ::hundun::flow::detail::DensityClosureBridge::set_post_assessment_fault(
      flow.impl_->material,
      fault == IdealGasPostAssessmentFault::final_transport_residual ? 0U : 1U,
      rank);
}
void test::IdealGasClosureTestAccess::set_attempt_layout_fault(
    FixedStepIdealGasFlow &flow, int rank) {
  if (!flow.impl_)
    throw runtime::Error("ideal-gas test flow has been moved from");
  flow.impl_->closure.set_attempt_layout_fault_for_test(rank);
}
void test::IdealGasClosureTestAccess::set_outlet_backflow_fault(
    FixedStepIdealGasFlow &flow) {
  if (!flow.impl_)
    throw runtime::Error("ideal-gas test flow has been moved from");
  flow.impl_->closure.set_outlet_backflow_fault_for_test();
}
void test::IdealGasClosureTestAccess::set_attempt_preparation_fault(
    FixedStepIdealGasFlow &flow, IdealGasAttemptPreparationFault fault,
    int rank) {
  if (!flow.impl_)
    throw runtime::Error("ideal-gas test flow has been moved from");
  const auto code = static_cast<std::uint8_t>(fault);
  if (fault >= IdealGasAttemptPreparationFault::post_density_views &&
      fault <= IdealGasAttemptPreparationFault::post_report_finalization) {
    ::hundun::flow::detail::DensityClosureBridge::set_post_assessment_fault(
        flow.impl_->material,
        static_cast<std::uint8_t>(
            code - static_cast<std::uint8_t>(
                       IdealGasAttemptPreparationFault::post_density_views) +
            2U),
        rank);
    return;
  }
  set_attempt_preparation_fault(flow.impl_->closure, fault, rank);
}
void test::IdealGasClosureTestAccess::set_controlled_allocation(
    FixedStepIdealGasFlow &flow, int rank) {
  if (!flow.impl_)
    throw runtime::Error("ideal-gas test flow has been moved from");
  set_controlled_allocation(flow.impl_->closure, rank);
}
bool test::IdealGasClosureTestAccess::allocation_observation_active(
    const FixedStepIdealGasFlow &flow) noexcept {
  return flow.impl_ && allocation_observation_active(flow.impl_->closure);
}
bool test::IdealGasClosureTestAccess::post_evidence_mutation_rejected(
    const IdealGasStepAttemptReport &report,
    IdealGasPostEvidenceMutation mutation) {
  return ::hundun::flow::detail::DensityClosureBridge::
      post_evidence_mutation_rejected(report.flow_,
                                      static_cast<std::uint8_t>(mutation));
}
#endif

} // namespace hundun::flow

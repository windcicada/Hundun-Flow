// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow/adaptive_time_control.hpp"

#include "adaptive_time_control_detail.hpp"
#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#include "adaptive_time_control_test_access.hpp"
#include "constant_density_piso_test_access.hpp"
#include "ideal_gas_closure_test_access.hpp"
#include "material_density_transport_test_access.hpp"
#endif

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstring>
#include <limits>
#include <new>
#include <sstream>
#include <type_traits>
#include <utility>
#include <vector>

namespace hundun::flow {
namespace {

constexpr std::uint64_t kFnvOffset = 14695981039346656037ULL;
constexpr std::uint64_t kFnvPrime = 1099511628211ULL;
constexpr char kSealDomain[] = "hundun-time-control-state-seal-v1";

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
std::atomic<test::TimeControlPostReturnMutation> post_return_mutation{
    test::TimeControlPostReturnMutation::none};
std::atomic<test::TimeControlPreflightFault> preflight_fault{
    test::TimeControlPreflightFault::none};
std::atomic<int> time_control_fault_rank{-1};
std::atomic<std::uint32_t> scheduled_recoverable_failures{};
std::atomic<StepFailureReason> scheduled_recoverable_reason{
    StepFailureReason::none};
std::atomic<std::uint64_t> post_return_iteration_value{
    std::numeric_limits<std::uint64_t>::max()};
std::atomic<std::size_t> time_control_raw_count{};
std::atomic<std::size_t> time_control_raw_fault_ordinal{};
std::atomic<int> time_control_raw_fault_rank{-1};
std::atomic<bool> trusted_tail_active{};
std::atomic<std::uint64_t> trusted_tail_allocation_attempts{};
std::atomic<std::uint64_t> trusted_tail_controller_collectives{};
std::atomic<std::uint64_t> trusted_tail_field_state_traversals{};
std::atomic<std::uint64_t> trusted_tail_callbacks_or_sinks{};
std::atomic<std::uint64_t> stability_reduction_calls{};
std::atomic<std::uint64_t> stability_reduced_scalars{};
#endif

std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

bool positive_zero(double value) noexcept { return bits(value) == 0U; }

void byte(std::uint64_t &hash, std::uint8_t value) noexcept {
  hash = (hash ^ value) * kFnvPrime;
}
template <class UInt>
void little_endian(std::uint64_t &hash, UInt value) noexcept {
  for (std::size_t i = 0; i < sizeof(UInt); ++i)
    byte(hash, static_cast<std::uint8_t>(value >> (8U * i)));
}
void fp64(std::uint64_t &hash, double value) noexcept {
  little_endian(hash, bits(value));
}

bool valid_time_mode(config::TimeMode mode) noexcept {
  return mode == config::TimeMode::fixed ||
         mode == config::TimeMode::adaptive;
}

bool valid_density_model(config::DensityModel model) noexcept {
  return model == config::DensityModel::constant ||
         model == config::DensityModel::material ||
         model == config::DensityModel::ideal_gas;
}

bool valid_config(const config::FlowTimeConfig &config) noexcept {
  return valid_time_mode(config.mode) && config.steps >= 0 &&
         std::isfinite(config.initial_dt_s) &&
         std::isfinite(config.min_dt_s) && std::isfinite(config.max_dt_s) &&
         config.initial_dt_s > 0.0 && config.min_dt_s > 0.0 &&
         config.max_dt_s >= config.min_dt_s &&
         config.initial_dt_s >= config.min_dt_s &&
         config.initial_dt_s <= config.max_dt_s &&
         config.cfl_target == 0.5 &&
         config.diffusion_number_target == 0.25 &&
         config.growth_factor == 1.25 && config.retry_factor == 0.5 &&
         config.max_retries == 8;
}

bool controls_valid_impl(const linear::SolveControl &control) noexcept {
  return std::isfinite(control.atol) && control.atol >= 0.0 &&
         std::isfinite(control.rtol) && control.rtol >= 0.0 &&
         control.max_iterations > 0U &&
         control.residual_recompute_interval > 0U;
}

bool valid_order(MomentumTimeOrder order) noexcept {
  return order == MomentumTimeOrder::backward_euler ||
         order == MomentumTimeOrder::bdf2;
}

runtime::CollectiveStatus
controller_collective_status(const runtime::MpiContext &mpi, bool local_ok,
                             std::string_view operation) {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  if (trusted_tail_active.load(std::memory_order_relaxed))
    trusted_tail_controller_collectives.fetch_add(1U,
                                                  std::memory_order_relaxed);
#endif
  return runtime::collective_status(mpi, local_ok, operation);
}

void record_controller_traversal() noexcept {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  if (trusted_tail_active.load(std::memory_order_relaxed))
    trusted_tail_field_state_traversals.fetch_add(
        1U, std::memory_order_relaxed);
#endif
}

void checked_time_control_mpi(int result, std::string_view operation) {
  runtime::check_mpi_result(result, operation);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  if (trusted_tail_active.load(std::memory_order_relaxed))
    trusted_tail_controller_collectives.fetch_add(1U,
                                                  std::memory_order_relaxed);
  const auto ordinal =
      time_control_raw_count.fetch_add(1U, std::memory_order_relaxed) + 1U;
  if (time_control_raw_fault_ordinal.load(std::memory_order_relaxed) ==
      ordinal)
    runtime::check_mpi_result(MPI_ERR_OTHER, operation);
#endif
}

struct AgreementFrame final {
  static constexpr std::size_t capacity = 64U;
  std::array<std::uint64_t, capacity> words{};
  std::size_t size{};
  bool valid{true};

  void add(std::uint64_t value) noexcept {
    if (size == capacity) {
      valid = false;
      return;
    }
    words[size++] = value;
  }
  void add(std::uint32_t value) noexcept { add(std::uint64_t{value}); }
  void add(bool value) noexcept { add(value ? std::uint64_t{1} : 0U); }
  void add(double value) noexcept { add(bits(value)); }
};

static_assert(AgreementFrame::capacity <=
              static_cast<std::size_t>(std::numeric_limits<int>::max()));
static_assert(21U <= AgreementFrame::capacity,
              "largest frozen Task 22 agreement frame must fit inline");

bool fixed_frame_agrees(const runtime::MpiContext &mpi,
                        const AgreementFrame &local,
                        std::string_view operation) {
  const auto prepared = controller_collective_status(
      mpi, local.valid && local.size <= AgreementFrame::capacity,
      "time-control.agreement.preparation");
  if (!prepared.ok)
    return false;
  std::uint64_t root_size = static_cast<std::uint64_t>(local.size);
  checked_time_control_mpi(
      MPI_Bcast(&root_size, 1, MPI_UINT64_T, 0, mpi.comm()), operation);
  const bool root_size_valid =
      root_size <= AgreementFrame::capacity &&
      root_size <=
          static_cast<std::uint64_t>(std::numeric_limits<int>::max());
  const auto root_prepared = controller_collective_status(
      mpi, root_size_valid, "time-control.agreement.root-size");
  if (!root_prepared.ok)
    return false;
  std::array<std::uint64_t, AgreementFrame::capacity> root{};
  if (mpi.rank() == 0)
    std::copy_n(local.words.begin(), local.size, root.begin());
  checked_time_control_mpi(
      MPI_Bcast(root.data(), static_cast<int>(root_size), MPI_UINT64_T, 0,
                mpi.comm()),
      operation);
  bool equal = static_cast<std::size_t>(root_size) == local.size;
  const auto common =
      std::min(static_cast<std::size_t>(root_size), local.size);
  for (std::size_t i = 0; i < common; ++i)
    equal = (root[i] == local.words[i]) && equal;
  return equal;
}

AgreementFrame config_frame(const config::FlowTimeConfig &config,
                            config::DensityModel model) {
  AgreementFrame frame;
  frame.add(static_cast<std::uint32_t>(config.mode));
  frame.add(static_cast<std::uint64_t>(config.steps));
  frame.add(config.initial_dt_s);
  frame.add(config.min_dt_s);
  frame.add(config.max_dt_s);
  frame.add(config.cfl_target);
  frame.add(config.diffusion_number_target);
  frame.add(config.growth_factor);
  frame.add(config.retry_factor);
  frame.add(static_cast<std::uint64_t>(config.max_retries));
  frame.add(static_cast<std::uint32_t>(model));
  return frame;
}

void append_control(AgreementFrame &frame,
                    const linear::SolveControl &control) {
  frame.add(control.atol);
  frame.add(control.rtol);
  frame.add(control.max_iterations);
  frame.add(control.residual_recompute_interval);
}

void append_metadata(AgreementFrame &frame,
                     const AcceptedStepMetadata &metadata) {
  frame.add(metadata.step);
  frame.add(metadata.time_s);
  frame.add(metadata.dt_s);
  frame.add(metadata.previous_dt_s);
  frame.add(static_cast<std::uint32_t>(metadata.order));
}

void append_state(AgreementFrame &frame, const TimeControlState &state) {
  frame.add(state.schema_version);
  frame.add(state.accepted_step);
  frame.add(state.proposed_next_dt_s);
  frame.add(state.last_accepted_dt_s);
  frame.add(static_cast<std::uint32_t>(state.last_accepted_order));
  frame.add(state.history_ready);
  frame.add(state.last_all_linear_solves_within_half_limit);
  frame.add(state.last_convective_rate_per_s);
  frame.add(state.last_diffusive_rate_per_s);
  frame.add(state.last_stability_metrics_available);
  frame.add(state.last_retry_count);
  frame.add(state.revision);
  frame.add(state.state_seal);
}

AgreementFrame assembly_frame(const mesh::MeshTopology &topology,
                              const mesh::MeshGeometry &geometry,
                              const runtime::MpiContext &mpi) {
  AgreementFrame frame;
  frame.add(static_cast<std::uint64_t>(mpi.size()));
  const auto extent = topology.global_extent();
  frame.add(static_cast<std::uint64_t>(extent.x));
  frame.add(static_cast<std::uint64_t>(extent.y));
  frame.add(static_cast<std::uint64_t>(extent.z));
  frame.add(topology.global_cell_count());
  frame.add(topology.global_face_count());
  frame.add(static_cast<std::uint32_t>(geometry.mapping_kind()));
  const auto origin = geometry.origin_m();
  const auto length = geometry.length_m();
  frame.add(origin.x);
  frame.add(origin.y);
  frame.add(origin.z);
  frame.add(length.x);
  frame.add(length.y);
  frame.add(length.z);
  return frame;
}

AgreementFrame state_frame(const TimeControlState &state,
                           const AcceptedStepMetadata &metadata) {
  AgreementFrame frame;
  append_state(frame, state);
  append_metadata(frame, metadata);
  return frame;
}

bool expected_proposal_is_finite(const config::FlowTimeConfig &config,
                                 const TimeControlState &state,
                                 double &value) noexcept {
  if (config.mode == config::TimeMode::fixed) {
    value = config.initial_dt_s;
    return std::isfinite(value) && value > 0.0;
  }
  const double growth =
      state.last_all_linear_solves_within_half_limit ? config.growth_factor
                                                     : 1.0;
  const double grown = state.last_accepted_dt_s * growth;
  if (!std::isfinite(grown) || grown <= 0.0)
    return false;
  value = std::min(config.max_dt_s, grown);
  if (state.last_convective_rate_per_s > 0.0) {
    const double candidate =
        config.cfl_target / state.last_convective_rate_per_s;
    if (!std::isfinite(candidate) || candidate <= 0.0)
      return false;
    value = std::min(value, candidate);
  }
  if (state.last_diffusive_rate_per_s > 0.0) {
    const double candidate =
        config.diffusion_number_target / state.last_diffusive_rate_per_s;
    if (!std::isfinite(candidate) || candidate <= 0.0)
      return false;
    value = std::min(value, candidate);
  }
  value = std::max(value, config.min_dt_s);
  return std::isfinite(value) && value > 0.0;
}

StepAttemptReport const &base_report(const StepAttemptReport &report) {
  return report;
}
StepAttemptReport const &
base_report(const MaterialDensityStepAttemptReport &report) {
  return report.flow();
}
StepAttemptReport const &base_report(const IdealGasStepAttemptReport &report) {
  return report.flow().flow();
}

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
template <class Report>
void apply_post_return_mutation(Report &report, int rank) noexcept {
  if (rank != time_control_fault_rank.load(std::memory_order_relaxed))
    return;
  auto &base = const_cast<StepAttemptReport &>(base_report(report));
  switch (post_return_mutation.load(std::memory_order_relaxed)) {
  case test::TimeControlPostReturnMutation::none:
    return;
  case test::TimeControlPostReturnMutation::attempted_dt:
    base.attempted_dt_s = std::nextafter(
        base.attempted_dt_s, std::numeric_limits<double>::infinity());
    return;
  case test::TimeControlPostReturnMutation::reason:
    base.reason = StepFailureReason::invalid_input;
    return;
  case test::TimeControlPostReturnMutation::momentum_x_iterations:
    base.momentum.components[0].iterations =
        post_return_iteration_value.load(std::memory_order_relaxed);
    return;
  case test::TimeControlPostReturnMutation::momentum_y_iterations:
    base.momentum.components[1].iterations =
        post_return_iteration_value.load(std::memory_order_relaxed);
    return;
  case test::TimeControlPostReturnMutation::momentum_z_iterations:
    base.momentum.components[2].iterations =
        post_return_iteration_value.load(std::memory_order_relaxed);
    return;
  case test::TimeControlPostReturnMutation::pressure_one_iterations:
    base.pressure[0].iterations =
        post_return_iteration_value.load(std::memory_order_relaxed);
    return;
  case test::TimeControlPostReturnMutation::pressure_two_iterations:
    base.pressure[1].iterations =
        post_return_iteration_value.load(std::memory_order_relaxed);
    return;
  }
}

bool preflight_fault_here(test::TimeControlPreflightFault fault,
                          int rank) noexcept {
  return rank == time_control_fault_rank.load(std::memory_order_relaxed) &&
         preflight_fault.load(std::memory_order_relaxed) == fault;
}

template <class Facade>
void configure_scheduled_attempt(Facade &facade, int rank,
                                 int local_rank) noexcept {
  test::ConstantDensityPisoTestAccess::reset();
  test::MaterialDensityTransportTestAccess::reset();
  const auto remaining =
      scheduled_recoverable_failures.load(std::memory_order_relaxed);
  if (remaining == 0U)
    return;
  const auto reason =
      scheduled_recoverable_reason.load(std::memory_order_relaxed);
  if constexpr (std::is_same_v<Facade, FixedStepIdealGasFlow>) {
    switch (reason) {
    case StepFailureReason::boundary_backflow:
      test::IdealGasClosureTestAccess::set_outlet_backflow_fault(facade);
      break;
    case StepFailureReason::transport_failure:
      test::IdealGasClosureTestAccess::set_post_assessment_fault(
          facade, test::IdealGasPostAssessmentFault::non_finite_state, rank);
      break;
    case StepFailureReason::final_transport_residual:
      test::IdealGasClosureTestAccess::set_post_assessment_fault(
          facade, test::IdealGasPostAssessmentFault::final_transport_residual,
          rank);
      break;
    case StepFailureReason::final_conservation_defect:
      test::IdealGasClosureTestAccess::set_post_assessment_fault(
          facade,
          test::IdealGasPostAssessmentFault::final_conservation_defect, rank);
      break;
    case StepFailureReason::none:
    case StepFailureReason::density_closure_failure:
      test::IdealGasClosureTestAccess::set_controlled_allocation(facade, rank);
      break;
    default:
      break;
    }
  } else if constexpr (std::is_same_v<Facade,
                                      FixedStepMaterialDensityFlow>) {
    if (reason == StepFailureReason::final_transport_residual)
      test::MaterialDensityTransportTestAccess::set_transport_residual(
          0U, 1.0, rank);
    else if (reason == StepFailureReason::none ||
             reason == StepFailureReason::final_continuity_residual)
      test::MaterialDensityTransportTestAccess::set_density_residual(1.0,
                                                                    rank);
  } else {
    switch (reason) {
    case StepFailureReason::final_momentum_residual:
      if (local_rank == rank)
        test::ConstantDensityPisoTestAccess::
            force_final_momentum_perturbation(0U, 1.0e-4);
      break;
    case StepFailureReason::final_transport_residual:
      if (local_rank == rank)
        test::ConstantDensityPisoTestAccess::
            force_final_transport_perturbation(0U, 1.0e-4);
      break;
    case StepFailureReason::final_conservation_defect:
      test::ConstantDensityPisoTestAccess::force_final_conservation_failure(
          local_rank == rank);
      break;
    case StepFailureReason::final_continuity_residual:
      test::ConstantDensityPisoTestAccess::force_final_continuity_failure(
          local_rank == rank);
      break;
    case StepFailureReason::final_pressure_residual:
      test::ConstantDensityPisoTestAccess::force_final_pressure_failure(
          local_rank == rank);
      break;
    case StepFailureReason::none:
    case StepFailureReason::non_finite_trial:
      test::ConstantDensityPisoTestAccess::set_attempt_failure_stage(
          test::AttemptFailureStage::after_begin);
      break;
    default:
      break;
    }
  }
}

void consume_scheduled_attempt() noexcept {
  auto remaining =
      scheduled_recoverable_failures.load(std::memory_order_relaxed);
  while (remaining != 0U &&
         !scheduled_recoverable_failures.compare_exchange_weak(
             remaining, remaining - 1U, std::memory_order_relaxed)) {
  }
  test::ConstantDensityPisoTestAccess::reset();
  test::MaterialDensityTransportTestAccess::reset();
}
#endif

bool recoverable_reason(StepFailureReason reason) noexcept {
  switch (reason) {
  case StepFailureReason::momentum_linear_solve:
  case StepFailureReason::pressure_linear_solve:
  case StepFailureReason::non_finite_trial:
  case StepFailureReason::boundary_backflow:
  case StepFailureReason::transport_failure:
  case StepFailureReason::final_momentum_residual:
  case StepFailureReason::final_transport_residual:
  case StepFailureReason::final_conservation_defect:
  case StepFailureReason::final_continuity_residual:
  case StepFailureReason::final_pressure_residual:
  case StepFailureReason::density_closure_failure:
    return true;
  default:
    return false;
  }
}

bool report_semantically_valid(const StepAttemptReport &report,
                               double attempted_dt) noexcept {
  const bool committed =
      report.disposition == StepAttemptDisposition::committed;
  const bool recoverable =
      report.disposition == StepAttemptDisposition::recoverable_failure;
  const bool non_retryable =
      report.disposition == StepAttemptDisposition::non_retryable_failure;
  if (!(committed || recoverable || non_retryable) ||
      bits(report.attempted_dt_s) != bits(attempted_dt))
    return false;
  if (committed)
    return report.reason == StepFailureReason::none &&
           report.lowest_failing_rank == -1;
  if (report.reason == StepFailureReason::none ||
      report.lowest_failing_rank < 0)
    return false;
  if (recoverable)
    return recoverable_reason(report.reason);
  return !recoverable_reason(report.reason) ||
         report.reason == StepFailureReason::collective_operation;
}

bool half_work_gate(const StepAttemptReport &report,
                    const linear::SolveControl &momentum,
                    const linear::SolveControl &pressure) noexcept {
  if (report.disposition != StepAttemptDisposition::committed)
    return false;
  const auto momentum_limit = momentum.max_iterations / 2U;
  const auto pressure_limit = pressure.max_iterations / 2U;
  return std::all_of(report.momentum.components.begin(),
                     report.momentum.components.end(),
                     [=](const auto &solve) {
                       return solve.iterations <= momentum_limit;
                     }) &&
         std::all_of(report.pressure.begin(), report.pressure.end(),
                     [=](const auto &solve) {
                       return solve.iterations <= pressure_limit;
                     });
}

std::string ranked_message(std::string_view prefix, int rank) {
  return std::string(prefix) + ": lowest failing rank " +
         std::to_string(rank);
}

} // namespace

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
extern "C" void hundun_task22_record_allocation_attempt() noexcept {
  if (trusted_tail_active.load(std::memory_order_relaxed))
    trusted_tail_allocation_attempts.fetch_add(1U,
                                               std::memory_order_relaxed);
}
#endif

bool detail::time_control_config_valid(
    const config::FlowTimeConfig &config) noexcept {
  return valid_config(config);
}

bool detail::density_model_valid(config::DensityModel model) noexcept {
  return valid_density_model(model);
}

bool detail::solve_control_valid(
    const linear::SolveControl &control) noexcept {
  return controls_valid_impl(control);
}

detail::TransportDiffusivityAuthority
detail::make_transport_diffusivity_authority(
    const std::vector<double> &coefficients) {
  TransportDiffusivityAuthority result;
  result.count = coefficients.size();
  for (double value : coefficients) {
    std::uint64_t value_bits = bits(value);
    for (std::size_t i = 0; i < sizeof(value_bits); ++i) {
      result.ordered_fingerprint ^=
          static_cast<std::uint8_t>(value_bits >> (8U * i));
      result.ordered_fingerprint *= kFnvPrime;
    }
    if (!std::isfinite(value) || value < 0.0)
      result.maximum = std::numeric_limits<double>::quiet_NaN();
    else if (std::isfinite(result.maximum))
      result.maximum = std::max(result.maximum, value);
  }
  return result;
}

std::uint64_t detail::TimeControlStateCodec::seal(
    const config::FlowTimeConfig &config, config::DensityModel model,
    const TimeControlState &state) noexcept {
  std::uint64_t hash = kFnvOffset;
  for (char ch : std::string_view(kSealDomain))
    byte(hash, static_cast<std::uint8_t>(ch));
  little_endian(hash, static_cast<std::uint32_t>(config.mode));
  fp64(hash, config.initial_dt_s);
  fp64(hash, config.min_dt_s);
  fp64(hash, config.max_dt_s);
  fp64(hash, config.cfl_target);
  fp64(hash, config.diffusion_number_target);
  fp64(hash, config.growth_factor);
  fp64(hash, config.retry_factor);
  little_endian(hash, static_cast<std::uint32_t>(config.max_retries));
  little_endian(hash, static_cast<std::uint32_t>(model));
  little_endian(hash, state.schema_version);
  little_endian(hash, state.accepted_step);
  fp64(hash, state.proposed_next_dt_s);
  fp64(hash, state.last_accepted_dt_s);
  byte(hash, static_cast<std::uint8_t>(state.last_accepted_order));
  byte(hash, state.history_ready ? 1U : 0U);
  byte(hash, state.last_all_linear_solves_within_half_limit ? 1U : 0U);
  fp64(hash, state.last_convective_rate_per_s);
  fp64(hash, state.last_diffusive_rate_per_s);
  byte(hash, state.last_stability_metrics_available ? 1U : 0U);
  little_endian(hash, state.last_retry_count);
  little_endian(hash, state.revision);
  return hash;
}

bool detail::TimeControlStateCodec::semantically_valid(
    const config::FlowTimeConfig &config, config::DensityModel model,
    const AcceptedStepMetadata &metadata,
    const TimeControlState &state) noexcept {
  if (!valid_config(config) || state.schema_version != 1U ||
      !valid_density_model(model) || !valid_order(state.last_accepted_order) ||
      !valid_order(metadata.order) ||
      state.accepted_step != metadata.step ||
      state.revision != state.accepted_step ||
      state.history_ready != (state.accepted_step != 0U) ||
      state.last_retry_count > 8U ||
      !std::isfinite(state.proposed_next_dt_s) ||
      state.proposed_next_dt_s < config.min_dt_s ||
      state.proposed_next_dt_s > config.max_dt_s ||
      !std::isfinite(metadata.time_s) || metadata.time_s < 0.0 ||
      state.state_seal != seal(config, model, state))
    return false;
  if (state.accepted_step == 0U)
    return positive_zero(metadata.time_s) && positive_zero(metadata.previous_dt_s) &&
           metadata.order == MomentumTimeOrder::backward_euler &&
           bits(metadata.dt_s) == bits(config.initial_dt_s) &&
           positive_zero(state.last_accepted_dt_s) &&
           state.last_accepted_order == MomentumTimeOrder::backward_euler &&
           !state.last_all_linear_solves_within_half_limit &&
           !state.last_stability_metrics_available &&
           positive_zero(state.last_convective_rate_per_s) &&
           positive_zero(state.last_diffusive_rate_per_s) &&
           state.last_retry_count == 0U &&
           bits(state.proposed_next_dt_s) == bits(config.initial_dt_s);
  if (!(metadata.time_s > 0.0) || !(metadata.dt_s > 0.0) ||
      !(metadata.previous_dt_s > 0.0) ||
      bits(state.last_accepted_dt_s) != bits(metadata.dt_s) ||
      state.last_accepted_order != metadata.order)
    return false;
  const double ratio = metadata.dt_s / metadata.previous_dt_s;
  const auto expected_order =
      state.accepted_step == 1U || ratio < 0.5 || ratio > 2.0
          ? MomentumTimeOrder::backward_euler
          : MomentumTimeOrder::bdf2;
  if (metadata.order != expected_order)
    return false;
  if (config.mode == config::TimeMode::fixed)
    return !state.last_stability_metrics_available &&
           positive_zero(state.last_convective_rate_per_s) &&
           positive_zero(state.last_diffusive_rate_per_s) &&
           bits(state.proposed_next_dt_s) == bits(config.initial_dt_s);
  double expected{};
  return state.last_stability_metrics_available &&
         std::isfinite(state.last_convective_rate_per_s) &&
         state.last_convective_rate_per_s >= 0.0 &&
         std::isfinite(state.last_diffusive_rate_per_s) &&
         state.last_diffusive_rate_per_s >= 0.0 &&
         expected_proposal_is_finite(config, state, expected) &&
         bits(state.proposed_next_dt_s) == bits(expected);
}

std::string detail::render_owned_cells(runtime::Box3 box) {
  return "cell.f64.owned." + std::to_string(box.begin.x) + "." +
         std::to_string(box.begin.y) + "." + std::to_string(box.begin.z) +
         "." + std::to_string(box.end.x) + "." + std::to_string(box.end.y) +
         "." + std::to_string(box.end.z);
}
std::string detail::render_global_cells(runtime::Int3 extent) {
  return "cell.f64.global." + std::to_string(extent.x) + "." +
         std::to_string(extent.y) + "." + std::to_string(extent.z);
}
std::string detail::render_owned_faces(std::size_t count) {
  return "face.f64.owned." + std::to_string(count);
}
std::string detail::render_global_faces(std::uint64_t count) {
  return "face.f64.global." + std::to_string(count);
}

TimeAdvanceReport::TimeAdvanceReport(ConstantReportTag) noexcept
    : final_attempt_(StepAttemptReport{}) {}
TimeAdvanceReport::TimeAdvanceReport(MaterialReportTag) noexcept
    : final_attempt_(MaterialDensityStepAttemptReport{}) {}
TimeAdvanceReport::TimeAdvanceReport(IdealGasReportTag) noexcept
    : final_attempt_(IdealGasStepAttemptReport{}) {}
bool TimeAdvanceReport::report_authenticated(
    const StepAttemptReport &) noexcept {
  return true;
}
bool TimeAdvanceReport::report_authenticated(
    const MaterialDensityStepAttemptReport &report) noexcept {
  return report.authenticated() && report.semantic_valid();
}
bool TimeAdvanceReport::report_authenticated(
    const IdealGasStepAttemptReport &report) noexcept {
  return report.authenticated() && report.semantic_valid();
}
TimeAdvanceReport::~TimeAdvanceReport() noexcept = default;
TimeAdvanceReport::TimeAdvanceReport(TimeAdvanceReport &&other) noexcept
    : attempts_(other.attempts_), attempt_count_(other.attempt_count_),
      final_attempt_(std::move(other.final_attempt_)),
      final_attempt_available_(other.final_attempt_available_),
      disposition_(other.disposition_), reason_(other.reason_),
      lowest_failing_rank_(other.lowest_failing_rank_),
      accepted_dt_s_(other.accepted_dt_s_),
      proposed_next_dt_s_(other.proposed_next_dt_s_),
      convective_rate_per_s_(other.convective_rate_per_s_),
      diffusive_rate_per_s_(other.diffusive_rate_per_s_),
      stability_metrics_available_(other.stability_metrics_available_),
      limited_by_min_dt_(other.limited_by_min_dt_),
      controller_identity_(other.controller_identity_),
      report_identity_(other.report_identity_),
      observed_flow_state_identity_(other.observed_flow_state_identity_),
      observed_step_(other.observed_step_),
      observed_time_s_(other.observed_time_s_),
      observed_metadata_(other.observed_metadata_),
      preflight_category_(other.preflight_category_) {
  other.reset_moved_from();
}
void TimeAdvanceReport::reset_moved_from() noexcept {
  attempts_ = {};
  attempt_count_ = 0U;
  final_attempt_available_ = false;
  disposition_ = TimeAdvanceDisposition::non_retryable_failure;
  reason_ = StepFailureReason::invalid_input;
  lowest_failing_rank_ = -1;
  accepted_dt_s_ = proposed_next_dt_s_ = convective_rate_per_s_ =
      diffusive_rate_per_s_ = 0.0;
  stability_metrics_available_ = limited_by_min_dt_ = false;
  controller_identity_ = report_identity_ = observed_flow_state_identity_ =
      observed_step_ = 0U;
  observed_time_s_ = 0.0;
  observed_metadata_ = {};
  preflight_category_ = PreflightCategory::none;
}
TimeAdvanceDisposition TimeAdvanceReport::disposition() const noexcept {
  return disposition_;
}
StepFailureReason TimeAdvanceReport::reason() const noexcept { return reason_; }
int TimeAdvanceReport::lowest_failing_rank() const noexcept {
  return lowest_failing_rank_;
}
std::size_t TimeAdvanceReport::attempt_count() const noexcept {
  return attempt_count_;
}
const TimeAttemptSummary &TimeAdvanceReport::attempt(std::size_t index) const {
  if (index >= attempt_count_)
    throw runtime::Error("time advance attempt index is out of range");
  return attempts_[index];
}
bool TimeAdvanceReport::final_attempt_available() const noexcept {
  return final_attempt_available_;
}
const DensityStepAttemptReport &TimeAdvanceReport::final_attempt() const {
  if (!final_attempt_available_)
    throw runtime::Error("time advance final attempt is unavailable");
  return final_attempt_;
}
double TimeAdvanceReport::accepted_dt_s() const noexcept {
  return accepted_dt_s_;
}
double TimeAdvanceReport::proposed_next_dt_s() const noexcept {
  return proposed_next_dt_s_;
}
double TimeAdvanceReport::convective_rate_per_s() const noexcept {
  return convective_rate_per_s_;
}
double TimeAdvanceReport::diffusive_rate_per_s() const noexcept {
  return diffusive_rate_per_s_;
}
bool TimeAdvanceReport::stability_metrics_available() const noexcept {
  return stability_metrics_available_;
}
bool TimeAdvanceReport::limited_by_min_dt() const noexcept {
  return limited_by_min_dt_;
}

TimeControlDiagnosticSource::TimeControlDiagnosticSource(
    std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
TimeControlDiagnosticSource::~TimeControlDiagnosticSource() noexcept = default;
TimeControlDiagnosticSource::TimeControlDiagnosticSource(
    TimeControlDiagnosticSource &&) noexcept = default;

struct Bdf2RetryController::Impl final {
  config::FlowTimeConfig config;
  config::DensityModel model;
  const mesh::MeshTopology *topology{};
  const mesh::MeshGeometry *geometry{};
  const runtime::MpiContext *mpi{};
  const FlowState *state{};
  std::uint64_t identity{1U};
  std::uint64_t report_identity{};
  bool active{};
};

namespace {

TimeControlState fresh_state(const config::FlowTimeConfig &config,
                             config::DensityModel model) noexcept {
  TimeControlState state;
  state.proposed_next_dt_s = config.initial_dt_s;
  state.state_seal = detail::TimeControlStateCodec::seal(config, model, state);
  return state;
}

void collective_factory_check(const runtime::MpiContext &mpi, bool local_ok,
                              std::string_view prefix) {
  const auto status = controller_collective_status(mpi, local_ok, prefix);
  if (!status.ok)
    throw runtime::Error(ranked_message(prefix, status.failing_rank));
}

std::uint64_t allocate_controller_identity(const runtime::MpiContext &mpi,
                                           std::string_view prefix) {
  static std::atomic<std::uint64_t> next{1U};
  std::uint64_t identity{};
  bool local_ok = true;
  if (mpi.rank() == 0) {
    identity = next.fetch_add(1U, std::memory_order_relaxed);
    local_ok =
        identity != 0U &&
        identity != std::numeric_limits<std::uint64_t>::max();
  }
  collective_factory_check(mpi, local_ok, prefix);
  checked_time_control_mpi(
      MPI_Bcast(&identity, 1, MPI_UINT64_T, 0, mpi.comm()),
      "MPI_Bcast(time-control controller identity)");
  if (identity == 0U ||
      identity == std::numeric_limits<std::uint64_t>::max())
    throw runtime::Error("time-control.create.identity: exhausted");
  return identity;
}

struct Stability final {
  double convection{};
  double diffusion{};
};

struct AdaptivePreparation final {
  std::vector<double> density;
  std::vector<double> flux;
  std::vector<double> flux_sum;
  std::vector<double> geometry_sum;
};

template <class ControllerImpl>
Stability stability_rates(const ControllerImpl &impl,
                          AdaptivePreparation &prepared,
                          double density_constant, bool variable_density,
                          double gamma) {
  record_controller_traversal();
  const auto cells = impl.topology->owned_cell_count();
  if (prepared.density.size() != cells ||
      prepared.flux.size() != impl.topology->local_face_count() ||
      prepared.flux_sum.size() != cells ||
      prepared.geometry_sum.size() != cells)
    throw runtime::Error("time-control stability layout is invalid");
  std::fill(prepared.flux_sum.begin(), prepared.flux_sum.end(), 0.0);
  std::fill(prepared.geometry_sum.begin(), prepared.geometry_sum.end(), 0.0);
  for (mesh::LocalFaceId face = 0; face < impl.topology->local_face_count();
       ++face) {
    const auto owner = impl.topology->owner(face);
    const double mass_flux = prepared.flux[face];
    const auto area = impl.geometry->face_area_m2(face);
    const auto displacement = impl.geometry->face_displacement_m(face);
    const auto area_vector =
        impl.geometry->face_area_vector_m2(face, mesh::FaceSide::owner);
    const double dot = area_vector.x * displacement.x +
                       area_vector.y * displacement.y +
                       area_vector.z * displacement.z;
    const double denominator = std::abs(dot);
    const double area_squared = area * area;
    if (!std::isfinite(mass_flux) || !std::isfinite(area) || area <= 0.0 ||
        !std::isfinite(dot) || !std::isfinite(denominator) ||
        denominator <= 0.0 || !std::isfinite(area_squared) ||
        area_squared <= 0.0)
      throw runtime::Error("time-control committed face state is invalid");
    const double factor = area_squared / denominator;
    if (!std::isfinite(factor) || factor <= 0.0)
      throw runtime::Error("time-control committed face state is invalid");
    const auto accumulate = [&](std::size_t cell) {
      const double next_flux =
          prepared.flux_sum[cell] + std::abs(mass_flux);
      const double next_geometry = prepared.geometry_sum[cell] + factor;
      if (!std::isfinite(next_flux) || !std::isfinite(next_geometry))
        throw runtime::Error("time-control stability accumulation overflow");
      prepared.flux_sum[cell] = next_flux;
      prepared.geometry_sum[cell] = next_geometry;
    };
    if (owner < cells) {
      accumulate(owner);
    }
    const auto neighbour = impl.topology->neighbour(face);
    if (neighbour && *neighbour < cells) {
      accumulate(*neighbour);
    }
  }
  Stability result;
  for (std::size_t cell = 0; cell < cells; ++cell) {
    const double rho =
        variable_density ? prepared.density[cell] : density_constant;
    const double volume = impl.geometry->cell_volume_m3(cell);
    if (!std::isfinite(rho) || rho <= 0.0 || !std::isfinite(volume) ||
        volume <= 0.0)
      throw runtime::Error("time-control committed cell state is invalid");
    const double rho_volume = rho * volume;
    const double denominator = 2.0 * rho_volume;
    const double diffusion_numerator = gamma * prepared.geometry_sum[cell];
    if (!std::isfinite(rho_volume) || rho_volume <= 0.0 ||
        !std::isfinite(denominator) || denominator <= 0.0 ||
        !std::isfinite(diffusion_numerator) ||
        diffusion_numerator < 0.0)
      throw runtime::Error("time-control stability arithmetic overflow");
    const double convection = prepared.flux_sum[cell] / denominator;
    const double diffusion = diffusion_numerator / denominator;
    if (!std::isfinite(convection) || convection < 0.0 ||
        !std::isfinite(diffusion) || diffusion < 0.0)
      throw runtime::Error("time-control stability rate is invalid");
    result.convection = std::max(result.convection, convection);
    result.diffusion = std::max(result.diffusion, diffusion);
  }
  return result;
}

} // namespace

struct detail::AdaptiveTimeControlEngine final {
template <class Report>
static void store_final(TimeAdvanceReport &destination,
                        Report &&source) noexcept {
  using Concrete = std::decay_t<Report>;
  std::get<Concrete>(destination.final_attempt_) =
      std::forward<Report>(source);
  destination.final_attempt_available_ = true;
}

template <class Facade, class Report, class Invoke>
static TimeAdvanceReport run(Bdf2RetryController &controller, FlowState &state,
                              Facade &facade, double density_constant,
                              bool variable_density,
                              config::DensityModel overload_model, double mu,
                              const linear::SolveControl &momentum,
                              const linear::SolveControl &pressure,
                              TimeAdvanceReport report, Invoke &&invoke) {
  auto &impl = *controller.impl_;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  trusted_tail_active.store(false, std::memory_order_relaxed);
  trusted_tail_allocation_attempts.store(0U, std::memory_order_relaxed);
  trusted_tail_controller_collectives.store(0U, std::memory_order_relaxed);
  trusted_tail_field_state_traversals.store(0U, std::memory_order_relaxed);
  trusted_tail_callbacks_or_sinks.store(0U, std::memory_order_relaxed);
#endif
  struct ActiveGuard final {
    bool &active;
    ~ActiveGuard() { active = false; }
  } guard{impl.active};
  report.controller_identity_ = impl.identity;
  report.report_identity_ =
      impl.report_identity == std::numeric_limits<std::uint64_t>::max()
          ? impl.report_identity
          : impl.report_identity + 1U;
  report.proposed_next_dt_s_ = controller.observer_state_.proposed_next_dt_s;
  const auto complete = [&]() noexcept {
    record_controller_traversal();
    const auto metadata = impl.state->metadata();
    report.observed_flow_state_identity_ =
        detail::AdaptiveTimeControlAccess::diagnostic_identity(*impl.state);
    report.observed_step_ = metadata.step;
    report.observed_time_s_ = metadata.time_s;
    report.observed_metadata_ = metadata;
    impl.report_identity = report.report_identity_;
  };
  const auto complete_prepared =
      [&](AcceptedStepMetadata metadata,
          std::uint64_t diagnostic_identity) noexcept {
        report.observed_flow_state_identity_ = diagnostic_identity;
        report.observed_step_ = metadata.step;
        report.observed_time_s_ = metadata.time_s;
        report.observed_metadata_ = metadata;
        impl.report_identity = report.report_identity_;
      };

  bool local_ok = valid_config(impl.config) &&
                  valid_density_model(impl.model) &&
                  controls_valid_impl(momentum) &&
                  controls_valid_impl(pressure) && std::isfinite(mu) &&
                  mu >= 0.0;
  if (!variable_density)
    local_ok = local_ok && std::isfinite(density_constant) &&
               density_constant > 0.0;
  auto status =
      controller_collective_status(*impl.mpi, local_ok,
                                   "time-control.preflight.config");
  if (!status.ok) {
    report.lowest_failing_rank_ = status.failing_rank;
    report.preflight_category_ = TimeAdvanceReport::PreflightCategory::config;
    complete();
    return report;
  }

  const auto facade_identity =
      detail::AdaptiveTimeControlAccess::assembly(facade);
  local_ok = overload_model == impl.model &&
             detail::AdaptiveTimeControlAccess::state_identity(state) ==
                 impl.state &&
             (!facade_identity.live ||
              (facade_identity.topology == impl.topology &&
               facade_identity.geometry == impl.geometry &&
               facade_identity.mpi == impl.mpi));
  status = controller_collective_status(
      *impl.mpi, local_ok, "time-control.preflight.identity");
  if (!status.ok) {
    report.lowest_failing_rank_ = status.failing_rank;
    report.preflight_category_ =
        TimeAdvanceReport::PreflightCategory::identity;
    complete();
    return report;
  }

  local_ok = impl.geometry->compatible(*impl.topology) &&
             detail::AdaptiveTimeControlAccess::state_layout_matches(
                 state, *impl.topology);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  local_ok =
      local_ok &&
      !preflight_fault_here(test::TimeControlPreflightFault::layout,
                            impl.mpi->rank());
#endif
  status = controller_collective_status(
      *impl.mpi, local_ok, "time-control.preflight.layout");
  if (!status.ok) {
    report.lowest_failing_rank_ = status.failing_rank;
    report.preflight_category_ =
        TimeAdvanceReport::PreflightCategory::layout;
    complete();
    return report;
  }

  local_ok = facade_identity.live &&
             detail::AdaptiveTimeControlAccess::state_live(state);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  local_ok =
      local_ok &&
      !preflight_fault_here(test::TimeControlPreflightFault::capability,
                            impl.mpi->rank());
#endif
  status = controller_collective_status(
      *impl.mpi, local_ok, "time-control.preflight.capability");
  if (!status.ok) {
    report.lowest_failing_rank_ = status.failing_rank;
    report.preflight_category_ =
        TimeAdvanceReport::PreflightCategory::capability;
    complete();
    return report;
  }

  AgreementFrame advance_frame = config_frame(impl.config, impl.model);
  advance_frame.add(mu);
  if (!variable_density)
    advance_frame.add(density_constant);
  append_control(advance_frame, momentum);
  append_control(advance_frame, pressure);
  const bool config_agrees = fixed_frame_agrees(
      *impl.mpi, advance_frame, "MPI_Bcast(time-control advance config)");
  status = controller_collective_status(
      *impl.mpi, config_agrees, "time-control.preflight.config");
  if (!status.ok) {
    report.lowest_failing_rank_ = status.failing_rank;
    report.preflight_category_ = TimeAdvanceReport::PreflightCategory::config;
    complete();
    return report;
  }

  const auto metadata_before = state.metadata();
  local_ok = detail::TimeControlStateCodec::semantically_valid(
                 impl.config, impl.model, metadata_before,
                 controller.observer_state_) &&
             controller.observer_state_.accepted_step ==
                 metadata_before.step &&
             controller.observer_state_.revision ==
                 controller.observer_state_.accepted_step;
  status = controller_collective_status(
      *impl.mpi, local_ok, "time-control.preflight.state");
  if (!status.ok) {
    report.lowest_failing_rank_ = status.failing_rank;
    report.preflight_category_ = TimeAdvanceReport::PreflightCategory::state;
    complete();
    return report;
  }
  const bool state_agrees = fixed_frame_agrees(
      *impl.mpi, state_frame(controller.observer_state_, metadata_before),
      "MPI_Bcast(time-control current state)");
  status = controller_collective_status(
      *impl.mpi, state_agrees, "time-control.preflight.state");
  if (!status.ok) {
    report.lowest_failing_rank_ = status.failing_rank;
    report.preflight_category_ = TimeAdvanceReport::PreflightCategory::state;
    complete();
    return report;
  }
  if (impl.report_identity == std::numeric_limits<std::uint64_t>::max() ||
      controller.observer_state_.revision ==
          std::numeric_limits<std::uint64_t>::max() ||
      state.metadata().step == std::numeric_limits<std::uint64_t>::max()) {
    report.lowest_failing_rank_ = 0;
    report.preflight_category_ = TimeAdvanceReport::PreflightCategory::state;
    complete();
    return report;
  }
  const auto &authority = detail::AdaptiveTimeControlAccess::authority(facade);
  const bool authority_ok =
      std::isfinite(authority.maximum) && authority.maximum >= 0.0;
  status = controller_collective_status(
      *impl.mpi, authority_ok,
      "time-control.preflight.transport-authority");
  if (!status.ok) {
    report.lowest_failing_rank_ = status.failing_rank;
    report.preflight_category_ =
        TimeAdvanceReport::PreflightCategory::transport_authority;
    complete();
    return report;
  }
  AgreementFrame authority_frame;
  authority_frame.add(authority.count);
  authority_frame.add(authority.ordered_fingerprint);
  authority_frame.add(authority.maximum);
  authority_frame.add(mu);
  if (!variable_density)
    authority_frame.add(density_constant);
  const bool authority_agrees = fixed_frame_agrees(
      *impl.mpi, authority_frame, "MPI_Bcast(time-control transport)");
  status = controller_collective_status(
      *impl.mpi, authority_agrees,
      "time-control.preflight.transport-authority");
  if (!status.ok) {
    report.lowest_failing_rank_ = status.failing_rank;
    report.preflight_category_ =
        TimeAdvanceReport::PreflightCategory::transport_authority;
    complete();
    return report;
  }

  AdaptivePreparation preparation;
  bool preparation_ok = true;
  bool capability_ok = true;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  preparation_ok =
      !preflight_fault_here(test::TimeControlPreflightFault::preparation,
                            impl.mpi->rank());
#endif
  if (impl.config.mode == config::TimeMode::adaptive) {
    try {
      if (!preparation_ok)
        throw std::bad_alloc();
      preparation.density =
          detail::AdaptiveTimeControlAccess::committed_density(state);
      preparation.flux =
          detail::AdaptiveTimeControlAccess::committed_face_mass_flux(state);
      preparation.flux_sum.resize(impl.topology->owned_cell_count());
      preparation.geometry_sum.resize(impl.topology->owned_cell_count());
    } catch (const std::bad_alloc &) {
      preparation_ok = false;
    } catch (...) {
      capability_ok = false;
    }
  }
  status = controller_collective_status(
      *impl.mpi, capability_ok, "time-control.preflight.capability");
  if (!status.ok) {
    report.lowest_failing_rank_ = status.failing_rank;
    report.preflight_category_ =
        TimeAdvanceReport::PreflightCategory::capability;
    complete();
    return report;
  }
  status = controller_collective_status(
      *impl.mpi, preparation_ok, "time-control.preflight.preparation");
  if (!status.ok) {
    report.lowest_failing_rank_ = status.failing_rank;
    report.preflight_category_ =
        TimeAdvanceReport::PreflightCategory::preparation;
    complete();
    return report;
  }

  Stability rates{};
  if (impl.config.mode == config::TimeMode::adaptive) {
    bool rates_ok = true;
    try {
      rates = stability_rates(impl, preparation, density_constant,
                              variable_density,
                              std::max(mu, authority.maximum));
    } catch (...) {
      rates_ok = false;
    }
    status = controller_collective_status(
        *impl.mpi, rates_ok, "time-control.preflight.state");
    if (!status.ok) {
      report.lowest_failing_rank_ = status.failing_rank;
      report.preflight_category_ = TimeAdvanceReport::PreflightCategory::state;
      complete();
      return report;
    }
    double pair[]{rates.convection, rates.diffusion};
    #ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    stability_reduction_calls.fetch_add(1U, std::memory_order_relaxed);
    stability_reduced_scalars.fetch_add(2U, std::memory_order_relaxed);
    #endif
    impl.mpi->allreduce_fp64_in_place(
        pair, 2U, runtime::Fp64ReductionOperation::maximum);
    rates = {pair[0], pair[1]};
    const bool reduced_rates_ok =
        std::isfinite(rates.convection) && rates.convection >= 0.0 &&
        std::isfinite(rates.diffusion) && rates.diffusion >= 0.0;
    status = controller_collective_status(
        *impl.mpi, reduced_rates_ok, "time-control.preflight.state");
    if (!status.ok) {
      report.lowest_failing_rank_ = status.failing_rank;
      report.preflight_category_ = TimeAdvanceReport::PreflightCategory::state;
      complete();
      return report;
    }
    report.stability_metrics_available_ = true;
    report.convective_rate_per_s_ = rates.convection;
    report.diffusive_rate_per_s_ = rates.diffusion;
  }

  double dt = impl.config.mode == config::TimeMode::fixed
                  ? impl.config.initial_dt_s
                  : std::clamp(controller.observer_state_.proposed_next_dt_s,
                               impl.config.min_dt_s, impl.config.max_dt_s);
  if (impl.config.mode == config::TimeMode::adaptive) {
    double candidate = dt;
    if (rates.convection > 0.0) {
      const double limited = impl.config.cfl_target / rates.convection;
      if (!std::isfinite(limited) || limited <= 0.0) {
        report.lowest_failing_rank_ = 0;
        report.preflight_category_ =
            TimeAdvanceReport::PreflightCategory::state;
        complete();
        return report;
      }
      candidate = std::min(candidate, limited);
    }
    if (rates.diffusion > 0.0) {
      const double limited =
          impl.config.diffusion_number_target / rates.diffusion;
      if (!std::isfinite(limited) || limited <= 0.0) {
        report.lowest_failing_rank_ = 0;
        report.preflight_category_ =
            TimeAdvanceReport::PreflightCategory::state;
        complete();
        return report;
      }
      candidate = std::min(candidate, limited);
    }
    if (!std::isfinite(candidate) || candidate <= 0.0) {
      report.lowest_failing_rank_ = 0;
      report.preflight_category_ =
          TimeAdvanceReport::PreflightCategory::state;
      complete();
      return report;
    }
    report.limited_by_min_dt_ = candidate <= impl.config.min_dt_s;
    dt = std::max(candidate, impl.config.min_dt_s);
  }

  std::uint32_t retry_count{};
  while (true) {
    const auto metadata = metadata_before;
    const double ratio = metadata.step == 0U ? 0.0 : dt / metadata.dt_s;
    const auto order =
        metadata.step != 0U && ratio >= 0.5 && ratio <= 2.0
            ? MomentumTimeOrder::bdf2
            : MomentumTimeOrder::backward_euler;
    const auto stencil =
        make_momentum_time_stencil(order, dt, metadata.dt_s);
    TimeControlState proposal_probe = controller.observer_state_;
    proposal_probe.last_accepted_dt_s = dt;
    proposal_probe.last_convective_rate_per_s = rates.convection;
    proposal_probe.last_diffusive_rate_per_s = rates.diffusion;
    proposal_probe.last_stability_metrics_available =
        impl.config.mode == config::TimeMode::adaptive;
    double slow_proposal{};
    double fast_proposal{};
    proposal_probe.last_all_linear_solves_within_half_limit = false;
    bool proposal_ok =
        expected_proposal_is_finite(impl.config, proposal_probe, slow_proposal);
    proposal_probe.last_all_linear_solves_within_half_limit = true;
    proposal_ok =
        expected_proposal_is_finite(impl.config, proposal_probe, fast_proposal) &&
        proposal_ok;
    const AcceptedStepMetadata accepted_metadata{
        metadata.step + 1U, metadata.time_s + dt, dt, metadata.dt_s, order};
    const auto current_diagnostic_identity =
        detail::AdaptiveTimeControlAccess::diagnostic_identity(state);
    proposal_ok =
        proposal_ok && std::isfinite(accepted_metadata.time_s) &&
        accepted_metadata.time_s > metadata.time_s &&
        current_diagnostic_identity !=
            std::numeric_limits<std::uint64_t>::max();
    TimeControlState slow_state = controller.observer_state_;
    slow_state.accepted_step = accepted_metadata.step;
    slow_state.last_accepted_dt_s = dt;
    slow_state.last_accepted_order = order;
    slow_state.history_ready = true;
    slow_state.last_all_linear_solves_within_half_limit = false;
    slow_state.last_convective_rate_per_s =
        impl.config.mode == config::TimeMode::adaptive ? rates.convection
                                                      : 0.0;
    slow_state.last_diffusive_rate_per_s =
        impl.config.mode == config::TimeMode::adaptive ? rates.diffusion
                                                      : 0.0;
    slow_state.last_stability_metrics_available =
        impl.config.mode == config::TimeMode::adaptive;
    slow_state.last_retry_count = retry_count;
    slow_state.revision = controller.observer_state_.revision + 1U;
    slow_state.proposed_next_dt_s = slow_proposal;
    slow_state.state_seal =
        detail::TimeControlStateCodec::seal(impl.config, impl.model,
                                            slow_state);
    TimeControlState fast_state = slow_state;
    fast_state.last_all_linear_solves_within_half_limit = true;
    fast_state.proposed_next_dt_s = fast_proposal;
    fast_state.state_seal =
        detail::TimeControlStateCodec::seal(impl.config, impl.model,
                                            fast_state);
    status = controller_collective_status(
        *impl.mpi, proposal_ok, "time-control.preflight.state");
    if (!status.ok) {
      report.lowest_failing_rank_ = status.failing_rank;
      report.preflight_category_ = TimeAdvanceReport::PreflightCategory::state;
      complete();
      return report;
    }
 #ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    configure_scheduled_attempt(facade,
                                time_control_fault_rank.load(
                                    std::memory_order_relaxed),
                                impl.mpi->rank());
 #endif
    Report attempt = invoke(stencil);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (base_report(attempt).disposition ==
        StepAttemptDisposition::committed)
      trusted_tail_active.store(true, std::memory_order_relaxed);
#endif
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    consume_scheduled_attempt();
    apply_post_return_mutation(attempt, impl.mpi->rank());
#endif
    const auto &base = base_report(attempt);
    const auto disposition = base.disposition;
    const auto reason = base.reason;
    const auto failing_rank = base.lowest_failing_rank;
    const bool work_gate = half_work_gate(base, momentum, pressure);
    report.attempts_[report.attempt_count_++] =
        {dt, order, disposition, reason, failing_rank,
         disposition == StepAttemptDisposition::committed && work_gate};

    if (disposition == StepAttemptDisposition::committed) {
      store_final(report, std::move(attempt));
      controller.observer_state_ = work_gate ? fast_state : slow_state;
      report.disposition_ = TimeAdvanceDisposition::committed;
      report.reason_ = StepFailureReason::none;
      report.lowest_failing_rank_ = -1;
      report.accepted_dt_s_ = dt;
      report.proposed_next_dt_s_ =
          controller.observer_state_.proposed_next_dt_s;
      complete_prepared(accepted_metadata, current_diagnostic_identity + 1U);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      trusted_tail_active.store(false, std::memory_order_relaxed);
#endif
      return report;
    }

    const bool normalized =
        TimeAdvanceReport::report_authenticated(attempt) &&
        report_semantically_valid(base, dt);
    status = controller_collective_status(
        *impl.mpi, normalized, "time-control.preflight.report");
    report.reason_ = reason;
    report.lowest_failing_rank_ = failing_rank;
    store_final(report, std::move(attempt));
    if (!status.ok) {
      report.disposition_ = TimeAdvanceDisposition::non_retryable_failure;
      report.reason_ = StepFailureReason::invalid_input;
      report.lowest_failing_rank_ = status.failing_rank;
      report.preflight_category_ = TimeAdvanceReport::PreflightCategory::report;
      complete();
      return report;
    }
    if (disposition != StepAttemptDisposition::recoverable_failure) {
      report.disposition_ = TimeAdvanceDisposition::non_retryable_failure;
      complete();
      return report;
    }
    if (bits(dt) == bits(impl.config.min_dt_s)) {
      report.disposition_ = TimeAdvanceDisposition::minimum_dt_failure;
      report.limited_by_min_dt_ = true;
      complete();
      return report;
    }
    if (retry_count == 8U) {
      report.disposition_ = TimeAdvanceDisposition::retry_limit_reached;
      complete();
      return report;
    }
    ++retry_count;
    const double candidate = impl.config.retry_factor * dt;
    if (!std::isfinite(candidate) || candidate <= 0.0) {
      report.disposition_ = TimeAdvanceDisposition::non_retryable_failure;
      report.reason_ = StepFailureReason::invalid_input;
      report.preflight_category_ = TimeAdvanceReport::PreflightCategory::state;
      complete();
      return report;
    }
    report.limited_by_min_dt_ = candidate <= impl.config.min_dt_s;
    dt = std::max(impl.config.min_dt_s, candidate);
  }
}
};

Bdf2RetryController Bdf2RetryController::create(
    const config::FlowTimeConfig &config, config::DensityModel model,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const runtime::MpiContext &mpi, const FlowState &state) {
  bool local_ok = valid_config(config);
  collective_factory_check(mpi, local_ok,
                           "time-control.create.invalid-config");
  collective_factory_check(mpi, valid_density_model(model),
                           "time-control.create.identity");
  local_ok = geometry.compatible(topology);
  collective_factory_check(mpi, local_ok, "time-control.create.layout");
  collective_factory_check(
      mpi, detail::AdaptiveTimeControlAccess::state_live(state),
      "time-control.create.state");
  auto observer = fresh_state(config, model);
  local_ok = detail::TimeControlStateCodec::semantically_valid(
      config, model, state.metadata(), observer);
  collective_factory_check(mpi, local_ok, "time-control.create.state");
  const auto metadata = state.metadata();
  bool agreement = fixed_frame_agrees(
      mpi, config_frame(config, model),
      "MPI_Bcast(time-control create config)");
  const bool assembly_agreement = fixed_frame_agrees(
      mpi, assembly_frame(topology, geometry, mpi),
      "MPI_Bcast(time-control create assembly)");
  agreement = assembly_agreement && agreement;
  const bool state_agreement = fixed_frame_agrees(
      mpi, state_frame(observer, metadata),
      "MPI_Bcast(time-control create state)");
  agreement = state_agreement && agreement;
  collective_factory_check(mpi, agreement, "time-control.create.agreement");
  const auto identity =
      allocate_controller_identity(mpi, "time-control.create.identity");
  std::unique_ptr<Impl> impl;
  bool construction_ok = true;
  try {
    impl = std::make_unique<Impl>();
  } catch (const std::bad_alloc &) {
    construction_ok = false;
  }
  collective_factory_check(mpi, construction_ok,
                           "time-control.create.state");
  impl->config = config;
  impl->model = model;
  impl->topology = &topology;
  impl->geometry = &geometry;
  impl->mpi = &mpi;
  impl->state = &state;
  impl->identity = identity;
  return Bdf2RetryController(std::move(impl), observer);
}

Bdf2RetryController Bdf2RetryController::restore(
    const config::FlowTimeConfig &config, config::DensityModel model,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const runtime::MpiContext &mpi, const FlowState &state,
    const TimeControlState &observer) {
  collective_factory_check(mpi, valid_config(config),
                           "time-control.restore.invalid-config");
  collective_factory_check(mpi, valid_density_model(model),
                           "time-control.restore.identity");
  collective_factory_check(mpi, geometry.compatible(topology),
                           "time-control.restore.layout");
  collective_factory_check(
      mpi, detail::AdaptiveTimeControlAccess::state_live(state),
      "time-control.restore.state");
  collective_factory_check(
      mpi, detail::TimeControlStateCodec::semantically_valid(
               config, model, state.metadata(), observer),
      "time-control.restore.state");
  const auto metadata = state.metadata();
  bool agreement = fixed_frame_agrees(
      mpi, config_frame(config, model),
      "MPI_Bcast(time-control restore config)");
  const bool assembly_agreement = fixed_frame_agrees(
      mpi, assembly_frame(topology, geometry, mpi),
      "MPI_Bcast(time-control restore assembly)");
  agreement = assembly_agreement && agreement;
  const bool state_agreement = fixed_frame_agrees(
      mpi, state_frame(observer, metadata),
      "MPI_Bcast(time-control restore state)");
  agreement = state_agreement && agreement;
  collective_factory_check(mpi, agreement, "time-control.restore.agreement");
  const auto identity =
      allocate_controller_identity(mpi, "time-control.restore.identity");
  std::unique_ptr<Impl> impl;
  bool construction_ok = true;
  try {
    impl = std::make_unique<Impl>();
  } catch (const std::bad_alloc &) {
    construction_ok = false;
  }
  collective_factory_check(mpi, construction_ok,
                           "time-control.restore.state");
  impl->config = config;
  impl->model = model;
  impl->topology = &topology;
  impl->geometry = &geometry;
  impl->mpi = &mpi;
  impl->state = &state;
  impl->identity = identity;
  return Bdf2RetryController(std::move(impl), observer);
}

Bdf2RetryController::Bdf2RetryController(
    std::unique_ptr<Impl> impl, TimeControlState observer_state) noexcept
    : impl_(std::move(impl)), observer_state_(observer_state) {}
Bdf2RetryController::~Bdf2RetryController() noexcept = default;
Bdf2RetryController::Bdf2RetryController(Bdf2RetryController &&other) noexcept
    : impl_(std::move(other.impl_)), observer_state_(other.observer_state_) {
  other.observer_state_ = {};
}
TimeControlState Bdf2RetryController::state() const noexcept {
  return observer_state_;
}

TimeAdvanceReport Bdf2RetryController::advance(
    FlowState &state, FixedStepConstantDensityFlow &facade, double rho_ref,
    double mu, const linear::SolveControl &momentum,
    const linear::SolveControl &pressure) {
  if (!impl_ || impl_->active)
    throw runtime::Error("time-control controller is moved-from or active");
  if (!detail::AdaptiveTimeControlAccess::state_live(*impl_->state) ||
      !detail::AdaptiveTimeControlAccess::state_live(state) ||
      !detail::AdaptiveTimeControlAccess::assembly(facade).live)
    throw runtime::Error("time-control bound object is moved-from");
  impl_->active = true;
  TimeAdvanceReport result(TimeAdvanceReport::ConstantReportTag{});
  return detail::AdaptiveTimeControlEngine::run<
      FixedStepConstantDensityFlow, StepAttemptReport>(
      *this, state, facade, rho_ref, false, config::DensityModel::constant, mu,
      momentum, pressure,
      std::move(result), [&](const auto &stencil) {
        return facade.attempt(state, rho_ref, mu, stencil, momentum, pressure);
      });
}
TimeAdvanceReport Bdf2RetryController::advance(
    FlowState &state, FixedStepMaterialDensityFlow &facade, double mu,
    const linear::SolveControl &momentum,
    const linear::SolveControl &pressure) {
  if (!impl_ || impl_->active)
    throw runtime::Error("time-control controller is moved-from or active");
  if (!detail::AdaptiveTimeControlAccess::state_live(*impl_->state) ||
      !detail::AdaptiveTimeControlAccess::state_live(state) ||
      !detail::AdaptiveTimeControlAccess::assembly(facade).live)
    throw runtime::Error("time-control bound object is moved-from");
  impl_->active = true;
  TimeAdvanceReport result(TimeAdvanceReport::MaterialReportTag{});
  return detail::AdaptiveTimeControlEngine::run<
      FixedStepMaterialDensityFlow, MaterialDensityStepAttemptReport>(
      *this, state, facade, 0.0, true, config::DensityModel::material, mu,
      momentum, pressure,
      std::move(result), [&](const auto &stencil) {
        return facade.attempt(state, mu, stencil, momentum, pressure);
      });
}
TimeAdvanceReport Bdf2RetryController::advance(
    FlowState &state, FixedStepIdealGasFlow &facade, double mu,
    const linear::SolveControl &momentum,
    const linear::SolveControl &pressure) {
  if (!impl_ || impl_->active)
    throw runtime::Error("time-control controller is moved-from or active");
  if (!detail::AdaptiveTimeControlAccess::state_live(*impl_->state) ||
      !detail::AdaptiveTimeControlAccess::state_live(state) ||
      !detail::AdaptiveTimeControlAccess::assembly(facade).live)
    throw runtime::Error("time-control bound object is moved-from");
  impl_->active = true;
  TimeAdvanceReport result(TimeAdvanceReport::IdealGasReportTag{});
  return detail::AdaptiveTimeControlEngine::run<
      FixedStepIdealGasFlow, IdealGasStepAttemptReport>(
      *this, state, facade, 0.0, true, config::DensityModel::ideal_gas, mu,
      momentum, pressure,
      std::move(result), [&](const auto &stencil) {
        return facade.attempt(state, mu, stencil, momentum, pressure);
      });
}

TimeControlDiagnosticSource Bdf2RetryController::diagnostic_source(
    const FlowState &state, const TimeAdvanceReport &report) const {
  if (!impl_)
    throw runtime::Error("time-control diagnostic source is stale");
  if (!detail::AdaptiveTimeControlAccess::state_live(*impl_->state) ||
      !detail::AdaptiveTimeControlAccess::state_live(state))
    throw runtime::Error("time-control diagnostic state is moved-from");
  const auto metadata = state.metadata();
  if (&state != impl_->state || report.controller_identity_ != impl_->identity ||
      report.report_identity_ == 0U ||
      report.report_identity_ != impl_->report_identity ||
      report.observed_flow_state_identity_ == 0U ||
      report.observed_flow_state_identity_ !=
          detail::AdaptiveTimeControlAccess::diagnostic_identity(state) ||
      report.observed_step_ != metadata.step ||
      bits(report.observed_time_s_) != bits(metadata.time_s) ||
      report.observed_metadata_.step != metadata.step ||
      bits(report.observed_metadata_.time_s) != bits(metadata.time_s) ||
      bits(report.observed_metadata_.dt_s) != bits(metadata.dt_s) ||
      bits(report.observed_metadata_.previous_dt_s) !=
          bits(metadata.previous_dt_s) ||
      report.observed_metadata_.order != metadata.order)
    throw runtime::Error("time-control diagnostic source identity mismatch");
  auto source = std::make_unique<TimeControlDiagnosticSource::Impl>();
  auto &snapshot = source->snapshot;
  snapshot.state = observer_state_;
  snapshot.attempts = report.attempts_;
  snapshot.attempt_count = report.attempt_count_;
  snapshot.disposition = report.disposition_;
  snapshot.reason = report.reason_;
  snapshot.lowest_failing_rank = report.lowest_failing_rank_;
  snapshot.accepted_dt_s = report.accepted_dt_s_;
  snapshot.proposed_next_dt_s = report.proposed_next_dt_s_;
  snapshot.convective_rate_per_s = report.convective_rate_per_s_;
  snapshot.diffusive_rate_per_s = report.diffusive_rate_per_s_;
  snapshot.stability_metrics_available = report.stability_metrics_available_;
  snapshot.limited_by_min_dt = report.limited_by_min_dt_;
  snapshot.preflight_category =
      static_cast<std::uint8_t>(report.preflight_category_);
  snapshot.config = impl_->config;
  snapshot.model = impl_->model;
  snapshot.controller_identity = impl_->identity;
  snapshot.report_identity = report.report_identity_;
  snapshot.flow_state_identity = report.observed_flow_state_identity_;
  snapshot.relative_rank = impl_->mpi->rank();
  snapshot.observed_step = report.observed_step_;
  snapshot.observed_time_s = report.observed_time_s_;
  snapshot.observed_metadata = report.observed_metadata_;
  snapshot.local_box = impl_->topology->owned_global_box();
  snapshot.global_extent = impl_->topology->global_extent();
  snapshot.local_faces = impl_->topology->local_face_count();
  snapshot.global_faces = impl_->topology->global_face_count();
  for (mesh::LocalFaceId face = 0; face < snapshot.local_faces; ++face)
    if (impl_->topology->cell_ownership(impl_->topology->owner(face)) ==
        mesh::EntityOwnership::owned)
      ++snapshot.canonical_owned_faces;
  snapshot.local_cell_layout = detail::render_owned_cells(snapshot.local_box);
  snapshot.global_cell_layout =
      detail::render_global_cells(snapshot.global_extent);
  snapshot.local_face_layout =
      detail::render_owned_faces(snapshot.canonical_owned_faces);
  snapshot.global_face_layout =
      detail::render_global_faces(snapshot.global_faces);
  return TimeControlDiagnosticSource(std::move(source));
}

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
void test::AdaptiveTimeControlTestAccess::set_post_return_mutation(
    test::TimeControlPostReturnMutation mutation, int rank) noexcept {
  post_return_mutation.store(mutation, std::memory_order_relaxed);
  time_control_fault_rank.store(rank, std::memory_order_relaxed);
}

void test::AdaptiveTimeControlTestAccess::set_preflight_fault(
    test::TimeControlPreflightFault fault, int rank) noexcept {
  preflight_fault.store(fault, std::memory_order_relaxed);
  time_control_fault_rank.store(rank, std::memory_order_relaxed);
}

void test::AdaptiveTimeControlTestAccess::set_recoverable_failures(
    std::uint32_t count, int rank) noexcept {
  scheduled_recoverable_failures.store(count, std::memory_order_relaxed);
  time_control_fault_rank.store(rank, std::memory_order_relaxed);
}

void test::AdaptiveTimeControlTestAccess::set_recoverable_failure_reason(
    StepFailureReason reason) noexcept {
  scheduled_recoverable_reason.store(reason, std::memory_order_relaxed);
}

void test::AdaptiveTimeControlTestAccess::set_post_return_iteration_value(
    std::uint64_t value) noexcept {
  post_return_iteration_value.store(value, std::memory_order_relaxed);
}

void test::AdaptiveTimeControlTestAccess::reset_faults() noexcept {
  post_return_mutation.store(test::TimeControlPostReturnMutation::none,
                             std::memory_order_relaxed);
  preflight_fault.store(test::TimeControlPreflightFault::none,
                        std::memory_order_relaxed);
  time_control_fault_rank.store(-1, std::memory_order_relaxed);
  trusted_tail_active.store(false, std::memory_order_relaxed);
  trusted_tail_allocation_attempts.store(0U, std::memory_order_relaxed);
  trusted_tail_controller_collectives.store(0U, std::memory_order_relaxed);
  trusted_tail_field_state_traversals.store(0U, std::memory_order_relaxed);
  trusted_tail_callbacks_or_sinks.store(0U, std::memory_order_relaxed);
  stability_reduction_calls.store(0U, std::memory_order_relaxed);
  stability_reduced_scalars.store(0U, std::memory_order_relaxed);
  time_control_raw_count.store(0U, std::memory_order_relaxed);
  time_control_raw_fault_ordinal.store(0U, std::memory_order_relaxed);
  time_control_raw_fault_rank.store(-1, std::memory_order_relaxed);
  scheduled_recoverable_failures.store(0U, std::memory_order_relaxed);
  scheduled_recoverable_reason.store(StepFailureReason::none,
                                     std::memory_order_relaxed);
  post_return_iteration_value.store(
      std::numeric_limits<std::uint64_t>::max(),
      std::memory_order_relaxed);
  test::ConstantDensityPisoTestAccess::reset();
  test::MaterialDensityTransportTestAccess::reset();
}

std::uint8_t test::AdaptiveTimeControlTestAccess::preflight_category(
    const TimeAdvanceReport &report) noexcept {
  return static_cast<std::uint8_t>(report.preflight_category_);
}

test::TimeControlTrustedTailObservation
test::AdaptiveTimeControlTestAccess::trusted_tail_observation() noexcept {
  return {
      trusted_tail_allocation_attempts.load(std::memory_order_relaxed),
      trusted_tail_controller_collectives.load(std::memory_order_relaxed),
      trusted_tail_field_state_traversals.load(std::memory_order_relaxed),
      trusted_tail_callbacks_or_sinks.load(std::memory_order_relaxed)};
}

void test::AdaptiveTimeControlTestAccess::set_raw_fault(
    std::size_t ordinal, int rank) noexcept {
  time_control_raw_count.store(0U, std::memory_order_relaxed);
  time_control_raw_fault_ordinal.store(ordinal, std::memory_order_relaxed);
  time_control_raw_fault_rank.store(rank, std::memory_order_relaxed);
}

std::size_t
test::AdaptiveTimeControlTestAccess::raw_operation_count() noexcept {
  return time_control_raw_count.load(std::memory_order_relaxed);
}
void test::AdaptiveTimeControlTestAccess::set_active(
    Bdf2RetryController &controller, bool active) noexcept {
  if (controller.impl_)
    controller.impl_->active = active;
}
std::array<double, 2> test::AdaptiveTimeControlTestAccess::stability_rates(
    Bdf2RetryController &controller, const FlowState &state,
    double density_constant, bool variable_density, double gamma) {
  AdaptivePreparation prepared;
  prepared.density = detail::AdaptiveTimeControlAccess::committed_density(state);
  prepared.flux =
      detail::AdaptiveTimeControlAccess::committed_face_mass_flux(state);
  prepared.flux_sum.resize(controller.impl_->topology->owned_cell_count());
  prepared.geometry_sum.resize(controller.impl_->topology->owned_cell_count());
  const auto result = ::hundun::flow::stability_rates(
      *controller.impl_, prepared, density_constant, variable_density, gamma);
  return {result.convection, result.diffusion};
}
std::array<std::uint64_t, 2>
test::AdaptiveTimeControlTestAccess::stability_reduction_observation()
    noexcept {
  return {stability_reduction_calls.load(std::memory_order_relaxed),
          stability_reduced_scalars.load(std::memory_order_relaxed)};
}
#endif

static_assert(
    std::is_nothrow_move_constructible_v<TimeAdvanceReport>);
static_assert(
    std::is_nothrow_move_constructible_v<TimeControlDiagnosticSource>);
static_assert(
    std::is_nothrow_move_constructible_v<Bdf2RetryController>);
static_assert(
    std::is_nothrow_move_assignable_v<MaterialDensityStepAttemptReport>);
static_assert(std::is_nothrow_move_assignable_v<IdealGasStepAttemptReport>);

} // namespace hundun::flow

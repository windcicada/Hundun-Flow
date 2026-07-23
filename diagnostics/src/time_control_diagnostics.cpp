// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/time_control_diagnostics.hpp"

#include "adaptive_time_control_detail.hpp"
#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
#include "time_control_diagnostics_test_access.hpp"
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace hundun::diagnostics {
namespace {

enum class FaultPoint : std::uint8_t {
  none,
  phase1_layout,
  phase2_request,
  projection_root_size,
  projection_payload,
  phase3_provider,
  phase4_payload,
  fingerprint_aggregation,
  eligible_aggregation,
  wire_root_size,
  wire_payload,
  phase4_wire,
  phase5_record,
  phase6_sink
};

#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
test::TimeControlDiagnosticFault diagnostic_fault{
    test::TimeControlDiagnosticFault::none};
int diagnostic_fault_rank{-1};
std::size_t diagnostic_phase_count{};
std::size_t diagnostic_raw_count{};
std::size_t diagnostic_submission_count{};
std::size_t diagnostic_raw_fault_ordinal{};
std::size_t request_mutation_offset{std::numeric_limits<std::size_t>::max()};
int request_mutation_rank{-1};
std::size_t provider_mutation_offset{std::numeric_limits<std::size_t>::max()};
int provider_mutation_rank{-1};
test::TimeControlWireMutation wire_mutation{
    test::TimeControlWireMutation::none};
int wire_mutation_rank{-1};

bool fault_here(FaultPoint point, int rank) noexcept {
  return static_cast<std::uint8_t>(diagnostic_fault) ==
             static_cast<std::uint8_t>(point) &&
         diagnostic_fault_rank == rank;
}
#else
bool fault_here(FaultPoint, int) noexcept { return false; }
#endif

void checked_mpi(int result, std::string_view operation) {
  runtime::check_mpi_result(result, operation);
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  ++diagnostic_raw_count;
  if (diagnostic_raw_fault_ordinal != 0U &&
      diagnostic_raw_count == diagnostic_raw_fault_ordinal)
    runtime::check_mpi_result(MPI_ERR_OTHER, operation);
#endif
}

constexpr DiagnosticCapabilityFlags kCapabilities =
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::summary) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::invariants) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::counters) |
    static_cast<DiagnosticCapabilityFlags>(
        DiagnosticCapability::bounded_state_sample) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::collective);

constexpr std::array<std::string_view, 7> kFields{
    "time-control.accepted-step", "time-control.history-ready",
    "time-control.last-accepted-dt", "time-control.last-order",
    "time-control.last-work-gate", "time-control.next-dt",
    "time-control.revision"};

struct Tuple final {
  std::string_view field;
  std::uint32_t component;
  std::string_view unit;
  double value;
};

std::array<Tuple, 9>
tuples(const flow::detail::TimeControlDiagnosticSnapshot &source) {
  const auto step = source.state.accepted_step;
  const auto revision = source.state.revision;
  return {{{kFields[0], 0U, "count",
            static_cast<double>(static_cast<std::uint32_t>(step))},
           {kFields[0], 1U, "count",
            static_cast<double>(static_cast<std::uint32_t>(step >> 32U))},
           {kFields[1], 0U, "1", source.state.history_ready ? 1.0 : 0.0},
           {kFields[2], 0U, "s", source.state.last_accepted_dt_s},
           {kFields[3], 0U, "count",
            static_cast<double>(
                static_cast<std::uint8_t>(source.state.last_accepted_order))},
           {kFields[4], 0U, "1",
            source.state.last_all_linear_solves_within_half_limit ? 1.0 : 0.0},
           {kFields[5], 0U, "s", source.state.proposed_next_dt_s},
           {kFields[6], 0U, "count",
            static_cast<double>(static_cast<std::uint32_t>(revision))},
           {kFields[6], 1U, "count",
            static_cast<double>(static_cast<std::uint32_t>(revision >> 32U))}}};
}

DiagnosticFailure failure_for(const flow::detail::TimeControlDiagnosticSnapshot &s,
                              bool collective) {
  if (s.disposition == flow::TimeAdvanceDisposition::committed)
    return {};
  DiagnosticFailure result;
  result.lowest_failing_rank = collective ? s.lowest_failing_rank : -1;
  if (s.attempt_count == 0U) {
    switch (s.preflight_category) {
    case 1U:
      return {DiagnosticFailureClass::invalid_input,
              "time-control.preflight.config", result.lowest_failing_rank};
    case 2U:
      return {DiagnosticFailureClass::invalid_input,
              "time-control.preflight.identity", result.lowest_failing_rank};
    case 3U:
      return {DiagnosticFailureClass::layout,
              "time-control.preflight.layout", result.lowest_failing_rank};
    case 4U:
      return {DiagnosticFailureClass::capability,
              "time-control.preflight.capability", result.lowest_failing_rank};
    case 5U:
      return {DiagnosticFailureClass::invalid_input,
              "time-control.preflight.state", result.lowest_failing_rank};
    case 6U:
      return {DiagnosticFailureClass::invalid_input,
              "time-control.preflight.transport-authority",
              result.lowest_failing_rank};
    case 7U:
      return {DiagnosticFailureClass::invalid_input,
              "time-control.preflight.preparation",
              result.lowest_failing_rank};
    case 8U:
      return {DiagnosticFailureClass::invalid_input,
              "time-control.preflight.report", result.lowest_failing_rank};
    default:
      return {DiagnosticFailureClass::invalid_input,
              "time-control.invalid-input", result.lowest_failing_rank};
    }
  }
  switch (s.reason) {
  case flow::StepFailureReason::momentum_linear_solve:
    result = {DiagnosticFailureClass::non_convergence,
              "time-control.momentum-linear-solve",
              result.lowest_failing_rank};
    break;
  case flow::StepFailureReason::pressure_linear_solve:
    result = {DiagnosticFailureClass::non_convergence,
              "time-control.pressure-linear-solve",
              result.lowest_failing_rank};
    break;
  case flow::StepFailureReason::non_finite_trial:
    result = {DiagnosticFailureClass::non_finite_state,
              "time-control.non-finite-trial", result.lowest_failing_rank};
    break;
  case flow::StepFailureReason::boundary_backflow:
    result = {DiagnosticFailureClass::boundary,
              "time-control.boundary-backflow", result.lowest_failing_rank};
    break;
  case flow::StepFailureReason::final_conservation_defect:
    result = {DiagnosticFailureClass::conservation,
              "time-control.final-conservation-defect",
              result.lowest_failing_rank};
    break;
  case flow::StepFailureReason::collective_operation:
    result = {DiagnosticFailureClass::collective_operation,
              "time-control.collective-operation", result.lowest_failing_rank};
    break;
  case flow::StepFailureReason::transport_failure:
  case flow::StepFailureReason::density_closure_failure:
    result = {DiagnosticFailureClass::numerical_breakdown,
              s.reason == flow::StepFailureReason::transport_failure
                  ? "time-control.transport-failure"
                  : "time-control.density-closure-failure",
              result.lowest_failing_rank};
    break;
  case flow::StepFailureReason::final_momentum_residual:
  case flow::StepFailureReason::final_transport_residual:
  case flow::StepFailureReason::final_continuity_residual:
  case flow::StepFailureReason::final_pressure_residual:
    result = {DiagnosticFailureClass::non_convergence,
              s.reason == flow::StepFailureReason::final_momentum_residual
                  ? "time-control.final-momentum-residual"
              : s.reason == flow::StepFailureReason::final_transport_residual
                  ? "time-control.final-transport-residual"
              : s.reason == flow::StepFailureReason::final_continuity_residual
                  ? "time-control.final-continuity-residual"
                  : "time-control.final-pressure-residual",
              result.lowest_failing_rank};
    break;
  case flow::StepFailureReason::invalid_input:
  default:
    result = {DiagnosticFailureClass::invalid_input,
              s.attempt_count == 0U ? "time-control.preflight.state"
                                    : "time-control.invalid-input",
              result.lowest_failing_rank};
    break;
  }
  return result;
}

bool layout_valid(const flow::detail::TimeControlDiagnosticSnapshot &s) {
  const auto &b = s.local_box;
  return b.begin.x >= 0 && b.begin.y >= 0 && b.begin.z >= 0 &&
         b.end.x >= b.begin.x && b.end.y >= b.begin.y &&
         b.end.z >= b.begin.z && s.global_extent.x > 0 &&
         s.global_extent.y > 0 && s.global_extent.z > 0 &&
         s.canonical_owned_faces <= s.local_faces &&
         s.local_cell_layout == flow::detail::render_owned_cells(s.local_box) &&
         s.global_cell_layout ==
             flow::detail::render_global_cells(s.global_extent) &&
         s.local_face_layout ==
             flow::detail::render_owned_faces(s.canonical_owned_faces) &&
         s.global_face_layout ==
             flow::detail::render_global_faces(s.global_faces);
}

bool selected_valid(const DiagnosticRequest &request) {
  if (!std::is_sorted(request.selected_fields.begin(),
                      request.selected_fields.end()) ||
      std::adjacent_find(request.selected_fields.begin(),
                         request.selected_fields.end()) !=
          request.selected_fields.end())
    return false;
  return std::all_of(request.selected_fields.begin(),
                     request.selected_fields.end(), [](auto item) {
                       return std::find(kFields.begin(), kFields.end(), item) !=
                              kFields.end();
                     });
}

bool selected(const DiagnosticRequest &request, std::string_view field) {
  return request.selected_fields.empty() ||
         std::binary_search(request.selected_fields.begin(),
                            request.selected_fields.end(), field);
}

bool supported_level(DiagnosticLevel level) noexcept {
  return level == DiagnosticLevel::summary ||
         level == DiagnosticLevel::invariants ||
         level == DiagnosticLevel::counters ||
         level == DiagnosticLevel::bounded_state_sample;
}

constexpr std::size_t kProjectionCapacity = 4096U;
constexpr std::size_t kWireCapacity = 512U;

template <std::size_t Capacity> struct FixedBytes final {
  std::array<unsigned char, Capacity> storage{};
  std::size_t length{};

  std::size_t size() const noexcept { return length; }
  bool empty() const noexcept { return length == 0U; }
  unsigned char *data() noexcept { return storage.data(); }
  const unsigned char *data() const noexcept { return storage.data(); }
  unsigned char &operator[](std::size_t index) noexcept {
    return storage[index];
  }
  const unsigned char &operator[](std::size_t index) const noexcept {
    return storage[index];
  }
  void push(unsigned char value) {
    if (length == Capacity)
      throw std::length_error("time-control canonical buffer is too large");
    storage[length++] = value;
  }
  void append(const unsigned char *source, std::size_t count) {
    if (count > Capacity - length)
      throw std::length_error("time-control canonical buffer is too large");
    std::copy_n(source, count, storage.data() + length);
    length += count;
  }
  void set_size(std::size_t count) {
    if (count > Capacity)
      throw std::length_error("time-control canonical buffer is too large");
    length = count;
  }
};

struct CanonicalWriter final {
  FixedBytes<kProjectionCapacity> bytes;

  void u8(std::uint8_t value) { bytes.push(value); }
  void u32(std::uint32_t value) {
    for (unsigned shift = 0; shift < 32U; shift += 8U)
      u8(static_cast<std::uint8_t>(value >> shift));
  }
  void u64(std::uint64_t value) {
    for (unsigned shift = 0; shift < 64U; shift += 8U)
      u8(static_cast<std::uint8_t>(value >> shift));
  }
  void fp64(double value) {
    std::uint64_t value_bits{};
    std::memcpy(&value_bits, &value, sizeof(value_bits));
    u64(value_bits);
  }
  void boolean(bool value) { u8(value ? 1U : 0U); }
  void string(std::string_view value) {
    u64(static_cast<std::uint64_t>(value.size()));
    bytes.append(reinterpret_cast<const unsigned char *>(value.data()),
                 value.size());
  }
};

FixedBytes<kProjectionCapacity>
request_projection(const DiagnosticRequest &request) {
  CanonicalWriter out;
  out.u8(static_cast<std::uint8_t>(request.level));
  out.u8(static_cast<std::uint8_t>(request.scope));
  out.u64(request.frame.step);
  out.fp64(request.frame.time_s);
  out.string(request.frame.phase);
  out.u64(static_cast<std::uint64_t>(request.selected_fields.size()));
  for (auto field : request.selected_fields)
    out.string(field);
  out.u64(static_cast<std::uint64_t>(request.sample_budget));
  return out.bytes;
}

void write_state(CanonicalWriter &out, const flow::TimeControlState &state) {
  out.u32(state.schema_version);
  out.u64(state.accepted_step);
  out.fp64(state.proposed_next_dt_s);
  out.fp64(state.last_accepted_dt_s);
  out.u8(static_cast<std::uint8_t>(state.last_accepted_order));
  out.boolean(state.history_ready);
  out.boolean(state.last_all_linear_solves_within_half_limit);
  out.fp64(state.last_convective_rate_per_s);
  out.fp64(state.last_diffusive_rate_per_s);
  out.boolean(state.last_stability_metrics_available);
  out.u32(state.last_retry_count);
  out.u64(state.revision);
  out.u64(state.state_seal);
}

void write_metadata(CanonicalWriter &out,
                    const flow::AcceptedStepMetadata &metadata) {
  out.u64(metadata.step);
  out.fp64(metadata.time_s);
  out.fp64(metadata.dt_s);
  out.fp64(metadata.previous_dt_s);
  out.u8(static_cast<std::uint8_t>(metadata.order));
}

FixedBytes<kProjectionCapacity> provider_projection(
    const flow::detail::TimeControlDiagnosticSnapshot &source,
    const DiagnosticRequest &request) {
  CanonicalWriter out;
  out.u32(kDiagnosticRecordSchemaV1);
  out.u32(static_cast<std::uint32_t>(DiagnosticModuleKind::time_control));
  out.string("hundun.flow.bdf2-retry-controller");
  out.string("primary");
  out.u32(kCapabilities);
  out.string("time-control.advance-result");
  const auto request_bytes = request_projection(request);
  out.u64(static_cast<std::uint64_t>(request_bytes.size()));
  out.bytes.append(request_bytes.data(), request_bytes.size());
  out.u32(static_cast<std::uint32_t>(source.global_extent.x));
  out.u32(static_cast<std::uint32_t>(source.global_extent.y));
  out.u32(static_cast<std::uint32_t>(source.global_extent.z));
  out.string(source.global_cell_layout);
  out.u64(source.global_faces);
  out.string(source.global_face_layout);
  write_state(out, source.state);
  out.u8(static_cast<std::uint8_t>(source.disposition));
  out.u8(static_cast<std::uint8_t>(source.reason));
  out.u32(static_cast<std::uint32_t>(source.lowest_failing_rank));
  out.u8(source.preflight_category);
  out.u64(static_cast<std::uint64_t>(source.attempt_count));
  for (std::size_t i = 0; i < source.attempt_count; ++i) {
    const auto &attempt = source.attempts[i];
    out.fp64(attempt.attempted_dt_s);
    out.u8(static_cast<std::uint8_t>(attempt.order));
    out.u8(static_cast<std::uint8_t>(attempt.disposition));
    out.u8(static_cast<std::uint8_t>(attempt.reason));
    out.u32(static_cast<std::uint32_t>(attempt.lowest_failing_rank));
    out.boolean(attempt.all_linear_solves_within_half_limit);
  }
  out.fp64(source.accepted_dt_s);
  out.fp64(source.proposed_next_dt_s);
  out.fp64(source.convective_rate_per_s);
  out.fp64(source.diffusive_rate_per_s);
  out.boolean(source.stability_metrics_available);
  out.boolean(source.limited_by_min_dt);
  out.u8(static_cast<std::uint8_t>(source.config.mode));
  out.u64(static_cast<std::uint64_t>(source.config.steps));
  out.fp64(source.config.initial_dt_s);
  out.fp64(source.config.min_dt_s);
  out.fp64(source.config.max_dt_s);
  out.fp64(source.config.cfl_target);
  out.fp64(source.config.diffusion_number_target);
  out.fp64(source.config.growth_factor);
  out.fp64(source.config.retry_factor);
  out.u64(static_cast<std::uint64_t>(source.config.max_retries));
  out.u8(static_cast<std::uint8_t>(source.model));
  out.u64(source.controller_identity);
  out.u64(source.report_identity);
  out.u64(source.flow_state_identity);
  out.u64(source.observed_step);
  out.fp64(source.observed_time_s);
  write_metadata(out, source.observed_metadata);
  return out.bytes;
}

bool exact_projection_agrees(const runtime::MpiContext &mpi,
                             const FixedBytes<kProjectionCapacity> &local,
                             std::string_view operation) {
  // Every rank owns the maximum receive capacity before the first raw
  // operation. Only the logical length is learned from rank zero.
  FixedBytes<kProjectionCapacity> root;
  std::uint64_t root_size = static_cast<std::uint64_t>(local.size());
  checked_mpi(
      MPI_Bcast(&root_size, 1, MPI_UINT64_T, 0, mpi.comm()), operation);
  if (root_size > kProjectionCapacity ||
      root_size >
          static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
    throw std::length_error("time-control projection is too large");
  root.set_size(static_cast<std::size_t>(root_size));
  if (mpi.rank() == 0 && root.size() == local.size())
    std::copy_n(local.data(), local.size(), root.data());
  checked_mpi(
      MPI_Bcast(root.data(), static_cast<int>(root.size()), MPI_BYTE, 0,
                mpi.comm()),
      operation);
  bool equal = root.size() == local.size();
  const auto common = std::min(root.size(), local.size());
  for (std::size_t i = 0; i < common; ++i)
    equal = (root[i] == local[i]) && equal;
  return equal;
}

struct CollectivePayload final {
  DiagnosticStateFingerprint fingerprint;
  std::uint64_t eligible_count{};
  std::vector<DiagnosticSample> samples;
};

struct PreparedCollectivePayload final {
  std::array<Tuple, 9> authority{};
  DiagnosticFingerprintParts fingerprint;
  std::uint64_t eligible_count{};
  FixedBytes<kWireCapacity> wire;
};

struct WireValidationError final {};

DiagnosticRecord
build_record(const flow::detail::TimeControlDiagnosticSnapshot &s,
             const DiagnosticRequest &request,
             const CollectivePayload *collective_payload = nullptr) {
  DiagnosticRecord record;
  record.module_kind = DiagnosticModuleKind::time_control;
  record.module_id = "hundun.flow.bdf2-retry-controller";
  record.instance_id = "primary";
  record.level = request.level;
  record.scope = request.scope;
  record.rank = request.frame.rank;
  record.step = request.frame.step;
  record.time_s = describe_fp64(request.frame.time_s);
  record.phase = "time-control.advance-result";
  record.status = s.disposition == flow::TimeAdvanceDisposition::committed
                      ? DiagnosticStatus::ok
                      : DiagnosticStatus::failed;
  record.failure =
      failure_for(s, request.scope == DiagnosticScope::collective);
  const auto local = request.scope == DiagnosticScope::local;
  record.identities = {
      {"flow-state.accepted.cells",
       local ? s.local_cell_layout : s.global_cell_layout,
       s.flow_state_identity, std::nullopt, std::nullopt},
      {"flow-state.accepted.faces",
       local ? s.local_face_layout : s.global_face_layout,
       s.flow_state_identity, std::nullopt, std::nullopt},
      {"time-control.controller", std::nullopt, s.state.revision, std::nullopt,
       std::nullopt},
      {"time-control.report", std::nullopt, s.report_identity, std::nullopt,
       std::nullopt}};

  if (request.level == DiagnosticLevel::summary) {
    const bool attempted = s.attempt_count != 0U;
    const double attempted_dt =
        attempted ? s.attempts[s.attempt_count - 1U].attempted_dt_s : 0.0;
    record.metrics = {
        {"time-control.accepted-dt", DiagnosticMetricKind::state_summary, "s",
         s.disposition == flow::TimeAdvanceDisposition::committed
             ? describe_fp64(s.accepted_dt_s)
             : DiagnosticFp64{}},
        {"time-control.attempted-dt", DiagnosticMetricKind::state_summary, "s",
         attempted ? describe_fp64(attempted_dt) : DiagnosticFp64{}},
        {"time-control.convective-number",
         DiagnosticMetricKind::state_summary, "1",
         attempted && s.stability_metrics_available
             ? describe_fp64(attempted_dt * s.convective_rate_per_s)
             : DiagnosticFp64{}},
        {"time-control.diffusion-number",
         DiagnosticMetricKind::state_summary, "1",
         attempted && s.stability_metrics_available
             ? describe_fp64(attempted_dt * s.diffusive_rate_per_s)
             : DiagnosticFp64{}},
        {"time-control.next-dt", DiagnosticMetricKind::state_summary, "s",
         describe_fp64(s.proposed_next_dt_s)}};
  } else if (request.level == DiagnosticLevel::invariants) {
    const auto invariant = [](std::string id, std::string unit, double observed,
                              double limit, InvariantRelation relation) {
      DiagnosticInvariant result{std::move(id), std::move(unit),
                                 describe_fp64(observed), describe_fp64(limit),
                                 relation, false};
      if (relation == InvariantRelation::finite ||
          relation == InvariantRelation::positive)
        result.limit = {};
      result.passed = evaluate_invariant(result);
      return result;
    };
    const bool state_valid =
        flow::detail::TimeControlStateCodec::semantically_valid(
            s.config, s.model, s.observed_metadata, s.state);
    bool adaptive_limit = s.config.mode == config::TimeMode::fixed ||
                          s.attempt_count == 0U ||
                          !s.stability_metrics_available ||
                          s.limited_by_min_dt;
    if (!adaptive_limit && s.attempt_count != 0U) {
      const double attempted =
          s.attempts[s.attempt_count - 1U].attempted_dt_s;
      const double cfl = attempted * s.convective_rate_per_s;
      const double diffusion = attempted * s.diffusive_rate_per_s;
      adaptive_limit =
          std::isfinite(cfl) && std::isfinite(diffusion) &&
          cfl <= s.config.cfl_target &&
          diffusion <= s.config.diffusion_number_target;
    }
    bool order_history =
        s.state.history_ready == (s.state.accepted_step != 0U) &&
        s.state.accepted_step == s.observed_metadata.step &&
        s.state.last_accepted_order == s.observed_metadata.order;
    if (s.state.accepted_step == 0U)
      order_history =
          order_history &&
          s.observed_metadata.order ==
              flow::MomentumTimeOrder::backward_euler;
    else if (s.state.accepted_step == 1U)
      order_history =
          order_history &&
          s.observed_metadata.order ==
              flow::MomentumTimeOrder::backward_euler;
    else {
      const double ratio =
          s.observed_metadata.dt_s / s.observed_metadata.previous_dt_s;
      const auto expected =
          std::isfinite(ratio) && ratio >= 0.5 && ratio <= 2.0
              ? flow::MomentumTimeOrder::bdf2
              : flow::MomentumTimeOrder::backward_euler;
      order_history = order_history && s.observed_metadata.order == expected;
    }
    record.invariants = {
        invariant("time-control.adaptive-limit-or-minimum", "1",
                  adaptive_limit ? 1.0 : 0.0, 1.0,
                  InvariantRelation::equal),
        invariant("time-control.attempt-count-bounded", "count",
                  static_cast<double>(s.attempt_count), 9.0,
                  InvariantRelation::less_equal),
        invariant("time-control.controller-state-valid", "1",
                  state_valid ? 1.0 : 0.0, 1.0, InvariantRelation::equal),
        invariant("time-control.next-dt-at-least-minimum", "s",
                  s.proposed_next_dt_s, s.config.min_dt_s,
                  InvariantRelation::greater_equal),
        invariant("time-control.next-dt-at-most-maximum", "s",
                  s.proposed_next_dt_s, s.config.max_dt_s,
                  InvariantRelation::less_equal),
        invariant("time-control.next-dt-positive", "s", s.proposed_next_dt_s,
                  0.0, InvariantRelation::positive),
        invariant("time-control.order-history-consistent", "1",
                  order_history ? 1.0 : 0.0, 1.0,
                  InvariantRelation::equal)};
  } else if (request.level == DiagnosticLevel::counters) {
    record.counters = {
        {"time-control.accepted-step", "count", s.state.accepted_step},
        {"time-control.attempt-count", "count", s.attempt_count},
        {"time-control.controller-revision", "count", s.state.revision},
        {"time-control.retry-count", "count",
         s.attempt_count == 0U ? 0U : s.attempt_count - 1U},
        {"time-control.stability-reductions", "count",
         s.stability_metrics_available ? 1U : 0U}};
  }

  DiagnosticFingerprintAccumulator fingerprint;
  const auto authority = tuples(s);
  if (collective_payload == nullptr) {
    if (s.relative_rank == 0) {
      for (const auto &tuple : authority)
        fingerprint.add(tuple.field, 0U, tuple.component,
                        describe_fp64(tuple.value));
    }
    record.state_fingerprint = fingerprint.finish();
  } else {
    record.state_fingerprint = collective_payload->fingerprint;
  }
  if (request.level == DiagnosticLevel::bounded_state_sample) {
    record.sample_budget = request.sample_budget;
    if (collective_payload != nullptr) {
      record.eligible_sample_count = collective_payload->eligible_count;
      record.samples = collective_payload->samples;
    } else if (s.relative_rank == 0) {
      record.eligible_sample_count =
          static_cast<std::uint64_t>(std::count_if(
              authority.begin(), authority.end(), [&](const auto &tuple) {
                return selected(request, tuple.field);
              }));
      for (const auto &tuple : authority) {
        if (!selected(request, tuple.field) ||
            record.samples.size() >= request.sample_budget)
          continue;
        record.samples.push_back(
            {std::string(tuple.field), 0U, tuple.component,
             std::string(tuple.unit), describe_fp64(tuple.value)});
      }
    }
    record.samples_truncated =
        record.samples.size() < record.eligible_sample_count;
  }
  return record;
}

void validate_request(const flow::detail::TimeControlDiagnosticSnapshot &s,
                      const DiagnosticRequest &request,
                      DiagnosticScope expected) {
  if (request.scope != expected || request.frame.rank != s.relative_rank ||
      request.frame.step != s.observed_step ||
      std::memcmp(&request.frame.time_s, &s.observed_time_s,
                  sizeof(double)) != 0 ||
      request.frame.phase != "time-control.advance-result")
    throw DiagnosticCollectionError(DiagnosticFailureClass::invalid_request,
                                    "time-control.diagnostics.frame", -1,
                                    "time-control diagnostic frame mismatch");
  if (!selected_valid(request))
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_request,
        "time-control.diagnostics.selected-field", -1,
        "time-control selected fields are invalid");
}

void require_phase(const runtime::MpiContext &mpi, bool local_ok,
                   DiagnosticFailureClass classification,
                   std::string_view code, std::string_view message) {
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
  ++diagnostic_phase_count;
#endif
  const auto status = runtime::collective_status(mpi, local_ok, code);
  if (!status.ok)
    throw DiagnosticCollectionError(classification, std::string(code),
                                    status.failing_rank,
                                    std::string(message));
}

std::uint32_t field_index(std::string_view field) {
  const auto found = std::find(kFields.begin(), kFields.end(), field);
  if (found == kFields.end())
    throw std::runtime_error("unknown time-control sample field");
  return static_cast<std::uint32_t>(found - kFields.begin());
}

std::uint32_t unit_code(std::string_view unit) {
  if (unit == "count")
    return 0U;
  if (unit == "s")
    return 1U;
  if (unit == "1")
    return 2U;
  throw std::runtime_error("unknown time-control sample unit");
}

std::string_view field_unit(std::uint32_t field) {
  if (field == 0U || field == 3U || field == 6U)
    return "count";
  if (field == 2U || field == 5U)
    return "s";
  return "1";
}

FixedBytes<kWireCapacity> sample_wire(
    const std::array<Tuple, 9> &authority, bool authority_rank,
    const DiagnosticRequest &request) {
  CanonicalWriter wire;
  wire.u32(1U);
  std::array<const Tuple *, 9> retained{};
  std::size_t retained_count{};
  if (authority_rank) {
    for (const auto &tuple : authority) {
      if (selected(request, tuple.field) &&
          retained_count < request.sample_budget)
        retained[retained_count++] = &tuple;
    }
  }
  wire.u64(static_cast<std::uint64_t>(retained_count));
  for (std::size_t index = 0; index < retained_count; ++index) {
    const auto *tuple = retained[index];
    wire.u32(field_index(tuple->field));
    wire.u64(0U);
    wire.u32(tuple->component);
    wire.u32(unit_code(tuple->unit));
    std::uint64_t value_bits{};
    std::memcpy(&value_bits, &tuple->value, sizeof(value_bits));
    wire.u64(value_bits);
  }
  FixedBytes<kWireCapacity> result;
  result.append(wire.bytes.data(), wire.bytes.size());
  return result;
}

std::uint32_t read_u32(const FixedBytes<kWireCapacity> &bytes,
                       std::size_t &offset) {
  if (offset > bytes.size() || bytes.size() - offset < 4U)
    throw std::runtime_error("short time-control sample wire");
  std::uint32_t result{};
  for (unsigned i = 0; i < 4U; ++i)
    result |= static_cast<std::uint32_t>(bytes[offset++]) << (8U * i);
  return result;
}

std::uint64_t read_u64(const FixedBytes<kWireCapacity> &bytes,
                       std::size_t &offset) {
  if (offset > bytes.size() || bytes.size() - offset < 8U)
    throw std::runtime_error("short time-control sample wire");
  std::uint64_t result{};
  for (unsigned i = 0; i < 8U; ++i)
    result |= static_cast<std::uint64_t>(bytes[offset++]) << (8U * i);
  return result;
}

PreparedCollectivePayload prepare_collective_payload(
    const flow::detail::TimeControlDiagnosticSnapshot &snapshot,
    int rank, const DiagnosticRequest &request) {
  PreparedCollectivePayload prepared;
  DiagnosticFingerprintAccumulator local_fingerprint;
  prepared.authority = tuples(snapshot);
  if (rank == 0) {
    for (const auto &tuple : prepared.authority)
      local_fingerprint.add(tuple.field, 0U, tuple.component,
                            describe_fp64(tuple.value));
  }
  prepared.fingerprint = local_fingerprint.parts();
  if (rank == 0 &&
      request.level == DiagnosticLevel::bounded_state_sample) {
    prepared.eligible_count = static_cast<std::uint64_t>(std::count_if(
        prepared.authority.begin(), prepared.authority.end(),
        [&](const auto &tuple) { return selected(request, tuple.field); }));
    prepared.wire = sample_wire(prepared.authority, true, request);
  }
  return prepared;
}

CollectivePayload collective_payload(
    const runtime::MpiContext &mpi, const DiagnosticRequest &request,
    const PreparedCollectivePayload &prepared) {
  auto parts = prepared.fingerprint;
  checked_mpi(
      MPI_Allreduce(MPI_IN_PLACE, &parts.xor64, 1, MPI_UINT64_T, MPI_BXOR,
                    mpi.comm()),
      "MPI_Allreduce(time-control fingerprint xor)");
  checked_mpi(
      MPI_Allreduce(MPI_IN_PLACE, &parts.sum64, 1, MPI_UINT64_T, MPI_SUM,
                    mpi.comm()),
      "MPI_Allreduce(time-control fingerprint sum)");

  std::uint64_t eligible = prepared.eligible_count;
  checked_mpi(
      MPI_Allreduce(MPI_IN_PLACE, &eligible, 1, MPI_UINT64_T, MPI_SUM,
                    mpi.comm()),
      "MPI_Allreduce(time-control eligible samples)");

  // The inline receive capacity was constructed during synchronized
  // preparation on every rank; raw operations only update its logical size.
  FixedBytes<kWireCapacity> wire = prepared.wire;
  std::uint64_t wire_size = static_cast<std::uint64_t>(wire.size());
  checked_mpi(
      MPI_Bcast(&wire_size, 1, MPI_UINT64_T, 0, mpi.comm()),
      "MPI_Bcast(time-control sample wire size)");
  if (wire_size > kWireCapacity ||
      wire_size >
          static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
    throw std::length_error("time-control sample wire is too large");
  if (mpi.rank() != 0)
    wire.set_size(static_cast<std::size_t>(wire_size));
  checked_mpi(
      MPI_Bcast(wire.data(), static_cast<int>(wire.size()), MPI_BYTE, 0,
                mpi.comm()),
      "MPI_Bcast(time-control sample wire)");

  CollectivePayload result;
  try {
    DiagnosticFingerprintAccumulator combined;
    combined.combine(parts);
    result.fingerprint = combined.finish();
  } catch (...) {
    throw;
  }
  result.eligible_count = eligible;
  if (request.level != DiagnosticLevel::bounded_state_sample) {
    if (!wire.empty() || eligible != 0U)
      throw std::runtime_error("unexpected non-sample wire");
    return result;
  }

  try {
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (mpi.rank() == wire_mutation_rank) {
      const auto mutation = wire_mutation;
      if (mutation == test::TimeControlWireMutation::short_payload &&
          !wire.empty())
        wire.set_size(wire.size() - 1U);
      else if (mutation == test::TimeControlWireMutation::trailing_payload)
        wire.push(0U);
      else if (wire.size() >= 68U) {
        std::size_t offset = 0U;
        switch (mutation) {
        case test::TimeControlWireMutation::count:
          offset = 4U;
          break;
        case test::TimeControlWireMutation::field:
          offset = 12U;
          break;
        case test::TimeControlWireMutation::global_id:
          offset = 16U;
          break;
        case test::TimeControlWireMutation::component:
          offset = 24U;
          break;
        case test::TimeControlWireMutation::unit:
          offset = 28U;
          break;
        case test::TimeControlWireMutation::value:
          offset = 32U;
          break;
        case test::TimeControlWireMutation::duplicate:
          std::copy_n(wire.data() + 12U, 28U, wire.data() + 40U);
          offset = wire.size();
          break;
        default:
          offset = wire.size();
          break;
        }
        if (offset < wire.size())
          wire[offset] ^= 0xffU;
      }
    }
#endif
    std::size_t offset{};
    if (read_u32(wire, offset) != 1U)
      throw std::runtime_error("unsupported time-control sample wire");
    const auto count = read_u64(wire, offset);
    std::array<const Tuple *, 9> expected{};
    std::size_t expected_count{};
    for (const auto &tuple : prepared.authority) {
      if (selected(request, tuple.field) &&
          expected_count < request.sample_budget)
        expected[expected_count++] = &tuple;
    }
    if (count != expected_count || count > eligible)
      throw std::runtime_error("invalid time-control sample count");
    std::string previous;
    std::uint32_t previous_component{};
    bool first = true;
    for (std::uint64_t index = 0; index < count; ++index) {
      const auto field = read_u32(wire, offset);
      const auto global_id = read_u64(wire, offset);
      const auto component = read_u32(wire, offset);
      const auto unit = read_u32(wire, offset);
      const auto value_bits = read_u64(wire, offset);
      const auto *expected_tuple = expected[static_cast<std::size_t>(index)];
      std::uint64_t expected_value_bits{};
      std::memcpy(&expected_value_bits, &expected_tuple->value,
                  sizeof(expected_value_bits));
      if (field >= kFields.size() || global_id != 0U ||
          !selected(request, kFields[field]) ||
          kFields[field] != expected_tuple->field ||
          component != expected_tuple->component ||
          unit != unit_code(expected_tuple->unit) ||
          value_bits != expected_value_bits)
        throw std::runtime_error("invalid time-control sample tuple");
      if (unit != unit_code(field_unit(field)))
        throw std::runtime_error("invalid time-control sample unit");
      if ((field == 0U || field == 6U) ? component > 1U : component != 0U)
        throw std::runtime_error("invalid time-control sample component");
      const std::string field_name(kFields[field]);
      if (!first &&
          (field_name < previous ||
           (field_name == previous && component <= previous_component)))
        throw std::runtime_error("noncanonical time-control sample tuple");
      double value{};
      std::memcpy(&value, &value_bits, sizeof(value));
      result.samples.push_back(
          {field_name, 0U, component, std::string(field_unit(field)),
           describe_fp64(value)});
      previous = field_name;
      previous_component = component;
      first = false;
    }
    if (offset != wire.size())
      throw std::runtime_error("trailing time-control sample wire bytes");
  } catch (...) {
    throw WireValidationError{};
  }
  return result;
}

} // namespace

namespace detail {
struct TimeControlAdapter final {
  static const flow::detail::TimeControlDiagnosticSnapshot &
  snapshot(const flow::TimeControlDiagnosticSource &source) {
    if (!source.impl_)
      throw DiagnosticCollectionError(
          DiagnosticFailureClass::invalid_input,
          "time-control.diagnostics.stale-source", -1,
          "time-control diagnostic source was moved from");
    return source.impl_->snapshot;
  }
};
} // namespace detail

DiagnosticDescriptor describe_time_control_diagnostics() noexcept {
  return {1U, DiagnosticModuleKind::time_control,
          "hundun.flow.bdf2-retry-controller", "primary", kCapabilities};
}
DiagnosticDescriptor
describe_diagnostics(const flow::TimeControlDiagnosticSource &) noexcept {
  return describe_time_control_diagnostics();
}
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const flow::TimeControlDiagnosticSource &) {
  return {kFields.begin(), kFields.end()};
}

void collect_diagnostics(const flow::TimeControlDiagnosticSource &source,
                         const DiagnosticRequest &request,
                         DiagnosticSink &sink) {
  const auto &snapshot = detail::TimeControlAdapter::snapshot(source);
  if (!layout_valid(snapshot))
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::layout,
        "time-control.diagnostics.local-layout", -1,
        "time-control diagnostic layout is invalid");
  if (!supported_level(request.level) ||
      request.scope != DiagnosticScope::local)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::capability,
        "time-control.diagnostics.capability", -1,
        "time-control diagnostic capability is unsupported");
  validate_request(snapshot, request, DiagnosticScope::local);
  auto record = build_record(snapshot, request);
  validate(record, describe_time_control_diagnostics(), request);
  try {
    sink.submit(record);
  } catch (...) {
    throw DiagnosticCollectionError(DiagnosticFailureClass::sink_failure,
                                    "diagnostics.sink.submit", -1,
                                    "diagnostic sink rejected record");
  }
}

void collect_diagnostics(const flow::TimeControlDiagnosticSource &source,
                         const runtime::MpiContext &mpi,
                         const DiagnosticRequest &request,
                         DiagnosticSink &sink) {
  const flow::detail::TimeControlDiagnosticSnapshot *snapshot{};
  bool phase_ok = true;
  try {
    snapshot = &detail::TimeControlAdapter::snapshot(source);
    phase_ok = layout_valid(*snapshot) &&
               !fault_here(FaultPoint::phase1_layout, mpi.rank());
  } catch (const runtime::MpiOperationError &) {
    throw;
  } catch (...) {
    phase_ok = false;
  }
  require_phase(mpi, phase_ok, DiagnosticFailureClass::layout,
                "time-control.diagnostics.local-layout",
                "time-control collective diagnostic layout is invalid");

  phase_ok = supported_level(request.level) &&
             request.scope == DiagnosticScope::collective;
  require_phase(mpi, phase_ok, DiagnosticFailureClass::capability,
                "time-control.diagnostics.capability",
                "time-control diagnostic capability is unsupported");

  FixedBytes<kProjectionCapacity> request_bytes;
  phase_ok = true;
  try {
    validate_request(*snapshot, request, DiagnosticScope::collective);
    validate(describe_time_control_diagnostics());
    validate(request, describe_time_control_diagnostics());
    request_bytes = request_projection(request);
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (mpi.rank() == request_mutation_rank && !request_bytes.empty() &&
        request_mutation_offset < request_bytes.size())
      request_bytes[request_mutation_offset] ^= 1U;
#endif
    if (fault_here(FaultPoint::phase2_request, mpi.rank()))
      throw std::runtime_error("injected request preparation failure");
    if (fault_here(FaultPoint::projection_root_size, mpi.rank()) ||
        fault_here(FaultPoint::projection_payload, mpi.rank()))
      throw std::runtime_error("injected projection preparation failure");
  } catch (const runtime::MpiOperationError &) {
    throw;
  } catch (...) {
    phase_ok = false;
  }
  require_phase(mpi, phase_ok, DiagnosticFailureClass::invalid_request,
                "time-control.diagnostics.request-preparation",
                "time-control diagnostic request preparation failed");
  const bool request_agrees = exact_projection_agrees(
      mpi, request_bytes, "MPI_Bcast(time-control diagnostic request)");
  require_phase(mpi, request_agrees,
                DiagnosticFailureClass::invalid_request,
                "time-control.diagnostics.request-agreement",
                "time-control collective diagnostic requests disagree");

  FixedBytes<kProjectionCapacity> provider_bytes;
  phase_ok = true;
  try {
    provider_bytes = provider_projection(*snapshot, request);
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    if (mpi.rank() == provider_mutation_rank && !provider_bytes.empty() &&
        provider_mutation_offset < provider_bytes.size())
      provider_bytes[provider_mutation_offset] ^= 1U;
#endif
    if (fault_here(FaultPoint::phase3_provider, mpi.rank()))
      throw std::runtime_error("injected provider preparation failure");
  } catch (const runtime::MpiOperationError &) {
    throw;
  } catch (...) {
    phase_ok = false;
  }
  require_phase(mpi, phase_ok, DiagnosticFailureClass::invalid_input,
                "time-control.diagnostics.provider-agreement",
                "time-control diagnostic provider preparation failed");
  const bool provider_agrees = exact_projection_agrees(
      mpi, provider_bytes, "MPI_Bcast(time-control diagnostic provider)");
  require_phase(mpi, provider_agrees, DiagnosticFailureClass::invalid_input,
                "time-control.diagnostics.provider-agreement",
                "time-control collective diagnostic providers disagree");

  CollectivePayload payload;
  PreparedCollectivePayload prepared_payload;
  phase_ok = true;
  try {
    prepared_payload =
        prepare_collective_payload(*snapshot, mpi.rank(), request);
    if (fault_here(FaultPoint::phase4_payload, mpi.rank()))
      throw std::runtime_error("injected payload preparation failure");
    if (fault_here(FaultPoint::wire_root_size, mpi.rank()) ||
        fault_here(FaultPoint::wire_payload, mpi.rank()))
      throw std::runtime_error("injected wire preparation failure");
  } catch (...) {
    phase_ok = false;
  }
  require_phase(mpi, phase_ok, DiagnosticFailureClass::layout,
                "time-control.diagnostics.sample-preparation",
                "time-control sample preparation failed");
  bool wire_ok = true;
  try {
    payload = collective_payload(mpi, request, prepared_payload);
    phase_ok =
        !fault_here(FaultPoint::fingerprint_aggregation, mpi.rank()) &&
        !fault_here(FaultPoint::eligible_aggregation, mpi.rank());
  } catch (const WireValidationError &) {
    wire_ok = false;
    phase_ok = true;
  } catch (const runtime::MpiOperationError &) {
    throw;
  } catch (...) {
    phase_ok = false;
  }
  require_phase(mpi, phase_ok, DiagnosticFailureClass::layout,
                "time-control.diagnostics.aggregation",
                "time-control diagnostic aggregation failed");

  wire_ok = wire_ok && !fault_here(FaultPoint::phase4_wire, mpi.rank());
  require_phase(mpi, wire_ok, DiagnosticFailureClass::layout,
                "time-control.diagnostics.sample-wire",
                "time-control sample wire is invalid");

  DiagnosticRecord record;
  phase_ok = true;
  try {
    record = build_record(*snapshot, request, &payload);
    validate(record, describe_time_control_diagnostics(), request);
    if (fault_here(FaultPoint::phase5_record, mpi.rank()))
      throw std::runtime_error("injected record validation failure");
  } catch (const runtime::MpiOperationError &) {
    throw;
  } catch (...) {
    phase_ok = false;
  }
  require_phase(mpi, phase_ok, DiagnosticFailureClass::invalid_input,
                "time-control.diagnostics.record",
                "time-control collective diagnostic record is invalid");

  bool sink_ok = true;
  try {
#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
    ++diagnostic_submission_count;
#endif
    if (fault_here(FaultPoint::phase6_sink, mpi.rank()))
      throw std::runtime_error("injected sink failure");
    sink.submit(record);
  } catch (...) {
    sink_ok = false;
  }
  require_phase(mpi, sink_ok, DiagnosticFailureClass::sink_failure,
                "diagnostics.sink.submit",
                "time-control collective diagnostic sink rejected record");
}

#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
void test::TimeControlDiagnosticsTestAccess::reset() noexcept {
  diagnostic_fault = test::TimeControlDiagnosticFault::none;
  diagnostic_fault_rank = -1;
  diagnostic_phase_count = 0U;
  diagnostic_raw_count = 0U;
  diagnostic_submission_count = 0U;
  diagnostic_raw_fault_ordinal = 0U;
  request_mutation_offset = std::numeric_limits<std::size_t>::max();
  request_mutation_rank = -1;
  provider_mutation_offset = std::numeric_limits<std::size_t>::max();
  provider_mutation_rank = -1;
  wire_mutation = test::TimeControlWireMutation::none;
  wire_mutation_rank = -1;
}
void test::TimeControlDiagnosticsTestAccess::set_raw_fault(
    std::size_t ordinal) noexcept {
  diagnostic_raw_fault_ordinal = ordinal;
}
void test::TimeControlDiagnosticsTestAccess::set_request_projection_mutation(
    std::size_t offset, int rank) noexcept {
  request_mutation_offset = offset;
  request_mutation_rank = rank;
}
void test::TimeControlDiagnosticsTestAccess::set_provider_projection_mutation(
    std::size_t offset, int rank) noexcept {
  provider_mutation_offset = offset;
  provider_mutation_rank = rank;
}
void test::TimeControlDiagnosticsTestAccess::set_wire_mutation(
    test::TimeControlWireMutation mutation, int rank) noexcept {
  wire_mutation = mutation;
  wire_mutation_rank = rank;
}
void test::TimeControlDiagnosticsTestAccess::set_fault(
    test::TimeControlDiagnosticFault fault, int rank) noexcept {
  diagnostic_fault = fault;
  diagnostic_fault_rank = rank;
}
std::size_t test::TimeControlDiagnosticsTestAccess::phase_count() noexcept {
  return diagnostic_phase_count;
}
std::size_t
test::TimeControlDiagnosticsTestAccess::raw_operation_count() noexcept {
  return diagnostic_raw_count;
}
std::size_t
test::TimeControlDiagnosticsTestAccess::submission_count() noexcept {
  return diagnostic_submission_count;
}
void test::TimeControlDiagnosticsTestAccess::set_state_counters(
    flow::TimeControlDiagnosticSource &source, std::uint64_t accepted_step,
    std::uint64_t revision) noexcept {
  if (!source.impl_)
    return;
  source.impl_->snapshot.state.accepted_step = accepted_step;
  source.impl_->snapshot.state.revision = revision;
}
void test::TimeControlDiagnosticsTestAccess::set_outcome(
    flow::TimeControlDiagnosticSource &source,
    flow::TimeAdvanceDisposition disposition, flow::StepFailureReason reason,
    int lowest_failing_rank, std::uint8_t preflight_category,
    std::size_t attempt_count) noexcept {
  if (!source.impl_)
    return;
  auto &snapshot = source.impl_->snapshot;
  snapshot.disposition = disposition;
  snapshot.reason = reason;
  snapshot.lowest_failing_rank = lowest_failing_rank;
  snapshot.preflight_category = preflight_category;
  snapshot.attempt_count =
      std::min(attempt_count, snapshot.attempts.size());
}
void test::TimeControlDiagnosticsTestAccess::set_local_mutation(
    flow::TimeControlDiagnosticSource &source,
    test::TimeControlLocalMutation mutation) {
  if (!source.impl_)
    return;
  auto &snapshot = source.impl_->snapshot;
  switch (mutation) {
  case test::TimeControlLocalMutation::local_box:
    snapshot.local_box.begin.x = -1;
    break;
  case test::TimeControlLocalMutation::canonical_owned_faces:
    snapshot.canonical_owned_faces = snapshot.local_faces + 1U;
    break;
  case test::TimeControlLocalMutation::local_faces:
    snapshot.local_faces = 0U;
    break;
  case test::TimeControlLocalMutation::local_cell_layout:
    snapshot.local_cell_layout += ".mutated";
    break;
  case test::TimeControlLocalMutation::global_cell_layout:
    snapshot.global_cell_layout += ".mutated";
    break;
  case test::TimeControlLocalMutation::local_face_layout:
    snapshot.local_face_layout += ".mutated";
    break;
  case test::TimeControlLocalMutation::global_face_layout:
    snapshot.global_face_layout += ".mutated";
    break;
  }
}
#endif

} // namespace hundun::diagnostics

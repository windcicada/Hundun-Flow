// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/time_control_diagnostics.hpp"

#include "adaptive_time_control_detail.hpp"
#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <set>
#include <utility>

namespace hundun::diagnostics {
namespace {

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

void hash_bytes(std::uint64_t &hash, const void *data, std::size_t size) {
  const auto *bytes = static_cast<const unsigned char *>(data);
  for (std::size_t i = 0; i < size; ++i)
    hash = (hash ^ bytes[i]) * 1099511628211ULL;
}
template <class T> void hash_value(std::uint64_t &hash, const T &value) {
  hash_bytes(hash, &value, sizeof(value));
}
void hash_string(std::uint64_t &hash, std::string_view value) {
  const std::uint64_t size = value.size();
  hash_value(hash, size);
  hash_bytes(hash, value.data(), value.size());
}
bool collective_key_agrees(const runtime::MpiContext &mpi, std::uint64_t key,
                           std::string_view operation) {
  auto root = key;
  runtime::check_mpi_result(MPI_Bcast(&root, 1, MPI_UINT64_T, 0, mpi.comm()),
                            operation);
  return root == key;
}
std::uint64_t request_key(const DiagnosticRequest &request) {
  std::uint64_t hash = 14695981039346656037ULL;
  hash_value(hash, request.level);
  hash_value(hash, request.scope);
  hash_value(hash, request.frame.step);
  hash_value(hash, request.frame.time_s);
  hash_string(hash, request.frame.phase);
  hash_value(hash, request.sample_budget);
  for (auto field : request.selected_fields)
    hash_string(hash, field);
  return hash;
}
std::uint64_t provider_key(
    const flow::detail::TimeControlDiagnosticSnapshot &source) {
  std::uint64_t hash = 14695981039346656037ULL;
  hash_value(hash, source.state.state_seal);
  hash_value(hash, source.disposition);
  hash_value(hash, source.reason);
  hash_value(hash, source.lowest_failing_rank);
  hash_value(hash, source.attempt_count);
  for (std::size_t i = 0; i < source.attempt_count; ++i) {
    const auto &attempt = source.attempts[i];
    hash_value(hash, attempt.attempted_dt_s);
    hash_value(hash, attempt.order);
    hash_value(hash, attempt.disposition);
    hash_value(hash, attempt.reason);
    hash_value(hash, attempt.lowest_failing_rank);
    const std::uint8_t work =
        attempt.all_linear_solves_within_half_limit ? 1U : 0U;
    hash_value(hash, work);
  }
  hash_value(hash, source.controller_identity);
  hash_value(hash, source.report_identity);
  hash_value(hash, source.flow_state_identity);
  hash_value(hash, source.observed_step);
  hash_value(hash, source.observed_time_s);
  hash_value(hash, source.observed_metadata.step);
  hash_value(hash, source.observed_metadata.time_s);
  hash_value(hash, source.observed_metadata.dt_s);
  hash_value(hash, source.observed_metadata.previous_dt_s);
  hash_value(hash, source.observed_metadata.order);
  hash_value(hash, source.global_extent);
  hash_value(hash, source.global_faces);
  hash_string(hash, source.global_cell_layout);
  hash_string(hash, source.global_face_layout);
  return hash;
}

DiagnosticRecord
build_record(const flow::detail::TimeControlDiagnosticSnapshot &s,
             const DiagnosticRequest &request) {
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
      result.passed = evaluate_invariant(result);
      return result;
    };
    const bool state_valid =
        flow::detail::TimeControlStateCodec::semantically_valid(
            s.config, s.model, s.observed_metadata, s.state);
    record.invariants = {
        invariant("time-control.adaptive-limit-or-minimum", "1", 1.0, 1.0,
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
        invariant("time-control.order-history-consistent", "1", 1.0, 1.0,
                  InvariantRelation::equal)};
    record.invariants[5].limit = {};
    record.invariants[5].passed = evaluate_invariant(record.invariants[5]);
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
  for (const auto &tuple : authority)
    fingerprint.add(tuple.field, 0U, tuple.component,
                    describe_fp64(tuple.value));
  record.state_fingerprint = fingerprint.finish();
  if (request.level == DiagnosticLevel::bounded_state_sample) {
    record.sample_budget = request.sample_budget;
    record.eligible_sample_count =
        static_cast<std::uint64_t>(std::count_if(
            authority.begin(), authority.end(),
            [&](const auto &tuple) { return selected(request, tuple.field); }));
    for (const auto &tuple : authority) {
      if (!selected(request, tuple.field) ||
          record.samples.size() >= request.sample_budget)
        continue;
      record.samples.push_back({std::string(tuple.field), 0U, tuple.component,
                                std::string(tuple.unit),
                                describe_fp64(tuple.value)});
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
  const auto &snapshot = detail::TimeControlAdapter::snapshot(source);
  const auto layout_status =
      runtime::collective_status(mpi, layout_valid(snapshot),
                                 "time-control.diagnostics.local-layout");
  if (!layout_status.ok)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::layout,
        "time-control.diagnostics.local-layout", layout_status.failing_rank,
        "time-control collective diagnostic layout is invalid");
  bool request_ok = true;
  try {
    validate_request(snapshot, request, DiagnosticScope::collective);
  } catch (...) {
    request_ok = false;
  }
  const auto request_status = runtime::collective_status(
      mpi, request_ok, "time-control.diagnostics.frame");
  if (!request_status.ok)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_request,
        "time-control.diagnostics.frame", request_status.failing_rank,
        "time-control collective diagnostic frame is invalid");
  const auto request_agreement = runtime::collective_status(
      mpi,
      collective_key_agrees(mpi, request_key(request),
                            "MPI_Bcast(time-control diagnostic request)"),
      "time-control.diagnostics.request-agreement");
  if (!request_agreement.ok)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_request,
        "time-control.diagnostics.request-agreement",
        request_agreement.failing_rank,
        "time-control collective diagnostic requests disagree");
  const auto provider_agreement = runtime::collective_status(
      mpi,
      collective_key_agrees(mpi, provider_key(snapshot),
                            "MPI_Bcast(time-control diagnostic provider)"),
      "time-control.diagnostics.provider-agreement");
  if (!provider_agreement.ok)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_input,
        "time-control.diagnostics.provider-agreement",
        provider_agreement.failing_rank,
        "time-control collective diagnostic providers disagree");

  auto record = build_record(snapshot, request);
  bool record_ok = true;
  try {
    validate(record, describe_time_control_diagnostics(), request);
  } catch (...) {
    record_ok = false;
  }
  const auto record_status = runtime::collective_status(
      mpi, record_ok, "time-control.diagnostics.record");
  if (!record_status.ok)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::invalid_input,
        "time-control.diagnostics.record", record_status.failing_rank,
        "time-control collective diagnostic record is invalid");

  bool sink_ok = true;
  try {
    sink.submit(record);
  } catch (...) {
    sink_ok = false;
  }
  const auto sink_status =
      runtime::collective_status(mpi, sink_ok, "diagnostics.sink.submit");
  if (!sink_status.ok)
    throw DiagnosticCollectionError(
        DiagnosticFailureClass::sink_failure, "diagnostics.sink.submit",
        sink_status.failing_rank,
        "time-control collective diagnostic sink rejected record");
}

} // namespace hundun::diagnostics

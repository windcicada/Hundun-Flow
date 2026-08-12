// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_les_wale.hpp"

#include <cstdint>
#include <string>

namespace hundun::diagnostics {
namespace {

constexpr std::string_view kModuleId = "hundun.les.wale";
constexpr std::string_view kInstanceId = "primary";
constexpr DiagnosticCapabilityFlags kCapabilities =
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::summary) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::counters);

void add_u64(DiagnosticFingerprintAccumulator &accumulator,
             std::string_view field_id, std::uint64_t value) {
  accumulator.add(field_id, 0U, 0U,
                  describe_fp64(static_cast<double>(
                      static_cast<std::uint32_t>(value))));
  accumulator.add(field_id, 0U, 1U,
                  describe_fp64(static_cast<double>(
                      static_cast<std::uint32_t>(value >> 32U))));
}

DiagnosticStateFingerprint fingerprint(const les::WaleSummary &summary) {
  DiagnosticFingerprintAccumulator accumulator;
  add_u64(accumulator, "identity", summary.identity.value);
  add_u64(accumulator, "nu-t.exact-zero-count", summary.exact_zero_count);
  accumulator.add("nu-t.l2", 0U, 0U,
                  describe_fp64(summary.l2_nu_t_m2_per_s));
  accumulator.add("nu-t.maximum", 0U, 0U,
                  describe_fp64(summary.maximum_nu_t_m2_per_s));
  accumulator.add("nu-t.minimum", 0U, 0U,
                  describe_fp64(summary.minimum_nu_t_m2_per_s));
  add_u64(accumulator, "owned-active-count", summary.owned_active_count);
  return accumulator.finish();
}

DiagnosticRecord build_record(const les::WaleSummary &summary,
                              const DiagnosticRequest &request) {
  DiagnosticRecord record;
  record.schema_version = kDiagnosticRecordSchemaV1;
  record.module_kind = DiagnosticModuleKind::les;
  record.module_id = std::string(kModuleId);
  record.instance_id = std::string(kInstanceId);
  record.level = request.level;
  record.scope = request.scope;
  record.rank = request.frame.rank;
  record.step = request.frame.step;
  record.time_s = describe_fp64(request.frame.time_s);
  record.phase = std::string(request.frame.phase);
  record.state_fingerprint = fingerprint(summary);
  if (request.level == DiagnosticLevel::summary) {
    record.metrics = {
        {"nu-t.l2", DiagnosticMetricKind::state_summary, "m2/s",
         describe_fp64(summary.l2_nu_t_m2_per_s)},
        {"nu-t.maximum", DiagnosticMetricKind::state_summary, "m2/s",
         describe_fp64(summary.maximum_nu_t_m2_per_s)},
        {"nu-t.minimum", DiagnosticMetricKind::state_summary, "m2/s",
         describe_fp64(summary.minimum_nu_t_m2_per_s)},
    };
  } else {
    record.counters = {
        {"identity", "count", summary.identity.value},
        {"nu-t.exact-zero-count", "count", summary.exact_zero_count},
        {"owned-active-count", "count", summary.owned_active_count},
    };
  }
  return record;
}

} // namespace

DiagnosticDescriptor describe_diagnostics(const les::WaleSummary &) noexcept {
  return {kDiagnosticRecordSchemaV1, DiagnosticModuleKind::les, kModuleId,
          kInstanceId, kCapabilities};
}

std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const les::WaleSummary &) {
  return {"identity", "nu-t.exact-zero-count", "nu-t.l2",
          "nu-t.maximum", "nu-t.minimum", "owned-active-count"};
}

void collect_diagnostics(const les::WaleSummary &summary,
                         const DiagnosticRequest &request,
                         DiagnosticSink &sink) {
  const auto descriptor = describe_diagnostics(summary);
  validate(request, descriptor);
  const auto record = build_record(summary, request);
  validate(record, descriptor, request);
  try {
    sink.submit(record);
  } catch (...) {
    throw DiagnosticCollectionError(DiagnosticFailureClass::sink_failure,
                                    "diagnostics.sink.submit", -1,
                                    "Diagnostic sink rejected a record");
  }
}

} // namespace hundun::diagnostics

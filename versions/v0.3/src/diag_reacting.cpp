// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/diag_reacting.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>

namespace hundun::diagnostics {
namespace {

constexpr std::string_view kModuleId = "hundun.flow.reacting";
constexpr std::string_view kInstanceId = "primary";
constexpr DiagnosticCapabilityFlags kCapabilities =
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::summary) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::invariants) |
    static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::counters);

bool valid(const ReactingDiagnosticsSnapshot &value) noexcept {
  if (value.revision == 0U || value.composition_fingerprint == 0U ||
      value.mechanism_sha256.size() != 64U ||
      !std::isfinite(value.final_residual) ||
      !std::isfinite(value.eos_relative_drift) ||
      !std::isfinite(value.enthalpy_budget_j) ||
      value.element_budget_kg.empty() || value.species_budget_kg.empty() ||
      (value.accepted_step &&
       (value.chemistry_calls != 2U ||
        value.pressure_corrector_calls != 2U ||
        value.failure_class != DiagnosticFailureClass::none)) ||
      (!value.accepted_step &&
       (value.failure_class == DiagnosticFailureClass::none ||
        value.failure_code == "none")))
    return false;
  const auto finite = [](double item) { return std::isfinite(item); };
  return std::all_of(value.element_budget_kg.begin(),
                     value.element_budget_kg.end(), finite) &&
         std::all_of(value.species_budget_kg.begin(),
                     value.species_budget_kg.end(), finite);
}

std::uint64_t hash(const ReactingDiagnosticsSnapshot &value) noexcept {
  std::uint64_t result = UINT64_C(14695981039346656037);
  const auto mix = [&](std::uint64_t item) {
    for (unsigned shift = 0; shift < 64U; shift += 8U) {
      result ^= (item >> shift) & 0xffU;
      result *= UINT64_C(1099511628211);
    }
  };
  mix(value.revision);
  mix(value.composition_fingerprint);
  mix(value.chemistry_calls);
  mix(value.chemistry_failures);
  mix(value.pressure_corrector_calls);
  mix(value.workspace_count);
  mix(value.rollback_count);
  for (char character : value.mechanism_sha256)
    mix(static_cast<unsigned char>(character));
  return result;
}

DiagnosticStateFingerprint fingerprint(
    const ReactingDiagnosticsSnapshot &value) {
  DiagnosticFingerprintAccumulator accumulator;
  accumulator.add("revision", 0U, 0U,
                  describe_fp64(static_cast<double>(value.revision)));
  accumulator.add("final-residual", 0U, 0U,
                  describe_fp64(value.final_residual));
  accumulator.add("eos-relative-drift", 0U, 0U,
                  describe_fp64(value.eos_relative_drift));
  accumulator.add("enthalpy-budget", 0U, 0U,
                  describe_fp64(value.enthalpy_budget_j));
  for (std::size_t index = 0; index < value.element_budget_kg.size(); ++index)
    accumulator.add("element-budget", index, 0U,
                    describe_fp64(value.element_budget_kg[index]));
  for (std::size_t index = 0; index < value.species_budget_kg.size(); ++index)
    accumulator.add("species-budget", index, 0U,
                    describe_fp64(value.species_budget_kg[index]));
  return accumulator.finish();
}

DiagnosticRecord record(const ReactingDiagnosticsSnapshot &value,
                        const DiagnosticRequest &request) {
  DiagnosticRecord result;
  result.schema_version = kDiagnosticRecordSchemaV1;
  result.module_kind = DiagnosticModuleKind::reacting;
  result.module_id = std::string(kModuleId);
  result.instance_id = std::string(kInstanceId);
  result.level = request.level;
  result.scope = request.scope;
  result.rank = request.frame.rank;
  result.step = request.frame.step;
  result.time_s = describe_fp64(request.frame.time_s);
  result.phase = std::string(request.frame.phase);
  result.status = value.accepted_step ? DiagnosticStatus::ok
                                      : DiagnosticStatus::failed;
  result.failure = value.accepted_step
                       ? DiagnosticFailure{}
                       : DiagnosticFailure{value.failure_class,
                                           value.failure_code, -1};
  result.state_fingerprint = fingerprint(value);
  result.identities = {
      {"composition", std::nullopt, value.composition_fingerprint,
       std::nullopt, std::nullopt},
      {"mechanism", value.mechanism_sha256, std::nullopt, std::nullopt,
       std::nullopt}};
  if (request.level == DiagnosticLevel::counters) {
    result.counters = {
        {"chemistry.calls", "count", value.chemistry_calls},
        {"chemistry.failures", "count", value.chemistry_failures},
        {"piso.calls", "count", value.pressure_corrector_calls},
        {"rollbacks", "count", value.rollback_count},
        {"workspaces", "count", value.workspace_count}};
  } else if (request.level == DiagnosticLevel::invariants) {
    result.invariants = {
        {"eos-drift.finite", "1", describe_fp64(value.eos_relative_drift),
         describe_fp64(0.0), InvariantRelation::finite,
         std::isfinite(value.eos_relative_drift)},
        {"final-residual.finite", "1", describe_fp64(value.final_residual),
         describe_fp64(0.0), InvariantRelation::finite,
         std::isfinite(value.final_residual)}};
  } else {
    result.metrics.clear();
    for (std::size_t index = 0; index < value.element_budget_kg.size(); ++index)
      result.metrics.push_back(
          {"element-budget." + std::to_string(index),
           DiagnosticMetricKind::conservation, "kg",
           describe_fp64(value.element_budget_kg[index])});
    result.metrics.push_back(
        {"enthalpy-budget", DiagnosticMetricKind::conservation, "J",
         describe_fp64(value.enthalpy_budget_j)});
    result.metrics.push_back(
        {"eos-relative-drift", DiagnosticMetricKind::residual, "1",
         describe_fp64(value.eos_relative_drift)});
    result.metrics.push_back(
        {"final-residual", DiagnosticMetricKind::residual, "1",
         describe_fp64(value.final_residual)});
    for (std::size_t index = 0; index < value.species_budget_kg.size(); ++index)
      result.metrics.push_back(
          {"species-budget." + std::to_string(index),
           DiagnosticMetricKind::conservation, "kg",
           describe_fp64(value.species_budget_kg[index])});
  }
  return result;
}

bool provider_validate(const void *source, std::string &message) noexcept {
  const auto *snapshot =
      static_cast<const ReactingDiagnosticsSnapshot *>(source);
  if (snapshot != nullptr && valid(*snapshot))
    return true;
  message = "reacting diagnostics snapshot is invalid";
  return false;
}
std::uint64_t provider_fingerprint(const void *source) noexcept {
  return hash(*static_cast<const ReactingDiagnosticsSnapshot *>(source));
}
void provider_collect(const void *source, const DiagnosticRequest &request,
                      DiagnosticSink &sink) {
  collect_diagnostics(
      *static_cast<const ReactingDiagnosticsSnapshot *>(source), request,
      sink);
}

} // namespace

DiagnosticDescriptor
describe_diagnostics(const ReactingDiagnosticsSnapshot &) noexcept {
  return {kDiagnosticRecordSchemaV1, DiagnosticModuleKind::reacting, kModuleId,
          kInstanceId, kCapabilities};
}

std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const ReactingDiagnosticsSnapshot &) {
  return {"element-budget", "enthalpy-budget", "eos-relative-drift",
          "final-residual", "revision", "species-budget"};
}

void collect_diagnostics(const ReactingDiagnosticsSnapshot &value,
                         const DiagnosticRequest &request,
                         DiagnosticSink &sink) {
  if (!valid(value))
    throw std::invalid_argument("reacting diagnostics snapshot is invalid");
  const auto descriptor = describe_diagnostics(value);
  validate(request, descriptor);
  const auto output = record(value, request);
  validate(output, descriptor, request);
  sink.submit(output);
}

DiagnosticProvider
make_reacting_diagnostic_provider(const ReactingDiagnosticsSnapshot &value) {
  return {DiagnosticModuleKind::reacting,
          "hundun.flow.reacting.primary",
          &value,
          &provider_validate,
          &provider_fingerprint,
          &provider_collect};
}

} // namespace hundun::diagnostics

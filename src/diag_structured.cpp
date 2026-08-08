// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_structured.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <tuple>

namespace hundun::diagnostics {
namespace {

[[noreturn]] void fail(DiagnosticFailureClass classification, std::string code,
                       std::string message) {
  throw DiagnosticCollectionError(classification, std::move(code), -1,
                                  std::move(message));
}

bool valid_utf8(std::string_view value) noexcept {
  std::size_t index = 0;
  while (index < value.size()) {
    const auto first = static_cast<unsigned char>(value[index]);
    std::size_t continuation = 0;
    std::uint32_t code = 0;
    if (first <= 0x7fU) {
      ++index;
      continue;
    }
    if ((first & 0xe0U) == 0xc0U) {
      continuation = 1;
      code = first & 0x1fU;
      if (code < 2U)
        return false;
    } else if ((first & 0xf0U) == 0xe0U) {
      continuation = 2;
      code = first & 0x0fU;
    } else if ((first & 0xf8U) == 0xf0U) {
      continuation = 3;
      code = first & 0x07U;
      if (code > 4U)
        return false;
    } else {
      return false;
    }
    if (index + continuation >= value.size())
      return false;
    for (std::size_t part = 1; part <= continuation; ++part) {
      const auto byte = static_cast<unsigned char>(value[index + part]);
      if ((byte & 0xc0U) != 0x80U)
        return false;
      code = (code << 6U) | (byte & 0x3fU);
    }
    if ((continuation == 2U && code < 0x800U) ||
        (continuation == 3U && code < 0x10000U) ||
        (code >= 0xd800U && code <= 0xdfffU) || code > 0x10ffffU)
      return false;
    index += continuation + 1U;
  }
  return true;
}

bool valid_id(std::string_view value) noexcept {
  if (value.empty() || value.size() > 128U || !valid_utf8(value))
    return false;
  return std::all_of(value.begin(), value.end(), [](char character) {
    return (character >= 'a' && character <= 'z') ||
           (character >= '0' && character <= '9') || character == '_' ||
           character == '.' || character == '-';
  });
}

void require_id(std::string_view value, const char *kind,
                DiagnosticFailureClass classification) {
  if (!valid_id(value))
    fail(classification, std::string("diagnostics.") + kind + ".invalid",
         std::string("invalid diagnostic ") + kind);
}

bool valid_unit(std::string_view unit) noexcept {
  constexpr std::array<std::string_view, 19> units{
      "1",    "s",     "m",      "m2",     "m3",   "m/s", "m2/s",
      "kg",   "kg/m3", "kg/s",   "kg*m/s", "Pa",   "J",   "J/m3",
      "J/kg", "K",     "degree", "byte",   "count"};
  return std::find(units.begin(), units.end(), unit) != units.end();
}

void require_unit(std::string_view unit) {
  if (!valid_unit(unit))
    fail(DiagnosticFailureClass::invalid_input, "diagnostics.unit.invalid",
         "invalid diagnostic unit");
}

DiagnosticValueStatus classify_bits(std::uint64_t bits) noexcept {
  const std::uint64_t exponent = (bits >> 52U) & UINT64_C(0x7ff);
  const std::uint64_t fraction = bits & UINT64_C(0x000fffffffffffff);
  if (exponent != UINT64_C(0x7ff))
    return DiagnosticValueStatus::finite;
  if (fraction == 0U)
    return (bits >> 63U) == 0U ? DiagnosticValueStatus::positive_infinity
                               : DiagnosticValueStatus::negative_infinity;
  return (fraction & UINT64_C(0x0008000000000000)) != 0U
             ? DiagnosticValueStatus::quiet_nan
             : DiagnosticValueStatus::signaling_nan;
}

void validate_fp64(DiagnosticFp64 value) {
  if (value.status == DiagnosticValueStatus::unavailable) {
    if (value.bits != 0U)
      fail(DiagnosticFailureClass::invalid_input,
           "diagnostics.fp64.unavailable-bits",
           "unavailable FP64 value has bits");
    return;
  }
  if (classify_bits(value.bits) != value.status)
    fail(DiagnosticFailureClass::invalid_input, "diagnostics.fp64.status-bits",
         "FP64 status and bits disagree");
}

double fp64_value(DiagnosticFp64 value) noexcept {
  double result{};
  std::memcpy(&result, &value.bits, sizeof(result));
  return result;
}

template <class T, class Id>
void require_strict_order(const std::vector<T> &values, Id id,
                          const char *code) {
  for (std::size_t index = 0; index < values.size(); ++index) {
    require_id(id(values[index]), code, DiagnosticFailureClass::invalid_input);
    if (index != 0U && !(id(values[index - 1U]) < id(values[index])))
      fail(DiagnosticFailureClass::invalid_input,
           std::string("diagnostics.") + code + ".order",
           "diagnostic vector is not strictly ordered");
  }
}

std::string_view module_kind_name(DiagnosticModuleKind value) {
  switch (value) {
  case DiagnosticModuleKind::runtime:
    return "runtime";
  case DiagnosticModuleKind::mpi:
    return "mpi";
  case DiagnosticModuleKind::mesh_topology:
    return "mesh_topology";
  case DiagnosticModuleKind::mesh_geometry:
    return "mesh_geometry";
  case DiagnosticModuleKind::field:
    return "field";
  case DiagnosticModuleKind::execution:
    return "execution";
  case DiagnosticModuleKind::linear_operator:
    return "linear_operator";
  case DiagnosticModuleKind::linear_solver:
    return "linear_solver";
  case DiagnosticModuleKind::halo:
    return "halo";
  case DiagnosticModuleKind::boundary:
    return "boundary";
  case DiagnosticModuleKind::finite_volume:
    return "finite_volume";
  case DiagnosticModuleKind::piso:
    return "piso";
  case DiagnosticModuleKind::density_transport:
    return "density_transport";
  case DiagnosticModuleKind::density_closure:
    return "density_closure";
  case DiagnosticModuleKind::time_control:
    return "time_control";
  case DiagnosticModuleKind::checkpoint:
    return "checkpoint";
  case DiagnosticModuleKind::flow_driver:
    return "flow_driver";
  case DiagnosticModuleKind::performance:
    return "performance";
  }
  fail(DiagnosticFailureClass::invalid_input, "diagnostics.enum.module-kind",
       "invalid module kind");
}

std::string_view level_name(DiagnosticLevel value) {
  switch (value) {
  case DiagnosticLevel::summary:
    return "summary";
  case DiagnosticLevel::invariants:
    return "invariants";
  case DiagnosticLevel::counters:
    return "counters";
  case DiagnosticLevel::bounded_state_sample:
    return "bounded_state_sample";
  }
  fail(DiagnosticFailureClass::invalid_input, "diagnostics.enum.level",
       "invalid diagnostic level");
}
std::string_view scope_name(DiagnosticScope value) {
  switch (value) {
  case DiagnosticScope::local:
    return "local";
  case DiagnosticScope::collective:
    return "collective";
  }
  fail(DiagnosticFailureClass::invalid_input, "diagnostics.enum.scope",
       "invalid diagnostic scope");
}
std::string_view status_name(DiagnosticStatus value) {
  switch (value) {
  case DiagnosticStatus::ok:
    return "ok";
  case DiagnosticStatus::warning:
    return "warning";
  case DiagnosticStatus::failed:
    return "failed";
  case DiagnosticStatus::unavailable:
    return "unavailable";
  }
  fail(DiagnosticFailureClass::invalid_input, "diagnostics.enum.status",
       "invalid diagnostic status");
}
std::string_view failure_name(DiagnosticFailureClass value) {
  switch (value) {
  case DiagnosticFailureClass::none:
    return "none";
  case DiagnosticFailureClass::invalid_request:
    return "invalid_request";
  case DiagnosticFailureClass::invalid_input:
    return "invalid_input";
  case DiagnosticFailureClass::layout:
    return "layout";
  case DiagnosticFailureClass::capability:
    return "capability";
  case DiagnosticFailureClass::non_finite_state:
    return "non_finite_state";
  case DiagnosticFailureClass::non_positive_state:
    return "non_positive_state";
  case DiagnosticFailureClass::numerical_breakdown:
    return "numerical_breakdown";
  case DiagnosticFailureClass::non_convergence:
    return "non_convergence";
  case DiagnosticFailureClass::conservation:
    return "conservation";
  case DiagnosticFailureClass::boundary:
    return "boundary";
  case DiagnosticFailureClass::file_integrity:
    return "file_integrity";
  case DiagnosticFailureClass::collective_operation:
    return "collective_operation";
  case DiagnosticFailureClass::sink_failure:
    return "sink_failure";
  case DiagnosticFailureClass::unavailable:
    return "unavailable";
  }
  fail(DiagnosticFailureClass::invalid_input, "diagnostics.enum.failure",
       "invalid failure class");
}
std::string_view value_status_name(DiagnosticValueStatus value) {
  switch (value) {
  case DiagnosticValueStatus::finite:
    return "finite";
  case DiagnosticValueStatus::positive_infinity:
    return "positive_infinity";
  case DiagnosticValueStatus::negative_infinity:
    return "negative_infinity";
  case DiagnosticValueStatus::quiet_nan:
    return "quiet_nan";
  case DiagnosticValueStatus::signaling_nan:
    return "signaling_nan";
  case DiagnosticValueStatus::unavailable:
    return "unavailable";
  }
  fail(DiagnosticFailureClass::invalid_input, "diagnostics.enum.value-status",
       "invalid value status");
}
std::string_view relation_name(InvariantRelation value) {
  switch (value) {
  case InvariantRelation::less_equal:
    return "less_equal";
  case InvariantRelation::greater_equal:
    return "greater_equal";
  case InvariantRelation::equal:
    return "equal";
  case InvariantRelation::finite:
    return "finite";
  case InvariantRelation::positive:
    return "positive";
  }
  fail(DiagnosticFailureClass::invalid_input, "diagnostics.enum.relation",
       "invalid invariant relation");
}
std::string_view metric_kind_name(DiagnosticMetricKind value) {
  switch (value) {
  case DiagnosticMetricKind::state_summary:
    return "state_summary";
  case DiagnosticMetricKind::residual:
    return "residual";
  case DiagnosticMetricKind::conservation:
    return "conservation";
  case DiagnosticMetricKind::performance:
    return "performance";
  }
  fail(DiagnosticFailureClass::invalid_input, "diagnostics.enum.metric-kind",
       "invalid metric kind");
}

std::string json_string(std::string_view value) {
  if (!valid_utf8(value))
    fail(DiagnosticFailureClass::invalid_input, "diagnostics.json.utf8",
         "invalid UTF-8 string");
  std::string result{"\""};
  constexpr char digits[] = "0123456789abcdef";
  for (const char raw_character : value) {
    const auto character = static_cast<unsigned char>(raw_character);
    switch (character) {
    case '"':
      result += "\\\"";
      break;
    case '\\':
      result += "\\\\";
      break;
    case '\b':
      result += "\\b";
      break;
    case '\f':
      result += "\\f";
      break;
    case '\n':
      result += "\\n";
      break;
    case '\r':
      result += "\\r";
      break;
    case '\t':
      result += "\\t";
      break;
    default:
      if (character < 0x20U) {
        result += "\\u00";
        result += digits[character >> 4U];
        result += digits[character & 0x0fU];
      } else {
        result += static_cast<char>(character);
      }
    }
  }
  result += '"';
  return result;
}

std::string hex64(std::uint64_t value) {
  std::ostringstream output;
  output.imbue(std::locale::classic());
  output << std::hex << std::nouppercase << std::setfill('0') << std::setw(16)
         << value;
  return output.str();
}

std::string finite_number(DiagnosticFp64 value) {
  std::ostringstream output;
  output.imbue(std::locale::classic());
  output << std::setprecision(std::numeric_limits<double>::max_digits10)
         << fp64_value(value);
  std::string result = output.str();
  if (result.find('.') == std::string::npos &&
      result.find('e') == std::string::npos &&
      result.find('E') == std::string::npos)
    result += ".0";
  return result;
}

std::string fp64_json(DiagnosticFp64 value) {
  validate_fp64(value);
  std::string result =
      "{\"status\":" + json_string(value_status_name(value.status));
  if (value.status == DiagnosticValueStatus::finite)
    result += ",\"value\":" + finite_number(value);
  if (value.status != DiagnosticValueStatus::unavailable)
    result += ",\"bits\":\"0x" + hex64(value.bits) + "\"";
  result += '}';
  return result;
}

std::uint64_t crc64(const unsigned char *data, std::size_t size) noexcept {
  std::uint64_t crc = 0U;
  constexpr std::uint64_t polynomial = UINT64_C(0x42f0e1eba9ea3693);
  for (std::size_t index = 0; index < size; ++index) {
    crc ^= static_cast<std::uint64_t>(data[index]) << 56U;
    for (int bit = 0; bit < 8; ++bit)
      crc = (crc & UINT64_C(0x8000000000000000)) != 0U
                ? (crc << 1U) ^ polynomial
                : crc << 1U;
  }
  return crc;
}

template <class Integer>
void append_little_endian(std::vector<unsigned char> &bytes, Integer value) {
  for (std::size_t index = 0; index < sizeof(Integer); ++index)
    bytes.push_back(static_cast<unsigned char>(
        (value >> static_cast<unsigned>(8U * index)) & Integer{0xff}));
}

DiagnosticCapability level_capability(DiagnosticLevel level) {
  switch (level) {
  case DiagnosticLevel::summary:
    return DiagnosticCapability::summary;
  case DiagnosticLevel::invariants:
    return DiagnosticCapability::invariants;
  case DiagnosticLevel::counters:
    return DiagnosticCapability::counters;
  case DiagnosticLevel::bounded_state_sample:
    return DiagnosticCapability::bounded_state_sample;
  }
  fail(DiagnosticFailureClass::invalid_request, "diagnostics.request.level",
       "invalid requested level");
}

} // namespace

#ifdef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
namespace detail {
std::uint64_t structured_crc64_ecma_for_test(std::string_view bytes) noexcept {
  return crc64(reinterpret_cast<const unsigned char *>(bytes.data()),
               bytes.size());
}
} // namespace detail
#endif

DiagnosticCollectionError::DiagnosticCollectionError(
    DiagnosticFailureClass classification, std::string code,
    int lowest_failing_rank, std::string message)
    : std::runtime_error(std::move(message)), classification_(classification),
      code_(std::move(code)), lowest_failing_rank_(lowest_failing_rank) {}
DiagnosticFailureClass
DiagnosticCollectionError::classification() const noexcept {
  return classification_;
}
std::string_view DiagnosticCollectionError::code() const noexcept {
  return code_;
}
int DiagnosticCollectionError::lowest_failing_rank() const noexcept {
  return lowest_failing_rank_;
}

DiagnosticFp64 describe_fp64(double value) noexcept {
  DiagnosticFp64 result;
  std::memcpy(&result.bits, &value, sizeof(value));
  result.status = classify_bits(result.bits);
  return result;
}

bool evaluate_invariant(const DiagnosticInvariant &invariant) {
  validate_fp64(invariant.observed);
  validate_fp64(invariant.limit);
  if (invariant.observed.status == DiagnosticValueStatus::unavailable)
    fail(DiagnosticFailureClass::invalid_input,
         "diagnostics.invariant.observed", "invariant observed unavailable");
  const bool unary = invariant.relation == InvariantRelation::finite ||
                     invariant.relation == InvariantRelation::positive;
  if (unary) {
    if (invariant.limit.status != DiagnosticValueStatus::unavailable)
      fail(DiagnosticFailureClass::invalid_input,
           "diagnostics.invariant.unary-limit", "unary invariant has limit");
    if (invariant.relation == InvariantRelation::finite)
      return invariant.observed.status == DiagnosticValueStatus::finite;
    return invariant.observed.status == DiagnosticValueStatus::finite &&
           fp64_value(invariant.observed) > 0.0;
  }
  if (invariant.limit.status != DiagnosticValueStatus::finite)
    fail(DiagnosticFailureClass::invalid_input,
         "diagnostics.invariant.binary-limit",
         "binary invariant limit is not finite");
  if (invariant.observed.status != DiagnosticValueStatus::finite)
    return false;
  const double observed = fp64_value(invariant.observed);
  const double limit = fp64_value(invariant.limit);
  switch (invariant.relation) {
  case InvariantRelation::less_equal:
    return observed <= limit;
  case InvariantRelation::greater_equal:
    return observed >= limit;
  case InvariantRelation::equal:
    return observed == limit;
  case InvariantRelation::finite:
  case InvariantRelation::positive:
    break;
  }
  fail(DiagnosticFailureClass::invalid_input, "diagnostics.invariant.relation",
       "invalid invariant relation");
}

bool has_capability(DiagnosticCapabilityFlags flags,
                    DiagnosticCapability capability) noexcept {
  return (flags & static_cast<DiagnosticCapabilityFlags>(capability)) != 0U;
}

void validate(const DiagnosticDescriptor &descriptor) {
  if (descriptor.schema_version != kDiagnosticRecordSchemaV1)
    fail(DiagnosticFailureClass::invalid_input, "diagnostics.descriptor.schema",
         "unsupported descriptor schema");
  static_cast<void>(module_kind_name(descriptor.module_kind));
  require_id(descriptor.module_id, "module-id",
             DiagnosticFailureClass::invalid_input);
  require_id(descriptor.instance_id, "instance-id",
             DiagnosticFailureClass::invalid_input);
  constexpr DiagnosticCapabilityFlags allowed =
      static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::summary) |
      static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::invariants) |
      static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::counters) |
      static_cast<DiagnosticCapabilityFlags>(
          DiagnosticCapability::bounded_state_sample) |
      static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::collective);
  if ((descriptor.capabilities & ~allowed) != 0U)
    fail(DiagnosticFailureClass::invalid_input,
         "diagnostics.descriptor.capabilities", "invalid capability bits");
}

void validate(const DiagnosticRequest &request) {
  static_cast<void>(level_name(request.level));
  static_cast<void>(scope_name(request.scope));
  if (request.frame.rank < 0 || !std::isfinite(request.frame.time_s) ||
      request.frame.time_s < 0.0)
    fail(DiagnosticFailureClass::invalid_request, "diagnostics.request.frame",
         "invalid diagnostic frame");
  require_id(request.frame.phase, "phase",
             DiagnosticFailureClass::invalid_request);
  for (std::size_t index = 0; index < request.selected_fields.size(); ++index) {
    require_id(request.selected_fields[index], "selected-field",
               DiagnosticFailureClass::invalid_request);
    if (index != 0U &&
        !(request.selected_fields[index - 1U] < request.selected_fields[index]))
      fail(DiagnosticFailureClass::invalid_request,
           "diagnostics.request.selected-fields",
           "selected fields are not strictly ordered");
  }
  if (request.level == DiagnosticLevel::bounded_state_sample) {
    if (request.sample_budget == 0U ||
        request.sample_budget > kMaximumStateSamplesV1)
      fail(DiagnosticFailureClass::invalid_request,
           "diagnostics.request.sample-budget", "invalid sample budget");
  } else if (request.sample_budget != 0U || !request.selected_fields.empty()) {
    fail(DiagnosticFailureClass::invalid_request,
         "diagnostics.request.non-sample-selection",
         "non-sample request has sample selection");
  }
}

void validate(const DiagnosticRequest &request,
              const DiagnosticDescriptor &descriptor) {
  validate(descriptor);
  validate(request);
  if (!has_capability(descriptor.capabilities,
                      level_capability(request.level)) ||
      (request.scope == DiagnosticScope::collective &&
       !has_capability(descriptor.capabilities,
                       DiagnosticCapability::collective)))
    fail(DiagnosticFailureClass::capability, "diagnostics.request.capability",
         "unsupported diagnostic request");
}

void validate(const DiagnosticRecord &record) {
  if (record.schema_version != kDiagnosticRecordSchemaV1)
    fail(DiagnosticFailureClass::invalid_input, "diagnostics.record.schema",
         "unsupported record schema");
  static_cast<void>(module_kind_name(record.module_kind));
  static_cast<void>(level_name(record.level));
  static_cast<void>(scope_name(record.scope));
  static_cast<void>(status_name(record.status));
  static_cast<void>(failure_name(record.failure.classification));
  require_id(record.module_id, "module-id",
             DiagnosticFailureClass::invalid_input);
  require_id(record.instance_id, "instance-id",
             DiagnosticFailureClass::invalid_input);
  require_id(record.phase, "phase", DiagnosticFailureClass::invalid_input);
  require_id(record.failure.code, "failure-code",
             DiagnosticFailureClass::invalid_input);
  if (record.rank < 0)
    fail(DiagnosticFailureClass::invalid_input, "diagnostics.record.rank",
         "invalid record rank");
  validate_fp64(record.time_s);
  if (record.time_s.status != DiagnosticValueStatus::finite ||
      fp64_value(record.time_s) < 0.0)
    fail(DiagnosticFailureClass::invalid_input, "diagnostics.record.time",
         "invalid record time");

  const bool normal = record.status == DiagnosticStatus::ok ||
                      record.status == DiagnosticStatus::warning;
  if (normal) {
    if (record.failure.classification != DiagnosticFailureClass::none ||
        record.failure.code != "none" ||
        record.failure.lowest_failing_rank != -1)
      fail(DiagnosticFailureClass::invalid_input,
           "diagnostics.record.failure-normal", "invalid normal failure data");
  } else if (record.status == DiagnosticStatus::failed) {
    if (record.failure.classification == DiagnosticFailureClass::none ||
        record.failure.classification == DiagnosticFailureClass::unavailable ||
        record.failure.code == "none")
      fail(DiagnosticFailureClass::invalid_input,
           "diagnostics.record.failure-failed", "invalid failed-state data");
  } else if (record.failure.classification !=
                 DiagnosticFailureClass::unavailable ||
             record.failure.code == "none") {
    fail(DiagnosticFailureClass::invalid_input,
         "diagnostics.record.failure-unavailable",
         "invalid unavailable-state data");
  }
  if (record.scope == DiagnosticScope::local &&
      record.failure.lowest_failing_rank != -1)
    fail(DiagnosticFailureClass::invalid_input,
         "diagnostics.record.local-failing-rank",
         "local record has a failing rank");
  if (record.scope == DiagnosticScope::collective) {
    const bool records_failure = record.status == DiagnosticStatus::failed ||
                                 record.status == DiagnosticStatus::unavailable;
    if ((records_failure && record.failure.lowest_failing_rank < 0) ||
        (!records_failure && record.failure.lowest_failing_rank != -1))
      fail(DiagnosticFailureClass::invalid_input,
           "diagnostics.record.collective-failing-rank",
           "collective failing rank is invalid");
  }

  require_strict_order(
      record.invariants,
      [](const auto &value) -> std::string_view { return value.id; },
      "invariant-id");
  require_strict_order(
      record.metrics,
      [](const auto &value) -> std::string_view { return value.id; },
      "metric-id");
  require_strict_order(
      record.counters,
      [](const auto &value) -> std::string_view { return value.id; },
      "counter-id");
  require_strict_order(
      record.identities,
      [](const auto &value) -> std::string_view { return value.subject_id; },
      "identity-id");
  for (const auto &invariant : record.invariants) {
    require_unit(invariant.unit);
    if (evaluate_invariant(invariant) != invariant.passed)
      fail(DiagnosticFailureClass::invalid_input,
           "diagnostics.invariant.passed", "invariant result disagrees");
  }
  for (const auto &metric : record.metrics) {
    static_cast<void>(metric_kind_name(metric.kind));
    require_unit(metric.unit);
    validate_fp64(metric.value);
  }
  for (const auto &counter : record.counters)
    require_unit(counter.unit);
  for (const auto &identity : record.identities) {
    if (identity.layout_fingerprint)
      require_id(*identity.layout_fingerprint, "layout-fingerprint",
                 DiagnosticFailureClass::invalid_input);
  }
  for (std::size_t index = 0; index < record.samples.size(); ++index) {
    const auto &sample = record.samples[index];
    require_id(sample.field_id, "sample-field",
               DiagnosticFailureClass::invalid_input);
    require_unit(sample.unit);
    validate_fp64(sample.value);
    if (sample.value.status == DiagnosticValueStatus::unavailable)
      fail(DiagnosticFailureClass::invalid_input,
           "diagnostics.sample.unavailable", "sample value unavailable");
    if (index != 0U) {
      const auto previous = std::tie(record.samples[index - 1U].field_id,
                                     record.samples[index - 1U].global_id,
                                     record.samples[index - 1U].component);
      const auto current =
          std::tie(sample.field_id, sample.global_id, sample.component);
      if (!(previous < current))
        fail(DiagnosticFailureClass::invalid_input, "diagnostics.sample.order",
             "samples not strictly ordered");
    }
  }

  const bool summary = record.level == DiagnosticLevel::summary;
  const bool invariants = record.level == DiagnosticLevel::invariants;
  const bool counters = record.level == DiagnosticLevel::counters;
  const bool samples = record.level == DiagnosticLevel::bounded_state_sample;
  if ((summary && record.metrics.empty()) ||
      (invariants && record.invariants.empty()) ||
      (counters && record.counters.empty()) ||
      (!summary && !record.metrics.empty()) ||
      (!invariants && !record.invariants.empty()) ||
      (!counters && !record.counters.empty()))
    fail(DiagnosticFailureClass::invalid_input, "diagnostics.record.payload",
         "record primary payload is invalid");
  if (!samples) {
    if (record.sample_budget != 0U || record.eligible_sample_count != 0U ||
        record.samples_truncated || !record.samples.empty())
      fail(DiagnosticFailureClass::invalid_input,
           "diagnostics.record.non-sample", "non-sample record has samples");
  } else {
    if (record.sample_budget == 0U ||
        record.sample_budget > kMaximumStateSamplesV1 ||
        record.samples.size() > record.sample_budget ||
        record.samples.size() > record.eligible_sample_count ||
        record.samples.size() !=
            std::min<std::uint64_t>(record.sample_budget,
                                    record.eligible_sample_count) ||
        record.samples_truncated !=
            (record.eligible_sample_count > record.samples.size()))
      fail(DiagnosticFailureClass::invalid_input,
           "diagnostics.record.sample-shape", "sample record shape invalid");
  }
  if (record.state_fingerprint.algorithm != kStateFingerprintAlgorithmV1 ||
      record.state_fingerprint.hex.size() != 32U ||
      !std::all_of(record.state_fingerprint.hex.begin(),
                   record.state_fingerprint.hex.end(), [](char character) {
                     return (character >= '0' && character <= '9') ||
                            (character >= 'a' && character <= 'f');
                   }))
    fail(DiagnosticFailureClass::invalid_input,
         "diagnostics.record.fingerprint", "invalid state fingerprint");
}

void validate(const DiagnosticRecord &record,
              const DiagnosticDescriptor &descriptor,
              const DiagnosticRequest &request) {
  validate(request, descriptor);
  validate(record);
  if (record.schema_version != descriptor.schema_version ||
      record.module_kind != descriptor.module_kind ||
      record.module_id != descriptor.module_id ||
      record.instance_id != descriptor.instance_id ||
      record.level != request.level || record.scope != request.scope ||
      record.rank != request.frame.rank || record.step != request.frame.step ||
      record.time_s.bits != describe_fp64(request.frame.time_s).bits ||
      record.phase != request.frame.phase ||
      record.sample_budget != request.sample_budget)
    fail(DiagnosticFailureClass::invalid_input,
         "diagnostics.record.cross-object", "record contract mismatch");
  if (request.level == DiagnosticLevel::bounded_state_sample &&
      !request.selected_fields.empty()) {
    for (const auto &sample : record.samples) {
      if (!std::binary_search(request.selected_fields.begin(),
                              request.selected_fields.end(),
                              std::string_view(sample.field_id)))
        fail(DiagnosticFailureClass::invalid_input,
             "diagnostics.record.sample-selection",
             "record sample is outside requested fields");
    }
  }
}

void DiagnosticFingerprintAccumulator::add(std::string_view field_id,
                                           std::uint64_t global_id,
                                           std::uint32_t component,
                                           DiagnosticFp64 value) {
  require_id(field_id, "fingerprint-field",
             DiagnosticFailureClass::invalid_input);
  validate_fp64(value);
  if (value.status == DiagnosticValueStatus::unavailable)
    fail(DiagnosticFailureClass::invalid_input,
         "diagnostics.fingerprint.unavailable",
         "fingerprint value unavailable");
  std::vector<unsigned char> bytes;
  bytes.reserve(field_id.size() + 20U);
  bytes.insert(bytes.end(), field_id.begin(), field_id.end());
  append_little_endian(bytes, global_id);
  append_little_endian(bytes, component);
  append_little_endian(bytes, value.bits);
  const std::uint64_t item = crc64(bytes.data(), bytes.size());
  const DiagnosticFingerprintParts candidate{parts_.xor64 ^ item,
                                             parts_.sum64 + item};
  parts_ = candidate;
}
void DiagnosticFingerprintAccumulator::combine(
    DiagnosticFingerprintParts other) noexcept {
  parts_.xor64 ^= other.xor64;
  parts_.sum64 += other.sum64;
}
DiagnosticFingerprintParts
DiagnosticFingerprintAccumulator::parts() const noexcept {
  return parts_;
}
DiagnosticStateFingerprint DiagnosticFingerprintAccumulator::finish() const {
  return {std::string(kStateFingerprintAlgorithmV1),
          hex64(parts_.xor64) + hex64(parts_.sum64)};
}

std::string to_canonical_json(const DiagnosticRecord &record) {
  validate(record);
  std::string output;
  output += "{\"schema_version\":" + std::to_string(record.schema_version);
  output +=
      ",\"module_kind\":" + json_string(module_kind_name(record.module_kind));
  output += ",\"module_id\":" + json_string(record.module_id);
  output += ",\"instance_id\":" + json_string(record.instance_id);
  output += ",\"level\":" + json_string(level_name(record.level));
  output += ",\"scope\":" + json_string(scope_name(record.scope));
  output += ",\"rank\":" + std::to_string(record.rank);
  output += ",\"step\":" + std::to_string(record.step);
  output += ",\"time\":" + fp64_json(record.time_s);
  output += ",\"phase\":" + json_string(record.phase);
  output += ",\"status\":" + json_string(status_name(record.status));
  output += ",\"failure\":{\"classification\":" +
            json_string(failure_name(record.failure.classification)) +
            ",\"code\":" + json_string(record.failure.code) +
            ",\"lowest_failing_rank\":" +
            std::to_string(record.failure.lowest_failing_rank) + "}";
  output += ",\"invariants\":[";
  for (std::size_t i = 0; i < record.invariants.size(); ++i) {
    if (i != 0U)
      output += ',';
    const auto &value = record.invariants[i];
    output += "{\"id\":" + json_string(value.id) +
              ",\"unit\":" + json_string(value.unit) +
              ",\"observed\":" + fp64_json(value.observed) +
              ",\"limit\":" + fp64_json(value.limit) +
              ",\"relation\":" + json_string(relation_name(value.relation)) +
              ",\"passed\":" + (value.passed ? "true" : "false") + "}";
  }
  output += "],\"metrics\":[";
  for (std::size_t i = 0; i < record.metrics.size(); ++i) {
    if (i != 0U)
      output += ',';
    const auto &value = record.metrics[i];
    output += "{\"id\":" + json_string(value.id) +
              ",\"kind\":" + json_string(metric_kind_name(value.kind)) +
              ",\"unit\":" + json_string(value.unit) +
              ",\"value\":" + fp64_json(value.value) + "}";
  }
  output += "],\"counters\":[";
  for (std::size_t i = 0; i < record.counters.size(); ++i) {
    if (i != 0U)
      output += ',';
    const auto &value = record.counters[i];
    output += "{\"id\":" + json_string(value.id) +
              ",\"unit\":" + json_string(value.unit) +
              ",\"value\":" + std::to_string(value.value) + "}";
  }
  output += "],\"identities\":[";
  for (std::size_t i = 0; i < record.identities.size(); ++i) {
    if (i != 0U)
      output += ',';
    const auto &value = record.identities[i];
    output +=
        "{\"subject_id\":" + json_string(value.subject_id) +
        ",\"layout_fingerprint\":" +
        (value.layout_fingerprint ? json_string(*value.layout_fingerprint)
                                  : "null") +
        ",\"revision\":" +
        (value.revision ? std::to_string(*value.revision) : "null") +
        ",\"generation\":" +
        (value.generation ? std::to_string(*value.generation) : "null") +
        ",\"allocation_identity\":" +
        (value.allocation_identity ? std::to_string(*value.allocation_identity)
                                   : "null") +
        "}";
  }
  output += "],\"state_fingerprint\":{\"algorithm\":" +
            json_string(record.state_fingerprint.algorithm) +
            ",\"hex\":" + json_string(record.state_fingerprint.hex) + "}";
  output += ",\"sample_budget\":" + std::to_string(record.sample_budget);
  output += ",\"eligible_sample_count\":" +
            std::to_string(record.eligible_sample_count);
  output += ",\"samples_truncated\":" +
            std::string(record.samples_truncated ? "true" : "false");
  output += ",\"samples\":[";
  for (std::size_t i = 0; i < record.samples.size(); ++i) {
    if (i != 0U)
      output += ',';
    const auto &value = record.samples[i];
    output += "{\"field_id\":" + json_string(value.field_id) +
              ",\"global_id\":" + std::to_string(value.global_id) +
              ",\"component\":" + std::to_string(value.component) +
              ",\"unit\":" + json_string(value.unit) +
              ",\"value\":" + fp64_json(value.value) + "}";
  }
  output += "]}";
  return output;
}

} // namespace hundun::diagnostics

// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_les_wale.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using namespace hundun;

std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

double from_bits(std::uint64_t value) noexcept {
  double result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

class RecordingSink final : public diagnostics::DiagnosticSink {
public:
  void submit(const diagnostics::DiagnosticRecord &record) override {
    records.push_back(record);
  }

  std::vector<diagnostics::DiagnosticRecord> records;
};

diagnostics::DiagnosticRequest request(diagnostics::DiagnosticLevel level) {
  diagnostics::DiagnosticRequest result;
  result.level = level;
  result.scope = diagnostics::DiagnosticScope::local;
  result.frame = {3, 17U, 0.125, "wale.attempt-summary"};
  return result;
}

les::WaleSummary summary() {
  return {{UINT64_C(0xfedcba9876543210)},
          0.0,
          0.03125,
          0.0125,
          UINT64_C(0x0020000000000003),
          UINT64_C(0x0020000000000011)};
}

bool same_summary(const les::WaleSummary &left,
                  const les::WaleSummary &right) noexcept {
  return left.identity == right.identity &&
         bits(left.minimum_nu_t_m2_per_s) ==
             bits(right.minimum_nu_t_m2_per_s) &&
         bits(left.maximum_nu_t_m2_per_s) ==
             bits(right.maximum_nu_t_m2_per_s) &&
         bits(left.l2_nu_t_m2_per_s) == bits(right.l2_nu_t_m2_per_s) &&
         left.exact_zero_count == right.exact_zero_count &&
         left.owned_active_count == right.owned_active_count;
}

diagnostics::DiagnosticRecord
collect_one(const les::WaleSummary &source,
            const diagnostics::DiagnosticRequest &diagnostic_request) {
  RecordingSink sink;
  try {
    diagnostics::collect_diagnostics(source, diagnostic_request, sink);
  } catch (const diagnostics::DiagnosticCollectionError &error) {
    HUNDUN_CHECK(error.code() != "wale.diagnostics.capability");
    throw;
  }
  HUNDUN_CHECK(sink.records.size() == 1U);
  return sink.records.front();
}

void check_descriptor_and_fields() {
  using namespace diagnostics;
  const auto source = summary();
  const auto descriptor = describe_diagnostics(source);
  HUNDUN_CHECK(descriptor.schema_version == kDiagnosticRecordSchemaV1);
  HUNDUN_CHECK(descriptor.module_kind == DiagnosticModuleKind::les);
  HUNDUN_CHECK(descriptor.module_id == "hundun.les.wale");
  HUNDUN_CHECK(descriptor.instance_id == "primary");
  HUNDUN_CHECK(has_capability(descriptor.capabilities,
                              DiagnosticCapability::summary));
  HUNDUN_CHECK(has_capability(descriptor.capabilities,
                              DiagnosticCapability::counters));
  HUNDUN_CHECK(!has_capability(descriptor.capabilities,
                               DiagnosticCapability::invariants));
  HUNDUN_CHECK(!has_capability(
      descriptor.capabilities, DiagnosticCapability::bounded_state_sample));
  HUNDUN_CHECK(!has_capability(descriptor.capabilities,
                               DiagnosticCapability::collective));
  validate(descriptor);

  const std::vector<std::string_view> expected{
      "identity", "nu-t.exact-zero-count", "nu-t.l2", "nu-t.maximum",
      "nu-t.minimum", "owned-active-count"};
  HUNDUN_CHECK(diagnostic_fingerprint_field_ids(source) == expected);
}

void check_added_module_kind_canonical_names() {
  using namespace diagnostics;
  auto record = collect_one(summary(), request(DiagnosticLevel::summary));
  constexpr std::array<std::pair<DiagnosticModuleKind, std::string_view>, 5>
      expected{{
          {DiagnosticModuleKind::immersed_surface, "immersed_surface"},
          {DiagnosticModuleKind::ghost_stencil, "ghost_stencil"},
          {DiagnosticModuleKind::local_flow_pattern, "local_flow_pattern"},
          {DiagnosticModuleKind::wall_force, "wall_force"},
          {DiagnosticModuleKind::les, "les"},
      }};
  constexpr std::string_view prefix = "\"module_kind\":\"";
  for (const auto &[kind, name] : expected) {
    record.module_kind = kind;
    const auto canonical = to_canonical_json(record);
    const auto begin = canonical.find(prefix);
    HUNDUN_CHECK(begin != std::string::npos);
    const auto value_begin = begin + prefix.size();
    const auto end = canonical.find('"', value_begin);
    HUNDUN_CHECK(end != std::string::npos);
    HUNDUN_CHECK(std::string_view(canonical).substr(value_begin,
                                                     end - value_begin) ==
                 name);
    HUNDUN_CHECK(canonical.find(prefix, value_begin) == std::string::npos);
  }
}

void check_summary_and_counters() {
  using namespace diagnostics;
  const auto source = summary();
  const auto summary_record =
      collect_one(source, request(DiagnosticLevel::summary));
  validate(summary_record, describe_diagnostics(source),
           request(DiagnosticLevel::summary));
  HUNDUN_CHECK(summary_record.module_kind == DiagnosticModuleKind::les);
  HUNDUN_CHECK(summary_record.module_id == "hundun.les.wale");
  HUNDUN_CHECK(summary_record.instance_id == "primary");
  HUNDUN_CHECK(summary_record.metrics.size() == 3U);
  HUNDUN_CHECK(summary_record.metrics[0U].id == "nu-t.l2");
  HUNDUN_CHECK(summary_record.metrics[0U].kind ==
               DiagnosticMetricKind::state_summary);
  HUNDUN_CHECK(summary_record.metrics[0U].unit == "m2/s");
  HUNDUN_CHECK(summary_record.metrics[0U].value.bits ==
               bits(source.l2_nu_t_m2_per_s));
  HUNDUN_CHECK(summary_record.metrics[1U].id == "nu-t.maximum");
  HUNDUN_CHECK(summary_record.metrics[1U].unit == "m2/s");
  HUNDUN_CHECK(summary_record.metrics[1U].value.bits ==
               bits(source.maximum_nu_t_m2_per_s));
  HUNDUN_CHECK(summary_record.metrics[2U].id == "nu-t.minimum");
  HUNDUN_CHECK(summary_record.metrics[2U].unit == "m2/s");
  HUNDUN_CHECK(summary_record.metrics[2U].value.bits ==
               bits(source.minimum_nu_t_m2_per_s));

  const auto counter_record =
      collect_one(source, request(DiagnosticLevel::counters));
  validate(counter_record, describe_diagnostics(source),
           request(DiagnosticLevel::counters));
  HUNDUN_CHECK(counter_record.counters.size() == 3U);
  HUNDUN_CHECK(counter_record.counters[0U].id == "identity");
  HUNDUN_CHECK(counter_record.counters[0U].unit == "count");
  HUNDUN_CHECK(counter_record.counters[0U].value == source.identity.value);
  HUNDUN_CHECK(counter_record.counters[1U].id ==
               "nu-t.exact-zero-count");
  HUNDUN_CHECK(counter_record.counters[1U].unit == "count");
  HUNDUN_CHECK(counter_record.counters[1U].value ==
               source.exact_zero_count);
  HUNDUN_CHECK(counter_record.counters[2U].id == "owned-active-count");
  HUNDUN_CHECK(counter_record.counters[2U].unit == "count");
  HUNDUN_CHECK(counter_record.counters[2U].value ==
               source.owned_active_count);
  HUNDUN_CHECK(counter_record.state_fingerprint.hex ==
               summary_record.state_fingerprint.hex);
}

void check_nonfinite_status_encoding() {
  using namespace diagnostics;
  auto source = summary();
  source.minimum_nu_t_m2_per_s =
      std::numeric_limits<double>::infinity();
  source.maximum_nu_t_m2_per_s =
      -std::numeric_limits<double>::infinity();
  source.l2_nu_t_m2_per_s = from_bits(UINT64_C(0x7ff8000000001234));
  const auto record = collect_one(source, request(DiagnosticLevel::summary));
  HUNDUN_CHECK(record.metrics[0U].value.status ==
               DiagnosticValueStatus::quiet_nan);
  HUNDUN_CHECK(record.metrics[0U].value.bits ==
               UINT64_C(0x7ff8000000001234));
  HUNDUN_CHECK(record.metrics[1U].value.status ==
               DiagnosticValueStatus::negative_infinity);
  HUNDUN_CHECK(record.metrics[2U].value.status ==
               DiagnosticValueStatus::positive_infinity);
  validate(record, describe_diagnostics(source),
           request(DiagnosticLevel::summary));
}

void check_repeat_collection_is_read_only_and_canonical() {
  using namespace diagnostics;
  const auto source = summary();
  const auto before = source;
  const auto first = collect_one(source, request(DiagnosticLevel::summary));
  const auto second = collect_one(source, request(DiagnosticLevel::summary));
  HUNDUN_CHECK(to_canonical_json(first) == to_canonical_json(second));
  HUNDUN_CHECK(same_summary(source, before));

  std::array<les::WaleSummary, 6> mutations{source, source, source,
                                             source, source, source};
  ++mutations[0U].identity.value;
  mutations[1U].minimum_nu_t_m2_per_s =
      std::nextafter(mutations[1U].minimum_nu_t_m2_per_s, 1.0);
  mutations[2U].maximum_nu_t_m2_per_s =
      std::nextafter(mutations[2U].maximum_nu_t_m2_per_s, 1.0);
  mutations[3U].l2_nu_t_m2_per_s =
      std::nextafter(mutations[3U].l2_nu_t_m2_per_s, 1.0);
  ++mutations[4U].exact_zero_count;
  ++mutations[5U].owned_active_count;
  for (const auto &changed : mutations) {
    const auto changed_record =
        collect_one(changed, request(DiagnosticLevel::summary));
    HUNDUN_CHECK(changed_record.state_fingerprint.hex !=
                 first.state_fingerprint.hex);
  }
  HUNDUN_CHECK(same_summary(source, before));
}

void check_absence_registers_nothing() {
  using namespace diagnostics;
  std::optional<les::WaleSummary> source;
  std::vector<DiagnosticDescriptor> descriptors;
  RecordingSink sink;
  const auto register_if_present = [&](const auto &candidate) {
    if (!candidate.has_value()) {
      return;
    }
    descriptors.push_back(describe_diagnostics(*candidate));
    collect_diagnostics(*candidate, request(DiagnosticLevel::summary), sink);
  };
  register_if_present(source);
  HUNDUN_CHECK(descriptors.empty());
  HUNDUN_CHECK(sink.records.empty());

  source = summary();
  register_if_present(source);
  HUNDUN_CHECK(descriptors.size() == 1U);
  HUNDUN_CHECK(sink.records.size() == 1U);
}

void check_unsupported_requests_fail_as_capability() {
  using namespace diagnostics;
  const auto source = summary();
  for (auto diagnostic_request : {
           request(DiagnosticLevel::invariants),
           request(DiagnosticLevel::summary),
       }) {
    if (diagnostic_request.level == DiagnosticLevel::summary) {
      diagnostic_request.scope = DiagnosticScope::collective;
    }
    RecordingSink sink;
    bool rejected = false;
    try {
      collect_diagnostics(source, diagnostic_request, sink);
    } catch (const DiagnosticCollectionError &error) {
      rejected = error.classification() == DiagnosticFailureClass::capability &&
                 error.code() != "wale.diagnostics.capability";
    }
    HUNDUN_CHECK(rejected);
    HUNDUN_CHECK(sink.records.empty());
  }
}

void run() {
  check_descriptor_and_fields();
  check_added_module_kind_canonical_names();
  check_summary_and_counters();
  check_nonfinite_status_encoding();
  check_repeat_collection_is_read_only_and_canonical();
  check_absence_registers_nothing();
  check_unsupported_requests_fail_as_capability();
}

} // namespace

int main() { return hundun::test::run(run); }

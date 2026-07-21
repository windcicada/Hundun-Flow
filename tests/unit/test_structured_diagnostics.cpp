// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/structured_diagnostics.hpp"
#include "tests/support/test_main.hpp"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

double from_bits(std::uint64_t bits) noexcept {
  double value{};
  std::memcpy(&value, &bits, sizeof(value));
  return value;
}

hundun::diagnostics::DiagnosticRecord summary_record() {
  using namespace hundun::diagnostics;
  DiagnosticRecord record;
  record.module_kind = DiagnosticModuleKind::density_transport;
  record.module_id = "hundun.flow.material_density_transport";
  record.instance_id = "primary";
  record.level = DiagnosticLevel::summary;
  record.scope = DiagnosticScope::local;
  record.rank = 0;
  record.step = 3;
  record.time_s = describe_fp64(-0.0);
  record.phase = "material.finalized-trial";
  record.metrics.push_back({"density.minimum",
                            DiagnosticMetricKind::state_summary, "kg/m3",
                            describe_fp64(0.8)});
  record.state_fingerprint = {std::string(kStateFingerprintAlgorithmV1),
                              "00000000000000000000000000000000"};
  return record;
}

template <class Function> void expect_rejected(Function &&function) {
  bool rejected = false;
  try {
    function();
  } catch (const hundun::diagnostics::DiagnosticCollectionError &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

void test_fp64() {
  using namespace hundun::diagnostics;
  HUNDUN_CHECK(describe_fp64(0.0).bits == UINT64_C(0));
  HUNDUN_CHECK(describe_fp64(-0.0).bits == UINT64_C(0x8000000000000000));
  HUNDUN_CHECK(describe_fp64(std::numeric_limits<double>::infinity()).status ==
               DiagnosticValueStatus::positive_infinity);
  HUNDUN_CHECK(describe_fp64(-std::numeric_limits<double>::infinity()).status ==
               DiagnosticValueStatus::negative_infinity);
  HUNDUN_CHECK(describe_fp64(from_bits(UINT64_C(0x7ff8000000001234))).status ==
               DiagnosticValueStatus::quiet_nan);
  HUNDUN_CHECK(describe_fp64(from_bits(UINT64_C(0x7ff0000000001234))).status ==
               DiagnosticValueStatus::signaling_nan);
}

void test_invariants() {
  using namespace hundun::diagnostics;
  DiagnosticInvariant equal{"signed-zero",
                            "1",
                            describe_fp64(+0.0),
                            describe_fp64(-0.0),
                            InvariantRelation::equal,
                            true};
  HUNDUN_CHECK(evaluate_invariant(equal));
  DiagnosticInvariant positive{
      "positive", "kg/m3", describe_fp64(0.1), {}, InvariantRelation::positive,
      true};
  HUNDUN_CHECK(evaluate_invariant(positive));
  positive.observed = describe_fp64(-0.0);
  positive.passed = false;
  HUNDUN_CHECK(!evaluate_invariant(positive));
  positive.passed = true;
  DiagnosticRecord record = summary_record();
  record.level = DiagnosticLevel::invariants;
  record.metrics.clear();
  record.invariants.push_back(positive);
  expect_rejected([&] { validate(record); });
}

void test_validation_and_json() {
  using namespace hundun::diagnostics;
  const DiagnosticDescriptor descriptor{
      1, DiagnosticModuleKind::density_transport,
      "hundun.flow.material_density_transport", "primary",
      static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::summary) |
          static_cast<DiagnosticCapabilityFlags>(
              DiagnosticCapability::collective)};
  DiagnosticRequest request;
  request.level = DiagnosticLevel::summary;
  request.scope = DiagnosticScope::local;
  request.frame = {0, 3, -0.0, "material.finalized-trial"};
  DiagnosticRecord record = summary_record();
  validate(record, descriptor, request);
  const std::string first = to_canonical_json(record);
  const std::string second = to_canonical_json(record);
  HUNDUN_CHECK(first == second);
  HUNDUN_CHECK(first.find("\"value\":-0.0") != std::string::npos);
  HUNDUN_CHECK(first.find("NaN") == std::string::npos);
  HUNDUN_CHECK(first.find("Infinity") == std::string::npos);

  auto bad = record;
  bad.metrics[0].unit = "J/m^3";
  expect_rejected([&] { validate(bad); });
  bad = record;
  bad.metrics.push_back(bad.metrics.front());
  expect_rejected([&] { validate(bad); });
  bad = record;
  bad.module_id = "Bad.Id";
  expect_rejected([&] { validate(bad); });
}

void test_fingerprint() {
  using namespace hundun::diagnostics;
  DiagnosticFingerprintAccumulator empty;
  HUNDUN_CHECK(empty.finish().hex == "00000000000000000000000000000000");
  DiagnosticFingerprintAccumulator value;
  value.add("rho", 7, 0, describe_fp64(1.25));
  const auto before = value.parts();
  expect_rejected([&] { value.add("Bad", 7, 0, describe_fp64(2.0)); });
  HUNDUN_CHECK(value.parts().xor64 == before.xor64);
  HUNDUN_CHECK(value.parts().sum64 == before.sum64);
  expect_rejected([&] { value.add("rho", 7, 0, {}); });
  HUNDUN_CHECK(value.parts().xor64 == before.xor64);
  HUNDUN_CHECK(value.parts().sum64 == before.sum64);
}

void test_all_status_level_shapes() {
  using namespace hundun::diagnostics;
  const DiagnosticCapabilityFlags all =
      static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::summary) |
      static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::invariants) |
      static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::counters) |
      static_cast<DiagnosticCapabilityFlags>(
          DiagnosticCapability::bounded_state_sample) |
      static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::collective);
  const DiagnosticDescriptor descriptor{
      1, DiagnosticModuleKind::density_transport,
      "hundun.flow.material_density_transport", "primary", all};
  constexpr DiagnosticStatus statuses[]{
      DiagnosticStatus::ok, DiagnosticStatus::warning, DiagnosticStatus::failed,
      DiagnosticStatus::unavailable};
  constexpr DiagnosticLevel levels[]{
      DiagnosticLevel::summary, DiagnosticLevel::invariants,
      DiagnosticLevel::counters, DiagnosticLevel::bounded_state_sample};
  constexpr DiagnosticScope scopes[]{DiagnosticScope::local,
                                     DiagnosticScope::collective};
  for (const auto scope : scopes) {
    for (const auto status : statuses) {
      for (const auto level : levels) {
        DiagnosticRequest request;
        request.level = level;
        request.scope = scope;
        request.frame = {0, 9U, -0.0, "material.finalized-trial"};
        if (level == DiagnosticLevel::bounded_state_sample)
          request.sample_budget = 1U;
        DiagnosticRecord record;
        record.module_kind = DiagnosticModuleKind::density_transport;
        record.module_id = "hundun.flow.material_density_transport";
        record.instance_id = "primary";
        record.level = level;
        record.scope = scope;
        record.rank = 0;
        record.step = 9U;
        record.time_s = describe_fp64(-0.0);
        record.phase = "material.finalized-trial";
        record.status = status;
        if (status == DiagnosticStatus::failed) {
          record.failure = {DiagnosticFailureClass::conservation,
                            "material.conservation",
                            scope == DiagnosticScope::collective ? 0 : -1};
        } else if (status == DiagnosticStatus::unavailable) {
          record.failure = {DiagnosticFailureClass::unavailable,
                            "material.unavailable",
                            scope == DiagnosticScope::collective ? 0 : -1};
        }
        record.state_fingerprint = {std::string(kStateFingerprintAlgorithmV1),
                                    "00000000000000000000000000000000"};
        if (level == DiagnosticLevel::summary) {
          record.metrics.push_back({"density.minimum",
                                    DiagnosticMetricKind::state_summary,
                                    "kg/m3", describe_fp64(1.0)});
        } else if (level == DiagnosticLevel::invariants) {
          DiagnosticInvariant invariant{
              "density.available",      "1",
              describe_fp64(0.0),       describe_fp64(1.0),
              InvariantRelation::equal, false};
          if (status != DiagnosticStatus::unavailable) {
            invariant.observed = describe_fp64(1.0);
            invariant.passed = true;
          }
          record.invariants.push_back(invariant);
        } else if (level == DiagnosticLevel::counters) {
          record.counters.push_back({"owned.cells", "count", 1U});
        } else {
          record.sample_budget = 1U;
          record.eligible_sample_count = 1U;
          record.samples.push_back(
              {"rho", 0U, 0U, "kg/m3", describe_fp64(1.0)});
        }
        validate(record, descriptor, request);
        HUNDUN_CHECK(!to_canonical_json(record).empty());
      }
    }
  }
}

} // namespace

int main() {
  return hundun::test::run([] {
    test_fp64();
    test_invariants();
    test_validation_and_json();
    test_fingerprint();
    test_all_status_level_shapes();
  });
}

// SPDX-License-Identifier: Apache-2.0

#include "diagnostics/src/structured_diagnostics_test_access.hpp"
#include "hundun/diagnostics/structured_diagnostics.hpp"
#include "tests/support/test_main.hpp"

#include <array>
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

class CountingSink final : public hundun::diagnostics::DiagnosticSink {
public:
  void submit(const hundun::diagnostics::DiagnosticRecord &) override {
    ++calls;
  }
  std::size_t calls{};
};

template <class Mutation>
void expect_record_mutation_rejected(
    const hundun::diagnostics::DiagnosticRecord &baseline,
    const hundun::diagnostics::DiagnosticDescriptor &descriptor,
    const hundun::diagnostics::DiagnosticRequest &request,
    Mutation &&mutation) {
  auto changed = baseline;
  mutation(changed);
  CountingSink sink;
  expect_rejected([&] {
    hundun::diagnostics::validate(changed, descriptor, request);
    sink.submit(changed);
  });
  HUNDUN_CHECK(sink.calls == 0U);
}

hundun::diagnostics::DiagnosticCapabilityFlags all_capabilities() {
  using namespace hundun::diagnostics;
  return static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::summary) |
         static_cast<DiagnosticCapabilityFlags>(
             DiagnosticCapability::invariants) |
         static_cast<DiagnosticCapabilityFlags>(
             DiagnosticCapability::counters) |
         static_cast<DiagnosticCapabilityFlags>(
             DiagnosticCapability::bounded_state_sample) |
         static_cast<DiagnosticCapabilityFlags>(
             DiagnosticCapability::collective);
}

hundun::diagnostics::DiagnosticDescriptor full_descriptor() {
  using namespace hundun::diagnostics;
  return {1, DiagnosticModuleKind::density_transport,
          "hundun.flow.material_density_transport", "primary",
          all_capabilities()};
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
  HUNDUN_CHECK(describe_fp64(from_bits(UINT64_C(0x7ff8000000001234))).bits ==
               UINT64_C(0x7ff8000000001234));

  for (const auto represented :
       {describe_fp64(1.0), describe_fp64(-0.0),
        describe_fp64(std::numeric_limits<double>::infinity()),
        describe_fp64(-std::numeric_limits<double>::infinity()),
        describe_fp64(from_bits(UINT64_C(0x7ff8000000001234))),
        describe_fp64(from_bits(UINT64_C(0x7ff0000000001234))),
        DiagnosticFp64{}}) {
    auto record = summary_record();
    record.metrics[0].value = represented;
    validate(record);
    const auto json = to_canonical_json(record);
    HUNDUN_CHECK(json.find(":NaN") == std::string::npos);
    HUNDUN_CHECK(json.find(":Infinity") == std::string::npos);
  }

  auto mismatched = summary_record();
  mismatched.metrics[0].value = {DiagnosticValueStatus::finite,
                                 UINT64_C(0x7ff0000000000000)};
  expect_rejected([&] { validate(mismatched); });
  mismatched = summary_record();
  mismatched.metrics[0].value = {DiagnosticValueStatus::unavailable,
                                 UINT64_C(1)};
  expect_rejected([&] { validate(mismatched); });
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

  const std::array<DiagnosticInvariant, 5> relations{
      DiagnosticInvariant{"less", "1", describe_fp64(1.0), describe_fp64(2.0),
                          InvariantRelation::less_equal, true},
      DiagnosticInvariant{"greater", "1", describe_fp64(2.0),
                          describe_fp64(1.0), InvariantRelation::greater_equal,
                          true},
      DiagnosticInvariant{"equal", "1", describe_fp64(+0.0),
                          describe_fp64(-0.0), InvariantRelation::equal, true},
      DiagnosticInvariant{"finite",
                          "1",
                          describe_fp64(1.0),
                          {},
                          InvariantRelation::finite,
                          true},
      DiagnosticInvariant{"positive",
                          "1",
                          describe_fp64(1.0),
                          {},
                          InvariantRelation::positive,
                          true}};
  for (const auto &invariant : relations)
    HUNDUN_CHECK(evaluate_invariant(invariant));

  auto invalid = relations[0];
  invalid.observed = {};
  expect_rejected([&] { static_cast<void>(evaluate_invariant(invalid)); });
  invalid = relations[0];
  invalid.limit = describe_fp64(std::numeric_limits<double>::infinity());
  expect_rejected([&] { static_cast<void>(evaluate_invariant(invalid)); });
  invalid = relations[3];
  invalid.limit = describe_fp64(0.0);
  expect_rejected([&] { static_cast<void>(evaluate_invariant(invalid)); });
  invalid = relations[0];
  invalid.observed = describe_fp64(from_bits(UINT64_C(0x7ff8000000001234)));
  HUNDUN_CHECK(!evaluate_invariant(invalid));
  invalid = relations[3];
  invalid.observed = describe_fp64(std::numeric_limits<double>::infinity());
  HUNDUN_CHECK(!evaluate_invariant(invalid));

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
  const auto module_pointer = record.module_id.data();
  const auto module_capacity = record.module_id.capacity();
  const auto phase_pointer = record.phase.data();
  const auto phase_capacity = record.phase.capacity();
  const auto metric_pointer = record.metrics.data();
  const auto metric_capacity = record.metrics.capacity();
  const auto metric_id_pointer = record.metrics[0].id.data();
  const auto metric_id_capacity = record.metrics[0].id.capacity();
  const auto original = record;
  const std::string first = to_canonical_json(record);
  const std::string second = to_canonical_json(record);
  HUNDUN_CHECK(first == second);
  const std::string exact =
      "{\"schema_version\":1,\"module_kind\":\"density_transport\","
      "\"module_id\":\"hundun.flow.material_density_transport\","
      "\"instance_id\":\"primary\",\"level\":\"summary\","
      "\"scope\":\"local\",\"rank\":0,\"step\":3,\"time\":{"
      "\"status\":\"finite\",\"value\":-0.0,"
      "\"bits\":\"0x8000000000000000\"},"
      "\"phase\":\"material.finalized-trial\",\"status\":\"ok\","
      "\"failure\":{\"classification\":\"none\",\"code\":\"none\","
      "\"lowest_failing_rank\":-1},\"invariants\":[],\"metrics\":[{"
      "\"id\":\"density.minimum\",\"kind\":\"state_summary\","
      "\"unit\":\"kg/m3\",\"value\":{\"status\":\"finite\","
      "\"value\":0.80000000000000004,"
      "\"bits\":\"0x3fe999999999999a\"}}],\"counters\":[],"
      "\"identities\":[],\"state_fingerprint\":{"
      "\"algorithm\":\"hundun-state-fp-v1\","
      "\"hex\":\"00000000000000000000000000000000\"},"
      "\"sample_budget\":0,\"eligible_sample_count\":0,"
      "\"samples_truncated\":false,\"samples\":[]}";
  HUNDUN_CHECK(first == exact);
  HUNDUN_CHECK(first.find("\"value\":-0.0") != std::string::npos);
  HUNDUN_CHECK(first.find("NaN") == std::string::npos);
  HUNDUN_CHECK(first.find("Infinity") == std::string::npos);
  HUNDUN_CHECK(record.module_id.data() == module_pointer);
  HUNDUN_CHECK(record.module_id.capacity() == module_capacity);
  HUNDUN_CHECK(record.phase.data() == phase_pointer);
  HUNDUN_CHECK(record.phase.capacity() == phase_capacity);
  HUNDUN_CHECK(record.metrics.data() == metric_pointer);
  HUNDUN_CHECK(record.metrics.capacity() == metric_capacity);
  HUNDUN_CHECK(record.metrics[0].id.data() == metric_id_pointer);
  HUNDUN_CHECK(record.metrics[0].id.capacity() == metric_id_capacity);
  HUNDUN_CHECK(std::memcmp(record.module_id.data(), original.module_id.data(),
                           record.module_id.size()) == 0);
  HUNDUN_CHECK(record.metrics[0].value.bits == original.metrics[0].value.bits);

  for (const std::string unit : {"degree", "J/m3"}) {
    auto accepted = record;
    accepted.metrics[0].unit = unit;
    validate(accepted);
  }
  for (const std::string unit :
       {"J/m^3", "J/m3 ", "J/m\xc2\xb3", "j/m3", "Degree", "DEGREE"}) {
    auto bad = record;
    bad.metrics[0].unit = unit;
    expect_rejected([&] { validate(bad); });
  }
  auto bad = record;
  bad.metrics.push_back(bad.metrics.front());
  expect_rejected([&] { validate(bad); });
  bad = record;
  bad.module_id = "Bad.Id";
  expect_rejected([&] { validate(bad); });
}

void test_fingerprint() {
  using namespace hundun::diagnostics;
  HUNDUN_CHECK(test::crc64_ecma("123456789") == UINT64_C(0x6c40df5f0b497347));
  DiagnosticFingerprintAccumulator empty;
  HUNDUN_CHECK(empty.finish().hex == "00000000000000000000000000000000");
  DiagnosticFingerprintAccumulator value;
  value.add("rho", 7, 0, describe_fp64(1.25));
  const auto before = value.parts();
  HUNDUN_CHECK(before.xor64 == UINT64_C(0x8bc2f00e0d369440));
  HUNDUN_CHECK(before.sum64 == UINT64_C(0x8bc2f00e0d369440));
  HUNDUN_CHECK(value.finish().hex == "8bc2f00e0d3694408bc2f00e0d369440");
  const std::string invalid_utf8{"\xc3\x28", 2U};
  const std::array<std::string, 5> invalid_ids{
      "", "Bad", "bad/id", std::string(129U, 'a'), invalid_utf8};
  for (const auto &id : invalid_ids) {
    expect_rejected([&] { value.add(id, 7, 0, describe_fp64(2.0)); });
    HUNDUN_CHECK(value.parts().xor64 == before.xor64);
    HUNDUN_CHECK(value.parts().sum64 == before.sum64);
  }
  for (const auto invalid_value :
       {DiagnosticFp64{}, DiagnosticFp64{DiagnosticValueStatus::finite,
                                         UINT64_C(0x7ff0000000000000)}}) {
    expect_rejected([&] { value.add("rho", 7, 0, invalid_value); });
    HUNDUN_CHECK(value.parts().xor64 == before.xor64);
    HUNDUN_CHECK(value.parts().sum64 == before.sum64);
  }
  DiagnosticFingerprintAccumulator nested;
  nested.add("rho_h", 9U, 1U, describe_fp64(-0.0));
  value.combine(nested.parts());
  HUNDUN_CHECK(value.parts().xor64 == UINT64_C(0xbaedc2ff3f806eae));
  HUNDUN_CHECK(value.parts().sum64 == UINT64_C(0xbcf222ff3fed8f2e));
  HUNDUN_CHECK(value.finish().hex == "baedc2ff3f806eaebcf222ff3fed8f2e");
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

void test_descriptor_request_record_negative_matrix() {
  using namespace hundun::diagnostics;
  const auto descriptor = full_descriptor();
  validate(descriptor);

  for (const std::string_view legal : {"0xdeadbeef", "cuda.stream"}) {
    auto accepted = descriptor;
    accepted.instance_id = legal;
    validate(accepted);
  }
  const std::array<std::string, 6> invalid_descriptor_ids{
      "",
      "Bad",
      "/tmp/value",
      "vendor detail",
      std::string(129U, 'a'),
      std::string("\xc3\x28", 2U)};
  for (const std::string &invalid : invalid_descriptor_ids) {
    auto bad = descriptor;
    bad.instance_id = invalid;
    expect_rejected([&] { validate(bad); });
  }
  {
    auto bad = descriptor;
    bad.schema_version = 2U;
    expect_rejected([&] { validate(bad); });
  }
  {
    auto bad = descriptor;
    bad.capabilities = UINT32_C(1) << 31U;
    expect_rejected([&] { validate(bad); });
  }

  DiagnosticRequest request;
  request.level = DiagnosticLevel::summary;
  request.scope = DiagnosticScope::local;
  request.frame = {0, 3U, -0.0, "material.finalized-trial"};
  validate(request, descriptor);
  for (const auto mutate : {0, 1, 2, 3, 4, 5}) {
    auto bad = request;
    switch (mutate) {
    case 0:
      bad.frame.rank = -1;
      break;
    case 1:
      bad.frame.time_s = -1.0;
      break;
    case 2:
      bad.frame.time_s = std::numeric_limits<double>::infinity();
      break;
    case 3:
      bad.frame.phase = "Bad";
      break;
    case 4:
      bad.sample_budget = 1U;
      break;
    default:
      bad.selected_fields = {"rho"};
      break;
    }
    expect_rejected([&] { validate(bad, descriptor); });
  }
  {
    auto bad = request;
    bad.level = DiagnosticLevel::bounded_state_sample;
    bad.sample_budget = 0U;
    expect_rejected([&] { validate(bad, descriptor); });
    bad.sample_budget = kMaximumStateSamplesV1 + 1U;
    expect_rejected([&] { validate(bad, descriptor); });
    bad.sample_budget = 1U;
    bad.selected_fields = {"rho", "rho"};
    expect_rejected([&] { validate(bad, descriptor); });
    bad.selected_fields = {"rho_h", "rho"};
    expect_rejected([&] { validate(bad, descriptor); });
  }
  {
    auto unsupported = descriptor;
    unsupported.capabilities =
        static_cast<DiagnosticCapabilityFlags>(DiagnosticCapability::summary);
    auto bad = request;
    bad.level = DiagnosticLevel::counters;
    expect_rejected([&] { validate(bad, unsupported); });
    bad = request;
    bad.scope = DiagnosticScope::collective;
    expect_rejected([&] { validate(bad, unsupported); });
  }

  const auto record = summary_record();
  const std::array mutations{
      +[](DiagnosticRecord &value) { value.schema_version = 2U; },
      +[](DiagnosticRecord &value) {
        value.module_kind = DiagnosticModuleKind::runtime;
      },
      +[](DiagnosticRecord &value) { value.module_id = "hundun.flow.other"; },
      +[](DiagnosticRecord &value) { value.instance_id = "secondary"; },
      +[](DiagnosticRecord &value) {
        value.level = DiagnosticLevel::counters;
        value.metrics.clear();
        value.counters.push_back({"owned.cells", "count", 1U});
      },
      +[](DiagnosticRecord &value) {
        value.scope = DiagnosticScope::collective;
      },
      +[](DiagnosticRecord &value) { value.rank = 1; },
      +[](DiagnosticRecord &value) { ++value.step; },
      +[](DiagnosticRecord &value) { value.time_s = describe_fp64(+0.0); },
      +[](DiagnosticRecord &value) { value.phase = "material.other"; }};
  for (const auto mutation : mutations)
    expect_record_mutation_rejected(record, descriptor, request, mutation);

  expect_record_mutation_rejected(
      record, descriptor, request,
      [](DiagnosticRecord &value) { value.metrics.clear(); });
  expect_record_mutation_rejected(
      record, descriptor, request, [](DiagnosticRecord &value) {
        value.invariants.push_back({"density.available", "1",
                                    describe_fp64(1.0), describe_fp64(1.0),
                                    InvariantRelation::equal, true});
      });
  expect_record_mutation_rejected(record, descriptor, request,
                                  [](DiagnosticRecord &value) {
                                    value.sample_budget = 1U;
                                    value.eligible_sample_count = 1U;
                                    value.samples_truncated = true;
                                  });
  expect_record_mutation_rejected(
      record, descriptor, request, [](DiagnosticRecord &value) {
        value.metrics.push_back(value.metrics.front());
      });
  expect_record_mutation_rejected(
      record, descriptor, request, [](DiagnosticRecord &value) {
        value.identities = {
            {"state", std::nullopt, 1U, std::nullopt, std::nullopt},
            {"state", std::nullopt, 2U, std::nullopt, std::nullopt}};
      });

  for (const auto scope :
       {DiagnosticScope::local, DiagnosticScope::collective}) {
    for (const auto status :
         {DiagnosticStatus::ok, DiagnosticStatus::warning,
          DiagnosticStatus::failed, DiagnosticStatus::unavailable}) {
      auto shaped = record;
      auto shaped_request = request;
      shaped.scope = scope;
      shaped_request.scope = scope;
      shaped.status = status;
      if (status == DiagnosticStatus::failed)
        shaped.failure = {DiagnosticFailureClass::conservation,
                          "material.conservation",
                          scope == DiagnosticScope::collective ? 0 : -1};
      else if (status == DiagnosticStatus::unavailable)
        shaped.failure = {DiagnosticFailureClass::unavailable,
                          "material.unavailable",
                          scope == DiagnosticScope::collective ? 0 : -1};
      validate(shaped, descriptor, shaped_request);
      if (status == DiagnosticStatus::ok ||
          status == DiagnosticStatus::warning) {
        expect_record_mutation_rejected(
            shaped, descriptor, shaped_request, [](DiagnosticRecord &value) {
              value.failure.classification =
                  DiagnosticFailureClass::conservation;
            });
        expect_record_mutation_rejected(shaped, descriptor, shaped_request,
                                        [](DiagnosticRecord &value) {
                                          value.failure.code = "material.bad";
                                        });
      } else {
        expect_record_mutation_rejected(
            shaped, descriptor, shaped_request, [](DiagnosticRecord &value) {
              value.failure.classification = DiagnosticFailureClass::none;
            });
        expect_record_mutation_rejected(
            shaped, descriptor, shaped_request,
            [](DiagnosticRecord &value) { value.failure.code = "none"; });
      }
      expect_record_mutation_rejected(
          shaped, descriptor, shaped_request,
          [scope, status](DiagnosticRecord &value) {
            const bool normal = status == DiagnosticStatus::ok ||
                                status == DiagnosticStatus::warning;
            value.failure.lowest_failing_rank =
                scope == DiagnosticScope::local ? 0 : (normal ? 0 : -1);
          });
    }
  }

  DiagnosticRequest sample_request;
  sample_request.level = DiagnosticLevel::bounded_state_sample;
  sample_request.scope = DiagnosticScope::local;
  sample_request.frame = {0, 9U, 0.0, "material.finalized-trial"};
  sample_request.sample_budget = 1U;
  DiagnosticRecord empty_sample;
  empty_sample.module_kind = DiagnosticModuleKind::density_transport;
  empty_sample.module_id = "hundun.flow.material_density_transport";
  empty_sample.instance_id = "primary";
  empty_sample.level = DiagnosticLevel::bounded_state_sample;
  empty_sample.scope = DiagnosticScope::local;
  empty_sample.rank = 0;
  empty_sample.step = 9U;
  empty_sample.time_s = describe_fp64(0.0);
  empty_sample.phase = "material.finalized-trial";
  empty_sample.state_fingerprint = {std::string(kStateFingerprintAlgorithmV1),
                                    "00000000000000000000000000000000"};
  empty_sample.sample_budget = 1U;
  validate(empty_sample, descriptor, sample_request);
  auto nonempty_sample = empty_sample;
  nonempty_sample.eligible_sample_count = 1U;
  nonempty_sample.samples.push_back(
      {"rho", 0U, 0U, "kg/m3", describe_fp64(1.0)});
  validate(nonempty_sample, descriptor, sample_request);
  expect_record_mutation_rejected(
      nonempty_sample, descriptor, sample_request,
      [](DiagnosticRecord &value) { value.sample_budget = 0U; });
  expect_record_mutation_rejected(
      nonempty_sample, descriptor, sample_request,
      [](DiagnosticRecord &value) { value.eligible_sample_count = 0U; });
  expect_record_mutation_rejected(
      nonempty_sample, descriptor, sample_request,
      [](DiagnosticRecord &value) { value.samples_truncated = true; });
  expect_record_mutation_rejected(
      nonempty_sample, descriptor, sample_request,
      [](DiagnosticRecord &value) { value.samples.clear(); });
  auto selected_request = sample_request;
  selected_request.selected_fields = {"rho_h"};
  expect_record_mutation_rejected(nonempty_sample, descriptor, selected_request,
                                  [](DiagnosticRecord &) {});
}

} // namespace

int main() {
  return hundun::test::run([] {
    test_fp64();
    test_invariants();
    test_validation_and_json();
    test_fingerprint();
    test_all_status_level_shapes();
    test_descriptor_request_record_negative_matrix();
  });
}

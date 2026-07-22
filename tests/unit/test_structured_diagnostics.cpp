// SPDX-License-Identifier: Apache-2.0

#include "diagnostics/src/structured_diagnostics_test_access.hpp"
#include "hundun/diagnostics/structured_diagnostics.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
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

struct DiagnosticFixture final {
  hundun::diagnostics::DiagnosticDescriptor descriptor;
  hundun::diagnostics::DiagnosticRequest request;
  hundun::diagnostics::DiagnosticRecord record;
};

DiagnosticFixture make_fixture(hundun::diagnostics::DiagnosticStatus status,
                               hundun::diagnostics::DiagnosticLevel level,
                               hundun::diagnostics::DiagnosticScope scope) {
  using namespace hundun::diagnostics;
  DiagnosticFixture fixture;
  fixture.descriptor = full_descriptor();
  fixture.request.level = level;
  fixture.request.scope = scope;
  fixture.request.frame = {0, 9U, -0.0, "material.finalized-trial"};
  fixture.record.module_kind = DiagnosticModuleKind::density_transport;
  fixture.record.module_id = "hundun.flow.material_density_transport";
  fixture.record.instance_id = "primary";
  fixture.record.level = level;
  fixture.record.scope = scope;
  fixture.record.rank = 0;
  fixture.record.step = 9U;
  fixture.record.time_s = describe_fp64(-0.0);
  fixture.record.phase = "material.finalized-trial";
  fixture.record.status = status;
  if (status == DiagnosticStatus::failed) {
    fixture.record.failure = {DiagnosticFailureClass::conservation,
                              "material.conservation",
                              scope == DiagnosticScope::collective ? 0 : -1};
  } else if (status == DiagnosticStatus::unavailable) {
    fixture.record.failure = {DiagnosticFailureClass::unavailable,
                              "material.unavailable",
                              scope == DiagnosticScope::collective ? 0 : -1};
  }
  fixture.record.identities = {
      {"field.rho", std::string("cell.f64.c1.g2plus.owned.0.0.0.1.1.1"),
       std::nullopt, std::nullopt, std::nullopt},
      {"layout.cells", std::string("cell.f64.c1.g2plus.owned.0.0.0.1.1.1"), 7U,
       std::nullopt, std::nullopt}};
  fixture.record.state_fingerprint = {std::string(kStateFingerprintAlgorithmV1),
                                      "00000000000000000000000000000000"};
  if (level == DiagnosticLevel::summary) {
    fixture.record.metrics = {
        {"density.maximum", DiagnosticMetricKind::state_summary, "kg/m3",
         describe_fp64(1.2)},
        {"density.minimum", DiagnosticMetricKind::state_summary, "kg/m3",
         describe_fp64(0.8)}};
  } else if (level == DiagnosticLevel::invariants) {
    fixture.record.invariants = {{"density.finite",
                                  "kg/m3",
                                  describe_fp64(0.8),
                                  {},
                                  InvariantRelation::finite,
                                  true},
                                 {"density.positive",
                                  "kg/m3",
                                  describe_fp64(0.8),
                                  {},
                                  InvariantRelation::positive,
                                  true}};
  } else if (level == DiagnosticLevel::counters) {
    fixture.record.counters = {{"owned.cells", "count", 1U},
                               {"transport.fields", "count", 2U}};
  } else {
    fixture.request.sample_budget = 2U;
    fixture.record.sample_budget = 2U;
    fixture.record.eligible_sample_count = 3U;
    fixture.record.samples_truncated = true;
    fixture.record.samples = {{"rho", 0U, 0U, "kg/m3", describe_fp64(1.0)},
                              {"rho_h", 0U, 0U, "J/m3", describe_fp64(3.0)}};
  }
  return fixture;
}

template <class Mutation>
void expect_fixture_mutation_rejected(const DiagnosticFixture &baseline,
                                      Mutation &&mutation) {
  static std::size_t mutation_index = 0U;
  const auto current_mutation = mutation_index++;
  auto changed = baseline;
  mutation(changed);
  CountingSink sink;
  bool rejected = false;
  try {
    hundun::diagnostics::validate(changed.record, changed.descriptor,
                                  changed.request);
    sink.submit(changed.record);
  } catch (const hundun::diagnostics::DiagnosticCollectionError &) {
    rejected = true;
  }
  if (!rejected) {
    std::string difference;
    const auto note = [&](bool condition, std::string_view label) {
      if (condition) {
        if (!difference.empty())
          difference += ',';
        difference += label;
      }
    };
    note(changed.descriptor.schema_version !=
             baseline.descriptor.schema_version,
         "descriptor.schema");
    note(changed.descriptor.module_kind != baseline.descriptor.module_kind,
         "descriptor.kind");
    note(changed.descriptor.module_id != baseline.descriptor.module_id,
         "descriptor.module");
    note(changed.descriptor.instance_id != baseline.descriptor.instance_id,
         "descriptor.instance");
    note(changed.descriptor.capabilities != baseline.descriptor.capabilities,
         "descriptor.capabilities");
    note(changed.request.level != baseline.request.level, "request.level");
    note(changed.request.scope != baseline.request.scope, "request.scope");
    note(changed.request.frame.rank != baseline.request.frame.rank,
         "request.rank");
    note(changed.request.frame.step != baseline.request.frame.step,
         "request.step");
    note(changed.request.frame.phase != baseline.request.frame.phase,
         "request.phase");
    note(changed.request.sample_budget != baseline.request.sample_budget,
         "request.budget");
    note(changed.request.selected_fields != baseline.request.selected_fields,
         "request.selection");
    note(changed.record.schema_version != baseline.record.schema_version,
         "record.schema");
    note(changed.record.module_kind != baseline.record.module_kind,
         "record.kind");
    note(changed.record.module_id != baseline.record.module_id,
         "record.module");
    note(changed.record.instance_id != baseline.record.instance_id,
         "record.instance");
    note(changed.record.level != baseline.record.level, "record.level");
    note(changed.record.scope != baseline.record.scope, "record.scope");
    note(changed.record.rank != baseline.record.rank, "record.rank");
    note(changed.record.step != baseline.record.step, "record.step");
    note(changed.record.time_s.status != baseline.record.time_s.status ||
             changed.record.time_s.bits != baseline.record.time_s.bits,
         "record.time");
    note(changed.record.phase != baseline.record.phase, "record.phase");
    note(changed.record.status != baseline.record.status, "record.status");
    note(changed.record.failure.classification !=
             baseline.record.failure.classification,
         "failure.class");
    note(changed.record.failure.code != baseline.record.failure.code,
         "failure.code");
    note(changed.record.failure.lowest_failing_rank !=
             baseline.record.failure.lowest_failing_rank,
         "failure.rank");
    note(changed.record.invariants.size() != baseline.record.invariants.size(),
         "invariants.size");
    note(changed.record.metrics.size() != baseline.record.metrics.size(),
         "metrics.size");
    note(changed.record.counters.size() != baseline.record.counters.size(),
         "counters.size");
    note(changed.record.identities.size() != baseline.record.identities.size(),
         "identities.size");
    note(changed.record.samples.size() != baseline.record.samples.size(),
         "samples.size");
    note(changed.record.sample_budget != baseline.record.sample_budget,
         "record.budget");
    note(changed.record.eligible_sample_count !=
             baseline.record.eligible_sample_count,
         "record.eligible");
    note(changed.record.samples_truncated != baseline.record.samples_truncated,
         "record.truncated");
    if (!changed.record.invariants.empty() &&
        !baseline.record.invariants.empty()) {
      const auto &left = changed.record.invariants.front();
      const auto &right = baseline.record.invariants.front();
      note(left.id != right.id, "invariant.id");
      note(left.observed.status != right.observed.status ||
               left.observed.bits != right.observed.bits,
           "invariant.observed");
      note(left.limit.status != right.limit.status ||
               left.limit.bits != right.limit.bits,
           "invariant.limit");
      note(left.relation != right.relation, "invariant.relation");
      note(left.passed != right.passed, "invariant.passed");
    }
    throw std::runtime_error(
        "diagnostic fixture mutation unexpectedly accepted: level=" +
        std::to_string(static_cast<int>(changed.record.level)) +
        " scope=" + std::to_string(static_cast<int>(changed.record.scope)) +
        " status=" + std::to_string(static_cast<int>(changed.record.status)) +
        " sink_calls=" + std::to_string(sink.calls) + " mutation=" +
        std::to_string(current_mutation) + " difference=" + difference);
  }
  HUNDUN_CHECK(sink.calls == 0U);
}

struct RecordStructuralSnapshot final {
  std::vector<std::uintptr_t> pointers;
  std::vector<std::size_t> sizes;
  std::vector<std::size_t> capacities;
  std::vector<std::string> strings;
  std::vector<std::uint64_t> words;
};

RecordStructuralSnapshot
snapshot_record(const hundun::diagnostics::DiagnosticRecord &record) {
  using namespace hundun::diagnostics;
  RecordStructuralSnapshot snapshot;
  const auto add_string = [&](const std::string &value) {
    snapshot.pointers.push_back(reinterpret_cast<std::uintptr_t>(value.data()));
    snapshot.sizes.push_back(value.size());
    snapshot.capacities.push_back(value.capacity());
    snapshot.strings.push_back(value);
  };
  const auto add_vector = [&](const auto &values) {
    snapshot.pointers.push_back(
        reinterpret_cast<std::uintptr_t>(values.data()));
    snapshot.sizes.push_back(values.size());
    snapshot.capacities.push_back(values.capacity());
  };
  const auto add_optional_u64 = [&](const std::optional<std::uint64_t> &value) {
    snapshot.words.push_back(value.has_value() ? 1U : 0U);
    snapshot.words.push_back(value.value_or(0U));
  };
  const auto add_fp64 = [&](DiagnosticFp64 value) {
    snapshot.words.push_back(static_cast<std::uint64_t>(value.status));
    snapshot.words.push_back(value.bits);
  };

  snapshot.words = {
      record.schema_version,
      static_cast<std::uint64_t>(record.module_kind),
      static_cast<std::uint64_t>(record.level),
      static_cast<std::uint64_t>(record.scope),
      static_cast<std::uint64_t>(static_cast<std::int64_t>(record.rank)),
      record.step,
      static_cast<std::uint64_t>(record.status),
      static_cast<std::uint64_t>(record.failure.classification),
      static_cast<std::uint64_t>(
          static_cast<std::int64_t>(record.failure.lowest_failing_rank)),
      static_cast<std::uint64_t>(record.sample_budget),
      record.eligible_sample_count,
      record.samples_truncated ? 1U : 0U};
  add_fp64(record.time_s);
  add_string(record.module_id);
  add_string(record.instance_id);
  add_string(record.phase);
  add_string(record.failure.code);
  add_string(record.state_fingerprint.algorithm);
  add_string(record.state_fingerprint.hex);

  add_vector(record.invariants);
  for (const auto &value : record.invariants) {
    add_string(value.id);
    add_string(value.unit);
    add_fp64(value.observed);
    add_fp64(value.limit);
    snapshot.words.push_back(static_cast<std::uint64_t>(value.relation));
    snapshot.words.push_back(value.passed ? 1U : 0U);
  }
  add_vector(record.metrics);
  for (const auto &value : record.metrics) {
    add_string(value.id);
    snapshot.words.push_back(static_cast<std::uint64_t>(value.kind));
    add_string(value.unit);
    add_fp64(value.value);
  }
  add_vector(record.counters);
  for (const auto &value : record.counters) {
    add_string(value.id);
    add_string(value.unit);
    snapshot.words.push_back(value.value);
  }
  add_vector(record.identities);
  for (const auto &value : record.identities) {
    add_string(value.subject_id);
    snapshot.words.push_back(value.layout_fingerprint.has_value() ? 1U : 0U);
    if (value.layout_fingerprint)
      add_string(*value.layout_fingerprint);
    add_optional_u64(value.revision);
    add_optional_u64(value.generation);
    add_optional_u64(value.allocation_identity);
  }
  add_vector(record.samples);
  for (const auto &value : record.samples) {
    add_string(value.field_id);
    snapshot.words.push_back(value.global_id);
    snapshot.words.push_back(value.component);
    add_string(value.unit);
    add_fp64(value.value);
  }
  return snapshot;
}

bool same_snapshot(const RecordStructuralSnapshot &left,
                   const RecordStructuralSnapshot &right) {
  return left.pointers == right.pointers && left.sizes == right.sizes &&
         left.capacities == right.capacities && left.strings == right.strings &&
         left.words == right.words;
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
        const auto fixture = make_fixture(status, level, scope);
        validate(fixture.record, fixture.descriptor, fixture.request);
        const auto before = snapshot_record(fixture.record);
        const auto first = to_canonical_json(fixture.record);
        const auto second = to_canonical_json(fixture.record);
        HUNDUN_CHECK(first == second);
        HUNDUN_CHECK(same_snapshot(before, snapshot_record(fixture.record)));

        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.descriptor.schema_version = 2U;
        });
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.descriptor.module_kind = DiagnosticModuleKind::runtime;
        });
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.descriptor.module_id = "hundun.flow.other";
        });
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.descriptor.instance_id = "secondary";
        });
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.descriptor.capabilities = UINT32_C(1) << 31U;
        });
        expect_fixture_mutation_rejected(fixture, [level](auto &changed) {
          changed.descriptor.capabilities &=
              ~static_cast<DiagnosticCapabilityFlags>(
                  level == DiagnosticLevel::summary
                      ? DiagnosticCapability::summary
                  : level == DiagnosticLevel::invariants
                      ? DiagnosticCapability::invariants
                  : level == DiagnosticLevel::counters
                      ? DiagnosticCapability::counters
                      : DiagnosticCapability::bounded_state_sample);
        });
        if (scope == DiagnosticScope::collective) {
          expect_fixture_mutation_rejected(fixture, [](auto &changed) {
            changed.descriptor.capabilities &=
                ~static_cast<DiagnosticCapabilityFlags>(
                    DiagnosticCapability::collective);
          });
        }
        expect_fixture_mutation_rejected(
            fixture, [](auto &changed) { changed.record.schema_version = 2U; });
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.record.module_kind = DiagnosticModuleKind::runtime;
        });
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.record.module_id = "hundun.flow.other";
        });
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.record.instance_id = "secondary";
        });
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.record.level = static_cast<DiagnosticLevel>(255);
        });
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.request.level = static_cast<DiagnosticLevel>(255);
        });
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.record.scope = changed.record.scope == DiagnosticScope::local
                                     ? DiagnosticScope::collective
                                     : DiagnosticScope::local;
        });
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.request.scope =
              changed.request.scope == DiagnosticScope::local
                  ? DiagnosticScope::collective
                  : DiagnosticScope::local;
        });
        expect_fixture_mutation_rejected(
            fixture, [](auto &changed) { ++changed.record.rank; });
        expect_fixture_mutation_rejected(
            fixture, [](auto &changed) { ++changed.request.frame.rank; });
        expect_fixture_mutation_rejected(
            fixture, [](auto &changed) { ++changed.record.step; });
        expect_fixture_mutation_rejected(
            fixture, [](auto &changed) { ++changed.request.frame.step; });
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.record.time_s = describe_fp64(+0.0);
        });
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.request.frame.time_s = +0.0;
        });
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.record.phase = "material.other";
        });
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.request.frame.phase = "material.other";
        });
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.record.status = static_cast<DiagnosticStatus>(255);
        });

        if (level == DiagnosticLevel::summary) {
          expect_fixture_mutation_rejected(
              fixture, [](auto &changed) { changed.record.metrics.clear(); });
        } else if (level == DiagnosticLevel::invariants) {
          expect_fixture_mutation_rejected(fixture, [](auto &changed) {
            changed.record.invariants.clear();
          });
        } else if (level == DiagnosticLevel::counters) {
          expect_fixture_mutation_rejected(
              fixture, [](auto &changed) { changed.record.counters.clear(); });
        } else {
          expect_fixture_mutation_rejected(
              fixture, [](auto &changed) { changed.record.samples.clear(); });
        }
        if (level != DiagnosticLevel::invariants) {
          expect_fixture_mutation_rejected(fixture, [](auto &changed) {
            changed.record.invariants.push_back({"forbidden.invariant",
                                                 "1",
                                                 describe_fp64(1.0),
                                                 {},
                                                 InvariantRelation::finite,
                                                 true});
          });
        }
        if (level != DiagnosticLevel::summary) {
          expect_fixture_mutation_rejected(fixture, [](auto &changed) {
            changed.record.metrics.push_back(
                {"forbidden.metric", DiagnosticMetricKind::state_summary, "1",
                 describe_fp64(1.0)});
          });
        }
        if (level != DiagnosticLevel::counters) {
          expect_fixture_mutation_rejected(fixture, [](auto &changed) {
            changed.record.counters.push_back(
                {"forbidden.counter", "count", 1U});
          });
        }
        if (level != DiagnosticLevel::bounded_state_sample) {
          expect_fixture_mutation_rejected(fixture, [](auto &changed) {
            changed.record.samples.push_back(
                {"rho", 0U, 0U, "kg/m3", describe_fp64(1.0)});
          });
          expect_fixture_mutation_rejected(fixture, [](auto &changed) {
            changed.request.selected_fields = {"rho"};
          });
          expect_fixture_mutation_rejected(fixture, [](auto &changed) {
            changed.request.sample_budget = 1U;
          });
        } else {
          expect_fixture_mutation_rejected(
              fixture, [](auto &changed) { ++changed.record.sample_budget; });
          expect_fixture_mutation_rejected(
              fixture, [](auto &changed) { ++changed.request.sample_budget; });
          expect_fixture_mutation_rejected(fixture, [](auto &changed) {
            changed.record.eligible_sample_count = 1U;
          });
          expect_fixture_mutation_rejected(fixture, [](auto &changed) {
            changed.record.samples_truncated = false;
          });
          expect_fixture_mutation_rejected(fixture, [](auto &changed) {
            changed.record.samples.front().value = {};
          });
          expect_fixture_mutation_rejected(fixture, [](auto &changed) {
            std::reverse(changed.record.samples.begin(),
                         changed.record.samples.end());
          });
          expect_fixture_mutation_rejected(fixture, [](auto &changed) {
            changed.request.selected_fields = {"rho_h"};
          });
        }

        if (status == DiagnosticStatus::ok ||
            status == DiagnosticStatus::warning) {
          expect_fixture_mutation_rejected(fixture, [](auto &changed) {
            changed.record.failure.classification =
                DiagnosticFailureClass::conservation;
          });
          expect_fixture_mutation_rejected(fixture, [](auto &changed) {
            changed.record.failure.code = "material.bad";
          });
        } else {
          expect_fixture_mutation_rejected(fixture, [](auto &changed) {
            changed.record.failure.classification =
                DiagnosticFailureClass::none;
          });
          expect_fixture_mutation_rejected(fixture, [](auto &changed) {
            changed.record.failure.code = "none";
          });
        }
        expect_fixture_mutation_rejected(
            fixture, [scope, status](auto &changed) {
              const bool normal = status == DiagnosticStatus::ok ||
                                  status == DiagnosticStatus::warning;
              changed.record.failure.lowest_failing_rank =
                  scope == DiagnosticScope::local ? 0 : (normal ? 0 : -1);
            });
      }
    }
  }
}

void test_snapshot_helper_mutation_oracle() {
  using namespace hundun::diagnostics;
  auto ordinary = make_fixture(DiagnosticStatus::ok, DiagnosticLevel::summary,
                               DiagnosticScope::local)
                      .record;
  const auto ordinary_before = snapshot_record(ordinary);
  ordinary.module_id.front() = 'x';
  HUNDUN_CHECK(!same_snapshot(ordinary_before, snapshot_record(ordinary)));

  auto nested =
      make_fixture(DiagnosticStatus::ok, DiagnosticLevel::bounded_state_sample,
                   DiagnosticScope::local)
          .record;
  const auto nested_before = snapshot_record(nested);
  nested.samples.back().value.bits ^= 1U;
  HUNDUN_CHECK(!same_snapshot(nested_before, snapshot_record(nested)));
}

void test_all_vector_order_contracts() {
  using namespace hundun::diagnostics;
  const auto status = DiagnosticStatus::ok;
  const auto scope = DiagnosticScope::local;
  for (const auto level :
       {DiagnosticLevel::summary, DiagnosticLevel::invariants,
        DiagnosticLevel::counters, DiagnosticLevel::bounded_state_sample}) {
    const auto fixture = make_fixture(status, level, scope);
    validate(fixture.record, fixture.descriptor, fixture.request);
    if (level == DiagnosticLevel::summary) {
      expect_fixture_mutation_rejected(fixture, [](auto &changed) {
        std::reverse(changed.record.metrics.begin(),
                     changed.record.metrics.end());
      });
      expect_fixture_mutation_rejected(fixture, [](auto &changed) {
        changed.record.metrics.back() = changed.record.metrics.front();
      });
    } else if (level == DiagnosticLevel::invariants) {
      expect_fixture_mutation_rejected(fixture, [](auto &changed) {
        std::reverse(changed.record.invariants.begin(),
                     changed.record.invariants.end());
      });
      expect_fixture_mutation_rejected(fixture, [](auto &changed) {
        changed.record.invariants.back() = changed.record.invariants.front();
      });
    } else if (level == DiagnosticLevel::counters) {
      expect_fixture_mutation_rejected(fixture, [](auto &changed) {
        std::reverse(changed.record.counters.begin(),
                     changed.record.counters.end());
      });
      expect_fixture_mutation_rejected(fixture, [](auto &changed) {
        changed.record.counters.back() = changed.record.counters.front();
      });
    } else {
      expect_fixture_mutation_rejected(fixture, [](auto &changed) {
        std::reverse(changed.record.samples.begin(),
                     changed.record.samples.end());
      });
      expect_fixture_mutation_rejected(fixture, [](auto &changed) {
        changed.record.samples.back() = changed.record.samples.front();
      });
    }
    expect_fixture_mutation_rejected(fixture, [](auto &changed) {
      std::reverse(changed.record.identities.begin(),
                   changed.record.identities.end());
    });
    expect_fixture_mutation_rejected(fixture, [](auto &changed) {
      changed.record.identities.back() = changed.record.identities.front();
    });
  }
}

hundun::diagnostics::DiagnosticInvariant
relation_fixture(hundun::diagnostics::InvariantRelation relation) {
  using namespace hundun::diagnostics;
  switch (relation) {
  case InvariantRelation::less_equal:
    return {"relation.primary", "1",      describe_fp64(1.0),
            describe_fp64(2.0), relation, true};
  case InvariantRelation::greater_equal:
    return {"relation.primary", "1",      describe_fp64(2.0),
            describe_fp64(1.0), relation, true};
  case InvariantRelation::equal:
    return {"relation.primary",  "1",      describe_fp64(+0.0),
            describe_fp64(-0.0), relation, true};
  case InvariantRelation::finite:
  case InvariantRelation::positive:
    return {"relation.primary", "1", describe_fp64(1.0), {}, relation, true};
  }
  throw std::runtime_error("unreachable invariant relation fixture");
}

void test_invariant_relation_status_scope_matrix() {
  using namespace hundun::diagnostics;
  constexpr DiagnosticStatus statuses[]{
      DiagnosticStatus::ok, DiagnosticStatus::warning, DiagnosticStatus::failed,
      DiagnosticStatus::unavailable};
  constexpr DiagnosticScope scopes[]{DiagnosticScope::local,
                                     DiagnosticScope::collective};
  constexpr InvariantRelation relations[]{
      InvariantRelation::less_equal, InvariantRelation::greater_equal,
      InvariantRelation::equal, InvariantRelation::finite,
      InvariantRelation::positive};
  constexpr std::uint64_t nonfinite_bits[]{
      UINT64_C(0x7ff0000000000000), UINT64_C(0xfff0000000000000),
      UINT64_C(0x7ff8000000000000), UINT64_C(0x7ff0000000001234),
      UINT64_C(0x7ff8000000005678)};

  for (const auto status : statuses) {
    for (const auto scope : scopes) {
      for (const auto relation : relations) {
        auto fixture = make_fixture(status, DiagnosticLevel::invariants, scope);
        fixture.record.invariants = {relation_fixture(relation),
                                     {"relation.secondary",
                                      "1",
                                      describe_fp64(1.0),
                                      {},
                                      InvariantRelation::finite,
                                      true}};
        validate(fixture.record, fixture.descriptor, fixture.request);
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.record.invariants.front().observed = {
              DiagnosticValueStatus::finite, UINT64_C(0x7ff0000000000000)};
        });
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.record.invariants.front().observed = {};
        });
        for (const auto bits : nonfinite_bits) {
          expect_fixture_mutation_rejected(fixture, [bits](auto &changed) {
            changed.record.invariants.front().observed =
                describe_fp64(from_bits(bits));
          });
        }
        expect_fixture_mutation_rejected(fixture, [relation](auto &changed) {
          changed.record.invariants.front().limit =
              relation == InvariantRelation::finite ||
                      relation == InvariantRelation::positive
                  ? describe_fp64(1.0)
                  : DiagnosticFp64{};
        });
        for (const auto bits : nonfinite_bits) {
          expect_fixture_mutation_rejected(fixture, [bits](auto &changed) {
            changed.record.invariants.front().limit =
                describe_fp64(from_bits(bits));
          });
        }
        expect_fixture_mutation_rejected(fixture, [relation](auto &changed) {
          auto &invariant = changed.record.invariants.front();
          switch (relation) {
          case InvariantRelation::less_equal:
            invariant.relation = InvariantRelation::greater_equal;
            break;
          case InvariantRelation::greater_equal:
            invariant.relation = InvariantRelation::less_equal;
            break;
          case InvariantRelation::equal:
            invariant.relation = InvariantRelation::positive;
            break;
          case InvariantRelation::finite:
          case InvariantRelation::positive:
            invariant.relation = InvariantRelation::less_equal;
            break;
          }
        });
        expect_fixture_mutation_rejected(fixture, [](auto &changed) {
          changed.record.invariants.front().passed = false;
        });
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
    test_snapshot_helper_mutation_oracle();
    test_all_vector_order_contracts();
    test_invariant_relation_status_scope_matrix();
    test_descriptor_request_record_negative_matrix();
  });
}

// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/diagnostics/structured_diagnostics.hpp"
#include "hundun/flow/checkpoint_v2.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <optional>
#include <string>
#include <string_view>

namespace hundun::test::checkpoint_v2_oracle {

inline std::string_view phase_name(flow::CheckpointV2Phase phase) {
  constexpr std::array names{
      std::string_view{"none"},           std::string_view{"preflight"},
      std::string_view{"transaction-entry"},
      std::string_view{"rank-payload"},
      std::string_view{"rank-temporary-file"},
      std::string_view{"rank-publish"},   std::string_view{"manifest"},
      std::string_view{"completed-marker"},
      std::string_view{"marker-read"},    std::string_view{"manifest-read"},
      std::string_view{"rank-read"},      std::string_view{"restore-prepare"},
      std::string_view{"restore-publish"}};
  const auto index = static_cast<std::size_t>(phase);
  HUNDUN_CHECK(index < names.size());
  return names[index];
}

inline diagnostics::DiagnosticFailureClass
failure_class(flow::CheckpointV2FailureReason reason) {
  using Reason = flow::CheckpointV2FailureReason;
  switch (reason) {
  case Reason::none:
    return diagnostics::DiagnosticFailureClass::none;
  case Reason::invalid_input:
  case Reason::state:
    return diagnostics::DiagnosticFailureClass::invalid_input;
  case Reason::layout:
    return diagnostics::DiagnosticFailureClass::layout;
  case Reason::file_integrity:
  case Reason::filesystem:
    return diagnostics::DiagnosticFailureClass::file_integrity;
  }
  HUNDUN_CHECK(false);
  return diagnostics::DiagnosticFailureClass::invalid_input;
}

inline std::string failure_code(const flow::CheckpointV2Report &report) {
  using Reason = flow::CheckpointV2FailureReason;
  if (report.reason() == Reason::none)
    return "none";
  constexpr std::array reason_names{
      std::string_view{"none"}, std::string_view{"invalid-input"},
      std::string_view{"layout"}, std::string_view{"state"},
      std::string_view{"file-integrity"}, std::string_view{"filesystem"}};
  const auto index = static_cast<std::size_t>(report.reason());
  HUNDUN_CHECK(index < reason_names.size());
  return std::string("checkpoint-v2.") +
         (report.operation() == flow::CheckpointV2Operation::write ? "write"
                                                                   : "read") +
         "." + std::string(phase_name(report.phase())) + "." +
         std::string(reason_names[index]);
}

inline diagnostics::DiagnosticStateFingerprint
state_fingerprint(const flow::CheckpointV2Report &report) {
  constexpr std::array<std::string_view, 26> ids{
      "checkpoint.collective-count",
      "checkpoint.crc-check-count",
      "checkpoint.disposition",
      "checkpoint.exact-size-eof-status",
      "checkpoint.file-count",
      "checkpoint.fingerprint-status",
      "checkpoint.global-actual-bytes",
      "checkpoint.global-logical-bytes",
      "checkpoint.local-actual-bytes",
      "checkpoint.local-crc64",
      "checkpoint.local-logical-bytes",
      "checkpoint.lowest-failing-rank",
      "checkpoint.manifest-crc-status",
      "checkpoint.manifest-crc64",
      "checkpoint.operation",
      "checkpoint.partition-status",
      "checkpoint.phase",
      "checkpoint.publication-status",
      "checkpoint.rank",
      "checkpoint.rank-crc-status",
      "checkpoint.reason",
      "checkpoint.rollback-status",
      "checkpoint.semantic-fingerprint",
      "checkpoint.step",
      "checkpoint.time",
      "checkpoint.transaction-entry-status"};
  diagnostics::DiagnosticFingerprintAccumulator accumulator;
  const auto limb = [&](std::size_t index, std::uint64_t value) {
    accumulator.add(
        ids[index], static_cast<std::uint64_t>(report.rank()), 0U,
        diagnostics::describe_fp64(
            static_cast<double>(static_cast<std::uint32_t>(value))));
    accumulator.add(
        ids[index], static_cast<std::uint64_t>(report.rank()), 1U,
        diagnostics::describe_fp64(
            static_cast<double>(static_cast<std::uint32_t>(value >> 32U))));
  };
  const auto count = [&](std::size_t index, std::int64_t value) {
    accumulator.add(ids[index], static_cast<std::uint64_t>(report.rank()), 0U,
                    diagnostics::describe_fp64(static_cast<double>(value)));
  };
  limb(0U, report.collective_count());
  limb(1U, report.crc_check_count());
  count(2U, static_cast<std::uint8_t>(report.disposition()));
  count(3U, static_cast<std::uint8_t>(report.exact_size_and_eof_status()));
  limb(4U, report.file_count());
  count(5U, static_cast<std::uint8_t>(report.fingerprint_status()));
  limb(6U, report.global_actual_bytes());
  limb(7U, report.global_logical_bytes());
  limb(8U, report.local_actual_bytes());
  limb(9U, report.local_crc64());
  limb(10U, report.local_logical_bytes());
  count(11U, report.lowest_failing_rank());
  count(12U, static_cast<std::uint8_t>(report.manifest_crc_status()));
  limb(13U, report.manifest_crc64());
  count(14U, static_cast<std::uint8_t>(report.operation()));
  count(15U, static_cast<std::uint8_t>(report.partition_status()));
  count(16U, static_cast<std::uint8_t>(report.phase()));
  count(17U, static_cast<std::uint8_t>(report.publication_status()));
  count(18U, report.rank());
  count(19U, static_cast<std::uint8_t>(report.rank_crc_status()));
  count(20U, static_cast<std::uint8_t>(report.reason()));
  count(21U, static_cast<std::uint8_t>(report.rollback_status()));
  limb(22U, report.semantic_fingerprint());
  limb(23U, report.step());
  accumulator.add(ids[24U], static_cast<std::uint64_t>(report.rank()), 0U,
                  diagnostics::describe_fp64(report.time_s()));
  count(25U, static_cast<std::uint8_t>(
                 report.transaction_entry_status()));
  return accumulator.finish();
}

inline void require_fp64(double expected,
                         const diagnostics::DiagnosticFp64 &actual) {
  HUNDUN_CHECK(actual.status == diagnostics::DiagnosticValueStatus::finite);
  std::uint64_t bits{};
  std::memcpy(&bits, &expected, sizeof(bits));
  HUNDUN_CHECK(actual.bits == bits);
}

inline void require_exact_product_diagnostic_record(
    const flow::CheckpointV2Report &report, diagnostics::DiagnosticLevel level,
    const diagnostics::DiagnosticRecord &record) {
  HUNDUN_CHECK(record.schema_version == 1U);
  HUNDUN_CHECK(record.module_kind ==
               diagnostics::DiagnosticModuleKind::checkpoint);
  HUNDUN_CHECK(record.module_id == "checkpoint-v2");
  HUNDUN_CHECK(record.instance_id == "checkpoint-v2");
  HUNDUN_CHECK(record.level == level);
  HUNDUN_CHECK(record.scope == diagnostics::DiagnosticScope::local);
  HUNDUN_CHECK(record.rank == report.rank());
  HUNDUN_CHECK(record.step == report.step());
  require_fp64(report.time_s(), record.time_s);
  HUNDUN_CHECK(record.phase == phase_name(report.phase()));
  HUNDUN_CHECK(
      record.status ==
      (report.disposition() == flow::CheckpointV2Disposition::completed
           ? diagnostics::DiagnosticStatus::ok
           : diagnostics::DiagnosticStatus::failed));
  HUNDUN_CHECK(record.failure.classification ==
               failure_class(report.reason()));
  HUNDUN_CHECK(record.failure.code == failure_code(report));
  HUNDUN_CHECK(record.failure.lowest_failing_rank == -1);
  HUNDUN_CHECK(record.identities.size() == 1U);
  HUNDUN_CHECK(record.identities.front().subject_id == "checkpoint-v2");
  std::array<char, 17> identity{};
  std::snprintf(identity.data(), identity.size(), "%016llx",
                static_cast<unsigned long long>(
                    report.semantic_fingerprint()));
  HUNDUN_CHECK(record.identities.front().layout_fingerprint ==
               std::optional<std::string>{identity.data()});
  HUNDUN_CHECK(!record.identities.front().revision);
  HUNDUN_CHECK(!record.identities.front().generation);
  HUNDUN_CHECK(!record.identities.front().allocation_identity);
  HUNDUN_CHECK(record.state_fingerprint.hex ==
               state_fingerprint(report).hex);
  HUNDUN_CHECK(record.state_fingerprint.algorithm ==
               diagnostics::kStateFingerprintAlgorithmV1);
  HUNDUN_CHECK(record.sample_budget == 0U);
  HUNDUN_CHECK(record.eligible_sample_count == 0U);
  HUNDUN_CHECK(!record.samples_truncated);
  HUNDUN_CHECK(record.samples.empty());

  if (level == diagnostics::DiagnosticLevel::summary) {
    struct Expected final {
      std::string_view id;
      std::string_view unit;
      double value;
    };
    const std::array expected{
        Expected{"checkpoint.disposition", "count",
                 static_cast<double>(
                     static_cast<std::uint8_t>(report.disposition()))},
        Expected{"checkpoint.global-actual-bytes-high", "byte",
                 static_cast<double>(static_cast<std::uint32_t>(
                     report.global_actual_bytes() >> 32U))},
        Expected{"checkpoint.global-actual-bytes-low", "byte",
                 static_cast<double>(static_cast<std::uint32_t>(
                     report.global_actual_bytes()))},
        Expected{"checkpoint.global-logical-bytes-high", "byte",
                 static_cast<double>(static_cast<std::uint32_t>(
                     report.global_logical_bytes() >> 32U))},
        Expected{"checkpoint.global-logical-bytes-low", "byte",
                 static_cast<double>(static_cast<std::uint32_t>(
                     report.global_logical_bytes()))},
        Expected{"checkpoint.local-actual-bytes-high", "byte",
                 static_cast<double>(static_cast<std::uint32_t>(
                     report.local_actual_bytes() >> 32U))},
        Expected{"checkpoint.local-actual-bytes-low", "byte",
                 static_cast<double>(static_cast<std::uint32_t>(
                     report.local_actual_bytes()))},
        Expected{"checkpoint.local-logical-bytes-high", "byte",
                 static_cast<double>(static_cast<std::uint32_t>(
                     report.local_logical_bytes() >> 32U))},
        Expected{"checkpoint.local-logical-bytes-low", "byte",
                 static_cast<double>(static_cast<std::uint32_t>(
                     report.local_logical_bytes()))},
        Expected{"checkpoint.operation", "count",
                 static_cast<double>(
                     static_cast<std::uint8_t>(report.operation()))},
        Expected{"checkpoint.reason", "count",
                 static_cast<double>(
                     static_cast<std::uint8_t>(report.reason()))}};
    HUNDUN_CHECK(record.metrics.size() == expected.size());
    HUNDUN_CHECK(record.invariants.empty() && record.counters.empty());
    for (std::size_t index = 0; index < expected.size(); ++index) {
      HUNDUN_CHECK(record.metrics[index].id == expected[index].id);
      HUNDUN_CHECK(record.metrics[index].kind ==
                   diagnostics::DiagnosticMetricKind::state_summary);
      HUNDUN_CHECK(record.metrics[index].unit == expected[index].unit);
      require_fp64(expected[index].value, record.metrics[index].value);
    }
    return;
  }
  if (level == diagnostics::DiagnosticLevel::invariants) {
    constexpr std::array<std::string_view, 8> ids{
        "checkpoint.exact-size-eof", "checkpoint.fingerprint",
        "checkpoint.manifest-crc",   "checkpoint.partition",
        "checkpoint.publication",    "checkpoint.rank-crc",
        "checkpoint.rollback",       "checkpoint.transaction-entry"};
    const std::array statuses{
        report.exact_size_and_eof_status(), report.fingerprint_status(),
        report.manifest_crc_status(),       report.partition_status(),
        report.publication_status(),        report.rank_crc_status(),
        report.rollback_status(),           report.transaction_entry_status()};
    HUNDUN_CHECK(record.invariants.size() == ids.size());
    HUNDUN_CHECK(record.metrics.empty() && record.counters.empty());
    for (std::size_t index = 0; index < ids.size(); ++index) {
      const auto observed =
          statuses[index] == flow::CheckpointV2CheckStatus::passed
              ? 1.0
              : statuses[index] == flow::CheckpointV2CheckStatus::failed
                    ? 0.0
                    : -1.0;
      HUNDUN_CHECK(record.invariants[index].id == ids[index]);
      HUNDUN_CHECK(record.invariants[index].unit == "1");
      require_fp64(observed, record.invariants[index].observed);
      require_fp64(1.0, record.invariants[index].limit);
      HUNDUN_CHECK(record.invariants[index].relation ==
                   diagnostics::InvariantRelation::equal);
      HUNDUN_CHECK(record.invariants[index].passed ==
                   (statuses[index] ==
                    flow::CheckpointV2CheckStatus::passed));
    }
    return;
  }
  struct ExpectedCounter final {
    std::string_view id;
    std::string_view unit;
    std::uint64_t value;
  };
  const std::array expected{
      ExpectedCounter{"checkpoint.collective-count", "count",
                      report.collective_count()},
      ExpectedCounter{"checkpoint.crc-check-count", "count",
                      report.crc_check_count()},
      ExpectedCounter{"checkpoint.file-count", "count", report.file_count()},
      ExpectedCounter{"checkpoint.global-actual-bytes", "byte",
                      report.global_actual_bytes()},
      ExpectedCounter{"checkpoint.global-logical-bytes", "byte",
                      report.global_logical_bytes()},
      ExpectedCounter{"checkpoint.local-actual-bytes", "byte",
                      report.local_actual_bytes()},
      ExpectedCounter{"checkpoint.local-crc64-high", "count",
                      report.local_crc64() >> 32U},
      ExpectedCounter{"checkpoint.local-crc64-low", "count",
                      static_cast<std::uint32_t>(report.local_crc64())},
      ExpectedCounter{"checkpoint.local-logical-bytes", "byte",
                      report.local_logical_bytes()},
      ExpectedCounter{"checkpoint.manifest-crc64-high", "count",
                      report.manifest_crc64() >> 32U},
      ExpectedCounter{"checkpoint.manifest-crc64-low", "count",
                      static_cast<std::uint32_t>(report.manifest_crc64())}};
  HUNDUN_CHECK(record.counters.size() == expected.size());
  HUNDUN_CHECK(record.metrics.empty() && record.invariants.empty());
  for (std::size_t index = 0; index < expected.size(); ++index) {
    HUNDUN_CHECK(record.counters[index].id == expected[index].id);
    HUNDUN_CHECK(record.counters[index].unit == expected[index].unit);
    HUNDUN_CHECK(record.counters[index].value == expected[index].value);
  }
}

} // namespace hundun::test::checkpoint_v2_oracle

// SPDX-License-Identifier: Apache-2.0

#include "flow/src/checkpoint_v2_detail.hpp"
#include "hundun/diagnostics/checkpoint_v2_diagnostics.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"
#include "tests/support/test_main.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <string>
#include <string_view>
#include <vector>

namespace {

class OneRecordSink final : public hundun::diagnostics::DiagnosticSink {
public:
  void submit(const hundun::diagnostics::DiagnosticRecord &value) override {
    record = value;
    ++calls;
  }
  hundun::diagnostics::DiagnosticRecord record;
  int calls{};
};

class DiagnosticErrorSink final : public hundun::diagnostics::DiagnosticSink {
public:
  void submit(const hundun::diagnostics::DiagnosticRecord &) override {
    throw hundun::diagnostics::DiagnosticCollectionError(
        hundun::diagnostics::DiagnosticFailureClass::invalid_request,
        "test.sink.original", -1, "deliberate sink rejection");
  }
};

template <class Function> bool rejects(Function &&function) {
  try {
    function();
  } catch (const hundun::diagnostics::DiagnosticCollectionError &) {
    return true;
  }
  return false;
}

hundun::flow::CheckpointV2Report report(int rank) {
  hundun::flow::detail::CheckpointV2ReportValues values;
  values.operation = hundun::flow::CheckpointV2Operation::write;
  values.disposition = hundun::flow::CheckpointV2Disposition::completed;
  values.reason = hundun::flow::CheckpointV2FailureReason::none;
  values.phase = hundun::flow::CheckpointV2Phase::completed_marker;
  values.rank = rank;
  values.step = 4U;
  values.time_s = 0.125;
  values.local_logical_bytes = 100U + static_cast<std::uint64_t>(rank);
  values.local_actual_bytes = 200U + static_cast<std::uint64_t>(rank);
  values.global_logical_bytes = 404U;
  values.global_actual_bytes = 900U;
  values.local_crc64 = 0x12340000U + static_cast<std::uint64_t>(rank);
  values.manifest_crc64 = 0x5678U;
  values.file_count = 6U;
  values.crc_check_count = 6U;
  values.collective_count = 12U;
  values.rank_crc = hundun::flow::CheckpointV2CheckStatus::passed;
  values.manifest_crc = hundun::flow::CheckpointV2CheckStatus::passed;
  values.exact_size_eof = hundun::flow::CheckpointV2CheckStatus::passed;
  values.fingerprint = hundun::flow::CheckpointV2CheckStatus::passed;
  values.partition = hundun::flow::CheckpointV2CheckStatus::passed;
  values.publication = hundun::flow::CheckpointV2CheckStatus::passed;
  return hundun::flow::detail::CheckpointV2Access::make(values);
}

hundun::diagnostics::DiagnosticRequest
request(int rank, hundun::diagnostics::DiagnosticLevel level,
        hundun::diagnostics::DiagnosticScope scope,
        std::string_view phase = "completed-marker") {
  return {level, scope, {rank, 4U, 0.125, phase}, {}, 0U};
}

std::string_view phase_name(hundun::flow::CheckpointV2Phase phase) {
  using Phase = hundun::flow::CheckpointV2Phase;
  switch (phase) {
  case Phase::preflight:
    return "preflight";
  case Phase::manifest:
    return "manifest";
  case Phase::rank_read:
    return "rank-read";
  case Phase::restore_prepare:
    return "restore-prepare";
  default:
    return "none";
  }
}

hundun::flow::CheckpointV2Report
failure_report(int rank, hundun::flow::CheckpointV2FailureReason reason,
               hundun::flow::CheckpointV2Phase phase) {
  hundun::flow::detail::CheckpointV2ReportValues values;
  values.operation =
      phase == hundun::flow::CheckpointV2Phase::rank_read ||
              phase == hundun::flow::CheckpointV2Phase::restore_prepare
          ? hundun::flow::CheckpointV2Operation::read
          : hundun::flow::CheckpointV2Operation::write;
  values.disposition = hundun::flow::CheckpointV2Disposition::failed;
  values.reason = reason;
  values.phase = phase;
  values.rank = rank;
  values.lowest_failing_rank = 0;
  values.step = 4U;
  values.time_s = 0.125;
  values.local_logical_bytes = 10U + static_cast<std::uint64_t>(rank);
  values.local_actual_bytes = 20U + static_cast<std::uint64_t>(rank);
  values.local_crc64 =
      UINT64_C(0x1020304050607000) + static_cast<std::uint64_t>(rank);
  values.manifest_crc = hundun::flow::CheckpointV2CheckStatus::not_checked;
  values.exact_size_eof = hundun::flow::CheckpointV2CheckStatus::failed;
  values.fingerprint = hundun::flow::CheckpointV2CheckStatus::not_checked;
  values.partition = hundun::flow::CheckpointV2CheckStatus::passed;
  values.transaction_entry = hundun::flow::CheckpointV2CheckStatus::passed;
  values.publication = hundun::flow::CheckpointV2CheckStatus::not_checked;
  values.rollback = hundun::flow::CheckpointV2CheckStatus::passed;
  values.rank_crc = rank == 0
                        ? hundun::flow::CheckpointV2CheckStatus::failed
                        : hundun::flow::CheckpointV2CheckStatus::not_checked;
  return hundun::flow::detail::CheckpointV2Access::make(values);
}

std::uint64_t independent_rank_crc_digest(int ranks) {
  std::vector<std::uint8_t> bytes;
  constexpr std::string_view domain =
      "hundun.checkpoint-v2.diagnostic-rank-crcs.v1";
  bytes.insert(bytes.end(), domain.begin(), domain.end());
  bytes.push_back(0U);
  const auto u32 = [&](std::uint32_t value) {
    for (unsigned shift = 0U; shift < 32U; shift += 8U)
      bytes.push_back(static_cast<std::uint8_t>(value >> shift));
  };
  const auto i32 = [&](std::int32_t value) {
    std::uint32_t wire{};
    std::memcpy(&wire, &value, sizeof(wire));
    u32(wire);
  };
  const auto u64 = [&](std::uint64_t value) {
    for (unsigned shift = 0U; shift < 64U; shift += 8U)
      bytes.push_back(static_cast<std::uint8_t>(value >> shift));
  };
  u32(1U);
  u32(static_cast<std::uint32_t>(ranks));
  for (int rank = 0; rank < ranks; ++rank) {
    i32(rank);
    u64(UINT64_C(0x1020304050607000) + static_cast<std::uint64_t>(rank));
  }
  constexpr std::uint64_t polynomial = UINT64_C(0x42F0E1EBA9EA3693);
  std::uint64_t crc{};
  for (const auto byte : bytes) {
    crc ^= static_cast<std::uint64_t>(byte) << 56U;
    for (unsigned bit = 0; bit < 8U; ++bit)
      crc = (crc & (UINT64_C(1) << 63U)) != 0U ? (crc << 1U) ^ polynomial
                                               : crc << 1U;
  }
  return crc;
}

double finite_value(const hundun::diagnostics::DiagnosticFp64 &value) {
  HUNDUN_CHECK(value.status ==
               hundun::diagnostics::DiagnosticValueStatus::finite);
  double result{};
  std::memcpy(&result, &value.bits, sizeof(result));
  return result;
}

hundun::diagnostics::DiagnosticStateFingerprint
expected_local_fingerprint(const hundun::flow::CheckpointV2Report &value) {
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
  hundun::diagnostics::DiagnosticFingerprintAccumulator accumulator;
  const auto limb = [&](std::size_t index, std::uint64_t item) {
    accumulator.add(ids[index], static_cast<std::uint64_t>(value.rank()), 0U,
        hundun::diagnostics::describe_fp64(
            static_cast<double>(static_cast<std::uint32_t>(item))));
    accumulator.add(ids[index], static_cast<std::uint64_t>(value.rank()), 1U,
                    hundun::diagnostics::describe_fp64(static_cast<double>(
                        static_cast<std::uint32_t>(item >> 32U))));
  };
  const auto count = [&](std::size_t index, std::int64_t item) {
    accumulator.add(
        ids[index], static_cast<std::uint64_t>(value.rank()), 0U,
        hundun::diagnostics::describe_fp64(static_cast<double>(item)));
  };
  limb(0U, value.collective_count());
  limb(1U, value.crc_check_count());
  count(2U, static_cast<std::uint8_t>(value.disposition()));
  count(3U, static_cast<std::uint8_t>(value.exact_size_and_eof_status()));
  limb(4U, value.file_count());
  count(5U, static_cast<std::uint8_t>(value.fingerprint_status()));
  limb(6U, value.global_actual_bytes());
  limb(7U, value.global_logical_bytes());
  limb(8U, value.local_actual_bytes());
  limb(9U, value.local_crc64());
  limb(10U, value.local_logical_bytes());
  count(11U, value.lowest_failing_rank());
  count(12U, static_cast<std::uint8_t>(value.manifest_crc_status()));
  limb(13U, value.manifest_crc64());
  count(14U, static_cast<std::uint8_t>(value.operation()));
  count(15U, static_cast<std::uint8_t>(value.partition_status()));
  count(16U, static_cast<std::uint8_t>(value.phase()));
  count(17U, static_cast<std::uint8_t>(value.publication_status()));
  count(18U, value.rank());
  count(19U, static_cast<std::uint8_t>(value.rank_crc_status()));
  count(20U, static_cast<std::uint8_t>(value.reason()));
  count(21U, static_cast<std::uint8_t>(value.rollback_status()));
  limb(22U, value.semantic_fingerprint());
  limb(23U, value.step());
  accumulator.add(ids[24U], static_cast<std::uint64_t>(value.rank()), 0U,
                  hundun::diagnostics::describe_fp64(value.time_s()));
  count(25U, static_cast<std::uint8_t>(value.transaction_entry_status()));
  return accumulator.finish();
}

void run(const hundun::runtime::MpiContext &mpi) {
  const auto local_report = report(mpi.rank());
  const auto expected_fingerprint = expected_local_fingerprint(local_report);
  auto source = hundun::flow::checkpoint_v2_diagnostic_source(local_report);
  const auto descriptor = hundun::diagnostics::describe_diagnostics(source);
  HUNDUN_CHECK(descriptor.module_kind ==
               hundun::diagnostics::DiagnosticModuleKind::checkpoint);
  HUNDUN_CHECK(!hundun::diagnostics::has_capability(
      descriptor.capabilities,
      hundun::diagnostics::DiagnosticCapability::bounded_state_sample));
  const auto ids =
      hundun::diagnostics::diagnostic_fingerprint_field_ids(source);
  HUNDUN_CHECK(ids.size() == 26U);
  HUNDUN_CHECK(std::is_sorted(ids.begin(), ids.end()));

  for (const auto level : {hundun::diagnostics::DiagnosticLevel::summary,
           hundun::diagnostics::DiagnosticLevel::invariants,
           hundun::diagnostics::DiagnosticLevel::counters}) {
    OneRecordSink first;
    const auto local =
        request(mpi.rank(), level, hundun::diagnostics::DiagnosticScope::local);
    hundun::diagnostics::collect_diagnostics(source, local, first);
    HUNDUN_CHECK(first.calls == 1);
    HUNDUN_CHECK(first.record.scope ==
                 hundun::diagnostics::DiagnosticScope::local);
    HUNDUN_CHECK(first.record.state_fingerprint.hex ==
                 expected_fingerprint.hex);
    const auto first_json =
        hundun::diagnostics::to_canonical_json(first.record);
    OneRecordSink second;
    hundun::diagnostics::collect_diagnostics(source, local, second);
    HUNDUN_CHECK(first_json ==
                 hundun::diagnostics::to_canonical_json(second.record));

    OneRecordSink collective;
    const auto global = request(
        mpi.rank(), level, hundun::diagnostics::DiagnosticScope::collective);
    hundun::diagnostics::collect_diagnostics(source, mpi, global, collective);
    HUNDUN_CHECK(collective.calls == 1);
    HUNDUN_CHECK(collective.record.scope ==
                 hundun::diagnostics::DiagnosticScope::collective);
    const auto expected_logical =
        static_cast<std::uint64_t>(mpi.size()) * 100U +
        static_cast<std::uint64_t>(mpi.size() * (mpi.size() - 1) / 2);
    const auto expected_actual =
        static_cast<std::uint64_t>(mpi.size()) * 200U +
        static_cast<std::uint64_t>(mpi.size() * (mpi.size() - 1) / 2);
    if (level == hundun::diagnostics::DiagnosticLevel::summary) {
      const auto find_metric = [&](std::string_view id) {
        return std::find_if(collective.record.metrics.begin(),
                            collective.record.metrics.end(),
                            [&](const auto &item) { return item.id == id; });
      };
      const auto logical = find_metric("checkpoint.local-logical-bytes-low");
      const auto actual = find_metric("checkpoint.local-actual-bytes-low");
      HUNDUN_CHECK(logical != collective.record.metrics.end());
      HUNDUN_CHECK(actual != collective.record.metrics.end());
      HUNDUN_CHECK(finite_value(logical->value) ==
                   static_cast<double>(expected_logical));
      HUNDUN_CHECK(finite_value(actual->value) ==
                   static_cast<double>(expected_actual));
    } else if (level == hundun::diagnostics::DiagnosticLevel::invariants) {
      const auto rank_crc = std::find_if(
          collective.record.invariants.begin(),
          collective.record.invariants.end(),
          [](const auto &item) { return item.id == "checkpoint.rank-crc"; });
      HUNDUN_CHECK(rank_crc != collective.record.invariants.end());
      HUNDUN_CHECK(rank_crc->passed);
    } else {
      const auto find_counter = [&](std::string_view id) {
        return std::find_if(collective.record.counters.begin(),
                            collective.record.counters.end(),
                            [&](const auto &item) { return item.id == id; });
      };
      const auto logical = find_counter("checkpoint.local-logical-bytes");
      const auto actual = find_counter("checkpoint.local-actual-bytes");
      HUNDUN_CHECK(logical != collective.record.counters.end());
      HUNDUN_CHECK(actual != collective.record.counters.end());
      HUNDUN_CHECK(logical->value == expected_logical);
      HUNDUN_CHECK(actual->value == expected_actual);
    }
    std::array<char, 33> root{};
    if (mpi.rank() == 0)
      std::copy(collective.record.state_fingerprint.hex.begin(),
                collective.record.state_fingerprint.hex.end(), root.begin());
    hundun::runtime::check_mpi_result(
        MPI_Bcast(root.data(), static_cast<int>(root.size()), MPI_CHAR, 0,
                  mpi.comm()),
        "MPI_Bcast(Checkpoint diagnostics fingerprint)");
    HUNDUN_CHECK(collective.record.state_fingerprint.hex ==
                 std::string(root.data()));
  }

  if (mpi.size() > 1) {
    auto disagreement =
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                hundun::diagnostics::DiagnosticScope::collective);
    if (mpi.rank() == mpi.size() - 1)
      ++disagreement.frame.step;
    bool rejected = false;
    try {
      OneRecordSink rejected_sink;
      hundun::diagnostics::collect_diagnostics(source, mpi, disagreement,
                                               rejected_sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      rejected = error.classification() ==
              hundun::diagnostics::DiagnosticFailureClass::
                  collective_operation &&
          error.lowest_failing_rank() == mpi.size() - 1;
    }
    HUNDUN_CHECK(rejected);
  }

  auto moved = std::move(source);
  HUNDUN_CHECK(rejects([&] {
    OneRecordSink sink;
    hundun::diagnostics::collect_diagnostics(
        source,
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                hundun::diagnostics::DiagnosticScope::local),
        sink);
  }));
  OneRecordSink sink;
  hundun::diagnostics::collect_diagnostics(
      moved,
      request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
              hundun::diagnostics::DiagnosticScope::local),
      sink);
  HUNDUN_CHECK(sink.calls == 1);

  {
    auto wrong_scope =
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                hundun::diagnostics::DiagnosticScope::collective);
    try {
      OneRecordSink rejected_sink;
      hundun::diagnostics::collect_diagnostics(moved, wrong_scope,
                                               rejected_sink);
      HUNDUN_CHECK(false);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      HUNDUN_CHECK(error.classification() ==
          hundun::diagnostics::DiagnosticFailureClass::capability);
    }
  }
  {
    DiagnosticErrorSink rejected_sink;
    try {
      hundun::diagnostics::collect_diagnostics(
          moved,
          request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::summary,
                  hundun::diagnostics::DiagnosticScope::local),
          rejected_sink);
      HUNDUN_CHECK(false);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      HUNDUN_CHECK(error.classification() ==
          hundun::diagnostics::DiagnosticFailureClass::sink_failure);
      HUNDUN_CHECK(error.code() == "diagnostics.sink.submit");
    }
  }

  using Reason = hundun::flow::CheckpointV2FailureReason;
  using Phase = hundun::flow::CheckpointV2Phase;
  struct FailureCase final {
    Reason reason;
    Phase phase;
    hundun::diagnostics::DiagnosticFailureClass classification;
    std::string_view suffix;
  };
  constexpr std::array failure_cases{
      FailureCase{Reason::invalid_input, Phase::preflight,
                  hundun::diagnostics::DiagnosticFailureClass::invalid_input,
                  "write.preflight.invalid-input"},
      FailureCase{Reason::layout, Phase::preflight,
                  hundun::diagnostics::DiagnosticFailureClass::layout,
                  "write.preflight.layout"},
      FailureCase{Reason::state, Phase::restore_prepare,
                  hundun::diagnostics::DiagnosticFailureClass::invalid_input,
                  "read.restore-prepare.state"},
      FailureCase{Reason::file_integrity, Phase::rank_read,
                  hundun::diagnostics::DiagnosticFailureClass::file_integrity,
                  "read.rank-read.file-integrity"},
      FailureCase{Reason::filesystem, Phase::manifest,
                  hundun::diagnostics::DiagnosticFailureClass::file_integrity,
                  "write.manifest.filesystem"}};
  for (const auto &item : failure_cases) {
    const auto failed = failure_report(mpi.rank(), item.reason, item.phase);
    auto failed_source = hundun::flow::checkpoint_v2_diagnostic_source(failed);
    OneRecordSink local_failure;
    hundun::diagnostics::collect_diagnostics(
        failed_source,
        request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::invariants,
                hundun::diagnostics::DiagnosticScope::local,
                phase_name(item.phase)),
        local_failure);
    HUNDUN_CHECK(local_failure.record.status ==
                 hundun::diagnostics::DiagnosticStatus::failed);
    HUNDUN_CHECK(local_failure.record.failure.classification ==
                 item.classification);
    HUNDUN_CHECK(local_failure.record.failure.code ==
                 std::string("checkpoint-v2.") + std::string(item.suffix));
    const auto exact = std::find_if(
        local_failure.record.invariants.begin(),
        local_failure.record.invariants.end(), [](const auto &invariant) {
          return invariant.id == "checkpoint.exact-size-eof";
        });
    const auto manifest = std::find_if(
        local_failure.record.invariants.begin(),
        local_failure.record.invariants.end(), [](const auto &invariant) {
          return invariant.id == "checkpoint.manifest-crc";
        });
    HUNDUN_CHECK(exact != local_failure.record.invariants.end());
    HUNDUN_CHECK(manifest != local_failure.record.invariants.end());
    HUNDUN_CHECK(finite_value(exact->observed) == 0.0);
    HUNDUN_CHECK(finite_value(manifest->observed) == -1.0);
  }

  const auto failed =
      failure_report(mpi.rank(), Reason::file_integrity, Phase::rank_read);
  auto failed_source = hundun::flow::checkpoint_v2_diagnostic_source(failed);
  OneRecordSink failed_invariants;
  hundun::diagnostics::collect_diagnostics(
      failed_source, mpi,
      request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::invariants,
              hundun::diagnostics::DiagnosticScope::collective, "rank-read"),
      failed_invariants);
  const auto rank_crc = std::find_if(
      failed_invariants.record.invariants.begin(),
      failed_invariants.record.invariants.end(), [](const auto &invariant) {
        return invariant.id == "checkpoint.rank-crc";
      });
  HUNDUN_CHECK(rank_crc != failed_invariants.record.invariants.end());
  HUNDUN_CHECK(!rank_crc->passed);
  HUNDUN_CHECK(finite_value(rank_crc->observed) == 0.0);

  OneRecordSink failed_counters;
  hundun::diagnostics::collect_diagnostics(
      failed_source, mpi,
      request(mpi.rank(), hundun::diagnostics::DiagnosticLevel::counters,
              hundun::diagnostics::DiagnosticScope::collective, "rank-read"),
      failed_counters);
  const auto counter = [&](std::string_view id) {
    return std::find_if(failed_counters.record.counters.begin(),
                        failed_counters.record.counters.end(),
                        [&](const auto &value) { return value.id == id; });
  };
  const auto high = counter("checkpoint.local-crc64-high");
  const auto low = counter("checkpoint.local-crc64-low");
  HUNDUN_CHECK(high != failed_counters.record.counters.end());
  HUNDUN_CHECK(low != failed_counters.record.counters.end());
  const auto digest = independent_rank_crc_digest(mpi.size());
  HUNDUN_CHECK(high->value == (digest >> 32U));
  HUNDUN_CHECK(low->value == static_cast<std::uint32_t>(digest));

  auto unsupported = request(
      mpi.rank(), hundun::diagnostics::DiagnosticLevel::bounded_state_sample,
      hundun::diagnostics::DiagnosticScope::local);
  unsupported.sample_budget = 1U;
  HUNDUN_CHECK(rejects([&] {
    OneRecordSink rejected;
    hundun::diagnostics::collect_diagnostics(moved, unsupported, rejected);
  }));
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
    run(mpi);
  });
}

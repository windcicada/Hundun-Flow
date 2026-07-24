// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/checkpoint_v2_diagnostics.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"
#include "flow/src/checkpoint_v2_detail.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <algorithm>
#include <cstring>
#include <string>

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

class DiagnosticErrorSink final
    : public hundun::diagnostics::DiagnosticSink {
public:
  void submit(const hundun::diagnostics::DiagnosticRecord &) override {
    throw hundun::diagnostics::DiagnosticCollectionError(
        hundun::diagnostics::DiagnosticFailureClass::invalid_request,
        "test.sink.original", -1, "deliberate sink rejection");
  }
};

template <class Function>
bool rejects(Function &&function) {
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

hundun::diagnostics::DiagnosticRequest request(
    int rank, hundun::diagnostics::DiagnosticLevel level,
    hundun::diagnostics::DiagnosticScope scope) {
  return {level, scope, {rank, 4U, 0.125, "completed-marker"}, {}, 0U};
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
    accumulator.add(
        ids[index], static_cast<std::uint64_t>(value.rank()), 0U,
        hundun::diagnostics::describe_fp64(
            static_cast<double>(static_cast<std::uint32_t>(item))));
    accumulator.add(
        ids[index], static_cast<std::uint64_t>(value.rank()), 1U,
        hundun::diagnostics::describe_fp64(
            static_cast<double>(static_cast<std::uint32_t>(item >> 32U))));
  };
  const auto count = [&](std::size_t index, std::int64_t item) {
    accumulator.add(ids[index],
                    static_cast<std::uint64_t>(value.rank()), 0U,
                    hundun::diagnostics::describe_fp64(
                        static_cast<double>(item)));
  };
  limb(0U, value.collective_count());
  limb(1U, value.crc_check_count());
  count(2U, static_cast<std::uint8_t>(value.disposition()));
  count(3U, static_cast<std::uint8_t>(
                value.exact_size_and_eof_status()));
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
  count(25U,
        static_cast<std::uint8_t>(value.transaction_entry_status()));
  return accumulator.finish();
}

void run(const hundun::runtime::MpiContext &mpi) {
  const auto local_report = report(mpi.rank());
  const auto expected_fingerprint =
      expected_local_fingerprint(local_report);
  auto source =
      hundun::flow::checkpoint_v2_diagnostic_source(local_report);
  const auto descriptor =
      hundun::diagnostics::describe_diagnostics(source);
  HUNDUN_CHECK(descriptor.module_kind ==
               hundun::diagnostics::DiagnosticModuleKind::checkpoint);
  HUNDUN_CHECK(!hundun::diagnostics::has_capability(
      descriptor.capabilities,
      hundun::diagnostics::DiagnosticCapability::bounded_state_sample));
  const auto ids =
      hundun::diagnostics::diagnostic_fingerprint_field_ids(source);
  HUNDUN_CHECK(ids.size() == 26U);
  HUNDUN_CHECK(std::is_sorted(ids.begin(), ids.end()));

  for (const auto level : {
           hundun::diagnostics::DiagnosticLevel::summary,
           hundun::diagnostics::DiagnosticLevel::invariants,
           hundun::diagnostics::DiagnosticLevel::counters}) {
    OneRecordSink first;
    const auto local =
        request(mpi.rank(), level,
                hundun::diagnostics::DiagnosticScope::local);
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
    const auto global =
        request(mpi.rank(), level,
                hundun::diagnostics::DiagnosticScope::collective);
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
      const auto logical =
          find_metric("checkpoint.local-logical-bytes-low");
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
      const auto logical =
          find_counter("checkpoint.local-logical-bytes");
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
      hundun::diagnostics::collect_diagnostics(
          source, mpi, disagreement, rejected_sink);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      rejected =
          error.classification() ==
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
      HUNDUN_CHECK(
          error.classification() ==
          hundun::diagnostics::DiagnosticFailureClass::capability);
    }
  }
  {
    DiagnosticErrorSink rejected_sink;
    try {
      hundun::diagnostics::collect_diagnostics(
          moved,
          request(mpi.rank(),
                  hundun::diagnostics::DiagnosticLevel::summary,
                  hundun::diagnostics::DiagnosticScope::local),
          rejected_sink);
      HUNDUN_CHECK(false);
    } catch (const hundun::diagnostics::DiagnosticCollectionError &error) {
      HUNDUN_CHECK(
          error.classification() ==
          hundun::diagnostics::DiagnosticFailureClass::sink_failure);
      HUNDUN_CHECK(error.code() == "diagnostics.sink.submit");
    }
  }
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

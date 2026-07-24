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

void run(const hundun::runtime::MpiContext &mpi) {
  auto source =
      hundun::flow::checkpoint_v2_diagnostic_source(report(mpi.rank()));
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

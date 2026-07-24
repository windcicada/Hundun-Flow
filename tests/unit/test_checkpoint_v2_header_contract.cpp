// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/checkpoint_v2_diagnostics.hpp"
#include "hundun/flow/checkpoint_v2.hpp"
#include "checkpoint_v2_detail.hpp"
#include "tests/support/test_main.hpp"

#include <cstring>
#include <type_traits>

static_assert(std::is_copy_constructible_v<hundun::flow::CheckpointV2Report>);
static_assert(
    std::is_copy_constructible_v<hundun::flow::CheckpointV2ReadResult>);
static_assert(std::is_move_constructible_v<
              hundun::flow::CheckpointV2DiagnosticSource>);
static_assert(!std::is_copy_constructible_v<
              hundun::flow::CheckpointV2DiagnosticSource>);
static_assert(!std::is_move_assignable_v<
              hundun::flow::CheckpointV2DiagnosticSource>);
static_assert(noexcept(hundun::flow::detail::CheckpointV2Access::make(
    hundun::flow::detail::CheckpointV2ReportValues{})));
static_assert(noexcept(hundun::flow::detail::CheckpointV2Access::failed(
    hundun::flow::CheckpointV2Operation::write, 0,
    hundun::flow::CheckpointV2FailureReason::state,
    hundun::flow::CheckpointV2Phase::preflight)));

namespace {

void test_report_fixture_and_semantic_seal() {
  using namespace hundun::flow;
  detail::CheckpointV2ReportValues values;
  values.operation = CheckpointV2Operation::read;
  values.disposition = CheckpointV2Disposition::completed;
  values.reason = CheckpointV2FailureReason::none;
  values.phase = CheckpointV2Phase::restore_publish;
  values.rank = 3;
  values.lowest_failing_rank = -1;
  values.step = 17U;
  values.time_s = -0.0;
  values.local_logical_bytes = 101U;
  values.local_actual_bytes = 102U;
  values.global_logical_bytes = 103U;
  values.global_actual_bytes = 104U;
  values.local_crc64 = 105U;
  values.manifest_crc64 = 106U;
  values.file_count = 107U;
  values.crc_check_count = 108U;
  values.collective_count = 109U;
  values.rank_crc = CheckpointV2CheckStatus::passed;
  values.manifest_crc = CheckpointV2CheckStatus::failed;
  values.exact_size_eof = CheckpointV2CheckStatus::passed;
  values.fingerprint = CheckpointV2CheckStatus::failed;
  values.partition = CheckpointV2CheckStatus::passed;
  values.transaction_entry = CheckpointV2CheckStatus::failed;
  values.publication = CheckpointV2CheckStatus::passed;
  values.rollback = CheckpointV2CheckStatus::failed;
  const auto report = detail::CheckpointV2Access::make(values);
  HUNDUN_CHECK(report.operation() == values.operation);
  HUNDUN_CHECK(report.disposition() == values.disposition);
  HUNDUN_CHECK(report.reason() == values.reason);
  HUNDUN_CHECK(report.phase() == values.phase);
  HUNDUN_CHECK(report.rank() == values.rank);
  HUNDUN_CHECK(report.lowest_failing_rank() ==
               values.lowest_failing_rank);
  HUNDUN_CHECK(report.step() == values.step);
  std::uint64_t report_time{};
  std::uint64_t expected_time{};
  const double actual_time = report.time_s();
  std::memcpy(&report_time, &actual_time, sizeof(report_time));
  std::memcpy(&expected_time, &values.time_s, sizeof(expected_time));
  HUNDUN_CHECK(report_time == expected_time);
  HUNDUN_CHECK(report.local_logical_bytes() ==
               values.local_logical_bytes);
  HUNDUN_CHECK(report.local_actual_bytes() == values.local_actual_bytes);
  HUNDUN_CHECK(report.global_logical_bytes() ==
               values.global_logical_bytes);
  HUNDUN_CHECK(report.global_actual_bytes() == values.global_actual_bytes);
  HUNDUN_CHECK(report.local_crc64() == values.local_crc64);
  HUNDUN_CHECK(report.manifest_crc64() == values.manifest_crc64);
  HUNDUN_CHECK(report.file_count() == values.file_count);
  HUNDUN_CHECK(report.crc_check_count() == values.crc_check_count);
  HUNDUN_CHECK(report.collective_count() == values.collective_count);
  HUNDUN_CHECK(report.rank_crc_status() == values.rank_crc);
  HUNDUN_CHECK(report.manifest_crc_status() == values.manifest_crc);
  HUNDUN_CHECK(report.exact_size_and_eof_status() ==
               values.exact_size_eof);
  HUNDUN_CHECK(report.fingerprint_status() == values.fingerprint);
  HUNDUN_CHECK(report.partition_status() == values.partition);
  HUNDUN_CHECK(report.transaction_entry_status() ==
               values.transaction_entry);
  HUNDUN_CHECK(report.publication_status() == values.publication);
  HUNDUN_CHECK(report.rollback_status() == values.rollback);
  HUNDUN_CHECK(report.semantic_fingerprint() != 0U);

  const auto seal_changes = [&](auto mutate) {
    auto candidate = values;
    mutate(candidate);
    return detail::CheckpointV2Access::make(candidate)
               .semantic_fingerprint() != report.semantic_fingerprint();
  };
  HUNDUN_CHECK(seal_changes(
      [](auto &v) { v.operation = CheckpointV2Operation::write; }));
  HUNDUN_CHECK(seal_changes(
      [](auto &v) { v.disposition = CheckpointV2Disposition::failed; }));
  HUNDUN_CHECK(seal_changes(
      [](auto &v) { v.reason = CheckpointV2FailureReason::state; }));
  HUNDUN_CHECK(seal_changes(
      [](auto &v) { v.phase = CheckpointV2Phase::rank_read; }));
  HUNDUN_CHECK(seal_changes([](auto &v) { --v.rank; }));
  HUNDUN_CHECK(seal_changes([](auto &v) { v.lowest_failing_rank = 2; }));
  HUNDUN_CHECK(seal_changes([](auto &v) { ++v.step; }));
  HUNDUN_CHECK(seal_changes([](auto &v) { v.time_s = 0.0; }));
  HUNDUN_CHECK(seal_changes([](auto &v) { ++v.local_logical_bytes; }));
  HUNDUN_CHECK(seal_changes([](auto &v) { ++v.local_actual_bytes; }));
  HUNDUN_CHECK(seal_changes([](auto &v) { ++v.global_logical_bytes; }));
  HUNDUN_CHECK(seal_changes([](auto &v) { ++v.global_actual_bytes; }));
  HUNDUN_CHECK(seal_changes([](auto &v) { ++v.local_crc64; }));
  HUNDUN_CHECK(seal_changes([](auto &v) { ++v.manifest_crc64; }));
  HUNDUN_CHECK(seal_changes([](auto &v) { ++v.file_count; }));
  HUNDUN_CHECK(seal_changes([](auto &v) { ++v.crc_check_count; }));
  HUNDUN_CHECK(seal_changes([](auto &v) { ++v.collective_count; }));
  HUNDUN_CHECK(seal_changes(
      [](auto &v) { v.rank_crc = CheckpointV2CheckStatus::failed; }));
  HUNDUN_CHECK(seal_changes(
      [](auto &v) { v.manifest_crc = CheckpointV2CheckStatus::passed; }));
  HUNDUN_CHECK(seal_changes(
      [](auto &v) { v.exact_size_eof = CheckpointV2CheckStatus::failed; }));
  HUNDUN_CHECK(seal_changes(
      [](auto &v) { v.fingerprint = CheckpointV2CheckStatus::passed; }));
  HUNDUN_CHECK(seal_changes(
      [](auto &v) { v.partition = CheckpointV2CheckStatus::failed; }));
  HUNDUN_CHECK(seal_changes([](auto &v) {
    v.transaction_entry = CheckpointV2CheckStatus::passed;
  }));
  HUNDUN_CHECK(seal_changes(
      [](auto &v) { v.publication = CheckpointV2CheckStatus::failed; }));
  HUNDUN_CHECK(seal_changes(
      [](auto &v) { v.rollback = CheckpointV2CheckStatus::passed; }));

  const auto failure = detail::CheckpointV2Access::failed(
      CheckpointV2Operation::write, 2,
      CheckpointV2FailureReason::filesystem,
      CheckpointV2Phase::rank_temporary_file);
  HUNDUN_CHECK(failure.disposition() ==
               CheckpointV2Disposition::failed);
  HUNDUN_CHECK(failure.reason() ==
               CheckpointV2FailureReason::filesystem);
  HUNDUN_CHECK(failure.lowest_failing_rank() == 2);
}

} // namespace

int main() {
  return hundun::test::run([] {
    const auto descriptor =
        hundun::diagnostics::describe_checkpoint_v2_diagnostics();
    HUNDUN_CHECK(descriptor.schema_version == 1U);
    test_report_fixture_and_semantic_seal();
  });
}

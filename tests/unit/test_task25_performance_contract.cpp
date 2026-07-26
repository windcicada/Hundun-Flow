// SPDX-License-Identifier: Apache-2.0

#include "applications/hundun/stage2_performance_timing.hpp"
#include "hundun/diagnostics/performance_artifact.hpp"
#include "hundun/diagnostics/performance_correctness.hpp"
#include "hundun/execution/execution.hpp"
#include "execution_test_access.hpp"
#include "tests/support/task25_counter_checks.hpp"
#include "tests/support/test_main.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace {

using hundun::diagnostics::ArtifactMetadata;
using hundun::diagnostics::ComparisonMode;
using hundun::diagnostics::ComparisonStatus;
using hundun::diagnostics::CompatibilityMetadata;
using hundun::diagnostics::PerformanceCorrectnessRecord;
using hundun::diagnostics::PerformanceWorkRecord;
using hundun::execution::AllocationCounters;
using hundun::execution::Buffer;
using hundun::execution::CpuReferenceContext;
using hundun::execution::allocation_counters;
using hundun::execution::test::ExecutionTestAccess;
using hundun::runtime::HaloPerformanceCounters;

CompatibilityMetadata compatibility() {
  CompatibilityMetadata value;
  value.hardware_identity = "cpu";
  value.node_identity = "node";
  value.mpi_identity = "mpi";
  value.compiler_identity = "compiler";
  value.compiler_version = "1";
  value.compiler_flags = "-O2";
  value.link_flags = "";
  value.build_type = "Release";
  value.cpu_affinity = "0";
  value.rank_placement = "compact";
  value.problem_fingerprint = "case";
  value.numerical_tolerance_contract = "stage2";
  value.measurement_method = "mpi-wtime-v1";
  value.warmup_steps = 1;
  value.measured_steps = 2;
  value.repetitions = 2;
  value.execution_backend = "cpu_reference";
  value.ranks = 1;
  value.threads = 1;
  value.process_grid = {1, 1, 1};
  value.global_owned_cell_extents = {8, 8, 4};
  value.per_rank_owned_cell_extents = {8, 8, 4};
  return value;
}

void check_delta(AllocationCounters before, AllocationCounters after,
                 std::uint64_t allocations, std::uint64_t allocated,
                 std::uint64_t deallocations, std::uint64_t deallocated,
                 std::uint64_t live) {
  HUNDUN_CHECK(after.allocation_events - before.allocation_events ==
               allocations);
  HUNDUN_CHECK(after.allocated_bytes - before.allocated_bytes == allocated);
  HUNDUN_CHECK(after.deallocation_events - before.deallocation_events ==
               deallocations);
  HUNDUN_CHECK(after.deallocated_bytes - before.deallocated_bytes ==
               deallocated);
  HUNDUN_CHECK(after.live_bytes == live);
  HUNDUN_CHECK(after.peak_live_bytes >= after.live_bytes);
}

void test_counter_equality_helper() {
  const HaloPerformanceCounters halo{
      1U, 2U, 3U, 4U, 5U, 6U, 7U, 8U, 9U, -0.0};
  auto halo_copy = halo;
  HUNDUN_CHECK(hundun::test::task25_counters_equal(halo, halo_copy));
  ++halo_copy.receive_messages;
  HUNDUN_CHECK(!hundun::test::task25_counters_equal(halo, halo_copy));
  halo_copy = halo;
  halo_copy.completed_wait_seconds = 0.0;
  HUNDUN_CHECK(!hundun::test::task25_counters_equal(halo, halo_copy));

  const AllocationCounters allocation{1U, 2U, 3U, 4U, 5U, 6U};
  auto allocation_copy = allocation;
  HUNDUN_CHECK(
      hundun::test::task25_counters_equal(allocation, allocation_copy));
  ++allocation_copy.peak_live_bytes;
  HUNDUN_CHECK(
      !hundun::test::task25_counters_equal(allocation, allocation_copy));
}

void test_measured_region_boundary() {
  const std::vector<double> clock_values{1.0, 2.0, 10.0,
                                         14.0, 20.0, 26.0};
  std::size_t clock_index = 0U;
  std::vector<std::string> events;
  const auto clock = [&] {
    events.emplace_back("clock");
    HUNDUN_CHECK(clock_index < clock_values.size());
    return clock_values[clock_index++];
  };

  double elapsed = 0.0;
  events.emplace_back("bookkeeping-before");
  elapsed += hundun::application::detail::measure_performance_region(
      clock, [&] { events.emplace_back("advance"); });
  events.emplace_back("retry-validation");
  elapsed += hundun::application::detail::measure_performance_region(
      clock, [&] { events.emplace_back("checkpoint-write"); });
  events.emplace_back("checkpoint-counters");
  elapsed += hundun::application::detail::measure_performance_region(
      clock, [&] { events.emplace_back("diagnostic-write"); });
  events.emplace_back("work-retention");

  HUNDUN_CHECK(elapsed == 11.0);
  HUNDUN_CHECK(clock_index == clock_values.size());
  const std::vector<std::string> expected{
      "bookkeeping-before", "clock", "advance", "clock",
      "retry-validation",   "clock", "checkpoint-write", "clock",
      "checkpoint-counters", "clock", "diagnostic-write", "clock",
      "work-retention"};
  HUNDUN_CHECK(events == expected);
}

void test_allocation_counters() {
  static_assert(std::is_trivially_copyable_v<AllocationCounters>);
  CpuReferenceContext context;
  const auto begin = allocation_counters();
  {
    Buffer zero(context, 0U);
    const auto after_zero = allocation_counters();
    check_delta(begin, after_zero, 1U, 0U, 0U, 0U, begin.live_bytes);

    Buffer value(context, 32U);
    const auto after_value = allocation_counters();
    check_delta(begin, after_value, 2U, 32U, 0U, 0U,
                begin.live_bytes + 32U);

    Buffer moved(std::move(value));
    HUNDUN_CHECK(allocation_counters().allocation_events ==
                 after_value.allocation_events);
    moved.reallocate(64U);
    const auto after_reallocate = allocation_counters();
    check_delta(begin, after_reallocate, 3U, 96U, 1U, 32U,
                begin.live_bytes + 64U);
  }
  const auto end = allocation_counters();
  check_delta(begin, end, 3U, 96U, 3U, 96U, begin.live_bytes);

  const auto unchanged = allocation_counters();
  const auto identity_before_unsupported_construction =
      ExecutionTestAccess::next_allocation_identity();
  bool failed = false;
  try {
    Buffer rejected(context, std::numeric_limits<std::size_t>::max());
  } catch (const hundun::runtime::Error& error) {
    failed = std::string(error.what()).find("byte size") != std::string::npos;
  }
  HUNDUN_CHECK(failed);
  HUNDUN_CHECK(hundun::test::task25_counters_equal(
      unchanged, allocation_counters()));
  HUNDUN_CHECK(ExecutionTestAccess::next_allocation_identity() ==
               identity_before_unsupported_construction);

  const auto identity_before_injected_construction =
      ExecutionTestAccess::next_allocation_identity();
  ExecutionTestAccess::fail_next_allocation();
  failed = false;
  try {
    Buffer rejected(context, 8U);
  } catch (const hundun::runtime::Error&) {
    failed = true;
  }
  HUNDUN_CHECK(failed);
  HUNDUN_CHECK(allocation_counters().allocation_events ==
               unchanged.allocation_events);
  HUNDUN_CHECK(ExecutionTestAccess::next_allocation_identity() ==
               identity_before_injected_construction);

  ExecutionTestAccess::set_allocation_counters_for_test(
      {std::numeric_limits<std::uint64_t>::max(), 0U, 0U, 0U, 0U, 0U});
  const auto identity_before_construction_counter_overflow =
      ExecutionTestAccess::next_allocation_identity();
  failed = false;
  try {
    Buffer rejected(context, 0U);
  } catch (const hundun::runtime::Error&) {
    failed = true;
  }
  HUNDUN_CHECK(failed);
  HUNDUN_CHECK(ExecutionTestAccess::next_allocation_identity() ==
               identity_before_construction_counter_overflow);
  ExecutionTestAccess::set_allocation_counters_for_test(unchanged);

  {
    const auto before_move_assignment = allocation_counters();
    Buffer destination(context, 8U);
    Buffer source(context, 16U);
    const auto before_move = allocation_counters();
    destination = std::move(source);
    const auto after_move = allocation_counters();
    HUNDUN_CHECK(after_move.allocation_events ==
                 before_move.allocation_events);
    HUNDUN_CHECK(after_move.deallocation_events ==
                 before_move.deallocation_events + 1U);
    HUNDUN_CHECK(after_move.deallocated_bytes ==
                 before_move.deallocated_bytes + 8U);
    HUNDUN_CHECK(after_move.live_bytes ==
                 before_move_assignment.live_bytes + 16U);
    HUNDUN_CHECK(source.allocation_identity() == 0U);
  }

  {
    Buffer stable(context, 32U);
    auto view = stable.view(0U, 1U);
    view[0] = 7.0;
    const auto identity = stable.allocation_identity();
    const auto epoch = stable.epoch();
    const auto bytes = stable.byte_size();
    const auto before_failure = allocation_counters();
    const auto identity_before_unsupported_reallocation =
        ExecutionTestAccess::next_allocation_identity();
    failed = false;
    try {
      stable.reallocate(std::numeric_limits<std::size_t>::max());
    } catch (const hundun::runtime::Error& error) {
      failed =
          std::string(error.what()).find("byte size") != std::string::npos;
    }
    HUNDUN_CHECK(failed);
    HUNDUN_CHECK(hundun::test::task25_counters_equal(
        before_failure, allocation_counters()));
    HUNDUN_CHECK(stable.allocation_identity() == identity);
    HUNDUN_CHECK(stable.epoch() == epoch);
    HUNDUN_CHECK(stable.byte_size() == bytes);
    HUNDUN_CHECK(view[0] == 7.0);
    HUNDUN_CHECK(ExecutionTestAccess::next_allocation_identity() ==
                 identity_before_unsupported_reallocation);

    const auto identity_before_injected_reallocation =
        ExecutionTestAccess::next_allocation_identity();
    ExecutionTestAccess::fail_next_allocation();
    failed = false;
    try {
      stable.reallocate(64U);
    } catch (const hundun::runtime::Error&) {
      failed = true;
    }
    HUNDUN_CHECK(failed);
    HUNDUN_CHECK(hundun::test::task25_counters_equal(
        before_failure, allocation_counters()));
    HUNDUN_CHECK(stable.allocation_identity() == identity);
    HUNDUN_CHECK(stable.epoch() == epoch);
    HUNDUN_CHECK(stable.byte_size() == bytes);
    HUNDUN_CHECK(view[0] == 7.0);
    HUNDUN_CHECK(ExecutionTestAccess::next_allocation_identity() ==
                 identity_before_injected_reallocation);

    auto injected = before_failure;
    injected.allocation_events =
        std::numeric_limits<std::uint64_t>::max();
    ExecutionTestAccess::set_allocation_counters_for_test(injected);
    const auto identity_before_counter_overflow =
        ExecutionTestAccess::next_allocation_identity();
    failed = false;
    try {
      stable.reallocate(64U);
    } catch (const hundun::runtime::Error&) {
      failed = true;
    }
    HUNDUN_CHECK(failed);
    HUNDUN_CHECK(hundun::test::task25_counters_equal(
        injected, allocation_counters()));
    HUNDUN_CHECK(stable.allocation_identity() == identity);
    HUNDUN_CHECK(stable.epoch() == epoch);
    HUNDUN_CHECK(stable.byte_size() == bytes);
    HUNDUN_CHECK(view[0] == 7.0);
    HUNDUN_CHECK(ExecutionTestAccess::next_allocation_identity() ==
                 identity_before_counter_overflow);
    ExecutionTestAccess::set_allocation_counters_for_test(before_failure);
  }

  const auto before_future_overflow = allocation_counters();
  auto future_event = before_future_overflow;
  future_event.deallocation_events =
      std::numeric_limits<std::uint64_t>::max();
  ExecutionTestAccess::set_allocation_counters_for_test(future_event);
  const auto identity_before_future_event_overflow =
      ExecutionTestAccess::next_allocation_identity();
  failed = false;
  try {
    Buffer rejected(context, 0U);
  } catch (const hundun::runtime::Error&) {
    failed = true;
  }
  HUNDUN_CHECK(failed);
  HUNDUN_CHECK(hundun::test::task25_counters_equal(
      future_event, allocation_counters()));
  HUNDUN_CHECK(ExecutionTestAccess::next_allocation_identity() ==
               identity_before_future_event_overflow);
  ExecutionTestAccess::set_allocation_counters_for_test(
      before_future_overflow);

  auto future_bytes = before_future_overflow;
  future_bytes.deallocated_bytes =
      std::numeric_limits<std::uint64_t>::max() - 3U;
  ExecutionTestAccess::set_allocation_counters_for_test(future_bytes);
  const auto identity_before_future_byte_overflow =
      ExecutionTestAccess::next_allocation_identity();
  failed = false;
  try {
    Buffer rejected(context, 8U);
  } catch (const hundun::runtime::Error&) {
    failed = true;
  }
  HUNDUN_CHECK(failed);
  HUNDUN_CHECK(hundun::test::task25_counters_equal(
      future_bytes, allocation_counters()));
  HUNDUN_CHECK(ExecutionTestAccess::next_allocation_identity() ==
               identity_before_future_byte_overflow);
  ExecutionTestAccess::set_allocation_counters_for_test(
      before_future_overflow);
}

void test_source_comparison() {
  ArtifactMetadata clean_a{"abc", true, "", compatibility()};
  ArtifactMetadata clean_b = clean_a;
  auto result = hundun::diagnostics::compare_artifact_metadata(
      clean_a, clean_b, ComparisonMode::identical);
  HUNDUN_CHECK(result.status == ComparisonStatus::comparable);

  clean_b.commit = "def";
  result = hundun::diagnostics::compare_artifact_metadata(
      clean_a, clean_b, ComparisonMode::identical);
  HUNDUN_CHECK(result.status == ComparisonStatus::incomparable);
  HUNDUN_CHECK(result.reasons.size() == 1U);
  HUNDUN_CHECK(result.reasons[0] == "source.commit.mismatch");

  clean_a.commit = "unavailable";
  clean_b.commit = "unavailable";
  result = hundun::diagnostics::compare_artifact_metadata(
      clean_a, clean_b, ComparisonMode::identical);
  HUNDUN_CHECK(result.reasons[0] == "source.commit.unavailable");

  ArtifactMetadata dirty_a{"abc", false, "diff:0123", compatibility()};
  ArtifactMetadata dirty_b = dirty_a;
  HUNDUN_CHECK(hundun::diagnostics::compare_artifact_metadata(
                   dirty_a, dirty_b, ComparisonMode::identical)
                   .status == ComparisonStatus::comparable);
  dirty_b.dirty_summary = "diff:4567";
  result = hundun::diagnostics::compare_artifact_metadata(
      dirty_a, dirty_b, ComparisonMode::identical);
  HUNDUN_CHECK(result.reasons[0] == "source.dirty-summary.mismatch");
  dirty_b = dirty_a;
  dirty_b.dirty_summary = "unavailable";
  result = hundun::diagnostics::compare_artifact_metadata(
      dirty_a, dirty_b, ComparisonMode::identical);
  HUNDUN_CHECK(result.reasons[0] ==
               "source.dirty-summary.unavailable");

  clean_a = {"abc", true, "", compatibility()};
  dirty_a = {"abc", false, "diff:0123", compatibility()};
  result = hundun::diagnostics::compare_artifact_metadata(
      clean_a, dirty_a, ComparisonMode::identical);
  HUNDUN_CHECK(result.reasons[0] == "source.clean.mismatch");

  clean_a = {"abc", true, "", compatibility()};
  clean_a.compatibility.mpi_identity = "unavailable";
  clean_b = clean_a;
  result = hundun::diagnostics::compare_artifact_metadata(
      clean_a, clean_b, ComparisonMode::identical);
  HUNDUN_CHECK(result.reasons[0] ==
               "platform.required-identity.unavailable");
}

void test_correctness_record() {
  PerformanceCorrectnessRecord value;
  value.passed = true;
  value.allocation_bytes_per_owned_cell = 8.0;
  value.peak_allocation_bytes_per_owned_cell = 16.0;
  value.repetitions = 2U;
  value.states = {{0U, "hundun-performance-state-fp-v1:0123456789abcdef"},
                  {1U, "hundun-performance-state-fp-v1:0123456789abcdef"}};
  for (std::uint64_t repetition = 0; repetition < 2U; ++repetition) {
    for (char phase : {'W', 'M'}) {
      const std::uint64_t steps = phase == 'W' ? 1U : 2U;
      for (std::uint64_t step = 0; step < steps; ++step) {
        for (std::uint64_t slot = 0; slot < 5U; ++slot) {
          value.work.push_back(PerformanceWorkRecord{
              repetition, phase, step, slot, "converged", 1U, 2U, 3U, 4U,
              UINT64_C(0x3ff0000000000000),
              UINT64_C(0x3fe0000000000000),
              UINT64_C(0x3fd0000000000000)});
        }
      }
    }
  }
  const std::string encoded =
      hundun::diagnostics::serialize_performance_correctness(value);
  HUNDUN_CHECK(encoded.find_first_of(" \n\r\t") == std::string::npos);
  const auto parsed =
      hundun::diagnostics::parse_performance_correctness(encoded);
  HUNDUN_CHECK(hundun::diagnostics::serialize_performance_correctness(parsed) ==
               encoded);

  const std::array malformed_values{
      encoded.substr(0U, encoded.size() - 1U),
      encoded + ";state=0:hundun-performance-state-fp-v1:0123456789abcdef",
      std::string("hundun-performance-correctness-v1;passed=2")};
  for (const std::string& malformed : malformed_values) {
    bool rejected = false;
    try {
      static_cast<void>(
          hundun::diagnostics::parse_performance_correctness(malformed));
    } catch (const std::invalid_argument&) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  }

  auto reordered = value;
  std::swap(reordered.work[0], reordered.work[1]);
  const auto reject_record = [](const PerformanceCorrectnessRecord& record) {
    bool rejected = false;
    try {
      static_cast<void>(
          hundun::diagnostics::serialize_performance_correctness(record));
    } catch (const std::invalid_argument&) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  };
  reject_record(reordered);

  auto peak_below_assembly = value;
  peak_below_assembly.peak_allocation_bytes_per_owned_cell = 4.0;
  reject_record(peak_below_assembly);

  auto failed_termination = value;
  failed_termination.work.front().termination = "maximum_iterations";
  reject_record(failed_termination);
  failed_termination.passed = false;
  static_cast<void>(
      hundun::diagnostics::serialize_performance_correctness(
          failed_termination));

  auto nan_residual = value;
  nan_residual.work.front().initial_residual_bits =
      UINT64_C(0x7ff8000000000000);
  reject_record(nan_residual);

  auto negative_residual = value;
  negative_residual.work.front().independent_final_residual_bits =
      UINT64_C(0xbff0000000000000);
  reject_record(negative_residual);

  auto changed_state = value;
  changed_state.states[1].second =
      "hundun-performance-state-fp-v1:1123456789abcdef";
  reject_record(changed_state);

  constexpr std::size_t second_repetition = 15U;
  auto changed_iteration = value;
  ++changed_iteration.work[second_repetition].iterations;
  reject_record(changed_iteration);

  auto changed_residual = value;
  changed_residual.work[second_repetition].recursive_residual_bits =
      UINT64_C(0x3fc0000000000000);
  reject_record(changed_residual);

  auto nonpassing_difference = changed_iteration;
  nonpassing_difference.passed = false;
  static_cast<void>(
      hundun::diagnostics::serialize_performance_correctness(
          nonpassing_difference));

  auto mixed_group = value;
  mixed_group.work[1].relative_step = 1U;
  reject_record(mixed_group);

  auto missing_slot = value;
  missing_slot.work.erase(missing_slot.work.begin() + 1);
  reject_record(missing_slot);

  auto missing_final_group = value;
  missing_final_group.work.erase(missing_final_group.work.end() - 5,
                                 missing_final_group.work.end());
  reject_record(missing_final_group);

  auto missing_phase_group = value;
  missing_phase_group.work.erase(missing_phase_group.work.begin() + 10,
                                 missing_phase_group.work.begin() + 15);
  missing_phase_group.work.erase(missing_phase_group.work.begin() + 20,
                                 missing_phase_group.work.begin() + 25);
  static_cast<void>(
      hundun::diagnostics::serialize_performance_correctness(
          missing_phase_group));
  bool rejected = false;
  try {
    hundun::diagnostics::validate_performance_correctness_coverage(
        missing_phase_group, 1U, 2U, 2U);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  auto missing_repetition = value;
  missing_repetition.states.pop_back();
  missing_repetition.repetitions = 1U;
  missing_repetition.work.erase(
      std::remove_if(missing_repetition.work.begin(),
                     missing_repetition.work.end(),
                     [](const auto& item) { return item.repetition == 1U; }),
      missing_repetition.work.end());
  static_cast<void>(
      hundun::diagnostics::serialize_performance_correctness(
          missing_repetition));
  rejected = false;
  try {
    hundun::diagnostics::validate_performance_correctness_coverage(
        missing_repetition, 1U, 2U, 2U);
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  std::string invalid_termination = encoded;
  const auto termination_offset = invalid_termination.find("converged");
  HUNDUN_CHECK(termination_offset != std::string::npos);
  invalid_termination.replace(termination_offset,
                              std::string("converged").size(), "unknown");
  rejected = false;
  try {
    static_cast<void>(hundun::diagnostics::parse_performance_correctness(
        invalid_termination));
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  std::string failed_termination_bits = encoded;
  const auto accepted_termination_offset =
      failed_termination_bits.find("converged");
  HUNDUN_CHECK(accepted_termination_offset != std::string::npos);
  failed_termination_bits.replace(
      accepted_termination_offset, std::string("converged").size(),
      "maximum_iterations");
  rejected = false;
  try {
    static_cast<void>(hundun::diagnostics::parse_performance_correctness(
        failed_termination_bits));
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  std::string nan_residual_bits = encoded;
  const auto nan_offset =
      nan_residual_bits.find("3ff0000000000000");
  HUNDUN_CHECK(nan_offset != std::string::npos);
  nan_residual_bits.replace(nan_offset, 16U, "7ff8000000000000");
  rejected = false;
  try {
    static_cast<void>(hundun::diagnostics::parse_performance_correctness(
        nan_residual_bits));
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  std::string negative_residual_bits = encoded;
  const auto negative_offset =
      negative_residual_bits.find("3fd0000000000000");
  HUNDUN_CHECK(negative_offset != std::string::npos);
  negative_residual_bits.replace(negative_offset, 16U,
                                 "bff0000000000000");
  rejected = false;
  try {
    static_cast<void>(hundun::diagnostics::parse_performance_correctness(
        negative_residual_bits));
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  std::string invalid_bits = encoded;
  const auto bits_offset = invalid_bits.find("3ff0000000000000");
  HUNDUN_CHECK(bits_offset != std::string::npos);
  invalid_bits[bits_offset] = 'g';
  rejected = false;
  try {
    static_cast<void>(
        hundun::diagnostics::parse_performance_correctness(invalid_bits));
  } catch (const std::invalid_argument&) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

}  // namespace

int main() {
  return hundun::test::run([] {
    test_counter_equality_helper();
    test_measured_region_boundary();
    test_allocation_counters();
    test_source_comparison();
    test_correctness_record();
  });
}

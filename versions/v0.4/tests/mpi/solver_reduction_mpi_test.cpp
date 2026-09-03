// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_linear.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <string_view>
#include <type_traits>
#include <utility>

namespace allocation_observer {

std::atomic<bool> enabled{false};
std::atomic<std::size_t> count{0U};

void observe() noexcept {
  if (enabled.load(std::memory_order_relaxed)) {
    count.fetch_add(1U, std::memory_order_relaxed);
  }
}

void* allocate(std::size_t size) {
  observe();
  void* const result = std::malloc(size == 0U ? 1U : size);
  if (result == nullptr) {
    throw std::bad_alloc{};
  }
  return result;
}

void* allocate_aligned(std::size_t size, std::size_t alignment) {
  observe();
  void* result = nullptr;
  const std::size_t requested = size == 0U ? alignment : size;
  if (posix_memalign(&result, alignment, requested) != 0 || result == nullptr) {
    throw std::bad_alloc{};
  }
  return result;
}

class Guard {
 public:
  Guard() noexcept {
    count.store(0U, std::memory_order_relaxed);
    enabled.store(true, std::memory_order_release);
  }
  ~Guard() { enabled.store(false, std::memory_order_release); }

  Guard(const Guard&) = delete;
  Guard& operator=(const Guard&) = delete;
};

}  // namespace allocation_observer

void* operator new(std::size_t size) {
  return allocation_observer::allocate(size);
}

void* operator new[](std::size_t size) {
  return allocation_observer::allocate(size);
}

void* operator new(std::size_t size, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}

void* operator new[](std::size_t size, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}

void* operator new(std::size_t size, const std::nothrow_t&) noexcept {
  try {
    return allocation_observer::allocate(size);
  } catch (...) {
    return nullptr;
  }
}

void* operator new[](std::size_t size, const std::nothrow_t&) noexcept {
  try {
    return allocation_observer::allocate(size);
  } catch (...) {
    return nullptr;
  }
}

void operator delete(void* pointer, const std::nothrow_t&) noexcept {
  std::free(pointer);
}

void operator delete[](void* pointer, const std::nothrow_t&) noexcept {
  std::free(pointer);
}

void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }
void operator delete(void* pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete[](void* pointer, std::size_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::size_t, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::size_t,
                       std::align_val_t) noexcept {
  std::free(pointer);
}

namespace {

using namespace hundun::v04;

constexpr std::uint32_t kCapacity = 8U;

static_assert(!std::is_copy_constructible_v<ReductionEngine>);
static_assert(!std::is_copy_assignable_v<ReductionEngine>);
static_assert(std::is_nothrow_move_constructible_v<ReductionEngine>);
static_assert(std::is_nothrow_move_assignable_v<ReductionEngine>);

bool expect(bool condition, int rank, std::string_view description) {
  if (!condition) {
    std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  }
  return condition;
}

bool all_true(bool local) {
  const int input = local ? 1 : 0;
  int output = 0;
  return MPI_Allreduce(&input, &output, 1, MPI_INT, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         output != 0;
}

std::uint64_t packed_status(Status status) noexcept {
  return (static_cast<std::uint64_t>(status.code) << 32U) |
         static_cast<std::uint64_t>(status.detail);
}

bool identical_u64(std::uint64_t value) {
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  return MPI_Allreduce(&value, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         MPI_Allreduce(&value, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         minimum == maximum;
}

std::uint64_t bits(double value) noexcept {
  std::uint64_t result = 0U;
  static_assert(sizeof(result) == sizeof(value));
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

bool same(LinearReductionCounters left,
          LinearReductionCounters right) noexcept {
  return left.calls == right.calls && left.scalars == right.scalars &&
         left.logical_bytes == right.logical_bytes &&
         left.tree_messages == right.tree_messages;
}

LinearReductionCounters delta(LinearReductionCounters before,
                              LinearReductionCounters after) noexcept {
  return LinearReductionCounters{
      after.calls - before.calls, after.scalars - before.scalars,
      after.logical_bytes - before.logical_bytes,
      after.tree_messages - before.tree_messages};
}

std::uint64_t expected_tree_messages(ReductionMode mode, int size,
                                     std::uint64_t calls = 1U) noexcept {
  return mode == ReductionMode::reproducible_tree
             ? calls * 2U * static_cast<std::uint64_t>(size - 1)
             : 0U;
}

bool unchanged(const std::array<double, 4U>& values,
               const std::array<double, 4U>& expected) noexcept {
  for (std::size_t index = 0U; index < values.size(); ++index) {
    if (bits(values[index]) != bits(expected[index])) {
      return false;
    }
  }
  return true;
}

std::array<double, 4U> exact_sum_oracle(int size) noexcept {
  switch (size) {
    case 1:
      return {1.0, 10.0, -2.0, 0.25};
    case 2:
      return {3.0, 19.0, 1.0, 0.75};
    case 4:
      return {10.0, 34.0, 2.0, 2.5};
    default:
      return {};
  }
}

std::array<double, 4U> exact_max_oracle(int size) noexcept {
  switch (size) {
    case 1:
      return {1.0, 10.0, -2.0, 0.25};
    case 2:
      return {2.0, 10.0, 3.0, 0.5};
    case 4:
      return {4.0, 10.0, 3.0, 1.0};
    default:
      return {};
  }
}

bool equal_exact(const std::array<double, 4U>& actual,
                 const std::array<double, 4U>& expected) noexcept {
  for (std::size_t index = 0U; index < actual.size(); ++index) {
    if (actual[index] != expected[index]) {
      return false;
    }
  }
  return true;
}

bool test_compile_contract_and_atomic_replacement(ReductionMode mode,
                                                  int rank) {
  ReductionEngine engine;
  bool passed = expect(
      static_cast<bool>(ReductionEngine::compile(MPI_COMM_WORLD, mode,
                                                  kCapacity, engine)),
      rank, "reduction engine compiles at registered capacity");
  passed = all_true(passed);
  if (!passed) {
    return false;
  }

  const std::uintptr_t send = engine.send_storage_address();
  const std::uintptr_t receive = engine.receive_storage_address();
  passed &= expect(engine.mode() == mode && engine.capacity() == kCapacity &&
                       send != 0U && receive != 0U &&
                       (mode != ReductionMode::reproducible_tree ||
                        send != receive) &&
                       engine.lowest_failing_rank() == -1 &&
                       same(engine.counters(), {}),
                   rank,
                   "compiled engine exposes two persistent buffers and zero counters");

  const Status zero_capacity =
      ReductionEngine::compile(MPI_COMM_WORLD, mode, 0U, engine);
  passed &= expect(zero_capacity.code == StatusCode::invalid_plan &&
                       engine.mode() == mode && engine.capacity() == kCapacity &&
                       engine.send_storage_address() == send &&
                       engine.receive_storage_address() == receive &&
                       same(engine.counters(), {}),
                   rank,
                   "zero-capacity replacement is rejected atomically");

  const Status null_communicator =
      ReductionEngine::compile(MPI_COMM_NULL, mode, kCapacity, engine);
  passed &= expect(null_communicator.code == StatusCode::invalid_plan &&
                       engine.capacity() == kCapacity &&
                       engine.send_storage_address() == send &&
                       engine.receive_storage_address() == receive &&
                       same(engine.counters(), {}),
                   rank,
                   "null-communicator replacement is rejected atomically");
  return all_true(passed);
}

bool test_packed_hand_oracles(ReductionMode mode, int rank, int size) {
  ReductionEngine engine;
  bool passed = expect(
      static_cast<bool>(ReductionEngine::compile(MPI_COMM_WORLD, mode,
                                                  kCapacity, engine)),
      rank, "oracle reduction engine compiles");
  passed = all_true(passed);
  if (!passed) {
    return false;
  }

  const std::array<double, 4U> local{
      static_cast<double>(rank + 1), static_cast<double>(10 - rank),
      rank % 2 == 0 ? -2.0 : 3.0, 0.25 * static_cast<double>(rank + 1)};
  std::array<double, 4U> global{-91.0, -92.0, -93.0, -94.0};
  const Status sum = engine.checked_sum(
      Span<const double>{local.data(), local.size()},
      Span<double>{global.data(), global.size()});
  passed &= expect(static_cast<bool>(sum) &&
                       equal_exact(global, exact_sum_oracle(size)) &&
                       engine.lowest_failing_rank() == -1 &&
                       same(engine.counters(),
                            LinearReductionCounters{
                                1U, 4U, 4U * sizeof(double),
                                expected_tree_messages(mode, size)}),
                   rank,
                   "one packed sum matches the hand-derived exact oracle and counters");

  global = {-81.0, -82.0, -83.0, -84.0};
  const Status maximum = engine.checked_max(
      Span<const double>{local.data(), local.size()},
      Span<double>{global.data(), global.size()});
  passed &= expect(static_cast<bool>(maximum) &&
                       equal_exact(global, exact_max_oracle(size)) &&
                       engine.lowest_failing_rank() == -1 &&
                       same(engine.counters(),
                            LinearReductionCounters{
                                2U, 8U, 8U * sizeof(double),
                                expected_tree_messages(mode, size, 2U)}),
                   rank,
                   "one packed maximum matches the hand-derived exact oracle and counters");
  return all_true(passed);
}

bool test_deterministic_maximum_locations(ReductionMode mode, int rank,
                                          int size) {
  ReductionEngine engine;
  bool passed = expect(
      static_cast<bool>(ReductionEngine::compile(MPI_COMM_WORLD, mode,
                                                  kCapacity, engine)),
      rank, "maximum-location reduction engine compiles");
  passed = all_true(passed);
  if (!passed) return false;

  const std::uint64_t tied_location =
      rank == 0 ? 7U : (rank == 1 || rank == 2 ? 3U : 9U);
  std::array<ReductionMaximumLocation, 2U> local{{
      {true, 10.0, tied_location, rank,
       {100.0 + rank, 101.0 + rank, 102.0 + rank,
        103.0 + rank, 104.0 + rank}},
      {true, static_cast<double>(rank + 1),
       static_cast<std::uint64_t>(100 - rank), rank,
       {200.0 + rank, 201.0 + rank, 202.0 + rank,
        203.0 + rank, 204.0 + rank}},
  }};
  std::array<ReductionMaximumLocation, 2U> global{};
  const Status reduced = engine.checked_max_locations(
      {local.data(), local.size()}, {global.data(), global.size()});
  const int tied_rank = size == 1 ? 0 : 1;
  const std::uint64_t expected_tied_location = size == 1 ? 7U : 3U;
  const LinearReductionCounters counters = engine.counters();
  passed &= expect(
      static_cast<bool>(reduced) && global[0U].valid &&
          global[0U].value == 10.0 &&
          global[0U].global_location == expected_tied_location &&
          global[0U].rank == tied_rank &&
          global[0U].payload[0U] == 100.0 + tied_rank &&
          global[1U].valid &&
          global[1U].value == static_cast<double>(size) &&
          global[1U].rank == size - 1 &&
          global[1U].global_location ==
              static_cast<std::uint64_t>(101 - size) &&
          global[1U].payload[4U] == 203.0 + size &&
          counters.calls == 1U && counters.scalars == 2U &&
          counters.logical_bytes ==
              2U * sizeof(ReductionMaximumLocation) &&
          counters.tree_messages == expected_tree_messages(mode, size) &&
          counters.blocking_operations == 2U,
      rank,
      "maximum-location tie uses global cell then rank and carries payload");

  const auto before_invalid = engine.counters();
  const auto unchanged_global = global;
  if (rank == size - 1) local[0U].rank = rank + 1;
  const Status invalid = engine.checked_max_locations(
      {local.data(), local.size()}, {global.data(), global.size()});
  passed &= expect(
      invalid.code == StatusCode::numerical_failure &&
          global[0U].global_location == unchanged_global[0U].global_location &&
          global[1U].global_location == unchanged_global[1U].global_location &&
          same(engine.counters(), before_invalid),
      rank,
      "rank-local maximum-location corruption fails collectively and atomically");
  return all_true(passed);
}

bool test_allreduce_roundoff_tolerance(int rank, int size) {
  ReductionEngine engine;
  bool passed = expect(
      static_cast<bool>(ReductionEngine::compile(
          MPI_COMM_WORLD, ReductionMode::mpi_allreduce, kCapacity, engine)),
      rank, "performance allreduce engine compiles");
  passed = all_true(passed);
  if (!passed) {
    return false;
  }

  const std::array<double, 2U> local{
      0.1 * static_cast<double>(rank + 1),
      (1.0 / 3.0) * static_cast<double>(rank + 1)};
  std::array<double, 2U> global{};
  const std::array<double, 2U> expected =
      size == 1 ? std::array<double, 2U>{0.1, 1.0 / 3.0}
                : (size == 2 ? std::array<double, 2U>{0.3, 1.0}
                             : std::array<double, 2U>{1.0, 10.0 / 3.0});
  const Status reduced = engine.checked_sum(
      Span<const double>{local.data(), local.size()},
      Span<double>{global.data(), global.size()});
  const double tolerance = 16.0 * std::numeric_limits<double>::epsilon();
  passed &= expect(static_cast<bool>(reduced) &&
                       std::abs(global[0] - expected[0]) <= tolerance &&
                       std::abs(global[1] - expected[1]) <=
                           tolerance * std::max(1.0, std::abs(expected[1])) &&
                       same(engine.counters(),
                            LinearReductionCounters{1U, 2U,
                                                    2U * sizeof(double), 0U}),
                   rank,
                   "MPI_Allreduce mode agrees with the FP64 roundoff oracle");
  return all_true(passed);
}

bool test_local_validation_is_collective_and_atomic(ReductionMode mode,
                                                    int rank, int size) {
  ReductionEngine engine;
  bool passed = expect(
      static_cast<bool>(ReductionEngine::compile(MPI_COMM_WORLD, mode,
                                                  kCapacity, engine)),
      rank, "validation reduction engine compiles");
  passed = all_true(passed);
  if (!passed) {
    return false;
  }

  const std::array<double, 4U> baseline{-101.0, -102.0, -103.0, -104.0};
  std::array<double, 4U> local{1.0, 2.0, 3.0, 4.0};
  std::array<double, 4U> global = baseline;
  const LinearReductionCounters zero = engine.counters();
  const std::uintptr_t send = engine.send_storage_address();
  const std::uintptr_t receive = engine.receive_storage_address();

  const Status empty = engine.checked_sum(Span<const double>{},
                                           Span<double>{});
  passed &= expect(empty.code == StatusCode::invalid_plan &&
                       same(engine.counters(), zero) &&
                       engine.send_storage_address() == send &&
                       engine.receive_storage_address() == receive,
                   rank, "empty reduction is rejected without publishing state");

  global = baseline;
  const Status mismatched = engine.checked_sum(
      Span<const double>{local.data(), 3U},
      Span<double>{global.data(), global.size()});
  passed &= expect(mismatched.code == StatusCode::invalid_plan &&
                       unchanged(global, baseline) &&
                       same(engine.counters(), zero),
                   rank,
                   "mismatched reduction spans preserve output and counters");

  std::array<double, kCapacity + 1U> oversized_local{};
  std::array<double, kCapacity + 1U> oversized_global{};
  oversized_global.fill(-205.0);
  const Status oversized = engine.checked_max(
      Span<const double>{oversized_local.data(), oversized_local.size()},
      Span<double>{oversized_global.data(), oversized_global.size()});
  bool oversized_unchanged = true;
  for (double value : oversized_global) {
    oversized_unchanged = oversized_unchanged && value == -205.0;
  }
  passed &= expect(oversized.code == StatusCode::invalid_plan &&
                       oversized_unchanged && same(engine.counters(), zero),
                   rank, "capacity overflow is rejected atomically");

  std::array<double, 5U> alias_values{7.0, 8.0, 9.0, 10.0, 11.0};
  const std::array<double, 5U> alias_before = alias_values;
  const Status exact_alias = engine.checked_sum(
      Span<const double>{alias_values.data(), 4U},
      Span<double>{alias_values.data(), 4U});
  passed &= expect(exact_alias.code == StatusCode::invalid_plan &&
                       alias_values == alias_before &&
                       same(engine.counters(), zero),
                   rank, "exact input/output alias is rejected atomically");

  const Status partial_alias = engine.checked_max(
      Span<const double>{alias_values.data(), 4U},
      Span<double>{alias_values.data() + 1U, 4U});
  passed &= expect(partial_alias.code == StatusCode::invalid_plan &&
                       alias_values == alias_before &&
                       same(engine.counters(), zero),
                   rank, "partial input/output overlap is rejected atomically");

  local = {1.0, 2.0, 3.0, 4.0};
  global = baseline;
  const std::size_t rank_local_input_size =
      rank == size - 1 ? local.size() - 1U : local.size();
  const Status rank_local_mismatch = engine.checked_sum(
      Span<const double>{local.data(), rank_local_input_size},
      Span<double>{global.data(), global.size()});
  const bool rank_local_mismatch_identical =
      identical_u64(packed_status(rank_local_mismatch));
  passed &= expect(rank_local_mismatch.code == StatusCode::invalid_plan &&
                       rank_local_mismatch.detail != 0U &&
                       rank_local_mismatch_identical &&
                       engine.lowest_failing_rank() == size - 1 &&
                       unchanged(global, baseline) &&
                       same(engine.counters(), zero),
                   rank,
                   "one-rank size mismatch reaches collective atomic failure");

  local = {1.0, 2.0, 3.0, 4.0};
  global = baseline;
  alias_values = alias_before;
  const Span<const double> rank_local_alias_input =
      rank == size - 1
          ? Span<const double>{alias_values.data(), 4U}
          : Span<const double>{local.data(), local.size()};
  const Span<double> rank_local_alias_output =
      rank == size - 1
          ? Span<double>{alias_values.data() + 1U, 4U}
          : Span<double>{global.data(), global.size()};
  const Status rank_local_alias = engine.checked_max(
      rank_local_alias_input, rank_local_alias_output);
  const bool rank_local_alias_identical =
      identical_u64(packed_status(rank_local_alias));
  const bool rank_local_output_unchanged =
      rank == size - 1 ? alias_values == alias_before
                       : unchanged(global, baseline);
  passed &= expect(rank_local_alias.code == StatusCode::invalid_plan &&
                       rank_local_alias.detail != 0U &&
                       rank_local_alias_identical &&
                       engine.lowest_failing_rank() == size - 1 &&
                       rank_local_output_unchanged &&
                       same(engine.counters(), zero),
                   rank,
                   "one-rank partial alias reaches collective atomic failure");

  local = {static_cast<double>(rank + 1), static_cast<double>(10 - rank),
           rank % 2 == 0 ? -2.0 : 3.0,
           0.25 * static_cast<double>(rank + 1)};
  global = baseline;
  if (rank == size - 1) {
    local[2] = std::numeric_limits<double>::quiet_NaN();
  }
  const Status nonfinite = engine.checked_sum(
      Span<const double>{local.data(), local.size()},
      Span<double>{global.data(), global.size()});
  const bool nonfinite_identical = identical_u64(packed_status(nonfinite));
  passed &= expect(nonfinite.code == StatusCode::numerical_failure &&
                       nonfinite.detail != 0U &&
                       nonfinite_identical &&
                       engine.lowest_failing_rank() == size - 1 &&
                       unchanged(global, baseline) &&
                       same(engine.counters(), zero) &&
                       engine.send_storage_address() == send &&
                       engine.receive_storage_address() == receive,
                   rank,
                   "one-rank NaN is classified identically and rejected atomically");

  local = {static_cast<double>(rank + 1), static_cast<double>(10 - rank),
           rank % 2 == 0 ? -2.0 : 3.0,
           0.25 * static_cast<double>(rank + 1)};
  global = baseline;
  const Status recovery = engine.checked_sum(
      Span<const double>{local.data(), local.size()},
      Span<double>{global.data(), global.size()});
  passed &= expect(static_cast<bool>(recovery) &&
                       engine.lowest_failing_rank() == -1 &&
                       equal_exact(global, exact_sum_oracle(size)) &&
                       same(engine.counters(),
                            LinearReductionCounters{
                                1U, 4U, 4U * sizeof(double),
                                expected_tree_messages(mode, size)}),
                   rank, "engine remains usable after validation failure");

  if (size > 1) {
    const std::array<double, 1U> finite_overflow_local{
        std::numeric_limits<double>::max()};
    const std::array<double, 1U> finite_overflow_baseline{-411.0};
    std::array<double, 1U> finite_overflow_global =
        finite_overflow_baseline;
    const LinearReductionCounters before_overflow = engine.counters();
    const Status overflow = engine.checked_sum(
        Span<const double>{finite_overflow_local.data(),
                           finite_overflow_local.size()},
        Span<double>{finite_overflow_global.data(),
                     finite_overflow_global.size()});
    const bool overflow_identical = identical_u64(packed_status(overflow));
    passed &= expect(overflow.code == StatusCode::numerical_failure &&
                         overflow.detail != 0U &&
                         overflow_identical &&
                         finite_overflow_global ==
                             finite_overflow_baseline &&
                         same(engine.counters(), before_overflow) &&
                         engine.send_storage_address() == send &&
                         engine.receive_storage_address() == receive,
                     rank,
                     "finite local DBL_MAX values reject a non-finite global sum atomically");
  }
  return all_true(passed);
}

bool test_lowest_rank_failure_propagation(ReductionMode mode, int rank,
                                          int size) {
  ReductionEngine engine;
  bool passed = expect(
      static_cast<bool>(ReductionEngine::compile(MPI_COMM_WORLD, mode,
                                                  kCapacity, engine)),
      rank, "failure-protocol reduction engine compiles");
  passed = all_true(passed);
  if (!passed) {
    return false;
  }

  const std::array<double, 2U> local{static_cast<double>(rank + 1), 2.0};
  const std::array<double, 2U> baseline{-301.0, -302.0};
  std::array<double, 2U> global = baseline;
  Status local_status{};
  if (rank == 0) {
    local_status = {StatusCode::rejected_step, 8101U};
  } else if (rank == size - 1) {
    local_status = {StatusCode::numerical_failure, 8199U};
  }

  const LinearReductionCounters before = engine.counters();
  const Status failure = engine.checked_max(
      Span<const double>{local.data(), local.size()},
      Span<double>{global.data(), global.size()}, local_status);
  const bool failure_identical = identical_u64(packed_status(failure));
  passed &= expect(failure.code == StatusCode::rejected_step &&
                       failure.detail == 8101U &&
                       failure_identical &&
                       engine.lowest_failing_rank() == 0 &&
                       global == baseline &&
                       same(engine.counters(), before),
                   rank,
                   "lowest failing rank supplies identical status/detail without publication");

  const Status recovery = engine.checked_max(
      Span<const double>{local.data(), local.size()},
      Span<double>{global.data(), global.size()});
  passed &= expect(static_cast<bool>(recovery) &&
                       engine.lowest_failing_rank() == -1 &&
                       global[0] == static_cast<double>(size) &&
                       global[1] == 2.0 &&
                       same(delta(before, engine.counters()),
                            LinearReductionCounters{
                                1U, 2U, 2U * sizeof(double),
                                expected_tree_messages(mode, size)}),
                   rank,
                   "collective local-status failure does not poison recovery");
  return all_true(passed);
}

bool test_fixed_tree_bitwise_repetition(int rank, int size) {
  ReductionEngine engine;
  bool passed = expect(
      static_cast<bool>(ReductionEngine::compile(
          MPI_COMM_WORLD, ReductionMode::reproducible_tree, kCapacity,
          engine)),
      rank, "fixed-tree engine compiles");
  passed = all_true(passed);
  if (!passed) {
    return false;
  }

  const std::array<double, 3U> local{
      rank == 0 ? 1.0e16
                : (rank == 1 ? 1.0 : (rank == 2 ? -1.0e16 : 3.0)),
      rank % 2 == 0 ? 1.0 / 10.0 : -1.0 / 7.0,
      static_cast<double>((rank + 1) * (rank + 1))};
  std::array<double, 3U> reference{};
  passed &= expect(static_cast<bool>(engine.checked_sum(
                       Span<const double>{local.data(), local.size()},
                       Span<double>{reference.data(), reference.size()})),
                   rank, "first fixed-tree reduction succeeds");
  const std::array<std::uint64_t, 3U> reference_bits{
      bits(reference[0]), bits(reference[1]), bits(reference[2])};
  for (std::uint64_t value : reference_bits) {
    passed &= identical_u64(value);
  }

  constexpr std::uint64_t repetitions = 31U;
  bool repeated_exactly = true;
  for (std::uint64_t repetition = 0U; repetition < repetitions;
       ++repetition) {
    std::array<double, 3U> result{};
    const Status status = engine.checked_sum(
        Span<const double>{local.data(), local.size()},
        Span<double>{result.data(), result.size()});
    repeated_exactly = repeated_exactly && static_cast<bool>(status) &&
                       bits(result[0]) == reference_bits[0] &&
                       bits(result[1]) == reference_bits[1] &&
                       bits(result[2]) == reference_bits[2];
  }
  passed &= expect(repeated_exactly &&
                       same(engine.counters(),
                            LinearReductionCounters{
                                repetitions + 1U, 3U * (repetitions + 1U),
                                3U * (repetitions + 1U) * sizeof(double),
                                expected_tree_messages(
                                    ReductionMode::reproducible_tree, size,
                                    repetitions + 1U)}),
                   rank,
                   "same communicator and rank order reproduce every FP64 bit");

  ReductionEngine second;
  passed &= expect(
      static_cast<bool>(ReductionEngine::compile(
          MPI_COMM_WORLD, ReductionMode::reproducible_tree, kCapacity,
          second)),
      rank, "second fixed-tree engine compiles on the same rank order");
  std::array<double, 3U> second_result{};
  const Status second_status = second.checked_sum(
      Span<const double>{local.data(), local.size()},
      Span<double>{second_result.data(), second_result.size()});
  passed &= expect(static_cast<bool>(second_status) &&
                       bits(second_result[0]) == reference_bits[0] &&
                       bits(second_result[1]) == reference_bits[1] &&
                       bits(second_result[2]) == reference_bits[2],
                   rank,
                   "independent fixed-tree engines preserve the same evaluation order");
  return all_true(passed);
}

bool test_hot_reuse_and_exact_counters(ReductionMode mode, int rank,
                                       int size) {
  ReductionEngine engine;
  bool passed = expect(
      static_cast<bool>(ReductionEngine::compile(MPI_COMM_WORLD, mode,
                                                  kCapacity, engine)),
      rank, "hot-path reduction engine compiles");
  passed = all_true(passed);
  if (!passed) {
    return false;
  }

  const std::array<double, 4U> local{
      static_cast<double>(rank + 1), static_cast<double>(10 - rank),
      rank % 2 == 0 ? -2.0 : 3.0, 0.25 * static_cast<double>(rank + 1)};
  std::array<double, 4U> global{};
  const std::uintptr_t send = engine.send_storage_address();
  const std::uintptr_t receive = engine.receive_storage_address();
  const LinearReductionCounters before = engine.counters();
  constexpr std::uint64_t repetitions = 100U;
  bool hot_ok = true;
  std::size_t hot_allocations = std::numeric_limits<std::size_t>::max();
  {
    allocation_observer::Guard guard;
    for (std::uint64_t repetition = 0U; repetition < repetitions;
         ++repetition) {
      const Status status = engine.checked_sum(
          Span<const double>{local.data(), local.size()},
          Span<double>{global.data(), global.size()});
      hot_ok = hot_ok && static_cast<bool>(status);
    }
    hot_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  const LinearReductionCounters actual_delta =
      delta(before, engine.counters());
  passed &= expect(hot_ok && hot_allocations == 0U &&
                       engine.send_storage_address() == send &&
                       engine.receive_storage_address() == receive &&
                       equal_exact(global, exact_sum_oracle(size)) &&
                       same(actual_delta,
                            LinearReductionCounters{
                                repetitions, 4U * repetitions,
                                4U * repetitions * sizeof(double),
                                expected_tree_messages(mode, size,
                                                       repetitions)}),
                   rank,
                   "100 hot reductions allocate nothing, reuse buffers, and count exactly");
  return all_true(passed);
}

bool test_move_ownership(int rank, int size) {
  ReductionEngine source;
  bool passed = expect(
      static_cast<bool>(ReductionEngine::compile(
          MPI_COMM_WORLD, ReductionMode::reproducible_tree, kCapacity,
          source)),
      rank, "move source engine compiles");
  passed = all_true(passed);
  if (!passed) {
    return false;
  }

  const std::array<double, 1U> local{static_cast<double>(rank + 1)};
  std::array<double, 1U> global{};
  passed &= expect(static_cast<bool>(source.checked_sum(
                       Span<const double>{local.data(), local.size()},
                       Span<double>{global.data(), global.size()})),
                   rank, "move source records live state");
  const std::uintptr_t send = source.send_storage_address();
  const std::uintptr_t receive = source.receive_storage_address();
  const LinearReductionCounters counters = source.counters();

  ReductionEngine moved(std::move(source));
  passed &= expect(source.capacity() == 0U &&
                       source.send_storage_address() == 0U &&
                       source.receive_storage_address() == 0U &&
                       moved.capacity() == kCapacity &&
                       moved.mode() == ReductionMode::reproducible_tree &&
                       moved.send_storage_address() == send &&
                       moved.receive_storage_address() == receive &&
                       same(moved.counters(), counters),
                   rank,
                   "move construction transfers communicator, buffers, and counters exactly");

  ReductionEngine destination;
  passed &= expect(
      static_cast<bool>(ReductionEngine::compile(
          MPI_COMM_WORLD, ReductionMode::mpi_allreduce, 2U, destination)),
      rank, "move-assignment destination initially owns another engine");
  destination = std::move(moved);
  passed &= expect(moved.capacity() == 0U &&
                       moved.send_storage_address() == 0U &&
                       moved.receive_storage_address() == 0U &&
                       destination.capacity() == kCapacity &&
                       destination.mode() ==
                           ReductionMode::reproducible_tree &&
                       destination.send_storage_address() == send &&
                       destination.receive_storage_address() == receive &&
                       same(destination.counters(), counters),
                   rank,
                   "move assignment releases old ownership and transfers live state");

  global[0] = -1.0;
  const Status recovery = destination.checked_sum(
      Span<const double>{local.data(), local.size()},
      Span<double>{global.data(), global.size()});
  passed &= expect(static_cast<bool>(recovery) &&
                       global[0] ==
                           0.5 * static_cast<double>(size) *
                               static_cast<double>(size + 1) &&
                       destination.counters().calls == counters.calls + 1U,
                   rank, "moved-to engine remains operational");
  return all_true(passed);
}

bool test_communicator_compatibility_is_local_and_read_only(int rank,
                                                            int size) {
  ReductionEngine engine;
  bool passed = expect(
      static_cast<bool>(ReductionEngine::compile(
          MPI_COMM_WORLD, ReductionMode::mpi_allreduce, kCapacity, engine)),
      rank, "communicator-contract reduction engine compiles");
  passed = all_true(passed);
  if (!passed) {
    return false;
  }

  MPI_Comm congruent = MPI_COMM_NULL;
  MPI_Comm reversed = MPI_COMM_NULL;
  passed &= expect(MPI_Comm_dup(MPI_COMM_WORLD, &congruent) == MPI_SUCCESS &&
                       MPI_Comm_split(MPI_COMM_WORLD, 0, size - 1 - rank,
                                      &reversed) == MPI_SUCCESS,
                   rank, "comparison communicators compile");
  if (!all_true(passed)) {
    if (congruent != MPI_COMM_NULL) (void)MPI_Comm_free(&congruent);
    if (reversed != MPI_COMM_NULL) (void)MPI_Comm_free(&reversed);
    return false;
  }

  const LinearReductionCounters before = engine.counters();
  const std::uintptr_t send = engine.send_storage_address();
  const std::uintptr_t receive = engine.receive_storage_address();
  const Status identical = engine.validate_communicator(MPI_COMM_WORLD);
  const Status compatible = engine.validate_communicator(congruent);
  const Status similar = engine.validate_communicator(reversed);
  const Status null = engine.validate_communicator(MPI_COMM_NULL);
  passed &= expect(static_cast<bool>(identical) &&
                       static_cast<bool>(compatible) &&
                       (size == 1 ? static_cast<bool>(similar)
                                  : similar.code == StatusCode::invalid_plan) &&
                       null.code == StatusCode::invalid_plan &&
                       same(engine.counters(), before) &&
                       engine.send_storage_address() == send &&
                       engine.receive_storage_address() == receive,
                   rank,
                   "CONGRUENT is accepted, multi-rank SIMILAR and null reject without mutation");

  const std::array<double, 1U> local{static_cast<double>(rank + 1)};
  std::array<double, 1U> global{-1.0};
  const Status recovered = engine.checked_sum(
      Span<const double>{local.data(), local.size()},
      Span<double>{global.data(), global.size()});
  passed &= expect(static_cast<bool>(recovered) &&
                       global[0] == 0.5 * static_cast<double>(size) *
                                        static_cast<double>(size + 1),
                   rank,
                   "rejected communicator query leaves reductions operational");

  (void)MPI_Comm_free(&reversed);
  (void)MPI_Comm_free(&congruent);
  return all_true(passed);
}

bool run_mode(ReductionMode mode, int rank, int size) {
  bool passed = test_compile_contract_and_atomic_replacement(mode, rank);
  passed &= test_packed_hand_oracles(mode, rank, size);
  passed &= test_deterministic_maximum_locations(mode, rank, size);
  passed &= test_local_validation_is_collective_and_atomic(mode, rank, size);
  passed &= test_lowest_rank_failure_propagation(mode, rank, size);
  passed &= test_hot_reuse_and_exact_counters(mode, rank, size);
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  bool passed = expect(size == 1 || size == 2 || size == 4, rank,
                       "reduction RED runs only at registered 1/2/4 ranks");
  passed &= run_mode(ReductionMode::reproducible_tree, rank, size);
  passed &= run_mode(ReductionMode::mpi_allreduce, rank, size);
  passed &= test_allreduce_roundoff_tolerance(rank, size);
  passed &= test_fixed_tree_bitwise_repetition(rank, size);
  passed &= test_move_ownership(rank, size);
  passed &= test_communicator_compatibility_is_local_and_read_only(rank, size);
  passed = all_true(passed);
  if (rank == 0 && !passed) {
    std::cerr << "v0.4 solver reduction MPI test failed for " << size
              << " ranks\n";
  }
  MPI_Finalize();
  return passed ? 0 : 1;
}

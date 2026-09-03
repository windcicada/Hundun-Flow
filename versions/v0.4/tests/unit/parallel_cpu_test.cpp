// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_parallel.hpp"

#include <mpi.h>
#include <sched.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <string_view>
#include <thread>
#include <type_traits>
#include <vector>

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
  void* result = std::malloc(size == 0U ? 1U : size);
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

namespace hundun::v04::detail {

void fail_cpu_thread_launch_after_for_test(
    std::size_t successful_launches) noexcept;
void reset_cpu_thread_launch_failure_for_test() noexcept;

}  // namespace hundun::v04::detail

namespace {

using namespace hundun::v04;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool same(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool finite_choice(const CpuExecutionPlan& plan) noexcept {
  switch (plan.kernel_variant()) {
    case CpuKernelVariant::scalar:
      return same(plan.tile_shape(), {32, 4, 4});
    case CpuKernelVariant::avx2:
      return same(plan.tile_shape(), {64, 4, 4});
    case CpuKernelVariant::avx512:
      return same(plan.tile_shape(), {128, 4, 2});
  }
  return false;
}

std::vector<std::int32_t> cores(const CpuExecutionPlan& plan) {
  const Span<const std::int32_t> selected = plan.core_ids();
  return {selected.data, selected.data + selected.size};
}

bool test_cold_default_and_explicit_plans() {
  static_assert(!std::is_copy_constructible_v<CpuExecutionPlan>);
  static_assert(!std::is_copy_assignable_v<CpuExecutionPlan>);
  static_assert(std::is_nothrow_move_constructible_v<CpuExecutionPlan>);
  static_assert(std::is_nothrow_move_assignable_v<CpuExecutionPlan>);

  CpuExecutionPlan plan;
  bool passed = expect(static_cast<bool>(CpuExecutionPlan::compile(
                           MPI_COMM_SELF, CpuExecutionRequest{}, plan)),
                       "default CPU plan compiles from cold discovery");
  const CpuExecutionDiagnostics diagnostics = plan.diagnostics();
  passed &= expect(plan.worker_count() > 0U &&
                       plan.worker_count() <= diagnostics.allowed_core_count,
                   "default workers stay inside the process affinity mask");
  passed &= expect(diagnostics.allowed_core_count > 0U &&
                       diagnostics.numa_node_count > 0U &&
                       diagnostics.local_rank_count == 1U &&
                       diagnostics.recommended_local_rank_count > 0U,
                   "cold diagnostics expose affinity NUMA and rank guidance");
  passed &= expect(!plan.pure_mpi() && finite_choice(plan),
                   "default plan freezes one finite ISA and tile choice");
  const std::vector<std::int32_t> discovered = cores(plan);
  passed &= expect(std::is_sorted(discovered.begin(), discovered.end()),
                   "automatic core placement is deterministic");

  CpuExecutionPlan repeated;
  passed &= expect(static_cast<bool>(CpuExecutionPlan::compile(
                       MPI_COMM_SELF, CpuExecutionRequest{}, repeated)) &&
                       cores(repeated) == discovered &&
                       repeated.kernel_variant() == plan.kernel_variant() &&
                       same(repeated.tile_shape(), plan.tile_shape()),
                   "repeated cold compilation freezes the same placement");

  const std::array<std::int32_t, 2U> duplicate{
      discovered.front(), discovered.front()};
  CpuExecutionRequest duplicate_request;
  duplicate_request.threads_per_rank = 2U;
  duplicate_request.core_ids = {duplicate.data(), duplicate.size()};
  CpuExecutionPlan explicit_plan;
  passed &= expect(static_cast<bool>(CpuExecutionPlan::compile(
                       MPI_COMM_SELF, duplicate_request, explicit_plan)) &&
                       explicit_plan.worker_count() == duplicate.size() &&
                       explicit_plan.diagnostics().duplicate_core_warning &&
                       cores(explicit_plan) ==
                           std::vector<std::int32_t>(duplicate.begin(),
                                                     duplicate.end()),
                   "duplicate explicit cores are retained with a warning");

  CpuExecutionRequest mismatch = duplicate_request;
  mismatch.threads_per_rank = 1U;
  const std::vector<std::int32_t> retained = cores(explicit_plan);
  const Status mismatch_status = CpuExecutionPlan::compile(
      MPI_COMM_SELF, mismatch, explicit_plan);
  passed &= expect(mismatch_status.code == StatusCode::invalid_plan &&
                       cores(explicit_plan) == retained,
                   "explicit thread/core mismatch rejects atomically");

  std::int32_t unavailable = std::numeric_limits<std::int32_t>::max();
  while (std::find(discovered.begin(), discovered.end(), unavailable) !=
         discovered.end()) {
    --unavailable;
  }
  CpuExecutionRequest invalid_core;
  invalid_core.threads_per_rank = 1U;
  invalid_core.core_ids = {&unavailable, 1U};
  const Status invalid_status = CpuExecutionPlan::compile(
      MPI_COMM_SELF, invalid_core, explicit_plan);
  passed &= expect(invalid_status.code == StatusCode::invalid_plan &&
                       cores(explicit_plan) == retained,
                   "core outside process affinity rejects atomically");

  const Status null_status = CpuExecutionPlan::compile(
      MPI_COMM_NULL, CpuExecutionRequest{}, explicit_plan);
  passed &= expect(null_status.code == StatusCode::invalid_plan &&
                       cores(explicit_plan) == retained,
                   "null communicator rejects without replacing the plan");
  return passed;
}

struct CountContext {
  std::atomic<std::uint64_t>* calls{};
  std::size_t count{};
};

void count_task(std::size_t worker, void* opaque) noexcept {
  auto& context = *static_cast<CountContext*>(opaque);
  if (worker < context.count) {
    context.calls[worker].fetch_add(1U, std::memory_order_relaxed);
  }
}

struct DispatchContext {
  std::atomic<std::uint64_t> calls{};
  std::atomic<std::uint64_t> worker_sum{};
  std::atomic<bool> entered{};
  std::atomic<bool>* release{};
};

void dispatch_task(std::size_t worker, void* opaque) noexcept {
  auto& context = *static_cast<DispatchContext*>(opaque);
  if (worker == 0U && context.release != nullptr) {
    context.entered.store(true, std::memory_order_release);
    while (!context.release->load(std::memory_order_acquire)) {
      std::this_thread::yield();
    }
  }
  context.worker_sum.fetch_add(static_cast<std::uint64_t>(worker) + 1U,
                               std::memory_order_relaxed);
  context.calls.fetch_add(1U, std::memory_order_release);
}

struct InlineContext {
  std::thread::id caller{};
  std::thread::id actual{};
  std::size_t worker{std::numeric_limits<std::size_t>::max()};
  std::uint64_t calls{};
};

void inline_task(std::size_t worker, void* opaque) noexcept {
  auto& context = *static_cast<InlineContext*>(opaque);
  context.actual = std::this_thread::get_id();
  context.worker = worker;
  ++context.calls;
}

bool test_inline_pure_mpi() {
  CpuExecutionRequest request;
  request.threads_per_rank = 999U;
  request.pure_mpi = true;
  request.bind_threads = true;
  CpuExecutionPlan plan;
  bool passed = expect(static_cast<bool>(CpuExecutionPlan::compile(
                           MPI_COMM_SELF, request, plan)) &&
                           plan.worker_count() == 1U && plan.pure_mpi() &&
                           !plan.binds_threads(),
                       "pure MPI overrides thread count and binding");

  CpuThreadTeam team;
  passed &= expect(static_cast<bool>(CpuThreadTeam::create(plan, team)) &&
                       team.worker_count() == 1U,
                   "pure MPI team creates without a worker thread");
  InlineContext context;
  context.caller = std::this_thread::get_id();
  passed &= expect(static_cast<bool>(team.run(inline_task, &context)) &&
                       context.calls == 1U && context.worker == 0U &&
                       context.actual == context.caller,
                   "pure MPI dispatch executes worker zero inline");
  return passed;
}

bool test_fixed_team_and_zero_allocation() {
  CpuExecutionPlan discovered;
  if (!CpuExecutionPlan::compile(MPI_COMM_SELF, CpuExecutionRequest{},
                                 discovered)) {
    return expect(false, "discovery plan exists for fixed-team test");
  }
  const std::vector<std::int32_t> available = cores(discovered);
  const std::size_t worker_count = std::min<std::size_t>(2U, available.size());
  CpuExecutionRequest request;
  request.threads_per_rank = static_cast<std::uint32_t>(worker_count);
  request.core_ids = {available.data(), worker_count};
  CpuExecutionPlan plan;
  bool passed = expect(static_cast<bool>(CpuExecutionPlan::compile(
                           MPI_COMM_SELF, request, plan)),
                       "fixed-team plan compiles");
  CpuThreadTeam team;
  passed &= expect(static_cast<bool>(CpuThreadTeam::create(plan, team)) &&
                       team.worker_count() == worker_count,
                   "fixed C++17 worker team starts once");

  std::array<std::atomic<std::uint64_t>, 2U> calls{};
  CountContext context{calls.data(), worker_count};
  passed &= expect(static_cast<bool>(team.run(count_task, &context)),
                   "fixed team completes an initial dispatch");

  std::size_t allocations = 0U;
  {
    allocation_observer::Guard guard;
    for (std::size_t iteration = 0; iteration < 1000U; ++iteration) {
      if (!team.run(count_task, &context)) {
        passed = false;
        break;
      }
    }
    allocations = allocation_observer::count.load(std::memory_order_acquire);
  }
  passed &= expect(allocations == 0U,
                   "repeated fixed-team dispatch allocates nothing");
  for (std::size_t worker = 0; worker < worker_count; ++worker) {
    passed &= expect(calls[worker].load(std::memory_order_relaxed) == 1001U,
                     "every fixed worker runs exactly once per dispatch");
  }


  std::atomic<bool> release_first{false};
  DispatchContext first;
  first.release = &release_first;
  DispatchContext second;
  Status first_status;
  Status second_status;
  std::thread first_caller([&] {
    first_status = team.run(dispatch_task, &first);
  });
  while (!first.entered.load(std::memory_order_acquire)) {
    std::this_thread::yield();
  }
  std::thread second_caller([&] {
    second_status = team.run(dispatch_task, &second);
  });
  for (std::size_t spin = 0U; spin < 10000U; ++spin) {
    std::this_thread::yield();
  }
  release_first.store(true, std::memory_order_release);
  first_caller.join();
  second_caller.join();
  const std::uint64_t expected_worker_sum =
      static_cast<std::uint64_t>(worker_count * (worker_count + 1U) / 2U);
  passed &= expect(first_status && second_status &&
                       first.calls.load(std::memory_order_acquire) ==
                           worker_count &&
                       second.calls.load(std::memory_order_acquire) ==
                           worker_count &&
                       first.worker_sum.load(std::memory_order_relaxed) ==
                           expected_worker_sum &&
                       second.worker_sum.load(std::memory_order_relaxed) ==
                           expected_worker_sum,
                   "concurrent callers serialize complete dispatch contexts");
  passed &= expect(team.run(nullptr, nullptr).code == StatusCode::invalid_plan,
                   "null CPU task is rejected");

  CpuThreadTeam moved = std::move(team);
  passed &= expect(team.worker_count() == 0U &&
                       moved.worker_count() == worker_count &&
                       static_cast<bool>(moved.run(count_task, &context)),
                   "fixed team ownership moves without rebuilding workers");
  return passed;
}

struct BindingContext {
  std::int32_t observed{-1};
};

void observe_cpu(std::size_t, void* opaque) noexcept {
  static_cast<BindingContext*>(opaque)->observed = sched_getcpu();
}

bool test_startup_binding() {
  CpuExecutionPlan discovered;
  if (!CpuExecutionPlan::compile(MPI_COMM_SELF, CpuExecutionRequest{},
                                 discovered)) {
    return expect(false, "discovery plan exists for binding test");
  }
  const std::int32_t selected = discovered.core_ids().data[0];
  CpuExecutionRequest request;
  request.threads_per_rank = 1U;
  request.bind_threads = true;
  request.core_ids = {&selected, 1U};
  CpuExecutionPlan plan;
  bool passed = expect(static_cast<bool>(CpuExecutionPlan::compile(
                           MPI_COMM_SELF, request, plan)) &&
                           plan.binds_threads(),
                       "explicit binding plan validates an allowed core");
  CpuThreadTeam team;
  passed &= expect(static_cast<bool>(CpuThreadTeam::create(plan, team)),
                   "bound worker completes startup binding");
  BindingContext context;
  passed &= expect(static_cast<bool>(team.run(observe_cpu, &context)) &&
                       context.observed == selected,
                   "persistent worker remains on its startup core");
  return passed;
}

bool test_partial_launch_failure_is_raii_safe() {
  CpuExecutionRequest inline_request;
  inline_request.pure_mpi = true;
  CpuExecutionPlan inline_plan;
  CpuThreadTeam retained;
  bool passed = expect(static_cast<bool>(CpuExecutionPlan::compile(
                           MPI_COMM_SELF, inline_request, inline_plan)) &&
                           static_cast<bool>(
                               CpuThreadTeam::create(inline_plan, retained)),
                       "retained inline team exists before injected failure");

  CpuExecutionPlan discovered;
  passed &= expect(static_cast<bool>(CpuExecutionPlan::compile(
                       MPI_COMM_SELF, CpuExecutionRequest{}, discovered)),
                   "discovery plan exists for launch injection");
  const std::vector<std::int32_t> available = cores(discovered);
  const std::size_t requested = std::min<std::size_t>(2U, available.size());
  CpuExecutionRequest thread_request;
  thread_request.threads_per_rank = static_cast<std::uint32_t>(requested);
  thread_request.core_ids = {available.data(), requested};
  CpuExecutionPlan thread_plan;
  passed &= expect(static_cast<bool>(CpuExecutionPlan::compile(
                       MPI_COMM_SELF, thread_request, thread_plan)),
                   "thread plan exists for launch injection");

  detail::fail_cpu_thread_launch_after_for_test(requested > 1U ? 1U : 0U);
  const Status injected = CpuThreadTeam::create(thread_plan, retained);
  detail::reset_cpu_thread_launch_failure_for_test();
  InlineContext context;
  context.caller = std::this_thread::get_id();
  passed &= expect(injected.code == StatusCode::allocation_failure &&
                       retained.worker_count() == 1U &&
                       static_cast<bool>(
                           retained.run(inline_task, &context)) &&
                       context.actual == context.caller,
                   "partial launch failure joins workers and preserves output");

  CpuExecutionPlan invalid_plan;
  const Status invalid = CpuThreadTeam::create(invalid_plan, retained);
  passed &= expect(invalid.code == StatusCode::invalid_plan &&
                       retained.worker_count() == 1U,
                   "invalid plan cannot replace a live team");
  return passed;
}

bool test_world_placement_and_failure_consensus() {
  CpuExecutionPlan retained;
  const Status baseline = CpuExecutionPlan::compile(
      MPI_COMM_WORLD, CpuExecutionRequest{}, retained);
  bool passed = expect(static_cast<bool>(baseline),
                       "world default CPU plan compiles collectively");
  if (!baseline) {
    return false;
  }
  const std::vector<std::int32_t> retained_cores = cores(retained);

  MPI_Comm local = MPI_COMM_NULL;
  passed &= expect(MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                                      MPI_INFO_NULL, &local) == MPI_SUCCESS,
                   "test obtains shared-memory communicator");
  if (local == MPI_COMM_NULL) {
    return false;
  }
  int local_size = 0;
  MPI_Comm_size(local, &local_size);
  const std::uint64_t local_count = retained_cores.size();
  std::vector<std::uint64_t> wide_counts(static_cast<std::size_t>(local_size));
  MPI_Allgather(&local_count, 1, MPI_UINT64_T, wide_counts.data(), 1,
                MPI_UINT64_T, local);
  std::vector<int> counts(static_cast<std::size_t>(local_size));
  std::vector<int> displacements(static_cast<std::size_t>(local_size));
  int total = 0;
  for (int rank = 0; rank < local_size; ++rank) {
    counts[static_cast<std::size_t>(rank)] =
        static_cast<int>(wide_counts[static_cast<std::size_t>(rank)]);
    displacements[static_cast<std::size_t>(rank)] = total;
    total += counts[static_cast<std::size_t>(rank)];
  }
  std::vector<std::int32_t> gathered(static_cast<std::size_t>(total));
  MPI_Allgatherv(retained_cores.data(), static_cast<int>(retained_cores.size()),
                 MPI_INT32_T, gathered.data(), counts.data(),
                 displacements.data(), MPI_INT32_T, local);
  std::sort(gathered.begin(), gathered.end());
  passed &= expect(std::adjacent_find(gathered.begin(), gathered.end()) ==
                           gathered.end() &&
                       retained.diagnostics().local_rank_count ==
                           static_cast<std::uint32_t>(local_size),
                   "default shared-rank core ownership never overlaps");
  MPI_Comm_free(&local);

  int world_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  const std::int32_t requested_core =
      world_rank == 0 ? std::numeric_limits<std::int32_t>::max()
                      : retained_cores.front();
  CpuExecutionRequest poisoned;
  poisoned.threads_per_rank = 1U;
  poisoned.core_ids = {&requested_core, 1U};
  const Status rejected = CpuExecutionPlan::compile(
      MPI_COMM_WORLD, poisoned, retained);
  const std::uint64_t packed =
      static_cast<std::uint64_t>(rejected.detail) << 16U |
      static_cast<std::uint64_t>(rejected.code);
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  MPI_Allreduce(&packed, &minimum, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&packed, &maximum, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  passed &= expect(rejected.code == StatusCode::invalid_plan &&
                       minimum == maximum && cores(retained) == retained_cores,
                   "one-rank invalid request returns one status and preserves all plans");

  if (local_size > 1) {
    cpu_set_t process_mask;
    CPU_ZERO(&process_mask);
    const bool affinity_ok = sched_getaffinity(0, sizeof(process_mask),
                                               &process_mask) == 0;
    std::array<int, CPU_SETSIZE> present{};
    for (int core = 0; core < CPU_SETSIZE; ++core) {
      present[static_cast<std::size_t>(core)] =
          affinity_ok && CPU_ISSET(core, &process_mask) != 0 ? 1 : 0;
    }
    std::array<int, CPU_SETSIZE> shared{};
    MPI_Allreduce(present.data(), shared.data(), CPU_SETSIZE, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    std::int32_t shared_core = -1;
    for (int core = 0; core < CPU_SETSIZE; ++core) {
      if (shared[static_cast<std::size_t>(core)] != 0) {
        shared_core = core;
        break;
      }
    }
    if (shared_core >= 0) {
      CpuExecutionRequest collision;
      collision.threads_per_rank = 1U;
      collision.core_ids = {&shared_core, 1U};
      CpuExecutionPlan collision_plan;
      const Status collision_status = CpuExecutionPlan::compile(
          MPI_COMM_WORLD, collision, collision_plan);
      passed &= expect(
          collision_status &&
              collision_plan.diagnostics().duplicate_core_warning,
          "cross-rank explicit core collision publishes a warning");
    }
  }
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  int provided = MPI_THREAD_SINGLE;
  if (MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided) !=
      MPI_SUCCESS) {
    return 2;
  }
  bool passed = expect(provided >= MPI_THREAD_FUNNELED,
                       "MPI supplies funneled thread support");
  passed &= test_cold_default_and_explicit_plans();
  passed &= test_inline_pure_mpi();
  passed &= test_fixed_team_and_zero_allocation();
  passed &= test_startup_binding();
  passed &= test_partial_launch_failure_is_raii_safe();
  passed &= test_world_placement_and_failure_consensus();
  MPI_Finalize();
  if (!passed) {
    return 1;
  }
  std::cout << "v0.4 CPU placement tests passed\n";
  return 0;
}

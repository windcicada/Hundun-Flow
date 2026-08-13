// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_parallel.hpp"

#include <mpi.h>
#include <pthread.h>
#include <sched.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <memory>
#include <mutex>
#include <new>
#include <string>
#include <thread>
#include <vector>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kCpuCommunicator = 601U;
constexpr std::uint32_t kCpuAffinity = 602U;
constexpr std::uint32_t kCpuRequest = 603U;
constexpr std::uint32_t kCpuCore = 604U;
constexpr std::uint32_t kCpuThread = 605U;
constexpr std::uint32_t kCpuBinding = 606U;

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
std::atomic<std::size_t> g_fail_thread_launch_after{
    std::numeric_limits<std::size_t>::max()};
#endif

Status invalid(std::uint32_t detail) noexcept {
  return {StatusCode::invalid_plan, detail};
}

bool valid_communicator(MPI_Comm communicator) noexcept {
  return communicator != MPI_COMM_NULL;
}

Status consensus_status(MPI_Comm communicator, Status local) noexcept {
  int rank = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kCpuCommunicator};
  }
  const int local_failure = local ? std::numeric_limits<int>::max() : rank;
  int first_failure = std::numeric_limits<int>::max();
  if (MPI_Allreduce(&local_failure, &first_failure, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kCpuCommunicator};
  }
  if (first_failure == std::numeric_limits<int>::max()) {
    return {};
  }
  std::uint32_t wire[2]{0U, 0U};
  if (rank == first_failure) {
    wire[0] = static_cast<std::uint32_t>(local.code);
    wire[1] = local.detail;
  }
  if (MPI_Bcast(wire, 2, MPI_UINT32_T, first_failure, communicator) !=
      MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kCpuCommunicator};
  }
  return {static_cast<StatusCode>(wire[0]), wire[1]};
}

Status discover_allowed_cores(std::vector<std::int32_t>& cores,
                              bool& constrained) noexcept {
  try {
    long configured = ::sysconf(_SC_NPROCESSORS_CONF);
    if (configured <= 0) {
      configured = CPU_SETSIZE;
    }
    std::size_t capacity = static_cast<std::size_t>(configured);
    capacity = std::max<std::size_t>(capacity, CPU_SETSIZE);
    constexpr std::size_t kMaximumCpuId = UINT32_C(1048576);

    for (;;) {
      if (capacity > kMaximumCpuId) {
        return invalid(kCpuAffinity);
      }
      const std::size_t bytes = CPU_ALLOC_SIZE(capacity);
      cpu_set_t* const mask = CPU_ALLOC(capacity);
      if (mask == nullptr) {
        return {StatusCode::allocation_failure, kCpuAffinity};
      }
      CPU_ZERO_S(bytes, mask);
      errno = 0;
      const int result = sched_getaffinity(0, bytes, mask);
      const int saved_errno = errno;
      if (result == 0) {
        cores.clear();
        cores.reserve(capacity);
        for (std::size_t cpu = 0; cpu < capacity; ++cpu) {
          if (CPU_ISSET_S(cpu, bytes, mask) != 0) {
            if (cpu > static_cast<std::size_t>(
                          std::numeric_limits<std::int32_t>::max())) {
              CPU_FREE(mask);
              return invalid(kCpuAffinity);
            }
            cores.push_back(static_cast<std::int32_t>(cpu));
          }
        }
        CPU_FREE(mask);
        if (cores.empty()) {
          return invalid(kCpuAffinity);
        }
        const long online = ::sysconf(_SC_NPROCESSORS_ONLN);
        constrained = online > 0 &&
                      cores.size() < static_cast<std::size_t>(online);
        return {};
      }
      CPU_FREE(mask);
      if (saved_errno != EINVAL) {
        return invalid(kCpuAffinity);
      }
      capacity *= 2U;
    }
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, kCpuAffinity};
  } catch (...) {
    return invalid(kCpuAffinity);
  }
}

bool parse_nonnegative_integer(std::string_view text,
                               std::int32_t& value) noexcept {
  if (text.empty()) {
    return false;
  }
  std::uint64_t parsed = 0U;
  for (const char digit : text) {
    if (digit < '0' || digit > '9') {
      return false;
    }
    parsed = parsed * 10U + static_cast<unsigned>(digit - '0');
    if (parsed > static_cast<std::uint64_t>(
                     std::numeric_limits<std::int32_t>::max())) {
      return false;
    }
  }
  value = static_cast<std::int32_t>(parsed);
  return true;
}

Status discover_core_numa_nodes(
    Span<const std::int32_t> allowed,
    std::vector<std::int32_t>& out) noexcept {
  try {
    std::vector<std::int32_t> nodes(allowed.size, 0);
    for (std::size_t index = 0; index < allowed.size; ++index) {
      const std::filesystem::path cpu_root =
          std::filesystem::path{"/sys/devices/system/cpu"} /
          ("cpu" + std::to_string(allowed.data[index]));
      std::error_code error;
      std::int32_t selected_node = std::numeric_limits<std::int32_t>::max();
      for (std::filesystem::directory_iterator iterator(cpu_root, error), end;
           !error && iterator != end; iterator.increment(error)) {
        const std::string name = iterator->path().filename().string();
        if (name.size() <= 4U || name.compare(0U, 4U, "node") != 0) {
          continue;
        }
        std::int32_t node = -1;
        if (parse_nonnegative_integer(
                std::string_view{name}.substr(4U), node)) {
          selected_node = std::min(selected_node, node);
        }
      }
      if (selected_node != std::numeric_limits<std::int32_t>::max()) {
        nodes[index] = selected_node;
      }
    }
    out = std::move(nodes);
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, kCpuAffinity};
  } catch (...) {
    return invalid(kCpuAffinity);
  }
}

struct SharedCorePlacement {
  std::uint32_t rank{};
  std::uint32_t count{};
  std::uint32_t numa_node_count{1U};
  std::vector<std::int32_t> default_cores;
};

struct CoreCandidate {
  std::int32_t core{};
  std::int32_t rank{};
};

Status discover_shared_core_placement(
    MPI_Comm communicator, Span<const std::int32_t> allowed,
    Span<const std::int32_t> local_numa_nodes,
    SharedCorePlacement& out) noexcept {
  if (allowed.size != local_numa_nodes.size) {
    return invalid(kCpuAffinity);
  }
  MPI_Comm local = MPI_COMM_NULL;
  if (MPI_Comm_split_type(communicator, MPI_COMM_TYPE_SHARED, 0,
                          MPI_INFO_NULL, &local) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kCpuCommunicator};
  }
  int local_rank = 0;
  int local_size = 0;
  const int rank_result = MPI_Comm_rank(local, &local_rank);
  const int size_result = MPI_Comm_size(local, &local_size);
  if (rank_result != MPI_SUCCESS || size_result != MPI_SUCCESS ||
      local_rank < 0 || local_size <= 0 ||
      local_rank >= local_size) {
    MPI_Comm_free(&local);
    return {StatusCode::mpi_failure, kCpuCommunicator};
  }

  SharedCorePlacement candidate;
  candidate.rank = static_cast<std::uint32_t>(local_rank);
  candidate.count = static_cast<std::uint32_t>(local_size);
  std::vector<std::uint64_t> wide_counts;
  Status local_status;
  try {
    wide_counts.resize(static_cast<std::size_t>(local_size));
  } catch (...) {
    local_status = {StatusCode::allocation_failure, kCpuAffinity};
  }
  Status stage_status = consensus_status(local, local_status);
  if (!stage_status) {
    MPI_Comm_free(&local);
    return stage_status;
  }

  const std::uint64_t local_count = allowed.size;
  if (MPI_Allgather(&local_count, 1, MPI_UINT64_T, wide_counts.data(), 1,
                    MPI_UINT64_T, local) != MPI_SUCCESS) {
    MPI_Comm_free(&local);
    return {StatusCode::mpi_failure, kCpuCommunicator};
  }

  std::uint64_t wide_total = 0U;
  bool valid_counts = true;
  for (const std::uint64_t count : wide_counts) {
    valid_counts = valid_counts &&
                   count <= static_cast<std::uint64_t>(
                                std::numeric_limits<int>::max()) &&
                   wide_total <=
                       static_cast<std::uint64_t>(
                           std::numeric_limits<int>::max()) - count;
    if (!valid_counts) {
      break;
    }
    wide_total += count;
  }
  local_status = valid_counts ? Status{} : invalid(kCpuAffinity);
  stage_status = consensus_status(local, local_status);
  if (!stage_status) {
    MPI_Comm_free(&local);
    return stage_status;
  }

  std::vector<int> counts;
  std::vector<int> displacements;
  std::vector<std::int32_t> gathered;
  std::vector<std::int32_t> gathered_numa;
  std::vector<CoreCandidate> candidates;
  std::vector<std::size_t> owner_load;
  local_status = {};
  try {
    const std::size_t rank_count = static_cast<std::size_t>(local_size);
    const std::size_t total = static_cast<std::size_t>(wide_total);
    counts.resize(rank_count);
    displacements.resize(rank_count);
    gathered.resize(total);
    gathered_numa.resize(total);
    candidates.resize(total);
    owner_load.assign(rank_count, 0U);
    candidate.default_cores.reserve(allowed.size);
  } catch (...) {
    local_status = {StatusCode::allocation_failure, kCpuAffinity};
  }
  stage_status = consensus_status(local, local_status);
  if (!stage_status) {
    MPI_Comm_free(&local);
    return stage_status;
  }

  int displacement = 0;
  for (int rank = 0; rank < local_size; ++rank) {
    counts[static_cast<std::size_t>(rank)] =
        static_cast<int>(wide_counts[static_cast<std::size_t>(rank)]);
    displacements[static_cast<std::size_t>(rank)] = displacement;
    displacement += counts[static_cast<std::size_t>(rank)];
  }
  if (MPI_Allgatherv(allowed.data, static_cast<int>(allowed.size),
                     MPI_INT32_T, gathered.data(), counts.data(),
                     displacements.data(), MPI_INT32_T, local) != MPI_SUCCESS) {
    MPI_Comm_free(&local);
    return {StatusCode::mpi_failure, kCpuCommunicator};
  }
  if (MPI_Allgatherv(local_numa_nodes.data,
                     static_cast<int>(local_numa_nodes.size), MPI_INT32_T,
                     gathered_numa.data(), counts.data(), displacements.data(),
                     MPI_INT32_T, local) != MPI_SUCCESS) {
    MPI_Comm_free(&local);
    return {StatusCode::mpi_failure, kCpuCommunicator};
  }

  std::size_t numa_union_count = 0U;
  for (std::size_t index = 0; index < gathered_numa.size(); ++index) {
    bool seen = false;
    for (std::size_t prior = 0; prior < index; ++prior) {
      if (gathered_numa[prior] == gathered_numa[index]) {
        seen = true;
        break;
      }
    }
    if (!seen) {
      ++numa_union_count;
    }
  }
  if (numa_union_count != 0U) {
    candidate.numa_node_count = static_cast<std::uint32_t>(numa_union_count);
  }

  for (int rank = 0; rank < local_size; ++rank) {
    const int begin = displacements[static_cast<std::size_t>(rank)];
    const int end = begin + counts[static_cast<std::size_t>(rank)];
    for (int offset = begin; offset < end; ++offset) {
      candidates[static_cast<std::size_t>(offset)] =
          {gathered[static_cast<std::size_t>(offset)], rank};
    }
  }

  bool identical_masks = true;
  for (int rank = 1; rank < local_size && identical_masks; ++rank) {
    const std::size_t current = static_cast<std::size_t>(rank);
    identical_masks = counts[current] == counts[0];
    for (int offset = 0; offset < counts[0] && identical_masks; ++offset) {
      identical_masks =
          gathered[static_cast<std::size_t>(displacements[current] + offset)] ==
          gathered[static_cast<std::size_t>(offset)];
    }
  }
  if (identical_masks) {
    const std::size_t core_count = static_cast<std::size_t>(counts[0]);
    for (int offset = 0; offset < counts[0]; ++offset) {
      candidates[static_cast<std::size_t>(offset)] = {
          gathered[static_cast<std::size_t>(offset)],
          gathered_numa[static_cast<std::size_t>(offset)]};
    }
    std::sort(candidates.begin(), candidates.begin() +
                                        static_cast<std::ptrdiff_t>(core_count),
              [](CoreCandidate left, CoreCandidate right) noexcept {
                return left.rank < right.rank ||
                       (left.rank == right.rank && left.core < right.core);
              });
    std::size_t node_count = 0U;
    for (std::size_t index = 0; index < core_count; ++index) {
      if (index == 0U || candidates[index].rank != candidates[index - 1U].rank) {
        ++node_count;
      }
    }
    const std::size_t rank_count = static_cast<std::size_t>(local_size);
    const std::size_t rank = static_cast<std::size_t>(local_rank);
    const std::size_t active_node_count = std::min(node_count, rank_count);
    const std::size_t node_quotient = rank_count / active_node_count;
    const std::size_t node_remainder = rank_count % active_node_count;
    std::size_t selected_node_index = 0U;
    std::size_t first_rank_for_node = 0U;
    std::size_t ranks_for_node = 0U;
    for (std::size_t node = 0U; node < active_node_count; ++node) {
      const std::size_t group =
          node_quotient + (node < node_remainder ? 1U : 0U);
      if (rank >= first_rank_for_node && rank < first_rank_for_node + group) {
        selected_node_index = node;
        ranks_for_node = group;
        break;
      }
      first_rank_for_node += group;
    }
    std::size_t node_begin = 0U;
    for (std::size_t node = 0U; node < selected_node_index; ++node) {
      const std::int32_t label = candidates[node_begin].rank;
      while (node_begin < core_count && candidates[node_begin].rank == label) {
        ++node_begin;
      }
    }
    std::size_t node_end = node_begin;
    const std::int32_t node_label = candidates[node_begin].rank;
    while (node_end < core_count && candidates[node_end].rank == node_label) {
      ++node_end;
    }
    const std::size_t node_core_count = node_end - node_begin;
    const std::size_t rank_in_node = rank - first_rank_for_node;
    const std::size_t quotient = node_core_count / ranks_for_node;
    const std::size_t remainder = node_core_count % ranks_for_node;
    const std::size_t count =
        quotient + (rank_in_node < remainder ? 1U : 0U);
    const std::size_t begin_index =
        node_begin + rank_in_node * quotient +
        std::min(rank_in_node, remainder);
    local_status = {};
    try {
      candidate.default_cores.reserve(count);
      for (std::size_t index = begin_index; index < begin_index + count;
           ++index) {
        candidate.default_cores.push_back(candidates[index].core);
      }
    } catch (...) {
      local_status = {StatusCode::allocation_failure, kCpuAffinity};
    }
    stage_status = consensus_status(local, local_status);
    if (!stage_status) {
      MPI_Comm_free(&local);
      return stage_status;
    }
    const int free_result = MPI_Comm_free(&local);
    if (free_result != MPI_SUCCESS) {
      return {StatusCode::mpi_failure, kCpuCommunicator};
    }
    out = std::move(candidate);
    return {};
  }

  local_status = {};
  try {
    std::sort(candidates.begin(), candidates.end(),
              [](CoreCandidate left, CoreCandidate right) noexcept {
                return left.core < right.core ||
                       (left.core == right.core && left.rank < right.rank);
              });
    std::size_t begin = 0U;
    while (begin < candidates.size()) {
      std::size_t end = begin + 1U;
      while (end < candidates.size() &&
             candidates[end].core == candidates[begin].core) {
        ++end;
      }
      std::int32_t owner = candidates[begin].rank;
      for (std::size_t index = begin + 1U; index < end; ++index) {
        const std::int32_t eligible = candidates[index].rank;
        const std::size_t owner_index = static_cast<std::size_t>(owner);
        const std::size_t eligible_index = static_cast<std::size_t>(eligible);
        if (owner_load[eligible_index] < owner_load[owner_index] ||
            (owner_load[eligible_index] == owner_load[owner_index] &&
             eligible < owner)) {
          owner = eligible;
        }
      }
      ++owner_load[static_cast<std::size_t>(owner)];
      if (owner == local_rank) {
        candidate.default_cores.push_back(candidates[begin].core);
      }
      begin = end;
    }
  } catch (...) {
    local_status = {StatusCode::allocation_failure, kCpuAffinity};
  }
  stage_status = consensus_status(local, local_status);
  if (!stage_status) {
    MPI_Comm_free(&local);
    return stage_status;
  }

  const int free_result = MPI_Comm_free(&local);
  if (free_result != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kCpuCommunicator};
  }
  out = std::move(candidate);
  return {};
}

CpuKernelVariant compile_kernel_variant() noexcept {
#if (defined(__x86_64__) || defined(__i386__)) && \
    (defined(__GNUC__) || defined(__clang__))
  // This is a cold compile-time decision for the execution plan.  Kernels
  // dispatch only on the frozen finite enum and never query CPUID themselves.
  __builtin_cpu_init();
  if (__builtin_cpu_supports("avx512f") != 0) {
    return CpuKernelVariant::avx512;
  }
  if (__builtin_cpu_supports("avx2") != 0) {
    return CpuKernelVariant::avx2;
  }
  return CpuKernelVariant::scalar;
#else
  return CpuKernelVariant::scalar;
#endif
}

Int3 tile_for(CpuKernelVariant variant) noexcept {
  switch (variant) {
    case CpuKernelVariant::scalar:
      return {32, 4, 4};
    case CpuKernelVariant::avx2:
      return {64, 4, 4};
    case CpuKernelVariant::avx512:
      return {128, 4, 2};
  }
  return {32, 4, 4};
}

bool contains_core(Span<const std::int32_t> allowed,
                   std::int32_t core) noexcept {
  return std::binary_search(allowed.data, allowed.data + allowed.size, core);
}

Status shared_duplicate_warning(MPI_Comm communicator,
                                Span<const std::int32_t> selected,
                                bool& warning) noexcept {
  MPI_Comm local = MPI_COMM_NULL;
  if (MPI_Comm_split_type(communicator, MPI_COMM_TYPE_SHARED, 0,
                          MPI_INFO_NULL, &local) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kCpuCommunicator};
  }
  int local_size = 0;
  if (MPI_Comm_size(local, &local_size) != MPI_SUCCESS || local_size <= 0) {
    MPI_Comm_free(&local);
    return {StatusCode::mpi_failure, kCpuCommunicator};
  }
  std::vector<std::uint64_t> wide_counts;
  Status local_status;
  try {
    wide_counts.resize(static_cast<std::size_t>(local_size));
  } catch (...) {
    local_status = {StatusCode::allocation_failure, kCpuRequest};
  }
  Status stage = consensus_status(local, local_status);
  if (!stage) {
    MPI_Comm_free(&local);
    return stage;
  }
  const std::uint64_t selected_count = selected.size;
  if (MPI_Allgather(&selected_count, 1, MPI_UINT64_T, wide_counts.data(), 1,
                    MPI_UINT64_T, local) != MPI_SUCCESS) {
    MPI_Comm_free(&local);
    return {StatusCode::mpi_failure, kCpuCommunicator};
  }
  std::uint64_t wide_total = 0U;
  for (const std::uint64_t count : wide_counts) {
    if (count > static_cast<std::uint64_t>(std::numeric_limits<int>::max()) ||
        wide_total > static_cast<std::uint64_t>(
                         std::numeric_limits<int>::max()) - count) {
      local_status = invalid(kCpuRequest);
      break;
    }
    wide_total += count;
  }
  stage = consensus_status(local, local_status);
  if (!stage) {
    MPI_Comm_free(&local);
    return stage;
  }
  std::vector<int> counts;
  std::vector<int> displacements;
  std::vector<std::int32_t> gathered;
  try {
    counts.resize(static_cast<std::size_t>(local_size));
    displacements.resize(static_cast<std::size_t>(local_size));
    gathered.resize(static_cast<std::size_t>(wide_total));
  } catch (...) {
    local_status = {StatusCode::allocation_failure, kCpuRequest};
  }
  stage = consensus_status(local, local_status);
  if (!stage) {
    MPI_Comm_free(&local);
    return stage;
  }
  int displacement = 0;
  for (int rank = 0; rank < local_size; ++rank) {
    counts[static_cast<std::size_t>(rank)] =
        static_cast<int>(wide_counts[static_cast<std::size_t>(rank)]);
    displacements[static_cast<std::size_t>(rank)] = displacement;
    displacement += counts[static_cast<std::size_t>(rank)];
  }
  if (MPI_Allgatherv(selected.data, static_cast<int>(selected.size),
                     MPI_INT32_T, gathered.data(), counts.data(),
                     displacements.data(), MPI_INT32_T, local) != MPI_SUCCESS) {
    MPI_Comm_free(&local);
    return {StatusCode::mpi_failure, kCpuCommunicator};
  }
  std::sort(gathered.begin(), gathered.end());
  warning = std::adjacent_find(gathered.begin(), gathered.end()) !=
            gathered.end();
  if (MPI_Comm_free(&local) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kCpuCommunicator};
  }
  return {};
}

Status detect_duplicate_cores(Span<const std::int32_t> cores,
                              bool& duplicate) noexcept {
  try {
    std::vector<std::int32_t> sorted(cores.data, cores.data + cores.size);
    std::sort(sorted.begin(), sorted.end());
    duplicate = std::adjacent_find(sorted.begin(), sorted.end()) !=
                sorted.end();
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, kCpuRequest};
  } catch (...) {
    return invalid(kCpuRequest);
  }
}

Status bind_this_thread(std::int32_t core) noexcept {
  if (core < 0) {
    return invalid(kCpuCore);
  }
  const long configured = ::sysconf(_SC_NPROCESSORS_CONF);
  const std::size_t configured_capacity =
      configured > 0 ? static_cast<std::size_t>(configured) : CPU_SETSIZE;
  const std::size_t capacity =
      std::max({static_cast<std::size_t>(core) + 1U,
                configured_capacity, static_cast<std::size_t>(CPU_SETSIZE)});
  const std::size_t bytes = CPU_ALLOC_SIZE(capacity);
  cpu_set_t* const mask = CPU_ALLOC(capacity);
  if (mask == nullptr) {
    return {StatusCode::allocation_failure, kCpuBinding};
  }
  CPU_ZERO_S(bytes, mask);
  CPU_SET_S(static_cast<std::size_t>(core), bytes, mask);
  const int result = pthread_setaffinity_np(pthread_self(), bytes, mask);
  CPU_FREE(mask);
  return result == 0 ? Status{} : invalid(kCpuBinding);
}

bool should_inject_thread_launch_failure(std::size_t successful) noexcept {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  return successful >=
         g_fail_thread_launch_after.load(std::memory_order_acquire);
#else
  static_cast<void>(successful);
  return false;
#endif
}

}  // namespace

Status CpuExecutionPlan::compile(MPI_Comm communicator,
                                 CpuExecutionRequest request,
                                 CpuExecutionPlan& out) noexcept {
  if (!valid_communicator(communicator)) {
    return invalid(kCpuCommunicator);
  }

  try {
    CpuExecutionPlan candidate;
    bool constrained = false;
    std::vector<std::int32_t> allowed;
    Status local_status = discover_allowed_cores(allowed, constrained);
    Status status = consensus_status(communicator, local_status);
    if (!status) {
      return status;
    }

    std::vector<std::int32_t> local_numa_nodes;
    local_status = discover_core_numa_nodes(
        {allowed.data(), allowed.size()}, local_numa_nodes);
    status = consensus_status(communicator, local_status);
    if (!status) {
      return status;
    }
    local_status = allowed.size() <= static_cast<std::size_t>(
                                        std::numeric_limits<std::uint32_t>::max())
                       ? Status{}
                       : invalid(kCpuAffinity);
    status = consensus_status(communicator, local_status);
    if (!status) {
      return status;
    }
    candidate.diagnostics_.allowed_core_count =
        static_cast<std::uint32_t>(allowed.size());
    candidate.diagnostics_.numa_node_count = 1U;
    candidate.diagnostics_.local_rank_count = 1U;
    candidate.diagnostics_.recommended_local_rank_count = 1U;
    candidate.diagnostics_.affinity_constrained = constrained;
    candidate.kernel_variant_ = compile_kernel_variant();
    candidate.tile_shape_ = tile_for(candidate.kernel_variant_);

    SharedCorePlacement shared_placement;
    status = discover_shared_core_placement(
        communicator, {allowed.data(), allowed.size()},
        {local_numa_nodes.data(), local_numa_nodes.size()}, shared_placement);
    status = consensus_status(communicator, status);
    if (!status) {
      return status;
    }
    candidate.diagnostics_.local_rank_count = shared_placement.count;
    candidate.diagnostics_.numa_node_count = shared_placement.numa_node_count;
    candidate.diagnostics_.recommended_local_rank_count =
        shared_placement.numa_node_count;

    try {
      if (request.pure_mpi) {
        if (request.core_ids.size != 0U) {
          local_status = invalid(kCpuRequest);
        } else {
          candidate.core_ids_.push_back(allowed.front());
          candidate.pure_mpi_ = true;
          candidate.bind_threads_ = false;
        }
      } else {
        const std::size_t requested_threads =
            request.threads_per_rank == 0U
                ? (request.core_ids.size == 0U
                       ? shared_placement.default_cores.size()
                       : request.core_ids.size)
                : static_cast<std::size_t>(request.threads_per_rank);
        if (requested_threads == 0U ||
            requested_threads > static_cast<std::size_t>(
                                    std::numeric_limits<std::uint32_t>::max())) {
          local_status = invalid(kCpuRequest);
        }

        if (local_status) {
          if (request.core_ids.size != 0U) {
            if (request.core_ids.data == nullptr ||
                request.core_ids.size != requested_threads) {
              local_status = invalid(kCpuRequest);
            }
            if (local_status) {
              for (std::size_t index = 0; index < request.core_ids.size;
                   ++index) {
                if (!contains_core({allowed.data(), allowed.size()},
                                   request.core_ids.data[index])) {
                  local_status = invalid(kCpuCore);
                  break;
                }
              }
            }
            if (local_status) {
              candidate.core_ids_.assign(
                  request.core_ids.data,
                  request.core_ids.data + request.core_ids.size);
            }
          } else {
            if (requested_threads > shared_placement.default_cores.size() ||
                shared_placement.default_cores.empty()) {
              local_status = invalid(kCpuRequest);
            } else {
              candidate.core_ids_.assign(
                  shared_placement.default_cores.begin(),
                  shared_placement.default_cores.begin() +
                      static_cast<std::ptrdiff_t>(requested_threads));
            }
          }
        }
        if (local_status) {
          local_status = detect_duplicate_cores(
              {candidate.core_ids_.data(), candidate.core_ids_.size()},
              candidate.diagnostics_.duplicate_core_warning);
        }
        if (local_status) {
          candidate.pure_mpi_ = false;
          candidate.bind_threads_ = request.bind_threads;
        }
      }
    } catch (const std::bad_alloc&) {
      local_status = {StatusCode::allocation_failure, kCpuRequest};
    } catch (...) {
      local_status = invalid(kCpuRequest);
    }

    status = consensus_status(communicator, local_status);
    if (!status) {
      return status;
    }
    bool shared_warning = false;
    status = shared_duplicate_warning(
        communicator, {candidate.core_ids_.data(), candidate.core_ids_.size()},
        shared_warning);
    status = consensus_status(communicator, status);
    if (!status) {
      return status;
    }
    candidate.diagnostics_.duplicate_core_warning =
        candidate.diagnostics_.duplicate_core_warning || shared_warning;
    if (candidate.diagnostics_.duplicate_core_warning) {
      std::fputs("HUNDUN-FLOW warning: duplicate CPU core placement was "
                 "requested; workers will share a core.\n",
                 stderr);
      std::fflush(stderr);
    }
    out = std::move(candidate);
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, kCpuRequest};
  } catch (...) {
    return invalid(kCpuRequest);
  }
}

struct CpuThreadTeam::Impl {
  struct Worker {
    std::size_t index{};
    std::int32_t core{};
    std::thread thread;
  };

  std::vector<Worker> workers;
  std::mutex dispatch_mutex;
  std::mutex mutex;
  std::condition_variable work_ready;
  std::condition_variable work_complete;
  std::condition_variable startup_complete;
  CpuTask task{};
  void* context{};
  std::uint64_t generation{};
  std::size_t completed{};
  std::size_t startup_count{};
  Status startup_status{};
  bool stopping{};
  bool inline_mode{};

  void worker_loop(std::size_t worker_index, std::int32_t core,
                   bool bind) noexcept {
    const Status binding = bind ? bind_this_thread(core) : Status{};
    std::uint64_t observed_generation = 0U;
    {
      std::lock_guard<std::mutex> lock(mutex);
      if (!binding && startup_status) {
        startup_status = binding;
      }
      ++startup_count;
    }
    startup_complete.notify_one();

    std::unique_lock<std::mutex> lock(mutex);
    for (;;) {
      work_ready.wait(lock, [&] {
        return stopping || generation != observed_generation;
      });
      if (stopping) {
        return;
      }
      const CpuTask current_task = task;
      void* const current_context = context;
      observed_generation = generation;
      lock.unlock();
      current_task(worker_index, current_context);
      lock.lock();
      ++completed;
      if (completed == workers.size()) {
        work_complete.notify_one();
      }
    }
  }

  void stop_and_join() noexcept {
    {
      std::lock_guard<std::mutex> lock(mutex);
      stopping = true;
    }
    work_ready.notify_all();
    for (Worker& worker : workers) {
      if (worker.thread.joinable()) {
        worker.thread.join();
      }
    }
  }
};

CpuThreadTeam::~CpuThreadTeam() noexcept { release(); }

CpuThreadTeam::CpuThreadTeam(CpuThreadTeam&& other) noexcept
    : implementation_(other.implementation_) {
  other.implementation_ = nullptr;
}

CpuThreadTeam& CpuThreadTeam::operator=(CpuThreadTeam&& other) noexcept {
  if (this != &other) {
    release();
    implementation_ = other.implementation_;
    other.implementation_ = nullptr;
  }
  return *this;
}

void CpuThreadTeam::release() noexcept {
  if (implementation_ != nullptr) {
    implementation_->stop_and_join();
    delete implementation_;
    implementation_ = nullptr;
  }
}

Status CpuThreadTeam::create(const CpuExecutionPlan& plan,
                             CpuThreadTeam& out) noexcept {
  if (plan.core_ids_.empty()) {
    return invalid(kCpuRequest);
  }

  std::unique_ptr<Impl> candidate;
  try {
    candidate = std::make_unique<Impl>();
    candidate->inline_mode = plan.pure_mpi_;
    if (!candidate->inline_mode) {
      candidate->workers.resize(plan.core_ids_.size());
      for (std::size_t index = 0; index < candidate->workers.size(); ++index) {
        candidate->workers[index].index = index;
        candidate->workers[index].core = plan.core_ids_[index];
      }
      std::size_t launched = 0U;
      for (; launched < candidate->workers.size(); ++launched) {
        if (should_inject_thread_launch_failure(launched)) {
          throw std::bad_alloc{};
        }
        Impl::Worker& worker = candidate->workers[launched];
        worker.thread = std::thread(
            [implementation = candidate.get(), index = worker.index,
             core = worker.core, bind = plan.bind_threads_] {
              implementation->worker_loop(index, core, bind);
            });
      }

      std::unique_lock<std::mutex> lock(candidate->mutex);
      candidate->startup_complete.wait(lock, [&] {
        return candidate->startup_count == candidate->workers.size();
      });
      if (!candidate->startup_status) {
        lock.unlock();
        candidate->stop_and_join();
        return candidate->startup_status;
      }
    }
  } catch (const std::bad_alloc&) {
    if (candidate != nullptr) {
      candidate->stop_and_join();
    }
    return {StatusCode::allocation_failure, kCpuThread};
  } catch (const std::system_error&) {
    if (candidate != nullptr) {
      candidate->stop_and_join();
    }
    return {StatusCode::allocation_failure, kCpuThread};
  } catch (...) {
    if (candidate != nullptr) {
      candidate->stop_and_join();
    }
    return invalid(kCpuThread);
  }

  CpuThreadTeam replacement(candidate.release());
  out = std::move(replacement);
  return {};
}

Status CpuThreadTeam::run(CpuTask task, void* context) noexcept {
  if (implementation_ == nullptr || task == nullptr) {
    return invalid(kCpuThread);
  }
  if (implementation_->inline_mode) {
    std::lock_guard<std::mutex> dispatch_lock(
        implementation_->dispatch_mutex);
    task(0U, context);
    return {};
  }

  std::lock_guard<std::mutex> dispatch_lock(implementation_->dispatch_mutex);
  std::unique_lock<std::mutex> lock(implementation_->mutex);
  if (implementation_->stopping) {
    return invalid(kCpuThread);
  }
  implementation_->task = task;
  implementation_->context = context;
  implementation_->completed = 0U;
  ++implementation_->generation;
  implementation_->work_ready.notify_all();
  implementation_->work_complete.wait(lock, [&] {
    return implementation_->completed == implementation_->workers.size();
  });
  return {};
}

std::size_t CpuThreadTeam::worker_count() const noexcept {
  if (implementation_ == nullptr) {
    return 0U;
  }
  return implementation_->inline_mode ? 1U
                                      : implementation_->workers.size();
}

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
namespace detail {

void fail_cpu_thread_launch_after_for_test(
    std::size_t successful_launches) noexcept {
  g_fail_thread_launch_after.store(successful_launches,
                                   std::memory_order_release);
}

void reset_cpu_thread_launch_failure_for_test() noexcept {
  g_fail_thread_launch_after.store(std::numeric_limits<std::size_t>::max(),
                                   std::memory_order_release);
}

}  // namespace detail
#endif

}  // namespace hundun::v04

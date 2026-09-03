// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_linear.hpp"

#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <new>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

constexpr std::size_t kPacketWords = 6U;
constexpr std::size_t kPacketSize = 0U;
constexpr std::size_t kReferenceRank = 1U;
constexpr std::size_t kReferenceCount = 2U;
constexpr std::size_t kFailureRank = 3U;
constexpr std::size_t kFailureCode = 4U;
constexpr std::size_t kFailureDetail = 5U;
constexpr int kGatherTag = 17041;
constexpr int kScatterTag = 17042;
constexpr int kMaximumLocationGatherTag = 17043;
constexpr int kMaximumLocationScatterTag = 17044;

constexpr std::uint32_t kInvalidEngine = 11001U;
constexpr std::uint32_t kInvalidCapacity = 11002U;
constexpr std::uint32_t kInvalidMode = 11003U;
constexpr std::uint32_t kCompileCollective = 11004U;
constexpr std::uint32_t kCompileMismatch = 11005U;
constexpr std::uint32_t kAllocation = 11006U;
constexpr std::uint32_t kMpiDatatype = 11007U;
constexpr std::uint32_t kMpiOperation = 11008U;
constexpr std::uint32_t kInvalidStatus = 11009U;
constexpr std::uint32_t kInvalidSpan = 11010U;
constexpr std::uint32_t kSpanMismatch = 11011U;
constexpr std::uint32_t kSpanCapacity = 11012U;
constexpr std::uint32_t kSpanOverlap = 11013U;
constexpr std::uint32_t kNonFiniteLocal = 11014U;
constexpr std::uint32_t kCounterOverflow = 11015U;
constexpr std::uint32_t kCollectiveFailure = 11016U;
constexpr std::uint32_t kPacketFailure = 11017U;
constexpr std::uint32_t kCountMismatch = 11018U;
constexpr std::uint32_t kNonFiniteGlobal = 11019U;
constexpr std::uint32_t kContractMismatch = 11020U;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
constexpr std::uint32_t kCheckedSumFaultForTest = 11021U;
std::uint64_t g_checked_sum_fault_ordinal = 0U;
std::uint64_t g_checked_sum_next_ordinal = 0U;
int g_checked_sum_fault_rank = -1;
bool g_checked_sum_fault_armed = false;
#endif

bool valid_mode(ReductionMode mode) noexcept {
  return mode == ReductionMode::reproducible_tree ||
         mode == ReductionMode::mpi_allreduce;
}

bool valid_status_code(StatusCode code) noexcept {
  return static_cast<std::uint16_t>(code) <=
         static_cast<std::uint16_t>(StatusCode::io_failure);
}

bool checked_add(std::uint64_t left, std::uint64_t right,
                 std::uint64_t& result) noexcept {
  if (left > std::numeric_limits<std::uint64_t>::max() - right) {
    return false;
  }
  result = left + right;
  return true;
}

bool checked_multiply(std::uint64_t left, std::uint64_t right,
                      std::uint64_t& result) noexcept {
  if (left != 0U &&
      right > std::numeric_limits<std::uint64_t>::max() / left) {
    return false;
  }
  result = left * right;
  return true;
}

bool spans_overlap(const double* input, double* output,
                   std::size_t count) noexcept {
  if (input == nullptr || output == nullptr || count == 0U ||
      count > std::numeric_limits<std::size_t>::max() / sizeof(double)) {
    return input == nullptr || output == nullptr;
  }
  const std::size_t bytes = count * sizeof(double);
  const auto input_begin = reinterpret_cast<std::uintptr_t>(input);
  const auto output_begin = reinterpret_cast<std::uintptr_t>(output);
  if (input_begin > std::numeric_limits<std::uintptr_t>::max() - bytes ||
      output_begin > std::numeric_limits<std::uintptr_t>::max() - bytes) {
    return true;
  }
  const std::uintptr_t input_end = input_begin + bytes;
  const std::uintptr_t output_end = output_begin + bytes;
  return input_begin < output_end && output_begin < input_end;
}

void select_failure(double rank, double code, double detail,
                    double* packet) noexcept {
  const double selected_rank = packet[kFailureRank];
  const bool earlier_rank = rank < selected_rank;
  const bool proposed_count_mismatch =
      code == static_cast<double>(StatusCode::invalid_plan) &&
      detail == static_cast<double>(kCountMismatch);
  const bool selected_count_mismatch =
      packet[kFailureCode] ==
          static_cast<double>(StatusCode::invalid_plan) &&
      packet[kFailureDetail] == static_cast<double>(kCountMismatch);
  const bool same_rank_better = rank == selected_rank &&
                                ((selected_count_mismatch &&
                                  !proposed_count_mismatch) ||
                                 (selected_count_mismatch ==
                                      proposed_count_mismatch &&
                                  (code < packet[kFailureCode] ||
                                   (code == packet[kFailureCode] &&
                                    detail < packet[kFailureDetail]))));
  if (earlier_rank || same_rank_better) {
    packet[kFailureRank] = rank;
    packet[kFailureCode] = code;
    packet[kFailureDetail] = detail;
  }
}

void combine_packets(const double* input, double* inout,
                     bool sum) noexcept {
  const std::size_t packet_size =
      static_cast<std::size_t>(inout[kPacketSize]);

  select_failure(input[kFailureRank], input[kFailureCode],
                 input[kFailureDetail], inout);

  const double input_rank = input[kReferenceRank];
  const double inout_rank = inout[kReferenceRank];
  const double input_count = input[kReferenceCount];
  const double inout_count = inout[kReferenceCount];
  if (input_rank < inout_rank) {
    if (input_count != inout_count) {
      select_failure(inout_rank,
                     static_cast<double>(StatusCode::invalid_plan),
                     static_cast<double>(kCountMismatch), inout);
    }
    inout[kReferenceRank] = input[kReferenceRank];
    inout[kReferenceCount] = input[kReferenceCount];
  } else if (input_count != inout_count) {
    select_failure(input_rank,
                   static_cast<double>(StatusCode::invalid_plan),
                   static_cast<double>(kCountMismatch), inout);
  }

  for (std::size_t index = kPacketWords; index < packet_size; ++index) {
    inout[index] = sum ? inout[index] + input[index]
                       : std::fmax(inout[index], input[index]);
  }
}

void sum_packet(void* input, void* inout, int* length,
                MPI_Datatype* datatype) noexcept {
  (void)datatype;
  if (length != nullptr && *length == 1 && input != nullptr &&
      inout != nullptr) {
    combine_packets(static_cast<const double*>(input),
                    static_cast<double*>(inout), true);
  }
}

void max_packet(void* input, void* inout, int* length,
                MPI_Datatype* datatype) noexcept {
  (void)datatype;
  if (length != nullptr && *length == 1 && input != nullptr &&
      inout != nullptr) {
    combine_packets(static_cast<const double*>(input),
                    static_cast<double*>(inout), false);
  }
}

bool better_maximum_location(const ReductionMaximumLocation& candidate,
                             const ReductionMaximumLocation& selected) {
  return candidate.valid &&
         (!selected.valid || candidate.value > selected.value ||
          (candidate.value == selected.value &&
           (candidate.global_location < selected.global_location ||
            (candidate.global_location == selected.global_location &&
             candidate.rank < selected.rank))));
}

void combine_maximum_locations(const ReductionMaximumLocation* input,
                               ReductionMaximumLocation* inout,
                               int count) noexcept {
  if (input == nullptr || inout == nullptr || count <= 0) return;
  for (int index = 0; index < count; ++index)
    if (better_maximum_location(input[index], inout[index]))
      inout[index] = input[index];
}

void maximum_location_op(void* input, void* inout, int* length,
                         MPI_Datatype* datatype) noexcept {
  (void)datatype;
  if (length != nullptr)
    combine_maximum_locations(
        static_cast<const ReductionMaximumLocation*>(input),
        static_cast<ReductionMaximumLocation*>(inout), *length);
}

Status raw_consensus(MPI_Comm communicator, int rank, int size, Status local,
                     int& lowest) noexcept {
  if (!valid_status_code(local.code)) {
    local = {StatusCode::invalid_plan, kInvalidStatus};
  }
  const int candidate = local ? size : rank;
  int selected = size;
  if (MPI_Allreduce(&candidate, &selected, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    lowest = -1;
    return {StatusCode::mpi_failure, kCompileCollective};
  }
  if (selected == size) {
    lowest = -1;
    return {};
  }
  std::uint64_t wire[2]{};
  if (rank == selected) {
    wire[0] = static_cast<std::uint64_t>(local.code);
    wire[1] = local.detail;
  }
  if (MPI_Bcast(wire, 2, MPI_UINT64_T, selected, communicator) !=
          MPI_SUCCESS ||
      wire[0] > static_cast<std::uint64_t>(StatusCode::io_failure) ||
      wire[1] > std::numeric_limits<std::uint32_t>::max()) {
    lowest = -1;
    return {StatusCode::mpi_failure, kCompileCollective};
  }
  lowest = selected;
  return {static_cast<StatusCode>(wire[0]),
          static_cast<std::uint32_t>(wire[1])};
}

bool mpi_is_live() noexcept {
  int initialized = 0;
  int finalized = 0;
  (void)MPI_Initialized(&initialized);
  if (initialized != 0) {
    (void)MPI_Finalized(&finalized);
  }
  return initialized != 0 && finalized == 0;
}

}  // namespace

struct ReductionEngine::Impl {
  MPI_Comm communicator{MPI_COMM_NULL};
  MPI_Datatype packet_type{MPI_DATATYPE_NULL};
  MPI_Op sum_operation{MPI_OP_NULL};
  MPI_Op max_operation{MPI_OP_NULL};
  MPI_Datatype maximum_location_type{MPI_DATATYPE_NULL};
  MPI_Op maximum_location_operation{MPI_OP_NULL};
  int rank{};
  int size{};
  ReductionMode reduction_mode{ReductionMode::mpi_allreduce};
  std::size_t maximum_scalars{};
  std::vector<double> send_storage;
  std::vector<double> receive_storage;
  std::vector<ReductionMaximumLocation> maximum_location_send_storage;
  std::vector<ReductionMaximumLocation> maximum_location_receive_storage;
  LinearReductionCounters reduction_counters{};
  int lowest_failing_rank{-1};
};

namespace {

template <class Implementation>
void destroy_impl(Implementation* implementation) noexcept {
  if (implementation == nullptr) {
    return;
  }
  if (mpi_is_live()) {
    if (implementation->sum_operation != MPI_OP_NULL) {
      (void)MPI_Op_free(&implementation->sum_operation);
    }
    if (implementation->max_operation != MPI_OP_NULL) {
      (void)MPI_Op_free(&implementation->max_operation);
    }
    if (implementation->maximum_location_operation != MPI_OP_NULL) {
      (void)MPI_Op_free(&implementation->maximum_location_operation);
    }
    if (implementation->maximum_location_type != MPI_DATATYPE_NULL) {
      (void)MPI_Type_free(&implementation->maximum_location_type);
    }
    if (implementation->packet_type != MPI_DATATYPE_NULL) {
      (void)MPI_Type_free(&implementation->packet_type);
    }
    if (implementation->communicator != MPI_COMM_NULL) {
      (void)MPI_Comm_free(&implementation->communicator);
    }
  }
  delete implementation;
}

template <class Implementation>
void initialize_packet(Implementation& implementation,
                       Span<const double> local, Status local_status) noexcept {
  double* const packet = implementation.send_storage.data();
  const std::size_t packet_size = implementation.send_storage.size();
  std::fill(packet, packet + packet_size, 0.0);
  packet[kPacketSize] = static_cast<double>(packet_size);
  packet[kReferenceRank] = static_cast<double>(implementation.rank);
  packet[kReferenceCount] = static_cast<double>(local.size);
  packet[kFailureRank] = static_cast<double>(implementation.size);
  packet[kFailureCode] = static_cast<double>(StatusCode::ok);
  packet[kFailureDetail] = 0.0;
  if (!local_status) {
    packet[kFailureRank] = static_cast<double>(implementation.rank);
    packet[kFailureCode] = static_cast<double>(local_status.code);
    packet[kFailureDetail] = static_cast<double>(local_status.detail);
  }
  if (local_status && local.data != nullptr &&
      local.size <= implementation.maximum_scalars) {
    std::copy(local.data, local.data + local.size, packet + kPacketWords);
  }
}

template <class Implementation>
Status packet_collective(Implementation& implementation,
                         bool sum) noexcept {
  const auto begin = std::chrono::steady_clock::now();
  int mpi_status = MPI_SUCCESS;
  if (implementation.reduction_mode == ReductionMode::mpi_allreduce) {
    mpi_status = MPI_Allreduce(
        implementation.send_storage.data(),
        implementation.receive_storage.data(), 1, implementation.packet_type,
        sum ? implementation.sum_operation : implementation.max_operation,
        implementation.communicator);
  } else if (implementation.rank == 0) {
    std::copy(implementation.send_storage.begin(),
              implementation.send_storage.end(),
              implementation.receive_storage.begin());
    for (int source = 1; source < implementation.size; ++source) {
      if (MPI_Recv(implementation.send_storage.data(), 1,
                   implementation.packet_type, source, kGatherTag,
                   implementation.communicator, MPI_STATUS_IGNORE) !=
          MPI_SUCCESS) {
        mpi_status = MPI_ERR_OTHER;
        break;
      }
      combine_packets(implementation.send_storage.data(),
                      implementation.receive_storage.data(), sum);
    }
    if (mpi_status == MPI_SUCCESS) {
      for (int destination = 1; destination < implementation.size;
           ++destination) {
        if (MPI_Send(implementation.receive_storage.data(), 1,
                     implementation.packet_type, destination, kScatterTag,
                     implementation.communicator) != MPI_SUCCESS) {
          mpi_status = MPI_ERR_OTHER;
          break;
        }
      }
    }
  } else {
    if (MPI_Send(implementation.send_storage.data(), 1,
                 implementation.packet_type, 0, kGatherTag,
                 implementation.communicator) != MPI_SUCCESS ||
        MPI_Recv(implementation.receive_storage.data(), 1,
                 implementation.packet_type, 0, kScatterTag,
                 implementation.communicator, MPI_STATUS_IGNORE) !=
            MPI_SUCCESS) {
      mpi_status = MPI_ERR_OTHER;
    }
  }
  const auto elapsed = static_cast<std::uint64_t>(
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          std::chrono::steady_clock::now() - begin)
          .count());
  LinearReductionCounters& counters = implementation.reduction_counters;
  counters.blocking_operations =
      counters.blocking_operations == UINT64_MAX
          ? UINT64_MAX
          : counters.blocking_operations + 1U;
  counters.wall_nanoseconds =
      elapsed > UINT64_MAX - counters.wall_nanoseconds
          ? UINT64_MAX
          : counters.wall_nanoseconds + elapsed;
  if (mpi_status != MPI_SUCCESS) {
    implementation.lowest_failing_rank = -1;
    return {StatusCode::mpi_failure, kCollectiveFailure};
  }
  const double* const packet = implementation.receive_storage.data();
  if (packet[kPacketSize] !=
          static_cast<double>(implementation.receive_storage.size()) ||
      packet[kFailureRank] < 0.0 ||
      packet[kFailureRank] > static_cast<double>(implementation.size) ||
      packet[kFailureCode] < 0.0 ||
      packet[kFailureCode] >
          static_cast<double>(StatusCode::io_failure) ||
      packet[kFailureDetail] < 0.0 ||
      packet[kFailureDetail] >
          static_cast<double>(std::numeric_limits<std::uint32_t>::max())) {
    implementation.lowest_failing_rank = -1;
    return {StatusCode::mpi_failure, kPacketFailure};
  }
  if (packet[kFailureRank] < static_cast<double>(implementation.size)) {
    implementation.lowest_failing_rank =
        static_cast<int>(packet[kFailureRank]);
    return {static_cast<StatusCode>(
                static_cast<std::uint16_t>(packet[kFailureCode])),
            static_cast<std::uint32_t>(packet[kFailureDetail])};
  }
  implementation.lowest_failing_rank = -1;
  return {};
}

template <class Implementation>
Status validate_reduction(Implementation& implementation,
                          Span<const double> local, Span<double> global,
                          Status supplied, bool sum) noexcept {
  Status local_status = supplied;
  if (!valid_status_code(local_status.code)) {
    local_status = {StatusCode::invalid_plan, kInvalidStatus};
  }
  if (local_status) {
    if (local.size == 0U || global.size == 0U || local.data == nullptr ||
        global.data == nullptr) {
      local_status = {StatusCode::invalid_plan, kInvalidSpan};
    } else if (local.size != global.size) {
      local_status = {StatusCode::invalid_plan, kSpanMismatch};
    } else if (local.size > implementation.maximum_scalars) {
      local_status = {StatusCode::invalid_plan, kSpanCapacity};
    } else if (spans_overlap(local.data, global.data, local.size)) {
      local_status = {StatusCode::invalid_plan, kSpanOverlap};
    } else {
      for (std::size_t index = 0U; index < local.size; ++index) {
        if (!std::isfinite(local.data[index])) {
          local_status = {StatusCode::numerical_failure, kNonFiniteLocal};
          break;
        }
      }
    }
  }

  std::uint64_t next_calls = 0U;
  std::uint64_t next_scalars = 0U;
  std::uint64_t next_bytes = 0U;
  std::uint64_t next_messages = 0U;
  std::uint64_t logical_bytes = 0U;
  std::uint64_t tree_messages = 0U;
  if (local_status) {
    const auto scalar_count = static_cast<std::uint64_t>(local.size);
    const auto peer_count = static_cast<std::uint64_t>(
        implementation.size > 0 ? implementation.size - 1 : 0);
    if (!checked_multiply(scalar_count, sizeof(double), logical_bytes) ||
        !checked_multiply(peer_count, 2U, tree_messages) ||
        !checked_add(implementation.reduction_counters.calls, 1U,
                     next_calls) ||
        !checked_add(implementation.reduction_counters.scalars, scalar_count,
                     next_scalars) ||
        !checked_add(implementation.reduction_counters.logical_bytes,
                     logical_bytes, next_bytes) ||
        !checked_add(implementation.reduction_counters.tree_messages,
                     implementation.reduction_mode ==
                             ReductionMode::reproducible_tree
                         ? tree_messages
                         : 0U,
                     next_messages)) {
      local_status = {StatusCode::invalid_plan, kCounterOverflow};
    }
  }

  initialize_packet(implementation, local, local_status);
  const Status collective = packet_collective(implementation, sum);
  if (!collective) {
    return collective;
  }

  const double* const values =
      implementation.receive_storage.data() + kPacketWords;
  for (std::size_t index = 0U; index < local.size; ++index) {
    if (!std::isfinite(values[index])) {
      implementation.lowest_failing_rank = -1;
      return {StatusCode::numerical_failure, kNonFiniteGlobal};
    }
  }

  std::copy(values, values + local.size, global.data);
  implementation.reduction_counters.calls = next_calls;
  implementation.reduction_counters.scalars = next_scalars;
  implementation.reduction_counters.logical_bytes = next_bytes;
  implementation.reduction_counters.tree_messages = next_messages;
  return {};
}

}  // namespace

ReductionEngine::~ReductionEngine() noexcept { release(); }

ReductionEngine::ReductionEngine(ReductionEngine&& other) noexcept
    : implementation_(std::exchange(other.implementation_, nullptr)) {}

ReductionEngine& ReductionEngine::operator=(ReductionEngine&& other) noexcept {
  if (this != &other) {
    release();
    implementation_ = std::exchange(other.implementation_, nullptr);
  }
  return *this;
}

void ReductionEngine::release() noexcept {
  destroy_impl(std::exchange(implementation_, nullptr));
}

Status ReductionEngine::compile(MPI_Comm communicator, ReductionMode mode,
                                std::size_t maximum_scalars,
                                ReductionEngine& out) noexcept {
  if (!mpi_is_live() || communicator == MPI_COMM_NULL) {
    return {StatusCode::invalid_plan, kInvalidEngine};
  }

  MPI_Comm duplicate = MPI_COMM_NULL;
  if (MPI_Comm_dup(communicator, &duplicate) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kCompileCollective};
  }
  int rank = 0;
  int size = 0;
  if (MPI_Comm_rank(duplicate, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(duplicate, &size) != MPI_SUCCESS || size <= 0) {
    (void)MPI_Comm_free(&duplicate);
    return {StatusCode::mpi_failure, kCompileCollective};
  }
  if (MPI_Comm_set_errhandler(duplicate, MPI_ERRORS_RETURN) != MPI_SUCCESS) {
    (void)MPI_Comm_free(&duplicate);
    return {StatusCode::mpi_failure, kCompileCollective};
  }

  Status local{};
  if (maximum_scalars == 0U ||
      maximum_scalars >
          static_cast<std::size_t>(std::numeric_limits<int>::max()) -
              kPacketWords) {
    local = {StatusCode::invalid_plan, kInvalidCapacity};
  } else if (!valid_mode(mode)) {
    local = {StatusCode::invalid_plan, kInvalidMode};
  }

  const std::uint64_t local_capacity = maximum_scalars;
  const std::uint64_t local_mode = static_cast<std::uint64_t>(mode);
  std::uint64_t minimums[2]{local_capacity, local_mode};
  std::uint64_t maximums[2]{local_capacity, local_mode};
  std::uint64_t global_minimums[2]{};
  std::uint64_t global_maximums[2]{};
  const int minimum_status = MPI_Allreduce(
      minimums, global_minimums, 2, MPI_UINT64_T, MPI_MIN, duplicate);
  const int maximum_status = MPI_Allreduce(
      maximums, global_maximums, 2, MPI_UINT64_T, MPI_MAX, duplicate);
  if (minimum_status != MPI_SUCCESS || maximum_status != MPI_SUCCESS) {
    (void)MPI_Comm_free(&duplicate);
    return {StatusCode::mpi_failure, kCompileCollective};
  }
  if (local && (global_minimums[0] != global_maximums[0] ||
                global_minimums[1] != global_maximums[1])) {
    local = {StatusCode::invalid_plan, kCompileMismatch};
  }
  int ignored_lowest = -1;
  Status consensus = raw_consensus(duplicate, rank, size, local,
                                   ignored_lowest);
  if (!consensus) {
    (void)MPI_Comm_free(&duplicate);
    return consensus;
  }

  Impl* candidate = new (std::nothrow) Impl;
  local = candidate == nullptr
              ? Status{StatusCode::allocation_failure, kAllocation}
              : Status{};
  consensus = raw_consensus(duplicate, rank, size, local, ignored_lowest);
  if (!consensus) {
    delete candidate;
    (void)MPI_Comm_free(&duplicate);
    return consensus;
  }
  candidate->communicator = duplicate;
  duplicate = MPI_COMM_NULL;
  candidate->rank = rank;
  candidate->size = size;
  candidate->reduction_mode = mode;
  candidate->maximum_scalars = maximum_scalars;

  try {
    candidate->send_storage.resize(maximum_scalars + kPacketWords);
    candidate->receive_storage.resize(maximum_scalars + kPacketWords);
    candidate->maximum_location_send_storage.resize(maximum_scalars);
    candidate->maximum_location_receive_storage.resize(maximum_scalars);
  } catch (...) {
    local = {StatusCode::allocation_failure, kAllocation};
  }
  consensus = raw_consensus(candidate->communicator, rank, size, local,
                            ignored_lowest);
  if (!consensus) {
    destroy_impl(candidate);
    return consensus;
  }

  const int packet_words = static_cast<int>(maximum_scalars + kPacketWords);
  local = MPI_Type_contiguous(packet_words, MPI_DOUBLE,
                              &candidate->packet_type) == MPI_SUCCESS &&
                  MPI_Type_commit(&candidate->packet_type) == MPI_SUCCESS &&
                  MPI_Type_contiguous(
                      static_cast<int>(sizeof(ReductionMaximumLocation)),
                      MPI_BYTE, &candidate->maximum_location_type) ==
                      MPI_SUCCESS &&
                  MPI_Type_commit(&candidate->maximum_location_type) ==
                      MPI_SUCCESS
              ? Status{}
              : Status{StatusCode::mpi_failure, kMpiDatatype};
  consensus = raw_consensus(candidate->communicator, rank, size, local,
                            ignored_lowest);
  if (!consensus) {
    destroy_impl(candidate);
    return consensus;
  }

  local = MPI_Op_create(&sum_packet, 1, &candidate->sum_operation) ==
                      MPI_SUCCESS &&
                  MPI_Op_create(&max_packet, 1,
                                &candidate->max_operation) == MPI_SUCCESS &&
                  MPI_Op_create(&maximum_location_op, 1,
                                &candidate->maximum_location_operation) ==
                      MPI_SUCCESS
              ? Status{}
              : Status{StatusCode::mpi_failure, kMpiOperation};
  consensus = raw_consensus(candidate->communicator, rank, size, local,
                            ignored_lowest);
  if (!consensus) {
    destroy_impl(candidate);
    return consensus;
  }

  out.release();
  out.implementation_ = candidate;
  return {};
}

Status ReductionEngine::checked_sum(Span<const double> local,
                                    Span<double> global,
                                    Status local_status) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kInvalidEngine};
  }
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  ++g_checked_sum_next_ordinal;
  if (g_checked_sum_fault_armed &&
      g_checked_sum_next_ordinal == g_checked_sum_fault_ordinal) {
    const bool inject = implementation_->rank == g_checked_sum_fault_rank;
    g_checked_sum_fault_armed = false;
    if (inject) {
      local_status = {StatusCode::numerical_failure,
                      kCheckedSumFaultForTest};
    }
  }
#endif
  return validate_reduction(*implementation_, local, global, local_status,
                            true);
}

Status ReductionEngine::checked_max(Span<const double> local,
                                    Span<double> global,
                                    Status local_status) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kInvalidEngine};
  }
  return validate_reduction(*implementation_, local, global, local_status,
                            false);
}

Status ReductionEngine::checked_max_locations(
    Span<const ReductionMaximumLocation> local,
    Span<ReductionMaximumLocation> global, Status local_status) noexcept {
  if (implementation_ == nullptr)
    return {StatusCode::invalid_plan, kInvalidEngine};
  Impl& implementation = *implementation_;
  if (!valid_status_code(local_status.code))
    local_status = {StatusCode::invalid_plan, kInvalidStatus};
  if (local_status) {
    if (local.data == nullptr || global.data == nullptr || local.size == 0U)
      local_status = {StatusCode::invalid_plan, kInvalidSpan};
    else if (local.size != global.size)
      local_status = {StatusCode::invalid_plan, kSpanMismatch};
    else if (local.size > implementation.maximum_scalars ||
             local.size > static_cast<std::size_t>(
                              std::numeric_limits<int>::max()))
      local_status = {StatusCode::invalid_plan, kSpanCapacity};
    else {
      const std::size_t bytes =
          local.size * sizeof(ReductionMaximumLocation);
      const auto input = reinterpret_cast<std::uintptr_t>(local.data);
      const auto output = reinterpret_cast<std::uintptr_t>(global.data);
      if (input > std::numeric_limits<std::uintptr_t>::max() - bytes ||
          output > std::numeric_limits<std::uintptr_t>::max() - bytes ||
          (input < output + bytes && output < input + bytes))
        local_status = {StatusCode::invalid_plan, kSpanOverlap};
    }
  }
  for (std::size_t index = 0U; index < local.size && local_status; ++index) {
    const ReductionMaximumLocation& value = local.data[index];
    if (!value.valid) continue;
    if (!std::isfinite(value.value) || value.value < 0.0 ||
        value.global_location == UINT64_MAX ||
        value.rank != implementation.rank) {
      local_status = {StatusCode::numerical_failure, kNonFiniteLocal};
      break;
    }
    for (double payload : value.payload)
      if (!std::isfinite(payload)) {
        local_status = {StatusCode::numerical_failure, kNonFiniteLocal};
        break;
      }
  }

  const auto begin = std::chrono::steady_clock::now();
  int ignored_lowest = -1;
  const Status agreed = raw_consensus(implementation.communicator,
                                      implementation.rank,
                                      implementation.size, local_status,
                                      ignored_lowest);
  if (!agreed) return agreed;
  std::copy(local.data, local.data + local.size,
            implementation.maximum_location_send_storage.data());
  int mpi_status = MPI_SUCCESS;
  const int count = static_cast<int>(local.size);
  if (implementation.reduction_mode == ReductionMode::mpi_allreduce) {
    mpi_status = MPI_Allreduce(
        implementation.maximum_location_send_storage.data(),
        implementation.maximum_location_receive_storage.data(), count,
        implementation.maximum_location_type,
        implementation.maximum_location_operation,
        implementation.communicator);
  } else if (implementation.rank == 0) {
    std::copy(implementation.maximum_location_send_storage.begin(),
              implementation.maximum_location_send_storage.begin() + count,
              implementation.maximum_location_receive_storage.begin());
    for (int source = 1; source < implementation.size; ++source) {
      if (MPI_Recv(implementation.maximum_location_send_storage.data(), count,
                   implementation.maximum_location_type, source,
                   kMaximumLocationGatherTag, implementation.communicator,
                   MPI_STATUS_IGNORE) != MPI_SUCCESS) {
        mpi_status = MPI_ERR_OTHER;
        break;
      }
      combine_maximum_locations(
          implementation.maximum_location_send_storage.data(),
          implementation.maximum_location_receive_storage.data(), count);
    }
    if (mpi_status == MPI_SUCCESS) {
      for (int destination = 1; destination < implementation.size;
           ++destination) {
        if (MPI_Send(
                implementation.maximum_location_receive_storage.data(), count,
                implementation.maximum_location_type, destination,
                kMaximumLocationScatterTag, implementation.communicator) !=
            MPI_SUCCESS) {
          mpi_status = MPI_ERR_OTHER;
          break;
        }
      }
    }
  } else if (MPI_Send(
                         implementation.maximum_location_send_storage.data(),
                         count, implementation.maximum_location_type, 0,
                         kMaximumLocationGatherTag,
                         implementation.communicator) != MPI_SUCCESS ||
             MPI_Recv(
                 implementation.maximum_location_receive_storage.data(),
                 count, implementation.maximum_location_type, 0,
                 kMaximumLocationScatterTag, implementation.communicator,
                 MPI_STATUS_IGNORE) != MPI_SUCCESS) {
    mpi_status = MPI_ERR_OTHER;
  }
  const auto elapsed = static_cast<std::uint64_t>(
      std::chrono::duration_cast<std::chrono::nanoseconds>(
          std::chrono::steady_clock::now() - begin)
          .count());
  if (mpi_status != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kCollectiveFailure};
  for (std::size_t index = 0U; index < local.size; ++index) {
    const ReductionMaximumLocation& value =
        implementation.maximum_location_receive_storage[index];
    if (!value.valid || !std::isfinite(value.value) || value.value < 0.0 ||
        value.global_location == UINT64_MAX || value.rank < 0 ||
        value.rank >= implementation.size) {
      return {StatusCode::numerical_failure, kNonFiniteGlobal};
    }
    for (double payload : value.payload)
      if (!std::isfinite(payload))
        return {StatusCode::numerical_failure, kNonFiniteGlobal};
  }

  LinearReductionCounters candidate = implementation.reduction_counters;
  const std::uint64_t count_u64 = static_cast<std::uint64_t>(local.size);
  const std::uint64_t bytes = count_u64 * sizeof(ReductionMaximumLocation);
  const std::uint64_t messages =
      implementation.reduction_mode == ReductionMode::reproducible_tree
          ? static_cast<std::uint64_t>(implementation.size - 1) * 2U
          : 0U;
  if (!checked_add(candidate.calls, 1U, candidate.calls) ||
      !checked_add(candidate.scalars, count_u64, candidate.scalars) ||
      !checked_add(candidate.logical_bytes, bytes, candidate.logical_bytes) ||
      !checked_add(candidate.tree_messages, messages,
                   candidate.tree_messages) ||
      !checked_add(candidate.blocking_operations, 2U,
                   candidate.blocking_operations) ||
      !checked_add(candidate.wall_nanoseconds, elapsed,
                   candidate.wall_nanoseconds)) {
    return {StatusCode::invalid_plan, kCounterOverflow};
  }
  std::copy(implementation.maximum_location_receive_storage.data(),
            implementation.maximum_location_receive_storage.data() +
                local.size,
            global.data);
  implementation.reduction_counters = candidate;
  return {};
}

Status ReductionEngine::consensus(Status local_status) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kInvalidEngine};
  }
  if (!valid_status_code(local_status.code)) {
    local_status = {StatusCode::invalid_plan, kInvalidStatus};
  }
  initialize_packet(*implementation_, {}, local_status);
  return packet_collective(*implementation_, true);
}

Status ReductionEngine::consensus_contract(
    PlanFingerprint local_fingerprint) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kInvalidEngine};
  }
  const auto begin = std::chrono::steady_clock::now();
  std::uint64_t operations = 0U;
  const auto record = [&]() noexcept {
    const auto elapsed = static_cast<std::uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - begin)
            .count());
    LinearReductionCounters& counters =
        implementation_->reduction_counters;
    counters.blocking_operations =
        operations > UINT64_MAX - counters.blocking_operations
            ? UINT64_MAX
            : counters.blocking_operations + operations;
    counters.wall_nanoseconds =
        elapsed > UINT64_MAX - counters.wall_nanoseconds
            ? UINT64_MAX
            : counters.wall_nanoseconds + elapsed;
  };
  const std::uint64_t local = local_fingerprint;
  std::uint64_t reference = local;
  ++operations;
  if (MPI_Bcast(&reference, 1, MPI_UINT64_T, 0,
                implementation_->communicator) != MPI_SUCCESS) {
    implementation_->lowest_failing_rank = -1;
    record();
    return {StatusCode::mpi_failure, kCollectiveFailure};
  }
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  operations += 2U;
  const int minimum_status = MPI_Allreduce(
      &local, &minimum, 1, MPI_UINT64_T, MPI_MIN,
      implementation_->communicator);
  const int maximum_status = MPI_Allreduce(
      &local, &maximum, 1, MPI_UINT64_T, MPI_MAX,
      implementation_->communicator);
  if (minimum_status != MPI_SUCCESS || maximum_status != MPI_SUCCESS) {
    implementation_->lowest_failing_rank = -1;
    record();
    return {StatusCode::mpi_failure, kCollectiveFailure};
  }
  if (local_fingerprint != 0U && minimum == maximum) {
    implementation_->lowest_failing_rank = -1;
    record();
    return {};
  }
  const int candidate = local_fingerprint == reference
                            ? implementation_->size
                            : implementation_->rank;
  int first_mismatch = implementation_->size;
  ++operations;
  if (MPI_Allreduce(&candidate, &first_mismatch, 1, MPI_INT, MPI_MIN,
                    implementation_->communicator) != MPI_SUCCESS) {
    implementation_->lowest_failing_rank = -1;
    record();
    return {StatusCode::mpi_failure, kCollectiveFailure};
  }
  implementation_->lowest_failing_rank =
      first_mismatch < implementation_->size ? first_mismatch : 0;
  record();
  return {StatusCode::invalid_plan, kContractMismatch};
}

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
Status ReductionEngine::arm_checked_sum_fault_for_test(
    std::uint64_t ordinal, int rank) noexcept {
  if (implementation_ == nullptr || ordinal == 0U || rank < 0 ||
      rank >= implementation_->size) {
    return {StatusCode::invalid_plan, kCheckedSumFaultForTest};
  }
  g_checked_sum_fault_ordinal = ordinal;
  g_checked_sum_next_ordinal = 0U;
  g_checked_sum_fault_rank = rank;
  g_checked_sum_fault_armed = true;
  return {};
}

void ReductionEngine::clear_checked_sum_fault_for_test() noexcept {
  g_checked_sum_fault_ordinal = 0U;
  g_checked_sum_next_ordinal = 0U;
  g_checked_sum_fault_rank = -1;
  g_checked_sum_fault_armed = false;
}
#endif

Status ReductionEngine::validate_communicator(
    MPI_Comm communicator) const noexcept {
  if (!mpi_is_live() || implementation_ == nullptr ||
      communicator == MPI_COMM_NULL) {
    return {StatusCode::invalid_plan, kInvalidEngine};
  }
  int relation = MPI_UNEQUAL;
  if (MPI_Comm_compare(implementation_->communicator, communicator,
                       &relation) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kCompileCollective};
  }
  return relation == MPI_IDENT || relation == MPI_CONGRUENT
             ? Status{}
             : Status{StatusCode::invalid_plan, kContractMismatch};
}

int ReductionEngine::lowest_failing_rank() const noexcept {
  return implementation_ == nullptr ? -1
                                    : implementation_->lowest_failing_rank;
}

LinearReductionCounters ReductionEngine::counters() const noexcept {
  return implementation_ == nullptr ? LinearReductionCounters{}
                                    : implementation_->reduction_counters;
}

std::uintptr_t ReductionEngine::send_storage_address() const noexcept {
  return implementation_ == nullptr
             ? 0U
             : reinterpret_cast<std::uintptr_t>(
                   implementation_->send_storage.data());
}

std::uintptr_t ReductionEngine::receive_storage_address() const noexcept {
  return implementation_ == nullptr
             ? 0U
             : reinterpret_cast<std::uintptr_t>(
                   implementation_->receive_storage.data());
}

std::size_t ReductionEngine::capacity() const noexcept {
  return implementation_ == nullptr ? 0U : implementation_->maximum_scalars;
}

ReductionMode ReductionEngine::mode() const noexcept {
  return implementation_ == nullptr ? ReductionMode::mpi_allreduce
                                    : implementation_->reduction_mode;
}

}  // namespace hundun::v04

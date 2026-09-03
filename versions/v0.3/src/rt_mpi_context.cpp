// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/rt_mpi_context.hpp"

#include "hundun/rt_error.hpp"
#include "rt_mpi_error_detail.hpp"

#include <climits>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>

namespace hundun::runtime {
namespace {

#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
thread_local int next_synchronous_fp64_allreduce_pre_call_error_for_test =
    MPI_SUCCESS;

int consume_synchronous_next_fp64_allreduce_pre_call_error_for_test() noexcept {
  return std::exchange(
      next_synchronous_fp64_allreduce_pre_call_error_for_test, MPI_SUCCESS);
}
#endif

}  // namespace

namespace detail {

#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
void inject_synchronous_next_fp64_allreduce_pre_call_error_raw(
    int mpi_error) {
  if (mpi_error == MPI_SUCCESS) {
    throw Error("synchronous FP64 allreduce pre-call test result must be an "
                "MPI error");
  }
  if (next_synchronous_fp64_allreduce_pre_call_error_for_test != MPI_SUCCESS) {
    throw Error(
        "synchronous FP64 allreduce pre-call test error is already pending");
  }
  next_synchronous_fp64_allreduce_pre_call_error_for_test = mpi_error;
}
#endif

}  // namespace detail

MpiContext MpiContext::duplicate(MPI_Comm source) {
  detail::require_mpi_active("duplicate MPI communicator");
  if (source == MPI_COMM_NULL) {
    throw Error("MPI context requires a valid intracommunicator");
  }

  int is_intercommunicator = 0;
  detail::check_mpi(MPI_Comm_test_inter(source, &is_intercommunicator),
                    "MPI_Comm_test_inter");
  if (is_intercommunicator != 0) {
    throw Error("MPI context requires an intracommunicator");
  }

  MPI_Comm duplicate = MPI_COMM_NULL;
  detail::check_mpi(MPI_Comm_dup(source, &duplicate), "MPI_Comm_dup");
  try {
    detail::check_mpi(
        MPI_Comm_set_errhandler(duplicate, MPI_ERRORS_RETURN),
        "MPI_Comm_set_errhandler");
    int rank = 0;
    int size = 0;
    int thread_level = MPI_THREAD_SINGLE;
    detail::check_mpi(MPI_Comm_rank(duplicate, &rank), "MPI_Comm_rank");
    detail::check_mpi(MPI_Comm_size(duplicate, &size), "MPI_Comm_size");
    detail::check_mpi(MPI_Query_thread(&thread_level), "MPI_Query_thread");
    return MpiContext(duplicate, rank, size, thread_level);
  } catch (...) {
    detail::free_communicator_without_throwing(duplicate);
    throw;
  }
}

MpiContext::MpiContext(MPI_Comm communicator, int rank, int size,
                       int thread_level) noexcept
    : communicator_(communicator),
      rank_(rank),
      size_(size),
      thread_level_(thread_level) {}

MpiContext::~MpiContext() noexcept {
  detail::free_communicator_without_throwing(communicator_);
}

MpiContext::MpiContext(MpiContext&& other) noexcept
    : communicator_(std::exchange(other.communicator_, MPI_COMM_NULL)),
      rank_(other.rank_),
      size_(other.size_),
      thread_level_(other.thread_level_),
      fp64_reduction_counters_(
          std::exchange(other.fp64_reduction_counters_, {})) {}

MpiContext& MpiContext::operator=(MpiContext&& other) noexcept {
  if (this != &other) {
    detail::free_communicator_without_throwing(communicator_);
    communicator_ = std::exchange(other.communicator_, MPI_COMM_NULL);
    rank_ = other.rank_;
    size_ = other.size_;
    thread_level_ = other.thread_level_;
    fp64_reduction_counters_ =
        std::exchange(other.fp64_reduction_counters_, {});
  }
  return *this;
}

void MpiContext::barrier() const {
  detail::require_mpi_active("enter MPI barrier");
  if (communicator_ == MPI_COMM_NULL) {
    throw Error("cannot enter MPI barrier with an empty MPI context");
  }
  detail::check_mpi(MPI_Barrier(communicator_), "MPI_Barrier");
}

void MpiContext::allreduce_fp64_in_place(
    double* values, std::size_t count,
    Fp64ReductionOperation operation) const {
  detail::require_mpi_active("reduce FP64 values");
  if (communicator_ == MPI_COMM_NULL) {
    throw Error("cannot reduce FP64 values with an empty MPI context");
  }

  MPI_Op mpi_operation = MPI_OP_NULL;
  switch (operation) {
    case Fp64ReductionOperation::sum:
      mpi_operation = MPI_SUM;
      break;
    case Fp64ReductionOperation::maximum:
      mpi_operation = MPI_MAX;
      break;
    default:
      throw Error("unsupported FP64 reduction operation");
  }

  if (count == 0U) {
    return;
  }
  if (values == nullptr) {
    throw Error("FP64 reduction requires a non-null value pointer");
  }
  if (count > static_cast<std::size_t>(INT_MAX)) {
    throw Error("FP64 reduction count exceeds MPI INT_MAX");
  }

  const auto checked_add = [](std::uint64_t current, std::uint64_t increment) {
    if (increment > std::numeric_limits<std::uint64_t>::max() - current) {
      throw Error("FP64 reduction counter would overflow");
    }
    return current + increment;
  };
  const std::uint64_t scalar_count = static_cast<std::uint64_t>(count);
  if (scalar_count >
      std::numeric_limits<std::uint64_t>::max() / sizeof(double)) {
    throw Error("FP64 reduction logical payload byte count would overflow");
  }
  const std::uint64_t logical_bytes =
      scalar_count * static_cast<std::uint64_t>(sizeof(double));
  const Fp64ReductionCounters updated{
      checked_add(fp64_reduction_counters_.collective_calls, 1U),
      checked_add(fp64_reduction_counters_.reduced_scalars, scalar_count),
      checked_add(fp64_reduction_counters_.logical_payload_bytes,
                  logical_bytes)};

  fp64_reduction_counters_ = updated;
#ifdef HUNDUN_RUNTIME_ENABLE_TEST_ACCESS
  const int injected_result =
      consume_synchronous_next_fp64_allreduce_pre_call_error_for_test();
  const int mpi_result =
      injected_result == MPI_SUCCESS
          ? MPI_Allreduce(MPI_IN_PLACE, values, static_cast<int>(count),
                          MPI_DOUBLE, mpi_operation, communicator_)
          : injected_result;
#else
  const int mpi_result =
      MPI_Allreduce(MPI_IN_PLACE, values, static_cast<int>(count), MPI_DOUBLE,
                    mpi_operation, communicator_);
#endif
  detail::check_mpi(
      mpi_result,
      "MPI_Allreduce FP64 values");
}

}  // namespace hundun::runtime

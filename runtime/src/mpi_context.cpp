// SPDX-License-Identifier: Apache-2.0

#include "hundun/runtime/mpi_context.hpp"

#include "hundun/runtime/error.hpp"
#include "mpi_error.hpp"
#include "mpi_context_test_seam.hpp"

#include <climits>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>

namespace hundun::runtime {
namespace {

thread_local int next_fp64_allreduce_result_for_test = MPI_SUCCESS;

int consume_next_fp64_allreduce_result_for_test() noexcept {
  return std::exchange(next_fp64_allreduce_result_for_test, MPI_SUCCESS);
}

}  // namespace

namespace detail {

void inject_next_fp64_allreduce_result_for_test(int mpi_result) {
  if (mpi_result == MPI_SUCCESS) {
    throw Error("FP64 allreduce test result must be an MPI error");
  }
  if (next_fp64_allreduce_result_for_test != MPI_SUCCESS) {
    throw Error("FP64 allreduce test result is already pending");
  }
  next_fp64_allreduce_result_for_test = mpi_result;
}

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
  const int injected_result =
      consume_next_fp64_allreduce_result_for_test();
  const int mpi_result =
      injected_result == MPI_SUCCESS
          ? MPI_Allreduce(MPI_IN_PLACE, values, static_cast<int>(count),
                          MPI_DOUBLE, mpi_operation, communicator_)
          : injected_result;
  detail::check_mpi(
      mpi_result,
      "MPI_Allreduce FP64 values");
}

}  // namespace hundun::runtime

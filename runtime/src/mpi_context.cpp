// SPDX-License-Identifier: Apache-2.0

#include "hundun/runtime/mpi_context.hpp"

#include "hundun/runtime/error.hpp"
#include "mpi_error.hpp"

#include <utility>

namespace hundun::runtime {

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
      thread_level_(other.thread_level_) {}

MpiContext& MpiContext::operator=(MpiContext&& other) noexcept {
  if (this != &other) {
    detail::free_communicator_without_throwing(communicator_);
    communicator_ = std::exchange(other.communicator_, MPI_COMM_NULL);
    rank_ = other.rank_;
    size_ = other.size_;
    thread_level_ = other.thread_level_;
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

}  // namespace hundun::runtime

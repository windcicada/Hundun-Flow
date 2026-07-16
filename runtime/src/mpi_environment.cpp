// SPDX-License-Identifier: Apache-2.0

#include "hundun/runtime/mpi_environment.hpp"

#include "hundun/runtime/error.hpp"

#include <string>

namespace hundun::runtime {
namespace {

void check_mpi(int result, const char* operation) {
  if (result != MPI_SUCCESS) {
    throw Error(std::string(operation) + " failed with MPI error " +
                std::to_string(result));
  }
}

void finalize_without_throwing() noexcept {
  int finalized = 0;
  if (MPI_Finalized(&finalized) == MPI_SUCCESS && finalized == 0) {
    (void)MPI_Finalize();
  }
}

}  // namespace

MpiEnvironment::MpiEnvironment(int& argc, char**& argv) {
  int finalized = 0;
  check_mpi(MPI_Finalized(&finalized), "MPI_Finalized");
  if (finalized != 0) {
    throw Error("cannot construct MPI environment after MPI_Finalize");
  }

  int initialized = 0;
  check_mpi(MPI_Initialized(&initialized), "MPI_Initialized");

  int provided = MPI_THREAD_SINGLE;
  if (initialized == 0) {
    check_mpi(MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided),
              "MPI_Init_thread");
    owns_mpi_ = true;
  } else {
    check_mpi(MPI_Query_thread(&provided), "MPI_Query_thread");
  }

  try {
    if (provided < MPI_THREAD_FUNNELED) {
      throw Error("MPI runtime does not provide MPI_THREAD_FUNNELED");
    }
    check_mpi(MPI_Comm_rank(MPI_COMM_WORLD, &rank_), "MPI_Comm_rank");
    check_mpi(MPI_Comm_size(MPI_COMM_WORLD, &size_), "MPI_Comm_size");
  } catch (...) {
    if (owns_mpi_) {
      finalize_without_throwing();
      owns_mpi_ = false;
    }
    throw;
  }
}

MpiEnvironment::~MpiEnvironment() noexcept {
  if (owns_mpi_) {
    finalize_without_throwing();
  }
}

void MpiEnvironment::barrier() const {
  int finalized = 0;
  check_mpi(MPI_Finalized(&finalized), "MPI_Finalized");
  if (finalized != 0) {
    throw Error("cannot enter MPI barrier after MPI_Finalize");
  }
  check_mpi(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier");
}

}  // namespace hundun::runtime

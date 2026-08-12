// SPDX-License-Identifier: Apache-2.0

#include "hundun/rt_mpi_environment.hpp"

#include "hundun/rt_error.hpp"
#include "rt_mpi_error_detail.hpp"

#include <mpi.h>

namespace hundun::runtime {
namespace {

void finalize_without_throwing() noexcept {
  if (detail::mpi_is_active()) {
    (void)MPI_Finalize();
  }
}

}  // namespace

MpiEnvironment::MpiEnvironment(int& argc, char**& argv) {
  int finalized = 0;
  detail::check_mpi(MPI_Finalized(&finalized), "MPI_Finalized");
  if (finalized != 0) {
    throw Error("cannot construct MPI environment after MPI_Finalize");
  }

  int initialized = 0;
  detail::check_mpi(MPI_Initialized(&initialized), "MPI_Initialized");

  int provided = MPI_THREAD_SINGLE;
  if (initialized == 0) {
    detail::check_mpi(
        MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided),
        "MPI_Init_thread");
    owns_mpi_ = true;
  } else {
    detail::check_mpi(MPI_Query_thread(&provided), "MPI_Query_thread");
  }

  try {
    if (provided < MPI_THREAD_FUNNELED) {
      throw Error("MPI runtime does not provide MPI_THREAD_FUNNELED");
    }
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

}  // namespace hundun::runtime

// SPDX-License-Identifier: Apache-2.0

#include "mpi_error.hpp"

#include "hundun/runtime/error.hpp"

#include <array>
#include <string>

namespace hundun::runtime::detail {

bool mpi_is_active() noexcept {
  int initialized = 0;
  if (MPI_Initialized(&initialized) != MPI_SUCCESS || initialized == 0) {
    return false;
  }
  int finalized = 0;
  return MPI_Finalized(&finalized) == MPI_SUCCESS && finalized == 0;
}

std::string mpi_error_message(std::string_view operation, int result,
                              MpiErrorStringFunction error_string) {
  std::string message(operation);
  message += " failed with MPI error ";
  message += std::to_string(result);

  if (!mpi_is_active() || error_string == nullptr) {
    return message;
  }

  std::array<char, MPI_MAX_ERROR_STRING> buffer{};
  int length = 0;
  if (error_string(result, buffer.data(), &length) != MPI_SUCCESS ||
      length <= 0 || length > static_cast<int>(buffer.size())) {
    return message;
  }
  message += ": ";
  message.append(buffer.data(), static_cast<std::size_t>(length));
  return message;
}

void check_mpi(int result, std::string_view operation) {
  if (result != MPI_SUCCESS) {
    throw Error(mpi_error_message(operation, result));
  }
}

void require_mpi_active(std::string_view operation) {
  int initialized = 0;
  check_mpi(MPI_Initialized(&initialized), "MPI_Initialized");
  if (initialized == 0) {
    throw Error("cannot " + std::string(operation) + " before MPI_Init");
  }

  int finalized = 0;
  check_mpi(MPI_Finalized(&finalized), "MPI_Finalized");
  if (finalized != 0) {
    throw Error("cannot " + std::string(operation) + " after MPI_Finalize");
  }
}

void free_communicator_without_throwing(MPI_Comm& communicator) noexcept {
  if (communicator != MPI_COMM_NULL && mpi_is_active()) {
    (void)MPI_Comm_free(&communicator);
  }
  communicator = MPI_COMM_NULL;
}

}  // namespace hundun::runtime::detail

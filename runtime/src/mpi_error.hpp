// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/runtime/error.hpp"

#include <mpi.h>

#include <string>
#include <string_view>

namespace hundun::runtime::detail {

class MpiOperationError final : public runtime::Error {
 public:
  explicit MpiOperationError(const std::string& message)
      : runtime::Error(message) {}
};

using MpiErrorStringFunction = int (*)(int, char*, int*);

std::string mpi_error_message(
    std::string_view operation, int result,
    MpiErrorStringFunction error_string = MPI_Error_string);
void check_mpi(int result, std::string_view operation);
void require_mpi_active(std::string_view operation);
bool mpi_is_active() noexcept;
void free_communicator_without_throwing(MPI_Comm& communicator) noexcept;

}  // namespace hundun::runtime::detail

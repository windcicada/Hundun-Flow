// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/runtime/error.hpp"

#include <string>
#include <string_view>

namespace hundun::runtime {

class MpiOperationError final : public Error {
public:
  explicit MpiOperationError(const std::string &message);
};

void check_mpi_result(int result, std::string_view operation);

} // namespace hundun::runtime

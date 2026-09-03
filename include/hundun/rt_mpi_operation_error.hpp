// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/rt_error.hpp"

#include <string>
#include <string_view>

namespace hundun::runtime {

class MpiOperationError final : public Error {
public:
  explicit MpiOperationError(const std::string &message);
};

void check_mpi_result(int result, std::string_view operation);

} // namespace hundun::runtime

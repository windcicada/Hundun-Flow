// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_types.hpp"

#include <cstdint>
#include <string_view>

namespace hundun::v04 {

enum class StatusCode : std::uint16_t {
  ok,
  invalid_case,
  invalid_plan,
  allocation_failure,
  mpi_failure,
  numerical_failure,
  rejected_step,
  io_failure
};

struct Status {
  StatusCode code{StatusCode::ok};
  std::uint32_t detail{};
  constexpr explicit operator bool() const noexcept {
    return code == StatusCode::ok;
  }
};

std::string_view status_message(Status status) noexcept;

}  // namespace hundun::v04

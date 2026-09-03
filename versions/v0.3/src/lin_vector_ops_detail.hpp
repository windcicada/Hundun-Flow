// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include <string>
#include <utility>

#include "hundun/rt_error.hpp"

namespace hundun::linear::detail {

class SynchronizedReductionError final : public runtime::Error {
 public:
  SynchronizedReductionError(std::string message,
                             bool has_local_source)
      : runtime::Error(std::move(message)),
        has_local_source_(has_local_source) {}

  bool has_local_source() const noexcept { return has_local_source_; }

 private:
  bool has_local_source_;
};

double bad_dot_product_sentinel() noexcept;

}  // namespace hundun::linear::detail

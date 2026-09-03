// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include <stdexcept>
#include <string>
#include <utility>

namespace hundun::runtime {

class Error : public std::runtime_error {
 public:
  explicit Error(const std::string& message) : std::runtime_error(message) {}
};

class ConfigError final : public Error {
 public:
  ConfigError(std::string pointer, std::string message)
      : Error(pointer.empty() ? message : pointer + ": " + message),
        pointer_(std::move(pointer)) {}

  const std::string& pointer() const noexcept { return pointer_; }

 private:
  std::string pointer_;
};

}  // namespace hundun::runtime

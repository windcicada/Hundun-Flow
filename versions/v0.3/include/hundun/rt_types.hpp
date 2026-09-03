// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include <cstdint>

namespace hundun::runtime {

struct Int3 {
  int x{};
  int y{};
  int z{};
};

struct Real3 {
  double x{};
  double y{};
  double z{};
};

struct Box3 {
  Int3 begin{};
  Int3 end{};
};

inline std::int64_t volume(Int3 extent) {
  return static_cast<std::int64_t>(extent.x) * extent.y * extent.z;
}

}  // namespace hundun::runtime

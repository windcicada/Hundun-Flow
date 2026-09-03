// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include <cstddef>
#include <cstdint>

namespace hundun::v04 {

struct Int3 {
  std::int32_t x{}, y{}, z{};
};

struct Real3 {
  double x{}, y{}, z{};
};

template <class T>
struct Span {
  T* data{};
  std::size_t size{};
};

using FieldId = std::uint16_t;
using StageId = std::uint16_t;
using RevisionToken = std::uint64_t;
using PlanFingerprint = std::uint64_t;

}  // namespace hundun::v04

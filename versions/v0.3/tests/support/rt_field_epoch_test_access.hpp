// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/rt_field_storage.hpp"

#include <cstdint>
#include <memory>

namespace hundun::runtime {

namespace detail {

class FieldEpochControl;

std::uint64_t field_epoch_generation_raw(
    const std::shared_ptr<FieldEpochControl> &epoch) noexcept;
void field_epoch_force_generation_raw(
    const std::shared_ptr<FieldEpochControl> &epoch,
    std::uint64_t generation) noexcept;

struct FieldEpochTestAccess final {
  static std::uint64_t generation(const FieldStorage &storage) noexcept {
    return field_epoch_generation_raw(storage.epoch_);
  }
  static void force_generation(FieldStorage &storage,
                               std::uint64_t generation) noexcept {
    field_epoch_force_generation_raw(storage.epoch_, generation);
  }
};

}  // namespace detail
}  // namespace hundun::runtime

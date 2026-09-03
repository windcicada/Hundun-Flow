// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include <cstddef>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "hundun/rt_field_descriptor.hpp"

namespace hundun::runtime {

class FieldRegistry final {
 public:
  FieldId declare_field(FieldDescriptor descriptor);
  void freeze() noexcept;
  bool frozen() const noexcept;
  std::size_t size() const noexcept;
  const FieldDescriptor &descriptor(FieldId id) const;
  FieldId field_id(std::string_view name) const;

 private:
  std::vector<FieldDescriptor> descriptors_;
  std::unordered_map<std::string, FieldId> field_ids_;
  bool frozen_{};
};

}  // namespace hundun::runtime

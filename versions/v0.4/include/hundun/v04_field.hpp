// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_status.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace hundun::v04 {

struct FieldDescriptor {
  FieldId id{};
  std::string stable_name;
  std::uint8_t components{};
  std::uint8_t ghost_width{};
};

class FieldSchema {
 public:
  using const_iterator = std::vector<FieldDescriptor>::const_iterator;

  std::size_t size() const noexcept { return fields_.size(); }
  bool empty() const noexcept { return fields_.empty(); }
  const FieldDescriptor& operator[](std::size_t index) const noexcept {
    return fields_[index];
  }
  const_iterator begin() const noexcept { return fields_.begin(); }
  const_iterator end() const noexcept { return fields_.end(); }

 private:
  friend class FieldRegistry;
  std::vector<FieldDescriptor> fields_;
};

class FieldRegistry {
 public:
  FieldRegistry() = default;
  FieldRegistry(const FieldRegistry&) = default;
  FieldRegistry& operator=(const FieldRegistry&) = default;
  FieldRegistry(FieldRegistry&&) noexcept = default;
  FieldRegistry& operator=(FieldRegistry&&) noexcept = default;

  Status declare_field(std::string_view stable_name, std::uint8_t components,
                       std::uint8_t ghost_width, FieldId& out);
  Status require_field(std::string_view stable_name, std::uint8_t components,
                       std::uint8_t minimum_ghost_width, FieldId& out);
  Status freeze_for_test(FieldSchema& out);
  PlanFingerprint fingerprint() const noexcept;

 private:
  std::vector<FieldDescriptor> fields_;
  std::unordered_set<std::string> names_;
  bool frozen_{};
};

}  // namespace hundun::v04

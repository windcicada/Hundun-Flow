// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_field.hpp"

#include <limits>

namespace hundun::v04 {

namespace {

constexpr std::uint8_t kMaxFieldComponents = 64U;
constexpr std::uint8_t kMaxGhostWidth = 16U;

}  // namespace

Status FieldRegistry::declare_field(std::string_view stable_name,
                                    std::uint8_t components,
                                    std::uint8_t ghost_width, FieldId& out) {
  if (frozen_ || stable_name.empty() ||
      components == 0 || components > kMaxFieldComponents ||
      ghost_width > kMaxGhostWidth ||
      fields_.size() > std::numeric_limits<FieldId>::max()) {
    return {StatusCode::invalid_plan, 1};
  }
  try {
    if (names_.find(std::string(stable_name)) != names_.end()) {
      return {StatusCode::invalid_plan, 2};
    }
    const FieldId id = static_cast<FieldId>(fields_.size());
    std::string owned_name(stable_name);
    const auto inserted = names_.insert(owned_name);
    if (!inserted.second) {
      return {StatusCode::invalid_plan, 2};
    }
    try {
      fields_.push_back(
          FieldDescriptor{id, std::move(owned_name), components, ghost_width});
    } catch (...) {
      names_.erase(inserted.first);
      throw;
    }
    out = id;
  } catch (...) {
    return {StatusCode::allocation_failure, 0};
  }
  return {};
}

Status FieldRegistry::freeze_for_test(FieldSchema& out) {
  if (frozen_) {
    return {StatusCode::invalid_plan, 3};
  }
  try {
    FieldSchema snapshot;
    snapshot.fields_ = fields_;
    out = std::move(snapshot);
  } catch (...) {
    return {StatusCode::allocation_failure, 0};
  }
  frozen_ = true;
  return {};
}

}  // namespace hundun::v04

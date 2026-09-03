// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_field.hpp"

#include <cstdint>
#include <limits>

namespace hundun::v04 {

namespace {

constexpr std::uint8_t kMaxGhostWidth = 16U;

std::uint64_t hash_mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= UINT64_C(1099511628211);
  return hash;
}

}  // namespace

Status FieldRegistry::declare_field(std::string_view stable_name,
                                    std::uint8_t components,
                                    std::uint8_t ghost_width, FieldId& out) {
  if (frozen_ || stable_name.empty() || components == 0U ||
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

Status FieldRegistry::require_field(std::string_view stable_name,
                                    std::uint8_t components,
                                    std::uint8_t minimum_ghost_width,
                                    FieldId& out) {
  if (frozen_ || stable_name.empty() || components == 0U ||
      minimum_ghost_width > kMaxGhostWidth) {
    return {StatusCode::invalid_plan, 1};
  }
  for (FieldDescriptor& field : fields_) {
    if (field.stable_name != stable_name) {
      continue;
    }
    if (field.components != components) {
      return {StatusCode::invalid_plan, 4};
    }
    if (field.ghost_width < minimum_ghost_width) {
      field.ghost_width = minimum_ghost_width;
    }
    out = field.id;
    return {};
  }
  return declare_field(stable_name, components, minimum_ghost_width, out);
}

Status FieldRegistry::freeze(FieldSchema& out) {
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

PlanFingerprint FieldRegistry::fingerprint() const noexcept {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  hash = hash_mix(hash, fields_.size());
  for (const FieldDescriptor& field : fields_) {
    hash = hash_mix(hash, field.id);
    hash = hash_mix(hash, field.components);
    hash = hash_mix(hash, field.ghost_width);
    hash = hash_mix(hash, field.stable_name.size());
    for (const unsigned char character : field.stable_name) {
      hash = hash_mix(hash, character);
    }
  }
  hash = hash_mix(hash, frozen_ ? 1U : 0U);
  return hash == 0U ? 1U : hash;
}

}  // namespace hundun::v04

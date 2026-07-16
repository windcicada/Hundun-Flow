// SPDX-License-Identifier: Apache-2.0

#include "hundun/runtime/field_registry.hpp"

#include <limits>
#include <string>
#include <utility>

#include "hundun/runtime/error.hpp"

namespace hundun::runtime {
namespace {

bool recognized(FunctionSpace space) {
  switch (space) {
    case FunctionSpace::cell_average:
    case FunctionSpace::face_value:
    case FunctionSpace::vertex_value:
    case FunctionSpace::element_dof:
    case FunctionSpace::quadrature_point:
    case FunctionSpace::particle:
      return true;
  }
  return false;
}

bool recognized(ScalarType scalar_type) {
  switch (scalar_type) {
    case ScalarType::float64:
    case ScalarType::int32:
    case ScalarType::uint8:
      return true;
  }
  return false;
}

bool recognized(RestartPolicy policy) {
  switch (policy) {
    case RestartPolicy::persistent:
    case RestartPolicy::transient:
      return true;
  }
  return false;
}

bool recognized(OutputPolicy policy) {
  switch (policy) {
    case OutputPolicy::never:
    case OutputPolicy::selected:
    case OutputPolicy::always:
      return true;
  }
  return false;
}

}  // namespace

FieldId FieldRegistry::declare_field(FieldDescriptor descriptor_value) {
  if (frozen_) {
    throw Error("cannot declare a field after the registry is frozen");
  }
  if (descriptor_value.name.empty()) {
    throw Error("field name must not be empty");
  }
  if (descriptor_value.unit.empty()) {
    throw Error("field unit must not be empty");
  }
  if (descriptor_value.owner.empty()) {
    throw Error("field owner must not be empty");
  }
  if (descriptor_value.components == 0U) {
    throw Error("field component count must be positive");
  }
  if (descriptor_value.ghost_width < 0) {
    throw Error("field ghost width must not be negative");
  }
  if (!recognized(descriptor_value.space) ||
      !recognized(descriptor_value.scalar_type) ||
      !recognized(descriptor_value.restart) ||
      !recognized(descriptor_value.output)) {
    throw Error("field descriptor contains an unrecognized enum value");
  }
  if (field_ids_.find(descriptor_value.name) != field_ids_.end()) {
    throw Error("field name is already declared: " + descriptor_value.name);
  }
  if (descriptors_.size() >
      static_cast<std::size_t>(std::numeric_limits<FieldId>::max())) {
    throw Error("field registry has exhausted sequential field IDs");
  }

  const auto id = static_cast<FieldId>(descriptors_.size());
  descriptors_.push_back(std::move(descriptor_value));
  try {
    const auto inserted = field_ids_.emplace(descriptors_.back().name, id);
    if (!inserted.second) {
      descriptors_.pop_back();
      throw Error("field name is already declared");
    }
  } catch (...) {
    if (descriptors_.size() == static_cast<std::size_t>(id) + 1U) {
      descriptors_.pop_back();
    }
    throw;
  }
  return id;
}

void FieldRegistry::freeze() noexcept { frozen_ = true; }

bool FieldRegistry::frozen() const noexcept { return frozen_; }

std::size_t FieldRegistry::size() const noexcept { return descriptors_.size(); }

const FieldDescriptor &FieldRegistry::descriptor(FieldId id) const {
  const auto index = static_cast<std::size_t>(id);
  if (index >= descriptors_.size()) {
    throw Error("field descriptor ID is out of bounds");
  }
  return descriptors_[index];
}

FieldId FieldRegistry::field_id(std::string_view name) const {
  const auto found = field_ids_.find(std::string(name));
  if (found == field_ids_.end()) {
    throw Error("field name is not registered: " + std::string(name));
  }
  return found->second;
}

}  // namespace hundun::runtime

// SPDX-License-Identifier: Apache-2.0

#include "hundun/rt_field_access_plan.hpp"

#include <algorithm>
#include <cstddef>

#include "hundun/rt_error.hpp"
#include "hundun/rt_field_registry.hpp"

namespace hundun::runtime {
namespace {

bool recognized_mode(AccessMode mode) noexcept {
  switch (mode) {
    case AccessMode::read:
    case AccessMode::write:
    case AccessMode::read_write:
      return true;
  }
  return false;
}

bool writes(AccessMode mode) noexcept {
  return mode == AccessMode::write || mode == AccessMode::read_write;
}

}  // namespace

FieldAccessPlan::FieldAccessPlan(const FieldRegistry &registry)
    : field_count_(registry.size()) {
  if (!registry.frozen()) {
    throw Error("field access plan requires a frozen registry");
  }
}

void FieldAccessPlan::declare_access(PhaseId phase, ActorId actor,
                                     FieldId field, AccessMode mode) {
  if (frozen_) {
    throw Error("field access plan is frozen");
  }
  if (static_cast<std::size_t>(field) >= field_count_) {
    throw Error("field access declaration references an unknown field");
  }
  if (!recognized_mode(mode)) {
    throw Error("field access declaration mode is unrecognized");
  }

  const auto same_triple = [phase, actor, field](const Declaration &entry) {
    return entry.phase == phase && entry.actor == actor &&
           entry.field == field;
  };
  if (std::find_if(declarations_.begin(), declarations_.end(), same_triple) !=
      declarations_.end()) {
    throw Error("field access declaration is duplicated");
  }

  if (writes(mode)) {
    const auto existing_writer =
        [phase, field](const Declaration &entry) {
          return entry.phase == phase && entry.field == field &&
                 writes(entry.mode);
        };
    if (std::find_if(declarations_.begin(), declarations_.end(),
                     existing_writer) != declarations_.end()) {
      throw Error(
          "field access plan already has a writer for this phase and field");
    }
  }

  declarations_.push_back(Declaration{phase, actor, field, mode});
}

void FieldAccessPlan::freeze() noexcept { frozen_ = true; }

bool FieldAccessPlan::frozen() const noexcept { return frozen_; }

const FieldAccessPlan::Declaration &FieldAccessPlan::require_declaration(
    PhaseId phase, ActorId actor, FieldId field,
    std::size_t storage_field_count) const {
  if (!frozen_) {
    throw Error("field access acquisition requires a frozen plan");
  }
  if (field_count_ != storage_field_count) {
    throw Error("field access plan field domain does not match storage");
  }
  const auto same_triple = [phase, actor, field](const Declaration &entry) {
    return entry.phase == phase && entry.actor == actor &&
           entry.field == field;
  };
  const auto declaration =
      std::find_if(declarations_.begin(), declarations_.end(), same_triple);
  if (declaration == declarations_.end()) {
    throw Error("field access was not declared");
  }
  return *declaration;
}

void FieldAccessPlan::require_read(PhaseId phase, ActorId actor, FieldId field,
                                   std::size_t storage_field_count) const {
  const auto &declaration =
      require_declaration(phase, actor, field, storage_field_count);
  if (declaration.mode == AccessMode::write) {
    throw Error("field access declaration does not permit reading");
  }
}

void FieldAccessPlan::require_write(PhaseId phase, ActorId actor,
                                    FieldId field,
                                    std::size_t storage_field_count) const {
  const auto &declaration =
      require_declaration(phase, actor, field, storage_field_count);
  if (declaration.mode == AccessMode::read) {
    throw Error("field access declaration does not permit writing");
  }
}

}  // namespace hundun::runtime

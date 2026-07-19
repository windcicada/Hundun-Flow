// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "hundun/runtime/field_descriptor.hpp"

namespace hundun::runtime {

class FieldRegistry;
class FieldStorage;

using ActorId = std::uint32_t;
using PhaseId = std::uint32_t;

enum class AccessMode { read, write, read_write };

class FieldAccessPlan final {
 public:
  explicit FieldAccessPlan(const FieldRegistry &registry);

  void declare_access(PhaseId phase, ActorId actor, FieldId field,
                      AccessMode mode);
  void freeze() noexcept;
  bool frozen() const noexcept;

 private:
  friend class FieldStorage;

  struct Declaration {
    PhaseId phase{};
    ActorId actor{};
    FieldId field{};
    AccessMode mode{};
  };

  void require_read(PhaseId phase, ActorId actor, FieldId field,
                    std::size_t storage_field_count) const;
  void require_write(PhaseId phase, ActorId actor, FieldId field,
                     std::size_t storage_field_count) const;
  const Declaration &require_declaration(PhaseId phase, ActorId actor,
                                         FieldId field,
                                         std::size_t storage_field_count) const;

  std::size_t field_count_{};
  std::vector<Declaration> declarations_;
  bool frozen_{};
};

}  // namespace hundun::runtime

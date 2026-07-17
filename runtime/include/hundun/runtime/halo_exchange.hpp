// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/runtime/exchange_plan.hpp"
#include "hundun/runtime/field_descriptor.hpp"

#include <memory>

namespace hundun::runtime {

class FieldStorage;
class StructuredDecomposition;

class HaloExchange final {
 public:
  static HaloExchange create(const StructuredDecomposition& decomposition,
                             ExchangePlan plan);
  ~HaloExchange() noexcept;

  HaloExchange(HaloExchange&&) noexcept;
  HaloExchange& operator=(HaloExchange&&) = delete;
  HaloExchange(const HaloExchange&) = delete;
  HaloExchange& operator=(const HaloExchange&) = delete;

  void exchange(FieldStorage& storage, FieldId id);
  void begin(const FieldStorage& storage, FieldId id);
  void wait(FieldStorage& storage, FieldId id);

 private:
  class Impl;

  explicit HaloExchange(std::unique_ptr<Impl> implementation) noexcept;

  std::unique_ptr<Impl> implementation_;
};

}  // namespace hundun::runtime

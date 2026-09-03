// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/rt_exchange_plan.hpp"
#include "hundun/rt_field_descriptor.hpp"
#include "hundun/rt_halo_performance_counters.hpp"

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

  // Returns the immutable held plan width without calling MPI. The value is
  // safe to query after MPI finalization; a moved-from object throws Error.
  [[nodiscard]] int ghost_width() const;
  // While MPI is active, reports whether the held communicator and immutable
  // plan match the supplied live decomposition. A moved-from Halo or
  // decomposition throws Error; rank reordering or a plan mismatch is false.
  [[nodiscard]] bool is_compatible_with(
      const StructuredDecomposition& decomposition) const;
  [[nodiscard]] HaloPerformanceCounters performance_counters() const;

  void exchange(FieldStorage& storage, FieldId id);
  void begin(const FieldStorage& storage, FieldId id);
  void wait(FieldStorage& storage, FieldId id);

 private:
  class Impl;

  explicit HaloExchange(std::unique_ptr<Impl> implementation) noexcept;

  std::unique_ptr<Impl> implementation_;
};

}  // namespace hundun::runtime

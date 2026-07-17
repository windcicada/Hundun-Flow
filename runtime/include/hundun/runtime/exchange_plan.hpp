// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/runtime/structured_decomposition.hpp"
#include "hundun/runtime/types.hpp"

#include <vector>

namespace hundun::runtime {

struct ExchangeRegion {
  Int3 offset;
  int neighbor_rank;
  Box3 send_box;
  Box3 receive_box;
};

class ExchangePlan {
 public:
  static ExchangePlan create(const StructuredDecomposition& decomposition,
                             Int3 local_extent, int ghost_width);
  const std::vector<ExchangeRegion>& regions() const noexcept;
  int ghost_width() const noexcept;

 private:
  friend class HaloExchange;

  ExchangePlan(Int3 local_extent, int ghost_width,
               std::vector<ExchangeRegion> regions) noexcept;

  Int3 local_extent_{};
  int ghost_width_{};
  std::vector<ExchangeRegion> regions_;
};

}  // namespace hundun::runtime

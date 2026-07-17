// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/runtime/types.hpp"

#include <array>
#include <cstdint>

namespace hundun::runtime::detail {

struct ExactProcessGridCost {
  std::uint64_t carry;
  std::uint64_t low;
};

ExactProcessGridCost process_grid_cost(
    Int3 cells, std::array<bool, 3> periodic, Int3 grid);
Int3 select_process_grid(int rank_count, Int3 cells,
                         std::array<bool, 3> periodic);
void validate_explicit_process_grid(Int3 grid, int rank_count, Int3 cells);

}  // namespace hundun::runtime::detail

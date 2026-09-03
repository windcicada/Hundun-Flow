// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_mesh.hpp"

#include "mesh_focus_detail.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <tuple>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kDecompositionInput = 201U;
constexpr std::uint32_t kDecompositionRanks = 202U;
constexpr std::uint32_t kTileInput = 203U;

std::int32_t component(Int3 value, int axis) noexcept {
  return axis == 0 ? value.x : (axis == 1 ? value.y : value.z);
}

void set_component(Int3& value, int axis, std::int32_t component_value) noexcept {
  if (axis == 0) {
    value.x = component_value;
  } else if (axis == 1) {
    value.y = component_value;
  } else {
    value.z = component_value;
  }
}

bool valid_cells(Int3 cells) noexcept {
  return cells.x > 0 && cells.y > 0 && cells.z > 0;
}

struct GridChoice {
  Int3 grid{};
  struct WideCost {
    std::uint64_t high{};
    std::uint64_t low{};
  } surface{std::numeric_limits<std::uint64_t>::max(),
            std::numeric_limits<std::uint64_t>::max()};
  std::uint64_t maximum_patch{std::numeric_limits<std::uint64_t>::max()};
  bool valid{};
};

bool lexicographically_less(Int3 left, Int3 right) noexcept {
  return std::tie(left.x, left.y, left.z) <
         std::tie(right.x, right.y, right.z);
}

bool less(GridChoice::WideCost left, GridChoice::WideCost right) noexcept {
  return left.high < right.high ||
         (left.high == right.high && left.low < right.low);
}

GridChoice::WideCost wide_multiply(std::uint64_t left,
                                   std::uint64_t right) noexcept {
  constexpr std::uint64_t mask = 0xffffffffULL;
  const std::uint64_t low_product = (left & mask) * right;
  const std::uint64_t high_product = (left >> 32U) * right;
  const std::uint64_t shifted_high = high_product << 32U;
  const std::uint64_t low = low_product + shifted_high;
  const std::uint64_t carry = low < low_product ? 1U : 0U;
  return {(high_product >> 32U) + carry, low};
}

GridChoice::WideCost add(GridChoice::WideCost left,
                         GridChoice::WideCost term) noexcept {
  const std::uint64_t previous = left.low;
  left.low += term.low;
  const std::uint64_t carry = left.low < previous ? 1U : 0U;
  const std::uint64_t high_increment = term.high + carry;
  left.high = high_increment < term.high ||
                      high_increment >
                          std::numeric_limits<std::uint64_t>::max() - left.high
                  ? std::numeric_limits<std::uint64_t>::max()
                  : left.high + high_increment;
  return left;
}

GridChoice::WideCost interface_surface(Int3 cells, Int3 grid) noexcept {
  const auto x = static_cast<std::uint64_t>(cells.x);
  const auto y = static_cast<std::uint64_t>(cells.y);
  const auto z = static_cast<std::uint64_t>(cells.z);
  const auto px = static_cast<std::uint64_t>(grid.x);
  const auto py = static_cast<std::uint64_t>(grid.y);
  const auto pz = static_cast<std::uint64_t>(grid.z);
  GridChoice::WideCost result{};
  result = add(result, wide_multiply(y * z, px - 1U));
  result = add(result, wide_multiply(x * z, py - 1U));
  return add(result, wide_multiply(x * y, pz - 1U));
}

void consider_grid(Int3 cells, Int3 candidate, GridChoice& best) noexcept {
  if (candidate.x > cells.x || candidate.y > cells.y ||
      candidate.z > cells.z) {
    return;
  }
  const GridChoice::WideCost surface = interface_surface(cells, candidate);
  const std::uint64_t maximum =
      detail::maximum_patch_cells(cells, candidate);
  if (maximum == 0U) {
    return;
  }
  if (!best.valid || less(surface, best.surface) ||
      (!less(best.surface, surface) && maximum < best.maximum_patch) ||
      (!less(surface, best.surface) && !less(best.surface, surface) &&
       maximum == best.maximum_patch &&
       lexicographically_less(candidate, best.grid))) {
    best = {candidate, surface, maximum, true};
  }
}

GridChoice choose_grid(Int3 cells, std::int32_t ranks) noexcept {
  GridChoice best;
  // Enumerate each unordered factor triple once, then evaluate its six axis
  // permutations.  This replaces an O(ranks) scan with divisor searches bound
  // by cbrt(ranks) and sqrt(remaining), including large prime rank counts.
  for (std::int32_t a = 1; a <= ranks / a / a; ++a) {
    if (ranks % a != 0) {
      continue;
    }
    const std::int32_t remaining = ranks / a;
    for (std::int32_t b = a; b <= remaining / b; ++b) {
      if (remaining % b != 0) {
        continue;
      }
      const std::int32_t c = remaining / b;
      if (b > c) {
        continue;
      }
      const std::array<Int3, 6U> permutations{
          Int3{a, b, c}, Int3{a, c, b}, Int3{b, a, c},
          Int3{b, c, a}, Int3{c, a, b}, Int3{c, b, a}};
      for (const Int3 candidate : permutations) {
        consider_grid(cells, candidate, best);
      }
    }
  }
  return best;
}

void split_axis(std::int32_t global, std::int32_t partitions,
                std::int32_t coordinate, std::int32_t& begin,
                std::int32_t& count) noexcept {
  const std::int32_t quotient = global / partitions;
  const std::int32_t remainder = global % partitions;
  count = quotient + (coordinate < remainder ? 1 : 0);
  const std::int64_t wide_begin =
      static_cast<std::int64_t>(coordinate) * quotient +
      std::min(coordinate, remainder);
  begin = static_cast<std::int32_t>(wide_begin);
}

}  // namespace

namespace detail {

std::uint64_t maximum_patch_cells(Int3 global_cells,
                                  Int3 process_grid) noexcept {
  if (!valid_cells(global_cells) || !valid_cells(process_grid)) {
    return 0U;
  }
  const std::uint64_t x = static_cast<std::uint64_t>(
      global_cells.x / process_grid.x +
      (global_cells.x % process_grid.x != 0 ? 1 : 0));
  const std::uint64_t y = static_cast<std::uint64_t>(
      global_cells.y / process_grid.y +
      (global_cells.y % process_grid.y != 0 ? 1 : 0));
  const std::uint64_t z = static_cast<std::uint64_t>(
      global_cells.z / process_grid.z +
      (global_cells.z % process_grid.z != 0 ? 1 : 0));
  if (x > std::numeric_limits<std::uint64_t>::max() / y ||
      x * y > std::numeric_limits<std::uint64_t>::max() / z) {
    return 0U;
  }
  return x * y * z;
}

Status make_mesh_patch(std::int32_t rank, std::int32_t rank_count,
                       Int3 global_cells, MeshPatch& out) noexcept {
  if (rank_count <= 0 || rank < 0 || rank >= rank_count ||
      !valid_cells(global_cells)) {
    return {StatusCode::invalid_plan, kDecompositionInput};
  }
  const GridChoice choice = choose_grid(global_cells, rank_count);
  if (!choice.valid) {
    return {StatusCode::invalid_plan, kDecompositionRanks};
  }
  MeshPatch candidate;
  candidate.process_grid = choice.grid;
  candidate.process_coord.x = rank % choice.grid.x;
  const std::int32_t yz_rank = rank / choice.grid.x;
  candidate.process_coord.y = yz_rank % choice.grid.y;
  candidate.process_coord.z = yz_rank / choice.grid.y;
  for (int axis = 0; axis < 3; ++axis) {
    std::int32_t begin = 0;
    std::int32_t count = 0;
    split_axis(component(global_cells, axis), component(choice.grid, axis),
               component(candidate.process_coord, axis), begin, count);
    set_component(candidate.begin, axis, begin);
    set_component(candidate.cells, axis, count);
  }
  out = candidate;
  return {};
}

}  // namespace detail

Status build_cpu_tiles(const MeshPatch& patch, Int3 target_cells,
                       std::vector<CpuTile>& out) noexcept {
  if (!valid_cells(patch.cells) || !valid_cells(target_cells) ||
      patch.begin.x < 0 || patch.begin.y < 0 || patch.begin.z < 0) {
    return {StatusCode::invalid_plan, kTileInput};
  }
  try {
    const auto patch_end_valid = [](std::int32_t begin,
                                    std::int32_t cells) noexcept {
      return static_cast<std::int64_t>(begin) + cells - 1 <=
             std::numeric_limits<std::int32_t>::max();
    };
    if (!patch_end_valid(patch.begin.x, patch.cells.x) ||
        !patch_end_valid(patch.begin.y, patch.cells.y) ||
        !patch_end_valid(patch.begin.z, patch.cells.z)) {
      return {StatusCode::invalid_plan, kTileInput};
    }
    const auto ceil_div = [](std::int32_t extent,
                             std::int32_t tile) noexcept -> std::uint64_t {
      const auto wide_extent = static_cast<std::uint64_t>(extent);
      const auto wide_tile = static_cast<std::uint64_t>(tile);
      return wide_extent / wide_tile +
             (wide_extent % wide_tile != 0U ? 1U : 0U);
    };
    const std::uint64_t nx = ceil_div(patch.cells.x, target_cells.x);
    const std::uint64_t ny = ceil_div(patch.cells.y, target_cells.y);
    const std::uint64_t nz = ceil_div(patch.cells.z, target_cells.z);
    if (nx > std::numeric_limits<std::size_t>::max() / ny ||
        nx * ny > std::numeric_limits<std::size_t>::max() / nz) {
      return {StatusCode::invalid_plan, kTileInput};
    }
    std::vector<CpuTile> candidate;
    candidate.reserve(static_cast<std::size_t>(nx * ny * nz));
    // z/y/x ordering makes successive tiles retain the solver's x-unit-stride
    // loop direction while leaving ownership independent from thread policy.
    for (std::uint64_t iz = 0U; iz < nz; ++iz) {
      const std::int64_t z = static_cast<std::int64_t>(iz) * target_cells.z;
      for (std::uint64_t iy = 0U; iy < ny; ++iy) {
        const std::int64_t y = static_cast<std::int64_t>(iy) * target_cells.y;
        for (std::uint64_t ix = 0U; ix < nx; ++ix) {
          const std::int64_t x =
              static_cast<std::int64_t>(ix) * target_cells.x;
          const std::int64_t global_x = patch.begin.x + x;
          const std::int64_t global_y = patch.begin.y + y;
          const std::int64_t global_z = patch.begin.z + z;
          if (global_x > std::numeric_limits<std::int32_t>::max() ||
              global_y > std::numeric_limits<std::int32_t>::max() ||
              global_z > std::numeric_limits<std::int32_t>::max()) {
            return {StatusCode::invalid_plan, kTileInput};
          }
          candidate.push_back(CpuTile{
              Int3{static_cast<std::int32_t>(global_x),
                   static_cast<std::int32_t>(global_y),
                   static_cast<std::int32_t>(global_z)},
              Int3{static_cast<std::int32_t>(std::min<std::int64_t>(
                       target_cells.x, patch.cells.x - x)),
                   static_cast<std::int32_t>(std::min<std::int64_t>(
                       target_cells.y, patch.cells.y - y)),
                   static_cast<std::int32_t>(std::min<std::int64_t>(
                       target_cells.z, patch.cells.z - z))}});
        }
      }
    }
    out = std::move(candidate);
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, 0U};
  } catch (...) {
    return {StatusCode::invalid_plan, kTileInput};
  }
}

}  // namespace hundun::v04

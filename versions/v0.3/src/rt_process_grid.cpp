// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "rt_process_grid_detail.hpp"

#include "hundun/rt_error.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <vector>

namespace hundun::runtime::detail {
namespace {

std::array<int, 3> components(Int3 value) {
  return {value.x, value.y, value.z};
}

void validate_cell_domain(Int3 cells) {
  const auto values = components(cells);
  if (std::any_of(values.begin(), values.end(),
                  [](int value) { return value <= 0; })) {
    throw Error("process-grid cell extents must be positive");
  }

  const std::uint64_t limit =
      static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max());
  std::uint64_t volume = static_cast<std::uint64_t>(cells.x);
  for (const int value : {cells.y, cells.z}) {
    const std::uint64_t factor = static_cast<std::uint64_t>(value);
    if (volume > limit / factor) {
      throw Error("process-grid cell count exceeds INT64_MAX");
    }
    volume *= factor;
  }
}

void validate_grid_geometry(Int3 grid, Int3 cells) {
  const auto grid_values = components(grid);
  const auto cell_values = components(cells);
  for (std::size_t axis = 0; axis < grid_values.size(); ++axis) {
    if (grid_values[axis] <= 0) {
      throw Error("process-grid dimensions must be positive");
    }
    if (grid_values[axis] > cell_values[axis]) {
      throw Error("process-grid dimension exceeds cell extent");
    }
  }
}

std::uint64_t checked_multiply(std::uint64_t lhs, std::uint64_t rhs,
                               const char* message) {
  const std::uint64_t limit = std::numeric_limits<std::uint64_t>::max();
  if (rhs != 0U && lhs > limit / rhs) {
    throw Error(message);
  }
  return lhs * rhs;
}

std::uint64_t cost_term(int partitions, bool periodic, int face_a,
                        int face_b) {
  const int wrap = periodic && partitions > 1 ? 1 : 0;
  const std::uint64_t coefficient =
      static_cast<std::uint64_t>(partitions - 1 + wrap);
  const std::uint64_t face = checked_multiply(
      static_cast<std::uint64_t>(face_a),
      static_cast<std::uint64_t>(face_b),
      "process-grid cost term overflow");
  return checked_multiply(coefficient, face,
                          "process-grid cost term overflow");
}

void add_limb(ExactProcessGridCost& cost, std::uint64_t term) {
  const std::uint64_t previous = cost.low;
  cost.low += term;
  if (cost.low < previous) {
    if (cost.carry == std::numeric_limits<std::uint64_t>::max()) {
      throw Error("process-grid exact cost carry overflow");
    }
    ++cost.carry;
  }
}

ExactProcessGridCost compute_cost(Int3 cells,
                                  std::array<bool, 3> periodic,
                                  Int3 grid) {
  ExactProcessGridCost cost{0U, 0U};
  add_limb(cost,
           cost_term(grid.x, periodic[0], cells.y, cells.z));
  add_limb(cost,
           cost_term(grid.y, periodic[1], cells.x, cells.z));
  add_limb(cost,
           cost_term(grid.z, periodic[2], cells.x, cells.y));
  return cost;
}

bool less_cost(ExactProcessGridCost lhs, ExactProcessGridCost rhs) {
  return lhs.carry < rhs.carry ||
         (lhs.carry == rhs.carry && lhs.low < rhs.low);
}

bool equal_cost(ExactProcessGridCost lhs, ExactProcessGridCost rhs) {
  return lhs.carry == rhs.carry && lhs.low == rhs.low;
}

bool lexicographically_less(Int3 lhs, Int3 rhs) {
  if (lhs.x != rhs.x) {
    return lhs.x < rhs.x;
  }
  if (lhs.y != rhs.y) {
    return lhs.y < rhs.y;
  }
  return lhs.z < rhs.z;
}

std::vector<int> divisors(int number) {
  std::vector<int> result;
  for (int candidate = 1; candidate <= number / candidate; ++candidate) {
    if (number % candidate != 0) {
      continue;
    }
    result.push_back(candidate);
    const int paired = number / candidate;
    if (paired != candidate) {
      result.push_back(paired);
    }
  }
  std::sort(result.begin(), result.end());
  return result;
}

}  // namespace

ExactProcessGridCost process_grid_cost(Int3 cells,
                                       std::array<bool, 3> periodic,
                                       Int3 grid) {
  validate_cell_domain(cells);
  validate_grid_geometry(grid, cells);
  return compute_cost(cells, periodic, grid);
}

Int3 select_process_grid(int rank_count, Int3 cells,
                         std::array<bool, 3> periodic) {
  if (rank_count <= 0) {
    throw Error("process-grid rank count must be positive");
  }
  validate_cell_domain(cells);

  bool found = false;
  Int3 selected{};
  ExactProcessGridCost selected_cost{};
  for (const int px : divisors(rank_count)) {
    const int remaining = rank_count / px;
    for (const int py : divisors(remaining)) {
      const int pz = remaining / py;
      const Int3 candidate{px, py, pz};
      if (candidate.x > cells.x || candidate.y > cells.y ||
          candidate.z > cells.z) {
        continue;
      }
      const ExactProcessGridCost candidate_cost =
          compute_cost(cells, periodic, candidate);
      if (!found || less_cost(candidate_cost, selected_cost) ||
          (equal_cost(candidate_cost, selected_cost) &&
           lexicographically_less(candidate, selected))) {
        found = true;
        selected = candidate;
        selected_cost = candidate_cost;
      }
    }
  }

  if (!found) {
    throw Error("no feasible process grid for rank count and cell shape");
  }
  return selected;
}

void validate_explicit_process_grid(Int3 grid, int rank_count, Int3 cells) {
  if (rank_count <= 0) {
    throw Error("process-grid rank count must be positive");
  }
  validate_cell_domain(cells);

  const auto grid_values = components(grid);
  if (std::any_of(grid_values.begin(), grid_values.end(),
                  [](int value) { return value <= 0; })) {
    throw Error("process-grid dimensions must be positive");
  }

  std::uint64_t product = 1U;
  for (const int value : grid_values) {
    product = checked_multiply(product, static_cast<std::uint64_t>(value),
                               "process-grid product overflow");
  }
  if (product >
      static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    throw Error("process-grid product exceeds MPI int rank domain");
  }
  if (static_cast<int>(product) != rank_count) {
    throw Error("process-grid product does not equal communicator size");
  }

  validate_grid_geometry(grid, cells);
}

}  // namespace hundun::runtime::detail

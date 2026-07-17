// SPDX-License-Identifier: Apache-2.0

#include "runtime/src/process_grid.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/types.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <climits>
#include <cstdint>
#include <string>

namespace {

using hundun::runtime::Int3;
using hundun::runtime::detail::ExactProcessGridCost;
using hundun::runtime::detail::process_grid_cost;
using hundun::runtime::detail::select_process_grid;
using hundun::runtime::detail::validate_explicit_process_grid;

bool same(Int3 lhs, Int3 rhs) {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

template <class Function>
void expect_runtime_error(Function&& function) {
  bool rejected = false;
  try {
    function();
  } catch (const hundun::runtime::Error& error) {
    HUNDUN_CHECK(!std::string(error.what()).empty());
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

void test_frozen_automatic_selections() {
  HUNDUN_CHECK(same(select_process_grid(
                        4, Int3{1, 2, 2}, {false, false, false}),
                    Int3{1, 2, 2}));
  HUNDUN_CHECK(same(select_process_grid(
                        12, Int3{1, 3, 4}, {false, false, false}),
                    Int3{1, 3, 4}));
  HUNDUN_CHECK(same(select_process_grid(
                        4, Int3{4, 4, 4}, {false, false, false}),
                    Int3{1, 2, 2}));
  HUNDUN_CHECK(same(select_process_grid(
                        2, Int3{8, 4, 2}, {true, false, false}),
                    Int3{1, 2, 1}));

  constexpr Int3 invariant_cells{17, 11, 7};
  constexpr std::array<bool, 3> invariant_periodic{true, false, true};
  HUNDUN_CHECK(same(select_process_grid(1, invariant_cells,
                                        invariant_periodic),
                    Int3{1, 1, 1}));
  HUNDUN_CHECK(same(select_process_grid(2, invariant_cells,
                                        invariant_periodic),
                    Int3{1, 2, 1}));
  HUNDUN_CHECK(same(select_process_grid(4, invariant_cells,
                                        invariant_periodic),
                    Int3{2, 2, 1}));
  HUNDUN_CHECK(same(select_process_grid(12, invariant_cells,
                                        invariant_periodic),
                    Int3{4, 3, 1}));
}

void test_frozen_exact_costs() {
  const ExactProcessGridCost ordinary = process_grid_cost(
      Int3{7, 11, 13}, {true, false, true}, Int3{2, 3, 1});
  HUNDUN_CHECK(ordinary.carry == 0U);
  HUNDUN_CHECK(ordinary.low == 468U);

  const ExactProcessGridCost carried = process_grid_cost(
      Int3{INT_MAX, INT_MAX, 2}, {true, true, true},
      Int3{INT_MAX, INT_MAX, 2});
  HUNDUN_CHECK(carried.carry == 1U);
  HUNDUN_CHECK(carried.low == UINT64_C(9223372011084972038));
}

void test_explicit_validation() {
  validate_explicit_process_grid(Int3{2, 2, 1}, 4, Int3{8, 4, 2});

  for (const Int3 grid : {Int3{0, 2, 1}, Int3{2, 0, 1}, Int3{2, 2, 0},
                          Int3{-1, 2, 1}, Int3{2, -1, 1},
                          Int3{2, 2, -1}}) {
    expect_runtime_error([&] {
      validate_explicit_process_grid(grid, 4, Int3{8, 4, 2});
    });
  }

  expect_runtime_error([] {
    validate_explicit_process_grid(
        Int3{INT_MAX, INT_MAX, INT_MAX}, 1,
        Int3{1, 1, 1});
  });
  expect_runtime_error([] {
    validate_explicit_process_grid(Int3{INT_MAX, 2, 1}, 1,
                                   Int3{INT_MAX, 2, 1});
  });
  expect_runtime_error([] {
    validate_explicit_process_grid(Int3{2, 2, 1}, 8, Int3{8, 4, 2});
  });
  expect_runtime_error([] {
    validate_explicit_process_grid(Int3{3, 1, 1}, 3, Int3{2, 2, 2});
  });
}

void test_automatic_rejections() {
  expect_runtime_error([] {
    static_cast<void>(select_process_grid(
        0, Int3{8, 4, 2}, {false, false, false}));
  });
  expect_runtime_error([] {
    static_cast<void>(select_process_grid(
        -1, Int3{8, 4, 2}, {false, false, false}));
  });
  expect_runtime_error([] {
    static_cast<void>(select_process_grid(
        6, Int3{2, 2, 2}, {false, false, false}));
  });
}

void test_cell_domain_rejections() {
  for (const Int3 cells : {Int3{0, 2, 2}, Int3{2, 0, 2}, Int3{2, 2, 0},
                           Int3{-1, 2, 2}, Int3{2, -1, 2},
                           Int3{2, 2, -1}}) {
    expect_runtime_error([&] {
      static_cast<void>(select_process_grid(
          1, cells, {false, false, false}));
    });
    expect_runtime_error([&] {
      static_cast<void>(process_grid_cost(
          cells, {false, false, false}, Int3{1, 1, 1}));
    });
    expect_runtime_error([&] {
      validate_explicit_process_grid(Int3{1, 1, 1}, 1, cells);
    });
  }

  expect_runtime_error([] {
    static_cast<void>(select_process_grid(
        1, Int3{INT_MAX, INT_MAX, INT_MAX}, {false, false, false}));
  });
  expect_runtime_error([] {
    static_cast<void>(process_grid_cost(
        Int3{INT_MAX, INT_MAX, INT_MAX}, {false, false, false},
        Int3{1, 1, 1}));
  });
  expect_runtime_error([] {
    validate_explicit_process_grid(
        Int3{1, 1, 1}, 1, Int3{INT_MAX, INT_MAX, INT_MAX});
  });
}

void run_all_tests() {
  test_frozen_automatic_selections();
  test_frozen_exact_costs();
  test_explicit_validation();
  test_automatic_rejections();
  test_cell_domain_rejections();
}

}  // namespace

int main() { return hundun::test::run(run_all_tests); }

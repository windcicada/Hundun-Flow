// SPDX-License-Identifier: Apache-2.0

#include "hundun/mesh/uniform_structured_mesh.hpp"

#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "hundun/runtime/types.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>

namespace {

using hundun::mesh::UniformStructuredMesh;
using hundun::runtime::Box3;
using hundun::runtime::DecompositionOptions;
using hundun::runtime::Error;
using hundun::runtime::Int3;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::Real3;
using hundun::runtime::StructuredDecomposition;

constexpr Int3 kGlobalCells{4, 2, 2};
constexpr Real3 kOrigin{-1.0, 0.0, 2.0};
constexpr Real3 kLength{2.0, 1.0, 4.0};
constexpr std::array<bool, 3> kNonperiodic{false, false, false};
constexpr double kTolerance = 1.0e-12;

static_assert(noexcept(std::declval<const UniformStructuredMesh&>().spacing_m()));
static_assert(
    noexcept(std::declval<const UniformStructuredMesh&>().cell_volume_m3()));
static_assert(
    noexcept(std::declval<const UniformStructuredMesh&>().local_extent()));
static_assert(noexcept(
    std::declval<const UniformStructuredMesh&>().owned_global_box()));

bool same(Int3 lhs, Int3 rhs) {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

bool same(Box3 lhs, Box3 rhs) {
  return same(lhs.begin, rhs.begin) && same(lhs.end, rhs.end);
}

void check_near(Real3 actual, Real3 expected) {
  HUNDUN_CHECK_NEAR(actual.x, expected.x, kTolerance);
  HUNDUN_CHECK_NEAR(actual.y, expected.y, kTolerance);
  HUNDUN_CHECK_NEAR(actual.z, expected.z, kTolerance);
}

Int3 with_component(Int3 value, int axis, int component) {
  if (axis == 0) {
    value.x = component;
  } else if (axis == 1) {
    value.y = component;
  } else {
    value.z = component;
  }
  return value;
}

Real3 with_component(Real3 value, int axis, double component) {
  if (axis == 0) {
    value.x = component;
  } else if (axis == 1) {
    value.y = component;
  } else {
    value.z = component;
  }
  return value;
}

template <class Function>
void expect_error(Function&& function) {
  bool threw = false;
  try {
    function();
  } catch (const Error& error) {
    threw = true;
    HUNDUN_CHECK(!std::string(error.what()).empty());
  }
  HUNDUN_CHECK(threw);
}

StructuredDecomposition make_decomposition(const MpiContext& mpi) {
  return StructuredDecomposition::create(
      mpi, kGlobalCells, kNonperiodic,
      DecompositionOptions{Int3{2, 1, 1}});
}

void test_geometry_and_ownership(const MpiContext& mpi) {
  auto decomposition = make_decomposition(mpi);
  const Int3 expected_extent = decomposition.local_extent();
  const Box3 expected_box = decomposition.owned_box();
  const UniformStructuredMesh mesh(kGlobalCells, kOrigin, kLength,
                                   decomposition);

  check_near(mesh.spacing_m(), Real3{0.5, 0.5, 2.0});
  HUNDUN_CHECK_NEAR(mesh.cell_volume_m3(), 0.5, kTolerance);
  HUNDUN_CHECK(same(mesh.local_extent(), expected_extent));
  HUNDUN_CHECK(same(mesh.owned_global_box(), expected_box));

  int local_first = 0;
  int local_last = 0;
  for (int k = 0; k < expected_extent.z; ++k) {
    for (int j = 0; j < expected_extent.y; ++j) {
      for (int i = 0; i < expected_extent.x; ++i) {
        const Int3 local{i, j, k};
        const Int3 global{expected_box.begin.x + i,
                          expected_box.begin.y + j,
                          expected_box.begin.z + k};
        const Real3 expected_center{
            kOrigin.x + (static_cast<double>(global.x) + 0.5) * 0.5,
            kOrigin.y + (static_cast<double>(global.y) + 0.5) * 0.5,
            kOrigin.z + (static_cast<double>(global.z) + 0.5) * 2.0};
        const Real3 actual = mesh.cell_center(local);
        check_near(actual, expected_center);

        if (same(global, Int3{0, 0, 0})) {
          check_near(actual, Real3{-0.75, 0.25, 3.0});
          local_first = 1;
        }
        if (same(global, Int3{3, 1, 1})) {
          check_near(actual, Real3{0.75, 0.75, 5.0});
          local_last = 1;
        }
      }
    }
  }

  const int local_nonzero_begin =
      expected_box.begin.x != 0 || expected_box.begin.y != 0 ||
              expected_box.begin.z != 0
          ? 1
          : 0;
  int global_first = 0;
  int global_last = 0;
  int ranks_with_nonzero_begin = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_first, &global_first, 1, MPI_INT, MPI_SUM,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(&local_last, &global_last, 1, MPI_INT, MPI_SUM,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(&local_nonzero_begin, &ranks_with_nonzero_begin,
                             1, MPI_INT, MPI_SUM,
                             mpi.comm()) == MPI_SUCCESS);
  HUNDUN_CHECK(global_first == 1);
  HUNDUN_CHECK(global_last == 1);
  HUNDUN_CHECK(ranks_with_nonzero_begin >= 1);

  for (int axis = 0; axis < 3; ++axis) {
    expect_error([&] {
      static_cast<void>(
          mesh.cell_center(with_component(Int3{0, 0, 0}, axis, -1)));
    });
    expect_error([&] {
      static_cast<void>(mesh.cell_center(
          with_component(Int3{0, 0, 0}, axis,
                         axis == 0 ? expected_extent.x
                                   : (axis == 1 ? expected_extent.y
                                                : expected_extent.z))));
    });
  }
}

void test_constructor_validation(const MpiContext& mpi) {
  auto decomposition = make_decomposition(mpi);

  for (int axis = 0; axis < 3; ++axis) {
    for (int invalid : {0, -1}) {
      expect_error([&] {
        static_cast<void>(UniformStructuredMesh(
            with_component(kGlobalCells, axis, invalid), kOrigin, kLength,
            decomposition));
      });
    }
  }

  expect_error([&] {
    static_cast<void>(UniformStructuredMesh(Int3{5, 2, 2}, kOrigin, kLength,
                                            decomposition));
  });

  for (int axis = 0; axis < 3; ++axis) {
    for (double invalid : {0.0, -1.0}) {
      expect_error([&] {
        static_cast<void>(UniformStructuredMesh(
            kGlobalCells, kOrigin, with_component(kLength, axis, invalid),
            decomposition));
      });
    }
  }

  const double infinity = std::numeric_limits<double>::infinity();
  const double quiet_nan = std::numeric_limits<double>::quiet_NaN();
  for (int axis = 0; axis < 3; ++axis) {
    for (double invalid : {infinity, quiet_nan}) {
      expect_error([&] {
        static_cast<void>(UniformStructuredMesh(
            kGlobalCells, with_component(kOrigin, axis, invalid), kLength,
            decomposition));
      });
      expect_error([&] {
        static_cast<void>(UniformStructuredMesh(
            kGlobalCells, kOrigin, with_component(kLength, axis, invalid),
            decomposition));
      });
    }
  }

  expect_error([&] {
    static_cast<void>(UniformStructuredMesh(
        kGlobalCells, kOrigin,
        Real3{std::numeric_limits<double>::denorm_min(), 1.0, 1.0},
        decomposition));
  });
  expect_error([&] {
    static_cast<void>(UniformStructuredMesh(
        kGlobalCells, kOrigin, Real3{1.0e-200, 1.0e-200, 1.0e-200},
        decomposition));
  });
  expect_error([&] {
    static_cast<void>(UniformStructuredMesh(
        kGlobalCells, kOrigin, Real3{1.0e200, 1.0e200, 1.0e200},
        decomposition));
  });
  expect_error([&] {
    const double maximum = std::numeric_limits<double>::max();
    static_cast<void>(UniformStructuredMesh(
        kGlobalCells, Real3{maximum, 0.0, 0.0},
        Real3{maximum, 1.0e-100, 1.0e-100}, decomposition));
  });
}

void test_value_state_lifetime(const MpiContext& mpi) {
  auto decomposition =
      std::make_unique<StructuredDecomposition>(make_decomposition(mpi));
  const Int3 expected_extent = decomposition->local_extent();
  const Box3 expected_box = decomposition->owned_box();
  const UniformStructuredMesh mesh(kGlobalCells, kOrigin, kLength,
                                   *decomposition);
  decomposition.reset();

  check_near(mesh.spacing_m(), Real3{0.5, 0.5, 2.0});
  HUNDUN_CHECK_NEAR(mesh.cell_volume_m3(), 0.5, kTolerance);
  HUNDUN_CHECK(same(mesh.local_extent(), expected_extent));
  HUNDUN_CHECK(same(mesh.owned_global_box(), expected_box));
  check_near(mesh.cell_center(Int3{0, 0, 0}),
             Real3{kOrigin.x +
                       (static_cast<double>(expected_box.begin.x) + 0.5) * 0.5,
                   kOrigin.y +
                       (static_cast<double>(expected_box.begin.y) + 0.5) * 0.5,
                   kOrigin.z +
                       (static_cast<double>(expected_box.begin.z) + 0.5) * 2.0});

  mpi.barrier();
}

void run_tests(const MpiContext& mpi) {
  HUNDUN_CHECK(mpi.size() == 2);
  test_geometry_and_ownership(mpi);
  test_constructor_validation(mpi);
  test_value_state_lifetime(mpi);
}

}  // namespace

int main(int argc, char** argv) {
  int result = EXIT_FAILURE;
  {
    MpiEnvironment environment(argc, argv);
    auto mpi = MpiContext::duplicate(MPI_COMM_WORLD);
    result = hundun::test::run([&] { run_tests(mpi); });
  }
  return result;
}

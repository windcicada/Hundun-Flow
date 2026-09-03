// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/ib_surface_query.hpp"

#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/ib_test_access.hpp"
#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

using hundun::immersed::CellRegion;
using hundun::immersed::ImmersedSurface;
using hundun::immersed::SurfaceQuery;
using hundun::runtime::MpiContext;
using hundun::runtime::Real3;
using hundun::test::Stage3TemporaryDirectory;

bool near(double actual, double expected, double scale = 1.0) {
  return std::abs(actual - expected) <=
         256.0 * std::numeric_limits<double>::epsilon() *
             std::max({1.0, std::abs(expected), scale});
}

void check_point(Real3 actual, Real3 expected) {
  HUNDUN_CHECK(near(actual.x, expected.x));
  HUNDUN_CHECK(near(actual.y, expected.y));
  HUNDUN_CHECK(near(actual.z, expected.z));
}

template <class Operation>
void expect_error(Operation operation, const std::string &marker) {
  bool threw = false;
  try {
    operation();
  } catch (const hundun::runtime::Error &error) {
    threw = true;
    HUNDUN_CHECK(std::string{error.what()}.find(marker) != std::string::npos);
  }
  HUNDUN_CHECK(threw);
}

SurfaceQuery make_owned_query(const std::filesystem::path &path,
                              const MpiContext &mpi) {
  const auto surface = ImmersedSurface::load_collective(path, 1.0, mpi, 0);
  return SurfaceQuery::create(surface);
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  auto mpi = MpiContext::duplicate(MPI_COMM_WORLD);
  HUNDUN_CHECK(mpi.size() == 1);
  Stage3TemporaryDirectory directory{"query-unit"};
  const auto path = directory.path() / "cube.stl";
  hundun::test::write_text(
      path, hundun::test::ascii_stl(hundun::test::outward_cube(), "cube"));

  const auto surface = ImmersedSurface::load_collective(path, 1.0, mpi, 0);
  const auto query = SurfaceQuery::create(surface);
  const std::uint64_t original_fingerprint = query.fingerprint();

  const auto closest = query.closest_point({-1.0, 0.25, 0.25});
  HUNDUN_CHECK(closest.triangle == 8U);
  check_point(closest.point_m, {0.0, 0.25, 0.25});
  check_point(closest.geometric_outward_normal, {-1.0, 0.0, 0.0});
  HUNDUN_CHECK(near(closest.squared_distance_m2, 1.0));

  const auto tied = query.closest_point({0.5, 0.5, 0.5});
  HUNDUN_CHECK(tied.triangle == 0U);
  HUNDUN_CHECK(near(tied.squared_distance_m2, 0.25));

  const auto edge_closest = query.closest_point({-1.0, -1.0, 0.5});
  HUNDUN_CHECK(edge_closest.triangle == 5U);
  check_point(edge_closest.point_m, {0.0, 0.0, 0.5});
  HUNDUN_CHECK(near(edge_closest.squared_distance_m2, 2.0));

  const auto vertex_closest = query.closest_point({-1.0, -1.0, -1.0});
  HUNDUN_CHECK(vertex_closest.triangle == 0U);
  check_point(vertex_closest.point_m, {0.0, 0.0, 0.0});
  HUNDUN_CHECK(near(vertex_closest.squared_distance_m2, 3.0));

  const auto intersections =
      query.segment_intersections({-1.0, 0.25, 0.25}, {2.0, 0.25, 0.25});
  HUNDUN_CHECK(intersections.size() == 2U);
  HUNDUN_CHECK(intersections[0].triangle == 8U);
  HUNDUN_CHECK(intersections[1].triangle == 10U);
  HUNDUN_CHECK(near(intersections[0].segment_fraction, 1.0 / 3.0));
  HUNDUN_CHECK(near(intersections[1].segment_fraction, 2.0 / 3.0));
  check_point(intersections[0].point_m, {0.0, 0.25, 0.25});
  check_point(intersections[1].point_m, {1.0, 0.25, 0.25});

  const auto reversed =
      query.segment_intersections({2.0, 0.25, 0.25}, {-1.0, 0.25, 0.25});
  HUNDUN_CHECK(reversed.size() == 2U);
  HUNDUN_CHECK(reversed[0].triangle == 10U);
  HUNDUN_CHECK(reversed[1].triangle == 8U);
  HUNDUN_CHECK(near(reversed[0].segment_fraction, 1.0 / 3.0));
  HUNDUN_CHECK(near(reversed[1].segment_fraction, 2.0 / 3.0));

  const auto edge_hits =
      query.segment_intersections({-1.0, -1.0, 0.5}, {2.0, 2.0, 0.5});
  HUNDUN_CHECK(edge_hits.size() == 2U);
  HUNDUN_CHECK(edge_hits[0].triangle == 5U);
  HUNDUN_CHECK(edge_hits[1].triangle == 7U);

  const auto vertex_hits =
      query.segment_intersections({-1.0, -1.0, -1.0}, {2.0, 2.0, 2.0});
  HUNDUN_CHECK(vertex_hits.size() == 2U);
  HUNDUN_CHECK(vertex_hits[0].triangle == 0U);
  HUNDUN_CHECK(vertex_hits[1].triangle == 2U);
  HUNDUN_CHECK(vertex_hits[0].segment_fraction <
               vertex_hits[1].segment_fraction);

  HUNDUN_CHECK(query.classify({0.5, 0.5, 0.5},
                              hundun::config::ImmersedFluidSide::outside) ==
               CellRegion::solid);
  HUNDUN_CHECK(query.classify({2.0, 2.0, 2.0},
                              hundun::config::ImmersedFluidSide::outside) ==
               CellRegion::fluid);
  HUNDUN_CHECK(query.classify({0.5, 0.5, 0.5},
                              hundun::config::ImmersedFluidSide::inside) ==
               CellRegion::fluid);
  HUNDUN_CHECK(query.classify({2.0, 2.0, 2.0},
                              hundun::config::ImmersedFluidSide::inside) ==
               CellRegion::solid);
  HUNDUN_CHECK(query.classify({-0.1, 0.5, 0.5},
                              hundun::config::ImmersedFluidSide::outside) ==
               CellRegion::fluid);
  HUNDUN_CHECK(query.classify({0.5, 1.1, 0.5},
                              hundun::config::ImmersedFluidSide::inside) ==
               CellRegion::solid);
  expect_error(
      [&] {
        const double near_boundary =
            -256.0 * std::numeric_limits<double>::epsilon();
        (void)query.classify({near_boundary, 0.5, 0.5},
                             hundun::config::ImmersedFluidSide::outside);
      },
      "boundary");

  const auto candidates =
      hundun::immersed::test::ImmersedTestAccess::bounded_candidates(
          query, {0.2, 0.2, -0.01}, {0.8, 0.8, 0.01});
  HUNDUN_CHECK(candidates ==
               std::vector<hundun::immersed::TriangleId>({0U, 1U}));
  HUNDUN_CHECK(candidates ==
               hundun::immersed::test::ImmersedTestAccess::bounded_candidates(
                   query, {0.2, 0.2, -0.01}, {0.8, 0.8, 0.01}));
  const auto all_candidates =
      hundun::immersed::test::ImmersedTestAccess::bounded_candidates(
          query, {-0.01, -0.01, -0.01}, {1.01, 1.01, 1.01});
  HUNDUN_CHECK(all_candidates.size() == surface.triangle_count());
  for (std::size_t index = 0U; index < all_candidates.size(); ++index) {
    HUNDUN_CHECK(all_candidates[index] == index);
  }
  const auto order =
      hundun::immersed::test::ImmersedTestAccess::bvh_order(query);
  HUNDUN_CHECK(order.size() == surface.triangle_count());
  auto sorted_order = order;
  std::sort(sorted_order.begin(), sorted_order.end());
  for (std::size_t index = 0; index < sorted_order.size(); ++index) {
    HUNDUN_CHECK(sorted_order[index] == index);
  }
  HUNDUN_CHECK(order == hundun::immersed::test::ImmersedTestAccess::bvh_order(
                            SurfaceQuery::create(surface)));
  std::vector<hundun::immersed::TriangleId> reversed_initial_order(
      surface.triangle_count());
  for (std::size_t index = 0U; index < reversed_initial_order.size(); ++index) {
    reversed_initial_order[index] = reversed_initial_order.size() - 1U - index;
  }
  const auto traversal_permuted_query =
      hundun::immersed::test::ImmersedTestAccess::
          create_query_with_initial_order(surface, reversed_initial_order);
  HUNDUN_CHECK(traversal_permuted_query.fingerprint() == query.fingerprint());
  HUNDUN_CHECK(hundun::immersed::test::ImmersedTestAccess::bvh_order(
                   traversal_permuted_query) == order);
  HUNDUN_CHECK(
      traversal_permuted_query.closest_point({-1.0, 0.25, 0.25}).triangle ==
      closest.triangle);
  HUNDUN_CHECK(hundun::immersed::test::ImmersedTestAccess::resolve_parity(
      {true, true, true}));
  HUNDUN_CHECK(!hundun::immersed::test::ImmersedTestAccess::resolve_parity(
      {false, false, false}));

  expect_error(
      [&] {
        (void)query.segment_intersections({0.2, 0.2, 0.0}, {0.8, 0.2, 0.0});
      },
      "coplanar");
  expect_error(
      [&] {
        (void)query.segment_intersections({0.5, 0.5, 0.5}, {0.5, 0.5, 0.5});
      },
      "length");
  expect_error(
      [&] {
        (void)query.closest_point(
            {std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0});
      },
      "non-finite");
  expect_error(
      [&] {
        (void)query.segment_intersections(
            {0.0, 0.0, 0.0},
            {std::numeric_limits<double>::infinity(), 0.0, 0.0});
      },
      "non-finite");
  expect_error(
      [&] {
        (void)query.classify({0.5, 0.5, 0.0},
                             hundun::config::ImmersedFluidSide::outside);
      },
      "boundary");
  expect_error(
      [&] {
        (void)hundun::immersed::test::ImmersedTestAccess::bounded_candidates(
            query, {1.0, 0.0, 0.0}, {0.0, 1.0, 1.0});
      },
      "bounds");
  expect_error(
      [&] {
        (void)hundun::immersed::test::ImmersedTestAccess::resolve_parity(
            {true, false, true});
      },
      "parity disagrees");
  expect_error(
      [&] {
        auto invalid_order = reversed_initial_order;
        invalid_order.front() = invalid_order.back();
        (void)hundun::immersed::test::ImmersedTestAccess::
            create_query_with_initial_order(surface, std::move(invalid_order));
      },
      "permutation");
  expect_error(
      [&] {
        (void)query.classify(
            {0.5, 0.5, 0.5},
            static_cast<hundun::config::ImmersedFluidSide>(255U));
      },
      "fluid side");

  HUNDUN_CHECK(query.fingerprint() == original_fingerprint);
  const auto owned_query = make_owned_query(path, mpi);
  HUNDUN_CHECK(owned_query.fingerprint() == query.fingerprint());
  HUNDUN_CHECK(
      owned_query.classify({0.5, 0.5, 0.5},
                           hundun::config::ImmersedFluidSide::outside) ==
      CellRegion::solid);

  auto permuted = hundun::test::outward_cube();
  std::rotate(permuted.begin(), permuted.begin() + 3, permuted.end());
  const auto permuted_path = directory.path() / "cube-permuted.stl";
  hundun::test::write_text(permuted_path,
                           hundun::test::ascii_stl(permuted, "permuted"));
  const auto permuted_surface =
      ImmersedSurface::load_collective(permuted_path, 1.0, mpi, 0);
  const auto permuted_query = SurfaceQuery::create(permuted_surface);
  HUNDUN_CHECK(permuted_query.fingerprint() != query.fingerprint());
  HUNDUN_CHECK(
      permuted_query.classify({0.5, 0.5, 0.5},
                              hundun::config::ImmersedFluidSide::outside) ==
      CellRegion::solid);

  const auto tiny_surface =
      ImmersedSurface::load_collective(path, 1.0e-6, mpi, 0);
  const auto tiny_query = SurfaceQuery::create(tiny_surface);
  const auto tiny_intersections = tiny_query.segment_intersections(
      {-1.0e-6, 0.25e-6, 0.25e-6}, {2.0e-6, 0.25e-6, 0.25e-6});
  HUNDUN_CHECK(tiny_intersections.size() == 2U);
  HUNDUN_CHECK(
      tiny_query.classify({0.5e-6, 0.5e-6, 0.5e-6},
                          hundun::config::ImmersedFluidSide::outside) ==
      CellRegion::solid);

  return 0;
}

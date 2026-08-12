// SPDX-License-Identifier: Apache-2.0

#include "hundun/ib_surface.hpp"

#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/ib_test_access.hpp"
#include "tests/support/stage3_stl_fixture.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

using hundun::immersed::ImmersedSurface;
using hundun::immersed::SurfaceTriangle;
using hundun::immersed::test::ImmersedTestAccess;
using hundun::runtime::MpiContext;
using hundun::runtime::Real3;
using hundun::test::Stage3TemporaryDirectory;
using hundun::test::StlFixtureTriangle;

std::uint64_t double_bits(double value) {
  std::uint64_t bits = 0U;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

bool near(double actual, double expected, double scale = 1.0) {
  return std::abs(actual - expected) <=
         256.0 * std::numeric_limits<double>::epsilon() *
             std::max({1.0, std::abs(expected), scale});
}

Real3 subtract(Real3 a, Real3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }

double dot(Real3 a, Real3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

struct TriangleEvidence final {
  std::uint64_t id{};
  std::array<Real3, 3> vertices{};
  Real3 normal{};
  double area{};
};

struct SurfaceEvidence final {
  std::size_t vertex_count{};
  Real3 minimum{};
  Real3 maximum{};
  double volume{};
  std::uint64_t fingerprint{};
  std::vector<TriangleEvidence> triangles;
};

bool same_real3(Real3 first, Real3 second) {
  return double_bits(first.x) == double_bits(second.x) &&
         double_bits(first.y) == double_bits(second.y) &&
         double_bits(first.z) == double_bits(second.z);
}

bool same_surface_evidence(const SurfaceEvidence &first,
                           const SurfaceEvidence &second) {
  if (first.vertex_count != second.vertex_count ||
      !same_real3(first.minimum, second.minimum) ||
      !same_real3(first.maximum, second.maximum) ||
      double_bits(first.volume) != double_bits(second.volume) ||
      first.fingerprint != second.fingerprint ||
      first.triangles.size() != second.triangles.size()) {
    return false;
  }
  for (std::size_t triangle = 0U; triangle < first.triangles.size();
       ++triangle) {
    const TriangleEvidence &a = first.triangles[triangle];
    const TriangleEvidence &b = second.triangles[triangle];
    if (a.id != b.id || a.vertices.size() != b.vertices.size() ||
        !same_real3(a.normal, b.normal) ||
        double_bits(a.area) != double_bits(b.area)) {
      return false;
    }
    for (std::size_t vertex = 0U; vertex < a.vertices.size(); ++vertex) {
      if (!same_real3(a.vertices[vertex], b.vertices[vertex])) {
        return false;
      }
    }
  }
  return true;
}

SurfaceEvidence capture_surface(const ImmersedSurface &surface) {
  SurfaceEvidence result;
  result.vertex_count = surface.vertex_count();
  result.minimum = surface.bounding_box_min_m();
  result.maximum = surface.bounding_box_max_m();
  result.volume = surface.closed_volume_m3();
  result.fingerprint = surface.fingerprint();
  result.triangles.reserve(surface.triangle_count());
  for (std::size_t index = 0U; index < surface.triangle_count(); ++index) {
    const auto &triangle = surface.triangle(index);
    result.triangles.push_back({triangle.id, triangle.vertices_m,
                                triangle.geometric_outward_normal,
                                triangle.area_m2});
  }
  return result;
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

std::filesystem::path
write_ascii_case(const Stage3TemporaryDirectory &directory,
                 const std::string &name,
                 const std::vector<StlFixtureTriangle> &triangles) {
  const auto path = directory.path() / (name + ".stl");
  hundun::test::write_text(path, hundun::test::ascii_stl(triangles, name));
  return path;
}

void reverse_all(std::vector<StlFixtureTriangle> &triangles) {
  for (StlFixtureTriangle &triangle : triangles) {
    std::swap(triangle.vertices[1], triangle.vertices[2]);
  }
}

void replace_vertex(std::vector<StlFixtureTriangle> &triangles, Real3 old_value,
                    Real3 new_value) {
  for (StlFixtureTriangle &triangle : triangles) {
    for (Real3 &vertex : triangle.vertices) {
      if (vertex.x == old_value.x && vertex.y == old_value.y &&
          vertex.z == old_value.z) {
        vertex = new_value;
      }
    }
  }
}

void check_tetrahedron(const ImmersedSurface &surface, double scale) {
  HUNDUN_CHECK(surface.vertex_count() == 4U);
  HUNDUN_CHECK(surface.triangle_count() == 4U);
  const Real3 minimum = surface.bounding_box_min_m();
  const Real3 maximum = surface.bounding_box_max_m();
  HUNDUN_CHECK(minimum.x == 0.0 && minimum.y == 0.0 && minimum.z == 0.0);
  HUNDUN_CHECK(maximum.x == scale && maximum.y == scale && maximum.z == scale);
  HUNDUN_CHECK(near(surface.closed_volume_m3(), scale * scale * scale / 6.0,
                    scale * scale * scale));

  const Real3 interior{0.25 * scale, 0.25 * scale, 0.25 * scale};
  for (std::size_t index = 0; index < surface.triangle_count(); ++index) {
    const auto &triangle = surface.triangle(index);
    HUNDUN_CHECK(triangle.id == index);
    const Real3 centre{(triangle.vertices_m[0].x + triangle.vertices_m[1].x +
                        triangle.vertices_m[2].x) /
                           3.0,
                       (triangle.vertices_m[0].y + triangle.vertices_m[1].y +
                        triangle.vertices_m[2].y) /
                           3.0,
                       (triangle.vertices_m[0].z + triangle.vertices_m[1].z +
                        triangle.vertices_m[2].z) /
                           3.0};
    HUNDUN_CHECK(dot(triangle.geometric_outward_normal,
                     subtract(centre, interior)) > 0.0);
    HUNDUN_CHECK(near(std::hypot(triangle.geometric_outward_normal.x,
                                 triangle.geometric_outward_normal.y,
                                 triangle.geometric_outward_normal.z),
                      1.0));
    HUNDUN_CHECK(triangle.area_m2 > 0.0);
  }
}

void check_coplanar_shared_vertex_tangent_cone() {
  const Real3 shared{0.47000000000000003, 0.26000000000000001, 0.62};
  const Real3 a1_outward{0.44, 0.20000000000000001,
                         0.56000000000000005};
  const SurfaceTriangle a1_target{
      4U,
      {shared,
       Real3{0.5, 0.28999999999999998, 0.62},
       Real3{0.5, 0.32000000000000001, 0.67999999999999994}},
      Real3{0.6666666666666663, -0.66666666666666663,
            0.33333333333333409},
      1.0};

  HUNDUN_CHECK(!ImmersedTestAccess::coplanar_contact_is_forbidden(
      a1_target, a1_outward, shared, false, true));
  HUNDUN_CHECK(ImmersedTestAccess::coplanar_contact_is_forbidden(
      a1_target, a1_outward, shared, false, false));
  HUNDUN_CHECK(ImmersedTestAccess::coplanar_contact_is_forbidden(
      a1_target, a1_outward, shared, true, false));
  HUNDUN_CHECK(ImmersedTestAccess::coplanar_contact_is_forbidden(
      a1_target, a1_outward, shared, true, true));
  const Real3 coordinate_near_nonshared{
      std::nextafter(shared.x, std::numeric_limits<double>::infinity()),
      shared.y, shared.z};
  HUNDUN_CHECK(ImmersedTestAccess::coplanar_contact_is_forbidden(
      a1_target, a1_outward, coordinate_near_nonshared, false, false));

  const SurfaceTriangle xy_target{
      0U,
      {Real3{0.0, 0.0, 0.0}, Real3{1.0, 0.0, 0.0},
       Real3{0.0, 1.0, 0.0}},
      Real3{0.0, 0.0, 1.0},
      0.5};
  constexpr double epsilon = std::numeric_limits<double>::epsilon();
  constexpr double rejected_parameter_window = 256.0 * epsilon;
  constexpr double below_window_width = 128.0 * epsilon;
  constexpr double above_window_width = 512.0 * epsilon;
  static_assert(below_window_width < rejected_parameter_window);
  static_assert(above_window_width > rejected_parameter_window);
  const auto inward_endpoint = [](double overlap_parameter_width) {
    const double extent = 1.0 / (2.0 * overlap_parameter_width);
    return Real3{extent, extent, 0.0};
  };
  HUNDUN_CHECK(ImmersedTestAccess::coplanar_contact_is_forbidden(
      xy_target, {0.0, 0.0, 0.0}, inward_endpoint(below_window_width), true,
      false));
  HUNDUN_CHECK(ImmersedTestAccess::coplanar_contact_is_forbidden(
      xy_target, {0.0, 0.0, 0.0}, inward_endpoint(above_window_width), true,
      false));

  HUNDUN_CHECK(ImmersedTestAccess::coplanar_contact_is_forbidden(
      xy_target, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, true, false));
  HUNDUN_CHECK(!ImmersedTestAccess::coplanar_contact_is_forbidden(
      xy_target, {0.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}, true, false));

  const SurfaceTriangle dropped_y_target{
      0U,
      {Real3{0.0, 0.0, 0.0}, Real3{1.0, 0.0, 0.0},
       Real3{0.0, 0.0, 1.0}},
      Real3{0.0, -1.0, 0.0},
      0.5};
  HUNDUN_CHECK(!ImmersedTestAccess::coplanar_contact_is_forbidden(
      dropped_y_target, {0.0, 0.0, 0.0}, {-1.0, 0.0, 0.0}, true, false));
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  auto mpi = MpiContext::duplicate(MPI_COMM_WORLD);
  HUNDUN_CHECK(mpi.size() == 1);
  Stage3TemporaryDirectory directory{"surface-unit"};
  check_coplanar_shared_vertex_tangent_cone();

  const auto tetrahedra = hundun::test::outward_tetrahedron();
  const auto ascii_path =
      write_ascii_case(directory, "tetra-ascii", tetrahedra);
  const auto ascii = ImmersedSurface::load_collective(ascii_path, 1.0, mpi, 0);
  check_tetrahedron(ascii, 1.0);
  const SurfaceEvidence ascii_evidence = capture_surface(ascii);
  const auto ascii_repeat =
      ImmersedSurface::load_collective(ascii_path, 1.0, mpi, 0);
  HUNDUN_CHECK(
      same_surface_evidence(capture_surface(ascii_repeat), ascii_evidence));
  SurfaceEvidence exact_copy = ascii_evidence;
  HUNDUN_CHECK(same_surface_evidence(ascii_evidence, exact_copy));
  exact_copy.volume = std::nextafter(exact_copy.volume, 1.0);
  HUNDUN_CHECK(!same_surface_evidence(ascii_evidence, exact_copy));
  exact_copy = ascii_evidence;
  exact_copy.triangles.front().vertices.front().x =
      std::nextafter(exact_copy.triangles.front().vertices.front().x, 1.0);
  HUNDUN_CHECK(!same_surface_evidence(ascii_evidence, exact_copy));

  const auto binary_path = directory.path() / "tetra-binary.stl";
  hundun::test::write_bytes(binary_path,
                            hundun::test::binary_stl(tetrahedra, true));
  const auto binary =
      ImmersedSurface::load_collective(binary_path, 1.0, mpi, 0);
  check_tetrahedron(binary, 1.0);
  HUNDUN_CHECK(same_surface_evidence(capture_surface(binary), ascii_evidence));

  const auto cube_ascii_path =
      write_ascii_case(directory, "cube-ascii", hundun::test::outward_cube());
  const auto cube_binary_path = directory.path() / "cube-binary.stl";
  hundun::test::write_bytes(
      cube_binary_path,
      hundun::test::binary_stl(hundun::test::outward_cube(), true));
  const auto cube_ascii =
      ImmersedSurface::load_collective(cube_ascii_path, 1.0, mpi, 0);
  const auto cube_binary =
      ImmersedSurface::load_collective(cube_binary_path, 1.0, mpi, 0);
  HUNDUN_CHECK(same_surface_evidence(capture_surface(cube_ascii),
                                     capture_surface(cube_binary)));

  const auto scaled = ImmersedSurface::load_collective(ascii_path, 2.0, mpi, 0);
  check_tetrahedron(scaled, 2.0);
  HUNDUN_CHECK(scaled.fingerprint() != ascii.fingerprint());

  auto inward_triangles = tetrahedra;
  reverse_all(inward_triangles);
  const auto inward_path =
      write_ascii_case(directory, "tetra-inward", inward_triangles);
  const auto inward =
      ImmersedSurface::load_collective(inward_path, 1.0, mpi, 0);
  check_tetrahedron(inward, 1.0);
  HUNDUN_CHECK(same_surface_evidence(capture_surface(inward), ascii_evidence));

  const auto translated =
      hundun::test::translated(tetrahedra, Real3{1.0e12, -1.0e12, 1.0e12});
  const auto translated_path =
      write_ascii_case(directory, "tetra-translated", translated);
  const auto translated_surface =
      ImmersedSurface::load_collective(translated_path, 1.0, mpi, 0);
  HUNDUN_CHECK(near(translated_surface.closed_volume_m3(), 1.0 / 6.0));

  auto permuted = tetrahedra;
  std::rotate(permuted.begin(), permuted.begin() + 1, permuted.end());
  const auto permuted_path =
      write_ascii_case(directory, "tetra-permuted", permuted);
  const auto permuted_surface =
      ImmersedSurface::load_collective(permuted_path, 1.0, mpi, 0);
  HUNDUN_CHECK(permuted_surface.fingerprint() != ascii.fingerprint());
  check_tetrahedron(permuted_surface, 1.0);

  auto one_ulp = tetrahedra;
  replace_vertex(one_ulp, {1.0, 0.0, 0.0},
                 {std::nextafter(1.0, 2.0), 0.0, 0.0});
  const auto one_ulp_path =
      write_ascii_case(directory, "tetra-one-ulp", one_ulp);
  const auto one_ulp_surface =
      ImmersedSurface::load_collective(one_ulp_path, 1.0, mpi, 0);
  HUNDUN_CHECK(one_ulp_surface.fingerprint() != ascii.fingerprint());

  auto within_weld = tetrahedra;
  within_weld[1].vertices[0].x = 1.0e-14;
  const auto within_weld_path =
      write_ascii_case(directory, "tetra-weld", within_weld);
  const auto within_weld_surface =
      ImmersedSurface::load_collective(within_weld_path, 1.0, mpi, 0);
  HUNDUN_CHECK(within_weld_surface.fingerprint() == ascii.fingerprint());

  auto within_weld_extreme = tetrahedra;
  within_weld_extreme[1].vertices[0].x = -1.0e-14;
  const auto within_weld_extreme_path =
      write_ascii_case(directory, "tetra-weld-extreme", within_weld_extreme);
  const auto within_weld_extreme_surface =
      ImmersedSurface::load_collective(within_weld_extreme_path, 1.0, mpi, 0);
  HUNDUN_CHECK(same_surface_evidence(
      capture_surface(within_weld_extreme_surface), ascii_evidence));

  auto outside_weld = tetrahedra;
  outside_weld[1].vertices[0].x = 1.0e-10;
  const auto outside_weld_path =
      write_ascii_case(directory, "tetra-no-weld", outside_weld);
  expect_error(
      [&] {
        (void)ImmersedSurface::load_collective(outside_weld_path, 1.0, mpi, 0);
      },
      "open edge");

  expect_error([&] { (void)ascii.triangle(ascii.triangle_count()); },
               "out of range");
  for (const double invalid :
       {0.0, -1.0, std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::quiet_NaN()}) {
    expect_error(
        [&] {
          (void)ImmersedSurface::load_collective(ascii_path, invalid, mpi, 0);
        },
        "length scale");
  }

  auto truncated_binary = hundun::test::binary_stl(tetrahedra);
  truncated_binary.pop_back();
  const auto truncated_path = directory.path() / "truncated.stl";
  hundun::test::write_bytes(truncated_path, truncated_binary);
  expect_error(
      [&] {
        (void)ImmersedSurface::load_collective(truncated_path, 1.0, mpi, 0);
      },
      "truncated");

  auto trailing_binary = hundun::test::binary_stl(tetrahedra);
  trailing_binary.push_back(0U);
  const auto trailing_binary_path = directory.path() / "trailing-binary.stl";
  hundun::test::write_bytes(trailing_binary_path, trailing_binary);
  expect_error(
      [&] {
        (void)ImmersedSurface::load_collective(trailing_binary_path, 1.0, mpi,
                                               0);
      },
      "trailing bytes");

  auto oversized_count = hundun::test::binary_stl(tetrahedra);
  oversized_count[80] = 0U;
  oversized_count[81] = 0U;
  oversized_count[82] = 0U;
  oversized_count[83] = 0x80U;
  const auto oversized_path = directory.path() / "oversized-count.stl";
  hundun::test::write_bytes(oversized_path, oversized_count);
  expect_error(
      [&] {
        (void)ImmersedSurface::load_collective(oversized_path, 1.0, mpi, 0);
      },
      "2^31-1");

  auto nonfinite_binary = hundun::test::binary_stl(tetrahedra);
  nonfinite_binary[84] = 0U;
  nonfinite_binary[85] = 0U;
  nonfinite_binary[86] = 0x80U;
  nonfinite_binary[87] = 0x7fU;
  const auto nonfinite_path = directory.path() / "nonfinite.stl";
  hundun::test::write_bytes(nonfinite_path, nonfinite_binary);
  expect_error(
      [&] {
        (void)ImmersedSurface::load_collective(nonfinite_path, 1.0, mpi, 0);
      },
      "non-finite file normal");

  auto nonfinite_coordinate = hundun::test::binary_stl(tetrahedra);
  nonfinite_coordinate[96] = 0U;
  nonfinite_coordinate[97] = 0U;
  nonfinite_coordinate[98] = 0x80U;
  nonfinite_coordinate[99] = 0x7fU;
  const auto nonfinite_coordinate_path =
      directory.path() / "nonfinite-coordinate.stl";
  hundun::test::write_bytes(nonfinite_coordinate_path, nonfinite_coordinate);
  expect_error(
      [&] {
        (void)ImmersedSurface::load_collective(nonfinite_coordinate_path, 1.0,
                                               mpi, 0);
      },
      "non-finite coordinate");

  for (const std::string token : {"nan", "inf"}) {
    std::string invalid_number = hundun::test::ascii_stl(tetrahedra);
    const auto vertex = invalid_number.find("vertex ");
    HUNDUN_CHECK(vertex != std::string::npos);
    const auto number_begin = vertex + std::string{"vertex "}.size();
    const auto number_end = invalid_number.find(' ', number_begin);
    HUNDUN_CHECK(number_end != std::string::npos);
    invalid_number.replace(number_begin, number_end - number_begin, token);
    const auto invalid_number_path =
        directory.path() / ("invalid-" + token + ".stl");
    hundun::test::write_text(invalid_number_path, invalid_number);
    expect_error(
        [&] {
          (void)ImmersedSurface::load_collective(invalid_number_path, 1.0, mpi,
                                                 0);
        },
        "invalid numeric token");
  }

  const auto malformed_path = directory.path() / "malformed-ascii.stl";
  std::string malformed = hundun::test::ascii_stl(tetrahedra);
  const auto facet = malformed.find("facet normal");
  HUNDUN_CHECK(facet != std::string::npos);
  malformed.replace(facet, std::string{"facet"}.size(), "face ");
  hundun::test::write_text(malformed_path, malformed);
  expect_error(
      [&] {
        (void)ImmersedSurface::load_collective(malformed_path, 1.0, mpi, 0);
      },
      "expected facet normal");

  const auto trailing_ascii_path = directory.path() / "trailing-ascii.stl";
  hundun::test::write_text(trailing_ascii_path,
                           hundun::test::ascii_stl(tetrahedra) +
                               "unexpected\n");
  expect_error(
      [&] {
        (void)ImmersedSurface::load_collective(trailing_ascii_path, 1.0, mpi,
                                               0);
      },
      "trailing tokens");

  auto zero_area = tetrahedra;
  zero_area[0].vertices[2] = zero_area[0].vertices[1];
  const auto zero_area_path =
      write_ascii_case(directory, "zero-area", zero_area);
  expect_error(
      [&] {
        (void)ImmersedSurface::load_collective(zero_area_path, 1.0, mpi, 0);
      },
      "degenerate");

  auto below_minimum_area = tetrahedra;
  below_minimum_area[0].vertices[1] = {0.0, 1.0e-13, 0.0};
  const auto below_minimum_area_path =
      write_ascii_case(directory, "below-minimum-area", below_minimum_area);
  expect_error(
      [&] {
        (void)ImmersedSurface::load_collective(below_minimum_area_path, 1.0,
                                               mpi, 0);
      },
      "area");

  auto collinear = tetrahedra;
  collinear[0].vertices[1] = {0.5, 0.0, 0.0};
  const auto collinear_path =
      write_ascii_case(directory, "collinear", collinear);
  expect_error(
      [&] {
        (void)ImmersedSurface::load_collective(collinear_path, 1.0, mpi, 0);
      },
      "area");

  auto open = tetrahedra;
  open.pop_back();
  const auto open_path = write_ascii_case(directory, "open", open);
  expect_error(
      [&] { (void)ImmersedSurface::load_collective(open_path, 1.0, mpi, 0); },
      "open edge");

  auto nonmanifold = tetrahedra;
  StlFixtureTriangle extra;
  extra.vertices = {Real3{0.0, 0.0, 0.0}, Real3{0.0, 1.0, 0.0},
                    Real3{-1.0, 0.5, 0.5}};
  nonmanifold.push_back(extra);
  const auto nonmanifold_path =
      write_ascii_case(directory, "nonmanifold", nonmanifold);
  expect_error(
      [&] {
        (void)ImmersedSurface::load_collective(nonmanifold_path, 1.0, mpi, 0);
      },
      "non-manifold edge");

  auto duplicate = tetrahedra;
  duplicate.push_back(duplicate.front());
  const auto duplicate_path =
      write_ascii_case(directory, "duplicate", duplicate);
  expect_error(
      [&] {
        (void)ImmersedSurface::load_collective(duplicate_path, 1.0, mpi, 0);
      },
      "duplicate triangle");

  auto reversed_face = tetrahedra;
  std::swap(reversed_face[0].vertices[1], reversed_face[0].vertices[2]);
  const auto reversed_face_path =
      write_ascii_case(directory, "reversed-face", reversed_face);
  expect_error(
      [&] {
        (void)ImmersedSurface::load_collective(reversed_face_path, 1.0, mpi, 0);
      },
      "orientation");

  auto zero_volume = tetrahedra;
  replace_vertex(zero_volume, {0.0, 0.0, 1.0}, {1.0, 1.0, 0.0});
  const auto zero_volume_path =
      write_ascii_case(directory, "zero-volume", zero_volume);
  expect_error(
      [&] {
        (void)ImmersedSurface::load_collective(zero_volume_path, 1.0, mpi, 0);
      },
      "zero enclosed volume");

  auto disconnected = tetrahedra;
  const auto second =
      hundun::test::translated(tetrahedra, Real3{3.0, 0.0, 0.0});
  disconnected.insert(disconnected.end(), second.begin(), second.end());
  const auto disconnected_path =
      write_ascii_case(directory, "disconnected", disconnected);
  expect_error(
      [&] {
        (void)ImmersedSurface::load_collective(disconnected_path, 1.0, mpi, 0);
      },
      "multiple connected components");

  auto self_intersecting = hundun::test::outward_cube();
  replace_vertex(self_intersecting, {1.0, 1.0, 1.0}, {0.5, 0.5, -0.5});
  const auto self_intersecting_path =
      write_ascii_case(directory, "self-intersecting", self_intersecting);
  expect_error(
      [&] {
        (void)ImmersedSurface::load_collective(self_intersecting_path, 1.0, mpi,
                                               0);
      },
      "self-intersection");

  auto shared_vertex_intersection = hundun::test::outward_cube();
  replace_vertex(shared_vertex_intersection, {1.0, 0.0, 1.0}, {0.5, 0.5, 0.0});
  const auto shared_vertex_intersection_path = write_ascii_case(
      directory, "shared-vertex-intersection", shared_vertex_intersection);
  expect_error(
      [&] {
        (void)ImmersedSurface::load_collective(shared_vertex_intersection_path,
                                               1.0, mpi, 0);
      },
      "self-intersection");

  return 0;
}

// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "ib_periodic_surface_window_detail.hpp"
#include "hundun/ib_surface.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/test_main.hpp"

#include <array>
#include <cmath>

using hundun::immersed::SurfaceTriangle;
using hundun::immersed::detail::PeriodicCellMapper;
using hundun::immersed::detail::PeriodicSurfaceWindow;
using hundun::runtime::Int3;
using hundun::runtime::Real3;

namespace {

SurfaceTriangle triangle(double z0, double z1, double z2) {
  SurfaceTriangle result;
  result.vertices_m = {Real3{0.5, 0.0, z0}, Real3{0.5, 0.1, z1},
                       Real3{0.5, 0.0, z2}};
  return result;
}

} // namespace

void run(const hundun::runtime::MpiContext &mpi) {
  const PeriodicSurfaceWindow window(
      {-5.0, -10.0, 0.0}, {15.0, 10.0, std::acos(-1.0)},
      {false, false, true});
  const double h = 0.05;

  HUNDUN_CHECK(window.bounding_box_is_admissible(
      {-0.5, -0.5, -1.0}, {0.5, 0.5, std::acos(-1.0) + 1.0}, h));
  HUNDUN_CHECK(!window.bounding_box_is_admissible(
      {-0.5, -0.5, 0.0}, {0.5, 0.5, std::acos(-1.0)}, h));
  HUNDUN_CHECK(!window.bounding_box_is_admissible(
      {-0.5, -0.5, -1.0}, {0.5, 0.5, std::acos(-1.0) - 0.1}, h));
  HUNDUN_CHECK(!window.bounding_box_is_admissible(
      {-5.0, -0.5, -1.0}, {0.5, 0.5, std::acos(-1.0) + 1.0}, h));

  const auto inside = triangle(0.0, 0.2, 0.1);
  const auto below = triangle(-1.0, -0.8, -0.9);
  const auto straddling = triangle(-0.1, 0.1, 0.0);
  HUNDUN_CHECK(window.triangle_is_split_at_periodic_planes(inside));
  HUNDUN_CHECK(window.triangle_is_split_at_periodic_planes(below));
  HUNDUN_CHECK(!window.triangle_is_split_at_periodic_planes(straddling));
  HUNDUN_CHECK(window.triangle_is_active(inside));
  HUNDUN_CHECK(!window.triangle_is_active(below));

  const PeriodicCellMapper cells({480, 480, 48}, {false, false, true},
                                 {20.0, 20.0, std::acos(-1.0)});
  const auto below_image = cells.image({113, 229, -1});
  HUNDUN_CHECK(below_image.has_value());
  HUNDUN_CHECK(below_image->canonical.z == 47);
  HUNDUN_CHECK(below_image->image.z == -1);
  HUNDUN_CHECK(std::abs(below_image->shift_m.z + std::acos(-1.0)) < 1.0e-14);
  const auto nearest = cells.nearest_image({113, 229, 47}, {113, 229, 0});
  HUNDUN_CHECK(nearest.canonical.z == 47);
  HUNDUN_CHECK(nearest.image.z == -1);
  const auto lower_box = cells.nearest_image_to_box(
      {113, 229, 47}, {{0, 0, 0}, {480, 480, 12}});
  HUNDUN_CHECK(lower_box.canonical.z == 47);
  HUNDUN_CHECK(lower_box.image.z == -1);
  const auto upper_box = cells.nearest_image_to_box(
      {113, 229, 0}, {{0, 0, 36}, {480, 480, 48}});
  HUNDUN_CHECK(upper_box.canonical.z == 0);
  HUNDUN_CHECK(upper_box.image.z == 48);
  const auto full_box = cells.nearest_image_to_box(
      {113, 229, 47}, {{0, 0, 0}, {480, 480, 48}});
  HUNDUN_CHECK(full_box.image.z == 47);
  HUNDUN_CHECK(!cells.image({-1, 229, 0}).has_value());

  const auto benchmark_surface =
      hundun::immersed::ImmersedSurface::load_collective(
          HUNDUN_CYLINDER_CASE_STL, 1.0, mpi, 0);
  const PeriodicSurfaceWindow benchmark_window(
      {-1.0, -1.0, 0.0}, {1.0, 1.0, std::acos(-1.0)},
      {false, true, true});
  HUNDUN_CHECK(benchmark_surface.bounding_box_min_m().z < 0.0);
  HUNDUN_CHECK(benchmark_surface.bounding_box_max_m().z > std::acos(-1.0));
  for (hundun::immersed::TriangleId triangle = 0U;
       triangle < benchmark_surface.triangle_count(); ++triangle)
    HUNDUN_CHECK(benchmark_window.triangle_is_split_at_periodic_planes(
        benchmark_surface.triangle(triangle)));
}

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  const auto mpi =
      hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  return hundun::test::run([&] { run(mpi); });
}

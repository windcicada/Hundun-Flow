// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/ib_wall_force.hpp"

#include <cstdint>
#include <vector>

namespace hundun::immersed::detail {

#ifdef HUNDUN_IMMERSED_ENABLE_TEST_ACCESS
struct WallForcePointSnapshot final {
  ImmersedLinkId link{};
  TriangleId triangle{};
  std::uint32_t point_index{};
  runtime::Real3 position_m{};
  runtime::Real3 solid_to_fluid_normal{};
  double weight_m2{};
  runtime::Real3 pressure_force_N{};
  runtime::Real3 viscous_force_N{};
  runtime::Real3 total_force_N{};
};

struct WallForceTrace final {
  WallForceSample reduced;
  std::vector<WallForcePointSnapshot> local_points;
};

WallForceTrace trace_wall_force_for_test(
    const WallForceIntegrator &,
    const runtime::FieldView<const double> &mechanical_pressure,
    const runtime::FieldView<const double> &velocity,
    const runtime::FieldView<const double> &velocity_gradient,
    const runtime::FieldView<const double> &mu_eff_by_cell);
#endif

struct WallPressureNormalGradient final {
  std::uint64_t link{};
  double value{};
};

WallForceSample integrate_with_wall_pressure_authority(
    const WallForceIntegrator &,
    const runtime::FieldView<const double> &mechanical_pressure,
    const runtime::FieldView<const double> &velocity,
    const runtime::FieldView<const double> &velocity_gradient,
    const runtime::FieldView<const double> &mu_eff_by_cell,
    const std::vector<WallPressureNormalGradient> &);

} // namespace hundun::immersed::detail

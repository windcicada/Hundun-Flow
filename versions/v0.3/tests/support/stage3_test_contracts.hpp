// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/rt_types.hpp"

#include <array>
#include <cstdint>

namespace hundun::test::stage3 {

enum class BodyKind : std::uint8_t {
  sphere,
  finite_cylinder,
  oblique_rectangular_prism,
  inside_sphere_cavity
};

struct BodySpec final {
  BodyKind kind{BodyKind::sphere};
  runtime::Real3 centre_m{0.5, 0.5, 0.5};
  double radius_m{};
  double length_m{};
  runtime::Real3 axis{};
  runtime::Real3 half_lengths_m{};
  double prism_mms_factor_multiplier{1.0};
};

inline constexpr double kReferenceLengthM = 1.0;
inline constexpr double kReferenceVelocityMPerS = 1.0;
inline constexpr double kReferenceDensityKgPerM3 = 1.0;
inline constexpr double kDynamicViscosityPaS = 0.01;
inline constexpr std::array<double, 3> kWarpAmplitude{0.02, -0.015, 0.01};
inline constexpr double kManufacturedSolveAtol = 1.0e-16;
inline constexpr double kManufacturedSolveRtol = 1.0e-15;
inline constexpr std::uint64_t kManufacturedSolveMaxIterations = 5000U;
inline constexpr std::uint64_t kManufacturedResidualRecomputeInterval = 5U;

BodySpec approved_body(BodyKind);
BodySpec translated_sphere();
BodySpec force_certified_oblique_prism();

} // namespace hundun::test::stage3

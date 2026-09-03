// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/ib_surface.hpp"
#include "hundun/rt_types.hpp"

#include <cstdint>

namespace hundun::flow::detail {

struct ImmersedWallPressureInput final {
  immersed::ImmersedLinkId link{};
  double rho_wall_kg_per_m3{};
  double effective_transformed_measure_m2{};
  runtime::Real3 solid_to_fluid_unit_normal{};
  runtime::Real3 momentum_velocity_correction_m3_s_per_kg{};
  double current_normal_gradient_pa_per_m{};
  double predictor_mass_flux_kg_per_s{};
};

struct ImmersedWallPressureCondition final {
  immersed::ImmersedLinkId link{};
  double predictor_mass_flux_kg_per_s{};
  double correction_coefficient{};
  double current_normal_gradient_pa_per_m{};
  double correction_normal_gradient_pa_per_m{};
  double corrected_normal_gradient_pa_per_m{};
  double corrected_mass_flux_kg_per_s{};
};

ImmersedWallPressureCondition
make_immersed_wall_pressure_condition(const ImmersedWallPressureInput &);

std::uint64_t make_immersed_pressure_revision(
    std::uint64_t previous_revision, std::uint64_t density_fingerprint,
    std::uint64_t momentum_diagonal_fingerprint,
    std::uint64_t geometry_fingerprint,
    std::uint64_t active_layout_fingerprint,
    std::uint64_t wall_condition_fingerprint);

} // namespace hundun::flow::detail

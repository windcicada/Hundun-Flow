// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "src/flow_immersed_piso_detail.hpp"

#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>

namespace {

using hundun::flow::detail::ImmersedWallPressureInput;
using hundun::flow::detail::make_immersed_pressure_revision;
using hundun::flow::detail::make_immersed_wall_pressure_condition;

void run() {
  auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
  const double rank_scale = 1.0 + static_cast<double>(mpi.rank());
  const ImmersedWallPressureInput input{
      7U, 2.0, 3.0, {0.6, 0.8, 0.0}, {0.25, 0.5, 0.75},
      -0.125 * rank_scale, 4.0 * rank_scale};
  const auto condition = make_immersed_wall_pressure_condition(input);
  const double mobility = 0.6 * 0.6 * 0.25 + 0.8 * 0.8 * 0.5;
  const double expected_d = 2.0 * 3.0 * mobility;
  HUNDUN_CHECK(condition.link == 7U);
  HUNDUN_CHECK(condition.correction_coefficient == expected_d);
  HUNDUN_CHECK(condition.correction_normal_gradient_pa_per_m ==
               input.predictor_mass_flux_kg_per_s / expected_d);
  HUNDUN_CHECK(condition.current_normal_gradient_pa_per_m ==
               input.current_normal_gradient_pa_per_m);
  HUNDUN_CHECK(condition.corrected_normal_gradient_pa_per_m ==
               input.current_normal_gradient_pa_per_m +
                   condition.correction_normal_gradient_pa_per_m);
  HUNDUN_CHECK(condition.corrected_mass_flux_kg_per_s == 0.0);
  HUNDUN_CHECK(!std::signbit(condition.corrected_mass_flux_kg_per_s));
  auto changed_density = input;
  changed_density.rho_wall_kg_per_m3 = 4.0;
  HUNDUN_CHECK(make_immersed_wall_pressure_condition(changed_density)
                   .correction_coefficient != condition.correction_coefficient);
  auto changed_measure = input;
  changed_measure.effective_transformed_measure_m2 = 5.0;
  HUNDUN_CHECK(make_immersed_wall_pressure_condition(changed_measure)
                   .correction_coefficient != condition.correction_coefficient);
  auto changed_diagonal = input;
  changed_diagonal.momentum_velocity_correction_m3_s_per_kg.x = 0.5;
  HUNDUN_CHECK(make_immersed_wall_pressure_condition(changed_diagonal)
                   .correction_coefficient != condition.correction_coefficient);

  const auto be =
      make_immersed_wall_pressure_condition(ImmersedWallPressureInput{
          input.link, input.rho_wall_kg_per_m3,
          input.effective_transformed_measure_m2,
          input.solid_to_fluid_unit_normal,
          input.momentum_velocity_correction_m3_s_per_kg,
          input.current_normal_gradient_pa_per_m, 1.25});
  const auto bdf2 =
      make_immersed_wall_pressure_condition(ImmersedWallPressureInput{
          input.link, input.rho_wall_kg_per_m3,
          input.effective_transformed_measure_m2,
          input.solid_to_fluid_unit_normal,
          input.momentum_velocity_correction_m3_s_per_kg,
          input.current_normal_gradient_pa_per_m, 1.75});
  HUNDUN_CHECK(be.correction_normal_gradient_pa_per_m !=
               bdf2.correction_normal_gradient_pa_per_m);

  constexpr std::uint64_t previous = 11U;
  const std::uint64_t revision =
      make_immersed_pressure_revision(previous, 12U, 13U, 14U, 15U, 16U);
  HUNDUN_CHECK(revision != make_immersed_pressure_revision(previous, 17U, 13U,
                                                           14U, 15U, 16U));
  HUNDUN_CHECK(revision != make_immersed_pressure_revision(previous, 12U, 17U,
                                                           14U, 15U, 16U));
  HUNDUN_CHECK(revision != make_immersed_pressure_revision(previous, 12U, 13U,
                                                           17U, 15U, 16U));
  HUNDUN_CHECK(revision != make_immersed_pressure_revision(previous, 12U, 13U,
                                                           14U, 17U, 16U));
  HUNDUN_CHECK(revision != make_immersed_pressure_revision(previous, 12U, 13U,
                                                           14U, 15U, 17U));
  HUNDUN_CHECK(revision !=
               make_immersed_pressure_revision(17U, 12U, 13U, 14U, 15U, 16U));

  for (const auto bad :
       std::array<double, 4>{0.0, -1.0, std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::quiet_NaN()}) {
    bool rejected = false;
    try {
      auto invalid = input;
      invalid.rho_wall_kg_per_m3 = bad;
      static_cast<void>(make_immersed_wall_pressure_condition(invalid));
    } catch (const hundun::runtime::Error &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  }
  for (const auto bad :
       std::array<double, 4>{0.0, -1.0, std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::quiet_NaN()}) {
    bool rejected = false;
    try {
      auto invalid = input;
      invalid.effective_transformed_measure_m2 = bad;
      static_cast<void>(make_immersed_wall_pressure_condition(invalid));
    } catch (const hundun::runtime::Error &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  }
  {
    auto invalid = input;
    invalid.momentum_velocity_correction_m3_s_per_kg = {-1.0, -1.0, -1.0};
    bool rejected = false;
    try {
      static_cast<void>(make_immersed_wall_pressure_condition(invalid));
    } catch (const hundun::runtime::Error &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  }
  {
    auto invalid = input;
    invalid.solid_to_fluid_unit_normal = {1.0, 1.0, 0.0};
    bool rejected = false;
    try {
      static_cast<void>(make_immersed_wall_pressure_condition(invalid));
    } catch (const hundun::runtime::Error &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  }
  for (const auto bad :
       std::array<double, 2>{std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::quiet_NaN()}) {
    bool rejected = false;
    try {
      auto invalid = input;
      invalid.current_normal_gradient_pa_per_m = bad;
      static_cast<void>(make_immersed_wall_pressure_condition(invalid));
    } catch (const hundun::runtime::Error &) {
      rejected = true;
    }
    HUNDUN_CHECK(rejected);
  }
}

} // namespace

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run(run);
}

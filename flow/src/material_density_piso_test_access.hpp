// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#error "material-density PISO test access requires tests-on build"
#endif

#include <cmath>
#include <cstdint>

namespace hundun::flow {
class FlowState;
class FixedStepMaterialDensityFlow;
}

namespace hundun::flow::test {

struct MaterialDensityVortexSource final {
  double x{};
  double y{};
  double z{};
};

class MaterialDensityPisoTestAccess final {
public:
  static constexpr std::uint32_t contract_version() noexcept { return 1U; }

  static MaterialDensityVortexSource vortex_source(double x, double y,
                                                    double mu) noexcept {
    const double density = 1.0 + 0.1 * std::sin(x) * std::sin(y);
    return {density * std::sin(x) * std::cos(x) +
                2.0 * mu * std::sin(x) * std::cos(y),
            density * std::sin(y) * std::cos(y) -
                2.0 * mu * std::cos(x) * std::sin(y),
            0.0};
  }

  static void force_state_diagnostic_identity(FlowState &,
                                               std::uint64_t) noexcept;
  static std::uint64_t
  state_diagnostic_identity(const FlowState &) noexcept;
  static void force_flow_attempt_identity(FixedStepMaterialDensityFlow &,
                                          std::uint64_t) noexcept;
  static void enable_vortex_source(FixedStepMaterialDensityFlow &,
                                   bool) noexcept;
  static double material_face_value(
      double predictor_velocity_interpolation,
      double predictor_momentum_interpolation, double face_density,
      double cell_momentum_n, double face_momentum_n,
      double cell_momentum_n_minus_1, double face_momentum_n_minus_1,
      double mobility, double dt_s, double alpha1, double alpha2,
      double pressure_correction) noexcept;
};

} // namespace hundun::flow::test

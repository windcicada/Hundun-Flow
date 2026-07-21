// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>

namespace hundun::flow::test {

class ConstantDensityPisoTestAccess final {
public:
  static void reset() noexcept;
  static void force_final_continuity_failure(bool enabled) noexcept;
  static void force_final_pressure_failure(bool enabled) noexcept;
  static void force_local_derived_failure(bool enabled) noexcept;
  static void set_provisional_transport_sentinel(bool enabled) noexcept;
  static void set_final_uniform_x_mass_flux(double value) noexcept;
  static double last_momentum_rhs(std::size_t component) noexcept;
  static double last_momentum_diagonal(std::size_t component) noexcept;
  static std::size_t provisional_transport_calls() noexcept;
  static std::size_t final_transport_calls() noexcept;
};

} // namespace hundun::flow::test

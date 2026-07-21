// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>

namespace hundun::flow::test {

enum class MomentumAssemblyMutation {
  none,
  omit_convection,
  omit_viscosity,
  omit_pressure,
  omit_alpha1,
  omit_alpha2,
  replace_history_flux
};

enum class TransportAssemblyMutation {
  none,
  omit_history_spatial,
  use_provisional_flux
};

enum class AttemptFailureStage {
  none,
  after_begin,
  after_momentum,
  after_face_predictor,
  after_corrector_1,
  after_provisional_transport,
  after_corrector_2,
  after_final_transport,
  before_commit
};

class ConstantDensityPisoTestAccess final {
public:
  static void reset() noexcept;
  static void force_final_continuity_failure(bool enabled) noexcept;
  static void force_final_pressure_failure(bool enabled) noexcept;
  static void force_local_derived_failure(bool enabled) noexcept;
  static void force_final_momentum_perturbation(std::size_t component,
                                                double delta) noexcept;
  static void force_final_transport_perturbation(std::size_t field_index,
                                                 double delta) noexcept;
  static void force_final_conservation_failure(bool enabled) noexcept;
  static void set_momentum_assembly_mutation(
      MomentumAssemblyMutation mutation) noexcept;
  static void set_transport_assembly_mutation(
      TransportAssemblyMutation mutation) noexcept;
  static void set_attempt_failure_stage(AttemptFailureStage stage) noexcept;
  static void set_provisional_transport_sentinel(bool enabled) noexcept;
  static void set_final_uniform_x_mass_flux(double value) noexcept;
  static void set_final_uniform_x_mass_flux_override(bool enabled) noexcept;
  static double last_momentum_rhs(std::size_t component) noexcept;
  static double last_momentum_diagonal(std::size_t component) noexcept;
  static std::size_t provisional_transport_calls() noexcept;
  static std::size_t final_transport_calls() noexcept;
  static int last_pressure_constraint_mode() noexcept;
  static int last_pressure_operator_mode() noexcept;
};

} // namespace hundun::flow::test

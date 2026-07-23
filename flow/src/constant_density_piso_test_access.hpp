// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace hundun::flow {

class FixedStepConstantDensityFlow;

namespace test {

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

struct ConservationDiagnostic final {
  double quantity_nm1{};
  double quantity_n{};
  double quantity_np1{};
  double signed_boundary_flux{};
  double absolute_boundary_flux{};
  double boundary_integral{};
  double history_correction{};
  double raw_defect{};
  double relative_defect{};
};

struct MeshWorkspaceSnapshot final {
  std::size_t vector_count{};
  std::size_t total_capacity{};
  std::uintptr_t data_identity{};
};

struct PressureOperatorSnapshot final {
  std::uintptr_t identity{};
  std::uint64_t revision{};
  bool present{};
};

struct ConstantWorkspaceCacheEntry final {
  std::uintptr_t identity{};
  std::size_t capacity{};
};

struct ConstantMomentumOperatorCacheEntry final {
  std::uintptr_t identity{};
  std::uint64_t revision{};
  std::vector<double> diagonal;
};

struct ConstantFacadeCacheSnapshot final {
  std::vector<ConstantWorkspaceCacheEntry> workspaces;
  std::array<ConstantMomentumOperatorCacheEntry, 3> operators{};
  std::size_t operator_count{};
  bool delegated{};
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
  static void set_final_mass_defect_perturbation(double delta) noexcept;
  static void set_final_momentum_norm_squares(
      std::size_t component, double residual_square,
      double scale_square) noexcept;
  static void set_final_transport_norm_squares(
      std::size_t field_index, double residual_square,
      double scale_square) noexcept;
  static void set_momentum_conservation_parts(
      std::size_t component, const std::array<double, 8> &values) noexcept;
  static void force_momentum_conservation_aggregate_overflow(
      std::size_t component, bool enabled) noexcept;
  static void force_transport_conservation_aggregate_overflow(
      std::size_t field_index, bool enabled) noexcept;
  static void set_momentum_assembly_mutation(
      MomentumAssemblyMutation mutation) noexcept;
  static void set_transport_assembly_mutation(
      TransportAssemblyMutation mutation) noexcept;
  static void set_attempt_failure_stage(AttemptFailureStage stage) noexcept;
  static void
  set_pressure_operator_construction_failure_rank(int rank) noexcept;
  static void set_pressure_operator_refresh_failure_rank(int rank) noexcept;
  static void set_provisional_transport_sentinel(bool enabled) noexcept;
  static void set_final_uniform_x_mass_flux(double value) noexcept;
  static void set_final_uniform_x_mass_flux_override(bool enabled) noexcept;
  static double last_momentum_rhs(std::size_t component) noexcept;
  static double last_momentum_diagonal(std::size_t component) noexcept;
  static std::size_t provisional_transport_calls() noexcept;
  static std::size_t final_transport_calls() noexcept;
  static int last_pressure_constraint_mode() noexcept;
  static int last_pressure_operator_mode() noexcept;
  static ConservationDiagnostic last_mass_conservation() noexcept;
  static ConservationDiagnostic
  last_momentum_conservation(std::size_t component) noexcept;
  static ConservationDiagnostic
  last_transport_conservation(std::size_t field_index) noexcept;
  static MeshWorkspaceSnapshot
  mesh_workspace_snapshot(const FixedStepConstantDensityFlow &flow) noexcept;
  static ConstantFacadeCacheSnapshot
  facade_cache_snapshot(const FixedStepConstantDensityFlow &flow);
  static PressureOperatorSnapshot pressure_operator_snapshot(
      const FixedStepConstantDensityFlow &flow) noexcept;
};

} // namespace test
} // namespace hundun::flow

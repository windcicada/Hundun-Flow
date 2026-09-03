// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "flow_adaptive_time_control_detail.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace hundun::flow {

class FixedStepConstantDensityFlow;

namespace detail {
void constant_density_reset_raw() noexcept;
void constant_density_force_final_continuity_failure_raw(bool) noexcept;
void constant_density_force_final_pressure_failure_raw(bool) noexcept;
void constant_density_force_local_derived_failure_raw(bool) noexcept;
void constant_density_force_final_momentum_perturbation_raw(
    std::size_t, double) noexcept;
void constant_density_force_final_transport_perturbation_raw(
    std::size_t, double) noexcept;
void constant_density_force_final_conservation_failure_raw(bool) noexcept;
void constant_density_set_final_mass_defect_perturbation_raw(double) noexcept;
void constant_density_set_final_momentum_norm_squares_raw(
    std::size_t, double, double) noexcept;
void constant_density_set_final_transport_norm_squares_raw(
    std::size_t, double, double) noexcept;
void constant_density_set_momentum_conservation_parts_raw(
    std::size_t, const std::array<double, 8> &) noexcept;
void constant_density_force_momentum_conservation_aggregate_overflow_raw(
    std::size_t, bool) noexcept;
void constant_density_force_transport_conservation_aggregate_overflow_raw(
    std::size_t, bool) noexcept;
void constant_density_set_momentum_assembly_mutation_raw(
    std::uint8_t) noexcept;
void constant_density_set_transport_assembly_mutation_raw(
    std::uint8_t) noexcept;
void constant_density_set_attempt_failure_stage_raw(std::uint8_t) noexcept;
void constant_density_set_pressure_operator_construction_failure_rank_raw(
    int) noexcept;
void constant_density_set_pressure_operator_refresh_failure_rank_raw(
    int) noexcept;
void constant_density_set_provisional_transport_sentinel_raw(bool) noexcept;
void constant_density_set_final_uniform_x_mass_flux_raw(double) noexcept;
void constant_density_set_final_uniform_x_mass_flux_override_raw(bool) noexcept;
double constant_density_last_momentum_rhs_raw(std::size_t) noexcept;
double constant_density_last_momentum_diagonal_raw(std::size_t) noexcept;
std::size_t constant_density_provisional_transport_calls_raw() noexcept;
std::size_t constant_density_final_transport_calls_raw() noexcept;
int constant_density_last_pressure_constraint_mode_raw() noexcept;
int constant_density_last_pressure_operator_mode_raw() noexcept;
void constant_density_last_mass_conservation_raw(
    std::array<double, 9> &) noexcept;
void constant_density_last_momentum_conservation_raw(
    std::size_t, std::array<double, 9> &) noexcept;
void constant_density_last_transport_conservation_raw(
    std::size_t, std::array<double, 9> &) noexcept;
} // namespace detail

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
  static void reset() noexcept { detail::constant_density_reset_raw(); }
  static void force_final_continuity_failure(bool enabled) noexcept {
    detail::constant_density_force_final_continuity_failure_raw(enabled);
  }
  static void force_final_pressure_failure(bool enabled) noexcept {
    detail::constant_density_force_final_pressure_failure_raw(enabled);
  }
  static void force_local_derived_failure(bool enabled) noexcept {
    detail::constant_density_force_local_derived_failure_raw(enabled);
  }
  static void force_final_momentum_perturbation(std::size_t component,
                                                double delta) noexcept {
    detail::constant_density_force_final_momentum_perturbation_raw(component,
                                                                   delta);
  }
  static void force_final_transport_perturbation(std::size_t field_index,
                                                 double delta) noexcept {
    detail::constant_density_force_final_transport_perturbation_raw(field_index,
                                                                    delta);
  }
  static void force_final_conservation_failure(bool enabled) noexcept {
    detail::constant_density_force_final_conservation_failure_raw(enabled);
  }
  static void set_final_mass_defect_perturbation(double delta) noexcept {
    detail::constant_density_set_final_mass_defect_perturbation_raw(delta);
  }
  static void set_final_momentum_norm_squares(
      std::size_t component, double residual_square,
      double scale_square) noexcept {
    detail::constant_density_set_final_momentum_norm_squares_raw(
        component, residual_square, scale_square);
  }
  static void set_final_transport_norm_squares(
      std::size_t field_index, double residual_square,
      double scale_square) noexcept {
    detail::constant_density_set_final_transport_norm_squares_raw(
        field_index, residual_square, scale_square);
  }
  static void set_momentum_conservation_parts(
      std::size_t component, const std::array<double, 8> &values) noexcept {
    detail::constant_density_set_momentum_conservation_parts_raw(component,
                                                                 values);
  }
  static void force_momentum_conservation_aggregate_overflow(
      std::size_t component, bool enabled) noexcept {
    detail::constant_density_force_momentum_conservation_aggregate_overflow_raw(
        component, enabled);
  }
  static void force_transport_conservation_aggregate_overflow(
      std::size_t field_index, bool enabled) noexcept {
    detail::constant_density_force_transport_conservation_aggregate_overflow_raw(
        field_index, enabled);
  }
  static void set_momentum_assembly_mutation(
      MomentumAssemblyMutation mutation) noexcept {
    detail::constant_density_set_momentum_assembly_mutation_raw(
        static_cast<std::uint8_t>(mutation));
  }
  static void set_transport_assembly_mutation(
      TransportAssemblyMutation mutation) noexcept {
    detail::constant_density_set_transport_assembly_mutation_raw(
        static_cast<std::uint8_t>(mutation));
  }
  static void set_attempt_failure_stage(AttemptFailureStage stage) noexcept {
    detail::constant_density_set_attempt_failure_stage_raw(
        static_cast<std::uint8_t>(stage));
  }
  static void
  set_pressure_operator_construction_failure_rank(int rank) noexcept {
    detail::constant_density_set_pressure_operator_construction_failure_rank_raw(
        rank);
  }
  static void set_pressure_operator_refresh_failure_rank(int rank) noexcept {
    detail::constant_density_set_pressure_operator_refresh_failure_rank_raw(
        rank);
  }
  static void set_provisional_transport_sentinel(bool enabled) noexcept {
    detail::constant_density_set_provisional_transport_sentinel_raw(enabled);
  }
  static void set_final_uniform_x_mass_flux(double value) noexcept {
    detail::constant_density_set_final_uniform_x_mass_flux_raw(value);
  }
  static void set_final_uniform_x_mass_flux_override(bool enabled) noexcept {
    detail::constant_density_set_final_uniform_x_mass_flux_override_raw(enabled);
  }
  static double last_momentum_rhs(std::size_t component) noexcept {
    return detail::constant_density_last_momentum_rhs_raw(component);
  }
  static double last_momentum_diagonal(std::size_t component) noexcept {
    return detail::constant_density_last_momentum_diagonal_raw(component);
  }
  static std::size_t provisional_transport_calls() noexcept {
    return detail::constant_density_provisional_transport_calls_raw();
  }
  static std::size_t final_transport_calls() noexcept {
    return detail::constant_density_final_transport_calls_raw();
  }
  static int last_pressure_constraint_mode() noexcept {
    return detail::constant_density_last_pressure_constraint_mode_raw();
  }
  static int last_pressure_operator_mode() noexcept {
    return detail::constant_density_last_pressure_operator_mode_raw();
  }
  static ConservationDiagnostic last_mass_conservation() noexcept {
    std::array<double, 9> values{};
    detail::constant_density_last_mass_conservation_raw(values);
    return {values[0], values[1], values[2], values[3], values[4],
            values[5], values[6], values[7], values[8]};
  }
  static ConservationDiagnostic
  last_momentum_conservation(std::size_t component) noexcept {
    std::array<double, 9> values{};
    detail::constant_density_last_momentum_conservation_raw(component, values);
    return {values[0], values[1], values[2], values[3], values[4],
            values[5], values[6], values[7], values[8]};
  }
  static ConservationDiagnostic
  last_transport_conservation(std::size_t field_index) noexcept {
    std::array<double, 9> values{};
    detail::constant_density_last_transport_conservation_raw(field_index,
                                                             values);
    return {values[0], values[1], values[2], values[3], values[4],
            values[5], values[6], values[7], values[8]};
  }
  static MeshWorkspaceSnapshot
  mesh_workspace_snapshot(const FixedStepConstantDensityFlow &flow) noexcept {
    MeshWorkspaceSnapshot result;
    detail::AdaptiveTimeControlAccess::constant_workspace_values_raw(
        flow, result.vector_count, result.total_capacity,
        result.data_identity);
    return result;
  }
  static ConstantFacadeCacheSnapshot
  facade_cache_snapshot(const FixedStepConstantDensityFlow &flow) {
    ConstantFacadeCacheSnapshot result;
    std::vector<std::uintptr_t> identities;
    std::vector<std::size_t> capacities;
    std::array<std::uintptr_t, 3> operator_ids{};
    std::array<std::uint64_t, 3> operator_revisions{};
    std::array<std::vector<double>, 3> diagonals{};
    detail::AdaptiveTimeControlAccess::constant_cache_values_raw(
        flow, identities, capacities, operator_ids, operator_revisions,
        diagonals, result.operator_count, result.delegated);
    result.workspaces.reserve(identities.size());
    for (std::size_t index = 0; index < identities.size(); ++index)
      result.workspaces.push_back({identities[index], capacities[index]});
    for (std::size_t index = 0; index < result.operator_count; ++index)
      result.operators[index] = {operator_ids[index], operator_revisions[index],
                                 std::move(diagonals[index])};
    return result;
  }
  static PressureOperatorSnapshot pressure_operator_snapshot(
      const FixedStepConstantDensityFlow &flow) noexcept {
    PressureOperatorSnapshot result;
    detail::AdaptiveTimeControlAccess::constant_pressure_operator_values_raw(
        flow, result.identity, result.revision, result.present);
    return result;
  }
};

} // namespace test
} // namespace hundun::flow

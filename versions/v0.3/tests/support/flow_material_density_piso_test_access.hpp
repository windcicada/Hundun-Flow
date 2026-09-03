// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#ifndef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#error "material-density PISO test access requires tests-on build"
#endif

#include "flow_adaptive_time_control_detail.hpp"
#include "flow_checkpoint_v2_detail.hpp"
#include "flow_density_closure_detail.hpp"
#include "hundun/rt_halo_performance_counters.hpp"
#include "hundun/rt_error.hpp"

#include <cmath>
#include <array>
#include <cstdint>
#include <cstddef>
#include <utility>
#include <vector>

namespace hundun::flow {
class FlowState;
class FixedStepMaterialDensityFlow;
class MaterialDensityStepAttemptReport;
class PisoCoupler;
enum class FlowLayer : std::uint8_t;
struct AcceptedStepMetadata;
}
namespace hundun::runtime {
class MpiContext;
}

namespace hundun::flow::detail {

double material_face_value_raw(
    double predictor_velocity_interpolation,
    double predictor_momentum_interpolation, double face_density,
    double cell_momentum_n, double face_momentum_n,
    double cell_momentum_n_minus_1, double face_momentum_n_minus_1,
    double mobility, double dt_s, double alpha1, double alpha2,
    double pressure_correction) noexcept;
void material_terminal_set_raw(std::uint8_t point,
                               std::uint8_t mode) noexcept;
void material_terminal_reset_raw() noexcept;
std::uint64_t material_terminal_calls_raw(std::uint8_t point) noexcept;
std::uint8_t material_terminal_reach_raw(std::uint8_t point) noexcept;
void material_face_flux_observation_begin_raw() noexcept;
void material_face_flux_observation_end_raw() noexcept;
bool material_face_flux_observation_active_raw() noexcept;
runtime::HaloPerformanceCounters material_combine_pressure_halo_counters_raw(
    runtime::HaloPerformanceCounters,
    const runtime::HaloPerformanceCounters &);

} // namespace hundun::flow::detail

namespace hundun::flow::test {

struct MaterialDensityVortexSource final {
  double x{};
  double y{};
  double z{};
};

struct MaterialMomentumConservationInput final {
  double momentum_n_minus_1{};
  double momentum_n{};
  double momentum_n_plus_1{};
  double boundary_n_minus_1{};
  double boundary_n{};
  double source_n_plus_1{};
  double momentum_abs_n_minus_1{};
  double momentum_abs_n{};
  double momentum_abs_n_plus_1{};
  double boundary_abs_n_minus_1{};
  double boundary_abs_n{};
  double source_abs_n_plus_1{};
  double dt_s{};
  double alpha0{};
  double alpha2{};
  bool bdf2{};
};

struct MaterialPhaseFailureForTest final {
  bool failed{};
  std::uint8_t reason{};
  bool recoverable{};
};

struct MaterialPhaseSelectionForTest final {
  std::uint8_t reason{};
  int lowest_failing_rank{-1};
  bool recoverable{};
};

struct MaterialPressureEvidenceForTest final {
  std::vector<double> rhs_raw;
  std::vector<double> rhs_solve;
  std::vector<double> correction;
  std::vector<double> final_face_density;
  double rhs_l2{};
  std::uint32_t corrector_ordinal{};
  bool token_available{};
  bool final_operator_available{};
  std::uint64_t final_operator_revision{};
};

struct MaterialPressureHaloCountersForTest final {
  bool ordinary_available{};
  runtime::HaloPerformanceCounters ordinary;
  bool material_final_available{};
  runtime::HaloPerformanceCounters material_final;
};

struct FacadeCacheSnapshot final {
  struct Workspace final {
    std::uintptr_t identity{};
    std::size_t capacity{};
  };
  struct MomentumOperator final {
    std::uintptr_t identity{};
    std::uint64_t revision{};
    std::vector<double> diagonal;
  };

  std::vector<Workspace> workspaces;
  std::array<MomentumOperator, 3> operators{};
  std::size_t operator_count{};
  bool delegated{};
};

enum class MaterialReportCorruptionForTest : std::uint8_t {
  success_corrector_count,
  success_provenance,
  success_shared_field,
  reliable_collective_rank,
  material_reason_mapping,
  parent_transport_size,
  material_transport_size,
  unavailable_numeric_value,
  material_count_zero,
  material_count_plus_two,
  material_count_five,
  parent_transport_residual_value,
  nested_transport_residual_value,
  parent_transport_conservation_value,
  nested_transport_conservation_value,
  nested_density_residual_availability,
  nested_transport_residual_availability,
  nested_mass_conservation_availability,
  nested_transport_conservation_availability,
  nested_minimum_density_availability,
  nested_attempt_identity,
  nested_finalization_identity,
  nested_shared_field,
  nested_provenance,
  nested_residual_outer_size,
  nested_conservation_outer_size
};

enum class MaterialTerminalModeForTest : std::uint8_t {
  none,
  returned_rankless,
  thrown_operation,
  returned_reliable
};

enum class MaterialTransportAuthorityMutation : std::uint8_t {
  omitted,
  reordered,
  same_maximum_different_sequence
};

enum class MaterialTerminalPointForTest : std::uint8_t {
  predictor_stage,
  momentum_x,
  momentum_y,
  momentum_z,
  pressure_corrector_one,
  provisional_stage,
  pressure_corrector_two,
  public_finalizer,
  final_continuity_reduction,
  final_continuity_status,
  final_pressure_entry,
  final_pressure_gamma_sum,
  final_pressure_gamma_count,
  final_pressure_residual_reduction,
  final_momentum_residual_reduction,
  final_momentum_conservation_reduction,
  final_momentum_status,
  final_conservation_status,
  count
};

class MaterialDensityPisoTestAccess final {
public:
  static constexpr std::uint32_t contract_version() noexcept { return 1U; }

  static double momentum_conservation_defect(
      const MaterialMomentumConservationInput &) noexcept;
  static bool report_authenticated(
      const MaterialDensityStepAttemptReport &) noexcept;
  static void corrupt_report(MaterialDensityStepAttemptReport &,
                             MaterialReportCorruptionForTest);
  static void set_preflight_allocation_failure_rank(int) noexcept;
  static void reset_preflight_allocation_failure() noexcept;
  static void set_terminal_fault(MaterialTerminalPointForTest point,
                                 MaterialTerminalModeForTest mode) noexcept {
    detail::material_terminal_set_raw(static_cast<std::uint8_t>(point),
                                      static_cast<std::uint8_t>(mode));
  }
  static void reset_terminal_fault() noexcept {
    detail::material_terminal_reset_raw();
  }
  static std::uint64_t
  terminal_point_calls(MaterialTerminalPointForTest point) noexcept {
    return detail::material_terminal_calls_raw(
        static_cast<std::uint8_t>(point));
  }
  static MaterialPhaseSelectionForTest select_phase_failure(
      const runtime::MpiContext &, MaterialPhaseFailureForTest);
  static void require_reliable_collective_result(std::uint8_t reason,
                                                 int lowest_failing_rank);
  static MaterialPressureEvidenceForTest
  material_pressure_evidence(const FixedStepMaterialDensityFlow &);
  static MaterialPressureHaloCountersForTest
  material_pressure_halo_counters(const FixedStepMaterialDensityFlow &);
  static MaterialPressureHaloCountersForTest
  material_pressure_halo_counters(const PisoCoupler &);
  static runtime::HaloPerformanceCounters combine_pressure_halo_counters(
      runtime::HaloPerformanceCounters,
      const runtime::HaloPerformanceCounters &);
  static FacadeCacheSnapshot
  facade_cache_snapshot(const FixedStepMaterialDensityFlow &);
  static MaterialPressureEvidenceForTest
  material_pressure_evidence(const PisoCoupler &);
  static const std::vector<double> &
  finalizer_flux_evidence(const FixedStepMaterialDensityFlow &);

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
  static bool state_attempt_active(const FlowState &) noexcept;
  static std::uint64_t state_attempt_identity(const FlowState &) noexcept;
  static std::uint64_t state_allocation_identity(const FlowState &,
                                                 FlowLayer);
  static void set_accepted_face_mass_flux(FlowState &, std::size_t face,
                                          double value);
  static void force_state_metadata(FlowState &,
                                   AcceptedStepMetadata) noexcept;
  static void prepare_state_commit(FlowState &, AcceptedStepMetadata);
  static void publish_state_commit(FlowState &) noexcept;
  static bool face_flux_path_observation_active() noexcept;
  static void force_flow_attempt_identity(FixedStepMaterialDensityFlow &,
                                          std::uint64_t) noexcept;
  static bool has_diagnostic_report(
      const FixedStepMaterialDensityFlow &) noexcept;
  static void enable_vortex_source(FixedStepMaterialDensityFlow &,
                                   bool) noexcept;
  static void mutate_transport_authority(
      FixedStepMaterialDensityFlow &,
      MaterialTransportAuthorityMutation) noexcept;
  static double material_face_value(
      double predictor_velocity_interpolation,
      double predictor_momentum_interpolation, double face_density,
      double cell_momentum_n, double face_momentum_n,
      double cell_momentum_n_minus_1, double face_momentum_n_minus_1,
      double mobility, double dt_s, double alpha1, double alpha2,
      double pressure_correction) noexcept {
    return detail::material_face_value_raw(
        predictor_velocity_interpolation, predictor_momentum_interpolation,
        face_density, cell_momentum_n, face_momentum_n,
        cell_momentum_n_minus_1, face_momentum_n_minus_1, mobility, dt_s,
        alpha1, alpha2, pressure_correction);
  }
};

inline double MaterialDensityPisoTestAccess::momentum_conservation_defect(
    const MaterialMomentumConservationInput &input) noexcept {
  return detail::material_momentum_conservation_defect_raw(
      input.momentum_n_minus_1, input.momentum_n, input.momentum_n_plus_1,
      input.boundary_n_minus_1, input.boundary_n, input.source_n_plus_1,
      input.momentum_abs_n_minus_1, input.momentum_abs_n,
      input.momentum_abs_n_plus_1, input.boundary_abs_n_minus_1,
      input.boundary_abs_n, input.source_abs_n_plus_1, input.dt_s,
      input.alpha0, input.alpha2, input.bdf2);
}

inline bool MaterialDensityPisoTestAccess::report_authenticated(
    const MaterialDensityStepAttemptReport &report) noexcept {
  return detail::DensityClosureBridge::material_report_authenticated_raw(
      report);
}

inline void MaterialDensityPisoTestAccess::corrupt_report(
    MaterialDensityStepAttemptReport &report,
    MaterialReportCorruptionForTest change) {
  detail::DensityClosureBridge::corrupt_material_report_raw(
      report, static_cast<std::uint8_t>(change));
}

inline void MaterialDensityPisoTestAccess::
    set_preflight_allocation_failure_rank(int rank) noexcept {
  detail::material_set_preflight_allocation_failure_rank_raw(rank);
}

inline void MaterialDensityPisoTestAccess::
    reset_preflight_allocation_failure() noexcept {
  detail::material_reset_preflight_allocation_failure_raw();
}

inline bool MaterialDensityPisoTestAccess::
    face_flux_path_observation_active() noexcept {
  return detail::material_face_flux_observation_active_raw();
}

inline MaterialPhaseSelectionForTest
MaterialDensityPisoTestAccess::select_phase_failure(
    const runtime::MpiContext &mpi, MaterialPhaseFailureForTest input) {
  MaterialPhaseSelectionForTest result;
  detail::material_select_phase_failure_raw(
      mpi, input.failed, input.reason, input.recoverable, result.reason,
      result.lowest_failing_rank, result.recoverable);
  return result;
}

inline void MaterialDensityPisoTestAccess::
    require_reliable_collective_result(std::uint8_t reason,
                                       int lowest_failing_rank) {
  detail::material_require_reliable_collective_result_raw(
      reason, lowest_failing_rank);
}

inline MaterialPressureEvidenceForTest
MaterialDensityPisoTestAccess::material_pressure_evidence(
    const FixedStepMaterialDensityFlow &flow) {
  MaterialPressureEvidenceForTest result;
  detail::AdaptiveTimeControlAccess::material_flow_pressure_values_raw(
      flow, result.rhs_raw, result.rhs_solve, result.correction,
      result.final_face_density, result.rhs_l2, result.corrector_ordinal,
      result.token_available, result.final_operator_available,
      result.final_operator_revision);
  return result;
}

inline MaterialPressureHaloCountersForTest
MaterialDensityPisoTestAccess::material_pressure_halo_counters(
    const FixedStepMaterialDensityFlow &flow) {
  MaterialPressureHaloCountersForTest result;
  detail::AdaptiveTimeControlAccess::material_flow_pressure_halo_values_raw(
      flow, result.ordinary_available, result.ordinary,
      result.material_final_available, result.material_final);
  return result;
}

inline FacadeCacheSnapshot MaterialDensityPisoTestAccess::facade_cache_snapshot(
    const FixedStepMaterialDensityFlow &flow) {
  FacadeCacheSnapshot result;
  std::vector<std::uintptr_t> identities;
  std::vector<std::size_t> capacities;
  std::array<std::uintptr_t, 3> operator_ids{};
  std::array<std::uint64_t, 3> operator_revisions{};
  std::array<std::vector<double>, 3> diagonals{};
  detail::AdaptiveTimeControlAccess::material_cache_values_raw(
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

inline const std::vector<double> &
MaterialDensityPisoTestAccess::finalizer_flux_evidence(
    const FixedStepMaterialDensityFlow &flow) {
  return detail::AdaptiveTimeControlAccess::material_finalizer_flux_values_raw(
      flow);
}

inline void MaterialDensityPisoTestAccess::force_flow_attempt_identity(
    FixedStepMaterialDensityFlow &flow, std::uint64_t value) noexcept {
  detail::AdaptiveTimeControlAccess::material_force_attempt_identity_raw(flow,
                                                                         value);
}

inline bool MaterialDensityPisoTestAccess::has_diagnostic_report(
    const FixedStepMaterialDensityFlow &flow) noexcept {
  return detail::AdaptiveTimeControlAccess::material_has_diagnostic_report_raw(
      flow);
}

inline void MaterialDensityPisoTestAccess::enable_vortex_source(
    FixedStepMaterialDensityFlow &flow, bool enabled) noexcept {
  detail::AdaptiveTimeControlAccess::material_enable_vortex_source_raw(
      flow, enabled);
}

inline void MaterialDensityPisoTestAccess::mutate_transport_authority(
    FixedStepMaterialDensityFlow &flow,
    MaterialTransportAuthorityMutation change) noexcept {
  detail::AdaptiveTimeControlAccess::material_change_transport_authority_raw(
      flow, static_cast<std::uint8_t>(change));
}

inline MaterialPressureEvidenceForTest
MaterialDensityPisoTestAccess::material_pressure_evidence(
    const PisoCoupler &coupler) {
  MaterialPressureEvidenceForTest result;
  detail::AdaptiveTimeControlAccess::material_pressure_values_raw(
      coupler, result.rhs_raw, result.rhs_solve, result.correction,
      result.final_face_density, result.rhs_l2, result.corrector_ordinal,
      result.token_available, result.final_operator_available,
      result.final_operator_revision);
  return result;
}

inline MaterialPressureHaloCountersForTest
MaterialDensityPisoTestAccess::material_pressure_halo_counters(
    const PisoCoupler &coupler) {
  MaterialPressureHaloCountersForTest result;
  detail::AdaptiveTimeControlAccess::material_pressure_halo_values_raw(
      coupler, result.ordinary_available, result.ordinary,
      result.material_final_available, result.material_final);
  return result;
}

inline runtime::HaloPerformanceCounters
MaterialDensityPisoTestAccess::combine_pressure_halo_counters(
    runtime::HaloPerformanceCounters left,
    const runtime::HaloPerformanceCounters &right) {
  return detail::material_combine_pressure_halo_counters_raw(left, right);
}

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
inline void MaterialDensityPisoTestAccess::force_state_diagnostic_identity(
    FlowState &state, std::uint64_t value) noexcept {
  state.impl_->diagnostic_mutation_identity = value;
}

inline std::uint64_t
MaterialDensityPisoTestAccess::state_diagnostic_identity(
    const FlowState &state) noexcept {
  return state.diagnostic_mutation_identity();
}

inline bool MaterialDensityPisoTestAccess::state_attempt_active(
    const FlowState &state) noexcept {
  return state.attempt_active();
}

inline std::uint64_t MaterialDensityPisoTestAccess::state_attempt_identity(
    const FlowState &state) noexcept {
  return state.attempt_identity();
}

inline std::uint64_t
MaterialDensityPisoTestAccess::state_allocation_identity(
    const FlowState &state, FlowLayer selected) {
  const runtime::FieldStorage *storage{};
  switch (selected) {
  case FlowLayer::history:
    storage = &state.impl_->history;
    break;
  case FlowLayer::committed:
    storage = &state.impl_->committed;
    break;
  case FlowLayer::trial:
    storage = &state.impl_->trial;
    break;
  }
  if (storage == nullptr)
    throw runtime::Error("flow-state test layer is invalid");
  const auto view = storage->acquire_read<double>(
      state.impl_->access, runtime::PhaseId{1800U}, runtime::ActorId{1800U},
      state.impl_->fields.density);
  return static_cast<std::uint64_t>(
      reinterpret_cast<std::uintptr_t>(&view(0, 0, 0, 0)));
}

inline void MaterialDensityPisoTestAccess::set_accepted_face_mass_flux(
    FlowState &state, std::size_t face, double value) {
  auto view = state.impl_->committed.acquire_face_write<double>(
      state.impl_->access, runtime::PhaseId{1800U}, runtime::ActorId{1800U},
      state.impl_->fields.face_mass_flux);
  if (face >= view.face_count())
    throw runtime::Error("flow-state test face is out of range");
  view(face, 0) = value;
}

inline void MaterialDensityPisoTestAccess::force_state_metadata(
    FlowState &state, AcceptedStepMetadata metadata) noexcept {
  state.impl_->metadata = metadata;
}

inline void MaterialDensityPisoTestAccess::prepare_state_commit(
    FlowState &state, AcceptedStepMetadata accepted) {
  state.prepare_commit_attempt(accepted);
}

inline void MaterialDensityPisoTestAccess::publish_state_commit(
    FlowState &state) noexcept {
  state.publish_commit_attempt();
}
#endif

} // namespace hundun::flow::test

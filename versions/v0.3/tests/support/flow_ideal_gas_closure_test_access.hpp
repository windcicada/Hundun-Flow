// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#error "ideal-gas closure test access is unavailable in tests-off builds"
#endif

#include "flow_adaptive_time_control_detail.hpp"
#include "flow_density_closure_detail.hpp"
#include "hundun/flow_ideal_gas_piso.hpp"
#include "hundun/rt_error.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace hundun::flow::test {

enum class IdealGasPostEvidenceMutation : std::uint8_t {
  post_report_begin = 0,
  post_authority_begin = 29,
  post_coherent_begin = 58,
  pre_report_begin = 87,
  pre_authority_begin = 116,
  pre_coherent_begin = 145,
  outer_begin = 174,
  count = 194
};

enum class IdealGasPrepareFault : std::uint8_t {
  state_prepare,
  closure_prepare
};

enum class IdealGasPostAssessmentFault : std::uint8_t {
  non_finite_state,
  non_positive_density,
  final_transport_residual,
  final_conservation_defect
};

enum class IdealGasOuterFailurePoint : std::uint8_t {
  momentum_after_predictor,
  pressure_after_first_corrector
};

enum class IdealGasCreateFault : std::uint8_t {
  mode_disagreement,
  ownership_gap,
  ownership_overlap,
  ownership_swap,
  local_preparation,
  preflight_workspace_allocation,
  construction_allocation,
  inlet_cp,
  inlet_gas_constant,
  inlet_pressure,
  sum_reduction,
  maximum_reduction
};

enum class IdealGasAttemptPreparationFault : std::uint8_t {
  predictor_local,
  predictor_write_capability,
  provisional_local,
  provisional_write_capability,
  final_local,
  final_write_capability,
  final_readback,
  post_density_views,
  post_geometry,
  post_transport_views,
  post_next_integrals,
  post_report_finalization,
  pre_authority_preparation,
  post_authority_preparation
};

enum class IdealGasMetricGateFault : std::uint8_t {
  eos,
  rho_remap,
  rho_h_remap,
  mass,
  enthalpy
};

enum class IdealGasStepReportCorruption : std::uint8_t {
  nested_material,
  closure_stage,
  closure_candidate_pressure,
  closure_seal,
  outer_attempt_identity,
  outer_closure_presence,
  outer_seal
};

struct IdealGasHaloTraceEntry final {
  IdealGasClosureStage stage{};
  runtime::FieldId density{};
  runtime::FieldId enthalpy_density{};
};

struct IdealGasFacadeCacheSnapshot final {
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

class IdealGasClosureTestAccess final {
public:
  static void set_uniform_enthalpy_rate(FixedStepIdealGasFlow &,
                                        double rate_J_per_kg_s);
  static bool report_authenticated(const IdealGasClosureReport &) noexcept;
  static bool candidate_pressure_mutation_rejected(
      const IdealGasClosureReport &, bool availability) noexcept;
  static bool report_authenticated(const IdealGasStepAttemptReport &) noexcept;
  static bool post_eos_evidence_authenticated(
      const MaterialDensityStepAttemptReport &) noexcept;
  static bool same_rank_reason_precedence_is_enum_order() noexcept;
  static bool post_store_rank_marker_is_collision_free(int ranks) noexcept;
  static bool final_gate_rank_marker_is_collision_free(int ranks) noexcept;
  static IdealGasClosure
  create(const mesh::MeshTopology &, const mesh::MeshGeometry &,
         const boundary::BoundaryRegistry &, const runtime::MpiContext &,
         const runtime::FieldRegistry &, const FlowFieldIds &,
         const FlowState &, IdealGasClosureSpec, IdealGasCreateFault, int rank);
  static IdealGasClosure
  create_immersed(const mesh::MeshTopology &, const mesh::MeshGeometry &,
                  const boundary::BoundaryRegistry &,
                  const immersed::ImmersedDomain &,
                  const runtime::MpiContext &,
                  const runtime::FieldRegistry &, const FlowFieldIds &,
                  const FlowState &, IdealGasClosureSpec);
  static int preflight_failure_rank(const runtime::Error &) noexcept;
  static std::uint64_t
  preflight_wire_exchange_count(const IdealGasClosure &) noexcept;
  static IdealGasClosureFailureReason
  create_validation_failure_reason(const runtime::Error &) noexcept;
  static void set_restore_preparation_fault(int rank) noexcept;
  static void set_restore_snapshot_shape_fault(int rank) noexcept;
  static void set_facade_create_fault(IdealGasClosure &, int rank);
  static void set_material_factory_create_fault(IdealGasClosure &, int rank);
  static bool consume_facade_create_fault(IdealGasClosure &, int rank) noexcept;
  static int
  consume_material_factory_create_fault(IdealGasClosure &) noexcept;
  static void begin_attempt(IdealGasClosure &, FlowState &,
                            std::uint64_t identity);
  static IdealGasClosureReport evaluate(IdealGasClosure &, FlowState &,
                                        IdealGasClosureStage);
  static int prepare_commit(IdealGasClosure &);
  static void publish_commit(IdealGasClosure &) noexcept;
  static void rollback(IdealGasClosure &) noexcept;
  static void set_stage_failure(IdealGasClosure &, IdealGasClosureStage,
                                IdealGasClosureFailureReason, int rank);
  static void set_metric_gate_failure(IdealGasClosure &,
                                      IdealGasMetricGateFault, int rank);
  static void set_post_store_corruption(IdealGasClosure &, int rank,
                                        bool enthalpy_density);
  static void exhaust_source_generation(FixedStepIdealGasFlow &);
  static void force_finalization_identity_wrap(FixedStepIdealGasFlow &);
  static std::vector<IdealGasHaloTraceEntry>
  halo_trace(const FixedStepIdealGasFlow &);
  static IdealGasFacadeCacheSnapshot
  facade_cache_snapshot(const FixedStepIdealGasFlow &);
  static IdealGasFacadeCacheSnapshot
  delegated_material_cache_snapshot(const FixedStepIdealGasFlow &);
  static std::uint64_t
  source_generation(const IdealGasClosureDiagnosticSource &);
  static void set_post_store_corruption(FixedStepIdealGasFlow &, int rank,
                                        bool enthalpy_density);
  static void set_candidate_precedence_fault(FixedStepIdealGasFlow &, int rank);
  static void set_stage_failure(FixedStepIdealGasFlow &, IdealGasClosureStage,
                                IdealGasClosureFailureReason, int rank);
  static void set_outer_failure(FixedStepIdealGasFlow &,
                                IdealGasOuterFailurePoint, int rank);
  static void set_prepare_fault(FixedStepIdealGasFlow &, IdealGasPrepareFault,
                                int rank);
  static void set_post_store_mpi_fault(FixedStepIdealGasFlow &, int rank);
  static void set_post_assessment_fault(FixedStepIdealGasFlow &,
                                        IdealGasPostAssessmentFault, int rank);
  static void set_attempt_layout_fault(FixedStepIdealGasFlow &, int rank);
  static void set_attempt_layout_fault(IdealGasClosure &, int rank);
  static void set_outlet_backflow_fault(FixedStepIdealGasFlow &);
  static void set_attempt_preparation_fault(
      FixedStepIdealGasFlow &, IdealGasAttemptPreparationFault, int rank);
  static void set_attempt_preparation_fault(
      IdealGasClosure &, IdealGasAttemptPreparationFault, int rank);
  static void set_controlled_allocation(FixedStepIdealGasFlow &, int rank);
  static void set_controlled_allocation(IdealGasClosure &, int rank);
  static bool
  allocation_observation_active(const FixedStepIdealGasFlow &) noexcept;
  static bool
  allocation_observation_active(const IdealGasClosure &) noexcept;
  static void begin_allocation_observation(IdealGasClosure &) noexcept;
  static void end_allocation_observation(IdealGasClosure &) noexcept;
  static void consume_attempt_preparation_fault(IdealGasClosure &);
  static bool post_evidence_mutation_rejected(
      const IdealGasStepAttemptReport &, IdealGasPostEvidenceMutation);
  static bool closure_field_mutation_rejected(
      const IdealGasStepAttemptReport &, std::uint8_t field);
  static void corrupt_report(IdealGasStepAttemptReport &,
                             IdealGasStepReportCorruption);
};

inline bool IdealGasClosureTestAccess::report_authenticated(
    const IdealGasClosureReport &report) noexcept {
  return detail::DensityClosureAdapter::closure_report_authenticated_raw(
      report);
}
inline bool IdealGasClosureTestAccess::candidate_pressure_mutation_rejected(
    const IdealGasClosureReport &report, bool availability) noexcept {
  return detail::DensityClosureAdapter::candidate_pressure_change_rejected_raw(
      report, availability);
}
inline bool IdealGasClosureTestAccess::
    same_rank_reason_precedence_is_enum_order() noexcept {
  return detail::DensityClosureAdapter::same_rank_reason_precedence_raw();
}
inline bool IdealGasClosureTestAccess::
    post_store_rank_marker_is_collision_free(int ranks) noexcept {
  return detail::DensityClosureAdapter::
      post_store_rank_marker_collision_free_raw(ranks);
}
inline bool IdealGasClosureTestAccess::
    final_gate_rank_marker_is_collision_free(int ranks) noexcept {
  return detail::DensityClosureAdapter::
      final_gate_rank_marker_collision_free_raw(ranks);
}
inline IdealGasClosure IdealGasClosureTestAccess::create(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const runtime::MpiContext &mpi, const runtime::FieldRegistry &registry,
    const FlowFieldIds &fields, const FlowState &state,
    IdealGasClosureSpec spec, IdealGasCreateFault fault, int rank) {
  return detail::DensityClosureAdapter::create_raw(
      topology, geometry, boundaries, mpi, registry, fields, state, spec,
      static_cast<std::uint8_t>(fault), rank);
}
inline IdealGasClosure IdealGasClosureTestAccess::create_immersed(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const immersed::ImmersedDomain &domain, const runtime::MpiContext &mpi,
    const runtime::FieldRegistry &registry, const FlowFieldIds &fields,
    const FlowState &state, IdealGasClosureSpec spec) {
  return detail::DensityClosureAdapter::create_immersed(
      topology, geometry, boundaries, domain, mpi, registry, fields, state,
      spec);
}
inline int IdealGasClosureTestAccess::preflight_failure_rank(
    const runtime::Error &error) noexcept {
  return detail::DensityClosureAdapter::preflight_failure_rank_raw(error);
}
inline std::uint64_t IdealGasClosureTestAccess::preflight_wire_exchange_count(
    const IdealGasClosure &closure) noexcept {
  return detail::DensityClosureAdapter::preflight_wire_exchange_count_raw(
      closure);
}
inline IdealGasClosureFailureReason
IdealGasClosureTestAccess::create_validation_failure_reason(
    const runtime::Error &error) noexcept {
  return detail::DensityClosureAdapter::create_validation_failure_reason_raw(
      error);
}
inline void IdealGasClosureTestAccess::set_restore_preparation_fault(
    int rank) noexcept {
  detail::DensityClosureAdapter::set_restore_preparation_fault_raw(rank);
}
inline void IdealGasClosureTestAccess::set_restore_snapshot_shape_fault(
    int rank) noexcept {
  detail::DensityClosureAdapter::set_restore_shape_fault_raw(rank);
}
inline void IdealGasClosureTestAccess::set_facade_create_fault(
    IdealGasClosure &closure, int rank) {
  detail::DensityClosureAdapter::set_facade_create_fault_raw(closure, rank);
}
inline void IdealGasClosureTestAccess::set_material_factory_create_fault(
    IdealGasClosure &closure, int rank) {
  detail::DensityClosureAdapter::set_material_factory_create_fault_raw(
      closure, rank);
}
inline bool IdealGasClosureTestAccess::consume_facade_create_fault(
    IdealGasClosure &closure, int rank) noexcept {
  return detail::DensityClosureAdapter::consume_facade_create_fault_raw(
      closure, rank);
}
inline int IdealGasClosureTestAccess::consume_material_factory_create_fault(
    IdealGasClosure &closure) noexcept {
  return detail::DensityClosureAdapter::
      consume_material_factory_create_fault_raw(closure);
}
inline void IdealGasClosureTestAccess::begin_attempt(
    IdealGasClosure &closure, FlowState &state, std::uint64_t identity) {
  detail::DensityClosureAdapter::begin(&closure, state, identity);
}
inline IdealGasClosureReport IdealGasClosureTestAccess::evaluate(
    IdealGasClosure &closure, FlowState &state, IdealGasClosureStage stage) {
  return detail::DensityClosureAdapter::evaluate_raw(closure, state, stage);
}
inline int IdealGasClosureTestAccess::prepare_commit(
    IdealGasClosure &closure) {
  return closure.prepare_commit();
}
inline void IdealGasClosureTestAccess::publish_commit(
    IdealGasClosure &closure) noexcept {
  closure.publish_commit();
}
inline void IdealGasClosureTestAccess::rollback(
    IdealGasClosure &closure) noexcept {
  detail::DensityClosureAdapter::rollback(&closure);
}
inline void IdealGasClosureTestAccess::set_stage_failure(
    IdealGasClosure &closure, IdealGasClosureStage stage,
    IdealGasClosureFailureReason reason, int rank) {
  detail::DensityClosureAdapter::set_stage_failure_raw(closure, stage, reason,
                                                       rank);
}
inline void IdealGasClosureTestAccess::set_metric_gate_failure(
    IdealGasClosure &closure, IdealGasMetricGateFault fault, int rank) {
  detail::DensityClosureAdapter::set_metric_gate_failure_raw(
      closure, static_cast<std::uint8_t>(fault), rank);
}
inline void IdealGasClosureTestAccess::set_post_store_corruption(
    IdealGasClosure &closure, int rank, bool enthalpy_density) {
  detail::DensityClosureAdapter::set_post_store_corruption_raw(
      closure, rank, enthalpy_density);
}
inline void IdealGasClosureTestAccess::set_attempt_preparation_fault(
    IdealGasClosure &closure, IdealGasAttemptPreparationFault fault,
    int rank) {
  detail::DensityClosureAdapter::set_attempt_preparation_fault_raw(
      closure, static_cast<std::uint8_t>(fault), rank);
}
inline void IdealGasClosureTestAccess::set_controlled_allocation(
    IdealGasClosure &closure, int rank) {
  detail::DensityClosureAdapter::set_controlled_allocation_raw(closure, rank);
}
inline bool IdealGasClosureTestAccess::allocation_observation_active(
    const IdealGasClosure &closure) noexcept {
  return detail::DensityClosureAdapter::allocation_observation_active_raw(
      closure);
}
inline void IdealGasClosureTestAccess::begin_allocation_observation(
    IdealGasClosure &closure) noexcept {
  detail::DensityClosureAdapter::begin_allocation_observation(&closure);
}
inline void IdealGasClosureTestAccess::end_allocation_observation(
    IdealGasClosure &closure) noexcept {
  detail::DensityClosureAdapter::end_allocation_observation(&closure);
}
inline void IdealGasClosureTestAccess::consume_attempt_preparation_fault(
    IdealGasClosure &closure) {
  detail::DensityClosureAdapter::consume_attempt_preparation_fault_raw(
      closure);
}

inline void IdealGasClosureTestAccess::set_uniform_enthalpy_rate(
    FixedStepIdealGasFlow &flow, double rate_J_per_kg_s) {
  detail::AdaptiveTimeControlAccess::ideal_set_enthalpy_rate_raw(
      flow, rate_J_per_kg_s);
}

inline bool IdealGasClosureTestAccess::report_authenticated(
    const IdealGasStepAttemptReport &report) noexcept {
  return report.authenticated();
}

inline bool IdealGasClosureTestAccess::post_eos_evidence_authenticated(
    const MaterialDensityStepAttemptReport &report) noexcept {
  return detail::DensityClosureBridge::post_eos_evidence_authenticated(report);
}

inline void IdealGasClosureTestAccess::exhaust_source_generation(
    FixedStepIdealGasFlow &flow) {
  detail::AdaptiveTimeControlAccess::ideal_exhaust_source_generation_raw(flow);
}

inline void IdealGasClosureTestAccess::force_finalization_identity_wrap(
    FixedStepIdealGasFlow &flow) {
  detail::AdaptiveTimeControlAccess::
      ideal_force_finalization_identity_wrap_raw(flow);
}

inline std::vector<IdealGasHaloTraceEntry>
IdealGasClosureTestAccess::halo_trace(const FixedStepIdealGasFlow &flow) {
  std::vector<IdealGasHaloTraceEntry> result;
  for (const auto &entry :
       detail::AdaptiveTimeControlAccess::ideal_halo_trace_raw(flow))
    result.push_back({static_cast<IdealGasClosureStage>(entry[0]),
                      static_cast<runtime::FieldId>(entry[1]),
                      static_cast<runtime::FieldId>(entry[2])});
  return result;
}

inline IdealGasFacadeCacheSnapshot ideal_cache_values(
    const FixedStepIdealGasFlow &flow, bool delegated) {
  IdealGasFacadeCacheSnapshot result;
  std::vector<std::uintptr_t> identities;
  std::vector<std::size_t> capacities;
  std::array<std::uintptr_t, 3> operator_ids{};
  std::array<std::uint64_t, 3> operator_revisions{};
  std::array<std::vector<double>, 3> diagonals{};
  detail::AdaptiveTimeControlAccess::ideal_cache_values_raw(
      flow, delegated, identities, capacities, operator_ids,
      operator_revisions, diagonals, result.operator_count,
      result.delegated);
  result.workspaces.reserve(identities.size());
  for (std::size_t index = 0; index < identities.size(); ++index)
    result.workspaces.push_back({identities[index], capacities[index]});
  for (std::size_t index = 0; index < result.operator_count; ++index)
    result.operators[index] = {operator_ids[index], operator_revisions[index],
                               std::move(diagonals[index])};
  return result;
}

inline IdealGasFacadeCacheSnapshot
IdealGasClosureTestAccess::facade_cache_snapshot(
    const FixedStepIdealGasFlow &flow) {
  return ideal_cache_values(flow, true);
}

inline IdealGasFacadeCacheSnapshot
IdealGasClosureTestAccess::delegated_material_cache_snapshot(
    const FixedStepIdealGasFlow &flow) {
  return ideal_cache_values(flow, false);
}

inline std::uint64_t IdealGasClosureTestAccess::source_generation(
    const IdealGasClosureDiagnosticSource &source) {
  return detail::DensityClosureDiagnosticAccess::source_generation_raw(source);
}

inline void IdealGasClosureTestAccess::set_post_store_corruption(
    FixedStepIdealGasFlow &flow, int rank, bool enthalpy_density) {
  detail::AdaptiveTimeControlAccess::ideal_set_post_store_corruption_raw(
      flow, rank, enthalpy_density);
}

inline void IdealGasClosureTestAccess::set_candidate_precedence_fault(
    FixedStepIdealGasFlow &flow, int rank) {
  detail::AdaptiveTimeControlAccess::
      ideal_set_candidate_precedence_fault_raw(flow, rank);
}

inline void IdealGasClosureTestAccess::set_stage_failure(
    FixedStepIdealGasFlow &flow, IdealGasClosureStage stage,
    IdealGasClosureFailureReason reason, int rank) {
  detail::AdaptiveTimeControlAccess::ideal_set_stage_failure_raw(
      flow, stage, reason, rank);
}

inline void IdealGasClosureTestAccess::set_outer_failure(
    FixedStepIdealGasFlow &flow, IdealGasOuterFailurePoint point, int rank) {
  detail::AdaptiveTimeControlAccess::ideal_set_outer_failure_raw(
      flow, static_cast<std::uint8_t>(point), rank);
}

inline void IdealGasClosureTestAccess::set_prepare_fault(
    FixedStepIdealGasFlow &flow, IdealGasPrepareFault fault, int rank) {
  detail::AdaptiveTimeControlAccess::ideal_set_prepare_fault_raw(
      flow, static_cast<std::uint8_t>(fault), rank);
}

inline void IdealGasClosureTestAccess::set_post_store_mpi_fault(
    FixedStepIdealGasFlow &flow, int rank) {
  detail::AdaptiveTimeControlAccess::ideal_set_post_store_mpi_fault_raw(flow,
                                                                        rank);
}

inline void IdealGasClosureTestAccess::set_post_assessment_fault(
    FixedStepIdealGasFlow &flow, IdealGasPostAssessmentFault fault, int rank) {
  detail::AdaptiveTimeControlAccess::ideal_set_post_assessment_fault_raw(
      flow, static_cast<std::uint8_t>(fault), rank);
}

inline void IdealGasClosureTestAccess::set_attempt_layout_fault(
    FixedStepIdealGasFlow &flow, int rank) {
  detail::AdaptiveTimeControlAccess::ideal_set_attempt_layout_fault_raw(flow,
                                                                        rank);
}

inline void IdealGasClosureTestAccess::set_attempt_layout_fault(
    IdealGasClosure &closure, int rank) {
  detail::DensityClosureAdapter::set_attempt_layout_fault_raw(closure, rank);
}

inline void IdealGasClosureTestAccess::set_outlet_backflow_fault(
    FixedStepIdealGasFlow &flow) {
  detail::AdaptiveTimeControlAccess::ideal_set_outlet_backflow_fault_raw(flow);
}

inline void IdealGasClosureTestAccess::set_attempt_preparation_fault(
    FixedStepIdealGasFlow &flow, IdealGasAttemptPreparationFault fault,
    int rank) {
  detail::AdaptiveTimeControlAccess::ideal_set_attempt_preparation_fault_raw(
      flow, static_cast<std::uint8_t>(fault), rank);
}

inline void IdealGasClosureTestAccess::set_controlled_allocation(
    FixedStepIdealGasFlow &flow, int rank) {
  detail::AdaptiveTimeControlAccess::ideal_set_controlled_allocation_raw(flow,
                                                                         rank);
}

inline bool IdealGasClosureTestAccess::allocation_observation_active(
    const FixedStepIdealGasFlow &flow) noexcept {
  return detail::AdaptiveTimeControlAccess::
      ideal_allocation_observation_active_raw(flow);
}

inline bool IdealGasClosureTestAccess::post_evidence_mutation_rejected(
    const IdealGasStepAttemptReport &report,
    IdealGasPostEvidenceMutation mutation) {
  return detail::DensityClosureBridge::post_evidence_mutation_rejected(
      report.flow_, static_cast<std::uint8_t>(mutation));
}

inline bool IdealGasClosureTestAccess::closure_field_mutation_rejected(
    const IdealGasStepAttemptReport &report, std::uint8_t field) {
  if (!report.closure_report_ || field >= 28U)
    return false;
  auto copy = report;
  auto &closure = *copy.closure_report_;
  const auto next = [](double value) noexcept {
    return std::nextafter(value, std::numeric_limits<double>::infinity());
  };
  switch (field) {
  case 0U:
    closure.disposition_ =
        closure.disposition_ == IdealGasClosureDisposition::closed
            ? IdealGasClosureDisposition::recoverable_failure
            : IdealGasClosureDisposition::closed;
    break;
  case 1U:
    closure.reason_ = closure.reason_ == IdealGasClosureFailureReason::none
                          ? IdealGasClosureFailureReason::invalid_input
                          : IdealGasClosureFailureReason::none;
    break;
  case 2U:
    closure.stage_ = closure.stage_ == IdealGasClosureStage::final
                         ? IdealGasClosureStage::predictor
                         : IdealGasClosureStage::final;
    break;
  case 3U:
    ++closure.lowest_failing_rank_;
    break;
  case 4U:
    ++closure.attempt_identity_;
    break;
  case 5U:
    ++closure.evaluation_count_;
    break;
  case 6U:
    ++closure.collective_count_;
    break;
  case 7U:
    closure.pressure_mode_ =
        closure.pressure_mode_ == IdealGasPressureMode::closed_dynamic
            ? IdealGasPressureMode::open_fixed
            : IdealGasPressureMode::closed_dynamic;
    break;
  case 8U:
    closure.configured_pressure_pa_ = next(closure.configured_pressure_pa_);
    break;
  case 9U:
    closure.candidate_pressure_available_ =
        !closure.candidate_pressure_available_;
    break;
  case 10U:
    closure.candidate_pressure_pa_ = next(closure.candidate_pressure_pa_);
    break;
  case 11U:
    if (closure.target_mass_kg_)
      closure.target_mass_kg_.reset();
    else
      closure.target_mass_kg_ = 1.0;
    break;
  case 12U:
    closure.target_mass_kg_ = next(closure.target_mass_kg_.value_or(1.0));
    break;
  case 13U:
    closure.final_metrics_available_ = !closure.final_metrics_available_;
    break;
  case 14U:
    closure.actual_mass_kg_ = next(closure.actual_mass_kg_);
    break;
  case 15U:
    closure.temperature_min_K_ = next(closure.temperature_min_K_);
    break;
  case 16U:
    closure.temperature_max_K_ = next(closure.temperature_max_K_);
    break;
  case 17U:
    closure.enthalpy_min_J_per_kg_ = next(closure.enthalpy_min_J_per_kg_);
    break;
  case 18U:
    closure.enthalpy_max_J_per_kg_ = next(closure.enthalpy_max_J_per_kg_);
    break;
  case 19U:
    closure.density_min_kg_per_m3_ = next(closure.density_min_kg_per_m3_);
    break;
  case 20U:
    closure.density_max_kg_per_m3_ = next(closure.density_max_kg_per_m3_);
    break;
  case 21U:
    closure.rho_remap_normalized_l2_ =
        next(closure.rho_remap_normalized_l2_);
    break;
  case 22U:
    closure.rho_h_remap_normalized_l2_ =
        next(closure.rho_h_remap_normalized_l2_);
    break;
  case 23U:
    closure.rho_remap_relative_conservation_defect_ =
        next(closure.rho_remap_relative_conservation_defect_);
    break;
  case 24U:
    closure.rho_h_remap_relative_conservation_defect_ =
        next(closure.rho_h_remap_relative_conservation_defect_);
    break;
  case 25U:
    closure.enthalpy_temperature_max_relative_error_ =
        next(closure.enthalpy_temperature_max_relative_error_);
    break;
  case 26U:
    closure.eos_max_relative_error_ = next(closure.eos_max_relative_error_);
    break;
  case 27U:
    closure.seal_ ^= 1U;
    break;
  }
  return !copy.authenticated();
}

inline void IdealGasClosureTestAccess::corrupt_report(
    IdealGasStepAttemptReport &report, IdealGasStepReportCorruption change) {
  switch (change) {
  case IdealGasStepReportCorruption::nested_material:
    detail::DensityClosureBridge::corrupt_material_report_raw(report.flow_,
                                                               0U);
    break;
  case IdealGasStepReportCorruption::closure_stage:
    report.closure_report_->stage_ = IdealGasClosureStage::predictor;
    break;
  case IdealGasStepReportCorruption::closure_candidate_pressure:
    report.closure_report_->candidate_pressure_pa_ = std::nextafter(
        report.closure_report_->candidate_pressure_pa_,
        std::numeric_limits<double>::infinity());
    break;
  case IdealGasStepReportCorruption::closure_seal:
    report.closure_report_->seal_ ^= 1U;
    break;
  case IdealGasStepReportCorruption::outer_attempt_identity:
    ++report.attempt_identity_;
    break;
  case IdealGasStepReportCorruption::outer_closure_presence:
    report.closure_report_.reset();
    break;
  case IdealGasStepReportCorruption::outer_seal:
    report.seal_ ^= 1U;
    break;
  }
}

} // namespace hundun::flow::test

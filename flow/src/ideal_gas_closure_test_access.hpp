// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#error "ideal-gas closure test access is unavailable in tests-off builds"
#endif

#include "hundun/flow/ideal_gas_piso.hpp"
#include "hundun/runtime/error.hpp"

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
  post_authority_preparation,
  pre_authority_publication,
  post_authority_publication
};

enum class IdealGasMetricGateFault : std::uint8_t {
  eos,
  rho_remap,
  rho_h_remap,
  mass,
  enthalpy
};

struct IdealGasHaloTraceEntry final {
  IdealGasClosureStage stage{};
  runtime::FieldId density{};
  runtime::FieldId enthalpy_density{};
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
  static int preflight_failure_rank(const runtime::Error &) noexcept;
  static IdealGasClosureFailureReason
  create_validation_failure_reason(const runtime::Error &) noexcept;
  static void set_facade_create_fault(IdealGasClosure &, int rank);
  static bool consume_facade_create_fault(IdealGasClosure &, int rank) noexcept;
  static void begin_attempt(IdealGasClosure &, FlowState &,
                            std::uint64_t identity);
  static IdealGasClosureReport evaluate(IdealGasClosure &, FlowState &,
                                        IdealGasClosureStage);
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
};

} // namespace hundun::flow::test

// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#error "ideal-gas closure test access is unavailable in tests-off builds"
#endif

#include "hundun/flow/ideal_gas_piso.hpp"

#include <vector>

namespace hundun::flow::test {

enum class IdealGasPostEvidenceMutation : std::uint8_t {
  attempt_identity,
  shared_field,
  provenance,
  stencil,
  finalization_zero,
  finalization_not_later,
  residual_outer_size,
  residual_element,
  disposition,
  reason,
  lowest_failing_rank,
  density_availability,
  density_value,
  residual_availability,
  conservation_outer_size,
  conservation_element,
  conservation_availability,
  mass_availability,
  mass_value,
  minimum_density_availability,
  minimum_density_value,
  pre_residual_element,
  pre_conservation_element,
  pre_mass_value,
  count
};

enum class IdealGasPrepareFault : std::uint8_t {
  state_prepare,
  closure_prepare
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
  sum_reduction,
  maximum_reduction
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
  static void set_create_fault(IdealGasCreateFault, int rank) noexcept;
  static void reset_create_fault() noexcept;
  static void begin_attempt(IdealGasClosure &, FlowState &,
                            std::uint64_t identity);
  static IdealGasClosureReport evaluate(IdealGasClosure &, FlowState &,
                                        IdealGasClosureStage);
  static void rollback(IdealGasClosure &) noexcept;
  static void set_stage_failure(IdealGasClosure &, IdealGasClosureStage,
                                IdealGasClosureFailureReason, int rank);
  static void set_metric_gate_failure(IdealGasClosure &, int rank);
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
  static void set_post_assessment_corruption(FixedStepIdealGasFlow &,
                                             int rank);
  static void set_attempt_layout_fault(FixedStepIdealGasFlow &, int rank);
  static void set_outlet_backflow_fault(FixedStepIdealGasFlow &);
  static bool post_evidence_mutation_rejected(
      const IdealGasStepAttemptReport &, IdealGasPostEvidenceMutation);
};

} // namespace hundun::flow::test

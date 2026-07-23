// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#error "ideal-gas closure test access is unavailable in tests-off builds"
#endif

#include "hundun/flow/ideal_gas_piso.hpp"

namespace hundun::flow::test {

enum class IdealGasPostEvidenceMutation : std::uint8_t {
  attempt_identity,
  shared_field,
  provenance,
  stencil,
  finalization_zero,
  finalization_not_later,
  residual_outer_size,
  residual_element
};

enum class IdealGasPrepareFault : std::uint8_t {
  state_prepare,
  closure_prepare
};

enum class IdealGasOuterFailurePoint : std::uint8_t {
  momentum_after_predictor,
  pressure_after_first_corrector
};

class IdealGasClosureTestAccess final {
public:
  static void set_uniform_enthalpy_rate(FixedStepIdealGasFlow &,
                                        double rate_J_per_kg_s);
  static bool report_authenticated(const IdealGasClosureReport &) noexcept;
  static bool report_authenticated(const IdealGasStepAttemptReport &) noexcept;
  static bool post_eos_evidence_authenticated(
      const MaterialDensityStepAttemptReport &) noexcept;
  static bool same_rank_reason_precedence_is_enum_order() noexcept;
  static void exhaust_source_generation(FixedStepIdealGasFlow &);
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
  static void set_outlet_backflow_fault(FixedStepIdealGasFlow &);
  static bool post_evidence_mutation_rejected(
      const IdealGasStepAttemptReport &, IdealGasPostEvidenceMutation);
};

} // namespace hundun::flow::test

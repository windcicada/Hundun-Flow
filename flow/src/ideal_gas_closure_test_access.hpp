// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#error "ideal-gas closure test access is unavailable in tests-off builds"
#endif

#include "hundun/flow/ideal_gas_piso.hpp"

namespace hundun::flow::test {

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
};

} // namespace hundun::flow::test

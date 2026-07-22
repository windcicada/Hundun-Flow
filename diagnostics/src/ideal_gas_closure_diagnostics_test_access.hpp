// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
#error "ideal-gas diagnostic test access is unavailable in tests-off builds"
#endif

#include "hundun/diagnostics/structured_diagnostics.hpp"

#include <cstdint>

namespace hundun::diagnostics::test {

enum class IdealGasClosureDiagnosticFault : std::uint8_t {
  none,
  provider_agreement,
  ownership_layout,
  sample_wire,
  record_validation,
  status_warning,
  status_failed,
  status_unavailable
};

struct IdealGasClosureDiagnosticWork final {
  std::uint64_t observations{};
  std::uint64_t fingerprint_items{};
  std::uint64_t summary_items{};
  std::uint64_t sample_items{};
};

class IdealGasClosureDiagnosticTestAccess final {
public:
  static void set_fault(IdealGasClosureDiagnosticFault,
                        int rank = -1) noexcept;
  static void reset() noexcept;
  static IdealGasClosureDiagnosticWork work() noexcept;
  static DiagnosticInvariant
  positive_invariant(std::string id, std::string unit, double observed);
};

} // namespace hundun::diagnostics::test

// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
#error "ideal-gas diagnostic test access is unavailable in tests-off builds"
#endif

#include "hundun/diagnostics/structured_diagnostics.hpp"
#include "hundun/flow/ideal_gas_piso.hpp"

#include <cstdint>

namespace hundun::diagnostics::test {

enum class IdealGasClosureDiagnosticFault : std::uint8_t {
  none,
  provider_agreement,
  ownership_layout,
  sample_wire,
  record_validation,
  request_preparation,
  sample_preparation,
  aggregation,
  raw_mpi
};

struct IdealGasClosureDiagnosticWork final {
  std::uint64_t observations{};
  std::uint64_t fingerprint_items{};
  std::uint64_t summary_items{};
  std::uint64_t invariant_items{};
  std::uint64_t sample_items{};
  std::uint64_t retained_sample_items{};
  std::uint64_t allocation_events{};
  std::uint64_t full_field_copy_attempts{};
  std::uint64_t collective_calls{};
};

class IdealGasClosureDiagnosticTestAccess final {
public:
  static void set_fault(const flow::IdealGasClosureDiagnosticSource &,
                        IdealGasClosureDiagnosticFault,
                        int rank = -1) noexcept;
  static void reset(const flow::IdealGasClosureDiagnosticSource &) noexcept;
  static IdealGasClosureDiagnosticWork
  work(const flow::IdealGasClosureDiagnosticSource &) noexcept;
  static DiagnosticInvariant
  positive_invariant(std::string id, std::string unit, double observed);
};

} // namespace hundun::diagnostics::test

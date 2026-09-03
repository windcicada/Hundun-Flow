// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/diag_structured.hpp"
#include "hundun/flow_ideal_gas_piso.hpp"

#include <array>
#include <cstdint>
#include <utility>

namespace hundun::diagnostics::detail {

void ideal_gas_closure_diagnostic_set_fault_raw(
    const flow::IdealGasClosureDiagnosticSource &, std::uint8_t,
    int) noexcept;
void ideal_gas_closure_diagnostic_reset_raw(
    const flow::IdealGasClosureDiagnosticSource &) noexcept;
std::array<std::uint64_t, 9>
ideal_gas_closure_diagnostic_work_raw(
    const flow::IdealGasClosureDiagnosticSource &) noexcept;
DiagnosticInvariant ideal_gas_closure_positive_invariant_raw(
    std::string, std::string, double);
void ideal_gas_closure_submit_record_contract_raw(
    const flow::IdealGasClosureDiagnosticSource &, const DiagnosticRequest &,
    const DiagnosticRecord &, const DiagnosticRecord &, DiagnosticSink &);

} // namespace hundun::diagnostics::detail

namespace hundun::diagnostics::test {

enum class IdealGasClosureDiagnosticFault : std::uint8_t {
  none,
  stale_source,
  capability,
  provider_agreement,
  ownership_layout,
  global_extent,
  sample_wire,
  record_validation,
  record_preparation,
  request_preparation,
  sample_preparation,
  aggregation,
  oversized_agreement,
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
  static void set_fault(const flow::IdealGasClosureDiagnosticSource &source,
                        IdealGasClosureDiagnosticFault fault,
                        int rank = -1) noexcept {
    detail::ideal_gas_closure_diagnostic_set_fault_raw(
        source, static_cast<std::uint8_t>(fault), rank);
  }
  static void
  reset(const flow::IdealGasClosureDiagnosticSource &source) noexcept {
    detail::ideal_gas_closure_diagnostic_reset_raw(source);
  }
  static IdealGasClosureDiagnosticWork
  work(const flow::IdealGasClosureDiagnosticSource &source) noexcept {
    const auto raw = detail::ideal_gas_closure_diagnostic_work_raw(source);
    return {raw[0], raw[1], raw[2], raw[3], raw[4],
            raw[5], raw[6], raw[7], raw[8]};
  }
  static DiagnosticInvariant
  positive_invariant(std::string id, std::string unit, double observed) {
    return detail::ideal_gas_closure_positive_invariant_raw(
        std::move(id), std::move(unit), observed);
  }
  static void
  submit_record_contract(const flow::IdealGasClosureDiagnosticSource &source,
                         const DiagnosticRequest &request,
                         const DiagnosticRecord &expected,
                         const DiagnosticRecord &candidate,
                         DiagnosticSink &sink) {
    detail::ideal_gas_closure_submit_record_contract_raw(
        source, request, expected, candidate, sink);
  }
};

} // namespace hundun::diagnostics::test

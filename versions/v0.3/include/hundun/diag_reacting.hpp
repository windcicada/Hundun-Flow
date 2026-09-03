// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/diag_provider_registry.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace hundun::diagnostics {

struct ReactingDiagnosticsSnapshot final {
  std::uint64_t revision{};
  std::uint64_t composition_fingerprint{};
  std::string mechanism_sha256;
  double final_residual{};
  double eos_relative_drift{};
  std::vector<double> element_budget_kg;
  std::vector<double> species_budget_kg;
  double enthalpy_budget_j{};
  std::uint64_t chemistry_calls{};
  std::uint64_t chemistry_failures{};
  std::uint64_t pressure_corrector_calls{};
  std::uint64_t workspace_count{};
  std::uint64_t rollback_count{};
  bool accepted_step{};
  DiagnosticFailureClass failure_class{DiagnosticFailureClass::none};
  std::string failure_code{"none"};
};

DiagnosticDescriptor
describe_diagnostics(const ReactingDiagnosticsSnapshot &) noexcept;
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const ReactingDiagnosticsSnapshot &);
void collect_diagnostics(const ReactingDiagnosticsSnapshot &,
                         const DiagnosticRequest &, DiagnosticSink &);
DiagnosticProvider
make_reacting_diagnostic_provider(const ReactingDiagnosticsSnapshot &);

} // namespace hundun::diagnostics

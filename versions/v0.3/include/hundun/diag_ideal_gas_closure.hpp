// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/diag_structured.hpp"
#include "hundun/flow_ideal_gas_piso.hpp"
#include "hundun/rt_mpi_context.hpp"

#include <vector>

namespace hundun::diagnostics {

DiagnosticDescriptor describe_ideal_gas_closure_diagnostics() noexcept;
DiagnosticDescriptor
describe_diagnostics(const flow::IdealGasClosureDiagnosticSource &) noexcept;
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const flow::IdealGasClosureDiagnosticSource &);
void collect_diagnostics(const flow::IdealGasClosureDiagnosticSource &,
                         const DiagnosticRequest &, DiagnosticSink &);
void collect_diagnostics(const flow::IdealGasClosureDiagnosticSource &,
                         const runtime::MpiContext &, const DiagnosticRequest &,
                         DiagnosticSink &);

} // namespace hundun::diagnostics

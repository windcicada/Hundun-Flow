// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/diag_structured.hpp"
#include "hundun/flow_adaptive_time_control.hpp"
#include "hundun/rt_mpi_context.hpp"

#include <vector>

namespace hundun::diagnostics {

DiagnosticDescriptor describe_time_control_diagnostics() noexcept;
DiagnosticDescriptor
describe_diagnostics(const flow::TimeControlDiagnosticSource &) noexcept;
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const flow::TimeControlDiagnosticSource &);
void collect_diagnostics(const flow::TimeControlDiagnosticSource &,
                         const DiagnosticRequest &, DiagnosticSink &);
void collect_diagnostics(const flow::TimeControlDiagnosticSource &,
                         const runtime::MpiContext &, const DiagnosticRequest &,
                         DiagnosticSink &);

} // namespace hundun::diagnostics

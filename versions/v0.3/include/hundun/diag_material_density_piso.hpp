// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/diag_structured.hpp"

#include <string_view>
#include <vector>

namespace hundun::flow {
class MaterialDensityFlowDiagnosticSource;
}
namespace hundun::runtime {
class MpiContext;
}

namespace hundun::diagnostics {

DiagnosticDescriptor describe_diagnostics(
    const flow::MaterialDensityFlowDiagnosticSource &) noexcept;

std::vector<std::string_view> diagnostic_fingerprint_field_ids(
    const flow::MaterialDensityFlowDiagnosticSource &);

void collect_diagnostics(const flow::MaterialDensityFlowDiagnosticSource &,
                         const DiagnosticRequest &, DiagnosticSink &);

void collect_diagnostics(const flow::MaterialDensityFlowDiagnosticSource &,
                         const runtime::MpiContext &,
                         const DiagnosticRequest &, DiagnosticSink &);

} // namespace hundun::diagnostics

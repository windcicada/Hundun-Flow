// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/diag_structured.hpp"

#include <string_view>
#include <vector>

namespace hundun::flow {
class MaterialDensityDiagnosticSource;
}
namespace hundun::runtime {
class MpiContext;
}

namespace hundun::diagnostics {

DiagnosticDescriptor
describe_diagnostics(const flow::MaterialDensityDiagnosticSource &) noexcept;

std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const flow::MaterialDensityDiagnosticSource &);

void collect_diagnostics(const flow::MaterialDensityDiagnosticSource &,
                         const DiagnosticRequest &, DiagnosticSink &);

void collect_diagnostics(const flow::MaterialDensityDiagnosticSource &,
                         const runtime::MpiContext &, const DiagnosticRequest &,
                         DiagnosticSink &);

} // namespace hundun::diagnostics

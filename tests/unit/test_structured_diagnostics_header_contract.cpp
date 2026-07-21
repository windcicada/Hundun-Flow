// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/structured_diagnostics.hpp"

#include <type_traits>

static_assert(
    std::is_same_v<decltype(hundun::diagnostics::kDiagnosticRecordSchemaV1),
                   const std::uint32_t>);
static_assert(hundun::diagnostics::kMaximumStateSamplesV1 == 256U);
static_assert(
    std::has_virtual_destructor_v<hundun::diagnostics::DiagnosticSink>);

int main() {
  hundun::diagnostics::DiagnosticRequest request;
  request.level = hundun::diagnostics::DiagnosticLevel::summary;
  return request.sample_budget == 0U ? 0 : 1;
}

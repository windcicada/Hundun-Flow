// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/diagnostics/structured_diagnostics.hpp"
#include "hundun/flow/checkpoint_v2.hpp"
#include "hundun/runtime/mpi_context.hpp"

#include <vector>

namespace hundun::diagnostics {

DiagnosticDescriptor describe_checkpoint_v2_diagnostics() noexcept;
DiagnosticDescriptor
describe_diagnostics(const flow::CheckpointV2DiagnosticSource &) noexcept;
std::vector<std::string_view> diagnostic_fingerprint_field_ids(
    const flow::CheckpointV2DiagnosticSource &);
void collect_diagnostics(const flow::CheckpointV2DiagnosticSource &,
                         const DiagnosticRequest &, DiagnosticSink &);
void collect_diagnostics(const flow::CheckpointV2DiagnosticSource &,
                         const runtime::MpiContext &, const DiagnosticRequest &,
                         DiagnosticSink &);

} // namespace hundun::diagnostics

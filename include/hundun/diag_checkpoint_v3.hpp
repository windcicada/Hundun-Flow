// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/diag_structured.hpp"
#include "hundun/flow_checkpoint_v3.hpp"
#include "hundun/rt_mpi_context.hpp"

#include <string_view>
#include <vector>

namespace hundun::diagnostics {

DiagnosticDescriptor
describe_diagnostics(const flow::CheckpointV3Report &) noexcept;
std::vector<std::string_view>
diagnostic_fingerprint_field_ids(const flow::CheckpointV3Report &);
std::vector<DiagnosticModuleKind>
stage3_added_provider_inventory(flow::CheckpointV3Presence);
std::vector<DiagnosticModuleKind>
stage3_added_provider_inventory(const flow::CheckpointV3Report &);
void collect_diagnostics(const flow::CheckpointV3Report &,
                         const DiagnosticRequest &, DiagnosticSink &);
void collect_diagnostics(const flow::CheckpointV3Report &,
                         const runtime::MpiContext &,
                         const DiagnosticRequest &, DiagnosticSink &);

} // namespace hundun::diagnostics

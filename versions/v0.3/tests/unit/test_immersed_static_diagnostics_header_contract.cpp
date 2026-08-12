// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_immersed_static.hpp"
#include "hundun/diag_immersed_module.hpp"
#include "hundun/diag_checkpoint_v3.hpp"

#include <type_traits>
#include <utility>

namespace {

using Summary = hundun::diagnostics::ImmersedStaticDiagnosticSummary;
using Kind = hundun::diagnostics::DiagnosticModuleKind;
using Request = hundun::diagnostics::DiagnosticRequest;
using Sink = hundun::diagnostics::DiagnosticSink;
using Mpi = hundun::runtime::MpiContext;

using Summarize = Summary (*)(
    const hundun::immersed::ImmersedSurface &,
    const hundun::immersed::SurfaceQuery &,
    const hundun::immersed::ImmersedDomain &,
    const hundun::immersed::GhostStencilPlan &,
    const hundun::immersed::WallQuadraturePlan &,
    const hundun::immersed::LocalFlowPatternTransform &);
using Inventory = std::vector<Kind> (*)(hundun::flow::CheckpointV3Presence);
using WithLocalFlowPattern = Summary (*)(
    Summary, const hundun::flow::ImmersedFlowDiagnosticSource &);
using ReportInventory = std::vector<Kind> (*)(
    const hundun::flow::CheckpointV3Report &);
using Describe = hundun::diagnostics::DiagnosticDescriptor (*)(
    const Summary &, Kind);
using Fields = std::vector<std::string_view> (*)(const Summary &, Kind);
using CollectLocal = void (*)(const Summary &, Kind, const Request &, Sink &);
using CollectCollective = void (*)(const Summary &, Kind, const Mpi &,
                                   const Request &, Sink &);
using AttemptDescribe = hundun::diagnostics::DiagnosticDescriptor (*)(
    const hundun::flow::ImmersedFlowDiagnosticSource &, Kind);
using AttemptFields = std::vector<std::string_view> (*)(
    const hundun::flow::ImmersedFlowDiagnosticSource &, Kind);
using AttemptCollectLocal = void (*)(
    const hundun::flow::ImmersedFlowDiagnosticSource &, Kind,
    const Request &, Sink &);
using AttemptCollectCollective = void (*)(
    const hundun::flow::ImmersedFlowDiagnosticSource &, Kind, const Mpi &,
    const Request &, Sink &);

static_assert(std::is_final_v<Summary>);
static_assert(std::is_standard_layout_v<Summary>);
static_assert(std::is_trivially_copyable_v<Summary>);
static_assert(std::is_same_v<
              decltype(&hundun::diagnostics::summarize_immersed_static),
              Summarize>);
static_assert(std::is_same_v<
              decltype(&hundun::diagnostics::with_local_flow_pattern_snapshot),
              WithLocalFlowPattern>);
static_assert(std::is_same_v<
              decltype(static_cast<Inventory>(
                  &hundun::diagnostics::stage3_added_provider_inventory)),
              Inventory>);
static_assert(std::is_same_v<
              decltype(static_cast<ReportInventory>(
                  &hundun::diagnostics::stage3_added_provider_inventory)),
              ReportInventory>);
static_assert(std::is_same_v<
              decltype(static_cast<Describe>(
                  &hundun::diagnostics::describe_diagnostics)),
              Describe>);
static_assert(std::is_same_v<
              decltype(static_cast<Fields>(
                  &hundun::diagnostics::diagnostic_fingerprint_field_ids)),
              Fields>);
static_assert(std::is_same_v<
              decltype(static_cast<CollectLocal>(
                  &hundun::diagnostics::collect_diagnostics)),
              CollectLocal>);
static_assert(std::is_same_v<
              decltype(static_cast<CollectCollective>(
                  &hundun::diagnostics::collect_diagnostics)),
              CollectCollective>);
static_assert(std::is_same_v<
              decltype(static_cast<AttemptDescribe>(
                  &hundun::diagnostics::describe_diagnostics)),
              AttemptDescribe>);
static_assert(std::is_same_v<
              decltype(static_cast<AttemptFields>(
                  &hundun::diagnostics::diagnostic_fingerprint_field_ids)),
              AttemptFields>);
static_assert(std::is_same_v<
              decltype(static_cast<AttemptCollectLocal>(
                  &hundun::diagnostics::collect_diagnostics)),
              AttemptCollectLocal>);
static_assert(std::is_same_v<
              decltype(static_cast<AttemptCollectCollective>(
                  &hundun::diagnostics::collect_diagnostics)),
              AttemptCollectCollective>);

} // namespace

int main() { return 0; }

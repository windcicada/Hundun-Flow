// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/diag_les_wale.hpp"

#include <cstdint>
#include <optional>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

using hundun::diagnostics::DiagnosticDescriptor;
using hundun::diagnostics::DiagnosticModuleKind;
using hundun::diagnostics::DiagnosticRequest;
using hundun::diagnostics::DiagnosticSink;
using hundun::les::WaleSummary;

static_assert(static_cast<std::uint16_t>(DiagnosticModuleKind::runtime) == 0U);
static_assert(static_cast<std::uint16_t>(DiagnosticModuleKind::mpi) == 1U);
static_assert(static_cast<std::uint16_t>(
                  DiagnosticModuleKind::mesh_topology) == 2U);
static_assert(static_cast<std::uint16_t>(
                  DiagnosticModuleKind::mesh_geometry) == 3U);
static_assert(static_cast<std::uint16_t>(DiagnosticModuleKind::field) == 4U);
static_assert(static_cast<std::uint16_t>(DiagnosticModuleKind::execution) ==
              5U);
static_assert(static_cast<std::uint16_t>(
                  DiagnosticModuleKind::linear_operator) == 6U);
static_assert(static_cast<std::uint16_t>(
                  DiagnosticModuleKind::linear_solver) == 7U);
static_assert(static_cast<std::uint16_t>(DiagnosticModuleKind::halo) == 8U);
static_assert(static_cast<std::uint16_t>(DiagnosticModuleKind::boundary) ==
              9U);
static_assert(static_cast<std::uint16_t>(
                  DiagnosticModuleKind::finite_volume) == 10U);
static_assert(static_cast<std::uint16_t>(DiagnosticModuleKind::piso) == 11U);
static_assert(static_cast<std::uint16_t>(
                  DiagnosticModuleKind::density_transport) == 12U);
static_assert(static_cast<std::uint16_t>(
                  DiagnosticModuleKind::density_closure) == 13U);
static_assert(static_cast<std::uint16_t>(
                  DiagnosticModuleKind::time_control) == 14U);
static_assert(static_cast<std::uint16_t>(DiagnosticModuleKind::checkpoint) ==
              15U);
static_assert(static_cast<std::uint16_t>(
                  DiagnosticModuleKind::flow_driver) == 16U);
static_assert(static_cast<std::uint16_t>(DiagnosticModuleKind::performance) ==
              17U);
static_assert(static_cast<std::uint16_t>(
                  DiagnosticModuleKind::immersed_surface) == 18U);
static_assert(static_cast<std::uint16_t>(
                  DiagnosticModuleKind::ghost_stencil) == 19U);
static_assert(static_cast<std::uint16_t>(
                  DiagnosticModuleKind::local_flow_pattern) == 20U);
static_assert(static_cast<std::uint16_t>(DiagnosticModuleKind::wall_force) ==
              21U);
static_assert(static_cast<std::uint16_t>(DiagnosticModuleKind::les) == 22U);

using DescribeSignature =
    DiagnosticDescriptor (*)(const WaleSummary &) noexcept;
using FieldsSignature =
    std::vector<std::string_view> (*)(const WaleSummary &);
using CollectSignature = void (*)(const WaleSummary &,
                                  const DiagnosticRequest &, DiagnosticSink &);

static_assert(std::is_same_v<
              decltype(static_cast<DescribeSignature>(
                  &hundun::diagnostics::describe_diagnostics)),
              DescribeSignature>);
static_assert(std::is_same_v<
              decltype(static_cast<FieldsSignature>(
                  &hundun::diagnostics::diagnostic_fingerprint_field_ids)),
              FieldsSignature>);
static_assert(std::is_same_v<
              decltype(static_cast<CollectSignature>(
                  &hundun::diagnostics::collect_diagnostics)),
              CollectSignature>);
static_assert(!std::is_invocable_v<DescribeSignature,
                                   const std::optional<WaleSummary> &>);
static_assert(!std::is_invocable_v<CollectSignature,
                                   const std::optional<WaleSummary> &,
                                   const DiagnosticRequest &, DiagnosticSink &>);

int main() { return 0; }

// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/material_density_transport_diagnostics.hpp"
#include "hundun/flow/material_density_transport.hpp"

#include <string_view>
#include <type_traits>
#include <utility>

static_assert(
    !std::is_copy_constructible_v<hundun::flow::MaterialFaceMassFlux>);
static_assert(std::is_move_constructible_v<hundun::flow::MaterialFaceMassFlux>);
static_assert(
    !std::is_copy_constructible_v<hundun::flow::MaterialDensityTransport>);
static_assert(
    !std::is_aggregate_v<hundun::flow::MaterialDensityTransportReport>);
static_assert(
    std::is_copy_constructible_v<hundun::flow::MaterialDensityTransportReport>);
static_assert(
    std::is_move_constructible_v<hundun::flow::MaterialDensityTransportReport>);
static_assert(
    std::is_same_v<
        decltype(std::declval<
                     const hundun::flow::MaterialDensityDiagnosticSource &>()
                     .owned_cell_layout_fingerprint()),
        std::string_view>);
static_assert(
    std::is_same_v<
        decltype(std::declval<
                     const hundun::flow::MaterialDensityDiagnosticSource &>()
                     .global_cell_layout_fingerprint()),
        std::string_view>);

int main() {
  const hundun::flow::MaterialDensityTransportSpec spec{};
  return spec.scalar_densities.empty() ? 0 : 1;
}

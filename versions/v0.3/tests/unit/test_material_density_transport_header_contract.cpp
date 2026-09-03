// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/diag_material_density_transport.hpp"
#include "hundun/flow_material_density_transport.hpp"

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
static_assert(
    std::is_same_v<
        decltype(std::declval<
                     const hundun::flow::MaterialDensityDiagnosticSource &>()
                     .global_cell_extent()),
        hundun::runtime::Int3>);
static_assert(
    std::is_same_v<
        decltype(std::declval<
                     const hundun::flow::MaterialDensityDiagnosticSource &>()
                     .owned_global_box()),
        hundun::runtime::Box3>);

int main() {
  const hundun::flow::MaterialDensityTransportSpec spec{};
  return spec.scalar_densities.empty() ? 0 : 1;
}

// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/material_density_transport_diagnostics.hpp"
#include "hundun/flow/material_density_transport.hpp"

#include <type_traits>

static_assert(
    !std::is_copy_constructible_v<hundun::flow::MaterialFaceMassFlux>);
static_assert(std::is_move_constructible_v<hundun::flow::MaterialFaceMassFlux>);
static_assert(
    !std::is_copy_constructible_v<hundun::flow::MaterialDensityTransport>);

int main() {
  const hundun::flow::MaterialDensityTransportSpec spec{};
  return spec.scalar_densities.empty() ? 0 : 1;
}

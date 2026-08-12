// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_material_density_piso.hpp"
#include "hundun/flow_material_density_piso.hpp"

#include <type_traits>

int main() {
  using Flow = hundun::flow::FixedStepMaterialDensityFlow;
  static_assert(!std::is_copy_constructible_v<Flow>);
  return 0;
}

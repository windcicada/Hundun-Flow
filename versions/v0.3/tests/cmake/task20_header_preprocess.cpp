// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/diag_material_density_piso.hpp"
#include "hundun/flow_material_density_piso.hpp"

#include <type_traits>

int main() {
  using Flow = hundun::flow::FixedStepMaterialDensityFlow;
  static_assert(!std::is_copy_constructible_v<Flow>);
  return 0;
}

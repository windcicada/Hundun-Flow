// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/diag_ideal_gas_closure.hpp"
#include "hundun/flow_ideal_gas_closure.hpp"
#include "hundun/flow_ideal_gas_piso.hpp"

#include <type_traits>

int main() {
  static_assert(
      !std::is_copy_constructible_v<hundun::flow::FixedStepIdealGasFlow>);
  static_assert(!std::is_copy_constructible_v<
                hundun::flow::IdealGasClosureDiagnosticSource>);
  return 0;
}

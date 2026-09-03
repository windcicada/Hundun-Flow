// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/diag_material_density_transport.hpp"
#include "hundun/diag_structured.hpp"
#include "hundun/flow_material_density_transport.hpp"
#include "hundun/rt_mpi_operation_error.hpp"

int main() {
  try {
    hundun::runtime::check_mpi_result(1, "typed-result-check");
  } catch (const hundun::runtime::MpiOperationError &) {
    return 0;
  }
  return 1;
}

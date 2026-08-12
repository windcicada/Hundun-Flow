// SPDX-License-Identifier: Apache-2.0

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

// SPDX-License-Identifier: Apache-2.0

#include "hundun/diagnostics/material_density_transport_diagnostics.hpp"
#include "hundun/diagnostics/structured_diagnostics.hpp"
#include "hundun/flow/material_density_transport.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"

int main() {
  try {
    hundun::runtime::check_mpi_result(1, "typed-result-check");
  } catch (const hundun::runtime::MpiOperationError &) {
    return 0;
  }
  return 1;
}

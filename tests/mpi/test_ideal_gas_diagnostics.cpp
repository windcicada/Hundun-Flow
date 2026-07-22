// SPDX-License-Identifier: Apache-2.0

#include "diagnostics/src/ideal_gas_closure_diagnostics_test_access.hpp"
#include "hundun/diagnostics/ideal_gas_closure_diagnostics.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "tests/support/test_main.hpp"

#include <string_view>

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(mpi.size() == 1 || mpi.size() == 2 || mpi.size() == 4);
    const auto descriptor =
        hundun::diagnostics::describe_ideal_gas_closure_diagnostics();
    HUNDUN_CHECK(descriptor.module_kind ==
                 hundun::diagnostics::DiagnosticModuleKind::density_closure);
    HUNDUN_CHECK(descriptor.module_id ==
                 std::string_view("flow.ideal-gas-closure"));
    HUNDUN_CHECK(descriptor.instance_id == std::string_view("primary"));
  });
}

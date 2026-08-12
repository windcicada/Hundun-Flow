// SPDX-License-Identifier: Apache-2.0

#include "src/app_flow_driver_detail.hpp"
#include "hundun/diag_session.hpp"
#include "hundun/diag_mesh_v2.hpp"
#include "hundun/diag_module.hpp"

#include "tests/support/test_main.hpp"

#include <type_traits>

int main() {
  return hundun::test::run([] {
    static_assert(
        std::is_move_constructible_v<hundun::diagnostics::DiagnosticSession>);
    static_assert(
        std::is_move_constructible_v<hundun::diagnostics::MeshDiagnosticV2>);
    HUNDUN_CHECK(true);
  });
}

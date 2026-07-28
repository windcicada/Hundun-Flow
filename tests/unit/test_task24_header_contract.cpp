// SPDX-License-Identifier: Apache-2.0

#include "applications/hundun/stage2_driver.hpp"
#include "hundun/diagnostics/diagnostic_session.hpp"
#include "hundun/diagnostics/mesh_diagnostic_v2.hpp"
#include "hundun/diagnostics/stage2_module_diagnostics.hpp"

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

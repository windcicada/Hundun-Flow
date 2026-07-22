// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
#error "ideal-gas diagnostic test access is unavailable in tests-off builds"
#endif

namespace hundun::diagnostics::test {

class IdealGasClosureDiagnosticTestAccess final {};

} // namespace hundun::diagnostics::test

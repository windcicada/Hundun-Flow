// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
#error "material-density PISO diagnostics test access requires tests-on build"
#endif

#include <cstdint>

namespace hundun::diagnostics::test {

enum class MaterialDensityPisoDiagnosticFault : std::uint8_t {
  none,
  provider_agreement,
  cell_exact_cover,
  face_exact_cover,
  sample_send_preparation,
  sample_receive_preparation,
  record_validation
};

class MaterialDensityPisoDiagnosticTestAccess final {
public:
  static void set_fault(MaterialDensityPisoDiagnosticFault) noexcept;
  static void reset() noexcept;
};

} // namespace hundun::diagnostics::test

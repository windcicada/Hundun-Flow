// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
#error "material-density PISO diagnostics test access requires tests-on build"
#endif

#include "hundun/diagnostics/structured_diagnostics.hpp"

#include <cstdint>
#include <vector>

namespace hundun::diagnostics::test {

enum class MaterialDensityPisoDiagnosticFault : std::uint8_t {
  none,
  provider_agreement,
  cell_exact_cover,
  face_exact_cover,
  sample_send_preparation,
  sample_receive_preparation,
  sample_wire_malformed,
  record_validation
};

struct MaterialDensityPisoSampleWireItem final {
  std::uint32_t field{};
  std::uint64_t global_id{};
  std::uint32_t component{};
  std::uint64_t value_bits{};
  DiagnosticValueStatus status{DiagnosticValueStatus::finite};
};

class MaterialDensityPisoDiagnosticTestAccess final {
public:
  static void set_fault(MaterialDensityPisoDiagnosticFault) noexcept;
  static void reset() noexcept;
  static std::vector<unsigned char> encode_sample_wire(
      const std::vector<MaterialDensityPisoSampleWireItem> &items);
  static std::vector<MaterialDensityPisoSampleWireItem>
  decode_sample_wire(const std::vector<unsigned char> &bytes);
};

} // namespace hundun::diagnostics::test

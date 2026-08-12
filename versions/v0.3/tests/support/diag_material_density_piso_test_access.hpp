// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/diag_structured.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace hundun::diagnostics::detail {

using MaterialDensityPisoSampleWireRawItem = std::array<std::uint64_t, 5>;

void material_density_piso_set_fault_for_test(std::uint8_t) noexcept;
void material_density_piso_reset_for_test() noexcept;
std::vector<unsigned char> material_density_piso_encode_sample_wire_for_test(
    const std::vector<MaterialDensityPisoSampleWireRawItem> &items);
std::vector<MaterialDensityPisoSampleWireRawItem>
material_density_piso_decode_sample_wire_for_test(
    const std::vector<unsigned char> &bytes);

} // namespace hundun::diagnostics::detail

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
  static void set_fault(MaterialDensityPisoDiagnosticFault fault) noexcept {
    detail::material_density_piso_set_fault_for_test(
        static_cast<std::uint8_t>(fault));
  }
  static void reset() noexcept {
    detail::material_density_piso_reset_for_test();
  }
  static std::vector<unsigned char> encode_sample_wire(
      const std::vector<MaterialDensityPisoSampleWireItem> &items) {
    std::vector<detail::MaterialDensityPisoSampleWireRawItem> raw;
    raw.reserve(items.size());
    for (const auto &item : items)
      raw.push_back({item.field, item.global_id, item.component,
                     item.value_bits, static_cast<std::uint8_t>(item.status)});
    return detail::material_density_piso_encode_sample_wire_for_test(raw);
  }
  static std::vector<MaterialDensityPisoSampleWireItem>
  decode_sample_wire(const std::vector<unsigned char> &bytes) {
    auto raw =
        detail::material_density_piso_decode_sample_wire_for_test(bytes);
    std::vector<MaterialDensityPisoSampleWireItem> result;
    result.reserve(raw.size());
    for (const auto &item : raw)
      result.push_back({static_cast<std::uint32_t>(item[0]), item[1],
                        static_cast<std::uint32_t>(item[2]), item[3],
                        static_cast<DiagnosticValueStatus>(item[4])});
    return result;
  }
};

} // namespace hundun::diagnostics::test

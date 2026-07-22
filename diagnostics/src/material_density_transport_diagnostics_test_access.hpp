// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
#error "material-density diagnostics test access is test-build only"
#endif

#include <cstddef>
#include <cstdint>
#include <vector>

namespace hundun::flow {
class MaterialDensityDiagnosticSource;
}

namespace hundun::diagnostics::test {

enum class MaterialDiagnosticAllocationPoint : std::uint8_t {
  summary_gather,
  transport_totals,
  owned_id_local,
  ownership_counts,
  ownership_gather,
  eligible_counts,
  local_sample_wire_and_size_counts,
  sample_exchange_buffers,
  decoded_and_retained_samples,
  provider_key,
  provider_reference
};

enum class MaterialDiagnosticRawCollectivePoint : std::uint8_t {
  none,
  other,
  preparation_summary_gather,
  preparation_transport_totals,
  preparation_owned_id_local,
  preparation_ownership_counts,
  preparation_ownership_gather,
  preparation_eligible_counts,
  preparation_local_sample_wire_and_size_counts,
  preparation_sample_exchange_buffers,
  preparation_decoded_and_retained_samples,
  sample_size_exchange,
  preparation_provider_key,
  preparation_provider_reference
};

struct MaterialDiagnosticWorkCounts final {
  std::uint64_t field_reads{};
  std::uint64_t volume_reads{};
  std::uint64_t summary_accumulations{};
  std::uint64_t sample_candidates{};
  std::uint64_t fingerprint_items{};
  std::uint64_t raw_collectives{};
  MaterialDiagnosticRawCollectivePoint last_raw_collective{
      MaterialDiagnosticRawCollectivePoint::none};
};

class MaterialDensityTransportDiagnosticsTestAccess final {
public:
  static void reset() noexcept;
  static MaterialDiagnosticWorkCounts work_counts() noexcept;
  static void override_global_id(std::size_t local_cell,
                                 std::uint64_t global_id,
                                 int rank = -1) noexcept;
  static void inject_provider_failure(int rank = -1) noexcept;
  static void inject_record_failure(int rank = -1) noexcept;
  static void inject_request_size_overflow(int rank = -1) noexcept;
  static void inject_sample_wire_overflow(int rank = -1) noexcept;
  static void inject_allocation_failure(MaterialDiagnosticAllocationPoint,
                                        int rank = -1) noexcept;
  static void override_reported_sample_wire_bytes(const std::uint64_t *counts,
                                                  std::size_t count);
  static void inject_provider_key_value_difference(int rank = -1) noexcept;
  static std::vector<unsigned char>
  provider_key_bytes(const flow::MaterialDensityDiagnosticSource &);
};

} // namespace hundun::diagnostics::test

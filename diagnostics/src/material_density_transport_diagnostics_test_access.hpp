// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_DIAGNOSTICS_ENABLE_TEST_ACCESS
#error "material-density diagnostics test access is test-build only"
#endif

#include <cstddef>
#include <cstdint>

namespace hundun::diagnostics::test {

struct MaterialDiagnosticWorkCounts final {
  std::uint64_t field_reads{};
  std::uint64_t volume_reads{};
  std::uint64_t summary_accumulations{};
  std::uint64_t sample_candidates{};
  std::uint64_t fingerprint_items{};
  std::uint64_t raw_collectives{};
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
};

} // namespace hundun::diagnostics::test

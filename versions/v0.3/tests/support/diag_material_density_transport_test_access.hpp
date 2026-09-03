// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/diag_structured.hpp"
#include "hundun/rt_types.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace hundun::flow {
class MaterialDensityDiagnosticSource;
}

namespace hundun::diagnostics::detail {

void material_diagnostic_reset_raw() noexcept;
std::array<std::uint64_t, 6> material_diagnostic_work_counts_raw() noexcept;
std::uint8_t material_diagnostic_last_raw_collective_raw() noexcept;
void material_diagnostic_override_global_id_raw(std::size_t, std::uint64_t,
                                                int) noexcept;
void material_diagnostic_override_owned_box_raw(runtime::Box3, int) noexcept;
void material_diagnostic_inject_flag_raw(std::uint8_t, int) noexcept;
void material_diagnostic_inject_allocation_failure_raw(std::uint8_t,
                                                       int) noexcept;
void material_diagnostic_override_request_wire_raw(const unsigned char *,
                                                   std::size_t, int);
void material_diagnostic_override_sample_wire_raw(const unsigned char *,
                                                  std::size_t, int);
std::vector<unsigned char>
material_diagnostic_request_wire_raw(const DiagnosticRequest &);
bool material_diagnostic_request_wire_valid_raw(
    const std::vector<unsigned char> &) noexcept;
std::vector<unsigned char> material_diagnostic_request_wire_round_trip_raw(
    const std::vector<unsigned char> &);
std::vector<unsigned char>
material_diagnostic_sample_wire_raw(const DiagnosticSample &);
bool material_diagnostic_sample_wire_valid_raw(
    const std::vector<unsigned char> &) noexcept;
DiagnosticSample material_diagnostic_decode_single_sample_raw(
    const std::vector<unsigned char> &);
std::vector<unsigned char> material_diagnostic_provider_key_raw(
    const flow::MaterialDensityDiagnosticSource &);

} // namespace hundun::diagnostics::detail

namespace hundun::diagnostics::test {

enum class MaterialDiagnosticAllocationPoint : std::uint8_t {
  summary_gather,
  transport_totals,
  ownership_box_proof,
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
  preparation_ownership_box_proof,
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
  static void reset() noexcept { detail::material_diagnostic_reset_raw(); }
  static MaterialDiagnosticWorkCounts work_counts() noexcept {
    const auto raw = detail::material_diagnostic_work_counts_raw();
    return {raw[0], raw[1], raw[2], raw[3], raw[4], raw[5],
            static_cast<MaterialDiagnosticRawCollectivePoint>(
                detail::material_diagnostic_last_raw_collective_raw())};
  }
  static void override_global_id(std::size_t local_cell,
                                 std::uint64_t global_id,
                                 int rank = -1) noexcept {
    detail::material_diagnostic_override_global_id_raw(local_cell, global_id,
                                                       rank);
  }
  static void override_owned_global_box(runtime::Box3 box,
                                        int rank = -1) noexcept {
    detail::material_diagnostic_override_owned_box_raw(box, rank);
  }
  static void inject_transport_total_count_overflow(int rank = -1) noexcept {
    detail::material_diagnostic_inject_flag_raw(0U, rank);
  }
  static void inject_provider_failure(int rank = -1) noexcept {
    detail::material_diagnostic_inject_flag_raw(1U, rank);
  }
  static void inject_record_failure(int rank = -1) noexcept {
    detail::material_diagnostic_inject_flag_raw(2U, rank);
  }
  static void inject_request_size_overflow(int rank = -1) noexcept {
    detail::material_diagnostic_inject_flag_raw(3U, rank);
  }
  static void inject_sample_wire_overflow(int rank = -1) noexcept {
    detail::material_diagnostic_inject_flag_raw(4U, rank);
  }
  static void inject_sample_wire_prefix_overflow() noexcept {
    detail::material_diagnostic_inject_flag_raw(5U, -1);
  }
  static void inject_allocation_failure(MaterialDiagnosticAllocationPoint,
                                        int rank = -1) noexcept;
  static void inject_provider_key_value_difference(int rank = -1) noexcept {
    detail::material_diagnostic_inject_flag_raw(6U, rank);
  }
  static void override_request_wire_bytes(const unsigned char *, std::size_t,
                                          int rank = -1);
  static void override_sample_wire_bytes(const unsigned char *, std::size_t,
                                         int rank = -1);
  static std::vector<unsigned char>
  request_wire_bytes(const DiagnosticRequest &request) {
    return detail::material_diagnostic_request_wire_raw(request);
  }
  static bool
  request_wire_is_valid(const std::vector<unsigned char> &bytes) noexcept {
    return detail::material_diagnostic_request_wire_valid_raw(bytes);
  }
  static std::vector<unsigned char>
  request_wire_round_trip(const std::vector<unsigned char> &bytes) {
    return detail::material_diagnostic_request_wire_round_trip_raw(bytes);
  }
  static std::vector<unsigned char>
  sample_wire_bytes(const DiagnosticSample &sample) {
    return detail::material_diagnostic_sample_wire_raw(sample);
  }
  static bool
  sample_wire_is_valid(const std::vector<unsigned char> &bytes) noexcept {
    return detail::material_diagnostic_sample_wire_valid_raw(bytes);
  }
  static DiagnosticSample
  decode_single_sample_wire(const std::vector<unsigned char> &bytes) {
    return detail::material_diagnostic_decode_single_sample_raw(bytes);
  }
  static std::vector<unsigned char>
  provider_key_bytes(const flow::MaterialDensityDiagnosticSource &source) {
    return detail::material_diagnostic_provider_key_raw(source);
  }
};

inline void MaterialDensityTransportDiagnosticsTestAccess::
    inject_allocation_failure(MaterialDiagnosticAllocationPoint point,
                              int rank) noexcept {
  detail::material_diagnostic_inject_allocation_failure_raw(
      static_cast<std::uint8_t>(point), rank);
}

inline void MaterialDensityTransportDiagnosticsTestAccess::
    override_request_wire_bytes(const unsigned char *bytes, std::size_t size,
                                int rank) {
  detail::material_diagnostic_override_request_wire_raw(bytes, size, rank);
}

inline void MaterialDensityTransportDiagnosticsTestAccess::
    override_sample_wire_bytes(const unsigned char *bytes, std::size_t size,
                               int rank) {
  detail::material_diagnostic_override_sample_wire_raw(bytes, size, rank);
}

} // namespace hundun::diagnostics::test

// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/fvm_immersed_reconstruction.hpp"
#include "hundun/mesh_topology.hpp"
#include "hundun/rt_error.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace hundun::finite_volume::detail {

std::vector<std::uint64_t>
immersed_reconstruction_snapshot_raw(const void *identity);
void immersed_reconstruction_reset_trace_raw();
std::vector<std::uint8_t> immersed_reconstruction_trace_raw();
std::size_t immersed_reconstruction_inactive_read_attempts_raw() noexcept;
void immersed_reconstruction_fail_on_inactive_read_raw(bool enabled) noexcept;
void immersed_reconstruction_force_momentum_fallback_raw(bool enabled) noexcept;

} // namespace hundun::finite_volume::detail

namespace hundun::finite_volume::test {

enum class ReconstructionExecutionStage : std::uint8_t {
  halo_begin,
  active_interior,
  halo_wait,
  remote_ghost_symbols,
  partition_boundary,
  immersed_interface
};

struct ImmersedReconstructionPlanSnapshot final {
  std::vector<mesh::GlobalCellId> interior_rows;
  std::vector<mesh::GlobalCellId> partition_rows;
  std::vector<mesh::GlobalCellId> interface_rows;
  std::vector<mesh::GlobalFaceId> interface_faces;
  std::size_t owned_active_count{};
  std::uint64_t fingerprint{};
};

class ImmersedReconstructionTestAccess final {
public:
  static ImmersedReconstructionPlanSnapshot
  snapshot(const ImmersedReconstruction &provider) {
    const auto encoded = detail::immersed_reconstruction_snapshot_raw(
        static_cast<const void *>(provider.impl_.get()));
    if (encoded.size() < 6U)
      throw runtime::Error("immersed reconstruction snapshot is unavailable");
    const std::size_t interior = static_cast<std::size_t>(encoded[0]);
    const std::size_t partition = static_cast<std::size_t>(encoded[1]);
    const std::size_t interface = static_cast<std::size_t>(encoded[2]);
    const std::size_t faces = static_cast<std::size_t>(encoded[3]);
    if (interior > encoded.size() - 6U ||
        partition > encoded.size() - 6U - interior ||
        interface > encoded.size() - 6U - interior - partition ||
        faces != encoded.size() - 6U - interior - partition - interface)
      throw runtime::Error("immersed reconstruction snapshot is malformed");

    ImmersedReconstructionPlanSnapshot result{};
    result.owned_active_count = static_cast<std::size_t>(encoded[4]);
    result.fingerprint = encoded[5];
    std::size_t offset = 6U;
    const auto append = [&](std::size_t count, auto &target) {
      target.reserve(count);
      for (std::size_t index = 0U; index < count; ++index)
        target.push_back(encoded[offset++]);
    };
    append(interior, result.interior_rows);
    append(partition, result.partition_rows);
    append(interface, result.interface_rows);
    append(faces, result.interface_faces);
    return result;
  }

  static void reset_trace() {
    detail::immersed_reconstruction_reset_trace_raw();
  }

  static std::vector<ReconstructionExecutionStage> trace() {
    const auto encoded = detail::immersed_reconstruction_trace_raw();
    std::vector<ReconstructionExecutionStage> result;
    result.reserve(encoded.size());
    for (const auto stage : encoded)
      result.push_back(static_cast<ReconstructionExecutionStage>(stage));
    return result;
  }

  static std::size_t inactive_read_attempts() noexcept {
    return detail::immersed_reconstruction_inactive_read_attempts_raw();
  }

  static void fail_on_inactive_read(bool enabled) noexcept {
    detail::immersed_reconstruction_fail_on_inactive_read_raw(enabled);
  }

  static void force_momentum_fallback(bool enabled) noexcept {
    detail::immersed_reconstruction_force_momentum_fallback_raw(enabled);
  }
};

} // namespace hundun::finite_volume::test

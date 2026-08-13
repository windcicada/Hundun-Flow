// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_linear.hpp"

#include <cstddef>
#include <cstdint>

namespace hundun::v04::detail {

constexpr std::size_t kMgMaximumLevels = 32U;
constexpr std::uint8_t kMgAxisX = 1U;
constexpr std::uint8_t kMgAxisY = 2U;
constexpr std::uint8_t kMgAxisZ = 4U;

// All coefficient blocks live in one allocation owned by the plan.  Keeping
// offsets rather than pointers makes a numeric refresh independent of vector
// relocation and gives the hot V-cycle a compact, immutable level table.
struct MgLevelStorage {
  MgLevelView view{};
  MeshPatch patch{};
  std::uint8_t coarsen_mask{};
  std::size_t cells{};
  std::size_t x_faces{};
  std::size_t y_faces{};
  std::size_t z_faces{};
  std::size_t diagonal_offset{};
  std::size_t x_offset{};
  std::size_t y_offset{};
  std::size_t z_offset{};
  std::size_t volume_offset{};
};

inline double* block(double* base, std::size_t offset) noexcept {
  return base + offset;
}

inline const double* block(const double* base, std::size_t offset) noexcept {
  return base + offset;
}

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
void set_mg_runtime_counters_for_test(
    NativeCartesianMgPlan& plan, MgPlanCounters counters,
    RevisionToken generation) noexcept;
#endif

}  // namespace hundun::v04::detail

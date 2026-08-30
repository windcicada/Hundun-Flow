// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_linear.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

namespace hundun::v04::detail {

constexpr std::size_t kMgMaximumLevels = 32U;
constexpr std::uint8_t kMgAxisX = 1U;
constexpr std::uint8_t kMgAxisY = 2U;
constexpr std::uint8_t kMgAxisZ = 4U;

// The point operator has exactly one bounded mode for each local patch face.
// These six values are compiled from the immutable patch/boundary contract;
// they are deliberately not a per-cell map or a persistent pointer table.
enum class MgPointBoundaryMode : std::uint8_t {
  padded_ghost,
  self,
  zero,
};

// All coefficient blocks live in one allocation owned by the plan.  Keeping
// offsets rather than pointers makes a numeric refresh independent of vector
// relocation and gives the hot V-cycle a compact, immutable level table.
struct MgLevelStorage {
  MgLevelView view{};
  MeshPatch patch{};
  std::uint8_t coarsen_mask{};
  MgPointSmootherKind point_smoother{MgPointSmootherKind::red_black};
  double chebyshev_lower_spectrum_fraction{};
  std::size_t cells{};
  std::size_t x_faces{};
  std::size_t y_faces{};
  std::size_t z_faces{};
  std::size_t diagonal_offset{};
  std::size_t x_offset{};
  std::size_t y_offset{};
  std::size_t z_offset{};
  std::size_t volume_offset{};
  std::array<MgPointBoundaryMode, 6U> point_boundary_modes{
      MgPointBoundaryMode::padded_ghost, MgPointBoundaryMode::padded_ghost,
      MgPointBoundaryMode::padded_ghost, MgPointBoundaryMode::padded_ghost,
      MgPointBoundaryMode::padded_ghost, MgPointBoundaryMode::padded_ghost};
};

inline double* block(double* base, std::size_t offset) noexcept {
  return base + offset;
}

inline const double* block(const double* base, std::size_t offset) noexcept {
  return base + offset;
}

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
enum class MgCycleFailurePhase : std::uint8_t {
  none,
  pre_smooth,
  restriction,
  first_coarse,
  second_coarse,
  prolongation,
  post_smooth,
  terminal
};

struct MgMatrixWorkCounters {
  std::array<std::uint64_t, kMgMaximumLevels> cycle_level_calls{};
  std::uint64_t cycle_restrictions{};
  std::uint64_t cycle_prolongations{};
  std::uint64_t full_actions{};
  std::uint64_t full_cell_visits{};
  std::uint64_t residual_writes{};
  std::uint64_t retained_final_full_actions{};
  std::uint64_t retained_final_defect_actions{};
  std::uint64_t defect_cell_visits{};
  std::uint64_t defect_writes{};
  std::uint64_t residual_finish_actions{};
  std::uint64_t residual_finish_cell_visits{};
  std::uint64_t fused_color_actions{};
  std::uint64_t fused_cell_visits{};
  std::uint64_t fused_updates{};
  std::uint64_t separated_color_actions{};
  std::uint64_t separated_color_updates{};
  std::uint64_t chebyshev_stages{};
  std::uint64_t chebyshev_exchange_actions{};
  std::uint64_t chebyshev_defect_actions{};
  std::uint64_t chebyshev_retained_final_defect_actions{};
  // Chebyshev's streamed path performs the same cell stencil and recurrence
  // work as the tests-only two-pass oracle.  The additional fields make the
  // lifecycle differences explicit: intermediate defect publications are
  // eliminated, while an odd number of stages requires an interior copyback.
  std::uint64_t chebyshev_stencil_evaluations{};
  std::uint64_t chebyshev_updates{};
  std::uint64_t chebyshev_elided_intermediate_residual_publications{};
  std::uint64_t chebyshev_copyback_cells{};
};

void set_mg_runtime_counters_for_test(
    NativeCartesianMgPlan& plan, MgPlanCounters counters,
    RevisionToken generation) noexcept;

void set_mg_reference_point_actions_for_test(
    NativeCartesianMgPlan& plan, bool enabled) noexcept;

void set_mg_reference_chebyshev_lifecycle_for_test(
    NativeCartesianMgPlan& plan, bool enabled) noexcept;

void set_mg_pre_sweeps_for_test(NativeCartesianMgPlan& plan,
                                std::uint8_t sweeps) noexcept;

void set_mg_skip_final_projection_for_test(NativeCartesianMgPlan& plan,
                                           bool enabled) noexcept;

void set_mg_force_chebyshev_invalid_for_test(
    NativeCartesianMgPlan& plan, bool enabled) noexcept;

void set_mg_cycle_failure_for_test(
    NativeCartesianMgPlan& plan, MgCycleFailurePhase phase,
    std::size_t level, int failing_rank, bool deferred) noexcept;

// Mode: 0 = automatic policy, 1 = force the distributed reference route,
// 2 = request the replicated route (eligibility rules still apply).
void set_mg_replicated_coarse_mode_for_test(
    NativeCartesianMgPlan& plan, std::uint8_t mode) noexcept;

bool mg_replicated_coarse_enabled_for_test(
    const NativeCartesianMgPlan& plan) noexcept;

MgMatrixWorkCounters mg_matrix_work_counters_for_test(
    const NativeCartesianMgPlan& plan) noexcept;

double mg_level_diagonal_for_test(const NativeCartesianMgPlan& plan,
                                  std::size_t level,
                                  std::size_t cell) noexcept;

#endif

}  // namespace hundun::v04::detail

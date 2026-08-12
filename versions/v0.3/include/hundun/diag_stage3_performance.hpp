// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/diag_performance_artifact.hpp"

#include <array>
#include <cstdint>
#include <string_view>

namespace hundun::diagnostics {

struct Stage3PerformanceCounters final {
  std::uint64_t init_surface_triangles{};
  std::uint64_t init_query_closest_calls{};
  std::uint64_t init_query_segment_calls{};
  std::uint64_t init_classification_cells{};
  std::uint64_t init_ghost_qr_plans{};
  std::uint64_t init_ghost_rejected_plans{};
  std::uint64_t init_ghost_donor_references{};
  std::uint64_t init_wall_points{};
  std::uint64_t step_ghost_constraints{};
  std::uint64_t step_lfp_transforms{};
  std::uint64_t step_immersed_rows{};
  std::uint64_t step_pressure_wall_constraints{};
  std::uint64_t step_wall_quadrature_evaluations{};
  std::uint64_t step_force_reductions{};
  std::uint64_t step_wale_gradient_cells{};
  std::uint64_t step_wale_evaluations{};
  std::uint64_t checkpoint_logical_io_bytes{};
};

inline constexpr std::array<std::string_view, 17> kStage3PerformanceCounterIds{
    "checkpoint.logical-io-bytes",
    "init.classification.cells",
    "init.ghost.donor-references",
    "init.ghost.qr-plans",
    "init.ghost.rejected-plans",
    "init.query.closest-calls",
    "init.query.segment-calls",
    "init.surface.triangles",
    "init.wall.points",
    "step.force.reductions",
    "step.ghost.constraints",
    "step.immersed.rows",
    "step.lfp.transforms",
    "step.pressure.wall-constraints",
    "step.wale.evaluations",
    "step.wale.gradient-cells",
    "step.wall.quadrature-evaluations"};

Stage3PerformanceCounters
stage3_performance_counter_delta(const Stage3PerformanceCounters &before,
                                 const Stage3PerformanceCounters &after);
Stage3PerformanceCounters
add_stage3_performance_counters(const Stage3PerformanceCounters &left,
                                const Stage3PerformanceCounters &right);
CounterMap
stage3_algorithmic_work_map(const Stage3PerformanceCounters &counters);

} // namespace hundun::diagnostics

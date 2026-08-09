// SPDX-License-Identifier: Apache-2.0

#include "hundun/diag_stage3_performance.hpp"

#include <limits>
#include <stdexcept>

namespace hundun::diagnostics {
namespace {

std::uint64_t checked_add(std::uint64_t left, std::uint64_t right) {
  if (left > std::numeric_limits<std::uint64_t>::max() - right)
    throw std::overflow_error("Stage 3 performance counter would overflow");
  return left + right;
}

std::uint64_t checked_delta(std::uint64_t before, std::uint64_t after) {
  if (after < before)
    throw std::invalid_argument("Stage 3 performance counter decreased");
  return after - before;
}

#define HUNDUN_STAGE3_COUNTER_MEMBERS(X)                                     \
  X(init_surface_triangles)                                                   \
  X(init_query_closest_calls)                                                 \
  X(init_query_segment_calls)                                                 \
  X(init_classification_cells)                                                \
  X(init_ghost_qr_plans)                                                      \
  X(init_ghost_rejected_plans)                                                \
  X(init_ghost_donor_references)                                              \
  X(init_wall_points)                                                         \
  X(step_ghost_constraints)                                                   \
  X(step_lfp_transforms)                                                      \
  X(step_immersed_rows)                                                       \
  X(step_pressure_wall_constraints)                                           \
  X(step_wall_quadrature_evaluations)                                         \
  X(step_force_reductions)                                                    \
  X(step_wale_gradient_cells)                                                 \
  X(step_wale_evaluations)                                                    \
  X(checkpoint_logical_io_bytes)

} // namespace

Stage3PerformanceCounters
stage3_performance_counter_delta(const Stage3PerformanceCounters &before,
                                 const Stage3PerformanceCounters &after) {
  Stage3PerformanceCounters result;
#define HUNDUN_DELTA(name) result.name = checked_delta(before.name, after.name);
  HUNDUN_STAGE3_COUNTER_MEMBERS(HUNDUN_DELTA)
#undef HUNDUN_DELTA
  return result;
}

Stage3PerformanceCounters
add_stage3_performance_counters(const Stage3PerformanceCounters &left,
                                const Stage3PerformanceCounters &right) {
  Stage3PerformanceCounters result;
#define HUNDUN_ADD(name) result.name = checked_add(left.name, right.name);
  HUNDUN_STAGE3_COUNTER_MEMBERS(HUNDUN_ADD)
#undef HUNDUN_ADD
  return result;
}

CounterMap
stage3_algorithmic_work_map(const Stage3PerformanceCounters &counters) {
  return {{"checkpoint.logical-io-bytes",
           counters.checkpoint_logical_io_bytes},
          {"init.classification.cells", counters.init_classification_cells},
          {"init.ghost.donor-references",
           counters.init_ghost_donor_references},
          {"init.ghost.qr-plans", counters.init_ghost_qr_plans},
          {"init.ghost.rejected-plans", counters.init_ghost_rejected_plans},
          {"init.query.closest-calls", counters.init_query_closest_calls},
          {"init.query.segment-calls", counters.init_query_segment_calls},
          {"init.surface.triangles", counters.init_surface_triangles},
          {"init.wall.points", counters.init_wall_points},
          {"step.force.reductions", counters.step_force_reductions},
          {"step.ghost.constraints", counters.step_ghost_constraints},
          {"step.immersed.rows", counters.step_immersed_rows},
          {"step.lfp.transforms", counters.step_lfp_transforms},
          {"step.pressure.wall-constraints",
           counters.step_pressure_wall_constraints},
          {"step.wale.evaluations", counters.step_wale_evaluations},
          {"step.wale.gradient-cells", counters.step_wale_gradient_cells},
          {"step.wall.quadrature-evaluations",
           counters.step_wall_quadrature_evaluations}};
}

} // namespace hundun::diagnostics

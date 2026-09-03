// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/diag_stage3_performance.hpp"
#include "tests/support/test_main.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

using hundun::diagnostics::Stage3PerformanceCounters;

void check_inventory_and_mapping() {
  const auto &ids = hundun::diagnostics::kStage3PerformanceCounterIds;
  HUNDUN_CHECK(std::is_sorted(ids.begin(), ids.end()));
  HUNDUN_CHECK(std::adjacent_find(ids.begin(), ids.end()) == ids.end());
  Stage3PerformanceCounters counters;
  counters.init_surface_triangles = 1U;
  counters.init_query_closest_calls = 2U;
  counters.init_query_segment_calls = 3U;
  counters.init_classification_cells = 4U;
  counters.init_ghost_qr_plans = 5U;
  counters.init_ghost_rejected_plans = 6U;
  counters.init_ghost_donor_references = 7U;
  counters.init_wall_points = 8U;
  counters.step_ghost_constraints = 9U;
  counters.step_lfp_transforms = 10U;
  counters.step_immersed_rows = 11U;
  counters.step_pressure_wall_constraints = 12U;
  counters.step_wall_quadrature_evaluations = 13U;
  counters.step_force_reductions = 14U;
  counters.step_wale_gradient_cells = 15U;
  counters.step_wale_evaluations = 16U;
  counters.checkpoint_logical_io_bytes = 17U;
  const auto mapped =
      hundun::diagnostics::stage3_algorithmic_work_map(counters);
  HUNDUN_CHECK(mapped.size() == ids.size());
  for (const auto id : ids)
    HUNDUN_CHECK(mapped.count(std::string(id)) == 1U);
  const hundun::diagnostics::CounterMap expected{
      {"checkpoint.logical-io-bytes", 17U},
      {"init.classification.cells", 4U},
      {"init.ghost.donor-references", 7U},
      {"init.ghost.qr-plans", 5U},
      {"init.ghost.rejected-plans", 6U},
      {"init.query.closest-calls", 2U},
      {"init.query.segment-calls", 3U},
      {"init.surface.triangles", 1U},
      {"init.wall.points", 8U},
      {"step.force.reductions", 14U},
      {"step.ghost.constraints", 9U},
      {"step.immersed.rows", 11U},
      {"step.lfp.transforms", 10U},
      {"step.pressure.wall-constraints", 12U},
      {"step.wale.evaluations", 16U},
      {"step.wale.gradient-cells", 15U},
      {"step.wall.quadrature-evaluations", 13U}};
  HUNDUN_CHECK(mapped == expected);
}

void check_checked_arithmetic_and_retry_semantics() {
  Stage3PerformanceCounters before;
  before.step_ghost_constraints = 7U;
  Stage3PerformanceCounters after = before;
  after.step_ghost_constraints = 12U;
  const auto delta =
      hundun::diagnostics::stage3_performance_counter_delta(before, after);
  HUNDUN_CHECK(delta.step_ghost_constraints == 5U);

  // A failed attempt changes actual work even though its numerical state is
  // rolled back. Counter snapshots therefore remain monotonic.
  const auto retained = hundun::diagnostics::add_stage3_performance_counters(
      before, delta);
  HUNDUN_CHECK(retained.step_ghost_constraints ==
               after.step_ghost_constraints);

  bool rejected = false;
  try {
    static_cast<void>(
        hundun::diagnostics::stage3_performance_counter_delta(after, before));
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);

  Stage3PerformanceCounters maximum;
  maximum.step_wale_evaluations =
      std::numeric_limits<std::uint64_t>::max();
  Stage3PerformanceCounters one;
  one.step_wale_evaluations = 1U;
  rejected = false;
  try {
    static_cast<void>(hundun::diagnostics::add_stage3_performance_counters(
        maximum, one));
  } catch (const std::overflow_error &) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

} // namespace

int main() {
  return hundun::test::run([] {
    check_inventory_and_mapping();
    check_checked_arithmetic_and_retry_semantics();
  });
}

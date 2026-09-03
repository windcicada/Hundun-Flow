// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/cfg_resolved_case_v3.hpp"
#include "hundun/diag_stage3_performance.hpp"
#include "hundun/rt_types.hpp"
#include "tests/support/stage3_stl_fixture.hpp"

#include <vector>

namespace hundun::runtime {
class MpiContext;
}

namespace hundun::test {

runtime::Int3 stage3_performance_process_grid(int ranks);
config::FlowCaseConfig stage3_performance_case(int ranks, int cells);
std::vector<StlFixtureTriangle> stage3_performance_body();

struct Stage3PerformanceEvidence final {
  diagnostics::Stage3PerformanceCounters after_initialization;
  diagnostics::Stage3PerformanceCounters before_measured_step;
  diagnostics::Stage3PerformanceCounters after_measured_step;
  diagnostics::Stage3PerformanceCounters measured_delta;
  diagnostics::Stage3PerformanceCounters checkpoint_delta;
  diagnostics::Stage3PerformanceCounters failed_attempt_delta;
  std::uint64_t immersed_link_count{};
  std::uint64_t owned_active_cell_count{};
  std::uint64_t local_wall_point_count{};
  std::uint32_t pressure_corrector_count{};
  std::uint64_t warmup_steps{};
  std::uint64_t measured_steps{};
  double elapsed_seconds{};
  std::uint64_t surface_fingerprint{};
  std::uint64_t classification_fingerprint{};
  bool committed{};
  bool failed_attempt_rolled_back{};
};

struct Stage3PerformanceRun final {
  int cells{8};
  std::uint64_t warmup_steps{};
  std::uint64_t measured_steps{1U};
  int repetitions{1};
  bool inject_failed_attempt{true};
};

Stage3PerformanceEvidence
run_stage3_performance_evidence(const runtime::MpiContext &,
                                const Stage3PerformanceRun &);

} // namespace hundun::test

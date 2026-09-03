// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/diag_structured.hpp"
#include "hundun/exec_execution.hpp"
#include "hundun/flow_immersed.hpp"
#include "hundun/rt_halo_performance_counters.hpp"
#include "hundun/rt_mpi_context.hpp"

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

namespace hundun::flow {

class ImmersedFlowDiagnosticSource final {
public:
  const ImmersedFlowStepAttemptReport &report() const;
  int rank() const;
  std::uint64_t committed_step() const;
  double committed_time_s() const;
  double maximum_wall_penetration_m_per_s() const;
  double mean_wall_penetration_m_per_s() const;
  std::uint64_t classified_cell_count() const;
  std::uint64_t active_cell_count() const;
  std::uint64_t immersed_link_count() const;
  std::uint64_t donor_reference_count() const;
  std::uint64_t wall_quadrature_point_count() const;
  config::DensityModel density_model() const;
  bool wall_force_available() const;
  const immersed::WallForceSample &wall_force_sample() const;
  bool local_flow_pattern_available() const;
  std::uint64_t local_flow_pattern_algorithm_fingerprint() const;
  std::uint64_t local_flow_pattern_row_fingerprint() const;
  std::uint64_t local_flow_pattern_replacement_group_count() const;
  std::uint64_t local_flow_pattern_algebraic_occurrence_count() const;
  double local_flow_pattern_replacement_coefficient_l2() const;
  std::uint64_t local_flow_pattern_limiting_case_status() const;
  std::uint64_t snapshot_seal() const;
  runtime::Fp64ReductionCounters reduction_counters() const;
  runtime::HaloPerformanceCounters halo_counters() const;
  execution::AllocationCounters allocation_counters() const;

private:
  void validate() const;

  const std::uint64_t *attempt_generation_{};
  std::uint64_t expected_attempt_generation_{};
  ImmersedFlowStepAttemptReport report_;
  int rank_{};
  std::uint64_t step_{};
  double time_s_{};
  double maximum_wall_penetration_m_per_s_{};
  double mean_wall_penetration_m_per_s_{};
  std::uint64_t classified_cell_count_{};
  std::uint64_t active_cell_count_{};
  std::uint64_t immersed_link_count_{};
  std::uint64_t donor_reference_count_{};
  std::uint64_t wall_quadrature_point_count_{};
  config::DensityModel density_model_{config::DensityModel::constant};
  bool wall_force_available_{};
  immersed::WallForceSample wall_force_sample_;
  bool local_flow_pattern_available_{};
  std::uint64_t local_flow_pattern_algorithm_fingerprint_{};
  std::uint64_t local_flow_pattern_row_fingerprint_{};
  std::uint64_t local_flow_pattern_replacement_group_count_{};
  std::uint64_t local_flow_pattern_algebraic_occurrence_count_{};
  double local_flow_pattern_replacement_coefficient_l2_{};
  std::uint64_t local_flow_pattern_limiting_case_status_{};
  std::uint64_t snapshot_seal_{};
  runtime::Fp64ReductionCounters reduction_counters_{};
  runtime::HaloPerformanceCounters halo_counters_{};
  execution::AllocationCounters allocation_counters_{};

  friend class FixedStepImmersedFlow;
};

} // namespace hundun::flow

namespace hundun::diagnostics {

DiagnosticDescriptor
describe_diagnostics(const flow::ImmersedFlowDiagnosticSource &) noexcept;
std::vector<std::string_view> diagnostic_fingerprint_field_ids(
    const flow::ImmersedFlowDiagnosticSource &);
void collect_diagnostics(const flow::ImmersedFlowDiagnosticSource &,
                         const DiagnosticRequest &, DiagnosticSink &);
void collect_diagnostics(const flow::ImmersedFlowDiagnosticSource &,
                         const runtime::MpiContext &,
                         const DiagnosticRequest &, DiagnosticSink &);

DiagnosticDescriptor describe_diagnostics(
    const flow::ImmersedFlowDiagnosticSource &, DiagnosticModuleKind);
std::vector<std::string_view> diagnostic_fingerprint_field_ids(
    const flow::ImmersedFlowDiagnosticSource &, DiagnosticModuleKind);
void collect_diagnostics(const flow::ImmersedFlowDiagnosticSource &,
                         DiagnosticModuleKind, const DiagnosticRequest &,
                         DiagnosticSink &);
void collect_diagnostics(const flow::ImmersedFlowDiagnosticSource &,
                         DiagnosticModuleKind, const runtime::MpiContext &,
                         const DiagnosticRequest &, DiagnosticSink &);

} // namespace hundun::diagnostics

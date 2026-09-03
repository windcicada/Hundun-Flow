// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/flow_adaptive_time_control.hpp"
#include "hundun/rt_types.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace hundun::flow::detail {

bool time_control_config_valid(const config::FlowTimeConfig &) noexcept;
bool density_model_valid(config::DensityModel) noexcept;
bool solve_control_valid(const linear::SolveControl &) noexcept;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
double material_momentum_conservation_defect_raw(
    double, double, double, double, double, double, double, double, double,
    double, double, double, double, double, double, bool) noexcept;
void material_set_preflight_allocation_failure_rank_raw(int) noexcept;
void material_reset_preflight_allocation_failure_raw() noexcept;
void material_select_phase_failure_raw(
    const runtime::MpiContext &, bool, std::uint8_t, bool, std::uint8_t &,
    int &, bool &);
void material_require_reliable_collective_result_raw(std::uint8_t, int);
#endif

struct TransportDiffusivityAuthority final {
  std::uint64_t count{};
  std::uint64_t ordered_fingerprint{14695981039346656037ULL};
  double maximum{};
};

struct FacadeAssemblyIdentity final {
  bool live{};
  const void *facade{};
  const mesh::MeshTopology *topology{};
  const mesh::MeshGeometry *geometry{};
  const runtime::MpiContext *mpi{};
};

TransportDiffusivityAuthority
make_transport_diffusivity_authority(const std::vector<double> &coefficients);

struct TimeControlStateCodec final {
  static std::uint64_t seal(const config::FlowTimeConfig &,
                            config::DensityModel,
                            const TimeControlState &) noexcept;
  static bool semantically_valid(const config::FlowTimeConfig &,
                                 config::DensityModel,
                                 const AcceptedStepMetadata &,
                                 const TimeControlState &) noexcept;
};

std::string render_owned_cells(runtime::Box3);
std::string render_global_cells(runtime::Int3);
std::string render_owned_faces(std::size_t);
std::string render_global_faces(std::uint64_t);

struct AdaptiveTimeControlAccess final {
  static FacadeAssemblyIdentity
  assembly(const FixedStepConstantDensityFlow &) noexcept;
  static FacadeAssemblyIdentity
  assembly(const FixedStepMaterialDensityFlow &) noexcept;
  static FacadeAssemblyIdentity
  assembly(const FixedStepIdealGasFlow &) noexcept;
  static const TransportDiffusivityAuthority &
  authority(const FixedStepConstantDensityFlow &) noexcept;
  static const TransportDiffusivityAuthority &
  authority(const FixedStepMaterialDensityFlow &) noexcept;
  static const TransportDiffusivityAuthority &
  authority(const FixedStepIdealGasFlow &) noexcept;
  static const FlowState *state_identity(const FlowState &) noexcept;
  static bool state_live(const FlowState &) noexcept;
  static bool state_layout_matches(const FlowState &,
                                   const mesh::MeshTopology &) noexcept;
  static std::uint64_t diagnostic_identity(const FlowState &) noexcept;
  static std::vector<double> committed_density(const FlowState &);
  static std::vector<double> committed_face_mass_flux(const FlowState &);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  static void constant_workspace_values_raw(
      const FixedStepConstantDensityFlow &, std::size_t &, std::size_t &,
      std::uintptr_t &) noexcept;
  static void constant_cache_values_raw(
      const FixedStepConstantDensityFlow &,
      std::vector<std::uintptr_t> &, std::vector<std::size_t> &,
      std::array<std::uintptr_t, 3> &, std::array<std::uint64_t, 3> &,
      std::array<std::vector<double>, 3> &, std::size_t &, bool &);
  static void constant_pressure_operator_values_raw(
      const FixedStepConstantDensityFlow &, std::uintptr_t &,
      std::uint64_t &, bool &) noexcept;
  static void material_pressure_values_raw(
      const PisoCoupler &, std::vector<double> &, std::vector<double> &,
      std::vector<double> &, std::vector<double> &, double &,
      std::uint32_t &, bool &, bool &, std::uint64_t &);
  static void material_pressure_halo_values_raw(
      const PisoCoupler &, bool &, runtime::HaloPerformanceCounters &,
      bool &, runtime::HaloPerformanceCounters &);
  static void material_cache_values_raw(
      const FixedStepMaterialDensityFlow &, std::vector<std::uintptr_t> &,
      std::vector<std::size_t> &, std::array<std::uintptr_t, 3> &,
      std::array<std::uint64_t, 3> &,
      std::array<std::vector<double>, 3> &, std::size_t &, bool &);
  static void material_flow_pressure_values_raw(
      const FixedStepMaterialDensityFlow &, std::vector<double> &,
      std::vector<double> &, std::vector<double> &, std::vector<double> &,
      double &, std::uint32_t &, bool &, bool &, std::uint64_t &);
  static void material_flow_pressure_halo_values_raw(
      const FixedStepMaterialDensityFlow &, bool &,
      runtime::HaloPerformanceCounters &, bool &,
      runtime::HaloPerformanceCounters &);
  static const std::vector<double> &material_finalizer_flux_values_raw(
      const FixedStepMaterialDensityFlow &);
  static void material_force_attempt_identity_raw(
      FixedStepMaterialDensityFlow &, std::uint64_t) noexcept;
  static bool material_has_diagnostic_report_raw(
      const FixedStepMaterialDensityFlow &) noexcept;
  static void material_enable_vortex_source_raw(
      FixedStepMaterialDensityFlow &, bool) noexcept;
  static void material_change_transport_authority_raw(
      FixedStepMaterialDensityFlow &, std::uint8_t) noexcept;
  static void ideal_set_enthalpy_rate_raw(FixedStepIdealGasFlow &, double);
  static void ideal_exhaust_source_generation_raw(FixedStepIdealGasFlow &);
  static void ideal_force_finalization_identity_wrap_raw(
      FixedStepIdealGasFlow &);
  static std::vector<std::array<std::uint64_t, 3>>
  ideal_halo_trace_raw(const FixedStepIdealGasFlow &);
  static void ideal_cache_values_raw(
      const FixedStepIdealGasFlow &, bool,
      std::vector<std::uintptr_t> &, std::vector<std::size_t> &,
      std::array<std::uintptr_t, 3> &, std::array<std::uint64_t, 3> &,
      std::array<std::vector<double>, 3> &, std::size_t &, bool &);
  static void ideal_set_post_store_corruption_raw(FixedStepIdealGasFlow &,
                                                  int, bool);
  static void ideal_set_candidate_precedence_fault_raw(FixedStepIdealGasFlow &,
                                                       int);
  static void ideal_set_stage_failure_raw(FixedStepIdealGasFlow &,
                                          IdealGasClosureStage,
                                          IdealGasClosureFailureReason, int);
  static void ideal_set_outer_failure_raw(FixedStepIdealGasFlow &,
                                          std::uint8_t, int);
  static void ideal_set_prepare_fault_raw(FixedStepIdealGasFlow &,
                                          std::uint8_t, int);
  static void ideal_set_post_store_mpi_fault_raw(FixedStepIdealGasFlow &, int);
  static void ideal_set_post_assessment_fault_raw(FixedStepIdealGasFlow &,
                                                  std::uint8_t, int);
  static void ideal_set_attempt_layout_fault_raw(FixedStepIdealGasFlow &, int);
  static void ideal_set_outlet_backflow_fault_raw(FixedStepIdealGasFlow &);
  static void ideal_set_attempt_preparation_fault_raw(FixedStepIdealGasFlow &,
                                                      std::uint8_t, int);
  static void ideal_set_controlled_allocation_raw(FixedStepIdealGasFlow &, int);
  static bool ideal_allocation_observation_active_raw(
      const FixedStepIdealGasFlow &) noexcept;
#endif
};

struct TimeControlDiagnosticSnapshot final {
  TimeControlState state;
  std::array<TimeAttemptSummary, 9> attempts{};
  std::size_t attempt_count{};
  TimeAdvanceDisposition disposition{};
  StepFailureReason reason{};
  int lowest_failing_rank{-1};
  double accepted_dt_s{};
  double proposed_next_dt_s{};
  double convective_rate_per_s{};
  double diffusive_rate_per_s{};
  bool stability_metrics_available{};
  bool limited_by_min_dt{};
  std::uint8_t preflight_category{};
  config::FlowTimeConfig config{};
  config::DensityModel model{};
  std::uint64_t controller_identity{};
  std::uint64_t report_identity{};
  std::uint64_t flow_state_identity{};
  int relative_rank{};
  std::uint64_t observed_step{};
  double observed_time_s{};
  AcceptedStepMetadata observed_metadata{};
  runtime::Box3 local_box{};
  runtime::Int3 global_extent{};
  std::size_t canonical_owned_faces{};
  std::size_t local_faces{};
  std::uint64_t global_faces{};
  std::string local_cell_layout;
  std::string global_cell_layout;
  std::string local_face_layout;
  std::string global_face_layout;
};

} // namespace hundun::flow::detail

struct hundun::flow::TimeControlDiagnosticSource::Impl final {
  hundun::flow::detail::TimeControlDiagnosticSnapshot snapshot;
};

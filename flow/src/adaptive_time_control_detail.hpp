// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/flow/adaptive_time_control.hpp"
#include "hundun/runtime/types.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace hundun::flow::detail {

struct TransportDiffusivityAuthority final {
  std::uint64_t count{};
  std::uint64_t ordered_fingerprint{14695981039346656037ULL};
  double maximum{};
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
  static const TransportDiffusivityAuthority &
  authority(const FixedStepConstantDensityFlow &) noexcept;
  static const TransportDiffusivityAuthority &
  authority(const FixedStepMaterialDensityFlow &) noexcept;
  static const TransportDiffusivityAuthority &
  authority(const FixedStepIdealGasFlow &) noexcept;
  static const FlowState *state_identity(const FlowState &) noexcept;
  static std::uint64_t diagnostic_identity(const FlowState &) noexcept;
  static std::vector<double> committed_density(const FlowState &);
  static std::vector<double> committed_face_mass_flux(const FlowState &);
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

// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_ibm.hpp"
#include "hundun/v04_mesh.hpp"
#include "hundun/v04_spray.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

namespace hundun::v04::spray::detail {

inline constexpr std::size_t kMaximumStencilEntries = 8U;
inline constexpr std::size_t kMaximumTrajectorySegments = 256U;
inline constexpr std::size_t kMaximumReboundEvents = 32U;
inline constexpr std::size_t kMaximumResolvedSegments =
    kMaximumTrajectorySegments + 2U * kMaximumReboundEvents;

enum class ParcelGridStencilStatus : std::uint8_t {
  success,
  invalid_input,
  outside_domain,
  non_finite_output
};

struct ParcelGridStencilEntry {
  std::uint64_t global_cell{};
  Int3 global_index{};
  std::int32_t owner_rank{-1};
  double weight{};
};

struct ParcelGridCouplingStencil {
  bool available{};
  ParcelGridStencilStatus status{ParcelGridStencilStatus::invalid_input};
  std::uint8_t entry_count{};
  bool boundary_clamped{};
  std::array<ParcelGridStencilEntry, kMaximumStencilEntries> entries{};

  [[nodiscard]] bool succeeded() const noexcept {
    return available && status == ParcelGridStencilStatus::success;
  }
};

[[nodiscard]] ParcelGridCouplingStencil build_parcel_grid_coupling_stencil(
    const CartesianGeometryPlan& geometry, const MeshPatch& patch,
    Vector3 position_m) noexcept;

enum class TrajectoryCandidateKind : std::uint8_t { predictor, corrector };

struct ParcelKinematicSample {
  Vector3 position_m{};
  Vector3 velocity_m_per_s{};
  double elapsed_time_s{};
};

struct ParcelAccelerationSample {
  bool available{};
  Vector3 acceleration_m_per_s2{};
};

using StatelessParcelAcceleration = ParcelAccelerationSample (*)(
    const ParcelKinematicSample& sample, const void* immutable_context) noexcept;

struct ParcelTrajectoryInput {
  SprayParcelState parcel{};
  TrajectoryCandidateKind kind{TrajectoryCandidateKind::predictor};
  double duration_s{};
  // A non-negative value stops this candidate at the relative event time.
  // The default negative value means that no external event time is known.
  double event_stop_time_s{-1.0};
  double characteristic_length_m{};
  double maximum_particle_cfl{};
  std::uint32_t maximum_substeps{
      static_cast<std::uint32_t>(kMaximumTrajectorySegments)};
  StatelessParcelAcceleration acceleration{};
  const void* immutable_acceleration_context{};
};

enum class ParcelTrajectoryStatus : std::uint8_t {
  success,
  invalid_input,
  acceleration_unavailable,
  substep_limit,
  non_finite_output
};

struct ParcelTrajectorySegment {
  std::uint32_t ordinal{};
  double begin_time_s{};
  double end_time_s{};
  Vector3 begin_position_m{};
  Vector3 end_position_m{};
  Vector3 begin_velocity_m_per_s{};
  Vector3 end_velocity_m_per_s{};
};

struct ParcelTrajectoryCandidate {
  bool available{};
  ParcelTrajectoryStatus status{ParcelTrajectoryStatus::invalid_input};
  TrajectoryCandidateKind kind{TrajectoryCandidateKind::predictor};
  double requested_duration_s{};
  double advanced_duration_s{};
  SprayParcelState parcel{};
  std::uint32_t segment_count{};
  std::array<ParcelTrajectorySegment, kMaximumTrajectorySegments> segments{};

  [[nodiscard]] bool succeeded() const noexcept {
    return available && status == ParcelTrajectoryStatus::success;
  }
};

[[nodiscard]] ParcelTrajectoryCandidate make_parcel_trajectory_candidate(
    const ParcelTrajectoryInput& input) noexcept;

struct StaticIbmReboundConfig {
  ImmersedFluidSide fluid_side{ImmersedFluidSide::outside};
  double normal_restitution{1.0};
  double tangential_restitution{1.0};
  std::uint32_t maximum_events{1U};
};

enum class StaticIbmReboundStatus : std::uint8_t {
  success,
  invalid_input,
  surface_query_failure,
  event_limit,
  segment_limit,
  non_finite_output
};

struct StaticIbmReboundEvent {
  std::uint32_t ordinal{};
  std::uint32_t source_segment{};
  SurfaceTriangleId triangle{kInvalidSurfaceTriangle};
  double event_time_s{};
  Vector3 position_m{};
  Vector3 fluid_side_normal{};
  Vector3 incoming_velocity_m_per_s{};
  Vector3 outgoing_velocity_m_per_s{};
};

struct StaticIbmReboundCandidate {
  bool available{};
  StaticIbmReboundStatus status{StaticIbmReboundStatus::invalid_input};
  TrajectoryCandidateKind kind{TrajectoryCandidateKind::predictor};
  SprayParcelState parcel{};
  std::uint32_t event_count{};
  std::array<StaticIbmReboundEvent, kMaximumReboundEvents> events{};
  std::uint32_t segment_count{};
  std::array<ParcelTrajectorySegment, kMaximumResolvedSegments> segments{};

  [[nodiscard]] bool succeeded() const noexcept {
    return available && status == StaticIbmReboundStatus::success;
  }
};

[[nodiscard]] StaticIbmReboundCandidate resolve_static_ibm_rebounds(
    const ParcelTrajectoryCandidate& trajectory,
    const ImmersedSurfacePlan& surface,
    StaticIbmReboundConfig config) noexcept;

enum class PhysicalDomainFace : std::uint8_t {
  x_min,
  x_max,
  y_min,
  y_max,
  z_min,
  z_max
};

struct ParcelOutletLedger {
  bool recorded{};
  bool requires_parcel_removal{};
  bool permits_gas_source_deposition{};
  ParcelId parcel_id{};
  PhysicalDomainFace face{PhysicalDomainFace::x_min};
  double event_time_s{};
  Vector3 position_m{};
  Vector3 velocity_m_per_s{};
  double represented_liquid_mass_kg{};
  Vector3 represented_momentum_kg_m_per_s{};
};

enum class PhysicalDomainOutletStatus : std::uint8_t {
  success,
  invalid_input,
  non_finite_output
};

struct PhysicalDomainOutletCandidate {
  bool available{};
  PhysicalDomainOutletStatus status{
      PhysicalDomainOutletStatus::invalid_input};
  bool parcel_retained{};
  SprayParcelState parcel{};
  ParcelOutletLedger ledger{};

  [[nodiscard]] bool succeeded() const noexcept {
    return available && status == PhysicalDomainOutletStatus::success;
  }
};

[[nodiscard]] PhysicalDomainOutletCandidate resolve_physical_domain_outlet(
    const ParcelTrajectoryCandidate& trajectory,
    const CartesianGeometryPlan& geometry) noexcept;

}  // namespace hundun::v04::spray::detail

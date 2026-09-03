// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/diag_stage3_performance.hpp"

#include <cstdint>
#include <filesystem>

namespace hundun::boundary {
class BoundaryRegistry;
}

namespace hundun::config {
struct ImmersedFlowCaseConfig;
}

namespace hundun::immersed {
class GhostStencilPlan;
class ImmersedDomain;
class ImmersedSurface;
class LocalFlowPatternTransform;
class SurfaceQuery;
class WallQuadraturePlan;
}

namespace hundun::mesh {
class MeshGeometry;
class MeshTopology;
}

namespace hundun::les {
class WaleModel;
}

namespace hundun::runtime {
class MpiContext;
class StructuredDecomposition;
}

namespace hundun::flow {

class FixedStepImmersedFlow;
class FlowState;
class IdealGasClosure;
namespace detail {
struct CheckpointV3Access;
}

enum class CheckpointV3Operation : std::uint8_t { write, read };
enum class CheckpointV3Disposition : std::uint8_t { completed, failed };
enum class CheckpointV3FailureReason : std::uint8_t {
  none,
  invalid_input,
  state,
  layout,
  fingerprint,
  presence,
  file_integrity,
  filesystem,
  collective_operation
};
enum class CheckpointV3Phase : std::uint8_t {
  preflight,
  rank_payload,
  manifest,
  completed_marker,
  restore_prepare,
  restore_publish
};
enum class CheckpointV3Presence : std::uint8_t {
  constant_static_ibm = 1,
  constant_body_fitted_wale = 2,
  constant_static_ibm_wale = 3,
  material_static_ibm = 4,
  material_body_fitted_wale = 5,
  material_static_ibm_wale = 6,
  ideal_gas_static_ibm = 7,
  ideal_gas_body_fitted_wale = 8,
  ideal_gas_static_ibm_wale = 9
};
enum class CheckpointV3CheckStatus : std::uint8_t {
  not_checked,
  passed,
  failed
};

struct CheckpointV3ControlState final {
  double proposed_next_dt_s{};
  std::uint32_t last_retry_count{};
};

struct CheckpointV3WriteModules final {
  CheckpointV3Presence presence{CheckpointV3Presence::constant_static_ibm};
  const immersed::ImmersedSurface *surface{};
  const immersed::SurfaceQuery *query{};
  const immersed::ImmersedDomain *domain{};
  const immersed::GhostStencilPlan *ghost_plan{};
  const immersed::WallQuadraturePlan *wall_plan{};
  const immersed::LocalFlowPatternTransform *transform{};
  const les::WaleModel *wale{};
  const FixedStepImmersedFlow *flow{};
  const IdealGasClosure *ideal_gas{};
};

struct CheckpointV3ReadModules final {
  CheckpointV3Presence presence{CheckpointV3Presence::constant_static_ibm};
  const immersed::ImmersedSurface *surface{};
  const immersed::SurfaceQuery *query{};
  const immersed::ImmersedDomain *domain{};
  const immersed::GhostStencilPlan *ghost_plan{};
  const immersed::WallQuadraturePlan *wall_plan{};
  const immersed::LocalFlowPatternTransform *transform{};
  const les::WaleModel *wale{};
  FixedStepImmersedFlow *flow{};
  IdealGasClosure *ideal_gas{};
};

class CheckpointV3Report final {
public:
  CheckpointV3Operation operation() const noexcept;
  CheckpointV3Disposition disposition() const noexcept;
  CheckpointV3FailureReason reason() const noexcept;
  CheckpointV3Phase phase() const noexcept;
  CheckpointV3Presence presence() const noexcept;
  int lowest_failing_rank() const noexcept;
  std::uint64_t step() const noexcept;
  double time_s() const noexcept;
  std::uint64_t manifest_crc64() const noexcept;
  CheckpointV3CheckStatus crc_status() const noexcept;
  CheckpointV3CheckStatus fingerprint_status() const noexcept;
  CheckpointV3CheckStatus partition_status() const noexcept;
  CheckpointV3CheckStatus rollback_status() const noexcept;

private:
  CheckpointV3Operation operation_{CheckpointV3Operation::write};
  CheckpointV3Disposition disposition_{CheckpointV3Disposition::failed};
  CheckpointV3FailureReason reason_{CheckpointV3FailureReason::invalid_input};
  CheckpointV3Phase phase_{CheckpointV3Phase::preflight};
  CheckpointV3Presence presence_{CheckpointV3Presence::constant_static_ibm};
  int lowest_failing_rank_{-1};
  std::uint64_t step_{};
  double time_s_{};
  std::uint64_t manifest_crc64_{};
  CheckpointV3CheckStatus crc_status_{CheckpointV3CheckStatus::not_checked};
  CheckpointV3CheckStatus fingerprint_status_{
      CheckpointV3CheckStatus::not_checked};
  CheckpointV3CheckStatus partition_status_{
      CheckpointV3CheckStatus::not_checked};
  CheckpointV3CheckStatus rollback_status_{
      CheckpointV3CheckStatus::not_checked};

  friend struct detail::CheckpointV3Access;
};

class CheckpointV3ReadResult final {
public:
  const CheckpointV3Report &report() const noexcept;
  bool restored() const noexcept;
  const CheckpointV3ControlState &control_state() const;

private:
  CheckpointV3Report report_;
  CheckpointV3ControlState control_;
  bool restored_{};

  friend struct detail::CheckpointV3Access;
};

CheckpointV3Report write_checkpoint_v3(
    const runtime::MpiContext &, const runtime::StructuredDecomposition &,
    const mesh::MeshTopology &, const mesh::MeshGeometry &,
    const boundary::BoundaryRegistry &, const config::ImmersedFlowCaseConfig &,
    const CheckpointV3WriteModules &, const FlowState &,
    CheckpointV3ControlState, const std::filesystem::path &directory);

CheckpointV3ReadResult read_checkpoint_v3(
    const runtime::MpiContext &, const runtime::StructuredDecomposition &,
    const mesh::MeshTopology &, const mesh::MeshGeometry &,
    const boundary::BoundaryRegistry &, const config::ImmersedFlowCaseConfig &,
    const CheckpointV3ReadModules &, FlowState &,
    const std::filesystem::path &directory);

CheckpointV3Report write_checkpoint_v3(
    const runtime::MpiContext &, const runtime::StructuredDecomposition &,
    const mesh::MeshTopology &, const mesh::MeshGeometry &,
    const boundary::BoundaryRegistry &, const config::ImmersedFlowCaseConfig &,
    const immersed::ImmersedSurface &, const immersed::SurfaceQuery &,
    const immersed::ImmersedDomain &, const immersed::GhostStencilPlan &,
    const immersed::WallQuadraturePlan &,
    const immersed::LocalFlowPatternTransform &, const FlowState &,
    const FixedStepImmersedFlow &, CheckpointV3ControlState,
    const std::filesystem::path &directory);

CheckpointV3ReadResult read_checkpoint_v3(
    const runtime::MpiContext &, const runtime::StructuredDecomposition &,
    const mesh::MeshTopology &, const mesh::MeshGeometry &,
    const boundary::BoundaryRegistry &, const config::ImmersedFlowCaseConfig &,
    const immersed::ImmersedSurface &, const immersed::SurfaceQuery &,
    const immersed::ImmersedDomain &, const immersed::GhostStencilPlan &,
    const immersed::WallQuadraturePlan &,
    const immersed::LocalFlowPatternTransform &, FlowState &,
    FixedStepImmersedFlow &, const std::filesystem::path &directory);

diagnostics::Stage3PerformanceCounters
checkpoint_v3_performance_counters() noexcept;

} // namespace hundun::flow

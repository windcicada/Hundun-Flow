// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_case.hpp"
#include "hundun/v04_execution.hpp"
#include "hundun/v04_mesh.hpp"
#include "hundun/v04_parallel.hpp"

#include <mpi.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace hundun::v04 {

enum class BoundaryStage : std::uint8_t {
  thermo,
  momentum,
  pressure,
  enthalpy,
  scalar,
  diagnostic
};

enum class BoundaryRelation : std::uint8_t {
  dirichlet,
  zero_gradient,
  normal_gradient,
  reflect_normal,
  convective
};

enum class BoundaryValueSource : std::uint8_t {
  none,
  compiled_scalar,
  compiled_vector,
  resolved_scalar,
  resolved_vector,
  resolved_normal_gradient
};

enum class BoundaryFlowKernel : std::uint8_t {
  ordinary,
  mass_flow_constraint,
  static_state,
  total_state,
  pressure_backflow,
  characteristic
};

enum class BoundaryParameterRole : std::uint8_t {
  flow,
  thermal,
  transported_scalar
};

enum class BoundaryScalarRole : std::uint8_t { none, species, passive_scalar };

inline constexpr std::uint32_t kInvalidBoundaryParameter = UINT32_MAX;

struct BoundaryIndexSpan {
  BoundaryStage stage{BoundaryStage::thermo};
  BoundaryRelation relation{BoundaryRelation::zero_gradient};
  CartesianFace face{CartesianFace::x_min};
  FieldId field{};
  std::uint8_t component_begin{};
  std::uint8_t component_count{1U};
  std::uint8_t ghost_layers{1U};
  std::uint32_t tangent_inner_count{};
  std::uint32_t tangent_outer_count{};
  std::uint32_t parameter{};
  std::uint32_t resolved_begin{kInvalidBoundaryParameter};
  BoundaryValueSource value_source{BoundaryValueSource::none};
  // Number of tangential face cells owned by this resolved slice. Resolved
  // values are flattened as [outer][inner] starting at resolved_begin.
  // Compiled/none sources must keep this zero.
  std::uint32_t resolved_stride{};
};

struct BoundaryKernelBatch {
  BoundaryStage stage{BoundaryStage::thermo};
  BoundaryRelation relation{BoundaryRelation::zero_gradient};
  std::uint32_t span_begin{};
  std::uint32_t span_count{};
};

struct BoundaryFacePlan {
  BoundaryKind flow_kind{BoundaryKind::symmetry};
  BoundaryKind thermal_kind{BoundaryKind::none};
  BoundaryFlowKernel flow_kernel{BoundaryFlowKernel::ordinary};
  bool local_owner{};
  bool periodic{};
  std::uint32_t flow_parameter{kInvalidBoundaryParameter};
  std::uint32_t thermal_parameter{kInvalidBoundaryParameter};
  std::uint32_t scalar_begin{};
  std::uint32_t scalar_count{};
};

struct BoundaryTransportedField {
  FieldId field{};
  BoundaryScalarRole role{BoundaryScalarRole::passive_scalar};
};

struct NscbcPrimitive;
struct NscbcNormalGradient;
struct NscbcTarget;
struct NscbcWaves;
struct NscbcRates;

class SchemePlan {
 public:
  ConvectionScheme momentum() const noexcept { return momentum_; }
  ConvectionScheme enthalpy() const noexcept { return enthalpy_; }
  ConvectionScheme species() const noexcept { return species_; }
  ConvectionScheme passive_scalar() const noexcept { return passive_scalar_; }
  DiffusionScheme diffusion() const noexcept { return diffusion_; }
  double limiter() const noexcept { return limiter_; }
  std::uint8_t required_ghost_width() const noexcept {
    return required_ghost_width_;
  }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }

 private:
  friend class BoundaryCompiler;
  ConvectionScheme momentum_{ConvectionScheme::limited_central2};
  ConvectionScheme enthalpy_{ConvectionScheme::limited_central2};
  ConvectionScheme species_{ConvectionScheme::tvd2};
  ConvectionScheme passive_scalar_{ConvectionScheme::tvd2};
  DiffusionScheme diffusion_{DiffusionScheme::central2};
  double limiter_{1.0};
  std::uint8_t required_ghost_width_{2U};
  PlanFingerprint fingerprint_{};
};

class BoundaryPlan {
 public:
  BoundaryPlan() noexcept = default;
  BoundaryPlan(const BoundaryPlan&) = delete;
  BoundaryPlan& operator=(const BoundaryPlan&) = delete;
  BoundaryPlan(BoundaryPlan&&) noexcept = default;
  BoundaryPlan& operator=(BoundaryPlan&&) noexcept = default;

  Status face(CartesianFace selected,
              const BoundaryFacePlan*& out) const noexcept;
  Span<const BoundaryIndexSpan> spans() const noexcept {
    return {spans_.data(), spans_.size()};
  }
  Span<const BoundaryKernelBatch> batches() const noexcept {
    return {batches_.data(), batches_.size()};
  }
  Span<const BoundaryTransportedField> transported_fields() const noexcept {
    return {transported_fields_.data(), transported_fields_.size()};
  }
  Span<const double> velocity_x() const noexcept {
    return {velocity_x_.data(), velocity_x_.size()};
  }
  Span<const double> velocity_y() const noexcept {
    return {velocity_y_.data(), velocity_y_.size()};
  }
  Span<const double> velocity_z() const noexcept {
    return {velocity_z_.data(), velocity_z_.size()};
  }
  Span<const double> direction_x() const noexcept {
    return {direction_x_.data(), direction_x_.size()};
  }
  Span<const double> direction_y() const noexcept {
    return {direction_y_.data(), direction_y_.size()};
  }
  Span<const double> direction_z() const noexcept {
    return {direction_z_.data(), direction_z_.size()};
  }
  Span<const double> backflow_velocity_x() const noexcept {
    return {backflow_velocity_x_.data(), backflow_velocity_x_.size()};
  }
  Span<const double> backflow_velocity_y() const noexcept {
    return {backflow_velocity_y_.data(), backflow_velocity_y_.size()};
  }
  Span<const double> backflow_velocity_z() const noexcept {
    return {backflow_velocity_z_.data(), backflow_velocity_z_.size()};
  }
  Span<const double> scalar_targets() const noexcept {
    return {scalar_targets_.data(), scalar_targets_.size()};
  }
  Span<const double> scalar_backflow_targets() const noexcept {
    return {scalar_backflow_targets_.data(), scalar_backflow_targets_.size()};
  }
  Span<const double> pressure_targets() const noexcept {
    return {pressure_targets_.data(), pressure_targets_.size()};
  }
  Span<const double> temperature_targets() const noexcept {
    return {temperature_targets_.data(), temperature_targets_.size()};
  }
  Span<const double> total_pressure_targets() const noexcept {
    return {total_pressure_targets_.data(), total_pressure_targets_.size()};
  }
  Span<const double> total_temperature_targets() const noexcept {
    return {total_temperature_targets_.data(),
            total_temperature_targets_.size()};
  }
  Span<const double> backflow_temperature_targets() const noexcept {
    return {backflow_temperature_targets_.data(),
            backflow_temperature_targets_.size()};
  }
  Span<const double> heat_flux_targets() const noexcept {
    return {heat_flux_targets_.data(), heat_flux_targets_.size()};
  }
  Span<const double> mass_flow_targets() const noexcept {
    return {mass_flow_targets_.data(), mass_flow_targets_.size()};
  }
  Span<const double> relaxation_rates() const noexcept {
    return {relaxation_rates_.data(), relaxation_rates_.size()};
  }
  Span<const double> mach_limits() const noexcept {
    return {mach_limits_.data(), mach_limits_.size()};
  }
  Span<const std::uint8_t> allow_backflow() const noexcept {
    return {allow_backflow_.data(), allow_backflow_.size()};
  }
  Span<const BoundaryParameterRole> parameter_roles() const noexcept {
    return {parameter_roles_.data(), parameter_roles_.size()};
  }
  Span<const ScalarBoundaryKind> scalar_kinds() const noexcept {
    return {scalar_kinds_.data(), scalar_kinds_.size()};
  }
  Span<const BoundaryScalarRole> scalar_roles() const noexcept {
    return {scalar_roles_.data(), scalar_roles_.size()};
  }
  Span<const double> normal_distance_1() const noexcept {
    return {normal_distance_1_.data(), normal_distance_1_.size()};
  }
  Span<const double> normal_distance_2() const noexcept {
    return {normal_distance_2_.data(), normal_distance_2_.size()};
  }
  const HaloTopology& halo_topology() const noexcept { return halo_topology_; }
  FieldId velocity_field() const noexcept { return velocity_field_; }
  FieldId pressure_field() const noexcept { return pressure_field_; }
  FieldId enthalpy_field() const noexcept { return enthalpy_field_; }
  PressureReferenceKind pressure_reference() const noexcept {
    return pressure_reference_;
  }
  Int3 local_cells() const noexcept { return local_cells_; }
  std::size_t parameter_count() const noexcept { return scalar_targets_.size(); }
  std::size_t resolved_scalar_count() const noexcept {
    return resolved_scalar_count_;
  }
  std::size_t resolved_vector_count() const noexcept {
    return resolved_vector_count_;
  }
  std::size_t resolved_normal_gradient_count() const noexcept {
    return resolved_normal_gradient_count_;
  }
  std::uint8_t required_ghost_width() const noexcept {
    return required_ghost_width_;
  }
  RevisionToken revision() const noexcept { return revision_; }
  PlanFingerprint semantic_fingerprint() const noexcept {
    return semantic_fingerprint_;
  }
  PlanFingerprint local_layout_fingerprint() const noexcept {
    return local_layout_fingerprint_;
  }
  Status pressure_perturbation_target(CartesianFace face, double p_ref,
                                      double& out) const noexcept;
  Status nscbc_target(CartesianFace face, NscbcTarget& out) const noexcept;
  Status evaluate_nscbc(CartesianFace face,
                        const NscbcPrimitive& primitive,
                        const NscbcNormalGradient& gradient,
                        NscbcWaves& waves, NscbcRates& rates) const noexcept;

 private:
  friend class BoundaryCompiler;
  friend Status apply_boundary_ghosts(
      BoundaryStage, const BoundaryPlan&, Span<FieldView>,
      struct BoundaryResolvedValues) noexcept;

  std::array<BoundaryFacePlan, 6U> faces_{};
  std::vector<BoundaryIndexSpan> spans_;
  std::vector<BoundaryKernelBatch> batches_;
  std::vector<BoundaryTransportedField> transported_fields_;
  std::vector<double> velocity_x_;
  std::vector<double> velocity_y_;
  std::vector<double> velocity_z_;
  std::vector<double> direction_x_;
  std::vector<double> direction_y_;
  std::vector<double> direction_z_;
  std::vector<double> backflow_velocity_x_;
  std::vector<double> backflow_velocity_y_;
  std::vector<double> backflow_velocity_z_;
  std::vector<double> scalar_targets_;
  std::vector<double> scalar_backflow_targets_;
  std::vector<double> pressure_targets_;
  std::vector<double> temperature_targets_;
  std::vector<double> total_pressure_targets_;
  std::vector<double> total_temperature_targets_;
  std::vector<double> backflow_temperature_targets_;
  std::vector<double> heat_flux_targets_;
  std::vector<double> mass_flow_targets_;
  std::vector<double> relaxation_rates_;
  std::vector<double> mach_limits_;
  std::vector<std::uint8_t> allow_backflow_;
  std::vector<BoundaryParameterRole> parameter_roles_;
  std::vector<ScalarBoundaryKind> scalar_kinds_;
  std::vector<BoundaryScalarRole> scalar_roles_;
  std::vector<double> normal_distance_1_;
  std::vector<double> normal_distance_2_;
  HaloTopology halo_topology_{};
  FieldId velocity_field_{};
  FieldId pressure_field_{};
  FieldId enthalpy_field_{};
  PressureReferenceKind pressure_reference_{
      PressureReferenceKind::boundary_absolute};
  Int3 local_cells_{};
  std::size_t resolved_scalar_count_{};
  std::size_t resolved_vector_count_{};
  std::size_t resolved_normal_gradient_count_{};
  std::uint8_t required_ghost_width_{};
  RevisionToken revision_{};
  PlanFingerprint semantic_fingerprint_{};
  PlanFingerprint local_layout_fingerprint_{};
};

struct BoundaryResolvedValues {
  Span<const double> scalar;
  Span<const Real3> vector;
  Span<const double> normal_gradient;
};

Status apply_boundary_ghosts(BoundaryStage stage, const BoundaryPlan& plan,
                             Span<FieldView> fields,
                             BoundaryResolvedValues values) noexcept;

enum class TimeLimit : std::uint8_t {
  fixed,
  convective,
  viscous,
  thermal,
  species,
  acoustic,
  growth,
  retry
};
enum class StepOrigin : std::uint8_t {
  fresh_start,
  restart,
  accepted,
  retry
};

struct LocalTimeLimits {
  double convective{};
  double viscous{};
  double thermal{};
  double species{};
  double acoustic{};
};

struct BdfCoefficients {
  double a0{};
  double a1{};
  double a2{};
  std::uint8_t order{1U};
};

struct StepTime {
  double time{};
  double dt{};
  std::uint64_t accepted_step{};
  std::uint32_t attempt{};
  StepOrigin origin{StepOrigin::fresh_start};
  BdfCoefficients bdf{};
  std::uint64_t generation{};
};

class TimeSchemePlan {
 public:
  static Status compile(const TimeControlSpec& spec,
                        TimeSchemePlan& out) noexcept;
  Status local_candidate(LocalTimeLimits limits, double& dt,
                         TimeLimit& active_limit) const noexcept;
  const TimeControlSpec& spec() const noexcept { return spec_; }
  bool acoustic_hard_limit() const noexcept {
    return spec_.control == TimeControlKind::adaptive_acoustic;
  }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }

 private:
  TimeControlSpec spec_{};
  PlanFingerprint fingerprint_{};
};

class TimeControllerState {
 public:
  static Status start(const TimeSchemePlan& plan, double start_time,
                      TimeControllerState& out) noexcept;
  static Status restart(const TimeSchemePlan& plan, double time,
                        double last_accepted_dt,
                        std::uint64_t accepted_step,
                        TimeControllerState& out) noexcept;
  Status propose(MPI_Comm communicator, LocalTimeLimits limits,
                 StepTime& out) noexcept;
  Status finish(MPI_Comm communicator, const StepTime& step,
                Status local_outcome, StepTime& next) noexcept;

  double time() const noexcept { return time_; }
  double last_accepted_dt() const noexcept { return last_accepted_dt_; }
  std::uint64_t accepted_step() const noexcept { return accepted_step_; }
  std::uint32_t retry_count() const noexcept { return retry_count_; }
  int lowest_failing_rank() const noexcept { return lowest_failing_rank_; }
  bool has_active_proposal() const noexcept { return active_; }

 private:
  Status coefficients(double dt, bool force_backward_euler,
                      BdfCoefficients& out) const noexcept;
  bool matches_active(const StepTime& step) const noexcept;

  TimeSchemePlan plan_{};
  StepTime active_step_{};
  double time_{};
  double last_accepted_dt_{};
  std::uint64_t accepted_step_{};
  std::uint64_t next_generation_{1U};
  std::uint32_t retry_count_{};
  int lowest_failing_rank_{-1};
  StepOrigin next_origin_{StepOrigin::fresh_start};
  bool force_backward_euler_{true};
  bool active_{};
};

struct NscbcPrimitive {
  double density{};
  double pressure{};
  double temperature{};
  double normal_velocity{};
  double tangent_velocity_1{};
  double tangent_velocity_2{};
  double sound_speed{};
  double gamma{};
};

struct NscbcNormalGradient {
  double density{};
  double pressure{};
  double normal_velocity{};
  double tangent_velocity_1{};
  double tangent_velocity_2{};
};

struct NscbcTarget {
  double pressure{};
  double temperature{};
  double normal_velocity{};
  double tangent_velocity_1{};
  double tangent_velocity_2{};
  double pressure_relaxation{};
  double velocity_relaxation{};
  double temperature_relaxation{};
  double mach_limit{0.95};
  double reversal_tolerance{1.0e-8};
  bool allow_backflow{};
  double backflow_temperature{};
  double backflow_normal_velocity{};
  double backflow_tangent_velocity_1{};
  double backflow_tangent_velocity_2{};
};

struct NscbcWaves {
  double acoustic_incoming{};
  double entropy{};
  double tangent_1{};
  double tangent_2{};
  double acoustic_outgoing{};
};

struct NscbcRates {
  double density{};
  double pressure{};
  double temperature{};
  double normal_velocity{};
  double tangent_velocity_1{};
  double tangent_velocity_2{};
  bool used_backflow{};
};

Status evaluate_nscbc(BoundaryKind kind, const NscbcPrimitive& primitive,
                      const NscbcNormalGradient& gradient,
                      const NscbcTarget& target, NscbcWaves& waves,
                      NscbcRates& rates) noexcept;

class BoundaryCompiler {
 public:
  static Status compile(MPI_Comm communicator, const ValidatedModel& model,
                        const CartesianGeometryPlan& geometry,
                        const MeshPatch& patch, FieldRegistry& registry,
                        BoundaryPlan& boundary, SchemePlan& schemes,
                        TimeSchemePlan& time,
                        struct BoundaryCompileDiagnostics* diagnostics =
                            nullptr) noexcept;
};

struct BoundaryCompileDiagnostics {
  int lowest_failing_rank{-1};
};

}  // namespace hundun::v04

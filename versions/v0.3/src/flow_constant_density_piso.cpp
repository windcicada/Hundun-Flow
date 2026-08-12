// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow_constant_density_piso.hpp"

#include "flow_adaptive_time_control_detail.hpp"
#include "flow_body_fitted_wale_detail.hpp"
#include "flow_fixed_step_detail.hpp"

#include "hundun/fvm_matrix_free_poisson.hpp"
#include "hundun/fvm_poisson_boundary_adapter.hpp"
#include "hundun/les_wale.hpp"
#include "hundun/lin_ghosted_vector.hpp"
#include "hundun/rt_collective_status.hpp"
#include "hundun/rt_error.hpp"
#include "rt_mpi_error_detail.hpp"
#include <algorithm>
#include <array>
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#include <atomic>
#endif
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace hundun::flow {

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
namespace detail {
std::uint8_t material_terminal_reach_raw(std::uint8_t point) noexcept;
void material_face_flux_observation_begin_raw() noexcept;
void material_face_flux_observation_end_raw() noexcept;
}
#endif

namespace {

using mesh::EntityOwnership;
using mesh::FaceSide;
using mesh::LocalCellId;
using mesh::LocalFaceId;
using runtime::Int3;
using runtime::Real3;

constexpr runtime::PhaseId kScratchPhase = 1810U;
constexpr runtime::ActorId kScratchActor = 1810U;
constexpr runtime::PhaseId kStatePhase = 1800U;
constexpr runtime::ActorId kStateActor = 1800U;
constexpr double kContinuityTolerance = 1.0e-10;
constexpr double kFinalEquationTolerance = 1.0e-9;
constexpr double kFinalConservationTolerance = 5.0e-11;

std::uint64_t checked_halo_add(std::uint64_t left, std::uint64_t right) {
  if (right > std::numeric_limits<std::uint64_t>::max() - left)
    throw runtime::Error("pressure Halo counter would overflow");
  return left + right;
}

runtime::HaloPerformanceCounters combine_halo_counters(
    runtime::HaloPerformanceCounters left,
    const runtime::HaloPerformanceCounters &right) {
  left.completed_exchanges =
      checked_halo_add(left.completed_exchanges, right.completed_exchanges);
  left.begin_calls = checked_halo_add(left.begin_calls, right.begin_calls);
  left.wait_calls = checked_halo_add(left.wait_calls, right.wait_calls);
  left.send_payload_bytes =
      checked_halo_add(left.send_payload_bytes, right.send_payload_bytes);
  left.receive_payload_bytes = checked_halo_add(
      left.receive_payload_bytes, right.receive_payload_bytes);
  left.pack_bytes = checked_halo_add(left.pack_bytes, right.pack_bytes);
  left.unpack_bytes = checked_halo_add(left.unpack_bytes, right.unpack_bytes);
  left.send_messages =
      checked_halo_add(left.send_messages, right.send_messages);
  left.receive_messages =
      checked_halo_add(left.receive_messages, right.receive_messages);
  left.completed_wait_seconds += right.completed_wait_seconds;
  if (!std::isfinite(left.completed_wait_seconds) ||
      left.completed_wait_seconds < 0.0)
    throw runtime::Error("pressure Halo wait counter is invalid");
  return left;
}

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
constexpr std::uint8_t kMomentumMutationNone = 0U;
constexpr std::uint8_t kMomentumMutationOmitConvection = 1U;
constexpr std::uint8_t kMomentumMutationOmitViscosity = 2U;
constexpr std::uint8_t kMomentumMutationOmitPressure = 3U;
constexpr std::uint8_t kMomentumMutationOmitAlpha1 = 4U;
constexpr std::uint8_t kMomentumMutationOmitAlpha2 = 5U;
constexpr std::uint8_t kMomentumMutationReplaceHistoryFlux = 6U;
constexpr std::uint8_t kTransportMutationNone = 0U;
constexpr std::uint8_t kTransportMutationOmitHistorySpatial = 1U;
constexpr std::uint8_t kTransportMutationUseProvisionalFlux = 2U;
constexpr std::uint8_t kAttemptFailureNone = 0U;
constexpr std::uint8_t kAttemptFailureAfterBegin = 1U;
constexpr std::uint8_t kAttemptFailureAfterMomentum = 2U;
constexpr std::uint8_t kAttemptFailureAfterFacePredictor = 3U;
constexpr std::uint8_t kAttemptFailureAfterCorrectorOne = 4U;
constexpr std::uint8_t kAttemptFailureAfterProvisionalTransport = 5U;
constexpr std::uint8_t kAttemptFailureAfterCorrectorTwo = 6U;
constexpr std::uint8_t kAttemptFailureAfterFinalTransport = 7U;
constexpr std::uint8_t kAttemptFailureBeforeCommit = 8U;
constexpr std::uint8_t kTerminalModeReturnedRankless = 1U;
constexpr std::uint8_t kTerminalModeThrownOperation = 2U;
constexpr std::uint8_t kTerminalModeReturnedReliable = 3U;
constexpr std::uint8_t kTerminalFinalPressureEntry = 10U;
constexpr std::uint8_t kTerminalFinalPressureGammaSum = 11U;
constexpr std::uint8_t kTerminalFinalPressureGammaCount = 12U;
constexpr std::uint8_t kTerminalFinalPressureResidualReduction = 13U;

struct ConservationDiagnosticRaw final {
  double quantity_nm1{};
  double quantity_n{};
  double quantity_np1{};
  double signed_boundary_flux{};
  double absolute_boundary_flux{};
  double boundary_integral{};
  double history_correction{};
  double raw_defect{};
  double relative_defect{};
};

class FaceFluxPathObservation final {
public:
  explicit FaceFluxPathObservation(bool enabled = true) noexcept
      : enabled_(enabled) {
    if (enabled_)
      detail::material_face_flux_observation_begin_raw();
  }
  ~FaceFluxPathObservation() noexcept {
    if (enabled_)
      detail::material_face_flux_observation_end_raw();
  }
  FaceFluxPathObservation(const FaceFluxPathObservation &) = delete;
  FaceFluxPathObservation &operator=(const FaceFluxPathObservation &) = delete;

private:
  bool enabled_{};
};

std::atomic<bool> force_final_continuity_failure{false};
std::atomic<bool> force_final_pressure_failure{false};
std::atomic<bool> force_local_derived_failure{false};
std::atomic<std::size_t> final_momentum_perturb_component{
    std::numeric_limits<std::size_t>::max()};
std::atomic<double> final_momentum_perturb_delta{0.0};
std::atomic<std::size_t> final_transport_perturb_index{
    std::numeric_limits<std::size_t>::max()};
std::atomic<double> final_transport_perturb_delta{0.0};
std::atomic<bool> force_final_conservation_failure{false};
std::atomic<double> final_mass_defect_perturbation{0.0};
std::array<std::atomic<bool>, 3> final_momentum_norm_armed{};
std::array<std::atomic<double>, 3> final_momentum_norm_residual_square{};
std::array<std::atomic<double>, 3> final_momentum_norm_scale_square{};
std::atomic<std::size_t> final_transport_norm_index{
    std::numeric_limits<std::size_t>::max()};
std::atomic<double> final_transport_norm_residual_square{0.0};
std::atomic<double> final_transport_norm_scale_square{0.0};
std::array<std::atomic<bool>, 3> momentum_conservation_parts_armed{};
std::array<std::array<std::atomic<double>, 8>, 3>
    momentum_conservation_parts_values{};
std::atomic<std::size_t> momentum_conservation_overflow_component{
    std::numeric_limits<std::size_t>::max()};
std::atomic<std::size_t> transport_conservation_overflow_index{
    std::numeric_limits<std::size_t>::max()};
ConservationDiagnosticRaw last_mass_conservation{};
std::array<ConservationDiagnosticRaw, 3> last_momentum_conservation{};
std::vector<ConservationDiagnosticRaw> last_transport_conservation;
std::atomic<std::uint8_t> momentum_assembly_mutation{kMomentumMutationNone};
std::atomic<std::uint8_t> transport_assembly_mutation{kTransportMutationNone};
std::atomic<std::uint8_t> attempt_failure_stage{kAttemptFailureNone};
std::atomic<int> last_pressure_constraint_mode{-1};
std::atomic<int> last_pressure_operator_mode{-1};
std::atomic<bool> provisional_transport_sentinel{false};
std::atomic<double> final_uniform_x_mass_flux{0.0};
std::atomic<bool> final_uniform_x_mass_flux_override{false};
std::array<std::atomic<double>, 3> last_momentum_rhs{};
std::array<std::atomic<double>, 3> last_momentum_diagonal{};
std::atomic<std::size_t> provisional_transport_call_count{0U};
std::atomic<std::size_t> final_transport_call_count{0U};
std::atomic<int> pressure_operator_construction_failure_rank{-1};
std::atomic<int> pressure_operator_refresh_failure_rank{-1};
#endif

bool finite(Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

Real3 add(Real3 left, Real3 right) noexcept {
  return {left.x + right.x, left.y + right.y, left.z + right.z};
}

Real3 subtract(Real3 left, Real3 right) noexcept {
  return {left.x - right.x, left.y - right.y, left.z - right.z};
}

Real3 multiply(double scale, Real3 value) noexcept {
  return {scale * value.x, scale * value.y, scale * value.z};
}

double dot(Real3 left, Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

double component(Real3 value, int index) noexcept {
  return index == 0 ? value.x : index == 1 ? value.y : value.z;
}

bool solve_success(linear::SolveTerminationReason reason) noexcept {
  return reason == linear::SolveTerminationReason::converged ||
         reason == linear::SolveTerminationReason::zero_right_hand_side;
}

std::size_t bytes_for(std::size_t count) {
  if (count > std::numeric_limits<std::size_t>::max() / sizeof(double)) {
    throw runtime::Error("vector byte count overflows");
  }
  return count * sizeof(double);
}

std::size_t multiplied_count(std::size_t count, std::size_t factor,
                             const char *subject) {
  if (factor != 0U &&
      count > std::numeric_limits<std::size_t>::max() / factor) {
    throw runtime::Error(std::string(subject) + " count overflows");
  }
  return count * factor;
}

struct StructuredIndex final {
  int i{};
  int j{};
  int k{};
};

StructuredIndex map_cell(Int3 global, runtime::Box3 owned, Int3 global_extent) {
  const Int3 local{owned.end.x - owned.begin.x, owned.end.y - owned.begin.y,
                   owned.end.z - owned.begin.z};
  const auto axis = [](int coordinate, int begin, int end, int global_n,
                       int local_n) {
    if (coordinate >= begin && coordinate < end)
      return coordinate - begin;
    if (coordinate == begin - 1 || (begin == 0 && coordinate == global_n - 1)) {
      return -1;
    }
    if (coordinate == end || (end == global_n && coordinate == 0)) {
      return local_n;
    }
    throw runtime::Error("cell has no structured field mapping");
  };
  return {axis(global.x, owned.begin.x, owned.end.x, global_extent.x, local.x),
          axis(global.y, owned.begin.y, owned.end.y, global_extent.y, local.y),
          axis(global.z, owned.begin.z, owned.end.z, global_extent.z, local.z)};
}

template <class T>
T &at(const runtime::FieldView<T> &view, StructuredIndex index,
      int field_component) {
  return view(index.i, index.j, index.k, field_component);
}

runtime::FieldDescriptor
cell_scratch(std::string name, std::uint32_t components, int ghost_width) {
  return {std::move(name),
          "1",
          "fixed_step_constant_density_flow",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          ghost_width,
          false,
          runtime::RestartPolicy::transient,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor face_scratch(std::string name,
                                      std::uint32_t components) {
  return {std::move(name),
          "1",
          "fixed_step_constant_density_flow",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::transient,
          runtime::OutputPolicy::never};
}

class ScratchFields final {
public:
  explicit ScratchFields(runtime::FieldLayoutSet layout) {
    lagged_velocity =
        registry.declare_field(cell_scratch("lagged_velocity", 3U, 2));
    velocity_gradient =
        registry.declare_field(cell_scratch("velocity_gradient", 9U, 2));
    pressure_gradient =
        registry.declare_field(cell_scratch("pressure_gradient", 3U, 2));
    momentum_face = registry.declare_field(face_scratch("momentum_face", 3U));
    momentum_residual =
        registry.declare_field(cell_scratch("momentum_residual", 3U, 0));
    actual_diagonal =
        registry.declare_field(cell_scratch("actual_diagonal", 3U, 2));
    scalar_gradient =
        registry.declare_field(cell_scratch("scalar_gradient", 3U, 2));
    scalar_face = registry.declare_field(face_scratch("scalar_face", 1U));
    scalar_gamma = registry.declare_field(face_scratch("scalar_gamma", 1U));
    sgs_dynamic_viscosity_cell = registry.declare_field(
        cell_scratch("sgs_dynamic_viscosity_cell", 1U, 2));
    sgs_dynamic_viscosity =
        registry.declare_field(face_scratch("sgs_dynamic_viscosity", 1U));
    effective_dynamic_viscosity = registry.declare_field(
        face_scratch("effective_dynamic_viscosity", 1U));
    scalar_residual =
        registry.declare_field(cell_scratch("scalar_residual", 1U, 0));
    pressure_correction =
        registry.declare_field(cell_scratch("pressure_correction", 1U, 2));
    mass_residual =
        registry.declare_field(cell_scratch("mass_residual", 1U, 0));
    registry.freeze();
    access = std::make_unique<runtime::FieldAccessPlan>(registry);
    for (runtime::FieldId field = 0U;
         field < static_cast<runtime::FieldId>(registry.size()); ++field) {
      access->declare_access(kScratchPhase, kScratchActor, field,
                             runtime::AccessMode::read_write);
    }
    access->freeze();
    storage = std::make_unique<runtime::FieldStorage>(registry, layout);
  }

  runtime::FieldRegistry registry;
  std::unique_ptr<runtime::FieldAccessPlan> access;
  std::unique_ptr<runtime::FieldStorage> storage;
  runtime::FieldId lagged_velocity{};
  runtime::FieldId velocity_gradient{};
  runtime::FieldId pressure_gradient{};
  runtime::FieldId momentum_face{};
  runtime::FieldId momentum_residual{};
  runtime::FieldId actual_diagonal{};
  runtime::FieldId scalar_gradient{};
  runtime::FieldId scalar_face{};
  runtime::FieldId scalar_gamma{};
  runtime::FieldId sgs_dynamic_viscosity_cell{};
  runtime::FieldId sgs_dynamic_viscosity{};
  runtime::FieldId effective_dynamic_viscosity{};
  runtime::FieldId scalar_residual{};
  runtime::FieldId pressure_correction{};
  runtime::FieldId mass_residual{};
};

class DiagonalMomentumOperator final : public linear::LinearOperator {
public:
  DiagonalMomentumOperator(execution::ExecutionContext &context,
                           linear::VectorLayout layout)
      : context_(&context), layout_(std::move(layout)),
        diagonal_(layout_.owned_count(), 1.0) {}

  void replace(const std::vector<double> &diagonal) {
    if (diagonal.size() != diagonal_.size() ||
        !std::all_of(diagonal.begin(), diagonal.end(), [](double value) {
          return value > 0.0 && std::isfinite(value);
        })) {
      throw runtime::Error("momentum diagonal is invalid");
    }
    if (revision_ == std::numeric_limits<std::uint64_t>::max()) {
      throw runtime::Error("momentum revision would wrap");
    }
    std::copy(diagonal.begin(), diagonal.end(), diagonal_.begin());
    ++revision_;
  }

  linear::VectorLayout domain_layout() const override { return layout_; }
  linear::VectorLayout range_layout() const override { return layout_; }
  const execution::ExecutionContext &context() const override {
    return *context_;
  }
  std::uint64_t revision() const override { return revision_; }
  execution::ExecutionEvent
  apply(execution::VectorView<const double> x,
        execution::VectorView<double> y) const override {
    if (x.size() != diagonal_.size() || y.size() != diagonal_.size() ||
        x.stride() != 1U || y.stride() != 1U || !y.writable() ||
        x.backend_identity() != context_->backend_identity() ||
        y.backend_identity() != context_->backend_identity()) {
      throw runtime::Error("momentum operator view is incompatible");
    }
    for (std::size_t cell = 0; cell < diagonal_.size(); ++cell) {
      y[cell] = diagonal_[cell] * x[cell];
    }
    return execution::ExecutionEvent::completed();
  }
  bool has_diagonal() const override { return true; }
  execution::ExecutionEvent
  diagonal(execution::VectorView<double> output) const override {
    if (output.size() != diagonal_.size() || output.stride() != 1U ||
        !output.writable() ||
        output.backend_identity() != context_->backend_identity()) {
      throw runtime::Error("momentum diagonal view is incompatible");
    }
    for (std::size_t cell = 0; cell < diagonal_.size(); ++cell) {
      output[cell] = diagonal_[cell];
    }
    return execution::ExecutionEvent::completed();
  }

private:
  execution::ExecutionContext *context_;
  linear::VectorLayout layout_;
  std::vector<double> diagonal_;
  std::uint64_t revision_{1U};
};

void zero_cell(runtime::FieldView<double> view) {
  const Int3 extent = view.interior_extent();
  for (int k = 0; k < extent.z; ++k) {
    for (int j = 0; j < extent.y; ++j) {
      for (int i = 0; i < extent.x; ++i) {
        for (std::uint32_t component_index = 0;
             component_index < view.components(); ++component_index) {
          view(i, j, k, static_cast<int>(component_index)) = 0.0;
        }
      }
    }
  }
}

double global_l2(const runtime::MpiContext &mpi,
                 execution::VectorView<const double> values) {
  double square = 0.0;
  for (std::size_t index = 0; index < values.size(); ++index) {
    square += values[index] * values[index];
  }
  mpi.allreduce_fp64_in_place(&square, 1U,
                              runtime::Fp64ReductionOperation::sum);
  return std::sqrt(square);
}

void validate_host_context(const execution::ExecutionContext &context) {
  if (context.backend_identity() == 0U ||
      context.space() != execution::ExecutionSpace::host ||
      !context.supports(execution::ExecutionCapability::host_access) ||
      !context.supports(execution::ExecutionCapability::buffer_allocation)) {
    throw runtime::Error("requires a host execution context");
  }
}

double face_factor(const mesh::MeshGeometry &geometry, LocalFaceId face) {
  const Real3 area = geometry.face_area_vector_m2(face, FaceSide::owner);
  const Real3 displacement = geometry.face_displacement_m(face);
  const double projection = dot(area, displacement);
  const double factor = dot(area, area) / projection;
  if (!finite(area) || !finite(displacement) || !(projection > 0.0) ||
      !(factor > 0.0) || !std::isfinite(factor)) {
    throw runtime::Error("pressure face metric is invalid");
  }
  return factor;
}

template <class PressureBoundaryValue>
void compute_pressure_gradient(const mesh::MeshTopology &topology,
                               const mesh::MeshGeometry &geometry,
                               const runtime::FieldView<const double> &pressure,
                               const runtime::FieldView<double> &gradient,
                               std::vector<Real3> &sums,
                               PressureBoundaryValue &&boundary_value) {
  if (sums.size() != topology.owned_cell_count()) {
    throw runtime::Error("pressure-gradient workspace is invalid");
  }
  std::fill(sums.begin(), sums.end(), Real3{});
  const auto owned = topology.owned_global_box();
  const Int3 global_extent = topology.global_extent();
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    const LocalCellId owner = topology.owner(face);
    const auto neighbour = topology.neighbour(face);
    const StructuredIndex owner_index =
        map_cell(topology.global_cell(owner), owned, global_extent);
    const double owner_value = at(pressure, owner_index, 0);
    double face_value = owner_value;
    if (neighbour.has_value()) {
      const StructuredIndex neighbour_index =
          map_cell(topology.global_cell(*neighbour), owned, global_extent);
      const double neighbour_value = at(pressure, neighbour_index, 0);
      face_value = owner_value == neighbour_value
                       ? owner_value
                       : 0.5 * (owner_value + neighbour_value);
    } else {
      face_value = boundary_value(face, owner_value);
    }
    if (!std::isfinite(face_value)) {
      throw runtime::Error("pressure face value is non-finite");
    }
    const Real3 owner_area =
        geometry.face_area_vector_m2(face, FaceSide::owner);
    if (topology.cell_ownership(owner) == EntityOwnership::owned) {
      sums[owner] = add(sums[owner], multiply(face_value, owner_area));
    }
    if (!topology.patch_id(face).has_value() && neighbour.has_value() &&
        topology.cell_ownership(*neighbour) == EntityOwnership::owned) {
      sums[*neighbour] =
          add(sums[*neighbour], multiply(-face_value, owner_area));
    }
  }
  for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const Real3 value =
        multiply(1.0 / geometry.cell_volume_m3(cell), sums[cell]);
    if (!finite(value)) {
      throw runtime::Error("pressure gradient is non-finite");
    }
    const Int3 index = topology.global_cell(cell);
    const StructuredIndex local = map_cell(index, owned, global_extent);
    gradient(local.i, local.j, local.k, 0) = value.x;
    gradient(local.i, local.j, local.k, 1) = value.y;
    gradient(local.i, local.j, local.k, 2) = value.z;
  }
}

StepAttemptReport base_report(double dt_s) {
  StepAttemptReport report;
  report.attempted_dt_s = dt_s;
  report.final_continuity_normalized_l2 =
      std::numeric_limits<double>::infinity();
  report.final_pressure_residual_l2 = std::numeric_limits<double>::infinity();
  report.final_momentum_normalized_l2.fill(
      std::numeric_limits<double>::infinity());
  report.final_mass_relative_conservation_defect =
      std::numeric_limits<double>::infinity();
  report.final_momentum_relative_conservation_defect.fill(
      std::numeric_limits<double>::infinity());
  return report;
}

StepAttemptReport numerical_failure(StepAttemptReport report,
                                    StepFailureReason reason, int rank) {
  report.disposition = StepAttemptDisposition::recoverable_failure;
  report.reason = reason;
  report.lowest_failing_rank = rank;
  report.suggested_dt_s = 0.5 * report.attempted_dt_s;
  return report;
}

StepAttemptReport fatal_failure(StepAttemptReport report,
                                StepFailureReason reason, int rank) {
  report.disposition = StepAttemptDisposition::non_retryable_failure;
  report.reason = reason;
  report.lowest_failing_rank = rank;
  report.suggested_dt_s = 0.0;
  return report;
}

PressureCorrectionReport correction_failure(
    PressureCorrectionReport report, PressureCorrectionDisposition disposition,
    StepFailureReason reason, int rank) noexcept {
  report.disposition = disposition;
  report.reason = reason;
  report.lowest_failing_rank = rank;
  report.accepted = false;
  return report;
}

double low_u32(std::uint64_t value) noexcept {
  return static_cast<double>(value & UINT64_C(0xffffffff));
}

double high_u32(std::uint64_t value) noexcept {
  return static_cast<double>(value >> 32U);
}

std::uint64_t fp64_bits(double value) noexcept {
  std::uint64_t bits{};
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

bool valid_solve_control(const linear::SolveControl &control) noexcept {
  return std::isfinite(control.atol) && control.atol >= 0.0 &&
         std::isfinite(control.rtol) && control.rtol >= 0.0 &&
         control.residual_recompute_interval != 0U;
}

runtime::CollectiveStatus
agree_solve_control(const runtime::MpiContext &mpi,
                    const linear::SolveControl &control,
                    const char *invalid_message,
                    const char *mismatch_message) {
  const auto validity =
      runtime::collective_status(mpi, valid_solve_control(control),
                                 invalid_message);
  if (!validity.ok)
    return validity;

  const std::array<std::uint64_t, 4> local_words{
      fp64_bits(control.atol), fp64_bits(control.rtol),
      control.max_iterations, control.residual_recompute_interval};
  std::array<double, 8> root_halves{};
  if (mpi.rank() == 0) {
    for (std::size_t word = 0; word < local_words.size(); ++word) {
      root_halves[word * 2U] = low_u32(local_words[word]);
      root_halves[word * 2U + 1U] = high_u32(local_words[word]);
    }
  }
  mpi.allreduce_fp64_in_place(root_halves.data(), root_halves.size(),
                              runtime::Fp64ReductionOperation::sum);
  bool matches = true;
  for (std::size_t word = 0; word < local_words.size(); ++word) {
    matches = matches && root_halves[word * 2U] == low_u32(local_words[word]) &&
              root_halves[word * 2U + 1U] == high_u32(local_words[word]);
  }
  return runtime::collective_status(mpi, matches, mismatch_message);
}

template <std::size_t Count>
runtime::CollectiveStatus
agree_fp64_inputs(const runtime::MpiContext &mpi,
                  const std::array<double, Count> &local, const char *message) {
  auto negated_minimum = local;
  for (double &value : negated_minimum)
    value = -value;
  mpi.allreduce_fp64_in_place(negated_minimum.data(), negated_minimum.size(),
                              runtime::Fp64ReductionOperation::maximum);
  bool agrees = true;
  for (std::size_t index = 0; index < local.size(); ++index) {
    agrees = agrees && local[index] == -negated_minimum[index];
  }
  return runtime::collective_status(mpi, agrees, message);
}

struct SynchronizedAttemptFailure final {
  StepFailureReason reason;
  int failing_rank;
  bool recoverable;
};

template <class Function>
void synchronized_local_phase(const runtime::MpiContext &mpi,
                              StepFailureReason reason, bool recoverable,
                              const char *message, Function &&function) {
  bool local_ok = true;
  try {
    std::forward<Function>(function)();
  } catch (const runtime::detail::MpiOperationError &) {
    throw;
  } catch (...) {
    local_ok = false;
  }
  const auto status = runtime::collective_status(mpi, local_ok, message);
  if (!status.ok) {
    throw SynchronizedAttemptFailure{reason, status.failing_rank, recoverable};
  }
}

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
void inject_attempt_stage_failure(const runtime::MpiContext &mpi,
                                  std::uint8_t stage) {
  synchronized_local_phase(
      mpi, StepFailureReason::non_finite_trial, true,
      "injected transaction-stage failure", [&] {
        if (attempt_failure_stage.load(std::memory_order_relaxed) == stage)
          throw runtime::Error("injected transaction-stage failure");
      });
}
#endif

} // namespace

struct PisoCoupler::Impl final {
  Impl(const runtime::StructuredDecomposition &supplied_decomposition,
       const mesh::MeshTopology &supplied_topology,
       const mesh::MeshGeometry &supplied_geometry,
       const boundary::BoundaryRegistry &supplied_boundaries,
       const runtime::MpiContext &supplied_mpi,
       execution::ExecutionContext &supplied_execution,
       runtime::HaloExchange &supplied_halo,
       const linear::LinearSolver &supplied_solver,
       linear::Preconditioner &supplied_preconditioner)
      : decomposition(&supplied_decomposition), topology(&supplied_topology),
        geometry(&supplied_geometry), boundaries(&supplied_boundaries),
        mpi(&supplied_mpi), execution(&supplied_execution),
        halo(&supplied_halo), solver(&supplied_solver),
        preconditioner(&supplied_preconditioner),
        fvm(finite_volume::CellCenteredFvmOperators::create(supplied_topology,
                                                            supplied_geometry)),
        prepared_face_flux(
            finite_volume::FaceMassFlux::prepare(supplied_topology)),
        scratch({supplied_decomposition.local_extent(),
                 supplied_topology.local_face_count()}),
        gamma(supplied_execution,
              bytes_for(supplied_topology.local_face_count())),
        rhs(supplied_execution,
            bytes_for(supplied_topology.owned_cell_count())),
        correction(supplied_execution,
                   bytes_for(supplied_topology.owned_cell_count())),
        residual(supplied_execution,
                 bytes_for(supplied_topology.owned_cell_count())),
        candidate_pressure(supplied_execution,
                           bytes_for(supplied_topology.owned_cell_count())),
        pressure_gradient_sums(supplied_topology.owned_cell_count()),
        velocity_candidate(multiplied_count(
            supplied_topology.owned_cell_count(), 3U,
            "corrected velocity workspace")),
        pressure_candidate(supplied_topology.owned_cell_count()),
        face_velocity_candidate(multiplied_count(
            supplied_topology.local_face_count(), 3U,
            "corrected face-velocity workspace")),
        mass_flux_candidate(supplied_topology.local_face_count()),
        material_face_density(supplied_topology.local_face_count()),
        material_rhs_raw(supplied_topology.owned_cell_count()),
        material_rhs_solve(supplied_topology.owned_cell_count()),
        material_correction(supplied_topology.owned_cell_count()),
        material_diagonal(multiplied_count(
            supplied_topology.local_cell_count(), 3U,
            "material pressure diagonal workspace")),
        material_periodic_gamma(supplied_topology.global_face_count()),
        material_periodic_count(supplied_topology.global_face_count()) {}

  const runtime::StructuredDecomposition *decomposition;
  const mesh::MeshTopology *topology;
  const mesh::MeshGeometry *geometry;
  const boundary::BoundaryRegistry *boundaries;
  const runtime::MpiContext *mpi;
  execution::ExecutionContext *execution;
  runtime::HaloExchange *halo;
  const linear::LinearSolver *solver;
  linear::Preconditioner *preconditioner;
  finite_volume::CellCenteredFvmOperators fvm;
  finite_volume::FaceMassFlux::PreparedStatePtr prepared_face_flux;
  ScratchFields scratch;
  execution::Buffer gamma;
  execution::Buffer rhs;
  execution::Buffer correction;
  execution::Buffer residual;
  execution::Buffer candidate_pressure;
  std::vector<Real3> pressure_gradient_sums;
  std::vector<double> velocity_candidate;
  std::vector<double> pressure_candidate;
  std::vector<double> face_velocity_candidate;
  std::vector<double> mass_flux_candidate;
  std::vector<double> material_face_density;
  std::vector<double> material_rhs_raw;
  std::vector<double> material_rhs_solve;
  std::vector<double> material_correction;
  std::vector<double> material_diagonal;
  std::vector<double> material_periodic_gamma;
  std::vector<double> material_periodic_count;
  const FlowState *material_state{};
  std::uint64_t material_attempt_identity{};
  std::uint32_t material_ordinal{};
  bool material_token_available{};
  double material_rhs_l2{};
  linear::SolveControl material_control{};
  std::optional<finite_volume::MatrixFreePoissonOperator> pressure_operator;
  std::optional<finite_volume::MatrixFreePoissonOperator>
      material_final_operator;
  std::optional<finite_volume::PoissonConstraint> constraint;
};

PisoCoupler
PisoCoupler::create(const runtime::StructuredDecomposition &decomposition,
                    const mesh::MeshTopology &topology,
                    const mesh::MeshGeometry &geometry,
                    const boundary::BoundaryRegistry &boundaries,
                    const runtime::MpiContext &mpi,
                    execution::ExecutionContext &execution_context,
                    runtime::HaloExchange &cell_halo,
                    const linear::LinearSolver &pressure_solver,
                    linear::Preconditioner &pressure_preconditioner) {
  validate_host_context(execution_context);
  geometry.require_compatible(topology);
  if (!cell_halo.is_compatible_with(decomposition) ||
      cell_halo.ghost_width() != 2) {
    throw runtime::Error("PISO requires a compatible width-two Halo");
  }
  return PisoCoupler(std::make_unique<Impl>(
      decomposition, topology, geometry, boundaries, mpi, execution_context,
      cell_halo, pressure_solver, pressure_preconditioner));
}

PisoCoupler::PisoCoupler(std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
PisoCoupler::~PisoCoupler() noexcept = default;
PisoCoupler::PisoCoupler(PisoCoupler &&) noexcept = default;

runtime::HaloPerformanceCounters
PisoCoupler::pressure_halo_performance_counters() const {
  if (!impl_)
    throw runtime::Error("PISO coupler has been moved from");
  runtime::HaloPerformanceCounters result{};
  if (impl_->pressure_operator)
    result = combine_halo_counters(
        result, impl_->pressure_operator->halo_performance_counters());
  if (impl_->material_final_operator)
    result = combine_halo_counters(
        result, impl_->material_final_operator->halo_performance_counters());
  return result;
}

void PisoCoupler::prepare_material_density_assessment() {
  if (!impl_)
    throw runtime::Error("material PISO coupler has been moved from");
  if (impl_->material_final_operator)
    return;
  auto gamma = impl_->gamma.view(0U, impl_->topology->local_face_count());
  for (std::size_t face = 0; face < impl_->topology->local_face_count(); ++face)
    gamma[face] = 1.0;
  impl_->material_final_operator.emplace(
      finite_volume::MatrixFreePoissonOperator::create(
          *impl_->decomposition, *impl_->topology, *impl_->geometry,
          *impl_->execution,
          static_cast<const execution::Buffer &>(impl_->gamma).view(
              0U, impl_->topology->local_face_count()),
          finite_volume::make_poisson_boundary_spec(*impl_->boundaries)));
}

PressureCorrectionReport PisoCoupler::correct_common_throwing(
    FlowState &state, double rho_ref,
    const MomentumTimeStencil *material_stencil,
    const runtime::FieldView<const double> &actual_momentum_diagonal,
    const linear::SolveControl &control,
    PressureCorrectionReport &result) const {
  const std::size_t count = impl_->topology->owned_cell_count();
  const Int3 local_extent = impl_->decomposition->local_extent();
  const auto owned = impl_->topology->owned_global_box();
  const Int3 global_extent = impl_->topology->global_extent();
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::invalid_input, false,
      "pressure-correction input validation failed", [&] {
        if (!state.attempt_active()) {
          throw runtime::Error("pressure correction requires an active trial");
        }
        if (material_stencil == nullptr &&
            (!(rho_ref > 0.0) || !std::isfinite(rho_ref))) {
          throw runtime::Error("pressure correction density is invalid");
        }
        if (material_stencil != nullptr) {
          const auto expected = make_momentum_time_stencil(
              material_stencil->order, material_stencil->dt_s,
              material_stencil->previous_dt_s);
          if (expected.alpha0 != material_stencil->alpha0 ||
              expected.alpha1 != material_stencil->alpha1 ||
              expected.alpha2 != material_stencil->alpha2) {
            throw runtime::Error(
                "material pressure correction stencil is inconsistent");
          }
          if (impl_->material_state != &state ||
              impl_->material_attempt_identity != state.attempt_identity()) {
            impl_->material_state = &state;
            impl_->material_attempt_identity = state.attempt_identity();
            impl_->material_ordinal = 0U;
            impl_->material_token_available = false;
          }
          if (impl_->material_ordinal >= 2U) {
            throw runtime::Error(
                "material pressure correction count exceeds two");
          }
        }
        if (actual_momentum_diagonal.interior_extent().x != local_extent.x ||
            actual_momentum_diagonal.interior_extent().y != local_extent.y ||
            actual_momentum_diagonal.interior_extent().z != local_extent.z ||
            actual_momentum_diagonal.components() != 3U ||
            actual_momentum_diagonal.ghost_width() < 2) {
          throw runtime::Error(
              "pressure correction diagonal layout is invalid");
        }
        for (LocalCellId cell = 0; cell < impl_->topology->local_cell_count();
             ++cell) {
          const StructuredIndex index = map_cell(
              impl_->topology->global_cell(cell), owned, global_extent);
          const double scalar_diagonal =
              at(actual_momentum_diagonal, index, 0);
          for (int component_index = 0; component_index < 3;
               ++component_index) {
            const double component_diagonal =
                at(actual_momentum_diagonal, index, component_index);
            if (!(component_diagonal > 0.0) ||
                !std::isfinite(component_diagonal) ||
                component_diagonal != scalar_diagonal) {
              throw runtime::Error(
                  "pressure correction requires a positive finite isotropic "
                  "momentum diagonal");
            }
          }
        }
      });

  auto &trial = state.solver_layer(FlowLayer::trial);
  const auto &access = state.solver_access_plan();
  const auto &registry = state.solver_registry();
  const auto fields = state.fields();
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::invalid_input, false,
      "pressure-system derivation failed", [&] {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
        FaceFluxPathObservation allocation_observation(
            material_stencil != nullptr);
#endif
        auto mass = material_stencil != nullptr
                        ? finite_volume::FaceMassFlux::bind_prepared(
                              *impl_->prepared_face_flux, registry, trial,
                              access, kStatePhase, kStateActor,
                              fields.face_mass_flux, *impl_->topology)
                        : finite_volume::FaceMassFlux::acquire(
                              registry, trial, access, kStatePhase,
                              kStateActor, fields.face_mass_flux,
                              *impl_->topology);
        auto mass_residual = impl_->scratch.storage->acquire_write<double>(
            *impl_->scratch.access, kScratchPhase, kScratchActor,
            impl_->scratch.mass_residual);
        zero_cell(mass_residual);
        impl_->fvm.accumulate_mass_residual(mass, mass_residual);

        if (material_stencil != nullptr) {
          const auto rho_trial = trial.acquire_read<double>(
              access, kStatePhase, kStateActor, fields.density);
          const auto rho_n = state.solver_layer(FlowLayer::committed)
                                 .acquire_read<double>(
                                     access, kStatePhase, kStateActor,
                                     fields.density);
          const auto rho_nm1 = state.solver_layer(FlowLayer::history)
                                   .acquire_read<double>(
                                       access, kStatePhase, kStateActor,
                                       fields.density);
          auto face_density =
              impl_->scratch.storage->acquire_face_write<double>(
                  *impl_->scratch.access, kScratchPhase, kScratchActor,
                  impl_->scratch.scalar_face);
          impl_->fvm.reconstruct_transport_faces(
              finite_volume::FiniteVolumeQuantity::density(),
              *impl_->boundaries, mass, rho_trial, face_density);
          for (LocalFaceId face = 0;
               face < impl_->topology->local_face_count(); ++face) {
            const double value = face_density(face, 0);
            if (!(value > 0.0) || !std::isfinite(value)) {
              throw runtime::Error(
                  "material pressure face density is invalid");
            }
            impl_->material_face_density[face] = value;
          }
          for (LocalCellId cell = 0; cell < count; ++cell) {
            const StructuredIndex index = map_cell(
                impl_->topology->global_cell(cell), owned, global_extent);
            const double temporal =
                (material_stencil->alpha0 *
                     at(rho_trial, index, 0) +
                 material_stencil->alpha1 * at(rho_n, index, 0) +
                 material_stencil->alpha2 * at(rho_nm1, index, 0)) *
                impl_->geometry->cell_volume_m3(cell) /
                material_stencil->dt_s;
            if (!std::isfinite(temporal)) {
              throw runtime::Error(
                  "material pressure temporal residual is non-finite");
            }
            mass_residual(index.i, index.j, index.k, 0) += temporal;
          }
        } else {
          std::fill(impl_->material_face_density.begin(),
                    impl_->material_face_density.end(), rho_ref);
        }

        auto gamma = impl_->gamma.view(0U, impl_->topology->local_face_count());
        for (LocalFaceId face = 0; face < impl_->topology->local_face_count();
             ++face) {
          const LocalCellId owner = impl_->topology->owner(face);
          const StructuredIndex owner_index = map_cell(
              impl_->topology->global_cell(owner), owned, global_extent);
          const double volume_owner = impl_->geometry->cell_volume_m3(owner);
          const double owner_rate =
              at(actual_momentum_diagonal, owner_index, 0) / volume_owner;
          double lambda = owner_rate;
          const auto neighbour = impl_->topology->neighbour(face);
          if (neighbour.has_value()) {
            const StructuredIndex neighbour_index = map_cell(
                impl_->topology->global_cell(*neighbour), owned, global_extent);
            const double volume_neighbour =
                impl_->geometry->cell_volume_m3(*neighbour);
            const double neighbour_rate =
                at(actual_momentum_diagonal, neighbour_index, 0) /
                volume_neighbour;
            for (int component_index = 0; component_index < 3;
                 ++component_index) {
              const double owner_component =
                  at(actual_momentum_diagonal, owner_index, component_index) /
                  volume_owner;
              const double neighbour_component =
                  at(actual_momentum_diagonal, neighbour_index,
                     component_index) /
                  volume_neighbour;
              if (owner_component != owner_rate ||
                  neighbour_component != neighbour_rate) {
                throw runtime::Error(
                    "pressure correction requires an isotropic momentum "
                    "diagonal");
              }
            }
            if (material_stencil != nullptr &&
                impl_->topology->periodic_pair(face).has_value()) {
              lambda = 0.5 * (owner_rate + neighbour_rate);
            } else {
              const Real3 face_center = impl_->geometry->face_center_m(face);
              const Real3 owner_center = impl_->geometry->cell_center_m(owner);
              const Real3 neighbour_center =
                  impl_->geometry->cell_center_m(*neighbour);
              const double owner_distance =
                  std::sqrt(dot(subtract(face_center, owner_center),
                                subtract(face_center, owner_center)));
              const double neighbour_distance =
                  std::sqrt(dot(subtract(neighbour_center, face_center),
                                subtract(neighbour_center, face_center)));
              const double total = owner_distance + neighbour_distance;
              if (!(total > 0.0) || !std::isfinite(total)) {
                throw runtime::Error(
                    "pressure correction face weights are invalid");
              }
              lambda = (neighbour_distance / total) * owner_rate +
                       (owner_distance / total) * neighbour_rate;
            }
          } else {
            for (int component_index = 0; component_index < 3;
                 ++component_index) {
              if (at(actual_momentum_diagonal, owner_index, component_index) /
                      volume_owner !=
                  owner_rate) {
                throw runtime::Error(
                    "pressure correction requires an isotropic momentum "
                    "diagonal");
              }
            }
          }
          if (!(lambda > 0.0) || !std::isfinite(lambda)) {
            throw runtime::Error(
                "pressure correction face mobility is invalid");
          }
          gamma[face] = impl_->material_face_density[face] / lambda;
        }
        if (material_stencil != nullptr) {
          for (LocalFaceId face = 0;
               face < impl_->topology->local_face_count(); ++face) {
            const auto pair_global = impl_->topology->periodic_pair(face);
            if (!pair_global ||
                impl_->topology->global_face_id(face) < *pair_global)
              continue;
            const auto pair = impl_->topology->find_local_face(*pair_global);
            if (pair)
              gamma[face] = gamma[*pair];
          }
        }
      });

  if (material_stencil != nullptr) {
    std::fill(impl_->material_periodic_gamma.begin(),
              impl_->material_periodic_gamma.end(), 0.0);
    std::fill(impl_->material_periodic_count.begin(),
              impl_->material_periodic_count.end(), 0.0);
    auto gamma = impl_->gamma.view(0U, impl_->topology->local_face_count());
    for (LocalFaceId face = 0; face < impl_->topology->local_face_count();
         ++face) {
      const auto pair = impl_->topology->periodic_pair(face);
      const auto global = impl_->topology->global_face_id(face);
      if (!pair || global > *pair ||
          impl_->topology->cell_ownership(impl_->topology->owner(face)) !=
              EntityOwnership::owned)
        continue;
      const std::size_t key = static_cast<std::size_t>(global);
      impl_->material_periodic_gamma[key] += gamma[face];
      impl_->material_periodic_count[key] += 1.0;
    }
    impl_->mpi->allreduce_fp64_in_place(
        impl_->material_periodic_gamma.data(),
        impl_->material_periodic_gamma.size(),
        runtime::Fp64ReductionOperation::sum);
    impl_->mpi->allreduce_fp64_in_place(
        impl_->material_periodic_count.data(),
        impl_->material_periodic_count.size(),
        runtime::Fp64ReductionOperation::sum);
    for (LocalFaceId face = 0; face < impl_->topology->local_face_count();
         ++face) {
      const auto pair = impl_->topology->periodic_pair(face);
      if (!pair)
        continue;
      const std::size_t key = static_cast<std::size_t>(
          std::min(impl_->topology->global_face_id(face), *pair));
      const double count_value = impl_->material_periodic_count[key];
      if (!(count_value > 0.0) || !std::isfinite(count_value))
        throw runtime::Error(
            "material periodic pressure coefficient is unavailable");
      gamma[face] = impl_->material_periodic_gamma[key] / count_value;
    }
  }

  finite_volume::PoissonBoundarySpec boundary_spec;
  synchronized_local_phase(*impl_->mpi, StepFailureReason::invalid_input, false,
                           "pressure-boundary derivation failed", [&] {
                             boundary_spec =
                                 finite_volume::make_poisson_boundary_spec(
                                     *impl_->boundaries);
                           });
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  last_pressure_constraint_mode.store(static_cast<int>(boundary_spec.mode),
                                      std::memory_order_relaxed);
#endif
  if (!impl_->constraint.has_value()) {
    std::optional<finite_volume::PoissonConstraint> candidate_constraint;
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::invalid_input, false,
        "pressure-constraint construction failed", [&] {
          candidate_constraint.emplace(finite_volume::PoissonConstraint::create(
              *impl_->topology, *impl_->geometry, *impl_->execution,
              *impl_->mpi, boundary_spec.mode));
        });
    impl_->constraint.emplace(std::move(*candidate_constraint));
  }
  if (!impl_->pressure_operator.has_value()) {
    try {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      const int injected_rank =
          pressure_operator_construction_failure_rank.load(
              std::memory_order_relaxed);
      if (impl_->mpi->rank() == injected_rank &&
          impl_->topology->local_face_count() != 0U) {
        impl_->gamma.view(0U, impl_->topology->local_face_count())[0] =
            std::numeric_limits<double>::quiet_NaN();
      }
#endif
      impl_->pressure_operator.emplace(
          finite_volume::MatrixFreePoissonOperator::create(
              *impl_->decomposition, *impl_->topology, *impl_->geometry,
              *impl_->execution,
              static_cast<const execution::Buffer &>(impl_->gamma)
                  .view(0U, impl_->topology->local_face_count()),
              boundary_spec));
    } catch (const finite_volume::PoissonConstructionError &error) {
      throw SynchronizedAttemptFailure{StepFailureReason::invalid_input,
                                       error.failing_rank(), false};
    }
  } else {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    const int injected_rank = pressure_operator_refresh_failure_rank.load(
        std::memory_order_relaxed);
    if (impl_->mpi->rank() == injected_rank &&
        impl_->topology->local_face_count() != 0U) {
      impl_->gamma.view(0U, impl_->topology->local_face_count())[0] =
          std::numeric_limits<double>::quiet_NaN();
    }
#endif
    const auto replacement =
        impl_->pressure_operator->collectively_replace_face_coefficients(
            static_cast<const execution::Buffer &>(impl_->gamma)
                .view(0U, impl_->topology->local_face_count()),
            *impl_->mpi);
    if (!replacement.accepted) {
      throw SynchronizedAttemptFailure{StepFailureReason::invalid_input,
                                       replacement.lowest_failing_rank,
                                       false};
    }
  }
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  last_pressure_constraint_mode.store(
      static_cast<int>(impl_->constraint->mode()), std::memory_order_relaxed);
  last_pressure_operator_mode.store(
      static_cast<int>(impl_->pressure_operator->constraint_mode()),
      std::memory_order_relaxed);
#endif

  std::optional<execution::VectorView<double>> rhs;
  std::optional<execution::VectorView<double>> correction;
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::pressure_linear_solve, true,
      "pressure right-hand-side assembly failed", [&] {
        rhs.emplace(impl_->rhs.view(0U, count));
        correction.emplace(impl_->correction.view(0U, count));
        const auto mass_residual = impl_->scratch.storage->acquire_read<double>(
            *impl_->scratch.access, kScratchPhase, kScratchActor,
            impl_->scratch.mass_residual);
        for (LocalCellId cell = 0; cell < count; ++cell) {
          const Int3 global = impl_->topology->global_cell(cell);
          const StructuredIndex index = map_cell(global, owned, global_extent);
          (*rhs)[cell] = -mass_residual(index.i, index.j, index.k, 0) /
                         impl_->geometry->cell_volume_m3(cell);
          if (material_stencil != nullptr)
            impl_->material_rhs_raw[cell] = (*rhs)[cell];
          (*correction)[cell] = 0.0;
          if (!std::isfinite((*rhs)[cell])) {
            throw runtime::Error("pressure RHS is non-finite");
          }
        }
      });
  impl_->constraint->project_rhs(*rhs);
  if (material_stencil != nullptr) {
    for (std::size_t cell = 0; cell < count; ++cell)
      impl_->material_rhs_solve[cell] = (*rhs)[cell];
  }
  impl_->constraint->normalize_solution(*correction);
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::invalid_input, false,
      "pressure preconditioner update failed", [&] {
        impl_->preconditioner->update(*impl_->pressure_operator,
                                      impl_->pressure_operator->revision());
      });
  result.solve =
      impl_->solver->solve(*impl_->pressure_operator, *impl_->preconditioner,
                           *rhs, *correction, control);
  if (material_stencil != nullptr &&
      result.solve.reason == linear::SolveTerminationReason::collective_failure) {
    const int failing_rank = result.solve.lowest_failing_rank;
    if (failing_rank < 0) {
      throw runtime::MpiOperationError(
          "material-density flow collective failure has no reliable rank");
    }
    return correction_failure(
        std::move(result),
        PressureCorrectionDisposition::non_retryable_failure,
        StepFailureReason::collective_operation, failing_rank);
  }
  impl_->constraint->normalize_solution(*correction);
  std::optional<execution::VectorView<double>> residual;
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::pressure_linear_solve, true,
      "pressure residual view acquisition failed",
      [&] { residual.emplace(impl_->residual.view(0U, count)); });
  impl_->pressure_operator->apply(*correction, *residual).wait();
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::pressure_linear_solve, true,
      "pressure residual assembly failed", [&] {
        for (std::size_t cell = 0; cell < count; ++cell) {
          (*residual)[cell] = (*rhs)[cell] - (*residual)[cell];
          if (!std::isfinite((*residual)[cell])) {
            throw runtime::Error("pressure residual is non-finite");
          }
        }
      });
  result.independent_residual_l2 = global_l2(*impl_->mpi, *residual);
  result.rhs_l2 = global_l2(*impl_->mpi, *rhs);
  const double threshold = std::max(control.atol, control.rtol * result.rhs_l2);
  const double agreement = 64.0 * std::numeric_limits<double>::epsilon() *
                           std::max(1.0, result.independent_residual_l2);
  result.accepted = solve_success(result.solve.reason) &&
                    std::isfinite(result.independent_residual_l2) &&
                    result.independent_residual_l2 <= threshold &&
                    std::isfinite(result.solve.final_residual) &&
                    std::abs(result.solve.final_residual -
                             result.independent_residual_l2) <= agreement;
  const auto acceptance = runtime::collective_status(
      *impl_->mpi, result.accepted,
      "pressure correction did not satisfy the residual contract");
  if (!acceptance.ok) {
    result.accepted = false;
    const bool solver_failed = !solve_success(result.solve.reason);
    const int selected_rank =
        solver_failed && result.solve.lowest_failing_rank >= 0
            ? result.solve.lowest_failing_rank
            : acceptance.failing_rank;
    if (result.solve.lowest_failing_rank < 0) {
      result.solve.lowest_failing_rank = acceptance.failing_rank;
    }
    return correction_failure(
        std::move(result), PressureCorrectionDisposition::recoverable_failure,
        StepFailureReason::pressure_linear_solve, selected_rank);
  }

  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::non_finite_trial, true,
      "pressure-correction staging failed", [&] {
        auto pi_prime = impl_->scratch.storage->acquire_write<double>(
            *impl_->scratch.access, kScratchPhase, kScratchActor,
            impl_->scratch.pressure_correction);
        zero_cell(pi_prime);
        for (LocalCellId cell = 0; cell < count; ++cell) {
          const StructuredIndex index = map_cell(
              impl_->topology->global_cell(cell), owned, global_extent);
          pi_prime(index.i, index.j, index.k, 0) = (*correction)[cell];
        }
      });
  impl_->halo->exchange(*impl_->scratch.storage,
                        impl_->scratch.pressure_correction);
  std::optional<runtime::FieldView<const double>> pi_prime_read;
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::non_finite_trial, true,
      "pressure-correction gradient failed", [&] {
        pi_prime_read.emplace(impl_->scratch.storage->acquire_read<double>(
            *impl_->scratch.access, kScratchPhase, kScratchActor,
            impl_->scratch.pressure_correction));
        auto prime_gradient = impl_->scratch.storage->acquire_write<double>(
            *impl_->scratch.access, kScratchPhase, kScratchActor,
            impl_->scratch.pressure_gradient);
        compute_pressure_gradient(
            *impl_->topology, *impl_->geometry, *pi_prime_read, prime_gradient,
            impl_->pressure_gradient_sums,
            [&](LocalFaceId face, double owner_value) {
              const auto patch = impl_->topology->patch_id(face);
              if (patch.has_value() &&
                  impl_->boundaries->patch(*patch).pressure_rule() ==
                      boundary::PressureRule::prescribed_value) {
                return 0.0;
              }
              return owner_value;
            });
      });

  auto &velocity_candidate = impl_->velocity_candidate;
  auto &pressure_candidate = impl_->pressure_candidate;
  auto &face_velocity_candidate = impl_->face_velocity_candidate;
  auto &mass_flux_candidate = impl_->mass_flux_candidate;
  std::optional<execution::VectorView<double>> candidate_pressure;
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::non_finite_trial, true,
      "corrected cell candidate derivation failed", [&] {
        candidate_pressure.emplace(impl_->candidate_pressure.view(0U, count));
        const auto trial_velocity = trial.acquire_read<double>(
            access, kStatePhase, kStateActor, fields.velocity);
        const auto trial_pressure = trial.acquire_read<double>(
            access, kStatePhase, kStateActor, fields.mechanical_pressure);
        const auto prime_gradient =
            impl_->scratch.storage->acquire_read<double>(
                *impl_->scratch.access, kScratchPhase, kScratchActor,
                impl_->scratch.pressure_gradient);
        for (LocalCellId cell = 0; cell < count; ++cell) {
          const StructuredIndex index = map_cell(
              impl_->topology->global_cell(cell), owned, global_extent);
          const double volume = impl_->geometry->cell_volume_m3(cell);
          for (int component_index = 0; component_index < 3;
               ++component_index) {
            const double diagonal =
                at(actual_momentum_diagonal, index, component_index);
            velocity_candidate[cell * 3U +
                               static_cast<std::size_t>(component_index)] =
                trial_velocity(index.i, index.j, index.k, component_index) -
                (volume / diagonal) *
                    prime_gradient(index.i, index.j, index.k, component_index);
          }
          (*candidate_pressure)[cell] =
              trial_pressure(index.i, index.j, index.k, 0) +
              (*correction)[cell];
        }
      });
  impl_->constraint->normalize_solution(*candidate_pressure);
  double roundoff_scales[3]{};
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::non_finite_trial, true,
      "corrected face candidate derivation failed", [&] {
        const auto trial_face_velocity = trial.acquire_face_read<double>(
            access, kStatePhase, kStateActor, fields.face_velocity);
        const auto trial_mass_flux = trial.acquire_face_read<double>(
            access, kStatePhase, kStateActor, fields.face_mass_flux);
        const auto gamma_read =
            static_cast<const execution::Buffer &>(impl_->gamma)
                .view(0U, impl_->topology->local_face_count());
        for (LocalCellId cell = 0; cell < count; ++cell) {
          pressure_candidate[cell] = (*candidate_pressure)[cell];
        }
        for (LocalFaceId face = 0; face < impl_->topology->local_face_count();
             ++face) {
          const LocalCellId owner = impl_->topology->owner(face);
          const StructuredIndex owner_index = map_cell(
              impl_->topology->global_cell(owner), owned, global_extent);
          const double owner_value = at(*pi_prime_read, owner_index, 0);
          double neighbour_value = owner_value;
          const auto neighbour = impl_->topology->neighbour(face);
          bool correct_face = neighbour.has_value();
          if (neighbour.has_value()) {
            const StructuredIndex neighbour_index = map_cell(
                impl_->topology->global_cell(*neighbour), owned, global_extent);
            neighbour_value = at(*pi_prime_read, neighbour_index, 0);
          } else {
            const auto patch = impl_->topology->patch_id(face);
            correct_face = patch.has_value() &&
                           impl_->boundaries->patch(*patch).pressure_rule() ==
                               boundary::PressureRule::prescribed_value;
            if (correct_face)
              neighbour_value = 0.0;
          }
          const double lambda_inverse =
              gamma_read[face] / impl_->material_face_density[face];
          double flux_delta = 0.0;
          Real3 velocity_delta{};
          if (correct_face) {
            flux_delta = gamma_read[face] *
                         face_factor(*impl_->geometry, face) *
                         (owner_value - neighbour_value);
            const Real3 displacement =
                impl_->geometry->face_displacement_m(face);
            const double distance_squared = dot(displacement, displacement);
            const Real3 face_gradient =
                multiply((neighbour_value - owner_value) / distance_squared,
                         displacement);
            velocity_delta = multiply(-lambda_inverse, face_gradient);
          }
          mass_flux_candidate[face] = trial_mass_flux(face, 0) + flux_delta;
          for (int component_index = 0; component_index < 3;
               ++component_index) {
            face_velocity_candidate[face * 3U +
                                    static_cast<std::size_t>(component_index)] =
                trial_face_velocity(face, component_index) +
                component(velocity_delta, component_index);
          }
        }
        for (LocalFaceId face = 0; face < impl_->topology->local_face_count();
             ++face) {
          const auto pair_global = impl_->topology->periodic_pair(face);
          if (!pair_global.has_value() ||
              impl_->topology->global_face_id(face) < *pair_global) {
            continue;
          }
          const auto pair = impl_->topology->find_local_face(*pair_global);
          if (!pair.has_value())
            continue;
          mass_flux_candidate[face] = -mass_flux_candidate[*pair];
          for (int component_index = 0; component_index < 3;
               ++component_index) {
            face_velocity_candidate[face * 3U +
                                    static_cast<std::size_t>(component_index)] =
                face_velocity_candidate[*pair * 3U + static_cast<std::size_t>(
                                                         component_index)];
          }
        }
        for (LocalFaceId face = 0; face < impl_->topology->local_face_count();
             ++face) {
          roundoff_scales[0] =
              std::max(roundoff_scales[0], std::abs(mass_flux_candidate[face]));
          roundoff_scales[1] =
              std::max(roundoff_scales[1], impl_->geometry->face_area_m2(face));
          for (int component_index = 0; component_index < 3;
               ++component_index) {
            roundoff_scales[2] = std::max(
                roundoff_scales[2],
                std::abs(face_velocity_candidate[face * 3U +
                                                 static_cast<std::size_t>(
                                                     component_index)]));
          }
        }
      });
  impl_->mpi->allreduce_fp64_in_place(roundoff_scales, 3U,
                                      runtime::Fp64ReductionOperation::maximum);
  double density_scale = rho_ref;
  if (material_stencil != nullptr) {
    density_scale = 0.0;
    for (const double value : impl_->material_face_density)
      density_scale = std::max(density_scale, value);
  }
  const double flux_roundoff = 256.0 * std::numeric_limits<double>::epsilon() *
                               density_scale *
                               std::max(1.0, roundoff_scales[1]);
  if (roundoff_scales[0] <= flux_roundoff &&
      roundoff_scales[2] <= 256.0 * std::numeric_limits<double>::epsilon()) {
    std::fill(mass_flux_candidate.begin(), mass_flux_candidate.end(), 0.0);
    std::fill(face_velocity_candidate.begin(), face_velocity_candidate.end(),
              0.0);
  }
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::non_finite_trial, true,
      "corrected trial write failed", [&] {
        const auto all_finite = [](const std::vector<double> &values) {
          return std::all_of(values.begin(), values.end(),
                             [](double value) { return std::isfinite(value); });
        };
        if (!all_finite(velocity_candidate) ||
            !all_finite(pressure_candidate) ||
            !all_finite(face_velocity_candidate) ||
            !all_finite(mass_flux_candidate)) {
          throw runtime::Error("corrected trial is non-finite");
        }
        auto velocity_write = trial.acquire_write<double>(
            access, kStatePhase, kStateActor, fields.velocity);
        auto pressure_write = trial.acquire_write<double>(
            access, kStatePhase, kStateActor, fields.mechanical_pressure);
        auto face_velocity_write = trial.acquire_face_write<double>(
            access, kStatePhase, kStateActor, fields.face_velocity);
        auto mass_flux_write = trial.acquire_face_write<double>(
            access, kStatePhase, kStateActor, fields.face_mass_flux);
        for (LocalCellId cell = 0; cell < count; ++cell) {
          const StructuredIndex index = map_cell(
              impl_->topology->global_cell(cell), owned, global_extent);
          pressure_write(index.i, index.j, index.k, 0) =
              pressure_candidate[cell];
          for (int component_index = 0; component_index < 3;
               ++component_index) {
            velocity_write(index.i, index.j, index.k, component_index) =
                velocity_candidate[cell * 3U +
                                   static_cast<std::size_t>(component_index)];
          }
        }
        for (LocalFaceId face = 0; face < impl_->topology->local_face_count();
             ++face) {
          mass_flux_write(face, 0) = mass_flux_candidate[face];
          for (int component_index = 0; component_index < 3;
               ++component_index) {
            face_velocity_write(face, component_index) =
                face_velocity_candidate[face * 3U + static_cast<std::size_t>(
                                                        component_index)];
          }
        }
      });
  result.disposition = PressureCorrectionDisposition::accepted;
  result.reason = StepFailureReason::none;
  result.lowest_failing_rank = -1;
  result.accepted = true;
  if (material_stencil != nullptr) {
    ++impl_->material_ordinal;
    impl_->material_token_available = impl_->material_ordinal == 2U;
    impl_->material_rhs_l2 = result.rhs_l2;
    impl_->material_control = control;
    for (std::size_t cell = 0; cell < count; ++cell) {
      impl_->material_correction[cell] = (*correction)[cell];
    }
    for (LocalCellId cell = 0; cell < impl_->topology->local_cell_count();
         ++cell) {
      const StructuredIndex index = map_cell(
          impl_->topology->global_cell(cell), owned, global_extent);
      for (int direction = 0; direction < 3; ++direction)
        impl_->material_diagonal[cell * 3U +
                                 static_cast<std::size_t>(direction)] =
            at(actual_momentum_diagonal, index, direction);
    }
  }
  return result;
}

PressureCorrectionReport PisoCoupler::correct_throwing(
    FlowState &state, double rho_ref,
    const runtime::FieldView<const double> &actual_momentum_diagonal,
    const linear::SolveControl &control,
    PressureCorrectionReport &result) const {
  return correct_common_throwing(state, rho_ref, nullptr,
                                 actual_momentum_diagonal, control, result);
}

PressureCorrectionReport PisoCoupler::correct(
    FlowState &state, double rho_ref,
    const runtime::FieldView<const double> &actual_momentum_diagonal,
    const linear::SolveControl &control) const {
  PressureCorrectionReport result;
  try {
    return correct_throwing(state, rho_ref, actual_momentum_diagonal, control,
                            result);
  } catch (const SynchronizedAttemptFailure &failure) {
    return correction_failure(
        std::move(result),
        failure.recoverable
            ? PressureCorrectionDisposition::recoverable_failure
            : PressureCorrectionDisposition::non_retryable_failure,
        failure.reason, failure.failing_rank);
  } catch (const runtime::detail::MpiOperationError &) {
    return correction_failure(
        std::move(result), PressureCorrectionDisposition::non_retryable_failure,
        StepFailureReason::collective_operation, -1);
  } catch (const runtime::Error &) {
    return correction_failure(
        std::move(result), PressureCorrectionDisposition::non_retryable_failure,
        StepFailureReason::invalid_input, -1);
  } catch (...) {
    return correction_failure(
        std::move(result), PressureCorrectionDisposition::non_retryable_failure,
        StepFailureReason::invalid_input, -1);
  }
}

PressureCorrectionReport PisoCoupler::correct_material_density(
    FlowState &state, const MomentumTimeStencil &stencil,
    const runtime::FieldView<const double> &actual_momentum_diagonal,
    const linear::SolveControl &control) const {
  PressureCorrectionReport result;
  if (!impl_)
    return correction_failure(
        std::move(result), PressureCorrectionDisposition::non_retryable_failure,
        StepFailureReason::invalid_input, -1);
  try {
    result = correct_common_throwing(state, 1.0, &stencil,
                                     actual_momentum_diagonal, control, result);
  } catch (const SynchronizedAttemptFailure &failure) {
    result = correction_failure(
        std::move(result),
        failure.recoverable
            ? PressureCorrectionDisposition::recoverable_failure
            : PressureCorrectionDisposition::non_retryable_failure,
        failure.reason, failure.failing_rank);
  } catch (const runtime::detail::MpiOperationError &) {
    throw;
  } catch (const runtime::Error &) {
    result = correction_failure(
        std::move(result), PressureCorrectionDisposition::non_retryable_failure,
        StepFailureReason::invalid_input, -1);
  } catch (...) {
    result = correction_failure(
        std::move(result), PressureCorrectionDisposition::non_retryable_failure,
        StepFailureReason::invalid_input, -1);
  }
  if (!result.accepted)
    impl_->material_token_available = false;
  return result;
}

PisoCoupler::MaterialPressureAssessment
PisoCoupler::assess_final_material_density_pressure(
    FlowState &state, const MomentumTimeStencil &stencil,
    const linear::SolveControl &control) const {
  MaterialPressureAssessment result;
  if (!impl_)
    return result;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  const auto terminal_result = [&](std::uint8_t point) {
    const auto mode = detail::material_terminal_reach_raw(point);
    if (mode == kTerminalModeThrownOperation)
      throw runtime::MpiOperationError(
          "injected material terminal operation failure");
    if (mode != kTerminalModeReturnedRankless &&
        mode != kTerminalModeReturnedReliable)
      return false;
    result.disposition = PressureCorrectionDisposition::non_retryable_failure;
    result.reason = StepFailureReason::collective_operation;
    result.lowest_failing_rank =
        mode == kTerminalModeReturnedReliable ? 0 : -1;
    return true;
  };
  if (terminal_result(kTerminalFinalPressureEntry)) {
    impl_->material_token_available = false;
    return result;
  }
#endif
  try {
    const std::size_t count = impl_->topology->owned_cell_count();
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::invalid_input, false,
        "material final pressure token validation failed", [&] {
          const auto expected = make_momentum_time_stencil(
              stencil.order, stencil.dt_s, stencil.previous_dt_s);
          if (!state.attempt_active() || impl_->material_state != &state ||
              impl_->material_attempt_identity != state.attempt_identity() ||
              impl_->material_ordinal != 2U ||
              !impl_->material_token_available ||
              expected.alpha0 != stencil.alpha0 ||
              expected.alpha1 != stencil.alpha1 ||
              expected.alpha2 != stencil.alpha2 ||
              control.atol != impl_->material_control.atol ||
              control.rtol != impl_->material_control.rtol ||
              control.max_iterations !=
                  impl_->material_control.max_iterations ||
              control.residual_recompute_interval !=
                  impl_->material_control.residual_recompute_interval) {
            throw runtime::Error("material final pressure token is stale");
          }
        });

    auto &trial = state.solver_layer(FlowLayer::trial);
    const auto &access = state.solver_access_plan();
    const auto fields = state.fields();
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::invalid_input, false,
        "material final pressure coefficient refresh failed", [&] {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
          FaceFluxPathObservation allocation_observation;
#endif
          const auto mass = finite_volume::FaceMassFlux::bind_prepared(
              *impl_->prepared_face_flux, state.solver_registry(), trial,
              access, kStatePhase, kStateActor, fields.face_mass_flux,
              *impl_->topology);
          const auto density = trial.acquire_read<double>(
              access, kStatePhase, kStateActor, fields.density);
          auto face_density =
              impl_->scratch.storage->acquire_face_write<double>(
                  *impl_->scratch.access, kScratchPhase, kScratchActor,
                  impl_->scratch.scalar_face);
          impl_->fvm.reconstruct_transport_faces(
              finite_volume::FiniteVolumeQuantity::density(),
              *impl_->boundaries, mass, density, face_density);
          auto gamma = impl_->gamma.view(0U,
                                         impl_->topology->local_face_count());
          for (LocalFaceId face = 0;
               face < impl_->topology->local_face_count(); ++face) {
            const double density_face = face_density(face, 0);
            const LocalCellId owner = impl_->topology->owner(face);
            const double owner_volume =
                impl_->geometry->cell_volume_m3(owner);
            const double owner_rate =
                impl_->material_diagonal[owner * 3U] / owner_volume;
            double lambda = owner_rate;
            const auto neighbour = impl_->topology->neighbour(face);
            if (neighbour.has_value()) {
              const double neighbour_volume =
                  impl_->geometry->cell_volume_m3(*neighbour);
              const double neighbour_rate =
                  impl_->material_diagonal[*neighbour * 3U] /
                  neighbour_volume;
              if (impl_->topology->periodic_pair(face).has_value()) {
                lambda = 0.5 * (owner_rate + neighbour_rate);
              } else {
                const Real3 face_center =
                    impl_->geometry->face_center_m(face);
                const Real3 owner_center =
                    impl_->geometry->cell_center_m(owner);
                const Real3 neighbour_center =
                    impl_->geometry->cell_center_m(*neighbour);
                const double owner_distance =
                    std::sqrt(dot(subtract(face_center, owner_center),
                                  subtract(face_center, owner_center)));
                const double neighbour_distance = std::sqrt(dot(
                    subtract(neighbour_center, face_center),
                    subtract(neighbour_center, face_center)));
                lambda =
                    (neighbour_distance * owner_rate +
                     owner_distance * neighbour_rate) /
                    (owner_distance + neighbour_distance);
              }
            }
            if (!(density_face > 0.0) || !(lambda > 0.0) ||
                !std::isfinite(density_face) || !std::isfinite(lambda)) {
              throw runtime::Error(
                  "material final pressure coefficient is invalid");
            }
            gamma[face] = density_face / lambda;
          }
        });
    std::fill(impl_->material_periodic_gamma.begin(),
              impl_->material_periodic_gamma.end(), 0.0);
    std::fill(impl_->material_periodic_count.begin(),
              impl_->material_periodic_count.end(), 0.0);
    auto final_gamma =
        impl_->gamma.view(0U, impl_->topology->local_face_count());
    for (LocalFaceId face = 0; face < impl_->topology->local_face_count();
         ++face) {
      const auto pair = impl_->topology->periodic_pair(face);
      const auto global = impl_->topology->global_face_id(face);
      if (!pair || global > *pair ||
          impl_->topology->cell_ownership(impl_->topology->owner(face)) !=
              EntityOwnership::owned)
        continue;
      const std::size_t key = static_cast<std::size_t>(global);
      impl_->material_periodic_gamma[key] += final_gamma[face];
      impl_->material_periodic_count[key] += 1.0;
    }
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (terminal_result(
            kTerminalFinalPressureGammaSum)) {
      impl_->material_token_available = false;
      return result;
    }
#endif
    impl_->mpi->allreduce_fp64_in_place(
        impl_->material_periodic_gamma.data(),
        impl_->material_periodic_gamma.size(),
        runtime::Fp64ReductionOperation::sum);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (terminal_result(
            kTerminalFinalPressureGammaCount)) {
      impl_->material_token_available = false;
      return result;
    }
#endif
    impl_->mpi->allreduce_fp64_in_place(
        impl_->material_periodic_count.data(),
        impl_->material_periodic_count.size(),
        runtime::Fp64ReductionOperation::sum);
    for (LocalFaceId face = 0; face < impl_->topology->local_face_count();
         ++face) {
      const auto pair = impl_->topology->periodic_pair(face);
      if (!pair)
        continue;
      const std::size_t key = static_cast<std::size_t>(
          std::min(impl_->topology->global_face_id(face), *pair));
      const double count_value = impl_->material_periodic_count[key];
      if (!(count_value > 0.0) || !std::isfinite(count_value))
        throw runtime::Error(
            "material final periodic pressure coefficient is unavailable");
      final_gamma[face] =
          impl_->material_periodic_gamma[key] / count_value;
    }
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::invalid_input, false,
        "material final pressure operator is unavailable", [&] {
          if (!impl_->material_final_operator)
            throw runtime::Error(
                "material final pressure operator was not prepared");
        });
    const auto replacement =
        impl_->material_final_operator->collectively_replace_face_coefficients(
            static_cast<const execution::Buffer &>(impl_->gamma).view(
                0U, impl_->topology->local_face_count()),
            *impl_->mpi);
    if (!replacement.accepted)
      throw SynchronizedAttemptFailure{StepFailureReason::invalid_input,
                                       replacement.lowest_failing_rank, false};
    auto rhs = impl_->rhs.view(0U, count);
    auto correction = impl_->correction.view(0U, count);
    auto residual = impl_->residual.view(0U, count);
    for (std::size_t cell = 0; cell < count; ++cell) {
      rhs[cell] = impl_->material_rhs_solve[cell];
      correction[cell] = impl_->material_correction[cell];
    }
    impl_->material_final_operator->apply(correction, residual).wait();
    for (std::size_t cell = 0; cell < count; ++cell)
      residual[cell] = rhs[cell] - residual[cell];
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (terminal_result(kTerminalFinalPressureResidualReduction)) {
      impl_->material_token_available = false;
      return result;
    }
#endif
    result.independent_residual_l2 = global_l2(*impl_->mpi, residual);
    result.rhs_l2 = impl_->material_rhs_l2;
    const double threshold =
        std::max(control.atol, control.rtol * result.rhs_l2);
    result.normalized_residual =
        threshold > 0.0
            ? result.independent_residual_l2 / threshold
            : (result.independent_residual_l2 == 0.0
                   ? 0.0
                   : std::numeric_limits<double>::infinity());
    result.residual_available = true;
    result.accepted = std::isfinite(result.independent_residual_l2) &&
                      result.independent_residual_l2 <= threshold;
    result.disposition =
        result.accepted ? PressureCorrectionDisposition::accepted
                        : PressureCorrectionDisposition::recoverable_failure;
    result.reason = result.accepted ? StepFailureReason::none
                                    : StepFailureReason::final_pressure_residual;
    if (result.accepted) {
      result.lowest_failing_rank = -1;
    } else {
      const auto status = runtime::collective_status(
          *impl_->mpi, false,
          "material final pressure residual exceeds tolerance");
      result.lowest_failing_rank = status.failing_rank;
    }
  } catch (const SynchronizedAttemptFailure &failure) {
    result.disposition =
        failure.recoverable
            ? PressureCorrectionDisposition::recoverable_failure
            : PressureCorrectionDisposition::non_retryable_failure;
    result.reason = failure.reason;
    result.lowest_failing_rank = failure.failing_rank;
  } catch (const runtime::detail::MpiOperationError &) {
    throw;
  } catch (...) {
    result.disposition = PressureCorrectionDisposition::non_retryable_failure;
    result.reason = StepFailureReason::invalid_input;
    result.lowest_failing_rank = -1;
  }
  impl_->material_token_available = false;
  return result;
}

namespace {

struct MomentumSpatialResidual final {
  std::vector<double> convection;
  std::vector<double> viscosity;
  std::vector<finite_volume::PhysicalBoundaryMomentumContribution> boundary;
};

struct TransportSpatialResidual final {
  std::vector<double> convection;
  std::vector<double> diffusion;
  std::vector<finite_volume::PhysicalBoundaryTransportContribution> boundary;
};

} // namespace

struct FixedStepConstantDensityFlow::Impl final {
  Impl(const runtime::StructuredDecomposition &supplied_decomposition,
       const mesh::MeshTopology &supplied_topology,
       const mesh::MeshGeometry &supplied_geometry,
       const boundary::BoundaryRegistry &supplied_boundaries,
       const runtime::MpiContext &supplied_mpi,
       execution::ExecutionContext &supplied_execution,
       runtime::HaloExchange &supplied_halo,
       const linear::LinearSolver &supplied_momentum_solver,
       std::array<linear::Preconditioner *, 3> supplied_preconditioners,
       PisoCoupler supplied_coupler,
       std::vector<ConstantDensityTransportSpec> supplied_transport,
       bool supplied_transport_specs_valid)
      : decomposition(&supplied_decomposition), topology(&supplied_topology),
        geometry(&supplied_geometry), boundaries(&supplied_boundaries),
        mpi(&supplied_mpi), execution(&supplied_execution),
        halo(&supplied_halo), momentum_predictor(supplied_momentum_solver),
        momentum_preconditioners(supplied_preconditioners),
        coupler(std::move(supplied_coupler)),
        transport(std::move(supplied_transport)),
        transport_specs_valid(supplied_transport_specs_valid),
        fvm(finite_volume::CellCenteredFvmOperators::create(supplied_topology,
                                                            supplied_geometry)),
        face_assembler(TimeConsistentFaceVelocity::create(supplied_topology,
                                                          supplied_geometry)),
        scratch({supplied_decomposition.local_extent(),
                 supplied_topology.local_face_count()}),
        rhs{execution::Buffer(supplied_execution,
                              bytes_for(supplied_topology.owned_cell_count())),
            execution::Buffer(supplied_execution,
                              bytes_for(supplied_topology.owned_cell_count())),
            execution::Buffer(supplied_execution,
                              bytes_for(supplied_topology.owned_cell_count()))},
        predictor{
            execution::Buffer(supplied_execution,
                              bytes_for(supplied_topology.owned_cell_count())),
            execution::Buffer(supplied_execution,
                              bytes_for(supplied_topology.owned_cell_count())),
            execution::Buffer(supplied_execution,
                              bytes_for(supplied_topology.owned_cell_count()))},
        diagonal{
            execution::Buffer(supplied_execution,
                              bytes_for(supplied_topology.owned_cell_count())),
            execution::Buffer(supplied_execution,
                              bytes_for(supplied_topology.owned_cell_count())),
            execution::Buffer(
                supplied_execution,
                bytes_for(supplied_topology.owned_cell_count()))} {
    std::vector<double> coefficients;
    coefficients.reserve(transport.size());
    for (const auto &item : transport)
      coefficients.push_back(item.diffusivity_kg_per_m_s);
    transport_authority =
        detail::make_transport_diffusivity_authority(coefficients);
    const std::size_t cell_count = supplied_topology.owned_cell_count();
    const std::size_t momentum_count = multiplied_count(
        cell_count, 3U, "momentum workspace");
    pressure_gradient_sums.resize(cell_count);
    continuity_absolute.resize(cell_count);
    for (auto &values : diagonal_values)
      values.resize(cell_count);
    momentum_n.convection.resize(momentum_count);
    momentum_n.viscosity.resize(momentum_count);
    momentum_nm1.convection.resize(momentum_count);
    momentum_nm1.viscosity.resize(momentum_count);
    momentum_n.boundary.reserve(supplied_topology.local_face_count());
    momentum_nm1.boundary.reserve(supplied_topology.local_face_count());
    pressure_boundary.reserve(supplied_topology.local_face_count());
    transport_n.resize(transport.size());
    transport_nm1.resize(transport.size());
    for (auto *set : {&transport_n, &transport_nm1}) {
      for (auto &values : *set) {
        values.convection.resize(cell_count);
        values.diffusion.resize(cell_count);
        values.boundary.reserve(supplied_topology.local_face_count());
      }
    }
    const auto layout = linear::VectorLayout::from_topology(supplied_topology);
    for (std::size_t component_index = 0; component_index < 3U;
         ++component_index) {
      operators[component_index] = std::make_unique<DiagonalMomentumOperator>(
          supplied_execution, layout);
    }
  }

  const runtime::StructuredDecomposition *decomposition;
  const mesh::MeshTopology *topology;
  const mesh::MeshGeometry *geometry;
  const boundary::BoundaryRegistry *boundaries;
  const runtime::MpiContext *mpi;
  execution::ExecutionContext *execution;
  runtime::HaloExchange *halo;
  MomentumPredictor momentum_predictor;
  std::array<linear::Preconditioner *, 3> momentum_preconditioners;
  PisoCoupler coupler;
  std::vector<ConstantDensityTransportSpec> transport;
  detail::TransportDiffusivityAuthority transport_authority;
  bool transport_specs_valid{};
  finite_volume::CellCenteredFvmOperators fvm;
  TimeConsistentFaceVelocity face_assembler;
  ScratchFields scratch;
  std::array<execution::Buffer, 3> rhs;
  std::array<execution::Buffer, 3> predictor;
  std::array<execution::Buffer, 3> diagonal;
  std::array<std::unique_ptr<DiagonalMomentumOperator>, 3> operators;
  std::vector<Real3> pressure_gradient_sums;
  std::vector<double> continuity_absolute;
  std::array<std::vector<double>, 3> diagonal_values;
  MomentumSpatialResidual momentum_n;
  MomentumSpatialResidual momentum_nm1;
  std::vector<finite_volume::PhysicalBoundaryPressureContribution>
      pressure_boundary;
  std::vector<TransportSpatialResidual> transport_n;
  std::vector<TransportSpatialResidual> transport_nm1;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  std::vector<double> provisional_face_mass_flux;
#endif
};

FixedStepConstantDensityFlow FixedStepConstantDensityFlow::create(
    const runtime::StructuredDecomposition &decomposition,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const boundary::BoundaryRegistry &boundaries,
    const runtime::MpiContext &mpi,
    execution::ExecutionContext &execution_context,
    runtime::HaloExchange &cell_halo,
    const linear::LinearSolver &momentum_solver,
    std::array<linear::Preconditioner *, 3> momentum_preconditioners,
    const linear::LinearSolver &pressure_solver,
    linear::Preconditioner &pressure_preconditioner,
    std::vector<ConstantDensityTransportSpec> transported_fields) {
  validate_host_context(execution_context);
  geometry.require_compatible(topology);
  if (!cell_halo.is_compatible_with(decomposition) ||
      cell_halo.ghost_width() != 2) {
    throw runtime::Error("fixed-step flow requires width-two Halo");
  }
  if (std::any_of(momentum_preconditioners.begin(),
                  momentum_preconditioners.end(),
                  [](const auto *item) { return item == nullptr; })) {
    throw runtime::Error("momentum preconditioner is null");
  }
  bool transport_specs_valid = true;
  bool has_enthalpy = false;
  std::set<std::size_t> scalar_indices;
  std::set<runtime::FieldId> unique_transport;
  for (const auto &item : transported_fields) {
    const bool unique_field = unique_transport.insert(item.field).second;
    bool valid_quantity = false;
    switch (item.quantity.kind) {
    case finite_volume::FiniteVolumeQuantityKind::enthalpy:
      valid_quantity = item.quantity.scalar_index == 0U && !has_enthalpy;
      has_enthalpy = true;
      break;
    case finite_volume::FiniteVolumeQuantityKind::scalar:
      valid_quantity =
          item.quantity.scalar_index < boundaries.scalar_count() &&
          scalar_indices.insert(item.quantity.scalar_index).second;
      break;
    case finite_volume::FiniteVolumeQuantityKind::density:
    case finite_volume::FiniteVolumeQuantityKind::velocity:
    default:
      valid_quantity = false;
      break;
    }
    transport_specs_valid =
        transport_specs_valid && valid_quantity && unique_field &&
        item.diffusivity_kg_per_m_s >= 0.0 &&
        std::isfinite(item.diffusivity_kg_per_m_s);
  }
  auto coupler = PisoCoupler::create(
      decomposition, topology, geometry, boundaries, mpi, execution_context,
      cell_halo, pressure_solver, pressure_preconditioner);
  return FixedStepConstantDensityFlow(std::make_unique<Impl>(
      decomposition, topology, geometry, boundaries, mpi, execution_context,
      cell_halo, momentum_solver, momentum_preconditioners, std::move(coupler),
      std::move(transported_fields), transport_specs_valid));
}

FixedStepConstantDensityFlow::FixedStepConstantDensityFlow(
    std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}

runtime::HaloPerformanceCounters
FixedStepConstantDensityFlow::pressure_halo_performance_counters() const {
  if (!impl_)
    throw runtime::Error("fixed-step flow has been moved from");
  return impl_->coupler.pressure_halo_performance_counters();
}
FixedStepConstantDensityFlow::~FixedStepConstantDensityFlow() noexcept =
    default;
FixedStepConstantDensityFlow::FixedStepConstantDensityFlow(
    FixedStepConstantDensityFlow &&) noexcept = default;

const detail::TransportDiffusivityAuthority &
detail::AdaptiveTimeControlAccess::authority(
    const FixedStepConstantDensityFlow &flow) noexcept {
  static const TransportDiffusivityAuthority empty{};
  return flow.impl_ ? flow.impl_->transport_authority : empty;
}

detail::FacadeAssemblyIdentity detail::AdaptiveTimeControlAccess::assembly(
    const FixedStepConstantDensityFlow &flow) noexcept {
  if (!flow.impl_)
    return {};
  return {true, &flow, flow.impl_->topology, flow.impl_->geometry,
          flow.impl_->mpi};
}

namespace {

template <class FlowImplementation>
void assemble_spatial_residual(
    FlowImplementation &impl, const runtime::FieldRegistry &registry,
    const runtime::FieldAccessPlan &state_access,
    runtime::FieldStorage &accepted, const FlowFieldIds &fields, double mu,
    const runtime::FaceFieldView<const double> *dynamic_viscosity_by_face,
    MomentumSpatialResidual &values,
    runtime::FieldStorage *flux_storage = nullptr) {
  synchronized_local_phase(
      *impl.mpi, StepFailureReason::non_finite_trial, true,
      "accepted momentum gradient failed", [&] {
        const auto velocity = accepted.acquire_read<double>(
            state_access, kStatePhase, kStateActor, fields.velocity);
        auto gradient = impl.scratch.storage->template acquire_write<double>(
            *impl.scratch.access, kScratchPhase, kScratchActor,
            impl.scratch.velocity_gradient);
        impl.fvm.compute_gradient(
            finite_volume::GradientScheme::green_gauss,
            finite_volume::FiniteVolumeQuantity::velocity(), *impl.boundaries,
            velocity, gradient);
      });
  impl.halo->exchange(*impl.scratch.storage, impl.scratch.velocity_gradient);
  synchronized_local_phase(
      *impl.mpi, StepFailureReason::non_finite_trial, true,
      "accepted momentum residual failed", [&] {
        const auto velocity = accepted.acquire_read<double>(
            state_access, kStatePhase, kStateActor, fields.velocity);
        const auto gradient_read =
            impl.scratch.storage->template acquire_read<double>(
                *impl.scratch.access, kScratchPhase, kScratchActor,
                impl.scratch.velocity_gradient);
        const auto flux = finite_volume::FaceMassFlux::acquire(
            registry, flux_storage == nullptr ? accepted : *flux_storage,
            state_access, kStatePhase, kStateActor,
            fields.face_mass_flux, *impl.topology);
        auto face = impl.scratch.storage->template acquire_face_write<double>(
            *impl.scratch.access, kScratchPhase, kScratchActor,
            impl.scratch.momentum_face);
        impl.fvm.reconstruct_momentum_faces(*impl.boundaries, flux, velocity,
                                            face);
        const auto face_read =
            impl.scratch.storage->template acquire_face_read<double>(
                *impl.scratch.access, kScratchPhase, kScratchActor,
                impl.scratch.momentum_face);
        auto residual = impl.scratch.storage->template acquire_write<double>(
            *impl.scratch.access, kScratchPhase, kScratchActor,
            impl.scratch.momentum_residual);
        zero_cell(residual);
        impl.fvm.accumulate_convective_residual(flux, face_read, residual);
        const auto owned = impl.topology->owned_global_box();
        for (LocalCellId cell = 0; cell < impl.topology->owned_cell_count();
             ++cell) {
          const StructuredIndex index =
              map_cell(impl.topology->global_cell(cell), owned,
                       impl.topology->global_extent());
          for (int component_index = 0; component_index < 3;
               ++component_index) {
            values.convection[cell * 3U +
                              static_cast<std::size_t>(component_index)] =
                residual(index.i, index.j, index.k, component_index);
          }
        }
        zero_cell(residual);
        if (dynamic_viscosity_by_face == nullptr) {
          impl.fvm.accumulate_viscous_residual(*impl.boundaries, velocity,
                                               gradient_read, mu, residual);
        } else {
          impl.fvm.accumulate_viscous_residual(
              *impl.boundaries, velocity, gradient_read,
              *dynamic_viscosity_by_face, residual);
        }
        for (LocalCellId cell = 0; cell < impl.topology->owned_cell_count();
             ++cell) {
          const StructuredIndex index =
              map_cell(impl.topology->global_cell(cell), owned,
                       impl.topology->global_extent());
          for (int component_index = 0; component_index < 3;
               ++component_index) {
            values.viscosity[cell * 3U +
                             static_cast<std::size_t>(component_index)] =
                residual(index.i, index.j, index.k, component_index);
          }
        }
        if (dynamic_viscosity_by_face == nullptr) {
          impl.fvm.physical_boundary_momentum_contributions(
              *impl.boundaries, flux, face_read, velocity, gradient_read, mu,
              values.boundary);
        } else {
          impl.fvm.physical_boundary_momentum_contributions(
              *impl.boundaries, flux, face_read, velocity, gradient_read,
              *dynamic_viscosity_by_face, values.boundary);
        }
      });
}

template <class FlowImplementation, class StateImpl>
void exchange_accepted_layers(FlowImplementation &impl, StateImpl &state) {
  for (FlowLayer layer : {FlowLayer::history, FlowLayer::committed}) {
    auto &storage = detail::FlowStateSolverAccess::layer(state, layer);
    impl.halo->exchange(storage, state.fields().density);
    impl.halo->exchange(storage, state.fields().velocity);
    impl.halo->exchange(storage, state.fields().mechanical_pressure);
    for (const auto field : state.fields().transported_cell_fields) {
      impl.halo->exchange(storage, field);
    }
  }
}

template <class FlowImplementation>
void assemble_transport_spatial_residual(
    FlowImplementation &impl, FlowState &state,
    runtime::FieldStorage &accepted,
    const ConstantDensityTransportSpec &item,
    const runtime::FaceFieldView<const double> *sgs_viscosity_by_face,
    const les::WaleControl *wale_control,
    TransportSpatialResidual &spatial) {
  auto &trial = detail::FlowStateSolverAccess::layer(state, FlowLayer::trial);
  const auto &registry = detail::FlowStateSolverAccess::registry(state);
  const auto &access = detail::FlowStateSolverAccess::access(state);
  synchronized_local_phase(
      *impl.mpi, StepFailureReason::transport_failure, true,
      "transport gradient failed", [&] {
        const auto values = accepted.acquire_read<double>(
            access, kStatePhase, kStateActor, item.field);
        auto gradient = impl.scratch.storage->template acquire_write<double>(
            *impl.scratch.access, kScratchPhase, kScratchActor,
            impl.scratch.scalar_gradient);
        impl.fvm.compute_gradient(finite_volume::GradientScheme::green_gauss,
                                  item.quantity, *impl.boundaries, values,
                                  gradient);
      });
  impl.halo->exchange(*impl.scratch.storage, impl.scratch.scalar_gradient);
  synchronized_local_phase(
      *impl.mpi, StepFailureReason::transport_failure, true,
      "transport spatial residual failed", [&] {
        const auto values = accepted.acquire_read<double>(
            access, kStatePhase, kStateActor, item.field);
        const auto flux = finite_volume::FaceMassFlux::acquire(
            registry, trial, access, kStatePhase, kStateActor,
            state.fields().face_mass_flux, *impl.topology);
        const auto gradient =
            impl.scratch.storage->template acquire_read<double>(
                *impl.scratch.access, kScratchPhase, kScratchActor,
                impl.scratch.scalar_gradient);
        auto face = impl.scratch.storage->template acquire_face_write<double>(
            *impl.scratch.access, kScratchPhase, kScratchActor,
            impl.scratch.scalar_face);
        impl.fvm.reconstruct_transport_faces(item.quantity, *impl.boundaries,
                                             flux, values, face);
        auto gamma = impl.scratch.storage->template acquire_face_write<double>(
            *impl.scratch.access, kScratchPhase, kScratchActor,
            impl.scratch.scalar_gamma);
        const double turbulent_divisor =
            item.quantity.kind ==
                    finite_volume::FiniteVolumeQuantityKind::enthalpy
                ? (wale_control == nullptr
                       ? 1.0
                       : wale_control->turbulent_prandtl)
                : (wale_control == nullptr
                       ? 1.0
                       : wale_control->turbulent_schmidt);
        for (LocalFaceId local_face = 0;
             local_face < impl.topology->local_face_count(); ++local_face) {
          gamma(local_face, 0) =
              item.diffusivity_kg_per_m_s +
              (sgs_viscosity_by_face == nullptr
                   ? 0.0
                   : (*sgs_viscosity_by_face)(local_face, 0) /
                         turbulent_divisor);
          if (!(gamma(local_face, 0) >= 0.0) ||
              !std::isfinite(gamma(local_face, 0)))
            throw runtime::Error(
                "WALE transport diffusion coefficient is invalid");
        }
        const auto face_values =
            impl.scratch.storage->template acquire_face_read<double>(
                *impl.scratch.access, kScratchPhase, kScratchActor,
                impl.scratch.scalar_face);
        const auto gamma_values =
            impl.scratch.storage->template acquire_face_read<double>(
                *impl.scratch.access, kScratchPhase, kScratchActor,
                impl.scratch.scalar_gamma);
        auto residual = impl.scratch.storage->template acquire_write<double>(
            *impl.scratch.access, kScratchPhase, kScratchActor,
            impl.scratch.scalar_residual);
        zero_cell(residual);
        impl.fvm.accumulate_convective_residual(flux, face_values, residual);
        const auto owned = impl.topology->owned_global_box();
        for (LocalCellId cell = 0; cell < impl.topology->owned_cell_count();
             ++cell) {
          const StructuredIndex index =
              map_cell(impl.topology->global_cell(cell), owned,
                       impl.topology->global_extent());
          spatial.convection[cell] = residual(index.i, index.j, index.k, 0);
        }
        zero_cell(residual);
        impl.fvm.accumulate_scalar_diffusive_residual(
            item.quantity, *impl.boundaries, values, gradient, gamma_values,
            residual);
        for (LocalCellId cell = 0; cell < impl.topology->owned_cell_count();
             ++cell) {
          const StructuredIndex index =
              map_cell(impl.topology->global_cell(cell), owned,
                       impl.topology->global_extent());
          spatial.diffusion[cell] = residual(index.i, index.j, index.k, 0);
        }
        impl.fvm.physical_boundary_transport_contributions(
            item.quantity, *impl.boundaries, flux, face_values, values,
            gradient, gamma_values, spatial.boundary);
      });
}

template <class FlowImplementation>
void recompute_transport(FlowImplementation &impl, FlowState &state,
                         double rho_ref, const MomentumTimeStencil &stencil,
                         const runtime::FaceFieldView<const double> *
                             sgs_viscosity_by_face,
                         const les::WaleControl *wale_control,
                         bool final_call) {
  if (impl.transport.empty())
    return;
  auto &committed =
      detail::FlowStateSolverAccess::layer(state, FlowLayer::committed);
  auto &history =
      detail::FlowStateSolverAccess::layer(state, FlowLayer::history);
  auto &trial = detail::FlowStateSolverAccess::layer(state, FlowLayer::trial);
  const auto &access = detail::FlowStateSolverAccess::access(state);
  for (std::size_t transport_index = 0;
       transport_index < impl.transport.size(); ++transport_index) {
    const auto &item = impl.transport[transport_index];
    auto &spatial_n = impl.transport_n[transport_index];
    auto &spatial_nm1 = impl.transport_nm1[transport_index];
    assemble_transport_spatial_residual(
        impl, state, committed, item, sgs_viscosity_by_face, wale_control,
        spatial_n);
    if (stencil.order == MomentumTimeOrder::bdf2) {
      assemble_transport_spatial_residual(
          impl, state, history, item, sgs_viscosity_by_face, wale_control,
          spatial_nm1);
    }
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (final_call &&
        transport_assembly_mutation.load(std::memory_order_relaxed) ==
            kTransportMutationOmitHistorySpatial) {
      std::fill(spatial_nm1.convection.begin(), spatial_nm1.convection.end(),
                0.0);
      std::fill(spatial_nm1.diffusion.begin(), spatial_nm1.diffusion.end(),
                0.0);
    }
#else
    static_cast<void>(final_call);
#endif
    synchronized_local_phase(
        *impl.mpi, StepFailureReason::transport_failure, true,
        "transport finalization failed", [&] {
          const auto current = committed.acquire_read<double>(
              access, kStatePhase, kStateActor, item.field);
          const auto previous = history.acquire_read<double>(
              access, kStatePhase, kStateActor, item.field);
          auto output = trial.acquire_write<double>(access, kStatePhase,
                                                    kStateActor, item.field);
          const auto owned = impl.topology->owned_global_box();
          for (LocalCellId cell = 0; cell < impl.topology->owned_cell_count();
               ++cell) {
            const StructuredIndex index =
                map_cell(impl.topology->global_cell(cell), owned,
                         impl.topology->global_extent());
            const double spatial =
                stencil.order == MomentumTimeOrder::backward_euler
                    ? spatial_n.convection[cell] + spatial_n.diffusion[cell]
                    : 2.0 * (spatial_n.convection[cell] +
                             spatial_n.diffusion[cell]) -
                          (spatial_nm1.convection[cell] +
                           spatial_nm1.diffusion[cell]);
            const double value =
                (-(stencil.alpha1 * current(index.i, index.j, index.k, 0) +
                   stencil.alpha2 * previous(index.i, index.j, index.k, 0)) -
                 stencil.dt_s * spatial /
                     (rho_ref * impl.geometry->cell_volume_m3(cell))) /
                stencil.alpha0;
            if (!std::isfinite(value)) {
              throw runtime::Error(
                  "transport finalization is non-finite");
            }
            output(index.i, index.j, index.k, 0) = value;
          }
        });
  }
}

template <class FlowImplementation>
std::pair<double, double> continuity_norms(FlowImplementation &impl,
                                           FlowState &state) {
  auto &trial = detail::FlowStateSolverAccess::layer(state, FlowLayer::trial);
  double sums[2]{};
  synchronized_local_phase(
      *impl.mpi, StepFailureReason::final_continuity_residual, true,
      "continuity derivation failed", [&] {
        const auto flux = finite_volume::FaceMassFlux::acquire(
            detail::FlowStateSolverAccess::registry(state), trial,
            detail::FlowStateSolverAccess::access(state), kStatePhase,
            kStateActor, state.fields().face_mass_flux, *impl.topology);
        auto raw = impl.scratch.storage->template acquire_write<double>(
            *impl.scratch.access, kScratchPhase, kScratchActor,
            impl.scratch.mass_residual);
        zero_cell(raw);
        impl.fvm.accumulate_mass_residual(flux, raw);
        if (impl.continuity_absolute.size() !=
            impl.topology->owned_cell_count()) {
          throw runtime::Error("continuity workspace is invalid");
        }
        std::fill(impl.continuity_absolute.begin(),
                  impl.continuity_absolute.end(), 0.0);
        const auto flux_view = trial.acquire_face_read<double>(
            detail::FlowStateSolverAccess::access(state), kStatePhase,
            kStateActor, state.fields().face_mass_flux);
        for (LocalFaceId face = 0; face < impl.topology->local_face_count();
             ++face) {
          const double magnitude = std::abs(flux_view(face, 0));
          const LocalCellId owner = impl.topology->owner(face);
          if (impl.topology->cell_ownership(owner) == EntityOwnership::owned) {
            impl.continuity_absolute[owner] += magnitude;
          }
          const auto neighbour = impl.topology->neighbour(face);
          if (!impl.topology->patch_id(face).has_value() &&
              neighbour.has_value() &&
              impl.topology->cell_ownership(*neighbour) ==
                  EntityOwnership::owned) {
            impl.continuity_absolute[*neighbour] += magnitude;
          }
        }
        const auto owned = impl.topology->owned_global_box();
        for (LocalCellId cell = 0; cell < impl.topology->owned_cell_count();
             ++cell) {
          const StructuredIndex index =
              map_cell(impl.topology->global_cell(cell), owned,
                       impl.topology->global_extent());
          const double value = raw(index.i, index.j, index.k, 0);
          sums[0] += value * value;
          sums[1] += impl.continuity_absolute[cell] *
                     impl.continuity_absolute[cell];
        }
        if (!std::isfinite(sums[0]) || !std::isfinite(sums[1])) {
          throw runtime::Error("continuity norm is non-finite");
        }
      });
  impl.mpi->allreduce_fp64_in_place(sums, 2U,
                                    runtime::Fp64ReductionOperation::sum);
  const double raw_l2 = std::sqrt(sums[0]);
  const double scale_l2 = std::sqrt(sums[1]);
  const double normalized =
      scale_l2 == 0.0
          ? (raw_l2 == 0.0 ? 0.0 : std::numeric_limits<double>::infinity())
          : raw_l2 / scale_l2;
  return {raw_l2, normalized};
}

double normalized_l2(double residual_square, double scale_square) noexcept {
  const double residual = std::sqrt(residual_square);
  const double scale = std::sqrt(scale_square);
  return scale == 0.0
             ? (residual == 0.0 ? 0.0
                                : std::numeric_limits<double>::infinity())
             : residual / scale;
}

class CompensatedSum final {
public:
  void add(double value) noexcept {
    const double candidate = sum_ + value;
    if (std::abs(sum_) >= std::abs(value)) {
      correction_ += (sum_ - candidate) + value;
    } else {
      correction_ += (value - candidate) + sum_;
    }
    sum_ = candidate;
  }

  double value() const noexcept { return sum_ + correction_; }

private:
  double sum_{};
  double correction_{};
};

struct ConservationParts final {
  CompensatedSum quantity_nm1;
  CompensatedSum quantity_n;
  CompensatedSum quantity_np1;
  CompensatedSum absolute_quantity_nm1;
  CompensatedSum absolute_quantity_n;
  CompensatedSum absolute_quantity_np1;
  CompensatedSum signed_boundary_flux;
  CompensatedSum absolute_boundary_flux;
};

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
void replace_conservation_parts(
    ConservationParts &parts,
    const std::array<std::atomic<double>, 8> &values) noexcept {
  parts = {};
  parts.quantity_nm1.add(values[0].load(std::memory_order_relaxed));
  parts.quantity_n.add(values[1].load(std::memory_order_relaxed));
  parts.quantity_np1.add(values[2].load(std::memory_order_relaxed));
  parts.absolute_quantity_nm1.add(
      values[3].load(std::memory_order_relaxed));
  parts.absolute_quantity_n.add(values[4].load(std::memory_order_relaxed));
  parts.absolute_quantity_np1.add(
      values[5].load(std::memory_order_relaxed));
  parts.signed_boundary_flux.add(values[6].load(std::memory_order_relaxed));
  parts.absolute_boundary_flux.add(
      values[7].load(std::memory_order_relaxed));
}

void inject_finite_terms_that_overflow(ConservationParts &parts) noexcept {
  const double maximum = std::numeric_limits<double>::max();
  parts.signed_boundary_flux.add(maximum);
  parts.signed_boundary_flux.add(maximum);
}
#endif

struct ConservationReduction final {
  runtime::CollectiveStatus status;
  std::array<double, 8> values{};
};

ConservationReduction
reduce_conservation_parts(const runtime::MpiContext &mpi,
                          const ConservationParts &parts) {
  std::array<double, 8> values{
      parts.quantity_nm1.value(), parts.quantity_n.value(),
      parts.quantity_np1.value(), parts.absolute_quantity_nm1.value(),
      parts.absolute_quantity_n.value(),
      parts.absolute_quantity_np1.value(),
      parts.signed_boundary_flux.value(),
      parts.absolute_boundary_flux.value()};
  const bool local_finite =
      std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
      });
  auto status = runtime::collective_status(
      mpi, local_finite,
      "local conservation aggregate is non-finite");
  if (!status.ok)
    return {std::move(status), {}};
  mpi.allreduce_fp64_in_place(values.data(), values.size(),
                              runtime::Fp64ReductionOperation::sum);
  const bool global_finite =
      std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
      });
  status = runtime::collective_status(
      mpi, global_finite,
      "global conservation aggregate is non-finite");
  return {std::move(status), values};
}

struct ConservationDiagnosticValues final {
  double quantity_nm1{};
  double quantity_n{};
  double quantity_np1{};
  double signed_boundary_flux{};
  double absolute_boundary_flux{};
  double boundary_integral{};
  double history_correction{};
  double raw_defect{};
  double relative{};
};

ConservationDiagnosticValues make_conservation_diagnostic(
    const std::array<double, 8> &values,
    const MomentumTimeStencil &stencil,
    bool enable_cancellation_normalization,
    double defect_perturbation = 0.0) {
  const double history_correction =
      stencil.alpha2 * (values[0] - values[1]) / stencil.alpha0;
  const double boundary_integral =
      stencil.dt_s * values[6] / stencil.alpha0;
  const double absolute_boundary_scale =
      stencil.dt_s * values[7] / stencil.alpha0;
  const double raw_defect = values[2] - values[1] + boundary_integral +
                            history_correction + defect_perturbation;
  double denominator =
      std::max({std::abs(values[1]), std::abs(values[2]),
                absolute_boundary_scale, std::abs(history_correction), 0.0});
  const double cancellation_scale =
      std::max({values[4], values[5],
                stencil.order == MomentumTimeOrder::bdf2 ? values[3] : 0.0});
  if (enable_cancellation_normalization && cancellation_scale > 0.0 &&
      denominator <= 64.0 * std::numeric_limits<double>::epsilon() *
                         cancellation_scale) {
    denominator = cancellation_scale;
  }
  denominator =
      std::max(denominator, std::numeric_limits<double>::min());
  return {values[0],
          values[1],
          values[2],
          values[6],
          values[7],
          boundary_integral,
          history_correction,
          raw_defect,
          std::abs(raw_defect) / denominator};
}

template <class FlowImplementation>
runtime::CollectiveStatus assess_final_momentum(
    FlowImplementation &impl, FlowState &state, double rho_ref, double mu,
    const runtime::FaceFieldView<const double> *dynamic_viscosity_by_face,
    const MomentumTimeStencil &stencil, StepAttemptReport &report) {
  auto &committed =
      detail::FlowStateSolverAccess::layer(state, FlowLayer::committed);
  auto &history =
      detail::FlowStateSolverAccess::layer(state, FlowLayer::history);
  auto &trial = detail::FlowStateSolverAccess::layer(state, FlowLayer::trial);
  const auto &registry = detail::FlowStateSolverAccess::registry(state);
  const auto &access = detail::FlowStateSolverAccess::access(state);
  const auto fields = state.fields();
  auto &residual_n = impl.momentum_n;
  auto &residual_nm1 = impl.momentum_nm1;
  assemble_spatial_residual(impl, registry, access, committed, fields, mu,
                            dynamic_viscosity_by_face, residual_n);
  if (stencil.order == MomentumTimeOrder::bdf2) {
    assemble_spatial_residual(impl, registry, access, history, fields, mu,
                              dynamic_viscosity_by_face, residual_nm1);
  }
  impl.halo->exchange(trial, fields.mechanical_pressure);

  std::array<double, 6> norm_sums{};
  std::array<ConservationParts, 3> conservation;
  synchronized_local_phase(
      *impl.mpi, StepFailureReason::final_momentum_residual, true,
      "final momentum residual reconstruction failed", [&] {
        const auto velocity_trial = trial.acquire_read<double>(
            access, kStatePhase, kStateActor, fields.velocity);
        const auto velocity_n = committed.acquire_read<double>(
            access, kStatePhase, kStateActor, fields.velocity);
        const auto velocity_nm1 = history.acquire_read<double>(
            access, kStatePhase, kStateActor, fields.velocity);
        const auto pressure_trial = trial.acquire_read<double>(
            access, kStatePhase, kStateActor, fields.mechanical_pressure);
        auto pressure_gradient = impl.scratch.storage->template acquire_write<
            double>(*impl.scratch.access, kScratchPhase, kScratchActor,
                    impl.scratch.pressure_gradient);
        compute_pressure_gradient(
            *impl.topology, *impl.geometry, pressure_trial, pressure_gradient,
            impl.pressure_gradient_sums,
            [&](LocalFaceId face, double owner_value) {
              const auto patch = impl.topology->patch_id(face);
              if (!patch.has_value())
                return owner_value;
              return impl.boundaries->evaluate_pressure(*patch, owner_value)
                  .face;
            });
        const auto pressure_gradient_read =
            impl.scratch.storage->template acquire_read<double>(
                *impl.scratch.access, kScratchPhase, kScratchActor,
                impl.scratch.pressure_gradient);
        impl.fvm.physical_boundary_pressure_contributions(
            *impl.boundaries, pressure_trial, impl.pressure_boundary);
        if ((stencil.order == MomentumTimeOrder::bdf2 &&
             residual_n.boundary.size() != residual_nm1.boundary.size()) ||
            residual_n.boundary.size() != impl.pressure_boundary.size()) {
          throw runtime::Error(
              "physical momentum boundary sets differ");
        }
        for (std::size_t boundary_index = 0;
             boundary_index < residual_n.boundary.size(); ++boundary_index) {
          const auto &current = residual_n.boundary[boundary_index];
          const auto *previous =
              stencil.order == MomentumTimeOrder::bdf2
                  ? &residual_nm1.boundary[boundary_index]
                  : nullptr;
          const auto &pressure = impl.pressure_boundary[boundary_index];
          if ((previous != nullptr &&
               current.global_face_id != previous->global_face_id) ||
              current.global_face_id != pressure.global_face_id) {
            throw runtime::Error(
                "physical momentum boundary identities differ");
          }
          for (std::size_t component_index = 0; component_index < 3U;
               ++component_index) {
            const double convection =
                previous == nullptr
                    ? current.convective[component_index]
                    : 2.0 * current.convective[component_index] -
                          previous->convective[component_index];
            const double viscosity =
                previous == nullptr
                    ? current.viscous[component_index]
                    : 2.0 * current.viscous[component_index] -
                          previous->viscous[component_index];
            const double pressure_value = pressure.pressure[component_index];
            const double effective_flux =
                convection + viscosity + pressure_value;
            conservation[component_index].signed_boundary_flux.add(
                effective_flux);
            conservation[component_index].absolute_boundary_flux.add(
                std::abs(effective_flux));
          }
        }
        const auto owned = impl.topology->owned_global_box();
        for (LocalCellId cell = 0; cell < impl.topology->owned_cell_count();
             ++cell) {
          const StructuredIndex index =
              map_cell(impl.topology->global_cell(cell), owned,
                       impl.topology->global_extent());
          const double volume = impl.geometry->cell_volume_m3(cell);
          for (std::size_t component_index = 0; component_index < 3U;
               ++component_index) {
            const double trial_time =
                (rho_ref * volume / stencil.dt_s) * stencil.alpha0 *
                velocity_trial(index.i, index.j, index.k,
                               static_cast<int>(component_index));
            const double current_time =
                (rho_ref * volume / stencil.dt_s) * stencil.alpha1 *
                velocity_n(index.i, index.j, index.k,
                           static_cast<int>(component_index));
            const double history_time =
                (rho_ref * volume / stencil.dt_s) * stencil.alpha2 *
                velocity_nm1(index.i, index.j, index.k,
                             static_cast<int>(component_index));
            const std::size_t offset = cell * 3U + component_index;
            const double convection_n =
                (stencil.order == MomentumTimeOrder::backward_euler ? 1.0
                                                                     : 2.0) *
                residual_n.convection[offset];
            const double convection_nm1 =
                stencil.order == MomentumTimeOrder::backward_euler
                    ? 0.0
                    : -residual_nm1.convection[offset];
            const double viscosity_n =
                (stencil.order == MomentumTimeOrder::backward_euler ? 1.0
                                                                     : 2.0) *
                residual_n.viscosity[offset];
            const double viscosity_nm1 =
                stencil.order == MomentumTimeOrder::backward_euler
                    ? 0.0
                    : -residual_nm1.viscosity[offset];
            const double pressure =
                volume * pressure_gradient_read(
                             index.i, index.j, index.k,
                             static_cast<int>(component_index));
            const std::array<double, 8> terms{
                trial_time,       current_time,    history_time,
                convection_n,     convection_nm1, viscosity_n,
                viscosity_nm1,    pressure};
            double raw = 0.0;
            double scale = 0.0;
            for (std::size_t term = 0; term < terms.size(); ++term) {
              raw += terms[term];
              scale += std::abs(terms[term]);
            }
            const double quantity_n = rho_ref * volume *
                                      velocity_n(index.i, index.j, index.k,
                                                 static_cast<int>(
                                                     component_index));
            const double quantity_nm1 =
                rho_ref * volume *
                velocity_nm1(index.i, index.j, index.k,
                             static_cast<int>(component_index));
            const double quantity_trial =
                rho_ref * volume *
                velocity_trial(index.i, index.j, index.k,
                               static_cast<int>(component_index));
            conservation[component_index].quantity_nm1.add(quantity_nm1);
            conservation[component_index].quantity_n.add(quantity_n);
            conservation[component_index].quantity_np1.add(quantity_trial);
            conservation[component_index].absolute_quantity_nm1.add(
                std::abs(quantity_nm1));
            conservation[component_index].absolute_quantity_n.add(
                std::abs(quantity_n));
            conservation[component_index].absolute_quantity_np1.add(
                std::abs(quantity_trial));
            norm_sums[component_index * 2U] += raw * raw;
            norm_sums[component_index * 2U + 1U] += scale * scale;
            if (!std::isfinite(raw) || !std::isfinite(scale) ||
                !std::isfinite(quantity_nm1) || !std::isfinite(quantity_n) ||
                !std::isfinite(quantity_trial)) {
              throw runtime::Error(
                  "final momentum assessment is non-finite");
            }
          }
        }
      });
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  for (std::size_t component_index = 0; component_index < 3U;
       ++component_index) {
    if (final_momentum_norm_armed[component_index].load(
            std::memory_order_relaxed)) {
      norm_sums[component_index * 2U] =
          final_momentum_norm_residual_square[component_index].load(
              std::memory_order_relaxed);
      norm_sums[component_index * 2U + 1U] =
          final_momentum_norm_scale_square[component_index].load(
              std::memory_order_relaxed);
    }
  }
#endif
  const bool local_norms_finite =
      std::all_of(norm_sums.begin(), norm_sums.end(), [](double value) {
        return value >= 0.0 && std::isfinite(value);
      });
  const auto local_norm_status = runtime::collective_status(
      *impl.mpi, local_norms_finite,
      "final momentum local norm is non-finite");
  if (!local_norm_status.ok)
    return local_norm_status;
  impl.mpi->allreduce_fp64_in_place(norm_sums.data(), norm_sums.size(),
                                    runtime::Fp64ReductionOperation::sum);
  bool global_ok = true;
  for (std::size_t component_index = 0; component_index < 3U;
       ++component_index) {
    report.final_momentum_normalized_l2[component_index] = normalized_l2(
        norm_sums[component_index * 2U],
        norm_sums[component_index * 2U + 1U]);
    global_ok = global_ok &&
                std::isfinite(
                    report.final_momentum_normalized_l2[component_index]) &&
                report.final_momentum_normalized_l2[component_index] <=
                    kFinalEquationTolerance;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (momentum_conservation_parts_armed[component_index].load(
            std::memory_order_relaxed)) {
      replace_conservation_parts(
          conservation[component_index],
          momentum_conservation_parts_values[component_index]);
    }
    if (momentum_conservation_overflow_component.load(
            std::memory_order_relaxed) == component_index) {
      inject_finite_terms_that_overflow(conservation[component_index]);
    }
#endif
    const auto reduced =
        reduce_conservation_parts(*impl.mpi, conservation[component_index]);
    if (!reduced.status.ok)
      return reduced.status;
    const auto diagnostic =
        make_conservation_diagnostic(reduced.values, stencil, true);
    report.final_momentum_relative_conservation_defect[component_index] =
        diagnostic.relative;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    last_momentum_conservation[component_index] = {
        diagnostic.quantity_nm1,
        diagnostic.quantity_n,
        diagnostic.quantity_np1,
        diagnostic.signed_boundary_flux,
        diagnostic.absolute_boundary_flux,
        diagnostic.boundary_integral,
        diagnostic.history_correction,
        diagnostic.raw_defect,
        diagnostic.relative};
#endif
  }
  return runtime::collective_status(
      *impl.mpi, global_ok,
      "final momentum residual exceeds tolerance");
}

template <class FlowImplementation>
runtime::CollectiveStatus assess_final_transport(
    FlowImplementation &impl, FlowState &state, double rho_ref,
    const runtime::FaceFieldView<const double> *sgs_viscosity_by_face,
    const les::WaleControl *wale_control,
    const MomentumTimeStencil &stencil, StepAttemptReport &report) {
  auto &committed =
      detail::FlowStateSolverAccess::layer(state, FlowLayer::committed);
  auto &history =
      detail::FlowStateSolverAccess::layer(state, FlowLayer::history);
  auto &trial = detail::FlowStateSolverAccess::layer(state, FlowLayer::trial);
  const auto &access = detail::FlowStateSolverAccess::access(state);
  report.final_transport_normalized_l2.assign(
      impl.transport.size(), std::numeric_limits<double>::infinity());
  report.final_transport_relative_conservation_defect.assign(
      impl.transport.size(), std::numeric_limits<double>::infinity());
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  last_transport_conservation.assign(impl.transport.size(), {});
#endif
  for (std::size_t transport_index = 0;
       transport_index < impl.transport.size(); ++transport_index) {
    const auto &item = impl.transport[transport_index];
    auto &spatial_n = impl.transport_n[transport_index];
    auto &spatial_nm1 = impl.transport_nm1[transport_index];
    assemble_transport_spatial_residual(
        impl, state, committed, item, sgs_viscosity_by_face, wale_control,
        spatial_n);
    if (stencil.order == MomentumTimeOrder::bdf2) {
      assemble_transport_spatial_residual(
          impl, state, history, item, sgs_viscosity_by_face, wale_control,
          spatial_nm1);
    }
    double norm_sums[2]{};
    ConservationParts conservation;
    synchronized_local_phase(
        *impl.mpi, StepFailureReason::final_transport_residual, true,
        "final transport residual reconstruction failed", [&] {
          const auto values_trial = trial.acquire_read<double>(
              access, kStatePhase, kStateActor, item.field);
          const auto values_n = committed.acquire_read<double>(
              access, kStatePhase, kStateActor, item.field);
          const auto values_nm1 = history.acquire_read<double>(
              access, kStatePhase, kStateActor, item.field);
          if (stencil.order == MomentumTimeOrder::bdf2 &&
              spatial_n.boundary.size() != spatial_nm1.boundary.size()) {
            throw runtime::Error(
                "physical transport boundary sets differ");
          }
          for (std::size_t boundary_index = 0;
               boundary_index < spatial_n.boundary.size(); ++boundary_index) {
            const auto &current = spatial_n.boundary[boundary_index];
            const auto *previous =
                stencil.order == MomentumTimeOrder::bdf2
                    ? &spatial_nm1.boundary[boundary_index]
                    : nullptr;
            if (previous != nullptr &&
                current.global_face_id != previous->global_face_id) {
              throw runtime::Error(
                  "physical transport boundary identities differ");
            }
            const double convection =
                previous == nullptr
                    ? current.convective
                    : 2.0 * current.convective - previous->convective;
            const double diffusion =
                previous == nullptr
                    ? current.diffusive
                    : 2.0 * current.diffusive - previous->diffusive;
            const double effective_flux = convection + diffusion;
            conservation.signed_boundary_flux.add(effective_flux);
            conservation.absolute_boundary_flux.add(
                std::abs(effective_flux));
          }
          const auto owned = impl.topology->owned_global_box();
          for (LocalCellId cell = 0; cell < impl.topology->owned_cell_count();
               ++cell) {
            const StructuredIndex index =
                map_cell(impl.topology->global_cell(cell), owned,
                         impl.topology->global_extent());
            const double mass =
                rho_ref * impl.geometry->cell_volume_m3(cell) / stencil.dt_s;
            const std::array<double, 7> terms{
                mass * stencil.alpha0 *
                    values_trial(index.i, index.j, index.k, 0),
                mass * stencil.alpha1 *
                    values_n(index.i, index.j, index.k, 0),
                mass * stencil.alpha2 *
                    values_nm1(index.i, index.j, index.k, 0),
                (stencil.order == MomentumTimeOrder::backward_euler ? 1.0
                                                                     : 2.0) *
                    spatial_n.convection[cell],
                stencil.order == MomentumTimeOrder::backward_euler
                    ? 0.0
                    : -spatial_nm1.convection[cell],
                (stencil.order == MomentumTimeOrder::backward_euler ? 1.0
                                                                     : 2.0) *
                    spatial_n.diffusion[cell],
                stencil.order == MomentumTimeOrder::backward_euler
                    ? 0.0
                    : -spatial_nm1.diffusion[cell]};
            double raw = 0.0;
            double scale = 0.0;
            for (std::size_t term = 0; term < terms.size(); ++term) {
              raw += terms[term];
              scale += std::abs(terms[term]);
            }
            const double quantity_nm1 =
                rho_ref * impl.geometry->cell_volume_m3(cell) *
                values_nm1(index.i, index.j, index.k, 0);
            const double quantity_n =
                rho_ref * impl.geometry->cell_volume_m3(cell) *
                values_n(index.i, index.j, index.k, 0);
            const double quantity_trial =
                rho_ref * impl.geometry->cell_volume_m3(cell) *
                values_trial(index.i, index.j, index.k, 0);
            conservation.quantity_nm1.add(quantity_nm1);
            conservation.quantity_n.add(quantity_n);
            conservation.quantity_np1.add(quantity_trial);
            conservation.absolute_quantity_nm1.add(std::abs(quantity_nm1));
            conservation.absolute_quantity_n.add(std::abs(quantity_n));
            conservation.absolute_quantity_np1.add(
                std::abs(quantity_trial));
            norm_sums[0] += raw * raw;
            norm_sums[1] += scale * scale;
            if (!std::isfinite(raw) || !std::isfinite(scale) ||
                !std::isfinite(quantity_nm1) || !std::isfinite(quantity_n) ||
                !std::isfinite(quantity_trial)) {
              throw runtime::Error(
                  "final transport assessment is non-finite");
            }
          }
        });
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (final_transport_norm_index.load(std::memory_order_relaxed) ==
        transport_index) {
      norm_sums[0] = final_transport_norm_residual_square.load(
          std::memory_order_relaxed);
      norm_sums[1] =
          final_transport_norm_scale_square.load(std::memory_order_relaxed);
    }
#endif
    const bool local_norms_finite =
        norm_sums[0] >= 0.0 && std::isfinite(norm_sums[0]) &&
        norm_sums[1] >= 0.0 && std::isfinite(norm_sums[1]);
    const auto local_norm_status = runtime::collective_status(
        *impl.mpi, local_norms_finite,
        "final transport local norm is non-finite");
    if (!local_norm_status.ok)
      return local_norm_status;
    impl.mpi->allreduce_fp64_in_place(
        norm_sums, 2U, runtime::Fp64ReductionOperation::sum);
    report.final_transport_normalized_l2[transport_index] =
        normalized_l2(norm_sums[0], norm_sums[1]);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (transport_conservation_overflow_index.load(std::memory_order_relaxed) ==
        transport_index) {
      inject_finite_terms_that_overflow(conservation);
    }
#endif
    const auto reduced = reduce_conservation_parts(*impl.mpi, conservation);
    if (!reduced.status.ok)
      return reduced.status;
    const auto diagnostic =
        make_conservation_diagnostic(reduced.values, stencil, true);
    report.final_transport_relative_conservation_defect[transport_index] =
        diagnostic.relative;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    last_transport_conservation[transport_index] = {
        diagnostic.quantity_nm1,
        diagnostic.quantity_n,
        diagnostic.quantity_np1,
        diagnostic.signed_boundary_flux,
        diagnostic.absolute_boundary_flux,
        diagnostic.boundary_integral,
        diagnostic.history_correction,
        diagnostic.raw_defect,
        diagnostic.relative};
#endif
    const bool global_field_ok =
        std::isfinite(report.final_transport_normalized_l2[transport_index]) &&
        report.final_transport_normalized_l2[transport_index] <=
            kFinalEquationTolerance;
    const auto status = runtime::collective_status(
        *impl.mpi, global_field_ok,
        "final transport residual exceeds tolerance");
    if (!status.ok)
      return status;
  }
  return runtime::collective_status(*impl.mpi, true,
                                    "final transport residual passed");
}

template <class FlowImplementation>
runtime::CollectiveStatus assess_final_conservation(
    FlowImplementation &impl, FlowState &state, double rho_ref,
    const MomentumTimeStencil &stencil, StepAttemptReport &report) {
  static_cast<void>(rho_ref);
  auto &trial = detail::FlowStateSolverAccess::layer(state, FlowLayer::trial);
  auto &committed =
      detail::FlowStateSolverAccess::layer(state, FlowLayer::committed);
  auto &history =
      detail::FlowStateSolverAccess::layer(state, FlowLayer::history);
  const auto &access = detail::FlowStateSolverAccess::access(state);
  const auto fields = state.fields();
  ConservationParts conservation;
  synchronized_local_phase(
      *impl.mpi, StepFailureReason::final_conservation_defect, true,
      "final mass conservation reconstruction failed", [&] {
        const auto density_nm1 = history.acquire_read<double>(
            access, kStatePhase, kStateActor, fields.density);
        const auto density_n = committed.acquire_read<double>(
            access, kStatePhase, kStateActor, fields.density);
        const auto density_trial = trial.acquire_read<double>(
            access, kStatePhase, kStateActor, fields.density);
        const auto owned = impl.topology->owned_global_box();
        for (LocalCellId cell = 0; cell < impl.topology->owned_cell_count();
             ++cell) {
          const StructuredIndex index =
              map_cell(impl.topology->global_cell(cell), owned,
                       impl.topology->global_extent());
          const double volume = impl.geometry->cell_volume_m3(cell);
          const double mass_nm1 =
              volume * density_nm1(index.i, index.j, index.k, 0);
          const double mass_n =
              volume * density_n(index.i, index.j, index.k, 0);
          const double mass_trial =
              volume * density_trial(index.i, index.j, index.k, 0);
          conservation.quantity_nm1.add(mass_nm1);
          conservation.quantity_n.add(mass_n);
          conservation.quantity_np1.add(mass_trial);
        }
        const auto flux = trial.acquire_face_read<double>(
            access, kStatePhase, kStateActor, fields.face_mass_flux);
        for (LocalFaceId face = 0; face < impl.topology->local_face_count();
             ++face) {
          const auto patch = impl.topology->patch_id(face);
          if (!patch.has_value() ||
              impl.boundaries->patch(*patch).kind() ==
                  boundary::BoundaryKind::periodic) {
            continue;
          }
          const LocalCellId owner = impl.topology->owner(face);
          if (impl.topology->cell_ownership(owner) != EntityOwnership::owned) {
            throw runtime::Error(
                "physical mass boundary face is not owner-owned");
          }
          conservation.signed_boundary_flux.add(flux(face, 0));
          conservation.absolute_boundary_flux.add(std::abs(flux(face, 0)));
        }
      });
  const auto reduced = reduce_conservation_parts(*impl.mpi, conservation);
  if (!reduced.status.ok)
    return reduced.status;
  double perturbation = 0.0;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  if (force_final_conservation_failure.load(std::memory_order_relaxed)) {
    perturbation += 1.0e-4;
  }
  perturbation +=
      final_mass_defect_perturbation.load(std::memory_order_relaxed);
#endif
  const auto diagnostic =
      make_conservation_diagnostic(reduced.values, stencil, false,
                                   perturbation);
  report.final_mass_relative_conservation_defect =
      diagnostic.relative;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  last_mass_conservation = {diagnostic.quantity_nm1,
                            diagnostic.quantity_n,
                            diagnostic.quantity_np1,
                            diagnostic.signed_boundary_flux,
                            diagnostic.absolute_boundary_flux,
                            diagnostic.boundary_integral,
                            diagnostic.history_correction,
                            diagnostic.raw_defect,
                            diagnostic.relative};
#endif
  bool local_ok =
      std::isfinite(report.final_mass_relative_conservation_defect) &&
      report.final_mass_relative_conservation_defect <=
          kFinalConservationTolerance;
  for (double value : report.final_momentum_relative_conservation_defect) {
    local_ok = local_ok && std::isfinite(value) &&
               value <= kFinalConservationTolerance;
  }
  for (double value : report.final_transport_relative_conservation_defect) {
    local_ok = local_ok && std::isfinite(value) &&
               value <= kFinalConservationTolerance;
  }
  return runtime::collective_status(
      *impl.mpi, local_ok,
      "final conservation defect exceeds tolerance");
}

} // namespace

StepAttemptReport FixedStepConstantDensityFlow::attempt(
    FlowState &state, double rho_ref, double mu,
    const MomentumTimeStencil &stencil,
    const linear::SolveControl &momentum_control,
    const linear::SolveControl &pressure_control) const {
  return attempt_impl(state, rho_ref, mu, stencil, momentum_control,
                      pressure_control, nullptr, nullptr);
}

StepAttemptReport FixedStepConstantDensityFlow::attempt_with_wale(
    FlowState &state, double rho_ref, double molecular_mu,
    const MomentumTimeStencil &stencil,
    const linear::SolveControl &momentum_control,
    const linear::SolveControl &pressure_control, const les::WaleModel &wale,
    les::WaleSummary &summary) const {
  summary = {};
  return attempt_impl(state, rho_ref, molecular_mu, stencil, momentum_control,
                      pressure_control, &wale, &summary);
}

StepAttemptReport FixedStepConstantDensityFlow::attempt_impl(
    FlowState &state, double rho_ref, double mu,
    const MomentumTimeStencil &stencil,
    const linear::SolveControl &momentum_control,
    const linear::SolveControl &pressure_control,
    const les::WaleModel *wale, les::WaleSummary *wale_summary) const {
  StepAttemptReport report = base_report(stencil.dt_s);
  bool active = false;
  try {
    bool local_valid = false;
    try {
      const auto expected = make_momentum_time_stencil(
          stencil.order, stencil.dt_s, stencil.previous_dt_s);
      const auto layout = state.layer(FlowLayer::committed).layout_set();
      const auto local_extent = impl_->decomposition->local_extent();
      const auto fields = state.fields();
      const auto metadata = state.metadata();
      local_valid =
          rho_ref > 0.0 && std::isfinite(rho_ref) && mu >= 0.0 &&
          std::isfinite(mu) &&
          ((wale == nullptr && wale_summary == nullptr) ||
           (wale != nullptr && wale_summary != nullptr)) &&
          impl_->transport_specs_valid &&
          !state.attempt_active() &&
          layout.cell_interior_extent.x == local_extent.x &&
          layout.cell_interior_extent.y == local_extent.y &&
          layout.cell_interior_extent.z == local_extent.z &&
          layout.face_count == impl_->topology->local_face_count() &&
          expected.alpha0 == stencil.alpha0 &&
          expected.alpha1 == stencil.alpha1 &&
          expected.alpha2 == stencil.alpha2 &&
          (stencil.order != MomentumTimeOrder::bdf2 ||
           stencil.previous_dt_s == metadata.dt_s) &&
          fields.transported_cell_fields.size() == impl_->transport.size();
      for (std::size_t index = 0;
           local_valid && index < impl_->transport.size(); ++index) {
        local_valid = impl_->transport[index].field ==
                      fields.transported_cell_fields[index];
        if (local_valid) {
          const auto &descriptor =
              state.solver_registry().descriptor(impl_->transport[index].field);
          local_valid =
              descriptor.space == runtime::FunctionSpace::cell_average &&
              descriptor.scalar_type == runtime::ScalarType::float64 &&
              descriptor.components == 1U && descriptor.ghost_width >= 2;
        }
      }
      if (local_valid) {
        const auto &state_access = state.solver_access_plan();
        const auto committed_density =
            state.layer(FlowLayer::committed).acquire_read<double>(
                state_access, kStatePhase, kStateActor, fields.density);
        const auto history_density =
            state.layer(FlowLayer::history).acquire_read<double>(
                state_access, kStatePhase, kStateActor, fields.density);
        for (int k = 0; local_valid && k < local_extent.z; ++k) {
          for (int j = 0; local_valid && j < local_extent.y; ++j) {
            for (int i = 0; local_valid && i < local_extent.x; ++i) {
              local_valid = committed_density(i, j, k, 0) == rho_ref &&
                            history_density(i, j, k, 0) == rho_ref;
            }
          }
        }
      }
    } catch (const runtime::detail::MpiOperationError &) {
      throw;
    } catch (...) {
      local_valid = false;
    }
    const auto validation = runtime::collective_status(
        *impl_->mpi, local_valid, "fixed-step input is invalid");
    if (!validation.ok) {
      return fatal_failure(report, StepFailureReason::invalid_input,
                           validation.failing_rank);
    }
    const auto metadata = state.metadata();
    const auto wale_control =
        wale == nullptr ? les::WaleControl{} : wale->control();
    const std::array<double, 18> core_inputs{
        rho_ref,
        mu,
        wale == nullptr ? 0.0 : 1.0,
        wale_control.coefficient,
        wale_control.turbulent_prandtl,
        wale_control.turbulent_schmidt,
        static_cast<double>(static_cast<std::uint8_t>(stencil.order)),
        stencil.dt_s,
        stencil.previous_dt_s,
        stencil.alpha0,
        stencil.alpha1,
        stencil.alpha2,
        low_u32(metadata.step),
        high_u32(metadata.step),
        metadata.time_s,
        metadata.dt_s,
        metadata.previous_dt_s,
        static_cast<double>(static_cast<std::uint8_t>(metadata.order))};
    const auto core_agreement =
        agree_fp64_inputs(*impl_->mpi, core_inputs,
                          "fixed-step rank inputs are inconsistent");
    if (!core_agreement.ok) {
      return fatal_failure(report, StepFailureReason::invalid_input,
                           core_agreement.failing_rank);
    }
    const auto momentum_control_agreement = agree_solve_control(
        *impl_->mpi, momentum_control,
        "momentum solve control is invalid",
        "momentum solve controls are inconsistent");
    if (!momentum_control_agreement.ok) {
      return fatal_failure(report, StepFailureReason::invalid_input,
                           momentum_control_agreement.failing_rank);
    }
    const auto pressure_control_agreement = agree_solve_control(
        *impl_->mpi, pressure_control,
        "pressure solve control is invalid",
        "pressure solve controls are inconsistent");
    if (!pressure_control_agreement.ok) {
      return fatal_failure(report, StepFailureReason::invalid_input,
                           pressure_control_agreement.failing_rank);
    }
    const auto transport_count =
        static_cast<std::uint64_t>(impl_->transport.size());
    const std::array<double, 2> transport_count_input{
        low_u32(transport_count), high_u32(transport_count)};
    const auto transport_count_agreement =
        agree_fp64_inputs(*impl_->mpi, transport_count_input,
                          "transported-field counts are inconsistent");
    if (!transport_count_agreement.ok) {
      return fatal_failure(report, StepFailureReason::invalid_input,
                           transport_count_agreement.failing_rank);
    }
    for (const auto &transport : impl_->transport) {
      const auto scalar_index =
          static_cast<std::uint64_t>(transport.quantity.scalar_index);
      const std::array<double, 5> transport_input{
          static_cast<double>(transport.field),
          static_cast<double>(
              static_cast<std::uint8_t>(transport.quantity.kind)),
          low_u32(scalar_index), high_u32(scalar_index),
          transport.diffusivity_kg_per_m_s};
      const auto transport_agreement = agree_fp64_inputs(
          *impl_->mpi, transport_input,
          "transported-field inputs are inconsistent");
      if (!transport_agreement.ok) {
        return fatal_failure(report, StepFailureReason::invalid_input,
                             transport_agreement.failing_rank);
      }
    }

    bool begin_ok = true;
    try {
      state.begin_attempt();
      active = true;
    } catch (const runtime::detail::MpiOperationError &) {
      throw;
    } catch (...) {
      begin_ok = false;
    }
    const auto begin_status = runtime::collective_status(
        *impl_->mpi, begin_ok, "trial initialization failed");
    if (!begin_status.ok) {
      if (active) {
        state.rollback_attempt();
        active = false;
      }
      return fatal_failure(report, StepFailureReason::invalid_input,
                           begin_status.failing_rank);
    }
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    inject_attempt_stage_failure(*impl_->mpi,
                                 kAttemptFailureAfterBegin);
#endif
    exchange_accepted_layers(*impl_, state);
    auto &committed = state.solver_layer(FlowLayer::committed);
    auto &history = state.solver_layer(FlowLayer::history);
    auto &trial = state.solver_layer(FlowLayer::trial);
    const auto &registry = state.solver_registry();
    const auto &access = state.solver_access_plan();
    const auto fields = state.fields();
    std::optional<detail::BodyFittedWaleAttemptData> wale_data;
    std::uint64_t wale_evaluation_count{};
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::non_finite_trial, true,
        "WALE attempt preparation failed", [&] {
          if (wale != nullptr) {
            const auto rho_attempt = committed.acquire_read<double>(
                access, kStatePhase, kStateActor, fields.density);
            wale_data.emplace(detail::prepare_body_fitted_wale_attempt(
                *impl_, state, *wale, mu, stencil, rho_attempt,
                wale_evaluation_count, kStatePhase, kStateActor,
                kScratchPhase, kScratchActor));
          }
        });
    const runtime::FaceFieldView<const double> *effective_viscosity =
        wale_data.has_value()
            ? &*wale_data->effective_dynamic_viscosity_by_face
            : nullptr;
    const runtime::FaceFieldView<const double> *sgs_viscosity =
        wale_data.has_value() ? &*wale_data->sgs_dynamic_viscosity_by_face
                              : nullptr;
    const les::WaleControl *active_wale_control =
        wale_data.has_value() ? &wale_control : nullptr;
    auto &residual_n = impl_->momentum_n;
    auto &residual_nm1 = impl_->momentum_nm1;
    assemble_spatial_residual(*impl_, registry, access, committed, fields, mu,
                              effective_viscosity, residual_n);
    runtime::FieldStorage *history_flux_storage = nullptr;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (momentum_assembly_mutation.load(std::memory_order_relaxed) ==
        kMomentumMutationReplaceHistoryFlux) {
      history_flux_storage = &trial;
    }
#endif
    if (stencil.order == MomentumTimeOrder::bdf2) {
      assemble_spatial_residual(*impl_, registry, access, history, fields, mu,
                                effective_viscosity, residual_nm1,
                                history_flux_storage);
    }
    std::optional<runtime::FieldView<const double>> velocity_n;
    std::optional<runtime::FieldView<const double>> velocity_nm1;
    std::optional<runtime::FieldView<const double>> pressure_n;
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::non_finite_trial, true,
        "accepted field acquisition failed", [&] {
          velocity_n.emplace(committed.acquire_read<double>(
              access, kStatePhase, kStateActor, fields.velocity));
          velocity_nm1.emplace(history.acquire_read<double>(
              access, kStatePhase, kStateActor, fields.velocity));
          pressure_n.emplace(committed.acquire_read<double>(
              access, kStatePhase, kStateActor, fields.mechanical_pressure));
        });
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::non_finite_trial, true,
        "pressure-gradient derivation failed", [&] {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
          if (force_local_derived_failure.load(std::memory_order_relaxed)) {
            throw runtime::Error(
                "injected local derived-value failure");
          }
#endif
          auto pressure_gradient =
              impl_->scratch.storage->acquire_write<double>(
                  *impl_->scratch.access, kScratchPhase, kScratchActor,
                  impl_->scratch.pressure_gradient);
          compute_pressure_gradient(
              *impl_->topology, *impl_->geometry, *pressure_n,
              pressure_gradient, impl_->pressure_gradient_sums,
              [&](LocalFaceId face, double owner_value) {
                const auto patch = impl_->topology->patch_id(face);
                if (!patch.has_value())
                  return owner_value;
                return impl_->boundaries->evaluate_pressure(*patch, owner_value)
                    .face;
              });
        });
    impl_->halo->exchange(*impl_->scratch.storage,
                          impl_->scratch.pressure_gradient);
    std::optional<runtime::FieldView<const double>> pressure_gradient_read;
    const std::size_t count = impl_->topology->owned_cell_count();
    const auto owned = impl_->topology->owned_global_box();
    auto &diagonal_values = impl_->diagonal_values;
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::momentum_linear_solve, true,
        "momentum equation preparation failed", [&] {
          pressure_gradient_read.emplace(
              impl_->scratch.storage->acquire_read<double>(
                  *impl_->scratch.access, kScratchPhase, kScratchActor,
                  impl_->scratch.pressure_gradient));
          for (std::size_t component_index = 0; component_index < 3U;
               ++component_index) {
            auto rhs = impl_->rhs[component_index].view(0U, count);
            auto predictor = impl_->predictor[component_index].view(0U, count);
            for (LocalCellId cell = 0; cell < count; ++cell) {
              const StructuredIndex index =
                  map_cell(impl_->topology->global_cell(cell), owned,
                           impl_->topology->global_extent());
              const double volume = impl_->geometry->cell_volume_m3(cell);
              const double diagonal =
                  stencil.alpha0 * rho_ref * volume / stencil.dt_s;
              double convection =
                  stencil.order == MomentumTimeOrder::backward_euler
                      ? residual_n.convection[cell * 3U + component_index]
                      : 2.0 * residual_n.convection[cell * 3U +
                                                    component_index] -
                            residual_nm1.convection[cell * 3U +
                                                    component_index];
              double viscosity =
                  stencil.order == MomentumTimeOrder::backward_euler
                      ? residual_n.viscosity[cell * 3U + component_index]
                      : 2.0 * residual_n.viscosity[cell * 3U +
                                                   component_index] -
                            residual_nm1.viscosity[cell * 3U +
                                                   component_index];
              double alpha1 = stencil.alpha1;
              double alpha2 = stencil.alpha2;
              double pressure = volume * (*pressure_gradient_read)(
                                               index.i, index.j, index.k,
                                               static_cast<int>(
                                                   component_index));
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
              const auto mutation =
                  momentum_assembly_mutation.load(std::memory_order_relaxed);
              if (mutation ==
                  kMomentumMutationOmitConvection) {
                convection = 0.0;
              } else if (mutation ==
                         kMomentumMutationOmitViscosity) {
                viscosity = 0.0;
              } else if (mutation ==
                         kMomentumMutationOmitPressure) {
                pressure = 0.0;
              } else if (mutation ==
                         kMomentumMutationOmitAlpha1) {
                alpha1 = 0.0;
              } else if (mutation ==
                         kMomentumMutationOmitAlpha2) {
                alpha2 = 0.0;
              }
#endif
              rhs[cell] =
                  -(rho_ref * volume / stencil.dt_s) *
                      (alpha1 *
                           (*velocity_n)(index.i, index.j, index.k,
                                         static_cast<int>(component_index)) +
                       alpha2 *
                           (*velocity_nm1)(index.i, index.j, index.k,
                                           static_cast<int>(component_index))) -
                  convection - viscosity - pressure;
              predictor[cell] = (*velocity_n)(
                  index.i, index.j, index.k, static_cast<int>(component_index));
              diagonal_values[component_index][cell] = diagonal;
              if (!std::isfinite(rhs[cell]) || !(diagonal > 0.0) ||
                  !std::isfinite(diagonal)) {
                throw runtime::Error("momentum equation is non-finite");
              }
            }
            impl_->operators[component_index]->replace(
                diagonal_values[component_index]);
            impl_->momentum_preconditioners[component_index]->update(
                *impl_->operators[component_index],
                impl_->operators[component_index]->revision());
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
            if (count != 0U) {
              last_momentum_rhs[component_index].store(
                  rhs[0], std::memory_order_relaxed);
              last_momentum_diagonal[component_index].store(
                  diagonal_values[component_index][0],
                  std::memory_order_relaxed);
            }
#endif
          }
        });
    std::optional<std::array<MomentumComponentEquation, 3>> equations;
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::momentum_linear_solve, true,
        "momentum solve view acquisition failed", [&] {
          equations.emplace(std::array<MomentumComponentEquation, 3>{{
              {impl_->operators[0].get(), impl_->momentum_preconditioners[0],
               static_cast<const execution::Buffer &>(impl_->rhs[0])
                   .view(0U, count),
               impl_->predictor[0].view(0U, count),
               impl_->diagonal[0].view(0U, count)},
              {impl_->operators[1].get(), impl_->momentum_preconditioners[1],
               static_cast<const execution::Buffer &>(impl_->rhs[1])
                   .view(0U, count),
               impl_->predictor[1].view(0U, count),
               impl_->diagonal[1].view(0U, count)},
              {impl_->operators[2].get(), impl_->momentum_preconditioners[2],
               static_cast<const execution::Buffer &>(impl_->rhs[2])
                   .view(0U, count),
               impl_->predictor[2].view(0U, count),
               impl_->diagonal[2].view(0U, count)},
          }});
        });
    report.momentum = impl_->momentum_predictor.solve(*impl_->mpi, *equations,
                                                      momentum_control);
    if (!report.momentum.all_converged()) {
      int failing_rank = -1;
      for (const auto &component_report : report.momentum.components) {
        if (!solve_success(component_report.reason) &&
            component_report.lowest_failing_rank >= 0 &&
            (failing_rank < 0 ||
             component_report.lowest_failing_rank < failing_rank)) {
          failing_rank = component_report.lowest_failing_rank;
        }
      }
      state.rollback_attempt();
      active = false;
      return numerical_failure(report, StepFailureReason::momentum_linear_solve,
                               failing_rank);
    }
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::non_finite_trial, true,
        "predictor staging failed", [&] {
          auto trial_velocity = trial.acquire_write<double>(
              access, kStatePhase, kStateActor, fields.velocity);
          auto actual_diagonal = impl_->scratch.storage->acquire_write<double>(
              *impl_->scratch.access, kScratchPhase, kScratchActor,
              impl_->scratch.actual_diagonal);
          for (LocalCellId cell = 0; cell < count; ++cell) {
            const StructuredIndex index =
                map_cell(impl_->topology->global_cell(cell), owned,
                         impl_->topology->global_extent());
            for (std::size_t component_index = 0; component_index < 3U;
                 ++component_index) {
              const double velocity =
                  impl_->predictor[component_index].view(0U, count)[cell];
              const double diagonal =
                  impl_->diagonal[component_index].view(0U, count)[cell];
              if (!std::isfinite(velocity) || !(diagonal > 0.0) ||
                  !std::isfinite(diagonal)) {
                throw runtime::Error("predictor staging is non-finite");
              }
              trial_velocity(index.i, index.j, index.k,
                             static_cast<int>(component_index)) = velocity;
              actual_diagonal(index.i, index.j, index.k,
                              static_cast<int>(component_index)) = diagonal;
            }
          }
        });
    impl_->halo->exchange(trial, fields.velocity);
    impl_->halo->exchange(*impl_->scratch.storage,
                          impl_->scratch.actual_diagonal);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    inject_attempt_stage_failure(*impl_->mpi,
                                 kAttemptFailureAfterMomentum);
#endif
    std::optional<runtime::FieldView<const double>> actual_diagonal_read;
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::non_finite_trial, true,
        "face predictor assembly failed", [&] {
          actual_diagonal_read.emplace(
              impl_->scratch.storage->acquire_read<double>(
                  *impl_->scratch.access, kScratchPhase, kScratchActor,
                  impl_->scratch.actual_diagonal));
          auto trial_face = trial.acquire_face_write<double>(
              access, kStatePhase, kStateActor, fields.face_velocity);
          const auto face_n = committed.acquire_face_read<double>(
              access, kStatePhase, kStateActor, fields.face_velocity);
          const auto face_nm1 = history.acquire_face_read<double>(
              access, kStatePhase, kStateActor, fields.face_velocity);
          MomentumFaceHistory face_history{
              *velocity_n, face_n,
              stencil.order == MomentumTimeOrder::bdf2 ? &*velocity_nm1
                                                       : nullptr,
              stencil.order == MomentumTimeOrder::bdf2 ? &face_nm1 : nullptr};
          const auto trial_velocity_read = trial.acquire_read<double>(
              access, kStatePhase, kStateActor, fields.velocity);
          impl_->face_assembler.assemble_constant_density(
              *impl_->boundaries, rho_ref, stencil, trial_velocity_read,
              *pressure_n, *pressure_gradient_read, *actual_diagonal_read,
              face_history, trial_face, registry, trial, access, kStatePhase,
              kStateActor, fields.face_mass_flux);
        });
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    inject_attempt_stage_failure(
        *impl_->mpi, kAttemptFailureAfterFacePredictor);
#endif
    const auto first = impl_->coupler.correct(
        state, rho_ref, *actual_diagonal_read, pressure_control);
    report.pressure[0] = first.solve;
    if (first.disposition != PressureCorrectionDisposition::accepted ||
        !first.accepted) {
      state.rollback_attempt();
      active = false;
      return first.disposition ==
                     PressureCorrectionDisposition::recoverable_failure
                 ? numerical_failure(report, first.reason,
                                     first.lowest_failing_rank)
                 : fatal_failure(report, first.reason,
                                 first.lowest_failing_rank);
    }
    report.pressure_corrector_count = 1U;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    inject_attempt_stage_failure(*impl_->mpi,
                                 kAttemptFailureAfterCorrectorOne);
#endif
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::transport_failure, true,
        "provisional flux provenance capture failed", [&] {
          const auto provisional_flux = trial.acquire_face_read<double>(
              access, kStatePhase, kStateActor, fields.face_mass_flux);
          impl_->provisional_face_mass_flux.resize(
              impl_->topology->local_face_count());
          for (LocalFaceId face = 0;
               face < impl_->topology->local_face_count(); ++face) {
            impl_->provisional_face_mass_flux[face] =
                provisional_flux(face, 0);
          }
        });
#endif
    recompute_transport(*impl_, state, rho_ref, stencil, sgs_viscosity,
                        active_wale_control, false);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::transport_failure, true,
        "provisional transport test staging failed", [&] {
          provisional_transport_call_count.fetch_add(1U,
                                                     std::memory_order_relaxed);
          if (provisional_transport_sentinel.load(std::memory_order_relaxed)) {
            for (const auto field : fields.transported_cell_fields) {
              auto values = trial.acquire_write<double>(access, kStatePhase,
                                                        kStateActor, field);
              const Int3 extent = values.interior_extent();
              for (int k = 0; k < extent.z; ++k) {
                for (int j = 0; j < extent.y; ++j) {
                  for (int i = 0; i < extent.x; ++i) {
                    values(i, j, k, 0) = 9.87654321e123;
                  }
                }
              }
            }
          }
        });
#endif
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    inject_attempt_stage_failure(
        *impl_->mpi, kAttemptFailureAfterProvisionalTransport);
#endif
    const auto second = impl_->coupler.correct(
        state, rho_ref, *actual_diagonal_read, pressure_control);
    report.pressure[1] = second.solve;
    if (second.disposition != PressureCorrectionDisposition::accepted ||
        !second.accepted) {
      state.rollback_attempt();
      active = false;
      return second.disposition ==
                     PressureCorrectionDisposition::recoverable_failure
                 ? numerical_failure(report, second.reason,
                                     second.lowest_failing_rank)
                 : fatal_failure(report, second.reason,
                                 second.lowest_failing_rank);
    }
    report.pressure_corrector_count = 2U;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    inject_attempt_stage_failure(*impl_->mpi,
                                 kAttemptFailureAfterCorrectorTwo);
#endif
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::transport_failure, true,
        "final flux test staging failed", [&] {
          const double injected_flux =
              final_uniform_x_mass_flux.load(std::memory_order_relaxed);
          if (injected_flux != 0.0) {
            auto final_flux = trial.acquire_face_write<double>(
                access, kStatePhase, kStateActor, fields.face_mass_flux);
            for (LocalFaceId face = 0;
                 face < impl_->topology->local_face_count(); ++face) {
              const double value =
                  injected_flux *
                  impl_->geometry->face_area_vector_m2(face, FaceSide::owner).x;
              if (final_uniform_x_mass_flux_override.load(
                      std::memory_order_relaxed)) {
                final_flux(face, 0) = value;
              } else {
                final_flux(face, 0) += value;
              }
            }
          }
        });
#endif
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    std::vector<double> authoritative_final_flux;
    if (transport_assembly_mutation.load(std::memory_order_relaxed) ==
        kTransportMutationUseProvisionalFlux) {
      synchronized_local_phase(
          *impl_->mpi, StepFailureReason::transport_failure, true,
          "final flux provenance mutation failed", [&] {
            const auto final_flux = trial.acquire_face_read<double>(
                access, kStatePhase, kStateActor, fields.face_mass_flux);
            authoritative_final_flux.resize(
                impl_->topology->local_face_count());
            for (LocalFaceId face = 0;
                 face < impl_->topology->local_face_count(); ++face) {
              authoritative_final_flux[face] = final_flux(face, 0);
            }
            if (impl_->provisional_face_mass_flux.size() !=
                authoritative_final_flux.size()) {
              throw runtime::Error(
                  "provisional flux provenance is incomplete");
            }
            auto mutation_flux = trial.acquire_face_write<double>(
                access, kStatePhase, kStateActor, fields.face_mass_flux);
            for (LocalFaceId face = 0;
                 face < impl_->topology->local_face_count(); ++face) {
              mutation_flux(face, 0) = impl_->provisional_face_mass_flux[face];
            }
          });
    }
#endif
    recompute_transport(*impl_, state, rho_ref, stencil, sgs_viscosity,
                        active_wale_control, true);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (!authoritative_final_flux.empty()) {
      synchronized_local_phase(
          *impl_->mpi, StepFailureReason::transport_failure, true,
          "authoritative final flux restoration failed", [&] {
            auto final_flux = trial.acquire_face_write<double>(
                access, kStatePhase, kStateActor, fields.face_mass_flux);
            for (LocalFaceId face = 0;
                 face < impl_->topology->local_face_count(); ++face) {
              final_flux(face, 0) = authoritative_final_flux[face];
            }
          });
    }
#endif
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    final_transport_call_count.fetch_add(1U, std::memory_order_relaxed);
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::non_finite_trial, true,
        "final equation mutation staging failed", [&] {
          const std::size_t momentum_component =
              final_momentum_perturb_component.load(std::memory_order_relaxed);
          const double momentum_delta =
              final_momentum_perturb_delta.load(std::memory_order_relaxed);
          if (momentum_component < 3U && momentum_delta != 0.0 && count != 0U) {
            auto values = trial.acquire_write<double>(
                access, kStatePhase, kStateActor, fields.velocity);
            const StructuredIndex index = map_cell(
                impl_->topology->global_cell(0U), owned,
                impl_->topology->global_extent());
            values(index.i, index.j, index.k,
                   static_cast<int>(momentum_component)) += momentum_delta;
          }
          const std::size_t transport_index =
              final_transport_perturb_index.load(std::memory_order_relaxed);
          const double transport_delta =
              final_transport_perturb_delta.load(std::memory_order_relaxed);
          if (transport_index < fields.transported_cell_fields.size() &&
              transport_delta != 0.0 && count != 0U) {
            auto values = trial.acquire_write<double>(
                access, kStatePhase, kStateActor,
                fields.transported_cell_fields[transport_index]);
            const StructuredIndex index = map_cell(
                impl_->topology->global_cell(0U), owned,
                impl_->topology->global_extent());
            values(index.i, index.j, index.k, 0) += transport_delta;
          }
        });
#endif
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    inject_attempt_stage_failure(
        *impl_->mpi, kAttemptFailureAfterFinalTransport);
#endif
    report.final_pressure_residual_l2 = second.independent_residual_l2;
    bool local_pressure = std::isfinite(second.independent_residual_l2) &&
                          second.independent_residual_l2 <=
                              std::max(pressure_control.atol,
                                       pressure_control.rtol * second.rhs_l2);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (force_final_pressure_failure.load(std::memory_order_relaxed)) {
      local_pressure = false;
    }
#endif
    const auto pressure_status = runtime::collective_status(
        *impl_->mpi, local_pressure,
        "final pressure residual exceeds tolerance");
    if (!pressure_status.ok) {
      state.rollback_attempt();
      active = false;
      return numerical_failure(report,
                               StepFailureReason::final_pressure_residual,
                               pressure_status.failing_rank);
    }
    const auto continuity = continuity_norms(*impl_, state);
    report.final_continuity_normalized_l2 = continuity.second;
    bool local_continuity = std::isfinite(continuity.second) &&
                            continuity.second <= kContinuityTolerance;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (force_final_continuity_failure.load(std::memory_order_relaxed)) {
      local_continuity = false;
    }
#endif
    const auto continuity_status = runtime::collective_status(
        *impl_->mpi, local_continuity,
        "final continuity residual exceeds tolerance");
    if (!continuity_status.ok) {
      state.rollback_attempt();
      active = false;
      return numerical_failure(report,
                               StepFailureReason::final_continuity_residual,
                               continuity_status.failing_rank);
    }
    const auto momentum_status = assess_final_momentum(
        *impl_, state, rho_ref, mu, effective_viscosity, stencil, report);
    if (!momentum_status.ok) {
      state.rollback_attempt();
      active = false;
      return numerical_failure(report,
                               StepFailureReason::final_momentum_residual,
                               momentum_status.failing_rank);
    }
    const auto transport_status = assess_final_transport(
        *impl_, state, rho_ref, sgs_viscosity, active_wale_control, stencil,
        report);
    if (!transport_status.ok) {
      state.rollback_attempt();
      active = false;
      return numerical_failure(report,
                               StepFailureReason::final_transport_residual,
                               transport_status.failing_rank);
    }
    const auto conservation_status =
        assess_final_conservation(*impl_, state, rho_ref, stencil, report);
    if (!conservation_status.ok) {
      state.rollback_attempt();
      active = false;
      return numerical_failure(report,
                               StepFailureReason::final_conservation_defect,
                               conservation_status.failing_rank);
    }
    std::optional<runtime::FaceFieldView<const double>> final_flux;
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::boundary_backflow, true,
        "final flux acquisition failed", [&] {
          final_flux.emplace(trial.acquire_face_read<double>(
              access, kStatePhase, kStateActor, fields.face_mass_flux));
        });
    const auto admissibility =
        impl_->boundaries->assess_final_pressure_outlet_flux(
            *impl_->topology, *impl_->mpi, *final_flux,
            state.metadata().step + 1U, state.metadata().time_s + stencil.dt_s);
    report.final_backflow_evidence = admissibility.evidence;
    if (admissibility.decision ==
        boundary::FinalFluxDecision::outlet_backflow) {
      const int rank = admissibility.evidence.has_value()
                           ? admissibility.evidence->lowest_failing_rank
                           : -1;
      state.rollback_attempt();
      active = false;
      return numerical_failure(report, StepFailureReason::boundary_backflow,
                               rank);
    }
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    inject_attempt_stage_failure(*impl_->mpi,
                                 kAttemptFailureBeforeCommit);
#endif
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::non_finite_trial, true,
        "collective commit preparation failed", [&] {
          auto density = trial.acquire_write<double>(
              access, kStatePhase, kStateActor, fields.density);
          for (LocalCellId cell = 0; cell < count; ++cell) {
            const StructuredIndex index =
                map_cell(impl_->topology->global_cell(cell), owned,
                         impl_->topology->global_extent());
            density(index.i, index.j, index.k, 0) = rho_ref;
          }
          const auto old = state.metadata();
          state.prepare_commit_attempt({old.step + 1U,
                                        old.time_s + stencil.dt_s, stencil.dt_s,
                                        old.dt_s, stencil.order});
        });
    state.publish_commit_attempt();
    active = false;
    report.disposition = StepAttemptDisposition::committed;
    report.reason = StepFailureReason::none;
    report.lowest_failing_rank = -1;
    report.suggested_dt_s = 0.0;
    if (wale_summary != nullptr)
      *wale_summary = wale_data->summary;
    return report;
  } catch (const SynchronizedAttemptFailure &failure) {
    if (active)
      state.rollback_attempt();
    return failure.recoverable
               ? numerical_failure(report, failure.reason, failure.failing_rank)
               : fatal_failure(report, failure.reason, failure.failing_rank);
  } catch (const runtime::detail::MpiOperationError &) {
    if (active)
      state.rollback_attempt();
    return fatal_failure(report, StepFailureReason::collective_operation,
                         -1);
  } catch (const runtime::Error &) {
    if (active)
      state.rollback_attempt();
    return fatal_failure(report, StepFailureReason::invalid_input, -1);
  } catch (...) {
    if (active)
      state.rollback_attempt();
    return fatal_failure(report, StepFailureReason::invalid_input, -1);
  }
}

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
namespace detail {
void constant_density_reset_raw() noexcept {
  ::hundun::flow::force_final_continuity_failure.store(
      false, std::memory_order_relaxed);
  ::hundun::flow::force_final_pressure_failure.store(false,
                                                     std::memory_order_relaxed);
  ::hundun::flow::force_local_derived_failure.store(false,
                                                    std::memory_order_relaxed);
  final_momentum_perturb_component.store(
      std::numeric_limits<std::size_t>::max(), std::memory_order_relaxed);
  final_momentum_perturb_delta.store(0.0, std::memory_order_relaxed);
  final_transport_perturb_index.store(
      std::numeric_limits<std::size_t>::max(), std::memory_order_relaxed);
  final_transport_perturb_delta.store(0.0, std::memory_order_relaxed);
  ::hundun::flow::force_final_conservation_failure.store(
      false, std::memory_order_relaxed);
  final_mass_defect_perturbation.store(0.0, std::memory_order_relaxed);
  for (auto &value : final_momentum_norm_armed)
    value.store(false, std::memory_order_relaxed);
  for (auto &value : final_momentum_norm_residual_square)
    value.store(0.0, std::memory_order_relaxed);
  for (auto &value : final_momentum_norm_scale_square)
    value.store(0.0, std::memory_order_relaxed);
  final_transport_norm_index.store(std::numeric_limits<std::size_t>::max(),
                                   std::memory_order_relaxed);
  final_transport_norm_residual_square.store(0.0, std::memory_order_relaxed);
  final_transport_norm_scale_square.store(0.0, std::memory_order_relaxed);
  for (auto &armed : momentum_conservation_parts_armed)
    armed.store(false, std::memory_order_relaxed);
  for (auto &component : momentum_conservation_parts_values) {
    for (auto &value : component)
      value.store(0.0, std::memory_order_relaxed);
  }
  momentum_conservation_overflow_component.store(
      std::numeric_limits<std::size_t>::max(), std::memory_order_relaxed);
  transport_conservation_overflow_index.store(
      std::numeric_limits<std::size_t>::max(), std::memory_order_relaxed);
  ::hundun::flow::last_mass_conservation = {};
  ::hundun::flow::last_momentum_conservation.fill({});
  ::hundun::flow::last_transport_conservation.clear();
  momentum_assembly_mutation.store(kMomentumMutationNone,
                                   std::memory_order_relaxed);
  transport_assembly_mutation.store(kTransportMutationNone,
                                    std::memory_order_relaxed);
  attempt_failure_stage.store(kAttemptFailureNone,
                              std::memory_order_relaxed);
  ::hundun::flow::last_pressure_constraint_mode.store(
      -1, std::memory_order_relaxed);
  ::hundun::flow::last_pressure_operator_mode.store(
      -1, std::memory_order_relaxed);
  provisional_transport_sentinel.store(false, std::memory_order_relaxed);
  final_uniform_x_mass_flux.store(0.0, std::memory_order_relaxed);
  final_uniform_x_mass_flux_override.store(false,
                                           std::memory_order_relaxed);
  for (auto &value : ::hundun::flow::last_momentum_rhs) {
    value.store(0.0, std::memory_order_relaxed);
  }
  for (auto &value : ::hundun::flow::last_momentum_diagonal) {
    value.store(0.0, std::memory_order_relaxed);
  }
  provisional_transport_call_count.store(0U, std::memory_order_relaxed);
  final_transport_call_count.store(0U, std::memory_order_relaxed);
  pressure_operator_construction_failure_rank.store(-1,
                                                    std::memory_order_relaxed);
  pressure_operator_refresh_failure_rank.store(-1,
                                               std::memory_order_relaxed);
}

void constant_density_force_final_continuity_failure_raw(
    bool enabled) noexcept {
  ::hundun::flow::force_final_continuity_failure.store(
      enabled, std::memory_order_relaxed);
}

void constant_density_force_final_pressure_failure_raw(
    bool enabled) noexcept {
  ::hundun::flow::force_final_pressure_failure.store(enabled,
                                                     std::memory_order_relaxed);
}

void constant_density_force_local_derived_failure_raw(
    bool enabled) noexcept {
  ::hundun::flow::force_local_derived_failure.store(enabled,
                                                    std::memory_order_relaxed);
}

void constant_density_force_final_momentum_perturbation_raw(
    std::size_t component, double delta) noexcept {
  final_momentum_perturb_component.store(component,
                                         std::memory_order_relaxed);
  final_momentum_perturb_delta.store(delta, std::memory_order_relaxed);
}

void constant_density_force_final_transport_perturbation_raw(
    std::size_t field_index, double delta) noexcept {
  final_transport_perturb_index.store(field_index,
                                      std::memory_order_relaxed);
  final_transport_perturb_delta.store(delta, std::memory_order_relaxed);
}

void constant_density_force_final_conservation_failure_raw(
    bool enabled) noexcept {
  ::hundun::flow::force_final_conservation_failure.store(
      enabled, std::memory_order_relaxed);
}

void constant_density_set_final_mass_defect_perturbation_raw(
    double delta) noexcept {
  final_mass_defect_perturbation.store(delta, std::memory_order_relaxed);
}

void constant_density_set_final_momentum_norm_squares_raw(
    std::size_t component, double residual_square,
    double scale_square) noexcept {
  if (component >= final_momentum_norm_residual_square.size())
    return;
  final_momentum_norm_armed[component].store(true, std::memory_order_relaxed);
  final_momentum_norm_residual_square[component].store(
      residual_square, std::memory_order_relaxed);
  final_momentum_norm_scale_square[component].store(
      scale_square, std::memory_order_relaxed);
}

void constant_density_set_final_transport_norm_squares_raw(
    std::size_t field_index, double residual_square,
    double scale_square) noexcept {
  final_transport_norm_index.store(field_index, std::memory_order_relaxed);
  final_transport_norm_residual_square.store(residual_square,
                                             std::memory_order_relaxed);
  final_transport_norm_scale_square.store(scale_square,
                                          std::memory_order_relaxed);
}

void constant_density_set_momentum_conservation_parts_raw(
    std::size_t component, const std::array<double, 8> &values) noexcept {
  if (component >= momentum_conservation_parts_values.size())
    return;
  for (std::size_t index = 0; index < values.size(); ++index) {
    momentum_conservation_parts_values[component][index].store(
        values[index], std::memory_order_relaxed);
  }
  momentum_conservation_parts_armed[component].store(true,
                                                     std::memory_order_relaxed);
}

void constant_density_force_momentum_conservation_aggregate_overflow_raw(
        std::size_t component, bool enabled) noexcept {
  momentum_conservation_overflow_component.store(
      enabled ? component : std::numeric_limits<std::size_t>::max(),
      std::memory_order_relaxed);
}

void constant_density_force_transport_conservation_aggregate_overflow_raw(
        std::size_t field_index, bool enabled) noexcept {
  transport_conservation_overflow_index.store(
      enabled ? field_index : std::numeric_limits<std::size_t>::max(),
      std::memory_order_relaxed);
}

void constant_density_set_momentum_assembly_mutation_raw(
    std::uint8_t mutation) noexcept {
  momentum_assembly_mutation.store(mutation, std::memory_order_relaxed);
}

void constant_density_set_transport_assembly_mutation_raw(
    std::uint8_t mutation) noexcept {
  transport_assembly_mutation.store(mutation, std::memory_order_relaxed);
}

void constant_density_set_attempt_failure_stage_raw(
    std::uint8_t stage) noexcept {
  attempt_failure_stage.store(stage, std::memory_order_relaxed);
}

void constant_density_set_pressure_operator_construction_failure_rank_raw(
    int rank) noexcept {
  pressure_operator_construction_failure_rank.store(rank,
                                                    std::memory_order_relaxed);
}

void constant_density_set_pressure_operator_refresh_failure_rank_raw(
    int rank) noexcept {
  pressure_operator_refresh_failure_rank.store(rank,
                                               std::memory_order_relaxed);
}

void constant_density_set_provisional_transport_sentinel_raw(
    bool enabled) noexcept {
  provisional_transport_sentinel.store(enabled, std::memory_order_relaxed);
}

void constant_density_set_final_uniform_x_mass_flux_raw(
    double value) noexcept {
  final_uniform_x_mass_flux.store(value, std::memory_order_relaxed);
}

void constant_density_set_final_uniform_x_mass_flux_override_raw(
    bool enabled) noexcept {
  final_uniform_x_mass_flux_override.store(enabled,
                                           std::memory_order_relaxed);
}

double constant_density_last_momentum_rhs_raw(
    std::size_t component) noexcept {
  return component < ::hundun::flow::last_momentum_rhs.size()
             ? ::hundun::flow::last_momentum_rhs[component].load(
                   std::memory_order_relaxed)
             : std::numeric_limits<double>::quiet_NaN();
}

double constant_density_last_momentum_diagonal_raw(
    std::size_t component) noexcept {
  return component < ::hundun::flow::last_momentum_diagonal.size()
             ? ::hundun::flow::last_momentum_diagonal[component].load(
                   std::memory_order_relaxed)
             : std::numeric_limits<double>::quiet_NaN();
}

std::size_t constant_density_provisional_transport_calls_raw() noexcept {
  return provisional_transport_call_count.load(std::memory_order_relaxed);
}

std::size_t constant_density_final_transport_calls_raw() noexcept {
  return final_transport_call_count.load(std::memory_order_relaxed);
}

int constant_density_last_pressure_constraint_mode_raw() noexcept {
  return ::hundun::flow::last_pressure_constraint_mode.load(
      std::memory_order_relaxed);
}

int constant_density_last_pressure_operator_mode_raw() noexcept {
  return ::hundun::flow::last_pressure_operator_mode.load(
      std::memory_order_relaxed);
}

void constant_density_last_mass_conservation_raw(
    std::array<double, 9> &values) noexcept {
  const auto &value = ::hundun::flow::last_mass_conservation;
  values = {value.quantity_nm1, value.quantity_n, value.quantity_np1,
            value.signed_boundary_flux, value.absolute_boundary_flux,
            value.boundary_integral, value.history_correction,
            value.raw_defect, value.relative_defect};
}

void constant_density_last_momentum_conservation_raw(
    std::size_t component, std::array<double, 9> &values) noexcept {
  if (component >= ::hundun::flow::last_momentum_conservation.size()) {
    values = {};
    return;
  }
  const auto &value = ::hundun::flow::last_momentum_conservation[component];
  values = {value.quantity_nm1, value.quantity_n, value.quantity_np1,
            value.signed_boundary_flux, value.absolute_boundary_flux,
            value.boundary_integral, value.history_correction,
            value.raw_defect, value.relative_defect};
}

void constant_density_last_transport_conservation_raw(
    std::size_t field_index, std::array<double, 9> &values) noexcept {
  if (field_index >= ::hundun::flow::last_transport_conservation.size()) {
    values = {};
    return;
  }
  const auto &value = ::hundun::flow::last_transport_conservation[field_index];
  values = {value.quantity_nm1, value.quantity_n, value.quantity_np1,
            value.signed_boundary_flux, value.absolute_boundary_flux,
            value.boundary_integral, value.history_correction,
            value.raw_defect, value.relative_defect};
}

void AdaptiveTimeControlAccess::constant_workspace_values_raw(
    const FixedStepConstantDensityFlow &flow, std::size_t &vector_count,
    std::size_t &total_capacity, std::uintptr_t &data_identity) noexcept {
  const auto add = [&](const auto &values) {
    ++vector_count;
    total_capacity += values.capacity();
    const auto pointer = reinterpret_cast<std::uintptr_t>(values.data());
    data_identity ^=
        pointer + UINT64_C(0x9e3779b97f4a7c15) +
        (data_identity << 6U) + (data_identity >> 2U);
  };
  const auto &impl = *flow.impl_;
  add(impl.pressure_gradient_sums);
  add(impl.continuity_absolute);
  for (const auto &values : impl.diagonal_values)
    add(values);
  add(impl.momentum_n.convection);
  add(impl.momentum_n.viscosity);
  add(impl.momentum_n.boundary);
  add(impl.momentum_nm1.convection);
  add(impl.momentum_nm1.viscosity);
  add(impl.momentum_nm1.boundary);
  add(impl.pressure_boundary);
  for (const auto *set : {&impl.transport_n, &impl.transport_nm1}) {
    for (const auto &values : *set) {
      add(values.convection);
      add(values.diffusion);
      add(values.boundary);
    }
  }
  const auto &coupler = *impl.coupler.impl_;
  add(coupler.pressure_gradient_sums);
  add(coupler.velocity_candidate);
  add(coupler.pressure_candidate);
  add(coupler.face_velocity_candidate);
  add(coupler.mass_flux_candidate);
}

void AdaptiveTimeControlAccess::constant_cache_values_raw(
    const FixedStepConstantDensityFlow &flow,
    std::vector<std::uintptr_t> &workspace_ids,
    std::vector<std::size_t> &workspace_capacities,
    std::array<std::uintptr_t, 3> &operator_ids,
    std::array<std::uint64_t, 3> &operator_revisions,
    std::array<std::vector<double>, 3> &operator_diagonals,
    std::size_t &operator_count, bool &delegated) {
  if (!flow.impl_)
    return;
  const auto add = [&](const auto &values) {
    workspace_ids.push_back(
        reinterpret_cast<std::uintptr_t>(values.data()));
    workspace_capacities.push_back(values.capacity());
  };
  const auto &impl = *flow.impl_;
  add(impl.pressure_gradient_sums);
  add(impl.continuity_absolute);
  for (const auto &values : impl.diagonal_values)
    add(values);
  add(impl.momentum_n.convection);
  add(impl.momentum_n.viscosity);
  add(impl.momentum_n.boundary);
  add(impl.momentum_nm1.convection);
  add(impl.momentum_nm1.viscosity);
  add(impl.momentum_nm1.boundary);
  add(impl.pressure_boundary);
  for (const auto *set : {&impl.transport_n, &impl.transport_nm1}) {
    for (const auto &values : *set) {
      add(values.convection);
      add(values.diffusion);
      add(values.boundary);
    }
  }
  const auto &coupler = *impl.coupler.impl_;
  add(coupler.pressure_gradient_sums);
  add(coupler.velocity_candidate);
  add(coupler.pressure_candidate);
  add(coupler.face_velocity_candidate);
  add(coupler.mass_flux_candidate);
  for (std::size_t component = 0; component < impl.operators.size();
       ++component) {
    if (!impl.operators[component])
      continue;
    const std::size_t entry = operator_count++;
    operator_ids[entry] =
        reinterpret_cast<std::uintptr_t>(impl.operators[component].get());
    operator_revisions[entry] = impl.operators[component]->revision();
    const auto count =
        impl.operators[component]->range_layout().owned_count();
    execution::Buffer diagonal(*impl.execution, bytes_for(count));
    impl.operators[component]->diagonal(diagonal.view(0U, count)).wait();
    const auto values =
        static_cast<const execution::Buffer &>(diagonal).view(0U, count);
    operator_diagonals[entry].reserve(count);
    for (std::size_t index = 0; index < count; ++index)
      operator_diagonals[entry].push_back(values[index]);
  }
  delegated = false;
}

void AdaptiveTimeControlAccess::constant_pressure_operator_values_raw(
    const FixedStepConstantDensityFlow &flow, std::uintptr_t &identity,
    std::uint64_t &revision, bool &present) noexcept {
  const auto &impl = *flow.impl_;
  const auto &candidate = impl.coupler.impl_->pressure_operator;
  present = candidate.has_value();
  if (present) {
    identity = reinterpret_cast<std::uintptr_t>(&*candidate);
    revision = candidate->revision();
  }
}

void AdaptiveTimeControlAccess::material_pressure_values_raw(
    const PisoCoupler &coupler, std::vector<double> &rhs_raw,
    std::vector<double> &rhs_solve, std::vector<double> &correction,
    std::vector<double> &final_face_density, double &rhs_l2,
    std::uint32_t &corrector_ordinal, bool &token_available,
    bool &final_operator_available,
    std::uint64_t &final_operator_revision) {
  if (!coupler.impl_)
    return;
  const auto &impl = *coupler.impl_;
  rhs_raw = impl.material_rhs_raw;
  rhs_solve = impl.material_rhs_solve;
  correction = impl.material_correction;
  final_face_density = impl.material_face_density;
  rhs_l2 = impl.material_rhs_l2;
  corrector_ordinal = impl.material_ordinal;
  token_available = impl.material_token_available;
  final_operator_available = impl.material_final_operator.has_value();
  if (final_operator_available)
    final_operator_revision = impl.material_final_operator->revision();
}

void AdaptiveTimeControlAccess::material_pressure_halo_values_raw(
    const PisoCoupler &coupler, bool &ordinary_available,
    runtime::HaloPerformanceCounters &ordinary,
    bool &material_final_available,
    runtime::HaloPerformanceCounters &material_final) {
  if (!coupler.impl_)
    return;
  const auto &impl = *coupler.impl_;
  ordinary_available = impl.pressure_operator.has_value();
  if (ordinary_available)
    ordinary = impl.pressure_operator->halo_performance_counters();
  material_final_available = impl.material_final_operator.has_value();
  if (material_final_available)
    material_final = impl.material_final_operator->halo_performance_counters();
}

runtime::HaloPerformanceCounters
material_combine_pressure_halo_counters_raw(
    runtime::HaloPerformanceCounters left,
    const runtime::HaloPerformanceCounters &right) {
  return ::hundun::flow::combine_halo_counters(left, right);
}
} // namespace detail
#endif

} // namespace hundun::flow

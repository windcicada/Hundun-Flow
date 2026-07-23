// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow/material_density_piso.hpp"

#include "density_closure_detail.hpp"
#include "fixed_step_flow_detail.hpp"
#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_operation_error.hpp"
#include <mpi.h>
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#include "material_density_piso_test_access.hpp"
#include "material_density_transport_test_access.hpp"
#endif

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstring>
#include <limits>
#include <new>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace hundun::flow {
namespace {

constexpr std::uint64_t kReportSeed = 0x6d6174657269616cULL;
constexpr runtime::PhaseId kStatePhase = 1800U;
constexpr runtime::ActorId kStateActor = 1800U;
constexpr runtime::PhaseId kScratchPhase = 1820U;
constexpr runtime::ActorId kScratchActor = 1820U;

void mix(std::uint64_t &hash, std::uint64_t value) noexcept {
  hash ^= value + 0x9e3779b97f4a7c15ULL + (hash << 6U) + (hash >> 2U);
}

std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

bool same(runtime::Int3 left, runtime::Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool same_stencil(const MomentumTimeStencil &left,
                  const MomentumTimeStencil &right) noexcept {
  return left.order == right.order && bits(left.dt_s) == bits(right.dt_s) &&
         bits(left.previous_dt_s) == bits(right.previous_dt_s) &&
         bits(left.alpha0) == bits(right.alpha0) &&
         bits(left.alpha1) == bits(right.alpha1) &&
         bits(left.alpha2) == bits(right.alpha2);
}

bool solve_success(linear::SolveTerminationReason reason) noexcept {
  return reason == linear::SolveTerminationReason::converged ||
         reason == linear::SolveTerminationReason::zero_right_hand_side;
}

class CompensatedSum final {
public:
  void add(double value) noexcept {
    const double adjusted = value - correction_;
    const double next = sum_ + adjusted;
    correction_ = (next - sum_) - adjusted;
    sum_ = next;
  }
  double value() const noexcept { return sum_; }

private:
  double sum_{};
  double correction_{};
};

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
std::atomic<int> preflight_allocation_failure_rank{-1};
using FaceFluxPathObservation = test::detail::FaceFluxPathObservation;
#endif

struct MomentumConservationTerms final {
  double momentum_n_minus_1{};
  double momentum_n{};
  double momentum_n_plus_1{};
  double boundary_n_minus_1{};
  double boundary_n{};
  double source_n_plus_1{};
  double momentum_abs_n_minus_1{};
  double momentum_abs_n{};
  double momentum_abs_n_plus_1{};
  double boundary_abs_n_minus_1{};
  double boundary_abs_n{};
  double source_abs_n_plus_1{};
  double dt_s{};
  double alpha0{};
  double alpha2{};
  bool bdf2{};
};

double momentum_conservation_defect(
    const MomentumConservationTerms &input) noexcept {
  const double history =
      (input.alpha2 / input.alpha0) *
      (input.momentum_n_minus_1 - input.momentum_n);
  const double boundary_effective =
      input.bdf2 ? 2.0 * input.boundary_n - input.boundary_n_minus_1
                 : input.boundary_n;
  const double boundary =
      (input.dt_s / input.alpha0) * boundary_effective;
  const double source =
      (input.dt_s / input.alpha0) * input.source_n_plus_1;
  const double raw = input.momentum_n_plus_1 - input.momentum_n + history +
                     boundary - source;
  const double boundary_abs_effective =
      input.bdf2 ? 2.0 * input.boundary_abs_n +
                       input.boundary_abs_n_minus_1
                 : input.boundary_abs_n;
  double denominator = std::max(
      {std::abs(input.momentum_n), std::abs(input.momentum_n_plus_1),
       std::abs(history),
       (input.dt_s / std::abs(input.alpha0)) * boundary_abs_effective,
       (input.dt_s / std::abs(input.alpha0)) *
           input.source_abs_n_plus_1,
       std::numeric_limits<double>::min()});
  const double cancellation = std::max(
      {input.momentum_abs_n, input.momentum_abs_n_plus_1,
       input.bdf2 ? input.momentum_abs_n_minus_1 : 0.0});
  if (denominator <=
      64.0 * std::numeric_limits<double>::epsilon() * cancellation)
    denominator = std::max(denominator, cancellation);
  return std::abs(raw) / denominator;
}

std::size_t bytes_for(std::size_t count) {
  if (count > std::numeric_limits<std::size_t>::max() / sizeof(double))
    throw runtime::Error("material flow vector byte count overflows");
  return count * sizeof(double);
}

std::size_t multiply_count(std::size_t count, std::size_t factor) {
  if (factor != 0U &&
      count > std::numeric_limits<std::size_t>::max() / factor)
    throw runtime::Error("material flow workspace count overflows");
  return count * factor;
}

struct StructuredIndex final {
  int i{};
  int j{};
  int k{};
};

StructuredIndex map_cell(runtime::Int3 global, runtime::Box3 owned,
                         runtime::Int3 global_extent) {
  const runtime::Int3 local{owned.end.x - owned.begin.x,
                            owned.end.y - owned.begin.y,
                            owned.end.z - owned.begin.z};
  const auto axis = [](int coordinate, int begin, int end, int global_n,
                       int local_n) {
    if (coordinate >= begin && coordinate < end)
      return coordinate - begin;
    if (coordinate == begin - 1 || (begin == 0 && coordinate == global_n - 1))
      return -1;
    if (coordinate == end || (end == global_n && coordinate == 0))
      return local_n;
    throw runtime::Error("material flow cell has no structured mapping");
  };
  return {axis(global.x, owned.begin.x, owned.end.x, global_extent.x, local.x),
          axis(global.y, owned.begin.y, owned.end.y, global_extent.y, local.y),
          axis(global.z, owned.begin.z, owned.end.z, global_extent.z, local.z)};
}

template <class T>
T &at(const runtime::FieldView<T> &view, StructuredIndex index,
      int component) {
  return view(index.i, index.j, index.k, component);
}

runtime::Real3 add(runtime::Real3 left, runtime::Real3 right) noexcept {
  return {left.x + right.x, left.y + right.y, left.z + right.z};
}

runtime::Real3 multiply(double scale, runtime::Real3 value) noexcept {
  return {scale * value.x, scale * value.y, scale * value.z};
}

double dot(runtime::Real3 left, runtime::Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
double component(runtime::Real3 value, int direction) noexcept {
  return direction == 0 ? value.x : direction == 1 ? value.y : value.z;
}
#endif

runtime::FieldDescriptor scratch_cell(std::string name,
                                      std::uint32_t components,
                                      int ghost_width) {
  return {std::move(name),
          "1",
          "fixed_step_material_density_flow",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          components,
          ghost_width,
          false,
          runtime::RestartPolicy::transient,
          runtime::OutputPolicy::never};
}

runtime::FieldDescriptor scratch_face(std::string name,
                                      std::uint32_t components) {
  return {std::move(name),
          "1",
          "fixed_step_material_density_flow",
          runtime::FunctionSpace::face_value,
          runtime::ScalarType::float64,
          components,
          0,
          false,
          runtime::RestartPolicy::transient,
          runtime::OutputPolicy::never};
}

class MaterialScratch final {
public:
  explicit MaterialScratch(runtime::FieldLayoutSet layout) {
    velocity_gradient = registry.declare_field(
        scratch_cell("material_velocity_gradient", 9U, 2));
    pressure_gradient = registry.declare_field(
        scratch_cell("material_pressure_gradient", 3U, 2));
    momentum_face =
        registry.declare_field(scratch_face("material_momentum_face", 3U));
    momentum_residual = registry.declare_field(
        scratch_cell("material_momentum_residual", 3U, 0));
    actual_diagonal = registry.declare_field(
        scratch_cell("material_actual_diagonal", 3U, 2));
    face_density_n =
        registry.declare_field(scratch_face("material_face_density_n", 1U));
    face_density_nm1 = registry.declare_field(
        scratch_face("material_face_density_nm1", 1U));
    face_density_trial = registry.declare_field(
        scratch_face("material_face_density_trial", 1U));
    mass_residual = registry.declare_field(
        scratch_cell("material_final_mass_residual", 1U, 0));
    registry.freeze();
    access = std::make_unique<runtime::FieldAccessPlan>(registry);
    for (runtime::FieldId field = 0U;
         field < static_cast<runtime::FieldId>(registry.size()); ++field)
      access->declare_access(kScratchPhase, kScratchActor, field,
                             runtime::AccessMode::read_write);
    access->freeze();
    storage = std::make_unique<runtime::FieldStorage>(registry, layout);
  }

  runtime::FieldRegistry registry;
  std::unique_ptr<runtime::FieldAccessPlan> access;
  std::unique_ptr<runtime::FieldStorage> storage;
  runtime::FieldId velocity_gradient{};
  runtime::FieldId pressure_gradient{};
  runtime::FieldId momentum_face{};
  runtime::FieldId momentum_residual{};
  runtime::FieldId actual_diagonal{};
  runtime::FieldId face_density_n{};
  runtime::FieldId face_density_nm1{};
  runtime::FieldId face_density_trial{};
  runtime::FieldId mass_residual{};
};

class MaterialDiagonalOperator final : public linear::LinearOperator {
public:
  MaterialDiagonalOperator(execution::ExecutionContext &context,
                           linear::VectorLayout layout)
      : context_(&context), layout_(std::move(layout)),
        diagonal_(layout_.owned_count(), 1.0) {}

  void replace(const std::vector<double> &values) {
    if (values.size() != diagonal_.size() ||
        !std::all_of(values.begin(), values.end(), [](double value) {
          return value > 0.0 && std::isfinite(value);
        }))
      throw runtime::Error("material momentum diagonal is invalid");
    if (revision_ == std::numeric_limits<std::uint64_t>::max())
      throw runtime::Error("material momentum revision would wrap");
    std::copy(values.begin(), values.end(), diagonal_.begin());
    ++revision_;
  }

  linear::VectorLayout domain_layout() const override { return layout_; }
  linear::VectorLayout range_layout() const override { return layout_; }
  const execution::ExecutionContext &context() const override {
    return *context_;
  }
  std::uint64_t revision() const override { return revision_; }
  execution::ExecutionEvent apply(execution::VectorView<const double> x,
                                  execution::VectorView<double> y) const override {
    if (x.size() != diagonal_.size() || y.size() != diagonal_.size())
      throw runtime::Error("material momentum operator view is incompatible");
    for (std::size_t cell = 0; cell < diagonal_.size(); ++cell)
      y[cell] = diagonal_[cell] * x[cell];
    return execution::ExecutionEvent::completed();
  }
  bool has_diagonal() const override { return true; }
  execution::ExecutionEvent
  diagonal(execution::VectorView<double> output) const override {
    if (output.size() != diagonal_.size())
      throw runtime::Error("material momentum diagonal view is incompatible");
    for (std::size_t cell = 0; cell < diagonal_.size(); ++cell)
      output[cell] = diagonal_[cell];
    return execution::ExecutionEvent::completed();
  }

private:
  execution::ExecutionContext *context_;
  linear::VectorLayout layout_;
  std::vector<double> diagonal_;
  std::uint64_t revision_{1U};
};

void zero(runtime::FieldView<double> view) {
  const auto extent = view.interior_extent();
  for (int k = 0; k < extent.z; ++k)
    for (int j = 0; j < extent.y; ++j)
      for (int i = 0; i < extent.x; ++i)
        for (std::uint32_t c = 0; c < view.components(); ++c)
          view(i, j, k, static_cast<int>(c)) = 0.0;
}

void validate_host_context(const execution::ExecutionContext &context) {
  if (context.backend_identity() == 0U ||
      context.space() != execution::ExecutionSpace::host ||
      !context.supports(execution::ExecutionCapability::host_access) ||
      !context.supports(execution::ExecutionCapability::buffer_allocation))
    throw runtime::Error("material flow requires a live host context");
}

void compute_pressure_gradient(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry,
    const runtime::FieldView<const double> &pressure,
    const runtime::FieldView<double> &gradient,
    std::vector<runtime::Real3> &sums) {
  std::fill(sums.begin(), sums.end(), runtime::Real3{});
  const auto owned = topology.owned_global_box();
  const auto global_extent = topology.global_extent();
  for (mesh::LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    const auto owner = topology.owner(face);
    const auto neighbour = topology.neighbour(face);
    const auto owner_index =
        map_cell(topology.global_cell(owner), owned, global_extent);
    const double owner_value = at(pressure, owner_index, 0);
    double face_value = owner_value;
    if (neighbour) {
      const auto neighbour_index =
          map_cell(topology.global_cell(*neighbour), owned, global_extent);
      const double neighbour_value = at(pressure, neighbour_index, 0);
      face_value = owner_value == neighbour_value
                       ? owner_value
                       : 0.5 * (owner_value + neighbour_value);
    }
    const auto area =
        geometry.face_area_vector_m2(face, mesh::FaceSide::owner);
    if (topology.cell_ownership(owner) == mesh::EntityOwnership::owned)
      sums[owner] = add(sums[owner], multiply(face_value, area));
    if (!topology.patch_id(face) && neighbour &&
        topology.cell_ownership(*neighbour) == mesh::EntityOwnership::owned)
      sums[*neighbour] = add(sums[*neighbour], multiply(-face_value, area));
  }
  for (mesh::LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
    const auto index =
        map_cell(topology.global_cell(cell), owned, global_extent);
    const auto value =
        multiply(1.0 / geometry.cell_volume_m3(cell), sums[cell]);
    gradient(index.i, index.j, index.k, 0) = value.x;
    gradient(index.i, index.j, index.k, 1) = value.y;
    gradient(index.i, index.j, index.k, 2) = value.z;
  }
}

struct MomentumResidual final {
  std::vector<double> convection;
  std::vector<double> viscosity;
  std::vector<finite_volume::PhysicalBoundaryMomentumContribution> boundary;
};

StepAttemptReport base_flow_report(double dt, std::size_t fields) {
  StepAttemptReport result;
  result.attempted_dt_s = dt;
  result.final_transport_normalized_l2.assign(fields, 0.0);
  result.final_transport_relative_conservation_defect.assign(fields, 0.0);
  return result;
}

void set_failure(StepAttemptReport &report, StepFailureReason reason, int rank,
                 bool recoverable) {
  report.disposition = recoverable
                           ? StepAttemptDisposition::recoverable_failure
                           : StepAttemptDisposition::non_retryable_failure;
  report.reason = reason;
  report.lowest_failing_rank = rank;
  report.suggested_dt_s = recoverable ? 0.5 * report.attempted_dt_s
                                      : report.attempted_dt_s;
}

StepFailureReason map_material_failure(MaterialTransportFailureReason reason) {
  switch (reason) {
  case MaterialTransportFailureReason::invalid_input:
    return StepFailureReason::invalid_input;
  case MaterialTransportFailureReason::non_finite_state:
  case MaterialTransportFailureReason::non_positive_density:
    return StepFailureReason::transport_failure;
  case MaterialTransportFailureReason::final_density_residual:
    return StepFailureReason::final_continuity_residual;
  case MaterialTransportFailureReason::final_transport_residual:
    return StepFailureReason::final_transport_residual;
  case MaterialTransportFailureReason::final_conservation_defect:
    return StepFailureReason::final_conservation_defect;
  case MaterialTransportFailureReason::collective_operation:
    return StepFailureReason::collective_operation;
  case MaterialTransportFailureReason::none:
    return StepFailureReason::none;
  }
  return StepFailureReason::invalid_input;
}

struct SynchronizedMaterialFailure final {
  StepFailureReason reason{StepFailureReason::invalid_input};
  int failing_rank{-1};
  bool recoverable{};
};

struct LocalMaterialFailure final {
  bool failed{};
  StepFailureReason reason{StepFailureReason::invalid_input};
  bool recoverable{};
};

SynchronizedMaterialFailure synchronize_local_failure(
    const runtime::MpiContext &mpi, LocalMaterialFailure local) {
  const int candidate = local.failed ? mpi.rank() : mpi.size();
  int lowest = mpi.size();
  runtime::check_mpi_result(
      MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, mpi.comm()),
      "MPI_Allreduce(material flow phase failure rank)");
  if (lowest == mpi.size())
    return {StepFailureReason::none, -1, false};
  int payload[2]{};
  if (mpi.rank() == lowest) {
    payload[0] = static_cast<int>(local.reason);
    payload[1] = local.recoverable ? 1 : 0;
  }
  runtime::check_mpi_result(
      MPI_Bcast(payload, 2, MPI_INT, lowest, mpi.comm()),
      "MPI_Bcast(material flow phase failure)");
  if (payload[0] < static_cast<int>(StepFailureReason::none) ||
      payload[0] >
          static_cast<int>(StepFailureReason::density_closure_failure) ||
      (payload[1] != 0 && payload[1] != 1))
    throw runtime::Error("material flow phase failure payload is invalid");
  return {static_cast<StepFailureReason>(payload[0]), lowest,
          payload[1] != 0};
}

void require_reliable_collective_result(StepFailureReason reason,
                                        int lowest_failing_rank) {
  if (reason == StepFailureReason::collective_operation &&
      lowest_failing_rank < 0)
    throw runtime::MpiOperationError(
        "material-density flow collective failure has no reliable rank");
}

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
void apply_terminal_test_point(test::MaterialTerminalPointForTest point) {
  const auto mode = test::detail::reach_material_terminal_point(point);
  if (mode == test::MaterialTerminalModeForTest::returned_rankless)
    throw runtime::MpiOperationError(
        "material-density flow collective failure has no reliable rank");
  if (mode == test::MaterialTerminalModeForTest::thrown_operation)
    throw runtime::MpiOperationError(
        "injected material terminal operation failure");
}

#endif

template <class Function>
void synchronized_local_phase(const runtime::MpiContext &mpi,
                              StepFailureReason reason, bool recoverable,
                              Function &&function) {
  LocalMaterialFailure local{};
  try {
    function();
  } catch (const runtime::MpiOperationError &) {
    throw;
  } catch (...) {
    local = {true, reason, recoverable};
  }
  const auto selected = synchronize_local_failure(mpi, local);
  if (selected.failing_rank >= 0)
    throw selected;
}

} // namespace

struct FixedStepMaterialDensityFlow::Impl final {
  Impl(const runtime::StructuredDecomposition &decomposition,
       const mesh::MeshTopology &topology,
       const mesh::MeshGeometry &geometry,
       const boundary::BoundaryRegistry &boundaries,
       const runtime::MpiContext &mpi,
       execution::ExecutionContext &execution,
       runtime::HaloExchange &halo,
       const linear::LinearSolver &momentum_solver,
       std::array<linear::Preconditioner *, 3> preconditioners,
       PisoCoupler pressure_coupler, MaterialDensityTransport transport,
       const runtime::FieldRegistry &registry, FlowFieldIds fields,
       MaterialDensityTransportSpec specification)
      : decomposition(&decomposition), topology(&topology), geometry(&geometry),
        boundaries(&boundaries), mpi(&mpi), execution(&execution), halo(&halo),
        momentum_solver(&momentum_solver),
        momentum_preconditioners(preconditioners),
        coupler(std::move(pressure_coupler)),
        material_transport(std::move(transport)), registry(&registry),
        fields(std::move(fields)), specification(std::move(specification)),
        fvm(finite_volume::CellCenteredFvmOperators::create(topology,
                                                            geometry)),
        prepared_face_flux(finite_volume::FaceMassFlux::prepare(topology)),
        prepared_material_flux(MaterialFaceMassFlux::prepare(topology)),
        face_assembler(TimeConsistentFaceVelocity::create(topology, geometry)),
        scratch({decomposition.local_extent(), topology.local_face_count()}),
        rhs{execution::Buffer(execution, bytes_for(topology.owned_cell_count())),
            execution::Buffer(execution, bytes_for(topology.owned_cell_count())),
            execution::Buffer(execution, bytes_for(topology.owned_cell_count()))},
        predictor{
            execution::Buffer(execution, bytes_for(topology.owned_cell_count())),
            execution::Buffer(execution, bytes_for(topology.owned_cell_count())),
            execution::Buffer(execution, bytes_for(topology.owned_cell_count()))},
        diagonal{
            execution::Buffer(execution, bytes_for(topology.owned_cell_count())),
            execution::Buffer(execution, bytes_for(topology.owned_cell_count())),
            execution::Buffer(execution, bytes_for(topology.owned_cell_count()))},
        pressure_gradient_sums(topology.owned_cell_count()),
        partition_face_values(
            multiply_count(topology.global_face_count(), 4U)),
        partition_face_counts(topology.global_face_count()) {
    const std::size_t cells = topology.owned_cell_count();
    const std::size_t momentum_values = multiply_count(cells, 3U);
    for (auto &values : diagonal_values)
      values.resize(cells);
    momentum_n.convection.resize(momentum_values);
    momentum_n.viscosity.resize(momentum_values);
    momentum_n.boundary.reserve(topology.local_face_count());
    momentum_nm1.convection.resize(momentum_values);
    momentum_nm1.viscosity.resize(momentum_values);
    momentum_nm1.boundary.reserve(topology.local_face_count());
    pressure_boundary.reserve(topology.local_face_count());
    continuity_absolute.resize(cells);
    const auto layout = linear::VectorLayout::from_topology(topology);
    for (std::size_t direction = 0; direction < 3U; ++direction)
      operators[direction] =
          std::make_unique<MaterialDiagonalOperator>(execution, layout);
    material_field_count =
        static_cast<std::uint64_t>(1U + this->specification.scalar_densities.size());
    for (mesh::LocalFaceId face = 0; face < topology.local_face_count(); ++face)
      if (topology.cell_ownership(topology.owner(face)) ==
          mesh::EntityOwnership::owned)
        canonical_faces.push_back(face);
    const auto box = topology.owned_global_box();
    const auto global = topology.global_extent();
    owned_cell_fingerprint =
        "cell.f64.owned." + std::to_string(box.begin.x) + "." +
        std::to_string(box.begin.y) + "." + std::to_string(box.begin.z) +
        "." + std::to_string(box.end.x) + "." +
        std::to_string(box.end.y) + "." + std::to_string(box.end.z);
    global_cell_fingerprint =
        "cell.f64.global." + std::to_string(global.x) + "." +
        std::to_string(global.y) + "." + std::to_string(global.z);
    owned_face_fingerprint =
        "face.f64.owned." + std::to_string(canonical_faces.size());
    global_face_fingerprint =
        "face.f64.global." + std::to_string(topology.global_face_count());
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    finalizer_flux_input.resize(topology.local_face_count());
#endif
  }

  const runtime::StructuredDecomposition *decomposition;
  const mesh::MeshTopology *topology;
  const mesh::MeshGeometry *geometry;
  const boundary::BoundaryRegistry *boundaries;
  const runtime::MpiContext *mpi;
  execution::ExecutionContext *execution;
  runtime::HaloExchange *halo;
  const linear::LinearSolver *momentum_solver;
  std::array<linear::Preconditioner *, 3> momentum_preconditioners;
  PisoCoupler coupler;
  MaterialDensityTransport material_transport;
  const runtime::FieldRegistry *registry;
  FlowFieldIds fields;
  MaterialDensityTransportSpec specification;
  finite_volume::CellCenteredFvmOperators fvm;
  finite_volume::FaceMassFlux::PreparedStatePtr prepared_face_flux;
  MaterialFaceMassFlux::PreparedStatePtr prepared_material_flux;
  TimeConsistentFaceVelocity face_assembler;
  MaterialScratch scratch;
  std::array<execution::Buffer, 3> rhs;
  std::array<execution::Buffer, 3> predictor;
  std::array<execution::Buffer, 3> diagonal;
  std::array<std::unique_ptr<MaterialDiagonalOperator>, 3> operators;
  std::array<std::vector<double>, 3> diagonal_values;
  std::vector<runtime::Real3> pressure_gradient_sums;
  std::vector<double> continuity_absolute;
  std::vector<double> partition_face_values;
  std::vector<double> partition_face_counts;
  MomentumResidual momentum_n;
  MomentumResidual momentum_nm1;
  std::vector<finite_volume::PhysicalBoundaryPressureContribution>
      pressure_boundary;
  std::uint64_t material_field_count{};
  std::uint64_t attempt_identity{};
  std::uint64_t source_generation{1U};
  const FlowState *last_state{};
  std::uint64_t last_report_seal{};
  std::uint64_t last_state_identity{};
  std::vector<mesh::LocalFaceId> canonical_faces;
  std::string owned_cell_fingerprint;
  std::string global_cell_fingerprint;
  std::string owned_face_fingerprint;
  std::string global_face_fingerprint;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
  bool vortex_source_enabled{};
  std::vector<double> finalizer_flux_input;
#endif
};

struct MaterialDensityFlowDiagnosticSource::Impl final {
  const FixedStepMaterialDensityFlow::Impl *flow{};
  const FlowState *state{};
  MaterialDensityStepAttemptReport report;
  std::uint64_t flow_identity{};
  std::uint64_t source_generation{};
  std::uint64_t state_identity{};
  std::uint64_t report_seal{};
  int rank{};
  std::uint64_t step{};
  double time_s{};
};

namespace {

template <class Implementation, class FluxBinder>
void assemble_spatial(Implementation &impl,
                      const runtime::FieldAccessPlan &access,
                      runtime::FieldStorage &storage, const FlowFieldIds &fields,
                      double mu, MomentumResidual &output,
                      FluxBinder &&bind_flux) {
  synchronized_local_phase(*impl.mpi, StepFailureReason::invalid_input, false,
                           [&] {
    const auto velocity = storage.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.velocity);
    auto gradient = impl.scratch.storage->template acquire_write<double>(
        *impl.scratch.access, kScratchPhase, kScratchActor,
        impl.scratch.velocity_gradient);
    impl.fvm.compute_gradient(finite_volume::GradientScheme::green_gauss,
                              finite_volume::FiniteVolumeQuantity::velocity(),
                              *impl.boundaries, velocity, gradient);
  });
  impl.halo->exchange(*impl.scratch.storage, impl.scratch.velocity_gradient);
  synchronized_local_phase(*impl.mpi, StepFailureReason::invalid_input, false,
                           [&] {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    FaceFluxPathObservation allocation_observation;
#endif
    const auto velocity = storage.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.velocity);
    const auto gradient_read =
        impl.scratch.storage->template acquire_read<double>(
            *impl.scratch.access, kScratchPhase, kScratchActor,
            impl.scratch.velocity_gradient);
    const auto flux = bind_flux(storage);
    auto face = impl.scratch.storage->template acquire_face_write<double>(
        *impl.scratch.access, kScratchPhase, kScratchActor,
        impl.scratch.momentum_face);
    impl.fvm.reconstruct_momentum_faces(*impl.boundaries, flux, velocity, face);
    const auto face_read =
        impl.scratch.storage->template acquire_face_read<double>(
            *impl.scratch.access, kScratchPhase, kScratchActor,
            impl.scratch.momentum_face);
    auto residual = impl.scratch.storage->template acquire_write<double>(
        *impl.scratch.access, kScratchPhase, kScratchActor,
        impl.scratch.momentum_residual);
    zero(residual);
    impl.fvm.accumulate_convective_residual(flux, face_read, residual);
    const auto owned = impl.topology->owned_global_box();
    for (mesh::LocalCellId cell = 0; cell < impl.topology->owned_cell_count();
         ++cell) {
      const auto index = map_cell(impl.topology->global_cell(cell), owned,
                                  impl.topology->global_extent());
      for (int direction = 0; direction < 3; ++direction)
        output.convection[cell * 3U +
                          static_cast<std::size_t>(direction)] =
            residual(index.i, index.j, index.k, direction);
    }
    zero(residual);
    impl.fvm.accumulate_viscous_residual(*impl.boundaries, velocity,
                                         gradient_read, mu, residual);
    for (mesh::LocalCellId cell = 0; cell < impl.topology->owned_cell_count();
         ++cell) {
      const auto index = map_cell(impl.topology->global_cell(cell), owned,
                                  impl.topology->global_extent());
      for (int direction = 0; direction < 3; ++direction)
        output.viscosity[cell * 3U +
                         static_cast<std::size_t>(direction)] =
            residual(index.i, index.j, index.k, direction);
    }
    impl.fvm.physical_boundary_momentum_contributions(
        *impl.boundaries, flux, face_read, velocity, gradient_read, mu,
        output.boundary);
  });
}

template <class Implementation, class FluxBinder>
void reconstruct_density_faces(Implementation &impl,
                               const runtime::FieldAccessPlan &access,
                               runtime::FieldStorage &storage,
                               const FlowFieldIds &fields,
                               runtime::FieldId output_field,
                               FluxBinder &&bind_flux) {
  synchronized_local_phase(*impl.mpi, StepFailureReason::invalid_input, false,
                           [&] {
    const auto flux = bind_flux(storage);
    const auto density = storage.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.density);
    auto face_density =
        impl.scratch.storage->template acquire_face_write<double>(
            *impl.scratch.access, kScratchPhase, kScratchActor, output_field);
    impl.fvm.reconstruct_transport_faces(
        finite_volume::FiniteVolumeQuantity::density(), *impl.boundaries, flux,
        density, face_density);
    for (mesh::LocalFaceId face = 0; face < impl.topology->local_face_count();
         ++face) {
      if (!(face_density(face, 0) > 0.0) ||
          !std::isfinite(face_density(face, 0)))
        throw runtime::Error("material reconstructed face density is invalid");
    }
  });
}

template <class Implementation, class FluxBinder>
void build_density_weighted_flux(Implementation &impl, FlowState &state,
                                 const runtime::FieldAccessPlan &access,
                                 const FlowFieldIds &fields,
                                 FluxBinder &&bind_flux) {
  auto &trial = detail::FlowStateSolverAccess::layer(state, FlowLayer::trial);
  synchronized_local_phase(*impl.mpi, StepFailureReason::non_finite_trial, true,
                           [&] {
    const auto velocity = trial.acquire_face_read<double>(
        access, kStatePhase, kStateActor, fields.face_velocity);
    auto direction = trial.acquire_face_write<double>(
        access, kStatePhase, kStateActor, fields.face_mass_flux);
    for (mesh::LocalFaceId face = 0;
         face < impl.topology->local_face_count(); ++face) {
      const runtime::Real3 value{velocity(face, 0), velocity(face, 1),
                                 velocity(face, 2)};
      direction(face, 0) = dot(
          value,
          impl.geometry->face_area_vector_m2(face, mesh::FaceSide::owner));
      if (!std::isfinite(direction(face, 0)))
        throw runtime::Error("material face-flow direction is non-finite");
    }
  });
  reconstruct_density_faces(impl, access, trial, fields,
                            impl.scratch.face_density_trial,
                            std::forward<FluxBinder>(bind_flux));
  synchronized_local_phase(*impl.mpi, StepFailureReason::non_finite_trial, true,
                           [&] {
    const auto density =
        impl.scratch.storage->template acquire_face_read<double>(
            *impl.scratch.access, kScratchPhase, kScratchActor,
            impl.scratch.face_density_trial);
    auto flux = trial.acquire_face_write<double>(
        access, kStatePhase, kStateActor, fields.face_mass_flux);
    for (mesh::LocalFaceId face = 0; face < impl.topology->local_face_count();
         ++face) {
      flux(face, 0) *= density(face, 0);
      if (!std::isfinite(flux(face, 0)))
        throw runtime::Error("material face mass flux is non-finite");
    }
  });
}

template <class Implementation>
void synchronize_partition_faces(Implementation &impl, FlowState &state,
                                 const runtime::FieldAccessPlan &access,
                                 const FlowFieldIds &fields) {
  auto &trial = detail::FlowStateSolverAccess::layer(state, FlowLayer::trial);
  synchronized_local_phase(*impl.mpi, StepFailureReason::invalid_input, false,
                           [&] {
    std::fill(impl.partition_face_values.begin(),
              impl.partition_face_values.end(), 0.0);
    std::fill(impl.partition_face_counts.begin(),
              impl.partition_face_counts.end(), 0.0);
    const auto flux = trial.acquire_face_read<double>(
        access, kStatePhase, kStateActor, fields.face_mass_flux);
    const auto velocity = trial.acquire_face_read<double>(
        access, kStatePhase, kStateActor, fields.face_velocity);
    for (mesh::LocalFaceId face = 0; face < impl.topology->local_face_count();
         ++face) {
      if (impl.topology->cell_ownership(impl.topology->owner(face)) !=
          mesh::EntityOwnership::owned)
        continue;
      const std::size_t global =
          static_cast<std::size_t>(impl.topology->global_face_id(face));
      if (global >= impl.partition_face_counts.size())
        throw runtime::Error("material partition face identity is invalid");
      impl.partition_face_values[global * 4U] = flux(face, 0);
      for (std::size_t component = 0; component < 3U; ++component)
        impl.partition_face_values[global * 4U + 1U + component] =
            velocity(face, static_cast<int>(component));
      impl.partition_face_counts[global] = 1.0;
    }
  });
  impl.mpi->allreduce_fp64_in_place(
      impl.partition_face_values.data(), impl.partition_face_values.size(),
      runtime::Fp64ReductionOperation::sum);
  impl.mpi->allreduce_fp64_in_place(
      impl.partition_face_counts.data(), impl.partition_face_counts.size(),
      runtime::Fp64ReductionOperation::sum);
  synchronized_local_phase(*impl.mpi, StepFailureReason::non_finite_trial, true,
                           [&] {
    auto flux_write = trial.acquire_face_write<double>(
        access, kStatePhase, kStateActor, fields.face_mass_flux);
    auto velocity_write = trial.acquire_face_write<double>(
        access, kStatePhase, kStateActor, fields.face_velocity);
    for (mesh::LocalFaceId face = 0; face < impl.topology->local_face_count();
         ++face) {
      const std::size_t global =
          static_cast<std::size_t>(impl.topology->global_face_id(face));
      const auto periodic = impl.topology->periodic_pair(face);
      const std::size_t representative =
          periodic.has_value()
              ? static_cast<std::size_t>(std::min(
                    impl.topology->global_face_id(face), *periodic))
              : global;
      if (global >= impl.partition_face_counts.size() ||
          representative >= impl.partition_face_counts.size() ||
          impl.partition_face_counts[representative] != 1.0)
        throw runtime::Error("material partition face owner is unavailable");
      double value = impl.partition_face_values[representative * 4U];
      if (periodic.has_value() &&
          impl.topology->global_face_id(face) > *periodic)
        value = -value;
      if (!std::isfinite(value))
        throw runtime::Error("material partition face flux is non-finite");
      flux_write(face, 0) = value;
      for (std::size_t component = 0; component < 3U; ++component) {
        const double item =
            impl.partition_face_values[representative * 4U + 1U + component];
        if (!std::isfinite(item))
          throw runtime::Error(
              "material partition face velocity is non-finite");
        velocity_write(face, static_cast<int>(component)) = item;
      }
    }
  });
}

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
runtime::Real3 analytic_source(double x, double y, double mu) noexcept {
  const auto value =
      test::MaterialDensityPisoTestAccess::vortex_source(x, y, mu);
  return {value.x, value.y, value.z};
}
#endif

} // namespace

const StepAttemptReport &MaterialDensityStepAttemptReport::flow() const
    noexcept {
  return flow_;
}

bool MaterialDensityStepAttemptReport::material_report_available() const
    noexcept {
  return material_report_.has_value();
}

const MaterialDensityTransportReport &
MaterialDensityStepAttemptReport::material_report() const {
  if (!material_report_) {
    throw runtime::Error("material-density report is unavailable");
  }
  return *material_report_;
}

MaterialTransportFailureReason
MaterialDensityStepAttemptReport::material_failure_reason() const noexcept {
  return material_failure_reason_;
}

std::uint64_t MaterialDensityStepAttemptReport::material_field_count() const
    noexcept {
  return material_field_count_;
}

runtime::FieldId
MaterialDensityStepAttemptReport::shared_face_mass_flux_field() const noexcept {
  return shared_face_mass_flux_field_;
}

MaterialFluxProvenance
MaterialDensityStepAttemptReport::flux_provenance() const noexcept {
  return flux_provenance_;
}

std::uint64_t
MaterialDensityStepAttemptReport::attempt_identity() const noexcept {
  return attempt_identity_;
}

bool MaterialDensityStepAttemptReport::final_continuity_residual_available()
    const noexcept {
  return final_continuity_residual_available_;
}

bool MaterialDensityStepAttemptReport::final_pressure_residual_available() const
    noexcept {
  return final_pressure_residual_available_;
}

double MaterialDensityStepAttemptReport::final_pressure_normalized_residual()
    const noexcept {
  return final_pressure_normalized_residual_;
}

const std::array<std::uint8_t, 3> &
MaterialDensityStepAttemptReport::final_momentum_residual_availability() const
    noexcept {
  return final_momentum_residual_available_;
}

bool MaterialDensityStepAttemptReport::mass_conservation_available() const
    noexcept {
  return mass_conservation_available_;
}

const std::array<std::uint8_t, 3> &
MaterialDensityStepAttemptReport::momentum_conservation_availability() const
    noexcept {
  return momentum_conservation_available_;
}

std::uint64_t MaterialDensityStepAttemptReport::compute_seal() const noexcept {
  std::uint64_t result = kReportSeed;
  const auto mix_solve = [&](const linear::SolveReport &solve) {
    mix(result, static_cast<std::uint64_t>(solve.reason));
    mix(result, solve.iterations);
    mix(result, bits(solve.initial_residual));
    mix(result, bits(solve.recursive_residual));
    mix(result, bits(solve.final_residual));
    mix(result, solve.matvec_count);
    mix(result, solve.preconditioner_apply_count);
    mix(result, solve.global_reduction_count);
    mix(result, static_cast<std::uint64_t>(solve.lowest_failing_rank + 1));
  };
  mix(result, static_cast<std::uint64_t>(flow_.disposition));
  mix(result, static_cast<std::uint64_t>(flow_.reason));
  mix(result, static_cast<std::uint64_t>(flow_.lowest_failing_rank + 1));
  mix(result, flow_.pressure_corrector_count);
  mix(result, bits(flow_.attempted_dt_s));
  mix(result, bits(flow_.suggested_dt_s));
  for (const auto &solve : flow_.momentum.components)
    mix_solve(solve);
  for (const auto &solve : flow_.pressure)
    mix_solve(solve);
  mix(result, bits(flow_.final_continuity_normalized_l2));
  mix(result, bits(flow_.final_pressure_residual_l2));
  for (const double value : flow_.final_momentum_normalized_l2)
    mix(result, bits(value));
  mix(result, flow_.final_transport_normalized_l2.size());
  for (const double value : flow_.final_transport_normalized_l2)
    mix(result, bits(value));
  mix(result, bits(flow_.final_mass_relative_conservation_defect));
  for (const double value : flow_.final_momentum_relative_conservation_defect)
    mix(result, bits(value));
  mix(result, flow_.final_transport_relative_conservation_defect.size());
  for (const double value : flow_.final_transport_relative_conservation_defect)
    mix(result, bits(value));
  mix(result, flow_.final_backflow_evidence.has_value() ? 1U : 0U);
  if (flow_.final_backflow_evidence) {
    const auto &evidence = *flow_.final_backflow_evidence;
    mix(result, evidence.patch_id);
    mix(result, evidence.step);
    mix(result, bits(evidence.time_s));
    mix(result, bits(evidence.minimum_outward_mass_flux_kg_per_s));
    mix(result, evidence.global_face_id);
    mix(result,
        static_cast<std::uint64_t>(evidence.lowest_failing_rank + 1));
  }
  mix(result, material_report_.has_value() ? 1U : 0U);
  if (material_report_) {
    const auto &material = *material_report_;
    mix(result, static_cast<std::uint64_t>(material.disposition()));
    mix(result, static_cast<std::uint64_t>(material.reason()));
    mix(result,
        static_cast<std::uint64_t>(material.lowest_failing_rank() + 1));
    const auto &stencil = material.stencil();
    mix(result, static_cast<std::uint64_t>(stencil.order));
    mix(result, bits(stencil.dt_s));
    mix(result, bits(stencil.previous_dt_s));
    mix(result, bits(stencil.alpha0));
    mix(result, bits(stencil.alpha1));
    mix(result, bits(stencil.alpha2));
    mix(result, static_cast<std::uint64_t>(material.flux_provenance()));
    mix(result, material.attempt_identity());
    mix(result, material.finalization_identity());
    mix(result,
        static_cast<std::uint64_t>(material.shared_face_mass_flux_field()));
    mix(result, material.density_residual_available() ? 1U : 0U);
    mix(result, bits(material.density_normalized_l2()));
    mix(result, material.transport_residual_availability().size());
    for (const auto value : material.transport_residual_availability())
      mix(result, value);
    mix(result, material.transport_normalized_l2().size());
    for (const double value : material.transport_normalized_l2())
      mix(result, bits(value));
    mix(result, material.mass_conservation_available() ? 1U : 0U);
    mix(result, bits(material.mass_relative_conservation_defect()));
    mix(result, material.transport_conservation_availability().size());
    for (const auto value : material.transport_conservation_availability())
      mix(result, value);
    mix(result, material.transport_relative_conservation_defect().size());
    for (const double value : material.transport_relative_conservation_defect())
      mix(result, bits(value));
    mix(result, material.minimum_density_available() ? 1U : 0U);
    mix(result, bits(material.minimum_density_kg_per_m3()));
    mix(result, material.minimum_density_global_cell());
    mix(result,
        static_cast<std::uint64_t>(material.minimum_density_rank() + 1));
  }
  mix(result, static_cast<std::uint64_t>(material_failure_reason_));
  mix(result, material_field_count_);
  mix(result, static_cast<std::uint64_t>(shared_face_mass_flux_field_));
  mix(result, static_cast<std::uint64_t>(flux_provenance_));
  mix(result, attempt_identity_);
  mix(result, material_attempt_identity_);
  mix(result, material_finalization_identity_);
  mix(result, final_continuity_residual_available_);
  mix(result, final_pressure_residual_available_);
  mix(result, bits(final_pressure_normalized_residual_));
  for (const auto value : final_momentum_residual_available_) {
    mix(result, value);
  }
  mix(result, mass_conservation_available_);
  for (const auto value : momentum_conservation_available_) {
    mix(result, value);
  }
  mix(result, closure_origin_ ? 1U : 0U);
  mix(result, pre_closure_authority_.has_value() ? 1U : 0U);
  if (pre_closure_authority_)
    mix(result, pre_closure_authority_->compute_seal());
  mix(result, post_closure_evidence_available_ ? 1U : 0U);
  mix(result, post_closure_report_.has_value() ? 1U : 0U);
  if (post_closure_report_)
    mix(result, post_closure_report_->compute_seal());
  mix(result, post_closure_authority_.has_value() ? 1U : 0U);
  if (post_closure_authority_)
    mix(result, post_closure_authority_->compute_seal());
  return result;
}

void MaterialDensityStepAttemptReport::seal() noexcept {
  seal_ = semantic_valid() ? compute_seal() : 0U;
}

bool MaterialDensityStepAttemptReport::semantic_valid() const noexcept {
  const auto positive_zero = [](double value) noexcept {
    return bits(value) == bits(0.0);
  };
  const auto byte = [](std::uint8_t value) noexcept { return value <= 1U; };
  const std::size_t fields = static_cast<std::size_t>(material_field_count_);
  const auto same_transport_report = [&](
                                         const MaterialDensityTransportReport &a,
                                         const MaterialDensityTransportReport &b) {
    const auto same_double = [](double left, double right) noexcept {
      return bits(left) == bits(right);
    };
    const auto same_doubles = [&](const std::vector<double> &left,
                                  const std::vector<double> &right) {
      return left.size() == right.size() &&
             std::equal(left.begin(), left.end(), right.begin(), same_double);
    };
    return a.authenticated() && b.authenticated() &&
           a.disposition() == b.disposition() && a.reason() == b.reason() &&
           a.lowest_failing_rank() == b.lowest_failing_rank() &&
           same_stencil(a.stencil(), b.stencil()) &&
           a.flux_provenance() == b.flux_provenance() &&
           a.attempt_identity() == b.attempt_identity() &&
           a.finalization_identity() == b.finalization_identity() &&
           a.shared_face_mass_flux_field() ==
               b.shared_face_mass_flux_field() &&
           a.density_residual_available() ==
               b.density_residual_available() &&
           same_double(a.density_normalized_l2(), b.density_normalized_l2()) &&
           a.transport_residual_availability() ==
               b.transport_residual_availability() &&
           same_doubles(a.transport_normalized_l2(),
                        b.transport_normalized_l2()) &&
           a.mass_conservation_available() ==
               b.mass_conservation_available() &&
           same_double(a.mass_relative_conservation_defect(),
                       b.mass_relative_conservation_defect()) &&
           a.transport_conservation_availability() ==
               b.transport_conservation_availability() &&
           same_doubles(a.transport_relative_conservation_defect(),
                        b.transport_relative_conservation_defect()) &&
           a.minimum_density_available() == b.minimum_density_available() &&
           same_double(a.minimum_density_kg_per_m3(),
                       b.minimum_density_kg_per_m3()) &&
           a.minimum_density_global_cell() ==
               b.minimum_density_global_cell() &&
           a.minimum_density_rank() == b.minimum_density_rank();
  };
  if (static_cast<std::uint8_t>(flow_.reason) >
          static_cast<std::uint8_t>(
              StepFailureReason::density_closure_failure) ||
      attempt_identity_ == 0U || material_field_count_ == 0U ||
      flow_.final_transport_normalized_l2.size() != fields ||
      flow_.final_transport_relative_conservation_defect.size() != fields ||
      !std::all_of(final_momentum_residual_available_.begin(),
                   final_momentum_residual_available_.end(), byte) ||
      !std::all_of(momentum_conservation_available_.begin(),
                   momentum_conservation_available_.end(), byte))
    return false;
  if (material_report_) {
    const auto &material = *material_report_;
    if (!material.authenticated() ||
        material.transport_residual_availability().size() != fields ||
        material.transport_normalized_l2().size() != fields ||
        material.transport_conservation_availability().size() != fields ||
        material.transport_relative_conservation_defect().size() != fields ||
        material.reason() != material_failure_reason_ ||
        material_attempt_identity_ == 0U ||
        material_finalization_identity_ == 0U ||
        material.attempt_identity() != material_attempt_identity_ ||
        material.finalization_identity() != material_finalization_identity_ ||
        material.shared_face_mass_flux_field() !=
            shared_face_mass_flux_field_ ||
        material.flux_provenance() != flux_provenance_ ||
        !std::all_of(material.transport_residual_availability().begin(),
                     material.transport_residual_availability().end(), byte) ||
        !std::all_of(material.transport_conservation_availability().begin(),
                     material.transport_conservation_availability().end(),
                     byte))
      return false;
    if (!material.density_residual_available() &&
        !positive_zero(material.density_normalized_l2()))
      return false;
    for (std::size_t field = 0; field < fields; ++field) {
      const bool residual_available =
          material.transport_residual_availability()[field] != 0U;
      if ((!closure_origin_ && residual_available &&
           bits(flow_.final_transport_normalized_l2[field]) !=
               bits(material.transport_normalized_l2()[field])) ||
          (!residual_available &&
           (!positive_zero(flow_.final_transport_normalized_l2[field]) ||
            !positive_zero(material.transport_normalized_l2()[field]))))
        return false;
      const bool conservation_available =
          material.transport_conservation_availability()[field] != 0U;
      if ((!closure_origin_ && conservation_available &&
           bits(flow_.final_transport_relative_conservation_defect[field]) !=
               bits(
                   material.transport_relative_conservation_defect()[field])) ||
          (!conservation_available &&
           (!positive_zero(
                flow_.final_transport_relative_conservation_defect[field]) ||
            !positive_zero(
                material.transport_relative_conservation_defect()[field]))))
        return false;
    }
    if ((!closure_origin_ &&
         mass_conservation_available_ !=
             material.mass_conservation_available()) ||
        (!closure_origin_ && mass_conservation_available_ &&
         bits(flow_.final_mass_relative_conservation_defect) !=
             bits(material.mass_relative_conservation_defect())) ||
        (!closure_origin_ && !mass_conservation_available_ &&
         (!positive_zero(flow_.final_mass_relative_conservation_defect) ||
          !positive_zero(material.mass_relative_conservation_defect()))))
      return false;
    if (!material.minimum_density_available()) {
      if (!positive_zero(material.minimum_density_kg_per_m3()) ||
          material.minimum_density_global_cell() != 0U ||
          material.minimum_density_rank() != -1)
        return false;
    } else if (!(material.minimum_density_kg_per_m3() > 0.0) ||
               !std::isfinite(material.minimum_density_kg_per_m3()) ||
               material.minimum_density_rank() < 0) {
      return false;
    }
  } else {
    if (material_attempt_identity_ != 0U ||
        material_finalization_identity_ != 0U || mass_conservation_available_ ||
        !positive_zero(flow_.final_mass_relative_conservation_defect) ||
        std::any_of(flow_.final_transport_normalized_l2.begin(),
                    flow_.final_transport_normalized_l2.end(),
                    [&](double value) { return !positive_zero(value); }) ||
        std::any_of(flow_.final_transport_relative_conservation_defect.begin(),
                    flow_.final_transport_relative_conservation_defect.end(),
                    [&](double value) { return !positive_zero(value); }))
      return false;
  }
  if (!final_continuity_residual_available_ &&
      !positive_zero(flow_.final_continuity_normalized_l2))
    return false;
  if (!final_pressure_residual_available_ &&
      (!positive_zero(flow_.final_pressure_residual_l2) ||
       !positive_zero(final_pressure_normalized_residual_)))
    return false;
  for (std::size_t component_index = 0; component_index < 3U;
       ++component_index) {
    if (final_momentum_residual_available_[component_index] == 0U &&
        !positive_zero(flow_.final_momentum_normalized_l2[component_index]))
      return false;
    if (momentum_conservation_available_[component_index] == 0U &&
        !positive_zero(
            flow_.final_momentum_relative_conservation_defect[component_index]))
      return false;
  }
  if (!mass_conservation_available_ &&
      !positive_zero(flow_.final_mass_relative_conservation_defect))
    return false;

  const bool committed = flow_.disposition == StepAttemptDisposition::committed;
  if (committed) {
    if (flow_.reason != StepFailureReason::none ||
        flow_.lowest_failing_rank != -1 ||
        flow_.pressure_corrector_count != 2U || !material_report_ ||
        material_report_->disposition() !=
            MaterialTransportDisposition::finalized ||
        material_report_->reason() != MaterialTransportFailureReason::none ||
        material_failure_reason_ != MaterialTransportFailureReason::none ||
        flux_provenance_ != MaterialFluxProvenance::final_corrected ||
        material_report_->flux_provenance() !=
            MaterialFluxProvenance::final_corrected ||
        !material_report_->density_residual_available() ||
        !std::all_of(
            material_report_->transport_residual_availability().begin(),
            material_report_->transport_residual_availability().end(),
            [](std::uint8_t value) { return value == 1U; }) ||
        !material_report_->mass_conservation_available() ||
        !std::all_of(
            material_report_->transport_conservation_availability().begin(),
            material_report_->transport_conservation_availability().end(),
            [](std::uint8_t value) { return value == 1U; }) ||
        !material_report_->minimum_density_available() ||
        !final_continuity_residual_available_ ||
        !final_pressure_residual_available_ || !mass_conservation_available_ ||
        !std::all_of(final_momentum_residual_available_.begin(),
                     final_momentum_residual_available_.end(),
                     [](std::uint8_t value) { return value == 1U; }) ||
        !std::all_of(momentum_conservation_available_.begin(),
                     momentum_conservation_available_.end(),
                     [](std::uint8_t value) { return value == 1U; }))
      return false;
  } else if (flow_.reason == StepFailureReason::none ||
             flow_.disposition == StepAttemptDisposition::committed ||
             flow_.lowest_failing_rank < 0) {
    return false;
  }

  if (!committed) {
    const bool non_retryable =
        flow_.reason == StepFailureReason::invalid_input ||
        flow_.reason == StepFailureReason::collective_operation;
    if (flow_.disposition !=
            (non_retryable ? StepAttemptDisposition::non_retryable_failure
                           : StepAttemptDisposition::recoverable_failure) ||
        bits(flow_.suggested_dt_s) != bits(non_retryable
                                               ? flow_.attempted_dt_s
                                               : 0.5 * flow_.attempted_dt_s) ||
        flow_.pressure_corrector_count > 2U)
      return false;
  }
  if (material_report_ &&
      (flow_.pressure_corrector_count != 2U ||
       flux_provenance_ != MaterialFluxProvenance::final_corrected ||
       material_report_->flux_provenance() !=
           MaterialFluxProvenance::final_corrected))
    return false;
  if (final_pressure_residual_available_ &&
      !final_continuity_residual_available_)
    return false;
  for (std::size_t component_index = 0; component_index < 3U; ++component_index)
    if (momentum_conservation_available_[component_index] != 0U &&
        final_momentum_residual_available_[component_index] == 0U)
      return false;

  if (flow_.reason == StepFailureReason::collective_operation &&
      flow_.lowest_failing_rank < 0)
    return false;
  const bool post_failure_authority =
      closure_origin_ && post_closure_report_ &&
      post_closure_report_->disposition() !=
          MaterialTransportDisposition::finalized &&
      map_material_failure(post_closure_report_->reason()) == flow_.reason &&
      post_closure_report_->lowest_failing_rank() ==
          flow_.lowest_failing_rank;
  if (flow_.reason == StepFailureReason::transport_failure &&
      material_failure_reason_ !=
          MaterialTransportFailureReason::non_finite_state &&
      material_failure_reason_ !=
          MaterialTransportFailureReason::non_positive_density &&
      !post_failure_authority)
    return false;
  if (material_failure_reason_ != MaterialTransportFailureReason::none) {
    if (map_material_failure(material_failure_reason_) != flow_.reason)
      return false;
  } else if ((flow_.reason == StepFailureReason::transport_failure &&
              !post_failure_authority) ||
             (material_report_ && material_report_->reason() !=
                                      MaterialTransportFailureReason::none)) {
    return false;
  }
  if (!closure_origin_ &&
      flow_.reason == StepFailureReason::density_closure_failure)
    return false;
  if (closure_origin_) {
    if (material_report_.has_value() != pre_closure_authority_.has_value() ||
        (material_report_ &&
         !same_transport_report(*material_report_,
                                *pre_closure_authority_)) ||
        post_closure_report_.has_value() !=
            post_closure_authority_.has_value())
      return false;
    if (post_closure_evidence_available_ !=
            post_closure_report_.has_value() ||
        (committed && !post_closure_evidence_available_))
      return false;
    if (post_closure_report_) {
      const auto &post = *post_closure_report_;
      if (!post.authenticated() ||
          !same_transport_report(post, *post_closure_authority_) ||
          post.attempt_identity() != attempt_identity_ ||
          post.shared_face_mass_flux_field() !=
              shared_face_mass_flux_field_ ||
          post.flux_provenance() != MaterialFluxProvenance::final_corrected ||
          !material_report_ ||
          !same_stencil(post.stencil(), material_report_->stencil()) ||
          post.finalization_identity() == 0U ||
          post.finalization_identity() <= material_finalization_identity_ ||
          post.transport_residual_availability().size() != fields ||
          post.transport_conservation_availability().size() != fields ||
          post.transport_normalized_l2().size() != fields ||
          post.transport_relative_conservation_defect().size() != fields ||
          bits(flow_.final_mass_relative_conservation_defect) !=
              bits(post.mass_relative_conservation_defect()) ||
          mass_conservation_available_ != post.mass_conservation_available())
        return false;
      if (post.disposition() != MaterialTransportDisposition::finalized &&
          (flow_.reason != map_material_failure(post.reason()) ||
           flow_.lowest_failing_rank != post.lowest_failing_rank()))
        return false;
      for (std::size_t field = 0; field < fields; ++field)
        if (bits(flow_.final_transport_normalized_l2[field]) !=
                bits(post.transport_normalized_l2()[field]) ||
            bits(flow_.final_transport_relative_conservation_defect[field]) !=
                bits(post.transport_relative_conservation_defect()[field]))
          return false;
      if (committed &&
          (post.disposition() != MaterialTransportDisposition::finalized ||
           post.reason() != MaterialTransportFailureReason::none ||
           !post.density_residual_available() ||
           !post.mass_conservation_available() ||
           !std::all_of(post.transport_residual_availability().begin(),
                        post.transport_residual_availability().end(),
                        [](std::uint8_t value) { return value == 1U; }) ||
           !std::all_of(post.transport_conservation_availability().begin(),
                        post.transport_conservation_availability().end(),
                        [](std::uint8_t value) { return value == 1U; })))
        return false;
    }
  } else if (pre_closure_authority_ || post_closure_evidence_available_ ||
             post_closure_report_ || post_closure_authority_) {
    return false;
  }
  return true;
}

bool MaterialDensityStepAttemptReport::authenticated() const noexcept {
  return seal_ != 0U && seal_ == compute_seal() && semantic_valid();
}

FixedStepMaterialDensityFlow::FixedStepMaterialDensityFlow(
    std::unique_ptr<Impl> implementation) noexcept
    : impl_(std::move(implementation)) {}

FixedStepMaterialDensityFlow FixedStepMaterialDensityFlow::create(
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
    const runtime::FieldRegistry &registry, FlowFieldIds fields,
    MaterialDensityTransportSpec specification) {
  validate_host_context(execution_context);
  geometry.require_compatible(topology);
  if (geometry.mapping_kind() != mesh::MappingKind::uniform_box)
    throw runtime::Error("material flow supports uniform geometry only");
  if (!cell_halo.is_compatible_with(decomposition) ||
      cell_halo.ghost_width() != 2)
    throw runtime::Error("material flow requires a compatible width-two Halo");
  if (std::any_of(momentum_preconditioners.begin(),
                  momentum_preconditioners.end(),
                  [](const auto *value) { return value == nullptr; }))
    throw runtime::Error("material momentum preconditioner is null");
  if (1U + specification.scalar_densities.size() >
      std::numeric_limits<std::uint64_t>::max())
    throw runtime::Error("material field count overflows");
  for (mesh::LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    const auto patch = topology.patch_id(face);
    if (patch && boundaries.patch(*patch).kind() !=
                     boundary::BoundaryKind::periodic)
      throw runtime::Error("material flow requires fully periodic boundaries");
  }
  auto coupler = PisoCoupler::create(
      decomposition, topology, geometry, boundaries, mpi, execution_context,
      cell_halo, pressure_solver, pressure_preconditioner);
  coupler.prepare_material_density_assessment();
  auto transport = MaterialDensityTransport::create(
      registry, decomposition, topology, geometry, boundaries, mpi, cell_halo,
      fields, specification);
  return FixedStepMaterialDensityFlow(std::make_unique<Impl>(
      decomposition, topology, geometry, boundaries, mpi, execution_context,
      cell_halo, momentum_solver, momentum_preconditioners,
      std::move(coupler), std::move(transport), registry, std::move(fields),
      std::move(specification)));
}

FixedStepMaterialDensityFlow FixedStepMaterialDensityFlow::create_open_capable(
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
    const runtime::FieldRegistry &registry, FlowFieldIds fields,
    MaterialDensityTransportSpec specification) {
  validate_host_context(execution_context);
  geometry.require_compatible(topology);
  if (geometry.mapping_kind() != mesh::MappingKind::uniform_box)
    throw runtime::Error("closure flow supports uniform geometry only");
  if (!cell_halo.is_compatible_with(decomposition) ||
      cell_halo.ghost_width() != 2)
    throw runtime::Error("closure flow requires a compatible width-two Halo");
  if (std::any_of(momentum_preconditioners.begin(),
                  momentum_preconditioners.end(),
                  [](const auto *value) { return value == nullptr; }))
    throw runtime::Error("closure momentum preconditioner is null");
  auto coupler = PisoCoupler::create(
      decomposition, topology, geometry, boundaries, mpi, execution_context,
      cell_halo, pressure_solver, pressure_preconditioner);
  coupler.prepare_material_density_assessment();
  auto transport = MaterialDensityTransport::create(
      registry, decomposition, topology, geometry, boundaries, mpi, cell_halo,
      fields, specification);
  return FixedStepMaterialDensityFlow(std::make_unique<Impl>(
      decomposition, topology, geometry, boundaries, mpi, execution_context,
      cell_halo, momentum_solver, momentum_preconditioners,
      std::move(coupler), std::move(transport), registry, std::move(fields),
      std::move(specification)));
}

FixedStepMaterialDensityFlow detail::DensityClosureBridge::create_open_capable(
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
    const runtime::FieldRegistry &registry, FlowFieldIds fields,
    MaterialDensityTransportSpec specification) {
  return FixedStepMaterialDensityFlow::create_open_capable(
      decomposition, topology, geometry, boundaries, mpi, execution_context,
      cell_halo, momentum_solver, momentum_preconditioners, pressure_solver,
      pressure_preconditioner, registry, std::move(fields),
      std::move(specification));
}

FixedStepMaterialDensityFlow::~FixedStepMaterialDensityFlow() noexcept =
    default;

FixedStepMaterialDensityFlow::FixedStepMaterialDensityFlow(
    FixedStepMaterialDensityFlow &&other) noexcept
    : impl_(std::move(other.impl_)) {
  if (impl_) {
    ++impl_->source_generation;
    if (impl_->source_generation == 0U)
      impl_->source_generation = 1U;
  }
}

MaterialDensityStepAttemptReport FixedStepMaterialDensityFlow::attempt(
    FlowState &state, double mu, const MomentumTimeStencil &stencil,
    const linear::SolveControl &momentum_control,
    const linear::SolveControl &pressure_control) const {
  return attempt_common(state, mu, stencil, momentum_control,
                        pressure_control, nullptr);
}

MaterialDensityStepAttemptReport
detail::DensityClosureBridge::attempt(
    const FixedStepMaterialDensityFlow &flow,
    FlowState &state, double mu, const MomentumTimeStencil &stencil,
    const linear::SolveControl &momentum_control,
    const linear::SolveControl &pressure_control,
    const DensityClosureHooks &hooks) {
  return flow.attempt_common(state, mu, stencil, momentum_control,
                             pressure_control, &hooks);
}

MaterialDensityStepAttemptReport FixedStepMaterialDensityFlow::attempt_common(
    FlowState &state, double mu, const MomentumTimeStencil &stencil,
    const linear::SolveControl &momentum_control,
    const linear::SolveControl &pressure_control,
    const detail::DensityClosureHooks *closure) const {
  if (!impl_)
    throw runtime::Error("material flow object has been moved from");
  impl_->last_state = nullptr;
  impl_->last_report_seal = 0U;
  impl_->last_state_identity = 0U;
  if (impl_->attempt_identity == std::numeric_limits<std::uint64_t>::max())
    throw runtime::Error("material flow attempt identity would wrap");
  ++impl_->attempt_identity;

  MaterialDensityStepAttemptReport result;
  result.flow_ = base_flow_report(
      stencil.dt_s, static_cast<std::size_t>(impl_->material_field_count));
  result.material_field_count_ = impl_->material_field_count;
  result.shared_face_mass_flux_field_ = impl_->fields.face_mass_flux;
  result.flux_provenance_ = MaterialFluxProvenance::predictor;
  result.attempt_identity_ = impl_->attempt_identity;
  result.closure_origin_ = closure != nullptr;
  const auto finish = [&]() {
    result.seal();
    impl_->last_state = &state;
    impl_->last_report_seal = result.seal_;
    impl_->last_state_identity = state.diagnostic_mutation_identity();
    return result;
  };
  bool preparation_ok = true;
  try {
    impl_->material_transport.prepare_task20_attempt();
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    if (preflight_allocation_failure_rank.load(std::memory_order_relaxed) ==
        impl_->mpi->rank())
      throw std::bad_alloc();
#endif
  } catch (const runtime::MpiOperationError &) {
    throw;
  } catch (...) {
    preparation_ok = false;
  }
  const auto preparation = runtime::collective_status(
      *impl_->mpi, preparation_ok,
      "material flow attempt workspace preparation failed");
  if (!preparation.ok) {
    set_failure(result.flow_, StepFailureReason::invalid_input,
                preparation.failing_rank, false);
    return finish();
  }
  bool active = false;
  const auto fail = [&](StepFailureReason reason, int rank, bool recoverable) {
    if (closure != nullptr)
      closure->rollback(closure->object);
    if (active) {
      state.rollback_attempt();
      active = false;
    }
    set_failure(result.flow_, reason, rank, recoverable);
    return finish();
  };
  const auto material_failure = [&](MaterialTransportFailureReason reason,
                                    int rank) {
    require_reliable_collective_result(map_material_failure(reason), rank);
    result.material_failure_reason_ = reason;
    const StepFailureReason step_reason = map_material_failure(reason);
    const bool recoverable =
        reason != MaterialTransportFailureReason::invalid_input &&
        reason != MaterialTransportFailureReason::collective_operation;
    return fail(step_reason, rank, recoverable);
  };

  try {
    bool valid = true;
    try {
      const auto expected = make_momentum_time_stencil(
          stencil.order, stencil.dt_s, stencil.previous_dt_s);
      const auto metadata = state.metadata();
      const auto layout = state.layer(FlowLayer::committed).layout_set();
      valid = mu >= 0.0 && std::isfinite(mu) && !state.attempt_active() &&
              &detail::FlowStateSolverAccess::registry(state) ==
                  impl_->registry &&
              state.fields().density == impl_->fields.density &&
              state.fields().velocity == impl_->fields.velocity &&
              state.fields().mechanical_pressure ==
                  impl_->fields.mechanical_pressure &&
              state.fields().face_velocity == impl_->fields.face_velocity &&
              state.fields().face_mass_flux ==
                  impl_->fields.face_mass_flux &&
              state.fields().transported_cell_fields ==
                  impl_->fields.transported_cell_fields &&
              same(layout.cell_interior_extent,
                   impl_->decomposition->local_extent()) &&
              layout.face_count == impl_->topology->local_face_count() &&
              expected.alpha0 == stencil.alpha0 &&
              expected.alpha1 == stencil.alpha1 &&
              expected.alpha2 == stencil.alpha2 &&
              (stencil.order != MomentumTimeOrder::bdf2 ||
               stencil.previous_dt_s == metadata.dt_s) &&
              std::isfinite(momentum_control.atol) &&
              momentum_control.atol >= 0.0 &&
              std::isfinite(momentum_control.rtol) &&
              momentum_control.rtol >= 0.0 &&
              momentum_control.residual_recompute_interval != 0U &&
              std::isfinite(pressure_control.atol) &&
              pressure_control.atol >= 0.0 &&
              std::isfinite(pressure_control.rtol) &&
              pressure_control.rtol >= 0.0 &&
              pressure_control.residual_recompute_interval != 0U &&
              (closure == nullptr ||
               std::isfinite(closure->enthalpy_rate_J_per_kg_s));
    } catch (const runtime::MpiOperationError &) {
      throw;
    } catch (...) {
      valid = false;
    }
    const auto validation = runtime::collective_status(
        *impl_->mpi, valid, "material flow input is invalid");
    if (!validation.ok)
      return fail(StepFailureReason::invalid_input, validation.failing_rank,
                  false);

    state.begin_attempt();
    active = true;
    if (closure != nullptr)
      closure->begin(closure->object, state, result.attempt_identity_);
    auto &history =
        detail::FlowStateSolverAccess::layer(state, FlowLayer::history);
    auto &committed =
        detail::FlowStateSolverAccess::layer(state, FlowLayer::committed);
    auto &trial = detail::FlowStateSolverAccess::layer(state, FlowLayer::trial);
    const auto &access = detail::FlowStateSolverAccess::access(state);
    const auto &fields = state.fields();
    const auto bind_face_flux = [&](const runtime::FieldStorage &storage) {
      return finite_volume::FaceMassFlux::bind_prepared(
          *impl_->prepared_face_flux, *impl_->registry, storage, access,
          kStatePhase, kStateActor, fields.face_mass_flux, *impl_->topology);
    };
    const auto bind_material_flux = [&](MaterialFluxProvenance provenance) {
      return MaterialFaceMassFlux::bind_prepared(
          *impl_->prepared_material_flux, *impl_->registry, trial, access,
          kStatePhase, kStateActor, fields.face_mass_flux, *impl_->topology,
          provenance);
    };

    for (FlowLayer layer : {FlowLayer::history, FlowLayer::committed}) {
      auto &storage = detail::FlowStateSolverAccess::layer(state, layer);
      impl_->halo->exchange(storage, fields.density);
      impl_->halo->exchange(storage, fields.velocity);
      impl_->halo->exchange(storage, fields.mechanical_pressure);
      for (const auto field : fields.transported_cell_fields)
        impl_->halo->exchange(storage, field);
    }

    synchronized_local_phase(*impl_->mpi, StepFailureReason::non_finite_trial,
                             true, [&] {
      const auto current = committed.acquire_face_read<double>(
          access, kStatePhase, kStateActor, fields.face_mass_flux);
      const auto previous = history.acquire_face_read<double>(
          access, kStatePhase, kStateActor, fields.face_mass_flux);
      auto predictor_flux = trial.acquire_face_write<double>(
          access, kStatePhase, kStateActor, fields.face_mass_flux);
      for (mesh::LocalFaceId face = 0;
           face < impl_->topology->local_face_count(); ++face) {
        predictor_flux(face, 0) =
            stencil.order == MomentumTimeOrder::backward_euler
                ? current(face, 0)
                : 2.0 * current(face, 0) - previous(face, 0);
        if (!std::isfinite(predictor_flux(face, 0)))
          throw runtime::Error("material predictor flux is non-finite");
      }
    });
    MaterialDensityTransport::StagingResult predictor_stage;
    {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      FaceFluxPathObservation allocation_observation;
#endif
      std::optional<MaterialFaceMassFlux> predictor_flux;
      synchronized_local_phase(
          *impl_->mpi, StepFailureReason::invalid_input, false, [&] {
            predictor_flux.emplace(
                bind_material_flux(MaterialFluxProvenance::predictor));
          });
      predictor_stage = impl_->material_transport.stage_trial(
          state, *predictor_flux, stencil,
          closure == nullptr ? 0.0 : closure->enthalpy_rate_J_per_kg_s);
      predictor_flux.reset();
    }
    if (predictor_stage.reason != MaterialTransportFailureReason::none) {
      require_reliable_collective_result(
          map_material_failure(predictor_stage.reason),
          predictor_stage.lowest_failing_rank);
      return material_failure(predictor_stage.reason,
                              predictor_stage.lowest_failing_rank);
    }
    if (closure != nullptr) {
      const auto closure_result = closure->evaluate(
          closure->object, state, detail::DensityClosureStage::predictor);
      if (!closure_result.accepted)
        return fail(StepFailureReason::density_closure_failure,
                    closure_result.lowest_failing_rank,
                    closure_result.recoverable);
      impl_->halo->exchange(trial, fields.density);
      impl_->halo->exchange(trial, closure->enthalpy_density);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      if (closure->after_halo != nullptr)
        closure->after_halo(closure->object,
                            detail::DensityClosureStage::predictor,
                            fields.density, closure->enthalpy_density);
#endif
    } else {
      impl_->halo->exchange(trial, fields.density);
    }

    assemble_spatial(*impl_, access, committed, fields, mu, impl_->momentum_n,
                     bind_face_flux);
    if (stencil.order == MomentumTimeOrder::bdf2)
      assemble_spatial(*impl_, access, history, fields, mu, impl_->momentum_nm1,
                       bind_face_flux);

    const auto rho_trial = trial.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.density);
    const auto rho_n = committed.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.density);
    const auto rho_nm1 = history.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.density);
    const auto velocity_n = committed.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.velocity);
    const auto velocity_nm1 = history.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.velocity);
    const auto pressure_n = committed.acquire_read<double>(
        access, kStatePhase, kStateActor, fields.mechanical_pressure);
    synchronized_local_phase(*impl_->mpi, StepFailureReason::invalid_input,
                             false, [&] {
      auto pressure_gradient = impl_->scratch.storage->acquire_write<double>(
          *impl_->scratch.access, kScratchPhase, kScratchActor,
          impl_->scratch.pressure_gradient);
      compute_pressure_gradient(*impl_->topology, *impl_->geometry, pressure_n,
                                pressure_gradient,
                                impl_->pressure_gradient_sums);
    });
    impl_->halo->exchange(*impl_->scratch.storage,
                          impl_->scratch.pressure_gradient);
    const auto pressure_gradient_read =
        impl_->scratch.storage->acquire_read<double>(
            *impl_->scratch.access, kScratchPhase, kScratchActor,
            impl_->scratch.pressure_gradient);

    const std::size_t count = impl_->topology->owned_cell_count();
    const auto owned = impl_->topology->owned_global_box();
    for (std::size_t direction = 0; direction < 3U; ++direction) {
      LocalMaterialFailure assembly_failure{};
      try {
        auto rhs = impl_->rhs[direction].view(0U, count);
        auto predictor = impl_->predictor[direction].view(0U, count);
        auto diagonal_output = impl_->diagonal[direction].view(0U, count);
        for (mesh::LocalCellId cell = 0; cell < count; ++cell) {
          const auto index = map_cell(impl_->topology->global_cell(cell), owned,
                                      impl_->topology->global_extent());
          const double volume = impl_->geometry->cell_volume_m3(cell);
          const double diagonal = stencil.alpha0 *
                                  at(rho_trial, index, 0) * volume /
                                  stencil.dt_s;
          const double spatial =
              stencil.order == MomentumTimeOrder::backward_euler
                  ? impl_->momentum_n.convection[cell * 3U + direction] +
                        impl_->momentum_n.viscosity[cell * 3U + direction]
                  : 2.0 *
                            (impl_->momentum_n
                                 .convection[cell * 3U + direction] +
                             impl_->momentum_n
                                 .viscosity[cell * 3U + direction]) -
                        (impl_->momentum_nm1
                             .convection[cell * 3U + direction] +
                         impl_->momentum_nm1
                             .viscosity[cell * 3U + direction]);
          double source = 0.0;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
          if (impl_->vortex_source_enabled) {
            const auto centre = impl_->geometry->cell_center_m(cell);
            source = component(analytic_source(centre.x, centre.y, mu),
                               static_cast<int>(direction));
          }
#endif
          if (!std::isfinite(source)) {
            assembly_failure = {true, StepFailureReason::non_finite_trial,
                                true};
            break;
          }
          if (!(diagonal > 0.0) || !std::isfinite(diagonal)) {
            assembly_failure = {true, StepFailureReason::invalid_input, false};
            break;
          }
          rhs[cell] =
              -(volume / stencil.dt_s) *
                  (stencil.alpha1 * at(rho_n, index, 0) *
                       at(velocity_n, index, static_cast<int>(direction)) +
                   stencil.alpha2 * at(rho_nm1, index, 0) *
                       at(velocity_nm1, index, static_cast<int>(direction))) -
              spatial - volume *
                            at(pressure_gradient_read, index,
                               static_cast<int>(direction)) +
              volume * source;
          predictor[cell] =
              at(velocity_n, index, static_cast<int>(direction));
          impl_->diagonal_values[direction][cell] = diagonal;
          diagonal_output[cell] = diagonal;
          if (!std::isfinite(rhs[cell]) || !std::isfinite(predictor[cell])) {
            assembly_failure = {true, StepFailureReason::non_finite_trial,
                                true};
            break;
          }
        }
      } catch (const runtime::MpiOperationError &) {
        throw;
      } catch (...) {
        assembly_failure = {true, StepFailureReason::invalid_input, false};
      }
      const auto selected =
          synchronize_local_failure(*impl_->mpi, assembly_failure);
      if (selected.failing_rank >= 0)
        throw selected;
      synchronized_local_phase(*impl_->mpi, StepFailureReason::invalid_input,
                               false, [&] {
        impl_->operators[direction]->replace(
            impl_->diagonal_values[direction]);
        impl_->momentum_preconditioners[direction]->update(
            *impl_->operators[direction],
            impl_->operators[direction]->revision());
      });

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      apply_terminal_test_point(
          static_cast<test::MaterialTerminalPointForTest>(
              static_cast<int>(
                  test::MaterialTerminalPointForTest::momentum_x) +
              static_cast<int>(direction)));
      if (closure != nullptr && direction == 0U &&
          closure->outer_failure != nullptr) {
        const int rank = closure->outer_failure(
            closure->object,
            detail::DensityClosureOuterPoint::momentum_after_predictor);
        if (rank >= 0)
          return fail(StepFailureReason::invalid_input, rank, false);
      }
#endif
      result.flow_.momentum.components[direction] =
          impl_->momentum_solver->solve(
              *impl_->operators[direction],
              *impl_->momentum_preconditioners[direction],
              static_cast<const execution::Buffer &>(impl_->rhs[direction])
                  .view(0U, count),
              impl_->predictor[direction].view(0U, count), momentum_control);
      const auto &solve = result.flow_.momentum.components[direction];
      if (solve.reason == linear::SolveTerminationReason::collective_failure) {
        require_reliable_collective_result(
            StepFailureReason::collective_operation,
            solve.lowest_failing_rank);
        return fail(StepFailureReason::collective_operation,
                    solve.lowest_failing_rank, false);
      }
      const auto solve_status = synchronize_local_failure(
          *impl_->mpi,
          {!solve_success(solve.reason),
           StepFailureReason::momentum_linear_solve, true});
      if (solve_status.failing_rank >= 0)
        return fail(solve_status.reason, solve_status.failing_rank,
                    solve_status.recoverable);
    }

    synchronized_local_phase(*impl_->mpi, StepFailureReason::non_finite_trial,
                             true, [&] {
      auto velocity = trial.acquire_write<double>(
          access, kStatePhase, kStateActor, fields.velocity);
      auto actual_diagonal = impl_->scratch.storage->acquire_write<double>(
          *impl_->scratch.access, kScratchPhase, kScratchActor,
          impl_->scratch.actual_diagonal);
      for (mesh::LocalCellId cell = 0; cell < count; ++cell) {
        const auto index = map_cell(impl_->topology->global_cell(cell), owned,
                                    impl_->topology->global_extent());
        for (std::size_t direction = 0; direction < 3U; ++direction) {
          velocity(index.i, index.j, index.k, static_cast<int>(direction)) =
              impl_->predictor[direction].view(0U, count)[cell];
          actual_diagonal(index.i, index.j, index.k,
                          static_cast<int>(direction)) =
              impl_->diagonal[direction].view(0U, count)[cell];
          if (!std::isfinite(
                  velocity(index.i, index.j, index.k,
                           static_cast<int>(direction))))
            throw runtime::Error("material predictor velocity is non-finite");
        }
      }
    });
    impl_->halo->exchange(trial, fields.velocity);
    impl_->halo->exchange(*impl_->scratch.storage,
                          impl_->scratch.actual_diagonal);
    {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      FaceFluxPathObservation allocation_observation;
#endif
      reconstruct_density_faces(*impl_, access, committed, fields,
                                impl_->scratch.face_density_n, bind_face_flux);
      if (stencil.order == MomentumTimeOrder::bdf2)
        reconstruct_density_faces(*impl_, access, history, fields,
                                  impl_->scratch.face_density_nm1,
                                  bind_face_flux);
    }
    const auto face_density_n =
        impl_->scratch.storage->acquire_face_read<double>(
            *impl_->scratch.access, kScratchPhase, kScratchActor,
            impl_->scratch.face_density_n);
    const auto face_density_nm1 =
        impl_->scratch.storage->acquire_face_read<double>(
            *impl_->scratch.access, kScratchPhase, kScratchActor,
            impl_->scratch.face_density_nm1);
    const auto face_velocity_n = committed.acquire_face_read<double>(
        access, kStatePhase, kStateActor, fields.face_velocity);
    const auto face_velocity_nm1 = history.acquire_face_read<double>(
        access, kStatePhase, kStateActor, fields.face_velocity);
    const auto actual_diagonal = impl_->scratch.storage->acquire_read<double>(
        *impl_->scratch.access, kScratchPhase, kScratchActor,
        impl_->scratch.actual_diagonal);
    synchronized_local_phase(*impl_->mpi, StepFailureReason::non_finite_trial,
                             true, [&] {
      auto face_velocity = trial.acquire_face_write<double>(
          access, kStatePhase, kStateActor, fields.face_velocity);
      const auto predictor_velocity = trial.acquire_read<double>(
          access, kStatePhase, kStateActor, fields.velocity);
      const MaterialMomentumFaceHistory face_history{
          rho_n,
          velocity_n,
          face_density_n,
          face_velocity_n,
          stencil.order == MomentumTimeOrder::bdf2 ? &rho_nm1 : nullptr,
          stencil.order == MomentumTimeOrder::bdf2 ? &velocity_nm1 : nullptr,
          stencil.order == MomentumTimeOrder::bdf2 ? &face_density_nm1
                                                   : nullptr,
          stencil.order == MomentumTimeOrder::bdf2 ? &face_velocity_nm1
                                                   : nullptr};
      impl_->face_assembler.assemble_material_density(
          *impl_->boundaries, stencil, predictor_velocity, pressure_n,
          pressure_gradient_read, actual_diagonal, face_history,
          face_velocity);
      for (mesh::LocalFaceId face = 0;
           face < impl_->topology->local_face_count(); ++face)
        for (int component_index = 0; component_index < 3;
             ++component_index)
          if (!std::isfinite(face_velocity(face, component_index)))
            throw runtime::Error("material face velocity is non-finite");
    });
    {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      FaceFluxPathObservation allocation_observation;
#endif
      build_density_weighted_flux(*impl_, state, access, fields,
                                  bind_face_flux);
    }
    synchronize_partition_faces(*impl_, state, access, fields);

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    apply_terminal_test_point(
        test::MaterialTerminalPointForTest::pressure_corrector_one);
#endif
    const auto first = impl_->coupler.correct_material_density(
        state, stencil, actual_diagonal, pressure_control);
    result.flow_.pressure[0] = first.solve;
    if (!first.accepted) {
      require_reliable_collective_result(first.reason,
                                         first.lowest_failing_rank);
      return fail(first.reason, first.lowest_failing_rank,
                  first.disposition ==
                      PressureCorrectionDisposition::recoverable_failure);
    }
    result.flow_.pressure_corrector_count = 1U;
    synchronize_partition_faces(*impl_, state, access, fields);

    MaterialDensityTransport::StagingResult provisional_stage;
    {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      FaceFluxPathObservation allocation_observation;
#endif
      std::optional<MaterialFaceMassFlux> provisional_flux;
      synchronized_local_phase(
          *impl_->mpi, StepFailureReason::invalid_input, false, [&] {
            provisional_flux.emplace(
                bind_material_flux(MaterialFluxProvenance::provisional));
          });
      provisional_stage = impl_->material_transport.stage_trial(
          state, *provisional_flux, stencil,
          closure == nullptr ? 0.0 : closure->enthalpy_rate_J_per_kg_s);
      provisional_flux.reset();
    }
    if (provisional_stage.reason != MaterialTransportFailureReason::none) {
      require_reliable_collective_result(
          map_material_failure(provisional_stage.reason),
          provisional_stage.lowest_failing_rank);
      return material_failure(provisional_stage.reason,
                              provisional_stage.lowest_failing_rank);
    }
    if (closure != nullptr) {
      const auto closure_result = closure->evaluate(
          closure->object, state, detail::DensityClosureStage::provisional);
      if (!closure_result.accepted)
        return fail(StepFailureReason::density_closure_failure,
                    closure_result.lowest_failing_rank,
                    closure_result.recoverable);
      impl_->halo->exchange(trial, fields.density);
      impl_->halo->exchange(trial, closure->enthalpy_density);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      if (closure->after_halo != nullptr)
        closure->after_halo(closure->object,
                            detail::DensityClosureStage::provisional,
                            fields.density, closure->enthalpy_density);
#endif
    } else {
      impl_->halo->exchange(trial, fields.density);
    }

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    apply_terminal_test_point(
        test::MaterialTerminalPointForTest::pressure_corrector_two);
    if (closure != nullptr && closure->outer_failure != nullptr) {
      const int rank = closure->outer_failure(
          closure->object,
          detail::DensityClosureOuterPoint::pressure_after_first_corrector);
      if (rank >= 0)
        return fail(StepFailureReason::invalid_input, rank, false);
    }
#endif
    const auto second = impl_->coupler.correct_material_density(
        state, stencil, actual_diagonal, pressure_control);
    result.flow_.pressure[1] = second.solve;
    if (!second.accepted) {
      require_reliable_collective_result(second.reason,
                                         second.lowest_failing_rank);
      return fail(second.reason, second.lowest_failing_rank,
                  second.disposition ==
                      PressureCorrectionDisposition::recoverable_failure);
    }
    result.flow_.pressure_corrector_count = 2U;
    synchronize_partition_faces(*impl_, state, access, fields);

    {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      FaceFluxPathObservation allocation_observation;
#endif
      std::optional<MaterialFaceMassFlux> final_material_flux;
      synchronized_local_phase(
          *impl_->mpi, StepFailureReason::invalid_input, false, [&] {
            final_material_flux.emplace(
                bind_material_flux(MaterialFluxProvenance::final_corrected));
          });
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      const auto flux_values = trial.acquire_face_read<double>(
          access, kStatePhase, kStateActor, fields.face_mass_flux);
      for (mesh::LocalFaceId face = 0;
           face < impl_->topology->local_face_count(); ++face)
        impl_->finalizer_flux_input[face] = flux_values(face, 0);
#endif
      result.material_report_.emplace(
          closure == nullptr
              ? impl_->material_transport.finalize_trial(
                    state, *final_material_flux, stencil)
              : impl_->material_transport.finalize_trial_with_source(
                    state, *final_material_flux, stencil,
                    closure->enthalpy_rate_J_per_kg_s));
      if (closure != nullptr) {
        // Ideal-origin nested evidence belongs to the facade attempt even
        // when earlier input rejections made FlowState's private sequence
        // differ from the facade sequence.
        result.material_report_->attempt_identity_ = result.attempt_identity_;
        result.material_report_->seal();
      }
      result.material_attempt_identity_ =
          result.material_report_->attempt_identity();
      result.material_finalization_identity_ =
          result.material_report_->finalization_identity();
      if (closure != nullptr)
        result.pre_closure_authority_ = result.material_report_;
      final_material_flux.reset();
    }
    result.flux_provenance_ = MaterialFluxProvenance::final_corrected;
    result.material_failure_reason_ = result.material_report_->reason();
    result.flow_.final_transport_normalized_l2 =
        result.material_report_->transport_normalized_l2();
    result.flow_.final_transport_relative_conservation_defect =
        result.material_report_->transport_relative_conservation_defect();
    result.flow_.final_mass_relative_conservation_defect =
        result.material_report_->mass_relative_conservation_defect();
    result.mass_conservation_available_ =
        result.material_report_->mass_conservation_available();
    if (result.material_report_->disposition() !=
        MaterialTransportDisposition::finalized) {
      require_reliable_collective_result(
          map_material_failure(result.material_report_->reason()),
          result.material_report_->lowest_failing_rank());
      return material_failure(result.material_report_->reason(),
                              result.material_report_->lowest_failing_rank());
    }

    if (closure != nullptr) {
      const auto closure_result = closure->evaluate(
          closure->object, state, detail::DensityClosureStage::final);
      if (!closure_result.accepted)
        return fail(StepFailureReason::density_closure_failure,
                    closure_result.lowest_failing_rank,
                    closure_result.recoverable);
      impl_->halo->exchange(trial, fields.density);
      impl_->halo->exchange(trial, closure->enthalpy_density);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      if (closure->after_halo != nullptr)
        closure->after_halo(closure->object,
                            detail::DensityClosureStage::final,
                            fields.density, closure->enthalpy_density);
#endif

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      if (closure->before_post_assessment != nullptr)
        closure->before_post_assessment(closure->object, state);
#endif

      std::optional<MaterialFaceMassFlux> post_closure_flux;
      synchronized_local_phase(
          *impl_->mpi, StepFailureReason::invalid_input, false, [&] {
            post_closure_flux.emplace(
                bind_material_flux(MaterialFluxProvenance::final_corrected));
          });
      result.post_closure_report_.emplace(
          impl_->material_transport.assess_trial_after_closure(
              state, *post_closure_flux, stencil,
              closure->enthalpy_rate_J_per_kg_s));
      result.post_closure_report_->attempt_identity_ = result.attempt_identity_;
      result.post_closure_report_->seal();
      result.post_closure_evidence_available_ = true;
      result.post_closure_authority_ = result.post_closure_report_;
      post_closure_flux.reset();
      result.flow_.final_transport_normalized_l2 =
          result.post_closure_report_->transport_normalized_l2();
      result.flow_.final_transport_relative_conservation_defect =
          result.post_closure_report_->transport_relative_conservation_defect();
      result.flow_.final_mass_relative_conservation_defect =
          result.post_closure_report_->mass_relative_conservation_defect();
      result.mass_conservation_available_ =
          result.post_closure_report_->mass_conservation_available();
      if (result.post_closure_report_->disposition() !=
          MaterialTransportDisposition::finalized) {
        const auto reason = result.post_closure_report_->reason();
        require_reliable_collective_result(
            map_material_failure(reason),
            result.post_closure_report_->lowest_failing_rank());
        return fail(map_material_failure(reason),
                    result.post_closure_report_->lowest_failing_rank(), true);
      }
    }

    if (closure != nullptr && impl_->boundaries->open_domain()) {
      closure->before_outlet(closure->object, state);
      const auto final_flux = trial.acquire_face_read<double>(
          access, kStatePhase, kStateActor, fields.face_mass_flux);
      const auto assessment =
          impl_->boundaries->assess_final_pressure_outlet_flux(
              *impl_->topology, *impl_->mpi, final_flux,
              state.metadata().step + 1U,
              state.metadata().time_s + stencil.dt_s);
      if (assessment.decision ==
          boundary::FinalFluxDecision::outlet_backflow) {
        result.flow_.final_backflow_evidence = assessment.evidence;
        return fail(StepFailureReason::boundary_backflow,
                    assessment.evidence->lowest_failing_rank, true);
      }
    }

    double continuity_sums[2]{};
    synchronized_local_phase(*impl_->mpi,
                             StepFailureReason::final_continuity_residual,
                             true, [&] {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
      FaceFluxPathObservation allocation_observation;
#endif
      const auto rho_final = trial.acquire_read<double>(
          access, kStatePhase, kStateActor, fields.density);
      const auto final_flux = bind_face_flux(trial);
      auto mass_residual = impl_->scratch.storage->acquire_write<double>(
          *impl_->scratch.access, kScratchPhase, kScratchActor,
          impl_->scratch.mass_residual);
      zero(mass_residual);
      impl_->fvm.accumulate_mass_residual(final_flux, mass_residual);
      std::fill(impl_->continuity_absolute.begin(),
                impl_->continuity_absolute.end(), 0.0);
      const auto final_flux_values = trial.acquire_face_read<double>(
          access, kStatePhase, kStateActor, fields.face_mass_flux);
      for (mesh::LocalFaceId face = 0;
           face < impl_->topology->local_face_count(); ++face) {
        const double magnitude = std::abs(final_flux_values(face, 0));
        const auto owner = impl_->topology->owner(face);
        if (impl_->topology->cell_ownership(owner) ==
            mesh::EntityOwnership::owned)
          impl_->continuity_absolute[owner] += magnitude;
        const auto neighbour = impl_->topology->neighbour(face);
        if (!impl_->topology->patch_id(face) && neighbour &&
            impl_->topology->cell_ownership(*neighbour) ==
                mesh::EntityOwnership::owned)
          impl_->continuity_absolute[*neighbour] += magnitude;
      }
      for (mesh::LocalCellId cell = 0; cell < count; ++cell) {
        const auto index = map_cell(impl_->topology->global_cell(cell), owned,
                                    impl_->topology->global_extent());
        const double volume = impl_->geometry->cell_volume_m3(cell);
        const double rho_final_value = at(rho_final, index, 0);
        const double rho_n_value = at(rho_n, index, 0);
        const double rho_nm1_value = at(rho_nm1, index, 0);
        const double temporal =
            (stencil.alpha0 * rho_final_value +
             stencil.alpha1 * rho_n_value +
             stencil.alpha2 * rho_nm1_value) *
            volume / stencil.dt_s;
        const double raw = temporal + at(mass_residual, index, 0);
        continuity_sums[0] += raw * raw;
        const double scale =
            (std::abs(stencil.alpha0 * rho_final_value) +
             std::abs(stencil.alpha1 * rho_n_value) +
             std::abs(stencil.alpha2 * rho_nm1_value)) *
                volume / stencil.dt_s +
            impl_->continuity_absolute[cell];
        continuity_sums[1] += scale * scale;
      }
    });
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    apply_terminal_test_point(
        test::MaterialTerminalPointForTest::final_continuity_reduction);
#endif
    impl_->mpi->allreduce_fp64_in_place(
        continuity_sums, 2U, runtime::Fp64ReductionOperation::sum);
    const double raw_l2 = std::sqrt(continuity_sums[0]);
    const double scale_l2 = std::sqrt(continuity_sums[1]);
    result.flow_.final_continuity_normalized_l2 =
        scale_l2 == 0.0
            ? (raw_l2 == 0.0 ? 0.0
                             : std::numeric_limits<double>::infinity())
            : raw_l2 / scale_l2;
    result.final_continuity_residual_available_ = true;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    apply_terminal_test_point(
        test::MaterialTerminalPointForTest::final_continuity_status);
#endif
    const auto continuity_status = runtime::collective_status(
        *impl_->mpi,
        std::isfinite(result.flow_.final_continuity_normalized_l2) &&
            result.flow_.final_continuity_normalized_l2 <= 1.0e-10,
        "material final continuity residual exceeds tolerance");
    if (!continuity_status.ok)
      return fail(StepFailureReason::final_continuity_residual,
                  continuity_status.failing_rank, true);

    const auto pressure = impl_->coupler.assess_final_material_density_pressure(
        state, stencil, pressure_control);
    result.flow_.final_pressure_residual_l2 =
        pressure.independent_residual_l2;
    result.final_pressure_normalized_residual_ = pressure.normalized_residual;
    result.final_pressure_residual_available_ = pressure.residual_available;
    if (!pressure.accepted) {
      require_reliable_collective_result(pressure.reason,
                                         pressure.lowest_failing_rank);
      return fail(pressure.reason, pressure.lowest_failing_rank,
                  pressure.disposition ==
                      PressureCorrectionDisposition::recoverable_failure);
    }

    impl_->halo->exchange(trial, fields.mechanical_pressure);
    std::array<double, 6> momentum_sums{};
    std::array<double, 36> conservation_sums{};
    synchronized_local_phase(*impl_->mpi,
                             StepFailureReason::final_momentum_residual, true,
                             [&] {
      const auto rho_final = trial.acquire_read<double>(
          access, kStatePhase, kStateActor, fields.density);
      const auto pressure_final = trial.acquire_read<double>(
          access, kStatePhase, kStateActor, fields.mechanical_pressure);
      auto final_gradient = impl_->scratch.storage->acquire_write<double>(
          *impl_->scratch.access, kScratchPhase, kScratchActor,
          impl_->scratch.pressure_gradient);
      compute_pressure_gradient(*impl_->topology, *impl_->geometry,
                                pressure_final, final_gradient,
                                impl_->pressure_gradient_sums);
      impl_->fvm.physical_boundary_pressure_contributions(
          *impl_->boundaries, pressure_final, impl_->pressure_boundary);
      if (impl_->momentum_n.boundary.size() !=
              impl_->pressure_boundary.size() ||
          (stencil.order == MomentumTimeOrder::bdf2 &&
           impl_->momentum_n.boundary.size() !=
               impl_->momentum_nm1.boundary.size()))
        throw runtime::Error(
            "material physical momentum boundary sets differ");
      std::array<CompensatedSum, 3> boundary_n;
      std::array<CompensatedSum, 3> boundary_nm1;
      std::array<CompensatedSum, 3> boundary_abs_n;
      std::array<CompensatedSum, 3> boundary_abs_nm1;
      for (std::size_t boundary_index = 0;
           boundary_index < impl_->momentum_n.boundary.size();
           ++boundary_index) {
        const auto &current = impl_->momentum_n.boundary[boundary_index];
        const auto &pressure_entry = impl_->pressure_boundary[boundary_index];
        const auto *history_entry =
            stencil.order == MomentumTimeOrder::bdf2
                ? &impl_->momentum_nm1.boundary[boundary_index]
                : nullptr;
        if (current.global_face_id != pressure_entry.global_face_id ||
            (history_entry != nullptr &&
             current.global_face_id != history_entry->global_face_id))
          throw runtime::Error(
              "material physical momentum boundary identities differ");
        for (std::size_t direction = 0; direction < 3U; ++direction) {
          const double current_value = current.convective[direction] +
                                       current.viscous[direction] +
                                       pressure_entry.pressure[direction];
          boundary_n[direction].add(current_value);
          boundary_abs_n[direction].add(std::abs(current_value));
          if (history_entry != nullptr) {
            const double history_value =
                history_entry->convective[direction] +
                history_entry->viscous[direction] +
                pressure_entry.pressure[direction];
            boundary_nm1[direction].add(history_value);
            boundary_abs_nm1[direction].add(std::abs(history_value));
          }
        }
      }
      const auto final_velocity = trial.acquire_read<double>(
          access, kStatePhase, kStateActor, fields.velocity);
      std::array<CompensatedSum, 36> local_conservation{};
      for (mesh::LocalCellId cell = 0; cell < count; ++cell) {
      const auto index = map_cell(impl_->topology->global_cell(cell), owned,
                                  impl_->topology->global_extent());
      const double volume = impl_->geometry->cell_volume_m3(cell);
      for (std::size_t direction = 0; direction < 3U; ++direction) {
        const double spatial =
            stencil.order == MomentumTimeOrder::backward_euler
                ? impl_->momentum_n.convection[cell * 3U + direction] +
                      impl_->momentum_n.viscosity[cell * 3U + direction]
                : 2.0 *
                          (impl_->momentum_n.convection[cell * 3U + direction] +
                           impl_->momentum_n.viscosity[cell * 3U + direction]) -
                      (impl_->momentum_nm1
                           .convection[cell * 3U + direction] +
                       impl_->momentum_nm1
                           .viscosity[cell * 3U + direction]);
        double source = 0.0;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
        if (impl_->vortex_source_enabled) {
          const auto centre = impl_->geometry->cell_center_m(cell);
          source = component(analytic_source(centre.x, centre.y, mu),
                             static_cast<int>(direction));
        }
#endif
        const double lhs =
            stencil.alpha0 * at(rho_final, index, 0) * volume /
                stencil.dt_s *
                at(final_velocity, index, static_cast<int>(direction)) +
            (volume / stencil.dt_s) *
                (stencil.alpha1 * at(rho_n, index, 0) *
                     at(velocity_n, index, static_cast<int>(direction)) +
                 stencil.alpha2 * at(rho_nm1, index, 0) *
                     at(velocity_nm1, index, static_cast<int>(direction))) +
            spatial +
            volume * at(final_gradient, index, static_cast<int>(direction)) -
            volume * source;
        const double scale =
            std::abs(stencil.alpha0 * at(rho_final, index, 0) * volume /
                     stencil.dt_s *
                     at(final_velocity, index, static_cast<int>(direction)));
        momentum_sums[direction * 2U] += lhs * lhs;
        momentum_sums[direction * 2U + 1U] += scale * scale;
        const double momentum_n =
            volume * at(rho_n, index, 0) *
            at(velocity_n, index, static_cast<int>(direction));
        const double momentum_nm1 =
            volume * at(rho_nm1, index, 0) *
            at(velocity_nm1, index, static_cast<int>(direction));
        const double momentum_next =
            volume * at(rho_final, index, 0) *
            at(final_velocity, index, static_cast<int>(direction));
        const double source_integral = volume * source;
        const std::size_t offset = direction * 12U;
        local_conservation[offset].add(momentum_nm1);
        local_conservation[offset + 1U].add(momentum_n);
        local_conservation[offset + 2U].add(momentum_next);
        local_conservation[offset + 5U].add(source_integral);
        local_conservation[offset + 6U].add(std::abs(momentum_nm1));
        local_conservation[offset + 7U].add(std::abs(momentum_n));
        local_conservation[offset + 8U].add(std::abs(momentum_next));
        local_conservation[offset + 11U].add(std::abs(source_integral));
        }
      }
      for (std::size_t direction = 0; direction < 3U; ++direction) {
        const std::size_t offset = direction * 12U;
        local_conservation[offset + 3U].add(boundary_nm1[direction].value());
        local_conservation[offset + 4U].add(boundary_n[direction].value());
        local_conservation[offset + 9U].add(
            boundary_abs_nm1[direction].value());
        local_conservation[offset + 10U].add(
            boundary_abs_n[direction].value());
      }
      std::transform(local_conservation.begin(), local_conservation.end(),
                     conservation_sums.begin(),
                     [](const CompensatedSum &sum) { return sum.value(); });
    });
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    apply_terminal_test_point(
        test::MaterialTerminalPointForTest::final_momentum_residual_reduction);
#endif
    impl_->mpi->allreduce_fp64_in_place(
        momentum_sums.data(), momentum_sums.size(),
        runtime::Fp64ReductionOperation::sum);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    apply_terminal_test_point(test::MaterialTerminalPointForTest::
                                  final_momentum_conservation_reduction);
#endif
    impl_->mpi->allreduce_fp64_in_place(
        conservation_sums.data(), conservation_sums.size(),
        runtime::Fp64ReductionOperation::sum);
    for (std::size_t direction = 0; direction < 3U; ++direction) {
      result.flow_.final_momentum_normalized_l2[direction] =
          std::sqrt(momentum_sums[direction * 2U] /
                    std::max(momentum_sums[direction * 2U + 1U],
                             std::numeric_limits<double>::min()));
      result.final_momentum_residual_available_[direction] = 1U;
      const std::size_t offset = direction * 12U;
      result.flow_.final_momentum_relative_conservation_defect[direction] =
          momentum_conservation_defect(
              {conservation_sums[offset],
               conservation_sums[offset + 1U],
               conservation_sums[offset + 2U],
               conservation_sums[offset + 3U],
               conservation_sums[offset + 4U],
               conservation_sums[offset + 5U],
               conservation_sums[offset + 6U],
               conservation_sums[offset + 7U],
               conservation_sums[offset + 8U],
               conservation_sums[offset + 9U],
               conservation_sums[offset + 10U],
               conservation_sums[offset + 11U],
               stencil.dt_s,
               stencil.alpha0,
               stencil.alpha2,
               stencil.order == MomentumTimeOrder::bdf2});
      result.momentum_conservation_available_[direction] = 1U;
    }
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    apply_terminal_test_point(
        test::MaterialTerminalPointForTest::final_momentum_status);
#endif
    const auto momentum_status = runtime::collective_status(
        *impl_->mpi,
        std::all_of(result.flow_.final_momentum_normalized_l2.begin(),
                    result.flow_.final_momentum_normalized_l2.end(),
                    [](double value) {
                      return std::isfinite(value) && value <= 1.0e-9;
                    }),
        "material final momentum residual exceeds tolerance");
    if (!momentum_status.ok)
      return fail(StepFailureReason::final_momentum_residual,
                  momentum_status.failing_rank, true);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    apply_terminal_test_point(
        test::MaterialTerminalPointForTest::final_conservation_status);
#endif
    const auto conservation_status = runtime::collective_status(
        *impl_->mpi,
        std::all_of(
            result.flow_.final_momentum_relative_conservation_defect.begin(),
            result.flow_.final_momentum_relative_conservation_defect.end(),
            [](double value) {
              return std::isfinite(value) && value <= 5.0e-11;
            }),
        "material final momentum conservation exceeds tolerance");
    if (!conservation_status.ok)
      return fail(StepFailureReason::final_conservation_defect,
                  conservation_status.failing_rank, true);

    const auto old = state.metadata();
    const AcceptedStepMetadata accepted{old.step + 1U,
                                        old.time_s + stencil.dt_s,
                                        stencil.dt_s,
                                        old.dt_s,
                                        stencil.order};
    if (closure != nullptr) {
      const int preparation_rank =
          closure->before_prepare(closure->object, state, accepted);
      if (preparation_rank >= 0)
        return fail(StepFailureReason::invalid_input, preparation_rank, false);
    }
    state.prepare_commit_attempt(accepted);
    if (closure != nullptr) {
      const int preparation_rank = closure->prepare(closure->object);
      if (preparation_rank >= 0)
        return fail(StepFailureReason::invalid_input, preparation_rank, false);
    }
    state.publish_commit_attempt();
    if (closure != nullptr)
      closure->publish(closure->object);
    active = false;
    result.flow_.disposition = StepAttemptDisposition::committed;
    result.flow_.reason = StepFailureReason::none;
    result.flow_.lowest_failing_rank = -1;
    result.flow_.suggested_dt_s = 0.0;
    return finish();
  } catch (const SynchronizedMaterialFailure &failure) {
    return fail(failure.reason, failure.failing_rank, failure.recoverable);
  } catch (const detail::DensityClosurePreflightFailure &failure) {
    return fail(StepFailureReason::invalid_input, failure.failing_rank(),
                false);
  } catch (const runtime::MpiOperationError &) {
    if (active)
      state.rollback_attempt();
    if (closure != nullptr)
      closure->rollback(closure->object);
    impl_->last_state = nullptr;
    impl_->last_report_seal = 0U;
    throw;
  } catch (const runtime::Error &) {
    const auto selected = synchronize_local_failure(
        *impl_->mpi, {true, StepFailureReason::invalid_input, false});
    return fail(selected.reason, selected.failing_rank, selected.recoverable);
  } catch (...) {
    const auto selected = synchronize_local_failure(
        *impl_->mpi, {true, StepFailureReason::invalid_input, false});
    return fail(selected.reason, selected.failing_rank, selected.recoverable);
  }
}

bool detail::DensityClosureBridge::report_authenticated(
    const MaterialDensityStepAttemptReport &report) noexcept {
  return report.authenticated();
}

std::uint64_t detail::DensityClosureBridge::report_seal(
    const MaterialDensityStepAttemptReport &report) noexcept {
  return report.seal_;
}

MaterialDensityStepAttemptReport detail::DensityClosureBridge::make_report() {
  return MaterialDensityStepAttemptReport{};
}

bool detail::DensityClosureBridge::post_eos_evidence_authenticated(
    const MaterialDensityStepAttemptReport &report) noexcept {
  return report.authenticated() && report.post_closure_evidence_available_ &&
         report.post_closure_report_.has_value();
}

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
void detail::DensityClosureBridge::force_finalization_identity_wrap(
    FixedStepMaterialDensityFlow &flow) {
  if (!flow.impl_)
    throw runtime::Error("material flow object has been moved from");
  test::MaterialDensityTransportTestAccess::force_finalization_identity_wrap(
      flow.impl_->material_transport);
}
#endif


MaterialDensityFlowDiagnosticSource
FixedStepMaterialDensityFlow::diagnostic_source(
    const FlowState &state,
    const MaterialDensityStepAttemptReport &report) const {
  if (!impl_)
    throw runtime::Error("material flow object has been moved from");
  const std::uint64_t state_identity = state.diagnostic_mutation_identity();
  if (state_identity == 0U || state.attempt_active() ||
      impl_->last_state != &state || !report.authenticated() ||
      report.attempt_identity_ != impl_->attempt_identity ||
      impl_->last_report_seal != report.seal_ ||
      impl_->last_state_identity != state_identity)
    throw runtime::Error("material flow diagnostic source is stale");
  auto source = std::make_unique<MaterialDensityFlowDiagnosticSource::Impl>();
  source->flow = impl_.get();
  source->state = &state;
  source->report = report;
  source->flow_identity = impl_->attempt_identity;
  source->source_generation = impl_->source_generation;
  source->state_identity = state_identity;
  source->report_seal = report.seal_;
  source->rank = impl_->mpi->rank();
  source->step = state.metadata().step;
  source->time_s = state.metadata().time_s;
  return MaterialDensityFlowDiagnosticSource(std::move(source));
}

namespace {

constexpr std::array<std::string_view, 5> kFingerprintIds{
    "face_mass_flux", "face_velocity", "pi", "rho", "velocity"};
constexpr std::array<std::string_view, 5> kFingerprintUnits{
    "kg/s", "m/s", "Pa", "kg/m3", "m/s"};
constexpr std::array<std::size_t, 5> kFingerprintComponents{1U, 3U, 1U, 1U,
                                                           3U};

} // namespace

void MaterialDensityFlowDiagnosticSource::validate() const {
  if (!impl_ || impl_->flow == nullptr || impl_->state == nullptr ||
      impl_->state->diagnostic_mutation_identity() != impl_->state_identity ||
      impl_->state->attempt_active() ||
      impl_->flow->attempt_identity != impl_->flow_identity ||
      impl_->flow->source_generation != impl_->source_generation ||
      impl_->flow->last_state != impl_->state ||
      impl_->flow->last_report_seal != impl_->report_seal ||
      impl_->report.attempt_identity() != impl_->flow_identity ||
      !impl_->report.authenticated())
    throw runtime::Error("material flow diagnostic source is stale");
}

std::size_t MaterialDensityFlowDiagnosticSource::fingerprint_field_count()
    const {
  validate();
  return kFingerprintIds.size();
}

std::string_view MaterialDensityFlowDiagnosticSource::fingerprint_field_id(
    std::size_t field) const {
  validate();
  if (field >= kFingerprintIds.size())
    throw runtime::Error("material flow diagnostic field index is invalid");
  return kFingerprintIds[field];
}

std::string_view
MaterialDensityFlowDiagnosticSource::field_unit(std::size_t field) const {
  validate();
  if (field >= kFingerprintUnits.size())
    throw runtime::Error("material flow diagnostic field index is invalid");
  return kFingerprintUnits[field];
}

MaterialDensityDiagnosticEntity
MaterialDensityFlowDiagnosticSource::field_entity(std::size_t field) const {
  validate();
  if (field >= kFingerprintIds.size())
    throw runtime::Error("material flow diagnostic field index is invalid");
  return field < 2U ? MaterialDensityDiagnosticEntity::face
                    : MaterialDensityDiagnosticEntity::cell;
}

std::size_t MaterialDensityFlowDiagnosticSource::field_component_count(
    std::size_t field) const {
  validate();
  if (field >= kFingerprintComponents.size())
    throw runtime::Error("material flow diagnostic field index is invalid");
  return kFingerprintComponents[field];
}

std::size_t MaterialDensityFlowDiagnosticSource::field_item_count(
    std::size_t field) const {
  validate();
  if (field >= kFingerprintIds.size())
    throw runtime::Error("material flow diagnostic field index is invalid");
  return field < 2U ? impl_->flow->canonical_faces.size()
                    : impl_->flow->topology->owned_cell_count();
}

std::uint64_t MaterialDensityFlowDiagnosticSource::field_global_id(
    std::size_t field, std::size_t item) const {
  validate();
  if (field >= kFingerprintIds.size())
    throw runtime::Error("material flow diagnostic field index is invalid");
  const std::size_t count = field < 2U
                                ? impl_->flow->canonical_faces.size()
                                : impl_->flow->topology->owned_cell_count();
  if (item >= count)
    throw runtime::Error("material flow diagnostic item index is invalid");
  return field < 2U
             ? impl_->flow->topology->global_face_id(
                   impl_->flow->canonical_faces[item])
             : impl_->flow->topology->global_cell_id(item);
}

double MaterialDensityFlowDiagnosticSource::field_value(
    std::size_t field, std::size_t item, std::size_t field_component) const {
  validate();
  if (field >= kFingerprintIds.size() ||
      field_component >=
          (field < kFingerprintComponents.size()
               ? kFingerprintComponents[field]
               : 0U))
    throw runtime::Error("material flow diagnostic field component is invalid");
  const std::size_t count = field < 2U
                                ? impl_->flow->canonical_faces.size()
                                : impl_->flow->topology->owned_cell_count();
  if (item >= count)
    throw runtime::Error("material flow diagnostic item index is invalid");
  const auto &access = detail::FlowStateSolverAccess::access(*impl_->state);
  const auto &storage = impl_->state->layer(FlowLayer::committed);
  runtime::FieldId field_id{};
  if (field == 0U)
    field_id = impl_->flow->fields.face_mass_flux;
  else if (field == 1U)
    field_id = impl_->flow->fields.face_velocity;
  else if (field == 2U)
    field_id = impl_->flow->fields.mechanical_pressure;
  else if (field == 3U)
    field_id = impl_->flow->fields.density;
  else
    field_id = impl_->flow->fields.velocity;
  if (field < 2U) {
    const auto values = storage.acquire_face_read<double>(
        access, kStatePhase, kStateActor, field_id);
    return values(impl_->flow->canonical_faces[item],
                  static_cast<int>(field_component));
  }
  const auto values = storage.acquire_read<double>(
      access, kStatePhase, kStateActor, field_id);
  const auto index = map_cell(impl_->flow->topology->global_cell(item),
                              impl_->flow->topology->owned_global_box(),
                              impl_->flow->topology->global_extent());
  return values(index.i, index.j, index.k,
                static_cast<int>(field_component));
}

std::string_view
MaterialDensityFlowDiagnosticSource::owned_cell_layout_fingerprint() const {
  validate();
  return impl_->flow->owned_cell_fingerprint;
}

std::string_view
MaterialDensityFlowDiagnosticSource::global_cell_layout_fingerprint() const {
  validate();
  return impl_->flow->global_cell_fingerprint;
}

std::string_view
MaterialDensityFlowDiagnosticSource::owned_face_layout_fingerprint() const {
  validate();
  return impl_->flow->owned_face_fingerprint;
}

std::string_view
MaterialDensityFlowDiagnosticSource::global_face_layout_fingerprint() const {
  validate();
  return impl_->flow->global_face_fingerprint;
}

int MaterialDensityFlowDiagnosticSource::relative_rank() const {
  validate();
  return impl_->rank;
}

std::uint64_t MaterialDensityFlowDiagnosticSource::committed_step() const {
  validate();
  return impl_->step;
}

double MaterialDensityFlowDiagnosticSource::committed_time_s() const {
  validate();
  return impl_->time_s;
}

runtime::Int3 MaterialDensityFlowDiagnosticSource::global_cell_extent() const {
  validate();
  return impl_->flow->topology->global_extent();
}

runtime::Box3 MaterialDensityFlowDiagnosticSource::owned_global_box() const {
  validate();
  return impl_->flow->topology->owned_global_box();
}

std::size_t MaterialDensityFlowDiagnosticSource::owned_cell_count() const {
  validate();
  return impl_->flow->topology->owned_cell_count();
}

std::size_t
MaterialDensityFlowDiagnosticSource::canonical_owned_face_count() const {
  validate();
  return impl_->flow->canonical_faces.size();
}

double MaterialDensityFlowDiagnosticSource::cell_volume_m3(
    std::size_t owned_cell) const {
  validate();
  if (owned_cell >= impl_->flow->topology->owned_cell_count())
    throw runtime::Error("material flow diagnostic cell index is invalid");
  return impl_->flow->geometry->cell_volume_m3(owned_cell);
}

const MaterialDensityStepAttemptReport &
MaterialDensityFlowDiagnosticSource::report() const {
  validate();
  return impl_->report;
}

MaterialDensityFlowDiagnosticSource::MaterialDensityFlowDiagnosticSource(
    std::unique_ptr<Impl> implementation) noexcept
    : impl_(std::move(implementation)) {}

MaterialDensityFlowDiagnosticSource::~MaterialDensityFlowDiagnosticSource()
    noexcept = default;

MaterialDensityFlowDiagnosticSource::MaterialDensityFlowDiagnosticSource(
    MaterialDensityFlowDiagnosticSource &&) noexcept = default;

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
double test::MaterialDensityPisoTestAccess::momentum_conservation_defect(
    const MaterialMomentumConservationInput &input) noexcept {
  return ::hundun::flow::momentum_conservation_defect(
      {input.momentum_n_minus_1,
       input.momentum_n,
       input.momentum_n_plus_1,
       input.boundary_n_minus_1,
       input.boundary_n,
       input.source_n_plus_1,
       input.momentum_abs_n_minus_1,
       input.momentum_abs_n,
       input.momentum_abs_n_plus_1,
       input.boundary_abs_n_minus_1,
       input.boundary_abs_n,
       input.source_abs_n_plus_1,
       input.dt_s,
       input.alpha0,
       input.alpha2,
       input.bdf2});
}

bool test::MaterialDensityPisoTestAccess::report_authenticated(
    const MaterialDensityStepAttemptReport &report) noexcept {
  return report.authenticated();
}

void test::MaterialDensityPisoTestAccess::corrupt_report(
    MaterialDensityStepAttemptReport &report,
    MaterialReportCorruptionForTest corruption) {
  switch (corruption) {
  case MaterialReportCorruptionForTest::success_corrector_count:
    ++report.flow_.pressure_corrector_count;
    break;
  case MaterialReportCorruptionForTest::success_provenance:
    report.flux_provenance_ = MaterialFluxProvenance::predictor;
    break;
  case MaterialReportCorruptionForTest::success_shared_field:
    ++report.shared_face_mass_flux_field_;
    break;
  case MaterialReportCorruptionForTest::reliable_collective_rank:
    report.flow_.disposition = StepAttemptDisposition::non_retryable_failure;
    report.flow_.reason = StepFailureReason::collective_operation;
    report.flow_.lowest_failing_rank = -1;
    break;
  case MaterialReportCorruptionForTest::material_reason_mapping:
    report.material_failure_reason_ =
        report.material_failure_reason_ == MaterialTransportFailureReason::none
            ? MaterialTransportFailureReason::non_positive_density
            : MaterialTransportFailureReason::invalid_input;
    break;
  case MaterialReportCorruptionForTest::parent_transport_size:
    report.flow_.final_transport_normalized_l2.push_back(0.0);
    break;
  case MaterialReportCorruptionForTest::material_transport_size:
    ++report.material_field_count_;
    break;
  case MaterialReportCorruptionForTest::unavailable_numeric_value:
    report.final_continuity_residual_available_ = false;
    report.flow_.final_continuity_normalized_l2 = 1.0;
    break;
  case MaterialReportCorruptionForTest::material_count_zero:
    report.material_field_count_ = 0U;
    break;
  case MaterialReportCorruptionForTest::material_count_plus_two:
    report.material_field_count_ += 2U;
    break;
  case MaterialReportCorruptionForTest::material_count_five:
    report.material_field_count_ = 5U;
    break;
  case MaterialReportCorruptionForTest::parent_transport_residual_value:
    report.flow_.final_transport_normalized_l2.front() =
        std::nextafter(report.flow_.final_transport_normalized_l2.front(),
                       std::numeric_limits<double>::infinity());
    break;
  case MaterialReportCorruptionForTest::nested_transport_residual_value:
    report.material_report_->transport_normalized_l2_.front() = std::nextafter(
        report.material_report_->transport_normalized_l2_.front(),
        std::numeric_limits<double>::infinity());
    break;
  case MaterialReportCorruptionForTest::parent_transport_conservation_value:
    report.flow_.final_transport_relative_conservation_defect.front() =
        std::nextafter(
            report.flow_.final_transport_relative_conservation_defect.front(),
            std::numeric_limits<double>::infinity());
    break;
  case MaterialReportCorruptionForTest::nested_transport_conservation_value:
    report.material_report_->transport_relative_conservation_defect_.front() =
        std::nextafter(report.material_report_
                           ->transport_relative_conservation_defect_.front(),
                       std::numeric_limits<double>::infinity());
    break;
  case MaterialReportCorruptionForTest::nested_density_residual_availability:
    report.material_report_->density_residual_available_ = false;
    break;
  case MaterialReportCorruptionForTest::nested_transport_residual_availability:
    report.material_report_->transport_residual_available_.front() = 0U;
    break;
  case MaterialReportCorruptionForTest::nested_mass_conservation_availability:
    report.material_report_->mass_conservation_available_ = false;
    break;
  case MaterialReportCorruptionForTest::
      nested_transport_conservation_availability:
    report.material_report_->transport_conservation_available_.front() = 0U;
    break;
  case MaterialReportCorruptionForTest::nested_minimum_density_availability:
    report.material_report_->minimum_density_available_ = false;
    break;
  case MaterialReportCorruptionForTest::nested_attempt_identity:
    ++report.material_report_->attempt_identity_;
    break;
  case MaterialReportCorruptionForTest::nested_finalization_identity:
    report.material_report_->finalization_identity_ = 0U;
    break;
  case MaterialReportCorruptionForTest::nested_shared_field:
    ++report.material_report_->shared_face_mass_flux_field_;
    break;
  case MaterialReportCorruptionForTest::nested_provenance:
    report.material_report_->flux_provenance_ =
        MaterialFluxProvenance::predictor;
    break;
  case MaterialReportCorruptionForTest::nested_residual_outer_size:
    report.material_report_->transport_residual_available_.pop_back();
    break;
  case MaterialReportCorruptionForTest::nested_conservation_outer_size:
    report.material_report_->transport_conservation_available_.pop_back();
    break;
  }
  if (report.material_report_)
    report.material_report_->seal();
  report.seal_ = report.compute_seal();
}

void test::MaterialDensityPisoTestAccess::set_preflight_allocation_failure_rank(
    int rank) noexcept {
  ::hundun::flow::preflight_allocation_failure_rank.store(
      rank, std::memory_order_relaxed);
}

void test::MaterialDensityPisoTestAccess::reset_preflight_allocation_failure()
    noexcept {
  ::hundun::flow::preflight_allocation_failure_rank.store(
      -1, std::memory_order_relaxed);
}

bool test::MaterialDensityPisoTestAccess::face_flux_path_observation_active()
    noexcept {
  return test::detail::face_flux_path_observation_depth.load(
             std::memory_order_relaxed) != 0U;
}

test::MaterialPhaseSelectionForTest
test::MaterialDensityPisoTestAccess::select_phase_failure(
    const runtime::MpiContext &mpi, MaterialPhaseFailureForTest input) {
  const auto selected = ::hundun::flow::synchronize_local_failure(
      mpi, {input.failed, static_cast<StepFailureReason>(input.reason),
            input.recoverable});
  return {static_cast<std::uint8_t>(selected.reason), selected.failing_rank,
          selected.recoverable};
}

void test::MaterialDensityPisoTestAccess::require_reliable_collective_result(
    std::uint8_t reason, int lowest_failing_rank) {
  ::hundun::flow::require_reliable_collective_result(
      static_cast<StepFailureReason>(reason), lowest_failing_rank);
}

test::MaterialPressureEvidenceForTest
test::MaterialDensityPisoTestAccess::material_pressure_evidence(
    const FixedStepMaterialDensityFlow &flow) {
  if (!flow.impl_)
    return {};
  return material_pressure_evidence(flow.impl_->coupler);
}

const std::vector<double> &
test::MaterialDensityPisoTestAccess::finalizer_flux_evidence(
    const FixedStepMaterialDensityFlow &flow) {
  if (!flow.impl_)
    throw runtime::Error("material flow object has been moved from");
  return flow.impl_->finalizer_flux_input;
}

void test::MaterialDensityPisoTestAccess::force_flow_attempt_identity(
    FixedStepMaterialDensityFlow &flow, std::uint64_t value) noexcept {
  if (flow.impl_)
    flow.impl_->attempt_identity = value;
}

bool test::MaterialDensityPisoTestAccess::has_diagnostic_report(
    const FixedStepMaterialDensityFlow &flow) noexcept {
  return flow.impl_ != nullptr && flow.impl_->last_state != nullptr &&
         flow.impl_->last_report_seal != 0U;
}

void test::MaterialDensityPisoTestAccess::enable_vortex_source(
    FixedStepMaterialDensityFlow &flow, bool enabled) noexcept {
  if (flow.impl_)
    flow.impl_->vortex_source_enabled = enabled;
}
#endif

} // namespace hundun::flow

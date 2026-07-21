// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow/constant_density_piso.hpp"

#include "hundun/finite_volume/matrix_free_poisson.hpp"
#include "hundun/finite_volume/poisson_boundary_adapter.hpp"
#include "hundun/linear/ghosted_vector.hpp"
#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"
#include "mpi_error.hpp"
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#include "constant_density_piso_test_access.hpp"
#endif

#include <algorithm>
#include <array>
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#include <atomic>
#endif
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <new>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace hundun::flow {

struct detail::FlowStateSolverAccess final {
  static const runtime::FieldRegistry &registry(const FlowState &state) {
    return state.solver_registry();
  }
  static const runtime::FieldAccessPlan &access(const FlowState &state) {
    return state.solver_access_plan();
  }
  static runtime::FieldStorage &layer(FlowState &state, FlowLayer selected) {
    return state.solver_layer(selected);
  }
};

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

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
std::atomic<bool> force_final_continuity_failure{false};
std::atomic<bool> force_final_pressure_failure{false};
std::atomic<bool> force_local_derived_failure{false};
std::atomic<bool> provisional_transport_sentinel{false};
std::atomic<double> final_uniform_x_mass_flux{0.0};
std::array<std::atomic<double>, 3> last_momentum_rhs{};
std::array<std::atomic<double>, 3> last_momentum_diagonal{};
std::atomic<std::size_t> provisional_transport_call_count{0U};
std::atomic<std::size_t> final_transport_call_count{0U};
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
    throw runtime::Error("Task 18 vector byte count overflows");
  }
  return count * sizeof(double);
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
    throw runtime::Error("Task 18 cell has no structured field mapping");
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
  runtime::FieldId velocity_gradient{};
  runtime::FieldId pressure_gradient{};
  runtime::FieldId momentum_face{};
  runtime::FieldId momentum_residual{};
  runtime::FieldId actual_diagonal{};
  runtime::FieldId scalar_gradient{};
  runtime::FieldId scalar_face{};
  runtime::FieldId scalar_gamma{};
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

  void replace(std::vector<double> diagonal) {
    if (diagonal.size() != diagonal_.size() ||
        !std::all_of(diagonal.begin(), diagonal.end(), [](double value) {
          return value > 0.0 && std::isfinite(value);
        })) {
      throw runtime::Error("Task 18 momentum diagonal is invalid");
    }
    if (revision_ == std::numeric_limits<std::uint64_t>::max()) {
      throw runtime::Error("Task 18 momentum revision would wrap");
    }
    diagonal_ = std::move(diagonal);
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
      throw runtime::Error("Task 18 momentum operator view is incompatible");
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
      throw runtime::Error("Task 18 momentum diagonal view is incompatible");
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
    throw runtime::Error("Task 18 requires a host execution context");
  }
}

double face_factor(const mesh::MeshGeometry &geometry, LocalFaceId face) {
  const Real3 area = geometry.face_area_vector_m2(face, FaceSide::owner);
  const Real3 displacement = geometry.face_displacement_m(face);
  const double projection = dot(area, displacement);
  const double factor = dot(area, area) / projection;
  if (!finite(area) || !finite(displacement) || !(projection > 0.0) ||
      !(factor > 0.0) || !std::isfinite(factor)) {
    throw runtime::Error("Task 18 pressure face metric is invalid");
  }
  return factor;
}

template <class PressureBoundaryValue>
void compute_pressure_gradient(const mesh::MeshTopology &topology,
                               const mesh::MeshGeometry &geometry,
                               const runtime::FieldView<const double> &pressure,
                               const runtime::FieldView<double> &gradient,
                               PressureBoundaryValue &&boundary_value) {
  std::vector<Real3> sums(topology.owned_cell_count());
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
      throw runtime::Error("Task 18 pressure face value is non-finite");
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
      throw runtime::Error("Task 18 pressure gradient is non-finite");
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

double low_u32(std::uint64_t value) noexcept {
  return static_cast<double>(value & UINT64_C(0xffffffff));
}

double high_u32(std::uint64_t value) noexcept {
  return static_cast<double>(value >> 32U);
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
                           bytes_for(supplied_topology.owned_cell_count())) {}

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
  ScratchFields scratch;
  execution::Buffer gamma;
  execution::Buffer rhs;
  execution::Buffer correction;
  execution::Buffer residual;
  execution::Buffer candidate_pressure;
  std::optional<finite_volume::MatrixFreePoissonOperator> pressure_operator;
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
  if (geometry.mapping_kind() != mesh::MappingKind::uniform_box) {
    throw runtime::Error("Task 18 PISO supports uniform geometry only");
  }
  if (!cell_halo.is_compatible_with(decomposition) ||
      cell_halo.ghost_width() != 2) {
    throw runtime::Error("Task 18 PISO requires a compatible width-two Halo");
  }
  return PisoCoupler(std::make_unique<Impl>(
      decomposition, topology, geometry, boundaries, mpi, execution_context,
      cell_halo, pressure_solver, pressure_preconditioner));
}

PisoCoupler::PisoCoupler(std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
PisoCoupler::~PisoCoupler() noexcept = default;
PisoCoupler::PisoCoupler(PisoCoupler &&) noexcept = default;

PressureCorrectionReport PisoCoupler::correct(
    FlowState &state, double rho_ref,
    const runtime::FieldView<const double> &actual_momentum_diagonal,
    const linear::SolveControl &control) const {
  const std::size_t count = impl_->topology->owned_cell_count();
  const Int3 local_extent = impl_->decomposition->local_extent();
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::pressure_linear_solve, true,
      "Task 18 pressure-correction input validation failed", [&] {
        if (!state.attempt_active()) {
          throw runtime::Error("pressure correction requires an active trial");
        }
        if (!(rho_ref > 0.0) || !std::isfinite(rho_ref)) {
          throw runtime::Error("pressure correction density is invalid");
        }
        if (actual_momentum_diagonal.interior_extent().x != local_extent.x ||
            actual_momentum_diagonal.interior_extent().y != local_extent.y ||
            actual_momentum_diagonal.interior_extent().z != local_extent.z ||
            actual_momentum_diagonal.components() != 3U ||
            actual_momentum_diagonal.ghost_width() < 2) {
          throw runtime::Error(
              "pressure correction diagonal layout is invalid");
        }
      });

  auto &trial = state.solver_layer(FlowLayer::trial);
  const auto &access = state.solver_access_plan();
  const auto &registry = state.solver_registry();
  const auto fields = state.fields();
  const auto owned = impl_->topology->owned_global_box();
  const Int3 global_extent = impl_->topology->global_extent();
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::pressure_linear_solve, true,
      "Task 18 pressure-system derivation failed", [&] {
        auto mass = finite_volume::FaceMassFlux::acquire(
            registry, trial, access, kStatePhase, kStateActor,
            fields.face_mass_flux, *impl_->topology);
        auto mass_residual = impl_->scratch.storage->acquire_write<double>(
            *impl_->scratch.access, kScratchPhase, kScratchActor,
            impl_->scratch.mass_residual);
        zero_cell(mass_residual);
        impl_->fvm.accumulate_mass_residual(mass, mass_residual);

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
          gamma[face] = rho_ref / lambda;
        }
      });

  finite_volume::PoissonBoundarySpec boundary_spec;
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::pressure_linear_solve, true,
      "Task 18 pressure-boundary derivation failed", [&] {
        boundary_spec =
            finite_volume::make_poisson_boundary_spec(*impl_->boundaries);
      });
  if (!impl_->constraint.has_value()) {
    std::optional<finite_volume::PoissonConstraint> candidate_constraint;
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::pressure_linear_solve, true,
        "Task 18 pressure-constraint construction failed", [&] {
          candidate_constraint.emplace(finite_volume::PoissonConstraint::create(
              *impl_->topology, *impl_->geometry, *impl_->execution,
              *impl_->mpi, boundary_spec.mode));
        });
    impl_->constraint.emplace(std::move(*candidate_constraint));
  }
  if (!impl_->pressure_operator.has_value()) {
    impl_->pressure_operator.emplace(
        finite_volume::MatrixFreePoissonOperator::create(
            *impl_->decomposition, *impl_->topology, *impl_->geometry,
            *impl_->execution,
            static_cast<const execution::Buffer &>(impl_->gamma)
                .view(0U, impl_->topology->local_face_count()),
            boundary_spec));
  } else {
    try {
      synchronized_local_phase(
          *impl_->mpi, StepFailureReason::pressure_linear_solve, true,
          "Task 18 pressure-operator update failed", [&] {
            impl_->pressure_operator->replace_face_coefficients(
                static_cast<const execution::Buffer &>(impl_->gamma)
                    .view(0U, impl_->topology->local_face_count()));
          });
    } catch (const SynchronizedAttemptFailure &) {
      impl_->pressure_operator.reset();
      throw;
    }
  }

  std::optional<execution::VectorView<double>> rhs;
  std::optional<execution::VectorView<double>> correction;
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::pressure_linear_solve, true,
      "Task 18 pressure right-hand-side assembly failed", [&] {
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
          (*correction)[cell] = 0.0;
          if (!std::isfinite((*rhs)[cell])) {
            throw runtime::Error("Task 18 pressure RHS is non-finite");
          }
        }
      });
  impl_->constraint->project_rhs(*rhs);
  impl_->constraint->normalize_solution(*correction);
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::pressure_linear_solve, true,
      "Task 18 pressure preconditioner update failed", [&] {
        impl_->preconditioner->update(*impl_->pressure_operator,
                                      impl_->pressure_operator->revision());
      });
  PressureCorrectionReport result;
  result.solve =
      impl_->solver->solve(*impl_->pressure_operator, *impl_->preconditioner,
                           *rhs, *correction, control);
  impl_->constraint->normalize_solution(*correction);
  std::optional<execution::VectorView<double>> residual;
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::pressure_linear_solve, true,
      "Task 18 pressure residual view acquisition failed",
      [&] { residual.emplace(impl_->residual.view(0U, count)); });
  impl_->pressure_operator->apply(*correction, *residual).wait();
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::pressure_linear_solve, true,
      "Task 18 pressure residual assembly failed", [&] {
        for (std::size_t cell = 0; cell < count; ++cell) {
          (*residual)[cell] = (*rhs)[cell] - (*residual)[cell];
          if (!std::isfinite((*residual)[cell])) {
            throw runtime::Error("Task 18 pressure residual is non-finite");
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
      "Task 18 pressure correction did not satisfy the residual contract");
  if (!acceptance.ok) {
    result.accepted = false;
    if (result.solve.lowest_failing_rank < 0) {
      result.solve.lowest_failing_rank = acceptance.failing_rank;
    }
    return result;
  }

  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::non_finite_trial, true,
      "Task 18 pressure-correction staging failed", [&] {
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
      "Task 18 pressure-correction gradient failed", [&] {
        pi_prime_read.emplace(impl_->scratch.storage->acquire_read<double>(
            *impl_->scratch.access, kScratchPhase, kScratchActor,
            impl_->scratch.pressure_correction));
        auto prime_gradient = impl_->scratch.storage->acquire_write<double>(
            *impl_->scratch.access, kScratchPhase, kScratchActor,
            impl_->scratch.pressure_gradient);
        compute_pressure_gradient(
            *impl_->topology, *impl_->geometry, *pi_prime_read, prime_gradient,
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

  std::vector<double> velocity_candidate;
  std::vector<double> pressure_candidate;
  std::vector<double> face_velocity_candidate;
  std::vector<double> mass_flux_candidate;
  std::optional<execution::VectorView<double>> candidate_pressure;
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::non_finite_trial, true,
      "Task 18 corrected cell candidate derivation failed", [&] {
        velocity_candidate.resize(count * 3U);
        pressure_candidate.resize(count);
        face_velocity_candidate.resize(impl_->topology->local_face_count() *
                                       3U);
        mass_flux_candidate.resize(impl_->topology->local_face_count());
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
      "Task 18 corrected face candidate derivation failed", [&] {
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
          const double lambda_inverse = gamma_read[face] / rho_ref;
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
  const double flux_roundoff = 256.0 * std::numeric_limits<double>::epsilon() *
                               rho_ref * std::max(1.0, roundoff_scales[1]);
  if (roundoff_scales[0] <= flux_roundoff &&
      roundoff_scales[2] <= 256.0 * std::numeric_limits<double>::epsilon()) {
    std::fill(mass_flux_candidate.begin(), mass_flux_candidate.end(), 0.0);
    std::fill(face_velocity_candidate.begin(), face_velocity_candidate.end(),
              0.0);
  }
  synchronized_local_phase(
      *impl_->mpi, StepFailureReason::non_finite_trial, true,
      "Task 18 corrected trial write failed", [&] {
        const auto all_finite = [](const std::vector<double> &values) {
          return std::all_of(values.begin(), values.end(),
                             [](double value) { return std::isfinite(value); });
        };
        if (!all_finite(velocity_candidate) ||
            !all_finite(pressure_candidate) ||
            !all_finite(face_velocity_candidate) ||
            !all_finite(mass_flux_candidate)) {
          throw runtime::Error("Task 18 corrected trial is non-finite");
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
  return result;
}

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
       std::vector<ConstantDensityTransportSpec> supplied_transport)
      : decomposition(&supplied_decomposition), topology(&supplied_topology),
        geometry(&supplied_geometry), boundaries(&supplied_boundaries),
        mpi(&supplied_mpi), execution(&supplied_execution),
        halo(&supplied_halo), momentum_predictor(supplied_momentum_solver),
        momentum_preconditioners(supplied_preconditioners),
        coupler(std::move(supplied_coupler)),
        transport(std::move(supplied_transport)),
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
  finite_volume::CellCenteredFvmOperators fvm;
  TimeConsistentFaceVelocity face_assembler;
  ScratchFields scratch;
  std::array<execution::Buffer, 3> rhs;
  std::array<execution::Buffer, 3> predictor;
  std::array<execution::Buffer, 3> diagonal;
  std::array<std::unique_ptr<DiagonalMomentumOperator>, 3> operators;
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
  if (geometry.mapping_kind() != mesh::MappingKind::uniform_box) {
    throw runtime::Error(
        "Task 18 fixed-step flow supports uniform geometry only");
  }
  if (!cell_halo.is_compatible_with(decomposition) ||
      cell_halo.ghost_width() != 2) {
    throw runtime::Error("Task 18 fixed-step flow requires width-two Halo");
  }
  if (std::any_of(momentum_preconditioners.begin(),
                  momentum_preconditioners.end(),
                  [](const auto *item) { return item == nullptr; })) {
    throw runtime::Error("Task 18 momentum preconditioner is null");
  }
  std::set<runtime::FieldId> unique_transport;
  for (const auto &item : transported_fields) {
    if (!(item.diffusivity_kg_per_m_s >= 0.0) ||
        !std::isfinite(item.diffusivity_kg_per_m_s) ||
        !unique_transport.insert(item.field).second) {
      throw runtime::Error(
          "Task 18 transported-field specification is invalid");
    }
  }
  auto coupler = PisoCoupler::create(
      decomposition, topology, geometry, boundaries, mpi, execution_context,
      cell_halo, pressure_solver, pressure_preconditioner);
  return FixedStepConstantDensityFlow(std::make_unique<Impl>(
      decomposition, topology, geometry, boundaries, mpi, execution_context,
      cell_halo, momentum_solver, momentum_preconditioners, std::move(coupler),
      std::move(transported_fields)));
}

FixedStepConstantDensityFlow::FixedStepConstantDensityFlow(
    std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
FixedStepConstantDensityFlow::~FixedStepConstantDensityFlow() noexcept =
    default;
FixedStepConstantDensityFlow::FixedStepConstantDensityFlow(
    FixedStepConstantDensityFlow &&) noexcept = default;

namespace {

template <class FlowImplementation>
std::vector<double> assemble_spatial_residual(
    FlowImplementation &impl, const runtime::FieldRegistry &registry,
    const runtime::FieldAccessPlan &state_access,
    runtime::FieldStorage &accepted, const FlowFieldIds &fields, double mu) {
  synchronized_local_phase(
      *impl.mpi, StepFailureReason::non_finite_trial, true,
      "Task 18 accepted momentum gradient failed", [&] {
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
  std::vector<double> values;
  synchronized_local_phase(
      *impl.mpi, StepFailureReason::non_finite_trial, true,
      "Task 18 accepted momentum residual failed", [&] {
        const auto velocity = accepted.acquire_read<double>(
            state_access, kStatePhase, kStateActor, fields.velocity);
        const auto gradient_read =
            impl.scratch.storage->template acquire_read<double>(
                *impl.scratch.access, kScratchPhase, kScratchActor,
                impl.scratch.velocity_gradient);
        const auto flux = finite_volume::FaceMassFlux::acquire(
            registry, accepted, state_access, kStatePhase, kStateActor,
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
        impl.fvm.accumulate_viscous_residual(*impl.boundaries, velocity,
                                             gradient_read, mu, residual);
        values.resize(impl.topology->owned_cell_count() * 3U);
        const auto owned = impl.topology->owned_global_box();
        for (LocalCellId cell = 0; cell < impl.topology->owned_cell_count();
             ++cell) {
          const StructuredIndex index =
              map_cell(impl.topology->global_cell(cell), owned,
                       impl.topology->global_extent());
          for (int component_index = 0; component_index < 3;
               ++component_index) {
            values[cell * 3U + static_cast<std::size_t>(component_index)] =
                residual(index.i, index.j, index.k, component_index);
          }
        }
      });
  return values;
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
void recompute_transport(FlowImplementation &impl, FlowState &state,
                         double rho_ref, const MomentumTimeStencil &stencil) {
  if (impl.transport.empty())
    return;
  auto &committed =
      detail::FlowStateSolverAccess::layer(state, FlowLayer::committed);
  auto &history =
      detail::FlowStateSolverAccess::layer(state, FlowLayer::history);
  auto &trial = detail::FlowStateSolverAccess::layer(state, FlowLayer::trial);
  const auto &registry = detail::FlowStateSolverAccess::registry(state);
  const auto &access = detail::FlowStateSolverAccess::access(state);
  for (const auto &item : impl.transport) {
    synchronized_local_phase(
        *impl.mpi, StepFailureReason::transport_failure, true,
        "Task 18 transport gradient failed", [&] {
          const auto current = committed.acquire_read<double>(
              access, kStatePhase, kStateActor, item.field);
          auto gradient = impl.scratch.storage->template acquire_write<double>(
              *impl.scratch.access, kScratchPhase, kScratchActor,
              impl.scratch.scalar_gradient);
          impl.fvm.compute_gradient(finite_volume::GradientScheme::green_gauss,
                                    item.quantity, *impl.boundaries, current,
                                    gradient);
        });
    impl.halo->exchange(*impl.scratch.storage, impl.scratch.scalar_gradient);
    synchronized_local_phase(
        *impl.mpi, StepFailureReason::transport_failure, true,
        "Task 18 transport finalization failed", [&] {
          const auto current = committed.acquire_read<double>(
              access, kStatePhase, kStateActor, item.field);
          const auto previous = history.acquire_read<double>(
              access, kStatePhase, kStateActor, item.field);
          const auto flux = finite_volume::FaceMassFlux::acquire(
              registry, trial, access, kStatePhase, kStateActor,
              state.fields().face_mass_flux, *impl.topology);
          const auto gradient_read =
              impl.scratch.storage->template acquire_read<double>(
                  *impl.scratch.access, kScratchPhase, kScratchActor,
                  impl.scratch.scalar_gradient);
          auto face = impl.scratch.storage->template acquire_face_write<double>(
              *impl.scratch.access, kScratchPhase, kScratchActor,
              impl.scratch.scalar_face);
          impl.fvm.reconstruct_transport_faces(item.quantity, *impl.boundaries,
                                               flux, current, face);
          auto gamma =
              impl.scratch.storage->template acquire_face_write<double>(
                  *impl.scratch.access, kScratchPhase, kScratchActor,
                  impl.scratch.scalar_gamma);
          for (LocalFaceId local_face = 0;
               local_face < impl.topology->local_face_count(); ++local_face) {
            gamma(local_face, 0) = item.diffusivity_kg_per_m_s;
          }
          const auto face_read =
              impl.scratch.storage->template acquire_face_read<double>(
                  *impl.scratch.access, kScratchPhase, kScratchActor,
                  impl.scratch.scalar_face);
          const auto gamma_read =
              impl.scratch.storage->template acquire_face_read<double>(
                  *impl.scratch.access, kScratchPhase, kScratchActor,
                  impl.scratch.scalar_gamma);
          auto residual = impl.scratch.storage->template acquire_write<double>(
              *impl.scratch.access, kScratchPhase, kScratchActor,
              impl.scratch.scalar_residual);
          zero_cell(residual);
          impl.fvm.accumulate_convective_residual(flux, face_read, residual);
          impl.fvm.accumulate_scalar_diffusive_residual(
              item.quantity, *impl.boundaries, current, gradient_read,
              gamma_read, residual);
          auto output = trial.acquire_write<double>(access, kStatePhase,
                                                    kStateActor, item.field);
          const auto owned = impl.topology->owned_global_box();
          for (LocalCellId cell = 0; cell < impl.topology->owned_cell_count();
               ++cell) {
            const StructuredIndex index =
                map_cell(impl.topology->global_cell(cell), owned,
                         impl.topology->global_extent());
            const double value =
                (-(stencil.alpha1 * current(index.i, index.j, index.k, 0) +
                   stencil.alpha2 * previous(index.i, index.j, index.k, 0)) -
                 stencil.dt_s * residual(index.i, index.j, index.k, 0) /
                     (rho_ref * impl.geometry->cell_volume_m3(cell))) /
                stencil.alpha0;
            if (!std::isfinite(value)) {
              throw runtime::Error(
                  "Task 18 transport finalization is non-finite");
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
      "Task 18 continuity derivation failed", [&] {
        const auto flux = finite_volume::FaceMassFlux::acquire(
            detail::FlowStateSolverAccess::registry(state), trial,
            detail::FlowStateSolverAccess::access(state), kStatePhase,
            kStateActor, state.fields().face_mass_flux, *impl.topology);
        auto raw = impl.scratch.storage->template acquire_write<double>(
            *impl.scratch.access, kScratchPhase, kScratchActor,
            impl.scratch.mass_residual);
        zero_cell(raw);
        impl.fvm.accumulate_mass_residual(flux, raw);
        std::vector<double> absolute(impl.topology->owned_cell_count(), 0.0);
        const auto flux_view = trial.acquire_face_read<double>(
            detail::FlowStateSolverAccess::access(state), kStatePhase,
            kStateActor, state.fields().face_mass_flux);
        for (LocalFaceId face = 0; face < impl.topology->local_face_count();
             ++face) {
          const double magnitude = std::abs(flux_view(face, 0));
          const LocalCellId owner = impl.topology->owner(face);
          if (impl.topology->cell_ownership(owner) == EntityOwnership::owned) {
            absolute[owner] += magnitude;
          }
          const auto neighbour = impl.topology->neighbour(face);
          if (!impl.topology->patch_id(face).has_value() &&
              neighbour.has_value() &&
              impl.topology->cell_ownership(*neighbour) ==
                  EntityOwnership::owned) {
            absolute[*neighbour] += magnitude;
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
          sums[1] += absolute[cell] * absolute[cell];
        }
        if (!std::isfinite(sums[0]) || !std::isfinite(sums[1])) {
          throw runtime::Error("Task 18 continuity norm is non-finite");
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

} // namespace

StepAttemptReport FixedStepConstantDensityFlow::attempt(
    FlowState &state, double rho_ref, double mu,
    const MomentumTimeStencil &stencil,
    const linear::SolveControl &momentum_control,
    const linear::SolveControl &pressure_control) const {
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
          impl_->geometry->mapping_kind() == mesh::MappingKind::uniform_box &&
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
      }
      const auto committed = state.snapshot(FlowLayer::committed);
      const auto history = state.snapshot(FlowLayer::history);
      local_valid =
          local_valid &&
          std::all_of(committed.density.begin(), committed.density.end(),
                      [&](double value) { return value == rho_ref; }) &&
          std::all_of(history.density.begin(), history.density.end(),
                      [&](double value) { return value == rho_ref; });
    } catch (const runtime::detail::MpiOperationError &) {
      throw;
    } catch (...) {
      local_valid = false;
    }
    const auto validation = runtime::collective_status(
        *impl_->mpi, local_valid, "Task 18 fixed-step input is invalid");
    if (!validation.ok) {
      return fatal_failure(report, StepFailureReason::invalid_input,
                           validation.failing_rank);
    }
    const auto metadata = state.metadata();
    const std::array<double, 14> core_inputs{
        rho_ref,
        mu,
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
                          "Task 18 fixed-step rank inputs are inconsistent");
    if (!core_agreement.ok) {
      return fatal_failure(report, StepFailureReason::invalid_input,
                           core_agreement.failing_rank);
    }
    const auto transport_count =
        static_cast<std::uint64_t>(impl_->transport.size());
    const std::array<double, 2> transport_count_input{
        low_u32(transport_count), high_u32(transport_count)};
    const auto transport_count_agreement =
        agree_fp64_inputs(*impl_->mpi, transport_count_input,
                          "Task 18 transported-field counts are inconsistent");
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
          "Task 18 transported-field inputs are inconsistent");
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
        *impl_->mpi, begin_ok, "Task 18 trial initialization failed");
    if (!begin_status.ok) {
      if (active) {
        state.rollback_attempt();
        active = false;
      }
      return fatal_failure(report, StepFailureReason::invalid_input,
                           begin_status.failing_rank);
    }
    exchange_accepted_layers(*impl_, state);
    auto &committed = state.solver_layer(FlowLayer::committed);
    auto &history = state.solver_layer(FlowLayer::history);
    const auto &registry = state.solver_registry();
    const auto &access = state.solver_access_plan();
    const auto fields = state.fields();
    const auto residual_n = assemble_spatial_residual(*impl_, registry, access,
                                                      committed, fields, mu);
    const auto residual_nm1 = assemble_spatial_residual(
        *impl_, registry, access, history, fields, mu);
    std::optional<runtime::FieldView<const double>> velocity_n;
    std::optional<runtime::FieldView<const double>> velocity_nm1;
    std::optional<runtime::FieldView<const double>> pressure_n;
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::non_finite_trial, true,
        "Task 18 accepted field acquisition failed", [&] {
          velocity_n.emplace(committed.acquire_read<double>(
              access, kStatePhase, kStateActor, fields.velocity));
          velocity_nm1.emplace(history.acquire_read<double>(
              access, kStatePhase, kStateActor, fields.velocity));
          pressure_n.emplace(committed.acquire_read<double>(
              access, kStatePhase, kStateActor, fields.mechanical_pressure));
        });
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::non_finite_trial, true,
        "Task 18 pressure-gradient derivation failed", [&] {
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
          if (force_local_derived_failure.load(std::memory_order_relaxed)) {
            throw runtime::Error(
                "injected Task 18 local derived-value failure");
          }
#endif
          auto pressure_gradient =
              impl_->scratch.storage->acquire_write<double>(
                  *impl_->scratch.access, kScratchPhase, kScratchActor,
                  impl_->scratch.pressure_gradient);
          compute_pressure_gradient(
              *impl_->topology, *impl_->geometry, *pressure_n,
              pressure_gradient, [&](LocalFaceId face, double owner_value) {
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
    std::array<std::vector<double>, 3> diagonal_values;
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::momentum_linear_solve, true,
        "Task 18 momentum equation preparation failed", [&] {
          pressure_gradient_read.emplace(
              impl_->scratch.storage->acquire_read<double>(
                  *impl_->scratch.access, kScratchPhase, kScratchActor,
                  impl_->scratch.pressure_gradient));
          for (auto &values : diagonal_values)
            values.resize(count);
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
              const double spatial =
                  stencil.order == MomentumTimeOrder::backward_euler
                      ? residual_n[cell * 3U + component_index]
                      : 2.0 * residual_n[cell * 3U + component_index] -
                            residual_nm1[cell * 3U + component_index];
              rhs[cell] =
                  -(rho_ref * volume / stencil.dt_s) *
                      (stencil.alpha1 *
                           (*velocity_n)(index.i, index.j, index.k,
                                         static_cast<int>(component_index)) +
                       stencil.alpha2 *
                           (*velocity_nm1)(index.i, index.j, index.k,
                                           static_cast<int>(component_index))) -
                  spatial -
                  volume * (*pressure_gradient_read)(
                               index.i, index.j, index.k,
                               static_cast<int>(component_index));
              predictor[cell] = (*velocity_n)(
                  index.i, index.j, index.k, static_cast<int>(component_index));
              diagonal_values[component_index][cell] = diagonal;
              if (!std::isfinite(rhs[cell]) || !(diagonal > 0.0) ||
                  !std::isfinite(diagonal)) {
                throw runtime::Error("Task 18 momentum equation is non-finite");
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
        "Task 18 momentum solve view acquisition failed", [&] {
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
    auto &trial = state.solver_layer(FlowLayer::trial);
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::non_finite_trial, true,
        "Task 18 predictor staging failed", [&] {
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
                throw runtime::Error("Task 18 predictor staging is non-finite");
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
    std::optional<runtime::FieldView<const double>> actual_diagonal_read;
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::non_finite_trial, true,
        "Task 18 face predictor assembly failed", [&] {
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
    const auto first = impl_->coupler.correct(
        state, rho_ref, *actual_diagonal_read, pressure_control);
    report.pressure[0] = first.solve;
    if (!first.accepted) {
      state.rollback_attempt();
      active = false;
      return numerical_failure(report, StepFailureReason::pressure_linear_solve,
                               first.solve.lowest_failing_rank);
    }
    report.pressure_corrector_count = 1U;
    recompute_transport(*impl_, state, rho_ref, stencil);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::transport_failure, true,
        "Task 18 provisional transport test staging failed", [&] {
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
    const auto second = impl_->coupler.correct(
        state, rho_ref, *actual_diagonal_read, pressure_control);
    report.pressure[1] = second.solve;
    if (!second.accepted) {
      state.rollback_attempt();
      active = false;
      return numerical_failure(report, StepFailureReason::pressure_linear_solve,
                               second.solve.lowest_failing_rank);
    }
    report.pressure_corrector_count = 2U;
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::transport_failure, true,
        "Task 18 final flux test staging failed", [&] {
          const double injected_flux =
              final_uniform_x_mass_flux.load(std::memory_order_relaxed);
          if (injected_flux != 0.0) {
            auto final_flux = trial.acquire_face_write<double>(
                access, kStatePhase, kStateActor, fields.face_mass_flux);
            for (LocalFaceId face = 0;
                 face < impl_->topology->local_face_count(); ++face) {
              final_flux(face, 0) +=
                  injected_flux *
                  impl_->geometry->face_area_vector_m2(face, FaceSide::owner).x;
            }
          }
        });
#endif
    recompute_transport(*impl_, state, rho_ref, stencil);
#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
    final_transport_call_count.fetch_add(1U, std::memory_order_relaxed);
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
        "Task 18 final pressure residual exceeds tolerance");
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
        "Task 18 final continuity residual exceeds tolerance");
    if (!continuity_status.ok) {
      state.rollback_attempt();
      active = false;
      return numerical_failure(report,
                               StepFailureReason::final_continuity_residual,
                               continuity_status.failing_rank);
    }
    std::optional<runtime::FaceFieldView<const double>> final_flux;
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::boundary_backflow, true,
        "Task 18 final flux acquisition failed", [&] {
          final_flux.emplace(trial.acquire_face_read<double>(
              access, kStatePhase, kStateActor, fields.face_mass_flux));
        });
    const auto admissibility =
        impl_->boundaries->assess_final_pressure_outlet_flux(
            *impl_->topology, *impl_->mpi, *final_flux,
            state.metadata().step + 1U, state.metadata().time_s + stencil.dt_s);
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
    synchronized_local_phase(
        *impl_->mpi, StepFailureReason::non_finite_trial, true,
        "Task 18 collective commit preparation failed", [&] {
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
                         impl_->mpi->rank());
  } catch (const runtime::Error &) {
    if (active) {
      state.rollback_attempt();
      return numerical_failure(report, StepFailureReason::non_finite_trial,
                               impl_->mpi->rank());
    }
    return fatal_failure(report, StepFailureReason::invalid_input,
                         impl_->mpi->rank());
  } catch (...) {
    if (active)
      state.rollback_attempt();
    return fatal_failure(report, StepFailureReason::invalid_input,
                         impl_->mpi->rank());
  }
}

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
void test::ConstantDensityPisoTestAccess::reset() noexcept {
  ::hundun::flow::force_final_continuity_failure.store(
      false, std::memory_order_relaxed);
  ::hundun::flow::force_final_pressure_failure.store(false,
                                                     std::memory_order_relaxed);
  ::hundun::flow::force_local_derived_failure.store(false,
                                                    std::memory_order_relaxed);
  provisional_transport_sentinel.store(false, std::memory_order_relaxed);
  final_uniform_x_mass_flux.store(0.0, std::memory_order_relaxed);
  for (auto &value : ::hundun::flow::last_momentum_rhs) {
    value.store(0.0, std::memory_order_relaxed);
  }
  for (auto &value : ::hundun::flow::last_momentum_diagonal) {
    value.store(0.0, std::memory_order_relaxed);
  }
  provisional_transport_call_count.store(0U, std::memory_order_relaxed);
  final_transport_call_count.store(0U, std::memory_order_relaxed);
}

void test::ConstantDensityPisoTestAccess::force_final_continuity_failure(
    bool enabled) noexcept {
  ::hundun::flow::force_final_continuity_failure.store(
      enabled, std::memory_order_relaxed);
}

void test::ConstantDensityPisoTestAccess::force_final_pressure_failure(
    bool enabled) noexcept {
  ::hundun::flow::force_final_pressure_failure.store(enabled,
                                                     std::memory_order_relaxed);
}

void test::ConstantDensityPisoTestAccess::force_local_derived_failure(
    bool enabled) noexcept {
  ::hundun::flow::force_local_derived_failure.store(enabled,
                                                    std::memory_order_relaxed);
}

void test::ConstantDensityPisoTestAccess::set_provisional_transport_sentinel(
    bool enabled) noexcept {
  provisional_transport_sentinel.store(enabled, std::memory_order_relaxed);
}

void test::ConstantDensityPisoTestAccess::set_final_uniform_x_mass_flux(
    double value) noexcept {
  final_uniform_x_mass_flux.store(value, std::memory_order_relaxed);
}

double test::ConstantDensityPisoTestAccess::last_momentum_rhs(
    std::size_t component) noexcept {
  return component < ::hundun::flow::last_momentum_rhs.size()
             ? ::hundun::flow::last_momentum_rhs[component].load(
                   std::memory_order_relaxed)
             : std::numeric_limits<double>::quiet_NaN();
}

double test::ConstantDensityPisoTestAccess::last_momentum_diagonal(
    std::size_t component) noexcept {
  return component < ::hundun::flow::last_momentum_diagonal.size()
             ? ::hundun::flow::last_momentum_diagonal[component].load(
                   std::memory_order_relaxed)
             : std::numeric_limits<double>::quiet_NaN();
}

std::size_t
test::ConstantDensityPisoTestAccess::provisional_transport_calls() noexcept {
  return provisional_transport_call_count.load(std::memory_order_relaxed);
}

std::size_t
test::ConstantDensityPisoTestAccess::final_transport_calls() noexcept {
  return final_transport_call_count.load(std::memory_order_relaxed);
}
#endif

} // namespace hundun::flow

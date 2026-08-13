// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_linear.hpp"
#include "hundun/v04_parallel.hpp"

#include "field_view_interval_detail.hpp"
#include "solver_mg_detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kMgPlan = 7101U;
constexpr std::uint32_t kMgCoefficient = 7102U;
constexpr std::uint32_t kMgCollective = 7103U;
constexpr std::uint32_t kMgApply = 7104U;
constexpr std::uint32_t kMgCounter = 7105U;

std::uint64_t mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= UINT64_C(1099511628211);
  return hash;
}

PlanFingerprint finish(std::uint64_t hash) noexcept {
  return hash == 0U ? PlanFingerprint{1U} : hash;
}

bool same_shape(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

std::size_t cell_count(Int3 shape) noexcept {
  if (shape.x <= 0 || shape.y <= 0 || shape.z <= 0) {
    return 0U;
  }
  const std::size_t x = static_cast<std::size_t>(shape.x);
  const std::size_t y = static_cast<std::size_t>(shape.y);
  const std::size_t z = static_cast<std::size_t>(shape.z);
  if (x > std::numeric_limits<std::size_t>::max() / y ||
      x * y > std::numeric_limits<std::size_t>::max() / z) {
    return 0U;
  }
  return x * y * z;
}

bool valid_identity(LinearIdentity identity) noexcept {
  return identity.symbolic != 0U && identity.numeric != 0U &&
         identity.hierarchy != 0U && identity.workspace != 0U &&
         identity.fingerprint != 0U;
}

bool valid_boundary_pair(MgBoundaryKind minimum,
                         MgBoundaryKind maximum) noexcept {
  const bool min_periodic = minimum == MgBoundaryKind::periodic;
  const bool max_periodic = maximum == MgBoundaryKind::periodic;
  return min_periodic == max_periodic;
}

HaloTopology halo_topology(MgBoundarySet boundaries) noexcept {
  return {
      boundaries.x_min == MgBoundaryKind::periodic,
      boundaries.y_min == MgBoundaryKind::periodic,
      boundaries.z_min == MgBoundaryKind::periodic,
  };
}

bool valid_policy(MgHierarchyPolicy policy) noexcept {
  return std::isfinite(policy.anisotropy_threshold) &&
         policy.anisotropy_threshold > 1.0 &&
         std::isfinite(policy.coefficient_change_rebuild_ratio) &&
         policy.coefficient_change_rebuild_ratio >= 0.0 &&
         policy.pre_sweeps != 0U && policy.post_sweeps != 0U &&
         policy.maximum_levels >= 2U &&
         policy.maximum_levels <= detail::kMgMaximumLevels &&
         policy.coarse_sweeps != 0U && policy.minimum_coarse_extent >= 2U &&
         policy.line_relaxation_maximum_extent >= 2U;
}

bool valid_field(ConstFieldView field, Int3 shape) noexcept {
  return field.base != nullptr && field.components == 1U &&
         same_shape(field.interior, shape) && field.ghosts.x >= 0 &&
         field.ghosts.y >= 0 && field.ghosts.z >= 0 &&
         field.stride_y >= static_cast<std::size_t>(shape.x) &&
         field.stride_z >=
             field.stride_y * static_cast<std::size_t>(shape.y) &&
         field.component_stride >=
             field.stride_z * static_cast<std::size_t>(shape.z) &&
         field.storage_identity != 0U && field.revision_domain != 0U;
}

bool valid_face(ConstFaceFieldView field, CartesianAxis axis,
                Int3 cells) noexcept {
  Int3 expected = cells;
  if (axis == CartesianAxis::x) {
    ++expected.x;
  } else if (axis == CartesianAxis::y) {
    ++expected.y;
  } else {
    ++expected.z;
  }
  return field.base != nullptr && field.axis == axis &&
         same_shape(field.extents, expected) &&
         field.stride_y >= static_cast<std::size_t>(expected.x) &&
         field.stride_z >=
             field.stride_y * static_cast<std::size_t>(expected.y) &&
         field.storage_identity != 0U && field.revision_domain != 0U;
}

bool finite_coefficients(MgCoefficientViews coefficients,
                         Int3 cells) noexcept {
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const double value = coefficients.diagonal.unchecked({i, j, k}, 0U);
        if (!(value > 0.0) || !std::isfinite(value)) {
          return false;
        }
      }
    }
  }
  const auto finite_face = [](ConstFaceFieldView field) noexcept {
    for (std::int32_t k = 0; k < field.extents.z; ++k) {
      for (std::int32_t j = 0; j < field.extents.y; ++j) {
        for (std::int32_t i = 0; i < field.extents.x; ++i) {
          const double value = field.unchecked({i, j, k});
          if (value < 0.0 || !std::isfinite(value)) {
            return false;
          }
        }
      }
    }
    return true;
  };
  return finite_face(coefficients.x) && finite_face(coefficients.y) &&
         finite_face(coefficients.z);
}

bool valid_coefficients(MgCoefficientViews coefficients,
                        Int3 cells) noexcept {
  return valid_field(coefficients.diagonal, cells) &&
         valid_face(coefficients.x, CartesianAxis::x, cells) &&
         valid_face(coefficients.y, CartesianAxis::y, cells) &&
         valid_face(coefficients.z, CartesianAxis::z, cells) &&
         finite_coefficients(coefficients, cells);
}

bool valid_services(MgRuntimeServices services,
                    LinearIdentity) noexcept {
  return services.finest_halo != nullptr && services.reductions != nullptr &&
         services.workspace != nullptr && services.finest_halo->ready() &&
         services.reductions->capacity() != 0U &&
         services.workspace->fingerprint() != 0U;
}

double axis_anisotropy(const AxisMetrics& metric) noexcept {
  const Span<const double> widths = metric.widths();
  if (widths.data == nullptr || widths.size == 0U) {
    return std::numeric_limits<double>::infinity();
  }
  double minimum = std::numeric_limits<double>::infinity();
  double maximum = 0.0;
  for (std::size_t index = 0U; index < widths.size; ++index) {
    minimum = std::min(minimum, widths.data[index]);
    maximum = std::max(maximum, widths.data[index]);
  }
  return minimum > 0.0 ? maximum / minimum
                       : std::numeric_limits<double>::infinity();
}

double axis_minimum_width(const AxisMetrics& metric) noexcept {
  const Span<const double> widths = metric.widths();
  double minimum = std::numeric_limits<double>::infinity();
  for (std::size_t index = 0U; index < widths.size; ++index) {
    minimum = std::min(minimum, widths.data[index]);
  }
  return minimum;
}

struct Strategy {
  CoarseningKind coarsening{CoarseningKind::full_xyz};
  std::uint8_t coarsen_mask{detail::kMgAxisX | detail::kMgAxisY |
                            detail::kMgAxisZ};
  std::uint8_t line_mask{};
};

Strategy choose_strategy(const CartesianGeometryPlan& geometry,
                         MgHierarchyPolicy policy) noexcept {
  const double hx = axis_minimum_width(geometry.x());
  const double hy = axis_minimum_width(geometry.y());
  const double hz = axis_minimum_width(geometry.z());
  const double largest_minimum = std::max({hx, hy, hz});
  const bool strong_x =
      axis_anisotropy(geometry.x()) >= policy.anisotropy_threshold ||
      largest_minimum / hx >= policy.anisotropy_threshold;
  const bool strong_y =
      axis_anisotropy(geometry.y()) >= policy.anisotropy_threshold ||
      largest_minimum / hy >= policy.anisotropy_threshold;
  const bool strong_z =
      axis_anisotropy(geometry.z()) >= policy.anisotropy_threshold ||
      largest_minimum / hz >= policy.anisotropy_threshold;
  const std::uint8_t strong =
      (strong_x ? detail::kMgAxisX : 0U) |
      (strong_y ? detail::kMgAxisY : 0U) |
      (strong_z ? detail::kMgAxisZ : 0U);
  switch (strong) {
    case detail::kMgAxisX:
      return {CoarseningKind::semi_yz,
              detail::kMgAxisY | detail::kMgAxisZ, strong};
    case detail::kMgAxisY:
      return {CoarseningKind::semi_xz,
              detail::kMgAxisX | detail::kMgAxisZ, strong};
    case detail::kMgAxisZ:
      return {CoarseningKind::semi_xy,
              detail::kMgAxisX | detail::kMgAxisY, strong};
    case detail::kMgAxisX | detail::kMgAxisY:
      return {CoarseningKind::z_only, detail::kMgAxisZ, strong};
    case detail::kMgAxisX | detail::kMgAxisZ:
      return {CoarseningKind::y_only, detail::kMgAxisY, strong};
    case detail::kMgAxisY | detail::kMgAxisZ:
      return {CoarseningKind::x_only, detail::kMgAxisX, strong};
    default:
      return {};
  }
}

PlanFingerprint structural_contract(const NativeCartesianMgSpec& spec,
                                    Strategy strategy) noexcept {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  hash = mix(hash, spec.geometry == nullptr ? 0U : spec.geometry->fingerprint());
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.x_min));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.x_max));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.y_min));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.y_max));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.z_min));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.z_max));
  hash = mix(hash, static_cast<std::uint64_t>(spec.null_space));
  hash = mix(hash, static_cast<std::uint64_t>(strategy.coarsening));
  hash = mix(hash, strategy.line_mask);
  hash = mix(hash, spec.policy.pre_sweeps);
  hash = mix(hash, spec.policy.post_sweeps);
  hash = mix(hash, spec.policy.maximum_levels);
  hash = mix(hash, spec.policy.coarse_sweeps);
  hash = mix(hash, spec.policy.minimum_coarse_extent);
  hash = mix(hash, spec.policy.line_relaxation_maximum_extent);
  std::uint64_t threshold_bits = 0U;
  std::uint64_t rebuild_bits = 0U;
  static_assert(sizeof(threshold_bits) == sizeof(double));
  std::memcpy(&threshold_bits, &spec.policy.anisotropy_threshold,
              sizeof(threshold_bits));
  std::memcpy(&rebuild_bits,
              &spec.policy.coefficient_change_rebuild_ratio,
              sizeof(rebuild_bits));
  hash = mix(hash, threshold_bits);
  hash = mix(hash, rebuild_bits);
  return finish(hash);
}

PlanFingerprint public_symbolic_fingerprint(
    const NativeCartesianMgSpec& spec, Strategy strategy) noexcept {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  hash = mix(hash, spec.geometry->fingerprint());
  hash = mix(hash, spec.identity.symbolic);
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.x_min));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.x_max));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.y_min));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.y_max));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.z_min));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.z_max));
  hash = mix(hash, static_cast<std::uint64_t>(strategy.coarsening));
  return finish(hash);
}

PlanFingerprint make_numeric_fingerprint(MgCoefficientIdentity identity,
                                         PlanFingerprint symbolic) noexcept {
  return finish(mix(mix(mix(UINT64_C(1469598103934665603), symbolic),
                        identity.revision),
                    identity.fingerprint));
}

PlanFingerprint make_hierarchy_fingerprint(
    PlanFingerprint symbolic, MgCoefficientIdentity identity,
    RevisionToken generation) noexcept {
  return finish(mix(mix(mix(UINT64_C(1469598103934665603), symbolic),
                        identity.fingerprint),
                    generation));
}

Status validate_compile(const NativeCartesianMgSpec& spec,
                        MgRuntimeServices services,
                        MgCoefficientViews coefficients,
                        Strategy& strategy) noexcept {
  if (spec.communicator == MPI_COMM_NULL || spec.geometry == nullptr ||
      !same_shape(spec.patch.cells, coefficients.diagonal.interior) ||
      !valid_identity(spec.identity) || spec.coefficients.revision == 0U ||
      spec.coefficients.fingerprint == 0U ||
      !std::isfinite(spec.coefficients.maximum_relative_change) ||
      spec.coefficients.maximum_relative_change < 0.0 ||
      !valid_policy(spec.policy) ||
      !valid_boundary_pair(spec.boundaries.x_min, spec.boundaries.x_max) ||
      !valid_boundary_pair(spec.boundaries.y_min, spec.boundaries.y_max) ||
      !valid_boundary_pair(spec.boundaries.z_min, spec.boundaries.z_max) ||
      !valid_services(services, spec.identity) ||
      !valid_coefficients(coefficients, spec.patch.cells)) {
    return {StatusCode::invalid_plan, kMgPlan};
  }
  strategy = choose_strategy(*spec.geometry, spec.policy);
  return {};
}

bool increment(std::uint64_t current, std::uint64_t& next) noexcept {
  if (current == std::numeric_limits<std::uint64_t>::max()) {
    return false;
  }
  next = current + 1U;
  return true;
}

bool valid_status_code(StatusCode code) noexcept {
  return static_cast<std::uint16_t>(code) <=
         static_cast<std::uint16_t>(StatusCode::io_failure);
}

Status raw_consensus(MPI_Comm communicator, Status local,
                     int* lowest = nullptr) noexcept {
  if (communicator == MPI_COMM_NULL) {
    return local ? Status{StatusCode::invalid_plan, kMgPlan} : local;
  }
  if (!valid_status_code(local.code)) {
    local = {StatusCode::invalid_plan, kMgPlan};
  }
  int rank = 0;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    if (lowest != nullptr) *lowest = -1;
    return {StatusCode::mpi_failure, kMgCollective};
  }
  const int candidate = local ? size : rank;
  int selected = size;
  if (MPI_Allreduce(&candidate, &selected, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    if (lowest != nullptr) *lowest = -1;
    return {StatusCode::mpi_failure, kMgCollective};
  }
  if (selected == size) {
    if (lowest != nullptr) *lowest = -1;
    return {};
  }
  std::uint64_t wire[2]{};
  if (rank == selected) {
    wire[0] = static_cast<std::uint64_t>(local.code);
    wire[1] = local.detail;
  }
  if (MPI_Bcast(wire, 2, MPI_UINT64_T, selected, communicator) !=
          MPI_SUCCESS ||
      wire[0] > static_cast<std::uint64_t>(StatusCode::io_failure) ||
      wire[1] > std::numeric_limits<std::uint32_t>::max()) {
    if (lowest != nullptr) *lowest = -1;
    return {StatusCode::mpi_failure, kMgCollective};
  }
  if (lowest != nullptr) *lowest = selected;
  return {static_cast<StatusCode>(wire[0]),
          static_cast<std::uint32_t>(wire[1])};
}

Status validate_runtime_communicators(const NativeCartesianMgSpec& spec,
                                      MgRuntimeServices services) noexcept {
  Status local = services.reductions == nullptr ||
                         services.workspace == nullptr ||
                         services.finest_halo == nullptr
                     ? Status{StatusCode::invalid_plan, kMgPlan}
                     : services.reductions->validate_communicator(
                           spec.communicator);
  if (!local) return local;
  const MgWorkspaceRequirements& requirements =
      services.workspace->requirements();
  if (requirements.level_count < 2U ||
      services.level_halos.size + 1U != requirements.level_count ||
      services.workspace->level(0U, MgWorkspaceSlot::solution).base ==
          nullptr) {
    return {StatusCode::invalid_plan, kMgPlan};
  }
  const FieldView workspace_field =
      services.workspace->level(0U, MgWorkspaceSlot::solution);
  const std::array<HaloFieldSpec, 1U> fields{{
      {workspace_field.field, 1U, 1U},
  }};
  const HaloTopology topology = halo_topology(spec.boundaries);
  local = services.finest_halo->validate_contract(
      spec.communicator, requirements.levels[0U].patch,
      {fields.data(), fields.size()}, topology);
  if (!local) return local;
  for (std::size_t level = 1U; level < requirements.level_count; ++level) {
    if (services.level_halos.data == nullptr ||
        services.level_halos.data[level - 1U] == nullptr) {
      return {StatusCode::invalid_plan, kMgPlan};
    }
    local = services.level_halos.data[level - 1U]->validate_contract(
        spec.communicator, requirements.levels[level].patch,
        {fields.data(), fields.size()}, topology);
    if (!local) return local;
  }
  return {};
}

}  // namespace

struct NativeCartesianMgPlan::Impl {
  MPI_Comm communicator{MPI_COMM_NULL};
  NativeCartesianMgSpec spec{};
  MgRuntimeServices services{};
  // Service objects remain caller-owned and must outlive the plan. Only this
  // pointer table is copied, so the caller's Span backing storage may change.
  std::vector<HaloEngine*> level_halo_table;
  std::vector<std::uintptr_t> level_halo_identities;
  HaloEngine* finest_halo_object{};
  ReductionEngine* reductions_object{};
  MgWorkspace* workspace_object{};
  std::uintptr_t finest_halo_identity{};
  std::uintptr_t reductions_identity{};
  RevisionToken workspace_binding_identity{};
  PlanFingerprint workspace_fingerprint{};
  PlanFingerprint workspace_collective_fingerprint{};
  std::uintptr_t workspace_storage{};
  std::size_t workspace_storage_doubles{};
  MgCoefficientViews coefficients{};
  std::vector<detail::MgLevelStorage> levels;
  std::vector<double> hierarchy_storage;
  std::vector<double> inactive_hierarchy_storage;
  // Compile/update-only storage for conservative tensor aggregation.  The
  // first two blocks are ping-pong fields; the final two are MPI planes.
  // Keeping it on the plan makes numeric refresh allocation-free.
  std::vector<double> coefficient_tensor_storage;
  std::size_t coefficient_tensor_values{};
  std::size_t coefficient_tensor_plane{};
  Strategy strategy{};
  MgPlanCounters runtime_counters{};
  PlanFingerprint symbolic{};
  PlanFingerprint numeric{};
  PlanFingerprint hierarchy{};
  RevisionToken generation{};
  double last_initial{};
  double last_final{};
  int rank{};
  int size{};
  int lowest{-1};
};

namespace {

template <class Implementation>
Status validate_borrowed_services(const Implementation& implementation) noexcept {
  // These are borrowed object addresses. The objects themselves must outlive
  // the plan; this check detects move/rebind while each object is still alive.
  if (implementation.services.finest_halo != implementation.finest_halo_object ||
      implementation.services.reductions != implementation.reductions_object ||
      implementation.services.workspace != implementation.workspace_object ||
      implementation.services.finest_halo == nullptr ||
      implementation.services.reductions == nullptr ||
      implementation.services.workspace == nullptr ||
      !implementation.services.finest_halo->ready() ||
      implementation.services.finest_halo->instance_identity() !=
          implementation.finest_halo_identity ||
      implementation.services.reductions->instance_identity() !=
          implementation.reductions_identity ||
      implementation.services.workspace->binding_identity() !=
          implementation.workspace_binding_identity ||
      implementation.services.workspace->fingerprint() !=
          implementation.workspace_fingerprint ||
      implementation.services.workspace->collective_fingerprint() !=
          implementation.workspace_collective_fingerprint ||
      implementation.services.workspace->storage_address() !=
          implementation.workspace_storage ||
      implementation.services.workspace->storage_doubles() !=
          implementation.workspace_storage_doubles ||
      implementation.level_halo_table.size() !=
          implementation.level_halo_identities.size()) {
    return {StatusCode::invalid_plan, kMgPlan};
  }
  for (std::size_t level = 0U;
       level < implementation.level_halo_table.size(); ++level) {
    const HaloEngine* const halo = implementation.level_halo_table[level];
    if (halo == nullptr || !halo->ready() ||
        halo->instance_identity() != implementation.level_halo_identities[level]) {
      return {StatusCode::invalid_plan, kMgPlan};
    }
  }
  return {};
}

bool mpi_live() noexcept {
  int initialized = 0;
  int finalized = 0;
  (void)MPI_Initialized(&initialized);
  if (initialized != 0) {
    (void)MPI_Finalized(&finalized);
  }
  return initialized != 0 && finalized == 0;
}

template <class Implementation>
void destroy(Implementation* implementation) noexcept {
  if (implementation == nullptr) {
    return;
  }
  if (mpi_live() && implementation->communicator != MPI_COMM_NULL) {
    (void)MPI_Comm_free(&implementation->communicator);
  }
  delete implementation;
}

template <class Implementation>
Status consensus(Implementation& implementation, Status local) noexcept {
  const Status result = implementation.services.reductions->consensus(local);
  implementation.lowest =
      implementation.services.reductions->lowest_failing_rank();
  return result;
}

template <class Implementation>
Status build_levels(Implementation& implementation) {
  const MgWorkspaceRequirements& requirements =
      implementation.services.workspace->requirements();
  if (requirements.level_count < 2U ||
      requirements.level_count > detail::kMgMaximumLevels) {
    return {StatusCode::invalid_plan, kMgPlan};
  }
  std::size_t total_hierarchy = 0U;
  for (std::size_t level = 0U; level < requirements.level_count; ++level) {
    const MgWorkspaceLevelRequirements& required = requirements.levels[level];
    const Int3 local = required.patch.cells;
    const std::size_t cells = cell_count(local);
    if (cells == 0U || cells >
                           (std::numeric_limits<std::size_t>::max() -
                            total_hierarchy) /
                               4U) {
      return {StatusCode::invalid_plan, kMgPlan};
    }
    detail::MgLevelStorage storage;
    storage.view = {required.global_shape, local, required.coarsening,
                    required.line_axis_mask};
    storage.patch = required.patch;
    storage.coarsen_mask = required.coarsen_axis_mask;
    storage.cells = cells;
    storage.diagonal_offset = total_hierarchy;
    storage.x_faces = cell_count({local.x + 1, local.y, local.z});
    storage.y_faces = cell_count({local.x, local.y + 1, local.z});
    storage.z_faces = cell_count({local.x, local.y, local.z + 1});
    storage.x_offset = storage.diagonal_offset + cells;
    storage.y_offset = storage.x_offset + storage.x_faces;
    storage.z_offset = storage.y_offset + storage.y_faces;
    storage.volume_offset = storage.z_offset + storage.z_faces;
    implementation.levels.push_back(storage);
    total_hierarchy = storage.volume_offset + cells;
  }
  const MgWorkspaceLevelRequirements& finest = requirements.levels[0U];
  if (!same_shape(finest.global_shape,
                  implementation.spec.geometry->global_cells()) ||
      !same_shape(finest.patch.cells, implementation.spec.patch.cells) ||
      !same_shape(finest.patch.begin, implementation.spec.patch.begin)) {
    return {StatusCode::invalid_plan, kMgPlan};
  }
  implementation.strategy.coarsening = finest.coarsening;
  implementation.strategy.coarsen_mask = finest.coarsen_axis_mask;
  implementation.strategy.line_mask = finest.line_axis_mask;
  implementation.hierarchy_storage.assign(total_hierarchy, 0.0);
  implementation.inactive_hierarchy_storage.assign(total_hierarchy, 0.0);
  std::size_t maximum_values = 0U;
  std::size_t maximum_plane = 0U;
  for (const detail::MgLevelStorage& level : implementation.levels) {
    maximum_values = std::max(
        {maximum_values, level.cells, level.x_faces, level.y_faces,
         level.z_faces});
    const Int3 cells = level.view.local_shape;
    const std::size_t nx = static_cast<std::size_t>(cells.x) + 1U;
    const std::size_t ny = static_cast<std::size_t>(cells.y) + 1U;
    const std::size_t nz = static_cast<std::size_t>(cells.z) + 1U;
    maximum_plane =
        std::max({maximum_plane, nx * ny, nx * nz, ny * nz});
  }
  if (maximum_values == 0U || maximum_plane == 0U ||
      maximum_values >
          (std::numeric_limits<std::size_t>::max() - 2U * maximum_plane) /
              2U) {
    return {StatusCode::invalid_plan, kMgPlan};
  }
  implementation.coefficient_tensor_values = maximum_values;
  implementation.coefficient_tensor_plane = maximum_plane;
  implementation.coefficient_tensor_storage.assign(
      2U * maximum_values + 2U * maximum_plane, 0.0);
  return {};
}

template <class Implementation>
Status project_constant(Implementation& implementation, std::size_t level,
                        FieldView field) noexcept {
  const detail::MgLevelStorage& selected = implementation.levels[level];
  const double* const volumes = detail::block(
      implementation.hierarchy_storage.data(), selected.volume_offset);
  const Int3 cells = selected.view.local_shape;
  double local_weighted = 0.0;
  double local_volume = 0.0;
  for (std::int32_t k = 0; k < field.interior.z; ++k) {
    for (std::int32_t j = 0; j < field.interior.y; ++j) {
      for (std::int32_t i = 0; i < field.interior.x; ++i) {
        const std::size_t index = static_cast<std::size_t>(i) +
                                  static_cast<std::size_t>(cells.x) *
                                      (static_cast<std::size_t>(j) +
                                       static_cast<std::size_t>(cells.y) *
                                           static_cast<std::size_t>(k));
        const double volume = volumes[index];
        local_weighted += volume * field.unchecked({i, j, k}, 0U);
        local_volume += volume;
      }
    }
  }
  double local[2]{local_weighted, local_volume};
  double total[2]{};
  Status status = implementation.services.reductions->checked_sum(
      {local, 2U}, {total, 2U});
  implementation.lowest =
      implementation.services.reductions->lowest_failing_rank();
  if (!status || !(total[1] > 0.0) || !std::isfinite(total[0]) ||
      !std::isfinite(total[1])) {
    return status ? Status{StatusCode::numerical_failure, kMgApply} : status;
  }
  const double mean = total[0] / total[1];
  for (std::int32_t k = 0; k < field.interior.z; ++k)
    for (std::int32_t j = 0; j < field.interior.y; ++j)
      for (std::int32_t i = 0; i < field.interior.x; ++i)
        field.unchecked({i, j, k}, 0U) -= mean;
  return {};
}

std::size_t flat(Int3 shape, Int3 cell) noexcept {
  return static_cast<std::size_t>(cell.x) +
         static_cast<std::size_t>(shape.x) *
             (static_cast<std::size_t>(cell.y) +
              static_cast<std::size_t>(shape.y) *
                  static_cast<std::size_t>(cell.z));
}

std::size_t face_flat(Int3 extents, Int3 face) noexcept {
  return flat(extents, face);
}

template <class Implementation>
void copy_finest_coefficients(const Implementation& implementation,
                              double* base) noexcept {
  const detail::MgLevelStorage& level = implementation.levels[0U];
  double* const diagonal = detail::block(base, level.diagonal_offset);
  double* const x = detail::block(base, level.x_offset);
  double* const y = detail::block(base, level.y_offset);
  double* const z = detail::block(base, level.z_offset);
  double* const volumes = detail::block(base, level.volume_offset);
  const Int3 cells = level.view.local_shape;
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 cell{i, j, k};
        diagonal[flat(cells, cell)] =
            implementation.coefficients.diagonal.unchecked(cell, 0U);
        const Int3 global{level.patch.begin.x + i, level.patch.begin.y + j,
                          level.patch.begin.z + k};
        volumes[flat(cells, cell)] =
            implementation.spec.geometry->x().widths().data[
                static_cast<std::size_t>(global.x)] *
            implementation.spec.geometry->y().widths().data[
                static_cast<std::size_t>(global.y)] *
            implementation.spec.geometry->z().widths().data[
                static_cast<std::size_t>(global.z)];
      }
    }
  }
  const Int3 x_extents{cells.x + 1, cells.y, cells.z};
  const Int3 y_extents{cells.x, cells.y + 1, cells.z};
  const Int3 z_extents{cells.x, cells.y, cells.z + 1};
  for (std::int32_t k = 0; k < x_extents.z; ++k) {
    for (std::int32_t j = 0; j < x_extents.y; ++j) {
      for (std::int32_t i = 0; i < x_extents.x; ++i) {
        const Int3 face{i, j, k};
        x[face_flat(x_extents, face)] =
            implementation.coefficients.x.unchecked(face);
      }
    }
  }
  for (std::int32_t k = 0; k < y_extents.z; ++k) {
    for (std::int32_t j = 0; j < y_extents.y; ++j) {
      for (std::int32_t i = 0; i < y_extents.x; ++i) {
        const Int3 face{i, j, k};
        y[face_flat(y_extents, face)] =
            implementation.coefficients.y.unchecked(face);
      }
    }
  }
  for (std::int32_t k = 0; k < z_extents.z; ++k) {
    for (std::int32_t j = 0; j < z_extents.y; ++j) {
      for (std::int32_t i = 0; i < z_extents.x; ++i) {
        const Int3 face{i, j, k};
        z[face_flat(z_extents, face)] =
            implementation.coefficients.z.unchecked(face);
      }
    }
  }
}

double aggregated_axis_width(const AxisMetrics& metrics,
                             std::int32_t selected_index,
                             std::int32_t selected_extent,
                             std::int32_t finest_extent) noexcept {
  const std::int32_t begin = static_cast<std::int32_t>(
      (static_cast<std::int64_t>(selected_index) * finest_extent) /
      selected_extent);
  const std::int32_t end = static_cast<std::int32_t>(
      (static_cast<std::int64_t>(selected_index + 1) * finest_extent) /
      selected_extent);
  double width = 0.0;
  for (std::int32_t index = begin; index < std::max(begin + 1, end); ++index) {
    width += metrics.widths().data[static_cast<std::size_t>(
        std::min(index, finest_extent - 1))];
  }
  return width;
}

double level_cell_volume(const CartesianGeometryPlan& geometry,
                         Int3 level_shape, Int3 global_cell) noexcept {
  const Int3 finest = geometry.global_cells();
  return aggregated_axis_width(geometry.x(), global_cell.x, level_shape.x,
                               finest.x) *
         aggregated_axis_width(geometry.y(), global_cell.y, level_shape.y,
                               finest.y) *
         aggregated_axis_width(geometry.z(), global_cell.z, level_shape.z,
                               finest.z);
}

// Rediscretise the conservative face-flux operator by summing all fine faces
// that tile a coarse face.  This preserves constants and flux conservation;
// the coarse diagonal is rebuilt from those faces, including the physical
// boundary contribution already encoded in the fine diagonal.

std::int32_t axis_component(Int3 value, int axis) noexcept {
  return axis == 0 ? value.x : (axis == 1 ? value.y : value.z);
}

void set_axis_component(Int3& value, int axis,
                        std::int32_t selected) noexcept {
  if (axis == 0) {
    value.x = selected;
  } else if (axis == 1) {
    value.y = selected;
  } else {
    value.z = selected;
  }
}

Int3 tensor_extents(Int3 cells, int face_axis) noexcept {
  if (face_axis >= 0) {
    set_axis_component(cells, face_axis,
                       axis_component(cells, face_axis) + 1);
  }
  return cells;
}

std::size_t tensor_plane_size(Int3 extents, int axis) noexcept {
  if (axis == 0) {
    return static_cast<std::size_t>(extents.y) *
           static_cast<std::size_t>(extents.z);
  }
  if (axis == 1) {
    return static_cast<std::size_t>(extents.x) *
           static_cast<std::size_t>(extents.z);
  }
  return static_cast<std::size_t>(extents.x) *
         static_cast<std::size_t>(extents.y);
}

std::size_t tensor_plane_offset(Int3 extents, int axis,
                                Int3 point) noexcept {
  if (axis == 0) {
    return static_cast<std::size_t>(point.y) +
           static_cast<std::size_t>(extents.y) *
               static_cast<std::size_t>(point.z);
  }
  if (axis == 1) {
    return static_cast<std::size_t>(point.x) +
           static_cast<std::size_t>(extents.x) *
               static_cast<std::size_t>(point.z);
  }
  return static_cast<std::size_t>(point.x) +
         static_cast<std::size_t>(extents.x) *
             static_cast<std::size_t>(point.y);
}

void pack_tensor_plane(const double* source, Int3 extents, int axis,
                       std::int32_t selected, double* wire) noexcept {
  if (axis == 0) {
    for (std::int32_t k = 0; k < extents.z; ++k)
      for (std::int32_t j = 0; j < extents.y; ++j)
        wire[static_cast<std::size_t>(j) +
             static_cast<std::size_t>(extents.y) * k] =
            source[flat(extents, {selected, j, k})];
  } else if (axis == 1) {
    for (std::int32_t k = 0; k < extents.z; ++k)
      for (std::int32_t i = 0; i < extents.x; ++i)
        wire[static_cast<std::size_t>(i) +
             static_cast<std::size_t>(extents.x) * k] =
            source[flat(extents, {i, selected, k})];
  } else {
    for (std::int32_t j = 0; j < extents.y; ++j)
      for (std::int32_t i = 0; i < extents.x; ++i)
        wire[static_cast<std::size_t>(i) +
             static_cast<std::size_t>(extents.x) * j] =
            source[flat(extents, {i, j, selected})];
  }
}

int tensor_rank_stride(const MeshPatch& patch, int axis) noexcept {
  return axis == 0
             ? 1
             : (axis == 1 ? patch.process_grid.x
                          : patch.process_grid.x * patch.process_grid.y);
}

// Aggregate one cell or face channel in x->y->z order.  At an odd ownership
// boundary, the upper rank sends an already-aggregated plane to the lower
// rank.  Consequently the y plane contains completed x aggregation and the z
// plane contains completed x/y aggregation, including edges and corners.
template <class Implementation, class Loader>
Status aggregate_tensor_channel(Implementation& implementation,
                                const detail::MgLevelStorage& fine,
                                const detail::MgLevelStorage& coarse,
                                int face_axis, int channel,
                                Loader&& load, double* target) noexcept {
  const std::size_t values = implementation.coefficient_tensor_values;
  const std::size_t plane_capacity =
      implementation.coefficient_tensor_plane;
  double* const first = implementation.coefficient_tensor_storage.data();
  double* const second = first + values;
  double* const send_wire = second + values;
  double* const receive_wire = send_wire + plane_capacity;
  Int3 global = fine.view.global_shape;
  MeshPatch patch = fine.patch;
  Int3 extents = tensor_extents(patch.cells, face_axis);
  const std::size_t initial_values = cell_count(extents);
  Status local = initial_values != 0U && initial_values <= values &&
                         target != nullptr
                     ? Status{}
                     : Status{StatusCode::invalid_plan, kMgPlan};
  Status agreed = consensus(implementation, local);
  if (!agreed) return agreed;
  for (std::int32_t k = 0; k < extents.z; ++k)
    for (std::int32_t j = 0; j < extents.y; ++j)
      for (std::int32_t i = 0; i < extents.x; ++i)
        first[flat(extents, {i, j, k})] = load(Int3{i, j, k});

  double* current = first;
  double* next = second;
  for (int axis = 0; axis < 3; ++axis) {
    const std::uint8_t bit = static_cast<std::uint8_t>(1U << axis);
    if ((fine.coarsen_mask & bit) == 0U) continue;
    Int3 next_global = global;
    MeshPatch next_patch = patch;
    const std::int32_t begin = axis_component(patch.begin, axis);
    const std::int32_t end =
        begin + axis_component(patch.cells, axis);
    const std::int32_t global_extent = axis_component(global, axis);
    set_axis_component(next_global, axis, (global_extent + 1) / 2);
    set_axis_component(next_patch.begin, axis, (begin + 1) / 2);
    set_axis_component(next_patch.cells, axis,
                       (end + 1) / 2 - (begin + 1) / 2);
    const Int3 next_extents = tensor_extents(next_patch.cells, face_axis);
    const std::size_t next_values = cell_count(next_extents);
    const std::size_t plane = tensor_plane_size(extents, axis);
    const std::int32_t coordinate =
        axis_component(patch.process_coord, axis);
    const std::int32_t process_extent =
        axis_component(patch.process_grid, axis);
    const bool send = coordinate > 0 && begin > 0 && (begin & 1) != 0;
    const bool receive = coordinate + 1 < process_extent &&
                         end < global_extent && (end & 1) != 0;
    const int stride = tensor_rank_stride(patch, axis);
    const std::int64_t expected_rank =
        static_cast<std::int64_t>(patch.process_coord.x) +
        static_cast<std::int64_t>(patch.process_grid.x) *
            (static_cast<std::int64_t>(patch.process_coord.y) +
             static_cast<std::int64_t>(patch.process_grid.y) *
                 patch.process_coord.z);
    local = next_values != 0U && next_values <= values && plane != 0U &&
                    plane <= plane_capacity &&
                    plane <= static_cast<std::size_t>(
                                 std::numeric_limits<int>::max()) &&
                    expected_rank == implementation.rank
                ? Status{}
                : Status{StatusCode::invalid_plan, kMgPlan};
    agreed = consensus(implementation, local);
    if (!agreed) return agreed;

    if (send) {
      // A normal-face channel needs the aligned face one fine interval inside
      // the upper patch.  Cell and tangential-face channels donate index zero.
      const std::int32_t selected = face_axis == axis ? 1 : 0;
      pack_tensor_plane(current, extents, axis, selected, send_wire);
    }
    const int mpi_status = MPI_Sendrecv(
        send_wire, send ? static_cast<int>(plane) : 0, MPI_DOUBLE,
        send ? implementation.rank - stride : MPI_PROC_NULL,
        23000 + static_cast<int>((&coarse - implementation.levels.data()) *
                                 32U) +
            channel * 4 + axis,
        receive_wire, receive ? static_cast<int>(plane) : 0, MPI_DOUBLE,
        receive ? implementation.rank + stride : MPI_PROC_NULL,
        23000 + static_cast<int>((&coarse - implementation.levels.data()) *
                                 32U) +
            channel * 4 + axis,
        implementation.communicator, MPI_STATUS_IGNORE);
    agreed = consensus(
        implementation,
        mpi_status == MPI_SUCCESS
            ? Status{}
            : Status{StatusCode::mpi_failure, kMgCollective});
    if (!agreed) return agreed;

    bool mapping_valid = true;
    for (std::int32_t k = 0; k < next_extents.z; ++k) {
      for (std::int32_t j = 0; j < next_extents.y; ++j) {
        for (std::int32_t i = 0; i < next_extents.x; ++i) {
          const Int3 output{i, j, k};
          const std::int32_t coarse_global =
              axis_component(next_patch.begin, axis) +
              axis_component(output, axis);
          double value = 0.0;
          if (face_axis == axis) {
            const std::int32_t fine_global =
                std::min(2 * coarse_global, global_extent);
            Int3 source = output;
            set_axis_component(source, axis, fine_global - begin);
            const std::int32_t selected = axis_component(source, axis);
            if (selected >= 0 && selected < axis_component(extents, axis)) {
              value = current[flat(extents, source)];
            } else if (receive && fine_global == end + 1) {
              value = receive_wire[tensor_plane_offset(extents, axis,
                                                       output)];
            } else {
              mapping_valid = false;
            }
          } else {
            const std::int32_t first_child = 2 * coarse_global;
            for (int child = 0; child < 2; ++child) {
              const std::int32_t fine_global = first_child + child;
              if (fine_global >= global_extent) continue;
              Int3 source = output;
              set_axis_component(source, axis, fine_global - begin);
              const std::int32_t selected = axis_component(source, axis);
              if (selected >= 0 && selected < axis_component(extents, axis)) {
                value += current[flat(extents, source)];
              } else if (receive && fine_global == end) {
                value += receive_wire[tensor_plane_offset(extents, axis,
                                                          output)];
              } else {
                mapping_valid = false;
              }
            }
          }
          next[flat(next_extents, output)] = value;
        }
      }
    }
    agreed = consensus(
        implementation,
        mapping_valid ? Status{} : Status{StatusCode::invalid_plan, kMgPlan});
    if (!agreed) return agreed;
    std::swap(current, next);
    global = next_global;
    patch = next_patch;
    extents = next_extents;
  }
  const Int3 target_extents =
      tensor_extents(coarse.view.local_shape, face_axis);
  local = same_shape(global, coarse.view.global_shape) &&
                  same_shape(patch.begin, coarse.patch.begin) &&
                  same_shape(patch.cells, coarse.patch.cells) &&
                  same_shape(extents, target_extents)
              ? Status{}
              : Status{StatusCode::invalid_plan, kMgPlan};
  agreed = consensus(implementation, local);
  if (!agreed) return agreed;
  std::copy_n(current, cell_count(extents), target);
  return {};
}

// Conservative cell quantities and face transmissibilities use the same
// tensor path.  Rebuilding the diagonal last gives both ranks bit-identical
// shared-face contributions and preserves an independently aggregated
// reaction/mass term.
template <class Implementation>
Status build_coarse_coefficients(Implementation& implementation,
                                 double* base) noexcept {
  copy_finest_coefficients(implementation, base);
  for (std::size_t level_index = 1U;
       level_index < implementation.levels.size(); ++level_index) {
    const detail::MgLevelStorage& fine =
        implementation.levels[level_index - 1U];
    const detail::MgLevelStorage& coarse = implementation.levels[level_index];
    const Int3 fs = fine.view.local_shape;
    const Int3 cs = coarse.view.local_shape;
    const double* const fd = detail::block(base, fine.diagonal_offset);
    const double* const fx = detail::block(base, fine.x_offset);
    const double* const fy = detail::block(base, fine.y_offset);
    const double* const fz = detail::block(base, fine.z_offset);
    const double* const fv = detail::block(base, fine.volume_offset);
    double* const cd = detail::block(base, coarse.diagonal_offset);
    double* const cx = detail::block(base, coarse.x_offset);
    double* const cy = detail::block(base, coarse.y_offset);
    double* const cz = detail::block(base, coarse.z_offset);
    double* const cv = detail::block(base, coarse.volume_offset);
    const Int3 fxe{fs.x + 1, fs.y, fs.z};
    const Int3 fye{fs.x, fs.y + 1, fs.z};
    const Int3 fze{fs.x, fs.y, fs.z + 1};

    Status status = aggregate_tensor_channel(
        implementation, fine, coarse, -1, 0,
        [&](Int3 cell) noexcept { return fv[flat(fs, cell)]; }, cv);
    if (!status) return status;
    status = aggregate_tensor_channel(
        implementation, fine, coarse, -1, 1,
        [&](Int3 cell) noexcept {
          const double faces = fx[face_flat(fxe, cell)] +
                               fx[face_flat(fxe, {cell.x + 1, cell.y,
                                                  cell.z})] +
                               fy[face_flat(fye, cell)] +
                               fy[face_flat(fye, {cell.x, cell.y + 1,
                                                  cell.z})] +
                               fz[face_flat(fze, cell)] +
                               fz[face_flat(fze, {cell.x, cell.y,
                                                  cell.z + 1})];
          return std::max(0.0, fd[flat(fs, cell)] - faces);
        },
        cd);
    if (!status) return status;
    status = aggregate_tensor_channel(
        implementation, fine, coarse, 0, 2,
        [&](Int3 face) noexcept { return fx[face_flat(fxe, face)]; }, cx);
    if (!status) return status;
    status = aggregate_tensor_channel(
        implementation, fine, coarse, 1, 3,
        [&](Int3 face) noexcept { return fy[face_flat(fye, face)]; }, cy);
    if (!status) return status;
    status = aggregate_tensor_channel(
        implementation, fine, coarse, 2, 4,
        [&](Int3 face) noexcept { return fz[face_flat(fze, face)]; }, cz);
    if (!status) return status;

    const Int3 cxe{cs.x + 1, cs.y, cs.z};
    const Int3 cye{cs.x, cs.y + 1, cs.z};
    const Int3 cze{cs.x, cs.y, cs.z + 1};
    Status local{};
    for (std::int32_t k = 0; k < cs.z; ++k) {
      for (std::int32_t j = 0; j < cs.y; ++j) {
        for (std::int32_t i = 0; i < cs.x; ++i) {
          const Int3 cell{i, j, k};
          const std::size_t index = flat(cs, cell);
          cd[index] += cx[face_flat(cxe, cell)] +
                       cx[face_flat(cxe, {i + 1, j, k})] +
                       cy[face_flat(cye, cell)] +
                       cy[face_flat(cye, {i, j + 1, k})] +
                       cz[face_flat(cze, cell)] +
                       cz[face_flat(cze, {i, j, k + 1})];
          if (!(cv[index] > 0.0) || !std::isfinite(cv[index]) ||
              !(cd[index] > 0.0) || !std::isfinite(cd[index])) {
            local = {StatusCode::numerical_failure, kMgCoefficient};
          }
        }
      }
    }
    status = consensus(implementation, local);
    if (!status) return status;
  }
  return {};
}

template <class Implementation>
HaloEngine* halo_for(Implementation& implementation,
                     std::size_t level) noexcept {
  return level == 0U ? implementation.services.finest_halo
                     : implementation.services.level_halos.data[level - 1U];
}

template <class Implementation>
Status exchange_solution(Implementation& implementation, std::size_t level,
                         FieldView& solution, StageId stage) noexcept {
  Status local = implementation.services.workspace->revise_level(
      level, MgWorkspaceSlot::solution);
  solution = implementation.services.workspace->level(
      level, MgWorkspaceSlot::solution);
  if (!local) return local;
  HaloTicket ticket;
  std::array<FieldView, 1U> views{solution};
  local = halo_for(implementation, level)->begin(
      stage, {views.data(), views.size()}, ticket);
  if (local) {
    local = halo_for(implementation, level)->finish(
        ticket, {views.data(), views.size()});
    solution = views[0U];
  }
  return local;
}

bool exterior(const MeshPatch& patch, Int3 global, CartesianAxis axis,
              bool minimum) noexcept {
  const std::int32_t begin = axis == CartesianAxis::x
                                 ? patch.begin.x
                                 : (axis == CartesianAxis::y ? patch.begin.y
                                                             : patch.begin.z);
  const std::int32_t count = axis == CartesianAxis::x
                                 ? patch.cells.x
                                 : (axis == CartesianAxis::y ? patch.cells.y
                                                             : patch.cells.z);
  const std::int32_t extent = axis == CartesianAxis::x
                                  ? global.x
                                  : (axis == CartesianAxis::y ? global.y
                                                              : global.z);
  return minimum ? begin == 0 : begin + count == extent;
}

double level_neighbor(ConstFieldView field, const MeshPatch& patch,
                      Int3 global, Int3 cell, CartesianAxis axis,
                      int direction, MgBoundaryKind minimum,
                      MgBoundaryKind maximum) noexcept {
  Int3 selected = cell;
  std::int32_t* index = axis == CartesianAxis::x
                            ? &selected.x
                            : (axis == CartesianAxis::y ? &selected.y
                                                        : &selected.z);
  const std::int32_t extent = axis == CartesianAxis::x
                                  ? field.interior.x
                                  : (axis == CartesianAxis::y
                                         ? field.interior.y
                                         : field.interior.z);
  *index += direction;
  if (*index >= 0 && *index < extent) {
    return field.unchecked(selected, 0U);
  }
  const bool at_global = exterior(patch, global, axis, direction < 0);
  if (!at_global) {
    return field.unchecked(selected, 0U);
  }
  const MgBoundaryKind boundary = direction < 0 ? minimum : maximum;
  if (boundary == MgBoundaryKind::periodic) {
    // Periodic topology is part of every reserved level halo, including the
    // single-rank local-copy case.  Reading that ghost keeps the operator
    // identical when a periodic axis is repartitioned.
    return field.unchecked(selected, 0U);
  }
  if (boundary == MgBoundaryKind::neumann) {
    *index = direction < 0 ? 0 : extent - 1;
    return field.unchecked(selected, 0U);
  }
  return 0.0;
}

template <class Implementation>
void apply_level_matrix(const Implementation& implementation,
                        std::size_t level_index, ConstFieldView input,
                        FieldView output) noexcept {
  const detail::MgLevelStorage& level = implementation.levels[level_index];
  const double* const base = implementation.hierarchy_storage.data();
  const double* const diagonal = detail::block(base, level.diagonal_offset);
  const double* const x = detail::block(base, level.x_offset);
  const double* const y = detail::block(base, level.y_offset);
  const double* const z = detail::block(base, level.z_offset);
  const Int3 cells = level.view.local_shape;
  const Int3 xe{cells.x + 1, cells.y, cells.z};
  const Int3 ye{cells.x, cells.y + 1, cells.z};
  const Int3 ze{cells.x, cells.y, cells.z + 1};
  const MgBoundarySet boundary = implementation.spec.boundaries;
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 cell{i, j, k};
        double value = diagonal[flat(cells, cell)] *
                       input.unchecked(cell, 0U);
        value -= x[face_flat(xe, {i, j, k})] *
                 level_neighbor(input, level.patch, level.view.global_shape,
                                cell, CartesianAxis::x, -1, boundary.x_min,
                                boundary.x_max);
        value -= x[face_flat(xe, {i + 1, j, k})] *
                 level_neighbor(input, level.patch, level.view.global_shape,
                                cell, CartesianAxis::x, 1, boundary.x_min,
                                boundary.x_max);
        value -= y[face_flat(ye, {i, j, k})] *
                 level_neighbor(input, level.patch, level.view.global_shape,
                                cell, CartesianAxis::y, -1, boundary.y_min,
                                boundary.y_max);
        value -= y[face_flat(ye, {i, j + 1, k})] *
                 level_neighbor(input, level.patch, level.view.global_shape,
                                cell, CartesianAxis::y, 1, boundary.y_min,
                                boundary.y_max);
        value -= z[face_flat(ze, {i, j, k})] *
                 level_neighbor(input, level.patch, level.view.global_shape,
                                cell, CartesianAxis::z, -1, boundary.z_min,
                                boundary.z_max);
        value -= z[face_flat(ze, {i, j, k + 1})] *
                 level_neighbor(input, level.patch, level.view.global_shape,
                                cell, CartesianAxis::z, 1, boundary.z_min,
                                boundary.z_max);
        output.unchecked(cell, 0U) = value;
      }
    }
  }
}

template <class Implementation>
Status point_smooth(Implementation& implementation, std::size_t level_index,
                    std::uint32_t sweeps, bool reverse,
                    StageId& stage) noexcept {
  FieldView x = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::solution);
  FieldView rhs = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::rhs);
  FieldView residual = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::residual);
  const detail::MgLevelStorage& level = implementation.levels[level_index];
  const double* const diagonal = detail::block(
      implementation.hierarchy_storage.data(), level.diagonal_offset);
  const Int3 cells = level.view.local_shape;
  constexpr double omega = 0.72;
  for (std::uint32_t sweep = 0U; sweep < sweeps; ++sweep) {
    Status status = exchange_solution(implementation, level_index, x,
                                      stage++);
    if (!status) return status;
    apply_level_matrix(implementation, level_index, as_const(x), residual);
    const int first = reverse ? 1 : 0;
    for (int pass = 0; pass < 2; ++pass) {
      const int color = (first + pass) & 1;
      if (reverse) {
        for (std::int32_t k = cells.z; k-- > 0;) {
          for (std::int32_t j = cells.y; j-- > 0;) {
            for (std::int32_t i = cells.x; i-- > 0;) {
              const Int3 cell{i, j, k};
              if (((i + j + k) & 1) == color) {
                x.unchecked(cell, 0U) +=
                    omega * (rhs.unchecked(cell, 0U) -
                             residual.unchecked(cell, 0U)) /
                    diagonal[flat(cells, cell)];
              }
            }
          }
        }
      } else {
        for (std::int32_t k = 0; k < cells.z; ++k) {
          for (std::int32_t j = 0; j < cells.y; ++j) {
            for (std::int32_t i = 0; i < cells.x; ++i) {
              const Int3 cell{i, j, k};
              if (((i + j + k) & 1) == color) {
                x.unchecked(cell, 0U) +=
                    omega * (rhs.unchecked(cell, 0U) -
                             residual.unchecked(cell, 0U)) /
                    diagonal[flat(cells, cell)];
              }
            }
          }
        }
      }
      status = exchange_solution(implementation, level_index, x, stage++);
      if (!status) return status;
      apply_level_matrix(implementation, level_index, as_const(x), residual);
    }
  }
  Status status = implementation.services.workspace->revise_level(
      level_index, MgWorkspaceSlot::residual);
  return status;
}

template <class Implementation>
Status line_smooth(Implementation& implementation, std::size_t level_index,
                   std::uint32_t sweeps, bool reverse,
                   StageId& stage,
                   Status& deferred,
                   std::uint8_t selected_mask = 0U) noexcept {
  const detail::MgLevelStorage& level = implementation.levels[level_index];
  const std::uint8_t mask = selected_mask == 0U
                                ? level.view.line_axis_mask
                                : selected_mask;
  if (mask == 0U || (mask & (mask - 1U)) != 0U) {
    return point_smooth(implementation, level_index, sweeps, reverse, stage);
  }
  const CartesianAxis axis =
      mask == detail::kMgAxisX
          ? CartesianAxis::x
          : (mask == detail::kMgAxisY ? CartesianAxis::y : CartesianAxis::z);
  const Int3 cells = level.view.local_shape;
  const std::int32_t extent = axis == CartesianAxis::x
                                  ? cells.x
                                  : (axis == CartesianAxis::y ? cells.y
                                                              : cells.z);
  if (extent < 2 || extent > static_cast<std::int32_t>(
                                 implementation.spec.policy
                                     .line_relaxation_maximum_extent)) {
    return point_smooth(implementation, level_index, sweeps, reverse, stage);
  }
  FieldView x = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::solution);
  FieldView rhs = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::rhs);
  FieldView work = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::temporary);
  const double* const base = implementation.hierarchy_storage.data();
  const double* const diagonal = detail::block(base, level.diagonal_offset);
  const double* const xf = detail::block(base, level.x_offset);
  const double* const yf = detail::block(base, level.y_offset);
  const double* const zf = detail::block(base, level.z_offset);
  const Int3 xe{cells.x + 1, cells.y, cells.z};
  const Int3 ye{cells.x, cells.y + 1, cells.z};
  const Int3 ze{cells.x, cells.y, cells.z + 1};
  const MgBoundarySet boundary = implementation.spec.boundaries;
  const auto solve_one = [&](std::int32_t fixed_a,
                             std::int32_t fixed_b) noexcept {
    auto cell_at = [&](std::int32_t along) noexcept {
      return axis == CartesianAxis::x
                 ? Int3{along, fixed_a, fixed_b}
                 : (axis == CartesianAxis::y
                        ? Int3{fixed_a, along, fixed_b}
                        : Int3{fixed_a, fixed_b, along});
    };
    auto low_face = [&](Int3 cell) noexcept {
      return axis == CartesianAxis::x
                 ? xf[face_flat(xe, cell)]
                 : (axis == CartesianAxis::y
                        ? yf[face_flat(ye, cell)]
                        : zf[face_flat(ze, cell)]);
    };
    auto high_face = [&](Int3 cell) noexcept {
      return axis == CartesianAxis::x
                 ? xf[face_flat(xe, {cell.x + 1, cell.y, cell.z})]
                 : (axis == CartesianAxis::y
                        ? yf[face_flat(ye, {cell.x, cell.y + 1, cell.z})]
                        : zf[face_flat(ze, {cell.x, cell.y, cell.z + 1})]);
    };
    auto transverse_rhs = [&](Int3 c) noexcept {
      double value = rhs.unchecked(c, 0U);
      if (axis != CartesianAxis::x) {
        value += xf[face_flat(xe, c)] *
                 level_neighbor(as_const(x), level.patch,
                                level.view.global_shape, c, CartesianAxis::x,
                                -1, boundary.x_min, boundary.x_max);
        value += xf[face_flat(xe, {c.x + 1, c.y, c.z})] *
                 level_neighbor(as_const(x), level.patch,
                                level.view.global_shape, c, CartesianAxis::x,
                                1, boundary.x_min, boundary.x_max);
      }
      if (axis != CartesianAxis::y) {
        value += yf[face_flat(ye, c)] *
                 level_neighbor(as_const(x), level.patch,
                                level.view.global_shape, c, CartesianAxis::y,
                                -1, boundary.y_min, boundary.y_max);
        value += yf[face_flat(ye, {c.x, c.y + 1, c.z})] *
                 level_neighbor(as_const(x), level.patch,
                                level.view.global_shape, c, CartesianAxis::y,
                                1, boundary.y_min, boundary.y_max);
      }
      if (axis != CartesianAxis::z) {
        value += zf[face_flat(ze, c)] *
                 level_neighbor(as_const(x), level.patch,
                                level.view.global_shape, c, CartesianAxis::z,
                                -1, boundary.z_min, boundary.z_max);
        value += zf[face_flat(ze, {c.x, c.y, c.z + 1})] *
                 level_neighbor(as_const(x), level.patch,
                                level.view.global_shape, c, CartesianAxis::z,
                                1, boundary.z_min, boundary.z_max);
      }
      return value;
    };
    // Thomas forward pass. temporary stores c-prime while residual stores
    // d-prime; both arrays are persistent MG workspace slots.
    FieldView dprime = implementation.services.workspace->level(
        level_index, MgWorkspaceSlot::residual);
    for (std::int32_t n = 0; n < extent; ++n) {
      const Int3 c = cell_at(n);
      double a = n == 0 ? 0.0 : -low_face(c);
      double b = diagonal[flat(cells, c)];
      const double cc = n + 1 == extent ? 0.0 : -high_face(c);
      double d = transverse_rhs(c);
      if (n == 0) {
        d += low_face(c) * level_neighbor(
                                    as_const(x), level.patch,
                                    level.view.global_shape, c, axis, -1,
                                    axis == CartesianAxis::x
                                        ? boundary.x_min
                                        : (axis == CartesianAxis::y
                                               ? boundary.y_min
                                               : boundary.z_min),
                                    axis == CartesianAxis::x
                                        ? boundary.x_max
                                        : (axis == CartesianAxis::y
                                               ? boundary.y_max
                                               : boundary.z_max));
      }
      if (n + 1 == extent) {
        d += high_face(c) * level_neighbor(
                                     as_const(x), level.patch,
                                     level.view.global_shape, c, axis, 1,
                                     axis == CartesianAxis::x
                                         ? boundary.x_min
                                         : (axis == CartesianAxis::y
                                                ? boundary.y_min
                                                : boundary.z_min),
                                     axis == CartesianAxis::x
                                         ? boundary.x_max
                                         : (axis == CartesianAxis::y
                                                ? boundary.y_max
                                                : boundary.z_max));
      }
      if (n > 0) {
        const Int3 previous = cell_at(n - 1);
        const double denominator = b - a * work.unchecked(previous, 0U);
        if (!(std::abs(denominator) >
              std::numeric_limits<double>::min())) return false;
        work.unchecked(c, 0U) = cc / denominator;
        dprime.unchecked(c, 0U) =
            (d - a * dprime.unchecked(previous, 0U)) / denominator;
      } else {
        if (!(std::abs(b) > std::numeric_limits<double>::min())) return false;
        work.unchecked(c, 0U) = cc / b;
        dprime.unchecked(c, 0U) = d / b;
      }
    }
    for (std::int32_t n = extent; n-- > 0;) {
      const Int3 c = cell_at(n);
      const double value = n + 1 == extent
                               ? dprime.unchecked(c, 0U)
                               : dprime.unchecked(c, 0U) -
                                     work.unchecked(c, 0U) *
                                         x.unchecked(cell_at(n + 1), 0U);
      x.unchecked(c, 0U) = value;
    }
    return true;
  };
  for (std::uint32_t sweep = 0U; sweep < sweeps; ++sweep) {
    Status status = exchange_solution(implementation, level_index, x,
                                      stage++);
    if (!status) return status;
    bool okay = true;
    const std::int32_t a_extent = axis == CartesianAxis::x ? cells.y : cells.x;
    const std::int32_t b_extent = axis == CartesianAxis::z ? cells.y : cells.z;
    if (reverse) {
      for (std::int32_t b = b_extent; b-- > 0;)
        for (std::int32_t a = a_extent; a-- > 0;) okay = solve_one(a, b) && okay;
    } else {
      for (std::int32_t b = 0; b < b_extent; ++b)
        for (std::int32_t a = 0; a < a_extent; ++a) okay = solve_one(a, b) && okay;
    }
    if (!okay && deferred) {
      deferred = {StatusCode::numerical_failure, kMgApply};
    }
  }
  Status status = implementation.services.workspace->revise_level(
      level_index, MgWorkspaceSlot::solution);
  if (status) status = implementation.services.workspace->revise_level(
      level_index, MgWorkspaceSlot::temporary);
  if (status) status = implementation.services.workspace->revise_level(
      level_index, MgWorkspaceSlot::residual);
  return status;
}

template <class Implementation>
Status smooth(Implementation& implementation, std::size_t level,
              std::uint32_t sweeps, bool reverse, StageId& stage,
              Status& deferred) noexcept {
  const std::uint8_t mask = implementation.levels[level].view.line_axis_mask;
  if (mask == 0U) {
    return point_smooth(implementation, level, sweeps, reverse, stage);
  }
  const std::uint8_t axes[3]{detail::kMgAxisX, detail::kMgAxisY,
                             detail::kMgAxisZ};
  if (reverse) {
    for (std::size_t index = 3U; index-- > 0U;) {
      if ((mask & axes[index]) != 0U) {
        const Status status = line_smooth(implementation, level, sweeps,
                                          true, stage, deferred, axes[index]);
        if (!status) return status;
      }
    }
  } else {
    for (const std::uint8_t axis : axes) {
      if ((mask & axis) != 0U) {
        const Status status = line_smooth(implementation, level, sweeps,
                                          false, stage, deferred, axis);
        if (!status) return status;
      }
    }
  }
  return {};
}

template <class Implementation>
Status compute_residual(Implementation& implementation,
                        std::size_t level_index, StageId& stage) noexcept {
  FieldView x = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::solution);
  FieldView rhs = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::rhs);
  FieldView residual = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::residual);
  Status status = exchange_solution(implementation, level_index, x, stage++);
  if (!status) return status;
  apply_level_matrix(implementation, level_index, as_const(x), residual);
  const Int3 cells = x.interior;
  for (std::int32_t k = 0; k < cells.z; ++k)
    for (std::int32_t j = 0; j < cells.y; ++j)
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 cell{i, j, k};
        residual.unchecked(cell, 0U) = rhs.unchecked(cell, 0U) -
                                       residual.unchecked(cell, 0U);
      }
  return implementation.services.workspace->revise_level(
      level_index, MgWorkspaceSlot::residual);
}

double ghost_or_local(ConstFieldView field, Int3 local) noexcept {
  return field.unchecked(local, 0U);
}

template <class Implementation>
Status restrict_residual(Implementation& implementation,
                         std::size_t fine_index, StageId& stage) noexcept {
  const detail::MgLevelStorage& fine = implementation.levels[fine_index];
  const detail::MgLevelStorage& coarse = implementation.levels[fine_index + 1U];
  FieldView fr = implementation.services.workspace->level(
      fine_index, MgWorkspaceSlot::residual);
  // A transfer may need one donor owned by a neighbouring fine patch when a
  // ceil-mapped partition boundary is odd.
  Status status = implementation.services.workspace->revise_level(
      fine_index, MgWorkspaceSlot::residual);
  fr = implementation.services.workspace->level(
      fine_index, MgWorkspaceSlot::residual);
  if (!status) return status;
  HaloTicket ticket;
  std::array<FieldView, 1U> views{fr};
  status = halo_for(implementation, fine_index)->begin(
      stage++, {views.data(), views.size()}, ticket);
  if (status) {
    status = halo_for(implementation, fine_index)->finish(
        ticket, {views.data(), views.size()});
    fr = views[0U];
  }
  if (!status) return status;
  FieldView crhs = implementation.services.workspace->level(
      fine_index + 1U, MgWorkspaceSlot::rhs);
  FieldView cx = implementation.services.workspace->level(
      fine_index + 1U, MgWorkspaceSlot::solution);
  const std::uint8_t mask = fine.coarsen_mask;
  const bool mx = (mask & detail::kMgAxisX) != 0U;
  const bool my = (mask & detail::kMgAxisY) != 0U;
  const bool mz = (mask & detail::kMgAxisZ) != 0U;
  const Int3 cs = coarse.view.local_shape;
  const double* const fine_volumes = detail::block(
      implementation.hierarchy_storage.data(), fine.volume_offset);
  for (std::int32_t ck = 0; ck < cs.z; ++ck) {
    for (std::int32_t cj = 0; cj < cs.y; ++cj) {
      for (std::int32_t ci = 0; ci < cs.x; ++ci) {
        const Int3 cg{coarse.patch.begin.x + ci,
                      coarse.patch.begin.y + cj,
                      coarse.patch.begin.z + ck};
        const std::int32_t xb = mx ? 2 * cg.x : cg.x;
        const std::int32_t yb = my ? 2 * cg.y : cg.y;
        const std::int32_t zb = mz ? 2 * cg.z : cg.z;
        double weighted_sum = 0.0;
        double volume_sum = 0.0;
        for (int dz = 0; dz < (mz ? 2 : 1); ++dz)
          for (int dy = 0; dy < (my ? 2 : 1); ++dy)
            for (int dx = 0; dx < (mx ? 2 : 1); ++dx) {
              const Int3 fg{xb + dx, yb + dy, zb + dz};
              if (fg.x >= fine.view.global_shape.x ||
                  fg.y >= fine.view.global_shape.y ||
                  fg.z >= fine.view.global_shape.z) continue;
              const Int3 local{fg.x - fine.patch.begin.x,
                               fg.y - fine.patch.begin.y,
                               fg.z - fine.patch.begin.z};
              const bool owned = local.x >= 0 && local.x < fine.view.local_shape.x &&
                                 local.y >= 0 && local.y < fine.view.local_shape.y &&
                                 local.z >= 0 && local.z < fine.view.local_shape.z;
              // The at-most-one ghost donor at an odd partition boundary has
              // no separate volume halo. Cartesian geometry is replicated,
              // so obtain its exact volume from global metrics instead.
              const double volume =
                  owned
                      ? fine_volumes[flat(fine.view.local_shape, local)]
                      : level_cell_volume(*implementation.spec.geometry,
                                          fine.view.global_shape, fg);
              weighted_sum += volume * ghost_or_local(as_const(fr), local);
              volume_sum += volume;
            }
        crhs.unchecked({ci, cj, ck}, 0U) =
            volume_sum == 0.0 ? 0.0 : weighted_sum / volume_sum;
        cx.unchecked({ci, cj, ck}, 0U) = 0.0;
      }
    }
  }
  status = implementation.services.workspace->revise_level(
      fine_index + 1U, MgWorkspaceSlot::rhs);
  if (status && implementation.spec.null_space == MgNullSpace::constant) {
    crhs = implementation.services.workspace->level(
        fine_index + 1U, MgWorkspaceSlot::rhs);
    status = project_constant(implementation, fine_index + 1U, crhs);
    if (status) {
      status = implementation.services.workspace->revise_level(
          fine_index + 1U, MgWorkspaceSlot::rhs);
    }
  }
  if (status) status = implementation.services.workspace->revise_level(
      fine_index + 1U, MgWorkspaceSlot::solution);
  return status;
}

template <class Implementation>
Status prolongate_add(Implementation& implementation,
                      std::size_t fine_index, StageId& stage) noexcept {
  const detail::MgLevelStorage& fine = implementation.levels[fine_index];
  const detail::MgLevelStorage& coarse = implementation.levels[fine_index + 1U];
  FieldView fx = implementation.services.workspace->level(
      fine_index, MgWorkspaceSlot::solution);
  FieldView coarse_solution = implementation.services.workspace->level(
      fine_index + 1U, MgWorkspaceSlot::solution);
  Status status = exchange_solution(implementation, fine_index + 1U,
                                    coarse_solution, stage++);
  if (!status) return status;
  const ConstFieldView cx = as_const(coarse_solution);
  const bool mx = (fine.coarsen_mask & detail::kMgAxisX) != 0U;
  const bool my = (fine.coarsen_mask & detail::kMgAxisY) != 0U;
  const bool mz = (fine.coarsen_mask & detail::kMgAxisZ) != 0U;
  const Int3 fs = fine.view.local_shape;
  for (std::int32_t k = 0; k < fs.z; ++k)
    for (std::int32_t j = 0; j < fs.y; ++j)
      for (std::int32_t i = 0; i < fs.x; ++i) {
        const Int3 fg{fine.patch.begin.x + i, fine.patch.begin.y + j,
                      fine.patch.begin.z + k};
        const auto interpolation = [](std::int32_t fine_global,
                                      bool coarsened,
                                      std::int32_t coarse_extent,
                                      std::int32_t& left,
                                      std::int32_t& right,
                                      double& right_weight) noexcept {
          if (!coarsened) {
            left = fine_global;
            right = fine_global;
            right_weight = 0.0;
            return;
          }
          const std::int32_t parent = fine_global / 2;
          if ((fine_global & 1) == 0) {
            left = parent - 1;
            right = parent;
            right_weight = 0.75;
          } else {
            left = parent;
            right = parent + 1;
            right_weight = 0.25;
          }
          if (left < 0) {
            left = 0;
            right = 0;
            right_weight = 0.0;
          } else if (right >= coarse_extent) {
            left = coarse_extent - 1;
            right = coarse_extent - 1;
            right_weight = 0.0;
          }
        };
        std::int32_t xl = 0, xr = 0, yl = 0, yr = 0, zl = 0, zr = 0;
        double xw = 0.0, yw = 0.0, zw = 0.0;
        interpolation(fg.x, mx, coarse.view.global_shape.x, xl, xr, xw);
        interpolation(fg.y, my, coarse.view.global_shape.y, yl, yr, yw);
        interpolation(fg.z, mz, coarse.view.global_shape.z, zl, zr, zw);
        double value = 0.0;
        for (int az = 0; az < (zl == zr ? 1 : 2); ++az) {
          const std::int32_t gz = az == 0 ? zl : zr;
          const double wz = az == 0 ? 1.0 - zw : zw;
          for (int ay = 0; ay < (yl == yr ? 1 : 2); ++ay) {
            const std::int32_t gy = ay == 0 ? yl : yr;
            const double wy = ay == 0 ? 1.0 - yw : yw;
            for (int ax = 0; ax < (xl == xr ? 1 : 2); ++ax) {
              const std::int32_t gx = ax == 0 ? xl : xr;
              const double wx = ax == 0 ? 1.0 - xw : xw;
              const Int3 local{gx - coarse.patch.begin.x,
                               gy - coarse.patch.begin.y,
                               gz - coarse.patch.begin.z};
              value += wx * wy * wz * cx.unchecked(local, 0U);
            }
          }
        }
        fx.unchecked({i, j, k}, 0U) += value;
      }
  return {};
}

template <class Implementation>
Status v_cycle(Implementation& implementation, std::size_t level,
               StageId& stage, Status& deferred) noexcept {
  const bool coarse = level + 1U == implementation.levels.size();
  if (coarse) {
    return smooth(implementation, level,
                  implementation.spec.policy.coarse_sweeps, false, stage,
                  deferred);
  }
  Status status = smooth(implementation, level,
                         implementation.spec.policy.pre_sweeps, false, stage,
                         deferred);
  if (status) status = compute_residual(implementation, level, stage);
  if (status) status = restrict_residual(implementation, level, stage);
  if (status) status = v_cycle(implementation, level + 1U, stage, deferred);
  if (!status) return status;
  status = prolongate_add(implementation, level, stage);
  if (status) status = implementation.services.workspace->revise_level(
      level, MgWorkspaceSlot::solution);
  if (status) {
    status = smooth(implementation, level,
                    implementation.spec.policy.post_sweeps, true, stage,
                    deferred);
  }
  return status;
}

}  // namespace

NativeCartesianMgPlan::~NativeCartesianMgPlan() noexcept { release(); }

NativeCartesianMgPlan::NativeCartesianMgPlan(
    NativeCartesianMgPlan&& other) noexcept
    : implementation_(std::exchange(other.implementation_, nullptr)) {}

NativeCartesianMgPlan& NativeCartesianMgPlan::operator=(
    NativeCartesianMgPlan&& other) noexcept {
  if (this != &other) {
    release();
    implementation_ = std::exchange(other.implementation_, nullptr);
  }
  return *this;
}

void NativeCartesianMgPlan::release() noexcept {
  destroy(std::exchange(implementation_, nullptr));
}

Status NativeCartesianMgPlan::compile(const NativeCartesianMgSpec& spec,
                                      MgRuntimeServices services,
                                      MgCoefficientViews coefficients,
                                      NativeCartesianMgPlan& out,
                                      MgPlanCounters* counters) noexcept {
  Strategy strategy{};
  Status local = validate_compile(spec, services, coefficients, strategy);
  if (local) {
    local = validate_runtime_communicators(spec, services);
  }
  int raw_lowest = -1;
  Status agreed = raw_consensus(spec.communicator, local, &raw_lowest);
  if (!agreed) {
    if (out.implementation_ != nullptr) {
      out.implementation_->lowest = raw_lowest;
    }
    return agreed;
  }
  agreed = services.reductions->consensus(local);
  if (!agreed) {
    if (out.implementation_ != nullptr) {
      out.implementation_->lowest =
          services.reductions->lowest_failing_rank();
    }
    return agreed;
  }
  if (local) {
    const MgWorkspaceRequirements& required =
        services.workspace->requirements();
    local = services.level_halos.size + 1U != required.level_count
                ? Status{StatusCode::invalid_plan, kMgPlan}
                : Status{};
    if (local) {
      for (std::size_t level = 0U; level < services.level_halos.size;
           ++level) {
        if (services.level_halos.data[level] == nullptr ||
            !services.level_halos.data[level]->ready()) {
          local = {StatusCode::invalid_plan, kMgPlan};
          break;
        }
      }
    }
    if (local) {
      const MgWorkspaceRequirements& required =
          services.workspace->requirements();
      local = required.execution_revision == 0U ||
                      required.collective_fingerprint == 0U ||
                      required.fingerprint == 0U
                  ? Status{StatusCode::invalid_plan, kMgPlan}
                  : Status{};
    }
  }
  agreed = services.reductions->consensus(local);
  if (!agreed) {
    return agreed;
  }
  const MgWorkspaceLevelRequirements& required_finest =
      services.workspace->requirements().levels[0U];
  strategy.coarsening = required_finest.coarsening;
  strategy.coarsen_mask = required_finest.coarsen_axis_mask;
  strategy.line_mask = required_finest.line_axis_mask;
  const PlanFingerprint structural = structural_contract(spec, strategy);
  agreed = services.reductions->consensus_contract(structural);
  if (!agreed) {
    if (out.implementation_ != nullptr) {
      out.implementation_->lowest =
          services.reductions->lowest_failing_rank();
    }
    return agreed;
  }
  std::uint64_t next_symbolic = 0U;
  std::uint64_t next_rebuild = 0U;
  local = counters != nullptr &&
                  (!increment(counters->symbolic_builds, next_symbolic) ||
                   !increment(counters->hierarchy_rebuilds, next_rebuild))
              ? Status{StatusCode::invalid_plan, kMgCounter}
              : Status{};
  agreed = services.reductions->consensus(local);
  if (!agreed) {
    return agreed;
  }
  Impl* candidate = new (std::nothrow) Impl;
  local = candidate == nullptr ? Status{StatusCode::allocation_failure, 0U}
                               : Status{};
  agreed = services.reductions->consensus(local);
  if (!agreed) {
    delete candidate;
    return agreed;
  }
  candidate->spec = spec;
  candidate->services = services;
  candidate->coefficients = coefficients;
  candidate->strategy = strategy;
  candidate->rank = 0;
  candidate->size = 1;
  try {
    candidate->level_halo_table.assign(
        services.level_halos.data,
        services.level_halos.data + services.level_halos.size);
    candidate->level_halo_identities.reserve(
        candidate->level_halo_table.size());
    for (const HaloEngine* halo : candidate->level_halo_table) {
      candidate->level_halo_identities.push_back(
          halo == nullptr ? 0U : halo->instance_identity());
    }
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, 0U};
  } catch (...) {
    local = {StatusCode::invalid_plan, kMgPlan};
  }
  agreed = services.reductions->consensus(local);
  if (!agreed) {
    destroy(candidate);
    return agreed;
  }
  candidate->services.level_halos = {
      candidate->level_halo_table.data(), candidate->level_halo_table.size()};
  candidate->finest_halo_object = services.finest_halo;
  candidate->reductions_object = services.reductions;
  candidate->workspace_object = services.workspace;
  candidate->finest_halo_identity = services.finest_halo->instance_identity();
  candidate->reductions_identity = services.reductions->instance_identity();
  candidate->workspace_binding_identity =
      services.workspace->binding_identity();
  candidate->workspace_fingerprint = services.workspace->fingerprint();
  candidate->workspace_collective_fingerprint =
      services.workspace->collective_fingerprint();
  candidate->workspace_storage = services.workspace->storage_address();
  candidate->workspace_storage_doubles =
      services.workspace->storage_doubles();
  const int rank_status =
      MPI_Comm_rank(spec.communicator, &candidate->rank);
  const int size_status =
      MPI_Comm_size(spec.communicator, &candidate->size);
  const int duplicate_status =
      MPI_Comm_dup(spec.communicator, &candidate->communicator);
  const int error_handler_status =
      duplicate_status == MPI_SUCCESS
          ? MPI_Comm_set_errhandler(candidate->communicator, MPI_ERRORS_RETURN)
          : MPI_SUCCESS;
  local = rank_status == MPI_SUCCESS && size_status == MPI_SUCCESS &&
                  duplicate_status == MPI_SUCCESS &&
                  error_handler_status == MPI_SUCCESS
              ? Status{}
              : Status{StatusCode::mpi_failure, kMgCollective};
  agreed = services.reductions->consensus(local);
  if (!agreed) {
    destroy(candidate);
    return agreed;
  }
  try {
    local = build_levels(*candidate);
    if (local) {
      local = build_coarse_coefficients(
          *candidate, candidate->hierarchy_storage.data());
    }
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, 0U};
  } catch (...) {
    local = {StatusCode::invalid_plan, kMgPlan};
  }
  agreed = services.reductions->consensus(local);
  if (!agreed) {
    destroy(candidate);
    return agreed;
  }
  candidate->generation = 1U;
  candidate->symbolic = public_symbolic_fingerprint(spec, strategy);
  candidate->numeric = make_numeric_fingerprint(spec.coefficients,
                                                candidate->symbolic);
  candidate->hierarchy = make_hierarchy_fingerprint(
      candidate->symbolic, spec.coefficients, candidate->generation);
  candidate->runtime_counters.symbolic_builds = 1U;
  candidate->runtime_counters.hierarchy_rebuilds = 1U;
  candidate->lowest = -1;
  out.release();
  out.implementation_ = candidate;
  if (counters != nullptr) {
    counters->symbolic_builds = next_symbolic;
    counters->hierarchy_rebuilds = next_rebuild;
  }
  return {};
}

Status NativeCartesianMgPlan::update_coefficients(
    LinearIdentity next_identity, MgCoefficientIdentity identity,
    MgCoefficientViews coefficients,
    MgPlanCounters* counters) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kMgPlan};
  }
  Impl& implementation = *implementation_;
  Status local = validate_borrowed_services(implementation);
  Status agreed = implementation.services.reductions->consensus(local);
  implementation.lowest =
      implementation.services.reductions->lowest_failing_rank();
  if (!agreed) {
    return agreed;
  }
  std::uint64_t contract = UINT64_C(1469598103934665603);
  contract = mix(contract, next_identity.symbolic);
  contract = mix(contract, next_identity.numeric);
  contract = mix(contract, next_identity.hierarchy);
  contract = mix(contract, next_identity.workspace);
  contract = mix(contract, next_identity.fingerprint);
  contract = mix(contract, identity.revision);
  contract = mix(contract, identity.fingerprint);
  agreed = implementation.services.reductions->consensus_contract(
      finish(contract));
  implementation.lowest =
      implementation.services.reductions->lowest_failing_rank();
  if (!agreed) {
    return agreed;
  }
  double global_change = 0.0;
  double local_change = identity.maximum_relative_change;
  local = !valid_identity(next_identity) ||
                         next_identity.symbolic !=
                             implementation.spec.identity.symbolic ||
                         next_identity.workspace !=
                             implementation.spec.identity.workspace ||
                         identity.revision == 0U ||
                         identity.fingerprint == 0U ||
                         !std::isfinite(identity.maximum_relative_change) ||
                         identity.maximum_relative_change < 0.0
                     ? Status{StatusCode::numerical_failure, kMgCoefficient}
                     : Status{};
  agreed = implementation.services.reductions->checked_max(
      {&local_change, 1U}, {&global_change, 1U}, local);
  implementation.lowest =
      implementation.services.reductions->lowest_failing_rank();
  if (!agreed) {
    return agreed;
  }
  const bool unchanged =
      identity.revision == implementation.spec.coefficients.revision &&
      identity.fingerprint == implementation.spec.coefficients.fingerprint &&
      next_identity.fingerprint == implementation.spec.identity.fingerprint;
  if (unchanged) {
    implementation.lowest = -1;
    return {};
  }
  local = valid_coefficients(coefficients, implementation.spec.patch.cells)
              ? Status{}
              : Status{StatusCode::numerical_failure, kMgCoefficient};
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  std::uint64_t next_refresh = 0U;
  std::uint64_t next_rebuild = implementation.runtime_counters.hierarchy_rebuilds;
  RevisionToken next_generation = implementation.generation;
  const bool rebuild =
      global_change >
      implementation.spec.policy.coefficient_change_rebuild_ratio;
  local = !increment(implementation.runtime_counters.numeric_refreshes,
                     next_refresh) ||
                  (rebuild &&
                   (!increment(next_rebuild, next_rebuild) ||
                    !increment(implementation.generation, next_generation)))
              ? Status{StatusCode::invalid_plan, kMgCounter}
              : Status{};
  if (counters != nullptr && local) {
    std::uint64_t ignored = 0U;
    if (!increment(counters->numeric_refreshes, ignored) ||
        (rebuild &&
         !increment(counters->hierarchy_rebuilds, ignored))) {
      local = {StatusCode::invalid_plan, kMgCounter};
    }
  }
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  const MgCoefficientViews previous_coefficients = implementation.coefficients;
  implementation.coefficients = coefficients;
  local = build_coarse_coefficients(
      implementation, implementation.inactive_hierarchy_storage.data());
  implementation.coefficients = previous_coefficients;
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  std::copy(implementation.inactive_hierarchy_storage.begin(),
            implementation.inactive_hierarchy_storage.end(),
            implementation.hierarchy_storage.begin());
  implementation.coefficients = coefficients;
  implementation.spec.identity = next_identity;
  implementation.spec.coefficients = identity;
  implementation.spec.coefficients.maximum_relative_change = global_change;
  implementation.numeric =
      make_numeric_fingerprint(identity, implementation.symbolic);
  implementation.runtime_counters.numeric_refreshes = next_refresh;
  if (rebuild) {
    implementation.generation = next_generation;
    implementation.hierarchy = make_hierarchy_fingerprint(
        implementation.symbolic, identity, implementation.generation);
    implementation.runtime_counters.hierarchy_rebuilds = next_rebuild;
  }
  if (counters != nullptr) {
    ++counters->numeric_refreshes;
    if (rebuild) {
      ++counters->hierarchy_rebuilds;
    }
  }
  implementation.lowest = -1;
  return {};
}

LinearPreconditionerCertificate NativeCartesianMgPlan::certificate()
    const noexcept {
  if (implementation_ == nullptr) {
    return {};
  }
  return {implementation_->spec.identity, implementation_->hierarchy,
          LinearPreconditionerClass::flexible};
}

Status NativeCartesianMgPlan::apply(ConstFieldView residual,
                                    FieldView correction,
                                    std::uint32_t) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kMgApply};
  }
  Impl& implementation = *implementation_;
  Status local = validate_borrowed_services(implementation);
  Status agreed = implementation.services.reductions->consensus(local);
  implementation.lowest =
      implementation.services.reductions->lowest_failing_rank();
  if (!agreed) {
    return agreed;
  }
  local = !valid_field(residual, implementation.spec.patch.cells) ||
                         !valid_field(as_const(correction),
                                      implementation.spec.patch.cells) ||
                         detail::field_views_overlap(residual, correction) ||
                         implementation.services.workspace->overlaps_storage(
                             residual) ||
                         implementation.services.workspace->overlaps_storage(
                             correction)
                     ? Status{StatusCode::invalid_plan, kMgApply}
                     : Status{};
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  const Int3 cells = residual.interior;
  FieldView candidate =
      implementation.services.workspace->level(0U, MgWorkspaceSlot::solution);
  FieldView recursive = implementation.services.workspace->level(
      0U, MgWorkspaceSlot::residual);
  FieldView rhs = implementation.services.workspace->level(
      0U, MgWorkspaceSlot::rhs);
  FieldView temporary = implementation.services.workspace->level(
      0U, MgWorkspaceSlot::temporary);
  local = candidate.base == nullptr || recursive.base == nullptr ||
                  rhs.base == nullptr || temporary.base == nullptr
              ? Status{StatusCode::invalid_plan, kMgApply}
              : Status{};
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        candidate.unchecked({i, j, k}, 0U) = 0.0;
        rhs.unchecked({i, j, k}, 0U) =
            residual.unchecked({i, j, k}, 0U);
      }
    }
  }
  local = implementation.services.workspace->revise_level(
      0U, MgWorkspaceSlot::solution);
  local = local ? implementation.services.workspace->revise_level(
                      0U, MgWorkspaceSlot::rhs)
                : local;
  candidate = implementation.services.workspace->level(
      0U, MgWorkspaceSlot::solution);
  rhs = implementation.services.workspace->level(0U, MgWorkspaceSlot::rhs);
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  if (implementation.spec.null_space == MgNullSpace::constant) {
    local = project_constant(implementation, 0U, rhs);
    if (local) {
      local = implementation.services.workspace->revise_level(
          0U, MgWorkspaceSlot::rhs);
      rhs = implementation.services.workspace->level(
          0U, MgWorkspaceSlot::rhs);
    }
    agreed = consensus(implementation, local);
    if (!agreed) return agreed;
  }
  StageId stage = 700U;
  Status deferred{};
  local = v_cycle(implementation, 0U, stage, deferred);
  if (local && !deferred) {
    local = deferred;
  }
  agreed = consensus(implementation, local);
  if (!agreed) return agreed;
  candidate = implementation.services.workspace->level(
      0U, MgWorkspaceSlot::solution);
  if (implementation.spec.null_space == MgNullSpace::constant) {
    local = project_constant(implementation, 0U, candidate);
    if (local) {
      local = implementation.services.workspace->revise_level(
          0U, MgWorkspaceSlot::solution);
    }
    candidate = implementation.services.workspace->level(
        0U, MgWorkspaceSlot::solution);
    agreed = consensus(implementation, local);
    if (!agreed) {
      return agreed;
    }
  }
  local = compute_residual(implementation, 0U, stage);
  agreed = consensus(implementation, local);
  if (!agreed) return agreed;
  recursive = implementation.services.workspace->level(
      0U, MgWorkspaceSlot::residual);
  double local_projection[3]{};
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 cell{i, j, k};
        const double r = residual.unchecked(cell, 0U);
        const double az = r - recursive.unchecked(cell, 0U);
        local_projection[0] += r * r;
        local_projection[1] += r * az;
        local_projection[2] += az * az;
      }
    }
  }
  double projection[3]{};
  local = implementation.services.reductions->checked_sum(
      {local_projection, 3U}, {projection, 3U});
  implementation.lowest =
      implementation.services.reductions->lowest_failing_rank();
  if (!local) return local;
  const bool finite_projection = std::isfinite(projection[0]) &&
                                 std::isfinite(projection[1]) &&
                                 std::isfinite(projection[2]) &&
                                 projection[0] >= 0.0 &&
                                 projection[2] >= 0.0;
  if (!finite_projection) {
    return {StatusCode::numerical_failure, kMgApply};
  }
  const double scale = std::max(1.0, projection[0]);
  const double tiny = std::numeric_limits<double>::epsilon() * scale;
  const double alpha = projection[1] > 0.0 && projection[2] > tiny
                           ? projection[1] / projection[2]
                           : 0.0;
  const double minimized_squared =
      alpha == 0.0
          ? projection[0]
          : std::max(0.0, projection[0] -
                              projection[1] * projection[1] /
                                  projection[2]);
  implementation.last_initial = std::sqrt(projection[0]);
  implementation.last_final = std::sqrt(minimized_squared);
  if (!std::isfinite(alpha) || !std::isfinite(implementation.last_initial) ||
      !std::isfinite(implementation.last_final)) {
    return {StatusCode::numerical_failure, kMgApply};
  }
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 cell{i, j, k};
        candidate.unchecked(cell, 0U) *= alpha;
        recursive.unchecked(cell, 0U) =
            residual.unchecked(cell, 0U) -
            alpha * (residual.unchecked(cell, 0U) -
                     recursive.unchecked(cell, 0U));
      }
    }
  }
  std::uint64_t next_applications = 0U;
  local = increment(implementation.runtime_counters.applications,
                    next_applications)
              ? Status{}
              : Status{StatusCode::invalid_plan, kMgCounter};
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        correction.unchecked({i, j, k}, 0U) =
            candidate.unchecked({i, j, k}, 0U);
      }
    }
  }
  implementation.runtime_counters.applications = next_applications;
  implementation.lowest = -1;
  return {};
}

std::size_t NativeCartesianMgPlan::level_count() const noexcept {
  return implementation_ == nullptr ? 0U : implementation_->levels.size();
}

Status NativeCartesianMgPlan::level(std::size_t index,
                                    MgLevelView& out) const noexcept {
  if (implementation_ == nullptr || index >= implementation_->levels.size()) {
    return {StatusCode::invalid_plan, kMgPlan};
  }
  out = implementation_->levels[index].view;
  return {};
}

CoarseningKind NativeCartesianMgPlan::finest_coarsening() const noexcept {
  return implementation_ == nullptr ? CoarseningKind::full_xyz
                                    : implementation_->strategy.coarsening;
}

std::uint8_t NativeCartesianMgPlan::line_axis_mask() const noexcept {
  return implementation_ == nullptr ? 0U : implementation_->strategy.line_mask;
}

PlanFingerprint NativeCartesianMgPlan::symbolic_fingerprint() const noexcept {
  return implementation_ == nullptr ? 0U : implementation_->symbolic;
}

PlanFingerprint NativeCartesianMgPlan::numeric_fingerprint() const noexcept {
  return implementation_ == nullptr ? 0U : implementation_->numeric;
}

PlanFingerprint NativeCartesianMgPlan::hierarchy_fingerprint() const noexcept {
  return implementation_ == nullptr ? 0U : implementation_->hierarchy;
}

RevisionToken NativeCartesianMgPlan::generation() const noexcept {
  return implementation_ == nullptr ? 0U : implementation_->generation;
}

std::uintptr_t NativeCartesianMgPlan::hierarchy_storage_address()
    const noexcept {
  return implementation_ == nullptr
             ? 0U
             : reinterpret_cast<std::uintptr_t>(
                   implementation_->hierarchy_storage.data());
}

std::uintptr_t NativeCartesianMgPlan::workspace_storage_address()
    const noexcept {
  return implementation_ == nullptr
             ? 0U
             : implementation_->services.workspace->storage_address();
}

std::size_t NativeCartesianMgPlan::workspace_doubles() const noexcept {
  return implementation_ == nullptr || implementation_->levels.empty()
             ? 0U
             : implementation_->services.workspace->storage_doubles();
}

double NativeCartesianMgPlan::last_cycle_initial_residual() const noexcept {
  return implementation_ == nullptr ? 0.0 : implementation_->last_initial;
}

double NativeCartesianMgPlan::last_cycle_final_residual() const noexcept {
  return implementation_ == nullptr ? 0.0 : implementation_->last_final;
}

MgPlanCounters NativeCartesianMgPlan::counters() const noexcept {
  return implementation_ == nullptr ? MgPlanCounters{}
                                    : implementation_->runtime_counters;
}

int NativeCartesianMgPlan::lowest_failing_rank() const noexcept {
  return implementation_ == nullptr ? -1 : implementation_->lowest;
}

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
namespace detail {

struct NativeMgTestAccess {
  static void set(NativeCartesianMgPlan& plan, MgPlanCounters counters,
                  RevisionToken generation) noexcept {
    if (plan.implementation_ != nullptr) {
      plan.implementation_->runtime_counters = counters;
      plan.implementation_->generation = generation;
    }
  }
};

void set_mg_runtime_counters_for_test(
    NativeCartesianMgPlan& plan, MgPlanCounters counters,
    RevisionToken generation) noexcept {
  NativeMgTestAccess::set(plan, counters, generation);
}

}  // namespace detail
#endif

}  // namespace hundun::v04

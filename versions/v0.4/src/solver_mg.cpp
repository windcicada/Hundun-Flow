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
constexpr std::size_t kReplicatedCoarseCellLimit = 4096U;
constexpr std::size_t kReplicatedOperatorWidth = 7U;

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
         policy.line_relaxation_maximum_extent >= 2U &&
         (policy.point_smoother == MgPointSmootherKind::red_black ||
          policy.point_smoother == MgPointSmootherKind::chebyshev_jacobi) &&
         (policy.cycle == MgCycleKind::v_cycle ||
          policy.cycle == MgCycleKind::f_cycle) &&
         std::isfinite(policy.chebyshev_lower_spectrum_fraction) &&
         policy.chebyshev_lower_spectrum_fraction > 0.0 &&
         policy.chebyshev_lower_spectrum_fraction < 1.0;
}

bool valid_operator_class(MgOperatorClass operator_class) noexcept {
  return operator_class == MgOperatorClass::general ||
         operator_class ==
             MgOperatorClass::symmetric_diagonally_dominant_m_matrix;
}

bool valid_correction_scaling(MgCorrectionScaling scaling) noexcept {
  return scaling == MgCorrectionScaling::residual_minimizing ||
         scaling == MgCorrectionScaling::unit_linear;
}

bool valid_smoother_configuration(const NativeCartesianMgSpec& spec) noexcept {
  return spec.policy.point_smoother != MgPointSmootherKind::chebyshev_jacobi ||
         spec.operator_class ==
             MgOperatorClass::symmetric_diagonally_dominant_m_matrix;
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

bool valid_activity(MgDomainActivityView activity, Int3 cells) noexcept {
  const bool empty = activity.cells.size == 0U &&
                     activity.x_faces.size == 0U &&
                     activity.y_faces.size == 0U &&
                     activity.z_faces.size == 0U;
  if (empty) {
    return activity.local_fingerprint == 0U &&
           activity.collective_fingerprint == 0U;
  }
  const std::size_t cell_values = cell_count(cells);
  const std::size_t x_values = cell_count({cells.x + 1, cells.y, cells.z});
  const std::size_t y_values = cell_count({cells.x, cells.y + 1, cells.z});
  const std::size_t z_values = cell_count({cells.x, cells.y, cells.z + 1});
  if (activity.cells.data == nullptr || activity.x_faces.data == nullptr ||
      activity.y_faces.data == nullptr || activity.z_faces.data == nullptr ||
      activity.cells.size != cell_values ||
      activity.x_faces.size != x_values ||
      activity.y_faces.size != y_values ||
      activity.z_faces.size != z_values ||
      activity.local_fingerprint == 0U ||
      activity.collective_fingerprint == 0U) {
    return false;
  }
  const auto binary = [](Span<const std::uint8_t> values) noexcept {
    for (std::size_t index = 0U; index < values.size; ++index) {
      if (values.data[index] > 1U) return false;
    }
    return true;
  };
  if (!binary(activity.cells) || !binary(activity.x_faces) ||
      !binary(activity.y_faces) || !binary(activity.z_faces)) {
    return false;
  }
  const Int3 xe{cells.x + 1, cells.y, cells.z};
  const Int3 ye{cells.x, cells.y + 1, cells.z};
  const Int3 ze{cells.x, cells.y, cells.z + 1};
  const auto offset = [](Int3 shape, Int3 value) noexcept {
    return static_cast<std::size_t>(value.x) +
           static_cast<std::size_t>(shape.x) *
               (static_cast<std::size_t>(value.y) +
                static_cast<std::size_t>(shape.y) *
                    static_cast<std::size_t>(value.z));
  };
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        if (activity.cells.data[offset(cells, {i, j, k})] != 0U) continue;
        if (activity.x_faces.data[offset(xe, {i, j, k})] != 0U ||
            activity.x_faces.data[offset(xe, {i + 1, j, k})] != 0U ||
            activity.y_faces.data[offset(ye, {i, j, k})] != 0U ||
            activity.y_faces.data[offset(ye, {i, j + 1, k})] != 0U ||
            activity.z_faces.data[offset(ze, {i, j, k})] != 0U ||
            activity.z_faces.data[offset(ze, {i, j, k + 1})] != 0U) {
          return false;
        }
      }
    }
  }
  return true;
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
  hash = mix(hash, static_cast<std::uint64_t>(spec.operator_class));
  if (spec.correction_scaling != MgCorrectionScaling::residual_minimizing) {
      hash = mix(hash, UINT64_C(0x6d675f756e69745f));
      hash = mix(hash, static_cast<std::uint64_t>(spec.correction_scaling));
  }
  hash = mix(hash, static_cast<std::uint64_t>(spec.policy.point_smoother));
  hash = mix(hash, static_cast<std::uint64_t>(spec.policy.cycle));
  hash = mix(hash, static_cast<std::uint64_t>(strategy.coarsening));
  hash = mix(hash, strategy.line_mask);
  hash = mix(hash, spec.policy.pre_sweeps);
  hash = mix(hash, spec.policy.post_sweeps);
  hash = mix(hash, spec.policy.maximum_levels);
  hash = mix(hash, spec.policy.coarse_sweeps);
  hash = mix(hash, spec.policy.minimum_coarse_extent);
  hash = mix(hash, spec.policy.line_relaxation_maximum_extent);
  hash = mix(hash, spec.activity.collective_fingerprint);
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
  std::uint64_t chebyshev_fraction_bits = 0U;
  std::memcpy(&chebyshev_fraction_bits,
              &spec.policy.chebyshev_lower_spectrum_fraction,
              sizeof(chebyshev_fraction_bits));
  hash = mix(hash, chebyshev_fraction_bits);
  return finish(hash);
}

PlanFingerprint public_symbolic_fingerprint(
    const NativeCartesianMgSpec& spec, Strategy strategy) noexcept {
  std::uint64_t hash = UINT64_C(1469598103934665603);
  hash = mix(hash, spec.geometry->fingerprint());
  hash = mix(hash, spec.identity.symbolic);
  hash = mix(hash, static_cast<std::uint64_t>(spec.operator_class));
  if (spec.correction_scaling != MgCorrectionScaling::residual_minimizing) {
    hash = mix(hash, UINT64_C(0x6d675f756e69745f));
    hash = mix(hash, static_cast<std::uint64_t>(spec.correction_scaling));
  }
  hash = mix(hash, static_cast<std::uint64_t>(spec.policy.point_smoother));
  hash = mix(hash, static_cast<std::uint64_t>(spec.policy.cycle));
  hash = mix(hash, spec.policy.pre_sweeps);
  hash = mix(hash, spec.policy.post_sweeps);
  std::uint64_t fraction_bits = 0U;
  std::memcpy(&fraction_bits,
              &spec.policy.chebyshev_lower_spectrum_fraction,
              sizeof(fraction_bits));
  hash = mix(hash, fraction_bits);
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.x_min));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.x_max));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.y_min));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.y_max));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.z_min));
  hash = mix(hash, static_cast<std::uint64_t>(spec.boundaries.z_max));
  hash = mix(hash, static_cast<std::uint64_t>(strategy.coarsening));
  hash = mix(hash, spec.activity.collective_fingerprint);
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

PlanFingerprint make_collective_preconditioner_fingerprint(
    PlanFingerprint symbolic, RevisionToken generation,
    std::uint64_t numeric_refreshes) noexcept {
  return finish(mix(mix(mix(UINT64_C(1469598103934665603), symbolic),
                        generation),
                    numeric_refreshes));
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
      !valid_operator_class(spec.operator_class) ||
      !valid_correction_scaling(spec.correction_scaling) ||
      !valid_smoother_configuration(spec) ||
      !valid_boundary_pair(spec.boundaries.x_min, spec.boundaries.x_max) ||
      !valid_boundary_pair(spec.boundaries.y_min, spec.boundaries.y_max) ||
      !valid_boundary_pair(spec.boundaries.z_min, spec.boundaries.z_max) ||
      !valid_services(services, spec.identity) ||
      !valid_coefficients(coefficients, spec.patch.cells) ||
      !valid_activity(spec.activity, spec.patch.cells)) {
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
  // Cold-owned positive-zero row storage used by physical Dirichlet point
  // rows.  Its capacity covers every local x extent, so the hot kernel can
  // use a resolved pointer without a per-cell mode branch or a literal
  // reconstruction.
  std::vector<double> point_boundary_zero_storage;
  std::vector<std::uint8_t> activity_cells;
  std::vector<std::uint8_t> activity_x_faces;
  std::vector<std::uint8_t> activity_y_faces;
  std::vector<std::uint8_t> activity_z_faces;
  // Compile/update-only storage for conservative tensor aggregation.  The
  // first two blocks are ping-pong fields; the final two are MPI planes.
  // Keeping it on the plan makes numeric refresh allocation-free.
  std::vector<double> coefficient_tensor_storage;
  std::size_t coefficient_tensor_values{};
  std::size_t coefficient_tensor_plane{};
  Strategy strategy{};
  MgPlanCounters runtime_counters{};
  // A small point-smoother coarse level may be replicated on every rank.
  // All descriptor maps, packed layouts, and numeric workspaces are cold-owned
  // here so that the eligible apply path has no ownership or allocation work.
  bool replicated_coarse{};
  Int3 replicated_global_shape{};
  std::size_t replicated_global_cells{};
  std::vector<std::int32_t> replicated_descriptors;
  std::vector<int> replicated_rhs_counts;
  std::vector<int> replicated_rhs_displacements;
  std::vector<int> replicated_operator_counts;
  std::vector<int> replicated_operator_displacements;
  std::vector<int> replicated_owner;
  std::vector<std::size_t> replicated_owner_local;
  std::vector<std::size_t> replicated_owner_order;
  std::vector<std::uint8_t> replicated_owner_parity;
  std::vector<std::size_t> replicated_neighbors;
  std::vector<double> replicated_rhs_local;
  std::vector<double> replicated_rhs_gather;
  std::vector<double> replicated_rhs_global;
  std::vector<double> replicated_solution;
  RevisionToken replicated_solution_revision{};
  std::vector<double> replicated_snapshot;
  std::vector<double> replicated_action;
  std::vector<double> replicated_operator_local;
  std::vector<double> replicated_operator_gather;
  std::vector<double> replicated_operator_active;
  std::vector<double> replicated_operator_inactive;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  bool reference_point_actions{};
  bool reference_chebyshev_lifecycle{};
  bool reference_replicated_coarse{};
  bool request_replicated_coarse{};
  bool skip_final_projection{};
  bool force_chebyshev_invalid{};
  detail::MgCycleFailurePhase cycle_failure_phase{
      detail::MgCycleFailurePhase::none};
  std::size_t cycle_failure_level{};
  int cycle_failure_rank{-1};
  bool cycle_failure_deferred{};
  bool cycle_failure_consumed{};
  mutable detail::MgMatrixWorkCounters matrix_work{};
#endif
  PlanFingerprint symbolic{};
  PlanFingerprint numeric{};
  PlanFingerprint hierarchy{};
  RevisionToken generation{};
  double last_initial{};
  double last_final{};
  int rank{};
  int size{};
  int lowest{-1};
  bool prepared_halo_epoch{};
  Status prepared_halo_deferred{};
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

template <class Implementation>
Status validate_prepared_batch(
    const Implementation& implementation,
    const LinearPreconditionerBatchDescriptor& descriptor) noexcept {
  const SolverWorkspace* const workspace = descriptor.workspace;
  const Int3 shape = implementation.spec.patch.cells;
  const Status borrowed = validate_borrowed_services(implementation);
  if (!borrowed) {
    return borrowed;
  }
  if (workspace == nullptr || workspace->fingerprint() == 0U ||
      (workspace->requirements().algorithm != LinearAlgorithm::fgmres &&
       workspace->requirements().algorithm != LinearAlgorithm::bicgstab) ||
      descriptor.slot_count == 0U || descriptor.maximum_applications == 0U ||
      !same_shape(descriptor.shape, shape) ||
      descriptor.input_slot_begin >= workspace->requirements().vector_slots ||
      descriptor.output_slot_begin >= workspace->requirements().vector_slots) {
    return {StatusCode::invalid_plan, kMgApply};
  }
  if (workspace->requirements().algorithm == LinearAlgorithm::bicgstab &&
      (descriptor.input_slot_begin != 3U ||
       descriptor.output_slot_begin != 5U || descriptor.slot_count != 2U ||
       descriptor.maximum_applications % 2U != 0U)) {
    return {StatusCode::invalid_plan, kMgApply};
  }
  const std::size_t input_end =
      static_cast<std::size_t>(descriptor.input_slot_begin) +
      static_cast<std::size_t>(descriptor.slot_count);
  const std::size_t output_end =
      static_cast<std::size_t>(descriptor.output_slot_begin) +
      static_cast<std::size_t>(descriptor.slot_count);
  if (input_end > workspace->requirements().vector_slots ||
      output_end > workspace->requirements().vector_slots ||
      implementation.runtime_counters.applications >
          std::numeric_limits<std::uint64_t>::max() -
              descriptor.maximum_applications) {
    return {StatusCode::invalid_plan, kMgCounter};
  }
  for (std::size_t index = 0U; index < descriptor.slot_count; ++index) {
    const std::uint8_t input_slot = static_cast<std::uint8_t>(
        static_cast<std::size_t>(descriptor.input_slot_begin) + index);
    const std::uint8_t output_slot = static_cast<std::uint8_t>(
        static_cast<std::size_t>(descriptor.output_slot_begin) + index);
    const FieldView input = workspace->vector(input_slot, descriptor.shape);
    const FieldView output = workspace->vector(output_slot, descriptor.shape);
    if (!valid_field(as_const(input), shape) ||
        !valid_field(as_const(output), shape) ||
        detail::field_views_overlap(input, output) ||
        implementation.services.workspace->overlaps_storage(as_const(input)) ||
        implementation.services.workspace->overlaps_storage(as_const(output))) {
      return {StatusCode::invalid_plan, kMgApply};
    }
  }
  const MgWorkspaceRequirements& requirements =
      implementation.services.workspace->requirements();
  if (requirements.level_count == 0U ||
      requirements.level_count > detail::kMgMaximumLevels) {
    return {StatusCode::invalid_plan, kMgApply};
  }
  for (std::size_t level = 0U; level < requirements.level_count; ++level) {
    const Int3 level_shape = requirements.levels[level].patch.cells;
    for (std::uint8_t slot = 0U; slot < 4U; ++slot) {
      const MgWorkspaceSlot mg_slot =
          static_cast<MgWorkspaceSlot>(slot);
      const FieldView field = implementation.services.workspace->level(
          level, mg_slot);
      if (!valid_field(as_const(field), level_shape) ||
          workspace->overlaps_storage(as_const(field))) {
        return {StatusCode::invalid_plan, kMgApply};
      }
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

detail::MgPointBoundaryMode cold_point_boundary_mode(
    const MeshPatch& patch, Int3 global, CartesianAxis axis, bool minimum,
    MgBoundaryKind boundary) noexcept {
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
                                  : (axis == CartesianAxis::y
                                         ? global.y
                                         : global.z);
  // A local rank interface, and a periodic global face, are both populated
  // by the reserved one-cell halo.  Only a physical non-periodic face needs
  // a self/zero mode in the point row.
  const bool at_global = minimum ? begin == 0 : begin + count == extent;
  if (!at_global || boundary == MgBoundaryKind::periodic) {
    return detail::MgPointBoundaryMode::padded_ghost;
  }
  return boundary == MgBoundaryKind::neumann
             ? detail::MgPointBoundaryMode::self
             : detail::MgPointBoundaryMode::zero;
}

template <class Implementation>
Status consensus(Implementation& implementation, Status local) noexcept {
  const int inherited_collective_failure = implementation.lowest;
  const Status result = implementation.services.reductions->consensus(local);
  const int reduced =
      implementation.services.reductions->lowest_failing_rank();
  // HaloEngine and other inner collective authorities already return the
  // same Status on every rank.  Re-consensus of that Status would select
  // rank zero and destroy their true source.  A successful consensus clears
  // provenance; a failed one preserves an inherited collective source.
  implementation.lowest =
      !result && inherited_collective_failure >= 0
          ? inherited_collective_failure
          : reduced;
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
    storage.point_boundary_modes = {
        cold_point_boundary_mode(required.patch, required.global_shape,
                                 CartesianAxis::x, true,
                                 implementation.spec.boundaries.x_min),
        cold_point_boundary_mode(required.patch, required.global_shape,
                                 CartesianAxis::x, false,
                                 implementation.spec.boundaries.x_max),
        cold_point_boundary_mode(required.patch, required.global_shape,
                                 CartesianAxis::y, true,
                                 implementation.spec.boundaries.y_min),
        cold_point_boundary_mode(required.patch, required.global_shape,
                                 CartesianAxis::y, false,
                                 implementation.spec.boundaries.y_max),
        cold_point_boundary_mode(required.patch, required.global_shape,
                                 CartesianAxis::z, true,
                                 implementation.spec.boundaries.z_min),
        cold_point_boundary_mode(required.patch, required.global_shape,
                                 CartesianAxis::z, false,
                                 implementation.spec.boundaries.z_max)};
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
  for (std::size_t level = 0U; level < implementation.levels.size();
       ++level) {
    detail::MgLevelStorage& selected = implementation.levels[level];
    selected.point_smoother =
        level + 1U < implementation.levels.size() &&
                selected.view.line_axis_mask == 0U &&
                implementation.spec.policy.point_smoother ==
                    MgPointSmootherKind::chebyshev_jacobi
            ? MgPointSmootherKind::chebyshev_jacobi
            : MgPointSmootherKind::red_black;
    selected.chebyshev_lower_spectrum_fraction =
        selected.point_smoother == MgPointSmootherKind::chebyshev_jacobi
            ? implementation.spec.policy.chebyshev_lower_spectrum_fraction
            : 0.0;
  }
  implementation.hierarchy_storage.assign(total_hierarchy, 0.0);
  implementation.inactive_hierarchy_storage.assign(total_hierarchy, 0.0);
  std::size_t maximum_values = 0U;
  std::size_t maximum_x = 0U;
  std::size_t maximum_plane = 0U;
  for (const detail::MgLevelStorage& level : implementation.levels) {
    maximum_values = std::max(
        {maximum_values, level.cells, level.x_faces, level.y_faces,
         level.z_faces});
    const Int3 cells = level.view.local_shape;
    maximum_x = std::max(maximum_x, static_cast<std::size_t>(cells.x));
    const std::size_t nx = static_cast<std::size_t>(cells.x) + 1U;
    const std::size_t ny = static_cast<std::size_t>(cells.y) + 1U;
    const std::size_t nz = static_cast<std::size_t>(cells.z) + 1U;
    maximum_plane =
        std::max({maximum_plane, nx * ny, nx * nz, ny * nz});
  }
  if (maximum_values == 0U || maximum_x == 0U || maximum_plane == 0U ||
      maximum_values >
          (std::numeric_limits<std::size_t>::max() - 2U * maximum_plane) /
              2U) {
    return {StatusCode::invalid_plan, kMgPlan};
  }
  implementation.coefficient_tensor_values = maximum_values;
  implementation.coefficient_tensor_plane = maximum_plane;
  implementation.point_boundary_zero_storage.assign(maximum_x, 0.0);
  implementation.coefficient_tensor_storage.assign(
      2U * maximum_values + 2U * maximum_plane, 0.0);
  return {};
}

std::size_t flat(Int3 shape, Int3 cell) noexcept;
std::size_t face_flat(Int3 extents, Int3 face) noexcept;

bool add_counter(std::uint64_t current, std::uint64_t amount,
                 std::uint64_t& next) noexcept {
  if (amount > std::numeric_limits<std::uint64_t>::max() - current) {
    return false;
  }
  next = current + amount;
  return true;
}

template <class Implementation>
Status record_replicated_collective(Implementation& implementation,
                                    std::size_t logical_bytes) noexcept {
  std::uint64_t next_calls = 0U;
  std::uint64_t next_bytes = 0U;
  if (!increment(implementation.runtime_counters.blocking_collectives,
                 next_calls) ||
      !add_counter(implementation.runtime_counters.collective_logical_bytes,
                   static_cast<std::uint64_t>(logical_bytes), next_bytes)) {
    return {StatusCode::invalid_plan, kMgCounter};
  }
  implementation.runtime_counters.blocking_collectives = next_calls;
  implementation.runtime_counters.collective_logical_bytes = next_bytes;
  return {};
}

template <class Implementation>
Status initialize_replicated_coarse(Implementation& implementation,
                                    MgPlanCounters* external,
                                    std::uint64_t& external_collectives,
                                    std::uint64_t& external_bytes) noexcept {
  if (implementation.size <= 1 || implementation.levels.empty()) {
    return {};
  }
  const detail::MgLevelStorage& coarse = implementation.levels.back();
  const std::size_t global_cells = cell_count(coarse.view.global_shape);
  const std::size_t local_cells = coarse.cells;
  if (coarse.view.line_axis_mask != 0U || global_cells == 0U ||
      global_cells > kReplicatedCoarseCellLimit || local_cells == 0U ||
      local_cells > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      global_cells > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      global_cells > std::numeric_limits<std::size_t>::max() /
                          kReplicatedOperatorWidth) {
    return {};
  }

  const std::size_t rank_count = static_cast<std::size_t>(implementation.size);
  Status local{};
  try {
    implementation.replicated_global_shape = coarse.view.global_shape;
    implementation.replicated_global_cells = global_cells;
    implementation.replicated_descriptors.assign(rank_count * 6U, 0);
    implementation.replicated_rhs_counts.assign(rank_count, 0);
    implementation.replicated_rhs_displacements.assign(rank_count, 0);
    implementation.replicated_operator_counts.assign(rank_count, 0);
    implementation.replicated_operator_displacements.assign(rank_count, 0);
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, 0U};
  } catch (...) {
    local = {StatusCode::invalid_plan, kMgPlan};
  }
  Status agreed = consensus(implementation, local);
  if (!agreed) return agreed;

  const std::array<std::int32_t, 6U> descriptor{
      coarse.patch.begin.x, coarse.patch.begin.y, coarse.patch.begin.z,
      coarse.patch.cells.x, coarse.patch.cells.y, coarse.patch.cells.z};
  const std::size_t descriptor_bytes = rank_count * 6U * sizeof(std::int32_t);
  std::uint64_t descriptor_collectives = external_collectives;
  std::uint64_t descriptor_bytes_total = external_bytes;
  local = external == nullptr ||
                  (add_counter(external_collectives, 1U,
                               descriptor_collectives) &&
                   add_counter(external_bytes,
                               static_cast<std::uint64_t>(descriptor_bytes),
                               descriptor_bytes_total))
              ? Status{}
              : Status{StatusCode::invalid_plan, kMgCounter};
  if (local && external != nullptr) {
    external_collectives = descriptor_collectives;
    external_bytes = descriptor_bytes_total;
  }
  agreed = consensus(implementation, local);
  if (!agreed) return agreed;
  local = record_replicated_collective(implementation, descriptor_bytes);
  agreed = consensus(implementation, local);
  if (!agreed) return agreed;
  const int mpi_status = MPI_Allgather(
      descriptor.data(), static_cast<int>(descriptor.size()), MPI_INT32_T,
      implementation.replicated_descriptors.data(),
      static_cast<int>(descriptor.size()), MPI_INT32_T,
      implementation.communicator);
  agreed = consensus(
      implementation,
      mpi_status == MPI_SUCCESS
          ? Status{}
          : Status{StatusCode::mpi_failure, kMgCollective});
  if (!agreed) return agreed;

  const Int3 global = coarse.view.global_shape;
  std::vector<int> owners;
  std::vector<std::size_t> owner_local;
  std::vector<std::size_t> owner_order;
  std::vector<std::uint8_t> owner_parity;
  try {
    owners.assign(global_cells, -1);
    owner_local.assign(global_cells, 0U);
    owner_order.assign(global_cells, 0U);
    owner_parity.assign(global_cells, 0U);
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, 0U};
  } catch (...) {
    local = {StatusCode::invalid_plan, kMgPlan};
  }
  agreed = consensus(implementation, local);
  if (!agreed) return agreed;

  bool exact_cover = true;
  int rhs_displacement = 0;
  int operator_displacement = 0;
  for (int rank = 0; rank < implementation.size && exact_cover; ++rank) {
    const std::size_t offset = static_cast<std::size_t>(rank) * 6U;
    const Int3 begin{implementation.replicated_descriptors[offset + 0U],
                     implementation.replicated_descriptors[offset + 1U],
                     implementation.replicated_descriptors[offset + 2U]};
    const Int3 cells{implementation.replicated_descriptors[offset + 3U],
                     implementation.replicated_descriptors[offset + 4U],
                     implementation.replicated_descriptors[offset + 5U]};
    const std::size_t owned_cells = cell_count(cells);
    if (owned_cells == 0U || begin.x < 0 || begin.y < 0 || begin.z < 0 ||
        begin.x > global.x || begin.y > global.y || begin.z > global.z ||
        cells.x > global.x - begin.x || cells.y > global.y - begin.y ||
        cells.z > global.z - begin.z ||
        owned_cells > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        owned_cells >
            (static_cast<std::size_t>(std::numeric_limits<int>::max()) /
             kReplicatedOperatorWidth) ||
        rhs_displacement > std::numeric_limits<int>::max() -
                               static_cast<int>(owned_cells) ||
        operator_displacement >
            std::numeric_limits<int>::max() -
                static_cast<int>(owned_cells * kReplicatedOperatorWidth)) {
      exact_cover = false;
      break;
    }
    implementation.replicated_rhs_counts[rank] = static_cast<int>(owned_cells);
    implementation.replicated_rhs_displacements[rank] = rhs_displacement;
    implementation.replicated_operator_counts[rank] =
        static_cast<int>(owned_cells * kReplicatedOperatorWidth);
    implementation.replicated_operator_displacements[rank] =
        operator_displacement;
    rhs_displacement += static_cast<int>(owned_cells);
    operator_displacement += static_cast<int>(owned_cells *
                                               kReplicatedOperatorWidth);
    for (std::int32_t k = 0; k < cells.z; ++k) {
      for (std::int32_t j = 0; j < cells.y; ++j) {
        for (std::int32_t i = 0; i < cells.x; ++i) {
          const Int3 global_cell{begin.x + i, begin.y + j, begin.z + k};
          const std::size_t global_index = flat(global, global_cell);
          if (global_index >= owners.size() || owners[global_index] != -1) {
            exact_cover = false;
            break;
          }
          owners[global_index] = rank;
          const std::size_t local_index = flat(cells, {i, j, k});
          owner_local[global_index] = local_index;
          owner_order[static_cast<std::size_t>(
              implementation.replicated_rhs_displacements[rank]) +
                      local_index] = global_index;
          owner_parity[global_index] =
              static_cast<std::uint8_t>((i + j + k) & 1);
        }
        if (!exact_cover) break;
      }
      if (!exact_cover) break;
    }
  }
  if (exact_cover) {
    exact_cover = std::all_of(
        owners.begin(), owners.end(), [](int owner) noexcept { return owner >= 0; });
  }
  if (!exact_cover || rhs_displacement != static_cast<int>(global_cells) ||
      operator_displacement !=
          static_cast<int>(global_cells * kReplicatedOperatorWidth)) {
    // Descriptor mismatch is an ordinary eligibility failure.  The existing
    // distributed route remains the complete plan and no replicated state is
    // published.
    implementation.replicated_descriptors.clear();
    implementation.replicated_rhs_counts.clear();
    implementation.replicated_rhs_displacements.clear();
    implementation.replicated_operator_counts.clear();
    implementation.replicated_operator_displacements.clear();
    return {};
  }

  const std::size_t operator_bytes =
      global_cells * kReplicatedOperatorWidth * sizeof(double);
  std::uint64_t operator_collectives = external_collectives;
  std::uint64_t operator_bytes_total = external_bytes;
  local = external == nullptr ||
                  (add_counter(external_collectives, 1U,
                               operator_collectives) &&
                   add_counter(external_bytes,
                               static_cast<std::uint64_t>(operator_bytes),
                               operator_bytes_total))
              ? Status{}
              : Status{StatusCode::invalid_plan, kMgCounter};
  if (local && external != nullptr) {
    external_collectives = operator_collectives;
    external_bytes = operator_bytes_total;
  }
  agreed = consensus(implementation, local);
  if (!agreed) return agreed;

  try {
    implementation.replicated_owner = std::move(owners);
    implementation.replicated_owner_local = std::move(owner_local);
    implementation.replicated_owner_order = std::move(owner_order);
    implementation.replicated_owner_parity = std::move(owner_parity);
    implementation.replicated_neighbors.assign(global_cells * 6U,
                                               global_cells);
    implementation.replicated_rhs_local.assign(local_cells, 0.0);
    implementation.replicated_rhs_gather.assign(global_cells, 0.0);
    implementation.replicated_rhs_global.assign(global_cells, 0.0);
    implementation.replicated_solution.assign(global_cells, 0.0);
    implementation.replicated_snapshot.assign(global_cells + 1U, 0.0);
    implementation.replicated_action.assign(global_cells, 0.0);
    implementation.replicated_operator_local.assign(
        local_cells * kReplicatedOperatorWidth, 0.0);
    implementation.replicated_operator_gather.assign(
        global_cells * kReplicatedOperatorWidth, 0.0);
    implementation.replicated_operator_active.assign(
        global_cells * kReplicatedOperatorWidth, 0.0);
    implementation.replicated_operator_inactive.assign(
        global_cells * kReplicatedOperatorWidth, 0.0);
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, 0U};
  } catch (...) {
    local = {StatusCode::invalid_plan, kMgPlan};
  }
  agreed = consensus(implementation, local);
  if (!agreed) return agreed;
  const auto neighbor_index = [&](Int3 cell, CartesianAxis axis,
                                  int direction) noexcept -> std::size_t {
    Int3 selected = cell;
    std::int32_t* coordinate = axis == CartesianAxis::x
                                   ? &selected.x
                                   : (axis == CartesianAxis::y ? &selected.y
                                                               : &selected.z);
    const std::int32_t extent = axis == CartesianAxis::x
                                    ? global.x
                                    : (axis == CartesianAxis::y ? global.y
                                                                : global.z);
    *coordinate += direction;
    if (*coordinate >= 0 && *coordinate < extent) {
      return flat(global, selected);
    }
    const MgBoundaryKind minimum =
        axis == CartesianAxis::x
            ? implementation.spec.boundaries.x_min
            : (axis == CartesianAxis::y ? implementation.spec.boundaries.y_min
                                        : implementation.spec.boundaries.z_min);
    const MgBoundaryKind maximum =
        axis == CartesianAxis::x
            ? implementation.spec.boundaries.x_max
            : (axis == CartesianAxis::y ? implementation.spec.boundaries.y_max
                                        : implementation.spec.boundaries.z_max);
    const MgBoundaryKind selected_boundary =
        direction < 0 ? minimum : maximum;
    if (selected_boundary == MgBoundaryKind::periodic) {
      *coordinate = direction < 0 ? extent - 1 : 0;
      return flat(global, selected);
    }
    if (selected_boundary == MgBoundaryKind::neumann) {
      *coordinate = direction < 0 ? 0 : extent - 1;
      return flat(global, selected);
    }
    return global_cells;
  };
  for (std::int32_t k = 0; k < global.z; ++k) {
    for (std::int32_t j = 0; j < global.y; ++j) {
      for (std::int32_t i = 0; i < global.x; ++i) {
        const Int3 cell{i, j, k};
        const std::size_t index = flat(global, cell) * 6U;
        implementation.replicated_neighbors[index + 0U] =
            neighbor_index(cell, CartesianAxis::x, -1);
        implementation.replicated_neighbors[index + 1U] =
            neighbor_index(cell, CartesianAxis::x, 1);
        implementation.replicated_neighbors[index + 2U] =
            neighbor_index(cell, CartesianAxis::y, -1);
        implementation.replicated_neighbors[index + 3U] =
            neighbor_index(cell, CartesianAxis::y, 1);
        implementation.replicated_neighbors[index + 4U] =
            neighbor_index(cell, CartesianAxis::z, -1);
        implementation.replicated_neighbors[index + 5U] =
            neighbor_index(cell, CartesianAxis::z, 1);
      }
    }
  }
  implementation.replicated_coarse = true;
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
        const std::size_t index = face_flat(x_extents, face);
        x[index] = implementation.coefficients.x.unchecked(face) *
                   (implementation.spec.activity.x_faces.size == 0U
                        ? 1.0
                        : implementation.spec.activity.x_faces.data[index]);
      }
    }
  }
  for (std::int32_t k = 0; k < y_extents.z; ++k) {
    for (std::int32_t j = 0; j < y_extents.y; ++j) {
      for (std::int32_t i = 0; i < y_extents.x; ++i) {
        const Int3 face{i, j, k};
        const std::size_t index = face_flat(y_extents, face);
        y[index] = implementation.coefficients.y.unchecked(face) *
                   (implementation.spec.activity.y_faces.size == 0U
                        ? 1.0
                        : implementation.spec.activity.y_faces.data[index]);
      }
    }
  }
  for (std::int32_t k = 0; k < z_extents.z; ++k) {
    for (std::int32_t j = 0; j < z_extents.y; ++j) {
      for (std::int32_t i = 0; i < z_extents.x; ++i) {
        const Int3 face{i, j, k};
        const std::size_t index = face_flat(z_extents, face);
        z[index] = implementation.coefficients.z.unchecked(face) *
                   (implementation.spec.activity.z_faces.size == 0U
                        ? 1.0
                        : implementation.spec.activity.z_faces.data[index]);
      }
    }
  }
  if (implementation.spec.activity.cells.size != 0U) {
    for (std::int32_t k = 0; k < cells.z; ++k) {
      for (std::int32_t j = 0; j < cells.y; ++j) {
        for (std::int32_t i = 0; i < cells.x; ++i) {
          const Int3 cell{i, j, k};
          const std::size_t index = flat(cells, cell);
          if (implementation.spec.activity.cells.data[index] == 0U) {
            diagonal[index] = 1.0;
            continue;
          }
          const double original_faces =
              implementation.coefficients.x.unchecked(cell) +
              implementation.coefficients.x.unchecked({i + 1, j, k}) +
              implementation.coefficients.y.unchecked(cell) +
              implementation.coefficients.y.unchecked({i, j + 1, k}) +
              implementation.coefficients.z.unchecked(cell) +
              implementation.coefficients.z.unchecked({i, j, k + 1});
          const double active_faces =
              x[face_flat(x_extents, cell)] +
              x[face_flat(x_extents, {i + 1, j, k})] +
              y[face_flat(y_extents, cell)] +
              y[face_flat(y_extents, {i, j + 1, k})] +
              z[face_flat(z_extents, cell)] +
              z[face_flat(z_extents, {i, j, k + 1})];
          diagonal[index] = diagonal[index] - original_faces + active_faces;
        }
      }
    }
  }
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
Status validate_certified_hierarchy(const Implementation& implementation,
                                    const double* base) noexcept {
  if (implementation.spec.operator_class !=
      MgOperatorClass::symmetric_diagonally_dominant_m_matrix) {
    return {};
  }
  if (base == nullptr || implementation.levels.empty()) {
    return {StatusCode::numerical_failure, kMgCoefficient};
  }
  for (std::size_t level_index = 0U;
       level_index < implementation.levels.size(); ++level_index) {
    const detail::MgLevelStorage& level = implementation.levels[level_index];
    const Int3 cells = level.view.local_shape;
    const Int3 xe{cells.x + 1, cells.y, cells.z};
    const Int3 ye{cells.x, cells.y + 1, cells.z};
    const Int3 ze{cells.x, cells.y, cells.z + 1};
    const double* const diagonal = detail::block(base, level.diagonal_offset);
    const double* const x = detail::block(base, level.x_offset);
    const double* const y = detail::block(base, level.y_offset);
    const double* const z = detail::block(base, level.z_offset);
    for (std::int32_t k = 0; k < cells.z; ++k) {
      for (std::int32_t j = 0; j < cells.y; ++j) {
        for (std::int32_t i = 0; i < cells.x; ++i) {
          const Int3 cell{i, j, k};
          const std::size_t index = flat(cells, cell);
          const double faces =
              x[face_flat(xe, cell)] +
              x[face_flat(xe, {i + 1, j, k})] +
              y[face_flat(ye, cell)] +
              y[face_flat(ye, {i, j + 1, k})] +
              z[face_flat(ze, cell)] +
              z[face_flat(ze, {i, j, k + 1})];
          const bool inactive =
              level_index == 0U &&
              implementation.spec.activity.cells.size != 0U &&
              implementation.spec.activity.cells.data[index] == 0U;
          if (!std::isfinite(diagonal[index]) ||
              !(diagonal[index] > 0.0) || !std::isfinite(faces) ||
              faces < 0.0) {
            return {StatusCode::numerical_failure, kMgCoefficient};
          }
          if (inactive) {
            if (diagonal[index] != 1.0 || faces != 0.0) {
              return {StatusCode::numerical_failure, kMgCoefficient};
            }
          } else {
            if (faces > diagonal[index]) {
              return {StatusCode::numerical_failure, kMgCoefficient};
            }
          }
        }
      }
    }
  }
  return {};
}

template <class Implementation>
Status build_replicated_operator(Implementation& implementation,
                                  const double* base,
                                  std::vector<double>& target) noexcept {
  if (!implementation.replicated_coarse || base == nullptr ||
      implementation.levels.empty()) {
    return {};
  }
  const detail::MgLevelStorage& coarse = implementation.levels.back();
  const Int3 cells = coarse.view.local_shape;
  const Int3 x_extents{cells.x + 1, cells.y, cells.z};
  const Int3 y_extents{cells.x, cells.y + 1, cells.z};
  const Int3 z_extents{cells.x, cells.y, cells.z + 1};
  const double* const diagonal = detail::block(base, coarse.diagonal_offset);
  const double* const x = detail::block(base, coarse.x_offset);
  const double* const y = detail::block(base, coarse.y_offset);
  const double* const z = detail::block(base, coarse.z_offset);
  Status local = coarse.cells == implementation.replicated_rhs_counts[
                                    implementation.rank] &&
                         target.size() == implementation.replicated_global_cells *
                                              kReplicatedOperatorWidth
                     ? Status{}
                     : Status{StatusCode::invalid_plan, kMgPlan};
  Status agreed = consensus(implementation, local);
  if (!agreed) return agreed;
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 cell{i, j, k};
        const std::size_t local_index = flat(cells, cell);
        const std::size_t row = local_index * kReplicatedOperatorWidth;
        implementation.replicated_operator_local[row + 0U] =
            diagonal[local_index];
        implementation.replicated_operator_local[row + 1U] =
            x[face_flat(x_extents, cell)];
        implementation.replicated_operator_local[row + 2U] =
            x[face_flat(x_extents, {i + 1, j, k})];
        implementation.replicated_operator_local[row + 3U] =
            y[face_flat(y_extents, cell)];
        implementation.replicated_operator_local[row + 4U] =
            y[face_flat(y_extents, {i, j + 1, k})];
        implementation.replicated_operator_local[row + 5U] =
            z[face_flat(z_extents, cell)];
        implementation.replicated_operator_local[row + 6U] =
            z[face_flat(z_extents, {i, j, k + 1})];
      }
    }
  }
  local = record_replicated_collective(
      implementation,
      implementation.replicated_global_cells * kReplicatedOperatorWidth *
          sizeof(double));
  agreed = consensus(implementation, local);
  if (!agreed) return agreed;
  const int mpi_status = MPI_Allgatherv(
      implementation.replicated_operator_local.data(),
      implementation.replicated_operator_counts[implementation.rank],
      MPI_DOUBLE, implementation.replicated_operator_gather.data(),
      implementation.replicated_operator_counts.data(),
      implementation.replicated_operator_displacements.data(), MPI_DOUBLE,
      implementation.communicator);
  agreed = consensus(
      implementation,
      mpi_status == MPI_SUCCESS
          ? Status{}
          : Status{StatusCode::mpi_failure, kMgCollective});
  if (!agreed) return agreed;

  std::fill(target.begin(), target.end(), 0.0);
  for (std::size_t global_index = 0U;
       global_index < implementation.replicated_global_cells;
       ++global_index) {
    const int owner = implementation.replicated_owner[global_index];
    const std::size_t local_index =
        implementation.replicated_owner_local[global_index];
    const std::size_t source =
        static_cast<std::size_t>(implementation.replicated_operator_displacements[
            owner]) +
        local_index * kReplicatedOperatorWidth;
    const std::size_t destination = global_index * kReplicatedOperatorWidth;
    for (std::size_t coefficient = 0U;
         coefficient < kReplicatedOperatorWidth; ++coefficient) {
      const double value =
          implementation.replicated_operator_gather[source + coefficient];
      if (!std::isfinite(value) || value < 0.0 ||
          (coefficient == 0U && !(value > 0.0))) {
        local = {StatusCode::numerical_failure, kMgCoefficient};
      }
      target[destination + coefficient] = value;
    }
  }
  agreed = consensus(implementation, local);
  return agreed;
}

template <class Implementation>
HaloEngine* halo_for(Implementation& implementation,
                     std::size_t level) noexcept {
  return level == 0U ? implementation.services.finest_halo
                     : implementation.services.level_halos.data[level - 1U];
}

template <class Implementation>
Status enter_prepared_halo_epoch(Implementation& implementation) noexcept {
  Status local{};
  const auto enter = [&](HaloEngine* halo) noexcept {
    if (halo == nullptr) {
      if (local) local = {StatusCode::invalid_plan, kMgPlan};
      return;
    }
    const Status entered = halo->enter_prepared_epoch();
    if (local && !entered) local = entered;
  };
  enter(implementation.services.finest_halo);
  for (HaloEngine* halo : implementation.level_halo_table) {
    enter(halo);
  }
  implementation.prepared_halo_epoch = static_cast<bool>(local);
  implementation.prepared_halo_deferred = {};
  return local;
}

template <class Implementation>
void close_prepared_halo_epoch(Implementation& implementation, bool publish,
                               int lowest_failing_rank) noexcept {
  if (implementation.services.finest_halo != nullptr) {
    implementation.services.finest_halo->close_prepared_epoch(
        publish, lowest_failing_rank);
  }
  for (HaloEngine* halo : implementation.level_halo_table) {
    if (halo != nullptr) {
      halo->close_prepared_epoch(publish, lowest_failing_rank);
    }
  }
  implementation.prepared_halo_epoch = false;
  implementation.prepared_halo_deferred = {};
}

template <class Implementation>
Status exchange_level_slot(Implementation& implementation, std::size_t level,
                           MgWorkspaceSlot slot, FieldView& field,
                           StageId stage) noexcept {
  Status local = implementation.services.workspace->revise_level(level, slot);
  field = implementation.services.workspace->level(level, slot);
  if (!local) return local;
  HaloTicket ticket;
  std::array<FieldView, 1U> views{field};
  HaloEngine* const halo = halo_for(implementation, level);
  local = implementation.prepared_halo_epoch
              ? halo->begin_prepared(stage, {views.data(), views.size()},
                                     implementation.prepared_halo_deferred,
                                     ticket)
              : halo->begin(stage, {views.data(), views.size()}, ticket);
  if (local) {
    local = implementation.prepared_halo_epoch
                ? halo->finish_prepared(
                      ticket, {views.data(), views.size()},
                      implementation.prepared_halo_deferred)
                : halo->finish(ticket, {views.data(), views.size()});
    field = views[0U];
  }
  if (!local) implementation.lowest = halo->lowest_failing_rank();
  return local;
}

template <class Implementation>
Status exchange_solution(Implementation& implementation, std::size_t level,
                         FieldView& solution, StageId stage) noexcept {
  return exchange_level_slot(implementation, level, MgWorkspaceSlot::solution,
                             solution, stage);
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

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
double level_matrix_value_reference(const detail::MgLevelStorage& level,
                          ConstFieldView input, Int3 cells, Int3 xe, Int3 ye,
                          Int3 ze, const double* diagonal, const double* x,
                          const double* y, const double* z,
                          const MgBoundarySet& boundary, Int3 cell) noexcept {
  const std::int32_t i = cell.x;
  const std::int32_t j = cell.y;
  const std::int32_t k = cell.z;
  double value = diagonal[flat(cells, cell)] * input.unchecked(cell, 0U);
  value -= x[face_flat(xe, {i, j, k})] *
           level_neighbor(input, level.patch, level.view.global_shape, cell,
                          CartesianAxis::x, -1, boundary.x_min,
                          boundary.x_max);
  value -= x[face_flat(xe, {i + 1, j, k})] *
           level_neighbor(input, level.patch, level.view.global_shape, cell,
                          CartesianAxis::x, 1, boundary.x_min,
                          boundary.x_max);
  value -= y[face_flat(ye, {i, j, k})] *
           level_neighbor(input, level.patch, level.view.global_shape, cell,
                          CartesianAxis::y, -1, boundary.y_min,
                          boundary.y_max);
  value -= y[face_flat(ye, {i, j + 1, k})] *
           level_neighbor(input, level.patch, level.view.global_shape, cell,
                          CartesianAxis::y, 1, boundary.y_min,
                          boundary.y_max);
  value -= z[face_flat(ze, {i, j, k})] *
           level_neighbor(input, level.patch, level.view.global_shape, cell,
                          CartesianAxis::z, -1, boundary.z_min,
                          boundary.z_max);
  value -= z[face_flat(ze, {i, j, k + 1})] *
           level_neighbor(input, level.patch, level.view.global_shape, cell,
                          CartesianAxis::z, 1, boundary.z_min,
                          boundary.z_max);
  return value;
}
#endif

struct PointRowView {
  const double* input{};
  const double* input_x_low{};
  const double* input_x_high{};
  const double* input_y_low{};
  const double* input_y_high{};
  const double* input_z_low{};
  const double* input_z_high{};
  const double* diagonal{};
  const double* x_coefficient{};
  const double* y_coefficient_low{};
  const double* y_coefficient_high{};
  const double* z_coefficient_low{};
  const double* z_coefficient_high{};
};

const double* point_row_boundary_pointer(
    const double* padded, const double* self, const double* zero,
    detail::MgPointBoundaryMode mode) noexcept {
  switch (mode) {
    case detail::MgPointBoundaryMode::padded_ghost:
      return padded;
    case detail::MgPointBoundaryMode::self:
      return self;
    case detail::MgPointBoundaryMode::zero:
      return zero;
  }
  return zero;
}

PointRowView point_row(ConstFieldView input, const double* diagonal,
                       const double* x, const double* y, const double* z,
                       Int3 cells, std::int32_t j, std::int32_t k,
                       const std::array<detail::MgPointBoundaryMode, 6U>&
                           boundary_modes,
                       const double* zero) noexcept {
  const std::size_t cell_row =
      (static_cast<std::size_t>(k) * static_cast<std::size_t>(cells.y) +
       static_cast<std::size_t>(j)) * static_cast<std::size_t>(cells.x);
  const std::size_t x_row =
      (static_cast<std::size_t>(k) * static_cast<std::size_t>(cells.y) +
       static_cast<std::size_t>(j)) * static_cast<std::size_t>(cells.x + 1);
  const std::size_t y_row =
      (static_cast<std::size_t>(k) *
           static_cast<std::size_t>(cells.y + 1) +
       static_cast<std::size_t>(j)) * static_cast<std::size_t>(cells.x);
  const std::size_t z_row = cell_row;
  const std::size_t z_next_row =
      (static_cast<std::size_t>(k + 1) * static_cast<std::size_t>(cells.y) +
       static_cast<std::size_t>(j)) * static_cast<std::size_t>(cells.x);
  const double* const input_row =
      input.base + static_cast<std::size_t>(k) * input.stride_z +
      static_cast<std::size_t>(j) * input.stride_y;
  return {
      input_row,
      point_row_boundary_pointer(
          input_row - static_cast<std::ptrdiff_t>(1U), input_row, zero,
          boundary_modes[0U]),
      point_row_boundary_pointer(
          input_row + static_cast<std::ptrdiff_t>(cells.x),
          input_row + static_cast<std::ptrdiff_t>(cells.x - 1), zero,
          boundary_modes[1U]),
      j > 0
          ? input_row - static_cast<std::ptrdiff_t>(input.stride_y)
          : point_row_boundary_pointer(
                input_row - static_cast<std::ptrdiff_t>(input.stride_y),
                input_row, zero, boundary_modes[2U]),
      j + 1 < cells.y
          ? input_row + static_cast<std::ptrdiff_t>(input.stride_y)
          : point_row_boundary_pointer(
                input_row + static_cast<std::ptrdiff_t>(input.stride_y),
                input_row, zero, boundary_modes[3U]),
      k > 0
          ? input_row - static_cast<std::ptrdiff_t>(input.stride_z)
          : point_row_boundary_pointer(
                input_row - static_cast<std::ptrdiff_t>(input.stride_z),
                input_row, zero, boundary_modes[4U]),
      k + 1 < cells.z
          ? input_row + static_cast<std::ptrdiff_t>(input.stride_z)
          : point_row_boundary_pointer(
                input_row + static_cast<std::ptrdiff_t>(input.stride_z),
                input_row, zero, boundary_modes[5U]),
      diagonal + cell_row,
      x + x_row,
      y + y_row,
      y + y_row + static_cast<std::size_t>(cells.x),
      z + z_row,
      z + z_next_row,
  };
}

double point_row_value(const PointRowView& row, std::int32_t nx,
                       std::int32_t i) noexcept {
  const std::size_t index = static_cast<std::size_t>(i);
  const double x_low =
      i > 0
          ? row.input[index - 1U]
          : row.input_x_low[0U];
  const double x_high =
      i + 1 < nx
          ? row.input[index + 1U]
          : row.input_x_high[0U];
  double value = row.diagonal[index] * row.input[index];
  value -= row.x_coefficient[index] * x_low;
  value -= row.x_coefficient[index + 1U] * x_high;
  value -= row.y_coefficient_low[index] * row.input_y_low[index];
  value -= row.y_coefficient_high[index] * row.input_y_high[index];
  value -= row.z_coefficient_low[index] * row.input_z_low[index];
  value -= row.z_coefficient_high[index] * row.input_z_high[index];
  return value;
}

template <bool Defect>
void point_row_bulk_kernel(const PointRowView& row, const double* rhs,
                           double* output, std::int32_t begin,
                           std::int32_t end) noexcept {
  // All neighbor rows are resolved before this loop.  The scalar callers
  // prove the input/rhs/output and four MG slots are disjoint; the kernel
  // deliberately does not use restrict because a Neumann neighbor row is
  // allowed to alias the current input row.
  for (std::int32_t i = begin; i < end; ++i) {
    const std::size_t index = static_cast<std::size_t>(i);
    double value = row.diagonal[index] * row.input[index];
    value -= row.x_coefficient[index] * row.input[index - 1U];
    value -= row.x_coefficient[index + 1U] * row.input[index + 1U];
    value -= row.y_coefficient_low[index] * row.input_y_low[index];
    value -= row.y_coefficient_high[index] * row.input_y_high[index];
    value -= row.z_coefficient_low[index] * row.input_z_low[index];
    value -= row.z_coefficient_high[index] * row.input_z_high[index];
    if constexpr (Defect) {
      output[index] = rhs[index] - value;
    } else {
      output[index] = value;
    }
  }
}

template <bool StageZero>
double chebyshev_stream_direction(const double* direction, std::size_t index,
                                  double previous_factor,
                                  double scaled) noexcept {
  if constexpr (StageZero) {
    static_cast<void>(direction);
    static_cast<void>(index);
    static_cast<void>(previous_factor);
    return scaled;
  }
  return previous_factor * direction[index] + scaled;
}

template <bool StageZero>
std::uint32_t chebyshev_stream_cell(
    const PointRowView& row, std::int32_t nx, const double* rhs,
    const double* direction, double* residual, double* next, std::size_t index,
    double previous_factor, double residual_factor) noexcept {
  const double value =
      point_row_value(row, nx, static_cast<std::int32_t>(index));
  const double defect = rhs[index] - value;
  const double scaled = residual_factor * defect / row.diagonal[index];
  const double next_direction = chebyshev_stream_direction<StageZero>(
      direction, index, previous_factor, scaled);
  const double next_solution = row.input[index] + next_direction;
  const std::uint32_t invalid_mask =
      static_cast<std::uint32_t>(!std::isfinite(scaled)) |
      static_cast<std::uint32_t>(!std::isfinite(next_direction)) |
      static_cast<std::uint32_t>(!std::isfinite(next_solution));
  residual[index] = next_direction;
  next[index] = next_solution;
  return invalid_mask;
}

template <bool StageZero>
std::uint32_t chebyshev_stream_bulk_kernel(
    const PointRowView& row, const double* rhs, const double* direction,
    double* residual, double* next, std::int32_t begin, std::int32_t end,
    double previous_factor, double residual_factor) noexcept {
  std::uint32_t invalid_mask = 0U;
  for (std::int32_t i = begin; i < end; ++i) {
    const std::size_t index = static_cast<std::size_t>(i);
    double value = row.diagonal[index] * row.input[index];
    value -= row.x_coefficient[index] * row.input[index - 1U];
    value -= row.x_coefficient[index + 1U] * row.input[index + 1U];
    value -= row.y_coefficient_low[index] * row.input_y_low[index];
    value -= row.y_coefficient_high[index] * row.input_y_high[index];
    value -= row.z_coefficient_low[index] * row.input_z_low[index];
    value -= row.z_coefficient_high[index] * row.input_z_high[index];
    const double defect = rhs[index] - value;
    const double scaled = residual_factor * defect / row.diagonal[index];
    const double next_direction = chebyshev_stream_direction<StageZero>(
        direction, index, previous_factor, scaled);
    const double next_solution = row.input[index] + next_direction;
    const std::uint32_t cell_invalid =
        static_cast<std::uint32_t>(!std::isfinite(scaled)) |
        static_cast<std::uint32_t>(!std::isfinite(next_direction)) |
        static_cast<std::uint32_t>(!std::isfinite(next_solution));
    residual[index] = next_direction;
    next[index] = next_solution;
    invalid_mask |= cell_invalid;
  }
  return invalid_mask;
}

template <bool StageZero>
std::uint32_t chebyshev_stream_row(
    const PointRowView& row, std::int32_t nx, const double* rhs,
    const double* direction, double* residual, double* next,
    double previous_factor, double residual_factor) noexcept {
  std::uint32_t invalid_mask = 0U;
  if (nx <= 2) {
    for (std::int32_t i = 0; i < nx; ++i) {
      invalid_mask |= chebyshev_stream_cell<StageZero>(
          row, nx, rhs, direction, residual, next,
          static_cast<std::size_t>(i), previous_factor, residual_factor);
    }
    return invalid_mask;
  }
  invalid_mask |= chebyshev_stream_cell<StageZero>(
      row, nx, rhs, direction, residual, next, 0U, previous_factor,
      residual_factor);
  invalid_mask |= chebyshev_stream_bulk_kernel<StageZero>(
      row, rhs, direction, residual, next, 1, nx - 1, previous_factor,
      residual_factor);
  invalid_mask |= chebyshev_stream_cell<StageZero>(
      row, nx, rhs, direction, residual, next,
      static_cast<std::size_t>(nx - 1), previous_factor, residual_factor);
  return invalid_mask;
}

template <bool StageZero, class Implementation>
std::uint32_t chebyshev_stream_stage(
    Implementation& implementation, const detail::MgLevelStorage& level,
    ConstFieldView current, ConstFieldView rhs, FieldView residual,
    FieldView next, const double* diagonal, const double* x_coefficient,
    const double* y_coefficient, const double* z_coefficient, Int3 cells,
    double previous_factor, double residual_factor) noexcept {
  std::uint32_t invalid_mask = 0U;
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      const PointRowView row =
          point_row(current, diagonal, x_coefficient, y_coefficient,
                    z_coefficient, cells, j, k, level.point_boundary_modes,
                    implementation.point_boundary_zero_storage.data());
      const double* const rhs_row =
          rhs.base + static_cast<std::size_t>(k) * rhs.stride_z +
          static_cast<std::size_t>(j) * rhs.stride_y;
      double* const next_row =
          next.base + static_cast<std::size_t>(k) * next.stride_z +
          static_cast<std::size_t>(j) * next.stride_y;
      double* const residual_row =
          residual.base + static_cast<std::size_t>(k) * residual.stride_z +
          static_cast<std::size_t>(j) * residual.stride_y;
      invalid_mask |= chebyshev_stream_row<StageZero>(
          row, cells.x, rhs_row, residual_row, residual_row, next_row,
          previous_factor, residual_factor);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
      const std::uint64_t row_cells = static_cast<std::uint64_t>(cells.x);
      implementation.matrix_work.defect_cell_visits += row_cells;
      implementation.matrix_work.chebyshev_stencil_evaluations += row_cells;
      implementation.matrix_work.chebyshev_updates += row_cells;
#endif
    }
  }
  return invalid_mask;
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
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  const Int3 xe{cells.x + 1, cells.y, cells.z};
  const Int3 ye{cells.x, cells.y + 1, cells.z};
  const Int3 ze{cells.x, cells.y, cells.z + 1};
#endif
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  ++implementation.matrix_work.full_actions;
  if (implementation.reference_point_actions) {
    const MgBoundarySet boundary = implementation.spec.boundaries;
    for (std::int32_t k = 0; k < cells.z; ++k) {
      for (std::int32_t j = 0; j < cells.y; ++j) {
        for (std::int32_t i = 0; i < cells.x; ++i) {
          const Int3 cell{i, j, k};
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
          ++implementation.matrix_work.full_cell_visits;
#endif
          output.unchecked(cell, 0U) = level_matrix_value_reference(
              level, input, cells, xe, ye, ze, diagonal, x, y, z, boundary,
              cell);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
          ++implementation.matrix_work.residual_writes;
#endif
        }
      }
    }
    return;
  }
#endif
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      const PointRowView row =
          point_row(input, diagonal, x, y, z, cells, j, k,
                    level.point_boundary_modes,
                    implementation.point_boundary_zero_storage.data());
      double* const output_row = output.base +
                                  static_cast<std::size_t>(k) * output.stride_z +
                                  static_cast<std::size_t>(j) * output.stride_y;
      if (cells.x <= 2) {
        for (std::int32_t i = 0; i < cells.x; ++i) {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
          ++implementation.matrix_work.full_cell_visits;
#endif
          output_row[static_cast<std::size_t>(i)] =
              point_row_value(row, cells.x, i);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
          ++implementation.matrix_work.residual_writes;
#endif
        }
      } else {
        output_row[0U] = point_row_value(row, cells.x, 0);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
        ++implementation.matrix_work.full_cell_visits;
        ++implementation.matrix_work.residual_writes;
#endif
        point_row_bulk_kernel<false>(row, nullptr, output_row, 1,
                                      cells.x - 1);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
        implementation.matrix_work.full_cell_visits +=
            static_cast<std::uint64_t>(cells.x - 2);
        implementation.matrix_work.residual_writes +=
            static_cast<std::uint64_t>(cells.x - 2);
#endif
        output_row[static_cast<std::size_t>(cells.x - 1)] =
            point_row_value(row, cells.x, cells.x - 1);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
        ++implementation.matrix_work.full_cell_visits;
        ++implementation.matrix_work.residual_writes;
#endif
      }
    }
  }
}

template <class Implementation>
void apply_level_defect(const Implementation& implementation,
                        std::size_t level_index, ConstFieldView input,
                        ConstFieldView rhs, FieldView output) noexcept {
  const detail::MgLevelStorage& level = implementation.levels[level_index];
  const double* const base = implementation.hierarchy_storage.data();
  const double* const diagonal = detail::block(base, level.diagonal_offset);
  const double* const x = detail::block(base, level.x_offset);
  const double* const y = detail::block(base, level.y_offset);
  const double* const z = detail::block(base, level.z_offset);
  const Int3 cells = level.view.local_shape;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  const Int3 xe{cells.x + 1, cells.y, cells.z};
  const Int3 ye{cells.x, cells.y + 1, cells.z};
  const Int3 ze{cells.x, cells.y, cells.z + 1};
  const bool chebyshev_stencil =
      level.point_smoother == MgPointSmootherKind::chebyshev_jacobi;
#endif
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (!chebyshev_stencil) {
    ++implementation.matrix_work.retained_final_defect_actions;
  }
  if (implementation.reference_point_actions) {
    const MgBoundarySet boundary = implementation.spec.boundaries;
    for (std::int32_t k = 0; k < cells.z; ++k) {
      for (std::int32_t j = 0; j < cells.y; ++j) {
        for (std::int32_t i = 0; i < cells.x; ++i) {
          const Int3 cell{i, j, k};
          const double value = level_matrix_value_reference(
              level, input, cells, xe, ye, ze, diagonal, x, y, z, boundary,
              cell);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
          ++implementation.matrix_work.defect_cell_visits;
          if (chebyshev_stencil) {
            ++implementation.matrix_work.chebyshev_stencil_evaluations;
          }
#endif
          output.unchecked(cell, 0U) = rhs.unchecked(cell, 0U) - value;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
          ++implementation.matrix_work.defect_writes;
#endif
        }
      }
    }
    return;
  }
#endif
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      const PointRowView row =
          point_row(input, diagonal, x, y, z, cells, j, k,
                    level.point_boundary_modes,
                    implementation.point_boundary_zero_storage.data());
      const double* const rhs_row = rhs.base +
                                     static_cast<std::size_t>(k) * rhs.stride_z +
                                     static_cast<std::size_t>(j) * rhs.stride_y;
      double* const output_row = output.base +
                                  static_cast<std::size_t>(k) * output.stride_z +
                                  static_cast<std::size_t>(j) * output.stride_y;
      if (cells.x <= 2) {
        for (std::int32_t i = 0; i < cells.x; ++i) {
          const std::size_t index = static_cast<std::size_t>(i);
          const double value = point_row_value(row, cells.x, i);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
          ++implementation.matrix_work.defect_cell_visits;
          if (chebyshev_stencil) {
            ++implementation.matrix_work.chebyshev_stencil_evaluations;
          }
#endif
          output_row[index] = rhs_row[index] - value;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
          ++implementation.matrix_work.defect_writes;
#endif
        }
      } else {
        const std::size_t first = 0U;
        output_row[first] = rhs_row[first] -
                            point_row_value(row, cells.x, 0);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
        ++implementation.matrix_work.defect_cell_visits;
        ++implementation.matrix_work.defect_writes;
        if (chebyshev_stencil) {
          ++implementation.matrix_work.chebyshev_stencil_evaluations;
        }
#endif
        point_row_bulk_kernel<true>(row, rhs_row, output_row, 1,
                                     cells.x - 1);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
        implementation.matrix_work.defect_cell_visits +=
            static_cast<std::uint64_t>(cells.x - 2);
        implementation.matrix_work.defect_writes +=
            static_cast<std::uint64_t>(cells.x - 2);
        if (chebyshev_stencil) {
          implementation.matrix_work.chebyshev_stencil_evaluations +=
              static_cast<std::uint64_t>(cells.x - 2);
        }
#endif
        const std::size_t last = static_cast<std::size_t>(cells.x - 1);
        output_row[last] = rhs_row[last] -
                            point_row_value(row, cells.x, cells.x - 1);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
        ++implementation.matrix_work.defect_cell_visits;
        ++implementation.matrix_work.defect_writes;
        if (chebyshev_stencil) {
          ++implementation.matrix_work.chebyshev_stencil_evaluations;
        }
#endif
      }
    }
  }
}

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
// The original defect-publication/update lifecycle is retained only as a
// test oracle.  Production objects do not contain this path.
template <class Implementation>
Status chebyshev_point_smooth_legacy(Implementation& implementation,
                                     std::size_t level_index,
                                     std::uint32_t stages,
                                     bool retain_final_defect,
                                     StageId& stage,
                                     Status& deferred) noexcept {
  FieldView x = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::solution);
  FieldView rhs = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::rhs);
  FieldView residual = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::residual);
  FieldView direction = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::temporary);
  const detail::MgLevelStorage& level = implementation.levels[level_index];
  const double* const diagonal = detail::block(
      implementation.hierarchy_storage.data(), level.diagonal_offset);
  const Int3 cells = level.view.local_shape;
  constexpr double beta = 2.0;
  const double alpha =
      level.chebyshev_lower_spectrum_fraction * beta;
  const double theta = (beta + alpha) * 0.5;
  const double delta = (beta - alpha) * 0.5;
  const double sigma = theta / delta;
  double previous_rho = 1.0 / sigma;
  bool finite = std::isfinite(alpha) && std::isfinite(theta) &&
                std::isfinite(delta) && delta > 0.0 &&
                std::isfinite(sigma) && sigma > 0.0 &&
                std::isfinite(previous_rho);
  const auto consume_stage = [&](StageId& ordinal) noexcept {
    if (stage == std::numeric_limits<StageId>::max()) {
      return false;
    }
    ordinal = stage;
    ++stage;
    return true;
  };
  const auto update_cells = [&](double previous_factor,
                                double residual_factor) noexcept {
    for (std::int32_t k = 0; k < cells.z; ++k) {
      for (std::int32_t j = 0; j < cells.y; ++j) {
        double* const x_row = x.base +
                              static_cast<std::size_t>(k) * x.stride_z +
                              static_cast<std::size_t>(j) * x.stride_y;
        double* const direction_row =
            direction.base + static_cast<std::size_t>(k) * direction.stride_z +
            static_cast<std::size_t>(j) * direction.stride_y;
        const double* const residual_row =
            residual.base + static_cast<std::size_t>(k) * residual.stride_z +
            static_cast<std::size_t>(j) * residual.stride_y;
        const double* const diagonal_row =
            diagonal +
            (static_cast<std::size_t>(k) * static_cast<std::size_t>(cells.y) +
             static_cast<std::size_t>(j)) * static_cast<std::size_t>(cells.x);
        for (std::int32_t i = 0; i < cells.x; ++i) {
          const double scaled = residual_factor *
                                residual_row[static_cast<std::size_t>(i)] /
                                diagonal_row[static_cast<std::size_t>(i)];
          const double next_direction =
              previous_factor * direction_row[static_cast<std::size_t>(i)] +
              scaled;
          const double next_solution = x_row[static_cast<std::size_t>(i)] +
                                       next_direction;
          finite = finite && std::isfinite(scaled) &&
                   std::isfinite(next_direction) &&
                   std::isfinite(next_solution);
          direction_row[static_cast<std::size_t>(i)] = next_direction;
          x_row[static_cast<std::size_t>(i)] = next_solution;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
          ++implementation.matrix_work.chebyshev_updates;
#endif
        }
      }
    }
  };

  if (stages == 0U) {
    return {StatusCode::invalid_plan, kMgApply};
  }

  for (std::uint32_t polynomial_stage = 0U;
       polynomial_stage < stages; ++polynomial_stage) {
    StageId ordinal = 0U;
    if (!consume_stage(ordinal)) {
      return {StatusCode::invalid_plan, kMgApply};
    }
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
    ++implementation.matrix_work.chebyshev_stages;
    ++implementation.matrix_work.chebyshev_exchange_actions;
#endif
    Status status = exchange_solution(implementation, level_index, x, ordinal);
    if (!status) {
      return status;
    }
    apply_level_defect(implementation, level_index, as_const(x),
                       as_const(rhs), residual);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
    ++implementation.matrix_work.chebyshev_defect_actions;
#endif
    if (polynomial_stage == 0U) {
      update_cells(0.0, 1.0 / theta);
    } else {
      const double rho = 1.0 / (2.0 * sigma - previous_rho);
      const double previous_factor = rho * previous_rho;
      const double residual_factor = 2.0 * rho / delta;
      finite = finite && std::isfinite(rho) && std::isfinite(previous_factor) &&
               std::isfinite(residual_factor);
      update_cells(previous_factor, residual_factor);
      previous_rho = rho;
    }
  }
  Status status{};
  if (retain_final_defect) {
    StageId ordinal = 0U;
    if (!consume_stage(ordinal)) {
      return {StatusCode::invalid_plan, kMgApply};
    }
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
    ++implementation.matrix_work.chebyshev_exchange_actions;
    ++implementation.matrix_work.chebyshev_retained_final_defect_actions;
#endif
    status = exchange_solution(implementation, level_index, x, ordinal);
    if (status) {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
      ++implementation.matrix_work.retained_final_defect_actions;
#endif
      apply_level_defect(implementation, level_index, as_const(x),
                         as_const(rhs), residual);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
      ++implementation.matrix_work.chebyshev_defect_actions;
#endif
    }
  }
  if (status) status = implementation.services.workspace->revise_level(
      level_index, MgWorkspaceSlot::solution);
  if (status) status = implementation.services.workspace->revise_level(
      level_index, MgWorkspaceSlot::temporary);
  if (status) status = implementation.services.workspace->revise_level(
      level_index, MgWorkspaceSlot::residual);
  if (!finite && status && deferred) {
    deferred = {StatusCode::numerical_failure, kMgApply};
  }
  return status;
}

#endif

template <class Implementation>
Status chebyshev_point_smooth_streamed(Implementation& implementation,
                                       std::size_t level_index,
                                       std::uint32_t stages,
                                       bool retain_final_defect,
                                       StageId& stage,
                                       Status& deferred) noexcept {
  FieldView solution = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::solution);
  FieldView temporary = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::temporary);
  FieldView rhs = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::rhs);
  FieldView residual = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::residual);
  const detail::MgLevelStorage& level = implementation.levels[level_index];
  const double* const base = implementation.hierarchy_storage.data();
  const double* const diagonal = detail::block(base, level.diagonal_offset);
  const double* const x_coefficient = detail::block(base, level.x_offset);
  const double* const y_coefficient = detail::block(base, level.y_offset);
  const double* const z_coefficient = detail::block(base, level.z_offset);
  const Int3 cells = level.view.local_shape;
  constexpr double beta = 2.0;
  const double alpha =
      level.chebyshev_lower_spectrum_fraction * beta;
  const double theta = (beta + alpha) * 0.5;
  const double delta = (beta - alpha) * 0.5;
  const double sigma = theta / delta;
  double previous_rho = 1.0 / sigma;
  // This integer mask is status-only.  It does not guard, reorder or alter
  // any FP operation in the streamed cell recurrence.
  std::uint32_t invalid_mask =
      (std::isfinite(alpha) && std::isfinite(theta) &&
       std::isfinite(delta) && delta > 0.0 && std::isfinite(sigma) &&
       sigma > 0.0 && std::isfinite(previous_rho))
          ? 0U
          : 1U;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  invalid_mask |=
      static_cast<std::uint32_t>(implementation.force_chebyshev_invalid);
#endif
  const auto consume_stage = [&](StageId& ordinal) noexcept {
    if (stage == std::numeric_limits<StageId>::max()) {
      return false;
    }
    ordinal = stage;
    ++stage;
    return true;
  };

  if (stages == 0U) {
    return {StatusCode::invalid_plan, kMgApply};
  }

  MgWorkspaceSlot current_slot = MgWorkspaceSlot::solution;
  FieldView current = solution;
  for (std::uint32_t polynomial_stage = 0U;
       polynomial_stage < stages; ++polynomial_stage) {
    StageId ordinal = 0U;
    if (!consume_stage(ordinal)) {
      return {StatusCode::invalid_plan, kMgApply};
    }
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
    ++implementation.matrix_work.chebyshev_stages;
    ++implementation.matrix_work.chebyshev_exchange_actions;
#endif
    Status status = exchange_level_slot(implementation, level_index,
                                        current_slot, current, ordinal);
    if (!status) {
      return status;
    }

    const MgWorkspaceSlot next_slot =
        current_slot == MgWorkspaceSlot::solution
            ? MgWorkspaceSlot::temporary
            : MgWorkspaceSlot::solution;
    FieldView next = implementation.services.workspace->level(level_index,
                                                               next_slot);
    const double rho = polynomial_stage == 0U
                           ? 0.0
                           : 1.0 / (2.0 * sigma - previous_rho);
    const double previous_factor =
        polynomial_stage == 0U ? 0.0 : rho * previous_rho;
    const double residual_factor =
        polynomial_stage == 0U ? 1.0 / theta : 2.0 * rho / delta;
    if (polynomial_stage != 0U) {
      invalid_mask |= static_cast<std::uint32_t>(!std::isfinite(rho));
      invalid_mask |=
          static_cast<std::uint32_t>(!std::isfinite(previous_factor));
      invalid_mask |=
          static_cast<std::uint32_t>(!std::isfinite(residual_factor));
      previous_rho = rho;
    }

    if (polynomial_stage == 0U) {
      invalid_mask |= chebyshev_stream_stage<true>(
          implementation, level, as_const(current), as_const(rhs), residual,
          next, diagonal, x_coefficient, y_coefficient, z_coefficient, cells,
          previous_factor, residual_factor);
    } else {
      invalid_mask |= chebyshev_stream_stage<false>(
          implementation, level, as_const(current), as_const(rhs), residual,
          next, diagonal, x_coefficient, y_coefficient, z_coefficient, cells,
          previous_factor, residual_factor);
    }
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
    // The old route published rhs-A*x as a complete residual field before
    // every polynomial update.  The streamed route consumes that value in
    // the same cell traversal and publishes only the recurrence direction.
    implementation.matrix_work.chebyshev_defect_actions++;
    implementation.matrix_work.chebyshev_elided_intermediate_residual_publications +=
        static_cast<std::uint64_t>(cells.x) *
        static_cast<std::uint64_t>(cells.y) *
        static_cast<std::uint64_t>(cells.z);
#endif
    current_slot = next_slot;
    current = next;
  }

  // The odd-degree result lives in temporary.  Only the interior is copied;
  // the authoritative solution halo is rebuilt by the retained exchange or
  // by the next smoother stage.
  if (current_slot == MgWorkspaceSlot::temporary) {
    for (std::int32_t k = 0; k < cells.z; ++k) {
      for (std::int32_t j = 0; j < cells.y; ++j) {
        const double* const temporary_row =
            temporary.base + static_cast<std::size_t>(k) * temporary.stride_z +
            static_cast<std::size_t>(j) * temporary.stride_y;
        double* const solution_row =
            solution.base + static_cast<std::size_t>(k) * solution.stride_z +
            static_cast<std::size_t>(j) * solution.stride_y;
        for (std::int32_t i = 0; i < cells.x; ++i) {
          solution_row[static_cast<std::size_t>(i)] =
              temporary_row[static_cast<std::size_t>(i)];
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
          ++implementation.matrix_work.chebyshev_copyback_cells;
#endif
        }
      }
    }
  }

  Status status{};
  if (retain_final_defect) {
    StageId ordinal = 0U;
    if (!consume_stage(ordinal)) {
      return {StatusCode::invalid_plan, kMgApply};
    }
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
    ++implementation.matrix_work.chebyshev_exchange_actions;
    ++implementation.matrix_work.chebyshev_retained_final_defect_actions;
#endif
    status = exchange_level_slot(implementation, level_index,
                                 MgWorkspaceSlot::solution, solution, ordinal);
    if (status) {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
      ++implementation.matrix_work.retained_final_defect_actions;
#endif
      apply_level_defect(implementation, level_index, as_const(solution),
                         as_const(rhs), residual);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
      ++implementation.matrix_work.chebyshev_defect_actions;
#endif
    }
  }
  if (status) status = implementation.services.workspace->revise_level(
      level_index, MgWorkspaceSlot::solution);
  if (status) status = implementation.services.workspace->revise_level(
      level_index, MgWorkspaceSlot::temporary);
  if (status) status = implementation.services.workspace->revise_level(
      level_index, MgWorkspaceSlot::residual);
  if (invalid_mask != 0U && status && deferred) {
    deferred = {StatusCode::numerical_failure, kMgApply};
  }
  return status;
}

template <class Implementation>
Status chebyshev_point_smooth(Implementation& implementation,
                              std::size_t level_index, std::uint32_t stages,
                              bool retain_final_defect,
                              StageId& stage,
                              Status& deferred) noexcept {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (implementation.reference_point_actions ||
      implementation.reference_chebyshev_lifecycle) {
    return chebyshev_point_smooth_legacy(implementation, level_index, stages,
                                         retain_final_defect, stage,
                                         deferred);
  }
#endif
  return chebyshev_point_smooth_streamed(implementation, level_index, stages,
                                         retain_final_defect, stage,
                                         deferred);
}

template <class Implementation>
bool fused_point_update_allowed(const Implementation& implementation,
                                std::size_t level_index) noexcept {
  const detail::MgLevelStorage& level = implementation.levels[level_index];
  const Int3 global = level.view.global_shape;
  const MgBoundarySet boundary = implementation.spec.boundaries;
  const bool x_periodic = boundary.x_min == MgBoundaryKind::periodic ||
                          boundary.x_max == MgBoundaryKind::periodic;
  const bool y_periodic = boundary.y_min == MgBoundaryKind::periodic ||
                          boundary.y_max == MgBoundaryKind::periodic;
  const bool z_periodic = boundary.z_min == MgBoundaryKind::periodic ||
                          boundary.z_max == MgBoundaryKind::periodic;
  return (!x_periodic || (global.x & 1) == 0) &&
         (!y_periodic || (global.y & 1) == 0) &&
         (!z_periodic || (global.z & 1) == 0);
}

template <class Implementation>
bool fused_point_actions_enabled(const Implementation& implementation,
                                 std::size_t level_index) noexcept {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (implementation.reference_point_actions) return false;
#endif
  return fused_point_update_allowed(implementation, level_index);
}

template <class Implementation>
bool point_retained_defect_enabled(
    const Implementation& implementation) noexcept {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  return !implementation.reference_point_actions;
#else
  static_cast<void>(implementation);
  return true;
#endif
}

template <class Implementation>
void fused_point_color_update(Implementation& implementation,
                              std::size_t level_index, FieldView x,
                              FieldView rhs, const double* diagonal,
                              Int3 cells, int color, bool reverse,
                              double omega) noexcept {
  const detail::MgLevelStorage& level = implementation.levels[level_index];
  const double* const base = implementation.hierarchy_storage.data();
  const double* const x_coefficient = detail::block(base, level.x_offset);
  const double* const y_coefficient = detail::block(base, level.y_offset);
  const double* const z_coefficient = detail::block(base, level.z_offset);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  const Int3 xe{cells.x + 1, cells.y, cells.z};
  const Int3 ye{cells.x, cells.y + 1, cells.z};
  const Int3 ze{cells.x, cells.y, cells.z + 1};
  const MgBoundarySet boundary = implementation.spec.boundaries;
#endif
  const ConstFieldView input = as_const(x);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  ++implementation.matrix_work.fused_color_actions;
#endif
  const auto update_row = [&](std::int32_t k, std::int32_t j,
                              bool reverse_i) noexcept {
    const PointRowView row =
        point_row(input, diagonal, x_coefficient, y_coefficient,
                  z_coefficient, cells, j, k, level.point_boundary_modes,
                  implementation.point_boundary_zero_storage.data());
    double* const x_row = x.base + static_cast<std::size_t>(k) * x.stride_z +
                          static_cast<std::size_t>(j) * x.stride_y;
    const double* const rhs_row =
        rhs.base + static_cast<std::size_t>(k) * rhs.stride_z +
        static_cast<std::size_t>(j) * rhs.stride_y;
    const auto update = [&](std::int32_t i) noexcept {
      const std::size_t index = static_cast<std::size_t>(i);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
      ++implementation.matrix_work.fused_cell_visits;
      const double value =
          implementation.reference_point_actions
              ? level_matrix_value_reference(
                    level, input, cells, xe, ye, ze, diagonal,
                    x_coefficient, y_coefficient, z_coefficient, boundary,
                    {i, j, k})
              : point_row_value(row, cells.x, i);
#else
      const double value = point_row_value(row, cells.x, i);
#endif
      x_row[index] += omega * (rhs_row[index] - value) / row.diagonal[index];
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
      ++implementation.matrix_work.fused_updates;
#endif
    };
    if (reverse_i) {
      for (std::int32_t i = cells.x; i-- > 0;) {
        if (((i + j + k) & 1) == color) update(i);
      }
    } else {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        if (((i + j + k) & 1) == color) update(i);
      }
    }
  };
  if (reverse) {
    for (std::int32_t k = cells.z; k-- > 0;) {
      for (std::int32_t j = cells.y; j-- > 0;) {
        update_row(k, j, true);
      }
    }
  } else {
    for (std::int32_t k = 0; k < cells.z; ++k) {
      for (std::int32_t j = 0; j < cells.y; ++j) {
        update_row(k, j, false);
      }
    }
  }
}

template <class Implementation>
Status point_smooth(Implementation& implementation, std::size_t level_index,
                    std::uint32_t sweeps, bool reverse,
                    bool retain_final_matrix_action,
                    StageId& stage,
                    Status& deferred) noexcept {
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
  if (level.point_smoother == MgPointSmootherKind::chebyshev_jacobi) {
    return chebyshev_point_smooth(implementation, level_index, sweeps,
                                  retain_final_matrix_action, stage,
                                  deferred);
  }
  constexpr double omega = 0.72;
  const bool fused = fused_point_actions_enabled(implementation, level_index);
  const bool direct_defect = retain_final_matrix_action &&
                             point_retained_defect_enabled(implementation);
  Status status = exchange_solution(implementation, level_index, x, stage++);
  if (!status) return status;
  if (!fused || sweeps == 0U) {
    if (direct_defect && sweeps == 0U) {
      apply_level_defect(implementation, level_index, as_const(x),
                         as_const(rhs), residual);
    } else {
      apply_level_matrix(implementation, level_index, as_const(x), residual);
    }
  }
  for (std::uint32_t sweep = 0U; sweep < sweeps; ++sweep) {
    const int first = reverse ? 1 : 0;
    for (int pass = 0; pass < 2; ++pass) {
      const int color = (first + pass) & 1;
      if (!fused) {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
        ++implementation.matrix_work.separated_color_actions;
#endif
        const auto update_row = [&](std::int32_t k, std::int32_t j,
                                    bool reverse_i) noexcept {
          double* const x_row =
              x.base + static_cast<std::size_t>(k) * x.stride_z +
              static_cast<std::size_t>(j) * x.stride_y;
          const double* const rhs_row =
              rhs.base + static_cast<std::size_t>(k) * rhs.stride_z +
              static_cast<std::size_t>(j) * rhs.stride_y;
          const double* const residual_row =
              residual.base + static_cast<std::size_t>(k) * residual.stride_z +
              static_cast<std::size_t>(j) * residual.stride_y;
          const double* const diagonal_row =
              diagonal +
              (static_cast<std::size_t>(k) *
                   static_cast<std::size_t>(cells.y) +
               static_cast<std::size_t>(j)) *
                  static_cast<std::size_t>(cells.x);
          const auto update = [&](std::int32_t i) noexcept {
            const std::size_t index = static_cast<std::size_t>(i);
            x_row[index] +=
                omega * (rhs_row[index] - residual_row[index]) /
                diagonal_row[index];
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
            ++implementation.matrix_work.separated_color_updates;
#endif
          };
          if (reverse_i) {
            for (std::int32_t i = cells.x; i-- > 0;) {
              if (((i + j + k) & 1) == color) update(i);
            }
          } else {
            for (std::int32_t i = 0; i < cells.x; ++i) {
              if (((i + j + k) & 1) == color) update(i);
            }
          }
        };
        if (reverse) {
          for (std::int32_t k = cells.z; k-- > 0;) {
            for (std::int32_t j = cells.y; j-- > 0;) {
              update_row(k, j, true);
            }
          }
        } else {
          for (std::int32_t k = 0; k < cells.z; ++k) {
            for (std::int32_t j = 0; j < cells.y; ++j) {
              update_row(k, j, false);
            }
          }
        }
      } else {
        fused_point_color_update(implementation, level_index, x, rhs,
                                 diagonal, cells, color, reverse, omega);
      }
      const bool final_update = sweep + 1U == sweeps && pass == 1;
      if (!final_update || retain_final_matrix_action) {
        status = exchange_solution(implementation, level_index, x, stage++);
        if (!status) return status;
        if (!fused || (final_update && retain_final_matrix_action)) {
          if (direct_defect && final_update) {
            apply_level_defect(implementation, level_index, as_const(x),
                               as_const(rhs), residual);
          } else {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
            if (final_update && retain_final_matrix_action) {
              ++implementation.matrix_work.retained_final_full_actions;
            }
#endif
            apply_level_matrix(implementation, level_index, as_const(x),
                               residual);
          }
        }
      }
    }
  }
  if (!retain_final_matrix_action) {
    status = implementation.services.workspace->revise_level(
        level_index, MgWorkspaceSlot::solution);
    if (!status) return status;
  }
  status = implementation.services.workspace->revise_level(
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
    return point_smooth(implementation, level_index, sweeps, reverse, false,
                        stage, deferred);
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
    return point_smooth(implementation, level_index, sweeps, reverse, false,
                        stage, deferred);
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
              std::uint32_t sweeps, bool reverse,
              bool retain_final_matrix_action, StageId& stage,
              Status& deferred) noexcept {
  const std::uint8_t mask = implementation.levels[level].view.line_axis_mask;
  if (mask == 0U) {
    return point_smooth(implementation, level, sweeps, reverse,
                        retain_final_matrix_action, stage, deferred);
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
Status finish_residual(Implementation& implementation,
                       std::size_t level_index) noexcept {
  FieldView rhs = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::rhs);
  FieldView residual = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::residual);
  const Int3 cells = residual.interior;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  ++implementation.matrix_work.residual_finish_actions;
#endif
  for (std::int32_t k = 0; k < cells.z; ++k)
    for (std::int32_t j = 0; j < cells.y; ++j)
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 cell{i, j, k};
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
        ++implementation.matrix_work.residual_finish_cell_visits;
#endif
        residual.unchecked(cell, 0U) = rhs.unchecked(cell, 0U) -
                                       residual.unchecked(cell, 0U);
      }
  return implementation.services.workspace->revise_level(
      level_index, MgWorkspaceSlot::residual);
}

template <class Implementation>
Status compute_residual(Implementation& implementation,
                        std::size_t level_index, StageId& stage) noexcept {
  FieldView x = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::solution);
  FieldView residual = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::residual);
  Status status = exchange_solution(implementation, level_index, x, stage++);
  if (!status) return status;
  apply_level_matrix(implementation, level_index, as_const(x), residual);
  return finish_residual(implementation, level_index);
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
  HaloEngine* const halo = halo_for(implementation, fine_index);
  status = implementation.prepared_halo_epoch
               ? halo->begin_prepared(
                     stage++, {views.data(), views.size()},
                     implementation.prepared_halo_deferred, ticket)
               : halo->begin(stage++, {views.data(), views.size()}, ticket);
  if (status) {
    status = implementation.prepared_halo_epoch
                 ? halo->finish_prepared(
                       ticket, {views.data(), views.size()},
                       implementation.prepared_halo_deferred)
                 : halo->finish(ticket, {views.data(), views.size()});
    fr = views[0U];
  }
  if (!status) implementation.lowest = halo->lowest_failing_rank();
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
  for (std::int32_t ck = 0; ck < cs.z; ++ck) {
    for (std::int32_t cj = 0; cj < cs.y; ++cj) {
      for (std::int32_t ci = 0; ci < cs.x; ++ci) {
        const Int3 cg{coarse.patch.begin.x + ci,
                      coarse.patch.begin.y + cj,
                      coarse.patch.begin.z + ck};
        const std::int32_t xb = mx ? 2 * cg.x : cg.x;
        const std::int32_t yb = my ? 2 * cg.y : cg.y;
        const std::int32_t zb = mz ? 2 * cg.z : cg.z;
        double residual_sum = 0.0;
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
              residual_sum += ghost_or_local(as_const(fr), local);
            }
        // Pressure/momentum coefficients and right-hand sides are stored as
        // finite-volume integrated equations.  The coarse conservative row
        // is the sum of its child rows, so its defect must use the same sum;
        // a volume average is dimensionally inconsistent with the summed
        // reaction and face coefficients and suppresses coarse correction.
        crhs.unchecked({ci, cj, ck}, 0U) = residual_sum;
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
Status replicated_coarse_solve(Implementation& implementation,
                               std::size_t level_index,
                               std::uint32_t sweeps, bool reverse,
                               StageId& stage, Status& deferred) noexcept {
  if (!implementation.replicated_coarse ||
      level_index + 1U != implementation.levels.size()) {
    return {StatusCode::invalid_plan, kMgPlan};
  }
  const detail::MgLevelStorage& coarse = implementation.levels[level_index];
  const Int3 local_cells = coarse.view.local_shape;
  FieldView rhs = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::rhs);
  FieldView solution = implementation.services.workspace->level(
      level_index, MgWorkspaceSlot::solution);
  Status local = rhs.base == nullptr || solution.base == nullptr ||
                         rhs.interior.x != local_cells.x ||
                         rhs.interior.y != local_cells.y ||
                         rhs.interior.z != local_cells.z
                     ? Status{StatusCode::invalid_plan, kMgApply}
                     : Status{};
  if (!local && deferred) deferred = local;
  const bool valid_workspace = static_cast<bool>(local);
  const bool continue_solution =
      valid_workspace && implementation.replicated_solution_revision != 0U &&
      solution.revision == implementation.replicated_solution_revision;

  for (std::int32_t k = 0; k < local_cells.z; ++k) {
    for (std::int32_t j = 0; j < local_cells.y; ++j) {
      for (std::int32_t i = 0; i < local_cells.x; ++i) {
        const std::size_t local_index = flat(local_cells, {i, j, k});
        implementation.replicated_rhs_local[local_index] =
            valid_workspace ? rhs.unchecked({i, j, k}, 0U) : 0.0;
      }
    }
  }
  local = record_replicated_collective(
      implementation,
      implementation.replicated_global_cells * sizeof(double));
  if (!local && deferred) deferred = local;
  const int mpi_status = MPI_Allgatherv(
      implementation.replicated_rhs_local.data(),
      implementation.replicated_rhs_counts[implementation.rank], MPI_DOUBLE,
      implementation.replicated_rhs_gather.data(),
      implementation.replicated_rhs_counts.data(),
      implementation.replicated_rhs_displacements.data(), MPI_DOUBLE,
      implementation.communicator);
  if (mpi_status != MPI_SUCCESS && deferred) {
    deferred = {StatusCode::mpi_failure, kMgCollective};
  }

  for (std::size_t global_index = 0U;
       global_index < implementation.replicated_global_cells;
       ++global_index) {
    const int owner = implementation.replicated_owner[global_index];
    const std::size_t local_index =
        implementation.replicated_owner_local[global_index];
    implementation.replicated_rhs_global[global_index] =
        mpi_status == MPI_SUCCESS
            ? implementation.replicated_rhs_gather[
                  static_cast<std::size_t>(
                      implementation.replicated_rhs_displacements[owner]) +
                  local_index]
            : 0.0;
    if (!continue_solution) {
      implementation.replicated_solution[global_index] = 0.0;
    }
  }

  const Int3 global = implementation.replicated_global_shape;
  const auto matrix_value = [&](std::size_t index,
                                const double* input) noexcept -> double {
    const double* const row = implementation.replicated_operator_active.data() +
                              index * kReplicatedOperatorWidth;
    const std::size_t* const neighbors =
        implementation.replicated_neighbors.data() + index * 6U;
    double value = row[0U] * input[index];
    value -= row[1U] * input[neighbors[0U]];
    value -= row[2U] * input[neighbors[1U]];
    value -= row[3U] * input[neighbors[2U]];
    value -= row[4U] * input[neighbors[3U]];
    value -= row[5U] * input[neighbors[4U]];
    value -= row[6U] * input[neighbors[5U]];
    return value;
  };

  local = {};
  for (std::uint32_t sweep = 0U; sweep < sweeps; ++sweep) {
    for (int pass = 0; pass < 2; ++pass) {
      const int color = (reverse ? 1 : 0) ^ (pass & 1);
      std::copy(implementation.replicated_solution.begin(),
                implementation.replicated_solution.end(),
                implementation.replicated_snapshot.begin());
      for (std::size_t order = 0U;
           order < implementation.replicated_global_cells; ++order) {
        const std::size_t global_index =
            implementation.replicated_owner_order[order];
        if (implementation.replicated_owner_parity[global_index] != color) {
          continue;
        }
        const double* const row =
            implementation.replicated_operator_active.data() +
            global_index * kReplicatedOperatorWidth;
        const double action = matrix_value(
            global_index, implementation.replicated_snapshot.data());
        implementation.replicated_action[global_index] = action;
        const double update =
            0.72 *
            (implementation.replicated_rhs_global[global_index] - action) /
            row[0U];
        implementation.replicated_solution[global_index] += update;
      }
    }
  }
  // Keep the same stage namespace as the reference coarse point smoother:
  // one initial exchange plus one exchange after every non-final color.
  const std::uint64_t exchanges =
      sweeps == 0U ? 1U : 1U + 2U * static_cast<std::uint64_t>(sweeps) - 1U;
  if (exchanges > std::numeric_limits<StageId>::max() - stage) {
    local = {StatusCode::invalid_plan, kMgApply};
  } else {
    stage = static_cast<StageId>(stage + exchanges);
  }
  if (!local && deferred) deferred = local;

  for (std::int32_t k = 0; k < local_cells.z && deferred; ++k) {
    for (std::int32_t j = 0; j < local_cells.y; ++j) {
      for (std::int32_t i = 0; i < local_cells.x; ++i) {
        const std::size_t local_index = flat(local_cells, {i, j, k});
        const std::size_t global_index = flat(
            global, {coarse.patch.begin.x + i, coarse.patch.begin.y + j,
                     coarse.patch.begin.z + k});
        if (valid_workspace) {
          solution.unchecked({i, j, k}, 0U) =
              implementation.replicated_solution[global_index];
        }
      }
    }
  }
  if (deferred && valid_workspace) {
    local = implementation.services.workspace->revise_level(
        level_index, MgWorkspaceSlot::solution);
    if (!local && deferred) {
      deferred = local;
    } else if (local) {
      solution = implementation.services.workspace->level(
          level_index, MgWorkspaceSlot::solution);
      implementation.replicated_solution_revision = solution.revision;
    }
  }
  if (!deferred || !valid_workspace) {
    implementation.replicated_solution_revision = 0U;
  }
  return {};
}

template <class Implementation>
bool replicated_coarse_route_enabled(
    const Implementation& implementation) noexcept {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (implementation.reference_replicated_coarse) return false;
  if (implementation.request_replicated_coarse) {
    return implementation.replicated_coarse;
  }
#endif
  return implementation.replicated_coarse;
}

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
template <class Implementation>
Status inject_cycle_failure(
    Implementation& implementation, detail::MgCycleFailurePhase phase,
    std::size_t level, Status& deferred) noexcept {
  if (implementation.cycle_failure_consumed ||
      implementation.cycle_failure_phase != phase ||
      implementation.cycle_failure_level != level) {
    return {};
  }
  implementation.cycle_failure_consumed = true;
  const Status local = implementation.rank == implementation.cycle_failure_rank
                           ? Status{implementation.cycle_failure_deferred
                                        ? StatusCode::numerical_failure
                                        : StatusCode::invalid_plan,
                                    kMgApply}
                           : Status{};
  if (implementation.cycle_failure_deferred) {
    if (!local && deferred) deferred = local;
    return {};
  }
  return consensus(implementation, local);
}
#endif

template <class Implementation>
Status mg_cycle(Implementation& implementation, std::size_t level,
                std::uint8_t cycle_counter, StageId& stage,
                Status& deferred) noexcept {
  if ((cycle_counter != 1U && cycle_counter != 2U) ||
      level >= implementation.levels.size()) {
    return {StatusCode::invalid_plan, kMgApply};
  }
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  ++implementation.matrix_work.cycle_level_calls[level];
  Status injected = inject_cycle_failure(
      implementation, detail::MgCycleFailurePhase::pre_smooth, level,
      deferred);
  if (!injected) return injected;
#endif
  const bool coarse = level + 1U == implementation.levels.size();
  if (coarse) {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
    injected = inject_cycle_failure(
        implementation, detail::MgCycleFailurePhase::terminal, level,
        deferred);
    if (!injected) return injected;
#endif
    if (replicated_coarse_route_enabled(implementation)) {
      return replicated_coarse_solve(
          implementation, level, implementation.spec.policy.coarse_sweeps,
          false, stage, deferred);
    }
    return smooth(implementation, level,
                  implementation.spec.policy.coarse_sweeps, false, false,
                  stage, deferred);
  }
  Status status = smooth(implementation, level,
                         implementation.spec.policy.pre_sweeps, false, true,
                         stage, deferred);
  if (status) {
    if (implementation.levels[level].view.line_axis_mask == 0U) {
      // Every production point path writes rhs - A*x in its retained final
      // stencil action.  The test-only legacy reference retains A*x and
      // exercises this conversion pass for the oracle.
      // Chebyshev's retained action is already a defect on both routes; the
      // legacy conversion is only for the red/black matrix oracle.
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
      const bool chebyshev =
          implementation.levels[level].point_smoother ==
          MgPointSmootherKind::chebyshev_jacobi;
#else
      constexpr bool chebyshev = false;
#endif
      if (!point_retained_defect_enabled(implementation) && !chebyshev) {
        status = finish_residual(implementation, level);
      }
    } else {
      status = compute_residual(implementation, level, stage);
    }
  }
  if (status) {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
    status = inject_cycle_failure(
        implementation, detail::MgCycleFailurePhase::restriction, level,
        deferred);
#endif
    if (status) status = restrict_residual(implementation, level, stage);
  }
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (status) ++implementation.matrix_work.cycle_restrictions;
#endif
  if (status) {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
    status = inject_cycle_failure(
        implementation, detail::MgCycleFailurePhase::first_coarse, level,
        deferred);
#endif
  }
  if (status) {
    status = mg_cycle(implementation, level + 1U, cycle_counter, stage,
                      deferred);
  }
  if (status && cycle_counter > 1U) {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
    status = inject_cycle_failure(
        implementation, detail::MgCycleFailurePhase::second_coarse, level,
        deferred);
#endif
  }
  if (status && cycle_counter > 1U) {
    status = mg_cycle(implementation, level + 1U,
                      static_cast<std::uint8_t>(cycle_counter - 1U), stage,
                      deferred);
  }
  if (!status) return status;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  status = inject_cycle_failure(
      implementation, detail::MgCycleFailurePhase::prolongation, level,
      deferred);
  if (!status) return status;
#endif
  status = prolongate_add(implementation, level, stage);
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  if (status) ++implementation.matrix_work.cycle_prolongations;
#endif
  if (status) status = implementation.services.workspace->revise_level(
      level, MgWorkspaceSlot::solution);
  if (status) {
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
    status = inject_cycle_failure(
        implementation, detail::MgCycleFailurePhase::post_smooth, level,
        deferred);
#endif
  }
  if (status) {
    status = smooth(implementation, level,
                    implementation.spec.policy.post_sweeps, true, false,
                    stage, deferred);
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
    if (spec.activity.cells.size != 0U) {
      candidate->activity_cells.assign(
          spec.activity.cells.data,
          spec.activity.cells.data + spec.activity.cells.size);
      candidate->activity_x_faces.assign(
          spec.activity.x_faces.data,
          spec.activity.x_faces.data + spec.activity.x_faces.size);
      candidate->activity_y_faces.assign(
          spec.activity.y_faces.data,
          spec.activity.y_faces.data + spec.activity.y_faces.size);
      candidate->activity_z_faces.assign(
          spec.activity.z_faces.data,
          spec.activity.z_faces.data + spec.activity.z_faces.size);
      candidate->spec.activity = {
          {candidate->activity_cells.data(), candidate->activity_cells.size()},
          {candidate->activity_x_faces.data(),
           candidate->activity_x_faces.size()},
          {candidate->activity_y_faces.data(),
           candidate->activity_y_faces.size()},
          {candidate->activity_z_faces.data(),
           candidate->activity_z_faces.size()},
          spec.activity.local_fingerprint,
          spec.activity.collective_fingerprint};
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
  std::uint64_t next_external_collectives =
      counters == nullptr ? 0U : counters->blocking_collectives;
  std::uint64_t next_external_bytes =
      counters == nullptr ? 0U : counters->collective_logical_bytes;
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
      local = initialize_replicated_coarse(
          *candidate, counters, next_external_collectives,
          next_external_bytes);
    }
    if (local) {
      local = build_coarse_coefficients(
          *candidate, candidate->hierarchy_storage.data());
    }
    if (local) {
      local = validate_certified_hierarchy(
          *candidate, candidate->hierarchy_storage.data());
    }
    if (local && candidate->replicated_coarse) {
      local = build_replicated_operator(
          *candidate, candidate->hierarchy_storage.data(),
          candidate->replicated_operator_active);
      if (local) {
        std::copy(candidate->replicated_operator_active.begin(),
                  candidate->replicated_operator_active.end(),
                  candidate->replicated_operator_inactive.begin());
      }
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
    counters->blocking_collectives = next_external_collectives;
    counters->collective_logical_bytes = next_external_bytes;
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
  const bool locally_unchanged =
      identity.revision == implementation.spec.coefficients.revision &&
      identity.fingerprint == implementation.spec.coefficients.fingerprint &&
      next_identity.fingerprint == implementation.spec.identity.fingerprint;
  std::uint64_t contract = UINT64_C(1469598103934665603);
  contract = mix(contract, next_identity.symbolic);
  contract = mix(contract, next_identity.hierarchy);
  // The workspace and composite identity are rank-local bindings: an odd or
  // otherwise non-divisible Cartesian decomposition legitimately gives
  // different local workspace shapes/fingerprints.  Their exact values are
  // validated below against each rank's compiled plan; the collective
  // contract can require presence, but must not require bitwise equality.
  contract = mix(contract, next_identity.workspace != 0U);
  contract = mix(contract, next_identity.fingerprint != 0U);
  contract = mix(contract, next_identity.numeric != 0U);
  contract = mix(contract, identity.revision != 0U);
  contract = mix(contract, identity.fingerprint != 0U);
  // Exact coefficient/storage identities are rank-local.  The collective
  // contract requires every rank to choose the same lifecycle branch while
  // local validation below retains exact identity matching.
  contract = mix(contract, locally_unchanged ? 1U : 0U);
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
  const bool unchanged = locally_unchanged;
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
  std::uint64_t next_external_collectives =
      counters == nullptr ? 0U : counters->blocking_collectives;
  std::uint64_t next_external_bytes =
      counters == nullptr ? 0U : counters->collective_logical_bytes;
  if (counters != nullptr && implementation.replicated_coarse) {
    const std::size_t logical_bytes =
        implementation.replicated_global_cells *
        kReplicatedOperatorWidth * sizeof(double);
    local = add_counter(counters->blocking_collectives, 1U,
                        next_external_collectives) &&
                    add_counter(counters->collective_logical_bytes,
                                static_cast<std::uint64_t>(logical_bytes),
                                next_external_bytes)
                ? Status{}
                : Status{StatusCode::invalid_plan, kMgCounter};
    agreed = consensus(implementation, local);
    if (!agreed) {
      return agreed;
    }
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
  local = validate_certified_hierarchy(
      implementation, implementation.inactive_hierarchy_storage.data());
  agreed = consensus(implementation, local);
  if (!agreed) {
    return agreed;
  }
  if (implementation.replicated_coarse) {
    local = build_replicated_operator(
        implementation, implementation.inactive_hierarchy_storage.data(),
        implementation.replicated_operator_inactive);
    agreed = consensus(implementation, local);
    if (!agreed) {
      return agreed;
    }
  }
  std::copy(implementation.inactive_hierarchy_storage.begin(),
            implementation.inactive_hierarchy_storage.end(),
            implementation.hierarchy_storage.begin());
  if (implementation.replicated_coarse) {
    std::copy(implementation.replicated_operator_inactive.begin(),
              implementation.replicated_operator_inactive.end(),
              implementation.replicated_operator_active.begin());
  }
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
    if (implementation.replicated_coarse) {
      counters->blocking_collectives = next_external_collectives;
      counters->collective_logical_bytes = next_external_bytes;
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
  return {implementation_->spec.identity,
          make_collective_preconditioner_fingerprint(
              implementation_->symbolic, implementation_->generation,
              implementation_->runtime_counters.numeric_refreshes),
          implementation_->spec.correction_scaling ==
                  MgCorrectionScaling::unit_linear
              ? LinearPreconditionerClass::fixed_general
              : LinearPreconditionerClass::flexible,
          LinearPreconditionerStatusScope::collective,
          LinearPreconditionerApplyLifecycle::prepared_batch};
}

Status NativeCartesianMgPlan::apply_impl(ConstFieldView residual,
                                          FieldView correction,
                                          bool prepared) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kMgApply};
  }
  Impl& implementation = *implementation_;
  Status local =
      prepared ? Status{} : validate_borrowed_services(implementation);
  Status agreed{};
  if (!prepared) {
    agreed = implementation.services.reductions->consensus(local);
    implementation.lowest =
        implementation.services.reductions->lowest_failing_rank();
    if (!agreed) {
      return agreed;
    }
  }
  if (!prepared) {
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
  if (!prepared) {
    local = candidate.base == nullptr || recursive.base == nullptr ||
                    rhs.base == nullptr || temporary.base == nullptr
                ? Status{StatusCode::invalid_plan, kMgApply}
                : Status{};
    agreed = consensus(implementation, local);
    if (!agreed) {
      return agreed;
    }
  }
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const std::size_t index = flat(cells, {i, j, k});
        const bool active = implementation.spec.activity.cells.size == 0U ||
                            implementation.spec.activity.cells.data[index] !=
                                0U;
        candidate.unchecked({i, j, k}, 0U) = 0.0;
        rhs.unchecked({i, j, k}, 0U) =
            active ? residual.unchecked({i, j, k}, 0U) : 0.0;
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
  if (!prepared) {
    agreed = consensus(implementation, local);
    if (!agreed) {
      return agreed;
    }
  } else if (!local) {
    // Every level/slot and the maximum application count were certified by
    // prepare_batch().  Preserve the fixed communication schedule if a
    // revision counter is nevertheless found out of contract.
    implementation.prepared_halo_deferred = local;
    local = {};
  }
  if (implementation.spec.null_space == MgNullSpace::constant) {
    local = project_constant(implementation, 0U, rhs);
    if (prepared && !local) {
      close_prepared_halo_epoch(implementation, false,
                                implementation.lowest);
      return local;
    }
    if (local) {
      local = implementation.services.workspace->revise_level(
          0U, MgWorkspaceSlot::rhs);
      rhs = implementation.services.workspace->level(
          0U, MgWorkspaceSlot::rhs);
    }
    if (!prepared) {
      agreed = consensus(implementation, local);
      if (!agreed) return agreed;
    } else if (!local) {
      implementation.prepared_halo_deferred = local;
      local = {};
    }
  }
  StageId stage = 700U;
  Status deferred{};
  implementation.replicated_solution_revision = 0U;
  const std::uint8_t cycle_counter =
      implementation.spec.policy.cycle == MgCycleKind::f_cycle ? 2U : 1U;
  local = mg_cycle(implementation, 0U, cycle_counter, stage, deferred);
  if (!prepared) {
    if (local && !deferred) {
      local = deferred;
    }
    agreed = consensus(implementation, local);
    if (!agreed) return agreed;
  } else if (!local) {
    close_prepared_halo_epoch(implementation, false,
                              implementation.lowest);
    return local;
  }
  candidate = implementation.services.workspace->level(
      0U, MgWorkspaceSlot::solution);
  if (implementation.spec.null_space == MgNullSpace::constant) {
    local = project_constant(implementation, 0U, candidate);
    if (prepared && !local) {
      close_prepared_halo_epoch(implementation, false,
                                implementation.lowest);
      return local;
    }
    if (local) {
      local = implementation.services.workspace->revise_level(
          0U, MgWorkspaceSlot::solution);
    }
    candidate = implementation.services.workspace->level(
        0U, MgWorkspaceSlot::solution);
    if (!prepared) {
      agreed = consensus(implementation, local);
      if (!agreed) {
        return agreed;
      }
    } else if (!local) {
      implementation.prepared_halo_deferred = local;
      local = {};
    }
  }
  local = compute_residual(implementation, 0U, stage);
  if (!prepared) {
    agreed = consensus(implementation, local);
    if (!agreed) return agreed;
  } else if (!local) {
    close_prepared_halo_epoch(implementation, false,
                              implementation.lowest);
    return local;
  }
  recursive = implementation.services.workspace->level(
      0U, MgWorkspaceSlot::residual);
  double local_projection[3]{};
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 cell{i, j, k};
        const std::size_t index = flat(cells, cell);
        const bool active = implementation.spec.activity.cells.size == 0U ||
                            implementation.spec.activity.cells.data[index] !=
                                0U;
        const double r = active ? residual.unchecked(cell, 0U) : 0.0;
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
  if (!local) {
    if (prepared) {
      close_prepared_halo_epoch(implementation, false,
                                implementation.lowest);
    }
    return local;
  }
  const bool finite_projection = std::isfinite(projection[0]) &&
                                 std::isfinite(projection[1]) &&
                                 std::isfinite(projection[2]) &&
                                 projection[0] >= 0.0 &&
                                 projection[2] >= 0.0;
  if (!finite_projection) {
    if (!prepared) {
      return {StatusCode::numerical_failure, kMgApply};
    }
    if (deferred) {
      deferred = {StatusCode::numerical_failure, kMgApply};
    }
    projection[0] = 0.0;
    projection[1] = 0.0;
    projection[2] = 0.0;
  }
  const double scale = std::max(1.0, projection[0]);
  const double tiny = std::numeric_limits<double>::epsilon() * scale;
#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
  const bool skip_final_projection = implementation.spec.correction_scaling ==
                                         MgCorrectionScaling::unit_linear ||
                                     implementation.skip_final_projection;
#else
  const bool skip_final_projection = implementation.spec.correction_scaling ==
                                     MgCorrectionScaling::unit_linear;
#endif
  const double alpha = skip_final_projection
                           ? 1.0
                           : (projection[1] > 0.0 && projection[2] > tiny
                                  ? projection[1] / projection[2]
                                  : 0.0);
  const bool unit_linear = implementation.spec.correction_scaling ==
                           MgCorrectionScaling::unit_linear;
  const double final_squared =
      unit_linear
          ? std::max(0.0, projection[0] - 2.0 * projection[1] + projection[2])
          : (skip_final_projection
                 ? projection[2]
                 : (alpha == 0.0
                        ? projection[0]
                        : std::max(0.0, projection[0] - projection[1] *
                                                            projection[1] /
                                                            projection[2])));
  implementation.last_initial = std::sqrt(projection[0]);
  implementation.last_final = std::sqrt(final_squared);
  if (!std::isfinite(alpha) || !std::isfinite(implementation.last_initial) ||
      !std::isfinite(implementation.last_final)) {
    if (!prepared) {
      return {StatusCode::numerical_failure, kMgApply};
    }
    if (deferred) {
      deferred = {StatusCode::numerical_failure, kMgApply};
    }
  }
  for (std::int32_t k = 0; k < cells.z; ++k) {
    for (std::int32_t j = 0; j < cells.y; ++j) {
      for (std::int32_t i = 0; i < cells.x; ++i) {
        const Int3 cell{i, j, k};
        const std::size_t index = flat(cells, cell);
        const bool active = implementation.spec.activity.cells.size == 0U ||
                            implementation.spec.activity.cells.data[index] !=
                                0U;
        candidate.unchecked(cell, 0U) =
            active ? candidate.unchecked(cell, 0U) * alpha : 0.0;
        const double r = active ? residual.unchecked(cell, 0U) : 0.0;
        recursive.unchecked(cell, 0U) =
            r - alpha * (r - recursive.unchecked(cell, 0U));
      }
    }
  }
  std::uint64_t next_applications = 0U;
  local = increment(implementation.runtime_counters.applications,
                    next_applications)
              ? Status{}
              : Status{StatusCode::invalid_plan, kMgCounter};
  if (!prepared) {
    agreed = consensus(implementation, local);
    if (!agreed) {
      return agreed;
    }
  } else {
    if (local && !deferred) {
      local = deferred;
    }
    if (local && !implementation.prepared_halo_deferred) {
      local = implementation.prepared_halo_deferred;
    }
    agreed = consensus(implementation, local);
    close_prepared_halo_epoch(implementation, static_cast<bool>(agreed),
                              implementation.lowest);
    if (!agreed) {
      return agreed;
    }
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

Status NativeCartesianMgPlan::apply(ConstFieldView residual,
                                    FieldView correction,
                                    std::uint32_t iteration) noexcept {
  (void)iteration;
  return apply_impl(residual, correction, false);
}

Status NativeCartesianMgPlan::prepare_batch(
    const LinearPreconditionerBatchDescriptor& descriptor,
    LinearPreconditionerBatchTicket& ticket) noexcept {
  ticket = {};
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kMgApply};
  }
  Impl& implementation = *implementation_;
  const Status valid = validate_prepared_batch(implementation, descriptor);
  if (!valid) {
    return valid;
  }
  const PlanFingerprint preconditioner_fingerprint =
      make_collective_preconditioner_fingerprint(
          implementation.symbolic, implementation.generation,
          implementation.runtime_counters.numeric_refreshes);
  issue_batch_ticket(ticket, this, descriptor, preconditioner_fingerprint,
                     implementation.generation);
  ticket.preconditioner_application_base =
      implementation.runtime_counters.applications;
  return {};
}

Status NativeCartesianMgPlan::apply_prepared(
    ConstFieldView residual, FieldView correction, std::uint32_t iteration,
    const LinearPreconditionerBatchTicket& ticket) noexcept {
  if (implementation_ == nullptr) {
    return {StatusCode::invalid_plan, kMgApply};
  }
  Impl& implementation = *implementation_;
  const std::uintptr_t owner = reinterpret_cast<std::uintptr_t>(
      static_cast<const LinearPreconditioner*>(this));
  const bool valid_ticket =
      ticket.owner == owner && ticket.workspace != nullptr &&
      ticket.workspace->fingerprint() != 0U && ticket.slot_count != 0U &&
      ticket.maximum_applications != 0U &&
      ticket.preconditioner_fingerprint ==
          make_collective_preconditioner_fingerprint(
              implementation.symbolic, implementation.generation,
              implementation.runtime_counters.numeric_refreshes) &&
      ticket.preconditioner_generation == implementation.generation &&
      ticket.preconditioner_application_base <=
          implementation.runtime_counters.applications &&
      ticket.preconditioner_application_base <=
          std::numeric_limits<std::uint64_t>::max() -
              ticket.maximum_applications &&
      implementation.runtime_counters.applications <
          ticket.preconditioner_application_base +
              ticket.maximum_applications &&
      iteration < ticket.maximum_applications;
  // prepare_batch() already cold-certified every borrowed service identity;
  // the ticket binds that proof to the exact plan generation and application
  // budget.  Do not rescan the hierarchy on every Krylov application.
  Status local = valid_ticket
                     ? Status{}
                     : Status{StatusCode::invalid_plan, kMgApply};
  if (local) {
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
  }
  const Status epoch = enter_prepared_halo_epoch(implementation);
  if (local && !epoch) local = epoch;
  implementation.lowest = -1;
  const Status agreed = consensus(implementation, local);
  if (!agreed) {
    close_prepared_halo_epoch(implementation, false,
                              implementation.lowest);
    return agreed;
  }
  return apply_impl(residual, correction, true);
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

  static void reference_point_actions(NativeCartesianMgPlan& plan,
                                      bool enabled) noexcept {
    if (plan.implementation_ != nullptr) {
      plan.implementation_->reference_point_actions = enabled;
      plan.implementation_->matrix_work = {};
    }
  }

  static void reference_chebyshev_lifecycle(NativeCartesianMgPlan& plan,
                                            bool enabled) noexcept {
    if (plan.implementation_ != nullptr) {
      plan.implementation_->reference_chebyshev_lifecycle = enabled;
      plan.implementation_->matrix_work = {};
    }
  }

  static void pre_sweeps(NativeCartesianMgPlan& plan,
                         std::uint8_t sweeps) noexcept {
    if (plan.implementation_ != nullptr) {
      plan.implementation_->spec.policy.pre_sweeps = sweeps;
    }
  }

  static void skip_final_projection(NativeCartesianMgPlan& plan,
                                    bool enabled) noexcept {
    if (plan.implementation_ != nullptr) {
      plan.implementation_->skip_final_projection = enabled;
    }
  }

  static void force_chebyshev_invalid(NativeCartesianMgPlan& plan,
                                      bool enabled) noexcept {
    if (plan.implementation_ != nullptr) {
      plan.implementation_->force_chebyshev_invalid = enabled;
    }
  }

  static void cycle_failure(NativeCartesianMgPlan& plan,
                            MgCycleFailurePhase phase, std::size_t level,
                            int failing_rank, bool deferred) noexcept {
    if (plan.implementation_ != nullptr) {
      plan.implementation_->cycle_failure_phase = phase;
      plan.implementation_->cycle_failure_level = level;
      plan.implementation_->cycle_failure_rank = failing_rank;
      plan.implementation_->cycle_failure_deferred = deferred;
      plan.implementation_->cycle_failure_consumed = false;
    }
  }

  static void replicated_coarse_mode(NativeCartesianMgPlan& plan,
                                     std::uint8_t mode) noexcept {
    if (plan.implementation_ != nullptr) {
      plan.implementation_->reference_replicated_coarse = mode == 1U;
      plan.implementation_->request_replicated_coarse = mode == 2U;
    }
  }

  static bool replicated_coarse_enabled(
      const NativeCartesianMgPlan& plan) noexcept {
    return plan.implementation_ != nullptr &&
           replicated_coarse_route_enabled(*plan.implementation_);
  }

  static MgMatrixWorkCounters matrix_work(
      const NativeCartesianMgPlan& plan) noexcept {
    return plan.implementation_ == nullptr
               ? MgMatrixWorkCounters{}
               : plan.implementation_->matrix_work;
  }

  static double level_diagonal(const NativeCartesianMgPlan& plan,
                               std::size_t level,
                               std::size_t cell) noexcept {
    if (plan.implementation_ == nullptr ||
        level >= plan.implementation_->levels.size()) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    const MgLevelStorage& storage = plan.implementation_->levels[level];
    if (cell >= cell_count(storage.view.local_shape)) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return block(plan.implementation_->hierarchy_storage.data(),
                 storage.diagonal_offset)[cell];
  }

};

void set_mg_runtime_counters_for_test(
    NativeCartesianMgPlan& plan, MgPlanCounters counters,
    RevisionToken generation) noexcept {
  NativeMgTestAccess::set(plan, counters, generation);
}

void set_mg_reference_point_actions_for_test(
    NativeCartesianMgPlan& plan, bool enabled) noexcept {
  NativeMgTestAccess::reference_point_actions(plan, enabled);
}

void set_mg_reference_chebyshev_lifecycle_for_test(
    NativeCartesianMgPlan& plan, bool enabled) noexcept {
  NativeMgTestAccess::reference_chebyshev_lifecycle(plan, enabled);
}

void set_mg_pre_sweeps_for_test(NativeCartesianMgPlan& plan,
                                std::uint8_t sweeps) noexcept {
  NativeMgTestAccess::pre_sweeps(plan, sweeps);
}

void set_mg_skip_final_projection_for_test(NativeCartesianMgPlan& plan,
                                           bool enabled) noexcept {
  NativeMgTestAccess::skip_final_projection(plan, enabled);
}

void set_mg_force_chebyshev_invalid_for_test(
    NativeCartesianMgPlan& plan, bool enabled) noexcept {
  NativeMgTestAccess::force_chebyshev_invalid(plan, enabled);
}

void set_mg_cycle_failure_for_test(
    NativeCartesianMgPlan& plan, MgCycleFailurePhase phase,
    std::size_t level, int failing_rank, bool deferred) noexcept {
  NativeMgTestAccess::cycle_failure(plan, phase, level, failing_rank,
                                    deferred);
}

void set_mg_replicated_coarse_mode_for_test(
    NativeCartesianMgPlan& plan, std::uint8_t mode) noexcept {
  NativeMgTestAccess::replicated_coarse_mode(plan, mode);
}

bool mg_replicated_coarse_enabled_for_test(
    const NativeCartesianMgPlan& plan) noexcept {
  return NativeMgTestAccess::replicated_coarse_enabled(plan);
}

MgMatrixWorkCounters mg_matrix_work_counters_for_test(
    const NativeCartesianMgPlan& plan) noexcept {
  return NativeMgTestAccess::matrix_work(plan);
}

double mg_level_diagonal_for_test(const NativeCartesianMgPlan& plan,
                                  std::size_t level,
                                  std::size_t cell) noexcept {
  return NativeMgTestAccess::level_diagonal(plan, level, cell);
}

}  // namespace detail
#endif

}  // namespace hundun::v04

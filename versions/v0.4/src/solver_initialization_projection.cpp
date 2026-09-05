// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_initialization.hpp"

#include "field_view_interval_detail.hpp"
#include "solver_cartesian_detail.hpp"

#include "hundun/v04_parallel.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
#include <unordered_map>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace detail {

enum class FreshProjectionFailurePoint : std::uint8_t {
  none,
  exact_operator_boundary,
  prepare_mg_boundary,
  audit_boundary
};

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)
namespace {
std::atomic<FreshProjectionFailurePoint> g_fresh_projection_failure_point{
    FreshProjectionFailurePoint::none};
std::atomic<int> g_fresh_projection_failure_rank{-1};
} // namespace

void set_fresh_projection_failure_for_test(FreshProjectionFailurePoint point,
                                           int failing_rank) noexcept {
  g_fresh_projection_failure_rank.store(failing_rank,
                                        std::memory_order_relaxed);
  g_fresh_projection_failure_point.store(point, std::memory_order_relaxed);
}

void clear_fresh_projection_failure_for_test() noexcept {
  g_fresh_projection_failure_point.store(FreshProjectionFailurePoint::none,
                                         std::memory_order_relaxed);
  g_fresh_projection_failure_rank.store(-1, std::memory_order_relaxed);
}

bool fresh_projection_failure_injected_for_test(
    FreshProjectionFailurePoint point, int rank) noexcept {
  if (g_fresh_projection_failure_rank.load(std::memory_order_relaxed) != rank)
    return false;
  FreshProjectionFailurePoint expected = point;
  return g_fresh_projection_failure_point.compare_exchange_strong(
      expected, FreshProjectionFailurePoint::none, std::memory_order_relaxed);
}
#else
bool fresh_projection_failure_injected_for_test(FreshProjectionFailurePoint,
                                                int) noexcept {
  return false;
}
#endif

} // namespace detail
namespace {

constexpr std::uint32_t kFreshProjectionPlan = 9721U;
constexpr std::uint32_t kFreshProjectionPrepare = 9722U;
constexpr std::uint32_t kFreshProjectionCompatibility = 9723U;
constexpr std::uint32_t kFreshProjectionSolve = 9724U;
constexpr std::uint32_t kFreshProjectionAudit = 9725U;
constexpr std::uint32_t kFreshProjectionCommit = 9726U;
constexpr std::uint32_t kFreshProjectionNumerical = 9727U;
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);
constexpr std::uint32_t kNoComponent = UINT32_MAX;

std::uint64_t mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kFnvPrime;
  return hash;
}

std::uint64_t finish(std::uint64_t value) noexcept {
  return value == 0U ? 1U : value;
}

std::uint64_t double_bits(double value) noexcept {
  std::uint64_t bits = 0U;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

bool same_cells(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

std::size_t cell_count(Int3 cells) noexcept {
  return static_cast<std::size_t>(cells.x) * static_cast<std::size_t>(cells.y) *
         static_cast<std::size_t>(cells.z);
}

std::size_t flat_cell(Int3 cells, Int3 cell) noexcept {
  return static_cast<std::size_t>(cell.x) +
         static_cast<std::size_t>(cells.x) *
             (static_cast<std::size_t>(cell.y) +
              static_cast<std::size_t>(cells.y) *
                  static_cast<std::size_t>(cell.z));
}

Int3 face_extents(Int3 cells, CartesianAxis axis) noexcept {
  if (axis == CartesianAxis::x)
    ++cells.x;
  if (axis == CartesianAxis::y)
    ++cells.y;
  if (axis == CartesianAxis::z)
    ++cells.z;
  return cells;
}

std::size_t flat_face(Int3 extents, Int3 face) noexcept {
  return static_cast<std::size_t>(face.x) +
         static_cast<std::size_t>(extents.x) *
             (static_cast<std::size_t>(face.y) +
              static_cast<std::size_t>(extents.y) *
                  static_cast<std::size_t>(face.z));
}

std::int32_t axis_value(Int3 value, CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? value.x
             : (axis == CartesianAxis::y ? value.y : value.z);
}

void set_axis_value(Int3 &value, CartesianAxis axis,
                    std::int32_t selected) noexcept {
  if (axis == CartesianAxis::x)
    value.x = selected;
  if (axis == CartesianAxis::y)
    value.y = selected;
  if (axis == CartesianAxis::z)
    value.z = selected;
}

std::uint64_t global_gid(Int3 global, Int3 cell) noexcept {
  return static_cast<std::uint64_t>(cell.x) +
         static_cast<std::uint64_t>(global.x) *
             (static_cast<std::uint64_t>(cell.y) +
              static_cast<std::uint64_t>(global.y) *
                  static_cast<std::uint64_t>(cell.z));
}

bool valid_face_view(ConstFaceFieldView face, CartesianAxis axis,
                     Int3 cells) noexcept {
  const Int3 expected = face_extents(cells, axis);
  return face.base != nullptr && face.axis == axis &&
         same_cells(face.extents, expected) &&
         face.stride_y >= static_cast<std::size_t>(expected.x) &&
         face.stride_z >=
             face.stride_y * static_cast<std::size_t>(expected.y) &&
         face.storage_identity != 0U && face.revision_domain != 0U;
}

bool valid_face_view(FaceFieldView face, CartesianAxis axis,
                     Int3 cells) noexcept {
  return valid_face_view(as_const(face), axis, cells);
}

ConstFaceFieldView select(ConstFaceFluxView flux, CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? flux.x
             : (axis == CartesianAxis::y ? flux.y : flux.z);
}

FaceFieldView select(FaceFluxView flux, CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? flux.x
             : (axis == CartesianAxis::y ? flux.y : flux.z);
}

ConstFaceFieldView select(ConstFaceFieldView x, ConstFaceFieldView y,
                          ConstFaceFieldView z, CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x ? x : (axis == CartesianAxis::y ? y : z);
}

FaceFieldView select(FaceFieldView x, FaceFieldView y, FaceFieldView z,
                     CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x ? x : (axis == CartesianAxis::y ? y : z);
}

bool empty_activity(MgDomainActivityView activity) noexcept {
  return activity.cells.data == nullptr && activity.cells.size == 0U &&
         activity.x_faces.data == nullptr && activity.x_faces.size == 0U &&
         activity.y_faces.data == nullptr && activity.y_faces.size == 0U &&
         activity.z_faces.data == nullptr && activity.z_faces.size == 0U &&
         activity.local_fingerprint == 0U &&
         activity.collective_fingerprint == 0U;
}

bool binary_span(Span<const std::uint8_t> values,
                 std::size_t expected) noexcept {
  if (values.data == nullptr || values.size != expected)
    return false;
  for (std::size_t index = 0U; index < values.size; ++index)
    if (values.data[index] > 1U)
      return false;
  return true;
}

bool valid_activity(MgDomainActivityView activity, Int3 cells) noexcept {
  if (empty_activity(activity))
    return true;
  return binary_span(activity.cells, cell_count(cells)) &&
         binary_span(activity.x_faces,
                     cell_count(face_extents(cells, CartesianAxis::x))) &&
         binary_span(activity.y_faces,
                     cell_count(face_extents(cells, CartesianAxis::y))) &&
         binary_span(activity.z_faces,
                     cell_count(face_extents(cells, CartesianAxis::z))) &&
         activity.local_fingerprint != 0U &&
         activity.collective_fingerprint != 0U;
}

bool active_cell(MgDomainActivityView activity, Int3 cells,
                 Int3 cell) noexcept {
  return empty_activity(activity) ||
         activity.cells.data[flat_cell(cells, cell)] != 0U;
}

bool active_face(MgDomainActivityView activity, Int3 cells, CartesianAxis axis,
                 Int3 face) noexcept {
  if (empty_activity(activity))
    return true;
  const Int3 extents = face_extents(cells, axis);
  const Span<const std::uint8_t> selected =
      axis == CartesianAxis::x
          ? activity.x_faces
          : (axis == CartesianAxis::y ? activity.y_faces : activity.z_faces);
  return selected.data[flat_face(extents, face)] != 0U;
}

std::uint64_t hash_field_values(std::uint64_t hash, ConstFieldView field,
                                Int3 cells, std::uint8_t components) noexcept {
  hash = mix(hash, reinterpret_cast<std::uintptr_t>(field.base));
  hash = mix(hash, static_cast<std::uint32_t>(field.interior.x));
  hash = mix(hash, static_cast<std::uint32_t>(field.interior.y));
  hash = mix(hash, static_cast<std::uint32_t>(field.interior.z));
  hash = mix(hash, static_cast<std::uint32_t>(field.ghosts.x));
  hash = mix(hash, static_cast<std::uint32_t>(field.ghosts.y));
  hash = mix(hash, static_cast<std::uint32_t>(field.ghosts.z));
  hash = mix(hash, field.components);
  hash = mix(hash, field.stride_y);
  hash = mix(hash, field.stride_z);
  hash = mix(hash, field.component_stride);
  hash = mix(hash, field.replica);
  hash = mix(hash, field.field);
  hash = mix(hash, field.revision);
  hash = mix(hash, field.storage_identity);
  hash = mix(hash, field.revision_domain);
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        for (std::uint8_t component = 0U; component < components; ++component)
          hash = mix(hash, double_bits(field.unchecked({x, y, z}, component)));
  return hash;
}

std::uint64_t hash_field_numeric_values(std::uint64_t hash,
                                        ConstFieldView field, Int3 cells,
                                        std::uint8_t components) noexcept {
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        for (std::uint8_t component = 0U; component < components; ++component)
          hash = mix(hash, double_bits(field.unchecked({x, y, z}, component)));
  return hash;
}

std::uint64_t hash_flux_values(std::uint64_t hash,
                               ConstFaceFluxView flux) noexcept {
  hash = mix(hash, flux.revision);
  for (CartesianAxis axis :
       {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
    const ConstFaceFieldView face = select(flux, axis);
    hash = mix(hash, reinterpret_cast<std::uintptr_t>(face.base));
    hash = mix(hash, static_cast<std::uint32_t>(face.extents.x));
    hash = mix(hash, static_cast<std::uint32_t>(face.extents.y));
    hash = mix(hash, static_cast<std::uint32_t>(face.extents.z));
    hash = mix(hash, face.stride_y);
    hash = mix(hash, face.stride_z);
    hash = mix(hash, static_cast<std::uint8_t>(face.axis));
    hash = mix(hash, face.storage_identity);
    hash = mix(hash, face.revision_domain);
    for (std::int32_t z = 0; z < face.extents.z; ++z)
      for (std::int32_t y = 0; y < face.extents.y; ++y)
        for (std::int32_t x = 0; x < face.extents.x; ++x)
          hash = mix(hash, double_bits(face.unchecked({x, y, z})));
  }
  return hash;
}

bool same_field_binding(ConstFieldView left, ConstFieldView right) noexcept {
  return left.base == right.base && same_cells(left.interior, right.interior) &&
         same_cells(left.ghosts, right.ghosts) &&
         left.components == right.components &&
         left.stride_y == right.stride_y && left.stride_z == right.stride_z &&
         left.component_stride == right.component_stride &&
         left.replica == right.replica && left.field == right.field &&
         left.revision == right.revision &&
         left.storage_identity == right.storage_identity &&
         left.revision_domain == right.revision_domain;
}

bool same_flux_binding(ConstFaceFluxView left,
                       ConstFaceFluxView right) noexcept {
  if (left.revision != right.revision)
    return false;
  for (CartesianAxis axis :
       {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
    const ConstFaceFieldView a = select(left, axis);
    const ConstFaceFieldView b = select(right, axis);
    if (a.base != b.base || !same_cells(a.extents, b.extents) ||
        a.stride_y != b.stride_y || a.stride_z != b.stride_z ||
        a.axis != b.axis || a.storage_identity != b.storage_identity ||
        a.revision_domain != b.revision_domain)
      return false;
  }
  return true;
}

template <std::size_t Count>
bool disjoint_fields(const std::array<ConstFieldView, Count> &fields) noexcept {
  for (std::size_t left = 0U; left < Count; ++left)
    for (std::size_t right = left + 1U; right < Count; ++right)
      if (detail::field_views_overlap(fields[left], fields[right]))
        return false;
  return true;
}

template <std::size_t Count>
bool disjoint_faces(
    const std::array<ConstFaceFieldView, Count> &faces) noexcept {
  for (std::size_t left = 0U; left < Count; ++left)
    for (std::size_t right = left + 1U; right < Count; ++right)
      if (detail::face_views_overlap(faces[left], faces[right]))
        return false;
  return true;
}

template <std::size_t FieldCount, std::size_t FaceCount>
bool disjoint_cells_and_faces(
    const std::array<ConstFieldView, FieldCount> &fields,
    const std::array<ConstFaceFieldView, FaceCount> &faces) noexcept {
  for (ConstFieldView field : fields)
    for (ConstFaceFieldView face : faces)
      if (detail::cell_face_views_overlap(field, face))
        return false;
  return true;
}

std::array<ConstFieldView, 4U> workspace_cell_views(
    const FreshStartKinematicProjectionWorkspace &workspace) noexcept {
  return {as_const(workspace.chi), as_const(workspace.rhs),
          as_const(workspace.diagonal), as_const(workspace.candidate_velocity)};
}

std::array<ConstFaceFieldView, 9U> workspace_face_views(
    const FreshStartKinematicProjectionWorkspace &workspace) noexcept {
  return {as_const(workspace.x_physical_mobility),
          as_const(workspace.y_physical_mobility),
          as_const(workspace.z_physical_mobility),
          as_const(workspace.x_solver_coefficient),
          as_const(workspace.y_solver_coefficient),
          as_const(workspace.z_solver_coefficient),
          as_const(workspace.candidate_mass_flux.x),
          as_const(workspace.candidate_mass_flux.y),
          as_const(workspace.candidate_mass_flux.z)};
}

bool workspace_alias_matrix_valid(
    const FreshStartKinematicProjectionWorkspace &workspace) noexcept {
  const auto fields = workspace_cell_views(workspace);
  const auto faces = workspace_face_views(workspace);
  return disjoint_fields(fields) && disjoint_faces(faces) &&
         disjoint_cells_and_faces(fields, faces);
}

bool face_overlaps_raw_storage(ConstFaceFieldView face, std::uintptr_t address,
                               std::size_t doubles) noexcept {
  detail::FieldStorageInterval face_interval{};
  if (!detail::face_storage_interval(face, face_interval) || address == 0U ||
      doubles == 0U ||
      doubles > std::numeric_limits<std::uintptr_t>::max() / sizeof(double))
    return true;
  const std::uintptr_t bytes =
      static_cast<std::uintptr_t>(doubles) * sizeof(double);
  if (bytes > std::numeric_limits<std::uintptr_t>::max() - address)
    return true;
  return detail::storage_intervals_overlap(face_interval,
                                           {address, address + bytes});
}

bool face_overlaps_solver_workspace(ConstFaceFieldView face,
                                    const SolverWorkspace &workspace) noexcept {
  const LinearWorkspaceRequirements &requirements = workspace.requirements();
  for (std::uint8_t slot = 0U; slot < requirements.vector_slots; ++slot)
    if (detail::cell_face_views_overlap(
            as_const(workspace.vector(slot, requirements.maximum_shape)), face))
      return true;
  const Span<double> scalars =
      workspace.scalars(0U, requirements.scalar_doubles);
  return face_overlaps_raw_storage(
      face, reinterpret_cast<std::uintptr_t>(scalars.data), scalars.size);
}

bool face_overlaps_mg_workspace(ConstFaceFieldView face,
                                const MgWorkspace &workspace) noexcept {
  return face_overlaps_raw_storage(face, workspace.storage_address(),
                                   workspace.storage_doubles());
}

bool solver_and_mg_workspaces_disjoint(const SolverWorkspace &solver,
                                       const MgWorkspace &mg) noexcept {
  const LinearWorkspaceRequirements &requirements = solver.requirements();
  for (std::uint8_t slot = 0U; slot < requirements.vector_slots; ++slot)
    if (mg.overlaps_storage(
            as_const(solver.vector(slot, requirements.maximum_shape))))
      return false;
  const Span<double> scalars = solver.scalars(0U, requirements.scalar_doubles);
  const std::uintptr_t scalar_address =
      reinterpret_cast<std::uintptr_t>(scalars.data);
  const std::uintptr_t mg_address = mg.storage_address();
  if (scalar_address == 0U || scalars.size == 0U || mg_address == 0U ||
      mg.storage_doubles() == 0U ||
      scalars.size >
          std::numeric_limits<std::uintptr_t>::max() / sizeof(double) ||
      mg.storage_doubles() >
          std::numeric_limits<std::uintptr_t>::max() / sizeof(double))
    return false;
  const std::uintptr_t scalar_bytes =
      static_cast<std::uintptr_t>(scalars.size) * sizeof(double);
  const std::uintptr_t mg_bytes =
      static_cast<std::uintptr_t>(mg.storage_doubles()) * sizeof(double);
  if (scalar_bytes >
          std::numeric_limits<std::uintptr_t>::max() - scalar_address ||
      mg_bytes > std::numeric_limits<std::uintptr_t>::max() - mg_address)
    return false;
  return !detail::storage_intervals_overlap(
      {scalar_address, scalar_address + scalar_bytes},
      {mg_address, mg_address + mg_bytes});
}

bool input_alias_matrix_valid(
    const FreshStartKinematicProjectionInput &input,
    const FreshStartKinematicProjectionWorkspace &workspace,
    const SolverWorkspace *solver_workspace,
    const MgWorkspace *mg_workspace) noexcept {
  std::array<ConstFieldView, 4U> inputs{input.density, input.velocity,
                                        input.velocity_accepted_n_minus_one,
                                        input.velocity_trial};
  constexpr std::size_t input_count = 4U;
  const std::array<ConstFaceFieldView, 3U> input_faces{
      input.mass_flux.x, input.mass_flux.y, input.mass_flux.z};
  const auto scratch_fields = workspace_cell_views(workspace);
  const auto scratch_faces = workspace_face_views(workspace);
  for (std::size_t left = 0U; left < input_count; ++left) {
    for (std::size_t right = left + 1U; right < input_count; ++right)
      if (detail::field_views_overlap(inputs[left], inputs[right]))
        return false;
    for (ConstFieldView scratch : scratch_fields)
      if (detail::field_views_overlap(inputs[left], scratch))
        return false;
    for (ConstFaceFieldView face : input_faces)
      if (detail::cell_face_views_overlap(inputs[left], face))
        return false;
    for (ConstFaceFieldView scratch : scratch_faces)
      if (detail::cell_face_views_overlap(inputs[left], scratch))
        return false;
    if ((solver_workspace != nullptr &&
         solver_workspace->overlaps_storage(inputs[left])) ||
        (mg_workspace != nullptr &&
         mg_workspace->overlaps_storage(inputs[left])))
      return false;
  }
  if (!disjoint_faces(input_faces))
    return false;
  for (ConstFaceFieldView input_face : input_faces) {
    for (ConstFieldView scratch : scratch_fields)
      if (detail::cell_face_views_overlap(scratch, input_face))
        return false;
    for (ConstFaceFieldView scratch : scratch_faces)
      if (detail::face_views_overlap(input_face, scratch))
        return false;
    if ((solver_workspace != nullptr &&
         face_overlaps_solver_workspace(input_face, *solver_workspace)) ||
        (mg_workspace != nullptr &&
         face_overlaps_mg_workspace(input_face, *mg_workspace)))
      return false;
  }
  return true;
}

std::uint64_t
hash_projection_input(std::uint64_t hash,
                      const FreshStartKinematicProjectionInput &input,
                      Int3 cells) noexcept {
  hash = mix(hash, input.state);
  hash = hash_field_values(hash, input.density, cells, 1U);
  hash = hash_field_values(hash, input.velocity, cells, 3U);
  hash =
      hash_field_values(hash, input.velocity_accepted_n_minus_one, cells, 3U);
  hash = hash_field_values(hash, input.velocity_trial, cells, 3U);
  hash = hash_flux_values(hash, input.mass_flux);
  return hash;
}

class DisjointSet {
public:
  explicit DisjointSet(std::size_t size) : parent_(size), rank_(size, 0U) {
    for (std::size_t index = 0U; index < size; ++index)
      parent_[index] = index;
  }
  std::size_t find(std::size_t value) noexcept {
    while (parent_[value] != value) {
      parent_[value] = parent_[parent_[value]];
      value = parent_[value];
    }
    return value;
  }
  void unite(std::size_t left, std::size_t right) noexcept {
    left = find(left);
    right = find(right);
    if (left == right)
      return;
    if (rank_[left] < rank_[right])
      std::swap(left, right);
    parent_[right] = left;
    if (rank_[left] == rank_[right])
      ++rank_[left];
  }

private:
  std::vector<std::size_t> parent_;
  std::vector<std::uint8_t> rank_;
};

} // namespace

struct FreshStartKinematicProjectionPlan::Impl {
  struct ExactOperator final : LinearOperator {
    const Impl *owner{};
    LinearOperatorCertificate certificate() const noexcept override {
      return owner == nullptr ? LinearOperatorCertificate{}
                              : owner->linear_certificate;
    }
    Status apply(FieldView x, FieldView y) const noexcept override;
    LinearOperatorFailureProvenance
    failure_provenance() const noexcept override {
      return owner == nullptr ? LinearOperatorFailureProvenance{}
                              : owner->operator_failure;
    }
  } exact_operator;

  struct JacobiPreconditioner final : LinearPreconditioner {
    const Impl *owner{};
    LinearPreconditionerCertificate certificate() const noexcept override {
      if (owner == nullptr || owner->linear_identity.fingerprint == 0U)
        return {};
      return {owner->linear_identity, owner->collective_operator,
              LinearPreconditionerClass::fixed_spd,
              LinearPreconditionerStatusScope::rank_local,
              LinearPreconditionerApplyLifecycle::per_call_checked};
    }
    Status apply(ConstFieldView input, FieldView output,
                 std::uint32_t) noexcept override;
  } jacobi;

  struct ContinuityAudit final : LinearConvergenceAudit {
    Impl *owner{};
    LinearConvergenceAuditCertificate certificate() const noexcept override {
      return {owner == nullptr ? 0U : owner->collective_operator};
    }
    Status evaluate(ConstFieldView solution, ConstFieldView,
                    ReductionEngine &reductions,
                    LinearConvergenceAuditResult &result) noexcept override;
  } continuity_audit;

  MPI_Comm communicator{MPI_COMM_NULL};
  int rank{-1};
  int size{};
  const CartesianGeometryPlan *geometry{};
  const CartesianKernelPlan *kernels{};
  const BoundaryPlan *boundary{};
  MeshPatch patch{};
  MgDomainActivityView source_activity{};
  FreshStartProjectionLinearRoute route{
      FreshStartProjectionLinearRoute::native_mg_fgmres};
  LinearSolveControl solve_control{};
  MgHierarchyPolicy mg_policy{};
  double compatibility_absolute{};
  double compatibility_relative{};
  double continuity_absolute{};
  double continuity_relative{};
  FreshStartKinematicProjectionServices services{};
  FreshStartKinematicProjectionWorkspace workspace{};
  PressureCorrectionBoundaryPlan pressure_boundary{};
  std::vector<std::uint8_t> source_cells;
  std::vector<std::uint8_t> source_x_faces;
  std::vector<std::uint8_t> source_y_faces;
  std::vector<std::uint8_t> source_z_faces;
  std::vector<std::uint8_t> mg_cells;
  std::vector<std::uint8_t> mg_x_faces;
  std::vector<std::uint8_t> mg_y_faces;
  std::vector<std::uint8_t> mg_z_faces;
  MgDomainActivityView mg_activity{};
  std::vector<std::uint32_t> local_components;
  std::vector<std::uint64_t> component_labels;
  std::vector<std::uint64_t> component_anchors;
  std::vector<std::uint8_t> component_dirichlet;
  std::vector<double> component_local_sums;
  std::vector<double> component_local_scales;
  std::vector<std::uint8_t> component_distributed;
  std::vector<std::uint32_t> component_route_local_index;
  std::vector<std::uint32_t> component_route_send_slot;
  std::vector<int> component_send_counts;
  std::vector<int> component_send_displacements;
  std::vector<int> component_receive_counts;
  std::vector<int> component_receive_displacements;
  std::vector<int> component_send_value_counts;
  std::vector<int> component_send_value_displacements;
  std::vector<int> component_receive_value_counts;
  std::vector<int> component_receive_value_displacements;
  std::vector<std::uint64_t> component_send_labels;
  std::vector<std::uint64_t> component_receive_labels;
  std::vector<std::uint32_t> component_receive_owner_slot;
  std::vector<std::uint64_t> component_owner_labels;
  std::vector<double> component_send_values;
  std::vector<double> component_receive_values;
  std::vector<double> component_response_send_values;
  std::vector<double> component_response_receive_values;
  std::vector<double> component_owner_sums;
  std::vector<double> component_owner_scales;
  std::vector<std::uint64_t> rank_hashes;
  FreshStartKinematicProjectionRedCertificate red{};
  FreshStartKinematicProjectionInput input{};
  FreshStartKinematicProjectionPreparedCertificate prepared{};
  FreshStartKinematicProjectionSolvedCertificate solved{};
  FreshStartKinematicProjectionCandidateCertificate candidate{};
  NativeCartesianMgPlan native_mg{};
  LinearIdentity linear_identity{};
  LinearOperatorCertificate linear_certificate{};
  PlanFingerprint collective_operator{};
  mutable LinearOperatorFailureProvenance operator_failure{};
  RevisionToken numeric_generation{};
  RevisionToken chi_generation{};
  PlanFingerprint base_input_hash{};
  PlanFingerprint solved_chi_owned_hash{};
  double initial_continuity_maximum{};
  bool bypass_no_immersed{};
  bool candidate_consumed{true};

  Impl() noexcept {
    exact_operator.owner = this;
    jacobi.owner = this;
    continuity_audit.owner = this;
  }

  Status consensus(Status local) const noexcept {
    return services.reductions == nullptr
               ? Status{StatusCode::invalid_plan, kFreshProjectionPlan}
               : services.reductions->consensus(local);
  }

  PlanFingerprint collective_hash(std::uint64_t local,
                                  std::uint64_t domain) noexcept {
    if (rank_hashes.size() != static_cast<std::size_t>(size) ||
        MPI_Allgather(&local, 1, MPI_UINT64_T, rank_hashes.data(), 1,
                      MPI_UINT64_T, communicator) != MPI_SUCCESS)
      return 0U;
    std::uint64_t hash = mix(kFnvOffset, domain);
    hash = mix(hash, static_cast<std::uint64_t>(size));
    for (std::uint64_t value : rank_hashes)
      hash = mix(hash, value);
    return finish(hash);
  }

  Status replay_base_input(std::uint32_t detail) noexcept {
    if (base_input_hash == 0U)
      return {StatusCode::invalid_plan, detail};
    const PlanFingerprint observed = collective_hash(
        finish(hash_projection_input(kFnvOffset, input, patch.cells)),
        UINT64_C(0x6672657368626173));
    if (observed == 0U)
      return {StatusCode::mpi_failure, detail};
    return consensus(observed == base_input_hash
                         ? Status{}
                         : Status{StatusCode::invalid_plan, detail});
  }

  bool cell_active(Int3 cell) const noexcept {
    return source_cells[flat_cell(patch.cells, cell)] != 0U;
  }
  bool face_active(CartesianAxis axis, Int3 face) const noexcept {
    const Int3 extents = face_extents(patch.cells, axis);
    const std::vector<std::uint8_t> &selected =
        axis == CartesianAxis::x
            ? source_x_faces
            : (axis == CartesianAxis::y ? source_y_faces : source_z_faces);
    return selected[flat_face(extents, face)] != 0U;
  }
  bool anchor_gid(std::uint64_t gid) const noexcept {
    return std::binary_search(component_anchors.begin(),
                              component_anchors.end(), gid);
  }
  std::uint64_t gid(Int3 local) const noexcept {
    return global_gid(geometry->global_cells(),
                      {patch.begin.x + local.x, patch.begin.y + local.y,
                       patch.begin.z + local.z});
  }
};

namespace {

template <class Implementation>
Status copy_activity(const FreshStartKinematicProjectionSpec &spec,
                     Implementation &impl) {
  const Int3 cells = spec.patch.cells;
  try {
    impl.source_cells.assign(cell_count(cells), 1U);
    impl.source_x_faces.assign(
        cell_count(face_extents(cells, CartesianAxis::x)), 1U);
    impl.source_y_faces.assign(
        cell_count(face_extents(cells, CartesianAxis::y)), 1U);
    impl.source_z_faces.assign(
        cell_count(face_extents(cells, CartesianAxis::z)), 1U);
    if (!empty_activity(spec.activity)) {
      std::copy_n(spec.activity.cells.data, spec.activity.cells.size,
                  impl.source_cells.data());
      std::copy_n(spec.activity.x_faces.data, spec.activity.x_faces.size,
                  impl.source_x_faces.data());
      std::copy_n(spec.activity.y_faces.data, spec.activity.y_faces.size,
                  impl.source_y_faces.data());
      std::copy_n(spec.activity.z_faces.data, spec.activity.z_faces.size,
                  impl.source_z_faces.data());
    }
  } catch (const std::bad_alloc &) {
    return {StatusCode::allocation_failure, kFreshProjectionPlan};
  }
  impl.source_activity = {
      {impl.source_cells.data(), impl.source_cells.size()},
      {impl.source_x_faces.data(), impl.source_x_faces.size()},
      {impl.source_y_faces.data(), impl.source_y_faces.size()},
      {impl.source_z_faces.data(), impl.source_z_faces.size()},
      spec.activity.local_fingerprint == 0U ? 1U
                                            : spec.activity.local_fingerprint,
      spec.activity.collective_fingerprint == 0U
          ? PlanFingerprint{1U}
          : spec.activity.collective_fingerprint};
  return {};
}

template <class Implementation>
Status compile_components_scalable(Implementation &impl) {
  const Int3 cells = impl.patch.cells;
  const Int3 global = impl.geometry->global_cells();
  try {
    const std::size_t local_cell_count = cell_count(cells);
    DisjointSet sets(local_cell_count);
    Status local{};
    std::uint64_t local_active_count = 0U;
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x)
          if (impl.cell_active({x, y, z}))
            ++local_active_count;

    // The supplied cell/face activity is itself part of the SPD proof.  An
    // active Cartesian face must have two active endpoints whenever both
    // endpoints are owned locally.  This catches either directed form of the
    // formerly asymmetric "active face into inactive cell" defect.
    for (CartesianAxis axis :
         {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
      const Int3 extents = face_extents(cells, axis);
      const std::int32_t extent = axis_value(cells, axis);
      for (std::int32_t z = 0; z < extents.z; ++z)
        for (std::int32_t y = 0; y < extents.y; ++y)
          for (std::int32_t x = 0; x < extents.x; ++x) {
            const Int3 face{x, y, z};
            const std::int32_t normal = axis_value(face, axis);
            if (normal <= 0 || normal >= extent ||
                !impl.face_active(axis, face))
              continue;
            Int3 left = face;
            Int3 right = face;
            set_axis_value(left, axis, normal - 1);
            set_axis_value(right, axis, normal);
            if (!impl.cell_active(left) || !impl.cell_active(right)) {
              local = {StatusCode::invalid_plan, kFreshProjectionPlan};
              continue;
            }
            sets.unite(flat_cell(cells, left), flat_cell(cells, right));
          }
    }
    local = impl.consensus(local);
    if (!local)
      return local;

    const auto boundary_area = [&](CartesianAxis axis) noexcept {
      return axis == CartesianAxis::x
                 ? static_cast<std::size_t>(cells.y) * cells.z
                 : (axis == CartesianAxis::y
                        ? static_cast<std::size_t>(cells.x) * cells.z
                        : static_cast<std::size_t>(cells.x) * cells.y);
    };
    const auto boundary_cell = [&](CartesianAxis axis, bool high,
                                   std::size_t index) noexcept {
      if (axis == CartesianAxis::x)
        return Int3{high ? cells.x - 1 : 0,
                    static_cast<std::int32_t>(index % cells.y),
                    static_cast<std::int32_t>(index / cells.y)};
      if (axis == CartesianAxis::y)
        return Int3{static_cast<std::int32_t>(index % cells.x),
                    high ? cells.y - 1 : 0,
                    static_cast<std::int32_t>(index / cells.x)};
      return Int3{static_cast<std::int32_t>(index % cells.x),
                  static_cast<std::int32_t>(index / cells.x),
                  high ? cells.z - 1 : 0};
    };
    const auto boundary_face = [&](CartesianAxis axis, bool high,
                                   std::size_t index) noexcept {
      Int3 face = boundary_cell(axis, high, index);
      set_axis_value(face, axis, high ? axis_value(cells, axis) : 0);
      return face;
    };
    const auto rank_from_coord = [&](Int3 coordinate) noexcept {
      return coordinate.x +
             impl.patch.process_grid.x *
                 (coordinate.y + impl.patch.process_grid.y * coordinate.z);
    };
    constexpr int kInvalidNeighbor = std::numeric_limits<int>::min();
    const auto neighbor = [&](CartesianAxis axis, bool high) noexcept {
      const Int3 face = boundary_face(axis, high, 0U);
      PressureCorrectionFaceRule rule;
      if (!impl.pressure_boundary.face_rule(axis, face, rule))
        return kInvalidNeighbor;
      Int3 coordinate = impl.patch.process_coord;
      std::int32_t selected = axis_value(coordinate, axis);
      const std::int32_t count = axis_value(impl.patch.process_grid, axis);
      if (!rule.physical) {
        selected += high ? 1 : -1;
        if (selected < 0 || selected >= count)
          return kInvalidNeighbor;
      } else if (rule.kind == PressureCorrectionFaceKind::periodic) {
        selected = high ? 0 : count - 1;
      } else {
        return MPI_PROC_NULL;
      }
      set_axis_value(coordinate, axis, selected);
      return rank_from_coord(coordinate);
    };

    // A one-rank periodic direction is still a two-sided interface.  Validate
    // both stored face masks and join its endpoint cells before local labels
    // are compiled.
    for (CartesianAxis axis :
         {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
      const int low_neighbor = neighbor(axis, false);
      const int high_neighbor = neighbor(axis, true);
      if (low_neighbor == kInvalidNeighbor ||
          high_neighbor == kInvalidNeighbor) {
        local = {StatusCode::invalid_plan, kFreshProjectionPlan};
        continue;
      }
      const std::size_t area = boundary_area(axis);
      if (low_neighbor == impl.rank && high_neighbor == impl.rank) {
        for (std::size_t index = 0U; index < area; ++index) {
          const Int3 low_face = boundary_face(axis, false, index);
          const Int3 high_face = boundary_face(axis, true, index);
          const bool low_active = impl.face_active(axis, low_face);
          const bool high_active = impl.face_active(axis, high_face);
          const Int3 low_cell = boundary_cell(axis, false, index);
          const Int3 high_cell = boundary_cell(axis, true, index);
          if (low_active != high_active ||
              (low_active &&
               (!impl.cell_active(low_cell) || !impl.cell_active(high_cell)))) {
            local = {StatusCode::invalid_plan, kFreshProjectionPlan};
          } else if (low_active) {
            sets.unite(flat_cell(cells, low_cell), flat_cell(cells, high_cell));
          }
        }
      }
      for (bool high : {false, true}) {
        if ((high ? high_neighbor : low_neighbor) != MPI_PROC_NULL)
          continue;
        for (std::size_t index = 0U; index < area; ++index) {
          const Int3 face = boundary_face(axis, high, index);
          if (impl.face_active(axis, face) &&
              !impl.cell_active(boundary_cell(axis, high, index)))
            local = {StatusCode::invalid_plan, kFreshProjectionPlan};
        }
      }
    }
    local = impl.consensus(local);
    if (!local)
      return local;

    std::vector<std::uint64_t> root_label(local_cell_count, UINT64_MAX);
    std::vector<std::uint8_t> root_dirichlet(local_cell_count, 0U);
    std::vector<std::uint8_t> root_distributed(local_cell_count, 0U);
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          if (!impl.cell_active(cell))
            continue;
          const std::size_t root = sets.find(flat_cell(cells, cell));
          root_label[root] = std::min(root_label[root], impl.gid(cell));
        }
    for (CartesianAxis axis :
         {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z})
      for (bool high : {false, true}) {
        if (neighbor(axis, high) != MPI_PROC_NULL)
          continue;
        const std::size_t area = boundary_area(axis);
        for (std::size_t index = 0U; index < area; ++index) {
          const Int3 face = boundary_face(axis, high, index);
          if (!impl.face_active(axis, face))
            continue;
          PressureCorrectionFaceRule rule;
          if (!impl.pressure_boundary.face_rule(axis, face, rule)) {
            local = {StatusCode::invalid_plan, kFreshProjectionPlan};
            continue;
          }
          if (rule.kind == PressureCorrectionFaceKind::homogeneous_dirichlet) {
            const Int3 cell = boundary_cell(axis, high, index);
            root_dirichlet[sets.find(flat_cell(cells, cell))] = 1U;
          }
        }
      }
    local = impl.consensus(local);
    if (!local)
      return local;

    std::vector<std::uint64_t> send_low;
    std::vector<std::uint64_t> send_high;
    std::vector<std::uint64_t> receive_from_plus;
    std::vector<std::uint64_t> receive_from_minus;
    const auto exchange_axis = [&](CartesianAxis axis, std::size_t width,
                                   int tag) -> Status {
      const std::size_t area = boundary_area(axis);
      if (area > static_cast<std::size_t>(INT_MAX) / width)
        return {StatusCode::invalid_plan, kFreshProjectionPlan};
      const int count = static_cast<int>(area * width);
      const int minus = neighbor(axis, false);
      const int plus = neighbor(axis, true);
      if (minus == kInvalidNeighbor || plus == kInvalidNeighbor)
        return {StatusCode::invalid_plan, kFreshProjectionPlan};
      receive_from_plus.assign(area * width, UINT64_MAX);
      receive_from_minus.assign(area * width, UINT64_MAX);
      if (MPI_Sendrecv(send_low.data(), count, MPI_UINT64_T, minus, tag,
                       receive_from_plus.data(), count, MPI_UINT64_T, plus, tag,
                       impl.communicator, MPI_STATUS_IGNORE) != MPI_SUCCESS ||
          MPI_Sendrecv(send_high.data(), count, MPI_UINT64_T, plus, tag + 1,
                       receive_from_minus.data(), count, MPI_UINT64_T, minus,
                       tag + 1, impl.communicator,
                       MPI_STATUS_IGNORE) != MPI_SUCCESS)
        return {StatusCode::mpi_failure, kFreshProjectionPlan};
      return {};
    };

    // Exchange only the six partition/periodic surfaces.  No owned interior
    // cell or purely local component is published.
    for (CartesianAxis axis :
         {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
      const std::size_t area = boundary_area(axis);
      if (neighbor(axis, false) == MPI_PROC_NULL &&
          neighbor(axis, true) == MPI_PROC_NULL)
        continue;
      send_low.assign(3U * area, UINT64_MAX);
      send_high.assign(3U * area, UINT64_MAX);
      for (std::size_t index = 0U; index < area; ++index)
        for (bool high : {false, true}) {
          std::vector<std::uint64_t> &selected = high ? send_high : send_low;
          const Int3 cell = boundary_cell(axis, high, index);
          const bool active = impl.cell_active(cell);
          selected[3U * index] = active ? 1U : 0U;
          selected[3U * index + 1U] =
              impl.face_active(axis, boundary_face(axis, high, index)) ? 1U
                                                                       : 0U;
          selected[3U * index + 2U] =
              active ? root_label[sets.find(flat_cell(cells, cell))]
                     : UINT64_MAX;
        }
      local = exchange_axis(axis, 3U, 7300 + 4 * static_cast<int>(axis));
      local = impl.consensus(local);
      if (!local)
        return local;
      const auto validate_side = [&](bool high,
                                     const std::vector<std::uint64_t> &remote,
                                     int remote_rank) noexcept {
        if (remote_rank == MPI_PROC_NULL)
          return;
        const std::vector<std::uint64_t> &own = high ? send_high : send_low;
        for (std::size_t index = 0U; index < area; ++index) {
          const bool own_cell = own[3U * index] != 0U;
          const bool own_face = own[3U * index + 1U] != 0U;
          const bool remote_cell = remote[3U * index] != 0U;
          const bool remote_face = remote[3U * index + 1U] != 0U;
          if (own_face != remote_face ||
              (own_face && (!own_cell || !remote_cell ||
                            remote[3U * index + 2U] == UINT64_MAX))) {
            local = {StatusCode::invalid_plan, kFreshProjectionPlan};
            continue;
          }
          if (own_face) {
            const Int3 cell = boundary_cell(axis, high, index);
            root_distributed[sets.find(flat_cell(cells, cell))] = 1U;
          }
        }
      };
      validate_side(true, receive_from_plus, neighbor(axis, true));
      validate_side(false, receive_from_minus, neighbor(axis, false));
      local = impl.consensus(local);
      if (!local)
        return local;
    }

    std::uint64_t global_active_count = 0U;
    if (MPI_Allreduce(&local_active_count, &global_active_count, 1,
                      MPI_UINT64_T, MPI_SUM, impl.communicator) != MPI_SUCCESS)
      return {StatusCode::mpi_failure, kFreshProjectionPlan};
    if (global_active_count == 0U)
      return {StatusCode::invalid_plan, kFreshProjectionPlan};

    bool global_changed = true;
    std::uint64_t propagation_steps = 0U;
    while (global_changed && propagation_steps <= global_active_count) {
      bool changed = false;
      for (CartesianAxis axis :
           {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
        const std::size_t area = boundary_area(axis);
        if (neighbor(axis, false) == MPI_PROC_NULL &&
            neighbor(axis, true) == MPI_PROC_NULL)
          continue;
        send_low.assign(2U * area, UINT64_MAX);
        send_high.assign(2U * area, UINT64_MAX);
        for (std::size_t index = 0U; index < area; ++index)
          for (bool high : {false, true}) {
            std::vector<std::uint64_t> &selected = high ? send_high : send_low;
            const Int3 cell = boundary_cell(axis, high, index);
            if (!impl.cell_active(cell) ||
                !impl.face_active(axis, boundary_face(axis, high, index)))
              continue;
            const std::size_t root = sets.find(flat_cell(cells, cell));
            selected[2U * index] = root_label[root];
            selected[2U * index + 1U] = root_dirichlet[root];
          }
        local = exchange_axis(axis, 2U, 7340 + 4 * static_cast<int>(axis));
        local = impl.consensus(local);
        if (!local)
          return local;
        const auto merge_side = [&](bool high,
                                    const std::vector<std::uint64_t> &remote,
                                    int remote_rank) noexcept {
          if (remote_rank == MPI_PROC_NULL)
            return;
          for (std::size_t index = 0U; index < area; ++index) {
            const Int3 face = boundary_face(axis, high, index);
            if (!impl.face_active(axis, face))
              continue;
            const Int3 cell = boundary_cell(axis, high, index);
            const std::size_t root = sets.find(flat_cell(cells, cell));
            const std::uint64_t remote_label = remote[2U * index];
            if (remote_label == UINT64_MAX) {
              local = {StatusCode::invalid_plan, kFreshProjectionPlan};
              continue;
            }
            const std::uint64_t next_label =
                std::min(root_label[root], remote_label);
            const std::uint8_t next_dirichlet = static_cast<std::uint8_t>(
                root_dirichlet[root] != 0U || remote[2U * index + 1U] != 0U);
            changed = changed || next_label != root_label[root] ||
                      next_dirichlet != root_dirichlet[root];
            root_label[root] = next_label;
            root_dirichlet[root] = next_dirichlet;
          }
        };
        merge_side(true, receive_from_plus, neighbor(axis, true));
        merge_side(false, receive_from_minus, neighbor(axis, false));
        local = impl.consensus(local);
        if (!local)
          return local;
      }
      int local_changed = changed ? 1 : 0;
      int collective_changed = 0;
      if (MPI_Allreduce(&local_changed, &collective_changed, 1, MPI_INT,
                        MPI_LOR, impl.communicator) != MPI_SUCCESS)
        return {StatusCode::mpi_failure, kFreshProjectionPlan};
      global_changed = collective_changed != 0;
      ++propagation_steps;
    }
    if (global_changed)
      return {StatusCode::invalid_plan, kFreshProjectionPlan};

    impl.component_labels.clear();
    for (std::size_t flat = 0U; flat < local_cell_count; ++flat)
      if (sets.find(flat) == flat && root_label[flat] != UINT64_MAX)
        impl.component_labels.push_back(root_label[flat]);
    std::sort(impl.component_labels.begin(), impl.component_labels.end());
    impl.component_labels.erase(
        std::unique(impl.component_labels.begin(), impl.component_labels.end()),
        impl.component_labels.end());
    if (impl.component_labels.empty())
      return {StatusCode::invalid_plan, kFreshProjectionPlan};
    impl.component_dirichlet.assign(impl.component_labels.size(), 0U);
    impl.component_distributed.assign(impl.component_labels.size(), 0U);
    std::unordered_map<std::uint64_t, std::uint32_t> component_index;
    component_index.reserve(2U * impl.component_labels.size() + 1U);
    for (std::size_t component = 0U; component < impl.component_labels.size();
         ++component)
      component_index.emplace(impl.component_labels[component],
                              static_cast<std::uint32_t>(component));
    for (std::size_t flat = 0U; flat < local_cell_count; ++flat) {
      if (sets.find(flat) != flat || root_label[flat] == UINT64_MAX)
        continue;
      const std::uint32_t component = component_index[root_label[flat]];
      impl.component_dirichlet[component] =
          static_cast<std::uint8_t>(impl.component_dirichlet[component] != 0U ||
                                    root_dirichlet[flat] != 0U);
      impl.component_distributed[component] = static_cast<std::uint8_t>(
          impl.component_distributed[component] != 0U ||
          root_distributed[flat] != 0U);
    }
    impl.local_components.assign(local_cell_count, kNoComponent);
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          if (!impl.cell_active(cell))
            continue;
          const std::size_t flat = flat_cell(cells, cell);
          impl.local_components[flat] =
              component_index[root_label[sets.find(flat)]];
        }
    impl.component_anchors.clear();
    for (std::size_t component = 0U; component < impl.component_labels.size();
         ++component)
      if (impl.component_dirichlet[component] == 0U)
        impl.component_anchors.push_back(impl.component_labels[component]);

    const auto owner_axis = [](std::int32_t index, std::int32_t extent,
                               std::int32_t partitions) noexcept {
      const std::int32_t quotient = extent / partitions;
      const std::int32_t remainder = extent % partitions;
      const std::int32_t threshold = (quotient + 1) * remainder;
      return index < threshold ? index / (quotient + 1)
                               : remainder + (index - threshold) / quotient;
    };
    const auto owner_rank = [&](std::uint64_t gid) noexcept {
      const std::uint64_t xy = static_cast<std::uint64_t>(global.x) * global.y;
      const std::int32_t z = static_cast<std::int32_t>(gid / xy);
      gid -= static_cast<std::uint64_t>(z) * xy;
      const std::int32_t y = static_cast<std::int32_t>(gid / global.x);
      const std::int32_t x = static_cast<std::int32_t>(
          gid - static_cast<std::uint64_t>(y) * global.x);
      return rank_from_coord(
          {owner_axis(x, global.x, impl.patch.process_grid.x),
           owner_axis(y, global.y, impl.patch.process_grid.y),
           owner_axis(z, global.z, impl.patch.process_grid.z)});
    };

    std::uint64_t local_owned_components = 0U;
    std::uint64_t local_owned_anchors = 0U;
    std::uint64_t local_component_xor = 0U;
    std::uint64_t local_component_sum = 0U;
    for (std::size_t component = 0U; component < impl.component_labels.size();
         ++component) {
      if (owner_rank(impl.component_labels[component]) != impl.rank)
        continue;
      ++local_owned_components;
      if (impl.component_dirichlet[component] == 0U)
        ++local_owned_anchors;
      const std::uint64_t item =
          finish(mix(impl.component_labels[component],
                     impl.component_dirichlet[component]));
      local_component_xor ^= item;
      local_component_sum += item;
    }
    std::array<std::uint64_t, 2U> local_counts{local_owned_components,
                                               local_owned_anchors};
    std::array<std::uint64_t, 2U> global_counts{};
    if (MPI_Allreduce(local_counts.data(), global_counts.data(), 2,
                      MPI_UINT64_T, MPI_SUM, impl.communicator) != MPI_SUCCESS)
      return {StatusCode::mpi_failure, kFreshProjectionPlan};
    if (global_counts[0U] == 0U || global_counts[0U] > UINT32_MAX ||
        global_counts[1U] > UINT32_MAX)
      return {StatusCode::invalid_plan, kFreshProjectionPlan};
    std::uint64_t global_component_xor = 0U;
    std::uint64_t global_component_sum = 0U;
    if (MPI_Allreduce(&local_component_xor, &global_component_xor, 1,
                      MPI_UINT64_T, MPI_BXOR,
                      impl.communicator) != MPI_SUCCESS ||
        MPI_Allreduce(&local_component_sum, &global_component_sum, 1,
                      MPI_UINT64_T, MPI_SUM, impl.communicator) != MPI_SUCCESS)
      return {StatusCode::mpi_failure, kFreshProjectionPlan};

    impl.component_send_counts.assign(static_cast<std::size_t>(impl.size), 0);
    impl.component_send_displacements.assign(
        static_cast<std::size_t>(impl.size), 0);
    impl.component_receive_counts.assign(static_cast<std::size_t>(impl.size),
                                         0);
    impl.component_receive_displacements.assign(
        static_cast<std::size_t>(impl.size), 0);
    for (std::size_t component = 0U; component < impl.component_labels.size();
         ++component)
      if (impl.component_distributed[component] != 0U &&
          impl.component_dirichlet[component] == 0U)
        ++impl.component_send_counts[static_cast<std::size_t>(
            owner_rank(impl.component_labels[component]))];
    int send_total = 0;
    for (int rank = 0; rank < impl.size; ++rank) {
      impl.component_send_displacements[static_cast<std::size_t>(rank)] =
          send_total;
      send_total += impl.component_send_counts[static_cast<std::size_t>(rank)];
    }
    if (MPI_Alltoall(impl.component_send_counts.data(), 1, MPI_INT,
                     impl.component_receive_counts.data(), 1, MPI_INT,
                     impl.communicator) != MPI_SUCCESS)
      return {StatusCode::mpi_failure, kFreshProjectionPlan};
    int receive_total = 0;
    for (int rank = 0; rank < impl.size; ++rank) {
      impl.component_receive_displacements[static_cast<std::size_t>(rank)] =
          receive_total;
      receive_total +=
          impl.component_receive_counts[static_cast<std::size_t>(rank)];
    }
    impl.component_send_labels.resize(static_cast<std::size_t>(send_total));
    impl.component_receive_labels.resize(
        static_cast<std::size_t>(receive_total));
    impl.component_route_local_index.clear();
    impl.component_route_send_slot.clear();
    impl.component_route_local_index.reserve(
        static_cast<std::size_t>(send_total));
    impl.component_route_send_slot.reserve(
        static_cast<std::size_t>(send_total));
    std::vector<int> cursor = impl.component_send_displacements;
    for (std::size_t component = 0U; component < impl.component_labels.size();
         ++component) {
      if (impl.component_distributed[component] == 0U ||
          impl.component_dirichlet[component] != 0U)
        continue;
      const int owner = owner_rank(impl.component_labels[component]);
      const std::uint32_t slot =
          static_cast<std::uint32_t>(cursor[static_cast<std::size_t>(owner)]++);
      impl.component_route_local_index.push_back(
          static_cast<std::uint32_t>(component));
      impl.component_route_send_slot.push_back(slot);
      impl.component_send_labels[slot] = impl.component_labels[component];
    }
    if (MPI_Alltoallv(impl.component_send_labels.data(),
                      impl.component_send_counts.data(),
                      impl.component_send_displacements.data(), MPI_UINT64_T,
                      impl.component_receive_labels.data(),
                      impl.component_receive_counts.data(),
                      impl.component_receive_displacements.data(), MPI_UINT64_T,
                      impl.communicator) != MPI_SUCCESS)
      return {StatusCode::mpi_failure, kFreshProjectionPlan};
    impl.component_owner_labels = impl.component_receive_labels;
    std::sort(impl.component_owner_labels.begin(),
              impl.component_owner_labels.end());
    impl.component_owner_labels.erase(
        std::unique(impl.component_owner_labels.begin(),
                    impl.component_owner_labels.end()),
        impl.component_owner_labels.end());
    impl.component_receive_owner_slot.resize(
        impl.component_receive_labels.size());
    for (std::size_t index = 0U; index < impl.component_receive_labels.size();
         ++index) {
      if (owner_rank(impl.component_receive_labels[index]) != impl.rank)
        return {StatusCode::invalid_plan, kFreshProjectionPlan};
      impl.component_receive_owner_slot[index] = static_cast<std::uint32_t>(
          std::lower_bound(impl.component_owner_labels.begin(),
                           impl.component_owner_labels.end(),
                           impl.component_receive_labels[index]) -
          impl.component_owner_labels.begin());
    }
    impl.component_send_value_counts.resize(
        static_cast<std::size_t>(impl.size));
    impl.component_send_value_displacements.resize(
        static_cast<std::size_t>(impl.size));
    impl.component_receive_value_counts.resize(
        static_cast<std::size_t>(impl.size));
    impl.component_receive_value_displacements.resize(
        static_cast<std::size_t>(impl.size));
    for (int rank = 0; rank < impl.size; ++rank) {
      const std::size_t index = static_cast<std::size_t>(rank);
      if (impl.component_send_counts[index] > INT_MAX / 2 ||
          impl.component_receive_counts[index] > INT_MAX / 2)
        return {StatusCode::invalid_plan, kFreshProjectionPlan};
      impl.component_send_value_counts[index] =
          2 * impl.component_send_counts[index];
      impl.component_send_value_displacements[index] =
          2 * impl.component_send_displacements[index];
      impl.component_receive_value_counts[index] =
          2 * impl.component_receive_counts[index];
      impl.component_receive_value_displacements[index] =
          2 * impl.component_receive_displacements[index];
    }
    impl.component_send_values.resize(2U * send_total);
    impl.component_receive_values.resize(2U * receive_total);
    impl.component_response_send_values.resize(2U * receive_total);
    impl.component_response_receive_values.resize(2U * send_total);
    impl.component_owner_sums.resize(impl.component_owner_labels.size());
    impl.component_owner_scales.resize(impl.component_owner_labels.size());

    // Exact byte count of every simultaneously live, module-owned component
    // communication buffer, using actual vector capacities rather than a
    // formula over logical records.  The local DSU/topology is deliberately
    // excluded: it is neither transmitted nor replicated on another rank.
    std::uint64_t communication_bytes = 0U;
    bool payload_fits = true;
    const auto account = [&](std::size_t capacity,
                             std::size_t element_bytes) noexcept {
      if (capacity > UINT64_MAX / element_bytes ||
          communication_bytes >
              UINT64_MAX -
                  static_cast<std::uint64_t>(capacity) * element_bytes) {
        payload_fits = false;
        return;
      }
      communication_bytes +=
          static_cast<std::uint64_t>(capacity) * element_bytes;
    };
    for (const auto *values :
         {&send_low, &send_high, &receive_from_plus, &receive_from_minus,
          &impl.component_send_labels, &impl.component_receive_labels,
          &impl.component_owner_labels, &impl.rank_hashes})
      account(values->capacity(), sizeof(std::uint64_t));
    for (const auto *values :
         {&impl.component_send_values, &impl.component_receive_values,
          &impl.component_response_send_values,
          &impl.component_response_receive_values, &impl.component_owner_sums,
          &impl.component_owner_scales})
      account(values->capacity(), sizeof(double));
    for (const auto *values :
         {&impl.component_route_local_index, &impl.component_route_send_slot,
          &impl.component_receive_owner_slot})
      account(values->capacity(), sizeof(std::uint32_t));
    for (const auto *values :
         {&cursor, &impl.component_send_counts,
          &impl.component_send_displacements, &impl.component_receive_counts,
          &impl.component_receive_displacements,
          &impl.component_send_value_counts,
          &impl.component_send_value_displacements,
          &impl.component_receive_value_counts,
          &impl.component_receive_value_displacements})
      account(values->capacity(), sizeof(int));
    if (!payload_fits || communication_bytes > UINT64_MAX - 7U)
      return {StatusCode::invalid_plan, kFreshProjectionPlan};
    const std::uint64_t local_peak_payload = (communication_bytes + 7U) / 8U;
    std::uint64_t maximum_peak_payload = 0U;
    std::uint64_t global_peak_payload = 0U;
    if (MPI_Allreduce(&local_peak_payload, &maximum_peak_payload, 1,
                      MPI_UINT64_T, MPI_MAX,
                      impl.communicator) != MPI_SUCCESS ||
        MPI_Allreduce(&local_peak_payload, &global_peak_payload, 1,
                      MPI_UINT64_T, MPI_SUM, impl.communicator) != MPI_SUCCESS)
      return {StatusCode::mpi_failure, kFreshProjectionPlan};

    std::uint64_t graph = mix(kFnvOffset, UINT64_C(0x6672657368636f6d));
    graph = mix(graph, global_active_count);
    graph = mix(graph, global_counts[0U]);
    graph = mix(graph, global_counts[1U]);
    graph = mix(graph, global_component_xor);
    graph = mix(graph, global_component_sum);
    impl.red.active_cells = global_active_count;
    impl.red.component_collective_payload_u64 = maximum_peak_payload;
    impl.red.component_collective_global_payload_u64 = global_peak_payload;
    impl.red.component_count = static_cast<std::uint32_t>(global_counts[0U]);
    impl.red.anchored_component_count =
        static_cast<std::uint32_t>(global_counts[1U]);
    impl.red.component_graph = finish(graph);
    impl.red.component_collective_volume_independent = true;
  } catch (const std::bad_alloc &) {
    return {StatusCode::allocation_failure, kFreshProjectionPlan};
  }
  return {};
}

template <class Implementation>
bool adjacent_global_cells(const Implementation &impl, CartesianAxis axis,
                           Int3 face, std::array<std::uint64_t, 2U> &gids,
                           std::size_t &count) noexcept {
  count = 0U;
  const Int3 global = impl.geometry->global_cells();
  const std::int32_t normal = axis_value(face, axis);
  const std::int32_t local_extent = axis_value(impl.patch.cells, axis);
  const std::int32_t global_face = axis_value(impl.patch.begin, axis) + normal;
  Int3 left_local = face;
  set_axis_value(left_local, axis, normal - 1);
  Int3 right_local = face;
  set_axis_value(right_local, axis, normal);
  const auto add_local = [&](Int3 local) {
    if (local.x < 0 || local.y < 0 || local.z < 0 ||
        local.x >= impl.patch.cells.x || local.y >= impl.patch.cells.y ||
        local.z >= impl.patch.cells.z)
      return;
    gids[count++] = impl.gid(local);
  };
  if (normal > 0)
    add_local(left_local);
  if (normal < local_extent)
    add_local(right_local);
  if (normal != 0 && normal != local_extent)
    return true;
  PressureCorrectionFaceRule rule;
  if (!impl.pressure_boundary.face_rule(axis, face, rule))
    return false;
  if (!rule.physical) {
    Int3 remote = {impl.patch.begin.x + face.x, impl.patch.begin.y + face.y,
                   impl.patch.begin.z + face.z};
    set_axis_value(remote, axis, global_face + (normal == 0 ? -1 : 0));
    gids[count++] = global_gid(global, remote);
  } else if (rule.kind == PressureCorrectionFaceKind::periodic) {
    Int3 remote = {impl.patch.begin.x + face.x, impl.patch.begin.y + face.y,
                   impl.patch.begin.z + face.z};
    set_axis_value(remote, axis,
                   global_face == 0 ? axis_value(global, axis) - 1 : 0);
    gids[count++] = global_gid(global, remote);
  }
  return count != 0U && count <= 2U;
}

template <class Implementation>
Status build_anchored_activity(Implementation &impl) {
  const Int3 cells = impl.patch.cells;
  try {
    impl.mg_cells.assign(cell_count(cells), 0U);
    impl.mg_x_faces.assign(cell_count(face_extents(cells, CartesianAxis::x)),
                           0U);
    impl.mg_y_faces.assign(cell_count(face_extents(cells, CartesianAxis::y)),
                           0U);
    impl.mg_z_faces.assign(cell_count(face_extents(cells, CartesianAxis::z)),
                           0U);
  } catch (const std::bad_alloc &) {
    return {StatusCode::allocation_failure, kFreshProjectionPlan};
  }
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        impl.mg_cells[flat_cell(cells, cell)] =
            impl.cell_active(cell) && !impl.anchor_gid(impl.gid(cell)) ? 1U
                                                                       : 0U;
      }
  for (CartesianAxis axis :
       {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
    const Int3 extents = face_extents(cells, axis);
    std::vector<std::uint8_t> &target =
        axis == CartesianAxis::x
            ? impl.mg_x_faces
            : (axis == CartesianAxis::y ? impl.mg_y_faces : impl.mg_z_faces);
    for (std::int32_t z = 0; z < extents.z; ++z)
      for (std::int32_t y = 0; y < extents.y; ++y)
        for (std::int32_t x = 0; x < extents.x; ++x) {
          const Int3 face{x, y, z};
          if (!impl.face_active(axis, face))
            continue;
          std::array<std::uint64_t, 2U> adjacent{};
          std::size_t adjacent_count = 0U;
          if (!adjacent_global_cells(impl, axis, face, adjacent,
                                     adjacent_count))
            return {StatusCode::invalid_plan, kFreshProjectionPlan};
          bool enabled = true;
          for (std::size_t index = 0U; index < adjacent_count; ++index)
            enabled = enabled && !impl.anchor_gid(adjacent[index]);
          target[flat_face(extents, face)] = enabled ? 1U : 0U;
        }
  }
  std::uint64_t local = mix(kFnvOffset, UINT64_C(0x66726573686d6761));
  for (std::uint8_t value : impl.mg_cells)
    local = mix(local, value);
  for (std::uint8_t value : impl.mg_x_faces)
    local = mix(local, value);
  for (std::uint8_t value : impl.mg_y_faces)
    local = mix(local, value);
  for (std::uint8_t value : impl.mg_z_faces)
    local = mix(local, value);
  const PlanFingerprint collective =
      impl.collective_hash(finish(local), UINT64_C(0x66726573686d6763));
  if (collective == 0U)
    return {StatusCode::mpi_failure, kFreshProjectionPlan};
  impl.mg_activity = {{impl.mg_cells.data(), impl.mg_cells.size()},
                      {impl.mg_x_faces.data(), impl.mg_x_faces.size()},
                      {impl.mg_y_faces.data(), impl.mg_y_faces.size()},
                      {impl.mg_z_faces.data(), impl.mg_z_faces.size()},
                      finish(local),
                      collective};
  return {};
}

template <class Implementation>
double face_mobility(const Implementation &impl, ConstFieldView density,
                     CartesianAxis axis, Int3 face) noexcept {
  if (!impl.face_active(axis, face))
    return 0.0;
  PressureCorrectionFaceRule rule;
  if (!impl.pressure_boundary.face_rule(axis, face, rule))
    return std::numeric_limits<double>::quiet_NaN();
  if (rule.kind == PressureCorrectionFaceKind::homogeneous_neumann)
    return 0.0;
  const Int3 cells = impl.patch.cells;
  const std::int32_t normal = axis_value(face, axis);
  const std::int32_t extent = axis_value(cells, axis);
  Int3 left = face;
  set_axis_value(left, axis, normal - 1);
  const bool periodic = rule.kind == PressureCorrectionFaceKind::periodic;
  if (rule.physical && !periodic) {
    const Int3 owner = normal == 0 ? face : left;
    const Int3 ghost = normal == 0 ? left : face;
    const double rho_owner = density.unchecked(owner, 0U);
    const double rho_ghost = density.unchecked(ghost, 0U);
    const double distance =
        normal == 0
            ? detail::centre_coordinate(*impl.kernels, axis, 0) -
                  detail::face_coordinate(*impl.kernels, axis, 0)
            : detail::face_coordinate(*impl.kernels, axis, extent) -
                  detail::centre_coordinate(*impl.kernels, axis, extent - 1);
    return 2.0 * detail::face_area(*impl.kernels, axis, face) /
           (distance / rho_owner + distance / rho_ghost);
  }
  const double rho_left = density.unchecked(left, 0U);
  const double rho_right = density.unchecked(face, 0U);
  double left_distance = 0.0;
  double right_distance = 0.0;
  if (rule.physical && periodic) {
    const Span<const double> widths = impl.geometry->axis(axis).widths();
    left_distance = 0.5 * widths.data[widths.size - 1U];
    right_distance = 0.5 * widths.data[0U];
  } else {
    const double location =
        detail::face_coordinate(*impl.kernels, axis, normal);
    left_distance =
        location - detail::centre_coordinate(*impl.kernels, axis, normal - 1);
    right_distance =
        detail::centre_coordinate(*impl.kernels, axis, normal) - location;
  }
  return detail::face_area(*impl.kernels, axis, face) /
         (left_distance / rho_left + right_distance / rho_right);
}

double divergence(ConstFaceFluxView flux, Int3 cell) noexcept {
  return flux.x.unchecked({cell.x + 1, cell.y, cell.z}) -
         flux.x.unchecked(cell) +
         flux.y.unchecked({cell.x, cell.y + 1, cell.z}) -
         flux.y.unchecked(cell) +
         flux.z.unchecked({cell.x, cell.y, cell.z + 1}) -
         flux.z.unchecked(cell);
}

} // namespace

Status FreshStartKinematicProjectionPlan::Impl::ExactOperator::apply(
    FieldView x, FieldView y) const noexcept {
  if (owner != nullptr)
    owner->operator_failure = {};
  if (owner == nullptr || owner->linear_identity.fingerprint == 0U ||
      x.field != owner->services.krylov_field ||
      !detail::valid_cell_view(as_const(x), owner->patch.cells, 0U, 1U, 1U) ||
      !detail::valid_cell_view(y, owner->patch.cells, 0U, 1U) ||
      detail::field_views_overlap(as_const(x), as_const(y)))
    return {StatusCode::invalid_plan, kFreshProjectionSolve};
  std::array<FieldView, 1U> fields{x};
  HaloTicket ticket;
  Status status = owner->services.operator_halo->begin(
      owner->services.operator_halo_stage, {fields.data(), fields.size()},
      ticket);
  if (status)
    status = owner->services.operator_halo->finish(
        ticket, {fields.data(), fields.size()});
  if (!status) {
    owner->operator_failure = {
        status, LinearOperatorStatusScope::collective,
        owner->services.operator_halo->lowest_failing_rank()};
    return status;
  }
  x = fields[0U];
  status = owner->pressure_boundary.fill_ghosts(x);
  if (status &&
      detail::fresh_projection_failure_injected_for_test(
          detail::FreshProjectionFailurePoint::exact_operator_boundary,
          owner->rank))
    status = {StatusCode::invalid_plan, kFreshProjectionSolve};
  status = owner->consensus(status);
  if (!status) {
    owner->operator_failure = {
        status, LinearOperatorStatusScope::collective,
        owner->services.reductions->lowest_failing_rank()};
    return status;
  }
  const ConstFieldView input = as_const(x);
  const auto &workspace = owner->workspace;
  const Int3 cells = owner->patch.cells;
  Status arithmetic{};
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y_index = 0; y_index < cells.y; ++y_index)
      for (std::int32_t x_index = 0; x_index < cells.x; ++x_index) {
        const Int3 cell{x_index, y_index, z};
        double value =
            workspace.diagonal.unchecked(cell, 0U) * input.unchecked(cell, 0U);
        for (CartesianAxis axis :
             {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
          const ConstFaceFieldView coefficient =
              select(as_const(workspace.x_solver_coefficient),
                     as_const(workspace.y_solver_coefficient),
                     as_const(workspace.z_solver_coefficient), axis);
          Int3 high = cell;
          if (axis == CartesianAxis::x)
            ++high.x;
          if (axis == CartesianAxis::y)
            ++high.y;
          if (axis == CartesianAxis::z)
            ++high.z;
          value -=
              coefficient.unchecked(cell) *
              owner->pressure_boundary.neighbor_value(input, cell, axis, -1);
          value -=
              coefficient.unchecked(high) *
              owner->pressure_boundary.neighbor_value(input, cell, axis, 1);
        }
        y.unchecked(cell, 0U) = value;
        if (!std::isfinite(value))
          arithmetic = {StatusCode::numerical_failure,
                        kFreshProjectionNumerical};
      }
  if (!arithmetic)
    owner->operator_failure = {arithmetic,
                               LinearOperatorStatusScope::rank_local, -1};
  return arithmetic;
}

Status FreshStartKinematicProjectionPlan::Impl::JacobiPreconditioner::apply(
    ConstFieldView input, FieldView output, std::uint32_t) noexcept {
  if (owner == nullptr ||
      !detail::valid_cell_view(input, owner->patch.cells, 0U, 1U, 0U) ||
      !detail::valid_cell_view(output, owner->patch.cells, 0U, 1U) ||
      detail::field_views_overlap(input, as_const(output)))
    return {StatusCode::invalid_plan, kFreshProjectionSolve};
  const Int3 cells = owner->patch.cells;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const double diagonal = owner->workspace.diagonal.unchecked(cell, 0U);
        const double value = input.unchecked(cell, 0U) / diagonal;
        if (!std::isfinite(value) || !(diagonal > 0.0))
          return {StatusCode::numerical_failure, kFreshProjectionNumerical};
        output.unchecked(cell, 0U) = value;
      }
  return {};
}

bool FreshStartKinematicProjectionPreparedCertificate::valid() const noexcept {
  return issuer_ != nullptr && lineage_ != 0U && numeric_ != 0U &&
         state_ != 0U && component_count_ != 0U &&
         anchor_count_ <= component_count_ &&
         std::isfinite(maximum_compatibility_defect_) &&
         maximum_compatibility_defect_ >= 0.0;
}

bool FreshStartKinematicProjectionSolvedCertificate::valid() const noexcept {
  return issuer_ != nullptr && prepared_lineage_ != 0U && lineage_ != 0U &&
         chi_revision_ != 0U && result_.status.code == StatusCode::ok &&
         (result_.termination == LinearTermination::converged ||
          result_.termination == LinearTermination::zero_rhs);
}

bool FreshStartKinematicProjectionCandidateCertificate::valid() const noexcept {
  return issuer_ != nullptr && solved_lineage_ != 0U &&
         candidate_state_ != 0U && candidate_flux_ != 0U && lineage_ != 0U &&
         std::isfinite(initial_continuity_maximum_) &&
         initial_continuity_maximum_ >= 0.0 &&
         std::isfinite(final_continuity_maximum_) &&
         final_continuity_maximum_ >= 0.0 &&
         std::isfinite(final_continuity_l2_) && final_continuity_l2_ >= 0.0 &&
         fixed_flux_bitwise_unchanged_ && joint_candidate_;
}

FreshStartKinematicProjectionPlan::
    ~FreshStartKinematicProjectionPlan() noexcept {
  release();
}

FreshStartKinematicProjectionPlan::FreshStartKinematicProjectionPlan(
    FreshStartKinematicProjectionPlan &&other) noexcept
    : implementation_(std::exchange(other.implementation_, nullptr)) {}

FreshStartKinematicProjectionPlan &FreshStartKinematicProjectionPlan::operator=(
    FreshStartKinematicProjectionPlan &&other) noexcept {
  if (this != &other) {
    release();
    implementation_ = std::exchange(other.implementation_, nullptr);
  }
  return *this;
}

void FreshStartKinematicProjectionPlan::release() noexcept {
  delete std::exchange(implementation_, nullptr);
}

Status FreshStartKinematicProjectionPlan::compile(
    const FreshStartKinematicProjectionSpec &spec,
    FreshStartKinematicProjectionServices services,
    FreshStartKinematicProjectionWorkspace workspace,
    FreshStartKinematicProjectionPlan &out) noexcept {
  int rank = -1;
  int size = 0;
  const Int3 cells = spec.patch.cells;
  const bool communicator =
      spec.communicator != MPI_COMM_NULL &&
      MPI_Comm_rank(spec.communicator, &rank) == MPI_SUCCESS &&
      MPI_Comm_size(spec.communicator, &size) == MPI_SUCCESS && rank >= 0 &&
      rank < size && size > 0;
  const bool route_control =
      (spec.route == FreshStartProjectionLinearRoute::native_mg_fgmres &&
       spec.solve.restart != 0U) ||
      (spec.route == FreshStartProjectionLinearRoute::jacobi_pcg_oracle &&
       spec.solve.restart == 0U);
  const bool tolerances =
      std::isfinite(spec.compatibility_absolute_tolerance) &&
      spec.compatibility_absolute_tolerance >= 0.0 &&
      std::isfinite(spec.compatibility_relative_tolerance) &&
      spec.compatibility_relative_tolerance >= 0.0 &&
      std::isfinite(spec.continuity_absolute_tolerance) &&
      spec.continuity_absolute_tolerance >= 0.0 &&
      std::isfinite(spec.continuity_relative_tolerance) &&
      spec.continuity_relative_tolerance >= 0.0;
  const LinearWorkspaceRequirements *requirements =
      services.solver_workspace == nullptr
          ? nullptr
          : &services.solver_workspace->requirements();
  const LinearAlgorithm required_algorithm =
      spec.route == FreshStartProjectionLinearRoute::native_mg_fgmres
          ? LinearAlgorithm::fgmres
          : LinearAlgorithm::pcg;
  const bool linear_services =
      services.operator_halo != nullptr && services.operator_halo->ready() &&
      services.operator_halo_stage != 0U &&
      services.correction_halo != nullptr &&
      services.correction_halo->ready() &&
      services.correction_halo_stage != 0U &&
      services.solver_workspace != nullptr && services.reductions != nullptr &&
      requirements != nullptr &&
      requirements->algorithm == required_algorithm &&
      requirements->ghost_width >= 1U &&
      requirements->maximum_shape.x >= cells.x &&
      requirements->maximum_shape.y >= cells.y &&
      requirements->maximum_shape.z >= cells.z;
  const bool mg_services =
      spec.route != FreshStartProjectionLinearRoute::native_mg_fgmres ||
      (services.mg.finest_halo != nullptr &&
       services.mg.reductions == services.reductions &&
       services.mg.workspace != nullptr);
  const bool valid_workspace =
      detail::valid_cell_view(as_const(workspace.chi), cells, 0U, 1U, 1U) &&
      detail::valid_cell_view(as_const(workspace.rhs), cells, 0U, 1U, 0U) &&
      detail::valid_cell_view(as_const(workspace.diagonal), cells, 0U, 1U,
                              0U) &&
      detail::valid_cell_view(as_const(workspace.candidate_velocity), cells, 0U,
                              3U, 0U) &&
      valid_face_view(workspace.x_physical_mobility, CartesianAxis::x, cells) &&
      valid_face_view(workspace.y_physical_mobility, CartesianAxis::y, cells) &&
      valid_face_view(workspace.z_physical_mobility, CartesianAxis::z, cells) &&
      valid_face_view(workspace.x_solver_coefficient, CartesianAxis::x,
                      cells) &&
      valid_face_view(workspace.y_solver_coefficient, CartesianAxis::y,
                      cells) &&
      valid_face_view(workspace.z_solver_coefficient, CartesianAxis::z,
                      cells) &&
      detail::valid_flux_view(workspace.candidate_mass_flux, cells) &&
      workspace_alias_matrix_valid(workspace);
  bool workspace_disjoint_from_linear =
      valid_workspace && services.solver_workspace != nullptr;
  if (workspace_disjoint_from_linear) {
    for (ConstFieldView field : workspace_cell_views(workspace)) {
      workspace_disjoint_from_linear =
          workspace_disjoint_from_linear &&
          !services.solver_workspace->overlaps_storage(field);
      if (services.mg.workspace != nullptr)
        workspace_disjoint_from_linear =
            workspace_disjoint_from_linear &&
            !services.mg.workspace->overlaps_storage(field);
    }
    for (ConstFaceFieldView face : workspace_face_views(workspace)) {
      workspace_disjoint_from_linear =
          workspace_disjoint_from_linear &&
          !face_overlaps_solver_workspace(face, *services.solver_workspace);
      if (services.mg.workspace != nullptr)
        workspace_disjoint_from_linear =
            workspace_disjoint_from_linear &&
            !face_overlaps_mg_workspace(face, *services.mg.workspace);
    }
    if (services.mg.workspace != nullptr)
      workspace_disjoint_from_linear =
          workspace_disjoint_from_linear &&
          solver_and_mg_workspaces_disjoint(*services.solver_workspace,
                                            *services.mg.workspace);
  }
  if (!communicator || spec.geometry == nullptr ||
      spec.geometry->fingerprint() == 0U || spec.kernels == nullptr ||
      spec.kernels->fingerprint() == 0U || spec.boundary == nullptr ||
      spec.boundary->semantic_fingerprint() == 0U ||
      !same_cells(spec.geometry->global_cells(),
                  {spec.geometry->global_cells().x,
                   spec.geometry->global_cells().y,
                   spec.geometry->global_cells().z}) ||
      !same_cells(spec.kernels->cells(), cells) ||
      !same_cells(spec.boundary->local_cells(), cells) ||
      !valid_activity(spec.activity, cells) || !route_control || !tolerances ||
      !linear_services || !mg_services || !valid_workspace ||
      !workspace_disjoint_from_linear)
    return {StatusCode::invalid_plan, kFreshProjectionPlan};

  const std::array<HaloFieldSpec, 1U> operator_fields{
      {{services.krylov_field, 1U, 1U}}};
  const std::array<HaloFieldSpec, 1U> correction_fields{
      {{workspace.chi.field, 1U, 1U}}};
  if (!services.operator_halo->validate_contract(
          spec.communicator, spec.patch,
          {operator_fields.data(), operator_fields.size()},
          spec.boundary->halo_topology()) ||
      !services.correction_halo->validate_contract(
          spec.communicator, spec.patch,
          {correction_fields.data(), correction_fields.size()},
          spec.boundary->halo_topology()))
    return {StatusCode::invalid_plan, kFreshProjectionPlan};

  auto candidate = std::unique_ptr<Impl>(new (std::nothrow) Impl);
  if (!candidate)
    return {StatusCode::allocation_failure, kFreshProjectionPlan};
  candidate->communicator = spec.communicator;
  candidate->rank = rank;
  candidate->size = size;
  candidate->geometry = spec.geometry;
  candidate->kernels = spec.kernels;
  candidate->boundary = spec.boundary;
  candidate->patch = spec.patch;
  candidate->route = spec.route;
  candidate->solve_control = spec.solve;
  candidate->mg_policy = spec.mg_policy;
  candidate->compatibility_absolute = spec.compatibility_absolute_tolerance;
  candidate->compatibility_relative = spec.compatibility_relative_tolerance;
  candidate->continuity_absolute = spec.continuity_absolute_tolerance;
  candidate->continuity_relative = spec.continuity_relative_tolerance;
  candidate->services = services;
  candidate->workspace = workspace;
  candidate->bypass_no_immersed = empty_activity(spec.activity);
  try {
    candidate->rank_hashes.resize(static_cast<std::size_t>(size));
  } catch (const std::bad_alloc &) {
    return {StatusCode::allocation_failure, kFreshProjectionPlan};
  }
  Status status = PressureCorrectionBoundaryPlan::compile(
      *spec.geometry, spec.patch, *spec.boundary, candidate->pressure_boundary);
  if (status)
    status = copy_activity(spec, *candidate);
  if (status)
    status = compile_components_scalable(*candidate);
  if (status)
    status = build_anchored_activity(*candidate);
  if (!status)
    return status;
  try {
    candidate->component_local_sums.resize(candidate->component_labels.size());
    candidate->component_local_scales.resize(
        candidate->component_labels.size());
  } catch (const std::bad_alloc &) {
    return {StatusCode::allocation_failure, kFreshProjectionPlan};
  }

  std::uint64_t plan_hash = mix(kFnvOffset, UINT64_C(0x6672657368706c6e));
  plan_hash = mix(plan_hash, spec.geometry->fingerprint());
  plan_hash =
      mix(plan_hash, candidate->pressure_boundary.certificate().semantic);
  plan_hash = mix(plan_hash, candidate->red.component_graph);
  plan_hash = mix(plan_hash, candidate->mg_activity.collective_fingerprint);
  plan_hash = mix(plan_hash, static_cast<std::uint8_t>(spec.route));
  candidate->red.plan = finish(plan_hash);
  candidate->red.boundary_semantics =
      candidate->pressure_boundary.certificate().semantic;
  candidate->red.activity_graph = candidate->mg_activity.collective_fingerprint;
  candidate->red.physical_flux_jacobian =
      finish(mix(plan_hash, UINT64_C(0x706879736d6f6269)));
  candidate->red.anchored_schur =
      finish(mix(plan_hash, UINT64_C(0x616e636873636875)));
  candidate->red.physical_mobility_separate = true;
  candidate->red.per_component_compatibility = true;
  candidate->red.global_minimum_gid_anchors = true;
  candidate->red.mg_null_space_none = true;
  candidate->red.ibm_gradient_from_face_activity = true;
  candidate->red.cut_face_zero_mobility = true;
  candidate->red.inactive_velocity_stationary_wall_zero = true;
  candidate->red.runtime_workspace_preallocated = true;
  candidate->red.no_immersed_bitwise_bypass = true;
  candidate->red.three_layer_joint_commit = true;
  candidate->red.chi_writes_pressure = false;
  if (!candidate->red.valid())
    return {StatusCode::invalid_plan, kFreshProjectionPlan};

  out.release();
  out.implementation_ = candidate.release();
  return {};
}

const FreshStartKinematicProjectionRedCertificate &
FreshStartKinematicProjectionPlan::red() const noexcept {
  static const FreshStartKinematicProjectionRedCertificate empty{};
  return implementation_ == nullptr ? empty : implementation_->red;
}

Status FreshStartKinematicProjectionPlan::prepare(
    const FreshStartKinematicProjectionInput &input,
    FreshStartKinematicProjectionPreparedCertificate &certificate) noexcept {
  certificate = {};
  if (implementation_ == nullptr)
    return {StatusCode::invalid_plan, kFreshProjectionPrepare};
  Impl &impl = *implementation_;
  const Int3 cells = impl.patch.cells;
  impl.input = {};
  impl.base_input_hash = 0U;
  impl.solved_chi_owned_hash = 0U;
  impl.prepared = {};
  impl.solved = {};
  impl.candidate = {};
  impl.candidate_consumed = true;
  Status local{};
  if (input.state == 0U ||
      !detail::valid_cell_view(input.density, cells, 0U, 1U, 1U) ||
      !detail::valid_cell_view(input.velocity, cells, 0U, 3U, 0U) ||
      !detail::valid_cell_view(input.velocity_accepted_n_minus_one, cells, 0U,
                               3U, 0U) ||
      !detail::valid_cell_view(input.velocity_trial, cells, 0U, 3U, 0U) ||
      !detail::valid_flux_view(input.mass_flux, cells,
                               input.mass_flux.revision) ||
      !input_alias_matrix_valid(
          input, impl.workspace, impl.services.solver_workspace,
          impl.route == FreshStartProjectionLinearRoute::native_mg_fgmres
              ? impl.services.mg.workspace
              : nullptr))
    local = {StatusCode::invalid_plan, kFreshProjectionPrepare};
  local = impl.consensus(local);
  if (!local)
    return local;

  const PlanFingerprint base_input_hash = impl.collective_hash(
      finish(hash_projection_input(kFnvOffset, input, cells)),
      UINT64_C(0x6672657368626173));
  if (base_input_hash == 0U)
    return {StatusCode::mpi_failure, kFreshProjectionPrepare};
  std::fill(impl.component_local_sums.begin(), impl.component_local_sums.end(),
            0.0);
  std::fill(impl.component_local_scales.begin(),
            impl.component_local_scales.end(), 0.0);
  double local_initial_maximum = 0.0;
  for (std::int32_t z = 0; z < cells.z && local; ++z)
    for (std::int32_t y = 0; y < cells.y && local; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        const std::size_t flat = flat_cell(cells, cell);
        if (!impl.cell_active(cell)) {
          impl.workspace.rhs.unchecked(cell, 0U) = 0.0;
          impl.workspace.diagonal.unchecked(cell, 0U) = 1.0;
          continue;
        }
        const double density = input.density.unchecked(cell, 0U);
        const double residual =
            impl.bypass_no_immersed ? 0.0 : divergence(input.mass_flux, cell);
        const double rhs = -residual;
        if (!std::isfinite(density) || !(density > 0.0) ||
            !std::isfinite(rhs)) {
          local = {StatusCode::numerical_failure, kFreshProjectionNumerical};
          break;
        }
        const std::uint32_t component = impl.local_components[flat];
        if (component == kNoComponent ||
            component >= impl.component_labels.size()) {
          local = {StatusCode::invalid_plan, kFreshProjectionPrepare};
          break;
        }
        impl.component_local_sums[component] += rhs;
        impl.component_local_scales[component] += std::abs(rhs);
        local_initial_maximum =
            std::max(local_initial_maximum, std::abs(residual));
        impl.workspace.rhs.unchecked(cell, 0U) =
            impl.anchor_gid(impl.gid(cell)) ? 0.0 : rhs;
      }
  local = impl.consensus(local);
  if (!local)
    return local;

  for (CartesianAxis axis :
       {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
    FaceFieldView physical = select(impl.workspace.x_physical_mobility,
                                    impl.workspace.y_physical_mobility,
                                    impl.workspace.z_physical_mobility, axis);
    FaceFieldView solver = select(impl.workspace.x_solver_coefficient,
                                  impl.workspace.y_solver_coefficient,
                                  impl.workspace.z_solver_coefficient, axis);
    const std::vector<std::uint8_t> &mg_faces =
        axis == CartesianAxis::x
            ? impl.mg_x_faces
            : (axis == CartesianAxis::y ? impl.mg_y_faces : impl.mg_z_faces);
    for (std::int32_t z = 0; z < physical.extents.z && local; ++z)
      for (std::int32_t y = 0; y < physical.extents.y && local; ++y)
        for (std::int32_t x = 0; x < physical.extents.x; ++x) {
          const Int3 face{x, y, z};
          const double value =
              impl.bypass_no_immersed
                  ? 0.0
                  : face_mobility(impl, input.density, axis, face);
          if (!std::isfinite(value) || value < 0.0) {
            local = {StatusCode::numerical_failure, kFreshProjectionNumerical};
            break;
          }
          physical.unchecked(face) = value;
          solver.unchecked(face) =
              mg_faces[flat_face(physical.extents, face)] != 0U ? value : 0.0;
        }
  }
  if (local) {
    for (std::int32_t z = 0; z < cells.z && local; ++z)
      for (std::int32_t y = 0; y < cells.y && local; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          if (!impl.cell_active(cell) || impl.anchor_gid(impl.gid(cell))) {
            impl.workspace.diagonal.unchecked(cell, 0U) = 1.0;
            continue;
          }
          double diagonal = impl.bypass_no_immersed ? 1.0 : 0.0;
          for (CartesianAxis axis :
               {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
            const ConstFaceFieldView physical =
                select(as_const(impl.workspace.x_physical_mobility),
                       as_const(impl.workspace.y_physical_mobility),
                       as_const(impl.workspace.z_physical_mobility), axis);
            Int3 high = cell;
            if (axis == CartesianAxis::x)
              ++high.x;
            if (axis == CartesianAxis::y)
              ++high.y;
            if (axis == CartesianAxis::z)
              ++high.z;
            if (!impl.bypass_no_immersed) {
              // Match the certified MG hierarchy's left-to-right face sum.
              // Pairwise grouping can round one ulp below the identical face
              // sum and spuriously violate the exact SDD check.
              diagonal += physical.unchecked(cell);
              diagonal += physical.unchecked(high);
            }
          }
          if (!std::isfinite(diagonal) || !(diagonal > 0.0)) {
            local = {StatusCode::numerical_failure, kFreshProjectionNumerical};
            break;
          }
          impl.workspace.diagonal.unchecked(cell, 0U) = diagonal;
        }
  }
  local = impl.consensus(local);
  if (!local)
    return local;

  std::fill(impl.component_send_values.begin(),
            impl.component_send_values.end(), 0.0);
  double local_maximum_defect = 0.0;
  bool incompatible_component = false;
  for (std::size_t route = 0U; route < impl.component_route_local_index.size();
       ++route) {
    const std::uint32_t component = impl.component_route_local_index[route];
    const std::uint32_t slot = impl.component_route_send_slot[route];
    if (slot == kNoComponent || 2U * static_cast<std::size_t>(slot) + 1U >=
                                    impl.component_send_values.size()) {
      local = {StatusCode::invalid_plan, kFreshProjectionPrepare};
      continue;
    }
    impl.component_send_values[2U * slot] =
        impl.component_local_sums[component];
    impl.component_send_values[2U * slot + 1U] =
        impl.component_local_scales[component];
  }
  for (std::size_t component = 0U; component < impl.component_labels.size();
       ++component) {
    if (impl.component_dirichlet[component] != 0U ||
        impl.component_distributed[component] != 0U)
      continue;
    const double sum = impl.component_local_sums[component];
    const double scale = impl.component_local_scales[component];
    const double limit =
        impl.compatibility_absolute + impl.compatibility_relative * scale;
    local_maximum_defect = std::max(local_maximum_defect, std::abs(sum));
    incompatible_component = incompatible_component || !std::isfinite(sum) ||
                             !std::isfinite(scale) || std::abs(sum) > limit;
  }
  local = impl.consensus(local);
  if (!local)
    return local;
  local = MPI_Alltoallv(impl.component_send_values.data(),
                        impl.component_send_value_counts.data(),
                        impl.component_send_value_displacements.data(),
                        MPI_DOUBLE, impl.component_receive_values.data(),
                        impl.component_receive_value_counts.data(),
                        impl.component_receive_value_displacements.data(),
                        MPI_DOUBLE, impl.communicator) == MPI_SUCCESS
              ? Status{}
              : Status{StatusCode::mpi_failure, kFreshProjectionCompatibility};
  local = impl.consensus(local);
  if (!local)
    return local;
  std::fill(impl.component_owner_sums.begin(), impl.component_owner_sums.end(),
            0.0);
  std::fill(impl.component_owner_scales.begin(),
            impl.component_owner_scales.end(), 0.0);
  for (std::size_t record = 0U; record < impl.component_receive_labels.size();
       ++record) {
    const std::uint32_t owner_slot = impl.component_receive_owner_slot[record];
    const double sum = impl.component_receive_values[2U * record];
    const double scale = impl.component_receive_values[2U * record + 1U];
    if (owner_slot >= impl.component_owner_sums.size() || !std::isfinite(sum) ||
        !std::isfinite(scale) || scale < 0.0) {
      local = {StatusCode::numerical_failure, kFreshProjectionCompatibility};
      continue;
    }
    impl.component_owner_sums[owner_slot] += sum;
    impl.component_owner_scales[owner_slot] += scale;
  }
  local = impl.consensus(local);
  if (!local)
    return local;
  for (std::size_t record = 0U; record < impl.component_receive_labels.size();
       ++record) {
    const std::uint32_t owner_slot = impl.component_receive_owner_slot[record];
    impl.component_response_send_values[2U * record] =
        impl.component_owner_sums[owner_slot];
    impl.component_response_send_values[2U * record + 1U] =
        impl.component_owner_scales[owner_slot];
  }
  local =
      MPI_Alltoallv(impl.component_response_send_values.data(),
                    impl.component_receive_value_counts.data(),
                    impl.component_receive_value_displacements.data(),
                    MPI_DOUBLE, impl.component_response_receive_values.data(),
                    impl.component_send_value_counts.data(),
                    impl.component_send_value_displacements.data(), MPI_DOUBLE,
                    impl.communicator) == MPI_SUCCESS
          ? Status{}
          : Status{StatusCode::mpi_failure, kFreshProjectionCompatibility};
  local = impl.consensus(local);
  if (!local)
    return local;
  for (std::size_t route = 0U; route < impl.component_route_local_index.size();
       ++route) {
    const std::uint32_t component = impl.component_route_local_index[route];
    const std::size_t slot = impl.component_route_send_slot[route];
    const double sum = impl.component_response_receive_values[2U * slot];
    const double scale = impl.component_response_receive_values[2U * slot + 1U];
    const double limit =
        impl.compatibility_absolute + impl.compatibility_relative * scale;
    local_maximum_defect = std::max(local_maximum_defect, std::abs(sum));
    incompatible_component = incompatible_component || !std::isfinite(sum) ||
                             !std::isfinite(scale) || std::abs(sum) > limit;
  }
  double maximum_defect = 0.0;
  local = impl.services.reductions->checked_max(
      {&local_maximum_defect, 1U}, {&maximum_defect, 1U},
      incompatible_component
          ? Status{StatusCode::rejected_step, kFreshProjectionCompatibility}
          : Status{});
  if (!local)
    return local;

  double global_initial_maximum = 0.0;
  local = impl.services.reductions->checked_max({&local_initial_maximum, 1U},
                                                {&global_initial_maximum, 1U});
  if (!local)
    return local;
  impl.initial_continuity_maximum = global_initial_maximum;

  ++impl.numeric_generation;
  if (impl.numeric_generation == 0U)
    ++impl.numeric_generation;
  std::uint64_t numeric_local = mix(kFnvOffset, input.state);
  numeric_local = hash_field_values(numeric_local, input.density, cells, 1U);
  numeric_local = hash_flux_values(numeric_local, input.mass_flux);
  numeric_local = hash_field_values(
      numeric_local, as_const(impl.workspace.diagonal), cells, 1U);
  const PlanFingerprint numeric =
      impl.collective_hash(finish(numeric_local), UINT64_C(0x66726573686e756d));
  if (numeric == 0U)
    return {StatusCode::mpi_failure, kFreshProjectionPrepare};
  impl.linear_identity.symbolic = impl.red.anchored_schur;
  impl.linear_identity.numeric = numeric;
  impl.linear_identity.hierarchy =
      finish(mix(impl.red.activity_graph, impl.numeric_generation));
  impl.linear_identity.workspace =
      impl.services.solver_workspace->fingerprint();
  std::uint64_t identity = mix(kFnvOffset, impl.linear_identity.symbolic);
  identity = mix(identity, impl.linear_identity.numeric);
  identity = mix(identity, impl.linear_identity.hierarchy);
  identity = mix(identity, impl.linear_identity.workspace);
  impl.linear_identity.fingerprint = finish(identity);
  impl.collective_operator = finish(mix(impl.red.anchored_schur, numeric));
  impl.linear_certificate = {impl.linear_identity, impl.collective_operator,
                             cells, LinearOperatorClass::spd};

  if (impl.route == FreshStartProjectionLinearRoute::native_mg_fgmres &&
      !impl.bypass_no_immersed) {
    MgBoundarySet boundaries;
    local = impl.pressure_boundary.mg_boundaries(boundaries);
    if (local && detail::fresh_projection_failure_injected_for_test(
                     detail::FreshProjectionFailurePoint::prepare_mg_boundary,
                     impl.rank))
      local = {StatusCode::invalid_plan, kFreshProjectionPrepare};
    local = impl.consensus(local);
    if (!local)
      return local;
    NativeCartesianMgSpec mg_spec;
    mg_spec.communicator = impl.communicator;
    mg_spec.geometry = impl.geometry;
    mg_spec.patch = impl.patch;
    mg_spec.boundaries = boundaries;
    mg_spec.null_space = MgNullSpace::none;
    mg_spec.operator_class =
        MgOperatorClass::symmetric_diagonally_dominant_m_matrix;
    mg_spec.policy = impl.mg_policy;
    mg_spec.correction_scaling = MgCorrectionScaling::residual_minimizing;
    mg_spec.identity = impl.linear_identity;
    mg_spec.coefficients = {impl.numeric_generation, numeric, 0.0};
    // Identity rows and zero incident coefficients already encode inactive
    // cells and anchors in the exact numeric graph.  Passing a second
    // activity mask makes Native MG eliminate those identity rows and changes
    // the preconditioned problem.  Empty activity therefore means "retain the
    // complete anchored coefficient graph", not "ignore IBM".
    mg_spec.activity = {};
    const MgCoefficientViews coefficients{
        as_const(impl.workspace.diagonal),
        as_const(impl.workspace.x_solver_coefficient),
        as_const(impl.workspace.y_solver_coefficient),
        as_const(impl.workspace.z_solver_coefficient)};
    local = NativeCartesianMgPlan::compile(mg_spec, impl.services.mg,
                                           coefficients, impl.native_mg);
    if (!local)
      return local;
  }

  std::uint64_t lineage = mix(impl.red.plan, numeric);
  lineage = mix(lineage, base_input_hash);
  lineage = mix(lineage, input.state);
  FreshStartKinematicProjectionPreparedCertificate issued;
  issued.issuer_ = &impl;
  issued.lineage_ = finish(lineage);
  issued.numeric_ = numeric;
  issued.state_ = input.state;
  issued.component_count_ = impl.red.component_count;
  issued.anchor_count_ = impl.red.anchored_component_count;
  issued.maximum_compatibility_defect_ = maximum_defect;
  if (!issued.valid())
    return {StatusCode::invalid_plan, kFreshProjectionPrepare};
  impl.input = input;
  impl.base_input_hash = base_input_hash;
  impl.prepared = issued;
  impl.solved = {};
  impl.candidate = {};
  impl.candidate_consumed = true;
  certificate = issued;
  return {};
}

Status FreshStartKinematicProjectionPlan::solve(
    const FreshStartKinematicProjectionPreparedCertificate &prepared,
    FreshStartKinematicProjectionSolvedCertificate &certificate,
    ResourceCounters *resources) noexcept {
  certificate = {};
  if (implementation_ == nullptr)
    return {StatusCode::invalid_plan, kFreshProjectionSolve};
  Impl &impl = *implementation_;
  Status local{};
  if (!prepared.valid() || prepared.issuer_ != &impl ||
      prepared.lineage_ != impl.prepared.lineage_ ||
      prepared.numeric_ != impl.prepared.numeric_ ||
      impl.linear_identity.fingerprint == 0U)
    local = {StatusCode::invalid_plan, kFreshProjectionSolve};
  local = impl.consensus(local);
  if (!local)
    return local;
  local = impl.replay_base_input(kFreshProjectionSolve);
  if (!local)
    return local;
  const Int3 cells = impl.patch.cells;
  for (std::int32_t z = -impl.workspace.chi.ghosts.z;
       z < cells.z + impl.workspace.chi.ghosts.z; ++z)
    for (std::int32_t y = -impl.workspace.chi.ghosts.y;
         y < cells.y + impl.workspace.chi.ghosts.y; ++y)
      for (std::int32_t x = -impl.workspace.chi.ghosts.x;
           x < cells.x + impl.workspace.chi.ghosts.x; ++x)
        impl.workspace.chi.unchecked({x, y, z}, 0U) = 0.0;
  LinearSolveResult result;
  if (impl.bypass_no_immersed) {
    result.status = {};
    result.termination = LinearTermination::zero_rhs;
    result.initial_true_residual = 0.0;
    result.final_true_residual = 0.0;
    result.recursive_residual = 0.0;
  } else {
    const LinearSolveInvocation invocation{
        as_const(impl.workspace.rhs), impl.workspace.chi, impl.linear_identity,
        impl.solve_control,
        impl.route == FreshStartProjectionLinearRoute::native_mg_fgmres &&
                impl.continuity_absolute + impl.continuity_relative *
                                               impl.initial_continuity_maximum >
                    0.0
            ? &impl.continuity_audit
            : nullptr};
    result = impl.route == FreshStartProjectionLinearRoute::native_mg_fgmres
                 ? solve_fgmres(impl.exact_operator, impl.native_mg, invocation,
                                *impl.services.solver_workspace,
                                *impl.services.reductions, resources)
                 : solve_pcg(impl.exact_operator, impl.jacobi, invocation,
                             *impl.services.solver_workspace,
                             *impl.services.reductions, resources);
  }
  if (!result.status || (result.termination != LinearTermination::converged &&
                         result.termination != LinearTermination::zero_rhs))
    return result.status
               ? Status{StatusCode::rejected_step, kFreshProjectionSolve}
               : result.status;
  ++impl.chi_generation;
  if (impl.chi_generation == 0U)
    ++impl.chi_generation;
  impl.workspace.chi.revision = impl.chi_generation;
  std::uint64_t local_hash =
      hash_field_values(mix(kFnvOffset, prepared.lineage_),
                        as_const(impl.workspace.chi), cells, 1U);
  const PlanFingerprint chi_hash =
      impl.collective_hash(finish(local_hash), UINT64_C(0x6672657368636869));
  const PlanFingerprint chi_owned_hash = impl.collective_hash(
      finish(hash_field_numeric_values(kFnvOffset, as_const(impl.workspace.chi),
                                       cells, 1U)),
      UINT64_C(0x667265736863686f));
  if (chi_hash == 0U || chi_owned_hash == 0U)
    return {StatusCode::mpi_failure, kFreshProjectionSolve};
  FreshStartKinematicProjectionSolvedCertificate issued;
  issued.issuer_ = &impl;
  issued.prepared_lineage_ = prepared.lineage_;
  issued.lineage_ = finish(mix(prepared.lineage_, chi_hash));
  issued.chi_revision_ = impl.workspace.chi.revision;
  issued.result_ = result;
  if (!issued.valid())
    return {StatusCode::invalid_plan, kFreshProjectionSolve};
  impl.solved_chi_owned_hash = chi_owned_hash;
  impl.solved = issued;
  impl.candidate = {};
  certificate = issued;
  return {};
}

namespace {

template <class Implementation>
Status form_projected_flux(Implementation &impl, ConstFieldView chi,
                           bool &fixed_flux_unchanged) noexcept {
  Status local;
  fixed_flux_unchanged = true;
  for (CartesianAxis axis :
       {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
    const ConstFaceFieldView input = select(impl.input.mass_flux, axis);
    const ConstFaceFieldView mobility =
        select(as_const(impl.workspace.x_physical_mobility),
               as_const(impl.workspace.y_physical_mobility),
               as_const(impl.workspace.z_physical_mobility), axis);
    const FaceFieldView output =
        select(impl.workspace.candidate_mass_flux, axis);
    for (std::int32_t z = 0; z < output.extents.z && local; ++z)
      for (std::int32_t y = 0; y < output.extents.y && local; ++y)
        for (std::int32_t x = 0; x < output.extents.x; ++x) {
          const Int3 face{x, y, z};
          const double base = input.unchecked(face);
          double value = base;
          const double coefficient = mobility.unchecked(face);
          if (coefficient != 0.0)
            value += impl.pressure_boundary.mass_flux_response(chi, axis, face,
                                                               coefficient);
          if (!std::isfinite(value)) {
            local = {StatusCode::numerical_failure, kFreshProjectionNumerical};
            break;
          }
          output.unchecked(face) = value;
          if (coefficient == 0.0)
            fixed_flux_unchanged =
                fixed_flux_unchanged &&
                std::memcmp(&output.unchecked(face), &base, sizeof(base)) == 0;
        }
  }
  return local;
}

template <class Implementation>
double masked_gradient_component(const Implementation &impl, ConstFieldView chi,
                                 Int3 cell, CartesianAxis axis) noexcept {
  Int3 low_face = cell;
  Int3 high_face = cell;
  if (axis == CartesianAxis::x)
    ++high_face.x;
  if (axis == CartesianAxis::y)
    ++high_face.y;
  if (axis == CartesianAxis::z)
    ++high_face.z;
  const bool low = impl.face_active(axis, low_face);
  const bool high = impl.face_active(axis, high_face);
  if (!low && !high)
    return 0.0;
  const double centre = chi.unchecked(cell, 0U);
  if (low && high) {
    Int3 minus = cell;
    Int3 plus = cell;
    if (axis == CartesianAxis::x) {
      --minus.x;
      ++plus.x;
    } else if (axis == CartesianAxis::y) {
      --minus.y;
      ++plus.y;
    } else {
      --minus.z;
      ++plus.z;
    }
    const detail::DerivativeWeights weights =
        detail::derivative_weights(*impl.kernels, axis, axis_value(cell, axis));
    return weights.minus * chi.unchecked(minus, 0U) + weights.centre * centre +
           weights.plus * chi.unchecked(plus, 0U);
  }
  const int direction = high ? 1 : -1;
  const double neighbour =
      impl.pressure_boundary.neighbor_value(chi, cell, axis, direction);
  const double centre_location =
      detail::centre_coordinate(*impl.kernels, axis, axis_value(cell, axis));
  const double neighbour_location = detail::centre_coordinate(
      *impl.kernels, axis, axis_value(cell, axis) + direction);
  const double distance = std::abs(neighbour_location - centre_location);
  return high ? (neighbour - centre) / distance
              : (centre - neighbour) / distance;
}

} // namespace

Status FreshStartKinematicProjectionPlan::Impl::ContinuityAudit::evaluate(
    ConstFieldView solution, ConstFieldView, ReductionEngine &reductions,
    LinearConvergenceAuditResult &result) noexcept {
  Impl &impl = *owner;
  const Int3 cells = impl.patch.cells;
  // The anchored algebraic L2 norm does not bound continuity at the removed
  // gauge row.  Test the actual unanchored flux response before Krylov exits.
  // This only stages chi/phi; the accepted U/phi transaction remains untouched.
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        impl.workspace.chi.unchecked({x, y, z}, 0U) =
            solution.unchecked({x, y, z}, 0U);
  std::array<FieldView, 1U> fields{impl.workspace.chi};
  HaloTicket ticket;
  Status status = impl.services.correction_halo->begin(
      impl.services.correction_halo_stage, {fields.data(), fields.size()},
      ticket);
  if (status)
    status = impl.services.correction_halo->finish(
        ticket, {fields.data(), fields.size()});
  if (!status) return status;
  impl.workspace.chi = fields[0U];
  status = impl.pressure_boundary.fill_ghosts(impl.workspace.chi);
  status = reductions.consensus(status);
  if (!status) return status;
  bool fixed_flux_unchanged = true;
  status = form_projected_flux(impl, as_const(impl.workspace.chi),
                               fixed_flux_unchanged);
  if (status && !fixed_flux_unchanged)
    status = {StatusCode::invalid_plan, kFreshProjectionAudit};
  double local_maximum = 0.0;
  for (std::int32_t z = 0; z < cells.z && status; ++z)
    for (std::int32_t y = 0; y < cells.y && status; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        if (!impl.cell_active(cell)) continue;
        const double residual =
            divergence(as_const(impl.workspace.candidate_mass_flux), cell);
        if (!std::isfinite(residual)) {
          status = {StatusCode::numerical_failure, kFreshProjectionNumerical};
          break;
        }
        local_maximum = std::max(local_maximum, std::abs(residual));
      }
  double maximum = 0.0;
  status = reductions.checked_max({&local_maximum, 1U}, {&maximum, 1U}, status);
  if (!status) return status;
  result = {};
  result.metric = maximum;
  result.unscaled_metric = maximum;
  result.limit = impl.continuity_absolute +
                 impl.continuity_relative * impl.initial_continuity_maximum;
  result.accepted = maximum <= result.limit;
  return {};
}

Status FreshStartKinematicProjectionPlan::audit(
    const FreshStartKinematicProjectionSolvedCertificate &solved,
    FreshStartKinematicProjectionCandidateCertificate &certificate) noexcept {
  certificate = {};
  if (implementation_ == nullptr)
    return {StatusCode::invalid_plan, kFreshProjectionAudit};
  Impl &impl = *implementation_;
  Status local{};
  if (!solved.valid() || solved.issuer_ != &impl ||
      solved.lineage_ != impl.solved.lineage_ ||
      solved.prepared_lineage_ != impl.prepared.lineage_ ||
      solved.chi_revision_ != impl.workspace.chi.revision)
    local = {StatusCode::invalid_plan, kFreshProjectionAudit};
  local = impl.consensus(local);
  if (!local)
    return local;
  local = impl.replay_base_input(kFreshProjectionAudit);
  if (!local)
    return local;
  const PlanFingerprint observed_chi = impl.collective_hash(
      finish(hash_field_values(mix(kFnvOffset, solved.prepared_lineage_),
                               as_const(impl.workspace.chi), impl.patch.cells,
                               1U)),
      UINT64_C(0x6672657368636869));
  local = observed_chi != 0U && finish(mix(solved.prepared_lineage_,
                                           observed_chi)) == solved.lineage_
              ? Status{}
              : Status{StatusCode::invalid_plan, kFreshProjectionAudit};
  local = impl.consensus(local);
  if (!local)
    return local;

  std::array<FieldView, 1U> fields{impl.workspace.chi};
  HaloTicket ticket;
  local = impl.services.correction_halo->begin(
      impl.services.correction_halo_stage, {fields.data(), fields.size()},
      ticket);
  if (local)
    local = impl.services.correction_halo->finish(
        ticket, {fields.data(), fields.size()});
  if (!local)
    return local;
  impl.workspace.chi = fields[0U];
  local = impl.pressure_boundary.fill_ghosts(impl.workspace.chi);
  if (local &&
      detail::fresh_projection_failure_injected_for_test(
          detail::FreshProjectionFailurePoint::audit_boundary, impl.rank))
    local = {StatusCode::invalid_plan, kFreshProjectionAudit};
  local = impl.consensus(local);
  if (!local)
    return local;
  const PlanFingerprint post_halo_chi = impl.collective_hash(
      finish(hash_field_numeric_values(kFnvOffset, as_const(impl.workspace.chi),
                                       impl.patch.cells, 1U)),
      UINT64_C(0x667265736863686f));
  local = post_halo_chi != 0U && post_halo_chi == impl.solved_chi_owned_hash
              ? Status{}
              : Status{StatusCode::invalid_plan, kFreshProjectionAudit};
  local = impl.consensus(local);
  if (!local)
    return local;
  const ConstFieldView chi = as_const(impl.workspace.chi);
  const Int3 cells = impl.patch.cells;

  bool fixed_flux_unchanged = true;
  local = form_projected_flux(impl, chi, fixed_flux_unchanged);
  for (std::int32_t z = 0; z < cells.z && local; ++z)
    for (std::int32_t y = 0; y < cells.y && local; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          double value = 0.0;
          if (impl.cell_active(cell)) {
            value = impl.input.velocity.unchecked(cell, component);
            value -= masked_gradient_component(
                impl, chi, cell, static_cast<CartesianAxis>(component));
          }
          if (!std::isfinite(value)) {
            local = {StatusCode::numerical_failure, kFreshProjectionNumerical};
            break;
          }
          impl.workspace.candidate_velocity.unchecked(cell, component) = value;
        }
      }
  local = impl.consensus(local);
  if (!local)
    return local;

  double local_maximum = 0.0;
  double local_square = 0.0;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x) {
        const Int3 cell{x, y, z};
        if (!impl.cell_active(cell) || impl.bypass_no_immersed)
          continue;
        const double residual =
            divergence(as_const(impl.workspace.candidate_mass_flux), cell);
        if (!std::isfinite(residual)) {
          local = {StatusCode::numerical_failure, kFreshProjectionNumerical};
          break;
        }
        local_maximum = std::max(local_maximum, std::abs(residual));
        local_square += residual * residual;
      }
  std::array<double, 2U> local_metrics{local_square, 0.0};
  std::array<double, 2U> global_metrics{};
  local = impl.services.reductions->checked_sum(
      {local_metrics.data(), 1U}, {global_metrics.data(), 1U}, local);
  if (!local)
    return local;
  local = impl.services.reductions->checked_max({&local_maximum, 1U},
                                                {&global_metrics[1U], 1U});
  if (!local)
    return local;
  const double final_l2 = std::sqrt(global_metrics[0U]);
  const double final_maximum = global_metrics[1U];
  const double limit =
      impl.continuity_absolute +
      impl.continuity_relative * impl.initial_continuity_maximum;
  if (!fixed_flux_unchanged || !std::isfinite(final_l2) ||
      !std::isfinite(final_maximum) || final_maximum > limit)
    local = {StatusCode::rejected_step, kFreshProjectionAudit};
  local = impl.consensus(local);
  if (!local)
    return local;

  const PlanFingerprint state_hash = impl.collective_hash(
      finish(hash_field_values(mix(kFnvOffset, solved.lineage_),
                               as_const(impl.workspace.candidate_velocity),
                               cells, 3U)),
      UINT64_C(0x6672657368637374));
  const PlanFingerprint flux_hash = impl.collective_hash(
      finish(hash_flux_values(mix(kFnvOffset, solved.lineage_),
                              as_const(impl.workspace.candidate_mass_flux))),
      UINT64_C(0x667265736863666c));
  if (state_hash == 0U || flux_hash == 0U)
    return {StatusCode::mpi_failure, kFreshProjectionAudit};
  FreshStartKinematicProjectionCandidateCertificate issued;
  issued.issuer_ = &impl;
  issued.solved_lineage_ = solved.lineage_;
  issued.candidate_state_ = state_hash;
  issued.candidate_flux_ = flux_hash;
  issued.lineage_ = finish(mix(mix(solved.lineage_, state_hash), flux_hash));
  issued.initial_continuity_maximum_ = impl.initial_continuity_maximum;
  issued.final_continuity_maximum_ = final_maximum;
  issued.final_continuity_l2_ = final_l2;
  issued.fixed_flux_bitwise_unchanged_ = fixed_flux_unchanged;
  issued.joint_candidate_ = true;
  if (!issued.valid())
    return {StatusCode::invalid_plan, kFreshProjectionAudit};
  impl.candidate = issued;
  impl.candidate_consumed = false;
  certificate = issued;
  return {};
}

Status FreshStartKinematicProjectionPlan::commit(
    const FreshStartKinematicProjectionCandidateCertificate &candidate,
    Span<FieldView> velocity_layers, FaceFluxView mass_flux) noexcept {
  if (implementation_ == nullptr)
    return {StatusCode::invalid_plan, kFreshProjectionCommit};
  Impl &impl = *implementation_;
  const Int3 cells = impl.patch.cells;
  Status local{};
  std::array<ConstFieldView, 3U> output_velocity{};
  std::array<ConstFaceFieldView, 3U> output_flux{};
  if (velocity_layers.data != nullptr && velocity_layers.size == 3U) {
    for (std::size_t layer = 0U; layer < output_velocity.size(); ++layer)
      output_velocity[layer] = as_const(velocity_layers.data[layer]);
  }
  output_flux = {as_const(mass_flux.x), as_const(mass_flux.y),
                 as_const(mass_flux.z)};
  bool valid_output =
      velocity_layers.data != nullptr && velocity_layers.size == 3U &&
      detail::valid_flux_view(as_const(mass_flux), cells, mass_flux.revision);
  for (ConstFieldView layer : output_velocity)
    valid_output =
        valid_output && detail::valid_cell_view(layer, cells, 0U, 3U, 0U);
  valid_output = valid_output && disjoint_fields(output_velocity) &&
                 disjoint_faces(output_flux) &&
                 disjoint_cells_and_faces(output_velocity, output_flux);
  const auto scratch_fields = workspace_cell_views(impl.workspace);
  const auto scratch_faces = workspace_face_views(impl.workspace);
  for (ConstFieldView layer : output_velocity) {
    for (ConstFieldView scratch : scratch_fields)
      valid_output =
          valid_output && !detail::field_views_overlap(layer, scratch);
    for (ConstFaceFieldView scratch : scratch_faces)
      valid_output =
          valid_output && !detail::cell_face_views_overlap(layer, scratch);
  }
  for (ConstFaceFieldView face : output_flux) {
    for (ConstFieldView scratch : scratch_fields)
      valid_output =
          valid_output && !detail::cell_face_views_overlap(scratch, face);
    for (ConstFaceFieldView scratch : scratch_faces)
      valid_output = valid_output && !detail::face_views_overlap(face, scratch);
  }
  if (!valid_output || impl.candidate_consumed || !candidate.valid() ||
      candidate.issuer_ != &impl ||
      candidate.lineage_ != impl.candidate.lineage_ ||
      candidate.solved_lineage_ != impl.solved.lineage_ ||
      !same_field_binding(output_velocity[0U], impl.input.velocity) ||
      !same_field_binding(output_velocity[1U],
                          impl.input.velocity_accepted_n_minus_one) ||
      !same_field_binding(output_velocity[2U], impl.input.velocity_trial) ||
      !same_flux_binding(as_const(mass_flux), impl.input.mass_flux))
    local = {StatusCode::invalid_plan, kFreshProjectionCommit};
  local = impl.consensus(local);
  if (!local)
    return local;
  local = impl.replay_base_input(kFreshProjectionCommit);
  if (!local)
    return local;
  const PlanFingerprint state_hash = impl.collective_hash(
      finish(hash_field_values(mix(kFnvOffset, impl.solved.lineage_),
                               as_const(impl.workspace.candidate_velocity),
                               cells, 3U)),
      UINT64_C(0x6672657368637374));
  const PlanFingerprint flux_hash = impl.collective_hash(
      finish(hash_flux_values(mix(kFnvOffset, impl.solved.lineage_),
                              as_const(impl.workspace.candidate_mass_flux))),
      UINT64_C(0x667265736863666c));
  if (candidate.candidate_state_ != state_hash ||
      candidate.candidate_flux_ != flux_hash || state_hash == 0U ||
      flux_hash == 0U)
    local = {StatusCode::invalid_plan, kFreshProjectionCommit};
  local = impl.consensus(local);
  if (!local)
    return local;

  if (impl.bypass_no_immersed) {
    impl.candidate_consumed = true;
    impl.candidate = {};
    impl.solved = {};
    impl.prepared = {};
    impl.input = {};
    impl.base_input_hash = 0U;
    impl.solved_chi_owned_hash = 0U;
    return {};
  }

  // No-fail region.  The collective validation above is the last operation
  // that can return a failure.  U and phi are now published as one joint
  // candidate; pressure is not addressable by this module.
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x)
        for (std::uint8_t component = 0U; component < 3U; ++component)
          for (std::size_t layer = 0U; layer < velocity_layers.size; ++layer)
            velocity_layers.data[layer].unchecked({x, y, z}, component) =
                impl.workspace.candidate_velocity.unchecked({x, y, z},
                                                            component);
  for (CartesianAxis axis :
       {CartesianAxis::x, CartesianAxis::y, CartesianAxis::z}) {
    FaceFieldView destination = select(mass_flux, axis);
    const ConstFaceFieldView source =
        select(as_const(impl.workspace.candidate_mass_flux), axis);
    for (std::int32_t z = 0; z < destination.extents.z; ++z)
      for (std::int32_t y = 0; y < destination.extents.y; ++y)
        for (std::int32_t x = 0; x < destination.extents.x; ++x)
          destination.unchecked({x, y, z}) = source.unchecked({x, y, z});
  }
  impl.candidate_consumed = true;
  impl.candidate = {};
  impl.solved = {};
  impl.prepared = {};
  impl.input = {};
  impl.base_input_hash = 0U;
  impl.solved_chi_owned_hash = 0U;
  return {};
}

} // namespace hundun::v04

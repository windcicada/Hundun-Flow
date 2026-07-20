// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow/momentum_predictor.hpp"

#include "hundun/finite_volume/cell_centered_fvm.hpp"
#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace hundun::flow {
namespace {

using mesh::EntityOwnership;
using mesh::FaceAxis;
using mesh::FaceSide;
using mesh::GlobalFaceId;
using mesh::LocalCellId;
using mesh::LocalFaceId;
using runtime::Error;
using runtime::Int3;
using runtime::Real3;

bool same(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool finite(Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

Real3 subtract(Real3 left, Real3 right) noexcept {
  return {left.x - right.x, left.y - right.y, left.z - right.z};
}

Real3 add(Real3 left, Real3 right) noexcept {
  return {left.x + right.x, left.y + right.y, left.z + right.z};
}

Real3 multiply(double scale, Real3 value) noexcept {
  return {scale * value.x, scale * value.y, scale * value.z};
}

double dot(Real3 left, Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

Real3 cross(Real3 left, Real3 right) noexcept {
  return {left.y * right.z - left.z * right.y,
          left.z * right.x - left.x * right.z,
          left.x * right.y - left.y * right.x};
}

double norm(Real3 value) noexcept {
  const double scale =
      std::max({std::abs(value.x), std::abs(value.y), std::abs(value.z)});
  if (scale == 0.0)
    return 0.0;
  if (!std::isfinite(scale))
    return std::numeric_limits<double>::infinity();
  const double x = value.x / scale;
  const double y = value.y / scale;
  const double z = value.z / scale;
  return scale * std::sqrt(x * x + y * y + z * z);
}

std::array<Int3, 4> face_vertices(mesh::LogicalFace face) {
  const Int3 c = face.coordinate;
  switch (face.axis) {
  case FaceAxis::x:
    return {Int3{c.x, c.y, c.z}, Int3{c.x, c.y + 1, c.z},
            Int3{c.x, c.y + 1, c.z + 1}, Int3{c.x, c.y, c.z + 1}};
  case FaceAxis::y:
    return {Int3{c.x, c.y, c.z}, Int3{c.x, c.y, c.z + 1},
            Int3{c.x + 1, c.y, c.z + 1}, Int3{c.x + 1, c.y, c.z}};
  case FaceAxis::z:
    return {Int3{c.x, c.y, c.z}, Int3{c.x + 1, c.y, c.z},
            Int3{c.x + 1, c.y + 1, c.z}, Int3{c.x, c.y + 1, c.z}};
  }
  throw Error("invalid momentum face axis");
}

std::array<Real3, 4>
periodic_canonical_vertices(const mesh::MeshGeometry &geometry,
                            mesh::LogicalFace logical) {
  if (logical.axis == FaceAxis::x)
    logical.coordinate.x = 0;
  if (logical.axis == FaceAxis::y)
    logical.coordinate.y = 0;
  if (logical.axis == FaceAxis::z)
    logical.coordinate.z = 0;
  const auto coordinates = face_vertices(logical);
  std::array<Real3, 4> result{};
  for (std::size_t index = 0; index < result.size(); ++index) {
    result[index] = geometry.vertex_position_m(coordinates[index]);
  }
  return result;
}

Real3 canonical_face_center(const std::array<Real3, 4> &vertices,
                            bool uniform) noexcept {
  const Real3 anchor = vertices[0];
  if (uniform)
    return add(anchor, multiply(0.5, subtract(vertices[2], anchor)));
  Real3 center = anchor;
  for (std::size_t index = 1; index < vertices.size(); ++index) {
    center = add(center, multiply(0.25, subtract(vertices[index], anchor)));
  }
  return center;
}

Real3 positive_axis_area(const std::array<Real3, 4> &vertices) noexcept {
  const Real3 first = multiply(0.5, cross(subtract(vertices[1], vertices[0]),
                                          subtract(vertices[2], vertices[0])));
  const Real3 second = multiply(0.5, cross(subtract(vertices[2], vertices[0]),
                                           subtract(vertices[3], vertices[0])));
  return add(first, second);
}

bool solve_success(linear::SolveTerminationReason reason) noexcept {
  return reason == linear::SolveTerminationReason::converged ||
         reason == linear::SolveTerminationReason::zero_right_hand_side;
}

void validate_host_context(const execution::ExecutionContext &context) {
  if (context.backend_identity() == 0U ||
      context.space() != execution::ExecutionSpace::host ||
      !context.supports(execution::ExecutionCapability::host_access) ||
      !context.supports(execution::ExecutionCapability::buffer_allocation)) {
    throw Error("momentum predictor requires a live host execution context");
  }
}

template <class T>
void validate_vector_view(execution::VectorView<T> view,
                          const execution::ExecutionContext &context,
                          std::size_t expected_size, bool require_writable,
                          const char *name) {
  if (view.size() != expected_size || view.stride() != 1U ||
      view.backend_identity() != context.backend_identity() ||
      view.space() != execution::ExecutionSpace::host ||
      (require_writable && !view.writable())) {
    throw Error(std::string("momentum ") + name + " view is incompatible");
  }
  static_cast<void>(view.data());
}

struct ViewRange final {
  execution::AllocationIdentity allocation{};
  std::size_t begin{};
  std::size_t end{};
  bool writable{};
};

template <class T> ViewRange view_range(execution::VectorView<T> view) {
  if (view.size() >
      (std::numeric_limits<std::size_t>::max() - view.offset_bytes()) /
          sizeof(double)) {
    throw Error("momentum vector byte range overflows");
  }
  return {view.allocation_identity(), view.offset_bytes(),
          view.offset_bytes() + view.size() * sizeof(double), view.writable()};
}

bool overlaps(const ViewRange &left, const ViewRange &right) noexcept {
  return left.allocation == right.allocation && left.begin < right.end &&
         right.begin < left.end;
}

struct StructuredIndex final {
  int i{};
  int j{};
  int k{};
};

StructuredIndex map_cell(Int3 global, runtime::Box3 owned, Int3 global_extent) {
  const Int3 local_extent{owned.end.x - owned.begin.x,
                          owned.end.y - owned.begin.y,
                          owned.end.z - owned.begin.z};
  const auto map_axis = [](int coordinate, int begin, int end, int global_n,
                           int local_n) {
    if (coordinate >= begin && coordinate < end)
      return coordinate - begin;
    if (coordinate == begin - 1 || (begin == 0 && coordinate == global_n - 1))
      return -1;
    if (coordinate == end || (end == global_n && coordinate == 0))
      return local_n;
    throw Error("momentum face cell has no structured field mapping");
  };
  return {map_axis(global.x, owned.begin.x, owned.end.x, global_extent.x,
                   local_extent.x),
          map_axis(global.y, owned.begin.y, owned.end.y, global_extent.y,
                   local_extent.y),
          map_axis(global.z, owned.begin.z, owned.end.z, global_extent.z,
                   local_extent.z)};
}

StructuredIndex shifted(StructuredIndex index, FaceAxis axis, int amount) {
  switch (axis) {
  case FaceAxis::x:
    index.i += amount;
    break;
  case FaceAxis::y:
    index.j += amount;
    break;
  case FaceAxis::z:
    index.k += amount;
    break;
  default:
    throw Error("invalid momentum face axis");
  }
  return index;
}

template <class T>
T &value_at(const runtime::FieldView<T> &view, StructuredIndex index,
            int component) {
  return view(index.i, index.j, index.k, component);
}

struct TopologySignature final {
  struct Cell final {
    mesh::GlobalCellId global{};
    Int3 logical{};
    EntityOwnership ownership{};
    bool operator==(const Cell &other) const noexcept {
      return global == other.global && same(logical, other.logical) &&
             ownership == other.ownership;
    }
  };
  struct Face final {
    GlobalFaceId global{};
    mesh::LogicalFace logical{};
    EntityOwnership ownership{};
    mesh::GlobalCellId owner{};
    std::optional<mesh::GlobalCellId> neighbour;
    std::optional<std::uint32_t> patch;
    std::optional<GlobalFaceId> periodic_pair;
    bool operator==(const Face &other) const noexcept {
      return global == other.global && logical.axis == other.logical.axis &&
             same(logical.coordinate, other.logical.coordinate) &&
             ownership == other.ownership && owner == other.owner &&
             neighbour == other.neighbour && patch == other.patch &&
             periodic_pair == other.periodic_pair;
    }
  };
  struct Patch final {
    std::uint32_t stable_id{};
    std::string name;
    mesh::PatchPairingKind pairing{};
    std::optional<std::uint32_t> paired;
    std::vector<LocalFaceId> faces;
    bool operator==(const Patch &other) const noexcept {
      return stable_id == other.stable_id && name == other.name &&
             pairing == other.pairing && paired == other.paired &&
             faces == other.faces;
    }
  };

  Int3 extent{};
  runtime::Box3 box{};
  std::size_t owned_cells{};
  std::size_t ghost_cells{};
  std::size_t owned_faces{};
  std::size_t ghost_faces{};
  std::vector<Cell> cells;
  std::vector<Face> faces;
  std::array<Patch, 6> patches;
};

TopologySignature make_signature(const mesh::MeshTopology &topology) {
  TopologySignature result{};
  result.extent = topology.global_extent();
  result.box = topology.owned_global_box();
  result.owned_cells = topology.owned_cell_count();
  result.ghost_cells = topology.ghost_cell_count();
  result.owned_faces = topology.owned_face_count();
  result.ghost_faces = topology.ghost_face_count();
  result.cells.reserve(topology.local_cell_count());
  for (LocalCellId cell = 0; cell < topology.local_cell_count(); ++cell) {
    result.cells.push_back({topology.global_cell_id(cell),
                            topology.global_cell(cell),
                            topology.cell_ownership(cell)});
  }
  result.faces.reserve(topology.local_face_count());
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    const auto neighbour = topology.neighbour(face);
    result.faces.push_back(
        {topology.global_face_id(face), topology.logical_face(face),
         topology.face_ownership(face),
         topology.global_cell_id(topology.owner(face)),
         neighbour.has_value()
             ? std::optional<mesh::GlobalCellId>{topology.global_cell_id(
                   *neighbour)}
             : std::nullopt,
         topology.patch_id(face), topology.periodic_pair(face)});
  }
  for (std::size_t index = 0; index < result.patches.size(); ++index) {
    const auto &patch = topology.patch(static_cast<std::uint32_t>(index));
    result.patches[index] = {patch.stable_id(), std::string(patch.name()),
                             patch.pairing_kind(), patch.paired_patch_id(),
                             patch.local_faces()};
  }
  return result;
}

struct FaceStencil final {
  LocalFaceId local{};
  GlobalFaceId global{};
  std::optional<GlobalFaceId> periodic_pair;
  std::optional<LocalFaceId> local_periodic_partner;
  std::optional<std::uint32_t> patch;
  StructuredIndex p{};
  StructuredIndex n{};
  bool has_neighbour{};
  bool reversed{};
  double volume_p{};
  double volume_n{};
  Real3 area{};
  Real3 unit_normal{};
  Real3 displacement{};
  double normal_distance{};
  double weight_p{};
  double weight_n{};
  StructuredIndex physical_owner{};
  Real3 physical_owner_area{};
};

struct OperationGuard final {
  explicit OperationGuard(bool &active) : active_(active) {
    if (active_)
      throw Error("momentum face operation is already active");
    active_ = true;
  }
  ~OperationGuard() noexcept { active_ = false; }
  bool &active_;
};

double component(Real3 value, int index) noexcept {
  return index == 0 ? value.x : index == 1 ? value.y : value.z;
}

double interpolate(double weight_p, double p, double weight_n,
                   double n) noexcept {
  return p == n ? p : weight_p * p + weight_n * n;
}

void validate_cell_layout(const runtime::FieldView<const double> &view,
                          Int3 extent, std::uint32_t components,
                          int minimum_ghost, const char *name) {
  if (!same(view.interior_extent(), extent) ||
      view.components() != components || view.ghost_width() < minimum_ghost) {
    throw Error(std::string("momentum ") + name + " field layout is invalid");
  }
}

void validate_face_layout(const runtime::FaceFieldView<const double> &view,
                          std::size_t faces, const char *name) {
  if (view.face_count() != faces || view.components() != 3U) {
    throw Error(std::string("momentum ") + name + " face layout is invalid");
  }
}

} // namespace

MomentumTimeStencil make_momentum_time_stencil(MomentumTimeOrder order,
                                               double dt_s,
                                               double previous_dt_s) {
  if (!(dt_s > 0.0) || !std::isfinite(dt_s)) {
    throw Error("momentum time step must be positive and finite");
  }
  switch (order) {
  case MomentumTimeOrder::backward_euler:
    return {order, dt_s, 0.0, 1.0, -1.0, 0.0};
  case MomentumTimeOrder::bdf2: {
    if (!(previous_dt_s > 0.0) || !std::isfinite(previous_dt_s)) {
      throw Error("BDF2 previous time step must be positive and finite");
    }
    const double ratio = dt_s / previous_dt_s;
    if (!std::isfinite(ratio) || ratio < 0.5 || ratio > 2.0) {
      throw Error("BDF2 time-step ratio is outside [0.5, 2.0]");
    }
    const double denominator = 1.0 + ratio;
    const double alpha0 = (1.0 + 2.0 * ratio) / denominator;
    const double alpha1 = -(1.0 + ratio);
    const double alpha2 = ratio * ratio / denominator;
    if (!std::isfinite(alpha0) || !std::isfinite(alpha1) ||
        !std::isfinite(alpha2)) {
      throw Error("BDF2 coefficients are non-finite");
    }
    return {order, dt_s, previous_dt_s, alpha0, alpha1, alpha2};
  }
  default:
    throw Error("invalid momentum time order");
  }
}

bool MomentumPredictorReport::all_converged() const noexcept {
  return std::all_of(
      components.begin(), components.end(),
      [](const auto &item) { return solve_success(item.reason); });
}

MomentumPredictor::MomentumPredictor(
    const linear::LinearSolver &solver) noexcept
    : solver_(&solver) {}

MomentumPredictorReport MomentumPredictor::solve(
    const runtime::MpiContext &mpi,
    const std::array<MomentumComponentEquation, 3> &equations,
    const linear::SolveControl &control) const {
  const execution::ExecutionContext *common_context = nullptr;
  std::optional<linear::VectorLayout> common_layout;
  std::array<std::optional<execution::Buffer>, 3> staging;
  std::array<ViewRange, 9> ranges{};
  std::size_t range_count = 0U;
  std::size_t count = 0U;
  std::string local_failure;
  try {
    for (std::size_t component_index = 0; component_index < equations.size();
         ++component_index) {
      const auto &equation = equations[component_index];
      if (equation.linear_operator == nullptr ||
          equation.preconditioner == nullptr) {
        throw Error("momentum equation has a null operator or preconditioner");
      }
      const auto &context = equation.linear_operator->context();
      validate_host_context(context);
      if (common_context == nullptr)
        common_context = &context;
      if (common_context != &context) {
        throw Error("momentum equations do not share one execution context");
      }
      const auto domain = equation.linear_operator->domain_layout();
      const auto range = equation.linear_operator->range_layout();
      if (domain != range)
        throw Error("momentum operator is not square");
      if (!common_layout.has_value())
        common_layout = domain;
      if (*common_layout != domain) {
        throw Error("momentum equation layouts are not identical");
      }
      count = domain.owned_count();
      validate_vector_view(equation.rhs, context, count, false, "RHS");
      validate_vector_view(equation.predictor, context, count, true,
                           "predictor");
      validate_vector_view(equation.actual_diagonal, context, count, true,
                           "diagonal output");
      ranges[range_count++] = view_range(equation.rhs);
      ranges[range_count++] = view_range(equation.predictor);
      ranges[range_count++] = view_range(equation.actual_diagonal);
      if (!equation.linear_operator->has_diagonal()) {
        throw Error("momentum operator lacks diagonal capability");
      }
    }
    for (std::size_t left = 0; left < range_count; ++left) {
      for (std::size_t right = left + 1U; right < range_count; ++right) {
        if (overlaps(ranges[left], ranges[right]) &&
            (ranges[left].writable || ranges[right].writable)) {
          throw Error("momentum input and output views alias");
        }
      }
    }

    count = common_layout->owned_count();
    if (count > std::numeric_limits<std::size_t>::max() / sizeof(double)) {
      throw Error("momentum diagonal staging size overflows");
    }
    for (std::size_t component_index = 0; component_index < equations.size();
         ++component_index) {
      staging[component_index].emplace(
          const_cast<execution::ExecutionContext &>(*common_context),
          count * sizeof(double));
      auto candidate = staging[component_index]->view(0U, count);
      auto event =
          equations[component_index].linear_operator->diagonal(candidate);
      event.wait();
      for (std::size_t index = 0; index < count; ++index) {
        if (!(candidate[index] > 0.0) || !std::isfinite(candidate[index])) {
          throw Error("momentum operator diagonal must be finite and positive");
        }
      }
    }
  } catch (const Error &error) {
    local_failure = error.what();
  } catch (const std::bad_alloc &) {
    local_failure = "momentum diagonal staging allocation failed";
  } catch (const std::length_error &) {
    local_failure = "momentum diagonal staging size is unsupported";
  } catch (const std::exception &error) {
    local_failure =
        std::string("momentum diagonal extraction failed: ") + error.what();
  } catch (...) {
    local_failure = "momentum diagonal extraction failed";
  }

  const runtime::CollectiveStatus preflight = runtime::collective_status(
      mpi, local_failure.empty(), local_failure);
  if (!preflight.ok) {
    throw Error("momentum predictor preflight failed on rank " +
                std::to_string(preflight.failing_rank) + ": " +
                preflight.message);
  }

  for (std::size_t component_index = 0; component_index < equations.size();
       ++component_index) {
    const auto source =
        static_cast<const execution::Buffer &>(*staging[component_index])
            .view(0U, count);
    for (std::size_t index = 0; index < count; ++index) {
      equations[component_index].actual_diagonal[index] = source[index];
    }
  }

  MomentumPredictorReport report{};
  for (std::size_t component_index = 0; component_index < equations.size();
       ++component_index) {
    const auto &equation = equations[component_index];
    report.components[component_index] =
        solver_->solve(*equation.linear_operator, *equation.preconditioner,
                       equation.rhs, equation.predictor, control);
  }
  return report;
}

struct TimeConsistentFaceVelocity::Impl final {
  TopologySignature signature;
  Int3 local_extent{};
  bool needs_ghost{};
  std::vector<FaceStencil> faces;
  mutable std::vector<double> velocity_scratch;
  mutable std::vector<double> flux_scratch;
  mutable bool active{};
};

TimeConsistentFaceVelocity
TimeConsistentFaceVelocity::create(const mesh::MeshTopology &topology,
                                   const mesh::MeshGeometry &geometry) {
  try {
    geometry.require_compatible(topology);
    auto impl = std::make_unique<Impl>();
    impl->signature = make_signature(topology);
    const auto box = topology.owned_global_box();
    impl->local_extent = {box.end.x - box.begin.x, box.end.y - box.begin.y,
                          box.end.z - box.begin.z};
    impl->faces.reserve(topology.local_face_count());
    impl->velocity_scratch.resize(topology.local_face_count() * 3U);
    impl->flux_scratch.resize(topology.local_face_count());

    for (LocalFaceId local = 0; local < topology.local_face_count(); ++local) {
      FaceStencil face{};
      face.local = local;
      face.global = topology.global_face_id(local);
      face.periodic_pair = topology.periodic_pair(local);
      face.patch = topology.patch_id(local);
      const LocalCellId owner = topology.owner(local);
      face.physical_owner =
          map_cell(topology.global_cell(owner), box, topology.global_extent());
      face.physical_owner_area =
          geometry.face_area_vector_m2(local, FaceSide::owner);
      const auto neighbour = topology.neighbour(local);
      if (!neighbour.has_value()) {
        if (!face.patch.has_value() || !finite(face.physical_owner_area)) {
          throw Error("physical momentum face stencil is invalid");
        }
        impl->faces.push_back(face);
        continue;
      }

      StructuredIndex owner_index = face.physical_owner;
      StructuredIndex neighbour_index = map_cell(
          topology.global_cell(*neighbour), box, topology.global_extent());
      if (face.periodic_pair.has_value()) {
        const auto logical = topology.logical_face(local);
        const int plane = logical.axis == FaceAxis::x   ? logical.coordinate.x
                          : logical.axis == FaceAxis::y ? logical.coordinate.y
                                                        : logical.coordinate.z;
        neighbour_index =
            shifted(owner_index, logical.axis, plane == 0 ? -1 : 1);
      }
      face.reversed =
          face.periodic_pair.has_value() && face.global > *face.periodic_pair;
      face.p = face.reversed ? neighbour_index : owner_index;
      face.n = face.reversed ? owner_index : neighbour_index;
      face.has_neighbour = true;
      const LocalCellId p_cell = face.reversed ? *neighbour : owner;
      const LocalCellId n_cell = face.reversed ? owner : *neighbour;
      face.volume_p = geometry.cell_volume_m3(p_cell);
      face.volume_n = geometry.cell_volume_m3(n_cell);
      Real3 p_to_face{};
      if (face.periodic_pair.has_value()) {
        const auto vertices =
            periodic_canonical_vertices(geometry, topology.logical_face(local));
        const Real3 face_center =
            canonical_face_center(vertices, geometry.mapping_kind() ==
                                                mesh::MappingKind::uniform_box);
        face.area = multiply(-1.0, positive_axis_area(vertices));
        const Real3 p_center = geometry.cell_center_m(p_cell);
        const Real3 n_center = geometry.cell_center_m(n_cell);
        face.displacement = subtract(n_center, p_center);
        const Real3 length = geometry.length_m();
        const FaceAxis axis = topology.logical_face(local).axis;
        if (axis == FaceAxis::x)
          face.displacement.x -= length.x;
        if (axis == FaceAxis::y)
          face.displacement.y -= length.y;
        if (axis == FaceAxis::z)
          face.displacement.z -= length.z;
        p_to_face = subtract(face_center, p_center);
      } else {
        face.area = face.physical_owner_area;
        face.displacement = geometry.face_displacement_m(local);
        p_to_face = subtract(geometry.face_center_m(local),
                             geometry.cell_center_m(owner));
      }
      const double area_magnitude = norm(face.area);
      const double a = norm(p_to_face);
      const double b = norm(subtract(face.displacement, p_to_face));
      const double distance_sum = a + b;
      if (!(face.volume_p > 0.0) || !(face.volume_n > 0.0) ||
          !(area_magnitude > 0.0) || !(distance_sum > 0.0) ||
          !std::isfinite(face.volume_p) || !std::isfinite(face.volume_n) ||
          !finite(face.area) || !finite(face.displacement) ||
          !std::isfinite(distance_sum)) {
        throw Error("momentum face geometry is invalid");
      }
      face.unit_normal = multiply(1.0 / area_magnitude, face.area);
      face.normal_distance = dot(face.displacement, face.unit_normal);
      const double area_projection = dot(face.area, face.displacement);
      if (!(face.normal_distance > 0.0) || !(area_projection > 0.0) ||
          !std::isfinite(face.normal_distance) ||
          !std::isfinite(area_projection)) {
        throw Error("momentum face orientation is invalid");
      }
      face.weight_p = b / distance_sum;
      face.weight_n = a / distance_sum;
      if (!(face.weight_p >= 0.0) || !(face.weight_n >= 0.0) ||
          !std::isfinite(face.weight_p) || !std::isfinite(face.weight_n)) {
        throw Error("momentum face interpolation weights are invalid");
      }
      if (face.periodic_pair.has_value() ||
          topology.cell_ownership(owner) == EntityOwnership::ghost ||
          topology.cell_ownership(*neighbour) == EntityOwnership::ghost) {
        impl->needs_ghost = true;
      }
      impl->faces.push_back(face);
    }

    std::vector<std::pair<GlobalFaceId, LocalFaceId>> local_faces_by_global;
    local_faces_by_global.reserve(impl->faces.size());
    for (const FaceStencil &face : impl->faces) {
      local_faces_by_global.emplace_back(face.global, face.local);
    }
    std::sort(local_faces_by_global.begin(), local_faces_by_global.end());
    for (FaceStencil &face : impl->faces) {
      if (!face.periodic_pair.has_value())
        continue;
      const auto found = std::lower_bound(
          local_faces_by_global.begin(), local_faces_by_global.end(),
          std::pair<GlobalFaceId, LocalFaceId>{*face.periodic_pair, 0U},
          [](const auto &left, const auto &right) {
            return left.first < right.first;
          });
      if (found != local_faces_by_global.end() &&
          found->first == *face.periodic_pair) {
        face.local_periodic_partner = found->second;
      }
    }
    return TimeConsistentFaceVelocity(std::move(impl));
  } catch (const Error &) {
    throw;
  } catch (const std::bad_alloc &) {
    throw Error("momentum face stencil allocation failed");
  } catch (const std::length_error &) {
    throw Error("momentum face stencil size is unsupported");
  }
}

TimeConsistentFaceVelocity::TimeConsistentFaceVelocity(
    std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
TimeConsistentFaceVelocity::~TimeConsistentFaceVelocity() noexcept = default;
TimeConsistentFaceVelocity::TimeConsistentFaceVelocity(
    TimeConsistentFaceVelocity &&) noexcept = default;

void TimeConsistentFaceVelocity::assemble_constant_density(
    const boundary::BoundaryRegistry &boundaries, double rho_ref,
    const MomentumTimeStencil &stencil,
    const runtime::FieldView<const double> &predictor_velocity,
    const runtime::FieldView<const double> &mechanical_pressure,
    const runtime::FieldView<const double> &pressure_gradient,
    const runtime::FieldView<const double> &actual_diagonal,
    const MomentumFaceHistory &history,
    const runtime::FaceFieldView<double> &trial_face_velocity,
    const runtime::FieldRegistry &registry, runtime::FieldStorage &storage,
    const runtime::FieldAccessPlan &access_plan, runtime::PhaseId phase,
    runtime::ActorId actor, runtime::FieldId face_mass_flux_field) const {
  if (!impl_)
    throw Error("momentum face object has been moved from");
  OperationGuard operation(impl_->active);
  finite_volume::require_face_mass_flux_field(registry, face_mass_flux_field);
  if (!(rho_ref > 0.0) || !std::isfinite(rho_ref)) {
    throw Error("constant momentum density must be positive and finite");
  }
  const auto expected = make_momentum_time_stencil(stencil.order, stencil.dt_s,
                                                   stencil.previous_dt_s);
  if (stencil.order != expected.order || stencil.dt_s != expected.dt_s ||
      stencil.previous_dt_s != expected.previous_dt_s ||
      stencil.alpha0 != expected.alpha0 || stencil.alpha1 != expected.alpha1 ||
      stencil.alpha2 != expected.alpha2) {
    throw Error("momentum time stencil coefficients are inconsistent");
  }

  const int ghost = impl_->needs_ghost ? 1 : 0;
  validate_cell_layout(predictor_velocity, impl_->local_extent, 3U, ghost,
                       "predictor velocity");
  validate_cell_layout(mechanical_pressure, impl_->local_extent, 1U, ghost,
                       "mechanical pressure");
  validate_cell_layout(pressure_gradient, impl_->local_extent, 3U, ghost,
                       "pressure gradient");
  validate_cell_layout(actual_diagonal, impl_->local_extent, 3U, ghost,
                       "actual diagonal");
  validate_cell_layout(history.velocity_n, impl_->local_extent, 3U, ghost,
                       "history n velocity");
  validate_face_layout(history.face_velocity_n, impl_->faces.size(),
                       "history n");
  if (stencil.order == MomentumTimeOrder::bdf2) {
    if (history.velocity_n_minus_1 == nullptr ||
        history.face_velocity_n_minus_1 == nullptr) {
      throw Error("BDF2 momentum interpolation requires two history layers");
    }
    validate_cell_layout(*history.velocity_n_minus_1, impl_->local_extent, 3U,
                         ghost, "history n-1 velocity");
    validate_face_layout(*history.face_velocity_n_minus_1, impl_->faces.size(),
                         "history n-1");
  }
  if (trial_face_velocity.face_count() != impl_->faces.size() ||
      trial_face_velocity.components() != 3U ||
      storage.layout_set().face_count != impl_->faces.size()) {
    throw Error("momentum face output layout is invalid");
  }
  if (!impl_->faces.empty()) {
    static_cast<void>(trial_face_velocity(0U, 0));
  }

  for (int k = 0; k < impl_->local_extent.z; ++k) {
    for (int j = 0; j < impl_->local_extent.y; ++j) {
      for (int i = 0; i < impl_->local_extent.x; ++i) {
        const StructuredIndex cell{i, j, k};
        const double pressure = value_at(mechanical_pressure, cell, 0);
        if (!std::isfinite(pressure)) {
          throw Error("momentum pressure input is non-finite");
        }
        for (int component_index = 0; component_index < 3; ++component_index) {
          const double predictor =
              value_at(predictor_velocity, cell, component_index);
          const double gradient =
              value_at(pressure_gradient, cell, component_index);
          const double diagonal =
              value_at(actual_diagonal, cell, component_index);
          const double current =
              value_at(history.velocity_n, cell, component_index);
          if (!std::isfinite(predictor) || !std::isfinite(gradient) ||
              !std::isfinite(current) || !(diagonal > 0.0) ||
              !std::isfinite(diagonal)) {
            throw Error("momentum owned cell input is invalid");
          }
          if (stencil.order == MomentumTimeOrder::bdf2 &&
              !std::isfinite(value_at(*history.velocity_n_minus_1, cell,
                                      component_index))) {
            throw Error("momentum BDF2 owned history is non-finite");
          }
        }
      }
    }
  }
  for (const FaceStencil &face : impl_->faces) {
    for (int component_index = 0; component_index < 3; ++component_index) {
      if (!std::isfinite(
              history.face_velocity_n(face.local, component_index))) {
        throw Error("momentum face history is non-finite");
      }
      if (stencil.order == MomentumTimeOrder::bdf2 &&
          !std::isfinite((*history.face_velocity_n_minus_1)(face.local,
                                                            component_index))) {
        throw Error("momentum BDF2 face history is non-finite");
      }
    }
  }

  for (const FaceStencil &face : impl_->faces) {
    if (face.local_periodic_partner.has_value() &&
        face.periodic_pair.has_value() && face.global > *face.periodic_pair) {
      continue;
    }
    if (!face.has_neighbour) {
      const Real3 interior{
          value_at(predictor_velocity, face.physical_owner, 0),
          value_at(predictor_velocity, face.physical_owner, 1),
          value_at(predictor_velocity, face.physical_owner, 2)};
      if (!finite(interior) || !face.patch.has_value()) {
        throw Error("physical momentum face state is invalid");
      }
      const auto values = boundaries.evaluate_velocity(
          *face.patch, interior, face.physical_owner_area);
      if (!finite(values.face)) {
        throw Error("physical momentum face velocity is non-finite");
      }
      const auto rule = boundaries.patch(*face.patch).mass_flux_rule();
      double flux = 0.0;
      if (rule == boundary::MassFluxRule::identically_zero) {
        flux = 0.0;
      } else if (rule == boundary::MassFluxRule::prescribed_inlet_state ||
                 rule == boundary::MassFluxRule::outflow_only) {
        flux = rho_ref * dot(values.face, face.physical_owner_area);
      } else {
        throw Error("physical momentum face has an invalid mass-flux rule");
      }
      if (!std::isfinite(flux)) {
        throw Error("physical momentum mass flux is non-finite");
      }
      impl_->velocity_scratch[face.local * 3U] = values.face.x;
      impl_->velocity_scratch[face.local * 3U + 1U] = values.face.y;
      impl_->velocity_scratch[face.local * 3U + 2U] = values.face.z;
      impl_->flux_scratch[face.local] = flux;
      continue;
    }
    if (face.periodic_pair.has_value()) {
      if (!face.patch.has_value() ||
          boundaries.patch(*face.patch).mass_flux_rule() !=
              boundary::MassFluxRule::periodic_pair) {
        throw Error("periodic momentum face boundary identity is invalid");
      }
    }

    const double pi_p = value_at(mechanical_pressure, face.p, 0);
    const double pi_n = value_at(mechanical_pressure, face.n, 0);
    if (!std::isfinite(pi_p) || !std::isfinite(pi_n)) {
      throw Error("momentum pressure input is non-finite");
    }
    Real3 gradient{};
    for (int direction = 0; direction < 3; ++direction) {
      const double p = value_at(pressure_gradient, face.p, direction);
      const double n = value_at(pressure_gradient, face.n, direction);
      const double value = interpolate(face.weight_p, p, face.weight_n, n);
      if (!std::isfinite(p) || !std::isfinite(n) || !std::isfinite(value)) {
        throw Error("momentum pressure gradient input is non-finite");
      }
      if (direction == 0)
        gradient.x = value;
      if (direction == 1)
        gradient.y = value;
      if (direction == 2)
        gradient.z = value;
    }
    const double pressure_defect =
        (pi_n - pi_p - dot(gradient, face.displacement)) / face.normal_distance;
    if (!std::isfinite(pressure_defect)) {
      throw Error("momentum pressure interpolation defect is non-finite");
    }

    Real3 trial{};
    for (int velocity_component = 0; velocity_component < 3;
         ++velocity_component) {
      const double diagonal_p =
          value_at(actual_diagonal, face.p, velocity_component);
      const double diagonal_n =
          value_at(actual_diagonal, face.n, velocity_component);
      if (!(diagonal_p > 0.0) || !(diagonal_n > 0.0) ||
          !std::isfinite(diagonal_p) || !std::isfinite(diagonal_n)) {
        throw Error("momentum actual diagonal is not finite and positive");
      }
      const double lambda = face.weight_p * (diagonal_p / face.volume_p) +
                            face.weight_n * (diagonal_n / face.volume_n);
      if (!(lambda > 0.0) || !std::isfinite(lambda)) {
        throw Error("momentum face diagonal density is invalid");
      }
      const double mobility = 1.0 / lambda;
      if (!(mobility > 0.0) || !std::isfinite(mobility)) {
        throw Error("momentum face mobility is invalid");
      }
      const double predictor_p =
          value_at(predictor_velocity, face.p, velocity_component);
      const double predictor_n =
          value_at(predictor_velocity, face.n, velocity_component);
      const double predictor_face =
          interpolate(face.weight_p, predictor_p, face.weight_n, predictor_n);
      const double velocity_n_p =
          value_at(history.velocity_n, face.p, velocity_component);
      const double velocity_n_n =
          value_at(history.velocity_n, face.n, velocity_component);
      const double cell_history_n =
          interpolate(face.weight_p, velocity_n_p, face.weight_n, velocity_n_n);
      const double face_history_n =
          history.face_velocity_n(face.local, velocity_component);
      if (!std::isfinite(predictor_p) || !std::isfinite(predictor_n) ||
          !std::isfinite(predictor_face) || !std::isfinite(velocity_n_p) ||
          !std::isfinite(velocity_n_n) || !std::isfinite(cell_history_n) ||
          !std::isfinite(face_history_n)) {
        throw Error("momentum velocity input is non-finite");
      }
      const double discrepancy_n = face_history_n - cell_history_n;
      double discrepancy_n_minus_1 = 0.0;
      if (stencil.order == MomentumTimeOrder::bdf2) {
        const double old_p =
            value_at(*history.velocity_n_minus_1, face.p, velocity_component);
        const double old_n =
            value_at(*history.velocity_n_minus_1, face.n, velocity_component);
        const double old_cell =
            interpolate(face.weight_p, old_p, face.weight_n, old_n);
        const double old_face =
            (*history.face_velocity_n_minus_1)(face.local, velocity_component);
        if (!std::isfinite(old_p) || !std::isfinite(old_n) ||
            !std::isfinite(old_cell) || !std::isfinite(old_face)) {
          throw Error("momentum BDF2 history input is non-finite");
        }
        discrepancy_n_minus_1 = old_face - old_cell;
      }
      const double pressure_correction =
          -mobility * pressure_defect *
          component(face.unit_normal, velocity_component);
      const double time_correction =
          mobility * (rho_ref / stencil.dt_s) *
          ((-stencil.alpha1) * discrepancy_n +
           (-stencil.alpha2) * discrepancy_n_minus_1);
      const double value =
          predictor_face + pressure_correction + time_correction;
      if (!std::isfinite(value)) {
        throw Error("time-consistent momentum face velocity is non-finite");
      }
      if (velocity_component == 0)
        trial.x = value;
      if (velocity_component == 1)
        trial.y = value;
      if (velocity_component == 2)
        trial.z = value;
    }
    double canonical_flux = rho_ref * dot(trial, face.area);
    if (!std::isfinite(canonical_flux)) {
      throw Error("time-consistent momentum mass flux is non-finite");
    }
    if (canonical_flux == 0.0)
      canonical_flux = 0.0;
    impl_->velocity_scratch[face.local * 3U] = trial.x;
    impl_->velocity_scratch[face.local * 3U + 1U] = trial.y;
    impl_->velocity_scratch[face.local * 3U + 2U] = trial.z;
    impl_->flux_scratch[face.local] =
        face.reversed ? -canonical_flux : canonical_flux;
  }

  for (const FaceStencil &face : impl_->faces) {
    if (!face.local_periodic_partner.has_value() ||
        !face.periodic_pair.has_value() || face.global < *face.periodic_pair) {
      continue;
    }
    const LocalFaceId representative = *face.local_periodic_partner;
    for (std::size_t component_index = 0; component_index < 3U;
         ++component_index) {
      impl_->velocity_scratch[face.local * 3U + component_index] =
          impl_->velocity_scratch[representative * 3U + component_index];
    }
    impl_->flux_scratch[face.local] = -impl_->flux_scratch[representative];
  }

  auto flux_writer = storage.acquire_face_write<double>(
      access_plan, phase, actor, face_mass_flux_field);
  if (flux_writer.face_count() != impl_->faces.size() ||
      flux_writer.components() != 1U) {
    throw Error("momentum mass-flux writer layout is invalid");
  }
  if (!impl_->faces.empty())
    static_cast<void>(flux_writer(0U, 0));
  for (const FaceStencil &face : impl_->faces) {
    for (int component_index = 0; component_index < 3; ++component_index) {
      trial_face_velocity(face.local, component_index) =
          impl_->velocity_scratch[face.local * 3U +
                                  static_cast<std::size_t>(component_index)];
    }
    flux_writer(face.local, 0) = impl_->flux_scratch[face.local];
  }
}

} // namespace hundun::flow

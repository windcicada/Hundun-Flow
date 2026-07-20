// SPDX-License-Identifier: Apache-2.0

#include "hundun/finite_volume/cell_centered_fvm.hpp"

#include "cell_centered_fvm_test_seam.hpp"
#include "hundun/runtime/error.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace hundun::finite_volume {
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

constexpr std::size_t kInvalidFace = std::numeric_limits<std::size_t>::max();

bool finite(Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

Real3 add(Real3 lhs, Real3 rhs) noexcept {
  return Real3{lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

Real3 subtract(Real3 lhs, Real3 rhs) noexcept {
  return Real3{lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

Real3 multiply(double scale, Real3 value) noexcept {
  return Real3{scale * value.x, scale * value.y, scale * value.z};
}

double dot(Real3 lhs, Real3 rhs) noexcept {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

double safe_norm(Real3 value) noexcept {
  const double scale =
      std::max({std::abs(value.x), std::abs(value.y), std::abs(value.z)});
  if (scale == 0.0) {
    return 0.0;
  }
  if (!std::isfinite(scale)) {
    return std::numeric_limits<double>::infinity();
  }
  const double x = value.x / scale;
  const double y = value.y / scale;
  const double z = value.z / scale;
  return scale * std::sqrt(x * x + y * y + z * z);
}

struct StructuredIndex final {
  int i{};
  int j{};
  int k{};
};

bool same(Int3 lhs, Int3 rhs) noexcept {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

bool same(runtime::Box3 lhs, runtime::Box3 rhs) noexcept {
  return same(lhs.begin, rhs.begin) && same(lhs.end, rhs.end);
}

std::size_t owned_linear(StructuredIndex index, Int3 extent) {
  if (index.i < 0 || index.i >= extent.x || index.j < 0 ||
      index.j >= extent.y || index.k < 0 || index.k >= extent.z) {
    throw Error("finite-volume owned cell mapping is invalid");
  }
  return (static_cast<std::size_t>(index.k) *
              static_cast<std::size_t>(extent.y) +
          static_cast<std::size_t>(index.j)) *
             static_cast<std::size_t>(extent.x) +
         static_cast<std::size_t>(index.i);
}

StructuredIndex map_cell(Int3 global, runtime::Box3 owned, Int3 global_extent) {
  const Int3 local_extent{owned.end.x - owned.begin.x,
                          owned.end.y - owned.begin.y,
                          owned.end.z - owned.begin.z};
  const auto map_axis = [](int coordinate, int begin, int end, int global_n,
                           int local_n) {
    if (coordinate >= begin && coordinate < end) {
      return coordinate - begin;
    }
    if (coordinate == begin - 1 || (begin == 0 && coordinate == global_n - 1)) {
      return -1;
    }
    if (coordinate == end || (end == global_n && coordinate == 0)) {
      return local_n;
    }
    throw Error("compact topology cell has no adjacent structured mapping");
  };
  return StructuredIndex{map_axis(global.x, owned.begin.x, owned.end.x,
                                  global_extent.x, local_extent.x),
                         map_axis(global.y, owned.begin.y, owned.end.y,
                                  global_extent.y, local_extent.y),
                         map_axis(global.z, owned.begin.z, owned.end.z,
                                  global_extent.z, local_extent.z)};
}

int coordinate(StructuredIndex index, FaceAxis axis) noexcept {
  switch (axis) {
  case FaceAxis::x:
    return index.i;
  case FaceAxis::y:
    return index.j;
  case FaceAxis::z:
    return index.k;
  }
  return 0;
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
    throw Error("invalid finite-volume face axis");
  }
  return index;
}

template <class T>
T &value_at(const runtime::FieldView<T> &view, StructuredIndex index,
            int component) {
  return view(index.i, index.j, index.k, component);
}

struct TopologySignature final {
  Int3 global_extent{};
  runtime::Box3 owned_box{};
  struct Cell final {
    LocalCellId local_id{};
    mesh::GlobalCellId global_id{};
    Int3 global_cell{};
    EntityOwnership ownership{};

    bool operator==(const Cell &other) const noexcept {
      return local_id == other.local_id && global_id == other.global_id &&
             same(global_cell, other.global_cell) &&
             ownership == other.ownership;
    }
  };
  struct Face final {
    LocalFaceId local_id{};
    GlobalFaceId global_id{};
    mesh::LogicalFace logical{};
    EntityOwnership ownership{};
    LocalCellId owner_local_id{};
    mesh::GlobalCellId owner_global_id{};
    EntityOwnership owner_ownership{};
    std::optional<LocalCellId> neighbour_local_id;
    std::optional<mesh::GlobalCellId> neighbour_global_id;
    std::optional<EntityOwnership> neighbour_ownership;
    std::optional<std::uint32_t> patch_id;
    std::optional<GlobalFaceId> periodic_pair;

    bool operator==(const Face &other) const noexcept {
      return local_id == other.local_id && global_id == other.global_id &&
             logical.axis == other.logical.axis &&
             same(logical.coordinate, other.logical.coordinate) &&
             ownership == other.ownership &&
             owner_local_id == other.owner_local_id &&
             owner_global_id == other.owner_global_id &&
             owner_ownership == other.owner_ownership &&
             neighbour_local_id == other.neighbour_local_id &&
             neighbour_global_id == other.neighbour_global_id &&
             neighbour_ownership == other.neighbour_ownership &&
             patch_id == other.patch_id && periodic_pair == other.periodic_pair;
    }
  };
  struct Patch final {
    std::uint32_t stable_id{};
    std::string name;
    mesh::PatchPairingKind pairing_kind{};
    std::optional<std::uint32_t> paired_patch_id;
    std::vector<LocalFaceId> local_faces;

    bool operator==(const Patch &other) const noexcept {
      return stable_id == other.stable_id && name == other.name &&
             pairing_kind == other.pairing_kind &&
             paired_patch_id == other.paired_patch_id &&
             local_faces == other.local_faces;
    }
  };

  std::size_t owned_cell_count{};
  std::size_t ghost_cell_count{};
  std::size_t owned_face_count{};
  std::size_t ghost_face_count{};
  std::vector<Cell> cells;
  std::vector<Face> faces;
  std::array<Patch, 6> patches;
};

thread_local std::optional<test::TopologySignatureMutationForTest>
    next_topology_signature_mutation;

template <class Integer> Integer different(Integer value) noexcept {
  return value == std::numeric_limits<Integer>::max() ? value - 1 : value + 1;
}

EntityOwnership other_ownership(EntityOwnership ownership) noexcept {
  return ownership == EntityOwnership::owned ? EntityOwnership::ghost
                                             : EntityOwnership::owned;
}

void mutate_signature(TopologySignature &signature,
                      test::TopologySignatureMutationForTest mutation) {
  using Mutation = test::TopologySignatureMutationForTest;
  if (signature.cells.empty() || signature.faces.empty()) {
    throw Error("topology signature mutation requires non-empty entities");
  }
  auto face_with_neighbour = std::find_if(
      signature.faces.begin(), signature.faces.end(),
      [](const auto &face) { return face.neighbour_local_id.has_value(); });
  auto face_with_patch =
      std::find_if(signature.faces.begin(), signature.faces.end(),
                   [](const auto &face) { return face.patch_id.has_value(); });
  auto face_with_pair = std::find_if(
      signature.faces.begin(), signature.faces.end(),
      [](const auto &face) { return face.periodic_pair.has_value(); });
  switch (mutation) {
  case Mutation::cell_global_id:
    signature.cells.front().global_id =
        different(signature.cells.front().global_id);
    return;
  case Mutation::cell_ownership:
    signature.cells.front().ownership =
        other_ownership(signature.cells.front().ownership);
    return;
  case Mutation::face_ownership:
    signature.faces.front().ownership =
        other_ownership(signature.faces.front().ownership);
    return;
  case Mutation::face_owner_local_id:
    signature.faces.front().owner_local_id =
        different(signature.faces.front().owner_local_id);
    return;
  case Mutation::face_owner_global_id:
    signature.faces.front().owner_global_id =
        different(signature.faces.front().owner_global_id);
    return;
  case Mutation::face_owner_ownership:
    signature.faces.front().owner_ownership =
        other_ownership(signature.faces.front().owner_ownership);
    return;
  case Mutation::face_neighbour_presence:
    if (face_with_neighbour == signature.faces.end())
      throw Error("topology signature has no neighbour to mutate");
    face_with_neighbour->neighbour_local_id.reset();
    face_with_neighbour->neighbour_global_id.reset();
    face_with_neighbour->neighbour_ownership.reset();
    return;
  case Mutation::face_neighbour_local_id:
    if (face_with_neighbour == signature.faces.end())
      throw Error("topology signature has no neighbour to mutate");
    *face_with_neighbour->neighbour_local_id =
        different(*face_with_neighbour->neighbour_local_id);
    return;
  case Mutation::face_neighbour_global_id:
    if (face_with_neighbour == signature.faces.end())
      throw Error("topology signature has no neighbour to mutate");
    *face_with_neighbour->neighbour_global_id =
        different(*face_with_neighbour->neighbour_global_id);
    return;
  case Mutation::face_neighbour_ownership:
    if (face_with_neighbour == signature.faces.end())
      throw Error("topology signature has no neighbour to mutate");
    *face_with_neighbour->neighbour_ownership =
        other_ownership(*face_with_neighbour->neighbour_ownership);
    return;
  case Mutation::logical_face:
    signature.faces.front().logical.coordinate.x =
        different(signature.faces.front().logical.coordinate.x);
    return;
  case Mutation::face_patch_membership:
    if (face_with_patch == signature.faces.end())
      throw Error("topology signature has no patch face to mutate");
    face_with_patch->patch_id.reset();
    return;
  case Mutation::periodic_pair:
    if (face_with_pair == signature.faces.end())
      throw Error("topology signature has no periodic pair to mutate");
    face_with_pair->periodic_pair.reset();
    return;
  case Mutation::patch_stable_id:
    signature.patches.front().stable_id =
        different(signature.patches.front().stable_id);
    return;
  case Mutation::patch_name:
    signature.patches.front().name += "_mutated";
    return;
  case Mutation::patch_pairing_kind:
    signature.patches.front().pairing_kind =
        signature.patches.front().pairing_kind == mesh::PatchPairingKind::none
            ? mesh::PatchPairingKind::periodic
            : mesh::PatchPairingKind::none;
    return;
  case Mutation::patch_paired_id:
    if (signature.patches.front().paired_patch_id.has_value()) {
      signature.patches.front().paired_patch_id.reset();
    } else {
      signature.patches.front().paired_patch_id = 1U;
    }
    return;
  case Mutation::patch_exact_membership: {
    auto patch = std::find_if(
        signature.patches.begin(), signature.patches.end(),
        [](const auto &candidate) { return !candidate.local_faces.empty(); });
    if (patch == signature.patches.end())
      throw Error("topology signature has no patch membership to mutate");
    patch->local_faces.pop_back();
    return;
  }
  }
  throw Error("invalid topology signature mutation");
}

TopologySignature make_signature(const mesh::MeshTopology &topology) {
  TopologySignature signature{};
  signature.global_extent = topology.global_extent();
  signature.owned_box = topology.owned_global_box();
  signature.owned_cell_count = topology.owned_cell_count();
  signature.ghost_cell_count = topology.ghost_cell_count();
  signature.owned_face_count = topology.owned_face_count();
  signature.ghost_face_count = topology.ghost_face_count();
  signature.cells.reserve(topology.local_cell_count());
  for (LocalCellId cell = 0; cell < topology.local_cell_count(); ++cell) {
    signature.cells.push_back(TopologySignature::Cell{
        cell, topology.global_cell_id(cell), topology.global_cell(cell),
        topology.cell_ownership(cell)});
  }
  signature.faces.reserve(topology.local_face_count());
  for (LocalFaceId face = 0; face < topology.local_face_count(); ++face) {
    const LocalCellId owner = topology.owner(face);
    const auto neighbour = topology.neighbour(face);
    signature.faces.push_back(TopologySignature::Face{
        face, topology.global_face_id(face), topology.logical_face(face),
        topology.face_ownership(face), owner, topology.global_cell_id(owner),
        topology.cell_ownership(owner), neighbour,
        neighbour.has_value()
            ? std::optional<mesh::GlobalCellId>{topology.global_cell_id(
                  *neighbour)}
            : std::nullopt,
        neighbour.has_value()
            ? std::optional<EntityOwnership>{topology.cell_ownership(
                  *neighbour)}
            : std::nullopt,
        topology.patch_id(face), topology.periodic_pair(face)});
  }
  for (std::size_t index = 0; index < signature.patches.size(); ++index) {
    const auto &patch = topology.patch(static_cast<std::uint32_t>(index));
    signature.patches[index] = TopologySignature::Patch{
        patch.stable_id(), std::string(patch.name()), patch.pairing_kind(),
        patch.paired_patch_id(), patch.local_faces()};
  }
  if (next_topology_signature_mutation.has_value()) {
    const auto mutation = *next_topology_signature_mutation;
    next_topology_signature_mutation.reset();
    mutate_signature(signature, mutation);
  }
  return signature;
}

bool same_signature(const TopologySignature &lhs,
                    const TopologySignature &rhs) noexcept {
  return same(lhs.global_extent, rhs.global_extent) &&
         same(lhs.owned_box, rhs.owned_box) &&
         lhs.owned_cell_count == rhs.owned_cell_count &&
         lhs.ghost_cell_count == rhs.ghost_cell_count &&
         lhs.owned_face_count == rhs.owned_face_count &&
         lhs.ghost_face_count == rhs.ghost_face_count &&
         lhs.cells == rhs.cells && lhs.faces == rhs.faces &&
         lhs.patches == rhs.patches;
}

struct FaceStencil final {
  LocalFaceId local_face{};
  GlobalFaceId global_face{};
  std::optional<GlobalFaceId> periodic_pair;
  FaceAxis axis{FaceAxis::x};
  StructuredIndex owner{};
  std::optional<StructuredIndex> neighbour;
  std::optional<std::size_t> owner_owned;
  std::optional<std::size_t> neighbour_owned;
  std::optional<std::uint32_t> patch;
  Real3 owner_area{};
  Real3 face_center{};
  Real3 owner_center{};
  std::optional<Real3> neighbour_center;
  Real3 displacement{};
  double skewness{};
  double non_orthogonality_degrees{};
};

struct CellFaceRef final {
  std::size_t face{kInvalidFace};
  bool as_owner{true};
};

bool is_periodic(const FaceStencil &face) noexcept {
  return face.periodic_pair.has_value();
}

bool is_physical_nonperiodic(const FaceStencil &face) noexcept {
  return face.patch.has_value() && !is_periodic(face);
}

void validate_quantity(FiniteVolumeQuantity quantity,
                       const boundary::BoundaryRegistry &boundaries,
                       bool allow_velocity) {
  switch (quantity.kind) {
  case FiniteVolumeQuantityKind::density:
  case FiniteVolumeQuantityKind::enthalpy:
    if (quantity.scalar_index != 0U) {
      throw Error("non-scalar finite-volume quantity has a scalar index");
    }
    return;
  case FiniteVolumeQuantityKind::scalar:
    if (quantity.scalar_index >= boundaries.scalar_count()) {
      throw Error("finite-volume scalar index is out of bounds");
    }
    return;
  case FiniteVolumeQuantityKind::velocity:
    if (!allow_velocity || quantity.scalar_index != 0U) {
      throw Error("finite-volume quantity is invalid for this operation");
    }
    return;
  default:
    throw Error("invalid finite-volume quantity kind");
  }
}

boundary::ScalarBoundaryValues
evaluate_scalar_boundary(const boundary::BoundaryRegistry &boundaries,
                         FiniteVolumeQuantity quantity, std::uint32_t patch,
                         double interior) {
  switch (quantity.kind) {
  case FiniteVolumeQuantityKind::density:
    return boundaries.evaluate_density(patch, interior);
  case FiniteVolumeQuantityKind::enthalpy:
    return boundaries.evaluate_enthalpy(patch, interior);
  case FiniteVolumeQuantityKind::scalar:
    return boundaries.evaluate_scalar(patch, quantity.scalar_index, interior);
  default:
    throw Error("scalar boundary evaluation received velocity");
  }
}

double boundary_face_value(const boundary::BoundaryRegistry &boundaries,
                           FiniteVolumeQuantity quantity,
                           const FaceStencil &face, double interior) {
  if (!face.patch.has_value()) {
    throw Error("finite-volume physical face has no boundary patch");
  }
  return evaluate_scalar_boundary(boundaries, quantity, *face.patch, interior)
      .face;
}

double boundary_exterior_value(const boundary::BoundaryRegistry &boundaries,
                               FiniteVolumeQuantity quantity,
                               const FaceStencil &face, Real3 velocity,
                               double scalar, int component) {
  if (!face.patch.has_value()) {
    throw Error("finite-volume physical face has no boundary patch");
  }
  if (quantity.kind == FiniteVolumeQuantityKind::velocity) {
    const auto values =
        boundaries.evaluate_velocity(*face.patch, velocity, face.owner_area);
    return component == 0   ? values.exterior.x
           : component == 1 ? values.exterior.y
                            : values.exterior.z;
  }
  return evaluate_scalar_boundary(boundaries, quantity, *face.patch, scalar)
      .exterior;
}

bool descriptor_matches(const runtime::FieldDescriptor &descriptor) noexcept {
  return descriptor.name == "face_mass_flux" && descriptor.unit == "kg/s" &&
         descriptor.owner == "flow" &&
         descriptor.space == runtime::FunctionSpace::face_value &&
         descriptor.scalar_type == runtime::ScalarType::float64 &&
         descriptor.components == 1U && descriptor.ghost_width == 0 &&
         descriptor.conservative &&
         descriptor.restart == runtime::RestartPolicy::persistent &&
         descriptor.output == runtime::OutputPolicy::never;
}

struct OperationGuard final {
  explicit OperationGuard(bool &active) : active_(active) {
    if (active_) {
      throw Error("finite-volume operator is already active");
    }
    active_ = true;
  }
  ~OperationGuard() noexcept { active_ = false; }
  bool &active_;
};

struct MetricOverride final {
  bool armed{};
  GlobalFaceId global_face{};
  double skewness{};
  double non_orthogonality_degrees{};
};

thread_local MetricOverride next_metric_override{};
thread_local bool force_singular_least_squares{};
thread_local std::optional<test::FaceMassFluxConstructionFailureForTest>
    next_face_mass_flux_construction_failure;

MetricOverride consume_metric_override() noexcept {
  const MetricOverride result = next_metric_override;
  next_metric_override = MetricOverride{};
  return result;
}

void inject_face_mass_flux_construction_failure() {
  if (!next_face_mass_flux_construction_failure.has_value()) {
    return;
  }
  const auto failure = *next_face_mass_flux_construction_failure;
  next_face_mass_flux_construction_failure.reset();
  switch (failure) {
  case test::FaceMassFluxConstructionFailureForTest::bad_alloc:
    throw std::bad_alloc{};
  case test::FaceMassFluxConstructionFailureForTest::length_error:
    throw std::length_error("injected face mass flux construction failure");
  }
}

struct CanonicalFace final {
  StructuredIndex p{};
  StructuredIndex n{};
  Real3 area{};
  int step{};
  bool reversed{};
};

CanonicalFace canonical_face(const FaceStencil &face) {
  if (!face.neighbour.has_value()) {
    throw Error("canonical finite-volume face has no neighbour");
  }
  const bool reversed =
      face.periodic_pair.has_value() && face.global_face > *face.periodic_pair;
  CanonicalFace result{};
  result.p = reversed ? *face.neighbour : face.owner;
  result.n = reversed ? face.owner : *face.neighbour;
  result.area = reversed ? multiply(-1.0, face.owner_area) : face.owner_area;
  result.reversed = reversed;
  result.step =
      coordinate(result.n, face.axis) - coordinate(result.p, face.axis);
  if (result.step != -1 && result.step != 1) {
    throw Error("finite-volume face is not adjacent in structured space");
  }
  return result;
}

double reconstruct_mc(const runtime::FieldView<const double> &values,
                      const FaceStencil &face, const CanonicalFace &canonical,
                      int component, double canonical_mass_flux) {
  const StructuredIndex pm1 = shifted(canonical.p, face.axis, -canonical.step);
  const StructuredIndex np1 = shifted(canonical.n, face.axis, canonical.step);
  const double q_pm1 = value_at(values, pm1, component);
  const double q_p = value_at(values, canonical.p, component);
  const double q_n = value_at(values, canonical.n, component);
  const double q_np1 = value_at(values, np1, component);
  if (!std::isfinite(q_pm1) || !std::isfinite(q_p) || !std::isfinite(q_n) ||
      !std::isfinite(q_np1)) {
    throw Error("MUSCL stencil contains a non-finite value");
  }
  const double slope_p = monotonized_central(q_p - q_pm1, q_n - q_p);
  const double slope_n = monotonized_central(q_n - q_p, q_np1 - q_n);
  const double left = q_p + 0.5 * slope_p;
  const double right = q_n - 0.5 * slope_n;
  const double result = canonical_mass_flux >= 0.0 ? left : right;
  if (!std::isfinite(result)) {
    throw Error("MUSCL reconstruction produced a non-finite value");
  }
  return result;
}

} // namespace

struct FaceMassFlux::State final {
  State(runtime::FieldId field_value,
        runtime::FaceFieldView<const double> view_value,
        TopologySignature signature_value)
      : field(field_value), view(std::move(view_value)),
        signature(std::move(signature_value)) {}

  runtime::FieldId field{};
  runtime::FaceFieldView<const double> view;
  TopologySignature signature;
};

struct CellCenteredFvmOperators::Impl final {
  TopologySignature signature;
  Int3 local_extent{};
  bool needs_remote_or_periodic{};
  std::vector<FaceStencil> faces;
  std::vector<std::array<CellFaceRef, 6>> cell_faces;
  std::vector<double> inverse_cell_volumes;
  mutable std::vector<double> scratch;
  mutable std::vector<double> secondary;
  mutable std::vector<double> face_scratch;
  mutable bool active{};
};

namespace test {

void mutate_next_topology_signature(TopologySignatureMutationForTest mutation) {
  if (next_topology_signature_mutation.has_value()) {
    throw Error("finite-volume topology signature mutation is already armed");
  }
  next_topology_signature_mutation = mutation;
}

void override_next_face_metrics(mesh::GlobalFaceId global_face, double skewness,
                                double non_orthogonality_degrees) {
  if (next_metric_override.armed) {
    throw Error("finite-volume metric override is already armed");
  }
  if (!std::isfinite(skewness) || !std::isfinite(non_orthogonality_degrees)) {
    throw Error("finite-volume metric override is non-finite");
  }
  next_metric_override =
      MetricOverride{true, global_face, skewness, non_orthogonality_degrees};
}

void force_next_least_squares_singular() {
  if (force_singular_least_squares) {
    throw Error("finite-volume singular least-squares seam is already armed");
  }
  force_singular_least_squares = true;
}

void fail_next_face_mass_flux_construction(
    FaceMassFluxConstructionFailureForTest failure) {
  if (next_face_mass_flux_construction_failure.has_value()) {
    throw Error("face mass flux construction failure is already armed");
  }
  next_face_mass_flux_construction_failure = failure;
}

} // namespace test

FiniteVolumeQuantity FiniteVolumeQuantity::density() noexcept {
  return FiniteVolumeQuantity{FiniteVolumeQuantityKind::density, 0U};
}

FiniteVolumeQuantity FiniteVolumeQuantity::enthalpy() noexcept {
  return FiniteVolumeQuantity{FiniteVolumeQuantityKind::enthalpy, 0U};
}

FiniteVolumeQuantity FiniteVolumeQuantity::scalar(std::size_t index) noexcept {
  return FiniteVolumeQuantity{FiniteVolumeQuantityKind::scalar, index};
}

FiniteVolumeQuantity FiniteVolumeQuantity::velocity() noexcept {
  return FiniteVolumeQuantity{FiniteVolumeQuantityKind::velocity, 0U};
}

double monotonized_central(double left, double right) noexcept {
  if (!std::isfinite(left) || !std::isfinite(right) || left == 0.0 ||
      right == 0.0 || std::signbit(left) != std::signbit(right)) {
    return 0.0;
  }
  const double a = std::abs(left);
  const double b = std::abs(right);
  constexpr double maximum = std::numeric_limits<double>::max();
  const double twice_a = a > maximum * 0.5 ? maximum : 2.0 * a;
  const double twice_b = b > maximum * 0.5 ? maximum : 2.0 * b;
  const double mean = 0.5 * a + 0.5 * b;
  const double magnitude = std::min({twice_a, mean, twice_b});
  return std::signbit(left) ? -magnitude : magnitude;
}

runtime::FieldDescriptor face_mass_flux_descriptor() {
  return runtime::FieldDescriptor{"face_mass_flux",
                                  "kg/s",
                                  "flow",
                                  runtime::FunctionSpace::face_value,
                                  runtime::ScalarType::float64,
                                  1U,
                                  0,
                                  true,
                                  runtime::RestartPolicy::persistent,
                                  runtime::OutputPolicy::never};
}

runtime::FieldId declare_face_mass_flux(runtime::FieldRegistry &registry) {
  return registry.declare_field(face_mass_flux_descriptor());
}

void require_face_mass_flux_field(const runtime::FieldRegistry &registry,
                                  runtime::FieldId field) {
  if (!registry.frozen()) {
    throw Error("face mass flux requires a frozen field registry");
  }
  if (!descriptor_matches(registry.descriptor(field))) {
    throw Error("face mass flux descriptor does not match the frozen contract");
  }
}

FaceMassFlux FaceMassFlux::acquire(const runtime::FieldRegistry &registry,
                                   const runtime::FieldStorage &storage,
                                   const runtime::FieldAccessPlan &access_plan,
                                   runtime::PhaseId phase,
                                   runtime::ActorId actor,
                                   runtime::FieldId field,
                                   const mesh::MeshTopology &topology) {
  try {
    require_face_mass_flux_field(registry, field);
    if (storage.layout_set().face_count != topology.local_face_count()) {
      throw Error("face mass flux storage and topology face counts differ");
    }
    auto view =
        storage.acquire_face_read<double>(access_plan, phase, actor, field);
    if (view.face_count() != topology.local_face_count() ||
        view.components() != 1U) {
      throw Error("face mass flux view layout is invalid");
    }
    inject_face_mass_flux_construction_failure();
    return FaceMassFlux(std::make_unique<State>(field, std::move(view),
                                                make_signature(topology)));
  } catch (const Error &) {
    throw;
  } catch (const std::bad_alloc &) {
    throw Error("face mass flux allocation failed");
  } catch (const std::length_error &) {
    throw Error("face mass flux allocation size is unsupported");
  }
}

FaceMassFlux::FaceMassFlux(std::unique_ptr<State> state) noexcept
    : field_(state->field), face_count_(state->view.face_count()),
      state_(std::move(state)) {}
FaceMassFlux::~FaceMassFlux() noexcept = default;
FaceMassFlux::FaceMassFlux(FaceMassFlux &&) noexcept = default;
runtime::FieldId FaceMassFlux::field_id() const noexcept { return field_; }
std::size_t FaceMassFlux::face_count() const noexcept { return face_count_; }

CellCenteredFvmOperators
CellCenteredFvmOperators::create(const mesh::MeshTopology &topology,
                                 const mesh::MeshGeometry &geometry) {
  try {
    geometry.require_compatible(topology);
    auto impl = std::make_unique<Impl>();
    impl->signature = make_signature(topology);
    const auto box = topology.owned_global_box();
    impl->local_extent = Int3{box.end.x - box.begin.x, box.end.y - box.begin.y,
                              box.end.z - box.begin.z};
    impl->faces.reserve(topology.local_face_count());
    impl->cell_faces.resize(topology.owned_cell_count());
    impl->scratch.resize(topology.owned_cell_count() * 9U);
    impl->secondary.resize(topology.owned_cell_count() * 9U);
    impl->inverse_cell_volumes.resize(topology.owned_cell_count());
    impl->face_scratch.resize(topology.local_face_count() * 3U);
    for (LocalCellId cell = 0; cell < topology.owned_cell_count(); ++cell) {
      const double volume = geometry.cell_volume_m3(cell);
      if (!(volume > 0.0) || !std::isfinite(volume)) {
        throw Error("finite-volume cell volume is invalid");
      }
      impl->inverse_cell_volumes[cell] = 1.0 / volume;
    }

    for (LocalFaceId local_face = 0; local_face < topology.local_face_count();
         ++local_face) {
      const LocalCellId owner = topology.owner(local_face);
      const auto neighbour = topology.neighbour(local_face);
      FaceStencil face{};
      face.local_face = local_face;
      face.global_face = topology.global_face_id(local_face);
      face.periodic_pair = topology.periodic_pair(local_face);
      face.axis = topology.logical_face(local_face).axis;
      face.owner =
          map_cell(topology.global_cell(owner), box, topology.global_extent());
      if (neighbour.has_value()) {
        face.neighbour = map_cell(topology.global_cell(*neighbour), box,
                                  topology.global_extent());
        if (face.periodic_pair.has_value()) {
          const auto logical = topology.logical_face(local_face);
          int plane = logical.coordinate.x;
          if (logical.axis == FaceAxis::y)
            plane = logical.coordinate.y;
          if (logical.axis == FaceAxis::z)
            plane = logical.coordinate.z;
          face.neighbour =
              shifted(face.owner, logical.axis, plane == 0 ? -1 : 1);
        }
      }
      if (topology.cell_ownership(owner) == EntityOwnership::owned) {
        face.owner_owned = owned_linear(face.owner, impl->local_extent);
      }
      if (neighbour.has_value() && !face.periodic_pair.has_value() &&
          topology.cell_ownership(*neighbour) == EntityOwnership::owned) {
        face.neighbour_owned =
            owned_linear(*face.neighbour, impl->local_extent);
      }
      face.patch = topology.patch_id(local_face);
      face.owner_area =
          geometry.face_area_vector_m2(local_face, FaceSide::owner);
      face.face_center = geometry.face_center_m(local_face);
      face.owner_center = geometry.cell_center_m(owner);
      if (neighbour.has_value()) {
        face.neighbour_center = geometry.cell_center_m(*neighbour);
      }
      face.displacement = geometry.face_displacement_m(local_face);
      face.skewness = geometry.face_skewness(local_face);
      face.non_orthogonality_degrees =
          geometry.face_non_orthogonality_degrees(local_face);
      if (!finite(face.owner_area) || !finite(face.face_center) ||
          !finite(face.owner_center) || !finite(face.displacement) ||
          (face.neighbour_center.has_value() &&
           !finite(*face.neighbour_center)) ||
          !std::isfinite(face.skewness) ||
          !std::isfinite(face.non_orthogonality_degrees)) {
        throw Error("finite-volume stencil geometry is non-finite");
      }
      if (face.neighbour.has_value() &&
          (topology.cell_ownership(owner) == EntityOwnership::ghost ||
           topology.cell_ownership(*neighbour) == EntityOwnership::ghost)) {
        impl->needs_remote_or_periodic = true;
      }
      if (face.periodic_pair.has_value()) {
        impl->needs_remote_or_periodic = true;
      }
      impl->faces.push_back(face);
    }

    const auto insert_ref = [&](std::size_t owned_cell, std::size_t face,
                                bool as_owner) {
      auto &refs = impl->cell_faces.at(owned_cell);
      const auto slot = std::find_if(refs.begin(), refs.end(), [](auto ref) {
        return ref.face == kInvalidFace;
      });
      if (slot == refs.end()) {
        throw Error("finite-volume cell has more than six incident faces");
      }
      *slot = CellFaceRef{face, as_owner};
    };
    for (std::size_t index = 0; index < impl->faces.size(); ++index) {
      const auto &face = impl->faces[index];
      if (face.owner_owned.has_value()) {
        insert_ref(*face.owner_owned, index, true);
      }
      if (!face.patch.has_value() && face.neighbour_owned.has_value()) {
        insert_ref(*face.neighbour_owned, index, false);
      }
    }
    for (const auto &refs : impl->cell_faces) {
      if (std::any_of(refs.begin(), refs.end(),
                      [](auto ref) { return ref.face == kInvalidFace; })) {
        throw Error("finite-volume cell does not have six incident faces");
      }
    }
    return CellCenteredFvmOperators(std::move(impl));
  } catch (const Error &) {
    throw;
  } catch (const std::bad_alloc &) {
    throw Error("finite-volume operator allocation failed");
  } catch (const std::length_error &) {
    throw Error("finite-volume operator allocation size is unsupported");
  }
}

CellCenteredFvmOperators::CellCenteredFvmOperators(
    std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
CellCenteredFvmOperators::~CellCenteredFvmOperators() noexcept = default;
CellCenteredFvmOperators::CellCenteredFvmOperators(
    CellCenteredFvmOperators &&) noexcept = default;

CellCenteredFvmOperators::Impl &CellCenteredFvmOperators::require_impl() const {
  if (!impl_) {
    throw Error("finite-volume operator has been moved from");
  }
  return *impl_;
}

void CellCenteredFvmOperators::compute_gradient(
    GradientScheme scheme, FiniteVolumeQuantity quantity,
    const boundary::BoundaryRegistry &boundaries,
    const runtime::FieldView<const double> &cell_values,
    const runtime::FieldView<double> &cell_gradients) const {
  OperationGuard operation(require_impl().active);
  validate_quantity(quantity, boundaries, true);
  const int components =
      quantity.kind == FiniteVolumeQuantityKind::velocity ? 3 : 1;
  const int gradient_components = components * 3;
  if (!same(cell_values.interior_extent(), impl_->local_extent) ||
      !same(cell_gradients.interior_extent(), impl_->local_extent) ||
      cell_values.components() != static_cast<std::uint32_t>(components) ||
      cell_gradients.components() !=
          static_cast<std::uint32_t>(gradient_components)) {
    throw Error("gradient field layout does not match the finite-volume mesh");
  }
  const int required_ghost = impl_->needs_remote_or_periodic ? 1 : 0;
  if (cell_values.ghost_width() < required_ghost ||
      cell_gradients.ghost_width() < required_ghost) {
    throw Error("gradient fields do not provide the required ghost width");
  }
  if (scheme != GradientScheme::green_gauss &&
      scheme != GradientScheme::weighted_least_squares) {
    throw Error("invalid finite-volume gradient scheme");
  }
  static_cast<void>(cell_gradients(0, 0, 0, 0));

  const bool force_singular = force_singular_least_squares;
  force_singular_least_squares = false;

  std::fill(impl_->scratch.begin(), impl_->scratch.end(), 0.0);
  for (std::size_t cell = 0; cell < impl_->cell_faces.size(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(impl_->local_extent.x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(impl_->local_extent.x);
    const int j =
        static_cast<int>(yz % static_cast<std::size_t>(impl_->local_extent.y));
    const int k =
        static_cast<int>(yz / static_cast<std::size_t>(impl_->local_extent.y));
    const StructuredIndex cell_index{i, j, k};
    for (int component = 0; component < components; ++component) {
      const double phi_p = value_at(cell_values, cell_index, component);
      if (!std::isfinite(phi_p)) {
        throw Error("gradient input contains a non-finite owned value");
      }
      Real3 result{};
      bool green_gauss_constant_stencil = true;
      std::array<double, 6> matrix{};
      Real3 right{};
      for (const CellFaceRef ref : impl_->cell_faces[cell]) {
        const FaceStencil &face = impl_->faces[ref.face];
        const bool physical = is_physical_nonperiodic(face);
        double phi_face = 0.0;
        double phi_q = 0.0;
        Real3 d{};
        Real3 area =
            ref.as_owner ? face.owner_area : multiply(-1.0, face.owner_area);
        if (physical) {
          if (!ref.as_owner) {
            throw Error("physical face is attached as a neighbour");
          }
          if (quantity.kind == FiniteVolumeQuantityKind::velocity) {
            const Real3 velocity{value_at(cell_values, cell_index, 0),
                                 value_at(cell_values, cell_index, 1),
                                 value_at(cell_values, cell_index, 2)};
            if (!finite(velocity)) {
              throw Error("gradient velocity input is non-finite");
            }
            const auto values = boundaries.evaluate_velocity(
                *face.patch, velocity, face.owner_area);
            phi_face = component == 0   ? values.face.x
                       : component == 1 ? values.face.y
                                        : values.face.z;
            phi_q = component == 0   ? values.exterior.x
                    : component == 1 ? values.exterior.y
                                     : values.exterior.z;
          } else {
            phi_face = boundary_face_value(boundaries, quantity, face, phi_p);
            phi_q = boundary_exterior_value(boundaries, quantity, face, Real3{},
                                            phi_p, component);
          }
          d = multiply(2.0, subtract(face.face_center, face.owner_center));
        } else {
          if (!face.neighbour.has_value()) {
            throw Error("finite-volume internal face has no neighbour");
          }
          const StructuredIndex other_index =
              ref.as_owner ? *face.neighbour : face.owner;
          phi_q = value_at(cell_values, other_index, component);
          if (!std::isfinite(phi_q)) {
            throw Error("gradient input contains a non-finite neighbour value");
          }
          d = ref.as_owner ? face.displacement
                           : multiply(-1.0, face.displacement);
          const Real3 owner_to_face =
              subtract(face.face_center, face.owner_center);
          const double a = safe_norm(owner_to_face);
          const double b =
              safe_norm(subtract(face.displacement, owner_to_face));
          const double denominator = a + b;
          if (!(a >= 0.0) || !(b >= 0.0) || !(denominator > 0.0) ||
              !std::isfinite(denominator)) {
            throw Error("gradient face interpolation distance is invalid");
          }
          const double owner_value =
              value_at(cell_values, face.owner, component);
          const double neighbour_value =
              value_at(cell_values, *face.neighbour, component);
          phi_face = owner_value == neighbour_value
                         ? owner_value
                         : (b / denominator) * owner_value +
                               (a / denominator) * neighbour_value;
        }
        if (!std::isfinite(phi_face) || !std::isfinite(phi_q) || !finite(d)) {
          throw Error("gradient stencil produced a non-finite value");
        }
        if (scheme == GradientScheme::green_gauss) {
          result = add(result, multiply(phi_face, area));
          green_gauss_constant_stencil =
              green_gauss_constant_stencil && phi_face == phi_p;
        } else {
          const double distance_squared = dot(d, d);
          if (!(distance_squared > 0.0) || !std::isfinite(distance_squared)) {
            throw Error("least-squares stencil distance is invalid");
          }
          const double weight = 1.0 / distance_squared;
          const double delta = phi_q - phi_p;
          matrix[0] += weight * d.x * d.x;
          matrix[1] += weight * d.x * d.y;
          matrix[2] += weight * d.x * d.z;
          matrix[3] += weight * d.y * d.y;
          matrix[4] += weight * d.y * d.z;
          matrix[5] += weight * d.z * d.z;
          right = add(right, multiply(weight * delta, d));
        }
      }
      if (scheme == GradientScheme::green_gauss) {
        if (green_gauss_constant_stencil) {
          result = Real3{};
        }
        // Volume scaling is applied in a fixed pass after all face sums.
      } else {
        if (force_singular && cell == 0U && component == 0) {
          matrix[5] = 0.0;
          matrix[2] = 0.0;
          matrix[4] = 0.0;
        }
        const double scale_x = std::sqrt(matrix[0]);
        const double scale_y = std::sqrt(matrix[3]);
        const double scale_z = std::sqrt(matrix[5]);
        if (!(scale_x > 0.0) || !(scale_y > 0.0) || !(scale_z > 0.0) ||
            !std::isfinite(scale_x) || !std::isfinite(scale_y) ||
            !std::isfinite(scale_z)) {
          throw Error("least-squares system is singular");
        }
        const double c10 = matrix[1] / (scale_x * scale_y);
        const double c20 = matrix[2] / (scale_x * scale_z);
        const double c21 = matrix[4] / (scale_y * scale_z);
        const double p11 = 1.0 - c10 * c10;
        constexpr double pivot_floor =
            128.0 * std::numeric_limits<double>::epsilon();
        if (!(p11 > pivot_floor) || !std::isfinite(p11)) {
          throw Error("least-squares system is singular");
        }
        const double l11 = std::sqrt(p11);
        const double l21 = (c21 - c20 * c10) / l11;
        const double p22 = 1.0 - c20 * c20 - l21 * l21;
        if (!(p22 > pivot_floor) || !std::isfinite(p22)) {
          throw Error("least-squares system is singular");
        }
        const double l22 = std::sqrt(p22);
        const double scaled_rhs_x = right.x / scale_x;
        const double scaled_rhs_y = right.y / scale_y;
        const double scaled_rhs_z = right.z / scale_z;
        const double y0 = scaled_rhs_x;
        const double y1 = (scaled_rhs_y - c10 * y0) / l11;
        const double y2 = (scaled_rhs_z - c20 * y0 - l21 * y1) / l22;
        const double scaled_result_z = y2 / l22;
        const double scaled_result_y = (y1 - l21 * scaled_result_z) / l11;
        const double scaled_result_x =
            y0 - c10 * scaled_result_y - c20 * scaled_result_z;
        result.x = scaled_result_x / scale_x;
        result.y = scaled_result_y / scale_y;
        result.z = scaled_result_z / scale_z;
      }
      const std::size_t base =
          cell * 9U + static_cast<std::size_t>(component) * 3U;
      impl_->scratch[base] = result.x;
      impl_->scratch[base + 1U] = result.y;
      impl_->scratch[base + 2U] = result.z;
    }
  }

  // Green--Gauss volume scaling is applied in a second fixed pass.
  if (scheme == GradientScheme::green_gauss) {
    for (std::size_t cell = 0; cell < impl_->cell_faces.size(); ++cell) {
      const double inverse_volume = impl_->inverse_cell_volumes[cell];
      for (int component = 0; component < components; ++component) {
        const std::size_t base =
            cell * 9U + static_cast<std::size_t>(component) * 3U;
        impl_->scratch[base] *= inverse_volume;
        impl_->scratch[base + 1U] *= inverse_volume;
        impl_->scratch[base + 2U] *= inverse_volume;
      }
    }
  }
  for (double value : impl_->scratch) {
    if (!std::isfinite(value)) {
      throw Error("gradient output is non-finite");
    }
  }
  for (std::size_t cell = 0; cell < impl_->cell_faces.size(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(impl_->local_extent.x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(impl_->local_extent.x);
    const int j =
        static_cast<int>(yz % static_cast<std::size_t>(impl_->local_extent.y));
    const int k =
        static_cast<int>(yz / static_cast<std::size_t>(impl_->local_extent.y));
    for (int component = 0; component < gradient_components; ++component) {
      cell_gradients(i, j, k, component) =
          impl_->scratch[cell * 9U + static_cast<std::size_t>(component)];
    }
  }
}

void CellCenteredFvmOperators::reconstruct_transport_faces(
    FiniteVolumeQuantity quantity, const boundary::BoundaryRegistry &boundaries,
    const FaceMassFlux &mass_flux,
    const runtime::FieldView<const double> &cell_values,
    const runtime::FaceFieldView<double> &face_values) const {
  OperationGuard operation(require_impl().active);
  validate_quantity(quantity, boundaries, false);
  if (!mass_flux.state_) {
    throw Error("face mass flux handle has been moved from");
  }
  if (!same_signature(impl_->signature, mass_flux.state_->signature) ||
      mass_flux.state_->view.face_count() != impl_->faces.size()) {
    throw Error("face mass flux topology identity does not match the operator");
  }
  if (!same(cell_values.interior_extent(), impl_->local_extent) ||
      cell_values.components() != 1U || cell_values.ghost_width() < 2 ||
      face_values.face_count() != impl_->faces.size() ||
      face_values.components() != 1U) {
    throw Error("transport reconstruction field layout is invalid");
  }
  static_cast<void>(face_values(0, 0));
  for (const auto &face : impl_->faces) {
    double reconstructed = 0.0;
    if (is_physical_nonperiodic(face)) {
      const double interior = value_at(cell_values, face.owner, 0);
      if (!std::isfinite(interior)) {
        throw Error("transport boundary interior value is non-finite");
      }
      reconstructed = boundary_face_value(boundaries, quantity, face, interior);
    } else {
      const CanonicalFace canonical = canonical_face(face);
      const double stored_mass_flux =
          mass_flux.state_->view(face.local_face, 0);
      if (!std::isfinite(stored_mass_flux)) {
        throw Error("face mass flux contains a non-finite value");
      }
      const double canonical_mass_flux =
          canonical.reversed ? -stored_mass_flux : stored_mass_flux;
      reconstructed =
          reconstruct_mc(cell_values, face, canonical, 0, canonical_mass_flux);
    }
    if (!std::isfinite(reconstructed)) {
      throw Error("transport face reconstruction is non-finite");
    }
    impl_->face_scratch[face.local_face * 3U] = reconstructed;
  }
  for (const auto &face : impl_->faces) {
    face_values(face.local_face, 0) = impl_->face_scratch[face.local_face * 3U];
  }
}

void CellCenteredFvmOperators::reconstruct_momentum_faces(
    const boundary::BoundaryRegistry &boundaries, const FaceMassFlux &mass_flux,
    const runtime::FieldView<const double> &velocity,
    const runtime::FaceFieldView<double> &face_velocity) const {
  OperationGuard operation(require_impl().active);
  if (!mass_flux.state_) {
    throw Error("face mass flux handle has been moved from");
  }
  if (!same_signature(impl_->signature, mass_flux.state_->signature) ||
      !same(velocity.interior_extent(), impl_->local_extent) ||
      velocity.components() != 3U || velocity.ghost_width() < 2 ||
      face_velocity.face_count() != impl_->faces.size() ||
      face_velocity.components() != 3U) {
    throw Error("momentum reconstruction field layout is invalid");
  }
  static_cast<void>(face_velocity(0, 0));
  const MetricOverride metric_override = consume_metric_override();
  bool override_observed = !metric_override.armed;
  for (const auto &face : impl_->faces) {
    if (is_physical_nonperiodic(face)) {
      const Real3 interior{value_at(velocity, face.owner, 0),
                           value_at(velocity, face.owner, 1),
                           value_at(velocity, face.owner, 2)};
      if (!finite(interior)) {
        throw Error("momentum boundary interior value is non-finite");
      }
      const auto values =
          boundaries.evaluate_velocity(*face.patch, interior, face.owner_area);
      impl_->face_scratch[face.local_face * 3U] = values.face.x;
      impl_->face_scratch[face.local_face * 3U + 1U] = values.face.y;
      impl_->face_scratch[face.local_face * 3U + 2U] = values.face.z;
      continue;
    }
    const CanonicalFace canonical = canonical_face(face);
    double skewness = face.skewness;
    double non_orthogonality = face.non_orthogonality_degrees;
    if (metric_override.armed &&
        (face.global_face == metric_override.global_face ||
         face.periodic_pair == metric_override.global_face)) {
      skewness = metric_override.skewness;
      non_orthogonality = metric_override.non_orthogonality_degrees;
      override_observed = true;
    }
    const bool use_mc = skewness > 0.25 || non_orthogonality > 70.0;
    const double stored_mass_flux = mass_flux.state_->view(face.local_face, 0);
    if (!std::isfinite(stored_mass_flux)) {
      throw Error("face mass flux contains a non-finite value");
    }
    const double canonical_mass_flux =
        canonical.reversed ? -stored_mass_flux : stored_mass_flux;
    for (int component = 0; component < 3; ++component) {
      double reconstructed = 0.0;
      if (use_mc) {
        reconstructed = reconstruct_mc(velocity, face, canonical, component,
                                       canonical_mass_flux);
      } else {
        const double p = value_at(velocity, canonical.p, component);
        const double n = value_at(velocity, canonical.n, component);
        reconstructed = 0.5 * p + 0.5 * n;
      }
      if (!std::isfinite(reconstructed)) {
        throw Error("momentum face reconstruction is non-finite");
      }
      impl_->face_scratch[face.local_face * 3U +
                          static_cast<std::size_t>(component)] = reconstructed;
    }
  }
  if (!override_observed) {
    throw Error("finite-volume metric override face was not observed");
  }
  for (const auto &face : impl_->faces) {
    for (int component = 0; component < 3; ++component) {
      face_velocity(face.local_face, component) =
          impl_->face_scratch[face.local_face * 3U +
                              static_cast<std::size_t>(component)];
    }
  }
}

void CellCenteredFvmOperators::assemble_provisional_mass_flux(
    const boundary::BoundaryRegistry &boundaries,
    const runtime::FieldView<const double> &density,
    const runtime::FieldView<const double> &velocity,
    const runtime::FieldRegistry &registry, runtime::FieldStorage &storage,
    const runtime::FieldAccessPlan &access_plan, runtime::PhaseId phase,
    runtime::ActorId actor, runtime::FieldId mass_flux_field) const {
  OperationGuard operation(require_impl().active);
  require_face_mass_flux_field(registry, mass_flux_field);
  if (!same(density.interior_extent(), impl_->local_extent) ||
      density.components() != 1U || density.ghost_width() < 2 ||
      !same(velocity.interior_extent(), impl_->local_extent) ||
      velocity.components() != 3U || velocity.ghost_width() < 1 ||
      storage.layout_set().face_count != impl_->faces.size()) {
    throw Error("provisional mass-flux field layout is invalid");
  }
  auto writer = storage.acquire_face_write<double>(access_plan, phase, actor,
                                                   mass_flux_field);
  if (writer.face_count() != impl_->faces.size() || writer.components() != 1U) {
    throw Error("provisional mass-flux writer layout is invalid");
  }
  for (const auto &face : impl_->faces) {
    double flux = 0.0;
    if (is_physical_nonperiodic(face)) {
      const auto &descriptor = boundaries.patch(*face.patch);
      switch (descriptor.mass_flux_rule()) {
      case boundary::MassFluxRule::identically_zero:
        flux = 0.0;
        break;
      case boundary::MassFluxRule::prescribed_inlet_state: {
        if (!descriptor.inlet_state().has_value()) {
          throw Error("velocity inlet has no materialized inlet state");
        }
        const auto &inlet = *descriptor.inlet_state();
        flux = inlet.density_kg_per_m3 *
               dot(inlet.velocity_m_per_s, face.owner_area);
        break;
      }
      case boundary::MassFluxRule::outflow_only: {
        const double rho_p = value_at(density, face.owner, 0);
        const Real3 u_p{value_at(velocity, face.owner, 0),
                        value_at(velocity, face.owner, 1),
                        value_at(velocity, face.owner, 2)};
        if (!std::isfinite(rho_p) || !finite(u_p)) {
          throw Error("pressure-outlet interior state is non-finite");
        }
        const double rho_f =
            boundaries.evaluate_density(*face.patch, rho_p).face;
        const Real3 u_f =
            boundaries.evaluate_velocity(*face.patch, u_p, face.owner_area)
                .face;
        if (!(rho_f > 0.0) || !std::isfinite(rho_f) || !finite(u_f)) {
          throw Error("pressure-outlet face state is invalid");
        }
        flux = rho_f * dot(u_f, face.owner_area);
        break;
      }
      default:
        throw Error("non-periodic boundary has an invalid mass-flux rule");
      }
    } else {
      const CanonicalFace canonical = canonical_face(face);
      Real3 u_face{};
      for (int component = 0; component < 3; ++component) {
        const double p = value_at(velocity, canonical.p, component);
        const double n = value_at(velocity, canonical.n, component);
        const double candidate = 0.5 * p + 0.5 * n;
        if (!std::isfinite(candidate)) {
          throw Error("provisional face velocity is non-finite");
        }
        if (component == 0)
          u_face.x = candidate;
        if (component == 1)
          u_face.y = candidate;
        if (component == 2)
          u_face.z = candidate;
      }
      const double volumetric_flux = dot(u_face, canonical.area);
      if (!std::isfinite(volumetric_flux)) {
        throw Error("provisional volumetric flux is non-finite");
      }
      const double rho_face =
          reconstruct_mc(density, face, canonical, 0, volumetric_flux);
      if (!(rho_face > 0.0)) {
        throw Error("provisional reconstructed density is not positive");
      }
      double canonical_flux = rho_face * volumetric_flux;
      if (canonical_flux == 0.0) {
        canonical_flux = 0.0;
      }
      flux = canonical.reversed ? -canonical_flux : canonical_flux;
    }
    if (!std::isfinite(flux)) {
      throw Error("provisional mass flux is non-finite");
    }
    impl_->face_scratch[face.local_face * 3U] = flux;
  }
  for (const auto &face : impl_->faces) {
    writer(face.local_face, 0) = impl_->face_scratch[face.local_face * 3U];
  }
}

void CellCenteredFvmOperators::accumulate_mass_residual(
    const FaceMassFlux &mass_flux,
    const runtime::FieldView<double> &raw_residual) const {
  OperationGuard operation(require_impl().active);
  if (!mass_flux.state_) {
    throw Error("face mass flux handle has been moved from");
  }
  if (!same_signature(impl_->signature, mass_flux.state_->signature) ||
      !same(raw_residual.interior_extent(), impl_->local_extent) ||
      raw_residual.components() != 1U || raw_residual.ghost_width() < 0) {
    throw Error("mass residual field layout is invalid");
  }
  for (std::size_t cell = 0; cell < impl_->cell_faces.size(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(impl_->local_extent.x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(impl_->local_extent.x);
    const int j =
        static_cast<int>(yz % static_cast<std::size_t>(impl_->local_extent.y));
    const int k =
        static_cast<int>(yz / static_cast<std::size_t>(impl_->local_extent.y));
    const double initial = raw_residual(i, j, k, 0);
    if (!std::isfinite(initial)) {
      throw Error("mass residual contains a non-finite input");
    }
    impl_->secondary[cell * 9U] = initial;
  }
  for (const auto &face : impl_->faces) {
    const double flux = mass_flux.state_->view(face.local_face, 0);
    if (!std::isfinite(flux)) {
      throw Error("face mass flux contains a non-finite value");
    }
    if (face.owner_owned.has_value()) {
      const std::size_t cell = *face.owner_owned;
      const double candidate = impl_->secondary[cell * 9U] + flux;
      if (!std::isfinite(candidate)) {
        throw Error("mass residual accumulation is non-finite");
      }
      impl_->secondary[cell * 9U] = candidate;
    }
    if (!face.patch.has_value() && face.neighbour_owned.has_value()) {
      const std::size_t cell = *face.neighbour_owned;
      const double candidate = impl_->secondary[cell * 9U] - flux;
      if (!std::isfinite(candidate)) {
        throw Error("mass residual accumulation is non-finite");
      }
      impl_->secondary[cell * 9U] = candidate;
    }
  }
  for (std::size_t cell = 0; cell < impl_->cell_faces.size(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(impl_->local_extent.x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(impl_->local_extent.x);
    const int j =
        static_cast<int>(yz % static_cast<std::size_t>(impl_->local_extent.y));
    const int k =
        static_cast<int>(yz / static_cast<std::size_t>(impl_->local_extent.y));
    raw_residual(i, j, k, 0) = impl_->secondary[cell * 9U];
  }
}

void CellCenteredFvmOperators::accumulate_convective_residual(
    const FaceMassFlux &mass_flux,
    const runtime::FaceFieldView<const double> &transported_face_values,
    const runtime::FieldView<double> &raw_residual) const {
  OperationGuard operation(require_impl().active);
  if (!mass_flux.state_) {
    throw Error("face mass flux handle has been moved from");
  }
  const std::uint32_t components = transported_face_values.components();
  if (!same_signature(impl_->signature, mass_flux.state_->signature) ||
      transported_face_values.face_count() != impl_->faces.size() ||
      (components != 1U && components != 3U) ||
      !same(raw_residual.interior_extent(), impl_->local_extent) ||
      raw_residual.components() != components ||
      raw_residual.ghost_width() < 0) {
    throw Error("convective residual field layout is invalid");
  }
  for (std::size_t cell = 0; cell < impl_->cell_faces.size(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(impl_->local_extent.x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(impl_->local_extent.x);
    const int j =
        static_cast<int>(yz % static_cast<std::size_t>(impl_->local_extent.y));
    const int k =
        static_cast<int>(yz / static_cast<std::size_t>(impl_->local_extent.y));
    for (std::uint32_t component = 0; component < components; ++component) {
      const double initial = raw_residual(i, j, k, static_cast<int>(component));
      if (!std::isfinite(initial)) {
        throw Error("convective residual contains a non-finite input");
      }
      impl_->secondary[cell * 9U + component] = initial;
    }
  }
  for (const auto &face : impl_->faces) {
    const double flux = mass_flux.state_->view(face.local_face, 0);
    if (!std::isfinite(flux)) {
      throw Error("face mass flux contains a non-finite value");
    }
    for (std::uint32_t component = 0; component < components; ++component) {
      const double transported =
          transported_face_values(face.local_face, static_cast<int>(component));
      const double contribution = flux * transported;
      if (!std::isfinite(transported) || !std::isfinite(contribution)) {
        throw Error("convective face contribution is non-finite");
      }
      if (face.owner_owned.has_value()) {
        const std::size_t offset = *face.owner_owned * 9U + component;
        const double candidate = impl_->secondary[offset] + contribution;
        if (!std::isfinite(candidate)) {
          throw Error("convective residual accumulation is non-finite");
        }
        impl_->secondary[offset] = candidate;
      }
      if (!face.patch.has_value() && face.neighbour_owned.has_value()) {
        const std::size_t offset = *face.neighbour_owned * 9U + component;
        const double candidate = impl_->secondary[offset] - contribution;
        if (!std::isfinite(candidate)) {
          throw Error("convective residual accumulation is non-finite");
        }
        impl_->secondary[offset] = candidate;
      }
    }
  }
  for (std::size_t cell = 0; cell < impl_->cell_faces.size(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(impl_->local_extent.x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(impl_->local_extent.x);
    const int j =
        static_cast<int>(yz % static_cast<std::size_t>(impl_->local_extent.y));
    const int k =
        static_cast<int>(yz / static_cast<std::size_t>(impl_->local_extent.y));
    for (std::uint32_t component = 0; component < components; ++component) {
      raw_residual(i, j, k, static_cast<int>(component)) =
          impl_->secondary[cell * 9U + component];
    }
  }
}

void CellCenteredFvmOperators::accumulate_scalar_diffusive_residual(
    FiniteVolumeQuantity quantity, const boundary::BoundaryRegistry &boundaries,
    const runtime::FieldView<const double> &cell_values,
    const runtime::FieldView<const double> &cell_gradients,
    const runtime::FaceFieldView<const double> &gamma_by_face,
    const runtime::FieldView<double> &raw_residual) const {
  OperationGuard operation(require_impl().active);
  validate_quantity(quantity, boundaries, false);
  const int required_ghost = impl_->needs_remote_or_periodic ? 1 : 0;
  if (!same(cell_values.interior_extent(), impl_->local_extent) ||
      cell_values.components() != 1U ||
      cell_values.ghost_width() < required_ghost ||
      !same(cell_gradients.interior_extent(), impl_->local_extent) ||
      cell_gradients.components() != 3U ||
      cell_gradients.ghost_width() < required_ghost ||
      gamma_by_face.face_count() != impl_->faces.size() ||
      gamma_by_face.components() != 1U ||
      !same(raw_residual.interior_extent(), impl_->local_extent) ||
      raw_residual.components() != 1U || raw_residual.ghost_width() < 0) {
    throw Error("scalar diffusion field layout is invalid");
  }
  for (std::size_t cell = 0; cell < impl_->cell_faces.size(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(impl_->local_extent.x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(impl_->local_extent.x);
    const int j =
        static_cast<int>(yz % static_cast<std::size_t>(impl_->local_extent.y));
    const int k =
        static_cast<int>(yz / static_cast<std::size_t>(impl_->local_extent.y));
    const double initial = raw_residual(i, j, k, 0);
    if (!std::isfinite(initial)) {
      throw Error("scalar diffusion residual input is non-finite");
    }
    impl_->secondary[cell * 9U] = initial;
  }
  for (const auto &face : impl_->faces) {
    const double gamma = gamma_by_face(face.local_face, 0);
    if (!(gamma >= 0.0) || !std::isfinite(gamma)) {
      throw Error("scalar diffusion coefficient is invalid");
    }
    double flux = 0.0;
    if (is_physical_nonperiodic(face)) {
      const auto &descriptor = boundaries.patch(*face.patch);
      boundary::TransportRule rule = descriptor.density_rule();
      if (quantity.kind == FiniteVolumeQuantityKind::enthalpy) {
        rule = descriptor.enthalpy_rule();
      } else if (quantity.kind == FiniteVolumeQuantityKind::scalar) {
        rule = descriptor.scalar_rule();
      }
      if (rule == boundary::TransportRule::prescribed_value) {
        const double q_p = value_at(cell_values, face.owner, 0);
        const double q_n = boundary_exterior_value(boundaries, quantity, face,
                                                   Real3{}, q_p, 0);
        const Real3 d =
            multiply(2.0, subtract(face.face_center, face.owner_center));
        const double sd = dot(face.owner_area, d);
        const double factor = dot(face.owner_area, face.owner_area) / sd;
        const Real3 nonorth = subtract(face.owner_area, multiply(factor, d));
        const Real3 gradient{value_at(cell_gradients, face.owner, 0),
                             value_at(cell_gradients, face.owner, 1),
                             value_at(cell_gradients, face.owner, 2)};
        if (!std::isfinite(q_p) || !std::isfinite(q_n) || !finite(d) ||
            !finite(gradient) || !(sd > 0.0) || !std::isfinite(factor)) {
          throw Error("prescribed scalar diffusion stencil is invalid");
        }
        flux = gamma * ((q_n - q_p) * factor + dot(gradient, nonorth));
      } else if (rule == boundary::TransportRule::zero_normal_diffusive_flux ||
                 rule == boundary::TransportRule::copy_interior ||
                 rule == boundary::TransportRule::pure_outflow) {
        flux = 0.0;
      } else {
        throw Error("physical scalar diffusion rule is invalid");
      }
    } else {
      const CanonicalFace canonical = canonical_face(face);
      const Real3 d = canonical.reversed ? multiply(-1.0, face.displacement)
                                         : face.displacement;
      const double sd = dot(canonical.area, d);
      const double factor = dot(canonical.area, canonical.area) / sd;
      if (!(sd > 0.0) || !std::isfinite(sd) || !std::isfinite(factor)) {
        throw Error("scalar diffusion orthogonal projection is invalid");
      }
      const Real3 owner_to_face = subtract(face.face_center, face.owner_center);
      double a = safe_norm(owner_to_face);
      double b = safe_norm(subtract(face.displacement, owner_to_face));
      if (canonical.reversed)
        std::swap(a, b);
      const double denominator = a + b;
      if (!(denominator > 0.0) || !std::isfinite(denominator)) {
        throw Error("scalar diffusion interpolation distance is invalid");
      }
      const double q_p = value_at(cell_values, canonical.p, 0);
      const double q_n = value_at(cell_values, canonical.n, 0);
      Real3 gradient{};
      for (int direction = 0; direction < 3; ++direction) {
        const double gp = value_at(cell_gradients, canonical.p, direction);
        const double gn = value_at(cell_gradients, canonical.n, direction);
        const double candidate =
            (b / denominator) * gp + (a / denominator) * gn;
        if (direction == 0)
          gradient.x = candidate;
        if (direction == 1)
          gradient.y = candidate;
        if (direction == 2)
          gradient.z = candidate;
      }
      const Real3 nonorth = subtract(canonical.area, multiply(factor, d));
      const double canonical_flux =
          gamma * ((q_n - q_p) * factor + dot(gradient, nonorth));
      flux = canonical.reversed ? -canonical_flux : canonical_flux;
    }
    if (!std::isfinite(flux)) {
      throw Error("scalar diffusion contribution is non-finite");
    }
    if (face.owner_owned.has_value()) {
      const std::size_t offset = *face.owner_owned * 9U;
      const double candidate = impl_->secondary[offset] - flux;
      if (!std::isfinite(candidate)) {
        throw Error("scalar diffusion residual accumulation is non-finite");
      }
      impl_->secondary[offset] = candidate;
    }
    if (!face.patch.has_value() && face.neighbour_owned.has_value()) {
      const std::size_t offset = *face.neighbour_owned * 9U;
      const double candidate = impl_->secondary[offset] + flux;
      if (!std::isfinite(candidate)) {
        throw Error("scalar diffusion residual accumulation is non-finite");
      }
      impl_->secondary[offset] = candidate;
    }
  }
  for (std::size_t cell = 0; cell < impl_->cell_faces.size(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(impl_->local_extent.x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(impl_->local_extent.x);
    const int j =
        static_cast<int>(yz % static_cast<std::size_t>(impl_->local_extent.y));
    const int k =
        static_cast<int>(yz / static_cast<std::size_t>(impl_->local_extent.y));
    raw_residual(i, j, k, 0) = impl_->secondary[cell * 9U];
  }
}

void CellCenteredFvmOperators::accumulate_viscous_residual(
    const boundary::BoundaryRegistry &boundaries,
    const runtime::FieldView<const double> &velocity,
    const runtime::FieldView<const double> &velocity_gradients,
    double dynamic_viscosity_pa_s,
    const runtime::FieldView<double> &raw_momentum_residual) const {
  OperationGuard operation(require_impl().active);
  const int required_ghost = impl_->needs_remote_or_periodic ? 1 : 0;
  if (!(dynamic_viscosity_pa_s >= 0.0) ||
      !std::isfinite(dynamic_viscosity_pa_s) ||
      !same(velocity.interior_extent(), impl_->local_extent) ||
      velocity.components() != 3U || velocity.ghost_width() < required_ghost ||
      !same(velocity_gradients.interior_extent(), impl_->local_extent) ||
      velocity_gradients.components() != 9U ||
      velocity_gradients.ghost_width() < required_ghost ||
      !same(raw_momentum_residual.interior_extent(), impl_->local_extent) ||
      raw_momentum_residual.components() != 3U ||
      raw_momentum_residual.ghost_width() < 0) {
    throw Error("viscous residual field layout is invalid");
  }
  for (std::size_t cell = 0; cell < impl_->cell_faces.size(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(impl_->local_extent.x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(impl_->local_extent.x);
    const int j =
        static_cast<int>(yz % static_cast<std::size_t>(impl_->local_extent.y));
    const int k =
        static_cast<int>(yz / static_cast<std::size_t>(impl_->local_extent.y));
    for (int component = 0; component < 3; ++component) {
      const double initial = raw_momentum_residual(i, j, k, component);
      if (!std::isfinite(initial)) {
        throw Error("viscous residual input is non-finite");
      }
      impl_->secondary[cell * 9U + static_cast<std::size_t>(component)] =
          initial;
    }
  }
  for (const auto &face : impl_->faces) {
    Real3 d{};
    Real3 area{};
    StructuredIndex p{};
    StructuredIndex n{};
    bool reversed = false;
    bool physical = is_physical_nonperiodic(face);
    Real3 u_p{};
    Real3 u_n{};
    std::array<double, 9> gradient{};
    if (physical) {
      p = face.owner;
      d = multiply(2.0, subtract(face.face_center, face.owner_center));
      area = face.owner_area;
      u_p = Real3{value_at(velocity, p, 0), value_at(velocity, p, 1),
                  value_at(velocity, p, 2)};
      u_n = boundaries.evaluate_velocity(*face.patch, u_p, area).exterior;
      for (int value = 0; value < 9; ++value) {
        gradient[static_cast<std::size_t>(value)] =
            value_at(velocity_gradients, p, value);
      }
    } else {
      const CanonicalFace canonical = canonical_face(face);
      p = canonical.p;
      n = canonical.n;
      area = canonical.area;
      reversed = canonical.reversed;
      d = reversed ? multiply(-1.0, face.displacement) : face.displacement;
      u_p = Real3{value_at(velocity, p, 0), value_at(velocity, p, 1),
                  value_at(velocity, p, 2)};
      u_n = Real3{value_at(velocity, n, 0), value_at(velocity, n, 1),
                  value_at(velocity, n, 2)};
      const Real3 owner_to_face = subtract(face.face_center, face.owner_center);
      double a = safe_norm(owner_to_face);
      double b = safe_norm(subtract(face.displacement, owner_to_face));
      if (reversed)
        std::swap(a, b);
      const double denominator = a + b;
      if (!(denominator > 0.0) || !std::isfinite(denominator)) {
        throw Error("viscous interpolation distance is invalid");
      }
      for (int value = 0; value < 9; ++value) {
        const double gp = value_at(velocity_gradients, p, value);
        const double gn = value_at(velocity_gradients, n, value);
        gradient[static_cast<std::size_t>(value)] =
            (b / denominator) * gp + (a / denominator) * gn;
      }
    }
    const double d2 = dot(d, d);
    if (!finite(d) || !finite(area) || !finite(u_p) || !finite(u_n) ||
        !(d2 > 0.0) || !std::isfinite(d2)) {
      throw Error("viscous face stencil is invalid");
    }
    const std::array<double, 3> delta{u_n.x - u_p.x, u_n.y - u_p.y,
                                      u_n.z - u_p.z};
    const std::array<double, 3> d_values{d.x, d.y, d.z};
    for (int component = 0; component < 3; ++component) {
      const std::size_t base = static_cast<std::size_t>(component) * 3U;
      const double projected = gradient[base] * d.x +
                               gradient[base + 1U] * d.y +
                               gradient[base + 2U] * d.z;
      const double correction =
          (delta[static_cast<std::size_t>(component)] - projected) / d2;
      for (int direction = 0; direction < 3; ++direction) {
        gradient[base + static_cast<std::size_t>(direction)] +=
            correction * d_values[static_cast<std::size_t>(direction)];
      }
    }
    const double divergence = gradient[0] + gradient[4] + gradient[8];
    const std::array<double, 3> area_values{area.x, area.y, area.z};
    std::array<double, 3> traction{};
    for (int component = 0; component < 3; ++component) {
      for (int direction = 0; direction < 3; ++direction) {
        double stress = dynamic_viscosity_pa_s *
                        (gradient[static_cast<std::size_t>(component) * 3U +
                                  static_cast<std::size_t>(direction)] +
                         gradient[static_cast<std::size_t>(direction) * 3U +
                                  static_cast<std::size_t>(component)]);
        if (component == direction) {
          stress -= dynamic_viscosity_pa_s * (2.0 / 3.0) * divergence;
        }
        traction[static_cast<std::size_t>(component)] +=
            stress * area_values[static_cast<std::size_t>(direction)];
      }
    }
    if (reversed) {
      for (double &component : traction)
        component = -component;
    }
    for (double component : traction) {
      if (!std::isfinite(component)) {
        throw Error("viscous traction is non-finite");
      }
    }
    if (face.owner_owned.has_value()) {
      for (std::size_t component = 0; component < 3U; ++component) {
        const std::size_t offset = *face.owner_owned * 9U + component;
        const double candidate = impl_->secondary[offset] - traction[component];
        if (!std::isfinite(candidate)) {
          throw Error("viscous residual accumulation is non-finite");
        }
        impl_->secondary[offset] = candidate;
      }
    }
    if (!face.patch.has_value() && face.neighbour_owned.has_value()) {
      for (std::size_t component = 0; component < 3U; ++component) {
        const std::size_t offset = *face.neighbour_owned * 9U + component;
        const double candidate = impl_->secondary[offset] + traction[component];
        if (!std::isfinite(candidate)) {
          throw Error("viscous residual accumulation is non-finite");
        }
        impl_->secondary[offset] = candidate;
      }
    }
  }
  for (std::size_t cell = 0; cell < impl_->cell_faces.size(); ++cell) {
    const int i = static_cast<int>(
        cell % static_cast<std::size_t>(impl_->local_extent.x));
    const std::size_t yz =
        cell / static_cast<std::size_t>(impl_->local_extent.x);
    const int j =
        static_cast<int>(yz % static_cast<std::size_t>(impl_->local_extent.y));
    const int k =
        static_cast<int>(yz / static_cast<std::size_t>(impl_->local_extent.y));
    for (int component = 0; component < 3; ++component) {
      raw_momentum_residual(i, j, k, component) =
          impl_->secondary[cell * 9U + static_cast<std::size_t>(component)];
    }
  }
}

} // namespace hundun::finite_volume

// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_ibm.hpp"
#include "hundun/v04_flow.hpp"

#include "mesh_focus_detail.hpp"
#include "solver_equation_detail.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <utility>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kIbmEquationPlan = 9241U;
constexpr std::uint32_t kIbmEquationApply = 9242U;
constexpr std::uint32_t kIbmEquationNumerical = 9243U;
constexpr std::uint32_t kIbmPhysicalBoundaryPlan = 9244U;
constexpr std::uint32_t kIbmPhysicalBoundaryApply = 9245U;
constexpr std::uint32_t kIbmPhysicalBoundaryNumerical = 9246U;
constexpr std::uint64_t kFnvOffset = UINT64_C(1469598103934665603);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

std::uint64_t mix(std::uint64_t hash, std::uint64_t value) noexcept {
  hash ^= value;
  hash *= kFnvPrime;
  return hash;
}

bool same_patch(MeshPatch left, MeshPatch right) noexcept {
  return left.begin.x == right.begin.x && left.begin.y == right.begin.y &&
         left.begin.z == right.begin.z && left.cells.x == right.cells.x &&
         left.cells.y == right.cells.y && left.cells.z == right.cells.z &&
         left.process_grid.x == right.process_grid.x &&
         left.process_grid.y == right.process_grid.y &&
         left.process_grid.z == right.process_grid.z &&
         left.process_coord.x == right.process_coord.x &&
         left.process_coord.y == right.process_coord.y &&
         left.process_coord.z == right.process_coord.z;
}

bool mpi_live() noexcept {
  int initialized = 0;
  int finalized = 0;
  return MPI_Initialized(&initialized) == MPI_SUCCESS && initialized != 0 &&
         MPI_Finalized(&finalized) == MPI_SUCCESS && finalized == 0;
}

struct InterfaceFace {
  CartesianAxis axis{CartesianAxis::x};
  Int3 index{};
};

InterfaceFace interface_face(const ImmersedLink& link) noexcept {
  InterfaceFace result;
  result.index = link.fluid_local_index;
  switch (link.direction) {
    case ImmersedFaceDirection::x_negative:
      result.axis = CartesianAxis::x;
      break;
    case ImmersedFaceDirection::x_positive:
      result.axis = CartesianAxis::x;
      ++result.index.x;
      break;
    case ImmersedFaceDirection::y_negative:
      result.axis = CartesianAxis::y;
      break;
    case ImmersedFaceDirection::y_positive:
      result.axis = CartesianAxis::y;
      ++result.index.y;
      break;
    case ImmersedFaceDirection::z_negative:
      result.axis = CartesianAxis::z;
      break;
    case ImmersedFaceDirection::z_positive:
      result.axis = CartesianAxis::z;
      ++result.index.z;
      break;
  }
  return result;
}

bool positive_face(ImmersedFaceDirection direction) noexcept {
  return direction == ImmersedFaceDirection::x_positive ||
         direction == ImmersedFaceDirection::y_positive ||
         direction == ImmersedFaceDirection::z_positive;
}

Int3 offset(Int3 value, CartesianAxis axis, std::int32_t amount) noexcept {
  if (axis == CartesianAxis::x)
    value.x += amount;
  else if (axis == CartesianAxis::y)
    value.y += amount;
  else
    value.z += amount;
  return value;
}

std::int32_t normal_index(Int3 face, CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? face.x
             : (axis == CartesianAxis::y ? face.y : face.z);
}

double regular_cross_traction(const CartesianKernelPlan& kernels,
                              ConstFieldView gradient,
                              ConstFieldView viscosity, InterfaceFace face,
                              std::uint8_t momentum_component) noexcept {
  const Int3 left = offset(face.index, face.axis, -1);
  const double mu_face = detail::interpolate_face(
      kernels, face.axis, normal_index(face.index, face.axis),
      viscosity.unchecked(left, 0U),
      viscosity.unchecked(face.index, 0U));
  const auto interpolate_gradient = [&](std::uint8_t component) {
    return detail::interpolate_face(
        kernels, face.axis, normal_index(face.index, face.axis),
        gradient.unchecked(left, component),
        gradient.unchecked(face.index, component));
  };
  const double divergence = interpolate_gradient(0U) +
                            interpolate_gradient(4U) +
                            interpolate_gradient(8U);
  const std::uint8_t axis = static_cast<std::uint8_t>(face.axis);
  const double transpose = interpolate_gradient(
      static_cast<std::uint8_t>(3U * axis + momentum_component));
  return mu_face *
         (transpose -
          (axis == momentum_component ? (2.0 / 3.0) * divergence : 0.0));
}

FaceFieldView select(FaceFluxView flux, CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? flux.x
             : (axis == CartesianAxis::y ? flux.y : flux.z);
}

ConstFaceFieldView select(ConstFaceFluxView flux,
                          CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? flux.x
             : (axis == CartesianAxis::y ? flux.y : flux.z);
}

bool valid_plan_inputs(const CartesianKernelPlan& kernels,
                       const EBTopology& topology,
                       const BoundaryStencilPlan& boundary,
                       const IbmInterfaceMetricPlan& metric) noexcept {
  const Int3 cells = kernels.cells();
  const std::size_t cell_count =
      cells.x > 0 && cells.y > 0 && cells.z > 0
          ? static_cast<std::size_t>(cells.x) * cells.y * cells.z
          : 0U;
  return kernels.fingerprint() != 0U && topology.fingerprint() != 0U &&
         boundary.fingerprint() != 0U && metric.fingerprint() != 0U &&
         metric.geometry_revision() == topology.geometry_revision() &&
         metric.geometry_fingerprint() == topology.geometry_fingerprint() &&
         metric.surface_fingerprint() == topology.surface_fingerprint() &&
         metric.physical_fingerprint() ==
             topology.interface_metric().physical_fingerprint() &&
         metric.fingerprint() == topology.interface_metric().fingerprint() &&
         cell_count != 0U &&
         topology.region().size == cell_count &&
         topology.links().size == boundary.links().size &&
         topology.links().size == metric.links().size &&
         boundary.reconstruction().fingerprint() != 0U;
}

Status validate_bound(const IbmEquationInterfacePlan& plan,
                      const CartesianKernelPlan* kernels,
                      const EBTopology* topology,
                      const BoundaryStencilPlan* boundary,
                      const IbmInterfaceMetricPlan* metric) noexcept {
  return plan.fingerprint() != 0U && kernels != nullptr &&
                 topology != nullptr && boundary != nullptr &&
                 metric != nullptr && metric->fingerprint() != 0U
             ? Status{}
             : Status{StatusCode::invalid_plan, kIbmEquationApply};
}

}  // namespace

struct IbmPhysicalBoundaryFluxAuthority::Impl {
  ~Impl() noexcept {
    if (communicator != MPI_COMM_NULL && mpi_live())
      MPI_Comm_free(&communicator);
  }

  MPI_Comm communicator{MPI_COMM_NULL};
  const CartesianGeometryPlan* geometry{};
  const EBTopology* topology{};
  const IbmEquationInterfacePlan* immersed_interface{};
  MeshPatch patch{};
  Int3 cells{};
  std::array<bool, 6U> local_physical{};
  std::array<std::vector<std::uint8_t>, 6U> activity;
  PlanFingerprint geometry_fingerprint{};
  PlanFingerprint topology_fingerprint{};
  PlanFingerprint interface_fingerprint{};
  PlanFingerprint local_fingerprint{};
  PlanFingerprint collective_fingerprint{};
};

IbmPhysicalBoundaryFluxAuthority::~IbmPhysicalBoundaryFluxAuthority()
    noexcept {
  release();
}

IbmPhysicalBoundaryFluxAuthority::IbmPhysicalBoundaryFluxAuthority(
    IbmPhysicalBoundaryFluxAuthority&& other) noexcept
    : implementation_(std::exchange(other.implementation_, nullptr)) {}

IbmPhysicalBoundaryFluxAuthority&
IbmPhysicalBoundaryFluxAuthority::operator=(
    IbmPhysicalBoundaryFluxAuthority&& other) noexcept {
  if (this != &other) {
    release();
    implementation_ = std::exchange(other.implementation_, nullptr);
  }
  return *this;
}

void IbmPhysicalBoundaryFluxAuthority::release() noexcept {
  delete std::exchange(implementation_, nullptr);
}

Status IbmPhysicalBoundaryFluxAuthority::compile(
    MPI_Comm communicator, const CartesianGeometryPlan& geometry,
    MeshPatch patch, const EBTopology& topology,
    const IbmEquationInterfacePlan& immersed_interface,
    IbmPhysicalBoundaryFluxAuthority& out) noexcept {
  int rank = -1;
  int size = 0;
  if (!mpi_live() || communicator == MPI_COMM_NULL ||
      MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || rank < 0 ||
      size <= 0)
    return {StatusCode::invalid_plan, kIbmPhysicalBoundaryPlan};

  const auto multiply = [](std::size_t left, std::size_t right,
                           std::size_t& value) noexcept {
    if (left != 0U &&
        right > std::numeric_limits<std::size_t>::max() / left)
      return false;
    value = left * right;
    return true;
  };
  const Int3 global = geometry.global_cells();
  MeshPatch canonical_patch;
  const Status canonical_patch_status =
      detail::make_mesh_patch(rank, size, global, canonical_patch);
  const bool positive_shape =
      patch.cells.x > 0 && patch.cells.y > 0 && patch.cells.z > 0 &&
      global.x > 0 && global.y > 0 && global.z > 0;
  const bool patch_in_domain =
      positive_shape && patch.begin.x >= 0 && patch.begin.y >= 0 &&
      patch.begin.z >= 0 && patch.begin.x <= global.x - patch.cells.x &&
      patch.begin.y <= global.y - patch.cells.y &&
      patch.begin.z <= global.z - patch.cells.z;
  std::size_t xy = 0U;
  std::size_t cell_count = 0U;
  std::size_t x_plane = 0U;
  std::size_t y_plane = 0U;
  std::size_t z_plane = 0U;
  bool local_valid =
      patch_in_domain && canonical_patch_status &&
      same_patch(patch, canonical_patch) &&
      same_patch(topology.patch_, patch) &&
      topology.global_cells_.x == global.x &&
      topology.global_cells_.y == global.y &&
      topology.global_cells_.z == global.z &&
      multiply(static_cast<std::size_t>(patch.cells.x),
               static_cast<std::size_t>(patch.cells.y), xy) &&
      multiply(xy, static_cast<std::size_t>(patch.cells.z), cell_count) &&
      multiply(static_cast<std::size_t>(patch.cells.y),
               static_cast<std::size_t>(patch.cells.z), x_plane) &&
      multiply(static_cast<std::size_t>(patch.cells.x),
               static_cast<std::size_t>(patch.cells.z), y_plane) &&
      multiply(static_cast<std::size_t>(patch.cells.x),
               static_cast<std::size_t>(patch.cells.y), z_plane);
  local_valid =
      local_valid && geometry.fingerprint() != 0U &&
      geometry.topology_revision() != 0U && topology.fingerprint() != 0U &&
      topology.geometry_revision() == geometry.topology_revision() &&
      topology.geometry_fingerprint() == geometry.fingerprint() &&
      topology.surface_fingerprint() != 0U &&
      topology.region().data != nullptr &&
      topology.region().size == cell_count &&
      immersed_interface.fingerprint_ != 0U &&
      immersed_interface.kernels_ != nullptr &&
      immersed_interface.topology_ == &topology &&
      immersed_interface.boundary_ != nullptr &&
      immersed_interface.metric_ != nullptr &&
      immersed_interface.metric_->geometry_revision() ==
          topology.geometry_revision() &&
      immersed_interface.metric_->geometry_fingerprint() ==
          topology.geometry_fingerprint() &&
      immersed_interface.metric_->surface_fingerprint() ==
          topology.surface_fingerprint();
  const Span<const std::uint8_t> region = topology.region();
  for (std::size_t index = 0U; local_valid && index < region.size; ++index)
    local_valid =
        region.data[index] == static_cast<std::uint8_t>(RegionFlag::fluid) ||
        region.data[index] == static_cast<std::uint8_t>(RegionFlag::solid);

  const int local_flag = local_valid ? 1 : 0;
  int global_flag = 0;
  if (MPI_Allreduce(&local_flag, &global_flag, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kIbmPhysicalBoundaryPlan};
  if (global_flag == 0)
    return {StatusCode::invalid_plan, kIbmPhysicalBoundaryPlan};

  const std::array<std::uint64_t, 9U> local_semantics{{
      geometry.fingerprint(),
      geometry.topology_revision(),
      topology.geometry_fingerprint(),
      topology.surface_fingerprint(),
      topology.geometry_revision(),
      topology.fingerprint(),
      static_cast<std::uint64_t>(global.x),
      static_cast<std::uint64_t>(global.y),
      static_cast<std::uint64_t>(global.z),
  }};
  std::array<std::uint64_t, 9U> minimum_semantics{};
  std::array<std::uint64_t, 9U> maximum_semantics{};
  const int minimum_semantics_status = MPI_Allreduce(
      local_semantics.data(), minimum_semantics.data(),
      static_cast<int>(local_semantics.size()), MPI_UINT64_T, MPI_MIN,
      communicator);
  const int maximum_semantics_status = MPI_Allreduce(
      local_semantics.data(), maximum_semantics.data(),
      static_cast<int>(local_semantics.size()), MPI_UINT64_T, MPI_MAX,
      communicator);
  if (minimum_semantics_status != MPI_SUCCESS ||
      maximum_semantics_status != MPI_SUCCESS)
    return {StatusCode::mpi_failure, kIbmPhysicalBoundaryPlan};
  if (minimum_semantics != maximum_semantics)
    return {StatusCode::invalid_plan, kIbmPhysicalBoundaryPlan};

  Impl* candidate = new (std::nothrow) Impl;
  int local_failure = candidate == nullptr ? 1 : 0;
  if (candidate != nullptr) {
    try {
      candidate->activity[0U].assign(x_plane, 1U);
      candidate->activity[1U].assign(x_plane, 1U);
      candidate->activity[2U].assign(y_plane, 1U);
      candidate->activity[3U].assign(y_plane, 1U);
      candidate->activity[4U].assign(z_plane, 1U);
      candidate->activity[5U].assign(z_plane, 1U);
    } catch (...) {
      local_failure = 1;
    }
  }
  int global_failure = 0;
  if (MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX,
                    communicator) != MPI_SUCCESS) {
    delete candidate;
    return {StatusCode::mpi_failure, kIbmPhysicalBoundaryPlan};
  }
  if (global_failure != 0) {
    delete candidate;
    return {StatusCode::allocation_failure, kIbmPhysicalBoundaryPlan};
  }

  candidate->geometry = &geometry;
  candidate->topology = &topology;
  candidate->immersed_interface = &immersed_interface;
  candidate->patch = patch;
  candidate->cells = patch.cells;
  candidate->geometry_fingerprint = geometry.fingerprint();
  candidate->topology_fingerprint = topology.fingerprint();
  candidate->interface_fingerprint = immersed_interface.fingerprint();
  candidate->local_physical = {{
      patch.begin.x == 0,
      patch.begin.x + patch.cells.x == global.x,
      patch.begin.y == 0,
      patch.begin.y + patch.cells.y == global.y,
      patch.begin.z == 0,
      patch.begin.z + patch.cells.z == global.z,
  }};
  const auto cell_offset = [&](Int3 cell) noexcept {
    return static_cast<std::size_t>(cell.x) +
           static_cast<std::size_t>(patch.cells.x) *
               (static_cast<std::size_t>(cell.y) +
                static_cast<std::size_t>(patch.cells.y) *
                    static_cast<std::size_t>(cell.z));
  };
  for (std::int32_t z = 0; z < patch.cells.z; ++z)
    for (std::int32_t y = 0; y < patch.cells.y; ++y) {
      const std::size_t plane =
          static_cast<std::size_t>(y) +
          static_cast<std::size_t>(patch.cells.y) * z;
      if (candidate->local_physical[0U])
        candidate->activity[0U][plane] =
            region.data[cell_offset({0, y, z})];
      if (candidate->local_physical[1U])
        candidate->activity[1U][plane] =
            region.data[cell_offset({patch.cells.x - 1, y, z})];
    }
  for (std::int32_t z = 0; z < patch.cells.z; ++z)
    for (std::int32_t x = 0; x < patch.cells.x; ++x) {
      const std::size_t plane =
          static_cast<std::size_t>(x) +
          static_cast<std::size_t>(patch.cells.x) * z;
      if (candidate->local_physical[2U])
        candidate->activity[2U][plane] =
            region.data[cell_offset({x, 0, z})];
      if (candidate->local_physical[3U])
        candidate->activity[3U][plane] =
            region.data[cell_offset({x, patch.cells.y - 1, z})];
    }
  for (std::int32_t y = 0; y < patch.cells.y; ++y)
    for (std::int32_t x = 0; x < patch.cells.x; ++x) {
      const std::size_t plane =
          static_cast<std::size_t>(x) +
          static_cast<std::size_t>(patch.cells.x) * y;
      if (candidate->local_physical[4U])
        candidate->activity[4U][plane] =
            region.data[cell_offset({x, y, 0})];
      if (candidate->local_physical[5U])
        candidate->activity[5U][plane] =
            region.data[cell_offset({x, y, patch.cells.z - 1})];
    }

  std::uint64_t local = mix(kFnvOffset, UINT64_C(0x69626d7062666163));
  local = mix(local, geometry.fingerprint());
  local = mix(local, topology.fingerprint());
  local = mix(local, immersed_interface.fingerprint());
  for (const Int3 value : {patch.begin, patch.cells, patch.process_grid,
                           patch.process_coord}) {
    local = mix(local, static_cast<std::uint32_t>(value.x));
    local = mix(local, static_cast<std::uint32_t>(value.y));
    local = mix(local, static_cast<std::uint32_t>(value.z));
  }
  for (std::size_t face = 0U; face < candidate->activity.size(); ++face) {
    local = mix(local, candidate->local_physical[face] ? 1U : 0U);
    local = mix(local, candidate->activity[face].size());
    for (const std::uint8_t value : candidate->activity[face])
      local = mix(local, value);
  }
  candidate->local_fingerprint = local == 0U ? 1U : local;

  std::vector<std::uint64_t> gathered;
  try {
    gathered.resize(static_cast<std::size_t>(size));
  } catch (...) {
    local_failure = 1;
  }
  if (MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX,
                    communicator) != MPI_SUCCESS) {
    delete candidate;
    return {StatusCode::mpi_failure, kIbmPhysicalBoundaryPlan};
  }
  if (global_failure != 0) {
    delete candidate;
    return {StatusCode::allocation_failure, kIbmPhysicalBoundaryPlan};
  }
  if (MPI_Allgather(&candidate->local_fingerprint, 1, MPI_UINT64_T,
                    gathered.data(), 1, MPI_UINT64_T,
                    communicator) != MPI_SUCCESS) {
    delete candidate;
    return {StatusCode::mpi_failure, kIbmPhysicalBoundaryPlan};
  }
  std::uint64_t collective =
      mix(kFnvOffset, UINT64_C(0x69626d706266636f));
  for (const std::uint64_t value : local_semantics)
    collective = mix(collective, value);
  collective = mix(collective, static_cast<std::uint64_t>(size));
  for (const std::uint64_t value : gathered)
    collective = mix(collective, value);
  candidate->collective_fingerprint =
      collective == 0U ? 1U : collective;
  if (MPI_Comm_dup(communicator, &candidate->communicator) != MPI_SUCCESS) {
    delete candidate;
    return {StatusCode::mpi_failure, kIbmPhysicalBoundaryPlan};
  }
  out.release();
  out.implementation_ = candidate;
  return {};
}

bool IbmPhysicalBoundaryFluxAuthority::ready() const noexcept {
  return implementation_ != nullptr &&
         implementation_->local_fingerprint != 0U &&
         implementation_->collective_fingerprint != 0U;
}

PlanFingerprint IbmPhysicalBoundaryFluxAuthority::local_fingerprint() const
    noexcept {
  return implementation_ == nullptr ? 0U
                                    : implementation_->local_fingerprint;
}

PlanFingerprint
IbmPhysicalBoundaryFluxAuthority::collective_fingerprint() const noexcept {
  return implementation_ == nullptr ? 0U
                                    : implementation_->collective_fingerprint;
}

bool IbmPhysicalBoundaryFluxAuthority::matches(
    MPI_Comm communicator, const CartesianGeometryPlan* geometry,
    MeshPatch patch,
    const IbmEquationInterfacePlan* immersed_interface) const noexcept {
  const Impl* impl = implementation_;
  int communicator_relation = MPI_UNEQUAL;
  const bool communicator_matches =
      impl != nullptr && mpi_live() && communicator != MPI_COMM_NULL &&
      impl->communicator != MPI_COMM_NULL &&
      MPI_Comm_compare(impl->communicator, communicator,
                       &communicator_relation) == MPI_SUCCESS &&
      (communicator_relation == MPI_IDENT ||
       communicator_relation == MPI_CONGRUENT);
  return communicator_matches && geometry != nullptr &&
         immersed_interface != nullptr && impl->geometry == geometry &&
         impl->immersed_interface == immersed_interface &&
         impl->geometry_fingerprint == geometry->fingerprint() &&
         impl->topology != nullptr &&
         impl->topology_fingerprint == impl->topology->fingerprint() &&
         same_patch(impl->topology->patch_, patch) &&
         immersed_interface->topology_ == impl->topology &&
         impl->interface_fingerprint == immersed_interface->fingerprint() &&
         same_patch(impl->patch, patch) && ready();
}

Status IbmPhysicalBoundaryFluxAuthority::physical_face_active(
    CartesianAxis axis, Int3 face, bool& active) const noexcept {
  active = false;
  const Impl* impl = implementation_;
  if (impl == nullptr || !ready())
    return {StatusCode::invalid_plan, kIbmPhysicalBoundaryApply};
  const std::size_t axis_index = static_cast<std::size_t>(axis);
  if (axis_index >= 3U)
    return {StatusCode::invalid_plan, kIbmPhysicalBoundaryApply};
  const std::int32_t normal = axis == CartesianAxis::x
                                  ? face.x
                                  : (axis == CartesianAxis::y ? face.y
                                                              : face.z);
  const std::int32_t extent = axis == CartesianAxis::x
                                  ? impl->cells.x
                                  : (axis == CartesianAxis::y ? impl->cells.y
                                                              : impl->cells.z);
  if (normal != 0 && normal != extent)
    return {StatusCode::invalid_plan, kIbmPhysicalBoundaryApply};
  const std::size_t slot = 2U * axis_index + (normal == extent ? 1U : 0U);
  if (!impl->local_physical[slot]) {
    active = true;
    return {};
  }
  std::size_t offset = 0U;
  if (axis == CartesianAxis::x) {
    if (face.y < 0 || face.y >= impl->cells.y || face.z < 0 ||
        face.z >= impl->cells.z)
      return {StatusCode::invalid_plan, kIbmPhysicalBoundaryApply};
    offset = static_cast<std::size_t>(face.y) +
             static_cast<std::size_t>(impl->cells.y) * face.z;
  } else if (axis == CartesianAxis::y) {
    if (face.x < 0 || face.x >= impl->cells.x || face.z < 0 ||
        face.z >= impl->cells.z)
      return {StatusCode::invalid_plan, kIbmPhysicalBoundaryApply};
    offset = static_cast<std::size_t>(face.x) +
             static_cast<std::size_t>(impl->cells.x) * face.z;
  } else {
    if (face.x < 0 || face.x >= impl->cells.x || face.y < 0 ||
        face.y >= impl->cells.y)
      return {StatusCode::invalid_plan, kIbmPhysicalBoundaryApply};
    offset = static_cast<std::size_t>(face.x) +
             static_cast<std::size_t>(impl->cells.x) * face.y;
  }
  if (offset >= impl->activity[slot].size())
    return {StatusCode::invalid_plan, kIbmPhysicalBoundaryApply};
  active = impl->activity[slot][offset] != 0U;
  return {};
}

Status IbmPhysicalBoundaryFluxAuthority::zero_inactive_physical_boundary_flux(
    FaceFluxView flux) const noexcept {
  const Impl* impl = implementation_;
  if (impl == nullptr || !ready() ||
      !detail::valid_flux_view(flux, impl->cells))
    return {StatusCode::invalid_plan, kIbmPhysicalBoundaryApply};
  const std::array<FaceFieldView, 3U> faces{flux.x, flux.y, flux.z};
  for (std::size_t axis_index = 0U; axis_index < faces.size(); ++axis_index) {
    const CartesianAxis axis = static_cast<CartesianAxis>(axis_index);
    const std::int32_t extent = axis == CartesianAxis::x
                                    ? impl->cells.x
                                    : (axis == CartesianAxis::y
                                           ? impl->cells.y
                                           : impl->cells.z);
    for (std::size_t side = 0U; side < 2U; ++side) {
      const std::size_t slot = 2U * axis_index + side;
      if (!impl->local_physical[slot]) continue;
      const std::int32_t normal = side == 0U ? 0 : extent;
      const auto zero = [&](Int3 face) noexcept {
        bool active = false;
        const Status status = physical_face_active(axis, face, active);
        if (!status) return status;
        if (!active) faces[axis_index].unchecked(face) = 0.0;
        return Status{};
      };
      if (axis == CartesianAxis::x) {
        for (std::int32_t z = 0; z < impl->cells.z; ++z)
          for (std::int32_t y = 0; y < impl->cells.y; ++y) {
            const Status status = zero({normal, y, z});
            if (!status) return status;
          }
      } else if (axis == CartesianAxis::y) {
        for (std::int32_t z = 0; z < impl->cells.z; ++z)
          for (std::int32_t x = 0; x < impl->cells.x; ++x) {
            const Status status = zero({x, normal, z});
            if (!status) return status;
          }
      } else {
        for (std::int32_t y = 0; y < impl->cells.y; ++y)
          for (std::int32_t x = 0; x < impl->cells.x; ++x) {
            const Status status = zero({x, y, normal});
            if (!status) return status;
          }
      }
    }
  }
  return {};
}

Status
IbmPhysicalBoundaryFluxAuthority::validate_inactive_physical_boundary_flux(
    ConstFaceFluxView flux, double absolute_tolerance) const noexcept {
  const Impl* impl = implementation_;
  if (impl == nullptr || !ready() || !std::isfinite(absolute_tolerance) ||
      absolute_tolerance < 0.0 ||
      !detail::valid_flux_view(flux, impl->cells, flux.revision))
    return {StatusCode::invalid_plan, kIbmPhysicalBoundaryApply};
  const std::array<ConstFaceFieldView, 3U> faces{flux.x, flux.y, flux.z};
  for (std::size_t axis_index = 0U; axis_index < faces.size(); ++axis_index) {
    const CartesianAxis axis = static_cast<CartesianAxis>(axis_index);
    const std::int32_t extent = axis == CartesianAxis::x
                                    ? impl->cells.x
                                    : (axis == CartesianAxis::y
                                           ? impl->cells.y
                                           : impl->cells.z);
    for (std::size_t side = 0U; side < 2U; ++side) {
      const std::size_t slot = 2U * axis_index + side;
      if (!impl->local_physical[slot]) continue;
      const std::int32_t normal = side == 0U ? 0 : extent;
      const auto validate = [&](Int3 face) noexcept {
        bool active = false;
        const Status status = physical_face_active(axis, face, active);
        if (!status) return status;
        if (!active) {
          const double value = faces[axis_index].unchecked(face);
          if (!std::isfinite(value) ||
              std::abs(value) > absolute_tolerance)
            return Status{StatusCode::numerical_failure,
                          kIbmPhysicalBoundaryNumerical};
        }
        return Status{};
      };
      if (axis == CartesianAxis::x) {
        for (std::int32_t z = 0; z < impl->cells.z; ++z)
          for (std::int32_t y = 0; y < impl->cells.y; ++y) {
            const Status status = validate({normal, y, z});
            if (!status) return status;
          }
      } else if (axis == CartesianAxis::y) {
        for (std::int32_t z = 0; z < impl->cells.z; ++z)
          for (std::int32_t x = 0; x < impl->cells.x; ++x) {
            const Status status = validate({x, normal, z});
            if (!status) return status;
          }
      } else {
        for (std::int32_t y = 0; y < impl->cells.y; ++y)
          for (std::int32_t x = 0; x < impl->cells.x; ++x) {
            const Status status = validate({x, y, normal});
            if (!status) return status;
          }
      }
    }
  }
  return {};
}

Status IbmEquationInterfacePlan::compile(
    const CartesianKernelPlan& kernels, const EBTopology& topology,
    const BoundaryStencilPlan& boundary,
    IbmEquationInterfacePlan& out) noexcept {
  return compile(kernels, topology, boundary, topology.interface_metric(),
                 out);
}

Status IbmEquationInterfacePlan::compile(
    const CartesianKernelPlan& kernels, const EBTopology& topology,
    const BoundaryStencilPlan& boundary,
    const IbmInterfaceMetricPlan& metric,
    IbmEquationInterfacePlan& out) noexcept {
  if (!valid_plan_inputs(kernels, topology, boundary, metric))
    return {StatusCode::invalid_plan, kIbmEquationPlan};
  std::uint64_t fingerprint = kFnvOffset;
  fingerprint = mix(fingerprint, kernels.fingerprint());
  fingerprint = mix(fingerprint, topology.fingerprint());
  fingerprint = mix(fingerprint, boundary.fingerprint());
  fingerprint = mix(fingerprint, metric.fingerprint());
  fingerprint = mix(fingerprint, topology.geometry_revision());
  if (fingerprint == 0U) fingerprint = 1U;
  IbmEquationInterfacePlan candidate;
  candidate.kernels_ = &kernels;
  candidate.topology_ = &topology;
  candidate.boundary_ = &boundary;
  candidate.metric_ = &metric;
  candidate.fingerprint_ = fingerprint;
  try {
    const Span<const BoundaryStencilLink> links = boundary.links();
    const Span<const QuadraticAffineRow> rows =
        boundary.reconstruction().rows();
    const Span<const QuadraticStencilGroup> groups =
        boundary.reconstruction().groups();
    const Span<const double> weights = boundary.reconstruction().weights();
    const Span<const GlobalCellId> donor_global_cells =
        boundary.reconstruction().donor_global_cells();
    const Span<const ImmersedLink> topology_links = topology.links();
    const Span<const IbmInterfaceLinkMetric> physical_links = metric.links();
    candidate.wall_linearization_.reserve(links.size);
    for (std::size_t index = 0U; index < links.size; ++index) {
      const std::uint32_t row_index =
          links.data[index].wall_normal_gradient_row;
      if (row_index >= rows.size || rows.data[row_index].group >= groups.size ||
          links.data[index].topology_link >= topology_links.size)
        return {StatusCode::invalid_plan, kIbmEquationPlan};
      const QuadraticAffineRow& row = rows.data[row_index];
      const QuadraticStencilGroup& group = groups.data[row.group];
      const std::size_t count = group.quality.donor_count;
      if (count == 0U || row.weight_begin > weights.size ||
          count > weights.size - row.weight_begin ||
          group.donor_begin > donor_global_cells.size ||
          count > donor_global_cells.size - group.donor_begin)
        return {StatusCode::invalid_plan, kIbmEquationPlan};
      const ImmersedLink& link =
          topology_links.data[links.data[index].topology_link];
      const IbmInterfaceLinkMetric& physical =
          physical_links.data[links.data[index].topology_link];
      if (physical.global_link != link.global_link ||
          physical.source_triangle != link.triangle ||
          !std::isfinite(link.cartesian_control_face_area) ||
          !(link.cartesian_control_face_area > 0.0) ||
          !std::isfinite(physical.physical_quadrature_area) ||
          !(physical.physical_quadrature_area > 0.0))
        return {StatusCode::invalid_plan, kIbmEquationPlan};
      for (const double value : physical.normal_second_moment)
        if (!std::isfinite(value))
          return {StatusCode::invalid_plan, kIbmEquationPlan};
      const Int3 cell = link.fluid_local_index;
      const double wall_dx =
          detail::centre_coordinate(kernels, CartesianAxis::x, cell.x) -
          link.wall_point.x;
      const double wall_dy =
          detail::centre_coordinate(kernels, CartesianAxis::y, cell.y) -
          link.wall_point.y;
      const double wall_dz =
          detail::centre_coordinate(kernels, CartesianAxis::z, cell.z) -
          link.wall_point.z;
      const double wall_distance =
          std::sqrt(wall_dx * wall_dx + wall_dy * wall_dy +
                    wall_dz * wall_dz);
      if (!std::isfinite(wall_distance) || wall_distance <= 0.0)
        return {StatusCode::invalid_plan, kIbmEquationPlan};
      const double inverse_wall_distance = 1.0 / wall_distance;
      std::size_t owner_count = 0U;
      double correction_l1 = 0.0;
      for (std::size_t donor = 0U; donor < count; ++donor) {
        const bool owner =
            donor_global_cells.data[group.donor_begin + donor] ==
            link.fluid_cell;
        correction_l1 +=
            std::abs(weights.data[row.weight_begin + donor] -
                     (owner ? inverse_wall_distance : 0.0));
        owner_count += owner ? 1U : 0U;
      }
      if (owner_count > 1U)
        return {StatusCode::invalid_plan, kIbmEquationPlan};
      if (owner_count == 0U) correction_l1 += inverse_wall_distance;
      const double gradient_majorant =
          inverse_wall_distance + correction_l1;
      const InterfaceFace face = interface_face(link);
      const std::uint8_t axis = static_cast<std::uint8_t>(face.axis);
      const std::int32_t normal =
          axis == 0U ? cell.x : (axis == 1U ? cell.y : cell.z);
      const detail::DerivativeWeights derivative =
          kernels.geometry_kind() == GeometryKind::uniform
              ? detail::metric_derivative_weights<true>(kernels, axis,
                                                        normal)
              : detail::metric_derivative_weights<false>(kernels, axis,
                                                         normal);
      const double solid_pressure_derivative_weight =
          positive_face(link.direction) ? derivative.plus : derivative.minus;
      if (!std::isfinite(gradient_majorant) || gradient_majorant <= 0.0 ||
          !std::isfinite(solid_pressure_derivative_weight) ||
          solid_pressure_derivative_weight == 0.0)
        return {StatusCode::invalid_plan, kIbmEquationPlan};
      candidate.wall_linearization_.push_back(
          {wall_distance, gradient_majorant,
           solid_pressure_derivative_weight});
    }
    out = std::move(candidate);
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, kIbmEquationPlan};
  } catch (...) {
    return {StatusCode::invalid_plan, kIbmEquationPlan};
  }
}

Status IbmEquationInterfacePlan::zero_interface_flux(
    FaceFluxView flux) const noexcept {
  Status status =
      validate_bound(*this, kernels_, topology_, boundary_, metric_);
  if (!status || !detail::valid_flux_view(flux, kernels_->cells()))
    return status ? Status{StatusCode::invalid_plan, kIbmEquationApply}
                  : status;
  const Span<const ImmersedLink> links = topology_->links();
  for (std::size_t index = 0U; index < links.size; ++index) {
    const InterfaceFace face = interface_face(links.data[index]);
    select(flux, face.axis).unchecked(face.index) = 0.0;
  }
  return {};
}

Status IbmEquationInterfacePlan::validate_interface_flux(
    ConstFaceFluxView flux, double absolute_tolerance) const noexcept {
  Status status =
      validate_bound(*this, kernels_, topology_, boundary_, metric_);
  if (!status || !std::isfinite(absolute_tolerance) ||
      absolute_tolerance < 0.0 ||
      !detail::valid_flux_view(flux, kernels_->cells(), flux.revision))
    return status ? Status{StatusCode::invalid_plan, kIbmEquationApply}
                  : status;
  const Span<const ImmersedLink> links = topology_->links();
  for (std::size_t index = 0U; index < links.size; ++index) {
    const InterfaceFace face = interface_face(links.data[index]);
    const double value = select(flux, face.axis).unchecked(face.index);
    if (!std::isfinite(value) || std::abs(value) > absolute_tolerance)
      return {StatusCode::numerical_failure, kIbmEquationNumerical};
  }
  return {};
}

Status IbmEquationInterfacePlan::constrain_pressure_predictor(
    FieldView h_by_a, FaceFluxView phi_h_by_a) const noexcept {
  Status status =
      validate_bound(*this, kernels_, topology_, boundary_, metric_);
  const Int3 cells = kernels_->cells();
  if (!status || !detail::valid_cell_view(h_by_a, cells, 0U, 3U) ||
      !detail::valid_flux_view(phi_h_by_a, cells))
    return status ? Status{StatusCode::invalid_plan, kIbmEquationApply}
                  : status;
  const Span<const std::uint8_t> region = topology_->region();
  std::size_t flat = 0U;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x, ++flat)
        if (region.data[flat] == static_cast<std::uint8_t>(RegionFlag::solid))
          for (std::uint8_t component = 0U; component < 3U; ++component)
            h_by_a.unchecked({x, y, z}, component) = 0.0;
  return zero_interface_flux(phi_h_by_a);
}

Status IbmEquationInterfacePlan::constrain_corrected_state(
    FieldView velocity, FaceFluxView flux) const noexcept {
  return constrain_pressure_predictor(velocity, flux);
}

Status IbmEquationInterfacePlan::constrain_momentum(
    ConstFieldView velocity, ConstFieldView velocity_gradient,
    ConstFieldView pressure_perturbation, ConstFieldView density,
    ConstFieldView molecular_viscosity, ConstFieldView effective_viscosity,
    const TurbulencePlan* wall_treatment,
    IbmCellEquationView system) const noexcept {
  Status status =
      validate_bound(*this, kernels_, topology_, boundary_, metric_);
  const Int3 cells = kernels_->cells();
  if (!status ||
      !detail::valid_cell_view(velocity, cells, 0U, 3U,
                              boundary_->maximum_halo_reach()) ||
      !detail::valid_cell_view(velocity_gradient, cells, 0U, 9U, 1U) ||
      !detail::valid_cell_view(pressure_perturbation, cells, 0U, 1U,
                              boundary_->maximum_halo_reach()) ||
      !detail::valid_cell_view(density, cells, 0U, 1U, 0U) ||
      !detail::valid_cell_view(molecular_viscosity, cells, 0U, 1U, 0U) ||
      !detail::valid_cell_view(effective_viscosity, cells, 0U, 1U,
                              boundary_->maximum_halo_reach()) ||
      !detail::valid_cell_view(system.diagonal, cells, 0U, 3U) ||
      !detail::valid_cell_view(system.rhs, cells, 0U, 3U) ||
      !detail::valid_cell_view(system.residual, cells, 0U, 3U))
    return status ? Status{StatusCode::invalid_plan, kIbmEquationApply}
                  : status;
  const Span<const ImmersedLink> links = topology_->links();
  const Span<const IbmInterfaceLinkMetric> physical_links = metric_->links();
  const Span<const BoundaryStencilLink> rows = boundary_->links();
  if (wall_linearization_.size() != rows.size ||
      physical_links.size != links.size)
    return {StatusCode::invalid_plan, kIbmEquationApply};
  for (std::size_t index = 0U; index < rows.size; ++index) {
    const BoundaryStencilLink& row = rows.data[index];
    if (row.topology_link >= links.size)
      return {StatusCode::invalid_plan, kIbmEquationApply};
    const ImmersedLink& link = links.data[row.topology_link];
    const IbmInterfaceLinkMetric& physical =
        physical_links.data[row.topology_link];
    const InterfaceFace face = interface_face(link);
    const double transmissibility = detail::positive_transmissibility(
        *kernels_, effective_viscosity, face.axis, face.index);
    double wall_viscosity = 0.0;
    status = evaluate_positive_bounded_quadratic_row(
        boundary_->reconstruction(), row.wall_value_row,
        effective_viscosity, 0U, wall_viscosity);
    const double area = physical.physical_quadrature_area;
    const Int3 cell = link.fluid_local_index;
    const WallLinearization linearization = wall_linearization_[index];
    double pressure_ghost = 0.0;
    if (status)
      status = evaluate_quadratic_row(
          boundary_->reconstruction(), row.zero_normal_value_row,
          pressure_perturbation, 0U, 0.0, 0.0, pressure_ghost);
    const double solid_pressure =
        pressure_perturbation.unchecked(link.solid_local_index, 0U);
    const double pressure_gradient_correction =
        linearization.solid_pressure_derivative_weight *
        (pressure_ghost - solid_pressure);
    const double pressure_force_correction =
        detail::cell_volume(*kernels_, cell) * pressure_gradient_correction;
    if (!status || !std::isfinite(transmissibility) ||
        transmissibility <= 0.0 || !std::isfinite(wall_viscosity) ||
        wall_viscosity <= 0.0 || !std::isfinite(area) || area <= 0.0 ||
        !std::isfinite(linearization.distance) ||
        linearization.distance <= 0.0 ||
        !std::isfinite(linearization.gradient_majorant) ||
        linearization.gradient_majorant <= 0.0 ||
        !std::isfinite(linearization.solid_pressure_derivative_weight) ||
        !std::isfinite(pressure_ghost) ||
        !std::isfinite(solid_pressure) ||
        !std::isfinite(pressure_gradient_correction) ||
        !std::isfinite(pressure_force_correction))
      return {StatusCode::numerical_failure, kIbmEquationNumerical};
    const double resolved_wall_transmissibility =
        wall_viscosity * area * linearization.gradient_majorant;
    if (!std::isfinite(resolved_wall_transmissibility) ||
        resolved_wall_transmissibility <= 0.0)
      return {StatusCode::numerical_failure, kIbmEquationNumerical};
    double normal_derivative[3]{};
    for (std::uint8_t component = 0U; component < 3U; ++component) {
      status = evaluate_quadratic_row(
          boundary_->reconstruction(), row.wall_normal_gradient_row,
          velocity, component, 0.0, 0.0,
          normal_derivative[component]);
      if (!status || !std::isfinite(normal_derivative[component]))
        return status ? Status{StatusCode::numerical_failure,
                               kIbmEquationNumerical}
                      : status;
    }
    const double normal_components[3]{link.solid_to_fluid_normal.x,
                                      link.solid_to_fluid_normal.y,
                                      link.solid_to_fluid_normal.z};
    const double normal_l1 =
        std::abs(normal_components[0U]) +
        std::abs(normal_components[1U]) +
        std::abs(normal_components[2U]);
    double desired[3]{};
    for (std::uint8_t component = 0U; component < 3U; ++component) {
      double normal_stress = 0.0;
      for (std::uint8_t derivative = 0U; derivative < 3U; ++derivative)
        normal_stress +=
            physical.normal_second_moment[3U * component + derivative] *
            normal_derivative[derivative];
      desired[component] =
          wall_viscosity *
          (area * normal_derivative[component] +
           (1.0 / 3.0) * normal_stress);
    }
    bool equilibrium_wall = false;
    double wall_drag_coefficient = 0.0;
    if (wall_treatment != nullptr &&
        wall_treatment->wall_treatment() ==
            WallTreatmentKind::equilibrium_wall_function) {
      WallFunctionSample sample;
      sample.surface = WallSurfaceKind::immersed;
      sample.solid_to_fluid_normal = link.solid_to_fluid_normal;
      sample.wall_distance = linearization.distance;
      sample.fluid_velocity = {velocity.unchecked(cell, 0U),
                               velocity.unchecked(cell, 1U),
                               velocity.unchecked(cell, 2U)};
      sample.wall_velocity = {};
      sample.density = density.unchecked(cell, 0U);
      sample.molecular_viscosity =
          molecular_viscosity.unchecked(cell, 0U);
      // Momentum only consumes shear. Finite neutral thermal/scalar values
      // satisfy the shared wall-law sample contract without creating a
      // second transport closure here.
      sample.heat_capacity = 1.0;
      sample.molecular_conductivity = 0.0;
      sample.fluid_temperature = 1.0;
      sample.wall_temperature = 1.0;
      sample.molecular_mass_diffusivity = 0.0;
      sample.fluid_scalar = 0.0;
      sample.wall_scalar = 0.0;
      WallFunctionResult wall;
      status = wall_treatment->evaluate_wall_function(sample, wall);
      if (!status) return status;
      const double normal_speed =
          sample.fluid_velocity.x * normal_components[0U] +
          sample.fluid_velocity.y * normal_components[1U] +
          sample.fluid_velocity.z * normal_components[2U];
      const double tangential_x =
          sample.fluid_velocity.x - normal_speed * normal_components[0U];
      const double tangential_y =
          sample.fluid_velocity.y - normal_speed * normal_components[1U];
      const double tangential_z =
          sample.fluid_velocity.z - normal_speed * normal_components[2U];
      const double tangential_speed = std::sqrt(
          tangential_x * tangential_x + tangential_y * tangential_y +
          tangential_z * tangential_z);
      wall_drag_coefficient =
          tangential_speed > 0.0
              ? sample.density * wall.friction_velocity *
                    wall.friction_velocity * area / tangential_speed
              : 0.0;
      if (!std::isfinite(wall_drag_coefficient) ||
          wall_drag_coefficient < 0.0)
        return {StatusCode::numerical_failure, kIbmEquationNumerical};
      equilibrium_wall = true;
      const double resolved_normal = desired[0U] * normal_components[0U] +
                                     desired[1U] * normal_components[1U] +
                                     desired[2U] * normal_components[2U];
      desired[0U] = resolved_normal * normal_components[0U] -
                    wall.shear_on_fluid.x * area;
      desired[1U] = resolved_normal * normal_components[1U] -
                    wall.shear_on_fluid.y * area;
      desired[2U] = resolved_normal * normal_components[2U] -
                    wall.shear_on_fluid.z * area;
    }
    for (std::uint8_t component = 0U; component < 3U; ++component) {
      const double fluid =
          velocity.unchecked(link.fluid_local_index, component);
      const double solid =
          velocity.unchecked(link.solid_local_index, component);
      const double regular_laplacian =
          transmissibility * (fluid - solid);
      const double cartesian_cross = regular_cross_traction(
          *kernels_, velocity_gradient, effective_viscosity, face,
          component);
      const double regular_cross_residual =
          (positive_face(link.direction) ? -1.0 : 1.0) *
          cartesian_cross * link.cartesian_control_face_area;
      // The active IBM pressure graph deletes this fluid--solid Cartesian
      // face. Replace the ordinary solid value by the same quadratic
      // zero-normal ghost consumed by both PISO corrections.
      const double pressure_correction =
          component == static_cast<std::uint8_t>(face.axis)
              ? pressure_force_correction
              : 0.0;
      const double correction = desired[component] - regular_laplacian -
                                regular_cross_residual + pressure_correction;
      const double normal_squared =
          normal_components[component] * normal_components[component];
      double physical_normal_row_l1 = 0.0;
      for (std::uint8_t derivative = 0U; derivative < 3U; ++derivative)
        physical_normal_row_l1 += std::abs(
            physical.normal_second_moment[3U * component + derivative]);
      physical_normal_row_l1 /= area;
      const double tangential_row_l1 =
          1.0 - normal_squared +
          std::abs(normal_components[component]) *
              (normal_l1 - std::abs(normal_components[component]));
      // The exact quadratic traction remains in the residual.  Replace the
      // removed Cartesian fluid-solid face in the HbyA diagonal with a
      // positive geometric low-order operator plus an L1 majorant of the
      // exact-quadratic-minus-low-order correction.  Use the absolute row sum
      // of the traction tensor so cross-component n_i*n_j terms are covered
      // by the scalar component solve; for the wall-law path, also cover the
      // full tangential projector row.  This changes only the residual
      // splitting, not the assembled equation at the trial state.
      const double replacement_diagonal =
          equilibrium_wall
              ? (4.0 / 3.0) * resolved_wall_transmissibility *
                        physical_normal_row_l1 +
                    wall_drag_coefficient * tangential_row_l1
              : resolved_wall_transmissibility *
                    (1.0 + (1.0 / 3.0) * physical_normal_row_l1);
      const double diagonal =
          system.diagonal.unchecked(cell, component) +
          replacement_diagonal - transmissibility;
      if (!std::isfinite(fluid) || !std::isfinite(solid) ||
          !std::isfinite(cartesian_cross) ||
          !std::isfinite(desired[component]) ||
          !std::isfinite(correction) ||
          !std::isfinite(replacement_diagonal) ||
          replacement_diagonal < 0.0 || !std::isfinite(diagonal) ||
          diagonal <= 0.0)
        return {StatusCode::numerical_failure, kIbmEquationNumerical};
      const double residual =
          system.residual.unchecked(cell, component) + correction;
      const double rhs = diagonal * fluid - residual;
      if (!std::isfinite(residual) || !std::isfinite(rhs))
        return {StatusCode::numerical_failure, kIbmEquationNumerical};
      system.diagonal.unchecked(cell, component) = diagonal;
      system.residual.unchecked(cell, component) = residual;
      system.rhs.unchecked(cell, component) = rhs;
    }
  }
  const Span<const std::uint8_t> region = topology_->region();
  std::size_t flat = 0U;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x, ++flat)
        if (region.data[flat] == static_cast<std::uint8_t>(RegionFlag::solid))
          for (std::uint8_t component = 0U; component < 3U; ++component) {
            const Int3 cell{x, y, z};
            system.diagonal.unchecked(cell, component) = 1.0;
            system.rhs.unchecked(cell, component) = 0.0;
            system.residual.unchecked(cell, component) =
                velocity.unchecked(cell, component);
          }
  return {};
}

Status IbmEquationInterfacePlan::correct_pressure_gradient(
    ConstFieldView pressure, FieldView gradient) const noexcept {
  Status status =
      validate_bound(*this, kernels_, topology_, boundary_, metric_);
  const Int3 cells = kernels_->cells();
  if (!status ||
      !detail::valid_cell_view(pressure, cells, 0U, 1U,
                              boundary_->maximum_halo_reach()) ||
      !detail::valid_cell_view(gradient, cells, 0U, 3U) ||
      detail::field_views_overlap(pressure, as_const(gradient)))
    return status ? Status{StatusCode::invalid_plan, kIbmEquationApply}
                  : status;
  const Span<const ImmersedLink> links = topology_->links();
  const Span<const BoundaryStencilLink> rows = boundary_->links();
  if (wall_linearization_.size() != rows.size)
    return {StatusCode::invalid_plan, kIbmEquationApply};
  for (std::size_t index = 0U; index < rows.size; ++index) {
    const BoundaryStencilLink& row = rows.data[index];
    if (row.topology_link >= links.size)
      return {StatusCode::invalid_plan, kIbmEquationApply};
    const ImmersedLink& link = links.data[row.topology_link];
    const WallLinearization linearization = wall_linearization_[index];
    const InterfaceFace face = interface_face(link);
    const std::uint8_t component = static_cast<std::uint8_t>(face.axis);
    double ghost = 0.0;
    status = evaluate_quadratic_row(
        boundary_->reconstruction(), row.zero_normal_value_row, pressure, 0U,
        0.0, 0.0, ghost);
    const double solid =
        pressure.unchecked(link.solid_local_index, 0U);
    const double correction =
        linearization.solid_pressure_derivative_weight *
        (ghost - solid);
    const double value =
        gradient.unchecked(link.fluid_local_index, component) + correction;
    if (!status ||
        !std::isfinite(linearization.solid_pressure_derivative_weight) ||
        !std::isfinite(ghost) || !std::isfinite(solid) ||
        !std::isfinite(correction) || !std::isfinite(value))
      return status ? Status{StatusCode::numerical_failure,
                             kIbmEquationNumerical}
                    : status;
    gradient.unchecked(link.fluid_local_index, component) = value;
  }
  return {};
}

Status IbmEquationInterfacePlan::correct_pressure_work(
    ConstFieldView pressure, ConstFieldView velocity, FieldView rate) const
    noexcept {
  Status status =
      validate_bound(*this, kernels_, topology_, boundary_, metric_);
  const Int3 cells = kernels_->cells();
  if (!status ||
      !detail::valid_cell_view(pressure, cells, 0U, 1U,
                               boundary_->maximum_halo_reach()) ||
      !detail::valid_cell_view(velocity, cells, 0U, 3U, 0U) ||
      !detail::valid_cell_view(rate, cells, 0U, 1U) ||
      detail::field_views_overlap(pressure, velocity) ||
      detail::field_views_overlap(pressure, as_const(rate)) ||
      detail::field_views_overlap(velocity, as_const(rate)))
    return status ? Status{StatusCode::invalid_plan, kIbmEquationApply}
                  : status;
  const Span<const ImmersedLink> links = topology_->links();
  const Span<const BoundaryStencilLink> rows = boundary_->links();
  if (wall_linearization_.size() != rows.size)
    return {StatusCode::invalid_plan, kIbmEquationApply};
  for (std::size_t index = 0U; index < rows.size; ++index) {
    const BoundaryStencilLink& row = rows.data[index];
    if (row.topology_link >= links.size)
      return {StatusCode::invalid_plan, kIbmEquationApply};
    const ImmersedLink& link = links.data[row.topology_link];
    const InterfaceFace face = interface_face(link);
    const std::uint8_t component = static_cast<std::uint8_t>(face.axis);
    const WallLinearization linearization = wall_linearization_[index];
    double ghost = 0.0;
    status = evaluate_quadratic_row(
        boundary_->reconstruction(), row.zero_normal_value_row, pressure, 0U,
        0.0, 0.0, ghost);
    const double solid = pressure.unchecked(link.solid_local_index, 0U);
    const double pressure_gradient_correction =
        linearization.solid_pressure_derivative_weight * (ghost - solid);
    const double velocity_normal =
        velocity.unchecked(link.fluid_local_index, component);
    const double correction = velocity_normal * pressure_gradient_correction;
    const double value = rate.unchecked(link.fluid_local_index, 0U) + correction;
    if (!status ||
        !std::isfinite(linearization.solid_pressure_derivative_weight) ||
        !std::isfinite(ghost) || !std::isfinite(solid) ||
        !std::isfinite(pressure_gradient_correction) ||
        !std::isfinite(velocity_normal) || !std::isfinite(correction) ||
        !std::isfinite(value))
      return status ? Status{StatusCode::numerical_failure,
                             kIbmEquationNumerical}
                    : status;
    rate.unchecked(link.fluid_local_index, 0U) = value;
  }
  return {};
}

Status IbmEquationInterfacePlan::correct_velocity_gradient(
    ConstFieldView velocity, FieldView velocity_gradient,
    Real3 wall_velocity) const noexcept {
  Status status =
      validate_bound(*this, kernels_, topology_, boundary_, metric_);
  if (!status) return status;
  const Int3 cells = kernels_->cells();
  if (!detail::valid_cell_view(velocity, cells, 0U, 3U,
                              boundary_->maximum_halo_reach()) ||
      !detail::valid_cell_view(velocity_gradient, cells, 0U, 9U) ||
      !std::isfinite(wall_velocity.x) || !std::isfinite(wall_velocity.y) ||
      !std::isfinite(wall_velocity.z))
    return {StatusCode::invalid_plan, kIbmEquationApply};
  const double wall_components[3U]{wall_velocity.x, wall_velocity.y,
                                   wall_velocity.z};
  const Span<const ImmersedLink> links = topology_->links();
  const Span<const BoundaryStencilLink> rows = boundary_->links();
  if (wall_linearization_.size() != rows.size)
    return {StatusCode::invalid_plan, kIbmEquationApply};
  for (std::size_t index = 0U; index < rows.size; ++index) {
    const BoundaryStencilLink& row = rows.data[index];
    if (row.topology_link >= links.size)
      return {StatusCode::invalid_plan, kIbmEquationApply};
    const ImmersedLink& link = links.data[row.topology_link];
    const InterfaceFace face = interface_face(link);
    const std::uint8_t derivative = static_cast<std::uint8_t>(face.axis);
    const double donor_weight =
        wall_linearization_[index].solid_pressure_derivative_weight;
    for (std::uint8_t component = 0U; component < 3U; ++component) {
      double ghost = 0.0;
      status = evaluate_quadratic_row(
          boundary_->reconstruction(), row.dirichlet_value_row, velocity,
          component, wall_components[component], 0.0, ghost);
      const double solid =
          velocity.unchecked(link.solid_local_index, component);
      const std::uint8_t gradient_component =
          static_cast<std::uint8_t>(3U * component + derivative);
      const double corrected =
          velocity_gradient.unchecked(link.fluid_local_index,
                                      gradient_component) +
          donor_weight * (ghost - solid);
      if (!status || !std::isfinite(donor_weight) ||
          !std::isfinite(ghost) || !std::isfinite(solid) ||
          !std::isfinite(corrected))
        return status ? Status{StatusCode::numerical_failure,
                               kIbmEquationNumerical}
                      : status;
      velocity_gradient.unchecked(link.fluid_local_index,
                                  gradient_component) = corrected;
    }
  }
  return {};
}

Status IbmEquationInterfacePlan::correct_zero_normal_diffusion(
    ConstFieldView transported, ConstFieldView diffusivity,
    FieldView rate) const noexcept {
  Status status =
      validate_bound(*this, kernels_, topology_, boundary_, metric_);
  const Int3 cells = kernels_->cells();
  if (!status ||
      !detail::valid_cell_view(transported, cells, 0U, 1U,
                              boundary_->maximum_halo_reach()) ||
      !detail::valid_cell_view(diffusivity, cells, 0U, 1U, 1U) ||
      !detail::valid_cell_view(rate, cells, 0U, 1U))
    return status ? Status{StatusCode::invalid_plan, kIbmEquationApply}
                  : status;
  const Span<const ImmersedLink> links = topology_->links();
  const Span<const BoundaryStencilLink> rows = boundary_->links();
  for (std::size_t index = 0U; index < rows.size; ++index) {
    const BoundaryStencilLink& row = rows.data[index];
    if (row.topology_link >= links.size)
      return {StatusCode::invalid_plan, kIbmEquationApply};
    const ImmersedLink& link = links.data[row.topology_link];
    const InterfaceFace face = interface_face(link);
    const double transmissibility = detail::positive_transmissibility(
        *kernels_, diffusivity, face.axis, face.index);
    const double volume = detail::cell_volume(*kernels_, link.fluid_local_index);
    double ghost = 0.0;
    status = evaluate_quadratic_row(
        boundary_->reconstruction(), row.zero_normal_value_row, transported,
        0U, 0.0, 0.0, ghost);
    const double solid = transported.unchecked(link.solid_local_index, 0U);
    const double correction =
        transmissibility * (ghost - solid) / volume;
    const double value = rate.unchecked(link.fluid_local_index, 0U) + correction;
    if (!status || !std::isfinite(transmissibility) ||
        transmissibility <= 0.0 || !std::isfinite(volume) || volume <= 0.0 ||
        !std::isfinite(ghost) || !std::isfinite(solid) ||
        !std::isfinite(value))
      return status ? Status{StatusCode::numerical_failure,
                             kIbmEquationNumerical}
                    : status;
    rate.unchecked(link.fluid_local_index, 0U) = value;
  }
  const Span<const std::uint8_t> region = topology_->region();
  std::size_t flat = 0U;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x, ++flat)
        if (region.data[flat] == static_cast<std::uint8_t>(RegionFlag::solid))
          rate.unchecked({x, y, z}, 0U) = 0.0;
  return {};
}

Status IbmEquationInterfacePlan::correct_positive_bounded_zero_normal_diffusion(
    ConstFieldView transported, ConstFieldView diffusivity,
    FieldView rate) const noexcept {
  Status status =
      validate_bound(*this, kernels_, topology_, boundary_, metric_);
  const Int3 cells = kernels_->cells();
  if (!status ||
      !detail::valid_cell_view(transported, cells, 0U, 1U,
                               boundary_->maximum_halo_reach()) ||
      !detail::valid_cell_view(diffusivity, cells, 0U, 1U, 1U) ||
      !detail::valid_cell_view(rate, cells, 0U, 1U))
    return status ? Status{StatusCode::invalid_plan, kIbmEquationApply}
                  : status;
  const Span<const ImmersedLink> links = topology_->links();
  const Span<const BoundaryStencilLink> rows = boundary_->links();
  for (std::size_t index = 0U; index < rows.size; ++index) {
    const BoundaryStencilLink& row = rows.data[index];
    if (row.topology_link >= links.size)
      return {StatusCode::invalid_plan, kIbmEquationApply};
    const ImmersedLink& link = links.data[row.topology_link];
    const InterfaceFace face = interface_face(link);
    const double transmissibility = detail::positive_transmissibility(
        *kernels_, diffusivity, face.axis, face.index);
    const double volume = detail::cell_volume(*kernels_, link.fluid_local_index);
    double ghost = 0.0;
    status = evaluate_positive_bounded_quadratic_row(
        boundary_->reconstruction(), row.zero_normal_value_row, transported,
        0U, ghost);
    const double solid = transported.unchecked(link.solid_local_index, 0U);
    const double current = rate.unchecked(link.fluid_local_index, 0U);
    const double correction =
        transmissibility * (ghost - solid) / volume;
    const double value = current + correction;
    if (!status || !std::isfinite(transmissibility) ||
        transmissibility <= 0.0 || !std::isfinite(volume) || volume <= 0.0 ||
        !std::isfinite(ghost) || ghost <= 0.0 || !std::isfinite(solid) ||
        !std::isfinite(current) || !std::isfinite(correction) ||
        !std::isfinite(value))
      return status ? Status{StatusCode::numerical_failure,
                             kIbmEquationNumerical}
                    : status;
    rate.unchecked(link.fluid_local_index, 0U) += correction;
  }
  const Span<const std::uint8_t> region = topology_->region();
  std::size_t flat = 0U;
  for (std::int32_t z = 0; z < cells.z; ++z)
    for (std::int32_t y = 0; y < cells.y; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x, ++flat)
        if (region.data[flat] == static_cast<std::uint8_t>(RegionFlag::solid))
          rate.unchecked({x, y, z}, 0U) = 0.0;
  return {};
}

RevisionToken IbmEquationInterfacePlan::constrain_certificate(
    RevisionToken prior, RevisionToken field,
    RevisionToken coefficient) const noexcept {
  if (fingerprint_ == 0U || prior == 0U || field == 0U || coefficient == 0U)
    return 0U;
  std::uint64_t result = kFnvOffset;
  result = mix(result, prior);
  result = mix(result, field);
  result = mix(result, coefficient);
  result = mix(result, fingerprint_);
  return result == 0U ? 1U : result;
}

}  // namespace hundun::v04

// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_ibm.hpp"
#include "hundun/v04_flow.hpp"

#include "mesh_focus_detail.hpp"
#include "solver_equation_detail.hpp"
#include "solver_viscous_detail.hpp"
#include "solver_conservative_energy_detail.hpp"
#include "solver_ibm_scalar_transport_detail.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
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

std::uint64_t double_bits(double value) noexcept {
  std::uint64_t bits = 0U;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
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

double fluid_material_conductance(const CartesianKernelPlan& kernels,
                                 ConstFieldView diffusivity,
                                 const ImmersedLink& link,
                                 InterfaceFace face) noexcept {
  const std::int32_t normal = face.axis == CartesianAxis::x ? face.index.x
      : (face.axis == CartesianAxis::y ? face.index.y : face.index.z);
  return diffusivity.unchecked(link.fluid_local_index, 0U) *
         detail::face_area(kernels, face.axis, face.index) /
         (detail::centre_coordinate(kernels, face.axis, normal) -
          detail::centre_coordinate(kernels, face.axis, normal - 1));
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

// COAST pressure_delta_core mirrors the adjacent fluid difference across
// a cut face. Sum the two half differences before converting to Cartesian
// coordinates. Apply each cut axis once, including a one-cell fluid gap.
double coast_pressure_link_correction(const CartesianKernelPlan& kernels,
    const EBTopology& topology, const ImmersedLink& link,
    ConstFieldView pressure) noexcept {
  const auto axis = interface_face(link).axis;
  const unsigned a = static_cast<unsigned>(axis);
  const Int3 cell = link.fluid_local_index;
  const Int3 minus = offset(cell, axis, -1), plus = offset(cell, axis, 1);
  const auto extent = topology.global_cells();
  const auto fluid = [&](Int3 global) {
    // Physical ghost values remain owned by the ordinary boundary plan.
    if (global.x < 0 || global.y < 0 || global.z < 0 ||
        global.x >= extent.x || global.y >= extent.y || global.z >= extent.z)
      return true;
    return topology.is_fluid_global(global);
  };
  const bool fm = fluid(offset(link.fluid_global_index, axis, -1));
  const bool fp = fluid(offset(link.fluid_global_index, axis, 1));
  if (!fm && !fp && !positive_face(link.direction)) return 0.0;
  const int n = normal_index(cell, axis);
  const double pm = pressure.unchecked(minus, 0),
               pc = pressure.unchecked(cell, 0),
               pp = pressure.unchecked(plus, 0);
  const auto weights = kernels.geometry_kind() == GeometryKind::uniform
      ? detail::metric_derivative_weights<true>(kernels, a, n)
      : detail::metric_derivative_weights<false>(kernels, a, n);
  const double original = weights.minus * pm + weights.centre * pc + weights.plus * pp;
  const double dm = fm ? pc - pm : fp ? pp - pc : 0.0;
  const double dp = fp ? pp - pc : fm ? pc - pm : 0.0;
  const double desired = (dm + dp) /
      (detail::centre_coordinate(kernels, axis, n + 1) -
       detail::centre_coordinate(kernels, axis, n - 1));
  return desired - original;
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

FaceFieldView select(FrozenConvectionFaceOutput values,
                     CartesianAxis axis) noexcept {
  return axis == CartesianAxis::x
             ? values.x
             : (axis == CartesianAxis::y ? values.y : values.z);
}

bool valid_source_face_values(FrozenConvectionFaceOutput values,
                              Int3 cells) noexcept {
  const auto valid = [&](FaceFieldView view, CartesianAxis axis) noexcept {
    Int3 expected = cells;
    if (axis == CartesianAxis::x)
      ++expected.x;
    else if (axis == CartesianAxis::y)
      ++expected.y;
    else
      ++expected.z;
    detail::FieldStorageInterval interval;
    return detail::face_storage_interval(view, interval) &&
           view.axis == axis && view.extents.x == expected.x &&
           view.extents.y == expected.y && view.extents.z == expected.z &&
           view.storage_identity != 0U && view.revision_domain != 0U;
  };
  return valid(values.x, CartesianAxis::x) &&
         valid(values.y, CartesianAxis::y) &&
         valid(values.z, CartesianAxis::z) &&
         !detail::face_views_overlap(values.x, values.y) &&
         !detail::face_views_overlap(values.x, values.z) &&
         !detail::face_views_overlap(values.y, values.z);
}

bool same_index(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool inside(KernelBox box, Int3 cell) noexcept {
  return cell.x >= box.begin.x && cell.y >= box.begin.y &&
         cell.z >= box.begin.z && cell.x < box.begin.x + box.cells.x &&
         cell.y < box.begin.y + box.cells.y &&
         cell.z < box.begin.z + box.cells.z;
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

// Topology builders order links by fluid cell, then direction. Verify that
// invariant while visiting each boundary row once; no hot-path allocation.
template<class Evaluate>
Status replace_cut_scalar_rows(Span<const ImmersedLink> links, Int3 cells,
    KernelBox box, FieldView rate, Evaluate evaluate) noexcept {
  for(bool write : {false,true}) {
    std::size_t begin=0U,previous=0U; bool first=true;
    while(begin<links.size) {
      const Int3 c=links.data[begin].fluid_local_index;
      if(c.x<0 || c.x>=cells.x || c.y<0 || c.y>=cells.y || c.z<0 || c.z>=cells.z)
        return {StatusCode::invalid_plan,kIbmEquationApply};
      const std::size_t flat=std::size_t(c.x)+std::size_t(cells.x)*(std::size_t(c.y)+std::size_t(cells.y)*c.z);
      if(!first && flat<=previous) return {StatusCode::invalid_plan,kIbmEquationApply};
      first=false; previous=flat;
      std::array<std::size_t,6U> cut; cut.fill(links.size);
      std::size_t end=begin;
      while(end<links.size && same_index(links.data[end].fluid_local_index,c)) {
        const auto direction=static_cast<std::size_t>(links.data[end].direction);
        if(direction>=cut.size() || cut[direction]!=links.size)
          return {StatusCode::invalid_plan,kIbmEquationApply};
        cut[direction]=end++;
      }
      if(c.x>=box.begin.x && c.x<box.begin.x+box.cells.x &&
         c.y>=box.begin.y && c.y<box.begin.y+box.cells.y &&
         c.z>=box.begin.z && c.z<box.begin.z+box.cells.z) {
        double value=0.0; const Status status=evaluate(c,cut,value);
        if(!status) return status;
        if(write) rate.unchecked(c,0U)=value;
      }
      begin=end;
    }
  }
  return {};
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
    IbmEquationInterfacePlan& out,
    IbmPressureGradientKind pressure_gradient) noexcept {
  return compile_sources(
      kernels, topology, boundary, topology.interface_metric(), {}, {}, 0U,
      false, out, pressure_gradient);
}

Status IbmEquationInterfacePlan::compile(
    const CartesianKernelPlan& kernels, const EBTopology& topology,
    const BoundaryStencilPlan& boundary,
    Span<const IbmInterfaceMassFluxSource> mass_flux_sources,
    IbmEquationInterfacePlan& out,
    IbmPressureGradientKind pressure_gradient) noexcept {
  return compile_sources(kernels, topology, boundary,
                         topology.interface_metric(), mass_flux_sources, {},
                         0U, false, out, pressure_gradient);
}

Status IbmEquationInterfacePlan::compile(
    const CartesianKernelPlan& kernels, const EBTopology& topology,
    const BoundaryStencilPlan& boundary,
    Span<const IbmInterfaceInletState> inlet_states,
    std::size_t independent_species_count,
    IbmEquationInterfacePlan& out,
    IbmPressureGradientKind pressure_gradient) noexcept {
  return compile_sources(kernels, topology, boundary,
                         topology.interface_metric(), {}, inlet_states,
                         independent_species_count, true, out, pressure_gradient);
}

Status IbmEquationInterfacePlan::compile(
    const CartesianKernelPlan& kernels, const EBTopology& topology,
    const BoundaryStencilPlan& boundary,
    const IbmInterfaceMetricPlan& metric,
    IbmEquationInterfacePlan& out,
    IbmPressureGradientKind pressure_gradient) noexcept {
  return compile_sources(kernels, topology, boundary, metric, {}, {}, 0U,
                         false, out, pressure_gradient);
}

Status IbmEquationInterfacePlan::compile(
    const CartesianKernelPlan& kernels, const EBTopology& topology,
    const BoundaryStencilPlan& boundary,
    const IbmInterfaceMetricPlan& metric,
    Span<const IbmInterfaceMassFluxSource> mass_flux_sources,
    IbmEquationInterfacePlan& out,
    IbmPressureGradientKind pressure_gradient) noexcept {
  return compile_sources(kernels, topology, boundary, metric,
                         mass_flux_sources, {}, 0U, false, out, pressure_gradient);
}

Status IbmEquationInterfacePlan::compile(
    const CartesianKernelPlan& kernels, const EBTopology& topology,
    const BoundaryStencilPlan& boundary,
    const IbmInterfaceMetricPlan& metric,
    Span<const IbmInterfaceInletState> inlet_states,
    std::size_t independent_species_count,
    IbmEquationInterfacePlan& out,
    IbmPressureGradientKind pressure_gradient) noexcept {
  return compile_sources(kernels, topology, boundary, metric, {},
                         inlet_states, independent_species_count, true, out, pressure_gradient);
}

Status IbmEquationInterfacePlan::compile_sources(
    const CartesianKernelPlan& kernels, const EBTopology& topology,
    const BoundaryStencilPlan& boundary,
    const IbmInterfaceMetricPlan& metric,
    Span<const IbmInterfaceMassFluxSource> mass_flux_sources,
    Span<const IbmInterfaceInletState> inlet_states,
    std::size_t independent_species_count, bool inlet_state_bound,
    IbmEquationInterfacePlan& out,
    IbmPressureGradientKind pressure_gradient) noexcept {
  if (pressure_gradient != IbmPressureGradientKind::quadratic_neumann &&
      pressure_gradient != IbmPressureGradientKind::coast_fluid_delta)
    return {StatusCode::invalid_plan, kIbmEquationPlan};
  if (!valid_plan_inputs(kernels, topology, boundary, metric) ||
      (mass_flux_sources.size != 0U && mass_flux_sources.data == nullptr) ||
      (inlet_states.size != 0U && inlet_states.data == nullptr) ||
      (inlet_state_bound
           ? mass_flux_sources.size != 0U
           : inlet_states.size != 0U || independent_species_count != 0U))
    return {StatusCode::invalid_plan, kIbmEquationPlan};
  std::uint64_t fingerprint = kFnvOffset;
  fingerprint = mix(fingerprint, kernels.fingerprint());
  fingerprint = mix(fingerprint, topology.fingerprint());
  fingerprint = mix(fingerprint, boundary.fingerprint());
  fingerprint = mix(fingerprint, metric.fingerprint());
  fingerprint = mix(fingerprint, topology.geometry_revision());
  if (pressure_gradient == IbmPressureGradientKind::coast_fluid_delta)
    fingerprint = mix(fingerprint, UINT64_C(0x434f415354504431));
  if (fingerprint == 0U) fingerprint = 1U;
  IbmEquationInterfacePlan candidate;
  candidate.kernels_ = &kernels;
  candidate.topology_ = &topology;
  candidate.boundary_ = &boundary;
  candidate.metric_ = &metric;
  candidate.independent_species_count_ = independent_species_count;
  candidate.inlet_state_bound_ = inlet_state_bound;
  candidate.pressure_gradient_ = pressure_gradient;
  try {
    const Int3 cells=kernels.cells();
    const auto region=topology.region();
    candidate.scalar_blocked_faces_.resize(region.size);
    for(std::size_t i=0U;i<region.size;++i)
      candidate.scalar_blocked_faces_[i]=region.data[i]==0U ? 63U : 0U;
    const auto scalar_links=topology.links();
    for(std::size_t l=0U;l<scalar_links.size;++l) {
      const auto& link=scalar_links.data[l];
      const Int3 c=link.fluid_local_index;
      const auto i=std::size_t(c.x)+std::size_t(cells.x)*(std::size_t(c.y)+std::size_t(cells.y)*c.z);
      candidate.scalar_blocked_faces_[i] |= std::uint8_t(1U << unsigned(link.direction));
    }
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
    const std::size_t source_count = inlet_state_bound
                                         ? inlet_states.size
                                         : mass_flux_sources.size;
    if (source_count != 0U) {
      for (std::size_t index = 1U; index < topology_links.size; ++index)
        if (topology_links.data[index - 1U].global_link >=
            topology_links.data[index].global_link)
          return {StatusCode::invalid_plan, kIbmEquationPlan};
    }
    const auto resolve_topology_link =
        [&](std::uint64_t global_link, std::size_t& resolved) noexcept {
          std::size_t begin = 0U;
          std::size_t end = topology_links.size;
          while (begin < end) {
            const std::size_t middle = begin + (end - begin) / 2U;
            if (topology_links.data[middle].global_link < global_link)
              begin = middle + 1U;
            else
              end = middle;
          }
          resolved = begin;
          return begin < topology_links.size && begin <= UINT32_MAX &&
                 topology_links.data[begin].global_link == global_link;
        };
    if (!inlet_state_bound && mass_flux_sources.size != 0U) {
      std::vector<IbmInterfaceMassFluxSource> ordered_sources(
          mass_flux_sources.data,
          mass_flux_sources.data + mass_flux_sources.size);
      std::sort(ordered_sources.begin(), ordered_sources.end(),
                [](const IbmInterfaceMassFluxSource& left,
                   const IbmInterfaceMassFluxSource& right) noexcept {
                  return left.global_link < right.global_link;
                });
      candidate.prescribed_interface_fluxes_.reserve(ordered_sources.size());
      fingerprint = mix(fingerprint, UINT64_C(0x69626d736f757263));
      fingerprint = mix(fingerprint, ordered_sources.size());
      for (std::size_t source_index = 0U;
           source_index < ordered_sources.size(); ++source_index) {
        IbmInterfaceMassFluxSource source = ordered_sources[source_index];
        if (!std::isfinite(source.face_mass_flux) ||
            (source_index != 0U &&
             ordered_sources[source_index - 1U].global_link ==
                 source.global_link))
          return {StatusCode::invalid_plan, kIbmEquationPlan};
        if (source.face_mass_flux == 0.0) source.face_mass_flux = 0.0;
        std::size_t topology_link = 0U;
        if (!resolve_topology_link(source.global_link, topology_link))
          return {StatusCode::invalid_plan, kIbmEquationPlan};
        candidate.prescribed_interface_fluxes_.push_back(
            {static_cast<std::uint32_t>(topology_link),
             source.face_mass_flux});
        fingerprint = mix(fingerprint, source.global_link);
        fingerprint = mix(fingerprint, double_bits(source.face_mass_flux));
      }
    } else if (inlet_state_bound) {
      std::vector<IbmInterfaceInletState> ordered_states;
      if (inlet_states.size != 0U)
        ordered_states.assign(inlet_states.data,
                              inlet_states.data + inlet_states.size);
      std::sort(ordered_states.begin(), ordered_states.end(),
                [](const IbmInterfaceInletState& left,
                   const IbmInterfaceInletState& right) noexcept {
                  return left.global_link < right.global_link;
                });
      candidate.prescribed_interface_fluxes_.reserve(ordered_states.size());
      if (independent_species_count != 0U &&
          ordered_states.size() >
              std::numeric_limits<std::size_t>::max() /
                  independent_species_count)
        return {StatusCode::invalid_plan, kIbmEquationPlan};
      candidate.prescribed_independent_species_.reserve(
          ordered_states.size() * independent_species_count);
      fingerprint = mix(fingerprint, UINT64_C(0x69626d696e6c6574));
      fingerprint = mix(fingerprint, ordered_states.size());
      fingerprint = mix(fingerprint, independent_species_count);
      for (std::size_t source_index = 0U;
           source_index < ordered_states.size(); ++source_index) {
        const IbmInterfaceInletState& source = ordered_states[source_index];
        if (!std::isfinite(source.face_mass_flux) ||
            !std::isfinite(source.velocity.x) ||
            !std::isfinite(source.velocity.y) ||
            !std::isfinite(source.velocity.z) ||
            !std::isfinite(source.enthalpy) ||
            source.independent_species.size != independent_species_count ||
            (independent_species_count != 0U &&
             source.independent_species.data == nullptr) ||
            (source_index != 0U &&
             ordered_states[source_index - 1U].global_link ==
                 source.global_link))
          return {StatusCode::invalid_plan, kIbmEquationPlan};
        std::size_t topology_link = 0U;
        if (!resolve_topology_link(source.global_link, topology_link))
          return {StatusCode::invalid_plan, kIbmEquationPlan};
        const ImmersedFaceDirection direction =
            topology_links.data[topology_link].direction;
        const bool fluid_on_positive_side =
            direction == ImmersedFaceDirection::x_negative ||
            direction == ImmersedFaceDirection::y_negative ||
            direction == ImmersedFaceDirection::z_negative;
        if (fluid_on_positive_side ? !(source.face_mass_flux > 0.0)
                                   : !(source.face_mass_flux < 0.0))
          return {StatusCode::invalid_plan, kIbmEquationPlan};
        const double normal_velocity =
            direction == ImmersedFaceDirection::x_negative ||
                    direction == ImmersedFaceDirection::x_positive
                ? source.velocity.x
                : (direction == ImmersedFaceDirection::y_negative ||
                           direction == ImmersedFaceDirection::y_positive
                       ? source.velocity.y
                       : source.velocity.z);
        // The Cartesian mass flux and velocity component use the same
        // negative-to-positive axis orientation.  A prescribed inlet must
        // therefore drive both quantities from the solid side into the fluid
        // side; accepting opposite signs would bind mutually inconsistent
        // convective and momentum states to one source face.
        if (source.face_mass_flux > 0.0 ? !(normal_velocity > 0.0)
                                        : !(normal_velocity < 0.0))
          return {StatusCode::invalid_plan, kIbmEquationPlan};
        const std::size_t species_begin =
            candidate.prescribed_independent_species_.size();
        for (std::size_t species = 0U;
             species < independent_species_count; ++species) {
          const double value = source.independent_species.data[species];
          if (!std::isfinite(value))
            return {StatusCode::invalid_plan, kIbmEquationPlan};
          candidate.prescribed_independent_species_.push_back(value);
        }
        candidate.prescribed_interface_fluxes_.push_back(
            {static_cast<std::uint32_t>(topology_link), source.face_mass_flux,
             source.velocity, source.enthalpy, species_begin, true});
        fingerprint = mix(fingerprint, source.global_link);
        fingerprint = mix(fingerprint, double_bits(source.face_mass_flux));
        fingerprint = mix(fingerprint, double_bits(source.velocity.x));
        fingerprint = mix(fingerprint, double_bits(source.velocity.y));
        fingerprint = mix(fingerprint, double_bits(source.velocity.z));
        fingerprint = mix(fingerprint, double_bits(source.enthalpy));
        for (std::size_t species = 0U;
             species < independent_species_count; ++species)
          fingerprint = mix(
              fingerprint,
              double_bits(source.independent_species.data[species]));
      }
    }
    if (inlet_state_bound && !candidate.prescribed_interface_fluxes_.empty()) {
      candidate.prescribed_source_faces_.reserve(
          candidate.prescribed_interface_fluxes_.size());
      for (const PrescribedInterfaceFlux& source :
           candidate.prescribed_interface_fluxes_) {
        if (source.topology_link >= topology_links.size)
          return {StatusCode::invalid_plan, kIbmEquationPlan};
        const InterfaceFace face =
            interface_face(topology_links.data[source.topology_link]);
        candidate.prescribed_source_faces_.push_back({face.axis, face.index});
      }
      std::sort(candidate.prescribed_source_faces_.begin(),
                candidate.prescribed_source_faces_.end(),
                [](const FrozenConvectionFixedFace& left,
                   const FrozenConvectionFixedFace& right) noexcept {
                  const auto left_axis =
                      static_cast<std::uint8_t>(left.axis);
                  const auto right_axis =
                      static_cast<std::uint8_t>(right.axis);
                  if (left_axis != right_axis) return left_axis < right_axis;
                  if (left.index.z != right.index.z)
                    return left.index.z < right.index.z;
                  if (left.index.y != right.index.y)
                    return left.index.y < right.index.y;
                  return left.index.x < right.index.x;
                });
      for (std::size_t index = 1U;
           index < candidate.prescribed_source_faces_.size(); ++index) {
        const FrozenConvectionFixedFace& prior =
            candidate.prescribed_source_faces_[index - 1U];
        const FrozenConvectionFixedFace& current =
            candidate.prescribed_source_faces_[index];
        if (prior.axis == current.axis && prior.index.x == current.index.x &&
            prior.index.y == current.index.y &&
            prior.index.z == current.index.z) {
          return {StatusCode::invalid_plan, kIbmEquationPlan};
        }
      }
    }
    if (fingerprint == 0U) fingerprint = 1U;
    candidate.fingerprint_ = fingerprint;
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

Status IbmEquationInterfacePlan::constrain_interface_flux(
    FaceFluxView flux) const noexcept {
  Status status =
      validate_bound(*this, kernels_, topology_, boundary_, metric_);
  if (!status || !detail::valid_flux_view(flux, kernels_->cells()))
    return status ? Status{StatusCode::invalid_plan, kIbmEquationApply}
                  : status;
  const Span<const ImmersedLink> links = topology_->links();
  for (std::size_t index = 0U;
       index < prescribed_interface_fluxes_.size(); ++index) {
    const PrescribedInterfaceFlux& source =
        prescribed_interface_fluxes_[index];
    if (source.topology_link >= links.size ||
        !std::isfinite(source.face_mass_flux) ||
        (index != 0U &&
         prescribed_interface_fluxes_[index - 1U].topology_link >=
             source.topology_link))
      return {StatusCode::invalid_plan, kIbmEquationApply};
  }
  for (std::size_t index = 0U; index < links.size; ++index) {
    const InterfaceFace face = interface_face(links.data[index]);
    select(flux, face.axis).unchecked(face.index) = 0.0;
  }
  // Inactive cells are storage placeholders, not fluid control volumes.
  // Retiring only cut faces leaves old solid-solid fluxes in BDF history.
  const Int3 cells = kernels_->cells();
  const auto region = topology_->region();
  std::size_t flat = 0U;
  for (int z = 0; z < cells.z; ++z)
    for (int y = 0; y < cells.y; ++y)
      for (int x = 0; x < cells.x; ++x, ++flat)
        if (region.data[flat] == 0U) {
          flux.x.unchecked({x, y, z}) = flux.x.unchecked({x + 1, y, z}) = 0.0;
          flux.y.unchecked({x, y, z}) = flux.y.unchecked({x, y + 1, z}) = 0.0;
          flux.z.unchecked({x, y, z}) = flux.z.unchecked({x, y, z + 1}) = 0.0;
        }
  // Solid-cell retirement above also touches every immersed source face.
  // Reapply prescribed Cartesian fluxes last so these values are the sole
  // continuity authority on their links.
  for (const PrescribedInterfaceFlux& source : prescribed_interface_fluxes_) {
    const InterfaceFace face = interface_face(links.data[source.topology_link]);
    select(flux, face.axis).unchecked(face.index) = source.face_mass_flux;
  }
  return {};
}

Status IbmEquationInterfacePlan::zero_interface_flux(
    FaceFluxView flux) const noexcept {
  return constrain_interface_flux(flux);
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
  for (std::size_t index = 0U;
       index < prescribed_interface_fluxes_.size(); ++index) {
    const PrescribedInterfaceFlux& source =
        prescribed_interface_fluxes_[index];
    if (source.topology_link >= links.size ||
        !std::isfinite(source.face_mass_flux) ||
        (index != 0U &&
         prescribed_interface_fluxes_[index - 1U].topology_link >=
             source.topology_link))
      return {StatusCode::invalid_plan, kIbmEquationApply};
  }
  std::size_t source_index = 0U;
  for (std::size_t index = 0U; index < links.size; ++index) {
    const InterfaceFace face = interface_face(links.data[index]);
    const double value = select(flux, face.axis).unchecked(face.index);
    double expected = 0.0;
    if (source_index < prescribed_interface_fluxes_.size() &&
        prescribed_interface_fluxes_[source_index].topology_link == index) {
      expected = prescribed_interface_fluxes_[source_index].face_mass_flux;
      ++source_index;
    }
    if (!std::isfinite(value) ||
        std::abs(value - expected) > absolute_tolerance)
      return {StatusCode::numerical_failure, kIbmEquationNumerical};
  }
  return {};
}

bool IbmEquationInterfacePlan::prescribed_face_flux(
    CartesianAxis axis, Int3 face, double& phi) const noexcept {
  if (!validate_bound(*this, kernels_, topology_, boundary_, metric_))
    return false;
  const Span<const ImmersedLink> links = topology_->links();
  for (const PrescribedInterfaceFlux& source :
       prescribed_interface_fluxes_) {
    if (source.topology_link >= links.size ||
        !std::isfinite(source.face_mass_flux))
      return false;
    const InterfaceFace candidate =
        interface_face(links.data[source.topology_link]);
    if (candidate.axis == axis && same_index(candidate.index, face)) {
      phi = source.face_mass_flux;
      return true;
    }
  }
  return false;
}

Status IbmEquationInterfacePlan::override_source_face_values(
    IbmInterfaceInletField field,
    IbmInterfaceInletEvaluation evaluation,
    FrozenConvectionFaceOutput values) const noexcept {
  const Status bound =
      validate_bound(*this, kernels_, topology_, boundary_, metric_);
  if (!bound) return bound;
  if (!valid_source_face_values(values, kernels_->cells()) ||
      (evaluation != IbmInterfaceInletEvaluation::value &&
       evaluation != IbmInterfaceInletEvaluation::fixed_state_variation))
    return {StatusCode::invalid_plan, kIbmEquationApply};
  if (prescribed_interface_fluxes_.empty()) return {};
  const bool valid_field = [&]() noexcept {
    switch (field.kind) {
      case IbmInterfaceInletFieldKind::velocity:
        return field.component < 3U;
      case IbmInterfaceInletFieldKind::enthalpy:
        return field.component == 0U;
      case IbmInterfaceInletFieldKind::independent_species:
        return field.component < independent_species_count_;
      case IbmInterfaceInletFieldKind::kinetic_energy:
        return field.component == 0U;
    }
    return false;
  }();
  if (!inlet_state_bound_ || !valid_field)
    return {StatusCode::invalid_plan, kIbmEquationApply};

  const Span<const ImmersedLink> links = topology_->links();
  for (std::size_t index = 0U;
       index < prescribed_interface_fluxes_.size(); ++index) {
    const PrescribedInterfaceFlux& source =
        prescribed_interface_fluxes_[index];
    if (source.topology_link >= links.size || !source.has_inlet_state ||
        !std::isfinite(source.velocity.x) ||
        !std::isfinite(source.velocity.y) ||
        !std::isfinite(source.velocity.z) ||
        !std::isfinite(source.enthalpy) ||
        source.independent_species_begin >
            prescribed_independent_species_.size() ||
        independent_species_count_ >
            prescribed_independent_species_.size() -
                source.independent_species_begin ||
        (index != 0U &&
         prescribed_interface_fluxes_[index - 1U].topology_link >=
             source.topology_link))
      return {StatusCode::invalid_plan, kIbmEquationApply};
  }

  for (const PrescribedInterfaceFlux& source :
       prescribed_interface_fluxes_) {
    double prescribed = 0.0;
    if (evaluation == IbmInterfaceInletEvaluation::value) {
      switch (field.kind) {
        case IbmInterfaceInletFieldKind::velocity:
          prescribed = field.component == 0U
                           ? source.velocity.x
                           : (field.component == 1U ? source.velocity.y
                                                   : source.velocity.z);
          break;
        case IbmInterfaceInletFieldKind::enthalpy:
          prescribed = source.enthalpy;
          break;
        case IbmInterfaceInletFieldKind::independent_species:
          prescribed = prescribed_independent_species_[
              source.independent_species_begin + field.component];
          break;
        case IbmInterfaceInletFieldKind::kinetic_energy:
          prescribed =
              0.5 * (source.velocity.x * source.velocity.x +
                     source.velocity.y * source.velocity.y +
                     source.velocity.z * source.velocity.z);
          break;
      }
    }
    const InterfaceFace face =
        interface_face(links.data[source.topology_link]);
    select(values, face.axis).unchecked(face.index) = prescribed;
  }
  return {};
}

Status IbmEquationInterfacePlan::freeze_source_convection_faces(
    IbmInterfaceInletField field, ConvectionScheme scheme,
    ConstFaceFluxView target_flux, ConstFieldView transported,
    std::uint8_t component, FrozenConvectionContext context,
    FrozenConvectionFaceOutput output,
    FrozenConvectionFaceField& frozen) const noexcept {
  frozen = {};
  const Status bound =
      validate_bound(*this, kernels_, topology_, boundary_, metric_);
  if (!bound) return bound;
  if (!has_inlet_sources() || prescribed_source_faces_.size() !=
                                  prescribed_interface_fluxes_.size()) {
    return {StatusCode::invalid_plan, kIbmEquationApply};
  }
  const std::uint8_t expected_component =
      field.kind == IbmInterfaceInletFieldKind::velocity
          ? static_cast<std::uint8_t>(field.component)
          : 0U;
  if (field.component > std::numeric_limits<std::uint8_t>::max() ||
      component != expected_component) {
    return {StatusCode::invalid_plan, kIbmEquationApply};
  }
  Status status = validate_interface_flux(target_flux);
  if (status) {
    status = freeze_cartesian_target_convection_faces(
        *kernels_, scheme, target_flux, transported, component, context,
        output, frozen);
  }
  if (status) {
    status = override_source_face_values(
        field, IbmInterfaceInletEvaluation::value, output);
  }
  if (status) {
    status = seal_fixed_cartesian_target_convection_faces(
        *kernels_, scheme, target_flux, transported, component, context,
        {prescribed_source_faces_.data(), prescribed_source_faces_.size()},
        fingerprint_, frozen);
  }
  if (status) status = validate_frozen_source_face_values(field, frozen);
  if (!status) frozen = {};
  return status;
}

Status IbmEquationInterfacePlan::validate_frozen_source_face_values(
    IbmInterfaceInletField field,
    const FrozenConvectionFaceField& frozen) const noexcept {
  const Status bound =
      validate_bound(*this, kernels_, topology_, boundary_, metric_);
  if (!bound) return bound;
  const bool valid_field = [&]() noexcept {
    switch (field.kind) {
      case IbmInterfaceInletFieldKind::velocity:
        return field.component < 3U;
      case IbmInterfaceInletFieldKind::enthalpy:
      case IbmInterfaceInletFieldKind::kinetic_energy:
        return field.component == 0U;
      case IbmInterfaceInletFieldKind::independent_species:
        return field.component < independent_species_count_;
    }
    return false;
  }();
  const auto valid_face = [&](ConstFaceFieldView view,
                              CartesianAxis axis) noexcept {
    Int3 expected = kernels_->cells();
    if (axis == CartesianAxis::x)
      ++expected.x;
    else if (axis == CartesianAxis::y)
      ++expected.y;
    else
      ++expected.z;
    detail::FieldStorageInterval interval;
    return detail::face_storage_interval(view, interval) &&
           view.axis == axis && same_index(view.extents, expected) &&
           view.storage_identity != 0U && view.revision_domain != 0U;
  };
  if (!has_inlet_sources() || !valid_field || !frozen.valid() ||
      !valid_face(frozen.x, CartesianAxis::x) ||
      !valid_face(frozen.y, CartesianAxis::y) ||
      !valid_face(frozen.z, CartesianAxis::z) ||
      frozen.fixed_face_authority != fingerprint_ ||
      frozen.fixed_faces.data != prescribed_source_faces_.data() ||
      frozen.fixed_faces.size != prescribed_source_faces_.size() ||
      prescribed_source_faces_.size() !=
          prescribed_interface_fluxes_.size()) {
    return {StatusCode::invalid_plan, kIbmEquationApply};
  }

  const Span<const ImmersedLink> links = topology_->links();
  for (const PrescribedInterfaceFlux& source :
       prescribed_interface_fluxes_) {
    if (source.topology_link >= links.size || !source.has_inlet_state ||
        source.independent_species_begin >
            prescribed_independent_species_.size() ||
        independent_species_count_ >
            prescribed_independent_species_.size() -
                source.independent_species_begin) {
      return {StatusCode::invalid_plan, kIbmEquationApply};
    }
    double expected = 0.0;
    switch (field.kind) {
      case IbmInterfaceInletFieldKind::velocity:
        expected = field.component == 0U
                       ? source.velocity.x
                       : (field.component == 1U ? source.velocity.y
                                                : source.velocity.z);
        break;
      case IbmInterfaceInletFieldKind::enthalpy:
        expected = source.enthalpy;
        break;
      case IbmInterfaceInletFieldKind::independent_species:
        expected = prescribed_independent_species_[
            source.independent_species_begin + field.component];
        break;
      case IbmInterfaceInletFieldKind::kinetic_energy:
        expected = 0.5 * (source.velocity.x * source.velocity.x +
                          source.velocity.y * source.velocity.y +
                          source.velocity.z * source.velocity.z);
        break;
    }
    const InterfaceFace face =
        interface_face(links.data[source.topology_link]);
    const ConstFaceFieldView values =
        face.axis == CartesianAxis::x
            ? frozen.x
            : (face.axis == CartesianAxis::y ? frozen.y : frozen.z);
    const double actual = values.unchecked(face.index);
    if (!std::isfinite(expected) || !std::isfinite(actual) ||
        double_bits(actual) != double_bits(expected)) {
      return {StatusCode::invalid_plan, kIbmEquationApply};
    }
  }
  return {};
}

Status IbmEquationInterfacePlan::add_source_convection_correction(
    IbmInterfaceInletField field, ConvectionScheme scheme,
    ConstFieldView transported, double scale, FieldView output,
    KernelBox box) const noexcept {
  return add_source_convection_correction_impl(
      field, &scheme, transported, scale, output, box);
}

Status IbmEquationInterfacePlan::add_source_first_order_upwind_correction(
    IbmInterfaceInletField field, ConstFieldView transported, double scale,
    FieldView output, KernelBox box) const noexcept {
  return add_source_convection_correction_impl(
      field, nullptr, transported, scale, output, box);
}

Status IbmEquationInterfacePlan::add_source_convection_correction_impl(
    IbmInterfaceInletField field, const ConvectionScheme* scheme,
    ConstFieldView transported, double scale, FieldView output,
    KernelBox box, bool kinetic_from_velocity) const noexcept {
  const Status bound =
      validate_bound(*this, kernels_, topology_, boundary_, metric_);
  if (!bound) return bound;
  const Int3 cells = kernels_->cells();
  if (box.begin.x == 0 && box.begin.y == 0 && box.begin.z == 0 &&
      box.cells.x == 0 && box.cells.y == 0 && box.cells.z == 0)
    box = {{0, 0, 0}, cells};
  if (!detail::valid_kernel_box(box, cells) || !std::isfinite(scale) ||
      (scheme != nullptr &&
       static_cast<std::uint8_t>(*scheme) >
           static_cast<std::uint8_t>(ConvectionScheme::tvd2)))
    return {StatusCode::invalid_plan, kIbmEquationApply};

  const std::uint8_t transported_component =
      field.kind == IbmInterfaceInletFieldKind::velocity
          ? static_cast<std::uint8_t>(field.component)
          : 0U;
  const std::uint8_t output_component = transported_component;
  const std::uint8_t required_ghost_width =
      scheme == nullptr || *scheme == ConvectionScheme::central2 ? 1U : 2U;
  if (!detail::valid_cell_view(transported, cells, transported_component,
                               kinetic_from_velocity ? 3U : 1U,
                               required_ghost_width) ||
      !detail::valid_cell_view(output, cells, output_component, 1U) ||
      detail::field_views_overlap(transported, output))
    return {StatusCode::invalid_plan, kIbmEquationApply};
  if (prescribed_interface_fluxes_.empty()) return {};

  const bool valid_field = [&]() noexcept {
    switch (field.kind) {
      case IbmInterfaceInletFieldKind::velocity:
        return field.component < 3U;
      case IbmInterfaceInletFieldKind::enthalpy:
      case IbmInterfaceInletFieldKind::kinetic_energy:
        return field.component == 0U;
      case IbmInterfaceInletFieldKind::independent_species:
        return field.component < independent_species_count_;
    }
    return false;
  }();
  if (!inlet_state_bound_ || !valid_field)
    return {StatusCode::invalid_plan, kIbmEquationApply};

  const Span<const ImmersedLink> links = topology_->links();
  const auto prescribed_value = [&](const PrescribedInterfaceFlux& source,
                                    double& value) noexcept {
    if (source.topology_link >= links.size || !source.has_inlet_state ||
        !std::isfinite(source.face_mass_flux) ||
        source.independent_species_begin >
            prescribed_independent_species_.size() ||
        independent_species_count_ >
            prescribed_independent_species_.size() -
                source.independent_species_begin)
      return false;
    switch (field.kind) {
      case IbmInterfaceInletFieldKind::velocity:
        value = field.component == 0U
                    ? source.velocity.x
                    : (field.component == 1U ? source.velocity.y
                                             : source.velocity.z);
        break;
      case IbmInterfaceInletFieldKind::enthalpy:
        value = source.enthalpy;
        break;
      case IbmInterfaceInletFieldKind::independent_species:
        value = prescribed_independent_species_[
            source.independent_species_begin + field.component];
        break;
      case IbmInterfaceInletFieldKind::kinetic_energy:
        value = 0.5 * (source.velocity.x * source.velocity.x +
                       source.velocity.y * source.velocity.y +
                       source.velocity.z * source.velocity.z);
        break;
    }
    return std::isfinite(value);
  };
  const auto correction = [&](const PrescribedInterfaceFlux& source,
                              double& value) noexcept -> Status {
    if (source.topology_link >= links.size)
      return {StatusCode::invalid_plan, kIbmEquationApply};
    const ImmersedLink& link = links.data[source.topology_link];
    double prescribed = 0.0;
    if (!prescribed_value(source, prescribed))
      return {StatusCode::invalid_plan, kIbmEquationApply};
    double ordinary = 0.0;
    if (kinetic_from_velocity) {
      const InterfaceFace face = interface_face(link);
      ordinary = detail::kinetic_convection_face(*kernels_, *scheme,
          transported, face.axis, face.index, source.face_mass_flux);
    } else if (scheme == nullptr) {
      ordinary = transported.unchecked(link.solid_local_index,
                                       transported_component);
    } else {
      const InterfaceFace face = interface_face(link);
      const Status reconstructed = reconstruct_cartesian_convection_face(
          *kernels_, *scheme, transported, transported_component, face.axis,
          face.index, source.face_mass_flux, ordinary);
      if (!reconstructed) return reconstructed;
    }
    const double volume = detail::cell_volume(*kernels_, link.fluid_local_index);
    const double sign = positive_face(link.direction) ? 1.0 : -1.0;
    value = scale * sign * source.face_mass_flux *
            (prescribed - ordinary) / volume;
    if (!std::isfinite(ordinary) || !std::isfinite(volume) || volume <= 0.0 ||
        !std::isfinite(value))
      return {StatusCode::numerical_failure, kIbmEquationNumerical};
    return {};
  };
  const auto first_for_cell = [&](std::size_t index) noexcept {
    const Int3 cell = links.data[prescribed_interface_fluxes_[index]
                                     .topology_link]
                          .fluid_local_index;
    for (std::size_t prior = 0U; prior < index; ++prior) {
      const PrescribedInterfaceFlux& candidate =
          prescribed_interface_fluxes_[prior];
      if (candidate.topology_link < links.size &&
          same_index(links.data[candidate.topology_link].fluid_local_index,
                     cell))
        return false;
    }
    return true;
  };
  const auto cell_correction = [&](std::size_t index,
                                   double& value) noexcept -> Status {
    const Int3 cell = links.data[prescribed_interface_fluxes_[index]
                                     .topology_link]
                          .fluid_local_index;
    long double accumulated = 0.0L;
    for (const PrescribedInterfaceFlux& source :
         prescribed_interface_fluxes_) {
      if (source.topology_link >= links.size)
        return {StatusCode::invalid_plan, kIbmEquationApply};
      if (!same_index(links.data[source.topology_link].fluid_local_index,
                      cell))
        continue;
      double term = 0.0;
      const Status evaluated = correction(source, term);
      if (!evaluated) return evaluated;
      accumulated += static_cast<long double>(term);
    }
    value = static_cast<double>(accumulated);
    return std::isfinite(value)
               ? Status{}
               : Status{StatusCode::numerical_failure,
                        kIbmEquationNumerical};
  };

  // Preflight complete cell aggregates so a multi-face inlet cannot partially
  // modify the caller's output on arithmetic failure.
  for (std::size_t index = 0U;
       index < prescribed_interface_fluxes_.size(); ++index) {
    const PrescribedInterfaceFlux& source =
        prescribed_interface_fluxes_[index];
    if (source.topology_link >= links.size)
      return {StatusCode::invalid_plan, kIbmEquationApply};
    const Int3 cell = links.data[source.topology_link].fluid_local_index;
    if (!inside(box, cell) || !first_for_cell(index)) continue;
    double change = 0.0;
    const Status evaluated = cell_correction(index, change);
    const double candidate = output.unchecked(cell, output_component) + change;
    if (!evaluated) return evaluated;
    if (!std::isfinite(output.unchecked(cell, output_component)) ||
        !std::isfinite(candidate))
      return {StatusCode::numerical_failure, kIbmEquationNumerical};
  }
  for (std::size_t index = 0U;
       index < prescribed_interface_fluxes_.size(); ++index) {
    const PrescribedInterfaceFlux& source =
        prescribed_interface_fluxes_[index];
    const Int3 cell = links.data[source.topology_link].fluid_local_index;
    if (!inside(box, cell) || !first_for_cell(index)) continue;
    double change = 0.0;
    const Status evaluated = cell_correction(index, change);
    if (!evaluated) return evaluated;
    output.unchecked(cell, output_component) += change;
  }
  return {};
}

Status IbmEquationInterfacePlan::constrain_pressure_predictor(
    FieldView h_by_a, FaceFluxView phi_h_by_a) const noexcept {
  Status status =
      validate_bound(*this, kernels_, topology_, boundary_, metric_);
  if (!status) return status;
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
  return constrain_interface_flux(phi_h_by_a);
}

Status IbmEquationInterfacePlan::constrain_corrected_state(
    FieldView velocity, FaceFluxView flux) const noexcept {
  return constrain_pressure_predictor(velocity, flux);
}

const IbmEquationInterfacePlan::PrescribedInterfaceFlux*
IbmEquationInterfacePlan::inlet_for_link(
    std::uint32_t topology_link) const noexcept {
  const auto found = std::lower_bound(
      prescribed_interface_fluxes_.begin(), prescribed_interface_fluxes_.end(),
      topology_link, [](const PrescribedInterfaceFlux& source,
                        std::uint32_t index) {
        return source.topology_link < index;
      });
  return found != prescribed_interface_fluxes_.end() &&
                 found->topology_link == topology_link && found->has_inlet_state
             ? &*found : nullptr;
}

// Shared resolved/wall-law traction for momentum, boundary work and heating.
Status IbmEquationInterfacePlan::viscous_boundary_traction(
    std::size_t index, ConstFieldView velocity, ConstFieldView density,
    ConstFieldView molecular_viscosity, ConstFieldView effective_viscosity,
    const TurbulencePlan* wall_treatment,
    ViscousBoundaryTraction& out) const noexcept {
  const auto rows = boundary_->links();
  const auto links = topology_->links();
  const auto physical_links = metric_->links();
  if (index >= rows.size || index >= wall_linearization_.size() ||
      rows.data[index].topology_link >= links.size ||
      physical_links.size != links.size)
    return {StatusCode::invalid_plan, kIbmEquationApply};
  const auto& row = rows.data[index];
  const auto& link = links.data[row.topology_link];
  const auto& physical = physical_links.data[row.topology_link];
  const auto cell = link.fluid_local_index;
  const auto linearization = wall_linearization_[index];
  const double area = physical.physical_quadrature_area;
  double wall_viscosity = 0.0;
  Status status = evaluate_positive_bounded_quadratic_row(
      boundary_->reconstruction(), row.wall_value_row,
      effective_viscosity, 0U, wall_viscosity);
  if (!status || !std::isfinite(wall_viscosity) || wall_viscosity <= 0.0 ||
      !std::isfinite(area) || area <= 0.0)
    return status ? Status{StatusCode::numerical_failure, kIbmEquationNumerical} : status;
  double normal_derivative[3]{};
  const auto* inlet = inlet_for_link(row.topology_link);
  const Real3 boundary_velocity = inlet != nullptr ? inlet->velocity : Real3{};
  const double boundary_components[3]{boundary_velocity.x,
                                      boundary_velocity.y,
                                      boundary_velocity.z};
  for (std::uint8_t component = 0U; component < 3U; ++component) {
    status = evaluate_quadratic_row(
        boundary_->reconstruction(), row.wall_normal_gradient_row,
        velocity, component, boundary_components[component], 0.0,
        normal_derivative[component]);
    if (!status || !std::isfinite(normal_derivative[component]))
      return status ? Status{StatusCode::numerical_failure,
                             kIbmEquationNumerical}
                    : status;
  }
  const double normal_components[3]{link.solid_to_fluid_normal.x,
                                    link.solid_to_fluid_normal.y,
                                    link.solid_to_fluid_normal.z};
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
  if (inlet == nullptr && wall_treatment != nullptr &&
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
  ViscousBoundaryTraction value;
  for (std::size_t c = 0U; c < 3U; ++c) {
    if (!std::isfinite(desired[c]))
      return {StatusCode::numerical_failure, kIbmEquationNumerical};
    value.residual[c] = desired[c];
  }
  value.boundary_velocity = boundary_velocity;
  value.viscosity = wall_viscosity;
  value.wall_drag_coefficient = wall_drag_coefficient;
  value.equilibrium_wall = equilibrium_wall;
  out = value;
  return {};
}

Status IbmEquationInterfacePlan::constrain_momentum(
    ConstFieldView velocity, ConstFieldView velocity_gradient,
    ConstFieldView pressure_perturbation, ConstFieldView density,
    ConstFieldView molecular_viscosity, ConstFieldView effective_viscosity,
    const TurbulencePlan* wall_treatment,
    IbmCellEquationView system) const noexcept {
  Status status =
      validate_bound(*this, kernels_, topology_, boundary_, metric_);
  if (!status) return status;
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
  const bool coast_delta = pressure_gradient_ == IbmPressureGradientKind::coast_fluid_delta;
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
    ViscousBoundaryTraction traction;
    status = viscous_boundary_traction(index, velocity, density,
        molecular_viscosity, effective_viscosity, wall_treatment, traction);
    const double wall_viscosity = traction.viscosity;
    const double area = physical.physical_quadrature_area;
    const Int3 cell = link.fluid_local_index;
    const WallLinearization linearization = wall_linearization_[index];
    double pressure_ghost = 0.0;
    if (status && !coast_delta)
      status = evaluate_quadratic_row(
          boundary_->reconstruction(), row.zero_normal_value_row,
          pressure_perturbation, 0U, 0.0, 0.0, pressure_ghost);
    const double solid_pressure =
        pressure_perturbation.unchecked(link.solid_local_index, 0U);
    const double pressure_gradient_correction =
        coast_delta ? coast_pressure_link_correction(*kernels_, *topology_, link,
                                                     pressure_perturbation) :
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
    const double normal_components[3]{link.solid_to_fluid_normal.x,
                                      link.solid_to_fluid_normal.y,
                                      link.solid_to_fluid_normal.z};
    const double normal_l1 = std::abs(normal_components[0U]) +
                             std::abs(normal_components[1U]) +
                             std::abs(normal_components[2U]);
    const auto& desired = traction.residual;
    const bool equilibrium_wall = traction.equilibrium_wall;
    const double wall_drag_coefficient = traction.wall_drag_coefficient;
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
  if (!status) return status;
  const Int3 cells = kernels_->cells();
  if (!status ||
      !detail::valid_cell_view(pressure, cells, 0U, 1U,
                              boundary_->maximum_halo_reach()) ||
      !detail::valid_cell_view(gradient, cells, 0U, 3U) ||
      detail::field_views_overlap(pressure, as_const(gradient)))
    return status ? Status{StatusCode::invalid_plan, kIbmEquationApply}
                  : status;
  const bool coast_delta = pressure_gradient_ == IbmPressureGradientKind::coast_fluid_delta;
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
    if (!coast_delta)
      status = evaluate_quadratic_row(
          boundary_->reconstruction(), row.zero_normal_value_row, pressure, 0U,
          0.0, 0.0, ghost);
    const double solid =
        pressure.unchecked(link.solid_local_index, 0U);
    const double correction =
        coast_delta ? coast_pressure_link_correction(*kernels_, *topology_, link, pressure) :
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
  if (!status) return status;
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
  const bool coast_delta = pressure_gradient_ == IbmPressureGradientKind::coast_fluid_delta;
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
    if (!coast_delta)
      status = evaluate_quadratic_row(
          boundary_->reconstruction(), row.zero_normal_value_row, pressure, 0U,
          0.0, 0.0, ghost);
    const double solid = pressure.unchecked(link.solid_local_index, 0U);
    const double pressure_gradient_correction =
        coast_delta ? coast_pressure_link_correction(*kernels_, *topology_, link, pressure) :
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
      detail::field_views_overlap(velocity, as_const(velocity_gradient)) ||
      !std::isfinite(wall_velocity.x) || !std::isfinite(wall_velocity.y) ||
      !std::isfinite(wall_velocity.z))
    return {StatusCode::invalid_plan, kIbmEquationApply};
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
    const auto* inlet = inlet_for_link(row.topology_link);
    const Real3 boundary_velocity = inlet != nullptr ? inlet->velocity
                                                    : wall_velocity;
    const double wall_components[3U]{boundary_velocity.x, boundary_velocity.y,
                                     boundary_velocity.z};
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

Status IbmEquationInterfacePlan::add_source_kinetic_convection_correction(
    ConvectionScheme scheme, ConstFieldView velocity, double scale,
    FieldView output, KernelBox box) const noexcept {
  return add_source_convection_correction_impl(
      {IbmInterfaceInletFieldKind::kinetic_energy,0U}, &scheme, velocity,
      scale, output, box, true);
}

Status IbmEquationInterfacePlan::correct_viscous_heating(
    ConstFieldView velocity, ConstFieldView velocity_gradient,
    ConstFieldView density, ConstFieldView molecular_viscosity,
    ConstFieldView effective_viscosity, const TurbulencePlan* wall_treatment,
    FieldView rate, KernelBox box, bool total_energy_work) const noexcept {
  Status status = validate_bound(*this, kernels_, topology_, boundary_, metric_);
  if (!status) return status;
  const auto cells = kernels_->cells();
  if (box.cells.x == 0 && box.cells.y == 0 && box.cells.z == 0) box = {{0,0,0}, cells};
  if (!detail::valid_kernel_box(box, cells) ||
      !detail::valid_cell_view(velocity, cells, 0U, 3U, boundary_->maximum_halo_reach()) ||
      !detail::valid_cell_view(velocity_gradient, cells, 0U, 9U, 1U) ||
      !detail::valid_cell_view(density, cells, 0U, 1U, 0U) ||
      !detail::valid_cell_view(molecular_viscosity, cells, 0U, 1U, 0U) ||
      !detail::valid_cell_view(effective_viscosity, cells, 0U, 1U, boundary_->maximum_halo_reach()) ||
      !detail::valid_cell_view(rate, cells, 0U, 1U))
    return {StatusCode::invalid_plan, kIbmEquationApply};
  for (auto input : {velocity, velocity_gradient, density, molecular_viscosity, effective_viscosity})
    if (detail::field_views_overlap(input, as_const(rate)))
      return {StatusCode::invalid_plan, kIbmEquationApply};
  const auto rows = boundary_->links();
  const auto links = topology_->links();
  for (std::size_t index = 0U; index < rows.size; ++index) {
    if (rows.data[index].topology_link >= links.size)
      return {StatusCode::invalid_plan, kIbmEquationApply};
    const auto& link = links.data[rows.data[index].topology_link];
    const auto cell = link.fluid_local_index;
    if (cell.x < box.begin.x || cell.x >= box.begin.x + box.cells.x ||
        cell.y < box.begin.y || cell.y >= box.begin.y + box.cells.y ||
        cell.z < box.begin.z || cell.z >= box.begin.z + box.cells.z) continue;
    ViscousBoundaryTraction traction;
    status = viscous_boundary_traction(index, velocity, density, molecular_viscosity,
                                       effective_viscosity, wall_treatment, traction);
    if (!status) return status;
    const auto face = interface_face(link);
    const auto left = offset(face.index, face.axis, -1);
    const double boundary_u[3]{traction.boundary_velocity.x,
                               traction.boundary_velocity.y,
                               traction.boundary_velocity.z};
    double correction = 0.0;
    for (std::uint8_t c = 0U; c < 3U; ++c) {
      const double u = total_energy_work ? 0.0 : velocity.unchecked(cell, c);
      const double u_face = detail::interpolate_face(*kernels_, face.axis,
          normal_index(face.index, face.axis), velocity.unchecked(left, c),
          velocity.unchecked(face.index, c));
      const double cartesian = detail::viscous_face_traction_area(
          *kernels_, velocity, velocity_gradient, effective_viscosity, face.axis, face.index, c);
      // desired is an outward residual (-force on fluid). The physical
      // interface work is U_boundary dot force, zero at a stationary wall.
      correction += (u - boundary_u[c]) * traction.residual[c] -
          (positive_face(link.direction) ? 1.0 : -1.0) * (u_face - u) * cartesian;
    }
    const double value = rate.unchecked(cell, 0U) +
                         correction / detail::cell_volume(*kernels_, cell);
    if (!std::isfinite(value)) return {StatusCode::numerical_failure, kIbmEquationNumerical};
    rate.unchecked(cell, 0U) = value;
  }
  return {};
}

Status IbmEquationInterfacePlan::inlet_viscous_work_input(
    ConstFieldView velocity, ConstFieldView effective_viscosity,
    double& input) const noexcept {
  Status status = validate_bound(*this, kernels_, topology_, boundary_, metric_);
  if (!status) return status;
  if (!has_inlet_sources()) {
    input = 0.0;
    return {};
  }
  const Int3 cells = kernels_->cells();
  if (!detail::valid_cell_view(velocity, cells, 0U, 3U,
                               boundary_->maximum_halo_reach()) ||
      !detail::valid_cell_view(effective_viscosity, cells, 0U, 1U,
                               boundary_->maximum_halo_reach()))
    return {StatusCode::invalid_plan, kIbmEquationApply};
  long double work = 0.0L;
  const auto rows = boundary_->links();
  for (std::size_t index = 0U; index < rows.size; ++index) {
    const auto& row = rows.data[index];
    const auto* inlet = inlet_for_link(row.topology_link);
    if (inlet == nullptr) continue;
    const double u[3]{inlet->velocity.x, inlet->velocity.y, inlet->velocity.z};
    ViscousBoundaryTraction traction;
    status = viscous_boundary_traction(index, velocity, {}, {}, effective_viscosity,
                                       nullptr, traction);
    if (!status) return status;
    for (std::uint8_t c = 0U; c < 3U; ++c) {
      // constrain_momentum inserts positive outward viscous residual;
      // physical traction/work into the fluid has the opposite sign.
      work -= u[c] * traction.residual[c];
    }
  }
  if (!std::isfinite(static_cast<double>(work)))
    return {StatusCode::numerical_failure, kIbmEquationNumerical};
  input = static_cast<double>(work);
  return {};
}

Status IbmEquationInterfacePlan::correct_zero_normal_diffusion(
    ConstFieldView transported, ConstFieldView diffusivity,
    FieldView rate) const noexcept {
  Status status =
      validate_bound(*this, kernels_, topology_, boundary_, metric_);
  if (!status) return status;
  const Int3 cells = kernels_->cells();
  if (!status ||
      !detail::valid_cell_view(transported, cells, 0U, 1U,
                              boundary_->maximum_halo_reach()) ||
      !detail::valid_cell_view(diffusivity, cells, 0U, 1U, 1U) ||
      !detail::valid_cell_view(rate, cells, 0U, 1U) ||
      detail::field_views_overlap(transported, as_const(rate)) ||
      detail::field_views_overlap(diffusivity, as_const(rate)))
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
    if (inlet_for_link(row.topology_link) != nullptr)
      ghost = transported.unchecked(link.fluid_local_index, 0U);
    else
      status = evaluate_quadratic_row(
          boundary_->reconstruction(), row.zero_normal_value_row, transported,
          0U, 0.0, 0.0, ghost);
    const double solid = transported.unchecked(link.solid_local_index, 0U);
    const double fluid = transported.unchecked(link.fluid_local_index, 0U);
    // The prescribed state sets advective h/Y. The inlet's diffusive thermal
    // and species flux is explicitly zero, not a fitted solid-wall ghost.
    const double fluid_conductance =
        fluid_material_conductance(*kernels_, diffusivity, link, face);
    // Remove the Cartesian solid-material face, then insert the existing
    // link reconstruction with fluid-side material authority. An adiabatic,
    // non-conjugate placeholder is not a second material in series.
    const double correction =
        (transmissibility * (fluid - solid) +
         fluid_conductance * (ghost - fluid)) / volume;
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

double detail::IbmScalarTransport::diffusion_diagonal(
    const IbmEquationInterfacePlan& plan,ConstFieldView gamma,Int3 cell) noexcept {
  const auto& kernels=*plan.kernels_;
  const Int3 n=kernels.cells();
  const auto i=std::size_t(cell.x)+std::size_t(n.x)*(std::size_t(cell.y)+std::size_t(n.y)*cell.z);
  const auto blocked=plan.scalar_blocked_faces_[i];
  if(blocked==0U) return detail::diffusion_diagonal(kernels,gamma,cell);
  double sum=0.0;
  for(unsigned d=0U;d<6U;++d) {
    if((blocked & (1U << d))!=0U) continue;
    const auto axis=static_cast<CartesianAxis>(d/2U);
    Int3 face=cell;
    if(d%2U) (axis==CartesianAxis::x ? face.x : axis==CartesianAxis::y ? face.y : face.z)++;
    sum+=detail::positive_transmissibility(kernels,gamma,axis,face);
  }
  return sum;
}

Status detail::IbmScalarTransport::constrain_rows(
    const IbmEquationInterfacePlan& plan,ConstFieldView q,KernelBox box,
    EquationSystemView system) noexcept {
  const auto bound=validate_bound(plan,plan.kernels_,plan.topology_,plan.boundary_,plan.metric_);
  if(!bound) return bound;
  const Int3 cells=plan.kernels_->cells();
  if(!detail::valid_kernel_box(box,cells) || !detail::valid_cell_view(q,cells,0U,1U,0U) ||
      !detail::valid_cell_view(system.diagonal,cells,0U,1U) ||
      !detail::valid_cell_view(system.rhs,cells,0U,1U) ||
      !detail::valid_cell_view(system.residual,cells,0U,1U))
    return {StatusCode::invalid_plan,kIbmEquationApply};
  const auto region=plan.topology_->region();
  const std::array<FaceFieldView,3U> faces{system.x_coefficient,system.y_coefficient,system.z_coefficient};
  for(int z=box.begin.z;z<box.begin.z+box.cells.z;++z)
    for(int y=box.begin.y;y<box.begin.y+box.cells.y;++y)
      for(int x=box.begin.x;x<box.begin.x+box.cells.x;++x) {
        const Int3 c{x,y,z};
        const auto i=std::size_t(x)+std::size_t(cells.x)*(std::size_t(y)+std::size_t(cells.y)*z);
        if(region.data[i]==0U) {
          system.diagonal.unchecked(c,0U)=1.0;
          system.rhs.unchecked(c,0U)=q.unchecked(c,0U);
          system.residual.unchecked(c,0U)=0.0;
        }
        const auto blocked=plan.scalar_blocked_faces_[i];
        for(unsigned d=0U;d<6U;++d) {
          if((blocked & (1U<<d))==0U || faces[d/2U].base==nullptr) continue;
          Int3 face=c;
          if(d%2U) (d/2U==0U ? face.x : d/2U==1U ? face.y : face.z)++;
          faces[d/2U].unchecked(face)=0.0;
        }
      }
  return {};
}

Status detail::IbmScalarTransport::convection(const IbmEquationInterfacePlan& plan,
    std::size_t species,ConvectionScheme scheme,ConstFieldView q,
    ConstFaceFluxView flux,KernelBox box,FieldView rate) noexcept {
  Status status=validate_bound(plan,plan.kernels_,plan.topology_,plan.boundary_,plan.metric_);
  if(!status) return status;
  const auto& kernels=*plan.kernels_; const Int3 cells=kernels.cells();
  if(!detail::valid_kernel_box(box,cells) ||
     static_cast<unsigned>(scheme)>static_cast<unsigned>(ConvectionScheme::tvd2) ||
     !detail::valid_cell_view(q,cells,0U,1U,scheme==ConvectionScheme::central2 ? 1U : 2U) ||
     !detail::valid_cell_view(rate,cells,0U,1U) || detail::field_views_overlap(q,rate) ||
     detail::cell_face_views_overlap(rate,flux.x) ||
     detail::cell_face_views_overlap(rate,flux.y) ||
     detail::cell_face_views_overlap(rate,flux.z))
    return {StatusCode::invalid_plan,kIbmEquationApply};
  status=plan.validate_interface_flux(flux);
  if(!status) return status;
  if(!plan.prescribed_interface_fluxes_.empty() &&
     (!plan.inlet_state_bound_ || species>=plan.independent_species_count_))
    return {StatusCode::invalid_plan,kIbmEquationApply};
  const auto links=plan.topology_->links();
  const auto physical_row=[&](Int3 c,const std::array<std::size_t,6U>& cut,double& value) noexcept -> Status {
    long double sum=0.0L;
    for(std::size_t d=0U;d<6U;++d) {
      const auto axis=static_cast<CartesianAxis>(d/2U); const bool positive=d%2U!=0U;
      Int3 face=c;
      if(positive) (axis==CartesianAxis::x ? face.x : axis==CartesianAxis::y ? face.y : face.z)++;
      const double mass=detail::select(flux,axis).unchecked(face);
      if(mass==0.0) continue;
      long double face_q=0.0L;
      if(cut[d]<links.size) {
        const auto* source=plan.inlet_for_link(static_cast<std::uint32_t>(cut[d]));
        if(source==nullptr || !source->has_inlet_state ||
           source->independent_species_begin>=plan.prescribed_independent_species_.size() ||
           species>=plan.prescribed_independent_species_.size()-source->independent_species_begin)
          return {StatusCode::invalid_plan,kIbmEquationApply};
        face_q=plan.prescribed_independent_species_[source->independent_species_begin+species];
      } else {
        face_q=detail::precise_scalar_convection_face(kernels,scheme,q,axis,face,mass);
      }
      sum+=(positive ? 1.0L : -1.0L)*mass*face_q;
    }
    value=static_cast<double>(sum/detail::cell_volume(kernels,c));
    return std::isfinite(value) ? Status{} : Status{StatusCode::numerical_failure,kIbmEquationNumerical};
  };
  return replace_cut_scalar_rows(links,cells,box,rate,physical_row);
}

Status detail::IbmScalarTransport::diffusion(const IbmEquationInterfacePlan& plan,
    ConstFieldView q,ConstFieldView gamma,KernelBox box,FieldView rate) noexcept {
  Status status=validate_bound(plan,plan.kernels_,plan.topology_,plan.boundary_,plan.metric_);
  if(!status) return status;
  const auto& kernels=*plan.kernels_; const Int3 cells=kernels.cells();
  if(!detail::valid_kernel_box(box,cells) || !detail::valid_cell_view(q,cells,0U,1U,1U) ||
     !detail::valid_cell_view(gamma,cells,0U,1U,1U) || !detail::valid_cell_view(rate,cells,0U,1U) ||
     detail::field_views_overlap(q,rate) || detail::field_views_overlap(gamma,rate))
    return {StatusCode::invalid_plan,kIbmEquationApply};
  const auto links=plan.topology_->links();
  const auto physical_row=[&](Int3 c,const std::array<std::size_t,6U>& cut,double& value) noexcept -> Status {
    long double sum=0.0L;
    for(std::size_t d=0U;d<6U;++d) {
      if(cut[d]<links.size) continue; // No solid placeholder enters this equation.
      const auto axis=static_cast<CartesianAxis>(d/2U); const bool positive=d%2U!=0U;
      Int3 neighbour=c;
      (axis==CartesianAxis::x ? neighbour.x : axis==CartesianAxis::y ? neighbour.y : neighbour.z)+=positive ? 1 : -1;
      const Int3 face=positive ? neighbour : c;
      const double transmissibility=detail::positive_transmissibility(kernels,gamma,axis,face);
      const double jump=q.unchecked(neighbour,0U)-q.unchecked(c,0U);
      if(!std::isfinite(transmissibility) || transmissibility<=0.0 || !std::isfinite(jump))
        return {StatusCode::numerical_failure,kIbmEquationNumerical};
      sum+=static_cast<long double>(transmissibility)*jump;
    }
    value=static_cast<double>(sum/detail::cell_volume(kernels,c));
    return std::isfinite(value) ? Status{} : Status{StatusCode::numerical_failure,kIbmEquationNumerical};
  };
  status=replace_cut_scalar_rows(links,cells,box,rate,physical_row);
  if(!status) return status;
  const auto region=plan.topology_->region();
  for(int z=box.begin.z;z<box.begin.z+box.cells.z;++z)
    for(int y=box.begin.y;y<box.begin.y+box.cells.y;++y)
      for(int x=box.begin.x;x<box.begin.x+box.cells.x;++x) {
        const std::size_t flat=std::size_t(x)+std::size_t(cells.x)*(std::size_t(y)+std::size_t(cells.y)*z);
        if(region.data[flat]==static_cast<std::uint8_t>(RegionFlag::solid)) rate.unchecked({x,y,z},0U)=0.0;
      }
  return {};
}

Status IbmEquationInterfacePlan::correct_impermeable_scalar_diffusion(
    ConstFieldView transported, ConstFieldView diffusivity,
    FieldView rate) const noexcept {
  Status status = validate_bound(*this, kernels_, topology_, boundary_, metric_);
  if (!status) return status;
  const Int3 cells = kernels_->cells();
  if (!detail::valid_cell_view(transported, cells, 0U, 1U, 1U) ||
      !detail::valid_cell_view(diffusivity, cells, 0U, 1U, 1U) ||
      !detail::valid_cell_view(rate, cells, 0U, 1U) ||
      detail::field_views_overlap(transported, as_const(rate)) ||
      detail::field_views_overlap(diffusivity, as_const(rate)))
    return {StatusCode::invalid_plan, kIbmEquationApply};
  const auto links = topology_->links();
  for (std::size_t index = 0U; index < links.size; ++index) {
    const auto& link = links.data[index];
    const InterfaceFace face = interface_face(link);
    const double transmissibility = detail::positive_transmissibility(
        *kernels_, diffusivity, face.axis, face.index);
    const double volume = detail::cell_volume(*kernels_, link.fluid_local_index);
    // Remove exactly the cut-face term included by cartesian_diffusion.
    // Fluid-fluid face pairs are untouched and cancel in the global ledger.
    const double jump =
        (transported.unchecked(link.fluid_local_index, 0U) -
         transported.unchecked(link.solid_local_index, 0U));
    const double product=transmissibility*jump;
    const double correction = jump!=0.0 && std::abs(product)<std::numeric_limits<double>::min()
        ? static_cast<double>(static_cast<long double>(transmissibility)*jump*(1.0/volume))
        : product/volume;
    const double value = rate.unchecked(link.fluid_local_index, 0U) + correction;
    if (!std::isfinite(transmissibility) || transmissibility <= 0.0 ||
        !std::isfinite(volume) || volume <= 0.0 || !std::isfinite(value))
      return {StatusCode::numerical_failure, kIbmEquationNumerical};
    rate.unchecked(link.fluid_local_index, 0U) = value;
  }
  const auto region = topology_->region();
  std::size_t flat = 0U;
  for (int z = 0; z < cells.z; ++z)
    for (int y = 0; y < cells.y; ++y)
      for (int x = 0; x < cells.x; ++x, ++flat)
        if (region.data[flat] == static_cast<std::uint8_t>(RegionFlag::solid))
          rate.unchecked({x,y,z},0U) = 0.0;
  return {};
}

Status IbmEquationInterfacePlan::correct_positive_bounded_zero_normal_diffusion(
    ConstFieldView transported, ConstFieldView diffusivity,
    FieldView rate) const noexcept {
  Status status =
      validate_bound(*this, kernels_, topology_, boundary_, metric_);
  if (!status) return status;
  const Int3 cells = kernels_->cells();
  if (!status ||
      !detail::valid_cell_view(transported, cells, 0U, 1U,
                               boundary_->maximum_halo_reach()) ||
      !detail::valid_cell_view(diffusivity, cells, 0U, 1U, 1U) ||
      !detail::valid_cell_view(rate, cells, 0U, 1U) ||
      detail::field_views_overlap(transported, as_const(rate)) ||
      detail::field_views_overlap(diffusivity, as_const(rate)))
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
    if (inlet_for_link(row.topology_link) != nullptr)
      ghost = transported.unchecked(link.fluid_local_index, 0U);
    else
      status = evaluate_positive_bounded_quadratic_row(
          boundary_->reconstruction(), row.zero_normal_value_row, transported,
          0U, ghost);
    const double solid = transported.unchecked(link.solid_local_index, 0U);
    const double current = rate.unchecked(link.fluid_local_index, 0U);
    const double fluid = transported.unchecked(link.fluid_local_index, 0U);
    const double fluid_conductance =
        fluid_material_conductance(*kernels_, diffusivity, link, face);
    const double correction =
        (transmissibility * (fluid - solid) +
         fluid_conductance * (ghost - fluid)) / volume;
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

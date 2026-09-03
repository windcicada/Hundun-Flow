// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_ibm.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <type_traits>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kTopologyInput = 13501U;
constexpr std::uint32_t kTopologyCollective = 13502U;
constexpr std::uint32_t kTopologyRegion = 13503U;
constexpr std::uint32_t kTopologyGeometry = 13504U;
constexpr std::uint32_t kTopologyLinks = 13505U;
constexpr std::uint32_t kTopologyAllocation = 13506U;
constexpr std::uint8_t kRegionHalo = 4U;
constexpr std::uint32_t kTopologySchemaRevision = 2U;
constexpr int kHaloTagBase = 17040;
constexpr std::uint64_t kFnvOffset = UINT64_C(14695981039346656037);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);

static_assert(std::is_nothrow_move_assignable_v<EBTopology>);

class Hash64 {
 public:
  template <class Integer>
  void integer(Integer value) noexcept {
    using Unsigned = std::make_unsigned_t<Integer>;
    const Unsigned bits = static_cast<Unsigned>(value);
    for (std::size_t shift = 0U; shift < sizeof(bits) * 8U; shift += 8U) {
      byte(static_cast<std::uint8_t>((bits >> shift) & Unsigned{0xffU}));
    }
  }
  void real(double value) noexcept {
    std::uint64_t bits{};
    std::memcpy(&bits, &value, sizeof(bits));
    integer(bits);
  }
  PlanFingerprint finish() const noexcept { return value_ == 0U ? 1U : value_; }

 private:
  void byte(std::uint8_t value) noexcept {
    value_ ^= value;
    value_ *= kFnvPrime;
  }
  std::uint64_t value_{kFnvOffset};
};

bool mpi_live() noexcept {
  int initialized = 0;
  int finalized = 0;
  (void)MPI_Initialized(&initialized);
  if (initialized != 0) {
    (void)MPI_Finalized(&finalized);
  }
  return initialized != 0 && finalized == 0;
}

Status consensus(MPI_Comm communicator, int rank, int size, Status local,
                 int& lowest) noexcept {
  const int candidate = local ? size : rank;
  int selected = size;
  if (MPI_Allreduce(&candidate, &selected, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    lowest = -1;
    return {StatusCode::mpi_failure, kTopologyCollective};
  }
  if (selected == size) {
    lowest = -1;
    return {};
  }
  std::array<std::uint64_t, 2U> wire{};
  if (rank == selected) {
    wire = {static_cast<std::uint64_t>(local.code), local.detail};
  }
  if (selected < 0 || selected >= size ||
      MPI_Bcast(wire.data(), static_cast<int>(wire.size()), MPI_UINT64_T,
                selected, communicator) != MPI_SUCCESS ||
      wire[0] > static_cast<std::uint64_t>(StatusCode::io_failure) ||
      wire[1] > std::numeric_limits<std::uint32_t>::max()) {
    lowest = -1;
    return {StatusCode::mpi_failure, kTopologyCollective};
  }
  lowest = selected;
  return {static_cast<StatusCode>(wire[0]),
          static_cast<std::uint32_t>(wire[1])};
}

bool checked_product(Int3 cells, std::uint64_t& out) noexcept {
  if (cells.x <= 0 || cells.y <= 0 || cells.z <= 0) {
    return false;
  }
  const std::uint64_t x = static_cast<std::uint64_t>(cells.x);
  const std::uint64_t y = static_cast<std::uint64_t>(cells.y);
  const std::uint64_t z = static_cast<std::uint64_t>(cells.z);
  if (x > std::numeric_limits<std::uint64_t>::max() / y ||
      x * y > std::numeric_limits<std::uint64_t>::max() / z) {
    return false;
  }
  out = x * y * z;
  return true;
}

bool checked_add(std::uint64_t left, std::uint64_t right,
                 std::uint64_t& out) noexcept {
  if (right > std::numeric_limits<std::uint64_t>::max() - left) {
    return false;
  }
  out = left + right;
  return true;
}

bool checked_multiply(std::uint64_t left, std::uint64_t right,
                      std::uint64_t& out) noexcept {
  if (left != 0U &&
      right > std::numeric_limits<std::uint64_t>::max() / left) {
    return false;
  }
  out = left * right;
  return true;
}

void split_axis(std::int32_t global, std::int32_t partitions,
                std::int32_t coordinate, std::int32_t& begin,
                std::int32_t& cells) noexcept {
  const std::int32_t quotient = global / partitions;
  const std::int32_t remainder = global % partitions;
  cells = quotient + (coordinate < remainder ? 1 : 0);
  begin = coordinate * quotient + std::min(coordinate, remainder);
}

bool valid_patch(Int3 global, const MeshPatch& patch, int rank,
                 int size) noexcept {
  if (patch.process_grid.x <= 0 || patch.process_grid.y <= 0 ||
      patch.process_grid.z <= 0 || patch.process_coord.x < 0 ||
      patch.process_coord.y < 0 || patch.process_coord.z < 0 ||
      patch.process_coord.x >= patch.process_grid.x ||
      patch.process_coord.y >= patch.process_grid.y ||
      patch.process_coord.z >= patch.process_grid.z) {
    return false;
  }
  std::uint64_t count{};
  if (!checked_product(patch.process_grid, count) ||
      count != static_cast<std::uint64_t>(size)) {
    return false;
  }
  const std::int64_t expected_rank =
      patch.process_coord.x +
      static_cast<std::int64_t>(patch.process_grid.x) *
          (patch.process_coord.y +
           static_cast<std::int64_t>(patch.process_grid.y) *
               patch.process_coord.z);
  if (expected_rank != rank) {
    return false;
  }
  Int3 begin{};
  Int3 cells{};
  split_axis(global.x, patch.process_grid.x, patch.process_coord.x, begin.x,
             cells.x);
  split_axis(global.y, patch.process_grid.y, patch.process_coord.y, begin.y,
             cells.y);
  split_axis(global.z, patch.process_grid.z, patch.process_coord.z, begin.z,
             cells.z);
  return patch.begin.x == begin.x && patch.begin.y == begin.y &&
         patch.begin.z == begin.z && patch.cells.x == cells.x &&
         patch.cells.y == cells.y && patch.cells.z == cells.z;
}

bool inside(Int3 cell, Int3 global) noexcept {
  return cell.x >= 0 && cell.y >= 0 && cell.z >= 0 &&
         cell.x < global.x && cell.y < global.y && cell.z < global.z;
}

std::size_t flat(Int3 shape, Int3 cell) noexcept {
  return static_cast<std::size_t>(cell.x) +
         static_cast<std::size_t>(shape.x) *
             (static_cast<std::size_t>(cell.y) +
              static_cast<std::size_t>(shape.y) *
                  static_cast<std::size_t>(cell.z));
}

GlobalCellId global_id(Int3 cell, Int3 global) noexcept {
  return static_cast<GlobalCellId>(cell.x) +
         static_cast<GlobalCellId>(global.x) *
             (static_cast<GlobalCellId>(cell.y) +
              static_cast<GlobalCellId>(global.y) *
                  static_cast<GlobalCellId>(cell.z));
}

Real3 cell_centre(const CartesianGeometryPlan& geometry, Int3 cell) noexcept {
  return {geometry.x().centres().data[static_cast<std::size_t>(cell.x)],
          geometry.y().centres().data[static_cast<std::size_t>(cell.y)],
          geometry.z().centres().data[static_cast<std::size_t>(cell.z)]};
}

double width(const CartesianGeometryPlan& geometry, Int3 cell,
             std::size_t axis) noexcept {
  const AxisMetrics& metrics = geometry.axis(static_cast<CartesianAxis>(axis));
  const std::int32_t index =
      axis == 0U ? cell.x : (axis == 1U ? cell.y : cell.z);
  return metrics.widths().data[static_cast<std::size_t>(index)];
}

bool selected_fluid(std::uint8_t raw, ImmersedFluidSide side) noexcept {
  const bool outside = raw == static_cast<std::uint8_t>(RegionFlag::fluid);
  return side == ImmersedFluidSide::outside ? outside : !outside;
}

std::uint8_t selected_region(std::uint8_t raw,
                             ImmersedFluidSide side) noexcept {
  return static_cast<std::uint8_t>(selected_fluid(raw, side)
                                       ? RegionFlag::fluid
                                       : RegionFlag::solid);
}

constexpr std::array<Int3, 6U> kDirectionDelta{{
    {-1, 0, 0}, {1, 0, 0}, {0, -1, 0},
    {0, 1, 0},  {0, 0, -1}, {0, 0, 1},
}};

double face_area(const CartesianGeometryPlan& geometry, Int3 fluid,
                 std::size_t direction) noexcept {
  const std::size_t axis = direction / 2U;
  if (axis == 0U) {
    return width(geometry, fluid, 1U) * width(geometry, fluid, 2U);
  }
  if (axis == 1U) {
    return width(geometry, fluid, 0U) * width(geometry, fluid, 2U);
  }
  return width(geometry, fluid, 0U) * width(geometry, fluid, 1U);
}

int neighbor_rank(const MeshPatch& patch, int axis, int sign) noexcept {
  Int3 coordinate = patch.process_coord;
  std::int32_t* selected =
      axis == 0 ? &coordinate.x : (axis == 1 ? &coordinate.y : &coordinate.z);
  const std::int32_t extent =
      axis == 0 ? patch.process_grid.x
                : (axis == 1 ? patch.process_grid.y : patch.process_grid.z);
  *selected += sign;
  if (*selected < 0 || *selected >= extent) {
    return MPI_PROC_NULL;
  }
  return coordinate.x + patch.process_grid.x *
                            (coordinate.y + patch.process_grid.y * coordinate.z);
}

bool exchange_buffer_elements(const MeshPatch& patch, Int3 global, int axis,
                              Int3 shape, std::uint64_t& out) noexcept {
  const int h = static_cast<int>(kRegionHalo);
  const std::int32_t global_extent =
      axis == 0 ? global.x : (axis == 1 ? global.y : global.z);
  const std::int32_t process_extent =
      axis == 0 ? patch.process_grid.x
                : (axis == 1 ? patch.process_grid.y : patch.process_grid.z);
  if (process_extent == 1) {
    out = 0U;
    return true;
  }
  Int3 count{};
  const bool narrow = global_extent / process_extent < h;
  if (axis == 0) {
    count = {narrow ? 1 : h, patch.cells.y, patch.cells.z};
  } else if (axis == 1) {
    count = {shape.x, narrow ? 1 : h, patch.cells.z};
  } else {
    count = {shape.x, shape.y, narrow ? 1 : h};
  }
  return checked_product(count, out) &&
         out <= static_cast<std::uint64_t>(std::numeric_limits<int>::max());
}

std::uint8_t& halo_at(std::vector<std::uint8_t>& halo, Int3 shape,
                      Int3 index) noexcept {
  return halo[flat(shape, index)];
}

Status exchange_axis(MPI_Comm communicator, const MeshPatch& patch,
                     Int3 global, int axis, Int3 shape,
                     std::vector<std::uint8_t>& halo,
                     std::vector<std::uint8_t>& send_lower,
                     std::vector<std::uint8_t>& send_upper,
                     std::vector<std::uint8_t>& receive_lower,
                     std::vector<std::uint8_t>& receive_upper) {
  const int h = static_cast<int>(kRegionHalo);
  const int lower = neighbor_rank(patch, axis, -1);
  const int upper = neighbor_rank(patch, axis, 1);
  const std::int32_t local_extent =
      axis == 0 ? patch.cells.x : (axis == 1 ? patch.cells.y : patch.cells.z);
  const std::int32_t global_extent =
      axis == 0 ? global.x : (axis == 1 ? global.y : global.z);
  const std::int32_t process_extent =
      axis == 0 ? patch.process_grid.x
                : (axis == 1 ? patch.process_grid.y : patch.process_grid.z);
  if (process_extent == 1) {
    return {};
  }
  Int3 begin{0, 0, 0};
  Int3 count = shape;

  // A width-four slab is the common fast path.  When a legal decomposition
  // owns fewer than four cells on this axis, exchange one plane at a time so
  // data can be forwarded across more than one neighboring rank.  The branch
  // is collective because every rank derives it from the same global shape
  // and process grid.
  if (global_extent / process_extent < h) {
    if (axis == 0) {
      begin.y = h;
      begin.z = h;
      count = {1, patch.cells.y, patch.cells.z};
    } else if (axis == 1) {
      begin.z = h;
      count = {shape.x, 1, patch.cells.z};
    } else {
      count = {shape.x, shape.y, 1};
    }
    std::uint64_t elements64{};
    if (!checked_product(count, elements64) ||
        elements64 >
            static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
      return {StatusCode::invalid_plan, kTopologyRegion};
    }
    const std::size_t elements = static_cast<std::size_t>(elements64);
    if (send_lower.size() < elements || send_upper.size() < elements ||
        receive_lower.size() < elements || receive_upper.size() < elements) {
      return {StatusCode::invalid_plan, kTopologyRegion};
    }
    const int count_i = static_cast<int>(elements);
    for (int depth = 0; depth < h; ++depth) {
      std::size_t write = 0U;
      for (std::int32_t z = 0; z < count.z; ++z) {
        for (std::int32_t y = 0; y < count.y; ++y) {
          for (std::int32_t x = 0; x < count.x; ++x) {
            Int3 lower_index{begin.x + x, begin.y + y, begin.z + z};
            Int3 upper_index = lower_index;
            if (axis == 0) {
              lower_index.x = h + depth;
              upper_index.x = h + local_extent - 1 - depth;
            } else if (axis == 1) {
              lower_index.y = h + depth;
              upper_index.y = h + local_extent - 1 - depth;
            } else {
              lower_index.z = h + depth;
              upper_index.z = h + local_extent - 1 - depth;
            }
            send_lower[write] = halo_at(halo, shape, lower_index);
            send_upper[write] = halo_at(halo, shape, upper_index);
            ++write;
          }
        }
      }
      const int tag = kHaloTagBase + 8 * axis + 2 * depth;
      const int lower_exchange = MPI_Sendrecv(
          send_lower.data(), count_i, MPI_UINT8_T, lower, tag,
          receive_upper.data(), count_i, MPI_UINT8_T, upper, tag, communicator,
          MPI_STATUS_IGNORE);
      const int upper_exchange = MPI_Sendrecv(
          send_upper.data(), count_i, MPI_UINT8_T, upper, tag + 1,
          receive_lower.data(), count_i, MPI_UINT8_T, lower, tag + 1,
          communicator, MPI_STATUS_IGNORE);
      if (lower_exchange != MPI_SUCCESS || upper_exchange != MPI_SUCCESS) {
        return {StatusCode::mpi_failure, kTopologyCollective};
      }
      write = 0U;
      for (std::int32_t z = 0; z < count.z; ++z) {
        for (std::int32_t y = 0; y < count.y; ++y) {
          for (std::int32_t x = 0; x < count.x; ++x) {
            Int3 lower_index{begin.x + x, begin.y + y, begin.z + z};
            Int3 upper_index = lower_index;
            if (axis == 0) {
              lower_index.x = h - 1 - depth;
              upper_index.x = h + local_extent + depth;
            } else if (axis == 1) {
              lower_index.y = h - 1 - depth;
              upper_index.y = h + local_extent + depth;
            } else {
              lower_index.z = h - 1 - depth;
              upper_index.z = h + local_extent + depth;
            }
            if (lower != MPI_PROC_NULL) {
              halo_at(halo, shape, lower_index) = receive_lower[write];
            }
            if (upper != MPI_PROC_NULL) {
              halo_at(halo, shape, upper_index) = receive_upper[write];
            }
            ++write;
          }
        }
      }
    }
    return {};
  }

  if (axis == 0) {
    begin.y = h;
    begin.z = h;
    count = {h, patch.cells.y, patch.cells.z};
  } else if (axis == 1) {
    begin.z = h;
    count = {shape.x, h, patch.cells.z};
  } else {
    count = {shape.x, shape.y, h};
  }
  std::uint64_t elements64{};
  if (!checked_product(count, elements64) ||
      elements64 > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
    return {StatusCode::invalid_plan, kTopologyRegion};
  }
  const std::size_t elements = static_cast<std::size_t>(elements64);
  if (send_lower.size() < elements || send_upper.size() < elements ||
      receive_lower.size() < elements || receive_upper.size() < elements) {
    return {StatusCode::invalid_plan, kTopologyRegion};
  }
  const auto plane_coordinate = [&](bool upper_plane, int depth) {
    return h + (upper_plane ? local_extent - h + depth : depth);
  };
  std::size_t write = 0U;
  for (std::int32_t z = 0; z < count.z; ++z) {
    for (std::int32_t y = 0; y < count.y; ++y) {
      for (std::int32_t x = 0; x < count.x; ++x) {
        Int3 lower_index{begin.x + x, begin.y + y, begin.z + z};
        Int3 upper_index = lower_index;
        const int depth = axis == 0 ? x : (axis == 1 ? y : z);
        if (axis == 0) {
          lower_index.x = plane_coordinate(false, depth);
          upper_index.x = plane_coordinate(true, depth);
        } else if (axis == 1) {
          lower_index.y = plane_coordinate(false, depth);
          upper_index.y = plane_coordinate(true, depth);
        } else {
          lower_index.z = plane_coordinate(false, depth);
          upper_index.z = plane_coordinate(true, depth);
        }
        send_lower[write] = halo_at(halo, shape, lower_index);
        send_upper[write] = halo_at(halo, shape, upper_index);
        ++write;
      }
    }
  }
  const int count_i = static_cast<int>(elements);
  const int lower_exchange = MPI_Sendrecv(
      send_lower.data(), count_i, MPI_UINT8_T, lower,
      kHaloTagBase + 2 * axis, receive_upper.data(), count_i, MPI_UINT8_T,
      upper, kHaloTagBase + 2 * axis, communicator, MPI_STATUS_IGNORE);
  const int upper_exchange = MPI_Sendrecv(
      send_upper.data(), count_i, MPI_UINT8_T, upper,
      kHaloTagBase + 2 * axis + 1, receive_lower.data(), count_i, MPI_UINT8_T,
      lower, kHaloTagBase + 2 * axis + 1, communicator, MPI_STATUS_IGNORE);
  if (lower_exchange != MPI_SUCCESS || upper_exchange != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kTopologyCollective};
  }
  write = 0U;
  for (std::int32_t z = 0; z < count.z; ++z) {
    for (std::int32_t y = 0; y < count.y; ++y) {
      for (std::int32_t x = 0; x < count.x; ++x) {
        Int3 lower_index{begin.x + x, begin.y + y, begin.z + z};
        Int3 upper_index = lower_index;
        const int depth = axis == 0 ? x : (axis == 1 ? y : z);
        if (axis == 0) {
          lower_index.x = depth;
          upper_index.x = h + local_extent + depth;
        } else if (axis == 1) {
          lower_index.y = depth;
          upper_index.y = h + local_extent + depth;
        } else {
          lower_index.z = depth;
          upper_index.z = h + local_extent + depth;
        }
        if (lower != MPI_PROC_NULL) {
          halo_at(halo, shape, lower_index) = receive_lower[write];
        }
        if (upper != MPI_PROC_NULL) {
          halo_at(halo, shape, upper_index) = receive_upper[write];
        }
        ++write;
      }
    }
  }
  return {};
}

std::uint64_t mix_cell(GlobalCellId id, std::uint8_t region) noexcept {
  std::uint64_t value = id ^ (static_cast<std::uint64_t>(region) << 61U) ^
                        UINT64_C(0x9e3779b97f4a7c15);
  value ^= value >> 30U;
  value *= UINT64_C(0xbf58476d1ce4e5b9);
  value ^= value >> 27U;
  value *= UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31U);
}

PlanFingerprint topology_fingerprint(const CartesianGeometryPlan& geometry,
                                     const ImmersedSurfacePlan& surface,
                                     ImmersedFluidSide side,
                                     std::uint64_t global_count,
                                     std::uint64_t region_xor,
                                     std::uint64_t region_sum) noexcept {
  Hash64 hash;
  hash.integer(geometry.fingerprint());
  hash.integer(surface.fingerprint());
  hash.integer(kTopologySchemaRevision);
  hash.integer(kRegionHalo);
  hash.integer(static_cast<std::uint8_t>(side));
  hash.integer(global_count);
  hash.integer(region_xor);
  hash.integer(region_sum);
  return hash.finish();
}

PlanFingerprint compile_contract(const CartesianGeometryPlan& geometry,
                                 const MeshPatch& patch,
                                 const ImmersedSurfacePlan& surface,
                                 ImmersedFluidSide side,
                                 ImmersedPlanLimits limits) noexcept {
  Hash64 hash;
  hash.integer(geometry.fingerprint());
  hash.integer(geometry.topology_revision());
  hash.integer(surface.fingerprint());
  hash.integer(static_cast<std::uint8_t>(side));
  hash.integer(patch.process_grid.x);
  hash.integer(patch.process_grid.y);
  hash.integer(patch.process_grid.z);
  hash.integer(limits.maximum_persistent_bytes_per_rank);
  hash.integer(limits.maximum_peak_bytes_per_rank);
  hash.integer(limits.maximum_local_links);
  return hash.finish();
}

}  // namespace

bool EBTopology::is_fluid_global(Int3 global_index) const noexcept {
  if (fingerprint_ == 0U || !inside(global_index, global_cells_)) {
    return false;
  }
  const int h = static_cast<int>(region_halo_width_);
  const Int3 local{global_index.x - patch_.begin.x,
                   global_index.y - patch_.begin.y,
                   global_index.z - patch_.begin.z};
  if (local.x < -h || local.y < -h || local.z < -h ||
      local.x >= patch_.cells.x + h || local.y >= patch_.cells.y + h ||
      local.z >= patch_.cells.z + h || halo_stride_y_ == 0U ||
      halo_stride_z_ == 0U) {
    return false;
  }
  const std::size_t offset =
      static_cast<std::size_t>(local.x + h) +
      halo_stride_y_ * static_cast<std::size_t>(local.y + h) +
      halo_stride_z_ * static_cast<std::size_t>(local.z + h);
  return offset < halo_region_.size() &&
         halo_region_[offset] == static_cast<std::uint8_t>(RegionFlag::fluid);
}

Status EBTopologyCompiler::compile(
    MPI_Comm communicator, const CartesianGeometryPlan& geometry,
    const MeshPatch& patch, const StlScanPlan& scan,
    const ImmersedSurfacePlan& surface, ImmersedFluidSide fluid_side,
    ImmersedPlanLimits limits, EBTopology& out) noexcept {
  if (!mpi_live() || communicator == MPI_COMM_NULL) {
    return {StatusCode::invalid_plan, kTopologyInput};
  }
  int rank = -1;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return {StatusCode::mpi_failure, kTopologyCollective};
  }
  int lowest = -1;
  const Int3 global = geometry.global_cells();
  std::uint64_t global_count{};
  std::uint64_t local_count{};
  Status local =
      geometry.fingerprint() == 0U || geometry.topology_revision() == 0U ||
              scan.fingerprint() == 0U || surface.fingerprint() == 0U ||
              scan.triangles().fingerprint() == 0U ||
              surface.source_triangle_fingerprint() !=
                  scan.triangles().fingerprint() ||
              !checked_product(global, global_count) ||
              global_count > (UINT64_MAX - 5U) / 6U ||
              !checked_product(patch.cells, local_count) ||
              !valid_patch(global, patch, rank, size) ||
              (fluid_side != ImmersedFluidSide::outside &&
               fluid_side != ImmersedFluidSide::inside) ||
              limits.maximum_persistent_bytes_per_rank == 0U ||
              limits.maximum_peak_bytes_per_rank == 0U ||
              limits.maximum_peak_bytes_per_rank <
                  limits.maximum_persistent_bytes_per_rank ||
              limits.maximum_local_links == 0U ||
              limits.maximum_local_links > UINT32_MAX ||
              local_count > UINT32_MAX ||
              local_count > std::numeric_limits<std::size_t>::max()
          ? Status{StatusCode::invalid_plan, kTopologyInput}
          : Status{};
  Status agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) {
    out.lowest_failing_rank_ = lowest;
    return agreed;
  }

  const PlanFingerprint contract =
      compile_contract(geometry, patch, surface, fluid_side, limits);
  PlanFingerprint root_contract = contract;
  if (MPI_Bcast(&root_contract, 1, MPI_UINT64_T, 0, communicator) !=
      MPI_SUCCESS) {
    local = {StatusCode::mpi_failure, kTopologyCollective};
  } else if (contract != root_contract) {
    local = {StatusCode::invalid_plan, kTopologyInput};
  }
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) {
    out.lowest_failing_rank_ = lowest;
    return agreed;
  }
  EBTopology candidate;
  const std::int64_t halo_x =
      static_cast<std::int64_t>(patch.cells.x) + 2 * kRegionHalo;
  const std::int64_t halo_y =
      static_cast<std::int64_t>(patch.cells.y) + 2 * kRegionHalo;
  const std::int64_t halo_z =
      static_cast<std::int64_t>(patch.cells.z) + 2 * kRegionHalo;
  Int3 halo_shape{};
  std::uint64_t halo_count{};
  std::uint64_t maximum_exchange_elements{};
  std::vector<std::uint8_t> send_lower, send_upper, receive_lower,
      receive_upper;
  if (halo_x > std::numeric_limits<std::int32_t>::max() ||
      halo_y > std::numeric_limits<std::int32_t>::max() ||
      halo_z > std::numeric_limits<std::int32_t>::max()) {
    local = {StatusCode::invalid_plan, kTopologyRegion};
  } else {
    halo_shape = {static_cast<std::int32_t>(halo_x),
                  static_cast<std::int32_t>(halo_y),
                  static_cast<std::int32_t>(halo_z)};
    for (int axis = 0; axis < 3; ++axis) {
      std::uint64_t axis_elements{};
      if (!exchange_buffer_elements(patch, global, axis, halo_shape,
                                    axis_elements)) {
        local = {StatusCode::invalid_plan, kTopologyRegion};
        break;
      }
      maximum_exchange_elements =
          std::max(maximum_exchange_elements, axis_elements);
    }
  }
  if (local && !checked_product(halo_shape, halo_count)) {
    local = {StatusCode::invalid_plan, kTopologyRegion};
  }
  std::uint64_t persistent_bytes{};
  std::uint64_t exchange_bytes{};
  std::uint64_t peak_bytes{};
  if (local &&
      (!checked_add(local_count, halo_count, persistent_bytes) ||
       !checked_multiply(maximum_exchange_elements, 4U, exchange_bytes) ||
       !checked_add(persistent_bytes, exchange_bytes, peak_bytes) ||
       persistent_bytes > limits.maximum_persistent_bytes_per_rank ||
       peak_bytes > limits.maximum_peak_bytes_per_rank ||
       halo_count > std::numeric_limits<std::size_t>::max() ||
       maximum_exchange_elements >
           std::numeric_limits<std::size_t>::max())) {
    local = {StatusCode::invalid_plan, kTopologyRegion};
  }
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) {
    out.lowest_failing_rank_ = lowest;
    return agreed;
  }
  try {
    candidate.region_.resize(static_cast<std::size_t>(local_count));
    candidate.halo_region_.assign(
        static_cast<std::size_t>(halo_count),
        static_cast<std::uint8_t>(RegionFlag::solid));
    const std::size_t exchange_elements =
        static_cast<std::size_t>(maximum_exchange_elements);
    send_lower.resize(exchange_elements);
    send_upper.resize(exchange_elements);
    receive_lower.resize(exchange_elements);
    receive_upper.resize(exchange_elements);
    candidate.halo_stride_y_ = static_cast<std::size_t>(halo_shape.x);
    candidate.halo_stride_z_ = candidate.halo_stride_y_ *
                               static_cast<std::size_t>(halo_shape.y);
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, kTopologyAllocation};
  } catch (...) {
    local = {StatusCode::invalid_plan, kTopologyAllocation};
  }
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) {
    out.lowest_failing_rank_ = lowest;
    return agreed;
  }

  local = scan.classify(
      geometry, patch,
      {candidate.region_.data(), candidate.region_.size()});
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) {
    out.lowest_failing_rank_ = lowest;
    return agreed;
  }
  const int h = static_cast<int>(kRegionHalo);
  std::uint64_t local_xor = 0U;
  std::uint64_t local_sum = 0U;
  for (std::int32_t z = 0; z < patch.cells.z; ++z) {
    for (std::int32_t y = 0; y < patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        const Int3 owned{x, y, z};
        const Int3 global_index{patch.begin.x + x, patch.begin.y + y,
                                patch.begin.z + z};
        const std::size_t owned_index = flat(patch.cells, owned);
        const std::uint8_t region =
            selected_region(candidate.region_[owned_index], fluid_side);
        candidate.region_[owned_index] = region;
        halo_at(candidate.halo_region_, halo_shape,
                {x + h, y + h, z + h}) = region;
        const std::uint64_t contribution =
            mix_cell(global_id(global_index, global), region);
        local_xor ^= contribution;
        local_sum += contribution;
      }
    }
  }

  try {
    for (int axis = 0; axis < 3 && local; ++axis) {
      local = exchange_axis(communicator, patch, global, axis, halo_shape,
                            candidate.halo_region_, send_lower, send_upper,
                            receive_lower, receive_upper);
    }
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, kTopologyAllocation};
  } catch (...) {
    local = {StatusCode::invalid_plan, kTopologyRegion};
  }
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) {
    out.lowest_failing_rank_ = lowest;
    return agreed;
  }
  std::vector<std::uint8_t>().swap(send_lower);
  std::vector<std::uint8_t>().swap(send_upper);
  std::vector<std::uint8_t>().swap(receive_lower);
  std::vector<std::uint8_t>().swap(receive_upper);

  std::uint64_t local_link_count = 0U;
  std::uint64_t local_interface_count = 0U;
  for (std::int32_t z = 0; z < patch.cells.z; ++z) {
    for (std::int32_t y = 0; y < patch.cells.y; ++y) {
      for (std::int32_t x = 0; x < patch.cells.x; ++x) {
        const Int3 fluid_local{x, y, z};
        if (candidate.region_[flat(patch.cells, fluid_local)] !=
            static_cast<std::uint8_t>(RegionFlag::fluid)) {
          continue;
        }
        bool interface_cell = false;
        const Int3 fluid{patch.begin.x + x, patch.begin.y + y,
                         patch.begin.z + z};
        for (std::size_t direction = 0U;
             direction < kDirectionDelta.size(); ++direction) {
          const Int3 delta = kDirectionDelta[direction];
          const Int3 solid{fluid.x + delta.x, fluid.y + delta.y,
                           fluid.z + delta.z};
          if (!inside(solid, global)) {
            continue;
          }
          const Int3 halo_index{x + h + delta.x, y + h + delta.y,
                                z + h + delta.z};
          if (halo_at(candidate.halo_region_, halo_shape, halo_index) ==
              static_cast<std::uint8_t>(RegionFlag::solid)) {
            interface_cell = true;
            ++local_link_count;
          }
        }
        local_interface_count += interface_cell ? 1U : 0U;
      }
    }
  }
  std::uint64_t link_bytes{};
  std::uint64_t interface_bytes{};
  std::uint64_t full_persistent_bytes{};
  std::uint64_t full_peak_bytes{};
  if (local_link_count > limits.maximum_local_links ||
      local_link_count > UINT32_MAX ||
      local_interface_count > UINT32_MAX ||
      !checked_multiply(local_link_count, sizeof(ImmersedLink), link_bytes) ||
      !checked_multiply(local_interface_count, sizeof(std::uint32_t),
                        interface_bytes) ||
      !checked_add(persistent_bytes, link_bytes, full_persistent_bytes) ||
      !checked_add(full_persistent_bytes, interface_bytes,
                   full_persistent_bytes) ||
      full_persistent_bytes > limits.maximum_persistent_bytes_per_rank ||
      local_link_count > std::numeric_limits<std::size_t>::max() ||
      local_interface_count > std::numeric_limits<std::size_t>::max()) {
    local = {StatusCode::invalid_plan, kTopologyLinks};
  } else {
    full_peak_bytes = std::max(peak_bytes, full_persistent_bytes);
    if (full_peak_bytes > limits.maximum_peak_bytes_per_rank) {
      local = {StatusCode::invalid_plan, kTopologyLinks};
    }
  }
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) {
    out.lowest_failing_rank_ = lowest;
    return agreed;
  }
  try {
    candidate.links_.reserve(static_cast<std::size_t>(local_link_count));
    candidate.interface_cells_.reserve(
        static_cast<std::size_t>(local_interface_count));
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, kTopologyAllocation};
  } catch (...) {
    local = {StatusCode::invalid_plan, kTopologyAllocation};
  }
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) {
    out.lowest_failing_rank_ = lowest;
    return agreed;
  }

  std::uint64_t global_xor{};
  std::uint64_t global_sum{};
  const int xor_reduction =
      MPI_Allreduce(&local_xor, &global_xor, 1, MPI_UINT64_T, MPI_BXOR,
                    communicator);
  const int sum_reduction =
      MPI_Allreduce(&local_sum, &global_sum, 1, MPI_UINT64_T, MPI_SUM,
                    communicator);
  if (xor_reduction != MPI_SUCCESS || sum_reduction != MPI_SUCCESS) {
    local = {StatusCode::mpi_failure, kTopologyCollective};
  }
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) {
    out.lowest_failing_rank_ = lowest;
    return agreed;
  }

  try {
    for (std::int32_t z = 0; z < patch.cells.z && local; ++z) {
      for (std::int32_t y = 0; y < patch.cells.y && local; ++y) {
        for (std::int32_t x = 0; x < patch.cells.x && local; ++x) {
          const Int3 fluid_local{x, y, z};
          const Int3 fluid{patch.begin.x + x, patch.begin.y + y,
                           patch.begin.z + z};
          if (candidate.region_[flat(patch.cells, fluid_local)] !=
              static_cast<std::uint8_t>(RegionFlag::fluid)) {
            continue;
          }
          bool interface_cell = false;
          for (std::size_t direction = 0U;
               direction < kDirectionDelta.size(); ++direction) {
            const Int3 delta = kDirectionDelta[direction];
            const Int3 solid{fluid.x + delta.x, fluid.y + delta.y,
                             fluid.z + delta.z};
            if (!inside(solid, global)) {
              continue;
            }
            const Int3 halo_index{x + h + delta.x, y + h + delta.y,
                                  z + h + delta.z};
            if (halo_at(candidate.halo_region_, halo_shape, halo_index) !=
                static_cast<std::uint8_t>(RegionFlag::solid)) {
              continue;
            }
            interface_cell = true;
            if (candidate.links_.size() >= limits.maximum_local_links) {
              local = {StatusCode::invalid_plan, kTopologyLinks};
              break;
            }
            const Real3 fluid_centre = cell_centre(geometry, fluid);
            const Real3 solid_centre = cell_centre(geometry, solid);
            SurfaceSegmentIntersection intersection;
            local = surface.first_segment_intersection(
                fluid_centre, solid_centre, intersection);
            if (!local || intersection.triangle == kInvalidSurfaceTriangle ||
                intersection.triangle >= surface.triangles().size) {
              if (local) {
                local = {StatusCode::invalid_plan, kTopologyGeometry};
              }
              break;
            }
            const Real3 outward =
                surface.triangles()
                    .data[static_cast<std::size_t>(intersection.triangle)]
                    .geometric_outward_normal;
            const Real3 normal =
                fluid_side == ImmersedFluidSide::outside
                    ? outward
                    : Real3{-outward.x, -outward.y, -outward.z};
            const double area = face_area(geometry, fluid, direction);
            ImmersedLink link;
            // Stable sparse key: no global prefix scan or full-domain pass.
            link.global_link =
                6U * global_id(fluid, global) + direction;
            link.fluid_cell = global_id(fluid, global);
            link.solid_cell = global_id(solid, global);
            link.fluid_global_index = fluid;
            link.solid_global_index = solid;
            link.fluid_local_index = fluid_local;
            link.solid_local_index = {solid.x - patch.begin.x,
                                      solid.y - patch.begin.y,
                                      solid.z - patch.begin.z};
            link.direction = static_cast<ImmersedFaceDirection>(direction);
            link.triangle = intersection.triangle;
            link.wall_point = intersection.point;
            link.solid_to_fluid_normal = normal;
            // Cartesian finite-volume control metric.  The physical STL
            // quadrature metric is compiled independently below.
            link.cartesian_control_face_area = area;
            link.surface_patch_centroid = intersection.point;
            candidate.links_.push_back(link);
          }
          if (interface_cell) {
            candidate.interface_cells_.push_back(
                static_cast<std::uint32_t>(flat(patch.cells, fluid_local)));
          }
        }
      }
    }
  } catch (const std::bad_alloc&) {
    local = {StatusCode::allocation_failure, kTopologyAllocation};
  } catch (...) {
    local = {StatusCode::invalid_plan, kTopologyGeometry};
  }
  if (local &&
      (candidate.links_.size() != local_link_count ||
       candidate.interface_cells_.size() != local_interface_count)) {
    local = {StatusCode::invalid_plan, kTopologyLinks};
  }
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) {
    out.lowest_failing_rank_ = lowest;
    return agreed;
  }

  candidate.global_cells_ = global;
  candidate.patch_ = patch;
  candidate.fluid_side_ = fluid_side;
  candidate.region_halo_width_ = kRegionHalo;
  candidate.geometry_revision_ = geometry.topology_revision();
  candidate.geometry_fingerprint_ = geometry.fingerprint();
  candidate.surface_fingerprint_ = surface.fingerprint();
  candidate.fingerprint_ = topology_fingerprint(
      geometry, surface, fluid_side, global_count, global_xor, global_sum);
  local = IbmInterfaceMetricCompiler::compile_with_resident_storage(
      communicator, geometry, patch, surface, candidate, limits,
      full_persistent_bytes, candidate.interface_metric_);
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) {
    out.lowest_failing_rank_ = lowest;
    return agreed;
  }
  const IbmInterfaceMetricResources metric_resources =
      candidate.interface_metric_.resources();
  std::uint64_t combined_persistent = 0U;
  std::uint64_t combined_peak = 0U;
  if (!checked_add(full_persistent_bytes,
                   metric_resources.persistent_bytes_per_rank,
                   combined_persistent) ||
      !checked_add(full_persistent_bytes,
                   metric_resources.peak_bytes_per_rank, combined_peak) ||
      combined_persistent > limits.maximum_persistent_bytes_per_rank ||
      combined_peak > limits.maximum_peak_bytes_per_rank) {
    local = {StatusCode::invalid_plan, kTopologyLinks};
  }
  agreed = consensus(communicator, rank, size, local, lowest);
  if (!agreed) {
    out.lowest_failing_rank_ = lowest;
    return agreed;
  }
  Hash64 sealed;
  sealed.integer(candidate.fingerprint_);
  sealed.integer(candidate.interface_metric_.physical_fingerprint());
  candidate.fingerprint_ = sealed.finish();
  candidate.lowest_failing_rank_ = -1;
  out = std::move(candidate);
  return {};
}

}  // namespace hundun::v04

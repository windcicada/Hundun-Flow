// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_mesh.hpp"

#include "mesh_focus_detail.hpp"

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

constexpr std::uint32_t kGeometryCommunicator = 301U;
constexpr std::uint32_t kGeometryCollective = 302U;
constexpr std::uint32_t kGeometryCellLimit = 303U;
constexpr std::uint32_t kGeometryMemoryLimit = 304U;
constexpr std::uint32_t kGeometryWire = 305U;
constexpr std::uint64_t kFnvOffset = 14695981039346656037ULL;
constexpr std::uint64_t kFnvPrime = 1099511628211ULL;
constexpr std::uint64_t kAxisPersistentArrays = 5U;
constexpr std::uint64_t kRootMetricPeakCopies = 10U;

static_assert(sizeof(double) == sizeof(std::uint64_t) &&
                  std::numeric_limits<double>::is_iec559,
              "geometry fingerprint requires IEEE-754 binary64");
static_assert(std::is_nothrow_move_assignable_v<CartesianGeometryPlan>,
              "geometry publication must not throw");

class Hash64 {
 public:
  void bytes(const void* data, std::size_t size) noexcept {
    const auto* input = static_cast<const unsigned char*>(data);
    for (std::size_t index = 0U; index < size; ++index) {
      value_ ^= static_cast<std::uint64_t>(input[index]);
      value_ *= kFnvPrime;
    }
  }

  template <class Integer>
  void integer(Integer value) noexcept {
    using Unsigned = std::make_unsigned_t<Integer>;
    Unsigned bits = static_cast<Unsigned>(value);
    for (std::size_t byte = 0U; byte < sizeof(bits); ++byte) {
      const auto part = static_cast<unsigned char>(
          (bits >> (byte * 8U)) & static_cast<Unsigned>(0xffU));
      bytes(&part, 1U);
    }
  }

  void real(double value) noexcept {
    std::uint64_t bits = 0U;
    std::memcpy(&bits, &value, sizeof(bits));
    integer(bits);
  }

  PlanFingerprint finish() const noexcept { return value_ == 0U ? 1U : value_; }

 private:
  std::uint64_t value_{kFnvOffset};
};

bool checked_product(Int3 cells, std::uint64_t& out) noexcept {
  if (cells.x <= 0 || cells.y <= 0 || cells.z <= 0) {
    return false;
  }
  const auto x = static_cast<std::uint64_t>(cells.x);
  const auto y = static_cast<std::uint64_t>(cells.y);
  const auto z = static_cast<std::uint64_t>(cells.z);
  if (x > std::numeric_limits<std::uint64_t>::max() / y) {
    return false;
  }
  const std::uint64_t xy = x * y;
  if (xy > std::numeric_limits<std::uint64_t>::max() / z) {
    return false;
  }
  out = xy * z;
  return true;
}

bool checked_add_u64(std::uint64_t left, std::uint64_t right,
                     std::uint64_t& out) noexcept {
  if (right > std::numeric_limits<std::uint64_t>::max() - left) {
    return false;
  }
  out = left + right;
  return true;
}

bool checked_multiply_u64(std::uint64_t left, std::uint64_t right,
                          std::uint64_t& out) noexcept {
  if (left != 0U && right > std::numeric_limits<std::uint64_t>::max() / left) {
    return false;
  }
  out = left * right;
  return true;
}

Status geometry_payload_limit(const CartesianMeshSpec& mesh,
                              GeometryBudget budget,
                              std::uint64_t& out) noexcept {
  if (mesh.limits.max_memory_bytes_per_rank == 0U ||
      budget.fixed_bytes_per_rank >
          mesh.limits.max_memory_bytes_per_rank) {
    return {StatusCode::invalid_plan, kGeometryMemoryLimit};
  }
  out = mesh.limits.max_memory_bytes_per_rank - budget.fixed_bytes_per_rank;
  return out == 0U ? Status{StatusCode::invalid_plan, kGeometryMemoryLimit}
                   : Status{};
}

Status preflight_exact_geometry(const CartesianMeshSpec& mesh,
                                std::uint64_t payload_limit) noexcept {
  if (!mesh.has_exact_cells || mesh.exact_cells.x <= 0 ||
      mesh.exact_cells.y <= 0 || mesh.exact_cells.z <= 0) {
    return {};
  }
  const std::array<std::uint64_t, 3U> counts{
      static_cast<std::uint64_t>(mesh.exact_cells.x),
      static_cast<std::uint64_t>(mesh.exact_cells.y),
      static_cast<std::uint64_t>(mesh.exact_cells.z)};
  std::uint64_t total_elements = 0U;
  for (const std::uint64_t count : counts) {
    std::uint64_t axis_elements = 0U;
    if (!checked_multiply_u64(count, kRootMetricPeakCopies, axis_elements) ||
        !checked_add_u64(axis_elements, 1U, axis_elements) ||
        !checked_add_u64(total_elements, axis_elements, total_elements)) {
      return {StatusCode::invalid_plan, kGeometryMemoryLimit};
    }
  }
  std::uint64_t bytes = 0U;
  if (!checked_multiply_u64(total_elements, sizeof(double), bytes) ||
      bytes > payload_limit) {
    return {StatusCode::invalid_plan, kGeometryMemoryLimit};
  }
  return {};
}

Status metric_peak_gate(const std::array<std::vector<double>, 3U>& faces,
                        std::uint64_t payload_limit) noexcept {
  std::uint64_t total_cells = 0U;
  std::uint64_t total_faces = 0U;
  for (const auto& axis : faces) {
    if (axis.size() < 2U ||
        !checked_add_u64(total_faces, axis.size(), total_faces) ||
        !checked_add_u64(total_cells, axis.size() - 1U, total_cells)) {
      return {StatusCode::invalid_plan, kGeometryMemoryLimit};
    }
  }
  // During finish_axis_metrics the current face vector coexists with its five
  // persistent arrays; completed axes retain five arrays, unprocessed axes one.
  std::uint64_t peak_elements = 0U;
  std::uint64_t metric_elements = 0U;
  if (!checked_multiply_u64(total_cells, kAxisPersistentArrays,
                            metric_elements) ||
      !checked_add_u64(total_faces, metric_elements, peak_elements)) {
    return {StatusCode::invalid_plan, kGeometryMemoryLimit};
  }
  std::uint64_t bytes = 0U;
  if (!checked_multiply_u64(peak_elements, sizeof(double), bytes) ||
      bytes > payload_limit) {
    return {StatusCode::invalid_plan, kGeometryMemoryLimit};
  }
  return {};
}

Status collective_status(Status local, MPI_Comm communicator) noexcept {
  std::array<std::uint64_t, 2U> wire{};
  wire[0] = static_cast<std::uint64_t>(local.code);
  wire[1] = local.detail;
  if (MPI_Bcast(wire.data(), static_cast<int>(wire.size()), MPI_UINT64_T, 0,
                communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kGeometryCollective};
  }
  if (wire[0] > static_cast<std::uint64_t>(StatusCode::io_failure) ||
      wire[1] > std::numeric_limits<std::uint32_t>::max()) {
    return {StatusCode::invalid_plan, kGeometryWire};
  }
  return {static_cast<StatusCode>(wire[0]),
          static_cast<std::uint32_t>(wire[1])};
}

Status broadcast_geometry_metadata(GeometryKind& kind, Real3& lower,
                                   Real3& upper, int rank,
                                   MPI_Comm communicator) noexcept {
  std::uint8_t wire_kind =
      rank == 0 ? static_cast<std::uint8_t>(kind) : std::uint8_t{0U};
  std::array<double, 6U> bounds{};
  if (rank == 0) {
    bounds = {lower.x, lower.y, lower.z, upper.x, upper.y, upper.z};
  }
  if (MPI_Bcast(&wire_kind, 1, MPI_UINT8_T, 0, communicator) != MPI_SUCCESS ||
      MPI_Bcast(bounds.data(), static_cast<int>(bounds.size()), MPI_DOUBLE, 0,
                communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kGeometryCollective};
  }
  if (wire_kind >
          static_cast<std::uint8_t>(GeometryKind::tensor_stretched) ||
      !std::all_of(bounds.begin(), bounds.end(),
                   [](double value) { return std::isfinite(value); }) ||
      !(bounds[0] < bounds[3]) || !(bounds[1] < bounds[4]) ||
      !(bounds[2] < bounds[5])) {
    return {StatusCode::invalid_plan, kGeometryWire};
  }
  kind = static_cast<GeometryKind>(wire_kind);
  lower = {bounds[0], bounds[1], bounds[2]};
  upper = {bounds[3], bounds[4], bounds[5]};
  return {};
}

Status broadcast_faces(std::array<std::vector<double>, 3U>& faces,
                       int rank, MPI_Comm communicator) noexcept {
  std::array<std::uint64_t, 4U> header{};
  if (rank == 0) {
    header[0] = 1U;
    for (std::size_t axis = 0U; axis < faces.size(); ++axis) {
      header[axis + 1U] = faces[axis].size();
      if (header[axis + 1U] < 2U ||
          header[axis + 1U] > static_cast<std::uint64_t>(
                             std::numeric_limits<int>::max())) {
        header[0] = 0U;
      }
    }
  }
  if (MPI_Bcast(header.data(), static_cast<int>(header.size()), MPI_UINT64_T, 0,
                communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kGeometryCollective};
  }
  if (header[0] == 0U) {
    return {StatusCode::invalid_plan, kGeometryWire};
  }
  int allocation_ok = 1;
  if (rank != 0) {
    try {
      for (std::size_t axis = 0U; axis < faces.size(); ++axis) {
        const std::uint64_t count = header[axis + 1U];
        if (count < 2U ||
            count > static_cast<std::uint64_t>(
                               std::numeric_limits<int>::max())) {
          allocation_ok = 0;
          break;
        }
        faces[axis].resize(static_cast<std::size_t>(count));
      }
    } catch (...) {
      allocation_ok = 0;
    }
  }
  int all_allocations_ok = 0;
  if (MPI_Allreduce(&allocation_ok, &all_allocations_ok, 1, MPI_INT, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kGeometryCollective};
  }
  if (all_allocations_ok == 0) {
    return {StatusCode::allocation_failure, kGeometryWire};
  }
  for (std::size_t axis = 0U; axis < faces.size(); ++axis) {
    if (MPI_Bcast(faces[axis].data(), static_cast<int>(faces[axis].size()),
                  MPI_DOUBLE, 0, communicator) != MPI_SUCCESS) {
      return {StatusCode::mpi_failure, kGeometryCollective};
    }
  }
  return {};
}

Status stage_consensus(Status local, int rank, int size,
                       MPI_Comm communicator) noexcept {
  const int local_ok = local ? 1 : 0;
  int all_ok = 0;
  if (MPI_Allreduce(&local_ok, &all_ok, 1, MPI_INT, MPI_MIN, communicator) !=
      MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kGeometryCollective};
  }
  if (all_ok != 0) {
    return {};
  }
  const int candidate = local ? size : rank;
  int lowest = size;
  if (MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, communicator) !=
      MPI_SUCCESS || lowest < 0 || lowest >= size) {
    return {StatusCode::mpi_failure, kGeometryCollective};
  }
  std::array<std::uint64_t, 2U> wire{};
  if (rank == lowest) {
    wire[0] = static_cast<std::uint64_t>(local.code);
    wire[1] = local.detail;
  }
  if (MPI_Bcast(wire.data(), static_cast<int>(wire.size()), MPI_UINT64_T,
                lowest, communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kGeometryCollective};
  }
  if (wire[0] > static_cast<std::uint64_t>(StatusCode::io_failure) ||
      wire[1] > std::numeric_limits<std::uint32_t>::max()) {
    return {StatusCode::invalid_plan, kGeometryWire};
  }
  return {static_cast<StatusCode>(wire[0]),
          static_cast<std::uint32_t>(wire[1])};
}

PlanFingerprint geometry_fingerprint(
    GeometryKind kind, const std::array<detail::AxisPayload, 3U>& axes) noexcept {
  Hash64 hash;
  hash.integer(static_cast<std::uint8_t>(kind));
  for (const detail::AxisPayload& axis : axes) {
    hash.integer(static_cast<std::uint64_t>(axis.faces.size()));
    for (const double coordinate : axis.faces) {
      hash.real(coordinate);
    }
  }
  return hash.finish();
}

}  // namespace

Span<const double> AxisMetrics::faces() const noexcept {
  return {faces_.data(), faces_.size()};
}

Span<const double> AxisMetrics::centres() const noexcept {
  return {centres_.data(), centres_.size()};
}

Span<const double> AxisMetrics::widths() const noexcept {
  return {widths_.data(), widths_.size()};
}

Span<const double> AxisMetrics::inverse_widths() const noexcept {
  return {inverse_widths_.data(), inverse_widths_.size()};
}

Span<const double> AxisMetrics::inverse_centre_distances() const noexcept {
  return {inverse_centre_distances_.data(), inverse_centre_distances_.size()};
}

const AxisMetrics& CartesianGeometryPlan::axis(
    CartesianAxis selected) const noexcept {
  const auto index = static_cast<std::uint8_t>(selected);
  return axes_[index < 3U ? index : 0U];
}

Status CartesianGeometryCompiler::compile(
    MPI_Comm communicator, const CartesianMeshSpec& mesh, GeometryBudget budget,
    CartesianGeometryPlan& geometry, MeshPatch& patch) noexcept {
  if (communicator == MPI_COMM_NULL) {
    return {StatusCode::invalid_plan, kGeometryCommunicator};
  }
  int rank = 0;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return {StatusCode::mpi_failure, kGeometryCommunicator};
  }

  std::array<std::vector<double>, 3U> faces;
  Status root_status{};
  std::uint64_t root_payload_limit = 0U;
  if (rank == 0) {
    root_status = geometry_payload_limit(mesh, budget, root_payload_limit);
    if (root_status) {
      root_status = preflight_exact_geometry(mesh, root_payload_limit);
    }
    if (root_status) {
      root_status =
          detail::generate_cartesian_faces(mesh, root_payload_limit, faces);
    }
    if (root_status) {
      const bool representable =
          std::all_of(faces.begin(), faces.end(), [](const auto& axis) {
            return axis.size() >= 2U &&
                   axis.size() - 1U <= static_cast<std::size_t>(
                                               std::numeric_limits<std::int32_t>::max());
          });
      std::uint64_t global_cells = 0U;
      if (!representable) {
        root_status = {StatusCode::invalid_plan, kGeometryCellLimit};
      } else {
        const Int3 cells{
            static_cast<std::int32_t>(faces[0].size() - 1U),
            static_cast<std::int32_t>(faces[1].size() - 1U),
            static_cast<std::int32_t>(faces[2].size() - 1U)};
        if (!checked_product(cells, global_cells)) {
          root_status = {StatusCode::invalid_plan, kGeometryCellLimit};
        }
      }
      if (root_status &&
          (mesh.limits.max_global_cells == 0U ||
           global_cells > mesh.limits.max_global_cells)) {
        root_status = {StatusCode::invalid_plan, kGeometryCellLimit};
      }
      if (root_status) {
        root_status = metric_peak_gate(faces, root_payload_limit);
      }
    }
  }
  const Status generated = collective_status(root_status, communicator);
  if (!generated) {
    return generated;
  }
  GeometryKind authoritative_kind = mesh.kind;
  Real3 authoritative_lower = mesh.lower;
  Real3 authoritative_upper = mesh.upper;
  const Status metadata_status = broadcast_geometry_metadata(
      authoritative_kind, authoritative_lower, authoritative_upper, rank,
      communicator);
  if (!metadata_status) {
    return metadata_status;
  }
  const Status broadcast = broadcast_faces(faces, rank, communicator);
  if (!broadcast) {
    return broadcast;
  }

  std::array<detail::AxisPayload, 3U> axes;
  Status local_status{};
  for (std::size_t axis = 0U; axis < axes.size(); ++axis) {
    const Status metrics =
        detail::finish_axis_metrics(std::move(faces[axis]), axes[axis]);
    if (!metrics) {
      local_status = metrics;
      break;
    }
  }
  const Status metric_consensus =
      stage_consensus(local_status, rank, size, communicator);
  if (!metric_consensus) {
    return metric_consensus;
  }
  const Int3 global_cells{
      static_cast<std::int32_t>(axes[0].widths.size()),
      static_cast<std::int32_t>(axes[1].widths.size()),
      static_cast<std::int32_t>(axes[2].widths.size())};
  MeshPatch candidate_patch;
  const Status decomposed =
      detail::make_mesh_patch(rank, size, global_cells, candidate_patch);
  const Status decomposition_consensus =
      stage_consensus(decomposed, rank, size, communicator);
  if (!decomposition_consensus) {
    return decomposition_consensus;
  }

  const std::uint64_t maximum_local =
      detail::maximum_patch_cells(global_cells, candidate_patch.process_grid);
  local_status = {};
  if (maximum_local == 0U ||
      (budget.bytes_per_owned_cell_upper_bound != 0U &&
       maximum_local >
           (std::numeric_limits<std::uint64_t>::max() -
            budget.fixed_bytes_per_rank) /
               budget.bytes_per_owned_cell_upper_bound)) {
    local_status = {StatusCode::invalid_plan, kGeometryMemoryLimit};
  }
  if (local_status) {
    const std::uint64_t bytes =
        budget.fixed_bytes_per_rank +
        maximum_local * budget.bytes_per_owned_cell_upper_bound;
    if (mesh.limits.max_memory_bytes_per_rank == 0U ||
        bytes > mesh.limits.max_memory_bytes_per_rank) {
      local_status = {StatusCode::invalid_plan, kGeometryMemoryLimit};
    }
  }
  const Status memory_consensus =
      stage_consensus(local_status, rank, size, communicator);
  if (!memory_consensus) {
    return memory_consensus;
  }

  CartesianGeometryPlan candidate;
  candidate.kind_ = authoritative_kind;
  candidate.global_cells_ = global_cells;
  candidate.lower_ = authoritative_lower;
  candidate.upper_ = authoritative_upper;
  candidate.topology_revision_ = 1U;
  candidate.fingerprint_ = geometry_fingerprint(authoritative_kind, axes);
  for (std::size_t axis = 0U; axis < axes.size(); ++axis) {
    candidate.axes_[axis].faces_ = std::move(axes[axis].faces);
    candidate.axes_[axis].centres_ = std::move(axes[axis].centres);
    candidate.axes_[axis].widths_ = std::move(axes[axis].widths);
    candidate.axes_[axis].inverse_widths_ =
        std::move(axes[axis].inverse_widths);
    candidate.axes_[axis].inverse_centre_distances_ =
        std::move(axes[axis].inverse_centre_distances);
    candidate.axes_[axis].uniform_ = axes[axis].uniform;
    candidate.axes_[axis].uniform_width_ = axes[axis].uniform_width;
    candidate.axes_[axis].uniform_inverse_width_ =
        axes[axis].uniform_inverse_width;
  }
  const Status publication_consensus =
      stage_consensus({}, rank, size, communicator);
  if (!publication_consensus) {
    return publication_consensus;
  }
  geometry = std::move(candidate);
  patch = candidate_patch;
  return {};
}

}  // namespace hundun::v04

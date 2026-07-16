// SPDX-License-Identifier: Apache-2.0

#include "hundun/runtime/structured_decomposition.hpp"

#include "hundun/runtime/error.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>

namespace hundun::runtime {
namespace {

void check_mpi(int result, const char* operation) {
  if (result != MPI_SUCCESS) {
    throw Error(std::string(operation) + " failed with MPI error " +
                std::to_string(result));
  }
}

void validate_global_extent(Int3 extent) {
  if (extent.x <= 0 || extent.y <= 0 || extent.z <= 0) {
    throw Error("structured decomposition extents must be positive");
  }

  constexpr std::int64_t limit = std::numeric_limits<std::int64_t>::max();
  std::int64_t cell_count = extent.x;
  if (cell_count > limit / extent.y) {
    throw Error("structured decomposition cell count exceeds INT64_MAX");
  }
  cell_count *= extent.y;
  if (cell_count > limit / extent.z) {
    throw Error("structured decomposition cell count exceeds INT64_MAX");
  }
}

int split_begin(int cells, int processes, int coordinate) noexcept {
  const int quotient = cells / processes;
  const int remainder = cells % processes;
  return coordinate * quotient + std::min(coordinate, remainder);
}

int split_end(int cells, int processes, int coordinate) noexcept {
  const int quotient = cells / processes;
  const int remainder = cells % processes;
  return split_begin(cells, processes, coordinate) + quotient +
         (coordinate < remainder ? 1 : 0);
}

void free_communicator_without_throwing(MPI_Comm& communicator) noexcept {
  if (communicator == MPI_COMM_NULL) {
    return;
  }

  int initialized = 0;
  if (MPI_Initialized(&initialized) != MPI_SUCCESS || initialized == 0) {
    return;
  }
  int finalized = 0;
  if (MPI_Finalized(&finalized) != MPI_SUCCESS || finalized != 0) {
    return;
  }
  (void)MPI_Comm_free(&communicator);
}

std::uint64_t checked_multiply(std::uint64_t lhs, std::uint64_t rhs) {
  constexpr std::uint64_t limit =
      std::numeric_limits<std::uint64_t>::max();
  if (rhs != 0U && lhs > limit / rhs) {
    throw Error("global cell ID multiplication overflow");
  }
  return lhs * rhs;
}

std::uint64_t checked_add(std::uint64_t lhs, std::uint64_t rhs) {
  constexpr std::uint64_t limit =
      std::numeric_limits<std::uint64_t>::max();
  if (lhs > limit - rhs) {
    throw Error("global cell ID addition overflow");
  }
  return lhs + rhs;
}

bool resolve_coordinate(int& coordinate, int processes, bool periodic) {
  if (coordinate >= 0 && coordinate < processes) {
    return true;
  }
  if (!periodic) {
    return false;
  }
  coordinate = coordinate < 0 ? processes - 1 : 0;
  return true;
}

}  // namespace

StructuredDecomposition StructuredDecomposition::create(
    MPI_Comm communicator, Int3 global_extent,
    std::array<bool, 3> periodic) {
  if (communicator == MPI_COMM_NULL) {
    throw Error("structured decomposition requires a valid communicator");
  }

  int is_intercommunicator = 0;
  check_mpi(MPI_Comm_test_inter(communicator, &is_intercommunicator),
            "MPI_Comm_test_inter");
  if (is_intercommunicator != 0) {
    throw Error("structured decomposition requires an intracommunicator");
  }

  validate_global_extent(global_extent);

  int process_count = 0;
  check_mpi(MPI_Comm_size(communicator, &process_count), "MPI_Comm_size");
  std::array<int, 3> dimensions{0, 0, 0};
  check_mpi(MPI_Dims_create(process_count, 3, dimensions.data()),
            "MPI_Dims_create");
  if (dimensions[0] > global_extent.x ||
      dimensions[1] > global_extent.y ||
      dimensions[2] > global_extent.z) {
    throw Error("structured decomposition process grid exceeds global extent");
  }

  const std::array<int, 3> periods{periodic[0] ? 1 : 0,
                                   periodic[1] ? 1 : 0,
                                   periodic[2] ? 1 : 0};
  MPI_Comm cartesian = MPI_COMM_NULL;
  check_mpi(MPI_Cart_create(communicator, 3, dimensions.data(), periods.data(),
                            0, &cartesian),
            "MPI_Cart_create");
  if (cartesian == MPI_COMM_NULL) {
    throw Error("MPI_Cart_create returned MPI_COMM_NULL");
  }

  try {
    check_mpi(MPI_Comm_set_errhandler(cartesian, MPI_ERRORS_RETURN),
              "MPI_Comm_set_errhandler");
    int rank = 0;
    check_mpi(MPI_Comm_rank(cartesian, &rank), "MPI_Comm_rank");
    std::array<int, 3> coordinates{};
    check_mpi(MPI_Cart_coords(cartesian, rank, 3, coordinates.data()),
              "MPI_Cart_coords");

    const Int3 process_grid{dimensions[0], dimensions[1], dimensions[2]};
    const Int3 process_coordinates{coordinates[0], coordinates[1],
                                   coordinates[2]};
    const Box3 owned_box{
        Int3{split_begin(global_extent.x, process_grid.x,
                         process_coordinates.x),
             split_begin(global_extent.y, process_grid.y,
                         process_coordinates.y),
             split_begin(global_extent.z, process_grid.z,
                         process_coordinates.z)},
        Int3{split_end(global_extent.x, process_grid.x,
                       process_coordinates.x),
             split_end(global_extent.y, process_grid.y,
                       process_coordinates.y),
             split_end(global_extent.z, process_grid.z,
                       process_coordinates.z)}};
    return StructuredDecomposition(cartesian, global_extent, process_grid,
                                   process_coordinates, owned_box, periodic);
  } catch (...) {
    free_communicator_without_throwing(cartesian);
    throw;
  }
}

StructuredDecomposition::StructuredDecomposition(
    MPI_Comm communicator, Int3 global_extent, Int3 process_grid,
    Int3 process_coordinates, Box3 owned_box,
    std::array<bool, 3> periodic) noexcept
    : communicator_(communicator),
      global_extent_(global_extent),
      process_grid_(process_grid),
      process_coordinates_(process_coordinates),
      owned_box_(owned_box),
      periodic_(periodic) {}

StructuredDecomposition::~StructuredDecomposition() {
  free_communicator_without_throwing(communicator_);
}

StructuredDecomposition::StructuredDecomposition(
    StructuredDecomposition&& other) noexcept
    : communicator_(std::exchange(other.communicator_, MPI_COMM_NULL)),
      global_extent_(other.global_extent_),
      process_grid_(other.process_grid_),
      process_coordinates_(other.process_coordinates_),
      owned_box_(other.owned_box_),
      periodic_(other.periodic_) {}

MPI_Comm StructuredDecomposition::comm() const noexcept {
  return communicator_;
}

Int3 StructuredDecomposition::global_extent() const noexcept {
  return global_extent_;
}

Int3 StructuredDecomposition::process_grid() const noexcept {
  return process_grid_;
}

Int3 StructuredDecomposition::process_coordinates() const noexcept {
  return process_coordinates_;
}

Box3 StructuredDecomposition::owned_box() const noexcept {
  return owned_box_;
}

Int3 StructuredDecomposition::local_extent() const noexcept {
  return Int3{owned_box_.end.x - owned_box_.begin.x,
              owned_box_.end.y - owned_box_.begin.y,
              owned_box_.end.z - owned_box_.begin.z};
}

int StructuredDecomposition::neighbor_rank(Int3 offset) const {
  const bool in_range = offset.x >= -1 && offset.x <= 1 &&
                        offset.y >= -1 && offset.y <= 1 &&
                        offset.z >= -1 && offset.z <= 1;
  const bool is_zero = offset.x == 0 && offset.y == 0 && offset.z == 0;
  if (!in_range || is_zero) {
    throw Error("neighbor offset must identify one of the 26 neighbors");
  }

  std::array<int, 3> coordinates{
      process_coordinates_.x + offset.x,
      process_coordinates_.y + offset.y,
      process_coordinates_.z + offset.z};
  const std::array<int, 3> dimensions{process_grid_.x, process_grid_.y,
                                      process_grid_.z};
  for (std::size_t axis = 0; axis < coordinates.size(); ++axis) {
    if (!resolve_coordinate(coordinates[axis], dimensions[axis],
                            periodic_[axis])) {
      return MPI_PROC_NULL;
    }
  }

  int rank = MPI_PROC_NULL;
  check_mpi(MPI_Cart_rank(communicator_, coordinates.data(), &rank),
            "MPI_Cart_rank");
  return rank;
}

Int3 StructuredDecomposition::global_cell(Int3 local_cell) const {
  const Int3 extent = local_extent();
  if (local_cell.x < 0 || local_cell.x >= extent.x || local_cell.y < 0 ||
      local_cell.y >= extent.y || local_cell.z < 0 ||
      local_cell.z >= extent.z) {
    throw Error("local cell is outside the owned box");
  }
  return Int3{owned_box_.begin.x + local_cell.x,
              owned_box_.begin.y + local_cell.y,
              owned_box_.begin.z + local_cell.z};
}

std::uint64_t StructuredDecomposition::global_cell_id(
    Int3 local_cell) const {
  const Int3 cell = global_cell(local_cell);
  const std::uint64_t k = static_cast<std::uint64_t>(cell.z);
  const std::uint64_t j = static_cast<std::uint64_t>(cell.y);
  const std::uint64_t i = static_cast<std::uint64_t>(cell.x);
  const std::uint64_t global_y =
      static_cast<std::uint64_t>(global_extent_.y);
  const std::uint64_t global_x =
      static_cast<std::uint64_t>(global_extent_.x);
  const std::uint64_t plane_row = checked_add(checked_multiply(k, global_y), j);
  return checked_add(checked_multiply(plane_row, global_x), i);
}

}  // namespace hundun::runtime

// SPDX-License-Identifier: Apache-2.0

#include "hundun/runtime/structured_decomposition.hpp"

#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"
#include "mpi_error.hpp"
#include "process_grid.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <exception>
#include <limits>
#include <string>
#include <utility>

namespace hundun::runtime {
namespace {

void require_collective_input_agreement(const MpiContext& context,
                                        Int3 global_extent,
                                        std::array<bool, 3> periodic,
                                        const DecompositionOptions& options) {
  const Int3 grid = options.process_grid.value_or(Int3{});
  const std::array<int, 10> local{
      global_extent.x,
      global_extent.y,
      global_extent.z,
      periodic[0] ? 1 : 0,
      periodic[1] ? 1 : 0,
      periodic[2] ? 1 : 0,
      options.process_grid.has_value() ? 1 : 0,
      grid.x,
      grid.y,
      grid.z};
  std::array<int, 10> minimum{};
  std::array<int, 10> maximum{};
  detail::check_mpi(
      MPI_Allreduce(local.data(), minimum.data(),
                    static_cast<int>(local.size()), MPI_INT, MPI_MIN,
                    context.comm()),
      "MPI_Allreduce");
  detail::check_mpi(
      MPI_Allreduce(local.data(), maximum.data(),
                    static_cast<int>(local.size()), MPI_INT, MPI_MAX,
                    context.comm()),
      "MPI_Allreduce");
  const bool inputs_agree = minimum == maximum;
  const CollectiveStatus status = collective_status(
      context, inputs_agree,
      inputs_agree ? "" : "structured decomposition inputs differ across ranks");
  if (!status.ok) {
    throw Error(status.message);
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
    const MpiContext& context, Int3 global_extent,
    std::array<bool, 3> periodic, DecompositionOptions options) {
  detail::require_mpi_active("create structured decomposition");
  const MPI_Comm communicator = context.comm();
  if (communicator == MPI_COMM_NULL) {
    throw Error("structured decomposition requires a valid MPI context");
  }

  require_collective_input_agreement(context, global_extent, periodic,
                                     options);

  const int process_count = context.size();
  Int3 selected_grid{};
  bool local_ok = true;
  std::string local_message;
  try {
    if (options.process_grid.has_value()) {
      detail::validate_explicit_process_grid(
          *options.process_grid, process_count, global_extent);
      selected_grid = *options.process_grid;
    } else {
      selected_grid =
          detail::select_process_grid(process_count, global_extent, periodic);
    }
  } catch (const std::exception& error) {
    local_ok = false;
    local_message = error.what();
  } catch (...) {
    local_ok = false;
    local_message = "unknown process-grid selection failure";
  }
  const CollectiveStatus selection_status =
      collective_status(context, local_ok, local_message);
  if (!selection_status.ok) {
    throw Error(selection_status.message);
  }

  const std::array<int, 3> dimensions{selected_grid.x, selected_grid.y,
                                      selected_grid.z};

  const std::array<int, 3> periods{periodic[0] ? 1 : 0,
                                   periodic[1] ? 1 : 0,
                                   periodic[2] ? 1 : 0};
  MPI_Comm cartesian = MPI_COMM_NULL;
  detail::check_mpi(MPI_Cart_create(communicator, 3, dimensions.data(),
                                    periods.data(), 0, &cartesian),
                    "MPI_Cart_create");
  if (cartesian == MPI_COMM_NULL) {
    throw Error("MPI_Cart_create returned MPI_COMM_NULL");
  }

  try {
    detail::check_mpi(
        MPI_Comm_set_errhandler(cartesian, MPI_ERRORS_RETURN),
        "MPI_Comm_set_errhandler");
    int rank = 0;
    detail::check_mpi(MPI_Comm_rank(cartesian, &rank), "MPI_Comm_rank");
    std::array<int, 3> coordinates{};
    detail::check_mpi(
        MPI_Cart_coords(cartesian, rank, 3, coordinates.data()),
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
    detail::free_communicator_without_throwing(cartesian);
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
  detail::free_communicator_without_throwing(communicator_);
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

std::array<bool, 3> StructuredDecomposition::periodic() const noexcept {
  return periodic_;
}

int StructuredDecomposition::neighbor_rank(Int3 offset) const {
  const bool in_range = offset.x >= -1 && offset.x <= 1 &&
                        offset.y >= -1 && offset.y <= 1 &&
                        offset.z >= -1 && offset.z <= 1;
  const bool is_zero = offset.x == 0 && offset.y == 0 && offset.z == 0;
  if (!in_range || is_zero) {
    throw Error("neighbor offset must identify one of the 26 neighbors");
  }

  detail::require_mpi_active("query a structured-decomposition neighbor");

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
  detail::check_mpi(MPI_Cart_rank(communicator_, coordinates.data(), &rank),
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

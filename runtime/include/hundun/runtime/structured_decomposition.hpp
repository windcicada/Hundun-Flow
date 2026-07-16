// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/runtime/types.hpp"

#include <mpi.h>

#include <array>
#include <cstdint>

namespace hundun::runtime {

class StructuredDecomposition final {
 public:
  static StructuredDecomposition create(
      MPI_Comm, Int3 global_extent, std::array<bool, 3> periodic);
  ~StructuredDecomposition();
  StructuredDecomposition(StructuredDecomposition&&) noexcept;
  StructuredDecomposition(const StructuredDecomposition&) = delete;
  MPI_Comm comm() const noexcept;
  Int3 global_extent() const noexcept;
  Int3 process_grid() const noexcept;
  Int3 process_coordinates() const noexcept;
  Box3 owned_box() const noexcept;
  Int3 local_extent() const noexcept;
  int neighbor_rank(Int3 offset) const;
  std::uint64_t global_cell_id(Int3 local_cell) const;
  Int3 global_cell(Int3 local_cell) const;

 private:
  StructuredDecomposition(MPI_Comm communicator, Int3 global_extent,
                          Int3 process_grid, Int3 process_coordinates,
                          Box3 owned_box, std::array<bool, 3> periodic) noexcept;

  MPI_Comm communicator_{MPI_COMM_NULL};
  Int3 global_extent_{};
  Int3 process_grid_{};
  Int3 process_coordinates_{};
  Box3 owned_box_{};
  std::array<bool, 3> periodic_{};
};

}  // namespace hundun::runtime

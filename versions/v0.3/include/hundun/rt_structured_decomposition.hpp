// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_types.hpp"

#include <mpi.h>

#include <array>
#include <cstdint>
#include <optional>

namespace hundun::runtime {

struct DecompositionOptions {
  std::optional<Int3> process_grid;
};

class StructuredDecomposition final {
 public:
  // Collective over the context communicator. All members must call create in
  // a compatible order; the returned decomposition owns its Cartesian
  // communicator.
  static StructuredDecomposition create(
      const MpiContext&, Int3 global_extent, std::array<bool, 3> periodic,
      DecompositionOptions options = {});
  // While MPI is active, destruction collectively frees the Cartesian
  // communicator, so member ranks must destroy decompositions in a compatible
  // order and normally before MpiEnvironment. After MPI_Finalize, destruction
  // only clears the local handle.
  ~StructuredDecomposition();
  StructuredDecomposition(StructuredDecomposition&&) noexcept;
  StructuredDecomposition(const StructuredDecomposition&) = delete;
  // Borrowed handle: never free it. Pass it to MPI only while MPI is active
  // and this decomposition has not been moved from.
  MPI_Comm comm() const noexcept;
  Int3 global_extent() const noexcept;
  Int3 process_grid() const noexcept;
  Int3 process_coordinates() const noexcept;
  Box3 owned_box() const noexcept;
  Int3 local_extent() const noexcept;
  std::array<bool, 3> periodic() const noexcept;
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

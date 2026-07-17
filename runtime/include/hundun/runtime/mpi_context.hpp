// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <mpi.h>

namespace hundun::runtime {

class MpiContext final {
 public:
  // Collective over source. All source members must call duplicate in a
  // compatible order. The returned context owns a duplicate; source is
  // borrowed and is never modified.
  static MpiContext duplicate(MPI_Comm source);

  // While MPI is active, destruction collectively frees the owned duplicate;
  // communicator members must destroy contexts in a compatible order.
  // After MPI_Finalize, destruction only clears the local handle.
  ~MpiContext() noexcept;
  MpiContext(MpiContext&& other) noexcept;
  // While MPI is active, replacing a non-null destination collectively frees
  // its old duplicate and therefore requires compatible member-rank order.
  // After MPI_Finalize, replacement only clears the local destination handle.
  MpiContext& operator=(MpiContext&& other) noexcept;
  MpiContext(const MpiContext&) = delete;
  MpiContext& operator=(const MpiContext&) = delete;

  // Borrowed handle: never free it. Pass it to MPI only while MPI is active
  // and this context has not been moved from.
  MPI_Comm comm() const noexcept { return communicator_; }
  int rank() const noexcept { return rank_; }
  int size() const noexcept { return size_; }
  int thread_level() const noexcept { return thread_level_; }
  void barrier() const;

 private:
  MpiContext(MPI_Comm communicator, int rank, int size,
             int thread_level) noexcept;

  MPI_Comm communicator_{MPI_COMM_NULL};
  int rank_{0};
  int size_{0};
  int thread_level_{MPI_THREAD_SINGLE};
};

// Normally every MpiContext is destroyed before the MpiEnvironment that keeps
// MPI active. Intentionally reversed destruction remains local after finalize.

}  // namespace hundun::runtime

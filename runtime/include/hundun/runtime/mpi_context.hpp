// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <mpi.h>

namespace hundun::runtime {

class MpiContext final {
 public:
  // Collective over source. The returned context owns a duplicate; source is
  // borrowed and is never modified.
  static MpiContext duplicate(MPI_Comm source);

  ~MpiContext() noexcept;
  MpiContext(MpiContext&& other) noexcept;
  MpiContext& operator=(MpiContext&& other) noexcept;
  MpiContext(const MpiContext&) = delete;
  MpiContext& operator=(const MpiContext&) = delete;

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

// In normal lexical use, MpiContext objects are destroyed before the
// MpiEnvironment that keeps MPI active. Destruction is nevertheless guarded
// so an intentionally reversed order never calls MPI_Comm_free after finalize.

}  // namespace hundun::runtime

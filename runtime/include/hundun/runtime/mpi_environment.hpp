// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <mpi.h>

namespace hundun::runtime {

class MpiEnvironment final {
 public:
  MpiEnvironment(int& argc, char**& argv);
  ~MpiEnvironment() noexcept;

  MpiEnvironment(const MpiEnvironment&) = delete;
  MpiEnvironment& operator=(const MpiEnvironment&) = delete;

  MPI_Comm comm() const noexcept { return MPI_COMM_WORLD; }
  int rank() const noexcept { return rank_; }
  int size() const noexcept { return size_; }
  void barrier() const;

 private:
  bool owns_mpi_{false};
  int rank_{0};
  int size_{1};
};

}  // namespace hundun::runtime

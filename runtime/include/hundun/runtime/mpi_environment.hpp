// SPDX-License-Identifier: Apache-2.0
#pragma once

namespace hundun::runtime {

class MpiEnvironment final {
 public:
  MpiEnvironment(int& argc, char**& argv);
  ~MpiEnvironment() noexcept;

  MpiEnvironment(const MpiEnvironment&) = delete;
  MpiEnvironment& operator=(const MpiEnvironment&) = delete;

 private:
  bool owns_mpi_{false};
};

}  // namespace hundun::runtime

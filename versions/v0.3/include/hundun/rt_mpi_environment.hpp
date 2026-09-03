// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
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

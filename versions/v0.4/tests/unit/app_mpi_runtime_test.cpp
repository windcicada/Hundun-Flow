// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_mpi_runtime.hpp"

#include <mpi.h>

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <optional>
#include <string>

int main() {
  constexpr const char* key =
      "OMPI_MCA_btl_vader_single_copy_mechanism";
  const char* current = std::getenv(key);
  const std::optional<std::string> saved =
      current == nullptr ? std::nullopt
                         : std::optional<std::string>{current};

  bool passed = ::unsetenv(key) == 0 &&
                hundun::v04::prepare_mpi_runtime_environment();
#if defined(OPEN_MPI) && OMPI_MAJOR_VERSION == 2 && OMPI_MINOR_VERSION == 1
  const char* automatic = std::getenv(key);
  passed = passed && automatic != nullptr &&
           std::strcmp(automatic, "none") == 0;
#else
  passed = passed && std::getenv(key) == nullptr;
#endif

  passed = passed && ::setenv(key, "cma", 1) == 0 &&
           hundun::v04::prepare_mpi_runtime_environment();
  const char* explicit_value = std::getenv(key);
  passed = passed && explicit_value != nullptr &&
           std::strcmp(explicit_value, "cma") == 0;

  if (saved.has_value())
    passed = passed && ::setenv(key, saved->c_str(), 1) == 0;
  else
    passed = passed && ::unsetenv(key) == 0;

  std::cout << (passed ? "PASS" : "FAIL")
            << " app_mpi_runtime_environment\n";
  return passed ? 0 : 1;
}

// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_mpi_runtime.hpp"

#include <mpi.h>

#include <cstdlib>

namespace hundun::v04 {

bool prepare_mpi_runtime_environment() noexcept {
#if defined(OPEN_MPI) && OMPI_MAJOR_VERSION == 2 && OMPI_MINOR_VERSION == 1
  constexpr const char* key =
      "OMPI_MCA_btl_vader_single_copy_mechanism";
  // Open MPI 2.1's vader RGET/CMA receive path retains 32 KiB PML free-list
  // slabs under repeated persistent halo traffic. The copy path avoids that
  // unbounded high-water mark. Do not override an explicit operator choice.
  if (std::getenv(key) != nullptr) return true;
  return ::setenv(key, "none", 0) == 0;
#else
  return true;
#endif
}

}  // namespace hundun::v04

// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/rt_mpi_operation_error.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

int main() {
  return hundun::test::run([] {
    hundun::runtime::check_mpi_result(MPI_SUCCESS, "success");
    bool typed = false;
    try {
      hundun::runtime::check_mpi_result(MPI_ERR_OTHER, "direct-test");
    } catch (const hundun::runtime::MpiOperationError &) {
      typed = true;
    }
    HUNDUN_CHECK(typed);
  });
}

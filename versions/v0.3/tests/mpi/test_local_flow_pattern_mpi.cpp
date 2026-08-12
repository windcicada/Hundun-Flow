// SPDX-License-Identifier: Apache-2.0

#include "hundun/ib_local_flow_pattern.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

namespace {

using hundun::immersed::LocalCoefficientRow;
using hundun::immersed::LocalFlowPatternTransform;
using hundun::runtime::Error;
using hundun::runtime::MpiEnvironment;

std::uint64_t bits(double value) {
  std::uint64_t result = 0;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

void run_parallel_evidence() {
  int rank = 0;
  int size = 0;
  HUNDUN_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &size) == MPI_SUCCESS);
  const LocalFlowPatternTransform transform;
  const LocalCoefficientRow row{{2.0, 3.0, 5.0, 0.0, 0.0, 0.0}, 7.0, 11.0};
  std::vector<std::uint64_t> links{10U, 20U, 30U};
  if ((rank % 2) != 0) {
    std::reverse(links.begin(), links.end());
  }
  const auto plan = transform.plan_row(42U, links, row);
  const double value =
      transform.evaluate_wall_replacement(plan, row, {1.0, 2.0, 4.0});
  std::uint64_t local[2]{plan.fingerprint, bits(value)};
  std::uint64_t minimum[2]{};
  std::uint64_t maximum[2]{};
  HUNDUN_CHECK(MPI_Allreduce(local, minimum, 2, MPI_UINT64_T, MPI_MIN,
                             MPI_COMM_WORLD) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Allreduce(local, maximum, 2, MPI_UINT64_T, MPI_MAX,
                             MPI_COMM_WORLD) == MPI_SUCCESS);
  HUNDUN_CHECK(minimum[0] == maximum[0]);
  HUNDUN_CHECK(minimum[1] == maximum[1]);

  int local_failed = 0;
  int local_classification = 0;
  if (size > 1 && rank == 1) {
    try {
      static_cast<void>(transform.plan_row(42U, {10U, 10U}, row));
    } catch (const Error &error) {
      HUNDUN_CHECK(std::string(error.what()).find("duplicate") !=
                   std::string::npos);
      local_failed = 1;
      local_classification = 17;
    }
  } else if (size == 1) {
    try {
      static_cast<void>(transform.plan_row(42U, {10U, 10U}, row));
    } catch (const Error &) {
      local_failed = 1;
      local_classification = 17;
    }
  }
  int failure_count = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_failed, &failure_count, 1, MPI_INT, MPI_SUM,
                             MPI_COMM_WORLD) == MPI_SUCCESS);
  HUNDUN_CHECK(failure_count == 1);
  int classification = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_classification, &classification, 1, MPI_INT,
                             MPI_MAX, MPI_COMM_WORLD) == MPI_SUCCESS);
  HUNDUN_CHECK(classification == 17);
  const int local_failure_rank =
      local_failed != 0 ? rank : std::numeric_limits<int>::max();
  int lowest_failure_rank = -1;
  HUNDUN_CHECK(MPI_Allreduce(&local_failure_rank, &lowest_failure_rank, 1,
                             MPI_INT, MPI_MIN, MPI_COMM_WORLD) == MPI_SUCCESS);
  HUNDUN_CHECK(lowest_failure_rank == (size == 1 ? 0 : 1));
}

} // namespace

int main(int argc, char **argv) {
  MpiEnvironment environment(argc, argv);
  return hundun::test::run(run_parallel_evidence);
}

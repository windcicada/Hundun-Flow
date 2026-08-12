// SPDX-License-Identifier: Apache-2.0

#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <cstdlib>

int main(int argc, char** argv) {
  int provided = MPI_THREAD_SINGLE;
  const int init_result =
      MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  if (init_result != MPI_SUCCESS) {
    return EXIT_FAILURE;
  }
  if (provided < MPI_THREAD_FUNNELED) {
    const int finalize_result = MPI_Finalize();
    if (finalize_result != MPI_SUCCESS) {
      return EXIT_FAILURE;
    }
    return EXIT_FAILURE;
  }

  const int borrowed_result = hundun::test::run([&] {
    {
      hundun::runtime::MpiEnvironment borrowed(argc, argv);
      auto context =
          hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
      HUNDUN_CHECK(context.rank() == 0);
      HUNDUN_CHECK(context.size() == 1);
      context.barrier();
    }

    int finalized = 0;
    HUNDUN_CHECK(MPI_Finalized(&finalized) == MPI_SUCCESS);
    HUNDUN_CHECK(finalized == 0);
  });

  const int finalize_result = MPI_Finalize();
  return borrowed_result == EXIT_SUCCESS && finalize_result == MPI_SUCCESS
             ? EXIT_SUCCESS
             : EXIT_FAILURE;
}

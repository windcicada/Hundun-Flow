// SPDX-License-Identifier: Apache-2.0

#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/mpi_environment.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <cstdlib>

int main(int argc, char** argv) {
  int provided = MPI_THREAD_SINGLE;
  if (MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided) !=
          MPI_SUCCESS ||
      provided < MPI_THREAD_FUNNELED) {
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

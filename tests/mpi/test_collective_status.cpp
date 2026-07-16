// SPDX-License-Identifier: Apache-2.0

#include "runtime/include/hundun/runtime/collective_status.hpp"
#include "runtime/include/hundun/runtime/error.hpp"
#include "runtime/include/hundun/runtime/mpi_environment.hpp"
#include "tests/support/test_main.hpp"

#include <cstdlib>
#include <optional>
#include <string>
#include <string_view>

namespace {

using hundun::runtime::CollectiveStatus;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::collective_status;
using hundun::runtime::require_expected_ranks;

constexpr std::string_view kIntercommunicatorError =
    "collective_status requires an intracommunicator";

void check_status(const CollectiveStatus& status, bool ok, int failing_rank,
                  const std::string& message) {
  HUNDUN_CHECK(status.ok == ok);
  HUNDUN_CHECK(status.failing_rank == failing_rank);
  HUNDUN_CHECK(status.message == message);
}

void expect_rank_error(MPI_Comm comm, std::optional<int> expected_ranks,
                       const std::string& expected_message) {
  bool threw = false;
  try {
    require_expected_ranks(comm, expected_ranks);
  } catch (const hundun::runtime::Error& error) {
    threw = true;
    HUNDUN_CHECK(std::string(error.what()) == expected_message);
  }
  HUNDUN_CHECK(threw);
}

template <class Function>
void expect_intercommunicator_error(Function&& function) {
  bool threw = false;
  try {
    function();
  } catch (const hundun::runtime::Error& error) {
    threw = true;
    HUNDUN_CHECK(std::string_view(error.what()) == kIntercommunicatorError);
  }
  HUNDUN_CHECK(threw);
}

void test_intercommunicator_rejection(MpiEnvironment& mpi) {
  MPI_Comm local_comm = MPI_COMM_NULL;
  HUNDUN_CHECK(MPI_Comm_split(mpi.comm(), mpi.rank(), 0, &local_comm) ==
               MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Comm_set_errhandler(local_comm, MPI_ERRORS_RETURN) ==
               MPI_SUCCESS);

  int local_size = 0;
  HUNDUN_CHECK(MPI_Comm_size(local_comm, &local_size) == MPI_SUCCESS);
  HUNDUN_CHECK(local_size == 1);

  MPI_Comm intercomm = MPI_COMM_NULL;
  const int remote_leader = 1 - mpi.rank();
  HUNDUN_CHECK(MPI_Intercomm_create(local_comm, 0, mpi.comm(), remote_leader,
                                    17, &intercomm) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Comm_set_errhandler(intercomm, MPI_ERRORS_RETURN) ==
               MPI_SUCCESS);

  int is_inter = 0;
  HUNDUN_CHECK(MPI_Comm_test_inter(intercomm, &is_inter) == MPI_SUCCESS);
  HUNDUN_CHECK(is_inter != 0);

  expect_intercommunicator_error([&] {
    static_cast<void>(collective_status(
        intercomm, mpi.rank() != 0,
        mpi.rank() == 0 ? "intercommunicator failure" : ""));
  });
  expect_intercommunicator_error(
      [&] { require_expected_ranks(intercomm, local_size); });

  HUNDUN_CHECK(MPI_Comm_free(&intercomm) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Comm_free(&local_comm) == MPI_SUCCESS);
}

void test_collectives(MpiEnvironment& mpi, int& argc, char**& argv) {
  HUNDUN_CHECK(mpi.size() == 2);

  int initialized = 0;
  HUNDUN_CHECK(MPI_Initialized(&initialized) == MPI_SUCCESS);
  HUNDUN_CHECK(initialized != 0);

  int provided = MPI_THREAD_SINGLE;
  HUNDUN_CHECK(MPI_Query_thread(&provided) == MPI_SUCCESS);
  HUNDUN_CHECK(provided >= MPI_THREAD_FUNNELED);

  {
    MpiEnvironment borrowed(argc, argv);
    HUNDUN_CHECK(borrowed.rank() == mpi.rank());
    HUNDUN_CHECK(borrowed.size() == mpi.size());
    borrowed.barrier();
  }
  int finalized = 0;
  HUNDUN_CHECK(MPI_Finalized(&finalized) == MPI_SUCCESS);
  HUNDUN_CHECK(finalized == 0);

  check_status(collective_status(mpi.comm(), true, "ignored"), true, -1, "");

  const bool rank_one_ok = mpi.rank() != 1;
  check_status(collective_status(mpi.comm(), rank_one_ok,
                                 rank_one_ok ? "" : "rank-one failure"),
               false, 1, "rank-one failure");

  check_status(collective_status(
                   mpi.comm(), false,
                   mpi.rank() == 0 ? "rank-zero failure" : "later failure"),
               false, 0, "rank-zero failure");

  require_expected_ranks(mpi.comm(), std::nullopt);
  require_expected_ranks(mpi.comm(), mpi.size());
  require_expected_ranks(
      mpi.comm(), mpi.rank() == 0 ? std::optional<int>{mpi.size()}
                                  : std::optional<int>{});

  const std::string one_too_many =
      "expected MPI rank count " + std::to_string(mpi.size() + 1) +
      ", got " + std::to_string(mpi.size());
  expect_rank_error(mpi.comm(), mpi.size() + 1, one_too_many);
  expect_rank_error(mpi.comm(),
                    mpi.rank() == 0
                        ? std::optional<int>{}
                        : std::optional<int>{mpi.size() + 1},
                    one_too_many);

  const std::string rank_zero_message =
      "expected MPI rank count " + std::to_string(mpi.size() + 2) +
      ", got " + std::to_string(mpi.size());
  expect_rank_error(mpi.comm(),
                    mpi.rank() == 0
                        ? std::optional<int>{mpi.size() + 2}
                        : std::optional<int>{mpi.size() + 1},
                    rank_zero_message);

  test_intercommunicator_rejection(mpi);

  mpi.barrier();
}

}  // namespace

int main(int argc, char** argv) {
  int collective_result = EXIT_FAILURE;
  {
    MpiEnvironment mpi(argc, argv);
    collective_result = hundun::test::run(
        [&] { test_collectives(mpi, argc, argv); });
  }

  const int lifecycle_result = hundun::test::run([&] {
    int finalized = 0;
    HUNDUN_CHECK(MPI_Finalized(&finalized) == MPI_SUCCESS);
    HUNDUN_CHECK(finalized != 0);

    bool threw = false;
    try {
      MpiEnvironment rejected(argc, argv);
    } catch (const hundun::runtime::Error&) {
      threw = true;
    }
    HUNDUN_CHECK(threw);
  });

  return collective_result == EXIT_SUCCESS ? lifecycle_result
                                            : collective_result;
}

// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/rt_collective_status.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/test_main.hpp"

#include <cstdlib>
#include <optional>
#include <string>
#include <string_view>

namespace {

using hundun::runtime::CollectiveStatus;
using hundun::runtime::MpiContext;
using hundun::runtime::MpiEnvironment;
using hundun::runtime::collective_status;
using hundun::runtime::require_expected_ranks;

void check_status(const CollectiveStatus& status, bool ok, int failing_rank,
                  const std::string& message) {
  HUNDUN_CHECK(status.ok == ok);
  HUNDUN_CHECK(status.failing_rank == failing_rank);
  HUNDUN_CHECK(status.message == message);
}

void expect_rank_error(const MpiContext& context,
                       std::optional<int> expected_ranks,
                       const std::string& expected_message) {
  bool threw = false;
  try {
    require_expected_ranks(context, expected_ranks);
  } catch (const hundun::runtime::Error& error) {
    threw = true;
    HUNDUN_CHECK(std::string(error.what()) == expected_message);
  }
  HUNDUN_CHECK(threw);
}

template <class Function>
void expect_error(Function&& function, const std::string& expected_text) {
  bool threw = false;
  try {
    function();
  } catch (const hundun::runtime::Error& error) {
    threw = true;
    HUNDUN_CHECK(std::string(error.what()).find(expected_text) !=
                 std::string::npos);
  }
  HUNDUN_CHECK(threw);
}

void test_collectives(MpiContext& mpi, int& argc, char**& argv) {
  HUNDUN_CHECK(mpi.size() == 2);

  int initialized = 0;
  HUNDUN_CHECK(MPI_Initialized(&initialized) == MPI_SUCCESS);
  HUNDUN_CHECK(initialized != 0);

  int provided = MPI_THREAD_SINGLE;
  HUNDUN_CHECK(MPI_Query_thread(&provided) == MPI_SUCCESS);
  HUNDUN_CHECK(provided >= MPI_THREAD_FUNNELED);

  {
    MpiEnvironment borrowed(argc, argv);
    auto borrowed_context = MpiContext::duplicate(MPI_COMM_WORLD);
    HUNDUN_CHECK(borrowed_context.rank() == mpi.rank());
    HUNDUN_CHECK(borrowed_context.size() == mpi.size());
    borrowed_context.barrier();
  }
  int finalized = 0;
  HUNDUN_CHECK(MPI_Finalized(&finalized) == MPI_SUCCESS);
  HUNDUN_CHECK(finalized == 0);

  check_status(collective_status(mpi, true, "ignored"), true, -1, "");

  const bool rank_one_ok = mpi.rank() != 1;
  check_status(collective_status(mpi, rank_one_ok,
                                 rank_one_ok ? "" : "rank-one failure"),
               false, 1, "rank-one failure");

  check_status(collective_status(
                   mpi, false,
                   mpi.rank() == 0 ? "rank-zero failure" : "later failure"),
               false, 0, "rank-zero failure");

  require_expected_ranks(mpi, std::nullopt);
  require_expected_ranks(mpi, mpi.size());
  require_expected_ranks(
      mpi, mpi.rank() == 0 ? std::optional<int>{mpi.size()}
                           : std::optional<int>{});

  const std::string one_too_many =
      "expected MPI rank count " + std::to_string(mpi.size() + 1) +
      ", got " + std::to_string(mpi.size());
  expect_rank_error(mpi, mpi.size() + 1, one_too_many);
  expect_rank_error(mpi,
                    mpi.rank() == 0
                        ? std::optional<int>{}
                        : std::optional<int>{mpi.size() + 1},
                    one_too_many);

  const std::string rank_zero_message =
      "expected MPI rank count " + std::to_string(mpi.size() + 2) +
      ", got " + std::to_string(mpi.size());
  expect_rank_error(mpi,
                    mpi.rank() == 0
                        ? std::optional<int>{mpi.size() + 2}
                        : std::optional<int>{mpi.size() + 1},
                    rank_zero_message);

  mpi.barrier();
}

}  // namespace

int main(int argc, char** argv) {
  std::optional<MpiContext> mpi;
  int collective_result = EXIT_FAILURE;
  {
    MpiEnvironment environment(argc, argv);
    mpi.emplace(MpiContext::duplicate(MPI_COMM_WORLD));
    collective_result = hundun::test::run(
        [&] { test_collectives(*mpi, argc, argv); });
  }

  const int lifecycle_result = hundun::test::run([&] {
    int finalized = 0;
    HUNDUN_CHECK(MPI_Finalized(&finalized) == MPI_SUCCESS);
    HUNDUN_CHECK(finalized != 0);

    expect_error(
        [&] { static_cast<void>(collective_status(*mpi, true, "")); },
        "after MPI_Finalize");
    expect_error([&] { require_expected_ranks(*mpi, mpi->size()); },
                 "after MPI_Finalize");
    mpi.reset();

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

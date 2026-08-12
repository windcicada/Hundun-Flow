// SPDX-License-Identifier: Apache-2.0

#include "hundun/rt_collective_status.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "src/flow_reacting_transaction_detail.hpp"
#include "tests/support/test_main.hpp"

#include <cstdlib>
#include <optional>

namespace {

void test_collective_rollback(const hundun::runtime::MpiContext &mpi) {
  using namespace hundun::flow::detail;
  const double rank = static_cast<double>(mpi.rank());
  const ReactingLayerState history{{1.0 + rank, 2.0 + rank},
                                   {10.0 + rank}, 101000.0 + rank};
  const ReactingLayerState committed{{3.0 + rank, 4.0 + rank},
                                     {20.0 + rank}, 102000.0 + rank};
  ReactingAttemptState state(1U, 2U, history, committed);
  const auto original_history = state.history();
  const auto original_committed = state.committed();

  state.begin_attempt();
  ReactingSourceTransaction transaction(state);
  const ReactingSourceIdentity chemistry{"rank-chemistry",
                                          ReactingSourceKind::chemistry};
  transaction.add_species(chemistry, 0U, 0U, -0.25);
  transaction.add_species(chemistry, 0U, 1U, 0.25);
  const bool local_ok = mpi.rank() != 1;
  const auto status = hundun::runtime::collective_status(
      mpi, local_ok, local_ok ? "" : "injected chemistry failure");
  HUNDUN_CHECK(!transaction.commit(status));
  HUNDUN_CHECK(state.history() == original_history);
  HUNDUN_CHECK(state.committed() == original_committed);
  HUNDUN_CHECK(!state.attempt_active());
  mpi.barrier();
}

} // namespace

int main(int argc, char **argv) {
  std::optional<hundun::runtime::MpiContext> mpi;
  int result = EXIT_FAILURE;
  {
    hundun::runtime::MpiEnvironment environment(argc, argv);
    mpi.emplace(hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD));
    result = hundun::test::run([&] {
      HUNDUN_CHECK(mpi->size() == 2);
      test_collective_rollback(*mpi);
    });
    mpi.reset();
  }
  return result;
}

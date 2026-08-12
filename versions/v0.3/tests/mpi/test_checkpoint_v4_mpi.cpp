// SPDX-License-Identifier: Apache-2.0

#include "flow_checkpoint_v4_detail.hpp"

#include "hundun/rt_collective_status.hpp"
#include "hundun/rt_mpi_context.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/test_main.hpp"

int main(int argc, char **argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  return hundun::test::run([&] {
    auto mpi = hundun::runtime::MpiContext::duplicate(MPI_COMM_WORLD);
    const double rank = static_cast<double>(mpi.rank());
    hundun::flow::detail::ReactingCheckpointV4Data source{
        77U,
        {"A", "B"},
        {0.7 + rank, 0.3, 0.6 + rank, 0.4},
        {0.65 + rank, 0.35, 0.55 + rank, 0.45},
        {10.0 + rank, 11.0},
        {12.0 + rank, 13.0},
        100000.0,
        101325.0,
        9U,
        0.09,
        0.01,
        2U,
        "cantera-3.2.0",
        std::string(64U, 'a'),
        "gas"};
    const auto sections =
        hundun::flow::detail::encode_reacting_checkpoint_v4(source);
    auto restored = source;
    restored.step = 0U;
    std::string message;
    const bool local_ok =
        hundun::flow::detail::restore_reacting_checkpoint_v4(
            sections, 77U, std::string(64U, 'a'), restored, message) &&
        restored.step == 9U &&
        restored.committed_rho_y_kg_per_m3 ==
            source.committed_rho_y_kg_per_m3;
    const auto status = hundun::runtime::collective_status(
        mpi, local_ok, local_ok ? "" : message);
    HUNDUN_CHECK(status.ok);
    restored.committed_rho_y_kg_per_m3[0] += 0.01;
    source.committed_rho_y_kg_per_m3[0] += 0.01;
    HUNDUN_CHECK(restored.committed_rho_y_kg_per_m3 ==
                 source.committed_rho_y_kg_per_m3);
  });
}

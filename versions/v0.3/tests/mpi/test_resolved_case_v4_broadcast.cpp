// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "src/app_resolved_case_v4_broadcast_detail.hpp"
#include "src/cfg_resolved_case_v4_loader_detail.hpp"

#include "hundun/cfg_resolved_case_v4_loader.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <functional>
#include <stdexcept>
#include <string>

namespace {

using hundun::config::ResolvedReactingCaseV4;
using hundun::runtime::Error;

ResolvedReactingCaseV4 root_case(int rank, int size) {
  ResolvedReactingCaseV4 value{};
  if (rank == 0) {
    const std::string json =
        R"({"schema_version":4,"case":{"name":"stage4-broadcast"},"simulation":{"type":"reacting_flow","density_model":"ideal_gas"},"resources":{},"mesh":{"cells":[4,4,4],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"mapping":"uniform_box"},"time":{"mode":"fixed","steps":2,"initial_dt_s":0.001,"min_dt_s":0.000125,"max_dt_s":0.001,"cfl_target":0.5,"diffusion_number_target":0.25,"growth_factor":1.25,"retry_factor":0.5,"max_retries":8},"physics":{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.001,"inlet_consistency_rtol":1e-12},"scalars":[],"boundaries":[{"patch":"x_min","type":"no_slip_wall","reacting":{"species":"non_catalytic_impermeable","thermal":{"mode":"adiabatic"}}},{"patch":"x_max","type":"no_slip_wall","reacting":{"species":"non_catalytic_impermeable","thermal":{"mode":"isothermal","temperature_k":450.0}}},{"patch":"y_min","type":"periodic"},{"patch":"y_max","type":"periodic"},{"patch":"z_min","type":"periodic"},{"patch":"z_max","type":"periodic"}],"restart":{"read":false,"write_directory":"checkpoints","write_interval":2},"diagnostics":{"directory":"diagnostics","write_interval":1,"write_mesh":true},"performance":{"enabled":false,"directory":"performance","warmup_steps":5,"measured_steps":20,"repetitions":5},"immersed_boundary":{"model":"none"},"les":{"model":"none"},"reacting":{"chemistry":{"backend":"cantera","mechanism":{"file":"mechanisms/h2.yaml","sha256":"0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef","phase":"gas"},"relative_tolerance":1e-8,"absolute_tolerance":1e-14,"maximum_internal_steps":5000},"thermodynamics":{"initial_p0_pa":101325.0,"initial_temperature_k":300.0,"initial_mass_fractions":[{"species":"H2","value":0.25},{"species":"O2","value":0.75}]},"transport":{"model":"mixture_averaged"},"pressure_constraint":{"mode":"open_fixed_p0"}}})";
    value = hundun::config::detail::parse_resolved_reacting_case_v4_json(
        json, "stage4-broadcast-case.json");
    value.common_flow.resources.expected_ranks = size;
    value.common_flow.resources.process_grid =
        hundun::runtime::Int3{size, 1, 1};
  }
  return value;
}

void expect_collective_error(const std::function<void()> &action,
                             const std::string &token) {
  bool rejected = false;
  try {
    action();
  } catch (const Error &error) {
    const std::string message = error.what();
    HUNDUN_CHECK(message.find(token) != std::string::npos);
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
}

void run_success(int rank, int size) {
  auto root = root_case(rank, size);
  const auto value = hundun::config::broadcast_resolved_reacting_case_v4(
      MPI_COMM_WORLD, 0, rank == 0 ? &root : nullptr);
  const std::string canonical =
      hundun::config::to_resolved_reacting_json_v4(value);
  std::uint64_t length = rank == 0 ? canonical.size() : 0U;
  HUNDUN_CHECK(MPI_Bcast(&length, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD) ==
               MPI_SUCCESS);
  std::string expected = rank == 0 ? canonical : std::string(length, '\0');
  HUNDUN_CHECK(MPI_Bcast(expected.data(), static_cast<int>(length), MPI_CHAR, 0,
                         MPI_COMM_WORLD) == MPI_SUCCESS);
  HUNDUN_CHECK(canonical == expected);
  HUNDUN_CHECK(value.species_names[0] == "H2");
  HUNDUN_CHECK(value.species_names[1] == "O2");
}

void run_failures(int rank, int size) {
  auto root = root_case(rank, size);
  expect_collective_error(
      [&] {
        static_cast<void>(hundun::config::broadcast_resolved_reacting_case_v4(
            MPI_COMM_WORLD, 0, nullptr));
      },
      "root requires");
  if (size > 1) {
    expect_collective_error(
        [&] {
          static_cast<void>(hundun::config::broadcast_resolved_reacting_case_v4(
              MPI_COMM_WORLD, 0, &root));
        },
        "non-root");
  }

#if defined(HUNDUN_APPLICATION_ENABLE_TEST_ACCESS)
  using Fault = hundun::config::detail::ResolvedReactingCaseV4BroadcastFault;
  expect_collective_error(
      [&] {
        static_cast<void>(
            hundun::config::detail::
                broadcast_resolved_reacting_case_v4_with_fault(
                    MPI_COMM_WORLD, 0, rank == 0 ? &root : nullptr,
                    Fault::reverse_species));
      },
      "species order");
  expect_collective_error(
      [&] {
        static_cast<void>(
            hundun::config::detail::
                broadcast_resolved_reacting_case_v4_with_fault(
                    MPI_COMM_WORLD, 0, rank == 0 ? &root : nullptr,
                    Fault::rank_local_error));
      },
      "injected rank-local");
#else
#error "v4 broadcast failure tests require test access"
#endif
}

} // namespace

int main(int argc, char **argv) {
  return hundun::test::run([&] {
    hundun::runtime::MpiEnvironment environment(argc, argv);
    int rank = 0;
    int size = 0;
    HUNDUN_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank) == MPI_SUCCESS);
    HUNDUN_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &size) == MPI_SUCCESS);
    HUNDUN_CHECK(argc == 1 || argc == 2);
    const std::string mode = argc == 1 ? "success" : argv[1];
    if (mode == "success") {
      run_success(rank, size);
    } else if (mode == "failures") {
      run_failures(rank, size);
    } else {
      throw std::runtime_error("unknown mode");
    }
  });
}

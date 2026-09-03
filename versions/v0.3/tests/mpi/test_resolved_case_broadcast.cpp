// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "src/app_case_config_broadcast_detail.hpp"

#include "hundun/cfg_resolved_case_loader.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>

#include <unistd.h>

namespace {

using hundun::config::CaseConfig;
using hundun::config::FlowCaseConfig;
using hundun::config::ResolvedCase;
using hundun::runtime::Error;

std::string v1_json() {
  return R"({"schema_version":1,"case":{"name":"broadcast_v1"},"resources":{"expected_ranks":1,"process_grid":[1,1,1]},"mesh":{"cells":[8,4,2],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"periodic":[true,true,true]},"time":{"dt_s":0.01,"steps":2},"transport":{"velocity_m_per_s":[1.0,0.0,0.0],"diffusivity_m2_per_s":0.0},"initial_condition":{"type":"sine_x"},"restart":{"read":false,"write_directory":"Restart"},"output":{"directory":"output","write_interval":1,"restart_interval":1}})";
}

std::string v2_json() {
  return R"({"performance":{"repetitions":5,"measured_steps":20,"warmup_steps":5,"directory":"performance","enabled":false},"diagnostics":{"write_mesh":true,"write_interval":1,"directory":"diagnostics"},"restart":{"write_interval":10,"write_directory":"checkpoints","read_directory":"previous/step00000010","read":true},"boundaries":[{"type":"periodic","patch":"z_max"},{"type":"symmetry","patch":"y_min"},{"type":"pressure_outlet","patch":"x_max","pressure_perturbation_pa":0.0},{"type":"periodic","patch":"z_min"},{"scalar_values":[{"value":0.75,"name":"zeta"},{"value":0.25,"name":"alpha"}],"density_kg_per_m3":1.176829268292683,"enthalpy_J_per_kg":301500.0,"temperature_K":300.0,"thermal_authority":"temperature","velocity_m_per_s":[1.0,0.0,0.0],"type":"velocity_inlet","patch":"x_min"},{"type":"symmetry","patch":"y_max"}],"scalars":[{"diffusivity_m2_per_s":0.0,"name":"zeta"},{"diffusivity_m2_per_s":0.001,"name":"alpha"}],"physics":{"thermodynamic_pressure_pa":101325.0,"gas_constant_J_per_kg_K":287.0,"cp_J_per_kg_K":1005.0,"inlet_consistency_rtol":1e-12,"dynamic_viscosity_pa_s":0.001,"rho_ref_kg_per_m3":1.0},"time":{"max_retries":8,"retry_factor":0.5,"growth_factor":1.25,"diffusion_number_target":0.25,"cfl_target":0.5,"max_dt_s":0.001,"min_dt_s":0.000125,"initial_dt_s":0.001,"steps":10,"mode":"fixed"},"mesh":{"warp_amplitude":[0.02,-0.015,0.01],"mapping":"analytic_warped_box","length_m":[1.0,1.0,1.0],"origin_m":[0.0,0.0,0.0],"cells":[8,8,4]},"simulation":{"density_model":"ideal_gas","type":"variable_density_flow"},"case":{"name":"broadcast_v2"},"schema_version":2})";
}

class RootCaseFile final {
 public:
  RootCaseFile(int rank, int size, const std::string& contents) : rank_(rank) {
    path_ = std::filesystem::temp_directory_path() /
            ("hundun-resolved-broadcast-" + std::to_string(size) + "-" +
             std::to_string(static_cast<long long>(::getpid())) + ".json");
    if (rank_ == 0) {
      std::ofstream stream(path_, std::ios::binary);
      HUNDUN_CHECK(static_cast<bool>(stream));
      stream << contents;
      HUNDUN_CHECK(static_cast<bool>(stream));
    }
    HUNDUN_CHECK(MPI_Barrier(MPI_COMM_WORLD) == MPI_SUCCESS);
  }

  ~RootCaseFile() {
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank_ == 0) {
      std::error_code error;
      std::filesystem::remove(path_, error);
    }
  }

  const std::filesystem::path& path() const noexcept { return path_; }

 private:
  int rank_;
  std::filesystem::path path_;
};

void broadcast_expected_string(MPI_Comm comm, int rank, std::string& value) {
  int length = rank == 0 ? static_cast<int>(value.size()) : 0;
  HUNDUN_CHECK(MPI_Bcast(&length, 1, MPI_INT, 0, comm) == MPI_SUCCESS);
  HUNDUN_CHECK(length >= 0);
  value.resize(static_cast<std::size_t>(length));
  HUNDUN_CHECK(MPI_Bcast(value.data(), length, MPI_BYTE, 0, comm) ==
               MPI_SUCCESS);
}

ResolvedCase make_v1_root_case(int rank, int size) {
  RootCaseFile file(rank, size, v1_json());
  if (rank != 0) {
    return ResolvedCase(CaseConfig{});
  }
  CaseConfig config = hundun::config::load_case_config(file.path());
  config.expected_ranks = size;
  config.process_grid = hundun::runtime::Int3{size, 1, 1};
  return ResolvedCase(std::move(config));
}

ResolvedCase make_v2_root_case(int rank, int size) {
  RootCaseFile file(rank, size, v2_json());
  if (rank != 0) {
    return ResolvedCase(CaseConfig{});
  }
  ResolvedCase resolved = hundun::config::load_resolved_case(file.path());
  auto& config = std::get<FlowCaseConfig>(resolved);
  config.resources.expected_ranks = size;
  config.resources.process_grid = hundun::runtime::Int3{size, 1, 1};
  return resolved;
}

void check_broadcast_case(MPI_Comm comm, int rank, ResolvedCase root_case) {
  std::string expected;
  if (rank == 0) {
    expected = hundun::config::to_resolved_json(root_case);
  }
  broadcast_expected_string(comm, rank, expected);
  const ResolvedCase received = hundun::application::broadcast_resolved_case(
      comm, 0, rank == 0 ? &root_case : nullptr);
  HUNDUN_CHECK(hundun::config::to_resolved_json(received) == expected);
}

template <class Operation>
void expect_collective_error(MPI_Comm comm, Operation&& operation,
                             const std::string& stable_message = {}) {
  bool caught = false;
  std::string message;
  try {
    std::forward<Operation>(operation)();
  } catch (const Error& error) {
    caught = true;
    message = error.what();
  }
  int local_caught = caught ? 1 : 0;
  int every_caught = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local_caught, &every_caught, 1, MPI_INT, MPI_MIN,
                             comm) == MPI_SUCCESS);
  HUNDUN_CHECK(every_caught == 1);

  int rank = 0;
  HUNDUN_CHECK(MPI_Comm_rank(comm, &rank) == MPI_SUCCESS);
  std::string expected = rank == 0 ? message : std::string{};
  broadcast_expected_string(comm, rank, expected);
  HUNDUN_CHECK(!expected.empty());
  HUNDUN_CHECK(message == expected);
  if (!stable_message.empty()) {
    HUNDUN_CHECK(message == stable_message);
  }
}

void run_full(int rank, int size) {
  check_broadcast_case(MPI_COMM_WORLD, rank, make_v1_root_case(rank, size));
  check_broadcast_case(MPI_COMM_WORLD, rank, make_v2_root_case(rank, size));
}

void run_errors(int rank, int size) {
  ResolvedCase valid = make_v2_root_case(rank, size);

  expect_collective_error(MPI_COMM_WORLD, [&] {
    static_cast<void>(hundun::application::broadcast_resolved_case(
        MPI_COMM_WORLD, 0, nullptr));
  });
  expect_collective_error(MPI_COMM_WORLD, [&] {
    static_cast<void>(hundun::application::broadcast_resolved_case(
        MPI_COMM_WORLD, 0, &valid));
  });
  expect_collective_error(MPI_COMM_WORLD, [&] {
    static_cast<void>(hundun::application::broadcast_resolved_case(
        MPI_COMM_WORLD, rank % 2, rank == 0 ? &valid : nullptr));
  });
  expect_collective_error(MPI_COMM_WORLD, [&] {
    static_cast<void>(hundun::application::broadcast_resolved_case(
        MPI_COMM_WORLD, size, rank == 0 ? &valid : nullptr));
  });

  ResolvedCase mismatch = valid;
  if (rank == 0) {
    auto& flow = std::get<FlowCaseConfig>(mismatch);
    flow.resources.expected_ranks = 1;
    flow.resources.process_grid = hundun::runtime::Int3{1, 1, 1};
  }
  expect_collective_error(MPI_COMM_WORLD, [&] {
    static_cast<void>(hundun::application::broadcast_resolved_case(
        MPI_COMM_WORLD, 0, rank == 0 ? &mismatch : nullptr));
  });

  ResolvedCase grid_only_mismatch = valid;
  if (rank == 0) {
    auto& flow = std::get<FlowCaseConfig>(grid_only_mismatch);
    flow.resources.expected_ranks.reset();
    flow.resources.process_grid = hundun::runtime::Int3{1, 1, 1};
  }
  expect_collective_error(
      MPI_COMM_WORLD,
      [&] {
        static_cast<void>(hundun::application::broadcast_resolved_case(
            MPI_COMM_WORLD, 0,
            rank == 0 ? &grid_only_mismatch : nullptr));
      },
      "resolved-case process-grid product does not equal communicator size");

  bool null_rejected = false;
  try {
    static_cast<void>(hundun::application::broadcast_resolved_case(
        MPI_COMM_NULL, 0, nullptr));
  } catch (const Error&) {
    null_rejected = true;
  }
  HUNDUN_CHECK(null_rejected);

  HUNDUN_CHECK(size % 2 == 0);
  const int half = size / 2;
  const int color = rank < half ? 0 : 1;
  MPI_Comm local = MPI_COMM_NULL;
  HUNDUN_CHECK(MPI_Comm_split(MPI_COMM_WORLD, color, rank, &local) ==
               MPI_SUCCESS);
  MPI_Comm inter = MPI_COMM_NULL;
  const int remote_leader = color == 0 ? half : 0;
  HUNDUN_CHECK(MPI_Intercomm_create(local, 0, MPI_COMM_WORLD, remote_leader,
                                    77, &inter) == MPI_SUCCESS);
  bool inter_rejected = false;
  try {
    static_cast<void>(hundun::application::broadcast_resolved_case(
        inter, 0, nullptr));
  } catch (const Error&) {
    inter_rejected = true;
  }
  HUNDUN_CHECK(inter_rejected);
  HUNDUN_CHECK(MPI_Comm_free(&inter) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Comm_free(&local) == MPI_SUCCESS);
}

void run_v1_expected_rank_mismatch(int rank, int size) {
  ResolvedCase mismatch = make_v1_root_case(rank, size);
  const int expected_ranks = size == 1 ? 2 : 1;
  if (rank == 0) {
    auto& stage1 = std::get<CaseConfig>(mismatch);
    stage1.expected_ranks = expected_ranks;
    stage1.process_grid.reset();
  }
  expect_collective_error(
      MPI_COMM_WORLD,
      [&] {
        static_cast<void>(hundun::application::broadcast_resolved_case(
            MPI_COMM_WORLD, 0, rank == 0 ? &mismatch : nullptr));
      },
      "expected MPI rank count " + std::to_string(expected_ranks) + ", got " +
          std::to_string(size));
}

int run_active(int argc, char** argv) {
  hundun::runtime::MpiEnvironment environment(argc, argv);
  int rank = 0;
  int size = 0;
  HUNDUN_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &size) == MPI_SUCCESS);
  HUNDUN_CHECK(argc == 2);
  const std::string mode(argv[1]);
  if (mode == "full") {
    run_full(rank, size);
  } else if (mode == "errors") {
    run_errors(rank, size);
  } else if (mode == "v1_mismatch") {
    run_v1_expected_rank_mismatch(rank, size);
  } else {
    throw std::runtime_error("unknown test mode");
  }
  return EXIT_SUCCESS;
}

int run_finalized(int argc, char** argv) {
  HUNDUN_CHECK(MPI_Init(&argc, &argv) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Finalize() == MPI_SUCCESS);
  bool rejected = false;
  try {
    static_cast<void>(hundun::application::broadcast_resolved_case(
        MPI_COMM_WORLD, 0, nullptr));
  } catch (const Error&) {
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
  return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char** argv) {
  return hundun::test::run([&] {
    HUNDUN_CHECK(argc == 2);
    if (std::string(argv[1]) == "finalized") {
      HUNDUN_CHECK(run_finalized(argc, argv) == EXIT_SUCCESS);
    } else {
      HUNDUN_CHECK(run_active(argc, argv) == EXIT_SUCCESS);
    }
  });
}

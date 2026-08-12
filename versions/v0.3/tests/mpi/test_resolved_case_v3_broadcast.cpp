// SPDX-License-Identifier: Apache-2.0

#include "src/app_resolved_case_v3_broadcast_detail.hpp"

#include "hundun/cfg_resolved_case_v3_loader.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>

namespace {

using hundun::config::ImmersedFlowCaseConfig;
using hundun::config::ResolvedCaseV3;
using hundun::runtime::Error;

std::string v1_json() {
  return R"({"schema_version":1,"case":{"name":"stage3-broadcast-v1"},"resources":{},"mesh":{"cells":[4,4,4],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"periodic":[true,true,true]},"time":{"dt_s":0.01,"steps":1},"transport":{"velocity_m_per_s":[1.0,0.0,0.0],"diffusivity_m2_per_s":0.0},"initial_condition":{"type":"sine_x"},"restart":{"read":false,"write_directory":"Restart"},"output":{"directory":"output","write_interval":1,"restart_interval":1}})";
}

std::string v3_json(const std::string &ibm, const std::string &les) {
  return std::string(
             R"({"schema_version":3,"case":{"name":"stage3-broadcast"},"simulation":{"type":"variable_density_flow","density_model":"constant"},"resources":{},"mesh":{"cells":[4,4,4],"origin_m":[0.0,0.0,0.0],"length_m":[1.0,1.0,1.0],"mapping":"uniform_box"},"time":{"mode":"fixed","steps":2,"initial_dt_s":0.001,"min_dt_s":0.000125,"max_dt_s":0.001,"cfl_target":0.5,"diffusion_number_target":0.25,"growth_factor":1.25,"retry_factor":0.5,"max_retries":8},"physics":{"rho_ref_kg_per_m3":1.0,"dynamic_viscosity_pa_s":0.001,"inlet_consistency_rtol":1e-12},"scalars":[],"boundaries":[{"patch":"x_min","type":"periodic"},{"patch":"x_max","type":"periodic"},{"patch":"y_min","type":"periodic"},{"patch":"y_max","type":"periodic"},{"patch":"z_min","type":"periodic"},{"patch":"z_max","type":"periodic"}],"restart":{"read":false,"write_directory":"checkpoints","write_interval":2},"diagnostics":{"directory":"diagnostics","write_interval":1,"write_mesh":true},"performance":{"enabled":false,"directory":"performance","warmup_steps":5,"measured_steps":20,"repetitions":5},)") +
         ibm + "," + les + "}";
}

std::string ibm_none() { return R"("immersed_boundary":{"model":"none"})"; }
std::string ibm_lfp() {
  return R"("immersed_boundary":{"model":"local_flow_pattern_ghost_cell","geometry":{"format":"stl","file":"geometry/body.stl","length_scale_to_m":1.0,"fluid_side":"inside"},"wall":{"velocity_m_per_s":[0.0,0.0,0.0],"enthalpy":"zero_normal_diffusive_flux","scalars":"zero_normal_diffusive_flux"}})";
}
std::string les_none() { return R"("les":{"model":"none"})"; }
std::string les_wale() {
  return R"("les":{"model":"wale","wale":{"coefficient":0.5},"turbulent_prandtl":0.9,"turbulent_schmidt":0.7})";
}

std::string v2_json() {
  std::string result = v3_json(ibm_none(), les_wale());
  const auto stage3 = result.find(",\"immersed_boundary\"");
  HUNDUN_CHECK(stage3 != std::string::npos);
  result.erase(stage3);
  result += '}';
  const auto version = result.find("\"schema_version\":3");
  HUNDUN_CHECK(version != std::string::npos);
  result.replace(version, std::string("\"schema_version\":3").size(),
                 "\"schema_version\":2");
  return result;
}

class RootCaseFile final {
public:
  RootCaseFile(int rank, const std::string &contents) : rank_(rank) {
    const auto stamp =
        std::chrono::steady_clock::now().time_since_epoch().count();
    path_ = std::filesystem::temp_directory_path() /
            ("hundun-v3-broadcast-" + std::to_string(stamp) + ".json");
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

  const std::filesystem::path &path() const noexcept { return path_; }

private:
  int rank_{};
  std::filesystem::path path_;
};

void broadcast_expected(int rank, std::string &text) {
  std::uint64_t length = rank == 0 ? text.size() : 0U;
  HUNDUN_CHECK(MPI_Bcast(&length, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD) ==
               MPI_SUCCESS);
  text.resize(static_cast<std::size_t>(length));
  HUNDUN_CHECK(MPI_Bcast(text.data(), static_cast<int>(length), MPI_BYTE, 0,
                         MPI_COMM_WORLD) == MPI_SUCCESS);
}

template <class Operation>
void expect_collective_error(Operation &&operation,
                             const std::string &token = {}) {
  bool caught = false;
  std::string message;
  try {
    std::forward<Operation>(operation)();
  } catch (const Error &error) {
    caught = true;
    message = error.what();
  }
  int local = caught ? 1 : 0;
  int all = 0;
  HUNDUN_CHECK(MPI_Allreduce(&local, &all, 1, MPI_INT, MPI_MIN,
                             MPI_COMM_WORLD) == MPI_SUCCESS);
  HUNDUN_CHECK(all == 1);
  int rank = 0;
  HUNDUN_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank) == MPI_SUCCESS);
  std::string expected = rank == 0 ? message : std::string{};
  broadcast_expected(rank, expected);
  HUNDUN_CHECK(message == expected);
  if (!token.empty()) {
    HUNDUN_CHECK(message.find(token) != std::string::npos);
  }
}

ResolvedCaseV3 load_root(int rank, const std::string &json) {
  RootCaseFile file(rank, json);
  if (rank == 0) {
    return hundun::config::load_resolved_case_v3(file.path());
  }
  return ResolvedCaseV3{};
}

void check_branch(int rank, int size, const std::string &json) {
  ResolvedCaseV3 root = load_root(rank, json);
  std::string expected;
  if (rank == 0) {
    auto &stage3 = std::get<ImmersedFlowCaseConfig>(root);
    stage3.common_flow.resources.expected_ranks = size;
    stage3.common_flow.resources.process_grid = {size, 1, 1};
    expected = hundun::config::to_resolved_json_v3(root);
  }
  broadcast_expected(rank, expected);
  const auto received = hundun::config::broadcast_resolved_case_v3(
      MPI_COMM_WORLD, 0, rank == 0 ? &root : nullptr);
  HUNDUN_CHECK(hundun::config::to_resolved_json_v3(received) == expected);
}

void check_legacy(int rank, int size, const std::string &json) {
  ResolvedCaseV3 root = load_root(rank, json);
  std::string expected;
  if (rank == 0) {
    if (auto *stage1 = std::get_if<hundun::config::CaseConfig>(&root)) {
      stage1->expected_ranks = size;
      stage1->process_grid = hundun::runtime::Int3{size, 1, 1};
    } else {
      auto &flow = std::get<hundun::config::FlowCaseConfig>(root);
      flow.resources.expected_ranks = size;
      flow.resources.process_grid = hundun::runtime::Int3{size, 1, 1};
    }
    expected = hundun::config::to_resolved_json_v3(root);
  }
  broadcast_expected(rank, expected);
  const auto received = hundun::config::broadcast_resolved_case_v3(
      MPI_COMM_WORLD, 0, rank == 0 ? &root : nullptr);
  HUNDUN_CHECK(hundun::config::to_resolved_json_v3(received) == expected);
}

void run_success(int rank, int size) {
  check_legacy(rank, size, v1_json());
  check_legacy(rank, size, v2_json());
  check_branch(rank, size, v3_json(ibm_none(), les_wale()));
  check_branch(rank, size, v3_json(ibm_lfp(), les_none()));
  check_branch(rank, size, v3_json(ibm_lfp(), les_wale()));
}

void run_failures(int rank, int size) {
  ResolvedCaseV3 root = load_root(rank, v3_json(ibm_lfp(), les_wale()));
  if (rank == 0) {
    auto &stage3 = std::get<ImmersedFlowCaseConfig>(root);
    stage3.common_flow.resources.expected_ranks = size;
    stage3.common_flow.resources.process_grid = {size, 1, 1};
  }

  expect_collective_error(
      [&] {
        static_cast<void>(hundun::config::broadcast_resolved_case_v3(
            MPI_COMM_WORLD, 0, nullptr));
      },
      "root requires");
  expect_collective_error(
      [&] {
        static_cast<void>(hundun::config::broadcast_resolved_case_v3(
            MPI_COMM_WORLD, 0, &root));
      },
      "non-root");
  expect_collective_error(
      [&] {
        static_cast<void>(hundun::config::broadcast_resolved_case_v3(
            MPI_COMM_WORLD, rank % 2, rank == 0 ? &root : nullptr));
      },
      "root differs");

  ResolvedCaseV3 mismatch = root;
  if (rank == 0) {
    auto &flow = std::get<ImmersedFlowCaseConfig>(mismatch).common_flow;
    flow.resources.expected_ranks = size + 1;
    flow.resources.process_grid.reset();
  }
  expect_collective_error(
      [&] {
        static_cast<void>(hundun::config::broadcast_resolved_case_v3(
            MPI_COMM_WORLD, 0, rank == 0 ? &mismatch : nullptr));
      },
      "expected_ranks");

  ResolvedCaseV3 invalid_root = root;
  if (rank == 0) {
    std::get<ImmersedFlowCaseConfig>(invalid_root).les.wale->coefficient =
        std::numeric_limits<double>::infinity();
  }
  expect_collective_error(
      [&] {
        static_cast<void>(hundun::config::broadcast_resolved_case_v3(
            MPI_COMM_WORLD, 0, rank == 0 ? &invalid_root : nullptr));
      },
      "/les/wale/coefficient");

#if defined(HUNDUN_APPLICATION_ENABLE_TEST_ACCESS)
  using Fault = hundun::config::detail::ResolvedCaseV3BroadcastFault;
  const auto expect_fault = [&](Fault fault, const std::string &token) {
    expect_collective_error(
        [&] {
          static_cast<void>(
              hundun::config::detail::broadcast_resolved_case_v3_with_fault(
                  MPI_COMM_WORLD, 0, rank == 0 ? &root : nullptr, fault));
        },
        token);
  };
  expect_fault(Fault::version_tag, "protocol version");
  expect_fault(Fault::variant_tag, "does not match");
  expect_fault(Fault::truncate_payload, "invalid JSON");
  expect_fault(Fault::mpi_operation, "MPI operation failure");
#else
#error "v3 broadcast failure tests require the test-access contract"
#endif

  bool null_rejected = false;
  try {
    static_cast<void>(
        hundun::config::broadcast_resolved_case_v3(MPI_COMM_NULL, 0, nullptr));
  } catch (const Error &) {
    null_rejected = true;
  }
  HUNDUN_CHECK(null_rejected);
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
      throw std::runtime_error("unknown test mode");
    }
  });
}

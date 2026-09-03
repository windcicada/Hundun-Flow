// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "src/app_case_config_broadcast_detail.hpp"
#include "tests/support/app_case_config_broadcast_test.hpp"

#include "hundun/cfg_case_config.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_mpi_environment.hpp"
#include "tests/support/test_main.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using hundun::application::broadcast_case_config;
using hundun::config::CaseConfig;
using hundun::runtime::Error;
using hundun::runtime::Int3;
namespace broadcast_detail = hundun::application::detail;

std::uint64_t double_bits(double value) {
  std::uint64_t result = 0;
  static_assert(sizeof(result) == sizeof(value));
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

bool same(Int3 lhs, Int3 rhs) {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

CaseConfig make_config(std::string name, bool with_expected_ranks,
                       bool with_process_grid, bool restart_read,
                       int communicator_size) {
  CaseConfig config{};
  config.schema_version = 1;
  config.case_name = std::move(name);
  if (with_expected_ranks) {
    config.expected_ranks = communicator_size;
  }
  if (with_process_grid) {
    config.process_grid = Int3{communicator_size, 1, 1};
  }
  config.mesh.cells = Int3{communicator_size + 7, 5, 3};
  config.mesh.origin_m = {0.125, -0.25, 0.375};
  config.mesh.length_m = {1.5, 2.5, 3.5};
  config.mesh.periodic = {true, false, true};
  config.time.dt_s = 0.03125;
  config.time.steps = 27;
  config.transport.velocity_m_per_s = {1.25, -2.5, 0.625};
  config.transport.diffusivity_m2_per_s = 0.015625;
  config.initial_condition = "sine_x";
  config.restart.read = restart_read;
  if (restart_read) {
    config.restart.read_directory = "Restart/step00000012";
  }
  config.restart.write_directory = restart_read ? "Restart.resumed" : "Restart";
  config.output.directory = restart_read ? "output.resumed" : "output";
  config.output.write_interval = 3;
  config.output.restart_interval = 4;
  return config;
}

void check_equal(const CaseConfig &actual, const CaseConfig &expected) {
  HUNDUN_CHECK(actual.schema_version == expected.schema_version);
  HUNDUN_CHECK(actual.case_name == expected.case_name);
  HUNDUN_CHECK(actual.expected_ranks == expected.expected_ranks);
  HUNDUN_CHECK(actual.process_grid.has_value() ==
               expected.process_grid.has_value());
  if (expected.process_grid.has_value()) {
    HUNDUN_CHECK(same(*actual.process_grid, *expected.process_grid));
  }
  HUNDUN_CHECK(same(actual.mesh.cells, expected.mesh.cells));
  HUNDUN_CHECK(double_bits(actual.mesh.origin_m.x) ==
               double_bits(expected.mesh.origin_m.x));
  HUNDUN_CHECK(double_bits(actual.mesh.origin_m.y) ==
               double_bits(expected.mesh.origin_m.y));
  HUNDUN_CHECK(double_bits(actual.mesh.origin_m.z) ==
               double_bits(expected.mesh.origin_m.z));
  HUNDUN_CHECK(double_bits(actual.mesh.length_m.x) ==
               double_bits(expected.mesh.length_m.x));
  HUNDUN_CHECK(double_bits(actual.mesh.length_m.y) ==
               double_bits(expected.mesh.length_m.y));
  HUNDUN_CHECK(double_bits(actual.mesh.length_m.z) ==
               double_bits(expected.mesh.length_m.z));
  HUNDUN_CHECK(actual.mesh.periodic == expected.mesh.periodic);
  HUNDUN_CHECK(double_bits(actual.time.dt_s) ==
               double_bits(expected.time.dt_s));
  HUNDUN_CHECK(actual.time.steps == expected.time.steps);
  HUNDUN_CHECK(double_bits(actual.transport.velocity_m_per_s.x) ==
               double_bits(expected.transport.velocity_m_per_s.x));
  HUNDUN_CHECK(double_bits(actual.transport.velocity_m_per_s.y) ==
               double_bits(expected.transport.velocity_m_per_s.y));
  HUNDUN_CHECK(double_bits(actual.transport.velocity_m_per_s.z) ==
               double_bits(expected.transport.velocity_m_per_s.z));
  HUNDUN_CHECK(double_bits(actual.transport.diffusivity_m2_per_s) ==
               double_bits(expected.transport.diffusivity_m2_per_s));
  HUNDUN_CHECK(actual.initial_condition == expected.initial_condition);
  HUNDUN_CHECK(actual.restart.read == expected.restart.read);
  HUNDUN_CHECK(actual.restart.read_directory ==
               expected.restart.read_directory);
  HUNDUN_CHECK(actual.restart.write_directory ==
               expected.restart.write_directory);
  HUNDUN_CHECK(actual.output.directory == expected.output.directory);
  HUNDUN_CHECK(actual.output.write_interval == expected.output.write_interval);
  HUNDUN_CHECK(actual.output.restart_interval ==
               expected.output.restart_interval);
}

void run_broadcast_case(MPI_Comm communicator, int root,
                        const CaseConfig &expected) {
  int rank = 0;
  HUNDUN_CHECK(MPI_Comm_rank(communicator, &rank) == MPI_SUCCESS);
  const CaseConfig actual = broadcast_case_config(
      communicator, root, rank == root ? &expected : nullptr);
  check_equal(actual, expected);
}

template <class Function>
void expect_uniform_error(MPI_Comm comparison_communicator,
                          Function &&function) {
  std::string message;
  try {
    function();
  } catch (const Error &error) {
    message = error.what();
  }
  HUNDUN_CHECK(!message.empty());
  HUNDUN_CHECK(message.size() < 512U);

  std::array<char, 512> local{};
  std::copy(message.begin(), message.end(), local.begin());
  int size = 0;
  HUNDUN_CHECK(MPI_Comm_size(comparison_communicator, &size) == MPI_SUCCESS);
  std::vector<char> gathered(static_cast<std::size_t>(size) * local.size());
  HUNDUN_CHECK(MPI_Allgather(local.data(), static_cast<int>(local.size()),
                             MPI_CHAR, gathered.data(),
                             static_cast<int>(local.size()), MPI_CHAR,
                             comparison_communicator) == MPI_SUCCESS);
  for (int process = 0; process < size; ++process) {
    const char *begin =
        gathered.data() + static_cast<std::size_t>(process) * local.size();
    HUNDUN_CHECK(std::string(begin) == message);
  }
}

void run_full() {
  int argc = 1;
  char program[] = "test_case_config_broadcast";
  char *arguments[] = {program, nullptr};
  char **argv = arguments;
  hundun::runtime::MpiEnvironment environment(argc, argv);

  int world_rank = 0;
  int world_size = 0;
  HUNDUN_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &world_rank) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &world_size) == MPI_SUCCESS);
  HUNDUN_CHECK(world_size == 4);

  const CaseConfig first = make_config("full-members", true, true, false, 4);
  run_broadcast_case(MPI_COMM_WORLD, 0, first);

  const CaseConfig second = make_config("", false, true, true, 4);
  run_broadcast_case(MPI_COMM_WORLD, 2, second);

  const CaseConfig automatic_grid =
      make_config("automatic-grid", true, false, false, 4);
  run_broadcast_case(MPI_COMM_WORLD, 3, automatic_grid);

  MPI_Comm split = MPI_COMM_NULL;
  HUNDUN_CHECK(MPI_Comm_split(MPI_COMM_WORLD, world_rank % 2, world_rank,
                              &split) == MPI_SUCCESS);
  int split_size = 0;
  HUNDUN_CHECK(MPI_Comm_size(split, &split_size) == MPI_SUCCESS);
  HUNDUN_CHECK(split_size == 2);
  const CaseConfig split_config =
      make_config("split-relative-ranks", true, true, true, split_size);
  run_broadcast_case(split, 1, split_config);
  HUNDUN_CHECK(MPI_Comm_free(&split) == MPI_SUCCESS);
}

void run_mismatch() {
  int argc = 1;
  char program[] = "test_case_config_broadcast";
  char *arguments[] = {program, nullptr};
  char **argv = arguments;
  hundun::runtime::MpiEnvironment environment(argc, argv);

  int rank = 0;
  int size = 0;
  HUNDUN_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &size) == MPI_SUCCESS);
  HUNDUN_CHECK(size == 4);
  const CaseConfig valid = make_config("mismatch", true, true, false, size);

  expect_uniform_error(MPI_COMM_WORLD, [&] {
    const int local_root = rank == 0 ? 0 : 1;
    static_cast<void>(broadcast_case_config(
        MPI_COMM_WORLD, local_root, rank == local_root ? &valid : nullptr));
  });
  expect_uniform_error(MPI_COMM_WORLD, [&] {
    static_cast<void>(broadcast_case_config(MPI_COMM_WORLD, size, nullptr));
  });
  expect_uniform_error(MPI_COMM_WORLD, [&] {
    static_cast<void>(broadcast_case_config(MPI_COMM_WORLD, 0, nullptr));
  });

  CaseConfig grid_mismatch =
      make_config("grid-mismatch", false, false, false, size);
  grid_mismatch.process_grid = Int3{1, 1, 1};
  expect_uniform_error(MPI_COMM_WORLD, [&] {
    static_cast<void>(broadcast_case_config(
        MPI_COMM_WORLD, 0, rank == 0 ? &grid_mismatch : nullptr));
  });

  expect_uniform_error(MPI_COMM_WORLD, [&] {
    static_cast<void>(broadcast_case_config(MPI_COMM_NULL, 0, nullptr));
  });

  MPI_Comm local = MPI_COMM_NULL;
  MPI_Comm inter = MPI_COMM_NULL;
  const int color = rank < 2 ? 0 : 1;
  HUNDUN_CHECK(MPI_Comm_split(MPI_COMM_WORLD, color, rank, &local) ==
               MPI_SUCCESS);
  const int remote_leader = color == 0 ? 2 : 0;
  HUNDUN_CHECK(MPI_Intercomm_create(local, 0, MPI_COMM_WORLD, remote_leader,
                                    917, &inter) == MPI_SUCCESS);
  expect_uniform_error(MPI_COMM_WORLD, [&] {
    static_cast<void>(broadcast_case_config(inter, 0, nullptr));
  });
  HUNDUN_CHECK(MPI_Comm_free(&inter) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Comm_free(&local) == MPI_SUCCESS);
}

void run_path_failures() {
  int argc = 1;
  char program[] = "test_case_config_broadcast";
  char *arguments[] = {program, nullptr};
  char **argv = arguments;
  hundun::runtime::MpiEnvironment environment(argc, argv);

  int rank = 0;
  int size = 0;
  HUNDUN_CHECK(MPI_Comm_rank(MPI_COMM_WORLD, &rank) == MPI_SUCCESS);
  HUNDUN_CHECK(MPI_Comm_size(MPI_COMM_WORLD, &size) == MPI_SUCCESS);
  HUNDUN_CHECK(size == 4);

  const CaseConfig without_read_directory =
      make_config("path-preparation-failure", true, true, false, size);
  expect_uniform_error(MPI_COMM_WORLD, [&] {
    static_cast<void>(broadcast_detail::broadcast_case_config_with_path_failure(
        MPI_COMM_WORLD, 2, rank == 2 ? &without_read_directory : nullptr,
        {broadcast_detail::CaseConfigPathTransfer::restart_write_directory,
         broadcast_detail::CaseConfigPathFailurePhase::preparation, 2}));
  });

  const CaseConfig with_read_directory =
      make_config("path-reconstruction-failure", true, true, true, size);
  expect_uniform_error(MPI_COMM_WORLD, [&] {
    static_cast<void>(broadcast_detail::broadcast_case_config_with_path_failure(
        MPI_COMM_WORLD, 0, rank == 0 ? &with_read_directory : nullptr,
        {broadcast_detail::CaseConfigPathTransfer::restart_read_directory,
         broadcast_detail::CaseConfigPathFailurePhase::reconstruction, 3}));
  });
}

int run_finalized(int argc, char **argv) {
  int provided = MPI_THREAD_SINGLE;
  HUNDUN_CHECK(MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided) ==
               MPI_SUCCESS);
  HUNDUN_CHECK(provided >= MPI_THREAD_FUNNELED);
  HUNDUN_CHECK(MPI_Finalize() == MPI_SUCCESS);
  bool rejected = false;
  try {
    static_cast<void>(broadcast_case_config(MPI_COMM_WORLD, 0, nullptr));
  } catch (const Error &error) {
    HUNDUN_CHECK(!std::string(error.what()).empty());
    rejected = true;
  }
  HUNDUN_CHECK(rejected);
  return 0;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2 || argv[1] == nullptr) {
    return 2;
  }
  const std::string_view mode(argv[1]);
  if (mode == "finalized") {
    return hundun::test::run(
        [&] { static_cast<void>(run_finalized(argc, argv)); });
  }
  if (mode == "full") {
    return hundun::test::run(run_full);
  }
  if (mode == "mismatch") {
    return hundun::test::run(run_mismatch);
  }
  if (mode == "path_failures") {
    return hundun::test::run(run_path_failures);
  }
  return 2;
}

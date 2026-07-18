// SPDX-License-Identifier: Apache-2.0

#include "applications/hundun/case_config_broadcast.hpp"

#include "hundun/config/case_config_loader.hpp"
#include "hundun/runtime/error.hpp"
#include "mpi_error.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <utility>

namespace hundun::application {
namespace {

using runtime::Error;

void broadcast_bytes(MPI_Comm comm, int root, char *bytes,
                     std::uint64_t byte_count, std::string_view operation) {
  constexpr std::uint64_t chunk_limit =
      static_cast<std::uint64_t>(std::numeric_limits<int>::max());
  std::uint64_t offset = 0;
  while (offset < byte_count) {
    const std::uint64_t remaining = byte_count - offset;
    const int count = static_cast<int>(std::min(remaining, chunk_limit));
    runtime::detail::check_mpi(
        MPI_Bcast(bytes + static_cast<std::size_t>(offset), count, MPI_BYTE,
                  root, comm),
        operation);
    offset += static_cast<std::uint64_t>(count);
  }
}

[[noreturn]] void throw_uniform_error(MPI_Comm comm, int rank, int size,
                                      bool local_ok,
                                      std::string_view local_message) {
  const int local_failure = local_ok ? size : rank;
  int failing_rank = size;
  runtime::detail::check_mpi(
      MPI_Allreduce(&local_failure, &failing_rank, 1, MPI_INT, MPI_MIN, comm),
      "MPI_Allreduce case-config failure rank");
  if (failing_rank == size) {
    throw Error("case-config collective failure state is inconsistent");
  }

  std::uint64_t length = 0;
  if (rank == failing_rank) {
    length = static_cast<std::uint64_t>(local_message.size());
  }
  runtime::detail::check_mpi(
      MPI_Bcast(&length, 1, MPI_UINT64_T, failing_rank, comm),
      "MPI_Bcast case-config error length");
  if (length == 0U || length > static_cast<std::uint64_t>(
                                   std::numeric_limits<std::size_t>::max())) {
    throw Error("case-config collective failure has an invalid message");
  }

  std::string message;
  bool allocation_ok = true;
  try {
    message.resize(static_cast<std::size_t>(length));
  } catch (...) {
    allocation_ok = false;
  }
  const int local_allocation_ok = allocation_ok ? 1 : 0;
  int every_allocation_ok = 0;
  runtime::detail::check_mpi(MPI_Allreduce(&local_allocation_ok,
                                           &every_allocation_ok, 1, MPI_INT,
                                           MPI_MIN, comm),
                             "MPI_Allreduce case-config error allocation");
  if (every_allocation_ok == 0) {
    throw Error("unable to allocate case-config collective error message");
  }
  if (rank == failing_rank) {
    std::copy(local_message.begin(), local_message.end(), message.begin());
  }
  broadcast_bytes(comm, failing_rank, message.data(), length,
                  "MPI_Bcast case-config error bytes");
  throw Error(std::move(message));
}

void converge_or_throw(MPI_Comm comm, int rank, int size, bool local_ok,
                       std::string_view local_message) {
  const int local_failure = local_ok ? size : rank;
  int failing_rank = size;
  runtime::detail::check_mpi(
      MPI_Allreduce(&local_failure, &failing_rank, 1, MPI_INT, MPI_MIN, comm),
      "MPI_Allreduce case-config status");
  if (failing_rank != size) {
    throw_uniform_error(comm, rank, size, local_ok, local_message);
  }
}

void broadcast_int(MPI_Comm comm, int root, int &value,
                   std::string_view operation) {
  runtime::detail::check_mpi(MPI_Bcast(&value, 1, MPI_INT, root, comm),
                             operation);
}

void broadcast_double(MPI_Comm comm, int root, double &value,
                      std::string_view operation) {
  runtime::detail::check_mpi(MPI_Bcast(&value, 1, MPI_DOUBLE, root, comm),
                             operation);
}

void broadcast_int3(MPI_Comm comm, int root, runtime::Int3 &value,
                    std::string_view operation) {
  std::array<int, 3> values{value.x, value.y, value.z};
  runtime::detail::check_mpi(MPI_Bcast(values.data(),
                                       static_cast<int>(values.size()), MPI_INT,
                                       root, comm),
                             operation);
  value = {values[0], values[1], values[2]};
}

void broadcast_real3(MPI_Comm comm, int root, runtime::Real3 &value,
                     std::string_view operation) {
  std::array<double, 3> values{value.x, value.y, value.z};
  runtime::detail::check_mpi(MPI_Bcast(values.data(),
                                       static_cast<int>(values.size()),
                                       MPI_DOUBLE, root, comm),
                             operation);
  value = {values[0], values[1], values[2]};
}

void broadcast_bool(MPI_Comm comm, int root, bool &value,
                    std::string_view operation) {
  runtime::detail::check_mpi(MPI_Bcast(&value, 1, MPI_CXX_BOOL, root, comm),
                             operation);
}

void broadcast_string(MPI_Comm comm, int root, int rank, std::string &value,
                      std::string_view operation) {
  std::uint64_t length = 0;
  if (rank == root) {
    if (value.size() > std::numeric_limits<std::uint64_t>::max()) {
      throw Error("case-config string exceeds the uint64 wire domain");
    }
    length = static_cast<std::uint64_t>(value.size());
  }
  runtime::detail::check_mpi(MPI_Bcast(&length, 1, MPI_UINT64_T, root, comm),
                             operation);

  bool allocation_ok = length <= static_cast<std::uint64_t>(
                                     std::numeric_limits<std::size_t>::max());
  if (allocation_ok) {
    try {
      value.resize(static_cast<std::size_t>(length));
    } catch (...) {
      allocation_ok = false;
    }
  }
  int size = 0;
  runtime::detail::check_mpi(MPI_Comm_size(comm, &size),
                             "MPI_Comm_size case-config string");
  converge_or_throw(comm, rank, size, allocation_ok,
                    "unable to allocate case-config string");
  broadcast_bytes(comm, root, value.data(), length, operation);
}

void broadcast_path(MPI_Comm comm, int root, int rank,
                    std::filesystem::path &path, std::string_view operation) {
  std::string text = rank == root ? path.generic_string() : std::string{};
  broadcast_string(comm, root, rank, text, operation);
  path = std::filesystem::path(std::move(text));
}

std::uint64_t checked_grid_product(runtime::Int3 grid) {
  const std::array<int, 3> values{grid.x, grid.y, grid.z};
  std::uint64_t product = 1U;
  for (const int value : values) {
    if (value <= 0) {
      throw Error("case-config process grid contains a nonpositive value");
    }
    const auto factor = static_cast<std::uint64_t>(value);
    if (product > std::numeric_limits<std::uint64_t>::max() / factor) {
      throw Error("case-config process-grid product overflows uint64");
    }
    product *= factor;
  }
  return product;
}

std::string exception_message_or_fallback(const std::exception &error,
                                          std::string_view fallback) {
  const char *message = error.what();
  return message != nullptr && *message != '\0' ? std::string(message)
                                                : std::string(fallback);
}

} // namespace

config::CaseConfig
broadcast_case_config(MPI_Comm comm, int root,
                      const config::CaseConfig *root_config) {
  runtime::detail::require_mpi_active("broadcast typed case configuration");
  if (comm == MPI_COMM_NULL) {
    throw Error("typed case configuration requires a valid intracommunicator");
  }
  int is_intercommunicator = 0;
  runtime::detail::check_mpi(MPI_Comm_test_inter(comm, &is_intercommunicator),
                             "MPI_Comm_test_inter case-config broadcast");
  if (is_intercommunicator != 0) {
    throw Error("typed case configuration requires an intracommunicator");
  }

  int rank = 0;
  int size = 0;
  runtime::detail::check_mpi(MPI_Comm_rank(comm, &rank),
                             "MPI_Comm_rank case-config broadcast");
  runtime::detail::check_mpi(MPI_Comm_size(comm, &size),
                             "MPI_Comm_size case-config broadcast");

  int minimum_root = 0;
  int maximum_root = 0;
  runtime::detail::check_mpi(
      MPI_Allreduce(&root, &minimum_root, 1, MPI_INT, MPI_MIN, comm),
      "MPI_Allreduce case-config minimum root");
  runtime::detail::check_mpi(
      MPI_Allreduce(&root, &maximum_root, 1, MPI_INT, MPI_MAX, comm),
      "MPI_Allreduce case-config maximum root");
  if (minimum_root != maximum_root) {
    throw Error("case-config broadcast root differs across communicator ranks");
  }
  root = minimum_root;
  if (root < 0 || root >= size) {
    throw Error("case-config broadcast root is outside the communicator");
  }

  config::CaseConfig result{};
  bool preparation_ok = true;
  std::string preparation_message;
  if (rank == root) {
    try {
      if (root_config == nullptr) {
        throw Error("case-config broadcast root requires a configuration");
      }
      config::validate_case_config(*root_config);
      result = *root_config;
    } catch (const std::exception &error) {
      preparation_ok = false;
      preparation_message = exception_message_or_fallback(
          error, "case-config root preparation failed");
    } catch (...) {
      preparation_ok = false;
      preparation_message = "case-config root preparation failed";
    }
  }
  converge_or_throw(comm, rank, size, preparation_ok, preparation_message);

  broadcast_int(comm, root, result.schema_version,
                "MPI_Bcast case-config schema version");
  broadcast_string(comm, root, rank, result.case_name,
                   "MPI_Bcast case-config case name");

  int has_expected_ranks = result.expected_ranks.has_value() ? 1 : 0;
  broadcast_int(comm, root, has_expected_ranks,
                "MPI_Bcast case-config expected-ranks flag");
  if (has_expected_ranks != 0 && has_expected_ranks != 1) {
    throw Error("case-config expected-ranks flag is invalid");
  }
  int expected_ranks = result.expected_ranks.value_or(0);
  if (has_expected_ranks != 0) {
    broadcast_int(comm, root, expected_ranks,
                  "MPI_Bcast case-config expected ranks");
    result.expected_ranks = expected_ranks;
  } else {
    result.expected_ranks.reset();
  }

  int has_process_grid = result.process_grid.has_value() ? 1 : 0;
  broadcast_int(comm, root, has_process_grid,
                "MPI_Bcast case-config process-grid flag");
  if (has_process_grid != 0 && has_process_grid != 1) {
    throw Error("case-config process-grid flag is invalid");
  }
  runtime::Int3 process_grid = result.process_grid.value_or(runtime::Int3{});
  if (has_process_grid != 0) {
    broadcast_int3(comm, root, process_grid,
                   "MPI_Bcast case-config process grid");
    result.process_grid = process_grid;
  } else {
    result.process_grid.reset();
  }

  broadcast_int3(comm, root, result.mesh.cells,
                 "MPI_Bcast case-config mesh cells");
  broadcast_real3(comm, root, result.mesh.origin_m,
                  "MPI_Bcast case-config mesh origin");
  broadcast_real3(comm, root, result.mesh.length_m,
                  "MPI_Bcast case-config mesh length");
  for (bool &periodic : result.mesh.periodic) {
    broadcast_bool(comm, root, periodic, "MPI_Bcast case-config periodic flag");
  }
  broadcast_double(comm, root, result.time.dt_s,
                   "MPI_Bcast case-config time step");
  broadcast_int(comm, root, result.time.steps,
                "MPI_Bcast case-config target step");
  broadcast_real3(comm, root, result.transport.velocity_m_per_s,
                  "MPI_Bcast case-config transport velocity");
  broadcast_double(comm, root, result.transport.diffusivity_m2_per_s,
                   "MPI_Bcast case-config diffusivity");
  broadcast_string(comm, root, rank, result.initial_condition,
                   "MPI_Bcast case-config initial condition");

  broadcast_bool(comm, root, result.restart.read,
                 "MPI_Bcast case-config restart read");
  int has_read_directory = result.restart.read_directory.has_value() ? 1 : 0;
  broadcast_int(comm, root, has_read_directory,
                "MPI_Bcast case-config read-directory flag");
  if (has_read_directory != 0 && has_read_directory != 1) {
    throw Error("case-config read-directory flag is invalid");
  }
  std::filesystem::path read_directory =
      result.restart.read_directory.value_or(std::filesystem::path{});
  if (has_read_directory != 0) {
    broadcast_path(comm, root, rank, read_directory,
                   "MPI_Bcast case-config read directory");
    result.restart.read_directory = std::move(read_directory);
  } else {
    result.restart.read_directory.reset();
  }
  broadcast_path(comm, root, rank, result.restart.write_directory,
                 "MPI_Bcast case-config write directory");
  broadcast_path(comm, root, rank, result.output.directory,
                 "MPI_Bcast case-config output directory");
  broadcast_int(comm, root, result.output.write_interval,
                "MPI_Bcast case-config output interval");
  broadcast_int(comm, root, result.output.restart_interval,
                "MPI_Bcast case-config restart interval");

  bool validation_ok = true;
  std::string validation_message;
  try {
    config::validate_case_config(result);
    if (result.process_grid.has_value() &&
        checked_grid_product(*result.process_grid) !=
            static_cast<std::uint64_t>(size)) {
      throw Error(
          "case-config process-grid product does not equal communicator size");
    }
  } catch (const std::exception &error) {
    validation_ok = false;
    validation_message = exception_message_or_fallback(
        error, "broadcast case-config validation failed");
  } catch (...) {
    validation_ok = false;
    validation_message = "broadcast case-config validation failed";
  }
  converge_or_throw(comm, rank, size, validation_ok, validation_message);
  return result;
}

} // namespace hundun::application

// SPDX-License-Identifier: Apache-2.0

#include "applications/hundun/case_config_broadcast.hpp"
#include "applications/hundun/detail/case_config_broadcast_test.hpp"

#include "hundun/config/case_config_loader.hpp"
#include "hundun/config/resolved_case_loader.hpp"
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
#include <vector>

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

bool should_inject_path_failure(
    const detail::CaseConfigPathFailureInjection *injection,
    detail::CaseConfigPathTransfer transfer,
    detail::CaseConfigPathFailurePhase phase, int rank) noexcept {
  return injection != nullptr && injection->transfer == transfer &&
         injection->phase == phase && injection->rank == rank;
}

template <class RootPathSource>
std::filesystem::path
broadcast_path(MPI_Comm comm, int root, int rank, int size,
               detail::CaseConfigPathTransfer transfer,
               RootPathSource &&root_path_source, std::string_view operation,
               const detail::CaseConfigPathFailureInjection *injection) {
  std::string text;
  bool preparation_ok = true;
  if (rank == root) {
    if (should_inject_path_failure(
            injection, transfer,
            detail::CaseConfigPathFailurePhase::preparation, rank)) {
      preparation_ok = false;
    } else {
      try {
        const std::filesystem::path prepared_path(root_path_source());
        text = prepared_path.generic_string();
        if (text.size() > std::numeric_limits<std::uint64_t>::max()) {
          throw Error("case-config path exceeds the uint64 wire domain");
        }
      } catch (...) {
        preparation_ok = false;
      }
    }
  }
  converge_or_throw(comm, rank, size, preparation_ok,
                    "unable to prepare case-config path");

  broadcast_string(comm, root, rank, text, operation);

  std::filesystem::path reconstructed_path;
  bool reconstruction_ok = true;
  if (should_inject_path_failure(
          injection, transfer,
          detail::CaseConfigPathFailurePhase::reconstruction, rank)) {
    reconstruction_ok = false;
  } else {
    try {
      reconstructed_path = std::filesystem::path(text);
    } catch (...) {
      reconstruction_ok = false;
    }
  }
  converge_or_throw(comm, rank, size, reconstruction_ok,
                    "unable to reconstruct case-config path");
  return reconstructed_path;
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

namespace {

config::CaseConfig broadcast_case_config_impl(
    MPI_Comm comm, int root, const config::CaseConfig *root_config,
    const detail::CaseConfigPathFailureInjection *path_failure_injection) {
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
  if (has_read_directory != 0) {
    std::filesystem::path read_directory = broadcast_path(
        comm, root, rank, size,
        detail::CaseConfigPathTransfer::restart_read_directory,
        [&result]() -> const std::filesystem::path & {
          return result.restart.read_directory.value();
        },
        "MPI_Bcast case-config read directory", path_failure_injection);
    result.restart.read_directory.emplace(std::move(read_directory));
  } else {
    result.restart.read_directory.reset();
  }
  std::filesystem::path write_directory = broadcast_path(
      comm, root, rank, size,
      detail::CaseConfigPathTransfer::restart_write_directory,
      [&result]() -> const std::filesystem::path & {
        return result.restart.write_directory;
      },
      "MPI_Bcast case-config write directory", path_failure_injection);
  result.restart.write_directory = std::move(write_directory);
  std::filesystem::path output_directory = broadcast_path(
      comm, root, rank, size, detail::CaseConfigPathTransfer::output_directory,
      [&result]() -> const std::filesystem::path & {
        return result.output.directory;
      },
      "MPI_Bcast case-config output directory", path_failure_injection);
  result.output.directory = std::move(output_directory);
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

} // namespace

config::CaseConfig
broadcast_case_config(MPI_Comm comm, int root,
                      const config::CaseConfig *root_config) {
  return broadcast_case_config_impl(comm, root, root_config, nullptr);
}

namespace {

struct CollectiveIdentity {
  int rank;
  int size;
  int root;
};

CollectiveIdentity validate_resolved_collective(MPI_Comm comm, int root) {
  runtime::detail::require_mpi_active("broadcast resolved case");
  if (comm == MPI_COMM_NULL) {
    throw Error("resolved-case broadcast requires a valid intracommunicator");
  }
  int is_intercommunicator = 0;
  runtime::detail::check_mpi(
      MPI_Comm_test_inter(comm, &is_intercommunicator),
      "MPI_Comm_test_inter resolved-case broadcast");
  if (is_intercommunicator != 0) {
    throw Error("resolved-case broadcast requires an intracommunicator");
  }

  CollectiveIdentity identity{};
  runtime::detail::check_mpi(MPI_Comm_rank(comm, &identity.rank),
                             "MPI_Comm_rank resolved-case broadcast");
  runtime::detail::check_mpi(MPI_Comm_size(comm, &identity.size),
                             "MPI_Comm_size resolved-case broadcast");
  int minimum_root = root;
  int maximum_root = root;
  runtime::detail::check_mpi(
      MPI_Allreduce(&root, &minimum_root, 1, MPI_INT, MPI_MIN, comm),
      "MPI_Allreduce resolved-case minimum root");
  runtime::detail::check_mpi(
      MPI_Allreduce(&root, &maximum_root, 1, MPI_INT, MPI_MAX, comm),
      "MPI_Allreduce resolved-case maximum root");
  if (minimum_root != maximum_root) {
    throw Error("resolved-case broadcast root differs across communicator ranks");
  }
  if (minimum_root < 0 || minimum_root >= identity.size) {
    throw Error("resolved-case broadcast root is outside the communicator");
  }
  identity.root = minimum_root;
  return identity;
}

void broadcast_uint64(MPI_Comm comm, int root, std::uint64_t &value,
                      std::string_view operation) {
  runtime::detail::check_mpi(MPI_Bcast(&value, 1, MPI_UINT64_T, root, comm),
                             operation);
}

template <class Enum>
void broadcast_enum(MPI_Comm comm, int root, Enum &value,
                    std::string_view operation) {
  int encoded = static_cast<int>(value);
  broadcast_int(comm, root, encoded, operation);
  value = static_cast<Enum>(encoded);
}

void broadcast_optional_double(MPI_Comm comm, int root,
                               std::optional<double> &value,
                               std::string_view flag_operation,
                               std::string_view value_operation) {
  int present = value.has_value() ? 1 : 0;
  broadcast_int(comm, root, present, flag_operation);
  if (present != 0 && present != 1) {
    throw Error("resolved-case optional-number flag is invalid");
  }
  if (present == 0) {
    value.reset();
    return;
  }
  double number = value.value_or(0.0);
  broadcast_double(comm, root, number, value_operation);
  value = number;
}

std::filesystem::path broadcast_flow_path(MPI_Comm comm, int root, int rank,
                                          int size,
                                          const std::filesystem::path &path,
                                          std::string_view operation) {
  std::string text;
  bool preparation_ok = true;
  std::string preparation_message;
  if (rank == root) {
    try {
      text = path.generic_string();
    } catch (const std::exception &error) {
      preparation_ok = false;
      preparation_message = exception_message_or_fallback(
          error, "unable to prepare resolved-case path");
    } catch (...) {
      preparation_ok = false;
      preparation_message = "unable to prepare resolved-case path";
    }
  }
  converge_or_throw(comm, rank, size, preparation_ok, preparation_message);
  broadcast_string(comm, root, rank, text, operation);

  std::filesystem::path result;
  bool reconstruction_ok = true;
  std::string reconstruction_message;
  try {
    result = std::filesystem::path(text);
  } catch (const std::exception &error) {
    reconstruction_ok = false;
    reconstruction_message = exception_message_or_fallback(
        error, "unable to reconstruct resolved-case path");
  } catch (...) {
    reconstruction_ok = false;
    reconstruction_message = "unable to reconstruct resolved-case path";
  }
  converge_or_throw(comm, rank, size, reconstruction_ok,
                    reconstruction_message);
  return result;
}

template <class Element>
void broadcast_vector_count(MPI_Comm comm, int root, int rank, int size,
                            std::vector<Element> &values,
                            std::string_view operation) {
  std::uint64_t count = 0U;
  if (rank == root) {
    count = static_cast<std::uint64_t>(values.size());
  }
  broadcast_uint64(comm, root, count, operation);
  bool allocation_ok =
      count <= static_cast<std::uint64_t>(std::numeric_limits<int>::max()) &&
      count <= static_cast<std::uint64_t>(
                   std::numeric_limits<std::size_t>::max());
  std::string allocation_message;
  if (allocation_ok) {
    try {
      values.resize(static_cast<std::size_t>(count));
    } catch (const std::exception &error) {
      allocation_ok = false;
      allocation_message = exception_message_or_fallback(
          error, "unable to allocate resolved-case vector");
    } catch (...) {
      allocation_ok = false;
      allocation_message = "unable to allocate resolved-case vector";
    }
  } else {
    allocation_message = "resolved-case vector count exceeds the schema domain";
  }
  converge_or_throw(comm, rank, size, allocation_ok, allocation_message);
}

void canonicalize_flow_order(config::FlowCaseConfig &flow) {
  std::sort(flow.scalars.begin(), flow.scalars.end(),
            [](const auto &left, const auto &right) {
              return left.name < right.name;
            });
  for (auto &boundary : flow.boundaries) {
    if (boundary.scalar_values.has_value()) {
      std::sort(boundary.scalar_values->begin(), boundary.scalar_values->end(),
                [](const auto &left, const auto &right) {
                  return left.name < right.name;
                });
    }
  }
  std::sort(flow.boundaries.begin(), flow.boundaries.end(),
            [](const auto &left, const auto &right) {
              return static_cast<int>(left.patch) <
                     static_cast<int>(right.patch);
            });
}

config::FlowCaseConfig broadcast_flow_case(
    MPI_Comm comm, const CollectiveIdentity &identity,
    const config::FlowCaseConfig *root_flow) {
  config::FlowCaseConfig flow{};
  bool preparation_ok = true;
  std::string preparation_message;
  if (identity.rank == identity.root) {
    try {
      flow = *root_flow;
      canonicalize_flow_order(flow);
    } catch (const std::exception &error) {
      preparation_ok = false;
      preparation_message = exception_message_or_fallback(
          error, "unable to prepare typed flow case");
    } catch (...) {
      preparation_ok = false;
      preparation_message = "unable to prepare typed flow case";
    }
  }
  converge_or_throw(comm, identity.rank, identity.size, preparation_ok,
                    preparation_message);

  broadcast_int(comm, identity.root, flow.schema_version,
                "MPI_Bcast flow schema version");
  broadcast_string(comm, identity.root, identity.rank, flow.case_name,
                   "MPI_Bcast flow case name");
  broadcast_enum(comm, identity.root, flow.simulation_type,
                 "MPI_Bcast flow simulation type");
  broadcast_enum(comm, identity.root, flow.density_model,
                 "MPI_Bcast flow density model");

  int has_expected_ranks = flow.resources.expected_ranks.has_value() ? 1 : 0;
  broadcast_int(comm, identity.root, has_expected_ranks,
                "MPI_Bcast flow expected-ranks flag");
  if (has_expected_ranks != 0 && has_expected_ranks != 1) {
    throw Error("flow expected-ranks flag is invalid");
  }
  int expected_ranks = flow.resources.expected_ranks.value_or(0);
  if (has_expected_ranks != 0) {
    broadcast_int(comm, identity.root, expected_ranks,
                  "MPI_Bcast flow expected ranks");
    flow.resources.expected_ranks = expected_ranks;
  } else {
    flow.resources.expected_ranks.reset();
  }

  int has_process_grid = flow.resources.process_grid.has_value() ? 1 : 0;
  broadcast_int(comm, identity.root, has_process_grid,
                "MPI_Bcast flow process-grid flag");
  if (has_process_grid != 0 && has_process_grid != 1) {
    throw Error("flow process-grid flag is invalid");
  }
  runtime::Int3 process_grid =
      flow.resources.process_grid.value_or(runtime::Int3{});
  if (has_process_grid != 0) {
    broadcast_int3(comm, identity.root, process_grid,
                   "MPI_Bcast flow process grid");
    flow.resources.process_grid = process_grid;
  } else {
    flow.resources.process_grid.reset();
  }

  broadcast_int3(comm, identity.root, flow.mesh.cells,
                 "MPI_Bcast flow mesh cells");
  broadcast_real3(comm, identity.root, flow.mesh.origin_m,
                  "MPI_Bcast flow mesh origin");
  broadcast_real3(comm, identity.root, flow.mesh.length_m,
                  "MPI_Bcast flow mesh length");
  broadcast_enum(comm, identity.root, flow.mesh.mapping,
                 "MPI_Bcast flow mesh mapping");
  int has_warp = flow.mesh.warp_amplitude.has_value() ? 1 : 0;
  broadcast_int(comm, identity.root, has_warp,
                "MPI_Bcast flow warp-amplitude flag");
  if (has_warp != 0 && has_warp != 1) {
    throw Error("flow warp-amplitude flag is invalid");
  }
  runtime::Real3 warp = flow.mesh.warp_amplitude.value_or(runtime::Real3{});
  if (has_warp != 0) {
    broadcast_real3(comm, identity.root, warp,
                    "MPI_Bcast flow warp amplitude");
    flow.mesh.warp_amplitude = warp;
  } else {
    flow.mesh.warp_amplitude.reset();
  }

  broadcast_enum(comm, identity.root, flow.time.mode,
                 "MPI_Bcast flow time mode");
  broadcast_int(comm, identity.root, flow.time.steps,
                "MPI_Bcast flow steps");
  broadcast_double(comm, identity.root, flow.time.initial_dt_s,
                   "MPI_Bcast flow initial dt");
  broadcast_double(comm, identity.root, flow.time.min_dt_s,
                   "MPI_Bcast flow minimum dt");
  broadcast_double(comm, identity.root, flow.time.max_dt_s,
                   "MPI_Bcast flow maximum dt");
  broadcast_double(comm, identity.root, flow.time.cfl_target,
                   "MPI_Bcast flow CFL target");
  broadcast_double(comm, identity.root, flow.time.diffusion_number_target,
                   "MPI_Bcast flow diffusion-number target");
  broadcast_double(comm, identity.root, flow.time.growth_factor,
                   "MPI_Bcast flow growth factor");
  broadcast_double(comm, identity.root, flow.time.retry_factor,
                   "MPI_Bcast flow retry factor");
  broadcast_int(comm, identity.root, flow.time.max_retries,
                "MPI_Bcast flow maximum retries");

  broadcast_double(comm, identity.root, flow.physics.rho_ref_kg_per_m3,
                   "MPI_Bcast flow reference density");
  broadcast_double(comm, identity.root,
                   flow.physics.dynamic_viscosity_pa_s,
                   "MPI_Bcast flow dynamic viscosity");
  broadcast_double(comm, identity.root,
                   flow.physics.inlet_consistency_rtol,
                   "MPI_Bcast flow inlet consistency tolerance");
  broadcast_optional_double(comm, identity.root,
                            flow.physics.cp_J_per_kg_K,
                            "MPI_Bcast flow cp flag", "MPI_Bcast flow cp");
  broadcast_optional_double(
      comm, identity.root, flow.physics.gas_constant_J_per_kg_K,
      "MPI_Bcast flow gas-constant flag", "MPI_Bcast flow gas constant");
  broadcast_optional_double(
      comm, identity.root, flow.physics.thermodynamic_pressure_pa,
      "MPI_Bcast flow thermodynamic-pressure flag",
      "MPI_Bcast flow thermodynamic pressure");

  broadcast_vector_count(comm, identity.root, identity.rank, identity.size,
                         flow.scalars, "MPI_Bcast flow scalar count");
  for (auto &scalar : flow.scalars) {
    broadcast_string(comm, identity.root, identity.rank, scalar.name,
                     "MPI_Bcast flow scalar name");
    broadcast_double(comm, identity.root, scalar.diffusivity_m2_per_s,
                     "MPI_Bcast flow scalar diffusivity");
  }

  for (auto &boundary : flow.boundaries) {
    broadcast_enum(comm, identity.root, boundary.patch,
                   "MPI_Bcast flow boundary patch");
    broadcast_enum(comm, identity.root, boundary.type,
                   "MPI_Bcast flow boundary type");

    int has_velocity = boundary.velocity_m_per_s.has_value() ? 1 : 0;
    broadcast_int(comm, identity.root, has_velocity,
                  "MPI_Bcast flow inlet-velocity flag");
    if (has_velocity != 0 && has_velocity != 1) {
      throw Error("flow inlet-velocity flag is invalid");
    }
    runtime::Real3 velocity =
        boundary.velocity_m_per_s.value_or(runtime::Real3{});
    if (has_velocity != 0) {
      broadcast_real3(comm, identity.root, velocity,
                      "MPI_Bcast flow inlet velocity");
      boundary.velocity_m_per_s = velocity;
    } else {
      boundary.velocity_m_per_s.reset();
    }

    int has_authority = boundary.thermal_authority.has_value() ? 1 : 0;
    broadcast_int(comm, identity.root, has_authority,
                  "MPI_Bcast flow thermal-authority flag");
    if (has_authority != 0 && has_authority != 1) {
      throw Error("flow thermal-authority flag is invalid");
    }
    config::InletThermalAuthority authority =
        boundary.thermal_authority.value_or(
            config::InletThermalAuthority::temperature);
    if (has_authority != 0) {
      broadcast_enum(comm, identity.root, authority,
                     "MPI_Bcast flow thermal authority");
      boundary.thermal_authority = authority;
    } else {
      boundary.thermal_authority.reset();
    }
    broadcast_optional_double(comm, identity.root, boundary.temperature_K,
                              "MPI_Bcast flow inlet-temperature flag",
                              "MPI_Bcast flow inlet temperature");
    broadcast_optional_double(comm, identity.root,
                              boundary.enthalpy_J_per_kg,
                              "MPI_Bcast flow inlet-enthalpy flag",
                              "MPI_Bcast flow inlet enthalpy");
    broadcast_optional_double(comm, identity.root,
                              boundary.density_kg_per_m3,
                              "MPI_Bcast flow inlet-density flag",
                              "MPI_Bcast flow inlet density");

    int has_scalar_values = boundary.scalar_values.has_value() ? 1 : 0;
    broadcast_int(comm, identity.root, has_scalar_values,
                  "MPI_Bcast flow inlet-scalar-values flag");
    if (has_scalar_values != 0 && has_scalar_values != 1) {
      throw Error("flow inlet-scalar-values flag is invalid");
    }
    if (has_scalar_values != 0) {
      if (!boundary.scalar_values.has_value()) {
        boundary.scalar_values.emplace();
      }
      broadcast_vector_count(comm, identity.root, identity.rank,
                             identity.size, *boundary.scalar_values,
                             "MPI_Bcast flow inlet-scalar count");
      for (auto &value : *boundary.scalar_values) {
        broadcast_string(comm, identity.root, identity.rank, value.name,
                         "MPI_Bcast flow inlet-scalar name");
        broadcast_double(comm, identity.root, value.value,
                         "MPI_Bcast flow inlet-scalar value");
      }
    } else {
      boundary.scalar_values.reset();
    }
    broadcast_optional_double(
        comm, identity.root, boundary.pressure_perturbation_pa,
        "MPI_Bcast flow outlet-pressure flag",
        "MPI_Bcast flow outlet pressure");
  }

  broadcast_bool(comm, identity.root, flow.restart.read,
                 "MPI_Bcast flow restart-read flag");
  int has_read_directory = flow.restart.read_directory.has_value() ? 1 : 0;
  broadcast_int(comm, identity.root, has_read_directory,
                "MPI_Bcast flow restart-read-directory flag");
  if (has_read_directory != 0 && has_read_directory != 1) {
    throw Error("flow restart-read-directory flag is invalid");
  }
  if (has_read_directory != 0) {
    flow.restart.read_directory = broadcast_flow_path(
        comm, identity.root, identity.rank, identity.size,
        flow.restart.read_directory.value_or(std::filesystem::path{}),
        "MPI_Bcast flow restart read directory");
  } else {
    flow.restart.read_directory.reset();
  }
  flow.restart.write_directory = broadcast_flow_path(
      comm, identity.root, identity.rank, identity.size,
      flow.restart.write_directory, "MPI_Bcast flow restart write directory");
  broadcast_int(comm, identity.root, flow.restart.write_interval,
                "MPI_Bcast flow restart write interval");

  flow.diagnostics.directory = broadcast_flow_path(
      comm, identity.root, identity.rank, identity.size,
      flow.diagnostics.directory, "MPI_Bcast flow diagnostics directory");
  broadcast_int(comm, identity.root, flow.diagnostics.write_interval,
                "MPI_Bcast flow diagnostics interval");
  broadcast_bool(comm, identity.root, flow.diagnostics.write_mesh,
                 "MPI_Bcast flow diagnostics mesh flag");

  broadcast_bool(comm, identity.root, flow.performance.enabled,
                 "MPI_Bcast flow performance flag");
  flow.performance.directory = broadcast_flow_path(
      comm, identity.root, identity.rank, identity.size,
      flow.performance.directory, "MPI_Bcast flow performance directory");
  broadcast_int(comm, identity.root, flow.performance.warmup_steps,
                "MPI_Bcast flow performance warmup steps");
  broadcast_int(comm, identity.root, flow.performance.measured_steps,
                "MPI_Bcast flow performance measured steps");
  broadcast_int(comm, identity.root, flow.performance.repetitions,
                "MPI_Bcast flow performance repetitions");

  canonicalize_flow_order(flow);
  return flow;
}

} // namespace

config::ResolvedCase broadcast_resolved_case(
    MPI_Comm comm, int root, const config::ResolvedCase *root_case) {
  const CollectiveIdentity identity =
      validate_resolved_collective(comm, root);

  const bool pointer_ok =
      identity.rank == identity.root ? root_case != nullptr
                                     : root_case == nullptr;
  const std::string pointer_message =
      identity.rank == identity.root
          ? "resolved-case broadcast root requires case data"
          : "resolved-case broadcast non-root data must be null";
  converge_or_throw(comm, identity.rank, identity.size, pointer_ok,
                    pointer_message);

  int discriminant = 0;
  bool preparation_ok = true;
  std::string preparation_message;
  config::ResolvedCase prepared;
  if (identity.rank == identity.root) {
    try {
      prepared = *root_case;
      static_cast<void>(config::to_resolved_json(prepared));
      discriminant = static_cast<int>(prepared.index());
      if (discriminant < 0 || discriminant > 1) {
        throw Error("resolved-case discriminant is invalid");
      }
    } catch (const std::exception &error) {
      preparation_ok = false;
      preparation_message = exception_message_or_fallback(
          error, "resolved-case root preparation failed");
    } catch (...) {
      preparation_ok = false;
      preparation_message = "resolved-case root preparation failed";
    }
  }
  converge_or_throw(comm, identity.rank, identity.size, preparation_ok,
                    preparation_message);
  broadcast_int(comm, identity.root, discriminant,
                "MPI_Bcast resolved-case discriminant");
  if (discriminant != 0 && discriminant != 1) {
    throw Error("resolved-case discriminant is invalid");
  }

  config::ResolvedCase result;
  if (discriminant == 0) {
    const config::CaseConfig *root_config =
        identity.rank == identity.root
            ? &std::get<config::CaseConfig>(prepared)
            : nullptr;
    result = broadcast_case_config(comm, identity.root, root_config);
  } else {
    const config::FlowCaseConfig *root_flow =
        identity.rank == identity.root
            ? &std::get<config::FlowCaseConfig>(prepared)
            : nullptr;
    result = broadcast_flow_case(comm, identity, root_flow);
  }

  bool validation_ok = true;
  std::string validation_message;
  try {
    static_cast<void>(config::to_resolved_json(result));
    if (const auto *flow = std::get_if<config::FlowCaseConfig>(&result)) {
      if (flow->resources.expected_ranks.has_value() &&
          *flow->resources.expected_ranks != identity.size) {
        throw Error(
            "resolved-case expected_ranks does not equal communicator size");
      }
      if (flow->resources.process_grid.has_value() &&
          checked_grid_product(*flow->resources.process_grid) !=
              static_cast<std::uint64_t>(identity.size)) {
        throw Error(
            "resolved-case process-grid product does not equal communicator size");
      }
    }
  } catch (const std::exception &error) {
    validation_ok = false;
    validation_message = exception_message_or_fallback(
        error, "broadcast resolved-case validation failed");
  } catch (...) {
    validation_ok = false;
    validation_message = "broadcast resolved-case validation failed";
  }
  converge_or_throw(comm, identity.rank, identity.size, validation_ok,
                    validation_message);
  return result;
}

namespace detail {

config::CaseConfig broadcast_case_config_with_path_failure(
    MPI_Comm comm, int root, const config::CaseConfig *root_config,
    CaseConfigPathFailureInjection injection) {
  return broadcast_case_config_impl(comm, root, root_config, &injection);
}

} // namespace detail

} // namespace hundun::application

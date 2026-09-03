// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "app_resolved_case_v3_broadcast_detail.hpp"

#include "cfg_resolved_case_v3_loader_detail.hpp"

#include "hundun/cfg_resolved_case_v3_loader.hpp"
#include "hundun/rt_error.hpp"
#include "rt_mpi_error_detail.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>
#include <string_view>
#include <utility>

namespace hundun::config {
namespace {

using runtime::Error;

enum class BroadcastFault {
  none,
#if defined(HUNDUN_APPLICATION_ENABLE_TEST_ACCESS)
  version_tag,
  variant_tag,
  truncate_payload,
  mpi_operation
#endif
};

struct CollectiveIdentity final {
  int rank{};
  int size{};
  int root{};
};

std::string exception_message_or_fallback(const std::exception &error,
                                          std::string_view fallback) {
  const char *message = error.what();
  return message != nullptr && *message != '\0' ? std::string(message)
                                                : std::string(fallback);
}

void broadcast_bytes(MPI_Comm comm, int root, char *bytes,
                     std::uint64_t byte_count, std::string_view operation) {
  constexpr std::uint64_t chunk_limit =
      static_cast<std::uint64_t>(std::numeric_limits<int>::max());
  std::uint64_t offset = 0U;
  while (offset < byte_count) {
    const int count =
        static_cast<int>(std::min(byte_count - offset, chunk_limit));
    runtime::detail::check_mpi(
        MPI_Bcast(bytes + static_cast<std::size_t>(offset), count, MPI_BYTE,
                  root, comm),
        operation);
    offset += static_cast<std::uint64_t>(count);
  }
}

[[noreturn]] void throw_uniform_error(MPI_Comm comm,
                                      const CollectiveIdentity &identity,
                                      bool local_ok,
                                      std::string_view local_message) {
  const int local_failure = local_ok ? identity.size : identity.rank;
  int failing_rank = identity.size;
  runtime::detail::check_mpi(
      MPI_Allreduce(&local_failure, &failing_rank, 1, MPI_INT, MPI_MIN, comm),
      "MPI_Allreduce immersed-flow resolved-case failure rank");
  if (failing_rank == identity.size) {
    throw Error(
        "immersed-flow resolved-case collective failure state is inconsistent");
  }

  std::uint64_t length = 0U;
  if (identity.rank == failing_rank) {
    length = static_cast<std::uint64_t>(local_message.size());
  }
  runtime::detail::check_mpi(
      MPI_Bcast(&length, 1, MPI_UINT64_T, failing_rank, comm),
      "MPI_Bcast immersed-flow resolved-case error length");
  if (length == 0U || length > static_cast<std::uint64_t>(
                                   std::numeric_limits<std::size_t>::max())) {
    throw Error(
        "immersed-flow resolved-case collective failure has an invalid message");
  }

  std::string message;
  bool allocation_ok = true;
  try {
    message.resize(static_cast<std::size_t>(length));
  } catch (...) {
    allocation_ok = false;
  }
  const int local_allocation = allocation_ok ? 1 : 0;
  int every_allocation = 0;
  runtime::detail::check_mpi(
      MPI_Allreduce(&local_allocation, &every_allocation, 1, MPI_INT, MPI_MIN,
                    comm),
      "MPI_Allreduce immersed-flow resolved-case error allocation");
  if (every_allocation == 0) {
    throw Error(
        "unable to allocate immersed-flow resolved-case collective error message");
  }
  if (identity.rank == failing_rank) {
    std::copy(local_message.begin(), local_message.end(), message.begin());
  }
  broadcast_bytes(comm, failing_rank, message.data(), length,
                  "MPI_Bcast immersed-flow resolved-case error bytes");
  throw Error(std::move(message));
}

void converge_or_throw(MPI_Comm comm, const CollectiveIdentity &identity,
                       bool local_ok, std::string_view local_message) {
  const int local_failure = local_ok ? identity.size : identity.rank;
  int failing_rank = identity.size;
  runtime::detail::check_mpi(
      MPI_Allreduce(&local_failure, &failing_rank, 1, MPI_INT, MPI_MIN, comm),
      "MPI_Allreduce immersed-flow resolved-case status");
  if (failing_rank != identity.size) {
    throw_uniform_error(comm, identity, local_ok, local_message);
  }
}

CollectiveIdentity validate_collective(MPI_Comm comm, int root) {
  runtime::detail::require_mpi_active("broadcast immersed-flow resolved case");
  if (comm == MPI_COMM_NULL) {
    throw Error(
        "immersed-flow resolved-case broadcast requires a valid intracommunicator");
  }
  int is_intercommunicator = 0;
  runtime::detail::check_mpi(
      MPI_Comm_test_inter(comm, &is_intercommunicator),
      "MPI_Comm_test_inter immersed-flow resolved-case broadcast");
  if (is_intercommunicator != 0) {
    throw Error(
        "immersed-flow resolved-case broadcast requires an intracommunicator");
  }

  CollectiveIdentity identity{};
  runtime::detail::check_mpi(MPI_Comm_rank(comm, &identity.rank),
                             "MPI_Comm_rank immersed-flow resolved-case broadcast");
  runtime::detail::check_mpi(MPI_Comm_size(comm, &identity.size),
                             "MPI_Comm_size immersed-flow resolved-case broadcast");
  int minimum_root = root;
  int maximum_root = root;
  runtime::detail::check_mpi(
      MPI_Allreduce(&root, &minimum_root, 1, MPI_INT, MPI_MIN, comm),
      "MPI_Allreduce immersed-flow resolved-case minimum root");
  runtime::detail::check_mpi(
      MPI_Allreduce(&root, &maximum_root, 1, MPI_INT, MPI_MAX, comm),
      "MPI_Allreduce immersed-flow resolved-case maximum root");
  if (minimum_root != maximum_root) {
    throw Error("immersed-flow resolved-case broadcast root differs across "
                "communicator ranks");
  }
  if (minimum_root < 0 || minimum_root >= identity.size) {
    throw Error(
        "immersed-flow resolved-case broadcast root is outside the communicator");
  }
  identity.root = minimum_root;
  return identity;
}

std::uint64_t checked_grid_product(runtime::Int3 grid) {
  const std::array<int, 3> values{grid.x, grid.y, grid.z};
  std::uint64_t product = 1U;
  for (const int value : values) {
    if (value <= 0) {
      throw Error(
          "immersed-flow resolved-case process grid contains a nonpositive value");
    }
    const auto factor = static_cast<std::uint64_t>(value);
    if (product > std::numeric_limits<std::uint64_t>::max() / factor) {
      throw Error(
          "immersed-flow resolved-case process-grid product overflows uint64");
    }
    product *= factor;
  }
  return product;
}

const FlowCaseConfig *common_flow(const ResolvedCaseV3 &value) {
  if (const auto *flow = std::get_if<FlowCaseConfig>(&value)) {
    return flow;
  }
  if (const auto *immersed_flow = std::get_if<ImmersedFlowCaseConfig>(&value)) {
    return &immersed_flow->common_flow;
  }
  return nullptr;
}

void validate_rank_contract(const ResolvedCaseV3 &value,
                            const CollectiveIdentity &identity) {
  if (const auto *passive_scalar = std::get_if<CaseConfig>(&value)) {
    if (passive_scalar->expected_ranks.has_value() &&
        *passive_scalar->expected_ranks != identity.size) {
      throw Error("expected MPI rank count " +
                  std::to_string(*passive_scalar->expected_ranks) + ", got " +
                  std::to_string(identity.size));
    }
    if (passive_scalar->process_grid.has_value() &&
        checked_grid_product(*passive_scalar->process_grid) !=
            static_cast<std::uint64_t>(identity.size)) {
      throw Error("immersed-flow resolved-case process-grid product does not equal "
                  "communicator size");
    }
    return;
  }
  const FlowCaseConfig *flow = common_flow(value);
  if (flow == nullptr) {
    throw Error("immersed-flow resolved-case variant is invalid");
  }
  if (flow->resources.expected_ranks.has_value() &&
      *flow->resources.expected_ranks != identity.size) {
    throw Error("immersed-flow resolved-case expected_ranks does not equal "
                "communicator size");
  }
  if (flow->resources.process_grid.has_value() &&
      checked_grid_product(*flow->resources.process_grid) !=
          static_cast<std::uint64_t>(identity.size)) {
    throw Error("immersed-flow resolved-case process-grid product does not equal "
                "communicator size");
  }
}

ResolvedCaseV3 broadcast_impl(MPI_Comm comm, int root,
                              const ResolvedCaseV3 *root_case,
                              BroadcastFault fault) {
  const CollectiveIdentity identity = validate_collective(comm, root);

  const bool pointer_ok = identity.rank == identity.root ? root_case != nullptr
                                                         : root_case == nullptr;
  const std::string pointer_message =
      identity.rank == identity.root
          ? "immersed-flow resolved-case broadcast root requires case data"
          : "immersed-flow resolved-case broadcast non-root data must be null";
  converge_or_throw(comm, identity, pointer_ok, pointer_message);

  constexpr int protocol_version = 1;
  int transmitted_version = protocol_version;
  int variant_tag = 0;
  std::string payload;
  bool preparation_ok = true;
  std::string preparation_message;
  if (identity.rank == identity.root) {
    try {
      variant_tag = static_cast<int>(root_case->index());
      if (variant_tag < 0 || variant_tag > 2) {
        throw Error("immersed-flow resolved-case variant tag is invalid");
      }
      payload = to_resolved_json_v3(*root_case);
      if (payload.size() >
          static_cast<std::size_t>(std::numeric_limits<std::uint64_t>::max())) {
        throw Error(
            "immersed-flow resolved-case payload exceeds the uint64 wire domain");
      }
    } catch (const std::exception &error) {
      preparation_ok = false;
      preparation_message = exception_message_or_fallback(
          error, "immersed-flow resolved-case root preparation failed");
    } catch (...) {
      preparation_ok = false;
      preparation_message = "immersed-flow resolved-case root preparation failed";
    }
  }
  converge_or_throw(comm, identity, preparation_ok, preparation_message);

#if defined(HUNDUN_APPLICATION_ENABLE_TEST_ACCESS)
  if (fault == BroadcastFault::mpi_operation) {
    const bool operation_ok = identity.rank != 0;
    converge_or_throw(comm, identity, operation_ok,
                      "injected immersed-flow resolved-case MPI operation failure");
  }
  if (identity.rank == identity.root) {
    if (fault == BroadcastFault::version_tag) {
      transmitted_version = protocol_version + 1;
    } else if (fault == BroadcastFault::variant_tag) {
      variant_tag = (variant_tag + 1) % 3;
    }
  }
#else
  static_cast<void>(fault);
#endif

  runtime::detail::check_mpi(
      MPI_Bcast(&transmitted_version, 1, MPI_INT, identity.root, comm),
      "MPI_Bcast immersed-flow resolved-case protocol version");
  runtime::detail::check_mpi(
      MPI_Bcast(&variant_tag, 1, MPI_INT, identity.root, comm),
      "MPI_Bcast immersed-flow resolved-case variant tag");

  bool header_ok = transmitted_version == protocol_version &&
                   variant_tag >= 0 && variant_tag <= 2;
  std::string header_message;
  if (transmitted_version != protocol_version) {
    header_message = "immersed-flow resolved-case protocol version is unsupported";
  } else if (variant_tag < 0 || variant_tag > 2) {
    header_message = "immersed-flow resolved-case variant tag is invalid";
  }
  converge_or_throw(comm, identity, header_ok, header_message);

  std::uint64_t length = 0U;
  if (identity.rank == identity.root) {
    length = static_cast<std::uint64_t>(payload.size());
#if defined(HUNDUN_APPLICATION_ENABLE_TEST_ACCESS)
    if (fault == BroadcastFault::truncate_payload && length > 0U) {
      --length;
    }
#endif
  }
  runtime::detail::check_mpi(
      MPI_Bcast(&length, 1, MPI_UINT64_T, identity.root, comm),
      "MPI_Bcast immersed-flow resolved-case payload length");

  bool allocation_ok =
      length > 0U && length <= static_cast<std::uint64_t>(
                                   std::numeric_limits<std::size_t>::max());
  std::string allocation_message;
  if (allocation_ok) {
    try {
      payload.resize(static_cast<std::size_t>(length));
    } catch (const std::exception &error) {
      allocation_ok = false;
      allocation_message = exception_message_or_fallback(
          error, "unable to allocate immersed-flow resolved-case payload");
    } catch (...) {
      allocation_ok = false;
      allocation_message = "unable to allocate immersed-flow resolved-case payload";
    }
  } else {
    allocation_message =
        "immersed-flow resolved-case payload length is outside the host domain";
  }
  converge_or_throw(comm, identity, allocation_ok, allocation_message);
  broadcast_bytes(comm, identity.root, payload.data(), length,
                  "MPI_Bcast immersed-flow resolved-case payload bytes");

  ResolvedCaseV3 result;
  bool parse_ok = true;
  std::string parse_message;
  try {
    result = hundun::config::detail::parse_resolved_case_v3_json(
        payload, std::filesystem::path("immersed-flow-broadcast-case.json"));
    if (static_cast<int>(result.index()) != variant_tag) {
      throw Error(
          "immersed-flow resolved-case payload does not match its variant tag");
    }
    validate_rank_contract(result, identity);
  } catch (const std::exception &error) {
    parse_ok = false;
    parse_message = exception_message_or_fallback(
        error, "immersed-flow resolved-case payload validation failed");
  } catch (...) {
    parse_ok = false;
    parse_message = "immersed-flow resolved-case payload validation failed";
  }
  converge_or_throw(comm, identity, parse_ok, parse_message);
  return result;
}

} // namespace

ResolvedCaseV3 broadcast_resolved_case_v3(MPI_Comm comm, int root,
                                          const ResolvedCaseV3 *root_case) {
  return broadcast_impl(comm, root, root_case, BroadcastFault::none);
}

#if defined(HUNDUN_APPLICATION_ENABLE_TEST_ACCESS)
namespace detail {

ResolvedCaseV3
broadcast_resolved_case_v3_with_fault(MPI_Comm comm, int root,
                                      const ResolvedCaseV3 *root_case,
                                      ResolvedCaseV3BroadcastFault fault) {
  BroadcastFault internal = BroadcastFault::none;
  switch (fault) {
  case ResolvedCaseV3BroadcastFault::none:
    internal = BroadcastFault::none;
    break;
  case ResolvedCaseV3BroadcastFault::version_tag:
    internal = BroadcastFault::version_tag;
    break;
  case ResolvedCaseV3BroadcastFault::variant_tag:
    internal = BroadcastFault::variant_tag;
    break;
  case ResolvedCaseV3BroadcastFault::truncate_payload:
    internal = BroadcastFault::truncate_payload;
    break;
  case ResolvedCaseV3BroadcastFault::mpi_operation:
    internal = BroadcastFault::mpi_operation;
    break;
  }
  return broadcast_impl(comm, root, root_case, internal);
}

} // namespace detail
#endif

} // namespace hundun::config

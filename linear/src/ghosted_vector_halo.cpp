// SPDX-License-Identifier: Apache-2.0

#include "hundun/linear/ghosted_vector_halo.hpp"

#include "ghosted_vector_halo_detail.hpp"
#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "mpi_error.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <new>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace hundun::linear {
namespace {

constexpr int kMetadataTag = 17;
constexpr int kPayloadTag = 23;

bool same(runtime::Int3 left, runtime::Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool same(runtime::Box3 left, runtime::Box3 right) noexcept {
  return same(left.begin, right.begin) && same(left.end, right.end);
}

std::size_t checked_add(std::size_t left, std::size_t right,
                        const char* subject) {
  if (left > std::numeric_limits<std::size_t>::max() - right) {
    throw runtime::Error(std::string(subject) + " addition overflow");
  }
  return left + right;
}

std::size_t checked_size(std::uint64_t value, const char* subject) {
  if (value >
      static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
    throw runtime::Error(std::string(subject) +
                         " exceeds the local size range");
  }
  return static_cast<std::size_t>(value);
}

void require_collective(const runtime::MpiContext& context, bool local_ok,
                        std::string_view local_message) {
  const runtime::CollectiveStatus status = runtime::collective_status(
      context, local_ok, local_ok ? std::string_view{} : local_message);
  if (!status.ok) {
    throw runtime::Error("compact Buffer Halo collective failure: rank=" +
                         std::to_string(status.failing_rank) + " " +
                         status.message);
  }
}

struct PlanPeer {
  int rank{};
  std::vector<std::size_t> send_indices;
  std::vector<std::size_t> receive_indices;
};

struct BuiltPlan {
  std::vector<PlanPeer> peers;
  std::size_t send_value_count{};
  std::size_t receive_value_count{};
};

struct FailureRecord {
  std::int64_t category{};
  std::int64_t rank{-1};
  std::int64_t operation{};
  std::int64_t result{};
  std::int64_t peer{-1};
  std::int64_t chunk_offset{};
  std::int64_t chunk_count{};
  std::int64_t tag{-1};
};

static_assert(std::is_trivially_copyable_v<FailureRecord>);
static_assert(sizeof(FailureRecord) == 8U * sizeof(std::int64_t));

bool failed(FailureRecord record) noexcept { return record.category != 0; }

bool converge_failure(const runtime::MpiContext& context,
                      FailureRecord local, FailureRecord& global) noexcept {
  const int local_rank = failed(local) ? context.rank() : context.size();
  int failing_rank = context.size();
  if (MPI_Allreduce(&local_rank, &failing_rank, 1, MPI_INT, MPI_MIN,
                    context.comm()) != MPI_SUCCESS) {
    return false;
  }
  if (failing_rank == context.size()) {
    global = FailureRecord{};
    return true;
  }
  global = context.rank() == failing_rank ? local : FailureRecord{};
  if (MPI_Bcast(&global, static_cast<int>(sizeof(global)), MPI_BYTE,
                failing_rank, context.comm()) != MPI_SUCCESS) {
    return false;
  }
  global.rank = failing_rank;
  return true;
}

struct InjectionSelection {
  bool local_selected{};
};

enum class InjectionPhase {
  metadata_post,
  metadata_completion,
  payload_post,
  payload_completion
};

enum class InjectionMessage { mismatch, invalid, insufficient };

std::string_view injection_message(InjectionPhase phase,
                                   InjectionMessage message) noexcept {
  if (phase == InjectionPhase::metadata_post) {
    if (message == InjectionMessage::mismatch) {
      return "compact Buffer Halo request-ID post injection rank differs "
             "across ranks";
    }
    if (message == InjectionMessage::invalid) {
      return "compact Buffer Halo request-ID post injection rank is invalid";
    }
    return "compact Buffer Halo request-ID post injection requires at least "
           "two request descriptors on the selected rank";
  }
  if (phase == InjectionPhase::metadata_completion) {
    if (message == InjectionMessage::mismatch) {
      return "compact Buffer Halo request-ID completion injection rank "
             "differs across ranks";
    }
    if (message == InjectionMessage::invalid) {
      return "compact Buffer Halo request-ID completion injection rank is "
             "invalid";
    }
    return "compact Buffer Halo request-ID completion injection requires at "
           "least two request descriptors on the selected rank";
  }
  if (phase == InjectionPhase::payload_post) {
    if (message == InjectionMessage::mismatch) {
      return "compact Buffer Halo payload post injection rank differs across "
             "ranks";
    }
    if (message == InjectionMessage::invalid) {
      return "compact Buffer Halo payload post injection rank is invalid";
    }
    return "compact Buffer Halo payload post injection requires at least two "
           "request descriptors on the selected rank";
  }
  if (message == InjectionMessage::mismatch) {
    return "compact Buffer Halo payload completion injection rank differs "
           "across ranks";
  }
  if (message == InjectionMessage::invalid) {
    return "compact Buffer Halo payload completion injection rank is invalid";
  }
  return "compact Buffer Halo payload completion injection requires at least "
         "two request descriptors on the selected rank";
}

InjectionSelection validate_injection_selection(
    const runtime::MpiContext& context, int configured_rank,
    std::size_t local_descriptor_count, InjectionPhase phase) {
  const bool disabled = configured_rank == -1;
  const bool all_eligible =
      configured_rank == detail::kInjectAllEligibleRanks;
  const bool legal = disabled || all_eligible ||
                     (configured_rank >= 0 &&
                      configured_rank < context.size());
  if (!legal) {
    throw runtime::Error(
        std::string(injection_message(phase, InjectionMessage::invalid)));
  }
  if (disabled) {
    return {};
  }

  const bool local_eligible = local_descriptor_count >= 2U;
  int local_selected = 0;
  if (all_eligible) {
    local_selected = local_eligible ? 1 : 0;
  } else if (context.rank() == configured_rank && local_eligible) {
    local_selected = 1;
  }
  int selected_count = 0;
  if (detail::current_vector_halo_test_options().observe) {
    detail::mutable_vector_halo_test_snapshot()
        .injection_selection_collectives +=
        detail::injection_selection_collective_count(configured_rank);
  }
  runtime::detail::check_mpi(
      MPI_Allreduce(&local_selected, &selected_count, 1, MPI_INT, MPI_SUM,
                    context.comm()),
      "MPI_Allreduce compact Halo injection eligible count");
  const bool enough = all_eligible ? selected_count >= 2 : selected_count == 1;
  if (!enough) {
    throw runtime::Error(std::string(
        injection_message(phase, InjectionMessage::insufficient)));
  }
  return InjectionSelection{local_selected != 0};
}

void terminate_if_uncertain_post_output(
    const runtime::MpiContext& context, bool local_uncertain) noexcept {
  const int local = local_uncertain ? 1 : 0;
  int any_uncertain = 0;
  if (MPI_Allreduce(&local, &any_uncertain, 1, MPI_INT, MPI_MAX,
                    context.comm()) != MPI_SUCCESS ||
      any_uncertain != 0) {
    std::terminate();
  }
}

std::string format_failure(FailureRecord record) {
  std::string message = "compact Buffer Halo ";
  message += record.category == 1 ? "post" : "completion";
  message += " failure: rank=" + std::to_string(record.rank);
  message += " operation=";
  switch (record.operation) {
    case 1:
      message += "MPI_Irecv";
      break;
    case 2:
      message += "MPI_Isend";
      break;
    case 3:
      message += "MPI_Wait";
      break;
    case 4:
      message += "MPI_Irecv_request_ID";
      break;
    case 5:
      message += "MPI_Isend_request_ID";
      break;
    case 6:
      message += "MPI_Wait_request_ID";
      break;
    default:
      message += "unknown_MPI_operation";
      break;
  }
  message += " result=" + std::to_string(record.result);
  message += " peer=" + std::to_string(record.peer);
  message += " chunk_offset=" + std::to_string(record.chunk_offset);
  message += " chunk_count=" + std::to_string(record.chunk_count);
  message += " tag=" + std::to_string(record.tag);
  return message;
}

void observe_failure(FailureRecord record) noexcept {
  if (!detail::current_vector_halo_test_options().observe) {
    return;
  }
  auto& destination = detail::mutable_vector_halo_test_snapshot().failure;
  destination.valid = failed(record);
  destination.category = static_cast<int>(record.category);
  destination.rank = static_cast<int>(record.rank);
  destination.operation = static_cast<int>(record.operation);
  destination.result = static_cast<int>(record.result);
  destination.peer = static_cast<int>(record.peer);
  destination.value_offset =
      record.chunk_offset < 0
          ? 0U
          : static_cast<std::size_t>(record.chunk_offset);
  destination.value_count = static_cast<int>(record.chunk_count);
  destination.tag = static_cast<int>(record.tag);
}

std::size_t non_null_request_count(
    const std::vector<MPI_Request>& requests) noexcept {
  return static_cast<std::size_t>(std::count_if(
      requests.begin(), requests.end(), [](MPI_Request request) {
        return request != MPI_REQUEST_NULL;
      }));
}

bool requests_are_null(const std::vector<MPI_Request>& requests) noexcept {
  return non_null_request_count(requests) == 0U;
}

bool drain_requests_noexcept(std::vector<MPI_Request>& requests,
                             bool cancel) noexcept {
  bool success = true;
  if (cancel) {
    for (MPI_Request& request : requests) {
      if (request != MPI_REQUEST_NULL && MPI_Cancel(&request) != MPI_SUCCESS) {
        success = false;
      }
    }
  }
  for (MPI_Request& request : requests) {
    if (request != MPI_REQUEST_NULL &&
        MPI_Wait(&request, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
      success = false;
    }
  }
  return success && requests_are_null(requests);
}

struct RequestDescriptor {
  bool receive{};
  std::size_t peer_index{};
  int peer{};
  std::size_t value_offset{};
  int value_count{};
  int tag{};
};

FailureRecord failure_for(const RequestDescriptor& descriptor,
                          std::int64_t category, int result) noexcept {
  return FailureRecord{category,
                       -1,
                       descriptor.receive ? 1 : 2,
                       result,
                       descriptor.peer,
                       static_cast<std::int64_t>(descriptor.value_offset),
                       descriptor.value_count,
                       descriptor.tag};
}

struct MetadataRequestDescriptor {
  bool receive{};
  int peer{};
  std::size_t storage_offset{};
  std::size_t value_offset{};
  int value_count{};
  int tag{kMetadataTag};
};

FailureRecord metadata_failure_for(
    const MetadataRequestDescriptor& descriptor, std::int64_t category,
    int result, bool completion) noexcept {
  return FailureRecord{category,
                       -1,
                       completion ? 6 : descriptor.receive ? 4 : 5,
                       result,
                       descriptor.peer,
                       static_cast<std::int64_t>(descriptor.value_offset),
                       descriptor.value_count,
                       descriptor.tag};
}

void quarantine_context_or_terminate(runtime::MpiContext& context,
                                     bool metadata) noexcept {
  try {
    runtime::MpiContext replacement =
        runtime::MpiContext::duplicate(context.comm());
    context = std::move(replacement);
    if (detail::current_vector_halo_test_options().observe && metadata) {
      ++detail::mutable_vector_halo_test_snapshot()
            .metadata_context_replacements;
    }
  } catch (...) {
    std::terminate();
  }
}

void validate_create_inputs(
    const runtime::MpiContext& context,
    const runtime::StructuredDecomposition& decomposition,
    const mesh::MeshTopology& topology,
    execution::ExecutionContext& execution_context) {
  bool local_ok = true;
  std::string local_message;
  try {
    if (decomposition.comm() == MPI_COMM_NULL) {
      throw runtime::Error(
          "compact Buffer Halo requires a live decomposition");
    }
    if (!same(topology.global_extent(), decomposition.global_extent()) ||
        !same(topology.owned_global_box(), decomposition.owned_box())) {
      throw runtime::Error(
          "compact Buffer Halo topology and decomposition differ");
    }
    if (topology.local_cell_count() !=
        checked_add(topology.owned_cell_count(), topology.ghost_cell_count(),
                    "compact topology local count")) {
      throw runtime::Error("compact topology cell counts disagree");
    }
    for (std::size_t local = 0; local < topology.local_cell_count(); ++local) {
      const mesh::EntityOwnership expected =
          local < topology.owned_cell_count()
              ? mesh::EntityOwnership::owned
              : mesh::EntityOwnership::ghost;
      if (topology.cell_ownership(local) != expected) {
        throw runtime::Error(
            "compact topology does not order owned cells before ghosts");
      }
    }
    detail::DevicePathCapabilities capabilities{};
    capabilities.space = execution_context.space();
    capabilities.context_host_access =
        execution_context.supports(execution::ExecutionCapability::host_access);
    capabilities.buffer_host_access = capabilities.context_host_access;
    capabilities.host_mpi = true;
    if (execution_context.backend_identity() == 0U) {
      throw runtime::Error(
          "compact Buffer Halo backend identity must be nonzero");
    }
    if (!execution_context.supports(
            execution::ExecutionCapability::buffer_allocation)) {
      throw runtime::Error(
          "compact Buffer Halo context cannot allocate a Buffer");
    }
    if (detail::select_buffer_halo_path(capabilities) !=
        BufferHaloPath::host_direct) {
      throw runtime::Error(
          "Task 9 production compact Halo requires host-direct");
    }
  } catch (const std::exception& error) {
    local_ok = false;
    local_message = error.what();
  } catch (...) {
    local_ok = false;
    local_message = "unknown compact Buffer Halo input failure";
  }
  require_collective(context, local_ok, local_message);

  const auto options = detail::current_vector_halo_test_options();
  const std::array<std::uint64_t, 6> local_contract{
      execution_context.backend_identity(),
      static_cast<std::uint64_t>(execution_context.space()),
      static_cast<std::uint64_t>(BufferHaloPath::host_direct),
      static_cast<std::uint64_t>(
          static_cast<std::int64_t>(options.inject_metadata_post_failure_rank)),
      static_cast<std::uint64_t>(static_cast<std::int64_t>(
          options.inject_metadata_completion_failure_rank)),
      options.observe ? 1U : 0U};
  std::array<std::uint64_t, 6> minimum{};
  std::array<std::uint64_t, 6> maximum{};
  runtime::detail::check_mpi(
      MPI_Allreduce(local_contract.data(), minimum.data(),
                    static_cast<int>(local_contract.size()), MPI_UINT64_T,
                    MPI_MIN, context.comm()),
      "MPI_Allreduce compact Halo contract minimum");
  runtime::detail::check_mpi(
      MPI_Allreduce(local_contract.data(), maximum.data(),
                    static_cast<int>(local_contract.size()), MPI_UINT64_T,
                    MPI_MAX, context.comm()),
      "MPI_Allreduce compact Halo contract maximum");
  if (minimum[3] != maximum[3]) {
    throw runtime::Error(std::string(injection_message(
        InjectionPhase::metadata_post, InjectionMessage::mismatch)));
  }
  if (minimum[4] != maximum[4]) {
    throw runtime::Error(std::string(injection_message(
        InjectionPhase::metadata_completion, InjectionMessage::mismatch)));
  }
  if (minimum[5] != maximum[5]) {
    throw runtime::Error(
        "compact Buffer Halo observation option differs across ranks");
  }
  require_collective(context,
                     std::equal(minimum.begin(), minimum.begin() + 3,
                                maximum.begin()),
                     "compact Buffer Halo context or path differs across ranks");

  const auto local_chunk = static_cast<std::uint64_t>(options.chunk_limit);
  std::uint64_t minimum_chunk = 0U;
  std::uint64_t maximum_chunk = 0U;
  runtime::detail::check_mpi(
      MPI_Allreduce(&local_chunk, &minimum_chunk, 1, MPI_UINT64_T, MPI_MIN,
                    context.comm()),
      "MPI_Allreduce compact Halo chunk minimum");
  runtime::detail::check_mpi(
      MPI_Allreduce(&local_chunk, &maximum_chunk, 1, MPI_UINT64_T, MPI_MAX,
                    context.comm()),
      "MPI_Allreduce compact Halo chunk maximum");
  bool chunk_ok = minimum_chunk == maximum_chunk;
  std::string chunk_message =
      "compact Buffer Halo chunk limit differs across ranks";
  if (chunk_ok) {
    try {
      static_cast<void>(detail::split_count_ranges(0U, options.chunk_limit));
    } catch (const std::exception& error) {
      chunk_ok = false;
      chunk_message = error.what();
    }
  }
  require_collective(context, chunk_ok, chunk_message);

  local_ok = true;
  local_message.clear();
  try {
    void* attribute = nullptr;
    int present = 0;
    runtime::detail::check_mpi(
        MPI_Comm_get_attr(context.comm(), MPI_TAG_UB, &attribute, &present),
        "MPI_Comm_get_attr MPI_TAG_UB");
    static_cast<void>(detail::effective_tag_upper_bound(
        present != 0, static_cast<const int*>(attribute), kPayloadTag));
    MPI_Errhandler handler = MPI_ERRHANDLER_NULL;
    runtime::detail::check_mpi(
        MPI_Comm_get_errhandler(context.comm(), &handler),
        "MPI_Comm_get_errhandler compact Halo");
    const bool returns_errors = handler == MPI_ERRORS_RETURN;
    if (handler != MPI_ERRHANDLER_NULL) {
      runtime::detail::check_mpi(MPI_Errhandler_free(&handler),
                                 "MPI_Errhandler_free compact Halo");
    }
    if (!returns_errors) {
      throw runtime::Error(
          "compact Buffer Halo communicator lacks MPI_ERRORS_RETURN");
    }
  } catch (const std::exception& error) {
    local_ok = false;
    local_message = error.what();
  } catch (...) {
    local_ok = false;
    local_message = "unknown compact Buffer Halo communicator failure";
  }
  require_collective(context, local_ok, local_message);
}

BuiltPlan build_plan(runtime::MpiContext& context,
                     const runtime::StructuredDecomposition& decomposition,
                     const mesh::MeshTopology& topology) {
  const auto options = detail::current_vector_halo_test_options();
  bool local_ok = true;
  std::string local_message;
  std::vector<detail::OrderedRequest> requests;
  try {
    requests.reserve(topology.ghost_cell_count());
    const runtime::Int3 global = decomposition.global_extent();
    const runtime::Int3 grid = decomposition.process_grid();
    for (std::size_t local = topology.owned_cell_count();
         local < topology.local_cell_count(); ++local) {
      const runtime::Int3 cell = topology.global_cell(local);
      std::array<int, 3> owner_coordinates{
          detail::owner_coordinate(global.x, grid.x, cell.x),
          detail::owner_coordinate(global.y, grid.y, cell.y),
          detail::owner_coordinate(global.z, grid.z, cell.z)};
      int owner_rank = MPI_PROC_NULL;
      runtime::detail::check_mpi(
          MPI_Cart_rank(decomposition.comm(), owner_coordinates.data(),
                        &owner_rank),
          "MPI_Cart_rank compact ghost owner");
      if (owner_rank == MPI_PROC_NULL || owner_rank == context.rank()) {
        throw runtime::Error(
            "compact ghost does not have exactly one remote owner rank");
      }
      requests.push_back(detail::OrderedRequest{
          owner_rank, topology.global_cell_id(local), local});
    }
    requests = detail::order_requests(std::move(requests));
    if (context.rank() == options.inject_request_id_mismatch_rank) {
      if (requests.empty()) {
        throw runtime::Error(
            "request-ID mismatch seam selected a rank with no request");
      }
      requests.front().global_id = topology.global_cell_count();
    }
  } catch (const std::exception& error) {
    local_ok = false;
    local_message = error.what();
  } catch (...) {
    local_ok = false;
    local_message = "unknown compact ghost owner-plan failure";
  }
  require_collective(context, local_ok, local_message);

  std::vector<std::uint64_t> outgoing_counts;
  std::vector<std::uint64_t> incoming_counts;
  std::vector<std::size_t> outgoing_offsets;
  std::vector<std::size_t> incoming_offsets;
  std::vector<mesh::GlobalCellId> outgoing_ids;
  std::vector<mesh::GlobalCellId> incoming_ids;
  local_ok = true;
  local_message.clear();
  try {
    const std::size_t ranks = static_cast<std::size_t>(context.size());
    outgoing_counts.assign(ranks, 0U);
    incoming_counts.assign(ranks, 0U);
    outgoing_offsets.assign(ranks + 1U, 0U);
    incoming_offsets.assign(ranks + 1U, 0U);
    for (const detail::OrderedRequest& request : requests) {
      auto& count = outgoing_counts[static_cast<std::size_t>(request.peer)];
      if (count == std::numeric_limits<std::uint64_t>::max()) {
        throw runtime::Error("compact request count overflow");
      }
      ++count;
    }
  } catch (const std::exception& error) {
    local_ok = false;
    local_message = error.what();
  } catch (...) {
    local_ok = false;
    local_message = "unknown compact request-count allocation failure";
  }
  require_collective(context, local_ok, local_message);

  runtime::detail::check_mpi(
      MPI_Alltoall(outgoing_counts.data(), 1, MPI_UINT64_T,
                   incoming_counts.data(), 1, MPI_UINT64_T, context.comm()),
      "MPI_Alltoall compact Halo request counts");

  local_ok = true;
  local_message.clear();
  try {
    for (std::size_t rank = 0; rank < outgoing_counts.size(); ++rank) {
      outgoing_offsets[rank + 1U] = checked_add(
          outgoing_offsets[rank], checked_size(outgoing_counts[rank],
                                               "outgoing compact count"),
          "outgoing compact offsets");
      incoming_offsets[rank + 1U] = checked_add(
          incoming_offsets[rank], checked_size(incoming_counts[rank],
                                               "incoming compact count"),
          "incoming compact offsets");
      if (rank == static_cast<std::size_t>(context.rank()) &&
          (outgoing_counts[rank] != 0U || incoming_counts[rank] != 0U)) {
        throw runtime::Error("compact Halo request unexpectedly targets self");
      }
    }
    outgoing_ids.resize(outgoing_offsets.back());
    incoming_ids.resize(incoming_offsets.back());
    std::vector<std::size_t> cursor = outgoing_offsets;
    for (const detail::OrderedRequest& request : requests) {
      outgoing_ids[cursor[static_cast<std::size_t>(request.peer)]++] =
          request.global_id;
    }
  } catch (const std::exception& error) {
    local_ok = false;
    local_message = error.what();
  } catch (...) {
    local_ok = false;
    local_message = "unknown compact request-ID storage failure";
  }
  require_collective(context, local_ok, local_message);

  std::vector<MPI_Request> traffic;
  std::vector<MetadataRequestDescriptor> metadata_descriptors;
  local_ok = true;
  local_message.clear();
  try {
    std::size_t request_count = 0U;
    for (std::size_t rank = 0; rank < outgoing_counts.size(); ++rank) {
      request_count = checked_add(
          request_count,
          detail::split_count_ranges(
              incoming_offsets[rank + 1U] - incoming_offsets[rank],
              options.chunk_limit)
              .size(),
          "compact request-ID message count");
      request_count = checked_add(
          request_count,
          detail::split_count_ranges(
              outgoing_offsets[rank + 1U] - outgoing_offsets[rank],
              options.chunk_limit)
              .size(),
          "compact request-ID message count");
    }
    traffic.assign(request_count, MPI_REQUEST_NULL);
    metadata_descriptors.reserve(request_count);
    for (std::size_t rank = 0; rank < incoming_counts.size(); ++rank) {
      const auto chunks = detail::split_count_ranges(
          incoming_offsets[rank + 1U] - incoming_offsets[rank],
          options.chunk_limit);
      for (const detail::CountRange chunk : chunks) {
        metadata_descriptors.push_back(MetadataRequestDescriptor{
            true, static_cast<int>(rank), incoming_offsets[rank] + chunk.offset,
            chunk.offset, chunk.count, kMetadataTag});
      }
    }
    for (std::size_t rank = 0; rank < outgoing_counts.size(); ++rank) {
      const auto chunks = detail::split_count_ranges(
          outgoing_offsets[rank + 1U] - outgoing_offsets[rank],
          options.chunk_limit);
      for (const detail::CountRange chunk : chunks) {
        metadata_descriptors.push_back(MetadataRequestDescriptor{
            false, static_cast<int>(rank),
            outgoing_offsets[rank] + chunk.offset, chunk.offset, chunk.count,
            kMetadataTag});
      }
    }
    if (metadata_descriptors.size() != traffic.size()) {
      throw runtime::Error(
          "compact request-ID descriptor count disagrees with slots");
    }
  } catch (const std::exception& error) {
    local_ok = false;
    local_message = error.what();
  } catch (...) {
    local_ok = false;
    local_message = "unknown compact request-ID traffic allocation failure";
  }
  require_collective(context, local_ok, local_message);

  const InjectionSelection metadata_post_injection =
      validate_injection_selection(
          context, options.inject_metadata_post_failure_rank,
          metadata_descriptors.size(), InjectionPhase::metadata_post);
  const InjectionSelection metadata_completion_injection =
      validate_injection_selection(
          context, options.inject_metadata_completion_failure_rank,
          metadata_descriptors.size(), InjectionPhase::metadata_completion);

  FailureRecord local_traffic_failure{};
  bool uncertain_post_output = false;
  std::size_t successful_posts = 0U;
  for (std::size_t slot = 0U; slot < metadata_descriptors.size(); ++slot) {
    const MetadataRequestDescriptor& descriptor = metadata_descriptors[slot];
    if (metadata_post_injection.local_selected && successful_posts == 1U) {
      local_traffic_failure = metadata_failure_for(
          descriptor, 1, MPI_ERR_OTHER, false);
      uncertain_post_output =
          detail::nonblocking_post_issue_action(
              detail::NonblockingPostIssueOrigin::synthetic_before_call) ==
          detail::NonblockingPostIssueAction::terminate_process;
      if (options.observe) {
        detail::mutable_vector_halo_test_snapshot()
            .metadata_posts_before_failure = successful_posts;
      }
      break;
    }
    if (options.observe) {
      ++detail::mutable_vector_halo_test_snapshot().metadata_post_calls;
    }
    const int result = descriptor.receive
                           ? MPI_Irecv(incoming_ids.data() +
                                           descriptor.storage_offset,
                                       descriptor.value_count, MPI_UINT64_T,
                                       descriptor.peer, descriptor.tag,
                                       context.comm(), &traffic[slot])
                           : MPI_Isend(outgoing_ids.data() +
                                           descriptor.storage_offset,
                                       descriptor.value_count, MPI_UINT64_T,
                                       descriptor.peer, descriptor.tag,
                                       context.comm(), &traffic[slot]);
    if (result != MPI_SUCCESS) {
      local_traffic_failure =
          metadata_failure_for(descriptor, 1, result, false);
      uncertain_post_output =
          detail::nonblocking_post_issue_action(
              detail::NonblockingPostIssueOrigin::mpi_call_error) ==
          detail::NonblockingPostIssueAction::terminate_process;
      break;
    }
    ++successful_posts;
  }
  terminate_if_uncertain_post_output(context, uncertain_post_output);
  if (failed(local_traffic_failure)) {
    local_traffic_failure.rank = context.rank();
  }
  FailureRecord global_traffic_failure{};
  if (!converge_failure(context, local_traffic_failure,
                        global_traffic_failure)) {
    std::terminate();
  }
  if (failed(global_traffic_failure)) {
    auto& snapshot = detail::mutable_vector_halo_test_snapshot();
    if (options.observe) {
      snapshot.metadata_non_null_before_cleanup =
          non_null_request_count(traffic);
    }
    if (!drain_requests_noexcept(traffic, true)) {
      std::terminate();
    }
    if (options.observe) {
      snapshot.metadata_non_null_after_cleanup =
          non_null_request_count(traffic);
    }
    observe_failure(global_traffic_failure);
    quarantine_context_or_terminate(context, true);
    throw runtime::Error(format_failure(global_traffic_failure));
  }

  FailureRecord local_completion_failure{};
  std::size_t completed_prefix = 0U;
  for (std::size_t slot = 0U; slot < metadata_descriptors.size(); ++slot) {
    const MetadataRequestDescriptor& descriptor = metadata_descriptors[slot];
    if (metadata_completion_injection.local_selected &&
        completed_prefix == 1U) {
      local_completion_failure = metadata_failure_for(
          descriptor, 2, MPI_ERR_OTHER, true);
      if (options.observe) {
        detail::mutable_vector_halo_test_snapshot()
            .metadata_completion_prefix = completed_prefix;
      }
      break;
    }
    if (options.observe) {
      ++detail::mutable_vector_halo_test_snapshot().metadata_wait_calls;
    }
    const int result = MPI_Wait(&traffic[slot], MPI_STATUS_IGNORE);
    if (result != MPI_SUCCESS) {
      local_completion_failure =
          metadata_failure_for(descriptor, 2, result, true);
      break;
    }
    ++completed_prefix;
  }
  if (failed(local_completion_failure)) {
    local_completion_failure.rank = context.rank();
  }
  global_traffic_failure = FailureRecord{};
  if (!converge_failure(context, local_completion_failure,
                        global_traffic_failure)) {
    std::terminate();
  }
  if (failed(global_traffic_failure)) {
    auto& snapshot = detail::mutable_vector_halo_test_snapshot();
    if (options.observe) {
      snapshot.metadata_non_null_before_cleanup =
          non_null_request_count(traffic);
    }
    if (!drain_requests_noexcept(traffic, true)) {
      std::terminate();
    }
    if (options.observe) {
      snapshot.metadata_non_null_after_cleanup =
          non_null_request_count(traffic);
    }
    observe_failure(global_traffic_failure);
    quarantine_context_or_terminate(context, true);
    throw runtime::Error(format_failure(global_traffic_failure));
  }
  if (!requests_are_null(traffic)) {
    std::terminate();
  }

  std::vector<std::vector<std::size_t>> sends_by_rank;
  local_ok = true;
  local_message.clear();
  try {
    sends_by_rank.resize(static_cast<std::size_t>(context.size()));
    for (std::size_t rank = 0; rank < incoming_counts.size(); ++rank) {
      auto& sends = sends_by_rank[rank];
      sends.reserve(incoming_offsets[rank + 1U] - incoming_offsets[rank]);
      mesh::GlobalCellId previous = 0U;
      bool have_previous = false;
      for (std::size_t offset = incoming_offsets[rank];
           offset < incoming_offsets[rank + 1U]; ++offset) {
        const mesh::GlobalCellId requested = incoming_ids[offset];
        if (have_previous && requested <= previous) {
          throw runtime::Error(
              "compact request-ID order or duplicate disagreement");
        }
        previous = requested;
        have_previous = true;
        const auto local = topology.find_local_cell(requested);
        if (!local.has_value() || *local >= topology.owned_cell_count() ||
            topology.cell_ownership(*local) != mesh::EntityOwnership::owned) {
          throw runtime::Error(
              "compact request ID is missing, non-owned, or wrong-rank");
        }
        sends.push_back(*local);
      }
    }
  } catch (const std::exception& error) {
    local_ok = false;
    local_message = std::string("compact request-ID validation failed: ") +
                    error.what();
  } catch (...) {
    local_ok = false;
    local_message = "unknown compact request-ID validation failure";
  }
  require_collective(context, local_ok, local_message);

  BuiltPlan result;
  local_ok = true;
  local_message.clear();
  try {
    result.receive_value_count = topology.ghost_cell_count();
    std::size_t request_cursor = 0U;
    for (int rank = 0; rank < context.size(); ++rank) {
      PlanPeer peer{};
      peer.rank = rank;
      while (request_cursor < requests.size() &&
             requests[request_cursor].peer < rank) {
        ++request_cursor;
      }
      std::size_t cursor = request_cursor;
      while (cursor < requests.size() && requests[cursor].peer == rank) {
        peer.receive_indices.push_back(requests[cursor].local_index);
        ++cursor;
      }
      request_cursor = cursor;
      peer.send_indices =
          std::move(sends_by_rank[static_cast<std::size_t>(rank)]);
      result.send_value_count = checked_add(
          result.send_value_count, peer.send_indices.size(),
          "compact Halo send value count");
      if (!peer.send_indices.empty() || !peer.receive_indices.empty()) {
        result.peers.push_back(std::move(peer));
      }
    }
  } catch (const std::exception& error) {
    local_ok = false;
    local_message = error.what();
  } catch (...) {
    local_ok = false;
    local_message = "unknown compact plan finalization failure";
  }
  require_collective(context, local_ok, local_message);
  return result;
}

enum class PublicOperation : int { exchange = 1, begin = 2, wait = 3 };

enum class VectorIssue {
  none,
  layout,
  backend,
  residence,
  allocation,
  identity,
  validation
};

std::string_view vector_issue_message(VectorIssue issue,
                                      bool wait_target) noexcept {
  switch (issue) {
    case VectorIssue::none:
      return {};
    case VectorIssue::layout:
      return wait_target ? "compact Buffer Halo wait target layout differs"
                         : "compact Buffer Halo vector layout differs";
    case VectorIssue::backend:
      return wait_target ? "compact Buffer Halo wait target backend differs"
                         : "compact Buffer Halo vector backend differs";
    case VectorIssue::residence:
      return wait_target
                 ? "compact Buffer Halo wait target residence differs"
                 : "compact Buffer Halo vector residence is not host";
    case VectorIssue::allocation:
      return wait_target
                 ? "compact Buffer Halo wait target allocation is not live"
                 : "compact Buffer Halo vector allocation is not live";
    case VectorIssue::identity:
      return "compact Buffer Halo wait target identity or epoch differs";
    case VectorIssue::validation:
      return wait_target ? "compact Buffer Halo wait target validation failed"
                         : "compact Buffer Halo vector validation failed";
  }
  return "unknown compact Buffer Halo vector validation failure";
}

}  // namespace

class GhostedVectorHalo::Impl final {
 public:
  struct PeerState {
    PeerState(execution::ExecutionContext& context, PlanPeer plan,
              std::size_t chunk_limit)
        : rank(plan.rank),
          send_indices(std::move(plan.send_indices)),
          receive_indices(std::move(plan.receive_indices)),
          send_buffer(context,
                      detail::checked_vector_bytes(send_indices.size())),
          receive_buffer(context,
                         detail::checked_vector_bytes(receive_indices.size())),
          send_chunks(
              detail::split_count_ranges(send_indices.size(), chunk_limit)),
          receive_chunks(detail::split_count_ranges(receive_indices.size(),
                                                    chunk_limit)) {}

    PeerState(PeerState&&) noexcept = default;
    PeerState& operator=(PeerState&&) noexcept = default;
    PeerState(const PeerState&) = delete;
    PeerState& operator=(const PeerState&) = delete;

    int rank{};
    std::vector<std::size_t> send_indices;
    std::vector<std::size_t> receive_indices;
    execution::Buffer send_buffer;
    execution::Buffer receive_buffer;
    std::vector<detail::CountRange> send_chunks;
    std::vector<detail::CountRange> receive_chunks;
  };

  Impl(runtime::MpiContext context, VectorLayout layout,
       execution::BackendIdentity backend_identity,
       std::vector<PeerState> peers, std::size_t send_value_count,
       std::size_t receive_value_count, std::vector<MPI_Request> requests,
       std::vector<RequestDescriptor> request_descriptors,
       std::optional<detail::VectorHaloTestSnapshot> observation) noexcept
      : context_(std::move(context)),
        layout_(std::move(layout)),
        backend_identity_(backend_identity),
        peers_(std::move(peers)),
        send_value_count_(send_value_count),
        receive_value_count_(receive_value_count),
        requests_(std::move(requests)),
        request_descriptors_(std::move(request_descriptors)),
        observation_(std::move(observation)) {
    if (observation_.has_value()) {
      detail::activate_vector_halo_test_observation(&*observation_);
    }
  }

  ~Impl() noexcept {
    if (observation_.has_value()) {
      detail::deactivate_vector_halo_test_observation(&*observation_);
    }
    if (!active_) {
      return;
    }
    if (!runtime::detail::mpi_is_active() || !drain_requests_noexcept()) {
      std::terminate();
    }
    clear_active();
  }

  BufferHaloPath path() const noexcept { return BufferHaloPath::host_direct; }
  std::size_t owned_count() const noexcept { return layout_.owned_count(); }
  std::size_t ghost_count() const noexcept { return layout_.ghost_count(); }
  std::size_t send_value_count() const noexcept { return send_value_count_; }
  std::size_t receive_value_count() const noexcept {
    return receive_value_count_;
  }

  void begin(const GhostedVector& vector) {
    require_public_operation(PublicOperation::begin);
    activate_observation_if_requested();
    begin_internal(vector);
  }

  void wait(GhostedVector& vector) {
    require_public_operation(PublicOperation::wait);
    activate_observation_if_requested();
    require_wait_target(vector);
    finish_internal(vector);
  }

  void exchange(GhostedVector& vector) {
    require_public_operation(PublicOperation::exchange);
    activate_observation_if_requested();
    begin_internal(vector);
    finish_internal(vector);
  }

 private:
  struct PostOutcome {
    FailureRecord failure;
    bool uncertain_output{};
  };

  void require_public_operation(PublicOperation operation) const {
    runtime::detail::require_mpi_active("enter compact Buffer Halo operation");
    const auto options = detail::current_vector_halo_test_options();
    const std::array<int, 5> local{
        static_cast<int>(operation), active_ ? 1 : 0,
        options.inject_post_failure_rank,
        options.inject_completion_failure_rank, options.observe ? 1 : 0};
    std::array<int, 5> minimum{};
    std::array<int, 5> maximum{};
    runtime::detail::check_mpi(
        MPI_Allreduce(local.data(), minimum.data(),
                      static_cast<int>(local.size()), MPI_INT, MPI_MIN,
                      context_.comm()),
        "MPI_Allreduce compact Halo operation minimum");
    runtime::detail::check_mpi(
        MPI_Allreduce(local.data(), maximum.data(),
                      static_cast<int>(local.size()), MPI_INT, MPI_MAX,
                      context_.comm()),
        "MPI_Allreduce compact Halo operation maximum");
    if (minimum[2] != maximum[2]) {
      throw runtime::Error(std::string(injection_message(
          InjectionPhase::payload_post, InjectionMessage::mismatch)));
    }
    if (minimum[3] != maximum[3]) {
      throw runtime::Error(std::string(injection_message(
          InjectionPhase::payload_completion, InjectionMessage::mismatch)));
    }
    if (minimum[4] != maximum[4]) {
      throw runtime::Error(
          "compact Buffer Halo observation option differs across ranks");
    }
    const bool contract_agrees = minimum[0] == maximum[0] &&
                                 minimum[1] == maximum[1];
    const bool local_is_reference = local[0] == minimum[0] &&
                                    local[1] == minimum[1];
    require_collective(
        context_, contract_agrees || local_is_reference,
        "compact Buffer Halo public operation or state mismatch");
    if (!contract_agrees) {
      throw runtime::Error(
          "compact Buffer Halo public operation or state mismatch");
    }
    const bool requires_active = operation == PublicOperation::wait;
    const bool state_ok = active_ == requires_active;
    require_collective(
        context_, state_ok,
        requires_active
            ? "compact Buffer Halo wait requires an active operation"
            : "compact Buffer Halo begin/exchange requires an idle operation");
  }

  void activate_observation_if_requested() {
    if (!detail::current_vector_halo_test_options().observe) {
      return;
    }
    if (!observation_.has_value()) {
      throw runtime::Error(
          "compact Buffer Halo observation was not enabled at creation");
    }
    detail::activate_vector_halo_test_observation(&*observation_);
  }

  VectorIssue inspect_vector(const GhostedVector& vector,
                             bool wait_target) const noexcept {
    try {
      if (vector.layout() != layout_) {
        return VectorIssue::layout;
      }
      if (vector.backend_identity() != backend_identity_) {
        return VectorIssue::backend;
      }
      if (vector.space() != execution::ExecutionSpace::host) {
        return VectorIssue::residence;
      }
      if (vector.allocation_identity() == 0U || vector.epoch() == 0U) {
        return VectorIssue::allocation;
      }
      static_cast<void>(vector.local_view().data());
      if (wait_target &&
          (vector.allocation_identity() != pending_allocation_ ||
           vector.epoch() != pending_epoch_)) {
        return VectorIssue::identity;
      }
    } catch (const std::exception&) {
      return VectorIssue::validation;
    } catch (...) {
      return VectorIssue::validation;
    }
    return VectorIssue::none;
  }

  void require_wait_target(const GhostedVector& vector) const {
    const VectorIssue issue = inspect_vector(vector, true);
    require_collective(context_, issue == VectorIssue::none,
                       vector_issue_message(issue, true));
  }

  void record_observation() {
    const auto options = detail::current_vector_halo_test_options();
    if (!options.observe) {
      return;
    }
    if (!observation_.has_value()) {
      throw runtime::Error(
          "compact Buffer Halo observation was not enabled at creation");
    }
    detail::activate_vector_halo_test_observation(&*observation_);
    auto& snapshot = detail::mutable_vector_halo_test_snapshot();
    snapshot.receives_preceded_sends = true;
    snapshot.chunk_offsets_ordered = true;
    snapshot.receive_posts = 0U;
    snapshot.send_posts = 0U;
    snapshot.request_capacity = requests_.capacity();
    snapshot.runtime_posts_before_failure = 0U;
    snapshot.runtime_completion_prefix = 0U;
    snapshot.runtime_non_null_before_cleanup = 0U;
    snapshot.runtime_non_null_after_cleanup = 0U;
    snapshot.failure = detail::FailureDiagnosticSnapshot{};
    if (snapshot.send_wire_identities.size() != peers_.size() ||
        snapshot.receive_wire_identities.size() != peers_.size() ||
        snapshot.post_events.size() != requests_.size()) {
      std::terminate();
    }
    for (std::size_t index = 0U; index < peers_.size(); ++index) {
      snapshot.send_wire_identities[index] =
          peers_[index].send_buffer.allocation_identity();
      snapshot.receive_wire_identities[index] =
          peers_[index].receive_buffer.allocation_identity();
    }
  }

  void pack(const GhostedVector& vector) {
    const auto source = vector.local_view();
    for (PeerState& peer : peers_) {
      auto destination = peer.send_buffer.view(0U, peer.send_indices.size());
      for (std::size_t index = 0; index < peer.send_indices.size(); ++index) {
        const std::size_t local = peer.send_indices[index];
        if (local >= layout_.owned_count()) {
          throw runtime::Error(
              "compact Halo send index is not in the owned prefix");
        }
        destination[index] = source[local];
      }
    }
  }

  void begin_internal(const GhostedVector& vector) {
    const VectorIssue issue = inspect_vector(vector, false);
    require_collective(context_, issue == VectorIssue::none,
                       vector_issue_message(issue, false));

    bool local_prepared = true;
    std::string preparation_message;
    try {
      pack(vector);
      record_observation();
    } catch (const std::exception& error) {
      local_prepared = false;
      preparation_message = error.what();
    } catch (...) {
      local_prepared = false;
      preparation_message = "unknown compact Halo pack failure";
    }
    require_collective(context_, local_prepared, preparation_message);

    const InjectionSelection post_injection = validate_injection_selection(
        context_,
        detail::current_vector_halo_test_options().inject_post_failure_rank,
        request_descriptors_.size(), InjectionPhase::payload_post);

    pending_allocation_ = vector.allocation_identity();
    pending_epoch_ = vector.epoch();
    active_ = true;
    std::fill(requests_.begin(), requests_.end(), MPI_REQUEST_NULL);
    const PostOutcome post_outcome = post_all(post_injection);
    terminate_if_uncertain_post_output(context_,
                                       post_outcome.uncertain_output);
    const FailureRecord local_failure = post_outcome.failure;
    FailureRecord global_failure;
    if (!converge_failure(context_, local_failure, global_failure)) {
      std::terminate();
    }
    if (failed(global_failure)) {
      auto& snapshot = detail::mutable_vector_halo_test_snapshot();
      if (detail::current_vector_halo_test_options().observe) {
        snapshot.runtime_non_null_before_cleanup =
            non_null_request_count(requests_);
      }
      if (!cancel_and_drain_noexcept()) {
        std::terminate();
      }
      if (detail::current_vector_halo_test_options().observe) {
        snapshot.runtime_non_null_after_cleanup =
            non_null_request_count(requests_);
      }
      observe_failure(global_failure);
      replace_context();
      clear_active();
      throw runtime::Error(format_failure(global_failure));
    }
  }

  PostOutcome post_all(InjectionSelection injection) noexcept {
    PostOutcome outcome{};
    auto& snapshot = detail::mutable_vector_halo_test_snapshot();
    const auto options = detail::current_vector_halo_test_options();
    const bool observe = options.observe;
    std::size_t successful_posts = 0U;
    for (std::size_t slot = 0U; slot < request_descriptors_.size(); ++slot) {
      const RequestDescriptor& descriptor = request_descriptors_[slot];
      if (injection.local_selected && successful_posts == 1U) {
        outcome.failure = failure_for(descriptor, 1, MPI_ERR_OTHER);
        outcome.uncertain_output =
            detail::nonblocking_post_issue_action(
                detail::NonblockingPostIssueOrigin::synthetic_before_call) ==
            detail::NonblockingPostIssueAction::terminate_process;
        if (observe) {
          snapshot.runtime_posts_before_failure = successful_posts;
        }
        break;
      }
      PeerState& peer = peers_[descriptor.peer_index];
      const int result =
          descriptor.receive
              ? MPI_Irecv(peer.receive_buffer
                              .view(descriptor.value_offset,
                                    static_cast<std::size_t>(
                                        descriptor.value_count))
                              .data(),
                          descriptor.value_count, MPI_DOUBLE, descriptor.peer,
                          descriptor.tag, context_.comm(), &requests_[slot])
              : MPI_Isend(peer.send_buffer
                              .view(descriptor.value_offset,
                                    static_cast<std::size_t>(
                                        descriptor.value_count))
                              .data(),
                          descriptor.value_count, MPI_DOUBLE, descriptor.peer,
                          descriptor.tag, context_.comm(), &requests_[slot]);
      if (observe) {
        if (descriptor.receive) {
          ++snapshot.receive_posts;
        } else {
          ++snapshot.send_posts;
        }
        snapshot.post_events[slot] = detail::WirePostEvent{
            descriptor.receive, descriptor.peer, descriptor.value_offset,
            descriptor.value_count, descriptor.tag};
      }
      if (result != MPI_SUCCESS) {
        outcome.failure = failure_for(descriptor, 1, result);
        outcome.uncertain_output =
            detail::nonblocking_post_issue_action(
                detail::NonblockingPostIssueOrigin::mpi_call_error) ==
            detail::NonblockingPostIssueAction::terminate_process;
        break;
      }
      ++successful_posts;
    }
    if (observe) {
      bool send_seen = false;
      snapshot.receives_preceded_sends = true;
      for (std::size_t slot = 0U; slot < successful_posts; ++slot) {
        if (!request_descriptors_[slot].receive) {
          send_seen = true;
        } else if (send_seen) {
          snapshot.receives_preceded_sends = false;
        }
      }
    }
    if (failed(outcome.failure)) {
      outcome.failure.rank = context_.rank();
    }
    return outcome;
  }

  FailureRecord complete_all(InjectionSelection injection) noexcept {
    FailureRecord failure{};
    const auto options = detail::current_vector_halo_test_options();
    std::size_t completed_prefix = 0U;
    for (std::size_t slot = 0U; slot < requests_.size(); ++slot) {
      const RequestDescriptor& descriptor = request_descriptors_[slot];
      if (injection.local_selected && completed_prefix == 1U) {
        failure = FailureRecord{
            2,
            context_.rank(),
            3,
            MPI_ERR_OTHER,
            descriptor.peer,
            static_cast<std::int64_t>(descriptor.value_offset),
            descriptor.value_count,
            descriptor.tag};
        if (options.observe) {
          detail::mutable_vector_halo_test_snapshot()
              .runtime_completion_prefix = completed_prefix;
        }
        break;
      }
      const int result = MPI_Wait(&requests_[slot], MPI_STATUS_IGNORE);
      if (result != MPI_SUCCESS) {
        failure = FailureRecord{
            2,
            context_.rank(),
            3,
            result,
            descriptor.peer,
            static_cast<std::int64_t>(descriptor.value_offset),
            descriptor.value_count,
            descriptor.tag};
        break;
      }
      ++completed_prefix;
    }
    return failure;
  }

  void finish_internal(GhostedVector& vector) {
    const InjectionSelection completion_injection =
        validate_injection_selection(
            context_, detail::current_vector_halo_test_options()
                          .inject_completion_failure_rank,
            request_descriptors_.size(), InjectionPhase::payload_completion);
    const FailureRecord local_failure = complete_all(completion_injection);
    FailureRecord global_failure;
    if (!converge_failure(context_, local_failure, global_failure)) {
      std::terminate();
    }
    if (failed(global_failure)) {
      auto& snapshot = detail::mutable_vector_halo_test_snapshot();
      if (detail::current_vector_halo_test_options().observe) {
        snapshot.runtime_non_null_before_cleanup =
            non_null_request_count(requests_);
      }
      if (!cancel_and_drain_noexcept()) {
        std::terminate();
      }
      if (detail::current_vector_halo_test_options().observe) {
        snapshot.runtime_non_null_after_cleanup =
            non_null_request_count(requests_);
      }
      observe_failure(global_failure);
      replace_context();
      clear_active();
      throw runtime::Error(format_failure(global_failure));
    }

    auto target = vector.local_view();
    for (PeerState& peer : peers_) {
      const auto source =
          static_cast<const execution::Buffer&>(peer.receive_buffer)
              .view(0U, peer.receive_indices.size());
      for (std::size_t index = 0; index < peer.receive_indices.size();
           ++index) {
        const std::size_t local = peer.receive_indices[index];
        if (local < layout_.owned_count() || local >= layout_.local_count()) {
          std::terminate();
        }
        target[local] = source[index];
      }
    }
    clear_active();
  }

  bool all_requests_null() const noexcept {
    return std::all_of(requests_.begin(), requests_.end(),
                       [](MPI_Request request) {
                         return request == MPI_REQUEST_NULL;
                       });
  }

  bool drain_requests_noexcept() noexcept {
    bool success = true;
    for (MPI_Request& request : requests_) {
      if (request != MPI_REQUEST_NULL &&
          MPI_Wait(&request, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
        success = false;
      }
    }
    return success && all_requests_null();
  }

  bool cancel_and_drain_noexcept() noexcept {
    bool success = true;
    for (MPI_Request& request : requests_) {
      if (request != MPI_REQUEST_NULL && MPI_Cancel(&request) != MPI_SUCCESS) {
        success = false;
      }
    }
    return drain_requests_noexcept() && success;
  }

  void replace_context() noexcept {
    quarantine_context_or_terminate(context_, false);
    if (detail::current_vector_halo_test_options().observe) {
      ++detail::mutable_vector_halo_test_snapshot().context_replacements;
    }
  }

  void clear_active() noexcept {
    active_ = false;
    pending_allocation_ = 0U;
    pending_epoch_ = 0U;
  }

  runtime::MpiContext context_;
  VectorLayout layout_;
  execution::BackendIdentity backend_identity_{};
  std::vector<PeerState> peers_;
  std::size_t send_value_count_{};
  std::size_t receive_value_count_{};
  std::vector<MPI_Request> requests_;
  const std::vector<RequestDescriptor> request_descriptors_;
  std::optional<detail::VectorHaloTestSnapshot> observation_;
  bool active_{};
  execution::AllocationIdentity pending_allocation_{};
  std::uint64_t pending_epoch_{};
};

GhostedVectorHalo GhostedVectorHalo::create(
    const runtime::StructuredDecomposition& decomposition,
    const mesh::MeshTopology& topology,
    execution::ExecutionContext& execution_context) {
  runtime::detail::require_mpi_active("create compact Buffer Halo");
  const bool observe =
      detail::current_vector_halo_test_options().observe;
  if (observe) {
    detail::begin_vector_halo_creation_observation();
  }
  runtime::MpiContext context =
      runtime::MpiContext::duplicate(decomposition.comm());
  validate_create_inputs(context, decomposition, topology, execution_context);
  BuiltPlan plan = build_plan(context, decomposition, topology);
  VectorLayout layout;
  std::vector<Impl::PeerState> peers;
  std::vector<MPI_Request> requests;
  std::vector<RequestDescriptor> request_descriptors;
  bool local_ok = true;
  std::string local_message;
  try {
    layout = VectorLayout::from_topology(topology);
    peers.reserve(plan.peers.size());
    std::size_t request_count = 0U;
    const std::size_t chunk_limit =
        detail::current_vector_halo_test_options().chunk_limit;
    for (PlanPeer& peer : plan.peers) {
      request_count = checked_add(
          request_count,
          detail::split_count_ranges(peer.send_indices.size(), chunk_limit)
              .size(),
          "compact Halo request slots");
      request_count = checked_add(
          request_count,
          detail::split_count_ranges(peer.receive_indices.size(), chunk_limit)
              .size(),
          "compact Halo request slots");
      peers.emplace_back(execution_context, std::move(peer), chunk_limit);
    }
    requests.assign(request_count, MPI_REQUEST_NULL);
    request_descriptors.reserve(request_count);
    for (std::size_t peer_index = 0U; peer_index < peers.size();
         ++peer_index) {
      const Impl::PeerState& peer = peers[peer_index];
      for (const detail::CountRange chunk : peer.receive_chunks) {
        request_descriptors.push_back(RequestDescriptor{
            true, peer_index, peer.rank, chunk.offset, chunk.count,
            kPayloadTag});
      }
    }
    for (std::size_t peer_index = 0U; peer_index < peers.size();
         ++peer_index) {
      const Impl::PeerState& peer = peers[peer_index];
      for (const detail::CountRange chunk : peer.send_chunks) {
        request_descriptors.push_back(RequestDescriptor{
            false, peer_index, peer.rank, chunk.offset, chunk.count,
            kPayloadTag});
      }
    }
    if (request_descriptors.size() != requests.size()) {
      throw runtime::Error(
          "compact Halo request descriptors disagree with slots");
    }
    if (observe) {
      detail::prepare_vector_halo_creation_observation(peers.size(),
                                                       requests.size());
    }
  } catch (const std::exception& error) {
    local_ok = false;
    local_message = error.what();
  } catch (...) {
    local_ok = false;
    local_message = "unknown compact Halo workspace allocation failure";
  }
  require_collective(context, local_ok, local_message);

  void* storage = nullptr;
  local_ok = true;
  local_message.clear();
  try {
    storage = ::operator new(sizeof(Impl));
  } catch (const std::exception& error) {
    local_ok = false;
    local_message = error.what();
  } catch (...) {
    local_ok = false;
    local_message = "unknown compact Halo object allocation failure";
  }
  try {
    require_collective(context, local_ok, local_message);
  } catch (...) {
    ::operator delete(storage);
    throw;
  }
  std::optional<detail::VectorHaloTestSnapshot> observation;
  if (observe) {
    observation.emplace(
        detail::take_vector_halo_creation_observation());
  }
  auto* implementation = ::new (storage) Impl(
      std::move(context), std::move(layout),
      execution_context.backend_identity(), std::move(peers),
      plan.send_value_count, plan.receive_value_count, std::move(requests),
      std::move(request_descriptors), std::move(observation));
  return GhostedVectorHalo(std::unique_ptr<Impl>(implementation));
}

GhostedVectorHalo::GhostedVectorHalo(
    std::unique_ptr<Impl> implementation) noexcept
    : implementation_(std::move(implementation)) {}

GhostedVectorHalo::~GhostedVectorHalo() noexcept = default;
GhostedVectorHalo::GhostedVectorHalo(GhostedVectorHalo&&) noexcept = default;

BufferHaloPath GhostedVectorHalo::path() const {
  if (!implementation_) {
    throw runtime::Error("moved-from compact Buffer Halo has no state");
  }
  return implementation_->path();
}

std::size_t GhostedVectorHalo::owned_count() const {
  if (!implementation_) {
    throw runtime::Error("moved-from compact Buffer Halo has no state");
  }
  return implementation_->owned_count();
}

std::size_t GhostedVectorHalo::ghost_count() const {
  if (!implementation_) {
    throw runtime::Error("moved-from compact Buffer Halo has no state");
  }
  return implementation_->ghost_count();
}

std::size_t GhostedVectorHalo::send_value_count() const {
  if (!implementation_) {
    throw runtime::Error("moved-from compact Buffer Halo has no state");
  }
  return implementation_->send_value_count();
}

std::size_t GhostedVectorHalo::receive_value_count() const {
  if (!implementation_) {
    throw runtime::Error("moved-from compact Buffer Halo has no state");
  }
  return implementation_->receive_value_count();
}

void GhostedVectorHalo::exchange(GhostedVector& vector) {
  if (!implementation_) {
    throw runtime::Error("moved-from compact Buffer Halo has no state");
  }
  implementation_->exchange(vector);
}

void GhostedVectorHalo::begin(const GhostedVector& vector) {
  if (!implementation_) {
    throw runtime::Error("moved-from compact Buffer Halo has no state");
  }
  implementation_->begin(vector);
}

void GhostedVectorHalo::wait(GhostedVector& vector) {
  if (!implementation_) {
    throw runtime::Error("moved-from compact Buffer Halo has no state");
  }
  implementation_->wait(vector);
}

}  // namespace hundun::linear

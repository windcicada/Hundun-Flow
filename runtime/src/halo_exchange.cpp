// SPDX-License-Identifier: Apache-2.0

#include "hundun/runtime/halo_exchange.hpp"

#include "hundun/runtime/collective_status.hpp"
#include "hundun/runtime/error.hpp"
#include "hundun/runtime/field_storage.hpp"
#include "hundun/runtime/mpi_context.hpp"
#include "hundun/runtime/structured_decomposition.hpp"
#include "halo_detail.hpp"
#include "halo_test_access.hpp"
#include "mpi_error.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <new>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace hundun::runtime {
static_assert(std::is_trivially_copyable_v<double>);
static_assert(std::is_trivially_copyable_v<std::int32_t>);
static_assert(std::is_trivially_copyable_v<std::uint8_t>);
static_assert(std::is_trivially_copyable_v<detail::HaloFailureRecord>);
static_assert(sizeof(detail::HaloFailureRecord) ==
              8U * sizeof(std::int64_t));

namespace detail {
namespace {

thread_local HaloTestOptions test_options;
thread_local HaloTestSnapshot test_snapshot;

}  // namespace

void set_halo_test_options(HaloTestOptions options) {
  test_options = options;
}

void reset_halo_test_observation() noexcept {
  test_snapshot = HaloTestSnapshot{};
}

HaloTestSnapshot halo_test_snapshot() noexcept { return test_snapshot; }

HaloTestOptions current_halo_test_options() noexcept { return test_options; }

HaloTestSnapshot& mutable_halo_test_snapshot() noexcept {
  return test_snapshot;
}

}  // namespace detail
namespace {

bool same(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool same(Box3 left, Box3 right) noexcept {
  return same(left.begin, right.begin) && same(left.end, right.end);
}

bool same(const ExchangeRegion& left, const ExchangeRegion& right) noexcept {
  return same(left.offset, right.offset) &&
         left.neighbor_rank == right.neighbor_rank &&
         same(left.send_box, right.send_box) &&
         same(left.receive_box, right.receive_box);
}

std::size_t scalar_size(ScalarType scalar_type) noexcept {
  switch (scalar_type) {
    case ScalarType::float64:
      return sizeof(double);
    case ScalarType::int32:
      return sizeof(std::int32_t);
    case ScalarType::uint8:
      return sizeof(std::uint8_t);
  }
  return 0U;
}

bool checked_add(std::size_t left, std::size_t right,
                 std::size_t& result) noexcept {
  constexpr std::size_t limit = std::numeric_limits<std::size_t>::max();
  if (left > limit - right) {
    return false;
  }
  result = left + right;
  return true;
}

bool checked_multiply(std::size_t left, std::size_t right,
                      std::size_t& result) noexcept {
  constexpr std::size_t limit = std::numeric_limits<std::size_t>::max();
  if (right != 0U && left > limit / right) {
    return false;
  }
  result = left * right;
  return true;
}

enum class HaloError : int {
  none = 0,
  moved_from = 1,
  begin_while_active = 2,
  wait_while_idle = 3,
  invalid_field = 4,
  incompatible_local_layout = 5,
  wire_layout_mismatch = 6,
  wait_target_mismatch = 7,
  preparation_failure = 8,
  agreement_collective_failure = 9
};

const char* halo_error_text(HaloError error) noexcept {
  switch (error) {
    case HaloError::none:
      return "";
    case HaloError::moved_from:
      return "halo exchange object has been moved from";
    case HaloError::begin_while_active:
      return "halo begin requires an idle exchange object";
    case HaloError::wait_while_idle:
      return "halo wait requires an active exchange object";
    case HaloError::invalid_field:
      return "halo field ID or storage is invalid";
    case HaloError::incompatible_local_layout:
      return "halo field layout is incompatible with the exchange plan";
    case HaloError::wire_layout_mismatch:
      return "halo wire layout differs across communicator ranks";
    case HaloError::wait_target_mismatch:
      return "halo wait target differs from the pending field";
    case HaloError::preparation_failure:
      return "halo buffer preparation failed before communication";
    case HaloError::agreement_collective_failure:
      return "halo agreement collective failed before communication";
  }
  return "halo exchange failed";
}

[[noreturn]] void throw_halo_error(HaloError error) {
  throw Error(halo_error_text(error));
}

HaloError converge_fixed_error(const MpiContext& context,
                               HaloError local_error) {
  const int size = context.size();
  const int rank = context.rank();
  const int local_failing_rank =
      local_error == HaloError::none ? size : rank;
  int failing_rank = size;
  detail::check_mpi(MPI_Allreduce(&local_failing_rank, &failing_rank, 1,
                                  MPI_INT, MPI_MIN, context.comm()),
                    "MPI_Allreduce");
  if (failing_rank == size) {
    return HaloError::none;
  }
  int error_code = rank == failing_rank ? static_cast<int>(local_error) : 0;
  detail::check_mpi(
      MPI_Bcast(&error_code, 1, MPI_INT, failing_rank, context.comm()),
      "MPI_Bcast");
  return static_cast<HaloError>(error_code);
}

bool all_requests_null(const std::vector<MPI_Request>& requests) noexcept {
  return std::all_of(requests.begin(), requests.end(),
                     [](MPI_Request request) {
                       return request == MPI_REQUEST_NULL;
                     });
}

bool has_failure(detail::HaloFailureRecord record) noexcept {
  return record.category !=
         static_cast<std::int64_t>(detail::HaloFailureCategory::none);
}

std::int64_t diagnostic_offset(std::size_t value) noexcept {
  const auto maximum =
      static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max());
  return value > maximum ? std::numeric_limits<std::int64_t>::max()
                         : static_cast<std::int64_t>(value);
}

detail::HaloFailureRecord make_failure_record(
    detail::HaloFailureCategory category, int rank,
    detail::HaloMpiOperation operation, int mpi_result, int region_index,
    std::size_t chunk_offset, int chunk_count, int tag) noexcept {
  return detail::HaloFailureRecord{
      static_cast<std::int64_t>(category),
      static_cast<std::int64_t>(rank),
      static_cast<std::int64_t>(operation),
      static_cast<std::int64_t>(mpi_result),
      static_cast<std::int64_t>(region_index),
      diagnostic_offset(chunk_offset),
      static_cast<std::int64_t>(chunk_count),
      static_cast<std::int64_t>(tag)};
}

bool converge_failure_record(const MpiContext& context,
                             detail::HaloFailureRecord local,
                             detail::HaloFailureRecord& global) noexcept {
  const int local_failing_rank =
      has_failure(local) ? context.rank() : context.size();
  int failing_rank = context.size();
  if (MPI_Allreduce(&local_failing_rank, &failing_rank, 1, MPI_INT, MPI_MIN,
                    context.comm()) != MPI_SUCCESS) {
    return false;
  }
  if (failing_rank == context.size()) {
    global = detail::HaloFailureRecord{};
    return true;
  }
  global = context.rank() == failing_rank ? local
                                           : detail::HaloFailureRecord{};
  if (MPI_Bcast(&global, static_cast<int>(sizeof(global)), MPI_BYTE,
                failing_rank, context.comm()) != MPI_SUCCESS) {
    return false;
  }
  global.failing_rank = static_cast<std::int64_t>(failing_rank);
  return true;
}

const char* operation_text(std::int64_t operation) noexcept {
  switch (static_cast<detail::HaloMpiOperation>(operation)) {
    case detail::HaloMpiOperation::none:
      return "none";
    case detail::HaloMpiOperation::irecv:
      return "MPI_Irecv";
    case detail::HaloMpiOperation::isend:
      return "MPI_Isend";
    case detail::HaloMpiOperation::waitall:
      return "MPI_Waitall";
    case detail::HaloMpiOperation::wait:
      return "MPI_Wait";
  }
  return "unknown";
}

std::string format_failure(detail::HaloFailureRecord record) {
  const char* category =
      record.category ==
              static_cast<std::int64_t>(detail::HaloFailureCategory::post)
          ? "post"
          : "completion";
  std::string message = "halo ";
  message += category;
  message += " failure: rank=";
  message += std::to_string(record.failing_rank);
  message += " operation=";
  message += operation_text(record.operation);
  message += " result=";
  message += std::to_string(record.mpi_result);
  message += " region=";
  message += std::to_string(record.region_index);
  message += " chunk_offset=";
  message += std::to_string(record.chunk_offset);
  message += " chunk_count=";
  message += std::to_string(record.chunk_count);
  message += " tag=";
  message += std::to_string(record.tag);
  return message;
}

}  // namespace

class HaloExchange::Impl final {
 public:
  Impl(ExchangePlan plan, MpiContext context) noexcept
      : plan_(std::move(plan)), context_(std::move(context)) {
    observe_context_generation();
  }

  ~Impl() noexcept {
    if (!active_) {
      return;
    }
    const bool already_complete = all_requests_null(requests_);
    const auto action = detail::active_destruction_action(
        true, detail::mpi_is_active(), already_complete);
    if (action == detail::ActiveDestructionAction::terminate_process) {
      std::terminate();
    }
    if (action == detail::ActiveDestructionAction::drain_requests) {
      if (detail::current_halo_test_options().observe) {
        ++detail::mutable_halo_test_snapshot().destructor_drains;
      }
      const detail::CompletionOutcome outcome = drain_requests_noexcept();
      if (!outcome.requests_proven_null) {
        std::terminate();
      }
    }
    clear_pending();
  }

  Impl(const Impl&) = delete;
  Impl& operator=(const Impl&) = delete;

  int ghost_width() const noexcept { return plan_.ghost_width_; }

  bool is_compatible_with(
      const StructuredDecomposition& decomposition) const {
    if (decomposition.comm() == MPI_COMM_NULL) {
      throw Error("halo compatibility requires a live decomposition");
    }

    int communicator_relation = MPI_UNEQUAL;
    detail::check_mpi(
        MPI_Comm_compare(context_.comm(), decomposition.comm(),
                         &communicator_relation),
        "MPI_Comm_compare halo compatibility");
    if (communicator_relation != MPI_IDENT &&
        communicator_relation != MPI_CONGRUENT) {
      return false;
    }
    if (!same(plan_.local_extent_, decomposition.local_extent())) {
      return false;
    }

    const ExchangePlan expected = ExchangePlan::create(
        decomposition, decomposition.local_extent(), plan_.ghost_width_);
    if (plan_.regions_.size() != expected.regions_.size()) {
      return false;
    }
    for (std::size_t index = 0; index < plan_.regions_.size(); ++index) {
      if (!same(plan_.regions_[index], expected.regions_[index])) {
        return false;
      }
    }
    return true;
  }

  static std::unique_ptr<Impl> create(
      const StructuredDecomposition& decomposition, ExchangePlan plan) {
    detail::require_mpi_active("create halo exchange");
    auto context = MpiContext::duplicate(decomposition.comm());

    bool local_ok = true;
    std::string local_message;
    // Phase 1 is deliberately local-only. No rank may enter a later
    // collective until every rank has converged this validation result.
    try {
      auto expected = ExchangePlan::create(
          decomposition, decomposition.local_extent(), plan.ghost_width());
      if (!same(plan.local_extent_, decomposition.local_extent()) ||
          plan.regions_.size() != expected.regions_.size()) {
        throw Error("halo plan does not describe the supplied decomposition");
      }
      for (std::size_t index = 0; index < plan.regions_.size(); ++index) {
        if (!same(plan.regions_[index], expected.regions_[index])) {
          throw Error(
              "halo plan region does not match the supplied decomposition");
        }
      }
    } catch (const std::exception& error) {
      local_ok = false;
      local_message = error.what();
    } catch (...) {
      local_ok = false;
      local_message = "unknown halo plan validation failure";
    }
    CollectiveStatus status =
        collective_status(context, local_ok, local_message);
    if (!status.ok) {
      throw Error(status.message);
    }

    // Phase 2 has a fixed two-collective schedule. The second reduction
    // carries any first-call failure to every rank, so no rank branches or
    // throws between the scheduled collectives.
    const auto width_options = detail::current_halo_test_options();
    int minimum_width = plan.ghost_width();
    int first_width_result = MPI_Allreduce(
        MPI_IN_PLACE, &minimum_width, 1, MPI_INT, MPI_MIN, context.comm());
    if (first_width_result == MPI_SUCCESS &&
        context.rank() ==
            width_options.inject_plan_width_first_collective_error_rank) {
      first_width_result = MPI_ERR_OTHER;
      if (width_options.observe) {
        ++detail::mutable_halo_test_snapshot()
              .plan_width_first_collective_errors_injected;
      }
    }
    const std::array<int, 2> width_local{
        plan.ghost_width(), first_width_result == MPI_SUCCESS ? 0 : 1};
    std::array<int, 2> width_maximum{};
    if (width_options.observe) {
      ++detail::mutable_halo_test_snapshot()
            .plan_width_second_collective_entries;
    }
    const int second_width_result = MPI_Allreduce(
        width_local.data(), width_maximum.data(),
        static_cast<int>(width_local.size()), MPI_INT, MPI_MAX,
        context.comm());
    if (second_width_result != MPI_SUCCESS) {
      std::terminate();
    }
    if (width_maximum[1] != 0) {
      throw_halo_error(HaloError::agreement_collective_failure);
    }
    if (minimum_width != width_maximum[0]) {
      throw Error("halo plan width differs across communicator ranks");
    }

    // Phase 3 contains rank-local communicator attribute checks and converges
    // them before construction succeeds.
    local_ok = true;
    local_message.clear();
    try {
      void* attribute = nullptr;
      int attribute_present = 0;
      detail::check_mpi(
          MPI_Comm_get_attr(context.comm(), MPI_TAG_UB, &attribute,
                            &attribute_present),
          "MPI_Comm_get_attr");
      static_cast<void>(detail::effective_halo_tag_upper_bound(
          attribute_present != 0, static_cast<const int*>(attribute)));

      MPI_Errhandler handler = MPI_ERRHANDLER_NULL;
      detail::check_mpi(MPI_Comm_get_errhandler(context.comm(), &handler),
                        "MPI_Comm_get_errhandler");
      const bool errors_return = handler == MPI_ERRORS_RETURN;
      if (handler != MPI_ERRHANDLER_NULL) {
        detail::check_mpi(MPI_Errhandler_free(&handler),
                          "MPI_Errhandler_free");
      }
      if (!errors_return) {
        throw Error("halo communicator does not use MPI_ERRORS_RETURN");
      }

      if (detail::current_halo_test_options().observe) {
        int comparison = MPI_UNEQUAL;
        detail::check_mpi(
            MPI_Comm_compare(decomposition.comm(), context.comm(),
                             &comparison),
            "MPI_Comm_compare");
        auto& snapshot = detail::mutable_halo_test_snapshot();
        snapshot.communicator_is_distinct_congruent =
            decomposition.comm() != context.comm() &&
            comparison == MPI_CONGRUENT;
        snapshot.communicator_uses_errors_return = errors_return;
      }
    } catch (const std::exception& error) {
      local_ok = false;
      local_message = error.what();
    } catch (...) {
      local_ok = false;
      local_message = "unknown halo communicator validation failure";
    }

    status = collective_status(context, local_ok, local_message);
    if (!status.ok) {
      throw Error(status.message);
    }

    void* implementation_storage = nullptr;
    local_ok = true;
    local_message.clear();
    try {
      implementation_storage = ::operator new(sizeof(Impl));
    } catch (const std::exception& error) {
      local_ok = false;
      local_message = error.what();
    } catch (...) {
      local_ok = false;
      local_message = "unknown halo object allocation failure";
    }
    status = collective_status(context, local_ok, local_message);
    if (!status.ok) {
      ::operator delete(implementation_storage);
      throw Error(status.message);
    }
    auto* implementation =
        ::new (implementation_storage)
            Impl(std::move(plan), std::move(context));
    return std::unique_ptr<Impl>(implementation);
  }

  void begin(const FieldStorage& storage, FieldId id) {
    detail::require_mpi_active("begin halo exchange");
    const HaloError state_error = converge_fixed_error(
        context_, active_ ? HaloError::begin_while_active : HaloError::none);
    if (state_error != HaloError::none) {
      throw_halo_error(state_error);
    }

    Layout layout;
    HaloError local_error = inspect_layout(storage, id, layout);
    HaloError error = converge_fixed_error(context_, local_error);
    if (error != HaloError::none) {
      throw_halo_error(error);
    }

    error = require_wire_agreement(layout, id);
    if (error != HaloError::none) {
      throw_halo_error(error);
    }

    bool prepared = true;
    try {
      prepare(storage, layout, id);
    } catch (...) {
      prepared = false;
    }
    const CollectiveStatus preparation_status = collective_status(
        context_, prepared,
        prepared ? "" : halo_error_text(HaloError::preparation_failure));
    if (!preparation_status.ok) {
      throw Error(preparation_status.message);
    }

    pending_ = layout;
    pending_id_ = id;
    active_ = true;
    const detail::HaloFailureRecord post_failure = post_all();
    detail::HaloFailureRecord global_post_failure;
    if (!converge_failure_record(context_, post_failure,
                                 global_post_failure)) {
      std::terminate();
    }
    if (has_failure(global_post_failure)) {
      const detail::CompletionOutcome cleanup =
          cancel_and_drain_noexcept();
      const bool local_proven = cleanup.requests_proven_null;
      const int local_value = local_proven ? 1 : 0;
      int every_proven = 0;
      if (MPI_Allreduce(&local_value, &every_proven, 1, MPI_INT, MPI_MIN,
                        context_.comm()) != MPI_SUCCESS ||
          every_proven == 0) {
        std::terminate();
      }
      if (detail::failure_recovery_action(
              true, cleanup.requests_proven_null, false, false) !=
          detail::FailureRecoveryAction::replace_context) {
        std::terminate();
      }
      replace_context_or_terminate();
      clear_pending();
      throw Error(format_failure(global_post_failure));
    }
  }

  void wait(FieldStorage& storage, FieldId id) {
    detail::require_mpi_active("wait for halo exchange");
    const HaloError state_error = converge_fixed_error(
        context_, active_ ? HaloError::none : HaloError::wait_while_idle);
    if (state_error != HaloError::none) {
      throw_halo_error(state_error);
    }

    Layout layout;
    HaloError local_error = inspect_layout(storage, id, layout);
    if (local_error == HaloError::none &&
        (!same_layout(layout, pending_) || id != pending_id_)) {
      local_error = HaloError::wait_target_mismatch;
    }
    const HaloError target_error =
        converge_fixed_error(context_, local_error);
    if (target_error != HaloError::none) {
      throw_halo_error(target_error);
    }

    CompletionAttempt completion =
        wait_all_noexcept(WaitInjection::explicit_completion);
    const bool local_completion_ok =
        detail::completion_succeeded(completion.outcome);
    if (!local_completion_ok && !has_failure(completion.failure)) {
      completion.failure = make_failure_record(
          detail::HaloFailureCategory::completion, context_.rank(),
          detail::HaloMpiOperation::none, MPI_ERR_OTHER, -1, 0U, 0, -1);
    }
    detail::HaloFailureRecord global_completion_failure;
    if (!converge_failure_record(context_, completion.failure,
                                 global_completion_failure)) {
      std::terminate();
    }
    if (has_failure(global_completion_failure)) {
      const bool local_proven = completion.outcome.requests_proven_null;
      const int local_value = local_proven ? 1 : 0;
      int every_proven = 0;
      if (MPI_Allreduce(&local_value, &every_proven, 1, MPI_INT, MPI_MIN,
                        context_.comm()) != MPI_SUCCESS ||
          every_proven == 0 ||
          detail::completion_failure_action(completion.outcome) ==
              detail::CompletionFailureAction::terminate_process) {
        std::terminate();
      }
      if (detail::failure_recovery_action(
              true, completion.outcome.requests_proven_null, false, false) !=
          detail::FailureRecoveryAction::replace_context) {
        std::terminate();
      }
      replace_context_or_terminate();
      clear_pending();
      throw Error(format_failure(global_completion_failure));
    }

    unpack(storage);
    clear_pending();
  }

 private:
  using Entry = FieldStorage::Entry;

  struct Layout {
    const void* storage_token{};
    ScalarType scalar_type{};
    std::size_t scalar_bytes{};
    std::uint32_t components{};
    int field_ghost_width{};
    Int3 interior_extent{};
    std::size_t x_stride{};
    std::size_t y_stride{};
    std::size_t z_stride{};
    std::size_t data_offset{};
    std::size_t storage_bytes{};
  };

  struct RegionBuffers {
    std::vector<std::byte> send;
    std::vector<std::byte> receive;
    std::vector<detail::CountRange> chunks;
  };

  struct CompletionAttempt {
    detail::CompletionOutcome outcome;
    detail::HaloFailureRecord failure;
  };

  HaloError inspect_layout(const FieldStorage& storage, FieldId id,
                           Layout& layout) const noexcept {
    const std::size_t index = static_cast<std::size_t>(id);
    if (!storage.entries_ || index >= storage.entries_->size()) {
      return HaloError::invalid_field;
    }
    const Entry& entry = (*storage.entries_)[index];
    const std::size_t bytes_per_scalar = scalar_size(entry.scalar_type);
    if (bytes_per_scalar == 0U || entry.components == 0U ||
        entry.ghost_width < plan_.ghost_width_ ||
        !same(storage.interior_extent_, plan_.local_extent_)) {
      return HaloError::incompatible_local_layout;
    }

    std::size_t twice_ghost = 0U;
    std::size_t padded_x = 0U;
    std::size_t padded_y = 0U;
    std::size_t padded_z = 0U;
    const std::size_t ghost = static_cast<std::size_t>(entry.ghost_width);
    if (!checked_multiply(ghost, 2U, twice_ghost) ||
        !checked_add(static_cast<std::size_t>(storage.interior_extent_.x),
                     twice_ghost, padded_x) ||
        !checked_add(static_cast<std::size_t>(storage.interior_extent_.y),
                     twice_ghost, padded_y) ||
        !checked_add(static_cast<std::size_t>(storage.interior_extent_.z),
                     twice_ghost, padded_z)) {
      return HaloError::incompatible_local_layout;
    }
    std::size_t expected_y_stride = 0U;
    std::size_t expected_z_stride = 0U;
    std::size_t element_count = 0U;
    std::size_t required_bytes = 0U;
    if (entry.x_stride != static_cast<std::size_t>(entry.components) ||
        !checked_multiply(padded_x, entry.x_stride, expected_y_stride) ||
        entry.y_stride != expected_y_stride ||
        !checked_multiply(padded_y, entry.y_stride, expected_z_stride) ||
        entry.z_stride != expected_z_stride ||
        !checked_multiply(padded_z, entry.z_stride, element_count) ||
        !checked_multiply(element_count, bytes_per_scalar, required_bytes) ||
        entry.data_offset > entry.bytes.size() ||
        required_bytes > entry.bytes.size() - entry.data_offset) {
      return HaloError::incompatible_local_layout;
    }

    layout = Layout{storage.entries_.get(),
                    entry.scalar_type,
                    bytes_per_scalar,
                    entry.components,
                    entry.ghost_width,
                    storage.interior_extent_,
                    entry.x_stride,
                    entry.y_stride,
                    entry.z_stride,
                    entry.data_offset,
                    entry.bytes.size()};
    return HaloError::none;
  }

  HaloError require_wire_agreement(const Layout& layout, FieldId id) const {
    const std::array<std::uint64_t, 5> local{
        static_cast<std::uint64_t>(id),
        static_cast<std::uint64_t>(layout.scalar_type),
        static_cast<std::uint64_t>(layout.components),
        static_cast<std::uint64_t>(layout.field_ghost_width),
        static_cast<std::uint64_t>(plan_.ghost_width_)};
    std::array<std::uint64_t, 5> minimum{};
    const auto options = detail::current_halo_test_options();
    int first_result = MPI_Allreduce(
        local.data(), minimum.data(), static_cast<int>(local.size()),
        MPI_UINT64_T, MPI_MIN, context_.comm());
    if (first_result == MPI_SUCCESS &&
        context_.rank() ==
            options.inject_wire_first_collective_error_rank) {
      first_result = MPI_ERR_OTHER;
      if (options.observe) {
        ++detail::mutable_halo_test_snapshot()
              .wire_first_collective_errors_injected;
      }
    }
    std::array<std::uint64_t, 6> maximum_and_error{};
    std::copy(local.begin(), local.end(), maximum_and_error.begin());
    maximum_and_error.back() = first_result == MPI_SUCCESS ? 0U : 1U;
    if (options.observe) {
      ++detail::mutable_halo_test_snapshot().wire_second_collective_entries;
    }
    if (MPI_Allreduce(MPI_IN_PLACE, maximum_and_error.data(),
                      static_cast<int>(maximum_and_error.size()),
                      MPI_UINT64_T, MPI_MAX, context_.comm()) != MPI_SUCCESS) {
      std::terminate();
    }
    if (maximum_and_error.back() != 0U) {
      return HaloError::agreement_collective_failure;
    }
    return std::equal(minimum.begin(), minimum.end(),
                      maximum_and_error.begin())
               ? HaloError::none
               : HaloError::wire_layout_mismatch;
  }

  void prepare(const FieldStorage& storage, const Layout& layout, FieldId id) {
    const auto options = detail::current_halo_test_options();
    if (options.observe) {
      auto& snapshot = detail::mutable_halo_test_snapshot();
      snapshot.pack_row_copy_events = 0U;
      snapshot.unpack_row_copy_events = 0U;
    }
    std::size_t request_count = 0U;
    for (std::size_t index = 0; index < plan_.regions_.size(); ++index) {
      const ExchangeRegion& region = plan_.regions_[index];
      RegionBuffers& buffers = regions_[index];
      const std::size_t send_bytes = detail::checked_region_payload_bytes(
          region.send_box, layout.components, layout.scalar_bytes);
      const std::size_t receive_bytes = detail::checked_region_payload_bytes(
          region.receive_box, layout.components, layout.scalar_bytes);
      if (send_bytes != receive_bytes) {
        throw Error("halo send and receive region payloads differ");
      }
      const bool active_region =
          region.neighbor_rank != MPI_PROC_NULL && send_bytes != 0U;
      buffers.send.resize(active_region ? send_bytes : 0U);
      buffers.receive.resize(active_region ? receive_bytes : 0U);
      detail::split_count_ranges_into(
          active_region ? send_bytes : 0U, options.chunk_limit,
          "message chunks", buffers.chunks);
      std::size_t region_requests = 0U;
      if (!checked_multiply(buffers.chunks.size(), 2U, region_requests) ||
          !checked_add(request_count, region_requests, request_count)) {
        throw Error("halo aggregate request count overflow");
      }
    }
    if (request_count > requests_.max_size()) {
      throw Error("halo request count exceeds local capacity");
    }
    requests_.assign(request_count, MPI_REQUEST_NULL);
    detail::split_count_ranges_into(request_count, options.waitall_limit,
                                    "MPI_Waitall batches", wait_batches_);

    const Entry& entry =
        (*storage.entries_)[static_cast<std::size_t>(id)];
    for (std::size_t index = 0; index < plan_.regions_.size(); ++index) {
      if (!regions_[index].send.empty()) {
        const std::size_t events = pack_region(
            entry, layout, plan_.regions_[index].send_box,
            regions_[index].send);
        if (options.observe) {
          detail::mutable_halo_test_snapshot().pack_row_copy_events += events;
        }
      }
    }

    if (options.observe) {
      auto& snapshot = detail::mutable_halo_test_snapshot();
      snapshot.context_generation = context_generation_;
      snapshot.all_receives_preceded_sends = true;
      snapshot.chunk_offsets_ordered = true;
      snapshot.receive_posts = 0U;
      snapshot.send_posts = 0U;
      snapshot.send_buffer_capacity = 0U;
      snapshot.receive_buffer_capacity = 0U;
      snapshot.chunk_metadata_capacity = 0U;
      for (const auto& buffers : regions_) {
        snapshot.send_buffer_capacity += buffers.send.capacity();
        snapshot.receive_buffer_capacity += buffers.receive.capacity();
        snapshot.chunk_metadata_capacity += buffers.chunks.capacity();
      }
      snapshot.request_capacity = requests_.capacity();
      snapshot.wait_batch_capacity = wait_batches_.capacity();
    }
  }

  static detail::HaloRowLayout row_layout(const Layout& layout) noexcept {
    return detail::HaloRowLayout{
        layout.field_ghost_width,
        layout.components,
        layout.scalar_bytes,
        layout.x_stride,
        layout.y_stride,
        layout.z_stride,
        layout.storage_bytes - layout.data_offset};
  }

  static std::size_t pack_region(const Entry& entry, const Layout& layout,
                                 Box3 box,
                                 std::vector<std::byte>& buffer) {
    const std::byte* source = entry.bytes.data() + entry.data_offset;
    return detail::pack_halo_region_rows(
        source, row_layout(layout), box, buffer.data(), buffer.size());
  }

  static std::size_t unpack_region(
      Entry& entry, const Layout& layout, Box3 box,
      const std::vector<std::byte>& buffer) noexcept {
    std::byte* destination = entry.bytes.data() + entry.data_offset;
    try {
      return detail::unpack_halo_region_rows(
          destination, row_layout(layout), box, buffer.data(), buffer.size());
    } catch (...) {
      std::terminate();
    }
  }

  detail::HaloFailureRecord post_all() noexcept {
    const auto options = detail::current_halo_test_options();
    bool injected = false;
    std::size_t request_index = 0U;
    detail::HaloFailureRecord failure;
    detail::PostEventState post_events;
    post_events.expected_receive_posts = requests_.size() / 2U;

    for (std::size_t region_index = 0; region_index < regions_.size();
         ++region_index) {
      RegionBuffers& buffers = regions_[region_index];
      const ExchangeRegion& region = plan_.regions_[region_index];
      std::size_t previous_offset = 0U;
      bool first = true;
      for (const auto chunk : buffers.chunks) {
        int mpi_result = MPI_Irecv(
            buffers.receive.data() + chunk.offset, chunk.count, MPI_BYTE,
            region.neighbor_rank, detail::halo_receive_tag(region.offset),
            context_.comm(), &requests_[request_index]);
        detail::observe_post_event(post_events, detail::PostEvent::receive);
        if (mpi_result == MPI_SUCCESS && !injected &&
            context_.rank() == options.inject_post_error_rank) {
          mpi_result = MPI_ERR_OTHER;
          injected = true;
          if (options.observe) {
            ++detail::mutable_halo_test_snapshot().post_errors_injected;
          }
        }
        if (mpi_result != MPI_SUCCESS) {
          if (!has_failure(failure)) {
            failure = make_failure_record(
                detail::HaloFailureCategory::post, context_.rank(),
                detail::HaloMpiOperation::irecv, mpi_result,
                static_cast<int>(region_index), chunk.offset, chunk.count,
                detail::halo_receive_tag(region.offset));
          }
        }
        if (options.observe) {
          auto& snapshot = detail::mutable_halo_test_snapshot();
          if (!first && chunk.offset <= previous_offset) {
            snapshot.chunk_offsets_ordered = false;
          }
        }
        previous_offset = chunk.offset;
        first = false;
        ++request_index;
      }
    }

    for (std::size_t region_index = 0; region_index < regions_.size();
         ++region_index) {
      RegionBuffers& buffers = regions_[region_index];
      const ExchangeRegion& region = plan_.regions_[region_index];
      std::size_t previous_offset = 0U;
      bool first = true;
      for (const auto chunk : buffers.chunks) {
        int mpi_result = MPI_Isend(
            buffers.send.data() + chunk.offset, chunk.count, MPI_BYTE,
            region.neighbor_rank, detail::halo_offset_code(region.offset),
            context_.comm(), &requests_[request_index]);
        detail::observe_post_event(post_events, detail::PostEvent::send);
        if (mpi_result == MPI_SUCCESS && !injected &&
            context_.rank() == options.inject_post_error_rank) {
          mpi_result = MPI_ERR_OTHER;
          injected = true;
          if (options.observe) {
            ++detail::mutable_halo_test_snapshot().post_errors_injected;
          }
        }
        if (mpi_result != MPI_SUCCESS) {
          if (!has_failure(failure)) {
            failure = make_failure_record(
                detail::HaloFailureCategory::post, context_.rank(),
                detail::HaloMpiOperation::isend, mpi_result,
                static_cast<int>(region_index), chunk.offset, chunk.count,
                detail::halo_offset_code(region.offset));
          }
        }
        if (options.observe) {
          auto& snapshot = detail::mutable_halo_test_snapshot();
          if (!first && chunk.offset <= previous_offset) {
            snapshot.chunk_offsets_ordered = false;
          }
        }
        previous_offset = chunk.offset;
        first = false;
        ++request_index;
      }
    }
    if (options.observe) {
      auto& snapshot = detail::mutable_halo_test_snapshot();
      snapshot.all_receives_preceded_sends =
          detail::post_event_sequence_valid(post_events);
      snapshot.receive_posts = post_events.receive_posts;
      snapshot.send_posts = post_events.send_posts;
      snapshot.first_send_sequence = post_events.first_send_sequence;
    }
    return failure;
  }

  enum class WaitInjection { explicit_completion, cleanup };

  CompletionAttempt wait_all_noexcept(
      WaitInjection injection) noexcept {
    const auto options = detail::current_halo_test_options();
    bool injected = false;
    bool mpi_error_seen = false;
    detail::HaloFailureRecord failure;
    int injection_rank = -1;
    if (injection == WaitInjection::explicit_completion) {
      injection_rank = options.inject_wait_error_rank;
    } else if (injection == WaitInjection::cleanup) {
      injection_rank = options.inject_cleanup_wait_error_rank;
    }
    for (const auto batch : wait_batches_) {
      int result = MPI_Waitall(
          batch.count, requests_.data() + batch.offset,
          MPI_STATUSES_IGNORE);
      if (result == MPI_SUCCESS && !injected &&
          context_.rank() == injection_rank) {
        result = MPI_ERR_OTHER;
        injected = true;
        if (injection == WaitInjection::cleanup && options.observe) {
          ++detail::mutable_halo_test_snapshot()
                .cleanup_wait_errors_injected;
        } else if (injection == WaitInjection::explicit_completion &&
                   options.observe) {
          ++detail::mutable_halo_test_snapshot().wait_errors_injected;
        }
      }
      if (result != MPI_SUCCESS) {
        mpi_error_seen = true;
        if (!has_failure(failure)) {
          failure = make_failure_record(
              detail::HaloFailureCategory::completion, context_.rank(),
              detail::HaloMpiOperation::waitall, result, -1, batch.offset,
              batch.count, -1);
        }
        for (int index = 0; index < batch.count; ++index) {
          MPI_Request& request =
              requests_[batch.offset + static_cast<std::size_t>(index)];
          if (request != MPI_REQUEST_NULL) {
            const int wait_result = MPI_Wait(&request, MPI_STATUS_IGNORE);
            if (wait_result != MPI_SUCCESS) {
              mpi_error_seen = true;
              if (!has_failure(failure)) {
                failure = make_failure_record(
                    detail::HaloFailureCategory::completion,
                    context_.rank(), detail::HaloMpiOperation::wait,
                    wait_result, -1,
                    batch.offset + static_cast<std::size_t>(index), 1, -1);
              }
            }
          }
        }
      }
    }
    return CompletionAttempt{
        detail::CompletionOutcome{mpi_error_seen,
                                  all_requests_null(requests_)},
        failure};
  }

  detail::CompletionOutcome drain_requests_noexcept() noexcept {
    return wait_all_noexcept(WaitInjection::cleanup).outcome;
  }

  detail::CompletionOutcome cancel_and_drain_noexcept() noexcept {
    const auto options = detail::current_halo_test_options();
    bool mpi_error_seen = false;
    bool cancellation_calls_succeeded = true;
    for (auto& request : requests_) {
      if (request != MPI_REQUEST_NULL) {
        const int result = MPI_Cancel(&request);
        if (options.observe) {
          ++detail::mutable_halo_test_snapshot().cancel_calls;
        }
        if (result != MPI_SUCCESS) {
          cancellation_calls_succeeded = false;
        }
      }
    }
    const int local_cancellation_ok = cancellation_calls_succeeded ? 1 : 0;
    int every_cancellation_ok = 0;
    if (MPI_Allreduce(&local_cancellation_ok, &every_cancellation_ok, 1,
                      MPI_INT, MPI_MIN, context_.comm()) != MPI_SUCCESS ||
        every_cancellation_ok == 0) {
      std::terminate();
    }
    for (auto& request : requests_) {
      if (request != MPI_REQUEST_NULL) {
        MPI_Status status{};
        int result = MPI_Wait(&request, &status);
        if (result == MPI_SUCCESS) {
          int cancelled = 0;
          result = MPI_Test_cancelled(&status, &cancelled);
          if (options.observe) {
            ++detail::mutable_halo_test_snapshot()
                  .cancellation_status_checks;
          }
          if (result == MPI_SUCCESS) {
            if (cancelled != 0) {
              if (options.observe) {
                ++detail::mutable_halo_test_snapshot().cancelled_requests;
              }
            } else {
              if (options.observe) {
                ++detail::mutable_halo_test_snapshot().completed_requests;
              }
            }
          }
        }
        if (result != MPI_SUCCESS) {
          mpi_error_seen = true;
        }
      }
    }
    return detail::CompletionOutcome{mpi_error_seen,
                                     all_requests_null(requests_)};
  }

  void observe_context_generation() const noexcept {
    if (detail::current_halo_test_options().observe) {
      detail::mutable_halo_test_snapshot().context_generation =
          context_generation_;
    }
  }

  void replace_context_or_terminate() noexcept {
    try {
      if (context_generation_ == std::numeric_limits<std::size_t>::max()) {
        std::terminate();
      }
      MpiContext replacement = MpiContext::duplicate(context_.comm());

      bool local_valid = replacement.comm() != context_.comm();
      int comparison = MPI_UNEQUAL;
      if (MPI_Comm_compare(context_.comm(), replacement.comm(), &comparison) !=
              MPI_SUCCESS ||
          comparison != MPI_CONGRUENT) {
        local_valid = false;
      }

      void* attribute = nullptr;
      int attribute_present = 0;
      if (MPI_Comm_get_attr(replacement.comm(), MPI_TAG_UB, &attribute,
                            &attribute_present) != MPI_SUCCESS) {
        local_valid = false;
      } else {
        try {
          static_cast<void>(detail::effective_halo_tag_upper_bound(
              attribute_present != 0,
              static_cast<const int*>(attribute)));
        } catch (...) {
          local_valid = false;
        }
      }

      MPI_Errhandler handler = MPI_ERRHANDLER_NULL;
      if (MPI_Comm_get_errhandler(replacement.comm(), &handler) !=
          MPI_SUCCESS) {
        local_valid = false;
      } else {
        if (handler != MPI_ERRORS_RETURN) {
          local_valid = false;
        }
        if (handler != MPI_ERRHANDLER_NULL &&
            MPI_Errhandler_free(&handler) != MPI_SUCCESS) {
          local_valid = false;
        }
      }

      const int local_value = local_valid ? 1 : 0;
      int every_valid = 0;
      if (MPI_Allreduce(&local_value, &every_valid, 1, MPI_INT, MPI_MIN,
                        replacement.comm()) != MPI_SUCCESS ||
          every_valid == 0) {
        std::terminate();
      }

      MpiContext retired(std::move(context_));
      context_ = std::move(replacement);
      ++context_generation_;
      if (detail::current_halo_test_options().observe) {
        auto& snapshot = detail::mutable_halo_test_snapshot();
        snapshot.context_generation = context_generation_;
        ++snapshot.context_replacements;
        snapshot.last_context_replacement_distinct_congruent = true;
      }
    } catch (...) {
      std::terminate();
    }
  }

  void unpack(FieldStorage& storage) noexcept {
    Entry& entry =
        (*storage.entries_)[static_cast<std::size_t>(pending_id_)];
    for (std::size_t index = 0; index < regions_.size(); ++index) {
      if (!regions_[index].receive.empty()) {
        const std::size_t events = unpack_region(
            entry, pending_, plan_.regions_[index].receive_box,
            regions_[index].receive);
        if (detail::current_halo_test_options().observe) {
          detail::mutable_halo_test_snapshot().unpack_row_copy_events +=
              events;
        }
      }
    }
  }

  static bool same_layout(const Layout& left, const Layout& right) noexcept {
    return left.storage_token == right.storage_token &&
           left.scalar_type == right.scalar_type &&
           left.scalar_bytes == right.scalar_bytes &&
           left.components == right.components &&
           left.field_ghost_width == right.field_ghost_width &&
           same(left.interior_extent, right.interior_extent) &&
           left.x_stride == right.x_stride &&
           left.y_stride == right.y_stride &&
           left.z_stride == right.z_stride &&
           left.data_offset == right.data_offset &&
           left.storage_bytes == right.storage_bytes;
  }

  void clear_pending() noexcept {
    active_ = false;
    pending_ = Layout{};
    pending_id_ = 0U;
    requests_.clear();
    wait_batches_.clear();
  }

  ExchangePlan plan_;
  MpiContext context_;
  std::array<RegionBuffers, 26> regions_;
  std::vector<MPI_Request> requests_;
  std::vector<detail::CountRange> wait_batches_;
  Layout pending_{};
  FieldId pending_id_{};
  bool active_{};
  std::size_t context_generation_{1U};
};

HaloExchange HaloExchange::create(
    const StructuredDecomposition& decomposition, ExchangePlan plan) {
  return HaloExchange(Impl::create(decomposition, std::move(plan)));
}

HaloExchange::HaloExchange(std::unique_ptr<Impl> implementation) noexcept
    : implementation_(std::move(implementation)) {}

HaloExchange::~HaloExchange() noexcept = default;

HaloExchange::HaloExchange(HaloExchange&&) noexcept = default;

int HaloExchange::ghost_width() const {
  if (!implementation_) {
    throw_halo_error(HaloError::moved_from);
  }
  return implementation_->ghost_width();
}

bool HaloExchange::is_compatible_with(
    const StructuredDecomposition& decomposition) const {
  detail::require_mpi_active("validate halo compatibility");
  if (!implementation_) {
    throw_halo_error(HaloError::moved_from);
  }
  return implementation_->is_compatible_with(decomposition);
}

void HaloExchange::exchange(FieldStorage& storage, FieldId id) {
  begin(static_cast<const FieldStorage&>(storage), id);
  wait(storage, id);
}

void HaloExchange::begin(const FieldStorage& storage, FieldId id) {
  if (!implementation_) {
    throw_halo_error(HaloError::moved_from);
  }
  implementation_->begin(storage, id);
}

void HaloExchange::wait(FieldStorage& storage, FieldId id) {
  if (!implementation_) {
    throw_halo_error(HaloError::moved_from);
  }
  implementation_->wait(storage, id);
}

}  // namespace hundun::runtime

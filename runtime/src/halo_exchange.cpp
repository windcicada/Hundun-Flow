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
#include <cstring>
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
  post_failure = 9,
  completion_failure = 10
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
    case HaloError::post_failure:
      return "halo nonblocking post failed and was drained";
    case HaloError::completion_failure:
      return "halo completion failed and received data was discarded";
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

}  // namespace

class HaloExchange::Impl final {
 public:
  Impl(ExchangePlan plan, MpiContext context) noexcept
      : plan_(std::move(plan)), context_(std::move(context)) {}

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
      if (!drain_requests_noexcept()) {
        std::terminate();
      }
    }
    clear_pending();
  }

  Impl(const Impl&) = delete;
  Impl& operator=(const Impl&) = delete;

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

    // Phase 2 is collective only after the local plan phase succeeded.
    local_ok = true;
    local_message.clear();
    try {
      int minimum_width = plan.ghost_width();
      int maximum_width = plan.ghost_width();
      detail::check_mpi(MPI_Allreduce(MPI_IN_PLACE, &minimum_width, 1,
                                      MPI_INT, MPI_MIN, context.comm()),
                        "MPI_Allreduce");
      detail::check_mpi(MPI_Allreduce(MPI_IN_PLACE, &maximum_width, 1,
                                      MPI_INT, MPI_MAX, context.comm()),
                        "MPI_Allreduce");
      if (minimum_width != maximum_width) {
        throw Error("halo plan width differs across communicator ranks");
      }
    } catch (const std::exception& error) {
      local_ok = false;
      local_message = error.what();
    } catch (...) {
      local_ok = false;
      local_message = "unknown halo width agreement failure";
    }
    status = collective_status(context, local_ok, local_message);
    if (!status.ok) {
      throw Error(status.message);
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
      const bool has_attribute = attribute_present != 0 && attribute != nullptr;
      const int upper_bound =
          has_attribute ? *static_cast<int*>(attribute) : 0;
      detail::validate_halo_tag_upper_bound(has_attribute, upper_bound);

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
    const HaloError post_error = post_all();
    const HaloError global_post_error =
        converge_fixed_error(context_, post_error);
    if (global_post_error != HaloError::none) {
      const bool local_proven = cancel_and_drain_noexcept();
      const int local_value = local_proven ? 1 : 0;
      int every_proven = 0;
      if (MPI_Allreduce(&local_value, &every_proven, 1, MPI_INT, MPI_MIN,
                        context_.comm()) != MPI_SUCCESS ||
          every_proven == 0) {
        std::terminate();
      }
      clear_pending();
      throw_halo_error(global_post_error);
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

    bool local_completion_ok = wait_all_noexcept(true);
    const HaloError completion_error = converge_fixed_error(
        context_, local_completion_ok ? HaloError::none
                                      : HaloError::completion_failure);
    if (completion_error != HaloError::none) {
      const bool local_proven = all_requests_null(requests_);
      const int local_value = local_proven ? 1 : 0;
      int every_proven = 0;
      if (MPI_Allreduce(&local_value, &every_proven, 1, MPI_INT, MPI_MIN,
                        context_.comm()) != MPI_SUCCESS ||
          every_proven == 0 ||
          detail::completion_failure_action(local_proven) ==
              detail::CompletionFailureAction::terminate_process) {
        std::terminate();
      }
      clear_pending();
      throw_halo_error(completion_error);
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
    std::array<std::uint64_t, 5> maximum{};
    detail::check_mpi(MPI_Allreduce(local.data(), minimum.data(),
                                    static_cast<int>(local.size()),
                                    MPI_UINT64_T, MPI_MIN, context_.comm()),
                      "MPI_Allreduce");
    detail::check_mpi(MPI_Allreduce(local.data(), maximum.data(),
                                    static_cast<int>(local.size()),
                                    MPI_UINT64_T, MPI_MAX, context_.comm()),
                      "MPI_Allreduce");
    return minimum == maximum ? HaloError::none
                              : HaloError::wire_layout_mismatch;
  }

  void prepare(const FieldStorage& storage, const Layout& layout, FieldId id) {
    const auto options = detail::current_halo_test_options();
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
        pack_region(entry, layout, plan_.regions_[index].send_box,
                    regions_[index].send);
      }
    }

    if (options.observe) {
      auto& snapshot = detail::mutable_halo_test_snapshot();
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

  static std::size_t element_byte_offset(const Layout& layout, int i, int j,
                                         int k, std::uint32_t component) {
    // inspect_layout proves the complete padded element/byte count fits in
    // size_t. Plan boxes are a subset of that validated coordinate domain, so
    // these products are bounded by the already-checked total allocation.
    const auto ghost = static_cast<std::int64_t>(layout.field_ghost_width);
    const auto ii = static_cast<std::size_t>(static_cast<std::int64_t>(i) + ghost);
    const auto jj = static_cast<std::size_t>(static_cast<std::int64_t>(j) + ghost);
    const auto kk = static_cast<std::size_t>(static_cast<std::int64_t>(k) + ghost);
    const std::size_t linear =
        kk * layout.z_stride + jj * layout.y_stride +
        ii * layout.x_stride + static_cast<std::size_t>(component);
    return linear * layout.scalar_bytes;
  }

  static void pack_region(const Entry& entry, const Layout& layout, Box3 box,
                          std::vector<std::byte>& buffer) noexcept {
    const std::byte* source = entry.bytes.data() + entry.data_offset;
    std::size_t cursor = 0U;
    for (int k = box.begin.z; k < box.end.z; ++k) {
      for (int j = box.begin.y; j < box.end.y; ++j) {
        for (int i = box.begin.x; i < box.end.x; ++i) {
          for (std::uint32_t component = 0U; component < layout.components;
               ++component) {
            const std::size_t source_offset =
                element_byte_offset(layout, i, j, k, component);
            std::memcpy(buffer.data() + cursor, source + source_offset,
                        layout.scalar_bytes);
            cursor += layout.scalar_bytes;
          }
        }
      }
    }
  }

  static void unpack_region(Entry& entry, const Layout& layout, Box3 box,
                            const std::vector<std::byte>& buffer) noexcept {
    std::byte* destination = entry.bytes.data() + entry.data_offset;
    std::size_t cursor = 0U;
    for (int k = box.begin.z; k < box.end.z; ++k) {
      for (int j = box.begin.y; j < box.end.y; ++j) {
        for (int i = box.begin.x; i < box.end.x; ++i) {
          for (std::uint32_t component = 0U; component < layout.components;
               ++component) {
            const std::size_t destination_offset =
                element_byte_offset(layout, i, j, k, component);
            std::memcpy(destination + destination_offset,
                        buffer.data() + cursor, layout.scalar_bytes);
            cursor += layout.scalar_bytes;
          }
        }
      }
    }
  }

  HaloError post_all() noexcept {
    const auto options = detail::current_halo_test_options();
    bool injected = false;
    bool send_seen = false;
    std::size_t request_index = 0U;
    HaloError result = HaloError::none;

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
        if (mpi_result == MPI_SUCCESS && !injected &&
            context_.rank() == options.inject_post_error_rank) {
          mpi_result = MPI_ERR_OTHER;
          injected = true;
        }
        if (mpi_result != MPI_SUCCESS) {
          result = HaloError::post_failure;
        }
        if (options.observe) {
          auto& snapshot = detail::mutable_halo_test_snapshot();
          ++snapshot.receive_posts;
          if (send_seen) {
            snapshot.all_receives_preceded_sends = false;
          }
          if (!first && chunk.offset <= previous_offset) {
            snapshot.chunk_offsets_ordered = false;
          }
        }
        previous_offset = chunk.offset;
        first = false;
        ++request_index;
      }
    }

    send_seen = true;
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
        if (mpi_result == MPI_SUCCESS && !injected &&
            context_.rank() == options.inject_post_error_rank) {
          mpi_result = MPI_ERR_OTHER;
          injected = true;
        }
        if (mpi_result != MPI_SUCCESS) {
          result = HaloError::post_failure;
        }
        if (options.observe) {
          auto& snapshot = detail::mutable_halo_test_snapshot();
          ++snapshot.send_posts;
          if (!first && chunk.offset <= previous_offset) {
            snapshot.chunk_offsets_ordered = false;
          }
        }
        previous_offset = chunk.offset;
        first = false;
        ++request_index;
      }
    }
    static_cast<void>(send_seen);
    return result;
  }

  bool wait_all_noexcept(bool allow_injection) noexcept {
    const auto options = detail::current_halo_test_options();
    bool injected = false;
    bool success = true;
    for (const auto batch : wait_batches_) {
      int result = MPI_Waitall(
          batch.count, requests_.data() + batch.offset,
          MPI_STATUSES_IGNORE);
      if (result == MPI_SUCCESS && allow_injection && !injected &&
          context_.rank() == options.inject_wait_error_rank) {
        result = MPI_ERR_OTHER;
        injected = true;
      }
      if (result != MPI_SUCCESS) {
        success = false;
        for (int index = 0; index < batch.count; ++index) {
          MPI_Request& request =
              requests_[batch.offset + static_cast<std::size_t>(index)];
          if (request != MPI_REQUEST_NULL &&
              MPI_Wait(&request, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
            success = false;
          }
        }
      }
    }
    return success && all_requests_null(requests_);
  }

  bool drain_requests_noexcept() noexcept {
    return wait_all_noexcept(false);
  }

  bool cancel_and_drain_noexcept() noexcept {
    bool success = true;
    for (auto& request : requests_) {
      if (request != MPI_REQUEST_NULL) {
        // A cancellation error does not by itself make the backing buffer
        // unsafe. The subsequent successful wait and null handle are the
        // completion proof that controls recovery.
        static_cast<void>(MPI_Cancel(&request));
      }
    }
    for (auto& request : requests_) {
      if (request != MPI_REQUEST_NULL &&
          MPI_Wait(&request, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
        success = false;
      }
    }
    return success && all_requests_null(requests_);
  }

  void unpack(FieldStorage& storage) noexcept {
    Entry& entry =
        (*storage.entries_)[static_cast<std::size_t>(pending_id_)];
    for (std::size_t index = 0; index < regions_.size(); ++index) {
      if (!regions_[index].receive.empty()) {
        unpack_region(entry, pending_, plan_.regions_[index].receive_box,
                      regions_[index].receive);
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
};

HaloExchange HaloExchange::create(
    const StructuredDecomposition& decomposition, ExchangePlan plan) {
  return HaloExchange(Impl::create(decomposition, std::move(plan)));
}

HaloExchange::HaloExchange(std::unique_ptr<Impl> implementation) noexcept
    : implementation_(std::move(implementation)) {}

HaloExchange::~HaloExchange() noexcept = default;

HaloExchange::HaloExchange(HaloExchange&&) noexcept = default;

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

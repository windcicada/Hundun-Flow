// SPDX-License-Identifier: Apache-2.0

#include "hundun/exec_execution.hpp"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <memory>
#include <mutex>
#include <new>
#include <string>
#include <utility>
#include <vector>

namespace hundun::execution {
namespace {

constexpr BackendIdentity kCpuReferenceIdentity =
    UINT64_C(0x48554e44554e4350);

std::mutex identity_mutex;
AllocationIdentity next_identity = 1;
#ifdef HUNDUN_EXECUTION_ENABLE_TEST_ACCESS
std::atomic<bool> inject_allocation_failure{false};
std::atomic<bool> inject_completed_event_allocation_failure{false};
#endif

struct AllocationCounterState final {
  AllocationCounters counters;
  std::uint64_t outstanding_events{};
  std::uint64_t outstanding_bytes{};
};

std::mutex allocation_counter_mutex;
AllocationCounterState allocation_counter_state;

std::uint64_t checked_counter_add(std::uint64_t left, std::uint64_t right,
                                  const char* message) {
  if (right > std::numeric_limits<std::uint64_t>::max() - left) {
    throw runtime::Error(message);
  }
  return left + right;
}

AllocationCounterState allocation_candidate(std::uint64_t bytes) {
  auto candidate = allocation_counter_state;
  candidate.counters.allocation_events = checked_counter_add(
      candidate.counters.allocation_events, 1U,
      "allocation event counter would overflow");
  candidate.counters.allocated_bytes = checked_counter_add(
      candidate.counters.allocated_bytes, bytes,
      "allocated byte counter would overflow");
  candidate.counters.live_bytes = checked_counter_add(
      candidate.counters.live_bytes, bytes,
      "live allocation byte counter would overflow");
  candidate.counters.peak_live_bytes =
      std::max(candidate.counters.peak_live_bytes,
               candidate.counters.live_bytes);
  candidate.outstanding_events = checked_counter_add(
      candidate.outstanding_events, 1U,
      "outstanding allocation counter would overflow");
  candidate.outstanding_bytes = checked_counter_add(
      candidate.outstanding_bytes, bytes,
      "outstanding allocation byte counter would overflow");
  static_cast<void>(checked_counter_add(
      candidate.counters.deallocation_events, candidate.outstanding_events,
      "future deallocation event counter would overflow"));
  static_cast<void>(checked_counter_add(
      candidate.counters.deallocated_bytes, candidate.outstanding_bytes,
      "future deallocation byte counter would overflow"));
  return candidate;
}

void preflight_allocation(std::uint64_t bytes) {
  const std::lock_guard<std::mutex> lock(allocation_counter_mutex);
  static_cast<void>(allocation_candidate(bytes));
}

void register_allocation(std::uint64_t bytes) {
  const std::lock_guard<std::mutex> lock(allocation_counter_mutex);
  auto candidate = allocation_candidate(bytes);
  allocation_counter_state = candidate;
}

void retire_allocation_noexcept(std::uint64_t bytes) noexcept {
  const std::lock_guard<std::mutex> lock(allocation_counter_mutex);
  auto& state = allocation_counter_state;
  if (state.outstanding_events == 0U || state.outstanding_bytes < bytes ||
      state.counters.live_bytes < bytes ||
      state.counters.deallocation_events ==
          std::numeric_limits<std::uint64_t>::max() ||
      state.counters.deallocated_bytes >
          std::numeric_limits<std::uint64_t>::max() - bytes) {
    std::terminate();
  }
  ++state.counters.deallocation_events;
  state.counters.deallocated_bytes += bytes;
  state.counters.live_bytes -= bytes;
  --state.outstanding_events;
  state.outstanding_bytes -= bytes;
}

std::size_t checked_multiply(std::size_t left, std::size_t right,
                             const char* message) {
  if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left) {
    throw runtime::Error(message);
  }
  return left * right;
}

std::size_t checked_add(std::size_t left, std::size_t right,
                        const char* message) {
  if (right > std::numeric_limits<std::size_t>::max() - left) {
    throw runtime::Error(message);
  }
  return left + right;
}

AllocationIdentity issue_allocation_identity() {
  const std::lock_guard<std::mutex> lock(identity_mutex);
  if (next_identity == 0 ||
      next_identity == std::numeric_limits<AllocationIdentity>::max()) {
    throw runtime::Error("allocation identity counter would wrap");
  }
  const AllocationIdentity issued = next_identity;
  ++next_identity;
  return issued;
}

struct AlignedDelete final {
  void operator()(std::byte* pointer) const noexcept {
    ::operator delete(pointer,
                      std::align_val_t{alignof(std::max_align_t)});
  }
};

void validate_supported_byte_size(std::size_t byte_size) {
  if (byte_size >
      static_cast<std::size_t>(std::numeric_limits<std::ptrdiff_t>::max())) {
    throw runtime::Error("buffer byte size exceeds the supported range");
  }
}

std::unique_ptr<std::byte, AlignedDelete> allocate_bytes(
    std::size_t byte_size) {
  validate_supported_byte_size(byte_size);
#ifdef HUNDUN_EXECUTION_ENABLE_TEST_ACCESS
  if (inject_allocation_failure.exchange(false)) {
    throw runtime::Error("buffer allocation failed by the test seam");
  }
#endif
  if (byte_size == 0) {
    return {};
  }
  try {
    auto storage = std::unique_ptr<std::byte, AlignedDelete>(
        static_cast<std::byte*>(::operator new(
            byte_size, std::align_val_t{alignof(std::max_align_t)})));
    for (std::size_t offset = 0; byte_size - offset >= sizeof(double);
         offset += sizeof(double)) {
      ::new (static_cast<void*>(storage.get() + offset)) double{};
    }
    return storage;
  } catch (const std::bad_alloc&) {
    throw runtime::Error("buffer allocation failed");
  }
}

enum class EventResult : unsigned char { pending, success, error };

}  // namespace

namespace detail {

struct AllocationControl final {
  AllocationControl(std::unique_ptr<std::byte, AlignedDelete> storage_value,
                    std::size_t byte_size_value,
                    AllocationIdentity allocation_identity_value,
                    std::uint64_t epoch_value,
                    BackendIdentity backend_identity_value,
                    ExecutionSpace space_value) noexcept
      : storage(std::move(storage_value)),
        byte_size(byte_size_value),
        allocation_identity(allocation_identity_value),
        epoch(epoch_value),
        backend_identity(backend_identity_value),
        space(space_value) {}

  ~AllocationControl() {
    if (counted) {
      retire_allocation_noexcept(static_cast<std::uint64_t>(byte_size));
    }
  }

  std::unique_ptr<std::byte, AlignedDelete> storage;
  std::size_t byte_size;
  AllocationIdentity allocation_identity;
  std::uint64_t epoch;
  BackendIdentity backend_identity;
  ExecutionSpace space;
  bool owner_live{true};
  bool counted{false};
};

struct EventState final {
  std::atomic<EventResult> result{EventResult::pending};
  std::mutex mutex;
  std::condition_variable condition;
  std::size_t waiter_count{0};
  std::string error_message;
  std::vector<std::shared_ptr<AllocationControl>> retained_allocations;
};

void validate_view_metadata(
    const std::shared_ptr<AllocationControl>& control,
    AllocationIdentity allocation_identity, std::uint64_t epoch,
    BackendIdentity backend_identity, ExecutionSpace space,
    std::size_t offset_bytes, std::size_t element_count,
    std::size_t stride_elements, bool require_host) {
  if (!control || !control->owner_live) {
    throw runtime::Error("vector view owner allocation is not live");
  }
  if (allocation_identity == 0 ||
      control->allocation_identity != allocation_identity ||
      control->epoch != epoch) {
    throw runtime::Error("vector view allocation identity or epoch is not live");
  }
  if (backend_identity == 0 ||
      control->backend_identity != backend_identity || control->space != space) {
    throw runtime::Error("vector view metadata does not match its allocation");
  }
  if (stride_elements == 0) {
    throw runtime::Error("vector view stride must be positive");
  }
  if (offset_bytes % alignof(double) != 0) {
    throw runtime::Error("vector view offset violates double alignment");
  }
  if (offset_bytes > control->byte_size) {
    throw runtime::Error("vector view range exceeds its allocation");
  }
  if (element_count != 0) {
    const std::size_t stride_bytes = checked_multiply(
        stride_elements, sizeof(double), "vector view stride overflow");
    const std::size_t final_displacement = checked_multiply(
        element_count - 1, stride_bytes, "vector view range overflow");
    const std::size_t final_begin = checked_add(
        offset_bytes, final_displacement, "vector view range overflow");
    const std::size_t final_end = checked_add(
        final_begin, sizeof(double), "vector view range overflow");
    if (final_end > control->byte_size) {
      throw runtime::Error("vector view range exceeds its allocation");
    }
  }
  if (require_host && space != ExecutionSpace::host) {
    throw runtime::Error("host access to a device vector view is forbidden");
  }
}

double* checked_view_element(
    const std::weak_ptr<AllocationControl>& weak_control,
    AllocationIdentity allocation_identity, std::uint64_t epoch,
    BackendIdentity backend_identity, ExecutionSpace space,
    std::size_t offset_bytes, std::size_t element_count,
    std::size_t stride_elements, bool mutable_access_allowed,
    std::size_t index, bool require_index) {
  const auto control = weak_control.lock();
  validate_view_metadata(control, allocation_identity, epoch, backend_identity,
                         space, offset_bytes, element_count, stride_elements,
                         true);
  if (!mutable_access_allowed) {
    throw runtime::Error(
        "mutable vector view access requires writable capability");
  }
  if (require_index && index >= element_count) {
    throw runtime::Error("vector view index is out of bounds");
  }
  if (element_count == 0) {
    return nullptr;
  }
  const std::size_t byte_stride = checked_multiply(
      stride_elements, sizeof(double), "vector view stride overflow");
  const std::size_t displacement = checked_multiply(
      index, byte_stride, "vector view index overflow");
  const std::size_t byte_offset = checked_add(
      offset_bytes, displacement, "vector view index overflow");
  return reinterpret_cast<double*>(control->storage.get() + byte_offset);
}

}  // namespace detail

std::string_view CpuReferenceContext::backend_name() const noexcept {
  return "cpu_reference";
}

AllocationCounters allocation_counters() noexcept {
  const std::lock_guard<std::mutex> lock(allocation_counter_mutex);
  return allocation_counter_state.counters;
}

BackendIdentity CpuReferenceContext::backend_identity() const noexcept {
  return kCpuReferenceIdentity;
}

ExecutionSpace CpuReferenceContext::space() const noexcept {
  return ExecutionSpace::host;
}

bool CpuReferenceContext::ordered() const noexcept { return true; }

bool CpuReferenceContext::supports(
    ExecutionCapability capability) const noexcept {
  switch (capability) {
    case ExecutionCapability::buffer_allocation:
    case ExecutionCapability::host_access:
    case ExecutionCapability::transfer:
      return true;
    case ExecutionCapability::asynchronous_event:
      return false;
  }
  return false;
}

namespace {

std::shared_ptr<detail::AllocationControl> make_host_control(
    BackendIdentity backend_identity, ExecutionSpace space,
    std::size_t byte_size, std::uint64_t epoch) {
  validate_supported_byte_size(byte_size);
  preflight_allocation(static_cast<std::uint64_t>(byte_size));
  auto storage = allocate_bytes(byte_size);
  const auto identity = issue_allocation_identity();
  try {
    auto control = std::make_shared<detail::AllocationControl>(
        std::move(storage), byte_size, identity, epoch, backend_identity,
        space);
    register_allocation(static_cast<std::uint64_t>(byte_size));
    control->counted = true;
    return control;
  } catch (const std::bad_alloc&) {
    throw runtime::Error("buffer allocation control creation failed");
  }
}

std::shared_ptr<detail::AllocationControl> make_control(
    ExecutionContext& context, std::size_t byte_size, std::uint64_t epoch) {
  if (context.backend_identity() == 0) {
    throw runtime::Error("buffer context backend identity must be nonzero");
  }
  if (context.space() != ExecutionSpace::host ||
      !context.supports(ExecutionCapability::host_access)) {
    throw runtime::Error("buffer allocation requires host access");
  }
  if (!context.supports(ExecutionCapability::buffer_allocation)) {
    throw runtime::Error("context does not support buffer allocation");
  }
  return make_host_control(context.backend_identity(), context.space(),
                           byte_size, epoch);
}

std::shared_ptr<detail::EventState> inline_success_state() {
  try {
    auto state = std::make_shared<detail::EventState>();
    state->result.store(EventResult::success, std::memory_order_release);
    return state;
  } catch (const std::bad_alloc&) {
    throw runtime::Error("execution event allocation failed");
  }
}

}  // namespace

Buffer::Buffer(ExecutionContext& context, std::size_t byte_size)
    : control_(make_control(context, byte_size, 1)) {}

Buffer::~Buffer() { invalidate(); }

Buffer::Buffer(Buffer&& other) noexcept : control_(std::move(other.control_)) {}

Buffer& Buffer::operator=(Buffer&& other) noexcept {
  if (this != &other) {
    invalidate();
    control_ = std::move(other.control_);
  }
  return *this;
}

std::size_t Buffer::byte_size() const noexcept {
  return control_ ? control_->byte_size : 0;
}

AllocationIdentity Buffer::allocation_identity() const noexcept {
  return control_ ? control_->allocation_identity : 0;
}

std::uint64_t Buffer::epoch() const noexcept {
  return control_ ? control_->epoch : 0;
}

BackendIdentity Buffer::backend_identity() const noexcept {
  return control_ ? control_->backend_identity : 0;
}

ExecutionSpace Buffer::space() const noexcept {
  return control_ ? control_->space : ExecutionSpace::host;
}

void Buffer::reallocate(std::size_t byte_size) {
  require_live();
  if (control_->epoch == std::numeric_limits<std::uint64_t>::max()) {
    throw runtime::Error("buffer epoch would wrap");
  }

  auto replacement = make_host_control(
      control_->backend_identity, control_->space, byte_size,
      control_->epoch + 1);
  control_->owner_live = false;
  control_ = std::move(replacement);
}

VectorView<double> Buffer::view(std::size_t offset_elements,
                                std::size_t element_count,
                                std::size_t stride_elements) {
  validate_view_range(offset_elements, element_count, stride_elements);
  const std::size_t offset_bytes = checked_multiply(
      offset_elements, sizeof(double), "vector view offset overflow");
  return VectorView<double>(control_, offset_bytes, element_count,
                            stride_elements, control_->allocation_identity,
                            control_->epoch, control_->backend_identity,
                            control_->space, true);
}

VectorView<const double> Buffer::view(std::size_t offset_elements,
                                      std::size_t element_count,
                                      std::size_t stride_elements) const {
  validate_view_range(offset_elements, element_count, stride_elements);
  const std::size_t offset_bytes = checked_multiply(
      offset_elements, sizeof(double), "vector view offset overflow");
  return VectorView<const double>(
      control_, offset_bytes, element_count, stride_elements,
      control_->allocation_identity, control_->epoch,
      control_->backend_identity, control_->space, false);
}

void Buffer::invalidate() noexcept {
  if (control_) {
    control_->owner_live = false;
    control_.reset();
  }
}

void Buffer::require_live() const {
  if (!control_ || !control_->owner_live) {
    throw runtime::Error("moved-from buffer has no live allocation");
  }
}

void Buffer::validate_view_range(std::size_t offset_elements,
                                 std::size_t element_count,
                                 std::size_t stride_elements) const {
  require_live();
  if (stride_elements == 0) {
    throw runtime::Error("vector view stride must be positive");
  }
  const std::size_t offset_bytes = checked_multiply(
      offset_elements, sizeof(double), "vector view offset overflow");
  detail::validate_view_metadata(
      control_, control_->allocation_identity, control_->epoch,
      control_->backend_identity, control_->space, offset_bytes, element_count,
      stride_elements, false);
}

bool ExecutionEvent::ready() const noexcept {
  return state_ &&
         state_->result.load(std::memory_order_acquire) != EventResult::pending;
}

ExecutionEvent ExecutionEvent::completed() {
#ifdef HUNDUN_EXECUTION_ENABLE_TEST_ACCESS
  if (inject_completed_event_allocation_failure.exchange(false)) {
    throw runtime::Error(
        "completed event allocation failed by the test seam");
  }
#endif
  return ExecutionEvent(inline_success_state());
}

void ExecutionEvent::wait() const {
  if (!state_) {
    throw runtime::Error("moved-from execution event has no state");
  }
  std::unique_lock<std::mutex> lock(state_->mutex);
  if (state_->result.load(std::memory_order_acquire) == EventResult::pending) {
    ++state_->waiter_count;
    state_->condition.notify_all();
    state_->condition.wait(lock, [&] {
      return state_->result.load(std::memory_order_acquire) !=
             EventResult::pending;
    });
    --state_->waiter_count;
  }
  const EventResult result = state_->result.load(std::memory_order_acquire);
  if (result == EventResult::error) {
    throw runtime::Error(state_->error_message);
  }
}

ExecutionEvent transfer(VectorView<const double> source,
                        VectorView<double> destination,
                        ExecutionContext& context) {
  const auto source_control = source.control_.lock();
  const auto destination_control = destination.control_.lock();
  detail::validate_view_metadata(
      source_control, source.allocation_identity_, source.epoch_,
      source.backend_identity_, source.space_, source.offset_bytes_,
      source.element_count_, source.stride_elements_, false);
  detail::validate_view_metadata(
      destination_control, destination.allocation_identity_, destination.epoch_,
      destination.backend_identity_, destination.space_,
      destination.offset_bytes_, destination.element_count_,
      destination.stride_elements_, false);

  if (source.element_count_ != destination.element_count_) {
    throw runtime::Error("transfer source and destination sizes differ");
  }
  if (!destination.writable_) {
    throw runtime::Error("transfer destination is not writable");
  }
  if (source.scalar_format() != ScalarFormat::float64 ||
      destination.scalar_format() != ScalarFormat::float64) {
    throw runtime::Error("transfer requires float64 vector views");
  }
  if (context.backend_identity() == 0) {
    throw runtime::Error("transfer context backend identity must be nonzero");
  }
  if (!context.supports(ExecutionCapability::transfer)) {
    throw runtime::Error("context does not support transfer");
  }
  if (!context.supports(ExecutionCapability::host_access) ||
      context.space() != ExecutionSpace::host) {
    throw runtime::Error("transfer requires a host context");
  }
  if (context.backend_identity() != source.backend_identity_ ||
      context.backend_identity() != destination.backend_identity_) {
    throw runtime::Error("transfer context backend identity does not match views");
  }
  if (source.space_ != ExecutionSpace::host ||
      destination.space_ != ExecutionSpace::host) {
    throw runtime::Error("transfer requires host vector views");
  }

  const bool same_allocation =
      source.allocation_identity_ == destination.allocation_identity_;
  const bool exact_same_view =
      same_allocation && source.offset_bytes_ == destination.offset_bytes_ &&
      source.element_count_ == destination.element_count_ &&
      source.stride_elements_ == destination.stride_elements_;
  if (exact_same_view) {
    return ExecutionEvent(inline_success_state());
  }

  if (same_allocation && source.element_count_ != 0) {
    const auto span_end = [](std::size_t offset, std::size_t count,
                             std::size_t stride) {
      const auto stride_bytes = checked_multiply(
          stride, sizeof(double), "transfer span overflow");
      const auto last = checked_add(
          offset,
          checked_multiply(count - 1, stride_bytes,
                           "transfer span overflow"),
          "transfer span overflow");
      return checked_add(last, sizeof(double), "transfer span overflow");
    };
    const std::size_t source_end = span_end(
        source.offset_bytes_, source.element_count_, source.stride_elements_);
    const std::size_t destination_end =
        span_end(destination.offset_bytes_, destination.element_count_,
                 destination.stride_elements_);
    if (source.offset_bytes_ < destination_end &&
        destination.offset_bytes_ < source_end) {
      throw runtime::Error("transfer byte spans overlap");
    }
  }

  auto completion = inline_success_state();
  double* source_pointer = detail::checked_view_element(
      source.control_, source.allocation_identity_, source.epoch_,
      source.backend_identity_, source.space_, source.offset_bytes_,
      source.element_count_, source.stride_elements_, true, 0, false);
  double* destination_pointer = detail::checked_view_element(
      destination.control_, destination.allocation_identity_,
      destination.epoch_, destination.backend_identity_, destination.space_,
      destination.offset_bytes_, destination.element_count_,
      destination.stride_elements_, destination.writable_, 0, false);
  for (std::size_t index = 0; index < source.element_count_; ++index) {
    destination_pointer[index * destination.stride_elements_] =
        source_pointer[index * source.stride_elements_];
  }
  return ExecutionEvent(std::move(completion));
}

#ifdef HUNDUN_EXECUTION_ENABLE_TEST_ACCESS
namespace detail {

AllocationIdentity execution_next_allocation_identity_raw() noexcept {
  const std::lock_guard<std::mutex> lock(identity_mutex);
  return next_identity;
}

void execution_set_allocation_counters_raw(
    AllocationCounters counters) noexcept {
  const std::lock_guard<std::mutex> lock(allocation_counter_mutex);
  if (counters.live_bytes != allocation_counter_state.outstanding_bytes)
    std::terminate();
  allocation_counter_state.counters = counters;
}

void execution_set_next_allocation_identity_raw(
    AllocationIdentity next) noexcept {
  const std::lock_guard<std::mutex> lock(identity_mutex);
  next_identity = next;
}

void execution_fail_next_allocation_raw() noexcept {
  inject_allocation_failure.store(true);
}

void execution_fail_next_completed_event_allocation_raw() noexcept {
  inject_completed_event_allocation_failure.store(true);
}

void execution_set_control_epoch_raw(
    const std::shared_ptr<AllocationControl> &control,
    std::uint64_t epoch) {
  if (!control || !control->owner_live)
    throw runtime::Error("moved-from buffer has no live allocation");
  if (epoch == 0)
    throw runtime::Error("buffer epoch must be nonzero");
  control->epoch = epoch;
}

std::shared_ptr<AllocationControl> execution_make_metadata_control_raw(
    std::size_t byte_size, AllocationIdentity allocation_identity,
    std::uint64_t epoch, BackendIdentity backend_identity,
    ExecutionSpace space) {
  if (allocation_identity == 0 || epoch == 0 || backend_identity == 0)
    throw runtime::Error(
        "metadata-only allocation identities and epoch must be nonzero");
  if (space != ExecutionSpace::device)
    throw runtime::Error(
        "metadata-only allocation test seam requires device space");
  return std::make_shared<AllocationControl>(
      std::unique_ptr<std::byte, AlignedDelete>{}, byte_size,
      allocation_identity, epoch, backend_identity, space);
}

bool execution_control_has_storage_raw(
    const std::shared_ptr<AllocationControl> &control) noexcept {
  return control && control->storage != nullptr;
}

bool execution_control_owner_live_raw(
    const std::shared_ptr<AllocationControl> &control) noexcept {
  return control && control->owner_live;
}

std::shared_ptr<EventState> execution_make_pending_event_raw(
    std::vector<std::shared_ptr<AllocationControl>> retained_allocations) {
  auto state = std::make_shared<EventState>();
  state->retained_allocations = std::move(retained_allocations);
  return state;
}

void execution_complete_success_raw(
    const std::shared_ptr<EventState> &state) {
  if (!state)
    throw runtime::Error("moved-from execution event has no state");
  {
    const std::lock_guard<std::mutex> lock(state->mutex);
    if (state->result.load(std::memory_order_relaxed) != EventResult::pending)
      throw runtime::Error("execution event is already complete");
    state->retained_allocations.clear();
    state->result.store(EventResult::success, std::memory_order_release);
  }
  state->condition.notify_all();
}

void execution_complete_error_raw(const std::shared_ptr<EventState> &state,
                                  std::string message) {
  if (!state)
    throw runtime::Error("moved-from execution event has no state");
  {
    const std::lock_guard<std::mutex> lock(state->mutex);
    if (state->result.load(std::memory_order_relaxed) != EventResult::pending)
      throw runtime::Error("execution event is already complete");
    state->error_message = std::move(message);
    state->retained_allocations.clear();
    state->result.store(EventResult::error, std::memory_order_release);
  }
  state->condition.notify_all();
}

void execution_wait_until_waiter_count_raw(
    const std::shared_ptr<EventState> &state, std::size_t expected) {
  if (!state)
    throw runtime::Error("moved-from execution event has no state");
  std::unique_lock<std::mutex> lock(state->mutex);
  state->condition.wait(lock, [&] {
    return state->waiter_count >= expected ||
           state->result.load(std::memory_order_acquire) !=
               EventResult::pending;
  });
  if (state->waiter_count < expected)
    throw runtime::Error(
        "execution event completed before the expected waiter arrived");
}

} // namespace detail
#endif

}  // namespace hundun::execution

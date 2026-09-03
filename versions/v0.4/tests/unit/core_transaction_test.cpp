// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_execution.hpp"
#include "hundun/v04_field.hpp"

#include <mpi.h>

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <string_view>
#include <type_traits>
#include <vector>

namespace allocation_observer {

std::atomic<bool> enabled{false};
std::atomic<std::size_t> count{0U};

void observe() noexcept {
  if (enabled.load(std::memory_order_relaxed)) {
    count.fetch_add(1U, std::memory_order_relaxed);
  }
}

void* allocate(std::size_t size) {
  observe();
  void* result = std::malloc(size == 0U ? 1U : size);
  if (result == nullptr) {
    throw std::bad_alloc{};
  }
  return result;
}

void* allocate_aligned(std::size_t size, std::size_t alignment) {
  observe();
  void* result = nullptr;
  const std::size_t requested = size == 0U ? alignment : size;
  if (posix_memalign(&result, alignment, requested) != 0 || result == nullptr) {
    throw std::bad_alloc{};
  }
  return result;
}

class Guard {
 public:
  Guard() noexcept {
    count.store(0U, std::memory_order_relaxed);
    enabled.store(true, std::memory_order_release);
  }
  ~Guard() { enabled.store(false, std::memory_order_release); }

  Guard(const Guard&) = delete;
  Guard& operator=(const Guard&) = delete;
};

}  // namespace allocation_observer

void* operator new(std::size_t size) {
  return allocation_observer::allocate(size);
}

void* operator new[](std::size_t size) {
  return allocation_observer::allocate(size);
}

void* operator new(std::size_t size, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}

void* operator new[](std::size_t size, std::align_val_t alignment) {
  return allocation_observer::allocate_aligned(
      size, static_cast<std::size_t>(alignment));
}

void operator delete(void* pointer) noexcept { std::free(pointer); }
void operator delete[](void* pointer) noexcept { std::free(pointer); }
void operator delete(void* pointer, std::size_t) noexcept { std::free(pointer); }
void operator delete[](void* pointer, std::size_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete(void* pointer, std::size_t, std::align_val_t) noexcept {
  std::free(pointer);
}
void operator delete[](void* pointer, std::size_t,
                       std::align_val_t) noexcept {
  std::free(pointer);
}

namespace {

using hundun::v04::ArenaFieldRequest;
using hundun::v04::ArenaLayout;
using hundun::v04::AttemptTransaction;
using hundun::v04::AttemptFinishDecision;
using hundun::v04::ConstFaceFluxView;
using hundun::v04::FaceFluxStorage;
using hundun::v04::FaceFluxView;
using hundun::v04::FieldId;
using hundun::v04::FieldLifetime;
using hundun::v04::FieldPlacement;
using hundun::v04::FieldRegistry;
using hundun::v04::FieldSchema;
using hundun::v04::FieldView;
using hundun::v04::Int3;
using hundun::v04::PendingCacheStamp;
using hundun::v04::PendingFaceFluxView;
using hundun::v04::PreparedAttemptFinish;
using hundun::v04::RevisionSlotId;
using hundun::v04::RevisionDependency;
using hundun::v04::RevisionSourceId;
using hundun::v04::RevisionToken;
using hundun::v04::Span;
using hundun::v04::StateLayers;
using hundun::v04::StateRole;
using hundun::v04::Status;
using hundun::v04::StatusCode;
using hundun::v04::FinalFaceFluxAuthority;
using hundun::v04::FinalFaceFluxWriter;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

std::uint64_t checksum(FieldView view) {
  std::uint64_t hash = 14695981039346656037ULL;
  for (std::uint8_t component = 0U; component < view.components; ++component) {
    for (std::int32_t z = -view.ghosts.z;
         z < view.interior.z + view.ghosts.z; ++z) {
      for (std::int32_t y = -view.ghosts.y;
           y < view.interior.y + view.ghosts.y; ++y) {
        for (std::int32_t x = -view.ghosts.x;
             x < view.interior.x + view.ghosts.x; ++x) {
          std::uint64_t bits = 0U;
          const double value = view.unchecked(Int3{x, y, z}, component);
          std::memcpy(&bits, &value, sizeof(bits));
          for (std::size_t byte = 0U; byte < sizeof(bits); ++byte) {
            hash ^= (bits >> (byte * 8U)) & 0xffU;
            hash *= 1099511628211ULL;
          }
        }
      }
    }
  }
  return hash;
}

void fill(FieldView view, double seed) {
  double value = seed;
  for (std::uint8_t component = 0U; component < view.components; ++component) {
    for (std::int32_t z = -view.ghosts.z;
         z < view.interior.z + view.ghosts.z; ++z) {
      for (std::int32_t y = -view.ghosts.y;
           y < view.interior.y + view.ghosts.y; ++y) {
        for (std::int32_t x = -view.ghosts.x;
             x < view.interior.x + view.ghosts.x; ++x) {
          view.unchecked(Int3{x, y, z}, component) = value;
          value += 0.125;
        }
      }
    }
  }
}

bool make_layers(StateLayers& layers, FieldId& velocity, FieldId& pressure) {
  FieldRegistry registry;
  FieldSchema schema;
  bool passed = true;
  passed &= expect(
      static_cast<bool>(registry.declare_field("velocity", 3U, 1U, velocity)),
      "velocity declaration succeeds");
  passed &= expect(
      static_cast<bool>(registry.declare_field("pressure", 1U, 1U, pressure)),
      "pressure declaration succeeds");
  passed &= expect(static_cast<bool>(registry.freeze(schema)),
                   "transaction schema freezes");
  const std::array requests{
      ArenaFieldRequest{velocity, Int3{7, 4, 3}, FieldPlacement{0U},
                        FieldLifetime::state_layer},
      ArenaFieldRequest{pressure, Int3{7, 4, 3}, FieldPlacement{0U},
                        FieldLifetime::state_layer},
  };
  ArenaLayout layout;
  passed &= expect(static_cast<bool>(ArenaLayout::compile(
                       schema,
                       Span<const ArenaFieldRequest>{requests.data(),
                                                     requests.size()},
                       layout)),
                   "transaction arena compiles");
  passed &= expect(static_cast<bool>(StateLayers::allocate(layout, layers)),
                   "state layers allocate once from the arena");
  return passed;
}

bool get(StateLayers& layers, StateRole role, FieldId id, FieldView& out,
         std::string_view description) {
  return expect(static_cast<bool>(layers.view(role, id, out)), description);
}

RevisionDependency depends_on(const AttemptTransaction& transaction,
                              FieldId field) {
  return RevisionDependency{
      AttemptTransaction::field_revision_source(field),
      transaction.trial_revision(field)};
}

bool overwrite_pending_flux(PendingFaceFluxView& pending, double value) {
  return static_cast<bool>(
      hundun::v04::detail::overwrite_pending_face_flux_for_test(pending,
                                                                value));
}

std::vector<std::uint64_t> flux_bits(ConstFaceFluxView flux) {
  std::vector<std::uint64_t> bits;
  const std::array faces{flux.x, flux.y, flux.z};
  for (const auto face : faces) {
    const std::size_t count = static_cast<std::size_t>(face.extents.x) *
                              static_cast<std::size_t>(face.extents.y) *
                              static_cast<std::size_t>(face.extents.z);
    bits.reserve(bits.size() + count);
    for (std::int32_t z = 0; z < face.extents.z; ++z)
      for (std::int32_t y = 0; y < face.extents.y; ++y)
        for (std::int32_t x = 0; x < face.extents.x; ++x) {
          std::uint64_t value_bits = 0U;
          const double value = face.unchecked({x, y, z});
          std::memcpy(&value_bits, &value, sizeof(value_bits));
          bits.push_back(value_bits);
        }
  }
  return bits;
}

bool test_commit_rollback_retry_and_hot_allocations() {
  FieldId velocity{};
  FieldId pressure{};
  StateLayers layers;
  bool passed = make_layers(layers, velocity, pressure);

  FieldView accepted_n;
  FieldView accepted_nm1;
  FieldView trial;
  passed &= get(layers, StateRole::accepted_n, velocity, accepted_n,
                "accepted n view is available");
  passed &= get(layers, StateRole::accepted_n_minus_one, velocity,
                accepted_nm1, "accepted n-1 view is available");
  passed &= get(layers, StateRole::trial, velocity, trial,
                "trial view is available");
  fill(accepted_n, 10.0);
  fill(accepted_nm1, 20.0);
  fill(trial, 30.0);

  const std::uint64_t accepted_crc = checksum(accepted_n);
  const std::uint64_t previous_crc = checksum(accepted_nm1);
  double* const accepted_base = accepted_n.base;
  double* const previous_base = accepted_nm1.base;
  double* const trial_base = trial.base;
  const std::size_t old_n = layers.handle(StateRole::accepted_n);
  const std::size_t old_nm1 = layers.handle(StateRole::accepted_n_minus_one);
  const std::size_t old_trial = layers.handle(StateRole::trial);
  passed &= expect(old_n != old_nm1 && old_n != old_trial &&
                       old_nm1 != old_trial,
                   "three state roles start on distinct replicas");

  AttemptTransaction transaction;
  passed &= expect(static_cast<bool>(AttemptTransaction::create(
                       layers.field_count(), 3U, layers.field_count(),
                       transaction)),
                   "transaction preallocates its pending-cache log");
  const auto counters_before = layers.counters();
  const RevisionToken trial_before_begin =
      layers.revision(StateRole::trial, velocity);

  Status begin_status;
  Status revise_status;
  Status cache_status;
  Status finish_status;
  std::size_t hot_allocations = std::numeric_limits<std::size_t>::max();
  {
    allocation_observer::Guard guard;
    begin_status = transaction.begin(layers);
    revise_status = transaction.revise_trial(velocity);
    if (revise_status) {
      revise_status = transaction.revise_trial(pressure);
    }
    const std::array dependencies{depends_on(transaction, pressure)};
    cache_status = transaction.publish_pending_cache(
        RevisionSlotId{1U},
        Span<const RevisionDependency>{dependencies.data(),
                                       dependencies.size()},
        PendingCacheStamp{77U});
    passed &= expect(!transaction.pending_cache_valid(RevisionSlotId{1U}),
                     "pending cache is invisible before collective commit");
    finish_status = transaction.collective_finish(MPI_COMM_SELF, Status{});
    hot_allocations = allocation_observer::count.load(std::memory_order_relaxed);
  }
  passed &= expect(static_cast<bool>(begin_status) &&
                       static_cast<bool>(revise_status) &&
                       static_cast<bool>(cache_status) &&
                       static_cast<bool>(finish_status),
                   "successful transaction hot operations return success");
  passed &= expect(hot_allocations == 0U,
                   "begin/work/cache/collective commit allocate zero times");
  passed &= expect(transaction.committed(),
                   "successful consensus marks the transaction committed");
  passed &= expect(transaction.lowest_failing_rank() == -1,
                   "successful consensus exposes no failing rank");
  passed &= expect(layers.counters().whole_field_copies ==
                       counters_before.whole_field_copies,
                   "successful begin/commit performs no whole-field copy");
  passed &= expect(layers.handle(StateRole::accepted_n) == old_trial &&
                       layers.handle(StateRole::accepted_n_minus_one) == old_n &&
                       layers.handle(StateRole::trial) == old_nm1,
                   "success rotates trial to n, n to n-1, and n-1 to trial");

  FieldView committed_n;
  FieldView committed_nm1;
  FieldView new_trial;
  passed &= get(layers, StateRole::accepted_n, velocity, committed_n,
                "committed n view is available");
  passed &= get(layers, StateRole::accepted_n_minus_one, velocity,
                committed_nm1, "committed n-1 view is available");
  passed &= get(layers, StateRole::trial, velocity, new_trial,
                "rotated trial view is available");
  passed &= expect(committed_n.base == trial_base &&
                       committed_nm1.base == accepted_base &&
                       new_trial.base == previous_base,
                   "commit is a zero-copy handle rotation");
  passed &= expect(checksum(committed_nm1) == accepted_crc,
                   "old accepted n bytes become n-1 unchanged");
  passed &= expect(checksum(new_trial) == previous_crc,
                   "old n-1 bytes become the reusable trial unchanged");
  passed &= expect(layers.revision(StateRole::accepted_n, velocity) !=
                       trial_before_begin,
                   "trial work publishes a fresh revision on commit");
  passed &= expect(transaction.pending_cache_valid(RevisionSlotId{1U}) &&
                       transaction.pending_cache(RevisionSlotId{1U}) ==
                           77U,
                   "pending cache metadata publishes only after commit");

  const std::size_t rollback_n = layers.handle(StateRole::accepted_n);
  const std::size_t rollback_nm1 =
      layers.handle(StateRole::accepted_n_minus_one);
  const std::size_t rollback_trial = layers.handle(StateRole::trial);
  FieldView rollback_n_view;
  FieldView rollback_nm1_view;
  get(layers, StateRole::accepted_n, velocity, rollback_n_view,
      "pre-rollback n view is available");
  get(layers, StateRole::accepted_n_minus_one, velocity, rollback_nm1_view,
      "pre-rollback n-1 view is available");
  const std::uint64_t rollback_n_crc = checksum(rollback_n_view);
  const std::uint64_t rollback_nm1_crc = checksum(rollback_nm1_view);

  Status rollback_begin;
  Status rollback_revise;
  Status rollback_cache;
  Status rollback_finish;
  FieldView stale_rejected_trial;
  hot_allocations = std::numeric_limits<std::size_t>::max();
  {
    allocation_observer::Guard guard;
    rollback_begin = transaction.begin(layers);
    const Status view_status = layers.view(StateRole::trial, velocity,
                                           stale_rejected_trial);
    if (view_status) {
      stale_rejected_trial.unchecked(Int3{0, 0, 0}, 0U) = -1234.0;
    }
    rollback_revise = transaction.revise_trial(velocity);
    const std::array dependencies{depends_on(transaction, pressure)};
    rollback_cache = transaction.publish_pending_cache(
        RevisionSlotId{2U},
        Span<const RevisionDependency>{dependencies.data(),
                                       dependencies.size()},
        PendingCacheStamp{999U});
    passed &= expect(transaction.pending_cache_valid(RevisionSlotId{1U}) &&
                         !transaction.pending_cache_valid(RevisionSlotId{2U}),
                     "an attempt exposes committed cache only, not pending cache");
    rollback_finish = transaction.collective_finish(
        MPI_COMM_SELF, Status{StatusCode::rejected_step, 41U});
    hot_allocations = allocation_observer::count.load(std::memory_order_relaxed);
  }
  passed &= expect(static_cast<bool>(rollback_begin) &&
                       static_cast<bool>(rollback_revise) &&
                       static_cast<bool>(rollback_cache),
                   "failure attempt begins and records bounded metadata");
  passed &= expect(rollback_finish.code == StatusCode::rejected_step &&
                       rollback_finish.detail == 41U,
                   "single-rank failure status is returned exactly");
  passed &= expect(hot_allocations == 0U,
                   "begin/work/cache/collective rollback allocate zero times");
  passed &= expect(!transaction.committed(),
                   "failed consensus marks the transaction uncommitted");
  passed &= expect(transaction.lowest_failing_rank() == 0,
                   "single-rank failure reports failing rank zero");
  passed &= expect(layers.handle(StateRole::accepted_n) == rollback_n &&
                       layers.handle(StateRole::accepted_n_minus_one) ==
                           rollback_nm1 &&
                       layers.handle(StateRole::trial) == rollback_trial,
                   "rollback preserves every state role handle");
  get(layers, StateRole::accepted_n, velocity, rollback_n_view,
      "post-rollback n view is available");
  get(layers, StateRole::accepted_n_minus_one, velocity, rollback_nm1_view,
      "post-rollback n-1 view is available");
  passed &= expect(checksum(rollback_n_view) == rollback_n_crc &&
                       checksum(rollback_nm1_view) == rollback_nm1_crc,
                   "rollback preserves accepted n and n-1 bytes");
  passed &= expect(transaction.pending_cache_valid(RevisionSlotId{1U}) &&
                       transaction.pending_cache(RevisionSlotId{1U}) != 0U &&
                       !transaction.pending_cache_valid(RevisionSlotId{2U}),
                   "rollback keeps committed cache and discards pending metadata");

  const RevisionToken rejected_trial_revision =
      layers.revision(StateRole::trial, velocity);
  passed &= expect(static_cast<bool>(transaction.begin(layers)),
                   "retry begins without rebuilding the transaction");
  passed &= expect(layers.revision(StateRole::trial, velocity) !=
                       rejected_trial_revision,
                   "retry invalidates the rejected trial revision");
  const RevisionToken retry_trial_revision =
      layers.revision(StateRole::trial, velocity);
  passed &= expect(
      layers.checked_ptr(StateRole::trial, stale_rejected_trial,
                         Int3{0, 0, 0}, 0U) == nullptr,
      "a rejected trial view is stale against retry's authoritative token");
  FieldView retry_trial_view;
  passed &= expect(
      static_cast<bool>(layers.view(StateRole::trial, velocity,
                                    retry_trial_view)) &&
          retry_trial_view.revision == retry_trial_revision &&
          layers.checked_ptr(StateRole::trial, retry_trial_view,
                             Int3{0, 0, 0}, 0U) != nullptr,
      "the retry view matches its authoritative role revision");
  passed &= expect(!transaction.pending_cache_valid(RevisionSlotId{2U}),
                   "retry cannot observe rejected pending metadata");
  passed &= expect(static_cast<bool>(transaction.revise_trial(velocity)) &&
                       static_cast<bool>(transaction.revise_trial(pressure)),
                   "retry revises every state-layer field before commit");
  const std::array retry_dependencies{depends_on(transaction, pressure)};
  passed &= expect(static_cast<bool>(transaction.publish_pending_cache(
                       RevisionSlotId{0U},
                       Span<const RevisionDependency>{
                           retry_dependencies.data(),
                           retry_dependencies.size()},
                       PendingCacheStamp{123U})),
                   "retry records fresh pending metadata");
  passed &= expect(static_cast<bool>(
                       transaction.collective_finish(MPI_COMM_SELF, Status{})) &&
                       transaction.committed(),
                   "retry commits normally");
  passed &= expect(transaction.pending_cache_valid(RevisionSlotId{0U}) &&
                       transaction.pending_cache(RevisionSlotId{0U}) ==
                           123U,
                   "retry publishes only its fresh cache metadata");
  passed &= expect(layers.counters().whole_field_copies == 0U,
                   "the complete commit/rollback/retry sequence copies no field");
  return passed;
}

bool test_incomplete_write_set_rolls_back() {
  FieldId velocity{};
  FieldId pressure{};
  StateLayers layers;
  bool passed = make_layers(layers, velocity, pressure);
  AttemptTransaction transaction;
  passed &= expect(static_cast<bool>(AttemptTransaction::create(
                       layers.field_count(), 2U, layers.field_count(),
                       transaction)),
                   "incomplete-write transaction preallocates");

  const std::size_t old_n = layers.handle(StateRole::accepted_n);
  const std::size_t old_nm1 = layers.handle(StateRole::accepted_n_minus_one);
  const std::size_t old_trial = layers.handle(StateRole::trial);
  passed &= expect(static_cast<bool>(transaction.begin(layers)) &&
                       static_cast<bool>(transaction.revise_trial(velocity)),
                   "incomplete attempt revises only one of two fields");
  const Status result = transaction.collective_finish(MPI_COMM_SELF, Status{});
  passed &= expect(result.code == StatusCode::invalid_plan,
                   "local success cannot commit an incomplete state write set");
  passed &= expect(!transaction.committed() &&
                       transaction.lowest_failing_rank() == 0,
                   "incomplete write set is an authoritative local failure");
  passed &= expect(layers.handle(StateRole::accepted_n) == old_n &&
                       layers.handle(StateRole::accepted_n_minus_one) ==
                           old_nm1 &&
                       layers.handle(StateRole::trial) == old_trial,
                   "incomplete write set rolls back every role handle");
  return passed;
}

bool test_collective_prepare_failure_discards_attempt() {
  FieldId velocity{};
  FieldId pressure{};
  StateLayers layers;
  bool passed = make_layers(layers, velocity, pressure);
  AttemptTransaction transaction;
  passed &= expect(static_cast<bool>(AttemptTransaction::create(
                       layers.field_count(), 1U, layers.field_count(),
                       transaction)),
                   "prepare-failure transaction preallocates");
  passed &= expect(static_cast<bool>(transaction.begin(layers)) &&
                       static_cast<bool>(transaction.revise_trial(velocity)) &&
                       static_cast<bool>(transaction.revise_trial(pressure)),
                   "prepare-failure attempt completes its local write set");
  const std::array dependencies{depends_on(transaction, pressure)};
  passed &= expect(static_cast<bool>(transaction.publish_pending_cache(
                       RevisionSlotId{0U},
                       Span<const RevisionDependency>{dependencies.data(),
                                                      dependencies.size()},
                       PendingCacheStamp{901U})),
                   "prepare-failure attempt stages side state");

  PreparedAttemptFinish unavailable;
  const Status failed = transaction.collective_prepare(
      MPI_COMM_NULL, Status{}, unavailable);
  passed &= expect(failed.code == StatusCode::invalid_plan &&
                       !unavailable.valid(),
                   "invalid communicator cannot mint a commit certificate");
  passed &= expect(!transaction.active() && transaction.finished() &&
                       !transaction.committed() &&
                       transaction.lowest_failing_rank() == -1 &&
                       !transaction.pending_cache_valid(RevisionSlotId{0U}),
                   "collective prepare failure closes the attempt and "
                   "discards staged side state");

  passed &= expect(static_cast<bool>(transaction.begin(layers)) &&
                       static_cast<bool>(transaction.revise_trial(velocity)) &&
                       static_cast<bool>(transaction.revise_trial(pressure)),
                   "a clean retry begins after collective prepare failure");
  const std::array retry_dependencies{depends_on(transaction, pressure)};
  passed &= expect(static_cast<bool>(transaction.publish_pending_cache(
                       RevisionSlotId{0U},
                       Span<const RevisionDependency>{
                           retry_dependencies.data(),
                           retry_dependencies.size()},
                       PendingCacheStamp{902U})) &&
                       static_cast<bool>(transaction.collective_finish(
                           MPI_COMM_SELF, Status{})) &&
                       transaction.committed() &&
                       transaction.pending_cache(RevisionSlotId{0U}) == 902U,
                   "retry commits only freshly staged side state");
  return passed;
}

bool test_ignored_hot_error_is_latched() {
  bool passed = true;

  {
    FieldId velocity{};
    FieldId pressure{};
    StateLayers layers;
    passed &= make_layers(layers, velocity, pressure);
    AttemptTransaction transaction;
    passed &= expect(
        static_cast<bool>(AttemptTransaction::create(
            layers.field_count(), 2U, layers.field_count(), transaction)),
        "duplicate-cache transaction preallocates");
    const std::size_t old_n = layers.handle(StateRole::accepted_n);
    passed &= expect(static_cast<bool>(transaction.begin(layers)) &&
                         static_cast<bool>(transaction.revise_trial(velocity)) &&
                         static_cast<bool>(transaction.revise_trial(pressure)),
                     "duplicate-cache attempt is valid before the duplicate");
    const std::array dependencies{depends_on(transaction, pressure)};
    passed &= expect(static_cast<bool>(transaction.publish_pending_cache(
                         RevisionSlotId{0U},
                         Span<const RevisionDependency>{dependencies.data(),
                                                        dependencies.size()},
                         PendingCacheStamp{11U})),
                     "first cache publication succeeds");
    const Status duplicate = transaction.publish_pending_cache(
        RevisionSlotId{0U},
        Span<const RevisionDependency>{dependencies.data(),
                                       dependencies.size()},
        PendingCacheStamp{12U});
    passed &= expect(duplicate.code == StatusCode::invalid_plan,
                     "duplicate pending-cache publication returns an error");
    const Status finish =
        transaction.collective_finish(MPI_COMM_SELF, Status{});
    passed &= expect(finish.code == StatusCode::invalid_plan &&
                         !transaction.committed() &&
                         transaction.lowest_failing_rank() == 0,
                     "ignored duplicate-cache error is latched into rollback");
    passed &= expect(layers.handle(StateRole::accepted_n) == old_n &&
                         !transaction.pending_cache_valid(RevisionSlotId{0U}),
                     "duplicate-cache failure rotates nothing and publishes nothing");
  }

  {
    FieldId velocity{};
    FieldId pressure{};
    StateLayers layers;
    passed &= make_layers(layers, velocity, pressure);
    AttemptTransaction transaction;
    passed &= expect(
        static_cast<bool>(AttemptTransaction::create(
            layers.field_count(), 2U, layers.field_count(), transaction)),
        "out-of-range-cache transaction preallocates");
    const std::size_t old_n = layers.handle(StateRole::accepted_n);
    passed &= expect(static_cast<bool>(transaction.begin(layers)) &&
                         static_cast<bool>(transaction.revise_trial(velocity)) &&
                         static_cast<bool>(transaction.revise_trial(pressure)),
                     "out-of-range-cache attempt completes its field writes");
    const std::array dependencies{depends_on(transaction, pressure)};
    const Status out_of_range = transaction.publish_pending_cache(
        RevisionSlotId{2U},
        Span<const RevisionDependency>{dependencies.data(),
                                       dependencies.size()},
        PendingCacheStamp{21U});
    passed &= expect(out_of_range.code == StatusCode::invalid_plan,
                     "out-of-range cache slot returns an error");
    const Status finish =
        transaction.collective_finish(MPI_COMM_SELF, Status{});
    passed &= expect(finish.code == StatusCode::invalid_plan &&
                         !transaction.committed() &&
                         transaction.lowest_failing_rank() == 0 &&
                         layers.handle(StateRole::accepted_n) == old_n,
                     "ignored out-of-range slot error is latched into rollback");
  }

  {
    FieldId velocity{};
    FieldId pressure{};
    StateLayers layers;
    passed &= make_layers(layers, velocity, pressure);
    AttemptTransaction transaction;
    passed &= expect(
        static_cast<bool>(AttemptTransaction::create(
            layers.field_count(), 2U, layers.field_count(), transaction)),
        "invalid-field transaction preallocates");
    const std::size_t old_n = layers.handle(StateRole::accepted_n);
    passed &= expect(static_cast<bool>(transaction.begin(layers)),
                     "invalid-field attempt begins");
    const Status invalid_field = transaction.revise_trial(
        std::numeric_limits<FieldId>::max());
    passed &= expect(invalid_field.code == StatusCode::invalid_plan,
                     "invalid FieldId returns an error");
    const Status finish =
        transaction.collective_finish(MPI_COMM_SELF, Status{});
    passed &= expect(finish.code == StatusCode::invalid_plan &&
                         !transaction.committed() &&
                         transaction.lowest_failing_rank() == 0 &&
                         layers.handle(StateRole::accepted_n) == old_n,
                     "ignored invalid-FieldId error is latched into rollback");
  }
  return passed;
}

bool test_cache_stamp_and_layer_binding() {
  FieldId velocity{};
  FieldId pressure{};
  StateLayers layers;
  bool passed = make_layers(layers, velocity, pressure);
  AttemptTransaction transaction;
  passed &= expect(static_cast<bool>(AttemptTransaction::create(
                       layers.field_count(), 2U, layers.field_count() + 1U,
                       transaction)),
                   "cache-stamp transaction preallocates");

  passed &= expect(static_cast<bool>(transaction.begin(layers)) &&
                       static_cast<bool>(transaction.revise_trial(velocity)) &&
                       static_cast<bool>(transaction.revise_trial(pressure)),
                   "cache-stamp attempt completes its declared write set");
  const std::array independent_dependencies{
      depends_on(transaction, velocity)};
  passed &= expect(static_cast<bool>(transaction.publish_pending_cache(
                       RevisionSlotId{0U},
                       Span<const RevisionDependency>{
                           independent_dependencies.data(),
                           independent_dependencies.size()},
                       PendingCacheStamp{500U})),
                   "cache publication accepts its exact field dependency");
  const RevisionSourceId geometry_source = 70001U;
  RevisionToken geometry_revision = 900U;
  passed &= expect(static_cast<bool>(transaction.bind_dependency(
                       RevisionDependency{geometry_source,
                                          geometry_revision})),
                   "a non-field geometry dependency joins the authority catalog");
  passed &= expect(static_cast<bool>(transaction.revise_trial(pressure)),
                   "an unrelated trial field can change after cache production");
  const Status independent_cache_finish =
      transaction.collective_finish(MPI_COMM_SELF, Status{});
  passed &= expect(static_cast<bool>(independent_cache_finish) &&
                       transaction.committed() &&
                       transaction.pending_cache(RevisionSlotId{0U}) == 500U,
                   "an unrelated field revision preserves a valid cache");

  const std::size_t stale_old_n = layers.handle(StateRole::accepted_n);
  passed &= expect(static_cast<bool>(transaction.begin(layers)) &&
                       static_cast<bool>(transaction.revise_trial(velocity)) &&
                       static_cast<bool>(transaction.revise_trial(pressure)),
                   "stale-cache attempt completes its declared write set");
  passed &= expect(static_cast<bool>(transaction.bind_dependency(
                       RevisionDependency{geometry_source,
                                          geometry_revision})),
                   "geometry dependency is rebound for the new attempt");
  const std::array stale_dependencies{
      depends_on(transaction, velocity),
      RevisionDependency{geometry_source, geometry_revision}};
  passed &= expect(static_cast<bool>(transaction.publish_pending_cache(
                       RevisionSlotId{0U},
                       Span<const RevisionDependency>{
                           stale_dependencies.data(),
                           stale_dependencies.size()},
                       PendingCacheStamp{501U})),
                   "cache publication accepts a multi-source dependency tuple");
  ++geometry_revision;
  passed &= expect(static_cast<bool>(transaction.bind_dependency(
                       RevisionDependency{geometry_source,
                                          geometry_revision})),
                   "geometry authority can advance without reallocating");
  const Status stale_cache_finish =
      transaction.collective_finish(MPI_COMM_SELF, Status{});
  passed &= expect(stale_cache_finish.code == StatusCode::invalid_plan &&
                       !transaction.committed() &&
                       layers.handle(StateRole::accepted_n) == stale_old_n &&
                       transaction.pending_cache_valid(RevisionSlotId{0U}) &&
                       transaction.pending_cache(RevisionSlotId{0U}) == 500U,
                   "a changed tuple dependency rolls back and preserves committed cache");

  passed &= expect(static_cast<bool>(transaction.begin(layers)) &&
                       static_cast<bool>(transaction.revise_trial(velocity)) &&
                       static_cast<bool>(transaction.revise_trial(pressure)),
                   "transaction remains reusable after stale-cache rollback");
  const std::array arbitrary_dependencies{
      RevisionDependency{AttemptTransaction::field_revision_source(velocity),
                         RevisionToken{1U}}};
  const Status arbitrary = transaction.publish_pending_cache(
      RevisionSlotId{0U},
      Span<const RevisionDependency>{arbitrary_dependencies.data(),
                                     arbitrary_dependencies.size()},
      PendingCacheStamp{600U});
  passed &= expect(arbitrary.code == StatusCode::invalid_plan,
                   "cache publication rejects an arbitrary nonzero token");
  const Status arbitrary_finish =
      transaction.collective_finish(MPI_COMM_SELF, Status{});
  passed &= expect(arbitrary_finish.code == StatusCode::invalid_plan &&
                       !transaction.committed(),
                   "an ignored arbitrary cache stamp remains latched");

  FieldId other_velocity{};
  FieldId other_pressure{};
  StateLayers other_layers;
  passed &= make_layers(other_layers, other_velocity, other_pressure);
  const std::size_t other_old_n =
      other_layers.handle(StateRole::accepted_n);
  const Status rebind = transaction.begin(other_layers);
  passed &= expect(rebind.code == StatusCode::invalid_plan,
                   "a transaction rejects a different state-layer identity");
  const Status rebind_finish =
      transaction.collective_finish(MPI_COMM_SELF, Status{});
  passed &= expect(rebind_finish.code == StatusCode::invalid_plan &&
                       !transaction.committed() &&
                       other_layers.handle(StateRole::accepted_n) ==
                           other_old_n,
                   "a rejected layer rebind can still finish and rotates nothing");
  return passed;
}

bool test_invalid_state_machine_calls() {
  FieldId velocity{};
  FieldId pressure{};
  StateLayers layers;
  bool passed = make_layers(layers, velocity, pressure);
  AttemptTransaction transaction;
  passed &= expect(static_cast<bool>(AttemptTransaction::create(
                       layers.field_count(), 1U, layers.field_count(),
                       transaction)),
                   "state-machine transaction preallocates");
  passed &= expect(
      transaction.revise_trial(velocity).code == StatusCode::invalid_plan &&
          transaction.publish_pending_cache(
              RevisionSlotId{0U}, Span<const RevisionDependency>{},
              PendingCacheStamp{1U}).code ==
              StatusCode::invalid_plan &&
          transaction.collective_finish(MPI_COMM_SELF, Status{}).code ==
              StatusCode::invalid_plan,
      "hot transaction methods reject calls before begin");
  passed &= expect(static_cast<bool>(transaction.begin(layers)),
                   "state-machine attempt begins");
  passed &= expect(transaction.begin(layers).code == StatusCode::invalid_plan,
                   "double begin is rejected and latches failure");
  const Status finish = transaction.collective_finish(MPI_COMM_SELF, Status{});
  passed &= expect(finish.code == StatusCode::invalid_plan &&
                       !transaction.committed() &&
                       transaction.lowest_failing_rank() == 0,
                   "double begin cannot be ignored to commit a partial attempt");
  passed &= expect(transaction.collective_finish(MPI_COMM_SELF, Status{}).code ==
                       StatusCode::invalid_plan,
                   "finish after finish is rejected locally");
  return passed;
}

bool test_prepare_then_commit_transaction_and_final_flux() {
  FieldId velocity{};
  FieldId pressure{};
  StateLayers layers;
  bool passed = make_layers(layers, velocity, pressure);
  AttemptTransaction transaction;
  FaceFluxStorage source_storage;
  FaceFluxStorage final_storage;
  FaceFluxView source;
  FinalFaceFluxAuthority authority;
  FinalFaceFluxWriter writer;
  constexpr Int3 cells{2, 1, 1};
  passed &= expect(
      static_cast<bool>(AttemptTransaction::create(
          layers.field_count(), 1U, layers.field_count(), transaction)) &&
          static_cast<bool>(FaceFluxStorage::allocate_workspace(
              cells, 1U, source_storage)) &&
          static_cast<bool>(source_storage.workspace_view(0U, 41U, source)) &&
          static_cast<bool>(FaceFluxStorage::allocate_final(cells,
                                                            final_storage)) &&
          static_cast<bool>(authority.claim(71U, 0U, transaction, writer)),
      "two-phase transaction and final-flux authority initialize");
  passed &= expect(final_storage.counters().replicas == 3U,
                   "final flux reserves active, previous, and pending replicas");
  const std::array source_faces{source.x, source.y, source.z};
  for (std::size_t axis = 0U; axis < source_faces.size(); ++axis) {
    const auto face = source_faces[axis];
    for (std::int32_t z = 0; z < face.extents.z; ++z)
      for (std::int32_t y = 0; y < face.extents.y; ++y)
        for (std::int32_t x = 0; x < face.extents.x; ++x)
          face.unchecked({x, y, z}) =
              10.0 * static_cast<double>(axis + 1U) + x + 0.1 * y + 0.01 * z;
  }
  passed &= expect(static_cast<bool>(writer.initialize_committed(
                       final_storage, hundun::v04::as_const(source))),
                   "final flux has an accepted baseline");
  ConstFaceFluxView fresh_active_flux;
  ConstFaceFluxView fresh_previous_flux;
  passed &= expect(
      static_cast<bool>(writer.committed(final_storage,
                                         fresh_active_flux)) &&
          static_cast<bool>(writer.committed_previous(
              final_storage, fresh_previous_flux)) &&
          fresh_active_flux.revision == fresh_previous_flux.revision &&
          fresh_active_flux.certificate ==
              fresh_previous_flux.certificate &&
          flux_bits(fresh_active_flux) == flux_bits(fresh_previous_flux),
      "fresh final flux publishes identical accepted and previous history");

  const auto publish_attempt = [&](double value) {
    if (!transaction.begin(layers) || !transaction.revise_trial(velocity) ||
        !transaction.revise_trial(pressure)) {
      return false;
    }
    PendingFaceFluxView pending;
    const std::array dependencies{depends_on(transaction, pressure)};
    if (!writer.begin_pending(transaction, final_storage, pending) ||
        !overwrite_pending_flux(pending, value)) {
      return false;
    }
    return static_cast<bool>(writer.publish_pending(
        Span<const RevisionDependency>{dependencies.data(),
                                       dependencies.size()},
        pending));
  };

  const std::size_t accept_n = layers.handle(StateRole::accepted_n);
  const std::size_t accept_nm1 =
      layers.handle(StateRole::accepted_n_minus_one);
  const std::size_t accept_trial = layers.handle(StateRole::trial);
  ConstFaceFluxView accepted_flux;
  passed &= expect(publish_attempt(1000.0) &&
                       static_cast<bool>(writer.committed(final_storage,
                                                          accepted_flux)),
                   "accept candidate and prior committed flux are available");
  const RevisionToken accepted_flux_revision = accepted_flux.revision;
  const auto accepted_flux_certificate = accepted_flux.certificate;
  PreparedAttemptFinish accept;
  const Status accept_prepared =
      transaction.collective_prepare(MPI_COMM_SELF, Status{}, accept);
  ConstFaceFluxView during_accept_prepare;
  passed &= expect(
      static_cast<bool>(accept_prepared) && accept.valid() &&
          accept.decision() == AttemptFinishDecision::accept &&
          static_cast<bool>(accept.outcome()) &&
          accept.lowest_failing_rank() == -1 && transaction.active() &&
          !transaction.finished() && !transaction.committed() &&
          layers.handle(StateRole::accepted_n) == accept_n &&
          layers.handle(StateRole::accepted_n_minus_one) == accept_nm1 &&
          layers.handle(StateRole::trial) == accept_trial &&
          !transaction.pending_cache_valid(RevisionSlotId{0U}) &&
          static_cast<bool>(writer.committed(final_storage,
                                             during_accept_prepare)) &&
          during_accept_prepare.revision == accepted_flux_revision &&
          during_accept_prepare.certificate == accepted_flux_certificate,
      "accept preparation changes no role, cache, or active final flux");
  transaction.commit_accept(accept);
  ConstFaceFluxView committed_flux;
  passed &= expect(!accept.valid() && transaction.finished() &&
                       transaction.committed() && !transaction.active() &&
                       layers.handle(StateRole::accepted_n) == accept_trial &&
                       layers.handle(StateRole::accepted_n_minus_one) == accept_n &&
                       layers.handle(StateRole::trial) == accept_nm1 &&
                       transaction.pending_cache_valid(RevisionSlotId{0U}) &&
                       static_cast<bool>(writer.committed(final_storage,
                                                          committed_flux)) &&
                       committed_flux.revision > accepted_flux_revision,
                   "accept commit alone rotates state, cache, and final flux");

  ConstFaceFluxView committed_previous;
  passed &= expect(static_cast<bool>(writer.committed_previous(
                       final_storage, committed_previous)),
                   "accepted final flux retains the previous time level");
  const std::vector<std::uint64_t> committed_flux_bits =
      flux_bits(committed_flux);
  const std::vector<std::uint64_t> committed_previous_bits =
      flux_bits(committed_previous);
  const RevisionToken committed_previous_revision =
      committed_previous.revision;
  const auto committed_flux_certificate = committed_flux.certificate;
  const auto committed_previous_certificate =
      committed_previous.certificate;

  const std::size_t reject_n = layers.handle(StateRole::accepted_n);
  const std::size_t reject_nm1 =
      layers.handle(StateRole::accepted_n_minus_one);
  const std::size_t reject_trial = layers.handle(StateRole::trial);
  const RevisionToken reject_flux_revision = committed_flux.revision;
  passed &= expect(publish_attempt(-1.0e200),
                   "reject candidate writes poison into pending flux");
  PreparedAttemptFinish reject;
  const Status reject_prepared = transaction.collective_prepare(
      MPI_COMM_SELF, Status{StatusCode::rejected_step, 601U}, reject);
  ConstFaceFluxView during_reject_prepare;
  passed &= expect(
      static_cast<bool>(reject_prepared) && reject.valid() &&
          reject.decision() == AttemptFinishDecision::reject &&
          reject.outcome().code == StatusCode::rejected_step &&
          reject.outcome().detail == 601U &&
          reject.lowest_failing_rank() == 0 && transaction.active() &&
          layers.handle(StateRole::accepted_n) == reject_n &&
          layers.handle(StateRole::accepted_n_minus_one) == reject_nm1 &&
          layers.handle(StateRole::trial) == reject_trial &&
          static_cast<bool>(writer.committed(final_storage,
                                             during_reject_prepare)) &&
          during_reject_prepare.revision == reject_flux_revision,
      "reject preparation leaves accepted roles and final flux unchanged");
  transaction.commit_reject(reject);
  ConstFaceFluxView after_reject;
  ConstFaceFluxView previous_after_reject;
  passed &= expect(!reject.valid() && transaction.finished() &&
                       !transaction.committed() && !transaction.active() &&
                       transaction.lowest_failing_rank() == 0 &&
                       layers.handle(StateRole::accepted_n) == reject_n &&
                       layers.handle(StateRole::accepted_n_minus_one) == reject_nm1 &&
                       layers.handle(StateRole::trial) == reject_trial &&
                       static_cast<bool>(writer.committed(final_storage,
                                                          after_reject)) &&
                       static_cast<bool>(writer.committed_previous(
                           final_storage, previous_after_reject)) &&
                       after_reject.revision == reject_flux_revision &&
                       previous_after_reject.revision ==
                           committed_previous_revision &&
                       after_reject.certificate ==
                           committed_flux_certificate &&
                       previous_after_reject.certificate ==
                           committed_previous_certificate &&
                       flux_bits(after_reject) == committed_flux_bits &&
                       flux_bits(previous_after_reject) ==
                           committed_previous_bits,
                   "fresh reject preserves active and previous final-flux "
                   "bytes, revisions, and certificates bitwise");

  passed &= expect(publish_attempt(-2.0e200),
                   "abort candidate publishes pending flux");
  PreparedAttemptFinish abort;
  passed &= expect(static_cast<bool>(transaction.collective_prepare(
                       MPI_COMM_SELF, Status{}, abort)) &&
                       abort.decision() == AttemptFinishDecision::accept,
                   "an otherwise acceptable attempt can be prepared");
  transaction.commit_reject(abort);
  passed &= expect(!abort.valid() && !transaction.committed() &&
                       layers.handle(StateRole::accepted_n) == reject_n &&
                       writer.committed(final_storage, after_reject) &&
                       after_reject.revision == reject_flux_revision,
                   "a prepared accept can be aborted when a peer module cannot commit");
  return passed;
}

bool test_final_flux_publish_preflight_is_read_only_and_allocation_free() {
  FieldId velocity{};
  FieldId pressure{};
  StateLayers layers;
  bool passed = make_layers(layers, velocity, pressure);
  AttemptTransaction transaction;
  FaceFluxStorage source_storage;
  FaceFluxStorage final_storage;
  FaceFluxView source;
  FinalFaceFluxAuthority authority;
  FinalFaceFluxWriter writer;
  constexpr Int3 cells{2, 1, 1};
  passed &= expect(
      static_cast<bool>(AttemptTransaction::create(
          layers.field_count(), 1U, layers.field_count(), transaction)) &&
          static_cast<bool>(FaceFluxStorage::allocate_workspace(
              cells, 1U, source_storage)) &&
          static_cast<bool>(source_storage.workspace_view(0U, 141U, source)) &&
          static_cast<bool>(FaceFluxStorage::allocate_final(cells,
                                                            final_storage)) &&
          static_cast<bool>(authority.claim(171U, 0U, transaction, writer)),
      "final-flux preflight fixture initializes");
  const std::array source_faces{source.x, source.y, source.z};
  for (std::size_t axis = 0U; axis < source_faces.size(); ++axis) {
    const auto face = source_faces[axis];
    for (std::int32_t z = 0; z < face.extents.z; ++z)
      for (std::int32_t y = 0; y < face.extents.y; ++y)
        for (std::int32_t x = 0; x < face.extents.x; ++x)
          face.unchecked({x, y, z}) =
              10.0 * static_cast<double>(axis + 1U) + x + 0.1 * y;
  }
  passed &= expect(static_cast<bool>(writer.initialize_committed(
                       final_storage, hundun::v04::as_const(source))),
                   "final-flux preflight fixture initializes baseline");
  ConstFaceFluxView baseline;
  passed &= expect(static_cast<bool>(writer.committed(final_storage, baseline)),
                   "final-flux preflight baseline is observable");
  const auto baseline_bits = flux_bits(baseline);
  const RevisionToken baseline_revision = baseline.revision;
  const auto storage_counters = final_storage.counters();

  passed &= expect(static_cast<bool>(transaction.begin(layers)) &&
                       static_cast<bool>(transaction.revise_trial(velocity)) &&
                       static_cast<bool>(transaction.revise_trial(pressure)),
                   "final-flux preflight attempt begins");
  PendingFaceFluxView pending;
  passed &= expect(static_cast<bool>(writer.begin_pending(
                       transaction, final_storage, pending)) &&
                       overwrite_pending_flux(pending, 1700.0),
                   "final-flux preflight candidate owns writable pending flux");
  const RevisionToken pending_revision = pending.revision();
  const RevisionDependency correct = depends_on(transaction, pressure);
  const std::array stale_dependencies{RevisionDependency{
      correct.source, static_cast<RevisionToken>(correct.revision + 1U)}};
  const std::array duplicate_dependencies{correct, correct};

  Status stale_preflight;
  Status duplicate_preflight;
  Status correct_preflight;
  Status published;
  std::size_t hot_allocations = std::numeric_limits<std::size_t>::max();
  {
    allocation_observer::Guard guard;
    stale_preflight = writer.preflight_publish_pending(
        Span<const RevisionDependency>{stale_dependencies.data(),
                                       stale_dependencies.size()},
        pending);
    duplicate_preflight = writer.preflight_publish_pending(
        Span<const RevisionDependency>{duplicate_dependencies.data(),
                                       duplicate_dependencies.size()},
        pending);
    correct_preflight = writer.preflight_publish_pending(
        Span<const RevisionDependency>{&correct, 1U}, pending);
    published = writer.publish_pending(
        Span<const RevisionDependency>{&correct, 1U}, pending);
    hot_allocations =
        allocation_observer::count.load(std::memory_order_relaxed);
  }
  ConstFaceFluxView after_preflight;
  passed &= expect(
      stale_preflight.code == StatusCode::invalid_plan &&
          duplicate_preflight.code == StatusCode::invalid_plan &&
          static_cast<bool>(correct_preflight) && static_cast<bool>(published) &&
          !pending.valid() && transaction.active(),
      "invalid dependency preflights are retryable and the corrected publish "
      "succeeds in the same attempt");
  passed &= expect(
      static_cast<bool>(writer.committed(final_storage, after_preflight)) &&
          after_preflight.revision == baseline_revision &&
          flux_bits(after_preflight) == baseline_bits,
      "preflight and pending publication write no committed final-flux byte");
  const auto counters_after = final_storage.counters();
  passed &= expect(
      hot_allocations == 0U &&
          counters_after.aligned_payload_allocations ==
              storage_counters.aligned_payload_allocations &&
          counters_after.aligned_payload_bytes ==
              storage_counters.aligned_payload_bytes &&
          counters_after.replicas == storage_counters.replicas &&
          counters_after.directional_blocks ==
              storage_counters.directional_blocks,
      "preflight and normal publish allocate no storage or C++ resource");
  passed &= expect(pending_revision == baseline_revision + 1U,
                   "preflight preserves the issued pending revision until publish");

  PreparedAttemptFinish accepted;
  passed &= expect(static_cast<bool>(transaction.collective_prepare(
                       MPI_COMM_SELF, Status{}, accepted)) &&
                       accepted.decision() == AttemptFinishDecision::accept,
                   "corrected final-flux publication prepares accept");
  transaction.commit_accept(accepted);
  ConstFaceFluxView committed;
  passed &= expect(static_cast<bool>(writer.committed(final_storage, committed)) &&
                       committed.revision == pending_revision,
                   "corrected final-flux publication commits the pending revision");
  const auto committed_bits = flux_bits(committed);
  const auto committed_certificate = committed.certificate;
  passed &= expect(static_cast<bool>(transaction.begin(layers)) &&
                       static_cast<bool>(transaction.revise_trial(velocity)) &&
                       static_cast<bool>(transaction.revise_trial(pressure)),
                   "direct-publish fail-closed attempt begins");
  PendingFaceFluxView rejected_pending;
  passed &= expect(static_cast<bool>(writer.begin_pending(
                       transaction, final_storage, rejected_pending)) &&
                       overwrite_pending_flux(rejected_pending, -1700.0),
                   "direct-publish fail-closed candidate owns pending flux");
  const RevisionDependency current = depends_on(transaction, pressure);
  const std::array invalid_direct_dependencies{RevisionDependency{
      current.source, static_cast<RevisionToken>(current.revision + 1U)}};
  const Status invalid_direct = writer.publish_pending(
      Span<const RevisionDependency>{invalid_direct_dependencies.data(),
                                     invalid_direct_dependencies.size()},
      rejected_pending);
  const std::array current_dependencies{current};
  const Status ignored_retry = writer.publish_pending(
      Span<const RevisionDependency>{current_dependencies.data(),
                                     current_dependencies.size()},
      rejected_pending);
  const Status rejected =
      transaction.collective_finish(MPI_COMM_SELF, Status{});
  ConstFaceFluxView after_direct_reject;
  passed &= expect(
      invalid_direct.code == StatusCode::invalid_plan &&
          ignored_retry.code == StatusCode::invalid_plan &&
          rejected.code == StatusCode::invalid_plan &&
          !transaction.committed() && !rejected_pending.valid() &&
          static_cast<bool>(writer.committed(final_storage,
                                             after_direct_reject)) &&
          after_direct_reject.revision == pending_revision &&
          after_direct_reject.certificate == committed_certificate &&
          flux_bits(after_direct_reject) == committed_bits,
      "direct invalid publication remains fail-closed and rolls back without "
      "changing committed final flux");
  return passed;
}

bool test_restored_final_flux_starts_with_restart_lineage() {
  FieldId velocity{};
  FieldId pressure{};
  StateLayers layers;
  bool passed = make_layers(layers, velocity, pressure);
  AttemptTransaction transaction;
  FaceFluxStorage source_storage;
  FaceFluxStorage final_storage;
  FaceFluxView source;
  FinalFaceFluxAuthority authority;
  FinalFaceFluxWriter writer;
  constexpr Int3 cells{2, 2, 1};
  passed &= expect(
      static_cast<bool>(AttemptTransaction::create(
          layers.field_count(), 1U, layers.field_count(), transaction)) &&
          static_cast<bool>(FaceFluxStorage::allocate_workspace(
              cells, 1U, source_storage)) &&
          static_cast<bool>(source_storage.workspace_view(0U, 91U, source)) &&
          static_cast<bool>(FaceFluxStorage::allocate_final(cells,
                                                            final_storage)) &&
          static_cast<bool>(authority.claim(81U, 0U, transaction, writer)),
      "restart final-flux authority initializes without a fresh image");
  const std::array source_faces{source.x, source.y, source.z};
  for (std::size_t axis = 0U; axis < source_faces.size(); ++axis) {
    const auto face = source_faces[axis];
    for (std::int32_t z = 0; z < face.extents.z; ++z)
      for (std::int32_t y = 0; y < face.extents.y; ++y)
        for (std::int32_t x = 0; x < face.extents.x; ++x)
          face.unchecked({x, y, z}) =
              100.0 * static_cast<double>(axis + 1U) +
              10.0 * x + y + 0.01 * z;
  }
  const double finite_source = source.x.unchecked({0, 0, 0});
  source.x.unchecked({0, 0, 0}) =
      std::numeric_limits<double>::quiet_NaN();
  ConstFaceFluxView unavailable;
  passed &= expect(
      writer.initialize_restored(final_storage, as_const(source)).code ==
              StatusCode::numerical_failure &&
          writer.committed(final_storage, unavailable).code ==
              StatusCode::invalid_plan,
      "restart flux rejects a non-finite image before publishing lineage");
  source.x.unchecked({0, 0, 0}) = finite_source;
  passed &= expect(static_cast<bool>(writer.initialize_restored(
                       final_storage, as_const(source))),
                   "restart flux initializes directly from restored bytes");

  ConstFaceFluxView active;
  ConstFaceFluxView previous;
  passed &= expect(static_cast<bool>(writer.committed(final_storage, active)) &&
                       static_cast<bool>(writer.committed_previous(
                           final_storage, previous)) &&
                       active.revision == 2U && previous.revision == 2U &&
                       active.certificate.valid() &&
                       previous.certificate.valid(),
                   "restored active and BE history share revision two");
  const std::array active_faces{active.x, active.y, active.z};
  const std::array previous_faces{previous.x, previous.y, previous.z};
  for (std::size_t axis = 0U; axis < source_faces.size(); ++axis) {
    const auto expected = as_const(source_faces[axis]);
    for (std::int32_t z = 0; z < expected.extents.z; ++z)
      for (std::int32_t y = 0; y < expected.extents.y; ++y)
        for (std::int32_t x = 0; x < expected.extents.x; ++x) {
          const Int3 face{x, y, z};
          passed &= expect(active_faces[axis].unchecked(face) ==
                                   expected.unchecked(face) &&
                               previous_faces[axis].unchecked(face) ==
                                   expected.unchecked(face),
                           "restored active/history flux bytes are exact");
        }
  }

  passed &= expect(
      writer.initialize_restored(final_storage, as_const(source)).code ==
          StatusCode::invalid_plan,
      "restart flux lineage is initialized exactly once");
  passed &= expect(static_cast<bool>(transaction.begin(layers)) &&
                       static_cast<bool>(transaction.revise_trial(velocity)) &&
                       static_cast<bool>(transaction.revise_trial(pressure)),
                   "post-restart attempt begins");
  PendingFaceFluxView pending;
  passed &= expect(static_cast<bool>(writer.begin_pending(
                       transaction, final_storage, pending)) &&
                       pending.revision() == 3U,
                   "first post-restart pending flux advances revision two");
  passed &= expect(overwrite_pending_flux(pending, 2000.0),
                   "first post-restart pending bytes are writable");
  const std::array dependencies{depends_on(transaction, pressure)};
  passed &= expect(static_cast<bool>(writer.publish_pending(
                       Span<const RevisionDependency>{dependencies.data(),
                                                      dependencies.size()},
                       pending)),
                   "first post-restart pending flux publishes");
  PreparedAttemptFinish accepted;
  passed &= expect(static_cast<bool>(transaction.collective_prepare(
                       MPI_COMM_SELF, Status{}, accepted)) &&
                       accepted.decision() == AttemptFinishDecision::accept,
                   "first post-restart pending flux prepares an accept");
  transaction.commit_accept(accepted);

  ConstFaceFluxView accepted_active;
  ConstFaceFluxView accepted_previous;
  passed &= expect(
      static_cast<bool>(writer.committed(final_storage, accepted_active)) &&
          static_cast<bool>(writer.committed_previous(final_storage,
                                                      accepted_previous)) &&
          accepted_active.revision == 3U &&
          accepted_previous.revision == 2U,
      "post-restart accept establishes distinct active and previous lineage");
  const std::vector<std::uint64_t> accepted_active_bits =
      flux_bits(accepted_active);
  const std::vector<std::uint64_t> accepted_previous_bits =
      flux_bits(accepted_previous);
  const auto accepted_active_certificate = accepted_active.certificate;
  const auto accepted_previous_certificate = accepted_previous.certificate;

  passed &= expect(static_cast<bool>(transaction.begin(layers)) &&
                       static_cast<bool>(transaction.revise_trial(velocity)) &&
                       static_cast<bool>(transaction.revise_trial(pressure)),
                   "second post-restart attempt begins");
  PendingFaceFluxView poison;
  passed &= expect(static_cast<bool>(writer.begin_pending(
                       transaction, final_storage, poison)) &&
                       poison.revision() == 4U,
                   "second post-restart pending flux advances revision three");
  passed &= expect(overwrite_pending_flux(poison, -3.0e200),
                   "second post-restart pending bytes accept poison");
  const std::array poison_dependencies{depends_on(transaction, pressure)};
  passed &= expect(static_cast<bool>(writer.publish_pending(
                       Span<const RevisionDependency>{
                           poison_dependencies.data(),
                           poison_dependencies.size()},
                       poison)),
                   "second post-restart attempt publishes poison only to "
                   "pending flux");
  PreparedAttemptFinish rejected;
  passed &= expect(static_cast<bool>(transaction.collective_prepare(
                       MPI_COMM_SELF,
                       Status{StatusCode::rejected_step, 602U}, rejected)) &&
                       rejected.decision() == AttemptFinishDecision::reject,
                   "post-restart poison prepares a reject");
  transaction.commit_reject(rejected);
  ConstFaceFluxView after_reject;
  ConstFaceFluxView previous_after_reject;
  passed &= expect(!poison.valid() &&
                       static_cast<bool>(writer.committed(final_storage,
                                                          after_reject)) &&
                       static_cast<bool>(writer.committed_previous(
                           final_storage, previous_after_reject)) &&
                       after_reject.revision == 3U &&
                       previous_after_reject.revision == 2U &&
                       after_reject.certificate ==
                           accepted_active_certificate &&
                       previous_after_reject.certificate ==
                           accepted_previous_certificate &&
                       flux_bits(after_reject) == accepted_active_bits &&
                       flux_bits(previous_after_reject) ==
                           accepted_previous_bits,
                   "restart reject preserves active and previous final-flux "
                   "bytes, revisions, and certificates bitwise");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  static_assert(!std::is_move_constructible_v<AttemptTransaction> &&
                    !std::is_move_assignable_v<AttemptTransaction>,
                "an active transaction cannot be moved to a second owner");
  static_assert(noexcept(std::declval<AttemptTransaction&>().begin(
      std::declval<StateLayers&>())));
  static_assert(noexcept(std::declval<AttemptTransaction&>().revise_trial(
      std::declval<FieldId>())));
  static_assert(noexcept(
      std::declval<AttemptTransaction&>().publish_pending_cache(
          std::declval<RevisionSlotId>(),
          std::declval<Span<const RevisionDependency>>(),
          std::declval<PendingCacheStamp>())));
  static_assert(noexcept(
      std::declval<const FinalFaceFluxWriter&>().preflight_publish_pending(
          std::declval<Span<const RevisionDependency>>(),
          std::declval<const PendingFaceFluxView&>())));
  static_assert(noexcept(std::declval<AttemptTransaction&>().collective_finish(
      std::declval<MPI_Comm>(), std::declval<Status>())));
  static_assert(!std::is_copy_constructible_v<PreparedAttemptFinish> &&
                !std::is_move_constructible_v<PreparedAttemptFinish>);
  static_assert(noexcept(std::declval<AttemptTransaction&>().collective_prepare(
      std::declval<MPI_Comm>(), std::declval<Status>(),
      std::declval<PreparedAttemptFinish&>())));
  static_assert(noexcept(std::declval<AttemptTransaction&>().commit_accept(
      std::declval<PreparedAttemptFinish&>())));
  static_assert(noexcept(std::declval<AttemptTransaction&>().commit_reject(
      std::declval<PreparedAttemptFinish&>())));

  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  bool passed = test_commit_rollback_retry_and_hot_allocations();
  passed &= test_incomplete_write_set_rolls_back();
  passed &= test_collective_prepare_failure_discards_attempt();
  passed &= test_ignored_hot_error_is_latched();
  passed &= test_cache_stamp_and_layer_binding();
  passed &= test_invalid_state_machine_calls();
  passed &= test_prepare_then_commit_transaction_and_final_flux();
  passed &= test_final_flux_publish_preflight_is_read_only_and_allocation_free();
  passed &= test_restored_final_flux_starts_with_restart_lineage();
  MPI_Finalize();
  return passed ? 0 : 1;
}

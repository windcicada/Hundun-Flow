// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_execution.hpp"

#include "core_arena_detail.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <new>
#include <utility>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kTransactionState = 1U;
constexpr std::uint32_t kTransactionCacheSlot = 2U;
constexpr std::uint32_t kTransactionCacheRevision = 3U;
constexpr std::uint32_t kTransactionMpi = 4U;
constexpr std::uint32_t kTransactionIncompleteState = 5U;

bool valid_status_code(StatusCode code) noexcept {
  switch (code) {
    case StatusCode::ok:
    case StatusCode::invalid_case:
    case StatusCode::invalid_plan:
    case StatusCode::allocation_failure:
    case StatusCode::mpi_failure:
    case StatusCode::numerical_failure:
    case StatusCode::rejected_step:
    case StatusCode::io_failure:
      return true;
  }
  return false;
}

std::uint64_t failing_rank_key(int rank, bool failed) noexcept {
  return failed ? static_cast<std::uint64_t>(rank)
                : std::numeric_limits<std::uint64_t>::max();
}

}  // namespace

Status AttemptTransaction::create(std::size_t field_capacity,
                                  std::size_t revision_slot_capacity,
                                  std::size_t revision_source_capacity,
                                  AttemptTransaction& out) {
  if (field_capacity == 0U ||
      field_capacity >
          static_cast<std::size_t>(std::numeric_limits<FieldId>::max()) + 1U ||
      revision_slot_capacity >
      static_cast<std::size_t>(std::numeric_limits<RevisionSlotId>::max()) +
          1U ||
      revision_source_capacity < field_capacity ||
      revision_source_capacity >
          static_cast<std::size_t>(
              std::numeric_limits<RevisionSourceId>::max())) {
    return {StatusCode::invalid_plan, kTransactionCacheSlot};
  }
  try {
    std::vector<RevisionToken> active(revision_slot_capacity,
                                      RevisionToken{0U});
    std::vector<RevisionToken> pending(revision_slot_capacity,
                                       RevisionToken{0U});
    std::vector<RevisionToken> pending_state(revision_slot_capacity,
                                             RevisionToken{0U});
    std::vector<RevisionSourceId> pending_sources(revision_slot_capacity);
    std::vector<std::size_t> pending_counts(revision_slot_capacity, 0U);
    std::size_t dependency_slots = 0U;
    if (revision_slot_capacity != 0U &&
        revision_source_capacity > std::numeric_limits<std::size_t>::max() /
                             revision_slot_capacity) {
      return {StatusCode::invalid_plan, kTransactionCacheSlot};
    }
    dependency_slots = revision_slot_capacity * revision_source_capacity;
    pending_state.assign(dependency_slots, RevisionToken{0U});
    pending_sources.assign(dependency_slots, RevisionSourceId{0U});
    std::vector<RevisionDependency> dependency_catalog(
        revision_source_capacity);
    std::vector<RevisionToken> commit_buffer(revision_slot_capacity,
                                             RevisionToken{0U});
    std::vector<std::uint8_t> revised(field_capacity, std::uint8_t{0U});
    if (out.active_) {
      return {StatusCode::invalid_plan, kTransactionState};
    }
    out.active_caches_.swap(active);
    out.pending_caches_.swap(pending);
    out.pending_cache_dependency_revisions_.swap(pending_state);
    out.pending_cache_dependency_sources_.swap(pending_sources);
    out.pending_cache_dependency_counts_.swap(pending_counts);
    out.dependency_catalog_.swap(dependency_catalog);
    out.dependency_count_ = 0U;
    out.dependency_capacity_per_slot_ = revision_source_capacity;
    out.cache_commit_buffer_.swap(commit_buffer);
    out.revised_fields_.swap(revised);
    out.layers_ = nullptr;
    out.bound_layers_identity_ = 0U;
    out.active_ = false;
    out.committed_ = false;
    out.lowest_failing_rank_ = -1;
    out.attempt_status_ = {};
    out.attempt_identity_ = 0U;
    out.finished_ = false;
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, 0U};
  } catch (...) {
    return {StatusCode::invalid_plan, kTransactionCacheSlot};
  }
}

Status AttemptTransaction::begin(StateLayers& layers) noexcept {
  if (active_) {
    attempt_status_ = {StatusCode::invalid_plan, kTransactionState};
    return attempt_status_;
  }
  const std::uint64_t candidate_attempt_identity = detail::issue_identity();
  if (candidate_attempt_identity == 0U) {
    attempt_status_ = {StatusCode::invalid_plan, kTransactionState};
    return attempt_status_;
  }
  for (std::size_t slot = 0U; slot < pending_caches_.size(); ++slot) {
    pending_caches_[slot] = RevisionToken{0U};
    pending_cache_dependency_counts_[slot] = 0U;
  }
  for (std::size_t field = 0U; field < revised_fields_.size(); ++field) {
    revised_fields_[field] = 0U;
  }
  dependency_count_ = 0U;
  layers_ = &layers;
  active_ = true;
  finished_ = false;
  committed_ = false;
  attempt_identity_ = candidate_attempt_identity;
  lowest_failing_rank_ = -1;
  attempt_status_ = {};
  if (revised_fields_.size() != layers.field_count_) {
    attempt_status_ = {StatusCode::invalid_plan, kTransactionState};
    return attempt_status_;
  }
  const StorageIdentity layers_identity = layers.storage_identity();
  if (layers_identity == 0U ||
      (bound_layers_identity_ != 0U &&
       bound_layers_identity_ != layers_identity)) {
    attempt_status_ = {StatusCode::invalid_plan, kTransactionState};
    return attempt_status_;
  }
  RevisionToken required = RevisionToken{0U};
  for (std::size_t field = 0U; field < layers.field_count_; ++field) {
    if (layers.state_fields_[field] != 0U &&
        required == std::numeric_limits<RevisionToken>::max()) {
      attempt_status_ = {StatusCode::invalid_plan, kTransactionState};
      return attempt_status_;
    }
    required += layers.state_fields_[field] != 0U ? RevisionToken{1U}
                                                  : RevisionToken{0U};
  }
  if (required == 0U || layers.next_revision_ == 0U ||
      required > std::numeric_limits<RevisionToken>::max() -
                     layers.next_revision_) {
    attempt_status_ = {StatusCode::invalid_plan, kTransactionState};
    return attempt_status_;
  }
  bound_layers_identity_ = layers_identity;
  const std::size_t trial_replica = layers.handle(StateRole::trial);
  for (std::size_t field = 0U; field < layers.field_count_; ++field) {
    if (layers.state_fields_[field] != 0U) {
      const Status revised =
          layers.issue_revision(trial_replica, static_cast<FieldId>(field));
      if (!revised) {
        attempt_status_ = revised;
        return attempt_status_;
      }
    }
  }
  for (std::size_t field = 0U; field < layers.field_count_; ++field) {
    if (layers.state_fields_[field] != 0U) {
      const Status bound = bind_dependency(RevisionDependency{
          field_revision_source(static_cast<FieldId>(field)),
          layers.revision(StateRole::trial, static_cast<FieldId>(field))});
      if (!bound) {
        return bound;
      }
    }
  }
  return {};
}

Status AttemptTransaction::revise_trial(FieldId field) noexcept {
  if (!active_ || layers_ == nullptr) {
    return {StatusCode::invalid_plan, kTransactionState};
  }
  if (attempt_status_.code != StatusCode::ok) {
    return attempt_status_;
  }
  const std::size_t index = static_cast<std::size_t>(field);
  if (index >= revised_fields_.size() ||
      layers_->state_fields_[index] == 0U) {
    attempt_status_ = {StatusCode::invalid_plan, kTransactionState};
    return attempt_status_;
  }
  const Status revised =
      layers_->issue_revision(layers_->handle(StateRole::trial), field);
  if (!revised) {
    attempt_status_ = revised;
    return revised;
  }
  revised_fields_[index] = 1U;
  return bind_dependency(RevisionDependency{
      field_revision_source(field),
      layers_->revision(StateRole::trial, field)});
}

RevisionToken AttemptTransaction::trial_revision(FieldId field) const noexcept {
  if (!active_ || layers_ == nullptr ||
      attempt_status_.code != StatusCode::ok) {
    return RevisionToken{0U};
  }
  return layers_->revision(StateRole::trial, field);
}

Status AttemptTransaction::bind_dependency(
    RevisionDependency dependency) noexcept {
  if (!active_ || layers_ == nullptr) {
    return {StatusCode::invalid_plan, kTransactionState};
  }
  if (attempt_status_.code != StatusCode::ok) {
    return attempt_status_;
  }
  if (dependency.source == 0U || dependency.revision == 0U) {
    attempt_status_ = {StatusCode::invalid_plan,
                       kTransactionCacheRevision};
    return attempt_status_;
  }
  for (std::size_t index = 0U; index < dependency_count_; ++index) {
    RevisionDependency& bound = dependency_catalog_[index];
    if (bound.source == dependency.source) {
      bound.revision = dependency.revision;
      return {};
    }
  }
  if (dependency_count_ >= dependency_catalog_.size()) {
    attempt_status_ = {StatusCode::invalid_plan,
                       kTransactionCacheRevision};
    return attempt_status_;
  }
  dependency_catalog_[dependency_count_] = dependency;
  ++dependency_count_;
  return {};
}

Status AttemptTransaction::publish_pending_cache(
    RevisionSlotId slot, Span<const RevisionDependency> dependencies,
    PendingCacheStamp stamp) noexcept {
  const std::size_t index = static_cast<std::size_t>(slot);
  if (!active_ || layers_ == nullptr) {
    return {StatusCode::invalid_plan, kTransactionState};
  }
  if (attempt_status_.code != StatusCode::ok) {
    return attempt_status_;
  }
  if (index >= pending_caches_.size()) {
    attempt_status_ = {StatusCode::invalid_plan, kTransactionCacheSlot};
    return attempt_status_;
  }
  if (stamp.cache_revision == 0U || dependencies.data == nullptr ||
      dependencies.size == 0U ||
      dependencies.size > dependency_capacity_per_slot_ ||
      pending_caches_[index] != 0U) {
    attempt_status_ = {StatusCode::invalid_plan, kTransactionCacheRevision};
    return attempt_status_;
  }
  const std::size_t base = index * dependency_capacity_per_slot_;
  for (std::size_t dependency = 0U; dependency < dependencies.size;
       ++dependency) {
    const RevisionDependency& certificate = dependencies.data[dependency];
    bool matches_authority = false;
    for (std::size_t bound_index = 0U; bound_index < dependency_count_;
         ++bound_index) {
      const RevisionDependency& bound = dependency_catalog_[bound_index];
      if (bound.source == certificate.source &&
          bound.revision == certificate.revision) {
        matches_authority = true;
        break;
      }
    }
    if (certificate.source == 0U || certificate.revision == 0U ||
        !matches_authority) {
      attempt_status_ = {StatusCode::invalid_plan,
                         kTransactionCacheRevision};
      return attempt_status_;
    }
    for (std::size_t prior = 0U; prior < dependency; ++prior) {
      if (pending_cache_dependency_sources_[base + prior] ==
          certificate.source) {
        attempt_status_ = {StatusCode::invalid_plan,
                           kTransactionCacheRevision};
        return attempt_status_;
      }
    }
    pending_cache_dependency_sources_[base + dependency] = certificate.source;
    pending_cache_dependency_revisions_[base + dependency] =
        certificate.revision;
  }
  pending_caches_[index] = stamp.cache_revision;
  pending_cache_dependency_counts_[index] = dependencies.size;
  return {};
}

Status AttemptTransaction::collective_finish(MPI_Comm communicator,
                                              Status local_status) noexcept {
  const bool locally_active = active_ && layers_ != nullptr;
  if (!locally_active || !valid_status_code(local_status.code)) {
    local_status = {StatusCode::invalid_plan, kTransactionState};
  }
  if (locally_active && attempt_status_.code != StatusCode::ok) {
    local_status = attempt_status_;
  } else if (locally_active && local_status.code == StatusCode::ok) {
    for (std::size_t field = 0U; field < revised_fields_.size(); ++field) {
      if (layers_->state_fields_[field] != 0U && revised_fields_[field] == 0U) {
        local_status = {StatusCode::invalid_plan,
                        kTransactionIncompleteState};
        break;
      }
    }
    if (local_status.code == StatusCode::ok) {
      for (std::size_t slot = 0U; slot < pending_caches_.size(); ++slot) {
        const std::size_t base = slot * dependency_capacity_per_slot_;
        for (std::size_t dependency = 0U;
             dependency < pending_cache_dependency_counts_[slot];
             ++dependency) {
          bool matches_authority = false;
          for (std::size_t bound_index = 0U;
               bound_index < dependency_count_; ++bound_index) {
            const RevisionDependency& bound =
                dependency_catalog_[bound_index];
            if (bound.source ==
                    pending_cache_dependency_sources_[base + dependency] &&
                bound.revision ==
                    pending_cache_dependency_revisions_[base + dependency]) {
              matches_authority = true;
              break;
            }
          }
          if (!matches_authority) {
            local_status = {StatusCode::invalid_plan,
                            kTransactionCacheRevision};
            break;
          }
        }
        if (local_status.code != StatusCode::ok) {
          break;
        }
      }
    }
  }
  if (communicator == MPI_COMM_NULL) {
    if (locally_active) {
      discard_attempt();
    }
    lowest_failing_rank_ = -1;
    return {StatusCode::invalid_plan, kTransactionMpi};
  }

  int rank = 0;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    discard_attempt();
    lowest_failing_rank_ = -1;
    return {StatusCode::mpi_failure, kTransactionMpi};
  }

  const std::uint64_t candidate =
      failing_rank_key(rank, local_status.code != StatusCode::ok);
  std::uint64_t selected = std::numeric_limits<std::uint64_t>::max();
  if (MPI_Allreduce(&candidate, &selected, 1, MPI_UINT64_T, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    discard_attempt();
    lowest_failing_rank_ = -1;
    return {StatusCode::mpi_failure, kTransactionMpi};
  }

  if (selected == std::numeric_limits<std::uint64_t>::max()) {
    for (std::size_t slot = 0U; slot < pending_caches_.size(); ++slot) {
      cache_commit_buffer_[slot] = pending_caches_[slot];
    }
    active_caches_.swap(cache_commit_buffer_);
    if (layers_ == nullptr) {
      lowest_failing_rank_ = -1;
      return {StatusCode::invalid_plan, kTransactionState};
    }
    layers_->rotate_commit();
    active_ = false;
    layers_ = nullptr;
    finished_ = true;
    committed_ = true;
    lowest_failing_rank_ = -1;
    return {};
  }
  if (selected >= static_cast<std::uint64_t>(size)) {
    discard_attempt();
    lowest_failing_rank_ = -1;
    return {StatusCode::mpi_failure, kTransactionMpi};
  }

  const int selected_rank = static_cast<int>(selected);
  std::array<std::uint64_t, 2U> wire{};
  if (rank == selected_rank) {
    wire[0] = static_cast<std::uint64_t>(local_status.code);
    wire[1] = static_cast<std::uint64_t>(local_status.detail);
  }
  if (MPI_Bcast(wire.data(), static_cast<int>(wire.size()), MPI_UINT64_T,
                selected_rank, communicator) != MPI_SUCCESS ||
      wire[0] > static_cast<std::uint64_t>(StatusCode::io_failure) ||
      wire[1] > std::numeric_limits<std::uint32_t>::max()) {
    discard_attempt();
    lowest_failing_rank_ = -1;
    return {StatusCode::mpi_failure, kTransactionMpi};
  }

  discard_attempt();
  lowest_failing_rank_ = selected_rank;
  return {static_cast<StatusCode>(wire[0]),
          static_cast<std::uint32_t>(wire[1])};
}

void AttemptTransaction::discard_attempt() noexcept {
  for (std::size_t slot = 0U; slot < pending_caches_.size(); ++slot) {
    pending_caches_[slot] = RevisionToken{0U};
    pending_cache_dependency_counts_[slot] = 0U;
  }
  active_ = false;
  layers_ = nullptr;
  finished_ = true;
  committed_ = false;
}

bool AttemptTransaction::pending_cache_valid(
    RevisionSlotId slot) const noexcept {
  const std::size_t index = static_cast<std::size_t>(slot);
  return index < active_caches_.size() && active_caches_[index] != 0U;
}

RevisionToken AttemptTransaction::pending_cache(
    RevisionSlotId slot) const noexcept {
  const std::size_t index = static_cast<std::size_t>(slot);
  return index < active_caches_.size() ? active_caches_[index]
                                      : RevisionToken{0U};
}

}  // namespace hundun::v04

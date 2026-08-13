// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_field.hpp"

#include <mpi.h>

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

namespace hundun::v04 {

using ArenaOwnerId = std::uint16_t;
using RevisionSlotId = std::uint16_t;
using StorageIdentity = std::uint64_t;
using RevisionDomainIdentity = std::uint64_t;
using RevisionSourceId = std::uint32_t;

enum class FieldLifetime : std::uint8_t {
  state_layer,
  pending_cache,
  persistent_workspace,
  step_scratch
};

struct FieldPlacement {
  ArenaOwnerId owner{};
};

struct ArenaFieldRequest {
  FieldId id{};
  Int3 interior{};
  FieldPlacement placement{};
  FieldLifetime lifetime{FieldLifetime::state_layer};
};

struct ArenaFieldLayout {
  FieldId id{};
  Int3 interior{};
  Int3 ghosts{};
  std::uint8_t components{};
  std::size_t stride_y{};
  std::size_t stride_z{};
  std::size_t component_stride{};
  std::size_t replicas{};
  std::size_t replica_stride_doubles{};
  ArenaOwnerId owner{};
  std::size_t owner_index{};
  std::size_t raw_offset_doubles{};
  std::size_t offset_doubles{};
  std::size_t span_doubles{};
  FieldLifetime lifetime{FieldLifetime::state_layer};
};

struct ArenaOwnerLayout {
  ArenaOwnerId owner{};
  std::size_t offset_doubles{};
  std::size_t span_doubles{};
};

class ArenaLayout {
 public:
  static Status compile(const FieldSchema& schema,
                        Span<const ArenaFieldRequest> requests,
                        ArenaLayout& out);

  std::size_t field_count() const noexcept { return fields_.size(); }
  std::size_t owner_count() const noexcept { return owners_.size(); }
  std::size_t total_doubles() const noexcept { return total_doubles_; }
  const ArenaFieldLayout* field(FieldId id) const noexcept;
  const ArenaOwnerLayout* owner(std::size_t index) const noexcept;

 private:
  friend class FieldStorage;
  friend class StateLayers;
  std::vector<ArenaFieldLayout> fields_;
  std::vector<ArenaOwnerLayout> owners_;
  std::size_t total_doubles_{};
};

class RevisionSet {
 public:
  RevisionSet() noexcept = default;
  RevisionSet(const RevisionSet&) = delete;
  RevisionSet& operator=(const RevisionSet&) = delete;
  RevisionSet(RevisionSet&&) noexcept = default;
  RevisionSet& operator=(RevisionSet&&) noexcept = default;

  static Status create(const FieldSchema& schema, RevisionSet& out);
  RevisionToken token(FieldId field) const noexcept;
  Status revise(FieldId field) noexcept;

 private:
  friend class FieldStorage;
  friend class StateLayers;
  friend class AttemptTransaction;
  std::vector<RevisionToken> tokens_;
  RevisionDomainIdentity identity_{};
};

template <class T>
struct BasicFieldView {
  static_assert(std::is_same_v<T, double> ||
                    std::is_same_v<T, const double>,
                "BasicFieldView supports FP64 fields only");

  T* base{};
  Int3 interior{};
  Int3 ghosts{};
  std::uint8_t components{};
  std::size_t stride_y{};
  std::size_t stride_z{};
  std::size_t component_stride{};
  std::size_t replica{};
  FieldId field{};
  RevisionToken revision{};
  StorageIdentity storage_identity{};
  RevisionDomainIdentity revision_domain{};

  T& unchecked(Int3 index, std::uint8_t component) const noexcept {
    const auto x = static_cast<std::ptrdiff_t>(index.x);
    const auto y = static_cast<std::ptrdiff_t>(index.y);
    const auto z = static_cast<std::ptrdiff_t>(index.z);
    const auto sy = static_cast<std::ptrdiff_t>(stride_y);
    const auto sz = static_cast<std::ptrdiff_t>(stride_z);
    const auto sc = static_cast<std::ptrdiff_t>(component_stride);
    const auto c = static_cast<std::ptrdiff_t>(component);
    return base[x + y * sy + z * sz + c * sc];
  }
};

using FieldView = BasicFieldView<double>;
using ConstFieldView = BasicFieldView<const double>;

ConstFieldView as_const(FieldView view) noexcept;

struct ExecutionCounters {
  std::uint64_t aligned_payload_allocations{};
  std::uint64_t aligned_payload_bytes{};
  std::uint64_t whole_field_copies{};
};

class FieldStorage {
 public:
  FieldStorage() noexcept = default;
  ~FieldStorage() noexcept;
  FieldStorage(const FieldStorage&) = delete;
  FieldStorage& operator=(const FieldStorage&) = delete;
  FieldStorage(FieldStorage&& other) noexcept;
  FieldStorage& operator=(FieldStorage&& other) noexcept;

  static Status allocate(const ArenaLayout& layout, FieldStorage& out);
  Status view(FieldId field, std::size_t replica,
              const RevisionSet& revisions, FieldView& out) noexcept;
  Status view(FieldId field, std::size_t replica,
              const RevisionSet& revisions, ConstFieldView& out) const noexcept;
  double* checked_ptr(const FieldView& view, Int3 index,
                      std::uint8_t component,
                      const RevisionSet& revisions) noexcept;
  const double* checked_ptr(const ConstFieldView& view, Int3 index,
                            std::uint8_t component,
                            const RevisionSet& revisions) const noexcept;
  ExecutionCounters counters() const noexcept { return counters_; }

 private:
  friend class StateLayers;
  struct OwnerAllocation {
    ArenaOwnerId owner{};
    double* data{};
    std::size_t doubles{};
  };

  void release() noexcept;
  const ArenaFieldLayout* field_layout(FieldId field) const noexcept;
  Status view(FieldId field, std::size_t replica, RevisionToken revision,
              RevisionDomainIdentity revision_domain,
              FieldView& out) noexcept;
  Status view(FieldId field, std::size_t replica, RevisionToken revision,
              RevisionDomainIdentity revision_domain,
              ConstFieldView& out) const noexcept;
  double* checked_ptr(const FieldView& view, Int3 index,
                      std::uint8_t component,
                      RevisionToken expected_revision,
                      RevisionDomainIdentity expected_domain) noexcept;
  const double* checked_ptr(const ConstFieldView& view, Int3 index,
                            std::uint8_t component,
                            RevisionToken expected_revision,
                            RevisionDomainIdentity expected_domain) const noexcept;
  template <class T>
  T* checked_ptr_impl(const BasicFieldView<T>& view, Int3 index,
                      std::uint8_t component,
                      RevisionToken expected_revision,
                      RevisionDomainIdentity expected_domain) const noexcept;

  std::vector<ArenaFieldLayout> fields_;
  std::vector<OwnerAllocation> owners_;
  ExecutionCounters counters_{};
  StorageIdentity identity_{};
};

enum class StateRole : std::uint8_t {
  accepted_n,
  accepted_n_minus_one,
  trial
};

struct PendingCacheStamp {
  RevisionToken cache_revision{};
};

struct RevisionDependency {
  RevisionSourceId source{};
  RevisionToken revision{};
};

class StateLayers {
 public:
  StateLayers() noexcept = default;
  StateLayers(const StateLayers&) = delete;
  StateLayers& operator=(const StateLayers&) = delete;
  StateLayers(StateLayers&& other) noexcept;
  StateLayers& operator=(StateLayers&& other) noexcept;

  static Status allocate(const ArenaLayout& layout, StateLayers& out);
  Status view(StateRole role, FieldId field, FieldView& out) noexcept;
  Status view(StateRole role, FieldId field, ConstFieldView& out) const noexcept;
  std::size_t handle(StateRole role) const noexcept;
  RevisionToken revision(StateRole role, FieldId field) const noexcept;
  RevisionToken state_revision(StateRole role) const noexcept;
  std::size_t field_count() const noexcept { return field_count_; }
  double* checked_ptr(StateRole role, const FieldView& view, Int3 index,
                      std::uint8_t component) noexcept;
  const double* checked_ptr(StateRole role, const ConstFieldView& view,
                            Int3 index,
                            std::uint8_t component) const noexcept;
  ExecutionCounters counters() const noexcept { return storage_.counters(); }

 private:
  friend class AttemptTransaction;
  static constexpr std::size_t kRoleCount = 3U;
  std::size_t role_index(StateRole role) const noexcept;
  StorageIdentity storage_identity() const noexcept {
    return storage_.identity_;
  }
  void rotate_commit() noexcept;
  Status issue_revision(std::size_t replica, FieldId field) noexcept;

  FieldStorage storage_;
  std::vector<RevisionToken> revisions_;
  RevisionToken state_revisions_[kRoleCount]{};
  std::vector<std::uint8_t> state_fields_;
  std::size_t field_count_{};
  std::size_t role_handles_[kRoleCount]{0U, 1U, 2U};
  RevisionToken next_revision_{1U};
};

class AttemptTransaction {
 public:
  AttemptTransaction() noexcept = default;
  AttemptTransaction(const AttemptTransaction&) = delete;
  AttemptTransaction& operator=(const AttemptTransaction&) = delete;
  AttemptTransaction(AttemptTransaction&&) = delete;
  AttemptTransaction& operator=(AttemptTransaction&&) = delete;

  static Status create(std::size_t field_capacity,
                       std::size_t revision_slot_capacity,
                       std::size_t revision_source_capacity,
                       AttemptTransaction& out);
  Status begin(StateLayers& layers) noexcept;
  Status revise_trial(FieldId field) noexcept;
  RevisionToken trial_revision(FieldId field) const noexcept;
  static RevisionSourceId field_revision_source(FieldId field) noexcept {
    return static_cast<RevisionSourceId>(field) + 1U;
  }
  Status bind_dependency(RevisionDependency dependency) noexcept;
  Status publish_pending_cache(RevisionSlotId slot,
                               Span<const RevisionDependency> dependencies,
                               PendingCacheStamp stamp) noexcept;
  Status collective_finish(MPI_Comm communicator,
                           Status local_status) noexcept;

  bool active() const noexcept { return active_; }
  bool finished() const noexcept { return finished_; }
  bool committed() const noexcept { return committed_; }
  std::uint64_t attempt_identity() const noexcept {
    return attempt_identity_;
  }
  int lowest_failing_rank() const noexcept { return lowest_failing_rank_; }
  bool pending_cache_valid(RevisionSlotId slot) const noexcept;
  RevisionToken pending_cache(RevisionSlotId slot) const noexcept;

 private:
  void discard_attempt() noexcept;
  StateLayers* layers_{};
  std::vector<RevisionToken> active_caches_;
  std::vector<RevisionToken> pending_caches_;
  std::vector<RevisionToken> pending_cache_dependency_revisions_;
  std::vector<RevisionSourceId> pending_cache_dependency_sources_;
  std::vector<std::size_t> pending_cache_dependency_counts_;
  std::vector<RevisionDependency> dependency_catalog_;
  std::size_t dependency_count_{};
  std::size_t dependency_capacity_per_slot_{};
  std::vector<RevisionToken> cache_commit_buffer_;
  std::vector<std::uint8_t> revised_fields_;
  StorageIdentity bound_layers_identity_{};
  Status attempt_status_{};
  std::uint64_t attempt_identity_{};
  bool active_{};
  bool finished_{};
  bool committed_{};
  int lowest_failing_rank_{-1};
};

}  // namespace hundun::v04

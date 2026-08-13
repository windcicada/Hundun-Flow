// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_execution.hpp"
#include "hundun/v04_field.hpp"

#include <mpi.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string_view>

namespace {

using hundun::v04::ArenaFieldRequest;
using hundun::v04::ArenaLayout;
using hundun::v04::AttemptTransaction;
using hundun::v04::FieldId;
using hundun::v04::FieldLifetime;
using hundun::v04::FieldPlacement;
using hundun::v04::FieldRegistry;
using hundun::v04::FieldSchema;
using hundun::v04::PendingCacheStamp;
using hundun::v04::RevisionDependency;
using hundun::v04::RevisionSlotId;
using hundun::v04::Span;
using hundun::v04::StateLayers;
using hundun::v04::StateRole;
using hundun::v04::Status;
using hundun::v04::StatusCode;

bool expect(bool condition, int rank, std::string_view description) {
  if (!condition) {
    std::cerr << "rank " << rank << " FAIL: " << description << '\n';
  }
  return condition;
}

bool all_true(bool local) {
  const int input = local ? 1 : 0;
  int output = 0;
  return MPI_Allreduce(&input, &output, 1, MPI_INT, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         output != 0;
}

std::uint64_t packed_status(Status status) {
  return (static_cast<std::uint64_t>(status.code) << 32U) |
         static_cast<std::uint64_t>(status.detail);
}

bool identical(std::uint64_t value) {
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  return MPI_Allreduce(&value, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         MPI_Allreduce(&value, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         minimum == maximum;
}

bool make_layers(StateLayers& layers, FieldId& field) {
  FieldRegistry registry;
  FieldSchema schema;
  if (!registry.declare_field("state", 1U, 1U, field) ||
      !registry.freeze(schema)) {
    return false;
  }
  const std::array requests{
      ArenaFieldRequest{field, hundun::v04::Int3{8, 4, 2},
                        FieldPlacement{0U}, FieldLifetime::state_layer},
  };
  ArenaLayout layout;
  return ArenaLayout::compile(
             schema,
             Span<const ArenaFieldRequest>{requests.data(), requests.size()},
             layout) &&
         StateLayers::allocate(layout, layers);
}

bool test_success(int rank) {
  StateLayers layers;
  FieldId field{};
  bool passed = expect(make_layers(layers, field), rank,
                       "success fixture allocates");
  AttemptTransaction transaction;
  passed &= expect(static_cast<bool>(AttemptTransaction::create(
                       layers.field_count(), 2U, layers.field_count(),
                       transaction)),
                   rank, "success transaction preallocates");
  const std::size_t old_n = layers.handle(StateRole::accepted_n);
  const std::size_t old_nm1 = layers.handle(StateRole::accepted_n_minus_one);
  const std::size_t old_trial = layers.handle(StateRole::trial);

  passed &= expect(static_cast<bool>(transaction.begin(layers)), rank,
                   "success attempt begins");
  passed &= expect(static_cast<bool>(transaction.revise_trial(field)), rank,
                   "success attempt revises its field");
  const std::array dependencies{
      RevisionDependency{AttemptTransaction::field_revision_source(field),
                         transaction.trial_revision(field)}};
  passed &= expect(static_cast<bool>(transaction.publish_pending_cache(
                       RevisionSlotId{0U},
                       Span<const RevisionDependency>{dependencies.data(),
                                                      dependencies.size()},
                       PendingCacheStamp{
                           static_cast<std::uint64_t>(1000 + rank)})),
                   rank, "success attempt records pending cache");
  const Status result = transaction.collective_finish(MPI_COMM_WORLD, Status{});

  passed &= expect(static_cast<bool>(result) && transaction.committed(), rank,
                   "all-success consensus commits everywhere");
  passed &= expect(transaction.lowest_failing_rank() == -1, rank,
                   "all-success exposes no failing rank everywhere");
  passed &= expect(identical(packed_status(result)), rank,
                   "all-success status is identical everywhere");
  passed &= expect(layers.handle(StateRole::accepted_n) == old_trial &&
                       layers.handle(StateRole::accepted_n_minus_one) == old_n &&
                       layers.handle(StateRole::trial) == old_nm1,
                   rank, "all-success rotates identical role handles everywhere");
  passed &= expect(transaction.pending_cache_valid(RevisionSlotId{0U}) &&
                       transaction.pending_cache(RevisionSlotId{0U}) ==
                           static_cast<std::uint64_t>(1000 + rank),
                   rank, "success publishes each rank's local pending cache");
  passed &= expect(layers.counters().whole_field_copies == 0U, rank,
                   "collective success copies no whole field");
  const std::uint64_t handle_signature =
      (static_cast<std::uint64_t>(layers.handle(StateRole::accepted_n)) << 16U) |
      (static_cast<std::uint64_t>(
           layers.handle(StateRole::accepted_n_minus_one))
       << 8U) |
      static_cast<std::uint64_t>(layers.handle(StateRole::trial));
  passed &= expect(identical(handle_signature), rank,
                   "all ranks finish with identical role assignment");
  return all_true(passed);
}

bool test_lowest_failing_rank(int rank, int size) {
  StateLayers layers;
  FieldId field{};
  bool passed = expect(make_layers(layers, field), rank,
                       "failure fixture allocates");
  AttemptTransaction transaction;
  passed &= expect(static_cast<bool>(AttemptTransaction::create(
                       layers.field_count(), 2U, layers.field_count(),
                       transaction)),
                   rank, "failure transaction preallocates");
  const std::size_t old_n = layers.handle(StateRole::accepted_n);
  const std::size_t old_nm1 = layers.handle(StateRole::accepted_n_minus_one);
  const std::size_t old_trial = layers.handle(StateRole::trial);

  passed &= expect(static_cast<bool>(transaction.begin(layers)), rank,
                   "failure attempt begins");
  passed &= expect(static_cast<bool>(transaction.revise_trial(field)), rank,
                   "failure attempt revises trial metadata");
  const std::array dependencies{
      RevisionDependency{AttemptTransaction::field_revision_source(field),
                         transaction.trial_revision(field)}};
  passed &= expect(static_cast<bool>(transaction.publish_pending_cache(
                       RevisionSlotId{1U},
                       Span<const RevisionDependency>{dependencies.data(),
                                                      dependencies.size()},
                       PendingCacheStamp{
                           static_cast<std::uint64_t>(2000 + rank)})),
                   rank, "failure attempt records pending cache");

  Status local{};
  Status expected{};
  if (size == 1) {
    local = Status{StatusCode::rejected_step, 701U};
    expected = local;
  } else {
    const int lowest_failing_rank = 0;
    const int higher_failing_rank = size - 1;
    if (rank == lowest_failing_rank) {
      local = Status{StatusCode::numerical_failure, 701U};
    } else if (rank == higher_failing_rank) {
      // This code/detail is numerically smaller on purpose. Consensus must
      // select by rank, not by StatusCode or a packed-status minimum.
      local = Status{StatusCode::invalid_plan, 3U};
    }
    expected = Status{StatusCode::numerical_failure, 701U};
  }

  const Status result = transaction.collective_finish(MPI_COMM_WORLD, local);
  passed &= expect(result.code == expected.code && result.detail == expected.detail,
                   rank, "the lowest failing rank's exact status wins");
  passed &= expect(identical(packed_status(result)), rank,
                   "failure status category/detail is identical everywhere");
  passed &= expect(!transaction.committed(), rank,
                   "any local failure rolls back every rank");
  passed &= expect(transaction.lowest_failing_rank() == 0, rank,
                   "one-rank and multi-rank failure select rank zero");
  passed &= expect(layers.handle(StateRole::accepted_n) == old_n &&
                       layers.handle(StateRole::accepted_n_minus_one) ==
                           old_nm1 &&
                       layers.handle(StateRole::trial) == old_trial,
                   rank, "failure preserves every accepted/trial role handle");
  passed &= expect(!transaction.pending_cache_valid(RevisionSlotId{1U}), rank,
                   "failure publishes no pending cache on any rank");
  passed &= expect(layers.counters().whole_field_copies == 0U, rank,
                   "collective rollback copies no whole field");

  const std::uint64_t handle_signature =
      (static_cast<std::uint64_t>(layers.handle(StateRole::accepted_n)) << 16U) |
      (static_cast<std::uint64_t>(
           layers.handle(StateRole::accepted_n_minus_one))
       << 8U) |
      static_cast<std::uint64_t>(layers.handle(StateRole::trial));
  passed &= expect(identical(handle_signature), rank,
                   "rollback role assignment is identical everywhere");
  return all_true(passed);
}

bool test_rank_local_begin_failure(int rank, int size) {
  StateLayers layers;
  FieldId field{};
  bool passed = expect(make_layers(layers, field), rank,
                       "prepare-failure fixture allocates");
  const int failing_rank = size == 1 ? 0 : 1;
  const std::size_t capacity =
      layers.field_count() + (rank == failing_rank ? 1U : 0U);
  AttemptTransaction transaction;
  passed &= expect(static_cast<bool>(AttemptTransaction::create(
                       capacity, 1U, capacity, transaction)),
                   rank, "prepare-failure transaction preallocates locally");
  const std::size_t old_n = layers.handle(StateRole::accepted_n);

  const Status begin = transaction.begin(layers);
  if (rank == failing_rank) {
    passed &= expect(begin.code == StatusCode::invalid_plan, rank,
                     "selected rank reports local begin failure");
  } else {
    passed &= expect(static_cast<bool>(begin) &&
                         static_cast<bool>(transaction.revise_trial(field)),
                     rank, "nonfailing rank completes local prepare work");
  }
  const Status result = transaction.collective_finish(MPI_COMM_WORLD, Status{});
  passed &= expect(result.code == StatusCode::invalid_plan &&
                       result.detail == 1U,
                   rank, "rank-local begin failure reaches collective consensus");
  passed &= expect(identical(packed_status(result)), rank,
                   "prepare-failure status is identical everywhere");
  passed &= expect(!transaction.committed() &&
                       transaction.lowest_failing_rank() == failing_rank &&
                       layers.handle(StateRole::accepted_n) == old_n,
                   rank, "prepare failure rolls back every rank without deadlock");
  return all_true(passed);
}

bool test_rank_local_inactive_finish(int rank, int size) {
  StateLayers layers;
  FieldId field{};
  bool passed = expect(make_layers(layers, field), rank,
                       "inactive-finish fixture allocates");
  AttemptTransaction transaction;
  passed &= expect(static_cast<bool>(AttemptTransaction::create(
                       layers.field_count(), 1U, layers.field_count(),
                       transaction)),
                   rank, "inactive-finish transaction preallocates");
  const int inactive_rank = size == 1 ? 0 : size - 1;
  if (rank != inactive_rank) {
    passed &= expect(static_cast<bool>(transaction.begin(layers)) &&
                         static_cast<bool>(transaction.revise_trial(field)),
                     rank, "active ranks prepare before protocol divergence");
  }
  const Status result = transaction.collective_finish(MPI_COMM_WORLD, Status{});
  passed &= expect(result.code == StatusCode::invalid_plan &&
                       result.detail == 1U,
                   rank, "an inactive rank still joins failure consensus");
  passed &= expect(identical(packed_status(result)), rank,
                   "inactive-rank failure status is identical everywhere");
  passed &= expect(!transaction.committed() &&
                       transaction.lowest_failing_rank() == inactive_rank,
                   rank, "inactive finish rolls back all participating ranks");
  return all_true(passed);
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  bool passed = test_success(rank);
  passed &= test_lowest_failing_rank(rank, size);
  passed &= test_rank_local_begin_failure(rank, size);
  passed &= test_rank_local_inactive_finish(rank, size);
  passed = all_true(passed);
  if (rank == 0 && !passed) {
    std::cerr << "v04 core transaction MPI test failed for " << size
              << " ranks\n";
  }
  MPI_Finalize();
  return passed ? 0 : 1;
}

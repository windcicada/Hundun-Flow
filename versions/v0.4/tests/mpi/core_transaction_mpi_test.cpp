// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_execution.hpp"
#include "hundun/v04_field.hpp"

#include <mpi.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string_view>
#include <vector>

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
using hundun::v04::FinalFaceFluxAuthority;
using hundun::v04::FinalFaceFluxWriter;
using hundun::v04::Int3;
using hundun::v04::PendingCacheStamp;
using hundun::v04::PendingFaceFluxView;
using hundun::v04::PreparedAttemptFinish;
using hundun::v04::RevisionDependency;
using hundun::v04::RevisionSlotId;
using hundun::v04::RevisionToken;
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

void fill_flux(FaceFluxView flux, double seed) {
  const std::array faces{flux.x, flux.y, flux.z};
  double value = seed;
  for (const auto face : faces) {
    for (std::int32_t z = 0; z < face.extents.z; ++z)
      for (std::int32_t y = 0; y < face.extents.y; ++y)
        for (std::int32_t x = 0; x < face.extents.x; ++x) {
          face.unchecked({x, y, z}) = value;
          value += 0.125;
        }
  }
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

bool test_final_flux_collective_rollback(int rank) {
  StateLayers layers;
  FieldId field{};
  bool passed = expect(make_layers(layers, field), rank,
                       "final-flux rollback fixture allocates state");
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
          static_cast<bool>(source_storage.workspace_view(0U, 51U, source)) &&
          static_cast<bool>(FaceFluxStorage::allocate_final(cells,
                                                            final_storage)) &&
          static_cast<bool>(authority.claim(91U, 0U, transaction, writer)),
      rank, "final-flux rollback fixture initializes authority");
  passed &= expect(final_storage.counters().replicas == 3U, rank,
                   "final flux has a rollback-safe third replica");
  fill_flux(source, 100.0 * static_cast<double>(rank + 1));
  passed &= expect(static_cast<bool>(writer.initialize_committed(
                       final_storage, hundun::v04::as_const(source))),
                   rank, "final-flux rollback fixture initializes fresh state");
  ConstFaceFluxView fresh_active;
  ConstFaceFluxView fresh_previous;
  passed &= expect(
      static_cast<bool>(writer.committed(final_storage, fresh_active)) &&
          static_cast<bool>(writer.committed_previous(final_storage,
                                                       fresh_previous)) &&
          fresh_active.revision == fresh_previous.revision &&
          fresh_active.certificate == fresh_previous.certificate &&
          flux_bits(fresh_active) == flux_bits(fresh_previous),
      rank,
      "fresh final flux publishes identical accepted and previous history");

  const auto publish_attempt = [&](double value) {
    if (!transaction.begin(layers) || !transaction.revise_trial(field)) {
      return false;
    }
    PendingFaceFluxView pending;
    if (!writer.begin_pending(transaction, final_storage, pending) ||
        !overwrite_pending_flux(pending, value)) {
      return false;
    }
    const std::array dependencies{RevisionDependency{
        AttemptTransaction::field_revision_source(field),
        transaction.trial_revision(field)}};
    return static_cast<bool>(writer.publish_pending(
        Span<const RevisionDependency>{dependencies.data(),
                                       dependencies.size()},
        pending));
  };

  passed &= expect(publish_attempt(1000.0 + rank), rank,
                   "first collective candidate writes pending final flux");
  PreparedAttemptFinish accepted;
  passed &= expect(static_cast<bool>(transaction.collective_prepare(
                       MPI_COMM_WORLD, Status{}, accepted)) &&
                       accepted.decision() == AttemptFinishDecision::accept,
                   rank, "first collective candidate prepares accept");
  transaction.commit_accept(accepted);

  ConstFaceFluxView active;
  ConstFaceFluxView previous;
  passed &= expect(static_cast<bool>(writer.committed(final_storage, active)) &&
                       static_cast<bool>(writer.committed_previous(
                           final_storage, previous)),
                   rank, "accepted final-flux history is observable");
  const auto active_bits = flux_bits(active);
  const auto previous_bits = flux_bits(previous);
  const auto active_certificate = active.certificate;
  const auto previous_certificate = previous.certificate;
  const auto active_revision = active.revision;
  const auto previous_revision = previous.revision;

  passed &= expect(publish_attempt(-1.0e200 - rank), rank,
                   "second collective candidate poisons pending final flux");
  const Status local = rank == 0
                           ? Status{StatusCode::rejected_step, 801U}
                           : Status{};
  PreparedAttemptFinish rejected;
  passed &= expect(static_cast<bool>(transaction.collective_prepare(
                       MPI_COMM_WORLD, local, rejected)) &&
                       rejected.decision() == AttemptFinishDecision::reject &&
                       rejected.lowest_failing_rank() == 0,
                   rank, "rank-zero failure collectively prepares rollback");
  transaction.commit_reject(rejected);

  ConstFaceFluxView active_after_reject;
  ConstFaceFluxView previous_after_reject;
  passed &= expect(
      static_cast<bool>(writer.committed(final_storage,
                                         active_after_reject)) &&
          static_cast<bool>(writer.committed_previous(
              final_storage, previous_after_reject)) &&
          active_after_reject.revision == active_revision &&
          previous_after_reject.revision == previous_revision &&
          active_after_reject.certificate == active_certificate &&
          previous_after_reject.certificate == previous_certificate &&
          flux_bits(active_after_reject) == active_bits &&
          flux_bits(previous_after_reject) == previous_bits,
      rank,
      "one-rank rejection preserves all ranks' active/previous flux bytes, "
      "revisions, and certificates bitwise");
  return all_true(passed);
}

bool test_rank_local_final_flux_preflight_is_atomic(int rank) {
  StateLayers layers;
  FieldId field{};
  bool passed = expect(make_layers(layers, field), rank,
                       "rank-local preflight fixture allocates state");
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
          static_cast<bool>(source_storage.workspace_view(0U, 151U, source)) &&
          static_cast<bool>(FaceFluxStorage::allocate_final(cells,
                                                            final_storage)) &&
          static_cast<bool>(authority.claim(191U, 0U, transaction, writer)),
      rank, "rank-local preflight fixture initializes authority");
  fill_flux(source, 100.0 * static_cast<double>(rank + 1));
  passed &= expect(static_cast<bool>(writer.initialize_committed(
                       final_storage, hundun::v04::as_const(source))),
                   rank, "rank-local preflight fixture initializes baseline");
  ConstFaceFluxView baseline;
  passed &= expect(static_cast<bool>(writer.committed(final_storage, baseline)),
                   rank, "rank-local preflight baseline is observable");
  const auto baseline_bits = flux_bits(baseline);
  const auto baseline_certificate = baseline.certificate;
  const RevisionToken baseline_revision = baseline.revision;

  passed &= expect(static_cast<bool>(transaction.begin(layers)) &&
                       static_cast<bool>(transaction.revise_trial(field)),
                   rank, "rank-local preflight attempt begins");
  PendingFaceFluxView pending;
  passed &= expect(static_cast<bool>(writer.begin_pending(
                       transaction, final_storage, pending)) &&
                       overwrite_pending_flux(pending, 1900.0 + rank),
                   rank, "rank-local preflight candidate owns pending flux");
  const RevisionToken pending_revision = pending.revision();
  const RevisionDependency correct{
      AttemptTransaction::field_revision_source(field),
      transaction.trial_revision(field)};
  const RevisionDependency local_dependency =
      rank == 0
          ? RevisionDependency{
                correct.source,
                static_cast<RevisionToken>(correct.revision + 1U)}
          : correct;
  const Status local_preflight = writer.preflight_publish_pending(
      Span<const RevisionDependency>{&local_dependency, 1U}, pending);
  const int local_ready = local_preflight ? 1 : 0;
  int globally_ready = 0;
  passed &= expect(MPI_Allreduce(&local_ready, &globally_ready, 1, MPI_INT,
                                 MPI_MIN, MPI_COMM_WORLD) == MPI_SUCCESS,
                   rank, "rank-local final-flux preflight reaches consensus");
  ConstFaceFluxView after_rejected_preflight;
  passed &= expect(
      (rank == 0 ? local_preflight.code == StatusCode::invalid_plan
                 : static_cast<bool>(local_preflight)) &&
          globally_ready == 0 && pending.valid() &&
          pending.revision() == pending_revision && transaction.active() &&
          static_cast<bool>(writer.committed(final_storage,
                                             after_rejected_preflight)) &&
          after_rejected_preflight.revision == baseline_revision &&
          after_rejected_preflight.certificate == baseline_certificate &&
          flux_bits(after_rejected_preflight) == baseline_bits,
      rank,
      "one-rank dependency failure performs zero write and preserves every "
      "rank's retryable pending lease");

  const Status corrected_preflight = writer.preflight_publish_pending(
      Span<const RevisionDependency>{&correct, 1U}, pending);
  passed &= expect(all_true(static_cast<bool>(corrected_preflight)), rank,
                   "corrected dependency preflights on every rank");
  passed &= expect(static_cast<bool>(writer.publish_pending(
                       Span<const RevisionDependency>{&correct, 1U}, pending)),
                   rank, "corrected dependency publishes on every rank");
  PreparedAttemptFinish accepted;
  passed &= expect(static_cast<bool>(transaction.collective_prepare(
                       MPI_COMM_WORLD, Status{}, accepted)) &&
                       accepted.decision() == AttemptFinishDecision::accept,
                   rank, "corrected final flux collectively prepares accept");
  transaction.commit_accept(accepted);
  ConstFaceFluxView committed;
  passed &= expect(static_cast<bool>(writer.committed(final_storage, committed)) &&
                       committed.revision == pending_revision,
                   rank,
                   "corrected final flux commits after the rejected preflight");
  return all_true(passed);
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

bool test_failure_severity_then_lowest_rank(int rank, int size) {
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
  int expected_rank = 0;
  if (size == 1) {
    local = Status{StatusCode::rejected_step, 701U};
    expected = local;
  } else {
    if (rank == 0) {
      local = Status{StatusCode::numerical_failure, 701U};
    } else if (rank == 1) {
      local = Status{StatusCode::invalid_plan, 702U};
    } else if (rank == size - 1) {
      local = Status{StatusCode::invalid_plan, 703U};
    }
    // Fatal outcomes dominate retryable numerical outcomes even when the
    // retryable failure is reported by a lower rank.  Rank one is also the
    // lowest reporter within the fatal severity class.
    expected = Status{StatusCode::invalid_plan, 702U};
    expected_rank = 1;
  }

  PreparedAttemptFinish prepared;
  const Status prepare_status =
      transaction.collective_prepare(MPI_COMM_WORLD, local, prepared);
  const Status result = prepared.outcome();
  passed &= expect(static_cast<bool>(prepare_status) && prepared.valid() &&
                       prepared.decision() == AttemptFinishDecision::reject &&
                       prepared.lowest_failing_rank() == expected_rank &&
                       transaction.active() && !transaction.finished() &&
                       layers.handle(StateRole::accepted_n) == old_n &&
                       layers.handle(StateRole::accepted_n_minus_one) == old_nm1 &&
                       layers.handle(StateRole::trial) == old_trial &&
                       !transaction.pending_cache_valid(RevisionSlotId{1U}),
                   rank,
                   "collective prepare selects severity before rank without "
                   "publishing local state");
  passed &= expect(result.code == expected.code && result.detail == expected.detail,
                   rank,
                   "the lowest rank within the strongest severity supplies "
                   "the exact status");
  passed &= expect(identical(packed_status(result)), rank,
                   "failure status category/detail is identical everywhere");
  const std::uint64_t decision_wire =
      (static_cast<std::uint64_t>(prepared.decision()) << 32U) |
      static_cast<std::uint32_t>(prepared.lowest_failing_rank());
  passed &= expect(identical(decision_wire), rank,
                   "prepared failure decision and lowest rank are identical");
  transaction.commit_reject(prepared);
  passed &= expect(!transaction.committed(), rank,
                   "any local failure rolls back every rank");
  passed &= expect(transaction.lowest_failing_rank() == expected_rank, rank,
                   "rollback records the selected severity-class rank");
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
  passed &= test_failure_severity_then_lowest_rank(rank, size);
  passed &= test_rank_local_begin_failure(rank, size);
  passed &= test_rank_local_inactive_finish(rank, size);
  passed &= test_final_flux_collective_rollback(rank);
  passed &= test_rank_local_final_flux_preflight_is_atomic(rank);
  passed = all_true(passed);
  if (rank == 0 && !passed) {
    std::cerr << "v04 core transaction MPI test failed for " << size
              << " ranks\n";
  }
  MPI_Finalize();
  return passed ? 0 : 1;
}

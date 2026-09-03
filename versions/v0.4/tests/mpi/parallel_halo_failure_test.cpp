// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_parallel.hpp"
#include "parallel_halo_detail.hpp"

#include <mpi.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string_view>
#include <vector>

namespace {

using namespace hundun::v04;

constexpr Int3 kInterior{4, 3, 2};
constexpr Int3 kGhosts{1, 1, 1};
constexpr FieldId kField = 0U;
constexpr RevisionToken kRevision = 77U;

struct OwnedField {
  std::vector<double> storage;
  FieldView view{};
};

struct FailureCase {
  detail::HaloFailurePoint point;
  StatusCode status_code;
  std::uint32_t detail_code;
  std::string_view name;
};

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

bool identical_status(Status status) {
  const std::uint64_t packed = packed_status(status);
  std::uint64_t minimum = 0U;
  std::uint64_t maximum = 0U;
  return MPI_Allreduce(&packed, &minimum, 1, MPI_UINT64_T, MPI_MIN,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         MPI_Allreduce(&packed, &maximum, 1, MPI_UINT64_T, MPI_MAX,
                       MPI_COMM_WORLD) == MPI_SUCCESS &&
         minimum == maximum;
}

OwnedField make_field(int rank) {
  OwnedField result;
  const std::size_t stride_y =
      static_cast<std::size_t>(kInterior.x + 2 * kGhosts.x);
  const std::size_t stride_z =
      stride_y * static_cast<std::size_t>(kInterior.y + 2 * kGhosts.y);
  const std::size_t component_stride =
      stride_z * static_cast<std::size_t>(kInterior.z + 2 * kGhosts.z);
  result.storage.resize(component_stride);
  for (std::size_t index = 0U; index < result.storage.size(); ++index) {
    result.storage[index] =
        static_cast<double>(rank * 1000000) - 500000.0 -
        static_cast<double>(index);
  }
  result.view.base = result.storage.data() +
                     static_cast<std::ptrdiff_t>(kGhosts.x) +
                     static_cast<std::ptrdiff_t>(kGhosts.y) *
                         static_cast<std::ptrdiff_t>(stride_y) +
                     static_cast<std::ptrdiff_t>(kGhosts.z) *
                         static_cast<std::ptrdiff_t>(stride_z);
  result.view.interior = kInterior;
  result.view.ghosts = kGhosts;
  result.view.components = 1U;
  result.view.stride_y = stride_y;
  result.view.stride_z = stride_z;
  result.view.component_stride = component_stride;
  result.view.field = kField;
  result.view.revision = kRevision;
  result.view.storage_identity = 8080U;
  result.view.revision_domain = 9090U;
  for (std::int32_t z = 0; z < kInterior.z; ++z) {
    for (std::int32_t y = 0; y < kInterior.y; ++y) {
      for (std::int32_t x = 0; x < kInterior.x; ++x) {
        result.view.unchecked(Int3{x, y, z}, 0U) =
            static_cast<double>(rank * 1000000 + z * 100 + y * 10 + x);
      }
    }
  }
  return result;
}

void reset_ghosts(OwnedField& field, double sentinel) {
  for (std::int32_t z = -1; z <= kInterior.z; ++z) {
    for (std::int32_t y = -1; y <= kInterior.y; ++y) {
      for (std::int32_t x = -1; x <= kInterior.x; ++x) {
        const bool interior = x >= 0 && x < kInterior.x && y >= 0 &&
                              y < kInterior.y && z >= 0 && z < kInterior.z;
        if (!interior) {
          field.view.unchecked(Int3{x, y, z}, 0U) = sentinel;
        }
      }
    }
  }
}

bool verify_recovered_x_ghosts(const OwnedField& field, int rank, int size) {
  bool passed = true;
  const int lower_rank = rank == 0 ? size - 1 : rank - 1;
  const int upper_rank = rank + 1 == size ? 0 : rank + 1;
  for (std::int32_t z = 0; z < kInterior.z; ++z) {
    for (std::int32_t y = 0; y < kInterior.y; ++y) {
      const double expected_lower = static_cast<double>(
          lower_rank * 1000000 + z * 100 + y * 10 + kInterior.x - 1);
      const double expected_upper =
          static_cast<double>(upper_rank * 1000000 + z * 100 + y * 10);
      passed &= expect(field.view.unchecked(Int3{-1, y, z}, 0U) ==
                           expected_lower,
                       rank, "recovered x-minus ghost is exact");
      passed &= expect(field.view.unchecked(Int3{kInterior.x, y, z}, 0U) ==
                           expected_upper,
                       rank, "recovered x-plus ghost is exact");
    }
  }
  return passed;
}

MeshPatch patch_for(int rank, int size) {
  return MeshPatch{Int3{rank * kInterior.x, 0, 0}, kInterior,
                   Int3{size, 1, 1}, Int3{rank, 0, 0}};
}

bool run_failure_case(const FailureCase& scenario, int rank, int size) {
  OwnedField field = make_field(rank);
  const std::array specs{HaloFieldSpec{kField, 1U, 1U}};
  std::array views{field.view};
  HaloEngine engine;
  bool passed = expect(!engine.ready(), rank,
                       "default halo engine is not ready");
  passed &= expect(static_cast<bool>(engine.reserve(
                           MPI_COMM_WORLD, patch_for(rank, size),
                           Span<const HaloFieldSpec>{specs.data(), specs.size()},
                           HaloTopology{true, false, false})),
                       rank, "failure fixture reserves a periodic halo plan");
  passed &= expect(engine.ready(), rank,
                   "reserved halo engine reports ready");
  passed &= expect(engine.ghost_revision(kField) == 0U, rank,
                   "failure fixture begins without a certificate");

  HaloTicket warmup_ticket;
  passed &= expect(static_cast<bool>(engine.begin(
                       49U,
                       Span<const FieldView>{views.data(), views.size()},
                       warmup_ticket)) &&
                       static_cast<bool>(engine.finish(
                           warmup_ticket,
                           Span<FieldView>{views.data(), views.size()})) &&
                       engine.ghost_revision(kField) == kRevision,
                   rank, "warmup publishes the prior ghost certificate");
  constexpr double sentinel = -424242.0;
  reset_ghosts(field, sentinel);
  const std::vector<double> snapshot = field.storage;

  const int failing_rank = size == 1 ? 0 : size - 1;
  detail::set_halo_failure_for_test(scenario.point, failing_rank);
  HaloTicket ticket;
  Status result = engine.begin(
      51U, Span<const FieldView>{views.data(), views.size()}, ticket);
  if (scenario.point == detail::HaloFailurePoint::completion ||
      scenario.point == detail::HaloFailurePoint::unpack) {
    passed &= expect(static_cast<bool>(result) && ticket.active(), rank,
                     "completion/unpack injection permits collective begin");
    passed &= expect(engine.ghost_revision(kField) == 0U, rank,
                     "a new begin invalidates the prior ghost certificate");
    result = engine.finish(ticket,
                           Span<FieldView>{views.data(), views.size()});
  } else {
    passed &= expect(!result && !ticket.active(), rank,
                     "pack/start injection rejects collectively in begin");
  }
  detail::clear_halo_failure_for_test();

  passed &= expect(result.code == scenario.status_code &&
                       result.detail == scenario.detail_code,
                   rank, scenario.name);
  passed &= expect(identical_status(result), rank,
                   "all ranks report one exact failure category/detail");
  passed &= expect(engine.lowest_failing_rank() == failing_rank, rank,
                   "failure consensus exposes the selected lowest rank");
  passed &= expect(engine.ghost_revision(kField) == 0U, rank,
                   "failed exchange publishes no ghost certificate");
  passed &= expect(field.storage == snapshot, rank,
                   "failed exchange mutates no interior or ghost cell");
  passed &= expect(engine.ready(), rank,
                   "recoverable exchange failure preserves ready state");

  HaloTicket recovery_ticket;
  const Status recovery_begin = engine.begin(
      52U, Span<const FieldView>{views.data(), views.size()}, recovery_ticket);
  passed &= expect(static_cast<bool>(recovery_begin) &&
                       recovery_ticket.active() &&
                       engine.ghost_revision(kField) == 0U,
                   rank,
                   "the same engine starts cleanly after injected failure");
  const Status recovery_finish = engine.finish(
      recovery_ticket, Span<FieldView>{views.data(), views.size()});
  passed &= expect(static_cast<bool>(recovery_finish) &&
                       !recovery_ticket.active() &&
                       engine.ghost_revision(kField) == kRevision,
                   rank, "successful recovery republishes the exact certificate");
  passed &= verify_recovered_x_ghosts(field, rank, size);
  return all_true(passed);
}

bool run_reserve_failure_case(detail::HaloFailurePoint point, int rank,
                              int size) {
  OwnedField field = make_field(rank);
  const std::array specs{HaloFieldSpec{kField, 1U, 1U}};
  std::array views{field.view};
  HaloEngine engine;
  bool passed = expect(static_cast<bool>(engine.reserve(
                           MPI_COMM_WORLD, patch_for(rank, size),
                           Span<const HaloFieldSpec>{specs.data(), specs.size()},
                           HaloTopology{true, false, false})),
                       rank, "reserve-failure fixture publishes an old plan");
  const HaloPlanStats before = engine.plan_stats();
  HaloTicket warmup;
  passed &= expect(static_cast<bool>(engine.begin(
                       71U, Span<const FieldView>{views.data(), views.size()},
                       warmup)) &&
                       static_cast<bool>(engine.finish(
                           warmup, Span<FieldView>{views.data(), views.size()})) &&
                       engine.ghost_revision(kField) == kRevision,
                   rank, "old plan exchanges before failed replacement");

  const int failing_rank = size == 1 ? 0 : size - 1;
  detail::set_halo_failure_for_test(point, failing_rank);
  const Status result = engine.reserve(
      MPI_COMM_WORLD, patch_for(rank, size),
      Span<const HaloFieldSpec>{specs.data(), specs.size()},
      HaloTopology{true, true, false});
  detail::clear_halo_failure_for_test();
  const HaloPlanStats after = engine.plan_stats();
  passed &= expect(result.code == StatusCode::allocation_failure &&
                       result.detail ==
                           detail::halo_detail_reserve_allocation &&
                       identical_status(result) &&
                       engine.lowest_failing_rank() == failing_rank,
                   rank,
                   "reserve allocation failure reaches exact consensus");
  passed &= expect(before.request_storage_address ==
                           after.request_storage_address &&
                       before.send_storage_address ==
                           after.send_storage_address &&
                       before.receive_storage_address ==
                           after.receive_storage_address &&
                       engine.ghost_revision(kField) == kRevision,
                   rank, "failed reserve preserves the published old engine");

  HaloTicket recovery;
  passed &= expect(static_cast<bool>(engine.begin(
                       72U, Span<const FieldView>{views.data(), views.size()},
                       recovery)) &&
                       static_cast<bool>(engine.finish(
                           recovery, Span<FieldView>{views.data(), views.size()})) &&
                       engine.ghost_revision(kField) == kRevision,
                   rank, "old engine exchanges after failed reserve");
  return all_true(passed);
}

bool run_prerequisite_failure_case(int rank, int size) {
  OwnedField field = make_field(rank);
  const std::array specs{HaloFieldSpec{kField, 1U, 1U}};
  std::array views{field.view};
  HaloEngine engine;
  bool passed = expect(static_cast<bool>(engine.reserve(
                           MPI_COMM_WORLD, patch_for(rank, size),
                           Span<const HaloFieldSpec>{specs.data(), specs.size()},
                           HaloTopology{true, false, false})),
                       rank, "prerequisite-failure fixture reserves");
  constexpr double sentinel = -515151.0;
  reset_ghosts(field, sentinel);
  const std::vector<double> snapshot = field.storage;
  const int failing_rank = size == 1 ? 0 : size - 1;
  const Status prerequisite =
      rank == failing_rank
          ? Status{StatusCode::numerical_failure, 0x505245U}
          : Status{};
  HaloTicket ticket;
  const Status result = engine.begin(
      81U, Span<const FieldView>{views.data(), views.size()}, prerequisite,
      ticket);
  const HaloRuntimeCounters counters = engine.runtime_counters();
  passed &= expect(result.code == StatusCode::numerical_failure &&
                       result.detail == 0x505245U && identical_status(result),
                   rank,
                   "rank-local prerequisite reaches exact halo consensus");
  passed &= expect(engine.lowest_failing_rank() == failing_rank &&
                       !ticket.active() && !engine.active(),
                   rank,
                   "prerequisite failure starts no exchange");
  passed &= expect(counters.begin_calls == 1U &&
                       counters.finish_calls == 0U &&
                       counters.messages_started == 0U &&
                       counters.bytes_packed == 0U &&
                       counters.bytes_unpacked == 0U,
                   rank,
                   "prerequisite failure performs no transport work");
  passed &= expect(field.storage == snapshot &&
                       engine.ghost_revision(kField) == 0U,
                   rank,
                   "prerequisite failure publishes no field mutation");
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

  const std::array scenarios{
      FailureCase{detail::HaloFailurePoint::pack, StatusCode::invalid_plan,
                  detail::halo_detail_pack_failure,
                  "pack failure reaches exact collective consensus"},
      FailureCase{detail::HaloFailurePoint::start, StatusCode::mpi_failure,
                  detail::halo_detail_start_failure,
                  "start failure reaches exact collective consensus"},
      FailureCase{detail::HaloFailurePoint::completion, StatusCode::mpi_failure,
                  detail::halo_detail_completion_failure,
                  "completion failure reaches exact collective consensus"},
      FailureCase{detail::HaloFailurePoint::unpack, StatusCode::invalid_plan,
                  detail::halo_detail_unpack_failure,
                  "unpack failure reaches exact collective consensus"},
  };

  bool passed = true;
  for (const FailureCase& scenario : scenarios) {
    passed &= run_failure_case(scenario, rank, size);
  }
  passed &= run_reserve_failure_case(
      detail::HaloFailurePoint::reserve_before_contract, rank, size);
  passed &= run_reserve_failure_case(
      detail::HaloFailurePoint::reserve_before_alltoall, rank, size);
  passed &= run_prerequisite_failure_case(rank, size);
  passed = all_true(passed);
  if (rank == 0 && !passed) {
    std::cerr << "v0.4 halo failure test failed for " << size << " ranks\n";
  }
  MPI_Finalize();
  return passed ? 0 : 1;
}

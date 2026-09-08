// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09

#include "../../src/mesh_focus_detail.hpp"
#include "../../src/models_spray_migration_detail.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <mpi.h>
#include <new>
#include <vector>

namespace {
bool count_hot_allocations = false;
std::size_t hot_allocations = 0U;
} // namespace
void *operator new(std::size_t bytes) {
  if (count_hot_allocations)
    ++hot_allocations;
  if (void *result = std::malloc(bytes ? bytes : 1U))
    return result;
  throw std::bad_alloc{};
}
void *operator new[](std::size_t bytes) { return ::operator new(bytes); }
void operator delete(void *pointer) noexcept { std::free(pointer); }
void operator delete[](void *pointer) noexcept { std::free(pointer); }
void operator delete(void *pointer, std::size_t) noexcept {
  std::free(pointer);
}
void operator delete[](void *pointer, std::size_t) noexcept {
  std::free(pointer);
}

using namespace hundun::v04;
using namespace hundun::v04::spray;
using namespace hundun::v04::spray::detail;

namespace {
std::uint64_t global_id(Int3 index, Int3 cells) {
  return static_cast<std::uint64_t>(index.x) +
         static_cast<std::uint64_t>(cells.x) *
             (static_cast<std::uint64_t>(index.y) +
              static_cast<std::uint64_t>(cells.y) * index.z);
}
ParcelMigrationValue value(int origin, int ordinal, int ranks, Int3 cells) {
  MeshPatch target;
  hundun::v04::detail::make_mesh_patch((origin + ordinal + 1) % ranks, ranks,
                                       cells, target);
  ParcelMigrationValue result;
  result.parcel.id = {static_cast<std::uint64_t>(origin + 1),
                      static_cast<std::uint64_t>(ordinal + 1)};
  result.parcel.position_m = {0.25 + origin, 0.5 + ordinal, -0.0};
  result.parcel.velocity_m_per_s = {1.0, -2.0, 3.0};
  result.parcel.droplet_mass_kg = 1.0e-9 * (origin + ordinal + 1);
  result.parcel.droplet_diameter_m = 1.0e-4;
  result.parcel.multiplicity = 5.0;
  result.parcel.temperature_k = 300.0 + origin;
  result.parcel.liquid_material_fingerprint = 0x125;
  result.parcel.owner_global_cell = global_id(target.begin, cells);
  result.breakup_ordinal =
      UINT64_C(0x123456789abcdef0) + static_cast<std::uint64_t>(ordinal);
  result.parcel.age_s = 0.5;
  result.tab_deformation = -0.25 - ordinal;
  result.tab_deformation_rate_per_s = 4.0 + origin;
  return result;
}
bool same(const ParcelMigrationValue &a, const ParcelMigrationValue &b) {
  return a.parcel.id == b.parcel.id &&
         std::memcmp(a.parcel.position_m.data(), b.parcel.position_m.data(),
                     3 * sizeof(double)) == 0 &&
         a.parcel.velocity_m_per_s == b.parcel.velocity_m_per_s &&
         a.parcel.droplet_mass_kg == b.parcel.droplet_mass_kg &&
         a.parcel.droplet_diameter_m == b.parcel.droplet_diameter_m &&
         a.parcel.multiplicity == b.parcel.multiplicity &&
         a.parcel.temperature_k == b.parcel.temperature_k &&
         a.parcel.liquid_material_fingerprint ==
             b.parcel.liquid_material_fingerprint &&
         a.parcel.owner_global_cell == b.parcel.owner_global_cell &&
         a.parcel.age_s == b.parcel.age_s &&
         a.tab_deformation == b.tab_deformation &&
         a.tab_deformation_rate_per_s == b.tab_deformation_rate_per_s &&
         a.breakup_ordinal == b.breakup_ordinal;
}
} // namespace

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank = 0, ranks = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  bool passed = true;
  auto expect = [&](bool condition, const char *message) {
    if (!condition)
      std::cerr << "rank " << rank << ": " << message << '\n';
    passed &= condition;
  };
  const Int3 cells{8, 8, 8};
  MeshPatch patch;
  expect(static_cast<bool>(
             hundun::v04::detail::make_mesh_patch(rank, ranks, cells, patch)),
         "partition fixture");
  ParcelMigrationPlan migration;
  expect(static_cast<bool>(
             migration.configure(MPI_COMM_WORLD, cells, patch, 8U, 8U * ranks)),
         "bounded migration plan");
  std::vector<ParcelMigrationValue> input{value(rank, 0, ranks, cells),
                                          value(rank, 1, ranks, cells)};
  const auto original = input;
  ParcelMigrationReport report;
  count_hot_allocations = true;
  const Status first_prepare =
      migration.prepare({input.data(), input.size()}, report);
  count_hot_allocations = false;
  expect(static_cast<bool>(first_prepare),
         "same-rank, neighbor and multi-rank staged exchange");
  expect(hot_allocations == 0U, "prepare uses only cold-reserved C++ buffers");
  expect(report.available && report.global_parcels == 2U * ranks,
         "global count preserved");
  auto candidates = migration.candidates();
  expect(candidates.size == 2U, "exactly two incoming parcels");
  for (std::size_t i = 0; i < candidates.size; ++i) {
    const auto &actual = candidates.data[i];
    const auto expected =
        value(static_cast<int>(actual.parcel.id.high - 1),
              static_cast<int>(actual.parcel.id.low - 1), ranks, cells);
    expect(same(actual, expected),
           "all values including TAB and signed zero survive");
    expect((static_cast<int>(actual.parcel.id.high + actual.parcel.id.low - 1) %
            ranks) == rank,
           "one authoritative destination owner");
    if (i)
      expect(candidates.data[i - 1].parcel.id < actual.parcel.id,
             "candidate order is independent of arrival order");
  }
  expect(same(input[0], original[0]) && same(input[1], original[1]),
         "prepare never mutates source state");
  migration.discard();
  expect(migration.candidates().size == 0U,
         "discard withdraws every staged owner");

  // Failure is introduced on one rank only. All ranks must leave no candidates.
  if (rank == ranks - 1)
    input[0].parcel.droplet_mass_kg = -1.0;
  expect(!migration.prepare({input.data(), input.size()}, report) &&
             !report.available && migration.candidates().size == 0U,
         "one-rank validation failure rejects all staged owners");
  input = original;
  expect(static_cast<bool>(
             migration.prepare({input.data(), input.size()}, report)),
         "retry recomputes migration after failure");
  migration.discard();

  // Duplicate identities may be sent to different owners; audit must be global.
  input[0].parcel.id = {71U, 89U};
  if (ranks == 1)
    input[1].parcel.id = input[0].parcel.id;
  expect(!migration.prepare({input.data(), input.size()}, report) &&
             report.status.detail ==
                 static_cast<std::uint32_t>(MigrationDetail::duplicate_id) &&
             migration.candidates().size == 0U,
         "cross-owner duplicate ID rejected");
  input = original;
  if (rank == 0)
    input[0].parcel.owner_global_cell = UINT64_MAX;
  expect(!migration.prepare({input.data(), input.size()}, report) &&
             migration.candidates().size == 0U,
         "invalid global cell is collective failure");

  expect(static_cast<bool>(migration.prepare({}, report)) && report.available &&
             report.global_parcels == 0U && migration.candidates().size == 0U,
         "zero-parcel collective succeeds");
  migration.discard();
  ParcelMigrationPlan limited;
  expect(static_cast<bool>(
             limited.configure(MPI_COMM_WORLD, cells, patch, 1U, 8U * ranks)),
         "small capacity fixture");
  input = original;
  expect(!limited.prepare({input.data(), input.size()}, report) &&
             !report.available,
         "capacity excess fails before payload publication");

  // Capacity can overflow on the receiver even with valid local input counts.
  if (ranks > 1) {
    std::vector<ParcelMigrationValue> converging{original[0]};
    converging[0].parcel.owner_global_cell = 0U;
    expect(!limited.prepare({converging.data(), converging.size()}, report) &&
               report.status.detail ==
                   static_cast<std::uint32_t>(
                       MigrationDetail::capacity_exceeded) &&
               report.lowest_failing_rank == 0 && !report.available &&
               limited.candidates().size == 0U,
           "incoming capacity excess rejects every owner before exchange");
  }
  ParcelMigrationPlan audit_limited;
  expect(static_cast<bool>(
             audit_limited.configure(MPI_COMM_WORLD, cells, patch, 8U, 1U)),
         "small ID auditor capacity fixture");
  input = original;
  for (std::size_t i = 0; i < input.size(); ++i) {
    const std::uint64_t id = static_cast<std::uint64_t>(2 * rank) + i + 1U;
    input[i].parcel.id = {id, id}; // unique full IDs, common auditor bucket
  }
  expect(!audit_limited.prepare({input.data(), input.size()}, report) &&
             report.status.detail == static_cast<std::uint32_t>(
                                         MigrationDetail::capacity_exceeded) &&
             !report.available && audit_limited.candidates().size == 0U,
         "ID audit capacity failure withdraws already received candidates");
  expect(static_cast<bool>(audit_limited.prepare({}, report)) &&
             report.available,
         "empty retry works after audit-stage failure");
  struct Location final : ParcelLocationProvider {
    bool stale{};
    Status locate(const Vector3 &position,
                  hundun::v04::portable::Revision revision,
                  ParcelLocation &out) const noexcept override {
      out = {revision, 0U, 0};
      out.global_cell = static_cast<std::uint64_t>(position[0]);
      if (stale)
        ++out.revision.input_revision;
      return {};
    }
  } locator;
  input = original;
  for (auto &v : input) {
    v.parcel.owner_global_cell = 0;
    v.parcel.position_m = {0.0, 0.0, 0.0};
  }
  const hundun::v04::portable::Revision revision{11, 37, 1};
  CartesianParcelLocationProvider cartesian(
      {revision, {0, 0, 0}, {1, 1, 1}, cells, ranks});
  ParcelLocation located;
  expect(static_cast<bool>(
             cartesian.locate({1.25, 2.25, 3.25}, revision, located)) &&
             located.global_cell == 209U,
         "Cartesian physical position maps to x-fast cell 209");
  expect(!cartesian.locate({8.0, 0, 0}, revision, located) &&
             located.owner_rank == -1,
         "upper-domain outlet is not silently clamped into a cell");
  hot_allocations = 0;
  count_hot_allocations = true;
  const auto located_status = migration.prepare_checked(
      {input.data(), input.size()}, cartesian, revision, report);
  count_hot_allocations = false;
  expect(static_cast<bool>(located_status) && hot_allocations == 0,
         "checked Cartesian migration validates endpoint position and owner "
         "without hot allocations");
  locator.stale = true;
  expect(!migration.prepare_checked({input.data(), input.size()}, locator,
                                    revision, report) &&
             !report.available && migration.candidates().size == 0 &&
             report.status.detail ==
                 static_cast<std::uint32_t>(MigrationDetail::stale_revision),
         "stale geometry response withdraws all candidates collectively");
  locator.stale = false;
  if (rank == 0)
    input[0].parcel.position_m[0] = 1.0;
  expect(!migration.prepare_checked({input.data(), input.size()}, locator,
                                    revision, report) &&
             !report.available && migration.candidates().size == 0,
         "wrong position/cell mapping cannot migrate under a plausible owner "
         "field");
  int local = passed ? 1 : 0, global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Finalize();
  return global ? 0 : 1;
}

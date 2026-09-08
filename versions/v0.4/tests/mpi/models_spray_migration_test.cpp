// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09

#include "../../src/core_spray_gas_detail.hpp"
#include "../../src/mesh_focus_detail.hpp"
#include "../../src/models_spray_migration_detail.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
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
bool native_gas_halo(int rank, int ranks) {
  const Int3 global{17, 11, 7};
  CartesianMeshSpec mesh;
  mesh.kind = GeometryKind::coast_runtime_axes_v1;
  mesh.axes_file = "manufactured-axes.dat";
  mesh.has_exact_cells = true;
  mesh.exact_cells = global;
  mesh.lower = {0, 0, 0};
  mesh.upper = {1, 1, 1};
  mesh.minimum_spacing = {1e-8, 1e-8, 1e-8};
  mesh.max_growth_ratio = 10;
  mesh.limits = {100000, UINT64_C(67108864)};
  const int sizes[]{global.x, global.y, global.z};
  for (unsigned d = 0; d < 3; ++d)
    for (int i = 0; i <= sizes[d]; ++i) {
      const double t = double(i) / sizes[d];
      mesh.coast_runtime_faces[d].push_back(.5 * t + .5 * t * t);
    }
  CartesianGeometryPlan geometry;
  MeshPatch patch;
  if (!CartesianGeometryCompiler::compile(MPI_COMM_WORLD, mesh, {}, geometry,
                                          patch))
    return false;
  std::array<std::vector<double>, 5> data;
  std::array<FieldView, 5> fields{};
  std::array<ConstFieldView, 2> species{};
  std::array<HaloFieldSpec, 5> halo_fields{};
  for (std::size_t f = 0; f < fields.size(); ++f) {
    const unsigned nc = f == 2 ? 3 : 1;
    const auto sy = std::size_t(patch.cells.x + 2),
               sz = sy * (patch.cells.y + 2);
    const auto stride = sz * (patch.cells.z + 2);
    data[f].assign(nc * stride, std::numeric_limits<double>::quiet_NaN());
    auto &v = fields[f];
    v.base = data[f].data() + 1 + sy + sz;
    v.interior = patch.cells;
    v.ghosts = {1, 1, 1};
    v.components = nc;
    v.stride_y = sy;
    v.stride_z = sz;
    v.component_stride = stride;
    v.field = f + 1;
    v.revision = 9;
    v.storage_identity = f + 100;
    v.revision_domain = 1000;
    halo_fields[f] = {v.field, 1, std::uint8_t(nc)};
    for (int z = 0; z < patch.cells.z; ++z)
      for (int y = 0; y < patch.cells.y; ++y)
        for (int x = 0; x < patch.cells.x; ++x) {
          const double value = (x + patch.begin.x) + 2 * (y + patch.begin.y) +
                               3 * (z + patch.begin.z);
          v.unchecked({x, y, z}, 0) = f == 0   ? 10 + value
                                      : f == 1 ? 300000 + 50 * value
                                      : f == 2 ? value
                                      : f == 3 ? .1 + .001 * value
                                               : .2 + .001 * value;
          if (f == 2) {
            v.unchecked({x, y, z}, 1) = -value;
            v.unchecked({x, y, z}, 2) = 0;
          }
        }
  }
  species = {as_const(fields[3]), as_const(fields[4])};
  HaloEngine halo;
  if (!halo.reserve(MPI_COMM_WORLD, patch,
                    {halo_fields.data(), halo_fields.size()},
                    {true, true, true}))
    return false;
  std::vector<GlobalCellId> donor_cells;
  std::vector<Int3> donor_targets;
  for (int z = -1; z <= patch.cells.z; ++z)
    for (int y = -1; y <= patch.cells.y; ++y)
      for (int x = -1; x <= patch.cells.x; ++x) {
        const int outside = (x < 0 || x >= patch.cells.x) +
                            (y < 0 || y >= patch.cells.y) +
                            (z < 0 || z >= patch.cells.z);
        if (outside < 2)
          continue;
        const Int3 canonical{(x + patch.begin.x + global.x) % global.x,
                             (y + patch.begin.y + global.y) % global.y,
                             (z + patch.begin.z + global.z) % global.z};
        donor_cells.push_back(global_id(canonical, global));
        donor_targets.push_back({x, y, z});
      }
  std::array<RemoteDonorFieldSpec, 5> donor_fields{};
  for (std::size_t i = 0; i < fields.size(); ++i)
    donor_fields[i] = {fields[i].field, fields[i].components};
  RemoteDonorExchangePlan corners;
  const RemoteDonorTargets targets{
      geometry.fingerprint(),
      1,
      {true, true, true},
      {donor_cells.data(), donor_cells.size()},
      {donor_targets.data(), donor_targets.size()}};
  if (!RemoteDonorExchangePlan::analyze_cells(
          MPI_COMM_WORLD, global, patch, targets,
          {donor_fields.data(), donor_fields.size()}, 20, corners) ||
      !corners.bind(MPI_COMM_WORLD))
    return false;
  if (ranks > 1) {
    auto mismatched = donor_fields;
    if (rank == 0)
      std::swap(mismatched[0], mismatched[1]);
    RemoteDonorExchangePlan invalid;
    if (RemoteDonorExchangePlan::analyze_cells(
            MPI_COMM_WORLD, global, patch, targets,
            {mismatched.data(), mismatched.size()}, 20, invalid) ||
        invalid.fingerprint()) {
      std::cerr << "rank " << rank
                << ": mismatched corner field order accepted\n";
      return false;
    }
  }
  hundun::v04::detail::ProductParcelGas sampler;
  const std::size_t mapping[]{2, 0};
  const portable::Revision revision{11, 71, 1};
  if (!sampler.configure(geometry, patch, {true, true, true}, 100, {mapping, 2},
                         1) ||
      !sampler.bind(revision, .1, 101325, as_const(fields[0]),
                    as_const(fields[1]), as_const(fields[2]),
                    {species.data(), species.size()}))
    return false;
  auto p = value(rank, 0, ranks, global).parcel;
  const int current[]{patch.begin.x, patch.begin.y, patch.begin.z};
  double expected = 0;
  for (unsigned d = 0; d < 3; ++d) {
    const auto centres = geometry.axis(static_cast<CartesianAxis>(d)).centres();
    const int before = (current[d] + sizes[d] - 1) % sizes[d];
    const double left = centres.data[before] - (current[d] == 0 ? 1.0 : 0.0);
    p.position_m[d] = .5 * (left + centres.data[current[d]]);
    expected += (d + 1) * .5 * (before + current[d]);
  }
  std::array<double, 3> y{-1, -1, -1};
  bool passed = true;
  if (ranks > 1)
    passed &=
        sampler.sample(p, .05, ParcelPass::corrector, revision, y.data(), 3)
            .status != portable::Status::success;
  HaloTicket ticket;
  if (!halo.begin(10, {fields.data(), fields.size()}, ticket) ||
      !halo.finish(ticket, {fields.data(), fields.size()}))
    return false;
  const auto ready =
      corners.preflight_exchange(20, {fields.data(), fields.size()});
  int local_ready = ready ? 1 : 0, all_ready = 0;
  MPI_Allreduce(&local_ready, &all_ready, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (!all_ready)
    return false;
  hot_allocations = 0;
  count_hot_allocations = true;
  const auto gathered = corners.exchange(20, {fields.data(), fields.size()});
  const auto sampled =
      sampler.sample(p, .05, ParcelPass::corrector, revision, y.data(), 3);
  ParcelLocation located;
  const auto location = sampler.locate(p.position_m, revision, located);
  count_hot_allocations = false;
  const auto near = [](double a, double b) {
    return std::isfinite(a) &&
           std::abs(a - b) < 2e-12 * std::max(1.0, std::abs(b));
  };
  passed &=
      bool(gathered) && sampled.status == portable::Status::success &&
      hot_allocations == 0 && near(sampled.pressure_pa, 101335 + expected) &&
      near(sampled.enthalpy_j_per_kg, 300000 + 50 * expected) &&
      near(sampled.velocity_m_per_s[0], expected) &&
      near(sampled.velocity_m_per_s[1], -expected) &&
      near(y[2], .1 + .001 * expected) && near(y[0], .2 + .001 * expected) &&
      near(y[1], .7 - .002 * expected) && bool(location) &&
      located.owner_rank >= 0 && located.owner_rank < ranks;
  if (!passed)
    std::cerr
        << "rank " << rank
        << ": native parcel gas samples actual periodic corner Halo: status="
        << unsigned(sampled.status) << " location=" << unsigned(location.code)
        << " allocations=" << hot_allocations << " p=" << sampled.pressure_pa
        << " expected=" << expected << " h=" << sampled.enthalpy_j_per_kg
        << " y=" << y[0] << "," << y[1] << "," << y[2]
        << " grid=" << patch.process_grid.x << "," << patch.process_grid.y
        << "," << patch.process_grid.z << '\n';
  return passed;
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
  passed &= native_gas_halo(rank, ranks);
  int local = passed ? 1 : 0, global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Finalize();
  return global ? 0 : 1;
}

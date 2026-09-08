// SPDX-License-Identifier: Apache-2.0
#include "mesh_focus_detail.hpp"
#include "models_chemistry_adapter_detail.hpp"
#include "models_portable_composition_detail.hpp"
#include <cmath>
#include <iostream>
#include <mpi.h>
using namespace hundun::v04;
class FailingChemistry final : public portable::GasAdvanceProvider {
public:
  portable::GasAdvanceProvider &backend;
  int calls{}, fail_at{};
  explicit FailingChemistry(portable::GasAdvanceProvider &p) : backend(p) {}
  const portable::GasIdentity &gas_identity() const noexcept override {
    return backend.gas_identity();
  }
  portable::Status
  advance_gas(const portable::GasAdvanceQuery &q,
              portable::GasAdvanceOutput &out) noexcept override {
    if (++calls == fail_at)
      return portable::Status::provider_failure;
    return backend.advance_gas(q, out);
  }
};
class BadGas final : public portable::GasQueryProvider {
public:
  chemistry::detail::AnalyticIsomerBackend &backend;
  explicit BadGas(chemistry::detail::AnalyticIsomerBackend &b) : backend(b) {}
  const portable::GasIdentity &gas_identity() const noexcept override {
    return backend.gas_identity();
  }
  portable::Status query_gas(const portable::GasQuery &q,
                             portable::GasQueryOutput &out) noexcept override {
    auto status = backend.query_gas(q, out);
    out.sample.density_kg_per_m3 = -1;
    return status;
  }
};
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank = 0, ranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  chemistry::detail::AnalyticIsomerBackend backend;
  spray::detail::LiquidAsset material;
  material.pack.material_fingerprint = 17;
  material.content_fingerprint = 17;
  material.gas_identity = backend.gas_identity();
  material.vapor_species_index = 0;
  material.vapor_species_name = "A";
  material.vapor_molecular_weight_kg_per_kmol = 28;
  spray::SprayParcelState mapped_parcel;
  mapped_parcel.liquid_material_fingerprint = 17;
  bool mapping_ok = portable::valid_vapor_mapping(mapped_parcel, material,
                                                  backend.gas_identity(), 0) &&
                    !portable::valid_vapor_mapping(mapped_parcel, material,
                                                   backend.gas_identity(), 1);
  double values[]{1, 0, 101850, 1, 0, 101850, 1, 0, 101850, 1, 0, 101850};
  portable::CompositionCellInput cell;
  cell.inventory = {static_cast<std::uint64_t>(rank), 1, 1, {}};
  cell.pressure_pa = 101325;
  cell.fields = values;
  cell.transport.mixing_time_s = 1;
  portable::CompositionInput input;
  input.duration_s = 0.5;
  input.field_count = 4;
  input.identity = &backend.closure_identity();
  input.gas = &backend;
  input.chemistry = &backend;
  input.cells = &cell;
  input.cell_count = 1;
  portable::CompositionWorkspace workspace(1, 1, 8, 2);
  auto r = workspace.prepare(MPI_COMM_WORLD, input);
  bool ok = mapping_ok && r.available && r.chemistry_calls == 8;
  if (ok)
    ok = std::abs(r.cells[0].fields[0] - 0.36787944117144233) < 1e-13 &&
         r.cells[0].fields[2] == 101850 && r.parcel_count == 0;
  const auto values_candidate = portable::prepare_accepted_values(
      MPI_COMM_WORLD, workspace, input, r, {1, 1, 1});
  ok &= values_candidate.available &&
        values_candidate.candidate.accepted_revision.accepted_step == 1;
  const auto restored = portable::restore_composition_values(
      values_candidate.candidate, backend.gas_identity(), {1, 1, 1}, 1, 1,
      nullptr);
  ok &= restored.available && restored.candidate.cells[0].fields[0] ==
                                  values_candidate.candidate.cells[0].fields[0];
  auto bad_rng = values_candidate.candidate;
  bad_rng.cells[0].random.accepted_step = 0;
  ok &= !portable::restore_composition_values(bad_rng, backend.gas_identity(),
                                              {1, 1, 1}, 1, 1, nullptr)
             .available;
  if (ranks > 1) {
    const auto unequal = portable::prepare_accepted_values(
        MPI_COMM_WORLD, workspace, input, r,
        {1, static_cast<std::uint64_t>(rank + 1), 1});
    ok &= !unequal.available;
  }
  FailingChemistry fault(backend);
  fault.fail_at = rank == 0 ? 4 : 0;
  input.chemistry = &fault;
  auto failure = workspace.prepare(MPI_COMM_WORLD, input);
  ok &= !failure.available && !failure.cells &&
        failure.lowest_failing_rank == 0 && failure.failure_module == 6 &&
        !workspace.current(r) && values[0] == 1 &&
        values_candidate.candidate.accepted_revision.accepted_step == 1;
  auto rejected = portable::prepare_accepted_values(MPI_COMM_WORLD, workspace,
                                                    input, failure, {1, 1, 1});
  ok &= !rejected.available && rejected.candidate.cells.empty();
  input.chemistry = &backend;
  r = workspace.prepare(MPI_COMM_WORLD, input);
  ok &= r.available &&
        r.cells[0].fields[0] == values_candidate.candidate.cells[0].fields[0];
  chemistry::detail::AnalyticIsomerBackend different_rate(3);
  input.chemistry = &different_rate;
  ok &= !workspace.prepare(MPI_COMM_WORLD, input).available;
  input.chemistry = &backend;
  BadGas bad(backend);
  input.gas = &bad;
  input.reaction_enabled = false;
  ok &= !workspace.prepare(MPI_COMM_WORLD, input).available;
  input.gas = &backend;
  input.reaction_enabled = true;
  if (ranks > 1) {
    input.start_time_s = rank == 0 ? 0 : 1;
    auto interval_mismatch = workspace.prepare(MPI_COMM_WORLD, input);
    ok &= !interval_mismatch.available;
    input.start_time_s = 0;
    chemistry::detail::AnalyticIsomerBackend other(3);
    input.gas = rank == 0 ? &backend : &other;
    input.chemistry = rank == 0 ? &backend : &other;
    input.identity =
        rank == 0 ? &backend.closure_identity() : &other.closure_identity();
    auto mechanism_mismatch = workspace.prepare(MPI_COMM_WORLD, input);
    ok &= !mechanism_mismatch.available;
    input.gas = &backend;
    input.chemistry = &backend;
    input.identity = &backend.closure_identity();
  }
  auto corrupted = values_candidate.candidate;
  corrupted.cells[0].fields[0] = -1;
  ok &= !portable::restore_composition_values(corrupted, backend.gas_identity(),
                                              {1, 1, 1}, 1, 1, nullptr)
             .available;
  spray::detail::ParcelMigrationPlan migration;
  MeshPatch patch;
  detail::make_mesh_patch(rank, ranks, {4, 4, 4}, patch);
  ok &= static_cast<bool>(
      migration.configure(MPI_COMM_WORLD, {4, 4, 4}, patch, 2, 8));
  spray::detail::CartesianParcelLocationProvider location(
      {{}, {}, {1, 1, 1}, {4, 4, 4}, ranks});
  input.migration = &migration;
  input.location = rank == 0 ? nullptr : &location;
  auto missing = workspace.prepare(MPI_COMM_WORLD, input);
  ok &= !missing.available && missing.failure_module == 7;
  input.location = &location;
  r = workspace.prepare(MPI_COMM_WORLD, input);
  ok &= r.available && r.parcel_count == 0;
  int local = ok ? 1 : 0, global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (!global)
    std::cerr << "portable no-parcel analytic chemistry failed module "
              << r.failure_module << " status " << int(r.status) << '\n';
  MPI_Finalize();
  return global ? 0 : 1;
}

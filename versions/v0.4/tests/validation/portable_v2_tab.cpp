// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09

// A bounded, prescribed zero-slip TAB/composition fixture, not a flow driver.
// One invocation selects one registered N/ranks configuration. Execution is
// gated by the coordinator after all implementation packages are gathered.
#include "mesh_focus_detail.hpp"
#include "models_chemistry_adapter_detail.hpp"
#include "models_portable_composition_detail.hpp"
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <string_view>
#include <vector>

namespace {
using namespace hundun::v04;
using namespace hundun::v04::spray;
using namespace hundun::v04::spray::detail;
namespace p = hundun::v04::portable;
constexpr double dt = 1e-6;
constexpr double mass = 4.1887902047863909846168578443726705e-10;
constexpr double surface_delta = 1.0471975511965977461542144610931676e-7;
constexpr ParcelId parent_id{8100, 42}, survivor_id{8100, 43};
int rank_id{}, rank_count{};
bool near(double a, double b, double absolute, double relative) {
  return std::isfinite(a) && std::isfinite(b) &&
         std::abs(a - b) <= absolute + relative * std::abs(b);
}
bool check(bool local, const char *message) {
  if (!local)
    std::cerr << "rank " << rank_id << " FAIL " << message << '\n';
  int in = local ? 1 : 0, out{};
  return MPI_Allreduce(&in, &out, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) ==
             MPI_SUCCESS &&
         out;
}
double sum(double value) {
  double total{};
  MPI_Allreduce(&value, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return total;
}

// Explicit analytic no-transfer interval. TAB is the real oscillator and
// representative-size closure; only the prescribed trajectory is analytic.
struct Ballistic final : ParcelIntervalProvider {
  ParcelIntervalReport advance(const SprayParcelState &old, double, double t,
                               ParcelPass,
                               p::Revision revision) const noexcept override {
    ParcelIntervalReport r;
    if (!std::isfinite(t) || t < 0)
      return r;
    r.revision = revision;
    r.parcel = old;
    for (unsigned d = 0; d < 3; ++d)
      r.parcel.position_m[d] += t * old.velocity_m_per_s[d];
    r.parcel.age_s += t;
    r.elapsed_duration_s = t;
    r.initial_liquid_absolute_enthalpy_j_per_kg = -100000;
    r.liquid_absolute_enthalpy_j_per_kg = -100000;
    r.exchange.available = true;
    r.available = true;
    return r;
  }
};
struct NoGeometryEvents final : ParcelEventGeometryProvider {
  ParcelEventQueryReport query(const SprayParcelState &,
                               const SprayParcelState &, double, double,
                               ParcelPass,
                               p::Revision revision) const noexcept override {
    ParcelEventQueryReport r;
    r.revision = revision;
    r.available = true;
    return r;
  }
};
struct FixedEnvironment final : ParcelTransferEnvironmentProvider {
  const LiquidPropertyService &liquid;
  explicit FixedEnvironment(const LiquidPropertyService &service)
      : liquid(service) {}
  ParcelTransferEnvironment
  sample(const SprayParcelState &, double, ParcelPass,
         p::Revision revision) const noexcept override {
    ParcelTransferEnvironment r;
    r.available = true;
    r.revision = revision;
    r.liquid = &liquid;
    r.liquid_absolute_enthalpy_j_per_kg = -100000;
    r.far_gas_density_kg_per_m3 = 1;
    r.gas.gas_velocity_m_per_s = {1, 0, 0};
    return r;
  }
};

bool same_parcel(const ParcelMigrationValue &a, const ParcelMigrationValue &b) {
  const auto &x = a.parcel;
  const auto &y = b.parcel;
  return x.id == y.id && x.position_m == y.position_m &&
         x.velocity_m_per_s == y.velocity_m_per_s && x.age_s == y.age_s &&
         x.droplet_mass_kg == y.droplet_mass_kg &&
         x.droplet_diameter_m == y.droplet_diameter_m &&
         x.multiplicity == y.multiplicity &&
         x.temperature_k == y.temperature_k &&
         x.liquid_material_fingerprint == y.liquid_material_fingerprint &&
         x.owner_global_cell == y.owner_global_cell &&
         a.tab_deformation == b.tab_deformation &&
         a.tab_deformation_rate_per_s == b.tab_deformation_rate_per_s &&
         a.breakup_ordinal == b.breakup_ordinal;
}
bool same(const p::CompositionSnapshot &a, const p::CompositionSnapshot &b) {
  if (a.accepted_revision != b.accepted_revision || a.version != b.version ||
      a.accepted_time_s != b.accepted_time_s ||
      a.field_count != b.field_count ||
      a.reaction_enabled != b.reaction_enabled ||
      a.parcel_rng_seed != b.parcel_rng_seed ||
      a.parcel_rng_algorithm_version != b.parcel_rng_algorithm_version ||
      a.cells.size() != b.cells.size() ||
      a.parcels.size() != b.parcels.size() || !a.injectors.empty() ||
      !b.injectors.empty())
    return false;
  for (std::size_t i = 0; i < a.parcels.size(); ++i)
    if (!same_parcel(a.parcels[i], b.parcels[i]))
      return false;
  for (std::size_t i = 0; i < a.cells.size(); ++i) {
    const auto &x = a.cells[i];
    const auto &y = b.cells[i];
    if (x.inventory.global_cell != y.inventory.global_cell ||
        x.inventory.volume_m3 != y.inventory.volume_m3 ||
        x.inventory.gas_mass_kg != y.inventory.gas_mass_kg ||
        x.inventory.gas_momentum_kg_m_per_s !=
            y.inventory.gas_momentum_kg_m_per_s ||
        x.pressure_pa != y.pressure_pa || x.fields != y.fields ||
        x.tcr_history.revision != y.tcr_history.revision ||
        x.tcr_history.initialized != y.tcr_history.initialized ||
        x.random.seed != y.random.seed ||
        x.random.accepted_step != y.random.accepted_step ||
        x.random.stochastic_stage != y.random.stochastic_stage ||
        x.random.field_pair != y.random.field_pair ||
        x.random.spatial_direction != y.random.spatial_direction ||
        x.random.purpose != y.random.purpose)
      return false;
  }
  return true;
}

struct Step {
  p::CompositionValueReport values;
  p::CompositionReport report;
};
Step advance(const p::CompositionSnapshot &old, MeshPatch patch,
             const LiquidAsset &asset,
             chemistry::detail::AnalyticIsomerBackend &gas,
             p::CompositionWorkspace &workspace) {
  Step out;
  CartesianParcelLocationProvider location(
      {old.accepted_revision, {0, 0, 0}, {1, 1, 1}, {4, 4, 4}, rank_count});
  ParcelMigrationPlan migration;
  if (!check(bool(migration.configure(MPI_COMM_WORLD, {4, 4, 4}, patch, 8, 32)),
             "bounded migration preparation"))
    return out;
  LiquidPropertyService liquid(&asset.pack, 1);
  FixedEnvironment environment(liquid);
  FixedTabEvolutionProvider evolution(environment);
  FixedTabEventBreakupProvider breakup(environment, 2);
  Ballistic interval;
  NoGeometryEvents geometry;
  std::vector<ParcelEventsInput> jobs;
  for (const auto &v : old.parcels) {
    ParcelEventsInput q;
    q.accepted_parcel = v.parcel;
    q.accepted_auxiliary = {v.tab_deformation, v.tab_deformation_rate_per_s,
                            v.breakup_ordinal};
    q.revision = old.accepted_revision;
    q.interval = &interval;
    q.geometry = &geometry;
    q.tab_evolution = &evolution;
    q.breakup = &breakup;
    q.duration_s = q.initial_substep_s = dt;
    q.maximum_segments = 16;
    q.maximum_events = 8;
    q.maximum_attempts = 128;
    jobs.push_back(q);
  }
  std::vector<std::size_t> vapor(jobs.size(), asset.vapor_species_index);
  std::vector<const LiquidAsset *> materials(jobs.size(), &asset);
  std::vector<p::CompositionCellInput> cells;
  for (const auto &v : old.cells) {
    p::CompositionCellInput q;
    q.inventory = v.inventory;
    q.pressure_pa = v.pressure_pa;
    q.fields = v.fields.data();
    q.tcr_history = v.tcr_history;
    q.tcr_trial.expected_revision = old.accepted_revision;
    q.tcr_trial.mode = tcr::detail::Mode::off;
    q.transport.mixing_time_s = 1;
    q.transport.random = v.random;
    cells.push_back(q);
  }
  p::CompositionInput q;
  q.revision = old.accepted_revision;
  q.start_time_s = old.accepted_time_s;
  q.duration_s = dt;
  q.field_count = old.field_count;
  q.identity = &gas.closure_identity();
  q.gas = &gas;
  q.chemistry = &gas;
  q.reaction_enabled = false;
  q.cells = cells.data();
  q.cell_count = cells.size();
  q.parcels = jobs.data();
  q.parcel_count = jobs.size();
  q.vapor_species_indices = vapor.data();
  q.parcel_materials = materials.data();
  q.parcel_rng_seed = old.parcel_rng_seed;
  q.migration = &migration;
  q.location = &location;
  out.report = workspace.prepare(MPI_COMM_WORLD, q);
  if (out.report.available)
    out.values = p::prepare_accepted_values(
        MPI_COMM_WORLD, workspace, q, out.report,
        {q.revision.accepted_step + 1, q.revision.input_revision + 1, 1});
  return out;
}

bool run(std::size_t fields, const char *path) {
  chemistry::detail::AnalyticIsomerBackend gas(2, false, 1000);
  auto loaded = load_liquid_asset(path, UINT64_C(6004043157121730787),
                                  gas.gas_identity());
  if (!check(loaded.available, "alpha asset complete identity"))
    return false;
  MeshPatch patch;
  if (!check(bool(hundun::v04::detail::make_mesh_patch(rank_id, rank_count,
                                                       {4, 4, 4}, patch)),
             "unique Cartesian cell owners"))
    return false;
  p::CompositionSnapshot initial;
  initial.accepted_revision = {0, 1, 1};
  initial.field_count = fields;
  initial.gas_identity = gas.gas_identity();
  initial.reaction_enabled = false;
  initial.parcel_rng_seed = 998877;
  for (int z = patch.begin.z; z < patch.begin.z + patch.cells.z; ++z)
    for (int y = patch.begin.y; y < patch.begin.y + patch.cells.y; ++y)
      for (int x = patch.begin.x; x < patch.begin.x + patch.cells.x; ++x) {
        p::CompositionCellValue cell;
        cell.inventory = {std::uint64_t(x + 4 * (y + 4 * z)), 1, 1, {0, 0, 0}};
        cell.pressure_pa = 118228.688979204018071428571429;
        cell.tcr_history.revision = initial.accepted_revision;
        cell.random.seed = 887766;
        for (std::size_t f = 0; f < fields; ++f)
          cell.fields.insert(cell.fields.end(), {.01, .99, 101000});
        initial.cells.push_back(cell);
      }
  CartesianParcelLocationProvider location(
      {initial.accepted_revision, {0, 0, 0}, {1, 1, 1}, {4, 4, 4}, rank_count});
  bool located_ok = true;
  for (int which = 0; which < 2; ++which) {
    const Vector3 position =
        which ? Vector3{3.25, 3.5, 3.5} : Vector3{.25, .5, .5};
    ParcelLocation owner;
    if (!location.locate(position, initial.accepted_revision, owner)) {
      located_ok = false;
      continue;
    }
    if (owner.owner_rank == rank_id) {
      SprayParcelState parcel{which ? survivor_id : parent_id,
                              position,
                              {1, 0, 0},
                              mass,
                              1e-4,
                              which ? 1. : 100.,
                              298.15,
                              loaded.asset.pack.material_fingerprint,
                              owner.global_cell,
                              .125};
      initial.parcels.push_back({parcel, which ? 0. : 1., 0, which ? 19U : 7U});
    }
  }
  if (!check(located_ok, "initial parcel ownership"))
    return false;
  p::CompositionWorkspace workspace(64, 8, 128, 2);
  auto first = advance(initial, patch, loaded.asset, gas, workspace);
  if (!check(first.report.available && first.values.available,
             "real TAB to common composition candidate"))
    return false;
  const auto &energy = first.report.breakup_energy;
  const double deformation = sum(energy.deformation_consumed_j);
  const double surface = sum(energy.surface_increase_j);
  const double bulk = sum(energy.bulk_kinetic_increase_j);
  const double residual = sum(energy.residual_j);
  const double events = sum(double(energy.event_count));
  bool valid =
      near(deformation, surface_delta, 1e-22, 5e-13) &&
      near(surface, surface_delta, 1e-22, 5e-13) && near(bulk, 0, 1e-22, 0) &&
      near(residual, 0, 1e-22, 0) && events == 1 &&
      first.report.model_id == "portable_source_transport_reaction_v1" &&
      first.report.chemistry_calls == 0 && !first.report.reaction_enabled;
  double local_children{}, local_survivors{}, local_mass{};
  for (const auto &v : first.values.candidate.parcels) {
    const auto &a = v.parcel;
    ParcelLocation owner;
    valid =
        valid &&
        bool(location.locate(a.position_m, initial.accepted_revision, owner)) &&
        owner.owner_rank == rank_id && owner.global_cell == a.owner_global_cell;
    valid = valid && !(a.id == parent_id) &&
            a.velocity_m_per_s == Vector3{1, 0, 0} &&
            near(a.age_s, .125001, 2e-16, 0) && a.temperature_k == 298.15 &&
            v.tab_deformation == 0 && v.tab_deformation_rate_per_s == 0;
    if (a.id == survivor_id) {
      ++local_survivors;
      valid = valid && v.breakup_ordinal == 19 &&
              near(a.position_m[0], 3.250001, 2e-15, 0);
    } else {
      ++local_children;
      local_mass += a.droplet_mass_kg * a.multiplicity;
      valid =
          valid && v.breakup_ordinal == 0 && a.owner_global_cell == 0 &&
          near(a.position_m[0], .250001, 2e-16, 0) &&
          near(a.droplet_diameter_m, 4.285714285714285714285714285714e-5, 0,
               5e-13) &&
          near(a.multiplicity, 635.185185185185185185185185185185, 0, 5e-13);
    }
  }
  const double children = sum(local_children), survivors = sum(local_survivors);
  const double child_mass = sum(local_mass);
  valid = valid && children == 2 && survivors == 1 &&
          near(child_mass, 4.1887902047863909846168578443726705e-8, 0, 5e-13);
  // A false transfer of the 1e-7 J deformation reservoir to gas must be seen:
  // the h tolerance here is much smaller than that energy / 1 kg cell mass.
  for (std::size_t i = 0; i < first.values.candidate.cells.size(); ++i) {
    const auto &a = first.values.candidate.cells[i];
    valid = valid && a.inventory.gas_mass_kg == 1 &&
            a.inventory.gas_momentum_kg_m_per_s == Vector3{0, 0, 0};
    for (std::size_t f = 0; f < fields; ++f)
      valid = valid && a.fields[f * 3] == .01 && a.fields[f * 3 + 1] == .99 &&
              near(a.fields[f * 3 + 2], 101000, 1e-10, 0);
  }
  if (!check(valid, "parent replacement, child remaining interval and separate "
                    "TAB energy"))
    return false;
  const auto accepted = first.values.candidate;
  CartesianParcelLocationProvider next_location({accepted.accepted_revision,
                                                 {0, 0, 0},
                                                 {1, 1, 1},
                                                 {4, 4, 4},
                                                 rank_count});
  auto restored = p::restore_composition_values(accepted, gas.gas_identity(),
                                                accepted.accepted_revision, 64,
                                                8, &next_location);
  if (!check(restored.available && same(accepted, restored.candidate),
             "ordinary-value restore preserves TAB, ordinal and clocks"))
    return false;
  auto continued = advance(accepted, patch, loaded.asset, gas, workspace);
  auto replay =
      advance(restored.candidate, patch, loaded.asset, gas, workspace);
  if (!check(continued.values.available && replay.values.available &&
                 continued.report.breakup_energy.event_count == 0 &&
                 replay.report.breakup_energy.event_count == 0 &&
                 same(continued.values.candidate, replay.values.candidate),
             "continuous and restored next interval agree without re-trigger"))
    return false;
  if (rank_id == 0)
    std::cout << std::setprecision(17) << "PASS portable_v2_tab N=" << fields
              << " ranks=" << rank_count << " children=" << children
              << " deformation_j=" << deformation << " surface_j=" << surface
              << " bulk_j=" << bulk << " residual_j=" << residual << '\n';
  return true;
}
} // namespace
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
  MPI_Comm_size(MPI_COMM_WORLD, &rank_count);
  const bool arguments =
      argc == 3 &&
      (std::string_view(argv[1]) == "2" || std::string_view(argv[1]) == "4") &&
      (rank_count == 1 || rank_count == 2 || rank_count == 4);
  const bool ok =
      check(arguments, "usage: portable_v2_tab N alpha.asset; ranks=1/2/4") &&
      run(std::string_view(argv[1]) == "2" ? 2 : 4, argv[2]);
  MPI_Finalize();
  return ok ? 0 : 1;
}

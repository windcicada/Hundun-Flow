// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09

// Prescribed gas / ordinary-value composition only, not a flow driver.
// One MPI invocation selects one manifest configuration. Do not run before
// P1--P8 are frozen and the coordinator opens the V2 execution gate.
#include "mesh_focus_detail.hpp"
#include "models_chemistry_adapter_detail.hpp"
#include "models_portable_composition_detail.hpp"
#include "models_spray_film_bridge_detail.hpp"
#include "models_spray_parcel_detail.hpp"
#include <algorithm>
#include <array>
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
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr std::uint64_t rng_seed = 224466;
int rank_id{}, rank_count{};
bool near(double a, double b, double absolute, double relative) {
  return std::isfinite(a) && std::isfinite(b) &&
         std::abs(a - b) <= absolute + relative * std::abs(b);
}
bool all(bool local) {
  int value = local ? 1 : 0, result = 0;
  return MPI_Allreduce(&value, &result, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) ==
             MPI_SUCCESS &&
         result != 0;
}
bool check_all(bool local, const char *what) {
  if (!local)
    std::cerr << "rank " << rank_id << " FAIL " << what << '\n';
  return all(local);
}
std::uint64_t global_cell(int x, int y, int z) {
  return std::uint64_t(x + 4 * (y + 4 * z));
}
double dot(Vector3 a, Vector3 b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

struct PrescribedGas final : ParcelGasStateProvider {
  std::uint64_t fingerprint{};
  double pressure_pa{100000}, enthalpy_j_per_kg{};
  bool fail{};
  mutable std::uint64_t queries{};
  ParcelGasSample sample(const SprayParcelState &, double elapsed, ParcelPass,
                         p::Revision revision, double *y,
                         std::size_t capacity) const noexcept override {
    ++queries;
    if (fail || capacity < 2 || !std::isfinite(elapsed) || elapsed < 0)
      return {};
    // Both packs place their explicitly mapped vapor at index zero; their
    // chemistry species meanings differ (alpha A/B, beta B/A).
    y[0] = .01;
    y[1] = .99;
    return {p::Status::success,
            revision,
            fingerprint,
            pressure_pa,
            enthalpy_j_per_kg,
            {0, 0, 0},
            2};
  }
};

// Geometry consumes only the resolved begin/end segment, independently of
// the interval model. No mesh storage or halo is fabricated or duplicated.
class CartesianEvents final : public ParcelEventGeometryProvider {
public:
  explicit CartesianEvents(double width) : width_(width) {}
  ParcelEventQueryReport query(const SprayParcelState &begin,
                               const SprayParcelState &end, double t0,
                               double t1, ParcelPass,
                               p::Revision revision) const noexcept override {
    ParcelEventQueryReport r;
    r.revision = revision;
    if (!(t1 >= t0) || begin.owner_global_cell >= 64)
      return r;
    int cell[3]{int(begin.owner_global_cell % 4),
                int((begin.owner_global_cell / 4) % 4),
                int(begin.owner_global_cell / 16)};
    for (unsigned d = 0; d < 3; ++d) {
      const double displacement = end.position_m[d] - begin.position_m[d];
      if (displacement == 0)
        continue;
      const int sign = displacement > 0 ? 1 : -1;
      const double face = (cell[d] + (sign > 0 ? 1 : 0)) * width_;
      const double fraction = (face - begin.position_m[d]) / displacement;
      if (fraction < 0 || fraction > 1)
        continue;
      int next[3]{cell[0], cell[1], cell[2]};
      next[d] += sign;
      ParcelEvent e;
      e.elapsed_time_s = t0 + (t1 - t0) * fraction;
      e.identity = 1 + begin.owner_global_cell * 6 + 2 * d + (sign > 0 ? 1 : 0);
      if (next[d] < 0 || next[d] >= 4)
        e.kind = ParcelEventKind::physical_outlet;
      else {
        e.kind = ParcelEventKind::internal_cell_crossing;
        e.next_global_cell = global_cell(next[0], next[1], next[2]);
      }
      r.events[r.count++] = e;
    }
    r.available = true;
    return r;
  }

private:
  double width_{};
};

// Explicit analytic INTERPHASE exchange fixture, not an A-S replacement.
// Two equal masses receive a prescribed drag acceleration -1 m/s2. Its work
// and impulse are exact interval integrals; no external force is introduced.
struct ConstantInterphaseAcceleration final : ParcelIntervalProvider {
  ParcelIntervalReport advance(const SprayParcelState &old, double, double dt,
                               ParcelPass,
                               p::Revision revision) const noexcept override {
    ParcelIntervalReport r;
    if (!std::isfinite(dt) || dt < 0)
      return r;
    r.revision = revision;
    r.elapsed_duration_s = dt;
    r.parcel = old;
    r.initial_liquid_absolute_enthalpy_j_per_kg = -100000;
    r.liquid_absolute_enthalpy_j_per_kg = -100000;
    for (unsigned d = 0; d < 3; ++d)
      r.parcel.position_m[d] += old.velocity_m_per_s[d] * dt;
    r.parcel.position_m[0] -= .5 * dt * dt;
    r.parcel.velocity_m_per_s[0] -= dt;
    r.parcel.age_s += dt;
    const double mass = old.droplet_mass_kg * old.multiplicity;
    auto &e = r.exchange;
    e.available = true;
    e.parcel_momentum_delta_kg_m_per_s[0] = -mass * dt;
    e.drag_momentum_to_parcel_kg_m_per_s = e.parcel_momentum_delta_kg_m_per_s;
    e.gas_momentum_delta_kg_m_per_s[0] = mass * dt;
    e.parcel_kinetic_energy_delta_j =
        .5 * mass *
        (dot(r.parcel.velocity_m_per_s, r.parcel.velocity_m_per_s) -
         dot(old.velocity_m_per_s, old.velocity_m_per_s));
    e.drag_work_to_parcel_j = e.parcel_kinetic_energy_delta_j;
    r.available = true;
    return r;
  }
};

bool same_parcel(const ParcelMigrationValue &a, const ParcelMigrationValue &b) {
  const auto &x = a.parcel;
  const auto &y = b.parcel;
  return x.id == y.id && x.position_m == y.position_m &&
         x.velocity_m_per_s == y.velocity_m_per_s &&
         x.droplet_mass_kg == y.droplet_mass_kg &&
         x.droplet_diameter_m == y.droplet_diameter_m &&
         x.multiplicity == y.multiplicity &&
         x.temperature_k == y.temperature_k &&
         x.liquid_material_fingerprint == y.liquid_material_fingerprint &&
         x.owner_global_cell == y.owner_global_cell && x.age_s == y.age_s &&
         a.tab_deformation == b.tab_deformation &&
         a.tab_deformation_rate_per_s == b.tab_deformation_rate_per_s;
}
bool same_snapshot(const p::CompositionSnapshot &a,
                   const p::CompositionSnapshot &b) {
  if (a.version != b.version || a.accepted_revision != b.accepted_revision ||
      a.field_count != b.field_count ||
      a.accepted_time_s != b.accepted_time_s ||
      a.reaction_enabled != b.reaction_enabled ||
      a.parcel_rng_seed != b.parcel_rng_seed ||
      a.parcel_rng_algorithm_version != b.parcel_rng_algorithm_version ||
      a.gas_identity.mechanism_sha256 != b.gas_identity.mechanism_sha256 ||
      a.gas_identity.phase != b.gas_identity.phase ||
      a.gas_identity.enthalpy_reference != b.gas_identity.enthalpy_reference ||
      a.gas_identity.species_names != b.gas_identity.species_names ||
      a.gas_identity.element_names != b.gas_identity.element_names ||
      a.gas_identity.element_counts != b.gas_identity.element_counts ||
      a.gas_identity.molecular_weights_kg_per_kmol !=
          b.gas_identity.molecular_weights_kg_per_kmol ||
      a.gas_identity.composition_fingerprint !=
          b.gas_identity.composition_fingerprint ||
      a.gas_identity.closure_fingerprint !=
          b.gas_identity.closure_fingerprint ||
      a.cells.size() != b.cells.size() ||
      a.parcels.size() != b.parcels.size() ||
      a.injectors.size() != b.injectors.size())
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
        x.tcr_history.branch_sign != y.tcr_history.branch_sign ||
        x.tcr_history.control != y.tcr_history.control ||
        x.tcr_history.initialized != y.tcr_history.initialized ||
        x.tcr_history.signed_root != y.tcr_history.signed_root ||
        x.tcr_history.input.eta != y.tcr_history.input.eta ||
        x.tcr_history.input.rate_ratio != y.tcr_history.input.rate_ratio ||
        x.tcr_history.mapping_identity != y.tcr_history.mapping_identity ||
        x.tcr_history.initialization_sign !=
            y.tcr_history.initialization_sign ||
        x.tcr_history.fold_count != y.tcr_history.fold_count ||
        x.tcr_history.last_fold.eta != y.tcr_history.last_fold.eta ||
        x.tcr_history.last_fold.rate_ratio !=
            y.tcr_history.last_fold.rate_ratio ||
        x.tcr_history.last_fold_base_revision !=
            y.tcr_history.last_fold_base_revision ||
        x.random.seed != y.random.seed ||
        x.random.accepted_step != y.random.accepted_step ||
        x.random.stochastic_stage != y.random.stochastic_stage ||
        x.random.field_pair != y.random.field_pair ||
        x.random.spatial_direction != y.random.spatial_direction ||
        x.random.purpose != y.random.purpose)
      return false;
  }
  for (std::size_t i = 0; i < a.parcels.size(); ++i)
    if (!same_parcel(a.parcels[i], b.parcels[i]))
      return false;
  for (std::size_t i = 0; i < a.injectors.size(); ++i) {
    const auto &x = a.injectors[i];
    const auto &y = b.injectors[i];
    if (x.accepted != y.accepted || x.spec.injector_id != y.spec.injector_id ||
        x.spec.seed != y.spec.seed || x.spec.origin_m != y.spec.origin_m ||
        x.spec.represented_mass_per_parcel_kg !=
            y.spec.represented_mass_per_parcel_kg ||
        x.spec.liquid_material_fingerprint !=
            y.spec.liquid_material_fingerprint ||
        x.spec.shape != y.spec.shape || x.spec.axis != y.spec.axis ||
        x.spec.cone_half_angle_rad != y.spec.cone_half_angle_rad ||
        x.spec.injection_speed_m_per_s != y.spec.injection_speed_m_per_s ||
        x.spec.mass_flow_rate_kg_per_s != y.spec.mass_flow_rate_kg_per_s ||
        x.spec.droplet_mass_kg != y.spec.droplet_mass_kg ||
        x.spec.droplet_diameter_m != y.spec.droplet_diameter_m ||
        x.spec.temperature_k != y.spec.temperature_k ||
        x.spec.owner_global_cell != y.spec.owner_global_cell)
      return false;
  }
  return true;
}

struct Fixture {
  bool beta{}, closed{};
  int fields{};
  double width{}, duration{};
  chemistry::detail::AnalyticIsomerBackend backend;
  LiquidAsset asset;
  MeshPatch patch{};
  p::CompositionSnapshot initial;
  std::array<DeterministicInjector, 3> injectors;
  std::vector<ParcelMigrationValue> injection;
  std::vector<InjectorValueSnapshot> injector_candidates;
  Fixture(bool use_beta, bool use_closed, int n)
      : beta(use_beta), closed(use_closed), fields(n),
        width(use_closed ? 1 : .01), duration(use_closed ? .25 : 1e-5),
        backend(use_beta ? 3 : 2, use_beta, use_beta ? 1200 : 1000) {}
  bool initialize(const char *asset_path) {
    auto loaded = load_liquid_asset(asset_path,
                                    beta ? UINT64_C(668675689539421851)
                                         : UINT64_C(6004043157121730787),
                                    backend.gas_identity());
    if (!check_all(loaded.available, "synthetic liquid asset identity"))
      return false;
    asset = std::move(loaded.asset);
    if (!check_all(bool(hundun::v04::detail::make_mesh_patch(
                       rank_id, rank_count, {4, 4, 4}, patch)),
                   "Cartesian ownership fixture"))
      return false;
    initial.accepted_revision = {0, 1, 1};
    initial.field_count = fields;
    initial.gas_identity = backend.gas_identity();
    initial.reaction_enabled = !closed;
    initial.parcel_rng_seed = rng_seed;
    const double h = closed ? 101000 : (beta ? 221220 : 102850);
    const double pressure = closed ? 118228.688979204018071428571429 : 100000;
    double y[]{.01, .99}, d[2]{}, sh[2]{}, w[2]{};
    p::GasQuery query{initial.accepted_revision,
                      backend.gas_identity().composition_fingerprint,
                      p::GasStateCoordinates::pressure_enthalpy,
                      pressure,
                      h,
                      0,
                      y,
                      2};
    p::GasQueryOutput gas{{}, d, sh, w, 2};
    if (!check_all(backend.query_gas(query, gas) == p::Status::success,
                   "initial gas EOS"))
      return false;
    const double volume = width * width * width;
    for (int z = patch.begin.z; z < patch.begin.z + patch.cells.z; ++z)
      for (int iy = patch.begin.y; iy < patch.begin.y + patch.cells.y; ++iy)
        for (int x = patch.begin.x; x < patch.begin.x + patch.cells.x; ++x) {
          p::CompositionCellValue c;
          c.inventory = {global_cell(x, iy, z),
                         volume,
                         closed ? 1 : gas.sample.density_kg_per_m3 * volume,
                         {0, 0, 0}};
          c.pressure_pa = pressure;
          c.tcr_history.revision = initial.accepted_revision;
          c.random.seed = rng_seed;
          c.random.stochastic_stage = 17;
          for (int f = 0; f < fields; ++f)
            c.fields.insert(c.fields.end(), {.01, .99, h});
          initial.cells.push_back(std::move(c));
        }
    CartesianParcelLocationProvider location({initial.accepted_revision,
                                              {0, 0, 0},
                                              {width, width, width},
                                              {4, 4, 4},
                                              rank_count});
    if (closed) {
      ParcelLocation located;
      bool ok = bool(
          location.locate({.25, .5, .5}, initial.accepted_revision, located));
      if (!check_all(ok, "closed-cell initial location"))
        return false;
      if (located.owner_rank == rank_id)
        for (unsigned i = 0; i < 2; ++i) {
          SprayParcelState parcel{{17, i + 1},
                                  {.25, .5, .5},
                                  {double(i + 1), 0, 0},
                                  .1,
                                  std::cbrt(.6 / (800 * pi)),
                                  1,
                                  298.15,
                                  asset.pack.material_fingerprint,
                                  0,
                                  0};
          initial.parcels.push_back({parcel, 0, 0});
        }
      return true;
    }
    bool ok = true;
    const double rho = beta ? 900 : 800, mass = rho * pi * 1e-12 / 6;
    for (unsigned i = 0; i < 3; ++i) {
      Vector3 origin{.5 * width, .5 * width, .5 * width};
      origin[i] = 2 * width - 1e-6;
      ParcelLocation located;
      if (!location.locate(origin, initial.accepted_revision, located)) {
        ok = false;
        continue;
      }
      if (located.owner_rank != rank_id)
        continue;
      InjectorSpec spec;
      spec.seed = rng_seed;
      spec.injector_id = 100 + i;
      spec.origin_m = origin;
      spec.axis = {0, 0, 0};
      spec.axis[i] = 1;
      spec.injection_speed_m_per_s = 1;
      spec.mass_flow_rate_kg_per_s = 1.25 * mass / duration;
      spec.represented_mass_per_parcel_kg = mass;
      spec.droplet_mass_kg = mass;
      spec.droplet_diameter_m = 1e-4;
      spec.temperature_k = 298.15;
      spec.liquid_material_fingerprint = asset.pack.material_fingerprint;
      spec.owner_global_cell = located.global_cell;
      auto &injector = injectors[i];
      if (!injector.reserve(1) || !injector.configure(spec)) {
        ok = false;
        continue;
      }
      initial.injectors.push_back({spec, injector.committed_state()});
      const auto report = injector.begin_trial(0, duration);
      SprayParcelState parcel;
      if (!report.succeeded() || report.parcel_count != 1 ||
          !injector.candidate_at(0, parcel) || !injector.preflight_commit()) {
        ok = false;
        continue;
      }
      injection.push_back({parcel, 0, 0});
      injector_candidates.push_back(
          {spec, {report.residual_mass_after_kg, report.next_ordinal}});
    }
    return check_all(ok,
                     "deterministic injection candidates without publication");
  }
};

struct StepResult {
  p::CompositionValueReport values;
  p::CompositionReport report;
  std::uint64_t sampler_queries{};
};
StepResult step(Fixture &fixture, const p::CompositionSnapshot &accepted,
                const std::vector<ParcelMigrationValue> &injected,
                const std::vector<InjectorValueSnapshot> &injector_values,
                p::CompositionWorkspace &workspace, bool inject_failure) {
  StepResult result;
  const auto revision = accepted.accepted_revision;
  CartesianParcelLocationProvider location(
      {revision,
       {0, 0, 0},
       {fixture.width, fixture.width, fixture.width},
       {4, 4, 4},
       rank_count});
  ParcelMigrationPlan migration;
  if (!check_all(bool(migration.configure(MPI_COMM_WORLD, {4, 4, 4},
                                          fixture.patch, 8, 32)),
                 "prepared migration capacity"))
    return result;
  PrescribedGas sampler;
  sampler.fingerprint = fixture.backend.gas_identity().composition_fingerprint;
  sampler.enthalpy_j_per_kg = fixture.beta ? 221220 : 102850;
  sampler.fail = inject_failure && rank_id == 0;
  FilmEnvironmentBridge bridge(fixture.asset, fixture.backend, sampler,
                               revision);
  FixedAsParcelIntervalProvider actual_interval(bridge);
  ConstantInterphaseAcceleration analytic_interval;
  const ParcelIntervalProvider *interval =
      fixture.closed
          ? static_cast<const ParcelIntervalProvider *>(&analytic_interval)
          : static_cast<const ParcelIntervalProvider *>(&actual_interval);
  CartesianEvents geometry(fixture.width);
  std::vector<ParcelEventsInput> jobs;
  std::vector<std::size_t> vapor;
  std::vector<const LiquidAsset *> materials;
  auto add_job = [&](const ParcelMigrationValue &parcel) {
    ParcelEventsInput q;
    q.accepted_parcel = parcel.parcel;
    q.accepted_auxiliary = {parcel.tab_deformation,
                            parcel.tab_deformation_rate_per_s,
                            parcel.breakup_ordinal};
    q.revision = revision;
    q.interval = interval;
    q.geometry = &geometry;
    q.duration_s = fixture.duration;
    q.initial_substep_s = fixture.duration;
    q.minimum_substep_s = 1e-12;
    q.relative_tolerance = 1e-6;
    q.maximum_segments = 64;
    q.maximum_events = 16;
    q.maximum_attempts = 512;
    jobs.push_back(q);
    vapor.push_back(fixture.asset.vapor_species_index);
    materials.push_back(&fixture.asset);
  };
  for (const auto &parcel : accepted.parcels)
    add_job(parcel);
  for (const auto &parcel : injected)
    add_job(parcel);
  std::vector<p::CompositionCellInput> cells;
  for (const auto &source : accepted.cells) {
    p::CompositionCellInput cell;
    cell.inventory = source.inventory;
    cell.pressure_pa = source.pressure_pa;
    cell.fields = source.fields.data();
    cell.tcr_history = source.tcr_history;
    cell.tcr_trial.expected_revision = revision;
    cell.tcr_trial.mode = tcr::detail::Mode::off;
    cell.transport.mixing_time_s = 1;
    cell.transport.random = source.random;
    cells.push_back(cell);
  }
  p::CompositionInput input;
  input.revision = revision;
  input.start_time_s = accepted.accepted_time_s;
  input.duration_s = fixture.duration;
  input.field_count = fixture.fields;
  input.identity = &fixture.backend.closure_identity();
  input.gas = &fixture.backend;
  input.chemistry = &fixture.backend;
  input.reaction_enabled = !fixture.closed;
  input.cells = cells.data();
  input.cell_count = cells.size();
  input.parcels = jobs.data();
  input.parcel_count = jobs.size();
  input.vapor_species_indices = vapor.data();
  input.parcel_materials = materials.data();
  input.candidate_injectors = injector_values.data();
  input.injector_count = injector_values.size();
  input.parcel_rng_seed = rng_seed;
  input.migration = &migration;
  input.location = &location;
  result.report = workspace.prepare(MPI_COMM_WORLD, input);
  result.sampler_queries = sampler.queries;
  if (result.report.available)
    result.values = p::prepare_accepted_values(
        MPI_COMM_WORLD, workspace, input, result.report,
        {revision.accepted_step + 1, revision.input_revision + 1, 1});
  return result;
}

std::array<double, 5>
budget(const p::CompositionSnapshot &snapshot, const LiquidAsset &asset,
       const std::vector<ParcelMigrationValue> &extra = {}) {
  std::array<double, 5> local{};
  for (const auto &cell : snapshot.cells) {
    const double m = cell.inventory.gas_mass_kg;
    double mean_h = 0;
    for (std::size_t f = 0; f < snapshot.field_count; ++f)
      mean_h += cell.fields[f * 3 + 2] / snapshot.field_count;
    local[0] += m;
    for (unsigned d = 0; d < 3; ++d)
      local[d + 1] += cell.inventory.gas_momentum_kg_m_per_s[d];
    local[4] += m * mean_h + dot(cell.inventory.gas_momentum_kg_m_per_s,
                                 cell.inventory.gas_momentum_kg_m_per_s) /
                                 (2 * m);
  }
  auto add = [&](const ParcelMigrationValue &value) {
    const auto &parcel = value.parcel;
    const auto h = evaluate_liquid_enthalpy(asset, parcel.temperature_k);
    const double m = parcel.droplet_mass_kg * parcel.multiplicity;
    local[0] += m;
    for (unsigned d = 0; d < 3; ++d)
      local[d + 1] += m * parcel.velocity_m_per_s[d];
    local[4] += m * h.liquid_enthalpy_j_per_kg +
                .5 * m * dot(parcel.velocity_m_per_s, parcel.velocity_m_per_s);
    if (!h.available)
      local[4] = NAN;
  };
  for (const auto &parcel : snapshot.parcels)
    add(parcel);
  for (const auto &parcel : extra)
    add(parcel);
  std::array<double, 5> global{};
  MPI_Allreduce(local.data(), global.data(), 5, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  return global;
}
bool common_source(const p::CompositionSnapshot &snapshot) {
  bool ok = true;
  for (const auto &cell : snapshot.cells)
    for (std::size_t f = 1; f < snapshot.field_count; ++f)
      for (unsigned s = 0; s < 3; ++s)
        ok &= cell.fields[s] == cell.fields[f * 3 + s];
  return check_all(
      ok, "one common physical source in every initially uniform field");
}
bool untouched_chemistry(const p::CompositionSnapshot &snapshot, bool beta,
                         bool second) {
  const double expected = beta ? (second ? .989940601781964360534593584864151
                                         : .989970300445495545033412299526002)
                               : (second ? .009999600007999893334399991466724
                                         : .009999800001999986666733333066668);
  bool ok = true;
  for (const auto &cell : snapshot.cells)
    if (cell.inventory.global_cell == 63) {
      ok &= near(cell.fields[beta ? 1 : 0], expected, 1e-14, 1e-12);
      ok &= near(cell.fields[2], beta ? 221220 : 102850, 1e-8, 1e-12);
    }
  return check_all(ok, "untouched cell independently follows exact chemistry "
                       "without parcel source");
}
bool budget_check(const std::array<double, 5> &before,
                  const std::array<double, 5> &after, bool closed) {
  bool ok = near(after[0], before[0], closed ? 1e-12 : 1e-18, 1e-12);
  for (unsigned d = 1; d < 4; ++d)
    ok &= near(after[d], before[d], closed ? 1e-13 : 1e-18, 1e-12);
  ok &=
      near(after[4], before[4], closed ? 5e-8 : 1e-11, closed ? 1e-14 : 1e-12);
  if (rank_id == 0) {
    std::cout << std::setprecision(17)
              << "mass_residual_kg=" << after[0] - before[0]
              << "\nenergy_residual_j=" << after[4] - before[4]
              << "\ntotal_mass_kg=" << after[0]
              << "\ntotal_h_plus_k_j=" << after[4] << '\n';
  }
  return check_all(ok, "global parcel+gas mass momentum H+K budget");
}
bool run_case(Fixture &fixture, bool faults) {
  p::CompositionWorkspace workspace(64, 8, 256, 2);
  const auto initial_copy = fixture.initial;
  const auto initial_budget =
      budget(fixture.initial, fixture.asset, fixture.injection);
  const double expected_initial_mass =
      fixture.closed
          ? 64.2
          : (fixture.beta ? .000053883428775835381969127471080501156
                          : .000053883271696202702479465547948331992);
  if (!check_all(near(initial_budget[0], expected_initial_mass,
                      fixture.closed ? 1e-12 : 1e-18, 1e-12),
                 "independent ideal-gas plus injected-liquid inventory"))
    return false;
  auto first = step(fixture, fixture.initial, fixture.injection,
                    fixture.injector_candidates, workspace, false);
  if (!check_all(first.report.available && first.values.available,
                 "first real composition accepted-value candidate"))
    return false;
  if (faults) {
    auto failed = step(fixture, fixture.initial, fixture.injection,
                       fixture.injector_candidates, workspace, true);
    if (!check_all(
            !failed.report.available && !failed.values.available &&
                !failed.report.cells && !failed.report.parcels &&
                failed.report.lowest_failing_rank == 0 &&
                failed.report.failure_module == 2 &&
                !workspace.current(first.report) &&
                same_snapshot(fixture.initial, initial_copy),
            "one-rank fault exposes no partial output or input mutation"))
      return false;
    auto retry = step(fixture, fixture.initial, fixture.injection,
                      fixture.injector_candidates, workspace, false);
    if (!check_all(
            retry.values.available &&
                same_snapshot(first.values.candidate, retry.values.candidate),
            "same-address retry equals clean candidate exactly"))
      return false;
  }
  bool unpublished = true;
  for (const auto &injector : fixture.injectors)
    if (injector.configured())
      unpublished &= injector.committed_state() == InjectorCommittedState{};
  if (!check_all(unpublished,
                 "stateful injectors have not been published by composition"))
    return false;
  auto accepted = std::move(first.values.candidate);
  std::uint64_t chemistry_calls = 0;
  MPI_Allreduce(&first.report.chemistry_calls, &chemistry_calls, 1,
                MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  if (!check_all(
          chemistry_calls ==
              std::uint64_t(fixture.closed ? 0 : 64 * 2 * fixture.fields),
          "explicit disabled chemistry or global 2N call topology"))
    return false;
  std::uint64_t local_count = accepted.parcels.size(), global_count = 0;
  MPI_Allreduce(&local_count, &global_count, 1, MPI_UINT64_T, MPI_SUM,
                MPI_COMM_WORLD);
  if (!check_all(global_count == std::uint64_t(fixture.closed ? 2 : 3),
                 "global physical parcel count after migration"))
    return false;
  bool clocks = accepted.parcel_rng_seed == rng_seed &&
                accepted.parcel_rng_algorithm_version == 1;
  for (const auto &cell : accepted.cells)
    clocks &=
        cell.random.seed == rng_seed && cell.random.stochastic_stage == 17 &&
        cell.random.accepted_step == accepted.accepted_revision.accepted_step;
  for (const auto &injector : accepted.injectors)
    clocks &=
        injector.accepted.next_ordinal == 1 &&
        near(injector.accepted.residual_mass_kg,
             .25 * injector.spec.represented_mass_per_parcel_kg, 1e-22, 1e-12);
  if (!check_all(clocks,
                 "accepted ESF/parcel RNG and injector reservoir clocks"))
    return false;
  if (!common_source(accepted) ||
      !budget_check(initial_budget, budget(accepted, fixture.asset),
                    fixture.closed))
    return false;
  if (fixture.closed) {
    bool exact = true;
    for (const auto &cell : accepted.cells)
      if (cell.inventory.global_cell == 0) {
        exact &=
            near(cell.inventory.gas_momentum_kg_m_per_s[0], .05, 1e-14, 1e-12);
        exact &= near(cell.fields[2], 101000.0675, 1e-8, 1e-13);
        exact &= near(dot(cell.inventory.gas_momentum_kg_m_per_s,
                          cell.inventory.gas_momentum_kg_m_per_s) /
                          (2 * cell.inventory.gas_mass_kg),
                      .00125, 1e-15, 1e-12);
      }
    return check_all(exact && first.report.chemistry_calls == 0,
                     "two-parcel aggregate nonlinear gas kinetic correction; "
                     "explicit reaction off");
  }
  if (!untouched_chemistry(accepted, fixture.beta, false))
    return false;
  // New owner and revision are checked using the prescribed partition, not
  // an imported product Restart encoding.
  CartesianParcelLocationProvider location(
      {accepted.accepted_revision,
       {0, 0, 0},
       {fixture.width, fixture.width, fixture.width},
       {4, 4, 4},
       rank_count});
  std::array<std::uint64_t, 10> origins{};
  origins[0] = fixture.injection.size();
  for (std::size_t i = 0; i < fixture.injection.size(); ++i) {
    origins[1 + i * 3] = fixture.injection[i].parcel.id.high;
    origins[2 + i * 3] = fixture.injection[i].parcel.id.low;
    origins[3 + i * 3] = std::uint64_t(rank_id);
  }
  std::array<std::uint64_t, 40> global_origins{};
  MPI_Allgather(origins.data(), 10, MPI_UINT64_T, global_origins.data(), 10,
                MPI_UINT64_T, MPI_COMM_WORLD);
  bool owners = true;
  std::uint64_t moved = 0;
  for (const auto &parcel : accepted.parcels) {
    ParcelLocation actual;
    owners &= bool(location.locate(parcel.parcel.position_m,
                                   accepted.accepted_revision, actual)) &&
              actual.owner_rank == rank_id &&
              actual.global_cell == parcel.parcel.owner_global_cell;
    bool found = false;
    for (int origin_rank = 0; origin_rank < rank_count; ++origin_rank)
      for (std::size_t i = 0; i < global_origins[origin_rank * 10]; ++i) {
        const auto *row = global_origins.data() + origin_rank * 10 + 1 + i * 3;
        if (parcel.parcel.id == ParcelId{row[0], row[1]}) {
          found = true;
          moved += row[2] != std::uint64_t(rank_id);
        }
      }
    owners &= found;
  }
  std::uint64_t global_moved = 0, queries = 0;
  MPI_Allreduce(&moved, &global_moved, 1, MPI_UINT64_T, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&first.sampler_queries, &queries, 1, MPI_UINT64_T, MPI_SUM,
                MPI_COMM_WORLD);
  if (!check_all(
          owners && (rank_count == 1 || global_moved > 0) && queries > 0,
          "actual cross-rank migration and full-composition film queries"))
    return false;
  if (rank_id == 0)
    std::cout << "first_step_sampler_queries=" << queries
              << "\nphysical_parcels_migrated=" << global_moved << '\n';
  auto restored = p::restore_composition_values(
      accepted, fixture.backend.gas_identity(), accepted.accepted_revision, 64,
      8, &location);
  if (!check_all(restored.available &&
                     same_snapshot(accepted, restored.candidate),
                 "full accepted-value snapshot restoration"))
    return false;
  const std::vector<ParcelMigrationValue> no_injection;
  auto continuous = step(fixture, accepted, no_injection, accepted.injectors,
                         workspace, false);
  auto resumed = step(fixture, restored.candidate, no_injection,
                      restored.candidate.injectors, workspace, false);
  if (!check_all(continuous.values.available && resumed.values.available &&
                     same_snapshot(continuous.values.candidate,
                                   resumed.values.candidate),
                 "continuous and restored next step agree exactly"))
    return false;
  return common_source(resumed.values.candidate) &&
         untouched_chemistry(resumed.values.candidate, fixture.beta, true) &&
         budget_check(initial_budget,
                      budget(resumed.values.candidate, fixture.asset), false);
}
} // namespace
int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);
  MPI_Comm_size(MPI_COMM_WORLD, &rank_count);
  const std::string_view mode = argc > 1 ? argv[1] : "",
                         pack = argc > 2 ? argv[2] : "",
                         number = argc > 3 ? argv[3] : "";
  const bool valid =
      argc == 5 && (mode == "spray" || mode == "closed" || mode == "fault") &&
      (pack == "alpha" || pack == "beta") && (number == "2" || number == "4") &&
      (rank_count == 1 || rank_count == 2 || rank_count == 4) &&
      (mode != "closed" || pack == "alpha");
  bool passed = false;
  if (check_all(valid, "select exactly one registered V2 configuration")) {
    Fixture fixture(pack == "beta", mode == "closed", number == "2" ? 2 : 4);
    if (fixture.initialize(argv[4]))
      passed = run_case(fixture, mode == "fault");
  }
  passed = all(passed);
  if (rank_id == 0)
    std::cout << (passed ? "V2_CASE_PASS" : "V2_CASE_FAIL") << '\n';
  MPI_Finalize();
  return passed ? 0 : 1;
}

// SPDX-License-Identifier: Apache-2.0
#include "models_spray_events_detail.hpp"
#include <cmath>
#include <iostream>
using namespace hundun::v04::spray;
using namespace hundun::v04::spray::detail;
namespace {
bool check(bool b, const char *text) {
  if (!b)
    std::cerr << text << '\n';
  return b;
}
struct Physics final : ParcelIntervalProvider {
  double evaporation{}, acceleration{};
  double spurious_kinetic_power{};
  bool first_order{}, bad_budget{};
  double failure_time{1e99};
  mutable std::size_t post_collision_queries{};
  ParcelIntervalReport
  advance(const SprayParcelState &s, double elapsed, double dt, ParcelPass,
          hundun::v04::portable::Revision revision) const noexcept override {
    ParcelIntervalReport r;
    r.revision = revision;
    r.parcel = s;
    r.available = true;
    if (elapsed >= failure_time)
      return {};
    r.elapsed_duration_s = dt;
    r.initial_liquid_absolute_enthalpy_j_per_kg = 100;
    r.liquid_absolute_enthalpy_j_per_kg = 100;
    if (s.velocity_m_per_s[0] < 0)
      ++post_collision_queries;
    if (evaporation > 0 && s.droplet_mass_kg / evaporation <= dt) {
      r.elapsed_duration_s = s.droplet_mass_kg / evaporation;
      r.complete_evaporation = true;
    }
    dt = r.elapsed_duration_s;
    r.parcel.droplet_mass_kg =
        r.complete_evaporation ? 0 : s.droplet_mass_kg - evaporation * dt;
    r.parcel.droplet_diameter_m =
        r.complete_evaporation
            ? 0
            : s.droplet_diameter_m *
                  std::cbrt(r.parcel.droplet_mass_kg / s.droplet_mass_kg);
    r.parcel.position_m[0] += s.velocity_m_per_s[0] * dt +
                              (first_order ? 0 : 0.5 * acceleration * dt * dt);
    r.parcel.velocity_m_per_s[0] += acceleration * dt;
    r.parcel.age_s += dt;
    auto &e = r.exchange;
    e.available = true;
    const double m0 = s.droplet_mass_kg * s.multiplicity,
                 m1 = r.parcel.droplet_mass_kg * s.multiplicity;
    e.parcel_liquid_mass_delta_kg = m1 - m0;
    e.vapor_mass_to_gas_kg = m0 - m1;
    e.gas_mass_delta_kg = m0 - m1;
    e.parcel_thermochemical_enthalpy_delta_j = (m1 - m0) * 100;
    e.thermal_exchange_to_gas_j = (m0 - m1) * 100;
    e.vapor_absolute_thermochemical_enthalpy_to_gas_j = (m0 - m1) * 100;
    if (bad_budget) {
      e.thermal_exchange_state_residual_j = 1;
      e.thermal_exchange_to_gas_j += 1;
    }
    for (int d = 0; d < 3; ++d) {
      e.parcel_momentum_delta_kg_m_per_s[d] =
          m1 * r.parcel.velocity_m_per_s[d] - m0 * s.velocity_m_per_s[d];
      e.gas_momentum_delta_kg_m_per_s[d] =
          -e.parcel_momentum_delta_kg_m_per_s[d];
      e.drag_momentum_to_parcel_kg_m_per_s[d] =
          e.parcel_momentum_delta_kg_m_per_s[d];
      const long double old_v = s.velocity_m_per_s[d],
                        new_v = r.parcel.velocity_m_per_s[d];
      e.parcel_kinetic_energy_delta_j +=
          static_cast<double>(.5L * (m1 - m0) * old_v * old_v +
                              .5L * m1 * (new_v - old_v) * (new_v + old_v));
    }
    e.parcel_kinetic_energy_delta_j += spurious_kinetic_power * dt;
    return r;
  }
};
struct EmptyGeometry final : ParcelEventGeometryProvider {
  ParcelEventQueryReport
  query(const SprayParcelState &, const SprayParcelState &, double, double,
        ParcelPass,
        hundun::v04::portable::Revision revision) const noexcept override {
    ParcelEventQueryReport r;
    r.available = true;
    r.revision = revision;
    return r;
  }
};
struct Geometry final : ParcelEventGeometryProvider {
  bool wall{}, outlet{}, crossing{}, breakup{};
  double boundary_time{.5};
  ParcelEventQueryReport
  query(const SprayParcelState &begin, const SprayParcelState &end,
        double start, double finish, ParcelPass,
        hundun::v04::portable::Revision revision) const noexcept override {
    ParcelEventQueryReport r;
    r.available = true;
    r.revision = revision;
    if (wall && begin.velocity_m_per_s[0] > 0 && begin.position_m[0] < 1 &&
        end.position_m[0] >= 1) {
      // Independent constant-a=1 wall root, not endpoint interpolation.
      const double tau =
          -begin.velocity_m_per_s[0] +
          std::sqrt(begin.velocity_m_per_s[0] * begin.velocity_m_per_s[0] +
                    2 * (1 - begin.position_m[0]));
      r.events[r.count++] = {ParcelEventKind::wall_collision,
                             start + tau,
                             7,
                             0,
                             {-1, 0, 0},
                             1,
                             false};
    }
    if (outlet && start < boundary_time && finish >= boundary_time)
      r.events[r.count++] = {
          ParcelEventKind::physical_outlet, boundary_time, 8, 0, {}, 1, false};
    if (crossing && begin.owner_global_cell == 0 && start < boundary_time &&
        finish >= boundary_time)
      r.events[r.count++] = {ParcelEventKind::internal_cell_crossing,
                             boundary_time,
                             9,
                             1,
                             {},
                             1,
                             false};
    if (breakup && start < boundary_time && finish >= boundary_time)
      r.events[r.count++] = {
          ParcelEventKind::breakup, boundary_time, 10, 0, {}, 1, false};
    return r;
  }
};
struct Children final : ParcelEventBreakupProvider {
  bool corrupt_count{};
  BreakupChildReport
  split(const SprayParcelState &s, const ParcelEvent &,
        const ParcelAuxiliaryState &, const TabBreakupReport *, ParcelPass,
        hundun::v04::portable::Revision revision) const noexcept override {
    BreakupChildReport r;
    r.status = BreakupChildStatus::success;
    r.parent_id = s.id;
    r.accepted_step = revision.accepted_step;
    r.candidate.available = true;
    r.candidate.child_count = 2;
    for (unsigned i = 0; i < 2; ++i) {
      auto &child = r.candidate.children[i];
      child = s;
      child.id = {s.id.high, s.id.low + i + 1};
      child.droplet_mass_kg /= 8;
      child.droplet_diameter_m /= 2;
      child.multiplicity *= 4;
    }
    if (corrupt_count)
      r.candidate.child_count = 100;
    return r;
  }
};
struct Environment final : ParcelTransferEnvironmentProvider {
  LiquidPropertyPack pack{};
  LiquidPropertyService liquid{&pack, 1};
  LiquidPropertyService other_liquid{&pack, 1};
  bool replace_final_liquid{};
  mutable std::size_t samples{};
  mutable std::array<Vector3, 3> sampled_positions{};
  Environment() {
    pack.material_fingerprint = 91;
    pack.minimum_temperature_k = 250;
    pack.maximum_temperature_k = 650;
    pack.density_kg_per_m3.c[0] = 750;
    pack.cp_j_per_kg_k.c[0] = 2200;
    pack.latent_heat_j_per_kg.c[0] = 250000;
    pack.surface_tension_n_per_m.c[0] = .025;
    pack.viscosity_pa_s.c[0] = .0008;
    pack.saturation_pressure.antoine_a = 4;
    pack.saturation_pressure.pressure_scale_pa = 1;
  }
  ParcelTransferEnvironment
  sample(const SprayParcelState &s, double, ParcelPass,
         hundun::v04::portable::Revision revision) const noexcept override {
    ++samples;
    ParcelTransferEnvironment r;
    r.available = true;
    r.revision = revision;
    r.liquid = &liquid;
    sampled_positions[(samples - 1) % 3] = s.position_m;
    r.far_gas_density_kg_per_m3 = 1.2;
    if (replace_final_liquid && samples % 3 == 0)
      r.liquid = &other_liquid;
    r.liquid_absolute_enthalpy_j_per_kg = 2200 * s.temperature_k;
    auto &g = r.gas;
    g.gas_temperature_k = 350;
    g.gas_vapor_mass_fraction = .1;
    g.thermodynamic_pressure_pa = 1e5;
    g.film_density_kg_per_m3 = .8;
    g.film_dynamic_viscosity_pa_s = 2e-5;
    g.film_thermal_conductivity_w_per_m_k = .05;
    g.film_vapor_diffusivity_m2_per_s = 2e-5;
    g.film_cp_j_per_kg_k = 1100;
    g.vapor_molecular_weight_kg_per_kmol = 28;
    g.carrier_molecular_weight_kg_per_kmol = 28;
    g.vapor_absolute_thermochemical_enthalpy_j_per_kg = 1020000;
    return r;
  }
};
SprayParcelState parcel() {
  SprayParcelState p;
  p.id = {1, 2};
  p.velocity_m_per_s = {2, 0, 0};
  p.droplet_mass_kg = 1;
  p.droplet_diameter_m = 1;
  p.multiplicity = 3;
  p.temperature_k = 300;
  p.liquid_material_fingerprint = 1;
  return p;
}
} // namespace
int main() {
  Physics physics;
  EmptyGeometry geometry;
  ParcelEventsInput in;
  in.accepted_parcel = parcel();
  in.interval = &physics;
  in.geometry = &geometry;
  in.duration_s = .5;
  in.initial_substep_s = .5;
  const auto r = integrate_parcel_events(in);
  bool ok = check(
      r.available && r.parcel.position_m[0] == 1 && r.predictor_passes == 1 &&
          r.corrector_passes == 1 && r.advanced_duration_s == .5,
      "ballistic endpoint and one predictor/corrector macro contract");
  Geometry events;
  events.wall = true;
  physics.acceleration = 1;
  in.geometry = &events;
  in.duration_s = 1;
  const auto rebound = integrate_parcel_events(in);
  if (!rebound.available)
    std::cerr << "rebound status " << static_cast<unsigned>(rebound.status)
              << '\n';
  ok &= check(
      rebound.available && rebound.event_count == 1 &&
          physics.post_collision_queries > 0 &&
          std::abs(rebound.parcel.position_m[0] + .196938456699069) < 1e-12 &&
          std::abs(rebound.parcel.velocity_m_per_s[0] + 1.898979485566356) <
              1e-12 &&
          std::abs(rebound.exchange.parcel_momentum_delta_kg_m_per_s[0] - 3) <
              1e-12,
      "wall event locates root, resamples force after reflection and separates "
      "wall impulse");
  events.wall = false;
  events.outlet = true;
  physics.acceleration = 0;
  const auto outlet = integrate_parcel_events(in);
  ok &=
      check(outlet.available && outlet.physical_outlet &&
                outlet.parent_removed && outlet.outlet_inventory.mass_kg == 3 &&
                outlet.outlet_inventory.momentum_kg_m_per_s[0] == 6 &&
                outlet.outlet_inventory.thermochemical_enthalpy_j == 300 &&
                outlet.outlet_inventory.kinetic_energy_j == 6,
            "physical outlet carries exactly one complete m/P/H/K inventory");
  physics.evaporation = 2;
  const auto vanished = integrate_parcel_events(in);
  ok &= check(
      vanished.available && vanished.complete_evaporation &&
          !vanished.physical_outlet && vanished.parcel.droplet_mass_kg == 0 &&
          vanished.advanced_duration_s == .5 &&
          vanished.exchange.gas_mass_delta_kg == 3 &&
          vanished.outlet_inventory.mass_kg == 0 && vanished.event_count == 2 &&
          vanished.events[0].kind == ParcelEventKind::complete_evaporation &&
          vanished.events[0].applied && !vanished.events[1].applied,
      "coincident evaporation wins outlet, reports both events and never "
      "double counts inventory");
  physics.evaporation = 0;
  events.outlet = false;
  events.crossing = true;
  const auto crossing = integrate_parcel_events(in);
  ok &= check(crossing.available && crossing.parcel.owner_global_cell == 1 &&
                  crossing.segment_count == 2 &&
                  crossing.segments[0].global_cell == 0 &&
                  crossing.segments[1].global_cell == 1,
              "cross-cell event tags interval exchanges with each source cell "
              "exactly once");
  events.crossing = false;
  physics.acceleration = 1;
  physics.first_order = true;
  in.relative_tolerance = 0;
  in.position_absolute_tolerance_m = 1e-4;
  in.velocity_absolute_tolerance_m_per_s = 1e-8;
  const auto refined = integrate_parcel_events(in);
  ok &= check(refined.available && refined.rejected_substeps > 0 &&
                  std::abs(refined.parcel.position_m[0] - 2.5) < .01,
              "whole versus two half steps refines first-order state error "
              "toward analytic x=2.5");
  physics.bad_budget = true;
  in.minimum_substep_s = .05;
  const auto budget_failure = integrate_parcel_events(in);
  ok &= check(!budget_failure.available && budget_failure.segment_count == 0 &&
                  !budget_failure.exchange.available,
              "thermal budget failure at minimum dt withdraws every segment "
              "and exchange");
  physics.bad_budget = false;
  physics.first_order = false;
  in.maximum_segments = 1;
  const auto capacity = integrate_parcel_events(in);
  ok &= check(!capacity.available && capacity.segment_count == 0,
              "capacity exhaustion has no partial trajectory");
  in.maximum_segments = 256;
  physics.failure_time = .5;
  const auto failed = integrate_parcel_events(in);
  ok &= check(
      !failed.available && failed.segment_count == 0 &&
          !failed.exchange.available,
      "provider failure after an earlier interval publishes no partial output");
  physics.failure_time = 1e99;
  physics.spurious_kinetic_power = 1;
  const auto false_kinetic = integrate_parcel_events(in);
  ok &= check(!false_kinetic.available && !false_kinetic.exchange.available,
              "unchanged endpoint kinetic energy cannot admit a fabricated "
              "kinetic source");
  physics.spurious_kinetic_power = 0;
  Children children;
  physics.failure_time = 1e99;
  physics.acceleration = 0;
  events.breakup = true;
  in.breakup = &children;
  const auto split = integrate_parcel_events(in);
  ok &= check(split.available && split.breakup_requested &&
                  split.children.child_count == 2 &&
                  split.advanced_duration_s == .5 &&
                  split.children_remaining_duration_s == .5,
              "breakup returns child candidates and explicit unintegrated "
              "remaining interval");
  children.corrupt_count = true;
  ok &= check(!integrate_parcel_events(in).available,
              "corrupt breakup child count is rejected before candidate arrays "
              "can be consumed");
  children.corrupt_count = false;
  events.breakup = false;
  in.accepted_parcel = split.children.children[0];
  in.elapsed_offset_s = split.advanced_duration_s;
  in.duration_s = split.children_remaining_duration_s;
  const auto continued = integrate_parcel_events(in);
  ok &= check(continued.available && continued.advanced_duration_s == .5 &&
                  continued.segments[0].begin_time_s == .5 &&
                  continued.segments[continued.segment_count - 1].end_time_s ==
                      1 &&
                  continued.parcel.position_m[0] == 2,
              "child continuation consumes remaining time using the common "
              "accepted-step time origin");
  Environment environment;
  FixedAsParcelIntervalProvider physical(environment);
  ParcelEventsInput physical_input;
  physical_input.interval = &physical;
  physical_input.geometry = &geometry;
  physical_input.accepted_parcel = parcel();
  physical_input.accepted_parcel.liquid_material_fingerprint = 91;
  physical_input.accepted_parcel.droplet_diameter_m = 1e-4;
  physical_input.accepted_parcel.droplet_mass_kg = 3.9269908169872415e-10;
  physical_input.accepted_parcel.temperature_k = 350;
  physical_input.accepted_parcel.velocity_m_per_s = {1e-6, 0, 0};
  physical_input.duration_s = 1e-4;
  physical_input.initial_substep_s = 1e-4;
  const auto physical_report = integrate_parcel_events(physical_input);
  if (!physical_report.available)
    std::cerr << "physical status "
              << static_cast<unsigned>(physical_report.status) << '\n';
  ok &= check(
      physical_report.available && environment.samples > 6 &&
          physical_report.parcel.droplet_mass_kg ==
              physical_input.accepted_parcel.droplet_mass_kg &&
          std::abs(physical_report.parcel.velocity_m_per_s[0] -
                   9.952115015903097e-7) < 1e-10 &&
          std::abs(physical_report.exchange.thermal_exchange_state_residual_j) <
              1e-16,
      "existing A-S/SN physical adapter reaches saturated Stokes limit with "
      "fresh gas queries");
  environment.samples = 0;
  const auto spatial = physical.advance(physical_input.accepted_parcel, 0, 1e-4,
                                        ParcelPass::corrector, {});
  ok &= check(spatial.available && environment.sampled_positions[1][0] > 0 &&
                  environment.sampled_positions[2] == spatial.parcel.position_m,
              "predicted and corrected gas queries sample their advected "
              "positions, not the segment origin");
  environment.samples = 0;
  environment.replace_final_liquid = true;
  ok &= check(
      !physical
           .advance(physical_input.accepted_parcel, 0, 1e-4,
                    ParcelPass::corrector, {})
           .available,
      "liquid service identity cannot change only at the corrected endpoint");
  environment.replace_final_liquid = false;
  physical_input.accepted_parcel.velocity_m_per_s = {0, 0, 0};
  physical_input.duration_s = 1e-6;
  physical_input.initial_substep_s = 1e-6;
  physical_input.accepted_auxiliary = {.2, 0, 7};
  const auto tab_off = integrate_parcel_events(physical_input);
  ok &= check(tab_off.available && !tab_off.tab_evolved &&
                  tab_off.auxiliary.tab_deformation == .2 &&
                  tab_off.auxiliary.breakup_ordinal == 7,
              "TAB off preserves accepted history without claiming evolution");
  FixedTabEvolutionProvider tab_evolution(environment);
  physical_input.tab_evolution = &tab_evolution;
  const auto tab_on = integrate_parcel_events(physical_input);
  // Constant-coefficient damped oscillator: y0=.2, v0=0, k=2.133333333e9,
  // damping=2133.333333/s. Independent analytic reference at 1 us.
  ok &= check(tab_on.available && tab_on.tab_evolved &&
                  tab_on.auxiliary.breakup_ordinal == 7 &&
                  std::abs(tab_on.auxiliary.tab_deformation -
                           .19978685618038158) < 1e-12 &&
                  std::abs(tab_on.auxiliary.tab_deformation_rate_per_s +
                           426.06035319131558) < 1e-8,
              "fresh-physics TAB evolution returns the damped-oscillator "
              "terminal history");
  physical_input.initial_substep_s = 5e-7;
  const auto tab_refined = integrate_parcel_events(physical_input);
  ok &= check(tab_refined.available &&
                  std::abs(tab_refined.auxiliary.tab_deformation -
                           tab_on.auxiliary.tab_deformation) < 1e-12,
              "piecewise constant TAB coefficients have the exact "
              "constant-state refinement limit");
  FixedTabEventBreakupProvider tab_children(environment, 2);
  physical_input.breakup = &tab_children;
  physical_input.accepted_auxiliary = {1, 0, 7};
  const auto dynamic_split = integrate_parcel_events(physical_input);
  if (!dynamic_split.available)
    std::cerr << "dynamic TAB status "
              << static_cast<unsigned>(dynamic_split.status) << '\n';
  ok &= check(dynamic_split.available && dynamic_split.breakup_requested &&
                  dynamic_split.advanced_duration_s == 0 &&
                  dynamic_split.children_remaining_duration_s ==
                      physical_input.duration_s &&
                  dynamic_split.segment_count == 0 &&
                  dynamic_split.children.child_count == 2,
              "accepted trigger feeds genuine TAB representative children "
              "without losing remaining time");
  ok &= check(
      dynamic_split.breakup_budget.supplied_deformation_energy_j > 0 &&
          std::abs(
              dynamic_split.breakup_budget.unassigned_deformation_energy_j) <
              1e-14,
      "TAB budget survives the event boundary independently of gas exchange");
  physical_input.accepted_auxiliary = {.99, 20000, 7};
  physical_input.initial_substep_s = 1e-6;
  const auto evolved_split = integrate_parcel_events(physical_input);
  if (!evolved_split.available)
    std::cerr << "evolved TAB status "
              << static_cast<unsigned>(evolved_split.status) << '\n';
  ok &= check(
      evolved_split.available && evolved_split.breakup_requested &&
          evolved_split.advanced_duration_s > 0 &&
          evolved_split.advanced_duration_s < physical_input.duration_s &&
          std::abs(evolved_split.auxiliary.tab_deformation - 1) < 1e-10 &&
          evolved_split.children_remaining_duration_s ==
              physical_input.duration_s - evolved_split.advanced_duration_s,
      "dynamically located TAB threshold carries trial history directly into "
      "conservative children");
  Physics boosted;
  boosted.acceleration = 0x1p-23;
  ParcelEventsInput boost;
  boost.accepted_parcel = parcel();
  boost.accepted_parcel.velocity_m_per_s = {1e8, 0, 0};
  boost.interval = &boosted;
  boost.geometry = &geometry;
  boost.duration_s = 1;
  boost.initial_substep_s = 1;
  const auto large_background = integrate_parcel_events(boost);
  // m=3 kg, u0=1e8 m/s, du=2^-23 m/s. The kinetic increment is
  // 35.762786865234396 J despite endpoint energies of order 1e16 J.
  ok &= check(
      large_background.available &&
          std::abs(large_background.exchange.parcel_kinetic_energy_delta_j -
                   35.762786865234396) < 1e-11,
      "kinetic consistency retains a small physical increment on a large "
      "velocity background");
  return ok ? 0 : 1;
}

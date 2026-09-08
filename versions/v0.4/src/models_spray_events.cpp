// SPDX-License-Identifier: Apache-2.0
#include "models_spray_events_detail.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace hundun::v04::spray::detail {
namespace {
bool finite3(const Vector3 &v) {
  return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
}
double dot(const Vector3 &a, const Vector3 &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
double kinetic_delta(const SprayParcelState &before,
                     const SprayParcelState &after) noexcept {
  const long double old_mass =
      static_cast<long double>(before.droplet_mass_kg) * before.multiplicity;
  const long double new_mass =
      static_cast<long double>(after.droplet_mass_kg) * after.multiplicity;
  long double delta = 0;
  for (std::size_t d = 0; d < 3; ++d) {
    const long double old_v = before.velocity_m_per_s[d],
                      new_v = after.velocity_m_per_s[d];
    // K1-K0 = (m1-m0)u0^2/2 + m1(u1-u0)(u1+u0)/2.
    // Do not subtract two rounded, large-background endpoint energies.
    delta += .5L * (new_mass - old_mass) * old_v * old_v +
             .5L * new_mass * (new_v - old_v) * (new_v + old_v);
  }
  return static_cast<double>(delta);
}
bool near(double a, double b, double absolute, double relative) {
  return std::isfinite(a) && std::isfinite(b) &&
         std::abs(a - b) <=
             absolute + relative * std::max(std::abs(a), std::abs(b));
}
void add(TransferExchangeCandidate &a, const TransferExchangeCandidate &b) {
#define ADD(field) a.field += b.field
  ADD(parcel_liquid_mass_delta_kg);
  ADD(parcel_thermochemical_enthalpy_delta_j);
  ADD(parcel_kinetic_energy_delta_j);
  ADD(vapor_mass_to_gas_kg);
  ADD(vapor_absolute_thermochemical_enthalpy_to_gas_j);
  ADD(convective_heat_to_parcel_j);
  ADD(drag_work_to_parcel_j);
  ADD(gas_mass_delta_kg);
  ADD(mass_closure_residual_kg);
  ADD(thermal_exchange_to_gas_j);
  ADD(thermal_exchange_state_residual_j);
#undef ADD
  for (std::size_t d = 0; d < 3; ++d) {
    a.parcel_momentum_delta_kg_m_per_s[d] +=
        b.parcel_momentum_delta_kg_m_per_s[d];
    a.vapor_momentum_to_gas_kg_m_per_s[d] +=
        b.vapor_momentum_to_gas_kg_m_per_s[d];
    a.drag_momentum_to_parcel_kg_m_per_s[d] +=
        b.drag_momentum_to_parcel_kg_m_per_s[d];
    a.gas_momentum_delta_kg_m_per_s[d] += b.gas_momentum_delta_kg_m_per_s[d];
    a.momentum_closure_residual_kg_m_per_s[d] +=
        b.momentum_closure_residual_kg_m_per_s[d];
    a.momentum_quadrature_residual_kg_m_per_s[d] +=
        b.momentum_quadrature_residual_kg_m_per_s[d];
  }
  a.available = true;
}
bool valid_interval(const ParcelIntervalReport &r, double dt,
                    portable::Revision revision) {
  if (!r.available || !r.exchange.available || r.revision != revision ||
      !std::isfinite(r.elapsed_duration_s) || r.elapsed_duration_s < 0 ||
      r.elapsed_duration_s > dt ||
      !std::isfinite(r.liquid_absolute_enthalpy_j_per_kg))
    return false;
  const auto &e = r.exchange;
  const double scalars[]{e.parcel_liquid_mass_delta_kg,
                         e.parcel_thermochemical_enthalpy_delta_j,
                         e.parcel_kinetic_energy_delta_j,
                         e.vapor_mass_to_gas_kg,
                         e.vapor_absolute_thermochemical_enthalpy_to_gas_j,
                         e.convective_heat_to_parcel_j,
                         e.drag_work_to_parcel_j,
                         e.gas_mass_delta_kg,
                         e.mass_closure_residual_kg,
                         e.thermal_exchange_to_gas_j,
                         e.thermal_exchange_state_residual_j};
  for (double v : scalars)
    if (!std::isfinite(v))
      return false;
  if (!finite3(e.parcel_momentum_delta_kg_m_per_s) ||
      !finite3(e.vapor_momentum_to_gas_kg_m_per_s) ||
      !finite3(e.drag_momentum_to_parcel_kg_m_per_s) ||
      !finite3(e.gas_momentum_delta_kg_m_per_s) ||
      !finite3(e.momentum_closure_residual_kg_m_per_s) ||
      !finite3(e.momentum_quadrature_residual_kg_m_per_s))
    return false;
  if (r.complete_evaporation)
    return r.parcel.droplet_mass_kg == 0 && r.parcel.droplet_diameter_m == 0 &&
           finite3(r.parcel.position_m) && finite3(r.parcel.velocity_m_per_s);
  return r.elapsed_duration_s == dt &&
         validate_parcel_state(r.parcel) == ParcelStateStatus::success;
}
bool accurate(const ParcelIntervalReport &a, const ParcelIntervalReport &b,
              const ParcelEventsInput &in) {
  if (a.complete_evaporation != b.complete_evaporation ||
      !near(a.elapsed_duration_s, b.elapsed_duration_s,
            in.event_time_tolerance_s, in.relative_tolerance) ||
      !near(a.parcel.droplet_mass_kg, b.parcel.droplet_mass_kg,
            in.mass_absolute_tolerance_kg, in.relative_tolerance) ||
      !near(a.parcel.temperature_k, b.parcel.temperature_k,
            in.temperature_absolute_tolerance_k, in.relative_tolerance))
    return false;
  for (std::size_t d = 0; d < 3; ++d)
    if (!near(a.parcel.position_m[d], b.parcel.position_m[d],
              in.position_absolute_tolerance_m, in.relative_tolerance) ||
        !near(a.parcel.velocity_m_per_s[d], b.parcel.velocity_m_per_s[d],
              in.velocity_absolute_tolerance_m_per_s, in.relative_tolerance) ||
        !near(a.exchange.parcel_momentum_delta_kg_m_per_s[d],
              b.exchange.parcel_momentum_delta_kg_m_per_s[d],
              in.momentum_absolute_tolerance_kg_m_per_s,
              in.relative_tolerance) ||
        !near(b.exchange.momentum_quadrature_residual_kg_m_per_s[d], 0,
              in.momentum_absolute_tolerance_kg_m_per_s +
                  in.relative_tolerance *
                      std::abs(b.exchange.parcel_momentum_delta_kg_m_per_s[d]),
              0))
      return false;
  return near(a.exchange.parcel_thermochemical_enthalpy_delta_j,
              b.exchange.parcel_thermochemical_enthalpy_delta_j,
              in.energy_absolute_tolerance_j, in.relative_tolerance) &&
         near(a.exchange.parcel_kinetic_energy_delta_j,
              b.exchange.parcel_kinetic_energy_delta_j,
              in.energy_absolute_tolerance_j, in.relative_tolerance) &&
         near(
             b.exchange.thermal_exchange_state_residual_j, 0,
             in.energy_absolute_tolerance_j +
                 in.relative_tolerance *
                     std::max(
                         std::abs(
                             b.exchange.parcel_thermochemical_enthalpy_delta_j),
                         std::abs(b.exchange.thermal_exchange_to_gas_j)),
             0);
}
bool consistent(const SprayParcelState &begin, const ParcelIntervalReport &r,
                const ParcelEventsInput &in) {
  const auto &end = r.parcel;
  const auto &e = r.exchange;
  if (end.id != begin.id || end.multiplicity != begin.multiplicity ||
      end.liquid_material_fingerprint != begin.liquid_material_fingerprint ||
      !std::isfinite(r.initial_liquid_absolute_enthalpy_j_per_kg))
    return false;
  const double old_mass = begin.droplet_mass_kg * begin.multiplicity,
               new_mass = end.droplet_mass_kg * end.multiplicity;
  if (!near(e.parcel_liquid_mass_delta_kg, new_mass - old_mass,
            in.mass_absolute_tolerance_kg, in.relative_tolerance) ||
      !near(e.parcel_kinetic_energy_delta_j, kinetic_delta(begin, end),
            in.energy_absolute_tolerance_j, in.relative_tolerance) ||
      !near(e.parcel_thermochemical_enthalpy_delta_j,
            new_mass * r.liquid_absolute_enthalpy_j_per_kg -
                old_mass * r.initial_liquid_absolute_enthalpy_j_per_kg,
            in.energy_absolute_tolerance_j, in.relative_tolerance) ||
      !near(e.thermal_exchange_state_residual_j,
            e.parcel_thermochemical_enthalpy_delta_j +
                e.thermal_exchange_to_gas_j,
            in.energy_absolute_tolerance_j, in.relative_tolerance))
    return false;
  for (std::size_t d = 0; d < 3; ++d)
    if (!near(e.parcel_momentum_delta_kg_m_per_s[d],
              new_mass * end.velocity_m_per_s[d] -
                  old_mass * begin.velocity_m_per_s[d],
              in.momentum_absolute_tolerance_kg_m_per_s, in.relative_tolerance))
      return false;
  return true;
}
ParcelEventsReport fail(ParcelEventsStatus status,
                        portable::Revision revision) {
  ParcelEventsReport r;
  r.status = status;
  r.revision = revision;
  return r;
}
bool valid_children(const BreakupChildReport &r, const SprayParcelState &parent,
                    portable::Revision revision) {
  if (!r.succeeded() || r.parent_id != parent.id ||
      r.accepted_step != revision.accepted_step ||
      r.candidate.child_count < 2 ||
      r.candidate.child_count > kMaximumBreakupChildParcels)
    return false;
  long double mass = 0;
  std::array<long double, 3> momentum{};
  for (std::size_t i = 0; i < r.candidate.child_count; ++i) {
    const auto &c = r.candidate.children[i];
    if (validate_parcel_state(c) != ParcelStateStatus::success ||
        c.id == parent.id ||
        c.liquid_material_fingerprint != parent.liquid_material_fingerprint ||
        c.position_m != parent.position_m ||
        c.owner_global_cell != parent.owner_global_cell)
      return false;
    for (std::size_t j = 0; j < i; ++j)
      if (c.id == r.candidate.children[j].id)
        return false;
    const long double m =
        static_cast<long double>(c.droplet_mass_kg) * c.multiplicity;
    mass += m;
    for (std::size_t d = 0; d < 3; ++d)
      momentum[d] += m * c.velocity_m_per_s[d];
  }
  const double m = parent.droplet_mass_kg * parent.multiplicity;
  if (!near(static_cast<double>(mass), m, 1e-24, 1e-12))
    return false;
  for (std::size_t d = 0; d < 3; ++d)
    if (!near(static_cast<double>(momentum[d]), m * parent.velocity_m_per_s[d],
              1e-22, 1e-12))
      return false;
  return true;
}
bool apply_wall(SprayParcelState &state, const ParcelEvent &event,
                portable::ExchangeDelta &ledger) noexcept {
  const double norm = std::sqrt(dot(event.wall_normal, event.wall_normal));
  if (!finite3(event.wall_normal) || !std::isfinite(norm) || norm <= 0 ||
      !std::isfinite(event.restitution) || event.restitution < 0 ||
      event.restitution > 1)
    return false;
  Vector3 normal = event.wall_normal;
  for (double &v : normal)
    v /= norm;
  const Vector3 old = state.velocity_m_per_s;
  const double vn = dot(old, normal),
               mass = state.droplet_mass_kg * state.multiplicity;
  for (std::size_t d = 0; d < 3; ++d) {
    state.velocity_m_per_s[d] -= (1 + event.restitution) * vn * normal[d];
    ledger.momentum_kg_m_per_s[d] +=
        mass * (old[d] - state.velocity_m_per_s[d]);
  }
  ledger.kinetic_energy_j +=
      .5 * mass *
      (dot(old, old) - dot(state.velocity_m_per_s, state.velocity_m_per_s));
  return true;
}
void outlet_inventory(const SprayParcelState &state, double enthalpy,
                      portable::ExchangeDelta &inventory) noexcept {
  inventory.mass_kg = state.droplet_mass_kg * state.multiplicity;
  for (std::size_t d = 0; d < 3; ++d)
    inventory.momentum_kg_m_per_s[d] =
        inventory.mass_kg * state.velocity_m_per_s[d];
  inventory.thermochemical_enthalpy_j = inventory.mass_kg * enthalpy;
  inventory.kinetic_energy_j =
      .5 * inventory.mass_kg *
      dot(state.velocity_m_per_s, state.velocity_m_per_s);
}
ParcelEventsReport pass(const ParcelEventsInput &in, ParcelPass which) {
  ParcelEventsReport result;
  result.revision = in.revision;
  SprayParcelState state = in.accepted_parcel;
  ParcelAuxiliaryState auxiliary = in.accepted_auxiliary;
  double time = in.elapsed_offset_s,
         step = std::min(in.initial_substep_s, in.duration_s);
  const double target_time = in.elapsed_offset_s + in.duration_s;
  auto reject = [&](ParcelEventsStatus status, portable::Revision revision) {
    auto out = fail(status, revision);
    out.failure_time_s = time;
    out.failure_segment = result.segment_count;
    out.failure_pass = which;
    return out;
  };
  std::array<ParcelEvent, 8> located_events{};
  std::size_t located_count = 0;
  TabBreakupReport located_tab{};
  while (time < target_time) {
    if (++result.attempts > in.maximum_attempts)
      return reject(ParcelEventsStatus::capacity_exceeded, in.revision);
    step = std::min(step, target_time - time);
    auto query = [&](const SprayParcelState &s, double t, double dt) {
      ++result.provider_queries;
      return in.interval->advance(s, t, dt, which, in.revision);
    };
    const auto whole = query(state, time, step);
    auto first = query(state, time, step / 2);
    if (!valid_interval(whole, step, in.revision) ||
        !valid_interval(first, step / 2, in.revision))
      return reject(ParcelEventsStatus::provider_failure, in.revision);
    if (!consistent(state, whole, in) || !consistent(state, first, in))
      return reject(ParcelEventsStatus::budget_failure, in.revision);
    ParcelIntervalReport refined = first;
    if (!first.complete_evaporation) {
      refined = query(first.parcel, time + step / 2, step / 2);
      if (!valid_interval(refined, step / 2, in.revision))
        return reject(ParcelEventsStatus::provider_failure, in.revision);
      if (!consistent(first.parcel, refined, in))
        return reject(ParcelEventsStatus::budget_failure, in.revision);
      refined.elapsed_duration_s += step / 2;
      refined.initial_liquid_absolute_enthalpy_j_per_kg =
          first.initial_liquid_absolute_enthalpy_j_per_kg;
      add(refined.exchange, first.exchange);
    }
    if (!accurate(whole, refined, in)) {
      ++result.rejected_substeps;
      if (step / 2 < in.minimum_substep_s)
        return reject(ParcelEventsStatus::minimum_step, in.revision);
      step /= 2;
      continue;
    }
    const double endpoint = time + refined.elapsed_duration_s;
    TabBreakupReport tab;
    if (in.tab_evolution) {
      const auto queried = in.tab_evolution->advance(state, auxiliary, time,
                                                     refined.elapsed_duration_s,
                                                     which, in.revision);
      ++result.provider_queries;
      if (queried.revision != in.revision)
        return reject(ParcelEventsStatus::stale_revision, in.revision);
      tab = queried.tab;
      if (!tab.succeeded())
        return reject(ParcelEventsStatus::provider_failure, in.revision);
      if (tab.breakup_requested &&
          (!std::isfinite(tab.event_time_s) || tab.event_time_s < 0 ||
           tab.event_time_s > refined.elapsed_duration_s))
        return reject(ParcelEventsStatus::event_failure, in.revision);
      // TAB coefficients are piecewise constant over an accepted substep.
      // Compare against two half coefficient queries when no threshold event
      // has interrupted the interval, then refine the same state interval.
      if (!tab.breakup_requested && !first.complete_evaporation) {
        auto first_tab = in.tab_evolution->advance(
            state, auxiliary, time, step / 2, which, in.revision);
        ++result.provider_queries;
        if (first_tab.revision != in.revision || !first_tab.tab.succeeded())
          return reject(ParcelEventsStatus::provider_failure, in.revision);
        auto half_aux = auxiliary;
        half_aux.tab_deformation = first_tab.tab.candidate.deformation;
        half_aux.tab_deformation_rate_per_s =
            first_tab.tab.candidate.deformation_rate_per_s;
        auto second_tab = in.tab_evolution->advance(
            first.parcel, half_aux, time + step / 2,
            refined.elapsed_duration_s - step / 2, which, in.revision);
        ++result.provider_queries;
        if (second_tab.revision != in.revision || !second_tab.tab.succeeded())
          return reject(ParcelEventsStatus::provider_failure, in.revision);
        if (first_tab.tab.breakup_requested ||
            second_tab.tab.breakup_requested ||
            !near(tab.candidate.deformation,
                  second_tab.tab.candidate.deformation, 1e-6,
                  in.relative_tolerance) ||
            !near(tab.candidate.deformation_rate_per_s,
                  second_tab.tab.candidate.deformation_rate_per_s, 1e-4,
                  in.relative_tolerance)) {
          ++result.rejected_substeps;
          if (step / 2 < in.minimum_substep_s)
            return reject(ParcelEventsStatus::minimum_step, in.revision);
          step /= 2;
          continue;
        }
      }
    }
    auto geometry = in.geometry->query(state, refined.parcel, time, endpoint,
                                       which, in.revision);
    if (!geometry.available)
      return reject(ParcelEventsStatus::provider_failure, in.revision);
    if (geometry.revision != in.revision)
      return reject(ParcelEventsStatus::stale_revision, in.revision);
    if (geometry.count > geometry.events.size())
      return reject(ParcelEventsStatus::capacity_exceeded, in.revision);
    std::array<ParcelEvent, 10> events{};
    std::size_t count = geometry.count;
    for (std::size_t i = 0; i < count; ++i) {
      events[i] = geometry.events[i];
      if (!std::isfinite(events[i].elapsed_time_s) ||
          events[i].elapsed_time_s < time ||
          events[i].elapsed_time_s > endpoint + in.event_time_tolerance_s ||
          static_cast<unsigned>(events[i].kind) > 4)
        return reject(ParcelEventsStatus::event_failure, in.revision);
    }
    // Preserve a provider-located event through reintegration. Floating point
    // endpoint coordinates need not lie bitwise on the surface, so requiring
    // the geometry to rediscover that same hit can miss a collision entirely.
    for (std::size_t i = 0; i < located_count; ++i) {
      if (std::abs(located_events[i].elapsed_time_s - endpoint) >
          in.event_time_tolerance_s)
        continue;
      bool duplicate = false;
      for (std::size_t j = 0; j < count; ++j)
        duplicate |= events[j].kind == located_events[i].kind &&
                     events[j].identity == located_events[i].identity;
      if (!duplicate) {
        if (count >= events.size() - 2)
          return reject(ParcelEventsStatus::capacity_exceeded, in.revision);
        events[count++] = located_events[i];
      }
    }
    if (in.tab_evolution && tab.breakup_requested) {
      const double event_time = time + tab.event_time_s;
      bool duplicate = false;
      for (std::size_t i = 0; i < count; ++i)
        duplicate |= events[i].kind == ParcelEventKind::breakup;
      if (!duplicate)
        events[count++] = {ParcelEventKind::breakup,
                           event_time,
                           auxiliary.breakup_ordinal,
                           0,
                           {},
                           1,
                           false};
    }
    if (refined.complete_evaporation)
      events[count++] = {
          ParcelEventKind::complete_evaporation, endpoint, 0, 0, {}, 1, false};
    std::sort(events.begin(), events.begin() + count,
              [](const auto &a, const auto &b) {
                if (a.elapsed_time_s != b.elapsed_time_s)
                  return a.elapsed_time_s < b.elapsed_time_s;
                if (a.kind != b.kind)
                  return a.kind < b.kind;
                return a.identity < b.identity;
              });
    if (count &&
        events[0].elapsed_time_s < endpoint - in.event_time_tolerance_s) {
      // Reintegrate to the provider's located event; never interpolate exchange
      // or continue using pre-collision force samples for the remainder.
      located_count = 0;
      for (std::size_t i = 0; i < count; ++i)
        if (std::abs(events[i].elapsed_time_s - events[0].elapsed_time_s) <=
            in.event_time_tolerance_s) {
          if (located_count == located_events.size())
            return reject(ParcelEventsStatus::capacity_exceeded, in.revision);
          located_events[located_count++] = events[i];
          if (events[i].kind == ParcelEventKind::breakup &&
              tab.breakup_requested)
            located_tab = tab;
        }
      step = events[0].elapsed_time_s - time;
      if (step <= 0) {
        const auto kind = events[0].kind;
        if (kind == ParcelEventKind::internal_cell_crossing ||
            kind == ParcelEventKind::wall_collision ||
            kind == ParcelEventKind::physical_outlet) {
          if (result.event_count + located_count > in.maximum_events)
            return reject(ParcelEventsStatus::capacity_exceeded, in.revision);
          located_events[0].applied = true;
          for (std::size_t i = 0; i < located_count; ++i)
            result.events[result.event_count++] = located_events[i];
          if (kind == ParcelEventKind::physical_outlet) {
            outlet_inventory(state,
                             whole.initial_liquid_absolute_enthalpy_j_per_kg,
                             result.outlet_inventory);
            result.physical_outlet = true;
            result.parent_removed = true;
            break;
          }
          if (kind == ParcelEventKind::wall_collision) {
            if (!apply_wall(state, events[0], result.wall_exchange))
              return reject(ParcelEventsStatus::event_failure, in.revision);
          } else {
            if (state.owner_global_cell == events[0].next_global_cell)
              return reject(ParcelEventsStatus::event_failure, in.revision);
            state.owner_global_cell = events[0].next_global_cell;
          }
          // Topology and wall impulses have no duration and no interphase
          // exchange. Retry from the updated state; the event cap bounds an
          // inconsistent provider that repeats a zero-time event forever.
          located_count = 0;
          step = std::min(in.initial_substep_s, target_time - time);
          continue;
        }
        // An accepted y already at its trigger (e.g. after a higher-priority
        // wall event) must retire at this time, not lose the remaining step.
        if (events[0].kind != ParcelEventKind::breakup || !in.breakup ||
            !tab.breakup_requested)
          return reject(ParcelEventsStatus::event_failure, in.revision);
        auxiliary.tab_deformation = tab.candidate.deformation;
        auxiliary.tab_deformation_rate_per_s =
            tab.candidate.deformation_rate_per_s;
        const auto split = in.breakup->split(state, events[0], auxiliary, &tab,
                                             which, in.revision);
        if (!valid_children(split, state, in.revision))
          return reject(ParcelEventsStatus::provider_failure, in.revision);
        if (result.event_count == in.maximum_events)
          return reject(ParcelEventsStatus::capacity_exceeded, in.revision);
        events[0].applied = true;
        result.events[result.event_count++] = events[0];
        result.children = split.candidate;
        result.breakup_budget = split.conservation;
        result.breakup_requested = true;
        result.parent_removed = true;
        break;
      }
      continue;
    }
    if (refined.elapsed_duration_s <= 0)
      return reject(ParcelEventsStatus::event_failure, in.revision);
    if (result.segment_count == in.maximum_segments)
      return reject(ParcelEventsStatus::capacity_exceeded, in.revision);
    auto &segment = result.segments[result.segment_count];
    segment.ordinal = static_cast<std::uint32_t>(result.segment_count++);
    segment.global_cell = state.owner_global_cell;
    segment.begin_time_s = time;
    segment.end_time_s = endpoint;
    segment.begin = state;
    segment.end = refined.parcel;
    segment.exchange = refined.exchange;
    add(result.exchange, refined.exchange);
    if (in.tab_evolution) {
      auxiliary.tab_deformation = tab.candidate.deformation;
      auxiliary.tab_deformation_rate_per_s =
          tab.candidate.deformation_rate_per_s;
      result.tab_evolved = true;
    }
    state = refined.parcel;
    time = endpoint;
    if (count) {
      located_count = 0;
      // All coincident candidates remain visible, but only the fixed-priority
      // event applies inventory. In particular evaporation plus outlet is not
      // counted as both vapor transfer and exported liquid inventory.
      std::sort(events.begin(), events.begin() + count,
                [](const auto &a, const auto &b) {
                  if (a.kind != b.kind)
                    return a.kind < b.kind;
                  return a.identity < b.identity;
                });
      if (result.event_count + count > in.maximum_events)
        return reject(ParcelEventsStatus::capacity_exceeded, in.revision);
      events[0].applied = true;
      for (std::size_t i = 0; i < count; ++i)
        result.events[result.event_count++] = events[i];
      const auto &event = events[0];
      if (event.kind == ParcelEventKind::complete_evaporation) {
        if (!refined.complete_evaporation)
          return reject(ParcelEventsStatus::event_failure, in.revision);
        result.complete_evaporation = true;
        result.parent_removed = true;
        break;
      }
      if (event.kind == ParcelEventKind::physical_outlet) {
        result.physical_outlet = true;
        result.parent_removed = true;
        outlet_inventory(state, refined.liquid_absolute_enthalpy_j_per_kg,
                         result.outlet_inventory);
        break;
      }
      if (event.kind == ParcelEventKind::wall_collision) {
        if (!apply_wall(state, event, result.wall_exchange))
          return reject(ParcelEventsStatus::event_failure, in.revision);
      } else if (event.kind == ParcelEventKind::breakup) {
        if (!in.breakup)
          return reject(ParcelEventsStatus::provider_failure, in.revision);
        const TabBreakupReport *trigger = nullptr;
        if (in.tab_evolution) {
          trigger =
              tab.breakup_requested
                  ? &tab
                  : (located_tab.breakup_requested ? &located_tab : nullptr);
          if (trigger) {
            auxiliary.tab_deformation = trigger->candidate.deformation;
            auxiliary.tab_deformation_rate_per_s =
                trigger->candidate.deformation_rate_per_s;
          }
        }
        const auto split = in.breakup->split(state, event, auxiliary, trigger,
                                             which, in.revision);
        if (!valid_children(split, state, in.revision))
          return reject(ParcelEventsStatus::provider_failure, in.revision);
        result.children = split.candidate;
        result.breakup_budget = split.conservation;
        result.breakup_requested = true;
        result.parent_removed = true;
        break;
      } else
        state.owner_global_cell = event.next_global_cell;
    }
    step = std::min(in.initial_substep_s, target_time - time);
  }
  result.status = ParcelEventsStatus::success;
  result.available = true;
  result.parcel = state;
  result.auxiliary = auxiliary;
  result.advanced_duration_s = time - in.elapsed_offset_s;
  if (result.breakup_requested)
    result.children_remaining_duration_s = target_time - time;
  result.exchange.available = true;
  return result;
}
} // namespace
ParcelEventsReport
integrate_parcel_events(const ParcelEventsInput &in) noexcept {
  const double positive[]{in.initial_substep_s,
                          in.minimum_substep_s,
                          in.position_absolute_tolerance_m,
                          in.velocity_absolute_tolerance_m_per_s,
                          in.mass_absolute_tolerance_kg,
                          in.temperature_absolute_tolerance_k,
                          in.momentum_absolute_tolerance_kg_m_per_s,
                          in.energy_absolute_tolerance_j,
                          in.event_time_tolerance_s};
  if (!in.interval || !in.geometry ||
      validate_parcel_state(in.accepted_parcel) != ParcelStateStatus::success ||
      !std::isfinite(in.accepted_auxiliary.tab_deformation) ||
      !std::isfinite(in.accepted_auxiliary.tab_deformation_rate_per_s) ||
      !std::isfinite(in.elapsed_offset_s) || in.elapsed_offset_s < 0 ||
      !std::isfinite(in.elapsed_offset_s + in.duration_s) ||
      (in.duration_s > 0 &&
       in.elapsed_offset_s + in.duration_s <= in.elapsed_offset_s) ||
      !std::isfinite(in.duration_s) || in.duration_s < 0 ||
      !std::isfinite(in.relative_tolerance) || in.relative_tolerance < 0 ||
      in.revision.algorithm_version != 1 || in.maximum_segments == 0 ||
      in.maximum_segments > kParcelEventSegmentCapacity ||
      in.maximum_events > kParcelEventCapacity || in.maximum_attempts == 0)
    return fail(ParcelEventsStatus::invalid_input, in.revision);
  for (double v : positive)
    if (!std::isfinite(v) || v <= 0)
      return fail(ParcelEventsStatus::invalid_input, in.revision);
  const auto predictor = pass(in, ParcelPass::predictor);
  if (!predictor.available)
    return predictor;
  auto corrector = pass(in, ParcelPass::corrector);
  if (!corrector.available)
    return corrector;
  corrector.predictor_passes = 1;
  corrector.corrector_passes = 1;
  corrector.attempts += predictor.attempts;
  corrector.rejected_substeps += predictor.rejected_substeps;
  corrector.provider_queries += predictor.provider_queries;
  return corrector;
}

ParcelIntervalReport FixedAsParcelIntervalProvider::advance(
    const SprayParcelState &state, double elapsed, double dt,
    ParcelPass pass_kind, portable::Revision revision) const noexcept {
  ParcelIntervalReport out;
  const auto before = provider_.sample(state, elapsed, pass_kind, revision);
  if (!before.available || before.revision != revision || !before.liquid ||
      !std::isfinite(before.liquid_absolute_enthalpy_j_per_kg))
    return out;
  FixedExchangeInput request;
  request.committed_parcel = state;
  request.liquid_properties = before.liquid;
  request.predictor_environment = before.gas;
  request.corrector_environment = before.gas;
  request.duration_s = dt;
  auto predicted = integrate_fixed_exchange(request);
  if (!predicted.succeeded())
    return out;
  auto advected = [&](SprayParcelState candidate, double duration) {
    for (std::size_t d = 0; d < 3; ++d)
      candidate.position_m[d] = static_cast<double>(
          static_cast<long double>(state.position_m[d]) +
          .5L *
              (static_cast<long double>(state.velocity_m_per_s[d]) +
               candidate.velocity_m_per_s[d]) *
              duration);
    return candidate;
  };
  if (!predicted.has_terminal_event) {
    const auto predicted_state =
        advected(predicted.candidate_parcel, predicted.advanced_duration_s);
    const auto after = provider_.sample(predicted_state,
                                        elapsed + predicted.advanced_duration_s,
                                        pass_kind, revision);
    if (!after.available || after.revision != revision ||
        after.liquid != before.liquid)
      return out;
    request.corrector_environment = after.gas;
  }
  const auto result = integrate_fixed_exchange(request);
  if (!result.succeeded())
    return out;
  const auto final_state =
      advected(result.candidate_parcel, result.advanced_duration_s);
  double h = before.liquid_absolute_enthalpy_j_per_kg;
  if (!result.has_terminal_event) {
    const auto final = provider_.sample(
        final_state, elapsed + result.advanced_duration_s, pass_kind, revision);
    if (!final.available || final.revision != revision ||
        final.liquid != before.liquid ||
        !std::isfinite(final.liquid_absolute_enthalpy_j_per_kg))
      return out;
    h = final.liquid_absolute_enthalpy_j_per_kg;
  }
  out.available = true;
  out.revision = revision;
  out.parcel = final_state;
  out.initial_liquid_absolute_enthalpy_j_per_kg =
      before.liquid_absolute_enthalpy_j_per_kg;
  out.elapsed_duration_s = result.advanced_duration_s;
  out.complete_evaporation = result.has_terminal_event;
  out.liquid_absolute_enthalpy_j_per_kg = h;
  out.exchange = result.exchange;
  out.exchange.parcel_kinetic_energy_delta_j = kinetic_delta(state, out.parcel);
  out.exchange.parcel_thermochemical_enthalpy_delta_j =
      state.multiplicity *
      (out.parcel.droplet_mass_kg * h -
       state.droplet_mass_kg * before.liquid_absolute_enthalpy_j_per_kg);
  out.exchange.thermal_exchange_state_residual_j =
      out.exchange.parcel_thermochemical_enthalpy_delta_j +
      out.exchange.thermal_exchange_to_gas_j;
  return out;
}
ParcelTabIntervalReport
FixedTabEvolutionProvider::advance(const SprayParcelState &state,
                                   const ParcelAuxiliaryState &auxiliary,
                                   double elapsed, double dt, ParcelPass which,
                                   portable::Revision revision) const noexcept {
  ParcelTabIntervalReport out;
  out.revision = revision;
  const auto environment = provider_.sample(state, elapsed, which, revision);
  if (!environment.available || environment.revision != revision ||
      !environment.liquid ||
      !std::isfinite(environment.far_gas_density_kg_per_m3) ||
      environment.far_gas_density_kg_per_m3 <= 0)
    return out;
  const auto liquid = environment.liquid->evaluate(
      {state.liquid_material_fingerprint, state.temperature_k});
  if (!liquid.succeeded() ||
      validate_liquid_property_report(
          {state.liquid_material_fingerprint, state.temperature_k}, liquid) !=
          LiquidPropertyStatus::success)
    return out;
  Vector3 slip{};
  for (std::size_t d = 0; d < 3; ++d)
    slip[d] =
        environment.gas.gas_velocity_m_per_s[d] - state.velocity_m_per_s[d];
  TabBreakupInput input;
  input.initial_deformation = auxiliary.tab_deformation;
  input.initial_deformation_rate_per_s = auxiliary.tab_deformation_rate_per_s;
  input.relative_speed_m_per_s = std::sqrt(dot(slip, slip));
  input.gas_density_kg_per_m3 = environment.far_gas_density_kg_per_m3;
  input.liquid_density_kg_per_m3 = liquid.properties.density_kg_per_m3;
  input.liquid_viscosity_pa_s = liquid.properties.viscosity_pa_s;
  input.surface_tension_n_per_m = liquid.properties.surface_tension_n_per_m;
  input.droplet_radius_m = state.droplet_diameter_m / 2;
  input.breakup_threshold = 1;
  input.duration_s = dt;
  input.coefficients = {2.0 / 3.0, .5, 8, 5};
  out.tab = evaluate_tab_breakup(input);
  return out;
}
BreakupChildReport FixedTabEventBreakupProvider::split(
    const SprayParcelState &state, const ParcelEvent &event,
    const ParcelAuxiliaryState &auxiliary, const TabBreakupReport *trigger,
    ParcelPass which, portable::Revision revision) const noexcept {
  if (!trigger || !trigger->succeeded() || !trigger->breakup_requested ||
      event.kind != ParcelEventKind::breakup)
    return {};
  const auto environment =
      provider_.sample(state, event.elapsed_time_s, which, revision);
  if (!environment.available || environment.revision != revision ||
      !environment.liquid)
    return {};
  const auto liquid = environment.liquid->evaluate(
      {state.liquid_material_fingerprint, state.temperature_k});
  if (!liquid.succeeded() ||
      validate_liquid_property_report(
          {state.liquid_material_fingerprint, state.temperature_k}, liquid) !=
          LiquidPropertyStatus::success)
    return {};
  // Re-evaluate the zero-duration trigger at the located state. This provides
  // the current radius/material stiffness to the representative-size closure,
  // without advancing or modifying the accepted deformation history.
  const auto local = FixedTabEvolutionProvider(provider_).advance(
      state, auxiliary, event.elapsed_time_s, 0, which, revision);
  if (!local.tab.succeeded() || !local.tab.breakup_requested)
    return {};
  TabRepresentativeSplitInput input;
  input.parent = state;
  input.tab_trigger = local.tab;
  input.accepted_step = revision.accepted_step;
  input.breakup_ordinal = auxiliary.breakup_ordinal;
  input.child_parcel_count = child_parcels_;
  input.liquid_density_kg_per_m3 = liquid.properties.density_kg_per_m3;
  input.surface_tension_n_per_m = liquid.properties.surface_tension_n_per_m;
  input.liquid_absolute_thermochemical_enthalpy_j_per_kg =
      environment.liquid_absolute_enthalpy_j_per_kg;
  for (std::size_t d = 0; d < 3; ++d)
    input.breakup_axis[d] =
        environment.gas.gas_velocity_m_per_s[d] - state.velocity_m_per_s[d];
  // A zero-slip oscillating drop has no aerodynamic axis. Its fixed x-axis is
  // an explicit deterministic orientation convention, not a new random model.
  if (dot(input.breakup_axis, input.breakup_axis) == 0)
    input.breakup_axis = {1, 0, 0};
  return generate_tab_representative_children(input).split;
}
} // namespace hundun::v04::spray::detail

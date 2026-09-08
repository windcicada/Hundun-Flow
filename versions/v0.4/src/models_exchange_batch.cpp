// SPDX-License-Identifier: Apache-2.0
#include "models_exchange_batch_detail.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>

namespace hundun::v04::portable {
namespace {
bool finite(const spray::detail::ParcelExchangeBudget &d) noexcept {
  return std::isfinite(d.mass_delta_kg) &&
         std::isfinite(d.momentum_delta_kg_m_per_s[0]) &&
         std::isfinite(d.momentum_delta_kg_m_per_s[1]) &&
         std::isfinite(d.momentum_delta_kg_m_per_s[2]) &&
         std::isfinite(d.thermochemical_enthalpy_delta_j) &&
         std::isfinite(d.kinetic_energy_delta_j) &&
         std::isfinite(d.thermal_exchange_to_gas_j);
}
void add(spray::detail::ParcelExchangeBudget &sum,
         const spray::detail::ParcelExchangeBudget &d, double w) noexcept {
  sum.mass_delta_kg += w * d.mass_delta_kg;
  for (int k = 0; k < 3; ++k)
    sum.momentum_delta_kg_m_per_s[k] += w * d.momentum_delta_kg_m_per_s[k];
  sum.thermochemical_enthalpy_delta_j += w * d.thermochemical_enthalpy_delta_j;
  sum.kinetic_energy_delta_j += w * d.kinetic_energy_delta_j;
  sum.thermal_exchange_to_gas_j += w * d.thermal_exchange_to_gas_j;
}
void add(ExchangeDelta &sum,
         const spray::detail::ParcelExchangeBudget &d) noexcept {
  sum.mass_kg += d.mass_delta_kg;
  for (int k = 0; k < 3; ++k)
    sum.momentum_kg_m_per_s[k] += d.momentum_delta_kg_m_per_s[k];
  sum.thermochemical_enthalpy_j += d.thermochemical_enthalpy_delta_j;
  sum.kinetic_energy_j += d.kinetic_energy_delta_j;
}
bool finite(const ExchangeDelta &d) noexcept {
  return std::isfinite(d.mass_kg) &&
         std::isfinite(d.thermochemical_enthalpy_j) &&
         std::isfinite(d.kinetic_energy_j) &&
         std::isfinite(d.momentum_kg_m_per_s[0]) &&
         std::isfinite(d.momentum_kg_m_per_s[1]) &&
         std::isfinite(d.momentum_kg_m_per_s[2]);
}
} // namespace
ExchangeWorkspace::ExchangeWorkspace(std::size_t cells, std::size_t segments,
                                     std::size_t species)
    : cells_(cells), order_(segments), species_(species) {
  if (!species || cells > std::numeric_limits<std::size_t>::max() / species)
    throw std::invalid_argument("invalid exchange species capacity");
  species_mass_.resize(cells * species);
}
bool ExchangeWorkspace::current(const ExchangeBatchReport &r) const noexcept {
  return r.available && r.generation == generation_ && r.cells == cells_.data();
}
ExchangeBatchReport
ExchangeWorkspace::evaluate(Revision revision, const ExchangeCell *cells,
                            std::size_t nc, const ExchangeSegment *segments,
                            std::size_t ns, double atol, double rtol) noexcept {
  ++generation_;
  auto fail = [&](Status s, std::size_t i) {
    ExchangeBatchReport r;
    r.status = s;
    r.failure_index = i;
    r.revision = revision;
    r.generation = generation_;
    return r;
  };
  if (nc > cells_.size() || ns > order_.size())
    return fail(Status::capacity_exceeded, 0);
  if ((nc && !cells) || (ns && !segments) || revision.algorithm_version != 1 ||
      !std::isfinite(atol) || atol < 0 || !std::isfinite(rtol) || rtol < 0 ||
      rtol >= 1)
    return fail(Status::invalid_input, 0);
  for (std::size_t i = 0; i < nc; ++i) {
    if (!(cells[i].volume_m3 > 0) || !std::isfinite(cells[i].volume_m3) ||
        (i && cells[i - 1].global_cell >= cells[i].global_cell))
      return fail(Status::invalid_input, i);
  }
  for (std::size_t i = 0; i < ns; ++i) {
    const auto &s = segments[i];
    if (s.revision != revision)
      return fail(Status::stale_revision, i);
    if ((s.parcel_id.high == 0 && s.parcel_id.low == 0) || !finite(s.delta) ||
        !std::isfinite(s.deposition_weight) || s.deposition_weight <= 0 ||
        s.deposition_weight > 1 || s.channel > ExchangeChannel::external)
      return fail(Status::invalid_input, i);
    if (s.channel != ExchangeChannel::interphase && s.deposition_weight != 1)
      return fail(Status::invalid_input, i);
    if (s.channel == ExchangeChannel::interphase) {
      if (s.vapor_species_index >= species_)
        return fail(Status::identity_mismatch, i);
      bool found = false;
      for (std::size_t j = 0; j < nc; ++j)
        found |= cells[j].global_cell == s.global_cell;
      if (!found)
        return fail(Status::invalid_input, i);
      const long double residual =
          static_cast<long double>(s.delta.thermochemical_enthalpy_delta_j) +
          s.delta.thermal_exchange_to_gas_j;
      const long double limit =
          atol + static_cast<long double>(rtol) *
                     (std::abs(s.delta.thermochemical_enthalpy_delta_j) +
                      std::abs(s.delta.thermal_exchange_to_gas_j));
      if (std::abs(residual) > limit)
        return fail(Status::conservation_failure, i);
    }
    order_[i] = i;
  }
  // Validate deposition once per physical segment, not once per field.
  for (std::size_t i = 0; i < ns; ++i) {
    const auto &a = segments[i];
    long double weight = 0;
    for (std::size_t j = 0; j < ns; ++j) {
      const auto &b = segments[j];
      if (a.parcel_id != b.parcel_id || a.segment_ordinal != b.segment_ordinal)
        continue;
      if (a.channel != b.channel ||
          a.vapor_species_index != b.vapor_species_index)
        return fail(Status::invalid_input, j);
      if (j != i && (a.global_cell == b.global_cell ||
                     a.channel != ExchangeChannel::interphase))
        return fail(Status::invalid_input, j);
      // All deposition rows refer to the same unsplit physical exchange.
      if (a.delta.mass_delta_kg != b.delta.mass_delta_kg ||
          a.delta.momentum_delta_kg_m_per_s !=
              b.delta.momentum_delta_kg_m_per_s ||
          a.delta.thermochemical_enthalpy_delta_j !=
              b.delta.thermochemical_enthalpy_delta_j ||
          a.delta.kinetic_energy_delta_j != b.delta.kinetic_energy_delta_j ||
          a.delta.thermal_exchange_to_gas_j !=
              b.delta.thermal_exchange_to_gas_j)
        return fail(Status::invalid_input, j);
      weight += b.deposition_weight;
    }
    if (std::abs(weight - 1.0L) > 32 * std::numeric_limits<double>::epsilon())
      return fail(Status::conservation_failure, i);
  }
  std::sort(order_.begin(), order_.begin() + ns,
            [&](std::size_t i, std::size_t j) {
              const auto &a = segments[i];
              const auto &b = segments[j];
              return std::tie(a.global_cell, a.parcel_id.high, a.parcel_id.low,
                              a.segment_ordinal) <
                     std::tie(b.global_cell, b.parcel_id.high, b.parcel_id.low,
                              b.segment_ordinal);
            });
  ExchangeBatchReport result;
  std::fill(species_mass_.begin(), species_mass_.end(), 0.0);
  for (std::size_t i = 0; i < nc; ++i) {
    spray::detail::GasCellExchangeInput in;
    in.gas_mass_kg = cells[i].gas_mass_kg;
    in.gas_momentum_kg_m_per_s = cells[i].gas_momentum_kg_m_per_s;
    in.thermal_absolute_tolerance_j = atol;
    in.thermal_relative_tolerance = rtol;
    for (std::size_t j = 0; j < ns; ++j) {
      const auto &s = segments[order_[j]];
      if (s.channel == ExchangeChannel::interphase &&
          s.global_cell == cells[i].global_cell) {
        add(in.parcel, s.delta, s.deposition_weight);
        species_mass_[i * species_ + s.vapor_species_index] -=
            s.deposition_weight * s.delta.mass_delta_kg;
      }
    }
    auto gas = spray::detail::make_gas_cell_exchange_candidate(in);
    if (!gas.succeeded())
      return fail(Status::conservation_failure, i);
    for (std::size_t s = 0; s < species_; ++s)
      if (!std::isfinite(species_mass_[i * species_ + s]))
        return fail(Status::conservation_failure, i);
    cells_[i] = {cells[i].global_cell, cells[i].volume_m3, gas,
                 species_mass_.data() + i * species_, species_};
  }
  for (std::size_t j = 0; j < ns; ++j) {
    const auto &s = segments[order_[j]];
    switch (s.channel) {
    case ExchangeChannel::wall:
      add(result.wall, s.delta);
      break;
    case ExchangeChannel::outlet:
      add(result.outlet, s.delta);
      break;
    case ExchangeChannel::external:
      add(result.external, s.delta);
      break;
    default:
      break;
    }
  }
  if (!finite(result.wall) || !finite(result.outlet) ||
      !finite(result.external))
    return fail(Status::conservation_failure, ns);
  result.status = Status::success;
  result.available = true;
  result.generation = generation_;
  result.revision = revision;
  result.cells = cells_.data();
  result.cell_count = nc;
  return result;
}
} // namespace hundun::v04::portable

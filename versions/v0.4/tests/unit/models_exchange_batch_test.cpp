// SPDX-License-Identifier: Apache-2.0
#include "models_exchange_batch_detail.hpp"
#include <cmath>
#include <iostream>
using namespace hundun::v04::portable;
int main() {
  // Two 1 kg droplets each slow from 2 to 1 m/s; one 2 kg gas cell
  // starts at rest. Independent endpoint K budget: liquid -3 J, gas +1 J,
  // gas thermal +2 J. Per-parcel gas corrections would give wrong +0.5 J.
  ExchangeCell cell{7, 2.0, 2.0, {}};
  ExchangeSegment segments[2]{};
  for (int i = 0; i < 2; ++i) {
    segments[i].global_cell = 7;
    segments[i].parcel_id = {1, static_cast<std::uint64_t>(2 - i)};
    segments[i].delta.momentum_delta_kg_m_per_s = {-1, 0, 0};
    segments[i].delta.kinetic_energy_delta_j = -1.5;
  }
  ExchangeWorkspace workspace(1, 2);
  auto result = workspace.evaluate({}, &cell, 1, segments, 2, 0, 0);
  if (!result.available || result.cell_count != 1 ||
      result.cells[0].gas.gas_kinetic_energy_delta_j != 1.0 ||
      result.cells[0].gas.gas_thermochemical_enthalpy_delta_j != 2.0) {
    std::cerr << "aggregate before nonlinear kinetic correction\n";
    return 1;
  }
  auto previous = result;
  std::swap(segments[0], segments[1]);
  result = workspace.evaluate({}, &cell, 1, segments, 2, 0, 0);
  if (!result.available || !workspace.current(result) ||
      workspace.current(previous) ||
      result.cells[0].gas.gas_thermochemical_enthalpy_delta_j != 2.0)
    return 2;
  segments[1].parcel_id = segments[0].parcel_id;
  result = workspace.evaluate({}, &cell, 1, segments, 2, 0, 0);
  if (result.available || result.cells || result.cell_count ||
      result.status != Status::invalid_input || workspace.current(previous))
    return 3;
  segments[1].parcel_id = {2, 3};
  segments[1].revision.input_revision = 1;
  result = workspace.evaluate({}, &cell, 1, segments, 2, 0, 0);
  if (result.available || result.status != Status::stale_revision)
    return 4;
  segments[1].revision = {};
  result = workspace.evaluate({}, &cell, 1, segments, 3, 0, 0);
  if (result.available || result.status != Status::capacity_exceeded)
    return 5;
  // Separate outlet budget is not an interphase gas source.
  segments[1].channel = ExchangeChannel::outlet;
  result = workspace.evaluate({}, &cell, 1, segments, 2, 0, 0);
  if (!result.available ||
      result.cells[0].gas.gas_kinetic_energy_delta_j != 0.25 ||
      result.cells[0].gas.gas_thermochemical_enthalpy_delta_j != 1.25 ||
      result.outlet.kinetic_energy_j != -1.5)
    return 6;
  result = workspace.evaluate({}, &cell, 1, nullptr, 0, 0, 0);
  if (!result.available || result.cells[0].gas.gas_mass_delta_kg != 0 ||
      result.cells[0].gas.gas_thermochemical_enthalpy_delta_j != 0)
    return 7;
  // Two half-weight destinations split one represented exchange once.
  ExchangeCell cells[2]{{7, 1, 1, {}}, {8, 1, 1, {}}};
  segments[1] = segments[0];
  segments[1].global_cell = 8;
  segments[0].deposition_weight = segments[1].deposition_weight = 0.5;
  ExchangeWorkspace two_cells(2, 2);
  result = two_cells.evaluate({}, cells, 2, segments, 2, 0, 0);
  if (!result.available ||
      result.cells[0].gas.gas_kinetic_energy_delta_j != 0.125 ||
      result.cells[1].gas.gas_thermochemical_enthalpy_delta_j != 0.625)
    return 8;
  segments[1].deposition_weight = 0.4;
  result = two_cells.evaluate({}, cells, 2, segments, 2, 0, 0);
  if (result.available || result.status != Status::conservation_failure)
    return 9;
  // Vapor mapping follows explicit index, not liquid or gas species names.
  ExchangeWorkspace mapped(2, 2, 2);
  segments[0].vapor_species_index = segments[1].vapor_species_index = 1;
  segments[1].deposition_weight = 0.5;
  for (auto &s : segments)
    s.delta.mass_delta_kg = -0.25;
  result = mapped.evaluate({}, cells, 2, segments, 2, 0, 0);
  if (!result.available || result.cells[0].gas_species_mass_delta_kg[0] != 0 ||
      result.cells[0].gas_species_mass_delta_kg[1] != 0.125)
    return 10;
  return 0;
}

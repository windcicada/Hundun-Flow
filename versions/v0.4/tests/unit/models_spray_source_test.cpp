// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "models_spray_source_detail.hpp"

#include <cmath>
#include <iostream>
#include <limits>

using namespace hundun::v04::spray::detail;

int main() {
  bool passed = true;
  auto expect = [&](bool condition, const char* message) {
    if (!condition) std::cerr << message << '\n';
    passed &= condition;
  };
  GasCellExchangeInput input;
  input.gas_mass_kg = 1.0;
  input.parcel.momentum_delta_kg_m_per_s = {-1.0, 0.0, 0.0};
  // One kg liquid changes velocity 2 -> 1 m/s. Gas starts at rest.
  input.parcel.kinetic_energy_delta_j = -1.5;
  auto result = make_gas_cell_exchange_candidate(input);
  expect(result.succeeded() && result.model_id == "closed_isobaric_exchange_v1",
         "candidate has deterministic physical identity");
  expect(result.gas_momentum_delta_kg_m_per_s[0] == 1.0 &&
         result.gas_kinetic_energy_delta_j == 0.5 &&
         result.gas_thermochemical_enthalpy_delta_j == 1.0,
         "drag transfers one joule of resolved kinetic energy into thermal energy");
  expect(result.total_energy_residual_j == 0.0,
         "gas plus liquid H+K closes, not only opposite thermal deltas");

  // A common +4 m/s boost must leave the thermal source unchanged.
  input.gas_momentum_kg_m_per_s = {4.0, 0.0, 0.0};
  input.parcel.kinetic_energy_delta_j = -5.5;
  result = make_gas_cell_exchange_candidate(input);
  expect(result.succeeded() && result.gas_kinetic_energy_delta_j == 4.5 &&
         result.gas_thermochemical_enthalpy_delta_j == 1.0,
         "exchange thermal source is Galilean invariant");

  input = {};
  input.gas_mass_kg = 1.0;
  input.parcel.mass_delta_kg = -0.5;
  input.parcel.momentum_delta_kg_m_per_s = {-1.0, 0.0, 0.0};
  input.parcel.thermochemical_enthalpy_delta_j = -5.0;
  input.parcel.thermal_exchange_to_gas_j = 5.0;
  input.parcel.kinetic_energy_delta_j = -1.0;
  result = make_gas_cell_exchange_candidate(input);
  expect(result.succeeded() && result.gas_mass_delta_kg == 0.5 &&
         std::abs(result.gas_thermochemical_enthalpy_delta_j - 17.0 / 3.0) < 1e-14 &&
         std::abs(result.total_energy_residual_j) < 1e-14,
         "evaporated mass carries momentum and energy into finite gas inventory");

  input.parcel.thermal_exchange_to_gas_j = 4.0;
  result = make_gas_cell_exchange_candidate(input);
  expect(!result.succeeded() && !result.available &&
         result.status == GasCellExchangeStatus::inconsistent_thermal_budget,
         "physical flux/state residual cannot be silently absorbed into gas enthalpy");

  input.parcel.thermal_exchange_to_gas_j = 5.0 + 1e-8;
  input.thermal_absolute_tolerance_j = 1e-7;
  result = make_gas_cell_exchange_candidate(input);
  expect(result.succeeded() && result.thermal_budget_residual_j > 0.0 &&
         result.thermal_budget_residual_j < input.thermal_absolute_tolerance_j,
         "accepted integration residual is explicitly bounded and reported");

  input = {};
  input.gas_mass_kg = 1.0;
  input.gas_momentum_kg_m_per_s = {1.0e100, -1.0e100, 1.0};
  result = make_gas_cell_exchange_candidate(input);
  expect(result.succeeded() && result.gas_kinetic_energy_delta_j == 0.0 &&
         result.gas_thermochemical_enthalpy_delta_j == 0.0,
         "zero exchange is exact even at large background kinetic energy");

  input.parcel.mass_delta_kg = 1.0;
  result = make_gas_cell_exchange_candidate(input);
  expect(!result.succeeded() && !result.available &&
         result.gas_mass_delta_kg == 0.0,
         "exhausting gas mass rejects all partial outputs");
  input.parcel.mass_delta_kg = 0.0;
  input.parcel.thermochemical_enthalpy_delta_j =
      std::numeric_limits<double>::quiet_NaN();
  result = make_gas_cell_exchange_candidate(input);
  expect(!result.succeeded() && !result.available &&
         result.gas_thermochemical_enthalpy_delta_j == 0.0,
         "non-finite exchange produces only explicit failure");
  return passed ? 0 : 1;
}

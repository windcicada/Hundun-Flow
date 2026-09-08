// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "models_spray_source_detail.hpp"

#include <cmath>

namespace hundun::v04::spray::detail {
namespace {
bool finite(Vector3 value) noexcept {
  return std::isfinite(value[0]) && std::isfinite(value[1]) &&
         std::isfinite(value[2]);
}
double scalar(long double value) noexcept {
  return value == 0.0L ? 0.0 : static_cast<double>(value);
}
GasCellExchangeCandidate failure(GasCellExchangeStatus status) noexcept {
  GasCellExchangeCandidate result;
  result.status = status;
  return result;
}
}  // namespace

GasCellExchangeCandidate make_gas_cell_exchange_candidate(
    const GasCellExchangeInput& input) noexcept {
  const auto& parcel = input.parcel;
  if (!std::isfinite(input.gas_mass_kg) || input.gas_mass_kg <= 0.0 ||
      !finite(input.gas_momentum_kg_m_per_s) ||
      !std::isfinite(parcel.mass_delta_kg) ||
      !finite(parcel.momentum_delta_kg_m_per_s) ||
      !std::isfinite(parcel.thermochemical_enthalpy_delta_j) ||
      !std::isfinite(parcel.kinetic_energy_delta_j) ||
      !std::isfinite(parcel.thermal_exchange_to_gas_j) ||
      !std::isfinite(input.thermal_absolute_tolerance_j) ||
      input.thermal_absolute_tolerance_j < 0.0 ||
      !std::isfinite(input.thermal_relative_tolerance) ||
      input.thermal_relative_tolerance < 0.0 ||
      input.thermal_relative_tolerance >= 1.0)
    return failure(GasCellExchangeStatus::invalid_input);

  const long double mass = input.gas_mass_kg;
  const long double delta_mass = -static_cast<long double>(parcel.mass_delta_kg);
  const long double new_mass = mass + delta_mass;
  if (new_mass <= 0.0L)
    return failure(GasCellExchangeStatus::exhausted_gas_mass);

  GasCellExchangeCandidate result;
  const long double thermal_residual =
      static_cast<long double>(parcel.thermochemical_enthalpy_delta_j) +
      static_cast<long double>(parcel.thermal_exchange_to_gas_j);
  const long double thermal_limit = input.thermal_absolute_tolerance_j +
      static_cast<long double>(input.thermal_relative_tolerance) *
      (std::abs(static_cast<long double>(parcel.thermochemical_enthalpy_delta_j)) +
       std::abs(static_cast<long double>(parcel.thermal_exchange_to_gas_j)));
  if (std::abs(thermal_residual) > thermal_limit)
    return failure(GasCellExchangeStatus::inconsistent_thermal_budget);
  result.thermal_budget_residual_j = scalar(thermal_residual);
  result.gas_mass_delta_kg = scalar(delta_mass);
  result.gas_mass_candidate_kg = scalar(new_mass);
  long double delta_k = 0.0L;
  for (std::size_t component = 0; component < 3U; ++component) {
    const long double p = input.gas_momentum_kg_m_per_s[component];
    const long double delta_p = -static_cast<long double>(
        parcel.momentum_delta_kg_m_per_s[component]);
    const long double velocity = p / mass;
    const long double relative_impulse = delta_p - delta_mass * velocity;
    // Algebraically K(m+dm,p+dp)-K(m,p), evaluated without subtracting
    // two large background kinetic energies. Exact zero input stays zero.
    delta_k += velocity * delta_p -
        0.5L * delta_mass * velocity * velocity +
        0.5L * relative_impulse * relative_impulse / new_mass;
    result.gas_momentum_delta_kg_m_per_s[component] = scalar(delta_p);
    result.gas_momentum_candidate_kg_m_per_s[component] = scalar(p + delta_p);
  }
  result.gas_kinetic_energy_delta_j = scalar(delta_k);
  result.gas_thermochemical_enthalpy_delta_j = scalar(
      -static_cast<long double>(parcel.thermochemical_enthalpy_delta_j) -
      static_cast<long double>(parcel.kinetic_energy_delta_j) - delta_k);
  result.total_energy_residual_j = scalar(
      static_cast<long double>(result.gas_thermochemical_enthalpy_delta_j) +
      static_cast<long double>(result.gas_kinetic_energy_delta_j) +
      static_cast<long double>(parcel.thermochemical_enthalpy_delta_j) +
      static_cast<long double>(parcel.kinetic_energy_delta_j));
  if (!std::isfinite(result.gas_mass_candidate_kg) ||
      !finite(result.gas_momentum_candidate_kg_m_per_s) ||
      !std::isfinite(result.gas_kinetic_energy_delta_j) ||
      !std::isfinite(result.gas_thermochemical_enthalpy_delta_j) ||
      !std::isfinite(result.thermal_budget_residual_j) ||
      !std::isfinite(result.total_energy_residual_j))
    return failure(GasCellExchangeStatus::non_finite_output);
  result.status = GasCellExchangeStatus::success;
  result.available = true;
  return result;
}

}  // namespace hundun::v04::spray::detail

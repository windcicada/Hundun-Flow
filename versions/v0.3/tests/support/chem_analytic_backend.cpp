// SPDX-License-Identifier: Apache-2.0

#include "tests/support/chem_analytic_backend.hpp"

#include <cmath>
#include <stdexcept>

namespace hundun::test {
namespace {

constexpr double cp_j_per_kg_k = 1000.0;
constexpr double formation_h0_j_per_kg = -400.0;
constexpr double formation_h1_j_per_kg = 20.0;
constexpr double density_kg_per_m3 = 2.0;
constexpr double reaction_rate_per_s = 2.0;
constexpr double equilibrium_y0 = 0.25;

bool valid_point(const chemistry::ThermochemicalPoint &point) {
  return std::isfinite(point.p0_pa) && point.p0_pa > 0.0 &&
         std::isfinite(point.h_tc_j_per_kg) &&
         point.mass_fractions.size() == 2U &&
         std::isfinite(point.mass_fractions[0]) &&
         std::isfinite(point.mass_fractions[1]) &&
         point.mass_fractions[0] >= 0.0 &&
         point.mass_fractions[1] >= 0.0 &&
         std::abs(point.mass_fractions[0] + point.mass_fractions[1] - 1.0) <=
             1.0e-12;
}

} // namespace

AnalyticReactingBackend::AnalyticReactingBackend() {
  composition_.element_names = {"A"};
  composition_.species = {{"A2", 2.0, {2}}, {"A", 1.0, {1}}};
  composition_.fingerprint =
      chemistry::composition_identity_fingerprint(composition_);
  chemistry::validate_composition_identity(composition_);
}

AnalyticReactingBackend::~AnalyticReactingBackend() = default;

const chemistry::CompositionIdentity &
AnalyticReactingBackend::composition() const noexcept {
  return composition_;
}

chemistry::ThermodynamicProperties AnalyticReactingBackend::evaluate(
    const chemistry::ThermochemicalPoint &point) const {
  if (!valid_point(point)) {
    throw std::invalid_argument("invalid analytic thermochemical point");
  }
  const double formation_enthalpy =
      point.mass_fractions[0] * formation_h0_j_per_kg +
      point.mass_fractions[1] * formation_h1_j_per_kg;
  const double temperature =
      (point.h_tc_j_per_kg - formation_enthalpy) / cp_j_per_kg_k;
  if (!std::isfinite(temperature) || temperature <= 0.0) {
    throw std::invalid_argument("invalid analytic temperature inversion");
  }
  return {temperature, density_kg_per_m3, cp_j_per_kg_k,
          1.0 / (point.mass_fractions[0] / 2.0 + point.mass_fractions[1])};
}

chemistry::TransportProperties AnalyticReactingBackend::evaluate(
    const chemistry::ThermochemicalPoint &point,
    const chemistry::ThermodynamicProperties &thermodynamics) const {
  static_cast<void>(evaluate(point));
  if (!std::isfinite(thermodynamics.temperature_k) ||
      thermodynamics.temperature_k <= 0.0) {
    throw std::invalid_argument("invalid analytic transport temperature");
  }
  return {1.0e-5, 0.02, {1.0e-5, 2.0e-5}};
}

chemistry::ChemistryIntervalReport AnalyticReactingBackend::integrate(
    const chemistry::ChemistryIntervalRequest &request) {
  chemistry::ChemistryIntervalReport report;
  report.final_state = request.state;
  report.integrated_rho_y_delta_kg_per_m3 = {0.0, 0.0};
  if (!valid_point(request.state) || !std::isfinite(request.start_time_s) ||
      !std::isfinite(request.duration_s) || request.duration_s < 0.0) {
    report.status = chemistry::ChemistryStatus::invalid_input;
    return report;
  }
  const double initial_y0 = request.state.mass_fractions[0];
  const double final_y0 = equilibrium_y0 +
                          (initial_y0 - equilibrium_y0) *
                              std::exp(-reaction_rate_per_s * request.duration_s);
  const double delta_y0 = final_y0 - initial_y0;
  report.final_state.mass_fractions = {final_y0, 1.0 - final_y0};
  report.integrated_rho_y_delta_kg_per_m3 =
      {density_kg_per_m3 * delta_y0, -density_kg_per_m3 * delta_y0};
  report.status = chemistry::ChemistryStatus::success;
  report.completed_duration_s = request.duration_s;
  report.internal_step_count = request.duration_s == 0.0 ? 0U : 1U;
  return report;
}

std::unique_ptr<AnalyticReactingBackend>
make_analytic_reacting_backend_for_tests() {
  return std::unique_ptr<AnalyticReactingBackend>(
      new AnalyticReactingBackend());
}

} // namespace hundun::test

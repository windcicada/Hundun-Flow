// SPDX-License-Identifier: Apache-2.0
#include "models_chemistry_adapter_detail.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace hundun::v04::chemistry::detail {
namespace {
bool fractions(const double *y, std::size_t n) noexcept {
  if (!y || !n)
    return false;
  double sum = 0;
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(y[i]) || y[i] < 0 || y[i] > 1)
      return false;
    sum += y[i];
  }
  return std::abs(sum - 1) <= 1e-12 * n;
}
bool near(double a, double b) noexcept {
  return std::isfinite(a) && std::isfinite(b) &&
         std::abs(a - b) <= 1e-10 + 1e-9 * std::max(std::abs(a), std::abs(b));
}
portable::Status validate_query(const portable::GasQuery &q,
                                const portable::GasIdentity &id,
                                const portable::GasQueryOutput &out) noexcept {
  if (q.composition_fingerprint != id.composition_fingerprint)
    return portable::Status::identity_mismatch;
  if (out.capacity < id.species_names.size())
    return portable::Status::capacity_exceeded;
  if (!out.diffusivities_m2_per_s || !out.species_enthalpies_j_per_kg ||
      !out.net_mass_rates_kg_per_m3_s ||
      q.species_count != id.species_names.size() ||
      !fractions(q.mass_fractions, q.species_count) ||
      !std::isfinite(q.pressure_pa) || q.pressure_pa <= 0 ||
      q.revision.algorithm_version != 1)
    return portable::Status::invalid_input;
  return portable::Status::success;
}
} // namespace

AnalyticIsomerBackend::AnalyticIsomerBackend(double rate, bool reverse,
                                             double cp)
    : rate_(rate), cp_(cp) {
  if (!std::isfinite(rate) || rate < 0 || (cp != 1000 && cp != 1200))
    throw std::invalid_argument("invalid analytic specification");
  if (reverse) {
    reactant_ = 1;
    product_ = 0;
  }
  composition_.element_names = {"X"};
  composition_.species = {{reverse ? "B" : "A", 28, {1}},
                          {reverse ? "A" : "B", 28, {1}}};
  composition_.fingerprint = composition_identity_fingerprint(composition_);
  closure_.element_count = 1;
  closure_.species = {{28, {1}, reverse ? 0. : 100000.},
                      {28, {1}, reverse ? 100000. : 0.}};
  closure_.fingerprint = combustion::chemistry_identity_fingerprint(closure_);
  // SHA256 of the canonical built-in equation specification, not a fuel file.
  gas_.mechanism_sha256 =
      "bd156f5dfa9a7a528525696097015cbf4afad6a5f27bf330c62516169b4dd1c1";
  if (cp == 1200)
    gas_.mechanism_sha256 =
        "7f3b7689de32e783e83d67f3e0dd78efba12ac768fb3b18922ae4aca29ad4ffe";
  gas_.phase =
      "synthetic-isomers-v1-ratebits-" + std::to_string(double_bits(rate));
  gas_.species_names = {composition_.species[0].name,
                        composition_.species[1].name};
  gas_.element_names = {"X"};
  gas_.element_counts = {1, 1};
  gas_.molecular_weights_kg_per_kmol = {28, 28};
  gas_.enthalpy_reference = "absolute-standard-formation-298.15K-v1";
  gas_.composition_fingerprint = composition_.fingerprint;
  gas_.closure_fingerprint = closure_.fingerprint;
}
const CompositionIdentity &AnalyticIsomerBackend::composition() const noexcept {
  return composition_;
}
const portable::GasIdentity &
AnalyticIsomerBackend::gas_identity() const noexcept {
  return gas_;
}
const combustion::ChemistryIdentity &
AnalyticIsomerBackend::closure_identity() const noexcept {
  return closure_;
}
portable::Status
AnalyticIsomerBackend::query_gas(const portable::GasQuery &q,
                                 portable::GasQueryOutput &out) noexcept {
  out.sample = {};
  const auto valid = validate_query(q, gas_, out);
  if (valid != portable::Status::success)
    return valid;
  double t = q.temperature_k;
  if (q.coordinates == portable::GasStateCoordinates::pressure_enthalpy)
    t = 298.15 +
        (q.enthalpy_j_per_kg - 100000 * q.mass_fractions[reactant_]) / cp_;
  else if (q.coordinates != portable::GasStateCoordinates::pressure_temperature)
    return portable::Status::invalid_input;
  if (!std::isfinite(t) || t < 100 || t > 5000)
    return portable::Status::invalid_input;
  portable::GasSample sample{q.revision,
                             gas_.composition_fingerprint,
                             q.pressure_pa,
                             t,
                             q.pressure_pa * 28 / (8314.46261815324 * t),
                             cp_ * (t - 298.15) +
                                 100000 * q.mass_fractions[reactant_],
                             cp_,
                             1e-5,
                             0.025};
  const double rate =
      rate_ * sample.density_kg_per_m3 * q.mass_fractions[reactant_];
  if (!std::isfinite(rate) || !std::isfinite(sample.density_kg_per_m3) ||
      sample.density_kg_per_m3 <= 0)
    return portable::Status::provider_failure;
  for (std::size_t i = 0; i < 2; ++i) {
    out.diffusivities_m2_per_s[i] = 2e-5;
    out.species_enthalpies_j_per_kg[i] =
        cp_ * (t - 298.15) + closure_.species[i].formation_enthalpy_j_per_kg;
    out.net_mass_rates_kg_per_m3_s[i] = i == reactant_ ? -rate : rate;
  }
  out.sample = sample;
  return portable::Status::success;
}
ChemistryIntervalReport
AnalyticIsomerBackend::integrate(const ChemistryIntervalRequest &r) {
  ChemistryIntervalReport out;
  double d[2], h[2], w[2];
  portable::GasQueryOutput sample{{}, d, h, w, 2};
  portable::GasQuery q{{},
                       composition_.fingerprint,
                       portable::GasStateCoordinates::pressure_enthalpy,
                       r.state.p0_pa,
                       r.state.h_tc_j_per_kg,
                       0,
                       r.state.mass_fractions.data(),
                       r.state.mass_fractions.size()};
  if (!std::isfinite(r.start_time_s) || r.start_time_s < 0 ||
      !std::isfinite(r.duration_s) || r.duration_s < 0 ||
      query_gas(q, sample) != portable::Status::success)
    return out;
  out.final_state = r.state;
  const double loss =
      -r.state.mass_fractions[reactant_] * std::expm1(-rate_ * r.duration_s);
  out.final_state.mass_fractions[reactant_] -= loss;
  out.final_state.mass_fractions[product_] += loss;
  out.integrated_rho_y_delta_kg_per_m3.resize(2);
  out.integrated_rho_y_delta_kg_per_m3[reactant_] =
      -sample.sample.density_kg_per_m3 * loss;
  out.integrated_rho_y_delta_kg_per_m3[product_] =
      sample.sample.density_kg_per_m3 * loss;
  out.status = ChemistryStatus::success;
  out.completed_duration_s = r.duration_s;
  out.internal_step_count = r.duration_s > 0 ? 1 : 0;
  return out;
}
portable::Status
AnalyticIsomerBackend::advance_gas(const portable::GasAdvanceQuery &r,
                                   portable::GasAdvanceOutput &out) noexcept {
  out.final_sample = {};
  out.completed_duration_s = 0;
  out.internal_step_count = 0;
  out.integrated_heat_release_j_per_m3 = 0;
  if (out.capacity < 2)
    return portable::Status::capacity_exceeded;
  if (!out.final_mass_fractions ||
      !out.integrated_species_density_delta_kg_per_m3 ||
      r.state.coordinates != portable::GasStateCoordinates::pressure_enthalpy ||
      !std::isfinite(r.start_time_s) || r.start_time_s < 0 ||
      !std::isfinite(r.duration_s) || r.duration_s < 0 ||
      !std::isfinite(r.start_time_s + r.duration_s))
    return portable::Status::invalid_input;
  double d[2], h[2], w[2], candidate[2];
  portable::GasQueryOutput sample{{}, d, h, w, 2};
  auto status = query_gas(r.state, sample);
  if (status != portable::Status::success)
    return status;
  const double rho = sample.sample.density_kg_per_m3;
  const double loss =
      -r.state.mass_fractions[reactant_] * std::expm1(-rate_ * r.duration_s);
  candidate[reactant_] = r.state.mass_fractions[reactant_] - loss;
  candidate[product_] = r.state.mass_fractions[product_] + loss;
  auto q = r.state;
  q.mass_fractions = candidate;
  status = query_gas(q, sample);
  if (status != portable::Status::success)
    return status;
  const double heat = rho * loss * 100000;
  if (!std::isfinite(heat))
    return portable::Status::provider_failure;
  for (std::size_t i = 0; i < 2; ++i) {
    out.final_mass_fractions[i] = candidate[i];
    out.integrated_species_density_delta_kg_per_m3[i] =
        i == reactant_ ? -rho * loss : rho * loss;
  }
  out.final_sample = sample.sample;
  out.completed_duration_s = r.duration_s;
  out.final_sample.enthalpy_j_per_kg = r.state.enthalpy_j_per_kg;
  out.internal_step_count = r.duration_s > 0 ? 1 : 0;
  out.integrated_heat_release_j_per_m3 = heat;
  return portable::Status::success;
}
BackendAdapter::BackendAdapter(ChemistryBackend &b,
                               portable::GasQueryProvider &g,
                               combustion::ChemistryIdentity id,
                               portable::Revision revision)
    : backend_(b), gas_(g), identity_(std::move(id)), revision_(revision) {
  const auto &c = b.composition();
  const auto &a = g.gas_identity();
  // Composition identifies species, not the mechanism or kinetic parameters.
  // Require explicit full-identity capability from the reacting backend too.
  const auto *backend_gas =
      dynamic_cast<const portable::GasQueryProvider *>(&b);
  if (!backend_gas)
    throw std::invalid_argument("chemistry backend lacks full gas identity");
  const auto &actual = backend_gas->gas_identity();
  if (actual.mechanism_sha256 != a.mechanism_sha256 ||
      actual.phase != a.phase || actual.species_names != a.species_names ||
      actual.element_names != a.element_names ||
      actual.element_counts != a.element_counts ||
      actual.molecular_weights_kg_per_kmol != a.molecular_weights_kg_per_kmol ||
      actual.enthalpy_reference != a.enthalpy_reference ||
      actual.composition_fingerprint != a.composition_fingerprint ||
      actual.closure_fingerprint != a.closure_fingerprint)
    throw std::invalid_argument(
        "chemistry backend and gas query asset identity mismatch");
  validate_composition_identity(c);
  const auto n = c.species.size(), ne = c.element_names.size();
  if (revision.algorithm_version != 1 || a.mechanism_sha256.size() != 64 ||
      a.mechanism_sha256.find_first_not_of("0123456789abcdef") !=
          std::string::npos ||
      a.phase.empty() || a.enthalpy_reference.empty() ||
      a.element_names != c.element_names || a.species_names.size() != n ||
      a.molecular_weights_kg_per_kmol.size() != n ||
      a.element_counts.size() != n * ne ||
      a.composition_fingerprint != c.fingerprint ||
      a.closure_fingerprint != identity_.fingerprint ||
      identity_.species.size() != n || identity_.element_count != ne ||
      identity_.fingerprint !=
          combustion::chemistry_identity_fingerprint(identity_))
    throw std::invalid_argument("chemistry identity mapping mismatch");
  for (std::size_t i = 0; i < n; ++i) {
    if (a.species_names[i] != c.species[i].name ||
        a.molecular_weights_kg_per_kmol[i] !=
            c.species[i].molecular_weight_kg_per_kmol ||
        identity_.species[i].molecular_weight_kg_per_kmol !=
            c.species[i].molecular_weight_kg_per_kmol ||
        identity_.species[i].element_counts.size() != ne ||
        !std::isfinite(identity_.species[i].formation_enthalpy_j_per_kg))
      throw std::invalid_argument("chemistry species mapping mismatch");
    for (std::size_t e = 0; e < ne; ++e)
      if (a.element_counts[i * ne + e] !=
              static_cast<std::uint32_t>(c.species[i].element_counts[e]) ||
          identity_.species[i].element_counts[e] !=
              c.species[i].element_counts[e])
        throw std::invalid_argument("chemistry element mapping mismatch");
  }
  diffusion_.resize(n);
  enthalpies_.resize(n);
  rates_.resize(n);
}
const combustion::ChemistryIdentity &BackendAdapter::identity() const noexcept {
  return identity_;
}
combustion::ChemistryRepresentationKind
BackendAdapter::representation() const noexcept {
  return combustion::ChemistryRepresentationKind::direct_cantera;
}
combustion::ChemistryAdvanceReport BackendAdapter::integrate(
    const combustion::ChemistryAdvanceRequest &r) noexcept {
  combustion::ChemistryAdvanceReport out;
  try {
    const auto n = identity_.species.size();
    portable::GasQuery q{revision_,
                         gas_.gas_identity().composition_fingerprint,
                         portable::GasStateCoordinates::pressure_enthalpy,
                         r.state.pressure_pa,
                         r.state.total_thermochemical_enthalpy_j_per_kg,
                         0,
                         r.state.mass_fractions.data(),
                         r.state.mass_fractions.size()};
    portable::GasQueryOutput sample{
        {}, diffusion_.data(), enthalpies_.data(), rates_.data(), n};
    if (gas_.query_gas(q, sample) != portable::Status::success ||
        sample.sample.revision != revision_ ||
        sample.sample.composition_fingerprint != q.composition_fingerprint ||
        !near(sample.sample.density_kg_per_m3, r.state.density_kg_per_m3))
      return out;
    const auto initial_density = sample.sample.density_kg_per_m3;
    ChemistryIntervalRequest request{
        {r.state.pressure_pa, r.state.total_thermochemical_enthalpy_j_per_kg,
         r.state.mass_fractions},
        r.start_time_s,
        r.duration_s};
    const auto advanced = backend_.integrate(request);
    if (advanced.status != ChemistryStatus::success) {
      out.status =
          advanced.status == ChemistryStatus::workspace_failure
              ? combustion::ChemistryAdvanceStatus::workspace_failure
              : combustion::ChemistryAdvanceStatus::integration_failure;
      return out;
    }
    if (advanced.final_state.mass_fractions.size() != n ||
        advanced.integrated_rho_y_delta_kg_per_m3.size() != n ||
        advanced.completed_duration_s != r.duration_s ||
        advanced.final_state.p0_pa != r.state.pressure_pa ||
        advanced.final_state.h_tc_j_per_kg !=
            r.state.total_thermochemical_enthalpy_j_per_kg) {
      out.status = combustion::ChemistryAdvanceStatus::non_finite_output;
      return out;
    }
    double heat = 0, mass = 0, scale = 0;
    for (std::size_t i = 0; i < n; ++i) {
      const double delta = advanced.integrated_rho_y_delta_kg_per_m3[i];
      if (!near(delta,
                initial_density * (advanced.final_state.mass_fractions[i] -
                                   r.state.mass_fractions[i]))) {
        out.status = combustion::ChemistryAdvanceStatus::non_finite_output;
        return out;
      }
      heat -= identity_.species[i].formation_enthalpy_j_per_kg * delta;
      mass += delta;
      scale += std::abs(delta);
    }
    if (!std::isfinite(heat) || std::abs(mass) > 1e-12 + 1e-9 * scale) {
      out.status = combustion::ChemistryAdvanceStatus::non_finite_output;
      return out;
    }
    for (std::size_t e = 0; e < identity_.element_count; ++e) {
      double residual = 0, element_scale = 0;
      for (std::size_t i = 0; i < n; ++i) {
        const double contribution =
            advanced.integrated_rho_y_delta_kg_per_m3[i] *
            identity_.species[i].element_counts[e] /
            identity_.species[i].molecular_weight_kg_per_kmol;
        residual += contribution;
        element_scale += std::abs(contribution);
      }
      if (std::abs(residual) > 1e-12 + 1e-9 * element_scale) {
        out.status = combustion::ChemistryAdvanceStatus::non_finite_output;
        return out;
      }
    }
    q.mass_fractions = advanced.final_state.mass_fractions.data();
    if (gas_.query_gas(q, sample) != portable::Status::success ||
        sample.sample.revision != revision_ ||
        sample.sample.composition_fingerprint != q.composition_fingerprint) {
      out.status = combustion::ChemistryAdvanceStatus::integration_failure;
      return out;
    }
    out.final_state = {r.state.pressure_pa, sample.sample.density_kg_per_m3,
                       r.state.total_thermochemical_enthalpy_j_per_kg,
                       advanced.final_state.mass_fractions};
    out.integrated_species_density_delta_kg_per_m3 =
        advanced.integrated_rho_y_delta_kg_per_m3;
    out.integrated_heat_release_j_per_m3 = heat;
    out.completed_duration_s = advanced.completed_duration_s;
    out.internal_step_count = advanced.internal_step_count;
    out.status = combustion::ChemistryAdvanceStatus::success;
    return out;
  } catch (const std::bad_alloc &) {
    out = {};
    out.status = combustion::ChemistryAdvanceStatus::workspace_failure;
  } catch (...) {
    out = {};
    out.status = combustion::ChemistryAdvanceStatus::integration_failure;
  }
  return out;
}
combustion::ChemicalRateQueryReport
BackendAdapter::query(const combustion::ChemicalRateQuery &r) noexcept {
  combustion::ChemicalRateQueryReport out;
  out.model = r.model;
  out.composition_fingerprint = r.composition_fingerprint;
  if (r.composition_fingerprint != identity_.fingerprint ||
      r.model != combustion::ChemicalTimescaleModel::reactant_depletion_l1_v1) {
    out.status = combustion::ChemicalRateQueryStatus::invalid_input;
    return out;
  }
  try {
    portable::GasQuery q{revision_,
                         gas_.gas_identity().composition_fingerprint,
                         portable::GasStateCoordinates::pressure_enthalpy,
                         r.state.pressure_pa,
                         r.state.total_thermochemical_enthalpy_j_per_kg,
                         0,
                         r.state.mass_fractions.data(),
                         r.state.mass_fractions.size()};
    portable::GasQueryOutput sample{{},
                                    diffusion_.data(),
                                    enthalpies_.data(),
                                    rates_.data(),
                                    rates_.size()};
    if (gas_.query_gas(q, sample) != portable::Status::success ||
        sample.sample.revision != revision_ ||
        sample.sample.composition_fingerprint != q.composition_fingerprint ||
        !near(sample.sample.density_kg_per_m3, r.state.density_kg_per_m3)) {
      out.status = combustion::ChemicalRateQueryStatus::provider_failure;
      return out;
    }
    out.net_species_mass_rates_kg_per_m3_s = rates_;
    out.status = combustion::ChemicalRateQueryStatus::success;
  } catch (...) {
    out.net_species_mass_rates_kg_per_m3_s.clear();
    out.status = combustion::ChemicalRateQueryStatus::provider_failure;
  }
  return out;
}
GasBatchWorkspace::GasBatchWorkspace(std::size_t states, std::size_t species)
    : species_(species) {
  if (!states || !species ||
      species > std::numeric_limits<std::size_t>::max() / 3 ||
      states > std::numeric_limits<std::size_t>::max() / (3 * species))
    throw std::invalid_argument("invalid bounded gas workspace capacity");
  storage_.resize(states * species * 3);
  outputs_.resize(states);
  for (std::size_t i = 0; i < states; ++i) {
    auto *p = storage_.data() + i * species * 3;
    outputs_[i] = {{}, p, p + species, p + 2 * species, species};
  }
}
GasBatchReport GasBatchWorkspace::query(portable::GasQueryProvider &provider,
                                        const portable::GasQuery *inputs,
                                        std::size_t count,
                                        portable::Revision revision) noexcept {
  available_ = false;
  GasBatchReport report;
  for (auto &output : outputs_)
    output.sample = {};
  if (count > outputs_.size()) {
    report.status = portable::Status::capacity_exceeded;
    return report;
  }
  if ((count && !inputs) ||
      provider.gas_identity().species_names.size() != species_)
    return report;
  for (std::size_t i = 0; i < count; ++i) {
    report.failure_index = i;
    if (inputs[i].revision != revision) {
      report.status = portable::Status::stale_revision;
      return report;
    }
    report.status =
        validate_query(inputs[i], provider.gas_identity(), outputs_[i]);
    if (report.status != portable::Status::success)
      return report;
  }
  for (std::size_t i = 0; i < count; ++i) {
    report.failure_index = i;
    report.status = provider.query_gas(inputs[i], outputs_[i]);
    if (report.status != portable::Status::success)
      return report;
    if (outputs_[i].sample.revision != revision ||
        outputs_[i].sample.composition_fingerprint !=
            inputs[i].composition_fingerprint) {
      report.status = portable::Status::provider_failure;
      return report;
    }
  }
  report.status = portable::Status::success;
  report.count = count;
  report.failure_index = count;
  report.available = true;
  available_ = true;
  return report;
}
const portable::GasQueryOutput *GasBatchWorkspace::results() const noexcept {
  return available_ ? outputs_.data() : nullptr;
}
} // namespace hundun::v04::chemistry::detail

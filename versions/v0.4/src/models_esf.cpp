// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09 Philox address fold and balanced Wiener
// explicitly ported from HUNDUN governance 8ffdf2b; no field storage, driver,
// schema or RNG cursor retained.
#include "models_esf_detail.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace hundun::v04::esf::detail {
namespace {
std::size_t checked_capacity(std::size_t species) {
  if (species == 0 || species > std::numeric_limits<std::size_t>::max() / 4 - 1)
    throw std::invalid_argument("invalid ESF prepared species capacity");
  return species;
}
std::uint32_t rot(std::uint32_t x, unsigned n) noexcept {
  return (x << n) | (x >> (32 - n));
}
bool fraction_tuple(const double *y, std::size_t n) noexcept {
  double sum = 0;
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(y[i]) || y[i] < 0 || y[i] > 1)
      return false;
    sum += y[i];
  }
  return std::abs(sum - 1) <= 2e-12;
}
} // namespace
std::array<std::uint32_t, 4>
philox_words(std::array<std::uint32_t, 4> c,
             std::array<std::uint32_t, 2> k) noexcept {
  for (unsigned r = 0; r < 10; ++r) {
    const std::uint64_t a = std::uint64_t(0xd2511f53U) * c[0],
                        b = std::uint64_t(0xcd9e8d57U) * c[2];
    c = {std::uint32_t(b >> 32) ^ c[1] ^ k[0], std::uint32_t(b),
         std::uint32_t(a >> 32) ^ c[3] ^ k[1], std::uint32_t(a)};
    k[0] += 0x9e3779b9U;
    k[1] += 0xbb67ae85U;
  }
  return c;
}
std::array<std::uint32_t, 4> philox(const CounterAddress &a,
                                    std::uint32_t lane) noexcept {
  return philox_words({std::uint32_t(a.accepted_step),
                       std::uint32_t(a.accepted_step >> 32),
                       a.field_pair ^ rot(a.spatial_direction, 11),
                       lane ^ rot(a.stochastic_stage, 17)},
                      {std::uint32_t(a.seed) ^ rot(a.purpose, 5),
                       std::uint32_t(a.seed >> 32) ^ rot(1, 23)});
}
WienerReport balanced_wiener(std::size_t n, double dt,
                             const CounterAddress &a) noexcept {
  WienerReport r;
  if ((n != 2 && n != 4) || !std::isfinite(dt) || dt < 0)
    return r;
  for (std::size_t p = 0; p < n / 2; ++p)
    for (unsigned d = 0; d < 3; ++d) {
      auto x = a;
      x.field_pair = std::uint32_t(p);
      x.spatial_direction = d;
      const double w = dt == 0
                           ? 0
                           : ((philox(x, 0x57494e31U)[0] & 1) ? -std::sqrt(dt)
                                                              : std::sqrt(dt));
      r.increments[2 * p][d] = w;
      r.increments[2 * p + 1][d] = -w;
    }
  r.status = portable::Status::success;
  return r;
}
portable::Status iem_factor(double dt, double tau, double control,
                            double &f) noexcept {
  if (!std::isfinite(dt) || dt < 0 || !std::isfinite(tau) || tau <= 0 ||
      !std::isfinite(control) || control < 0)
    return portable::Status::invalid_input;
  f = (dt == 0 || control == 0)
          ? 1
          : std::exp(-0.5 * std::cbrt(control) * (dt / tau));
  return std::isfinite(f) ? portable::Status::success
                          : portable::Status::invalid_input;
}
portable::Status correct_species_flux(const double *y, const double *raw,
                                      std::size_t n, double *out) noexcept {
  if (!y || !raw || !out || n == 0 || !fraction_tuple(y, n))
    return portable::Status::invalid_input;
  double sum = 0;
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(raw[i]))
      return portable::Status::invalid_input;
    sum += raw[i];
  }
  if (!std::isfinite(sum))
    return portable::Status::invalid_input;
  // Preflight before touching caller output, including overflow.
  for (std::size_t i = 0; i < n; ++i)
    if (!std::isfinite(raw[i] - y[i] * sum))
      return portable::Status::invalid_input;
  for (std::size_t i = 0; i < n; ++i)
    out[i] = raw[i] - y[i] * sum;
  return portable::Status::success;
}
Workspace::Workspace(std::size_t ns)
    : capacity_(checked_capacity(ns)), candidate_(4 * (ns + 1)),
      transported_(4 * (ns + 1)), means_(ns + 1), variances_(ns + 1),
      next_y_(ns), species_delta_(ns), mean_species_delta_(ns) {}
Report Workspace::advance(const Request &q) noexcept {
  Report r;
  const auto &a = q.accepted;
  const auto ns = a.species, n = a.fields, c = ns + 1;
  const bool valid_borrow = !a.owner || a.owner->valid(a);
  ++generation_;
  if (!valid_borrow) {
    r.status = portable::Status::stale_revision;
    return r;
  }
  if (a.revision != q.expected_revision || a.revision.algorithm_version != 1 ||
      q.random.accepted_step != a.revision.accepted_step) {
    r.status = portable::Status::stale_revision;
    return r;
  }
  if (ns > capacity_) {
    r.status = portable::Status::capacity_exceeded;
    return r;
  }
  if (ns == 0 || (n != 2 && n != 4) || !a.values || !q.identity ||
      q.identity->species.size() != ns ||
      a.composition_fingerprint != q.identity->fingerprint ||
      !std::isfinite(q.turbulent_diffusivity_m2_s) ||
      q.turbulent_diffusivity_m2_s < 0 ||
      iem_factor(q.dt_s, q.mixing_time_s, q.tcr_control, r.relaxation_factor) !=
          portable::Status::success)
    return r;
  for (const auto &s : q.identity->species)
    if (!std::isfinite(s.molecular_weight_kg_per_kmol) ||
        s.molecular_weight_kg_per_kmol <= 0 ||
        s.element_counts.size() != q.identity->element_count)
      return r;
  const auto w = balanced_wiener(n, q.dt_s, q.random);
  if (w.status != portable::Status::success)
    return r;
  const double amplitude = std::sqrt(2 * q.turbulent_diffusivity_m2_s);
  for (std::size_t f = 0; f < n; ++f) {
    r.failure_field = f;
    if (!fraction_tuple(a.values + f * c, ns))
      return r;
    for (std::size_t j = 0; j < c; ++j) {
      r.failure_component = j;
      const auto k = f * c + j;
      if (!std::isfinite(a.values[k]))
        return r;
      const double rate = q.deterministic_rates ? q.deterministic_rates[k] : 0;
      if (!std::isfinite(rate))
        return r;
      double value = a.values[k] + q.dt_s * rate;
      if (amplitude != 0 && q.gradients)
        for (unsigned d = 0; d < 3; ++d) {
          const double g = q.gradients[3 * k + d];
          if (!std::isfinite(g))
            return r;
          value += amplitude * g * w.increments[f][d];
        }
      if (!std::isfinite(value))
        return r;
      transported_[k] = value;
    }
    if (!fraction_tuple(transported_.data() + f * c, ns)) {
      r.status = portable::Status::conservation_failure;
      return r;
    }
  }
  for (std::size_t j = 0; j < c; ++j) {
    double mean = 0;
    for (std::size_t f = 0; f < n; ++f)
      mean += transported_[f * c + j] / double(n);
    if (!std::isfinite(mean))
      return r;
    means_[j] = mean;
    double sum = 0, var = 0;
    for (std::size_t f = 0; f < n; ++f) {
      const auto k = f * c + j;
      const double old = transported_[k];
      const double v = r.relaxation_factor == 1
                           ? old
                           : std::fma(r.relaxation_factor, old - mean, mean);
      if (!std::isfinite(v))
        return r;
      candidate_[k] = v;
      sum += v / double(n);
      var += (v - mean) * (v - mean) / double(n);
    }
    r.max_mean_residual = std::max(r.max_mean_residual, std::abs(sum - mean));
    variances_[j] = var;
    if (!std::isfinite(var) ||
        std::abs(sum - mean) > 2e-12 + 2e-12 * std::abs(mean)) {
      r.status = portable::Status::conservation_failure;
      return r;
    }
  }
  for (std::size_t f = 0; f < n; ++f)
    if (!fraction_tuple(candidate_.data() + f * c, ns)) {
      r.status = portable::Status::conservation_failure;
      return r;
    }
  // Check ensemble elements across mixing, not individual-field elements.
  for (std::size_t e = 0; e < q.identity->element_count; ++e) {
    double before = 0, after = 0;
    for (std::size_t f = 0; f < n; ++f)
      for (std::size_t s = 0; s < ns; ++s) {
        const auto &id = q.identity->species[s];
        const double z = double(id.element_counts[e]) /
                         id.molecular_weight_kg_per_kmol / double(n);
        before += transported_[f * c + s] * z;
        after += candidate_[f * c + s] * z;
      }
    r.max_element_residual =
        std::max(r.max_element_residual, std::abs(after - before));
    if (!std::isfinite(after) ||
        std::abs(after - before) > 2e-12 + 2e-12 * std::abs(before)) {
      r.status = portable::Status::conservation_failure;
      return r;
    }
  }
  r.status = portable::Status::success;
  r.candidate = a;
  r.candidate.values = candidate_.data();
  r.candidate.owner = this;
  r.candidate.generation = generation_;
  r.means = means_.data();
  r.variances = variances_.data();
  return r;
}
Report Workspace::react(const ReactionRequest &q,
                        portable::GasAdvanceProvider &provider) noexcept {
  Report r;
  r.model_identity = 0x4553465245410001ULL;
  const auto &a = q.accepted;
  const auto ns = a.species, n = a.fields, c = ns + 1;
  const bool valid_borrow = !a.owner || a.owner->valid(a);
  ++generation_;
  if (!valid_borrow || a.revision != q.expected_revision ||
      a.revision.algorithm_version != 1) {
    r.status = portable::Status::stale_revision;
    return r;
  }
  if (ns > capacity_) {
    r.status = portable::Status::capacity_exceeded;
    return r;
  }
  if (ns == 0 || (n != 2 && n != 4) || !a.values || !q.chemistry_identity ||
      !q.gas_identity || !q.pressures_pa || !q.initial_densities_kg_per_m3 ||
      !std::isfinite(q.start_time_s) || q.start_time_s < 0 ||
      !std::isfinite(q.duration_s) || q.duration_s <= 0 ||
      !std::isfinite(q.start_time_s + q.duration_s) || q.duration_s / 2 <= 0 ||
      q.start_time_s + q.duration_s / 2 <= q.start_time_s)
    return r;
  const auto &id = *q.chemistry_identity;
  const auto &gas = *q.gas_identity;
  if (!portable::same_gas_identity(gas, provider.gas_identity())) {
    r.status = portable::Status::identity_mismatch;
    return r;
  }
  if (a.composition_fingerprint != id.fingerprint ||
      gas.closure_fingerprint != id.fingerprint ||
      gas.composition_fingerprint == 0 || gas.mechanism_sha256.empty() ||
      gas.phase.empty() || gas.enthalpy_reference.empty() ||
      id.species.size() != ns || gas.species_names.size() != ns ||
      gas.molecular_weights_kg_per_kmol.size() != ns ||
      gas.element_names.size() != id.element_count ||
      gas.element_counts.size() != ns * id.element_count) {
    r.status = portable::Status::identity_mismatch;
    return r;
  }
  for (std::size_t s = 0; s < ns; ++s) {
    const auto &si = id.species[s];
    if (!std::isfinite(si.molecular_weight_kg_per_kmol) ||
        si.molecular_weight_kg_per_kmol <= 0 ||
        !std::isfinite(si.formation_enthalpy_j_per_kg) ||
        si.element_counts.size() != id.element_count ||
        gas.molecular_weights_kg_per_kmol[s] !=
            si.molecular_weight_kg_per_kmol) {
      r.status = portable::Status::identity_mismatch;
      return r;
    }
    for (std::size_t e = 0; e < id.element_count; ++e)
      if (gas.element_counts[s * id.element_count + e] !=
          si.element_counts[e]) {
        r.status = portable::Status::identity_mismatch;
        return r;
      }
  }
  const auto near = [](double x, double y) {
    return std::isfinite(x) && std::isfinite(y) &&
           std::abs(x - y) <= 1e-10 + 2e-12 * std::abs(y);
  };
  std::fill(mean_species_delta_.begin(), mean_species_delta_.end(), 0.0);
  for (std::size_t f = 0; f < n; ++f) {
    r.failure_field = f;
    if (!fraction_tuple(a.values + f * c, ns) ||
        !std::isfinite(a.values[f * c + ns]) ||
        !std::isfinite(q.pressures_pa[f]) || q.pressures_pa[f] <= 0 ||
        !std::isfinite(q.initial_densities_kg_per_m3[f]) ||
        q.initial_densities_kg_per_m3[f] <= 0)
      return r;
    for (std::size_t j = 0; j < c; ++j)
      transported_[f * c + j] = a.values[f * c + j];
  }
  for (std::size_t f = 0; f < n; ++f) {
    r.failure_field = f;
    double density = q.initial_densities_kg_per_m3[f];
    for (unsigned half = 0; half < 2; ++half) {
      portable::GasAdvanceQuery advance{
          {a.revision, gas.composition_fingerprint,
           portable::GasStateCoordinates::pressure_enthalpy, q.pressures_pa[f],
           transported_[f * c + ns], 0, transported_.data() + f * c, ns},
          q.start_time_s + double(half) * q.duration_s / 2,
          q.duration_s / 2};
      portable::GasAdvanceOutput out{
          {}, next_y_.data(), species_delta_.data(), ns};
      ++r.chemistry_call_count;
      const auto status = provider.advance_gas(advance, out);
      if (status != portable::Status::success) {
        r.status = status;
        return r;
      }
      const auto &sample = out.final_sample;
      if (sample.revision != a.revision ||
          sample.composition_fingerprint != gas.composition_fingerprint ||
          !near(sample.pressure_pa, q.pressures_pa[f]) ||
          !near(sample.enthalpy_j_per_kg, advance.state.enthalpy_j_per_kg) ||
          !std::isfinite(sample.density_kg_per_m3) ||
          sample.density_kg_per_m3 <= 0 ||
          !std::isfinite(sample.temperature_k) || sample.temperature_k <= 0 ||
          !std::isfinite(sample.cp_j_per_kg_k) || sample.cp_j_per_kg_k <= 0 ||
          out.completed_duration_s != advance.duration_s ||
          out.final_mass_fractions != next_y_.data() ||
          out.integrated_species_density_delta_kg_per_m3 !=
              species_delta_.data() ||
          out.capacity != ns || !fraction_tuple(next_y_.data(), ns)) {
        r.status = portable::Status::provider_failure;
        return r;
      }
      double heat = 0, mass = 0;
      for (std::size_t s = 0; s < ns; ++s) {
        r.failure_component = s;
        const double expected =
            density * (next_y_[s] - transported_[f * c + s]);
        if (!near(species_delta_[s], expected)) {
          r.status = portable::Status::conservation_failure;
          return r;
        }
        heat -= id.species[s].formation_enthalpy_j_per_kg * species_delta_[s];
        mass += species_delta_[s];
      }
      if (!near(mass, 0) || !near(out.integrated_heat_release_j_per_m3, heat)) {
        r.status = portable::Status::conservation_failure;
        return r;
      }
      for (std::size_t e = 0; e < id.element_count; ++e) {
        double delta = 0;
        for (std::size_t s = 0; s < ns; ++s)
          delta += species_delta_[s] * double(id.species[s].element_counts[e]) /
                   id.species[s].molecular_weight_kg_per_kmol;
        if (!near(delta, 0)) {
          r.status = portable::Status::conservation_failure;
          return r;
        }
      }
      for (std::size_t s = 0; s < ns; ++s) {
        mean_species_delta_[s] += species_delta_[s] / double(n);
        if (!std::isfinite(mean_species_delta_[s])) {
          r.status = portable::Status::conservation_failure;
          return r;
        }
        transported_[f * c + s] = next_y_[s];
      }
      // Closed chemistry conserves its PH enthalpy coordinate exactly; the
      // validated EOS round-trip must not inject floating-point drift.
      density = sample.density_kg_per_m3;
      r.ensemble_heat_release_j_per_m3 += heat / double(n);
    }
    r.final_densities_kg_per_m3[f] = density;
  }
  for (std::size_t j = 0; j < c; ++j) {
    double mean = 0, var = 0;
    for (std::size_t f = 0; f < n; ++f)
      mean += transported_[f * c + j] / double(n);
    for (std::size_t f = 0; f < n; ++f) {
      const double d = transported_[f * c + j] - mean;
      var += d * d / double(n);
      candidate_[f * c + j] = transported_[f * c + j];
    }
    if (!std::isfinite(mean) || !std::isfinite(var) ||
        !std::isfinite(r.ensemble_heat_release_j_per_m3))
      return r;
    means_[j] = mean;
    variances_[j] = var;
  }
  r.status = portable::Status::success;
  r.candidate = a;
  r.candidate.values = candidate_.data();
  r.candidate.owner = this;
  r.candidate.generation = generation_;
  r.means = means_.data();
  r.variances = variances_.data();
  r.mean_integrated_species_density_delta_kg_per_m3 = mean_species_delta_.data();
  return r;
}
} // namespace hundun::v04::esf::detail

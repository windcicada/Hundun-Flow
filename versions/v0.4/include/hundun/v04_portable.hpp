// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace hundun::v04::portable {

// Values only: no field, product transaction, schema or persistence format.
struct Revision {
  std::uint64_t accepted_step{};
  std::uint64_t input_revision{};
  std::uint32_t algorithm_version{1};
};
inline bool operator==(const Revision& a, const Revision& b) noexcept {
  return a.accepted_step == b.accepted_step &&
         a.input_revision == b.input_revision &&
         a.algorithm_version == b.algorithm_version;
}
inline bool operator!=(const Revision& a, const Revision& b) noexcept {
  return !(a == b);
}

enum class Status : std::uint8_t {
  success, invalid_input, identity_mismatch, stale_revision,
  capacity_exceeded, provider_failure, unavailable, conservation_failure
};

// Construct/validate identity in preparation. Existing composition and closure
// fingerprints retain their meanings; neither substitutes for asset identity.
struct GasIdentity {
  std::string mechanism_sha256;
  std::string phase;
  std::vector<std::string> species_names;
  std::vector<std::string> element_names;
  std::vector<std::uint32_t> element_counts; // species-major [Ns, Ne]
  std::vector<double> molecular_weights_kg_per_kmol;
  std::string enthalpy_reference;
  std::uint64_t composition_fingerprint{};
  std::uint64_t closure_fingerprint{};
};

enum class GasStateCoordinates : std::uint8_t { pressure_enthalpy, pressure_temperature };
struct GasQuery {
  Revision revision{};
  std::uint64_t composition_fingerprint{};
  GasStateCoordinates coordinates{GasStateCoordinates::pressure_enthalpy};
  double pressure_pa{};
  double enthalpy_j_per_kg{};
  double temperature_k{};
  const double* mass_fractions{};
  std::size_t species_count{};
};
struct GasSample {
  Revision revision{};
  std::uint64_t composition_fingerprint{};
  double pressure_pa{}, temperature_k{}, density_kg_per_m3{};
  double enthalpy_j_per_kg{}, cp_j_per_kg_k{};
  double viscosity_pa_s{}, conductivity_w_per_m_k{};
};

// Caller-owned outputs of exactly Ns elements. Prepare all HUNDUN capacity
// before querying. On failure scalar sample is unavailable and arrays must not
// be consumed. Providers report third-party allocations separately.
struct GasQueryOutput {
  GasSample sample{};
  double* diffusivities_m2_per_s{};
  double* species_enthalpies_j_per_kg{};
  double* net_mass_rates_kg_per_m3_s{};
  std::size_t capacity{};
};
class GasQueryProvider {
public:
  virtual ~GasQueryProvider() = default;
  virtual const GasIdentity& gas_identity() const noexcept = 0;
  virtual Status query_gas(const GasQuery&, GasQueryOutput&) noexcept = 0;
};

// All exchange values are interval-integrated extensive totals, new minus old.
// A parcel exchange already includes multiplicity. Deposition weight is applied
// exactly once by the cell reduction. Rates and densities are different types
// in the consuming modules and must be converted explicitly there.
struct ExchangeDelta {
  double mass_kg{};
  double momentum_kg_m_per_s[3]{};
  double thermochemical_enthalpy_j{};
  double kinetic_energy_j{};
};

} // namespace hundun::v04::portable

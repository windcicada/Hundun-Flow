// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace hundun::chemistry {

struct SpeciesIdentity final {
  std::string name;
  double molecular_weight_kg_per_kmol{};
  std::vector<std::int32_t> element_counts;
};

struct CompositionIdentity final {
  std::vector<std::string> element_names;
  std::vector<SpeciesIdentity> species;
  std::uint64_t fingerprint{};
};

struct ThermochemicalPoint final {
  double p0_pa{};
  double h_tc_j_per_kg{};
  std::vector<double> mass_fractions;
};

struct ThermodynamicProperties final {
  double temperature_k{};
  double density_kg_per_m3{};
  double cp_j_per_kg_k{};
  double mixture_molecular_weight_kg_per_kmol{};
};

struct TransportProperties final {
  double viscosity_pa_s{};
  double conductivity_w_per_m_k{};
  std::vector<double> mixture_diffusivity_m2_per_s;
};

namespace detail {

inline constexpr std::uint64_t chemistry_fnv_offset =
    UINT64_C(14695981039346656037);
inline constexpr std::uint64_t chemistry_fnv_prime = UINT64_C(1099511628211);

inline void fingerprint_byte(std::uint64_t &hash, std::uint8_t value) noexcept {
  hash ^= static_cast<std::uint64_t>(value);
  hash *= chemistry_fnv_prime;
}

inline void fingerprint_u64(std::uint64_t &hash, std::uint64_t value) noexcept {
  for (std::uint32_t shift = 0U; shift < 64U; shift += 8U) {
    fingerprint_byte(hash,
                     static_cast<std::uint8_t>((value >> shift) & UINT64_C(0xff)));
  }
}

inline void fingerprint_string(std::uint64_t &hash,
                               std::string_view value) noexcept {
  fingerprint_u64(hash, static_cast<std::uint64_t>(value.size()));
  for (const char character : value) {
    fingerprint_byte(
        hash, static_cast<std::uint8_t>(static_cast<unsigned char>(character)));
  }
}

inline std::uint64_t double_bits(double value) noexcept {
  std::uint64_t bits{};
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

} // namespace detail

inline std::uint64_t composition_identity_fingerprint(
    const CompositionIdentity &identity) noexcept {
  std::uint64_t hash = detail::chemistry_fnv_offset;
  detail::fingerprint_u64(
      hash, static_cast<std::uint64_t>(identity.element_names.size()));
  for (const auto &element : identity.element_names) {
    detail::fingerprint_string(hash, element);
  }
  detail::fingerprint_u64(hash,
                          static_cast<std::uint64_t>(identity.species.size()));
  for (const auto &species : identity.species) {
    detail::fingerprint_string(hash, species.name);
    detail::fingerprint_u64(
        hash, detail::double_bits(species.molecular_weight_kg_per_kmol));
    detail::fingerprint_u64(
        hash, static_cast<std::uint64_t>(species.element_counts.size()));
    for (const auto count : species.element_counts) {
      detail::fingerprint_u64(
          hash, static_cast<std::uint64_t>(static_cast<std::int64_t>(count)));
    }
  }
  return hash;
}

inline void validate_composition_identity(const CompositionIdentity &identity) {
  if (identity.element_names.empty()) {
    throw std::invalid_argument("composition identity requires elements");
  }
  if (identity.species.empty()) {
    throw std::invalid_argument("composition identity requires species");
  }
  for (std::size_t element = 0; element < identity.element_names.size();
       ++element) {
    if (identity.element_names[element].empty()) {
      throw std::invalid_argument("composition element name must be nonempty");
    }
    if (std::find(identity.element_names.begin(),
                  identity.element_names.begin() +
                      static_cast<std::ptrdiff_t>(element),
                  identity.element_names[element]) !=
        identity.element_names.begin() +
            static_cast<std::ptrdiff_t>(element)) {
      throw std::invalid_argument("composition element names must be unique");
    }
  }
  for (std::size_t species_index = 0; species_index < identity.species.size();
       ++species_index) {
    const auto &species = identity.species[species_index];
    if (species.name.empty()) {
      throw std::invalid_argument("composition species name must be nonempty");
    }
    const auto duplicate = std::find_if(
        identity.species.begin(),
        identity.species.begin() + static_cast<std::ptrdiff_t>(species_index),
        [&](const SpeciesIdentity &candidate) {
          return candidate.name == species.name;
        });
    if (duplicate != identity.species.begin() +
                         static_cast<std::ptrdiff_t>(species_index)) {
      throw std::invalid_argument("composition species names must be unique");
    }
    if (!std::isfinite(species.molecular_weight_kg_per_kmol) ||
        species.molecular_weight_kg_per_kmol <= 0.0) {
      throw std::invalid_argument(
          "composition molecular weight must be finite and positive");
    }
    if (species.element_counts.size() != identity.element_names.size()) {
      throw std::invalid_argument("composition element matrix shape mismatch");
    }
    if (std::any_of(species.element_counts.begin(), species.element_counts.end(),
                    [](std::int32_t count) { return count < 0; })) {
      throw std::invalid_argument(
          "composition element counts must be nonnegative");
    }
  }
  if (identity.fingerprint == 0U ||
      identity.fingerprint != composition_identity_fingerprint(identity)) {
    throw std::invalid_argument("composition fingerprint mismatch");
  }
}

inline double mixture_total_thermochemical_enthalpy_j_per_kg(
    const std::vector<double> &mass_fractions,
    const std::vector<double> &total_species_enthalpy_j_per_kg) {
  if (mass_fractions.empty() ||
      mass_fractions.size() != total_species_enthalpy_j_per_kg.size()) {
    throw std::invalid_argument(
        "thermochemical enthalpy inputs must have one value per species");
  }
  double mass_fraction_sum = 0.0;
  double mixture_enthalpy = 0.0;
  for (std::size_t species = 0; species < mass_fractions.size(); ++species) {
    const double mass_fraction = mass_fractions[species];
    const double species_enthalpy = total_species_enthalpy_j_per_kg[species];
    if (!std::isfinite(mass_fraction) || mass_fraction < 0.0 ||
        !std::isfinite(species_enthalpy)) {
      throw std::invalid_argument(
          "thermochemical enthalpy inputs must be finite and nonnegative in Y");
    }
    mass_fraction_sum += mass_fraction;
    mixture_enthalpy += mass_fraction * species_enthalpy;
  }
  const double tolerance =
      1.0e-12 * static_cast<double>(mass_fractions.size());
  if (!std::isfinite(mass_fraction_sum) ||
      std::abs(mass_fraction_sum - 1.0) > tolerance ||
      !std::isfinite(mixture_enthalpy)) {
    throw std::invalid_argument(
        "thermochemical enthalpy mass fractions must sum to one");
  }
  return mixture_enthalpy;
}

} // namespace hundun::chemistry

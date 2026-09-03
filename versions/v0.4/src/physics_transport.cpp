// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_physics.hpp"

#include "physics_input_detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <type_traits>
#include <utility>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kTransportInput = 2101U;
constexpr std::uint32_t kTransportSpecies = 2102U;
constexpr std::uint32_t kTransportComposition = 2103U;
constexpr std::uint32_t kTransportTemperature = 2104U;
constexpr std::uint32_t kTransportNumerical = 2105U;
constexpr std::size_t kMaximumSpecies = 65U;
constexpr std::uint64_t kFnvOffset = UINT64_C(14695981039346656037);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);
constexpr double kCoastNativeAirPrandtl = 0.70;
constexpr std::array<double, 2U> kCoastNativeAirMoleFraction{0.21, 0.79};
constexpr std::array<double, 2U> kCoastNativeAirSpeciesMolecularWeight{
    31.9988, 28.0134};
constexpr std::array<double, 2U> kCoastNativeAirCriticalPressure{49.7, 33.5};
constexpr std::array<double, 2U> kCoastNativeAirCriticalTemperature{154.4,
                                                                    126.2};

static_assert(std::is_nothrow_move_assignable_v<TransportPlan>,
              "transport plan publication must not throw");
static_assert(sizeof(double) == sizeof(std::uint64_t) &&
                  std::numeric_limits<double>::is_iec559,
              "transport identity requires IEEE-754 binary64");

class Hash64 {
 public:
  void bytes(const void* data, std::size_t size) noexcept {
    const auto* input = static_cast<const unsigned char*>(data);
    for (std::size_t index = 0U; index < size; ++index) {
      value_ ^= static_cast<std::uint64_t>(input[index]);
      value_ *= kFnvPrime;
    }
  }

  template <class Integer>
  void integer(Integer value) noexcept {
    using Unsigned = std::make_unsigned_t<Integer>;
    Unsigned bits = static_cast<Unsigned>(value);
    for (std::size_t byte = 0U; byte < sizeof(bits); ++byte) {
      const auto part = static_cast<unsigned char>(
          (bits >> (byte * 8U)) & static_cast<Unsigned>(0xffU));
      bytes(&part, 1U);
    }
  }

  void real(double value) noexcept {
    std::uint64_t bits = 0U;
    std::memcpy(&bits, &value, sizeof(bits));
    integer(bits);
  }

  PlanFingerprint finish() const noexcept { return value_ == 0U ? 1U : value_; }

 private:
  std::uint64_t value_{kFnvOffset};
};

bool finite_positive(double value) noexcept {
  return std::isfinite(value) && value > 0.0;
}

bool composition(Span<const double> independent,
                 Span<const std::uint16_t> mapping,
                 std::size_t dependent_species, std::size_t species_count,
                 std::array<double, kMaximumSpecies>& mass_fraction) noexcept {
  if (species_count == 0U || species_count > kMaximumSpecies ||
      independent.size != mapping.size ||
      (independent.size != 0U && independent.data == nullptr)) {
    return false;
  }
  mass_fraction.fill(0.0);
  double sum = 0.0;
  for (std::size_t index = 0U; index < independent.size; ++index) {
    const double value = independent.data[index];
    const std::size_t species = mapping.data[index];
    if (!std::isfinite(value) || value < 0.0 || species >= species_count ||
        species == dependent_species || value > 1.0 - sum) {
      return false;
    }
    sum += value;
    mass_fraction[species] = value;
  }
  const double dependent = 1.0 - sum;
  if (!std::isfinite(dependent) || dependent < 0.0) {
    return false;
  }
  mass_fraction[dependent_species] = dependent;
  return true;
}

double species_cp(const std::array<std::vector<double>, 7U>& nasa_low,
                  const std::array<std::vector<double>, 7U>& nasa_high,
                  Span<const double> temperature_switch,
                  Span<const double> molecular_weight, std::size_t species,
                  double temperature,
                  double universal_gas_constant =
                      kUniversalGasConstant) noexcept {
  const auto& coefficients = temperature <= temperature_switch.data[species]
                                 ? nasa_low
                                 : nasa_high;
  double cp_over_r = coefficients[4U][species];
  for (std::size_t coefficient = 4U; coefficient-- > 0U;) {
    cp_over_r = coefficients[coefficient][species] + cp_over_r * temperature;
  }
  return universal_gas_constant * cp_over_r /
         molecular_weight.data[species];
}

// Independent implementation of the published Yoon--Thodos corresponding-
// states low-pressure gas-viscosity correlation, followed by Wilke gas-
// mixture weighting.  Tc is in K, Pc in atm, M in g/mol, the
// correlation result is in micropoise, and 1 micropoise = 1e-7 Pa s.
double coast_native_air_species_viscosity(std::size_t species,
                                          double temperature) noexcept {
  const double critical_temperature =
      kCoastNativeAirCriticalTemperature[species];
  const double critical_pressure = kCoastNativeAirCriticalPressure[species];
  const double molecular_weight =
      kCoastNativeAirSpeciesMolecularWeight[species];
  const double c1 =
      std::pow(critical_temperature, 1.0 / 6.0) /
      (std::sqrt(molecular_weight) *
       std::pow(critical_pressure, 2.0 / 3.0));
  const double reduced_temperature = temperature / critical_temperature;
  const double viscosity_micro_poises =
      4.610 * std::pow(reduced_temperature, 0.618) -
      2.04 * std::exp(-0.449 * reduced_temperature) +
      1.94 * std::exp(-4.058 * reduced_temperature) + 0.1;
  return viscosity_micro_poises / c1 * 1.0e-7;
}

double coast_native_air_viscosity(double temperature) noexcept {
  std::array<double, 2U> viscosity{};
  for (std::size_t species = 0U; species < viscosity.size(); ++species) {
    viscosity[species] =
        coast_native_air_species_viscosity(species, temperature);
  }
  double mixture = 0.0;
  for (std::size_t first = 0U; first < viscosity.size(); ++first) {
    double denominator = 0.0;
    for (std::size_t second = 0U; second < viscosity.size(); ++second) {
      const double phi =
          (1.0 / std::sqrt(8.0)) /
          std::sqrt(1.0 + kCoastNativeAirSpeciesMolecularWeight[first] /
                              kCoastNativeAirSpeciesMolecularWeight[second]) *
          std::pow(1.0 +
                       std::sqrt(viscosity[first] / viscosity[second]) *
                           std::pow(
                               kCoastNativeAirSpeciesMolecularWeight[second] /
                                   kCoastNativeAirSpeciesMolecularWeight[first],
                               0.25),
                   2.0);
      denominator += kCoastNativeAirMoleFraction[second] * phi;
    }
    mixture += kCoastNativeAirMoleFraction[first] * viscosity[first] /
               denominator;
  }
  return mixture;
}

double species_viscosity(Span<const double> viscosity_reference,
                         Span<const double> reference_temperature,
                         Span<const double> sutherland_temperature,
                         std::size_t species, double temperature) noexcept {
  if (reference_temperature.data[species] == 0.0) {
    return viscosity_reference.data[species];
  }
  const double reference = reference_temperature.data[species];
  const double sutherland = sutherland_temperature.data[species];
  const double ratio = temperature / reference;
  return viscosity_reference.data[species] * ratio * std::sqrt(ratio) *
         (reference + sutherland) / (temperature + sutherland);
}

double species_conductivity(
    const std::array<std::vector<double>, 7U>& nasa_low,
    const std::array<std::vector<double>, 7U>& nasa_high,
    Span<const double> temperature_switch,
    Span<const double> molecular_weight,
    Span<const double> reference_temperature, Span<const double> prandtl,
    Span<const double> constant_conductivity, std::size_t species,
    double temperature, double viscosity) noexcept {
  if (reference_temperature.data[species] == 0.0) {
    return constant_conductivity.data[species];
  }
  return viscosity *
         species_cp(nasa_low, nasa_high, temperature_switch,
                    molecular_weight, species, temperature) /
         prandtl.data[species];
}

}  // namespace

Status TransportPlan::compile(const ThermophysicalSpec& spec,
                              const ThermodynamicsPlan& thermodynamics,
                              TransportPlan& out) noexcept {
  const std::size_t count = spec.species.size();
  if (count == 0U || count > kMaximumSpecies ||
      thermodynamics.fingerprint_ == 0U ||
      thermodynamics.source_fingerprint_ == 0U ||
      detail::thermophysical_spec_fingerprint(spec) !=
          thermodynamics.source_fingerprint_ ||
      thermodynamics.inverse_molecular_weight_.size() != count ||
      thermodynamics.temperature_switch_.size() != count ||
      thermodynamics.dependent_species_ >= count ||
      thermodynamics.independent_to_species_.size() + 1U != count) {
    return {StatusCode::invalid_plan, kTransportInput};
  }
  ThermophysicalSpec canonical_spec;
  try {
    canonical_spec = spec;
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, kTransportInput};
  } catch (...) {
    return {StatusCode::invalid_plan, kTransportInput};
  }
  if (!detail::canonicalize_thermophysical_spec(canonical_spec)) {
    return {StatusCode::invalid_plan, kTransportInput};
  }
  std::array<std::uint8_t, kMaximumSpecies> seen{};
  seen[thermodynamics.dependent_species_] = 1U;
  for (const std::uint16_t mapped :
       thermodynamics.independent_to_species_) {
    const std::size_t index = static_cast<std::size_t>(mapped);
    if (index >= count || seen[index] != 0U) {
      return {StatusCode::invalid_plan, kTransportInput};
    }
    seen[index] = 1U;
  }
  if (!std::all_of(seen.begin(), seen.begin() + count,
                   [](std::uint8_t value) { return value == 1U; })) {
    return {StatusCode::invalid_plan, kTransportInput};
  }

  TransportPlan candidate;
  try {
    candidate.molecular_weight_.resize(count);
    candidate.viscosity_reference_.resize(count);
    candidate.reference_temperature_.resize(count);
    candidate.sutherland_temperature_.resize(count);
    candidate.prandtl_.resize(count);
    candidate.conductivity_.resize(count);
    candidate.temperature_switch_.resize(count);
    candidate.molecular_weight_ratio_quarter_.resize(count * count);
    candidate.wilke_denominator_reciprocal_.resize(count * count);
    candidate.independent_to_species_ =
        thermodynamics.independent_to_species_;
    for (auto& coefficients : candidate.nasa_low_) {
      coefficients.resize(count);
    }
    for (auto& coefficients : candidate.nasa_high_) {
      coefficients.resize(count);
    }
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, kTransportInput};
  } catch (...) {
    return {StatusCode::invalid_plan, kTransportInput};
  }

  bool any_sutherland = false;
  bool coast_native_air = false;
  for (std::size_t species = 0U; species < count; ++species) {
    const SpeciesThermophysicalSpec& input =
        canonical_spec.species[species];
    if (!finite_positive(input.molecular_weight) ||
        (input.transport_law != TransportLaw::coast_native_air &&
         !finite_positive(input.viscosity_reference)) ||
        !finite_positive(input.temperature_switch) ||
        std::abs(input.molecular_weight *
                     thermodynamics.inverse_molecular_weight_[species] -
                 1.0) >
            4.0 * std::numeric_limits<double>::epsilon() ||
        input.temperature_switch !=
            thermodynamics.temperature_switch_[species]) {
      return {StatusCode::invalid_plan, kTransportSpecies};
    }
    candidate.molecular_weight_[species] = input.molecular_weight;
    candidate.viscosity_reference_[species] = input.viscosity_reference;
    candidate.temperature_switch_[species] = input.temperature_switch;
    for (std::size_t coefficient = 0U; coefficient < 7U; ++coefficient) {
      const double low = input.nasa7_low[coefficient];
      const double high = input.nasa7_high[coefficient];
      if (!std::isfinite(low) || !std::isfinite(high)) {
        return {StatusCode::invalid_plan, kTransportSpecies};
      }
      candidate.nasa_low_[coefficient][species] = low;
      candidate.nasa_high_[coefficient][species] = high;
    }

    switch (input.transport_law) {
      case TransportLaw::constant:
        if (!finite_positive(input.conductivity)) {
          return {StatusCode::invalid_plan, kTransportSpecies};
        }
        candidate.reference_temperature_[species] = 0.0;
        candidate.sutherland_temperature_[species] = 0.0;
        candidate.prandtl_[species] = 0.0;
        candidate.conductivity_[species] = input.conductivity;
        break;
      case TransportLaw::sutherland:
        if (!finite_positive(input.transport_reference_temperature) ||
            !std::isfinite(input.sutherland_temperature) ||
            input.sutherland_temperature < 0.0 ||
            !finite_positive(input.prandtl)) {
          return {StatusCode::invalid_plan, kTransportSpecies};
        }
        candidate.reference_temperature_[species] =
            input.transport_reference_temperature;
        candidate.sutherland_temperature_[species] =
            input.sutherland_temperature;
        candidate.prandtl_[species] = input.prandtl;
        candidate.conductivity_[species] = 0.0;
        any_sutherland = true;
        break;
      case TransportLaw::coast_native_air:
        if (count != 1U || input.viscosity_reference != 0.0 ||
            input.transport_reference_temperature != 0.0 ||
            input.sutherland_temperature != 0.0 || input.prandtl != 0.0 ||
            input.conductivity != 0.0) {
          return {StatusCode::invalid_plan, kTransportSpecies};
        }
        candidate.reference_temperature_[species] = 0.0;
        candidate.sutherland_temperature_[species] = 0.0;
        candidate.prandtl_[species] = 0.0;
        candidate.conductivity_[species] = 0.0;
        coast_native_air = true;
        break;
      default:
        return {StatusCode::invalid_plan, kTransportSpecies};
    }
  }

  candidate.dependent_species_ = thermodynamics.dependent_species_;
  candidate.minimum_temperature_ = thermodynamics.minimum_temperature_;
  candidate.maximum_temperature_ = thermodynamics.maximum_temperature_;
  candidate.kernel_ = coast_native_air
                          ? TransportKernel::coast_native_air
                          : (any_sutherland
                                 ? TransportKernel::sutherland_wilke
                                 : TransportKernel::constant);
  candidate.enthalpy_prandtl_ =
      coast_native_air ? kCoastNativeAirPrandtl : 0.0;
  if (candidate.kernel_ == TransportKernel::constant) {
    try {
      candidate.constant_wilke_phi_.resize(count * count);
    } catch (const std::bad_alloc&) {
      return {StatusCode::allocation_failure, kTransportInput};
    } catch (...) {
      return {StatusCode::invalid_plan, kTransportInput};
    }
  }
  for (std::size_t first = 0U; first < count; ++first) {
    for (std::size_t second = 0U; second < count; ++second) {
      const std::size_t interaction = first * count + second;
      const double molecular_weight_ratio =
          candidate.molecular_weight_[second] /
          candidate.molecular_weight_[first];
      const double ratio_quarter =
          std::sqrt(std::sqrt(molecular_weight_ratio));
      const double denominator_reciprocal =
          1.0 / std::sqrt(8.0 * (1.0 + 1.0 / molecular_weight_ratio));
      if (!finite_positive(ratio_quarter) ||
          !finite_positive(denominator_reciprocal)) {
        return {StatusCode::invalid_plan, kTransportSpecies};
      }
      candidate.molecular_weight_ratio_quarter_[interaction] = ratio_quarter;
      candidate.wilke_denominator_reciprocal_[interaction] =
          denominator_reciprocal;
      if (candidate.kernel_ == TransportKernel::constant) {
        const double numerator =
            1.0 +
            std::sqrt(candidate.viscosity_reference_[first] /
                      candidate.viscosity_reference_[second]) *
                ratio_quarter;
        const double phi =
            numerator * numerator * denominator_reciprocal;
        if (!finite_positive(phi)) {
          return {StatusCode::invalid_plan, kTransportSpecies};
        }
        candidate.constant_wilke_phi_[interaction] = phi;
      }
    }
  }
  Hash64 hash;
  hash.integer(static_cast<std::uint8_t>(candidate.kernel_));
  hash.integer(thermodynamics.fingerprint_);
  hash.integer(static_cast<std::uint64_t>(count));
  hash.integer(static_cast<std::uint64_t>(candidate.dependent_species_));
  hash.real(candidate.minimum_temperature_);
  hash.real(candidate.maximum_temperature_);
  hash.real(candidate.enthalpy_prandtl_);
  hash.integer(static_cast<std::uint64_t>(
      candidate.constant_wilke_phi_.size()));
  for (std::size_t first = 0U; first < count; ++first) {
    for (std::size_t second = 0U; second < count; ++second) {
      const std::size_t interaction = first * count + second;
      hash.real(candidate.molecular_weight_ratio_quarter_[interaction]);
      hash.real(candidate.wilke_denominator_reciprocal_[interaction]);
      if (candidate.kernel_ == TransportKernel::constant) {
        hash.real(candidate.constant_wilke_phi_[interaction]);
      }
    }
  }
  for (std::size_t species = 0U; species < count; ++species) {
    hash.real(candidate.molecular_weight_[species]);
    hash.real(candidate.viscosity_reference_[species]);
    hash.real(candidate.reference_temperature_[species]);
    hash.real(candidate.sutherland_temperature_[species]);
    hash.real(candidate.prandtl_[species]);
    hash.real(candidate.conductivity_[species]);
    hash.real(candidate.temperature_switch_[species]);
    for (std::size_t coefficient = 0U; coefficient < 7U; ++coefficient) {
      hash.real(candidate.nasa_low_[coefficient][species]);
      hash.real(candidate.nasa_high_[coefficient][species]);
    }
  }
  for (const std::uint16_t mapped : candidate.independent_to_species_) {
    hash.integer(mapped);
  }
  candidate.fingerprint_ = hash.finish();
  out = std::move(candidate);
  return {};
}

Status TransportPlan::evaluate(
    double temperature, Span<const double> independent_mass_fractions,
    MolecularTransportState& out) const noexcept {
  const std::size_t count = molecular_weight_.size();
  if (!finite_positive(temperature) || temperature < minimum_temperature_ ||
      temperature > maximum_temperature_) {
    return {StatusCode::numerical_failure, kTransportTemperature};
  }
  if (count == 0U || count > kMaximumSpecies ||
      viscosity_reference_.size() != count ||
      reference_temperature_.size() != count ||
      sutherland_temperature_.size() != count || prandtl_.size() != count ||
      conductivity_.size() != count || temperature_switch_.size() != count ||
      molecular_weight_ratio_quarter_.size() != count * count ||
      wilke_denominator_reciprocal_.size() != count * count ||
      constant_wilke_phi_.size() !=
          (kernel_ == TransportKernel::constant ? count * count : 0U) ||
      independent_to_species_.size() + 1U != count ||
      dependent_species_ >= count) {
    return {StatusCode::invalid_plan, kTransportInput};
  }

  std::array<double, kMaximumSpecies> mass_fraction{};
  if (!composition(independent_mass_fractions,
                   {independent_to_species_.data(),
                    independent_to_species_.size()},
                   dependent_species_, count, mass_fraction)) {
    return {StatusCode::numerical_failure, kTransportComposition};
  }

  if (kernel_ == TransportKernel::coast_native_air) {
    const double viscosity = coast_native_air_viscosity(temperature);
    const double cp = species_cp(
        nasa_low_, nasa_high_,
        {temperature_switch_.data(), temperature_switch_.size()},
        {molecular_weight_.data(), molecular_weight_.size()}, 0U,
        temperature, kCoastNativeAirUniversalGasConstant);
    const double conductivity =
        cp * viscosity / kCoastNativeAirPrandtl;
    if (!finite_positive(viscosity) || !finite_positive(cp) ||
        !finite_positive(conductivity)) {
      return {StatusCode::numerical_failure, kTransportNumerical};
    }
    out = {viscosity, conductivity};
    return {};
  }

  std::array<double, kMaximumSpecies> mole_fraction{};
  std::array<double, kMaximumSpecies> viscosity{};
  std::array<double, kMaximumSpecies> conductivity{};
  double mole_sum = 0.0;
  for (std::size_t species = 0U; species < count; ++species) {
    const double moles = mass_fraction[species] / molecular_weight_[species];
    const double mu = species_viscosity(
        {viscosity_reference_.data(), viscosity_reference_.size()},
        {reference_temperature_.data(), reference_temperature_.size()},
        {sutherland_temperature_.data(), sutherland_temperature_.size()},
        species, temperature);
    const double lambda = species_conductivity(
        nasa_low_, nasa_high_,
        {temperature_switch_.data(), temperature_switch_.size()},
        {molecular_weight_.data(), molecular_weight_.size()},
        {reference_temperature_.data(), reference_temperature_.size()},
        {prandtl_.data(), prandtl_.size()},
        {conductivity_.data(), conductivity_.size()}, species, temperature,
        mu);
    if (!std::isfinite(moles) || moles < 0.0 || !finite_positive(mu) ||
        !finite_positive(lambda)) {
      return {StatusCode::numerical_failure, kTransportNumerical};
    }
    mole_fraction[species] = moles;
    viscosity[species] = mu;
    conductivity[species] = lambda;
    mole_sum += moles;
  }
  if (!finite_positive(mole_sum)) {
    return {StatusCode::numerical_failure, kTransportComposition};
  }
  for (std::size_t species = 0U; species < count; ++species) {
    mole_fraction[species] /= mole_sum;
  }

  double mixture_viscosity = 0.0;
  double mixture_conductivity = 0.0;
  if (kernel_ == TransportKernel::constant) {
    for (std::size_t first = 0U; first < count; ++first) {
      double wilke_denominator = 0.0;
      for (std::size_t second = 0U; second < count; ++second) {
        wilke_denominator +=
            mole_fraction[second] *
            constant_wilke_phi_[first * count + second];
      }
      if (!finite_positive(wilke_denominator)) {
        return {StatusCode::numerical_failure, kTransportNumerical};
      }
      mixture_viscosity +=
          mole_fraction[first] * viscosity[first] / wilke_denominator;
      mixture_conductivity +=
          mole_fraction[first] * conductivity[first] / wilke_denominator;
    }
  } else {
    for (std::size_t first = 0U; first < count; ++first) {
      double wilke_denominator = 0.0;
      for (std::size_t second = 0U; second < count; ++second) {
        const std::size_t interaction = first * count + second;
        const double numerator =
            1.0 + std::sqrt(viscosity[first] / viscosity[second]) *
                      molecular_weight_ratio_quarter_[interaction];
        const double phi = numerator * numerator *
                           wilke_denominator_reciprocal_[interaction];
        wilke_denominator += mole_fraction[second] * phi;
      }
      if (!finite_positive(wilke_denominator)) {
        return {StatusCode::numerical_failure, kTransportNumerical};
      }
      mixture_viscosity +=
          mole_fraction[first] * viscosity[first] / wilke_denominator;
      // Thermal conductivity uses the same finite Wilke interaction weights.
      // This keeps both molecular coefficients smooth at vanishing species
      // fractions without adding a second O(N^2) interaction table.
      mixture_conductivity +=
          mole_fraction[first] * conductivity[first] / wilke_denominator;
    }
  }
  if (!finite_positive(mixture_viscosity) ||
      !finite_positive(mixture_conductivity)) {
    return {StatusCode::numerical_failure, kTransportNumerical};
  }
  out = {mixture_viscosity, mixture_conductivity};
  return {};
}

Status TransportPlan::effective_enthalpy_transport(
    double molecular_viscosity, double effective_viscosity,
    double heat_capacity, double& conductivity,
    double& enthalpy_diffusivity) const noexcept {
  if (kernel_ != TransportKernel::coast_native_air ||
      enthalpy_prandtl_ != kCoastNativeAirPrandtl) {
    return {StatusCode::invalid_plan, kTransportInput};
  }
  if (!finite_positive(molecular_viscosity) ||
      !finite_positive(effective_viscosity) ||
      effective_viscosity < molecular_viscosity ||
      !finite_positive(heat_capacity)) {
    return {StatusCode::numerical_failure, kTransportNumerical};
  }
  const double candidate_diffusivity =
      effective_viscosity / enthalpy_prandtl_;
  const double candidate_conductivity =
      heat_capacity * candidate_diffusivity;
  if (!finite_positive(candidate_diffusivity) ||
      !finite_positive(candidate_conductivity)) {
    return {StatusCode::numerical_failure, kTransportNumerical};
  }
  conductivity = candidate_conductivity;
  enthalpy_diffusivity = candidate_diffusivity;
  return {};
}

}  // namespace hundun::v04

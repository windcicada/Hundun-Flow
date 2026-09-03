// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "hundun/v04_physics.hpp"

#include "field_view_interval_detail.hpp"
#include "physics_input_detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <string_view>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kThermoInput = 801U;
constexpr std::uint32_t kThermoComposition = 802U;
constexpr std::uint32_t kThermoRange = 803U;
constexpr std::uint32_t kThermoInversion = 804U;
constexpr std::uint32_t kThermoSpeciesBounds = 805U;
constexpr std::uint32_t kClosedMassInput = 821U;
constexpr std::uint32_t kClosedMassCollective = 822U;
constexpr std::uint32_t kClosedMassConvergence = 823U;
constexpr std::size_t kMaximumSpecies = 65U;
constexpr std::uint64_t kFnvOffset = 14695981039346656037ULL;
constexpr std::uint64_t kFnvPrime = 1099511628211ULL;
constexpr double kCoastNativeAirEnthalpyReferenceTemperature = 273.15;

bool finite(double value) noexcept { return std::isfinite(value); }

bool finite(Real3 value) noexcept {
  return finite(value.x) && finite(value.y) && finite(value.z);
}

bool same_tuple(ThermalRevisionTuple left,
                ThermalRevisionTuple right) noexcept {
  return left.enthalpy == right.enthalpy &&
         left.composition == right.composition &&
         left.enthalpy_storage == right.enthalpy_storage &&
         left.composition_storage == right.composition_storage &&
         left.enthalpy_revision_domain == right.enthalpy_revision_domain &&
         left.composition_revision_domain ==
             right.composition_revision_domain &&
         left.patch_identity == right.patch_identity;
}

bool valid_tuple(ThermalRevisionTuple value) noexcept {
  return value.enthalpy != 0U && value.composition != 0U &&
         value.enthalpy_storage != 0U &&
         value.composition_storage != 0U &&
         value.enthalpy_revision_domain != 0U &&
         value.composition_revision_domain != 0U &&
         value.patch_identity != 0U;
}

bool finite_coefficients(const std::array<double, 7U>& values) noexcept {
  for (const double value : values) {
    if (!finite(value)) {
      return false;
    }
  }
  return true;
}

class CompensatedSum {
 public:
  void add(double value) noexcept {
    const double tentative = sum_ + value;
    if (std::abs(sum_) >= std::abs(value)) {
      correction_ += (sum_ - tentative) + value;
    } else {
      correction_ += (value - tentative) + sum_;
    }
    sum_ = tentative;
  }

  double value() const noexcept { return sum_ + correction_; }

 private:
  double sum_{};
  double correction_{};
};

class Hash64 {
 public:
  void byte(std::uint8_t value) noexcept {
    value_ ^= value;
    value_ *= kFnvPrime;
  }

  template <class Integer>
  void integer(Integer value) noexcept {
    using Unsigned = std::make_unsigned_t<Integer>;
    Unsigned bits = static_cast<Unsigned>(value);
    for (std::size_t index = 0U; index < sizeof(bits); ++index) {
      byte(static_cast<std::uint8_t>((bits >> (8U * index)) & 0xffU));
    }
  }

  void real(double value) noexcept {
    if (value == 0.0) {
      value = 0.0;
    }
    std::uint64_t bits = 0U;
    std::memcpy(&bits, &value, sizeof(bits));
    integer(bits);
  }

  void text(std::string_view value) noexcept {
    integer(value.size());
    for (const char character : value) {
      byte(static_cast<std::uint8_t>(static_cast<unsigned char>(character)));
    }
  }

  PlanFingerprint finish() const noexcept {
    return value_ == 0U ? 1U : value_;
  }

 private:
  std::uint64_t value_{kFnvOffset};
};

bool valid_spec_header(const ThermophysicalSpec& spec) noexcept {
  return !spec.data_file.empty() && finite(spec.minimum_temperature) &&
         finite(spec.maximum_temperature) &&
         spec.minimum_temperature > 0.0 &&
         spec.maximum_temperature > spec.minimum_temperature &&
         finite(spec.temperature_relative_tolerance) &&
         spec.temperature_relative_tolerance > 0.0 &&
         spec.temperature_relative_tolerance < 1.0 &&
         spec.maximum_temperature_iterations > 0U &&
         spec.species.size() >= 1U &&
         spec.species.size() <= kMaximumSpecies;
}

double species_cp(double universal_gas_constant, double temperature,
                  double inverse_molecular_weight,
                  const std::array<double, 7U>& coefficients) noexcept {
  const double polynomial =
      ((((coefficients[4U] * temperature + coefficients[3U]) * temperature +
          coefficients[2U]) *
             temperature +
         coefficients[1U]) *
            temperature +
        coefficients[0U]);
  return universal_gas_constant * inverse_molecular_weight * polynomial;
}

double dimensionless_h(
    double temperature,
    const std::array<double, 7U>& coefficients) noexcept {
  const double t2 = temperature * temperature;
  const double t3 = t2 * temperature;
  const double t4 = t3 * temperature;
  return coefficients[0U] * temperature +
         coefficients[1U] * t2 * 0.5 + coefficients[2U] * t3 / 3.0 +
         coefficients[3U] * t4 * 0.25 +
         coefficients[4U] * t4 * temperature * 0.2 + coefficients[5U];
}

double species_h(double universal_gas_constant, double temperature,
                 double inverse_molecular_weight,
                 const std::array<double, 7U>& coefficients) noexcept {
  return universal_gas_constant * inverse_molecular_weight *
         dimensionless_h(temperature, coefficients);
}

Status collective_failure(MPI_Comm communicator, Status local,
                          int& lowest) noexcept {
  int rank = 0;
  int size = 0;
  if (communicator == MPI_COMM_NULL ||
      MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    lowest = -1;
    return {StatusCode::mpi_failure, kClosedMassCollective};
  }
  const std::uint64_t candidate =
      local.code == StatusCode::ok
          ? std::numeric_limits<std::uint64_t>::max()
          : static_cast<std::uint64_t>(rank);
  std::uint64_t selected = std::numeric_limits<std::uint64_t>::max();
  if (MPI_Allreduce(&candidate, &selected, 1, MPI_UINT64_T, MPI_MIN,
                    communicator) != MPI_SUCCESS) {
    lowest = -1;
    return {StatusCode::mpi_failure, kClosedMassCollective};
  }
  if (selected == std::numeric_limits<std::uint64_t>::max()) {
    lowest = -1;
    return {};
  }
  if (selected >= static_cast<std::uint64_t>(size)) {
    lowest = -1;
    return {StatusCode::mpi_failure, kClosedMassCollective};
  }
  const int failing_rank = static_cast<int>(selected);
  std::array<std::uint64_t, 2U> wire{};
  if (rank == failing_rank) {
    wire[0U] = static_cast<std::uint64_t>(local.code);
    wire[1U] = local.detail;
  }
  if (MPI_Bcast(wire.data(), static_cast<int>(wire.size()), MPI_UINT64_T,
                failing_rank, communicator) != MPI_SUCCESS ||
      wire[0U] > static_cast<std::uint64_t>(StatusCode::io_failure) ||
      wire[1U] > std::numeric_limits<std::uint32_t>::max()) {
    lowest = -1;
    return {StatusCode::mpi_failure, kClosedMassCollective};
  }
  lowest = failing_rank;
  return {static_cast<StatusCode>(wire[0U]),
          static_cast<std::uint32_t>(wire[1U])};
}

}  // namespace

Status ThermodynamicsPlan::compile(
    const ThermophysicalSpec& spec,
    Span<const TransportedScalarSpec> scalar_catalog,
    ThermodynamicsPlan& out) noexcept {
  if (!valid_spec_header(spec) ||
      !detail::valid_thermophysical_spec(spec) ||
      (scalar_catalog.size != 0U && scalar_catalog.data == nullptr)) {
    return {StatusCode::invalid_plan, kThermoInput};
  }
  std::size_t independent_count = 0U;
  for (std::size_t index = 0U; index < scalar_catalog.size; ++index) {
    if (scalar_catalog.data[index].role == TransportedScalarRole::species) {
      ++independent_count;
    }
  }
  if (spec.species.size() != independent_count + 1U) {
    return {StatusCode::invalid_plan, kThermoInput};
  }
  for (std::size_t species = 0U; species < spec.species.size(); ++species) {
    const SpeciesThermophysicalSpec& value = spec.species[species];
    if (value.stable_name.empty() || !finite(value.molecular_weight) ||
        value.molecular_weight <= 0.0 ||
        !finite(value.temperature_switch) ||
        value.temperature_switch < spec.minimum_temperature ||
        value.temperature_switch > spec.maximum_temperature ||
        !finite_coefficients(value.nasa7_low) ||
        !finite_coefficients(value.nasa7_high)) {
      return {StatusCode::invalid_plan, kThermoInput};
    }
    for (std::size_t prior = 0U; prior < species; ++prior) {
      if (spec.species[prior].stable_name == value.stable_name) {
        return {StatusCode::invalid_plan, kThermoInput};
      }
    }
  }
  std::array<std::uint8_t, kMaximumSpecies> matched{};
  for (std::size_t catalog = 0U; catalog < scalar_catalog.size; ++catalog) {
    if (scalar_catalog.data[catalog].role != TransportedScalarRole::species) {
      continue;
    }
    std::size_t match = spec.species.size();
    for (std::size_t species = 0U; species < spec.species.size(); ++species) {
      if (spec.species[species].stable_name ==
          scalar_catalog.data[catalog].stable_name) {
        match = species;
        break;
      }
    }
    if (match == spec.species.size() || matched[match] != 0U) {
      return {StatusCode::invalid_plan, kThermoInput};
    }
    matched[match] = 1U;
  }
  std::size_t dependent_species = spec.species.size();
  for (std::size_t species = 0U; species < spec.species.size(); ++species) {
    if (matched[species] == 0U) {
      if (dependent_species != spec.species.size()) {
        return {StatusCode::invalid_plan, kThermoInput};
      }
      dependent_species = species;
    }
  }
  if (dependent_species == spec.species.size()) {
    return {StatusCode::invalid_plan, kThermoInput};
  }

  try {
    ThermophysicalSpec canonical_spec = spec;
    if (!detail::canonicalize_thermophysical_spec(canonical_spec)) {
      return {StatusCode::invalid_plan, kThermoInput};
    }
    ThermodynamicsPlan candidate;
    const std::size_t count = canonical_spec.species.size();
    candidate.inverse_molecular_weight_.resize(count);
    candidate.temperature_switch_.resize(count);
    candidate.species_enthalpy_minimum_.resize(count);
    candidate.species_enthalpy_maximum_.resize(count);
    for (std::vector<double>& coefficients : candidate.nasa_low_) {
      coefficients.resize(count);
    }
    for (std::vector<double>& coefficients : candidate.nasa_high_) {
      coefficients.resize(count);
    }
    candidate.independent_to_species_.reserve(independent_count);
    for (std::size_t catalog = 0U; catalog < scalar_catalog.size; ++catalog) {
      if (scalar_catalog.data[catalog].role == TransportedScalarRole::species) {
        for (std::size_t species = 0U; species < count; ++species) {
          if (spec.species[species].stable_name ==
              scalar_catalog.data[catalog].stable_name) {
            candidate.independent_to_species_.push_back(
                static_cast<std::uint16_t>(species));
            break;
          }
        }
      }
    }
    candidate.dependent_species_ = dependent_species;
    candidate.minimum_temperature_ = canonical_spec.minimum_temperature;
    candidate.maximum_temperature_ = canonical_spec.maximum_temperature;
    candidate.relative_tolerance_ =
        canonical_spec.temperature_relative_tolerance;
    candidate.maximum_iterations_ =
        canonical_spec.maximum_temperature_iterations;
    const bool coast_native_air =
        canonical_spec.species.size() == 1U &&
        canonical_spec.species.front().transport_law ==
            TransportLaw::coast_native_air;
    candidate.universal_gas_constant_ =
        coast_native_air ? kCoastNativeAirUniversalGasConstant
                         : kUniversalGasConstant;
    candidate.source_fingerprint_ =
        detail::thermophysical_spec_fingerprint(canonical_spec);
    if (candidate.source_fingerprint_ == 0U) {
      return {StatusCode::invalid_plan, kThermoInput};
    }

    bool constant_cp = true;
    bool common_enthalpy_offset = true;
    double enthalpy_offset = 0.0;
    Hash64 hash;
    hash.text("HUNDUN-FLOW-v0.4-thermodynamics-v1");
    hash.integer(count);
    hash.real(candidate.minimum_temperature_);
    hash.real(candidate.maximum_temperature_);
    hash.real(candidate.relative_tolerance_);
    hash.integer(candidate.maximum_iterations_);
    hash.real(candidate.universal_gas_constant_);
    hash.integer(static_cast<std::uint16_t>(candidate.dependent_species_));
    for (const std::uint16_t mapped : candidate.independent_to_species_) {
      hash.integer(mapped);
    }
    for (std::size_t species = 0U; species < count; ++species) {
      const SpeciesThermophysicalSpec& source =
          canonical_spec.species[species];
      std::array<double, 7U> low = source.nasa7_low;
      std::array<double, 7U> high = source.nasa7_high;
      if (coast_native_air) {
        // Ordinary non-reacting COAST stores primary h as h(T)-h(273.15 K).
        const double reference = dimensionless_h(
            kCoastNativeAirEnthalpyReferenceTemperature, low);
        low[5U] -= reference;
        high[5U] -= reference;
      }
      candidate.inverse_molecular_weight_[species] =
          1.0 / source.molecular_weight;
      candidate.temperature_switch_[species] = source.temperature_switch;
      hash.text(source.stable_name);
      hash.real(source.molecular_weight);
      hash.real(source.temperature_switch);
      for (std::size_t coefficient = 0U; coefficient < 7U; ++coefficient) {
        candidate.nasa_low_[coefficient][species] = low[coefficient];
        candidate.nasa_high_[coefficient][species] = high[coefficient];
        hash.real(low[coefficient]);
        hash.real(high[coefficient]);
      }
      candidate.species_enthalpy_minimum_[species] = species_h(
          candidate.universal_gas_constant_,
          candidate.minimum_temperature_,
          candidate.inverse_molecular_weight_[species],
          candidate.minimum_temperature_ <= candidate.temperature_switch_[species]
              ? low
              : high);
      candidate.species_enthalpy_maximum_[species] = species_h(
          candidate.universal_gas_constant_,
          candidate.maximum_temperature_,
          candidate.inverse_molecular_weight_[species],
          candidate.maximum_temperature_ <= candidate.temperature_switch_[species]
              ? low
              : high);
      if (!finite(candidate.species_enthalpy_minimum_[species]) ||
          !finite(candidate.species_enthalpy_maximum_[species])) {
        return {StatusCode::invalid_plan, kThermoSpeciesBounds};
      }
      for (std::size_t coefficient = 1U; coefficient <= 4U; ++coefficient) {
        constant_cp &= low[coefficient] == 0.0 &&
                       high[coefficient] == 0.0;
      }
      constant_cp &= low[0U] == high[0U];
      constant_cp &= low[5U] == high[5U];
      const double source_offset = candidate.universal_gas_constant_ *
                                   candidate.inverse_molecular_weight_[species] *
                                   low[5U];
      if (species == 0U) {
        enthalpy_offset = source_offset;
      } else {
        common_enthalpy_offset &= source_offset == enthalpy_offset;
      }
    }
    constant_cp &= common_enthalpy_offset;
    candidate.kernel_ = constant_cp ? ThermodynamicsKernel::constant_cp
                                    : ThermodynamicsKernel::nasa7;
    candidate.constant_enthalpy_offset_ =
        constant_cp ? enthalpy_offset : 0.0;
    hash.integer(static_cast<std::uint8_t>(candidate.kernel_));
    hash.integer(candidate.source_fingerprint_);
    candidate.fingerprint_ = hash.finish();

    // The shared cold validator above proves cp>R over each complete NASA7
    // interval and enforces cp/h continuity.  Do not duplicate it with a
    // weaker sample-only pass here.
    out = std::move(candidate);
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, kThermoInput};
  } catch (...) {
    return {StatusCode::invalid_plan, kThermoInput};
  }
}

Status ThermodynamicsPlan::composition(
    Span<const double> independent_mass_fractions,
    double& dependent) const noexcept {
  if (fingerprint_ == 0U ||
      independent_mass_fractions.size != independent_to_species_.size() ||
      (independent_mass_fractions.size != 0U &&
       independent_mass_fractions.data == nullptr)) {
    return {StatusCode::invalid_plan, kThermoComposition};
  }
  double sum = 0.0;
  for (std::size_t index = 0U; index < independent_mass_fractions.size;
       ++index) {
    const double value = independent_mass_fractions.data[index];
    if (!finite(value) || value < 0.0 || value > 1.0) {
      return {StatusCode::numerical_failure, kThermoComposition};
    }
    sum += value;
    if (!finite(sum) || sum > 1.0) {
      return {StatusCode::numerical_failure, kThermoComposition};
    }
  }
  dependent = 1.0 - sum;
  return finite(dependent) && dependent >= 0.0
             ? Status{}
             : Status{StatusCode::numerical_failure, kThermoComposition};
}

Status ThermodynamicsPlan::mixture_properties(
    double temperature, Span<const double> independent_mass_fractions,
    double dependent, double& enthalpy, double& cp,
    double& gas_constant) const noexcept {
  if (!finite(temperature) || temperature < minimum_temperature_ ||
      temperature > maximum_temperature_ || !finite(dependent) ||
      dependent < 0.0 || inverse_molecular_weight_.empty() ||
      dependent_species_ >= inverse_molecular_weight_.size()) {
    return {StatusCode::numerical_failure, kThermoRange};
  }
  double mixture_h = 0.0;
  double mixture_cp = 0.0;
  double mixture_r = 0.0;
  for (std::size_t independent = 0U;
       independent < independent_mass_fractions.size; ++independent) {
    const std::size_t species = independent_to_species_[independent];
    const double mass_fraction = independent_mass_fractions.data[independent];
    const bool low = temperature <= temperature_switch_[species];
    std::array<double, 7U> coefficients{};
    for (std::size_t index = 0U; index < 7U; ++index) {
      coefficients[index] =
          (low ? nasa_low_[index] : nasa_high_[index])[species];
    }
    mixture_h += mass_fraction *
                 species_h(universal_gas_constant_, temperature,
                           inverse_molecular_weight_[species],
                           coefficients);
    mixture_cp += mass_fraction *
                  species_cp(universal_gas_constant_, temperature,
                             inverse_molecular_weight_[species],
                             coefficients);
    mixture_r += mass_fraction * universal_gas_constant_ *
                 inverse_molecular_weight_[species];
  }
  {
    const std::size_t species = dependent_species_;
    const bool low = temperature <= temperature_switch_[species];
    std::array<double, 7U> coefficients{};
    for (std::size_t index = 0U; index < 7U; ++index) {
      coefficients[index] =
          (low ? nasa_low_[index] : nasa_high_[index])[species];
    }
    mixture_h += dependent *
                 species_h(universal_gas_constant_, temperature,
                           inverse_molecular_weight_[species],
                           coefficients);
    mixture_cp += dependent *
                  species_cp(universal_gas_constant_, temperature,
                             inverse_molecular_weight_[species],
                             coefficients);
    mixture_r += dependent * universal_gas_constant_ *
                 inverse_molecular_weight_[species];
  }
  if (!finite(mixture_h) || !finite(mixture_cp) || !finite(mixture_r) ||
      mixture_r <= 0.0 || mixture_cp <= mixture_r) {
    return {StatusCode::numerical_failure, kThermoRange};
  }
  enthalpy = mixture_h;
  cp = mixture_cp;
  gas_constant = mixture_r;
  return {};
}

Status ThermodynamicsPlan::mixture_enthalpy(
    double temperature, Span<const double> independent_mass_fractions,
    double& enthalpy, double& cp, double& gas_constant) const noexcept {
  double dependent = 0.0;
  const Status valid = composition(independent_mass_fractions, dependent);
  if (!valid) {
    return valid;
  }
  double candidate_h = 0.0;
  double candidate_cp = 0.0;
  double candidate_r = 0.0;
  const Status evaluated = mixture_properties(
      temperature, independent_mass_fractions, dependent, candidate_h,
      candidate_cp, candidate_r);
  if (!evaluated) {
    return evaluated;
  }
  enthalpy = candidate_h;
  cp = candidate_cp;
  gas_constant = candidate_r;
  return {};
}

Status ThermodynamicsPlan::species_enthalpy_bounds(
    std::size_t full_species_index, double& at_minimum,
    double& at_maximum) const noexcept {
  if (fingerprint_ == 0U ||
      full_species_index >= species_enthalpy_minimum_.size() ||
      species_enthalpy_minimum_.size() !=
          species_enthalpy_maximum_.size()) {
    return {StatusCode::invalid_plan, kThermoSpeciesBounds};
  }
  const double minimum = species_enthalpy_minimum_[full_species_index];
  const double maximum = species_enthalpy_maximum_[full_species_index];
  if (!finite(minimum) || !finite(maximum)) {
    return {StatusCode::numerical_failure, kThermoSpeciesBounds};
  }
  at_minimum = minimum;
  at_maximum = maximum;
  return {};
}

Status ThermodynamicsPlan::independent_species_enthalpy_bounds(
    std::size_t independent_index, double& at_minimum,
    double& at_maximum) const noexcept {
  if (fingerprint_ == 0U ||
      independent_index >= independent_to_species_.size()) {
    return {StatusCode::invalid_plan, kThermoSpeciesBounds};
  }
  return species_enthalpy_bounds(independent_to_species_[independent_index],
                                 at_minimum, at_maximum);
}

Status ThermodynamicsPlan::evaluate(
    double p_abs, double h,
    Span<const double> independent_mass_fractions, Real3 velocity,
    ThermoState& out, double temperature_hint) const noexcept {
  ThermalState thermal;
  const ThermalRevisionTuple one_shot{};
  const Status evaluated = evaluate_thermal_impl(
      h, independent_mass_fractions, one_shot, false, thermal,
      temperature_hint);
  if (!evaluated) {
    return evaluated;
  }
  return complete_state_impl(p_abs, thermal, one_shot, false, velocity, out);
}

Status ThermodynamicsPlan::evaluate_thermal(
    double h, Span<const double> independent_mass_fractions,
    ThermalRevisionTuple revisions, ThermalState& out,
    double temperature_hint) const noexcept {
  return evaluate_thermal_impl(h, independent_mass_fractions, revisions, true,
                               out, temperature_hint);
}

Status ThermodynamicsPlan::evaluate_thermal_impl(
    double h, Span<const double> independent_mass_fractions,
    ThermalRevisionTuple revisions, bool require_certificate,
    ThermalState& out, double temperature_hint) const noexcept {
  if (!finite(h) || (require_certificate && !valid_tuple(revisions))) {
    return {StatusCode::numerical_failure, kThermoInput};
  }
  double dependent = 0.0;
  const Status valid = composition(independent_mass_fractions, dependent);
  if (!valid) {
    return valid;
  }

  if (kernel_ == ThermodynamicsKernel::constant_cp) {
    double cp = 0.0;
    double gas = 0.0;
    for (std::size_t independent = 0U;
         independent < independent_mass_fractions.size; ++independent) {
      const std::size_t species = independent_to_species_[independent];
      const double fraction = independent_mass_fractions.data[independent];
      const double species_gas =
          universal_gas_constant_ * inverse_molecular_weight_[species];
      gas += fraction * species_gas;
      cp += fraction * species_gas * nasa_low_[0U][species];
    }
    const double dependent_gas =
        universal_gas_constant_ *
        inverse_molecular_weight_[dependent_species_];
    gas += dependent * dependent_gas;
    cp += dependent * dependent_gas *
          nasa_low_[0U][dependent_species_];
    if (!finite(cp) || !finite(gas) || gas <= 0.0 || cp <= gas) {
      return {StatusCode::numerical_failure, kThermoRange};
    }
    const double temperature = (h - constant_enthalpy_offset_) / cp;
    if (!finite(temperature) || temperature < minimum_temperature_ ||
        temperature > maximum_temperature_) {
      return {StatusCode::numerical_failure, kThermoInversion};
    }
    const double gamma = cp / (cp - gas);
    if (!finite(gamma) || gamma <= 1.0) {
      return {StatusCode::numerical_failure, kThermoRange};
    }
    const double psi = 1.0 / (gas * temperature);
    const double dpsi_dh = -psi / (temperature * cp);
    const double sound_squared = gamma * gas * temperature;
    if (!finite(psi) || psi <= 0.0 || !finite(dpsi_dh) ||
        dpsi_dh >= 0.0 || !finite(sound_squared) ||
        sound_squared <= 0.0) {
      return {StatusCode::numerical_failure, kThermoRange};
    }
    ThermalState candidate;
    candidate.temperature_ = temperature;
    candidate.cp_ = cp;
    candidate.gas_constant_ = gas;
    candidate.gamma_ = gamma;
    candidate.drho_dp_hY_ = psi;
    candidate.dpressure_compressibility_dh_hY_ = dpsi_dh;
    candidate.sound_speed_ = std::sqrt(sound_squared);
    candidate.thermodynamics_ = fingerprint_;
    candidate.revisions_ = revisions;
    candidate.certified_ = true;
    out = candidate;
    return {};
  }

  double low_h = 0.0;
  double low_cp = 0.0;
  double low_r = 0.0;
  double high_h = 0.0;
  double high_cp = 0.0;
  double high_r = 0.0;
  Status status = mixture_properties(
      minimum_temperature_, independent_mass_fractions, dependent, low_h,
      low_cp, low_r);
  if (status) {
    status = mixture_properties(maximum_temperature_,
                                independent_mass_fractions, dependent, high_h,
                                high_cp, high_r);
  }
  if (!status || h < low_h || h > high_h) {
    return {StatusCode::numerical_failure, kThermoInversion};
  }

  double lower = minimum_temperature_;
  double upper = maximum_temperature_;
  double temperature = 0.0;
  double evaluated_h = 0.0;
  double cp = 0.0;
  double gas = 0.0;
  bool converged = false;
  if (h == low_h) {
    temperature = lower;
    evaluated_h = low_h;
    cp = low_cp;
    gas = low_r;
    converged = true;
  } else if (h == high_h) {
    temperature = upper;
    evaluated_h = high_h;
    cp = high_cp;
    gas = high_r;
    converged = true;
  } else {
    temperature =
        finite(temperature_hint) && temperature_hint > lower &&
                temperature_hint < upper
            ? temperature_hint
            : lower + (upper - lower) * (h - low_h) / (high_h - low_h);
    if (!finite(temperature) || temperature <= lower ||
        temperature >= upper) {
      temperature = 0.5 * (lower + upper);
    }
    for (std::uint32_t iteration = 0U; iteration < maximum_iterations_;
         ++iteration) {
      status = mixture_properties(temperature, independent_mass_fractions,
                                  dependent, evaluated_h, cp, gas);
      if (!status) {
        return status;
      }
      const double residual = evaluated_h - h;
      const double scale =
          std::max({1.0, std::abs(h), std::abs(evaluated_h)});
      if (std::abs(residual) <= relative_tolerance_ * scale) {
        converged = true;
        break;
      }
      if (residual > 0.0) {
        upper = temperature;
      } else {
        lower = temperature;
      }
      double next = temperature - residual / cp;
      if (!finite(next) || next <= lower || next >= upper) {
        next = 0.5 * (lower + upper);
      }
      temperature = next;
    }
  }
  if (!converged) {
    return {StatusCode::numerical_failure, kThermoInversion};
  }

  const double gamma = cp / (cp - gas);
  if (!finite(gamma) || gamma <= 1.0) {
    return {StatusCode::numerical_failure, kThermoRange};
  }
  const double psi = 1.0 / (gas * temperature);
  const double dpsi_dh = -psi / (temperature * cp);
  const double sound_squared = gamma * gas * temperature;
  if (!finite(psi) || psi <= 0.0 || !finite(dpsi_dh) ||
      dpsi_dh >= 0.0 || !finite(sound_squared) ||
      sound_squared <= 0.0) {
    return {StatusCode::numerical_failure, kThermoRange};
  }
  ThermalState candidate;
  candidate.temperature_ = temperature;
  candidate.cp_ = cp;
  candidate.gas_constant_ = gas;
  candidate.gamma_ = gamma;
  candidate.drho_dp_hY_ = psi;
  candidate.dpressure_compressibility_dh_hY_ = dpsi_dh;
  candidate.sound_speed_ = std::sqrt(sound_squared);
  candidate.thermodynamics_ = fingerprint_;
  candidate.revisions_ = revisions;
  candidate.certified_ = true;
  out = candidate;
  return {};
}

Status ThermodynamicsPlan::evaluate_pressure(
    double p_abs, const ThermalState& thermal,
    ThermalRevisionTuple revisions, PressureThermoState& out) const noexcept {
  return evaluate_pressure_impl(p_abs, thermal, revisions, true, out);
}

Status ThermodynamicsPlan::evaluate_pressure_impl(
    double p_abs, const ThermalState& thermal,
    ThermalRevisionTuple revisions, bool require_certificate,
    PressureThermoState& out) const noexcept {
  if (!finite(p_abs) || p_abs <= 0.0 ||
      !thermal.certified_ || thermal.thermodynamics_ != fingerprint_ ||
      (require_certificate &&
       (!valid_tuple(revisions) ||
        !same_tuple(thermal.revisions_, revisions))) ||
      !finite(thermal.drho_dp_hY_) || thermal.drho_dp_hY_ <= 0.0 ||
      !finite(thermal.dpressure_compressibility_dh_hY_) ||
      thermal.dpressure_compressibility_dh_hY_ >= 0.0 ||
      !finite(thermal.temperature_) ||
      thermal.temperature_ < minimum_temperature_ ||
      thermal.temperature_ > maximum_temperature_) {
    return {StatusCode::numerical_failure, kThermoRange};
  }
  const double density = p_abs * thermal.drho_dp_hY_;
  const double drho_dh =
      p_abs * thermal.dpressure_compressibility_dh_hY_;
  if (!finite(density) || density <= 0.0 || !finite(drho_dh) ||
      drho_dh >= 0.0) {
    return {StatusCode::numerical_failure, kThermoRange};
  }
  out = {density, thermal.drho_dp_hY_, drho_dh};
  return {};
}

Status ThermodynamicsPlan::complete_state(double p_abs,
                                          const ThermalState& thermal,
                                          ThermalRevisionTuple revisions,
                                          Real3 velocity,
                                          ThermoState& out) const noexcept {
  return complete_state_impl(p_abs, thermal, revisions, true, velocity, out);
}

Status ThermodynamicsPlan::complete_state_impl(
    double p_abs, const ThermalState& thermal,
    ThermalRevisionTuple revisions, bool require_certificate, Real3 velocity,
    ThermoState& out) const noexcept {
  if (!finite(p_abs) || p_abs <= 0.0 || !finite(velocity) ||
      !thermal.certified_ || thermal.thermodynamics_ != fingerprint_ ||
      (require_certificate &&
       (!valid_tuple(revisions) ||
        !same_tuple(thermal.revisions_, revisions))) ||
      !finite(thermal.temperature_) ||
      thermal.temperature_ < minimum_temperature_ ||
      thermal.temperature_ > maximum_temperature_ || !finite(thermal.cp_) ||
      !finite(thermal.gas_constant_) || !finite(thermal.gamma_) ||
      thermal.gas_constant_ <= 0.0 ||
      thermal.cp_ <= thermal.gas_constant_ || thermal.gamma_ <= 1.0) {
    return {StatusCode::numerical_failure, kThermoRange};
  }
  PressureThermoState pressure;
  const Status pressure_status = evaluate_pressure_impl(
      p_abs, thermal, revisions, require_certificate, pressure);
  if (!pressure_status) {
    return pressure_status;
  }
  const double velocity_squared = velocity.x * velocity.x +
                                  velocity.y * velocity.y +
                                  velocity.z * velocity.z;
  if (!finite(thermal.sound_speed_) || thermal.sound_speed_ <= 0.0 ||
      !finite(velocity_squared) ||
      velocity_squared < 0.0) {
    return {StatusCode::numerical_failure, kThermoRange};
  }
  ThermoState candidate;
  candidate.rho = pressure.rho;
  candidate.temperature = thermal.temperature_;
  candidate.cp = thermal.cp_;
  candidate.gas_constant = thermal.gas_constant_;
  candidate.gamma = thermal.gamma_;
  candidate.drho_dp_hY = pressure.drho_dp_hY;
  candidate.drho_dh_pY = pressure.drho_dh_pY;
  candidate.sound_speed = thermal.sound_speed_;
  candidate.mach = std::sqrt(velocity_squared) / thermal.sound_speed_;
  if (!finite(candidate.drho_dp_hY) || candidate.drho_dp_hY <= 0.0 ||
      !finite(candidate.drho_dh_pY) || candidate.drho_dh_pY >= 0.0 ||
      !finite(candidate.mach)) {
    return {StatusCode::numerical_failure, kThermoRange};
  }
  out = candidate;
  return {};
}

Status ThermodynamicsPlan::evaluate_from_reference_pressure(
    double p_ref, double pi, double h,
    Span<const double> independent_mass_fractions, Real3 velocity,
    ThermoState& out, double temperature_hint) const noexcept {
  const double absolute_pressure = p_ref + pi;
  if (!finite(p_ref) || !finite(pi) || !finite(absolute_pressure) ||
      absolute_pressure <= 0.0) {
    return {StatusCode::numerical_failure, kThermoInput};
  }
  return evaluate(absolute_pressure, h, independent_mass_fractions, velocity,
                  out, temperature_hint);
}

Status ThermodynamicsPlan::evaluate_from_density(
    double density, double h,
    Span<const double> independent_mass_fractions, Real3 velocity,
    double& pressure_absolute, ThermoState& out,
    double temperature_hint) const noexcept {
  if (!finite(density) || density <= 0.0) {
    return {StatusCode::numerical_failure, kThermoInput};
  }
  ThermalState thermal;
  const ThermalRevisionTuple one_shot{};
  Status status = evaluate_thermal_impl(
      h, independent_mass_fractions, one_shot, false, thermal,
      temperature_hint);
  if (!status) return status;
  const double p_abs = density / thermal.drho_dp_hY_;
  if (!finite(p_abs) || p_abs <= 0.0) {
    return {StatusCode::numerical_failure, kThermoRange};
  }
  ThermoState candidate;
  status = complete_state_impl(p_abs, thermal, one_shot, false, velocity,
                               candidate);
  if (!status) return status;
  const double oracle_density = p_abs * candidate.drho_dp_hY;
  const double scale = std::max({1.0, density, std::abs(oracle_density)});
  if (!finite(oracle_density) || oracle_density <= 0.0 ||
      std::abs(oracle_density - density) >
          32.0 * std::numeric_limits<double>::epsilon() * scale) {
    return {StatusCode::numerical_failure, kThermoRange};
  }
  // rho* is the predictor authority.  Keep the independently formed EOS
  // product available to the caller while preventing a divide/multiply
  // roundoff from silently replacing the conservative state.
  candidate.rho = density;
  pressure_absolute = p_abs;
  out = candidate;
  return {};
}

Status ClosedMassPlan::compile(PressureReferenceKind authority,
                               const ThermophysicalSpec& spec,
                               ClosedMassPlan& out) noexcept {
  if (static_cast<std::uint8_t>(authority) >
          static_cast<std::uint8_t>(PressureReferenceKind::closed_mass) ||
      !finite(spec.closed_mass_relative_tolerance) ||
      spec.closed_mass_relative_tolerance <= 0.0 ||
      spec.closed_mass_relative_tolerance >= 1.0 ||
      !finite(spec.maximum_closed_mass_relative_step) ||
      spec.maximum_closed_mass_relative_step <= 0.0 ||
      spec.maximum_closed_mass_relative_step >= 1.0 ||
      spec.maximum_closed_mass_iterations == 0U ||
      spec.maximum_closed_mass_iterations > 256U) {
    return {StatusCode::invalid_plan, kClosedMassInput};
  }
  ClosedMassPlan candidate;
  candidate.authority_ = authority;
  candidate.relative_tolerance_ = spec.closed_mass_relative_tolerance;
  candidate.maximum_relative_step_ =
      spec.maximum_closed_mass_relative_step;
  candidate.maximum_iterations_ = spec.maximum_closed_mass_iterations;
  Hash64 hash;
  hash.text("HUNDUN-FLOW-v0.4-closed-mass-v1");
  hash.integer(static_cast<std::uint8_t>(authority));
  hash.real(candidate.relative_tolerance_);
  hash.real(candidate.maximum_relative_step_);
  hash.integer(candidate.maximum_iterations_);
  candidate.fingerprint_ = hash.finish();
  out = candidate;
  return {};
}

Status ClosedMassPlan::solve(MPI_Comm communicator,
                             const ThermodynamicsPlan& thermodynamics,
                             const ClosedMassCellView& cells,
                             double target_mass,
                             double current_pressure_reference,
                             ClosedMassResult& out) const noexcept {
  int size = 0;
  if (communicator == MPI_COMM_NULL ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return {StatusCode::mpi_failure, kClosedMassCollective};
  }
  const std::size_t count = cells.enthalpy.size;
  const std::size_t independent = thermodynamics.independent_species_count();
  const bool scalar_shape =
      cells.pressure_perturbation.size == count && cells.volume.size == count &&
      cells.active.size == count && count > 0U;
  const bool pointer_shape =
      cells.pressure_perturbation.data != nullptr &&
      cells.enthalpy.data != nullptr && cells.volume.data != nullptr &&
      cells.active.data != nullptr &&
      (independent == 0U || cells.independent_mass_fractions.data != nullptr);
  const bool composition_shape =
      independent == 0U
          ? cells.independent_mass_fractions.size == 0U
          : count <= std::numeric_limits<std::size_t>::max() / independent &&
                cells.independent_mass_fractions.size == count * independent;
  Status local{};
  if (authority_ != PressureReferenceKind::closed_mass ||
      thermodynamics.fingerprint() == 0U || !finite(target_mass) ||
      target_mass <= 0.0 || !finite(current_pressure_reference) ||
      current_pressure_reference <= 0.0 || !scalar_shape || !pointer_shape ||
      !composition_shape) {
    local = {StatusCode::invalid_plan, kClosedMassInput};
  }
  int lowest = -1;
  Status collective = collective_failure(communicator, local, lowest);
  if (!collective) {
    return collective;
  }
  std::uint64_t target_bits = 0U;
  std::uint64_t pressure_bits = 0U;
  std::memcpy(&target_bits, &target_mass, sizeof(target_bits));
  std::memcpy(&pressure_bits, &current_pressure_reference,
              sizeof(pressure_bits));
  const std::array<std::uint64_t, 4U> identity{
      thermodynamics.fingerprint(), fingerprint_, target_bits, pressure_bits};
  std::array<std::uint64_t, 4U> minimum_identity{};
  std::array<std::uint64_t, 4U> maximum_identity{};
  const int minimum_status = MPI_Allreduce(
      identity.data(), minimum_identity.data(),
      static_cast<int>(identity.size()), MPI_UINT64_T, MPI_MIN, communicator);
  const int maximum_status = MPI_Allreduce(
      identity.data(), maximum_identity.data(),
      static_cast<int>(identity.size()), MPI_UINT64_T, MPI_MAX, communicator);
  if (minimum_status != MPI_SUCCESS || maximum_status != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kClosedMassCollective};
  }
  if (minimum_identity != maximum_identity) {
    return {StatusCode::invalid_plan, kClosedMassCollective};
  }

  // At fixed h, Y, and pi the ideal-gas mass is affine in p_ref:
  //   M(p_ref) = p_ref * sum(V/(Rmix*T))
  //              + sum(pi*V/(Rmix*T)).
  // Invert h(T,Y) once per cell, then retain only these two local sums.  The
  // bounded Newton loop below is consequently scalar-only and performs no
  // repeated thermodynamic work or communication.
  CompensatedSum local_slope_sum;
  CompensatedSum local_intercept_sum;
  std::array<double, kMaximumSpecies - 1U> composition{};
  local = {};
  for (std::size_t cell = 0U; cell < count; ++cell) {
    if (cells.active.data[cell] > 1U) {
      local = {StatusCode::invalid_plan, kClosedMassInput};
      break;
    }
    if (cells.active.data[cell] == 0U) {
      continue;
    }
    const double volume = cells.volume.data[cell];
    const double enthalpy = cells.enthalpy.data[cell];
    const double pi = cells.pressure_perturbation.data[cell];
    const double initial_absolute_pressure =
        current_pressure_reference + pi;
    if (!finite(volume) || volume <= 0.0 || !finite(enthalpy) ||
        !finite(pi) || !finite(initial_absolute_pressure) ||
        initial_absolute_pressure <= 0.0) {
      local = {StatusCode::numerical_failure, kClosedMassInput};
      break;
    }
    for (std::size_t species = 0U; species < independent; ++species) {
      composition[species] =
          cells.independent_mass_fractions.data[species * count + cell];
    }
    ThermalState thermal;
    const Status evaluated = thermodynamics.evaluate_thermal_impl(
        enthalpy, Span<const double>{composition.data(), independent}, {},
        false, thermal, std::numeric_limits<double>::quiet_NaN());
    if (!evaluated) {
      local = evaluated;
      break;
    }
    const double slope = volume * thermal.drho_dp_hY_;
    const double intercept = pi * slope;
    if (!finite(slope) || slope <= 0.0 || !finite(intercept)) {
      local = {StatusCode::numerical_failure, kClosedMassInput};
      break;
    }
    local_slope_sum.add(slope);
    local_intercept_sum.add(intercept);
  }
  collective = collective_failure(communicator, local, lowest);
  if (!collective) {
    return collective;
  }
  const std::array<double, 2U> local_coefficients{
      local_slope_sum.value(), local_intercept_sum.value()};
  std::array<double, 2U> global_coefficients{};
  if (MPI_Allreduce(local_coefficients.data(), global_coefficients.data(),
                    static_cast<int>(global_coefficients.size()), MPI_DOUBLE,
                    MPI_SUM, communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kClosedMassCollective};
  }
  const double global_derivative = global_coefficients[0U];
  const double global_intercept = global_coefficients[1U];
  if (!finite(global_derivative) || global_derivative <= 0.0 ||
      !finite(global_intercept)) {
    return {StatusCode::numerical_failure, kClosedMassInput};
  }

  double pressure = current_pressure_reference;
  double global_mass = 0.0;
  std::uint32_t iterations = 0U;
  bool converged = false;
  for (std::uint32_t iteration = 0U; iteration < maximum_iterations_;
       ++iteration) {
    global_mass = pressure * global_derivative + global_intercept;
    const double residual = global_mass - target_mass;
    iterations = iteration + 1U;
    if (!finite(global_mass) || global_mass <= 0.0 ||
        !finite(global_derivative) || global_derivative <= 0.0 ||
        !finite(residual)) {
      return {StatusCode::numerical_failure, kClosedMassInput};
    }
    if (std::abs(residual) <= relative_tolerance_ * target_mass) {
      converged = true;
      break;
    }
    double step = -residual / global_derivative;
    const double limit = maximum_relative_step_ * pressure;
    step = std::max(-limit, std::min(limit, step));
    const double next = pressure + step;
    if (!finite(next) || next <= 0.0 || step == 0.0) {
      return {StatusCode::numerical_failure, kClosedMassConvergence};
    }
    pressure = next;
  }
  if (!converged) {
    return {StatusCode::numerical_failure, kClosedMassConvergence};
  }
  local = {};
  for (std::size_t cell = 0U; cell < count; ++cell) {
    if (cells.active.data[cell] != 0U &&
        (!finite(pressure + cells.pressure_perturbation.data[cell]) ||
         pressure + cells.pressure_perturbation.data[cell] <= 0.0)) {
      local = {StatusCode::numerical_failure, kClosedMassConvergence};
      break;
    }
  }
  collective = collective_failure(communicator, local, lowest);
  if (!collective) {
    return collective;
  }
  ClosedMassResult candidate;
  candidate.pressure_reference = pressure;
  candidate.mass = global_mass;
  candidate.residual = global_mass - target_mass;
  candidate.iterations = iterations;
  candidate.lowest_failing_rank = -1;
  out = candidate;
  return {};
}

Status ClosedMassPlan::solve_fields(
    MPI_Comm communicator, const ThermodynamicsPlan& thermodynamics,
    const ClosedMassFieldView& cells, double target_mass,
    double current_pressure_reference, ClosedMassResult& out) const noexcept {
  int size = 0;
  if (communicator == MPI_COMM_NULL ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return {StatusCode::mpi_failure, kClosedMassCollective};
  }
  const std::size_t independent = thermodynamics.independent_species_count();
  const Int3 shape = cells.patch.cells;
  const Int3 global = cells.geometry == nullptr
                          ? Int3{}
                          : cells.geometry->global_cells();
  std::size_t count = 0U;
  const bool count_valid =
      shape.x > 0 && shape.y > 0 && shape.z > 0 &&
      static_cast<std::size_t>(shape.x) <=
          std::numeric_limits<std::size_t>::max() /
              static_cast<std::size_t>(shape.y) &&
      static_cast<std::size_t>(shape.x) *
              static_cast<std::size_t>(shape.y) <=
          std::numeric_limits<std::size_t>::max() /
              static_cast<std::size_t>(shape.z);
  if (count_valid) {
    count = static_cast<std::size_t>(shape.x) *
            static_cast<std::size_t>(shape.y) *
            static_cast<std::size_t>(shape.z);
  }
  const auto valid_scalar_field = [&](ConstFieldView field) {
    detail::FieldStorageInterval interval{};
    return field.base != nullptr && field.interior.x == shape.x &&
           field.interior.y == shape.y && field.interior.z == shape.z &&
           field.components == 1U && field.revision != 0U &&
           field.storage_identity != 0U && field.revision_domain != 0U &&
           detail::field_storage_interval(field, interval);
  };
  bool species_valid =
      cells.independent_mass_fractions.size == independent &&
      (independent == 0U ||
       cells.independent_mass_fractions.data != nullptr);
  for (std::size_t species = 0U;
       species < cells.independent_mass_fractions.size && species_valid;
       ++species) {
    species_valid =
        valid_scalar_field(cells.independent_mass_fractions.data[species]);
  }
  const bool patch_valid =
      cells.geometry != nullptr && cells.geometry->fingerprint() != 0U &&
      cells.patch.begin.x >= 0 && cells.patch.begin.y >= 0 &&
      cells.patch.begin.z >= 0 && cells.patch.begin.x <= global.x - shape.x &&
      cells.patch.begin.y <= global.y - shape.y &&
      cells.patch.begin.z <= global.z - shape.z;
  Status local{};
  if (authority_ != PressureReferenceKind::closed_mass ||
      thermodynamics.fingerprint() == 0U || !finite(target_mass) ||
      target_mass <= 0.0 || !finite(current_pressure_reference) ||
      current_pressure_reference <= 0.0 || !count_valid || count == 0U ||
      !patch_valid || !valid_scalar_field(cells.pressure_perturbation) ||
      !valid_scalar_field(cells.enthalpy) || !species_valid ||
      (cells.active.size != 0U &&
       (cells.active.data == nullptr || cells.active.size != count))) {
    local = {StatusCode::invalid_plan, kClosedMassInput};
  }
  int lowest = -1;
  Status collective = collective_failure(communicator, local, lowest);
  if (!collective) return collective;

  std::uint64_t target_bits = 0U;
  std::uint64_t pressure_bits = 0U;
  std::memcpy(&target_bits, &target_mass, sizeof(target_bits));
  std::memcpy(&pressure_bits, &current_pressure_reference,
              sizeof(pressure_bits));
  const std::array<std::uint64_t, 5U> identity{
      thermodynamics.fingerprint(), fingerprint_, cells.geometry->fingerprint(),
      target_bits, pressure_bits};
  std::array<std::uint64_t, identity.size()> minimum_identity{};
  std::array<std::uint64_t, identity.size()> maximum_identity{};
  if (MPI_Allreduce(identity.data(), minimum_identity.data(),
                    static_cast<int>(identity.size()), MPI_UINT64_T, MPI_MIN,
                    communicator) != MPI_SUCCESS ||
      MPI_Allreduce(identity.data(), maximum_identity.data(),
                    static_cast<int>(identity.size()), MPI_UINT64_T, MPI_MAX,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kClosedMassCollective};
  }
  if (minimum_identity != maximum_identity)
    return {StatusCode::invalid_plan, kClosedMassCollective};

  CompensatedSum local_slope_sum;
  CompensatedSum local_intercept_sum;
  std::array<double, kMaximumSpecies - 1U> composition{};
  local = {};
  std::size_t flat = 0U;
  for (std::int32_t z = 0; z < shape.z && local; ++z) {
    for (std::int32_t y = 0; y < shape.y && local; ++y) {
      for (std::int32_t x = 0; x < shape.x; ++x, ++flat) {
        if (cells.active.size != 0U && cells.active.data[flat] > 1U) {
          local = {StatusCode::invalid_plan, kClosedMassInput};
          break;
        }
        if (cells.active.size != 0U && cells.active.data[flat] == 0U)
          continue;
        const Int3 cell{x, y, z};
        const double volume =
            cells.geometry->x().widths().data[
                static_cast<std::size_t>(cells.patch.begin.x + x)] *
            cells.geometry->y().widths().data[
                static_cast<std::size_t>(cells.patch.begin.y + y)] *
            cells.geometry->z().widths().data[
                static_cast<std::size_t>(cells.patch.begin.z + z)];
        const double enthalpy = cells.enthalpy.unchecked(cell, 0U);
        const double pi = cells.pressure_perturbation.unchecked(cell, 0U);
        const double initial_absolute_pressure =
            current_pressure_reference + pi;
        if (!finite(volume) || volume <= 0.0 || !finite(enthalpy) ||
            !finite(pi) || !finite(initial_absolute_pressure) ||
            initial_absolute_pressure <= 0.0) {
          local = {StatusCode::numerical_failure, kClosedMassInput};
          break;
        }
        for (std::size_t species = 0U; species < independent; ++species) {
          composition[species] = cells.independent_mass_fractions.data[species]
                                     .unchecked(cell, 0U);
        }
        ThermalState thermal;
        const Status evaluated = thermodynamics.evaluate_thermal_impl(
            enthalpy, {composition.data(), independent}, {}, false, thermal,
            std::numeric_limits<double>::quiet_NaN());
        if (!evaluated) {
          local = evaluated;
          break;
        }
        const double slope = volume * thermal.drho_dp_hY_;
        const double intercept = pi * slope;
        if (!finite(slope) || slope <= 0.0 || !finite(intercept)) {
          local = {StatusCode::numerical_failure, kClosedMassInput};
          break;
        }
        local_slope_sum.add(slope);
        local_intercept_sum.add(intercept);
      }
    }
  }
  collective = collective_failure(communicator, local, lowest);
  if (!collective) return collective;
  const std::array<double, 2U> local_coefficients{
      local_slope_sum.value(), local_intercept_sum.value()};
  std::array<double, 2U> global_coefficients{};
  if (MPI_Allreduce(local_coefficients.data(), global_coefficients.data(), 2,
                    MPI_DOUBLE, MPI_SUM, communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kClosedMassCollective};
  }
  const double global_derivative = global_coefficients[0U];
  const double global_intercept = global_coefficients[1U];
  if (!finite(global_derivative) || global_derivative <= 0.0 ||
      !finite(global_intercept)) {
    return {StatusCode::numerical_failure, kClosedMassInput};
  }

  double pressure = current_pressure_reference;
  double global_mass = 0.0;
  std::uint32_t iterations = 0U;
  bool converged = false;
  for (std::uint32_t iteration = 0U; iteration < maximum_iterations_;
       ++iteration) {
    global_mass = pressure * global_derivative + global_intercept;
    const double residual = global_mass - target_mass;
    iterations = iteration + 1U;
    if (!finite(global_mass) || global_mass <= 0.0 || !finite(residual))
      return {StatusCode::numerical_failure, kClosedMassInput};
    if (std::abs(residual) <= relative_tolerance_ * target_mass) {
      converged = true;
      break;
    }
    double step = -residual / global_derivative;
    const double limit = maximum_relative_step_ * pressure;
    step = std::max(-limit, std::min(limit, step));
    const double next = pressure + step;
    if (!finite(next) || next <= 0.0 || step == 0.0)
      return {StatusCode::numerical_failure, kClosedMassConvergence};
    pressure = next;
  }
  if (!converged)
    return {StatusCode::numerical_failure, kClosedMassConvergence};

  local = {};
  flat = 0U;
  for (std::int32_t z = 0; z < shape.z && local; ++z) {
    for (std::int32_t y = 0; y < shape.y && local; ++y) {
      for (std::int32_t x = 0; x < shape.x; ++x, ++flat) {
        if (cells.active.size != 0U && cells.active.data[flat] == 0U)
          continue;
        const double absolute =
            pressure + cells.pressure_perturbation.unchecked({x, y, z}, 0U);
        if (!finite(absolute) || absolute <= 0.0) {
          local = {StatusCode::numerical_failure, kClosedMassConvergence};
          break;
        }
      }
    }
  }
  collective = collective_failure(communicator, local, lowest);
  if (!collective) return collective;
  const ClosedMassResult candidate{pressure, global_mass,
                                   global_mass - target_mass, iterations, -1};
  out = candidate;
  return {};
}

Status ClosedMassPlan::certify_density_fields(
    MPI_Comm communicator, const ClosedMassDensityFieldView& cells,
    double target_mass, double current_pressure_reference,
    ClosedMassResult& out) const noexcept {
  int size = 0;
  if (communicator == MPI_COMM_NULL ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return {StatusCode::mpi_failure, kClosedMassCollective};
  }
  const Int3 shape = cells.patch.cells;
  const Int3 global = cells.geometry == nullptr
                          ? Int3{}
                          : cells.geometry->global_cells();
  std::size_t count = 0U;
  const bool count_valid =
      shape.x > 0 && shape.y > 0 && shape.z > 0 &&
      static_cast<std::size_t>(shape.x) <=
          std::numeric_limits<std::size_t>::max() /
              static_cast<std::size_t>(shape.y) &&
      static_cast<std::size_t>(shape.x) *
              static_cast<std::size_t>(shape.y) <=
          std::numeric_limits<std::size_t>::max() /
              static_cast<std::size_t>(shape.z);
  if (count_valid) {
    count = static_cast<std::size_t>(shape.x) *
            static_cast<std::size_t>(shape.y) *
            static_cast<std::size_t>(shape.z);
  }
  const auto valid_scalar_field = [&](ConstFieldView field) {
    detail::FieldStorageInterval interval{};
    return field.base != nullptr && field.interior.x == shape.x &&
           field.interior.y == shape.y && field.interior.z == shape.z &&
           field.components == 1U && field.revision != 0U &&
           field.storage_identity != 0U && field.revision_domain != 0U &&
           detail::field_storage_interval(field, interval);
  };
  const bool patch_valid =
      cells.geometry != nullptr && cells.geometry->fingerprint() != 0U &&
      cells.patch.begin.x >= 0 && cells.patch.begin.y >= 0 &&
      cells.patch.begin.z >= 0 && cells.patch.begin.x <= global.x - shape.x &&
      cells.patch.begin.y <= global.y - shape.y &&
      cells.patch.begin.z <= global.z - shape.z;
  Status local{};
  if (authority_ != PressureReferenceKind::closed_mass ||
      fingerprint_ == 0U || cells.predictor_state == 0U ||
      !finite(target_mass) || target_mass <= 0.0 ||
      !finite(current_pressure_reference) ||
      current_pressure_reference <= 0.0 || !count_valid || count == 0U ||
      !patch_valid || !valid_scalar_field(cells.pressure_perturbation) ||
      !valid_scalar_field(cells.density) ||
      !valid_scalar_field(cells.pressure_compressibility) ||
      (cells.active.size != 0U &&
       (cells.active.data == nullptr || cells.active.size != count))) {
    local = {StatusCode::invalid_plan, kClosedMassInput};
  }
  int lowest = -1;
  Status collective = collective_failure(communicator, local, lowest);
  if (!collective) return collective;

  std::uint64_t target_bits = 0U;
  std::uint64_t pressure_bits = 0U;
  std::memcpy(&target_bits, &target_mass, sizeof(target_bits));
  std::memcpy(&pressure_bits, &current_pressure_reference,
              sizeof(pressure_bits));
  // The predictor certificate is validated locally by the coupling wrapper.
  // Its state intentionally binds rank-local storage identities, revision
  // domains, and face addresses, so it is not a collective identity token.
  // Compare only the semantic inputs that must be identical on every rank.
  const std::array<std::uint64_t, 4U> identity{
      fingerprint_, cells.geometry->fingerprint(), target_bits, pressure_bits};
  std::array<std::uint64_t, identity.size()> minimum_identity{};
  std::array<std::uint64_t, identity.size()> maximum_identity{};
  if (MPI_Allreduce(identity.data(), minimum_identity.data(),
                    static_cast<int>(identity.size()), MPI_UINT64_T, MPI_MIN,
                    communicator) != MPI_SUCCESS ||
      MPI_Allreduce(identity.data(), maximum_identity.data(),
                    static_cast<int>(identity.size()), MPI_UINT64_T, MPI_MAX,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kClosedMassCollective};
  }
  if (minimum_identity != maximum_identity)
    return {StatusCode::invalid_plan, kClosedMassCollective};

  CompensatedSum local_mass_sum;
  local = {};
  std::size_t flat = 0U;
  for (std::int32_t z = 0; z < shape.z && local; ++z) {
    for (std::int32_t y = 0; y < shape.y && local; ++y) {
      for (std::int32_t x = 0; x < shape.x; ++x, ++flat) {
        if (cells.active.size != 0U && cells.active.data[flat] > 1U) {
          local = {StatusCode::invalid_plan, kClosedMassInput};
          break;
        }
        if (cells.active.size != 0U && cells.active.data[flat] == 0U)
          continue;
        const Int3 cell{x, y, z};
        const double volume =
            cells.geometry->x().widths().data[
                static_cast<std::size_t>(cells.patch.begin.x + x)] *
            cells.geometry->y().widths().data[
                static_cast<std::size_t>(cells.patch.begin.y + y)] *
            cells.geometry->z().widths().data[
                static_cast<std::size_t>(cells.patch.begin.z + z)];
        const double pi = cells.pressure_perturbation.unchecked(cell, 0U);
        const double density = cells.density.unchecked(cell, 0U);
        const double psi =
            cells.pressure_compressibility.unchecked(cell, 0U);
        const double absolute_pressure = current_pressure_reference + pi;
        const double oracle_density = absolute_pressure * psi;
        const double scale =
            std::max({1.0, density, std::abs(oracle_density)});
        if (!finite(volume) || volume <= 0.0 || !finite(pi) ||
            !finite(density) || density <= 0.0 || !finite(psi) ||
            psi <= 0.0 || !finite(absolute_pressure) ||
            absolute_pressure <= 0.0 || !finite(oracle_density) ||
            oracle_density <= 0.0 ||
            std::abs(oracle_density - density) >
                128.0 * std::numeric_limits<double>::epsilon() * scale) {
          local = {StatusCode::numerical_failure, kClosedMassInput};
          break;
        }
        local_mass_sum.add(volume * density);
      }
    }
  }
  collective = collective_failure(communicator, local, lowest);
  if (!collective) return collective;
  const double local_mass = local_mass_sum.value();
  double global_mass = 0.0;
  if (MPI_Allreduce(&local_mass, &global_mass, 1, MPI_DOUBLE, MPI_SUM,
                    communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kClosedMassCollective};
  }
  const double residual = global_mass - target_mass;
  if (!finite(global_mass) || global_mass <= 0.0 || !finite(residual) ||
      std::abs(residual) > relative_tolerance_ * target_mass) {
    return {StatusCode::numerical_failure, kClosedMassConvergence};
  }
  out = {current_pressure_reference, global_mass, residual, 1U, -1};
  return {};
}

}  // namespace hundun::v04

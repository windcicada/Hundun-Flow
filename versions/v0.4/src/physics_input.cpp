// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_physics.hpp"

#include "physics_input_detail.hpp"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <locale>
#include <new>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kInput = 801U;
constexpr std::uint32_t kSyntax = 802U;
constexpr std::uint32_t kValue = 803U;
constexpr std::uint32_t kCollective = 806U;
constexpr std::size_t kMaximumSpecies = 65U;
constexpr std::size_t kMaximumTextBytes = 4U * 1024U * 1024U;
constexpr std::size_t kMaximumNameBytes = 255U;
constexpr std::uint32_t kMaximumTemperatureIterations = 256U;
constexpr std::uint32_t kMaximumClosedMassIterations = 256U;
constexpr unsigned kMaximumBernsteinDepth = 48U;
constexpr std::size_t kMaximumBernsteinNodes = 4096U;
constexpr long double kContinuityRelativeTolerance = 1.0e-8L;
constexpr std::uint64_t kFnvOffset = 14695981039346656037ULL;
constexpr std::uint64_t kFnvPrime = 1099511628211ULL;

Status invalid(std::uint32_t detail) noexcept {
  return {StatusCode::invalid_case, detail};
}

bool valid_name(std::string_view value) noexcept {
  if (value.empty() || value.size() > kMaximumNameBytes) {
    return false;
  }
  for (const unsigned char character : value) {
    if (!((character >= 'a' && character <= 'z') ||
          (character >= 'A' && character <= 'Z') ||
          (character >= '0' && character <= '9') || character == '_' ||
          character == '-')) {
      return false;
    }
  }
  return true;
}

bool finite_positive(double value) noexcept {
  return std::isfinite(value) && value > 0.0;
}

bool parse_real(std::string_view token, double& out) noexcept {
  if (token.empty() || token.size() > 127U ||
      token.find('\0') != std::string_view::npos) {
    return false;
  }
  try {
    std::istringstream stream{std::string(token)};
    stream.imbue(std::locale::classic());
    double candidate = 0.0;
    stream >> std::noskipws >> candidate;
    if (stream.fail() ||
        stream.peek() != std::char_traits<char>::eof() ||
        !std::isfinite(candidate)) {
      return false;
    }
    out = candidate;
    return true;
  } catch (...) {
    return false;
  }
}

bool parse_u32(std::string_view token, std::uint32_t& out) noexcept {
  if (token.empty()) {
    return false;
  }
  const char* const begin = token.data();
  const char* const end = token.data() + token.size();
  const auto parsed = std::from_chars(begin, end, out, 10);
  return parsed.ec == std::errc{} && parsed.ptr == end;
}

Status tokenize(std::string_view text, std::vector<std::string>& out) {
  try {
    out.clear();
    out.reserve(std::min<std::size_t>(text.size() / 4U + 1U, 4096U));
    std::size_t position = 0U;
    while (position < text.size()) {
      while (position < text.size()) {
        const char current = text[position];
        if (current == '#') {
          while (position < text.size() && text[position] != '\n') {
            ++position;
          }
        } else if (current == ' ' || current == '\t' || current == '\r' ||
                   current == '\n') {
          ++position;
        } else {
          break;
        }
      }
      if (position == text.size()) {
        break;
      }
      const std::size_t begin = position;
      while (position < text.size()) {
        const char current = text[position];
        if (current == '#' || current == ' ' || current == '\t' ||
            current == '\r' || current == '\n') {
          break;
        }
        ++position;
      }
      if (position == begin) {
        return invalid(kSyntax);
      }
      out.emplace_back(text.substr(begin, position - begin));
    }
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, kInput};
  } catch (...) {
    return invalid(kSyntax);
  }
}

class Tokens {
 public:
  explicit Tokens(const std::vector<std::string>& tokens) noexcept
      : tokens_(tokens) {}

  bool exact(std::string_view expected) noexcept {
    if (position_ >= tokens_.size() || tokens_[position_] != expected) {
      return false;
    }
    ++position_;
    return true;
  }
  bool text(std::string& out) {
    if (position_ >= tokens_.size()) {
      return false;
    }
    out = tokens_[position_++];
    return true;
  }
  bool real(double& out) noexcept {
    return position_ < tokens_.size() &&
           parse_real(tokens_[position_++], out);
  }
  bool u32(std::uint32_t& out) noexcept {
    return position_ < tokens_.size() &&
           parse_u32(tokens_[position_++], out);
  }
  bool finished() const noexcept { return position_ == tokens_.size(); }

 private:
  const std::vector<std::string>& tokens_;
  std::size_t position_{};
};

using BernsteinQuartic = std::array<long double, 5U>;

long double dimensionless_cp(
    const std::array<double, 7U>& coefficients,
    long double temperature) noexcept {
  return (((static_cast<long double>(coefficients[4U]) * temperature +
            static_cast<long double>(coefficients[3U])) *
               temperature +
           static_cast<long double>(coefficients[2U])) *
              temperature +
          static_cast<long double>(coefficients[1U])) *
             temperature +
         static_cast<long double>(coefficients[0U]);
}

long double dimensionless_h(
    const std::array<double, 7U>& coefficients,
    long double temperature) noexcept {
  return temperature *
             (static_cast<long double>(coefficients[0U]) +
              temperature *
                  (0.5L * static_cast<long double>(coefficients[1U]) +
                   temperature *
                       (static_cast<long double>(coefficients[2U]) / 3.0L +
                        temperature *
                            (0.25L *
                                 static_cast<long double>(coefficients[3U]) +
                             temperature *
                                 (0.2L * static_cast<long double>(
                                             coefficients[4U])))))) +
         static_cast<long double>(coefficients[5U]);
}

bool finite_bernstein(const BernsteinQuartic& coefficients) noexcept {
  for (const long double coefficient : coefficients) {
    if (!std::isfinite(coefficient)) {
      return false;
    }
  }
  return true;
}

bool prove_positive_bernstein(const BernsteinQuartic& coefficients,
                              unsigned depth,
                              std::size_t& visited_nodes) noexcept {
  if (!finite_bernstein(coefficients) ||
      ++visited_nodes > kMaximumBernsteinNodes) {
    return false;
  }
  long double minimum = coefficients[0U];
  long double maximum = coefficients[0U];
  long double scale = 1.0L;
  for (const long double coefficient : coefficients) {
    minimum = std::min(minimum, coefficient);
    maximum = std::max(maximum, coefficient);
    scale = std::max(scale, std::abs(coefficient));
  }
  const long double margin =
      64.0L * std::numeric_limits<long double>::epsilon() * scale;
  if (minimum > margin) {
    return true;
  }
  if (maximum <= margin || depth == kMaximumBernsteinDepth) {
    return false;
  }

  BernsteinQuartic work = coefficients;
  BernsteinQuartic left{};
  BernsteinQuartic right{};
  left[0U] = work[0U];
  right[4U] = work[4U];
  for (std::size_t level = 1U; level <= 4U; ++level) {
    const std::size_t remaining = 5U - level;
    for (std::size_t index = 0U; index < remaining; ++index) {
      work[index] = 0.5L * (work[index] + work[index + 1U]);
    }
    left[level] = work[0U];
    right[4U - level] = work[remaining - 1U];
  }
  return prove_positive_bernstein(left, depth + 1U, visited_nodes) &&
         prove_positive_bernstein(right, depth + 1U, visited_nodes);
}

bool prove_cp_greater_than_r(
    const std::array<double, 7U>& coefficients, double lower_temperature,
    double upper_temperature) noexcept {
  const long double lower = lower_temperature;
  const long double delta =
      static_cast<long double>(upper_temperature) - lower;
  if (!(delta > 0.0L) || !std::isfinite(delta)) {
    return false;
  }

  const long double a0 =
      static_cast<long double>(coefficients[0U]) - 1.0L;
  const long double a1 = coefficients[1U];
  const long double a2 = coefficients[2U];
  const long double a3 = coefficients[3U];
  const long double a4 = coefficients[4U];
  const long double delta2 = delta * delta;
  const long double delta3 = delta2 * delta;
  const long double delta4 = delta3 * delta;
  const long double lower2 = lower * lower;
  const std::array<long double, 5U> power{
      (((a4 * lower + a3) * lower + a2) * lower + a1) * lower + a0,
      delta * (((4.0L * a4 * lower + 3.0L * a3) * lower +
                 2.0L * a2) *
                    lower +
                a1),
      delta2 * (6.0L * a4 * lower2 + 3.0L * a3 * lower + a2),
      delta3 * (4.0L * a4 * lower + a3), delta4 * a4};
  const BernsteinQuartic bernstein{
      power[0U],
      power[0U] + 0.25L * power[1U],
      power[0U] + 0.5L * power[1U] + power[2U] / 6.0L,
      power[0U] + 0.75L * power[1U] + 0.5L * power[2U] +
          0.25L * power[3U],
      power[0U] + power[1U] + power[2U] + power[3U] + power[4U]};
  std::size_t visited_nodes = 0U;
  return prove_positive_bernstein(bernstein, 0U, visited_nodes);
}

bool continuous(long double left, long double right) noexcept {
  if (!std::isfinite(left) || !std::isfinite(right)) {
    return false;
  }
  const long double scale =
      std::max({1.0L, std::abs(left), std::abs(right)});
  return std::abs(left - right) <= kContinuityRelativeTolerance * scale;
}

double canonical_high_enthalpy_constant(
    const SpeciesThermophysicalSpec& species) noexcept {
  const double temperature = species.temperature_switch;
  const double t2 = temperature * temperature;
  const double t3 = t2 * temperature;
  const double t4 = t3 * temperature;
  const auto without_constant =
      [&](const std::array<double, 7U>& coefficients) noexcept {
        return coefficients[0U] * temperature +
               coefficients[1U] * t2 * 0.5 +
               coefficients[2U] * t3 / 3.0 +
               coefficients[3U] * t4 * 0.25 +
               coefficients[4U] * t4 * temperature * 0.2;
      };
  const double low_at_switch =
      without_constant(species.nasa7_low) + species.nasa7_low[5U];
  return low_at_switch - without_constant(species.nasa7_high);
}

bool valid_spec(const ThermophysicalSpec& spec) noexcept {
  if (!finite_positive(spec.minimum_temperature) ||
      !finite_positive(spec.maximum_temperature) ||
      spec.minimum_temperature >= spec.maximum_temperature ||
      !finite_positive(spec.temperature_relative_tolerance) ||
      spec.temperature_relative_tolerance >= 1.0 ||
      spec.maximum_temperature_iterations == 0U ||
      spec.maximum_temperature_iterations > kMaximumTemperatureIterations ||
      !finite_positive(spec.closed_mass_relative_tolerance) ||
      spec.closed_mass_relative_tolerance >= 1.0 ||
      spec.maximum_closed_mass_iterations == 0U ||
      spec.maximum_closed_mass_iterations > kMaximumClosedMassIterations ||
      !finite_positive(spec.maximum_closed_mass_relative_step) ||
      spec.maximum_closed_mass_relative_step >= 1.0 || spec.species.empty() ||
      spec.species.size() > kMaximumSpecies) {
    return false;
  }
  for (std::size_t species_index = 0U;
       species_index < spec.species.size(); ++species_index) {
    const SpeciesThermophysicalSpec& species = spec.species[species_index];
    if (!valid_name(species.stable_name) ||
        !finite_positive(species.molecular_weight) ||
        !finite_positive(species.temperature_switch) ||
        species.temperature_switch <= spec.minimum_temperature ||
        species.temperature_switch >= spec.maximum_temperature ||
        !finite_positive(species.viscosity_reference)) {
      return false;
    }
    for (std::size_t prior = 0U; prior < species_index; ++prior) {
      if (spec.species[prior].stable_name == species.stable_name) {
        return false;
      }
    }
    for (const double coefficient : species.nasa7_low) {
      if (!std::isfinite(coefficient)) {
        return false;
      }
    }
    for (const double coefficient : species.nasa7_high) {
      if (!std::isfinite(coefficient)) {
        return false;
      }
    }
    if (species.transport_law == TransportLaw::constant) {
      if (!finite_positive(species.conductivity) ||
          species.transport_reference_temperature != 0.0 ||
          species.sutherland_temperature != 0.0 || species.prandtl != 0.0) {
        return false;
      }
    } else if (species.transport_law == TransportLaw::sutherland) {
      if (!finite_positive(species.transport_reference_temperature) ||
          species.sutherland_temperature < 0.0 ||
          !std::isfinite(species.sutherland_temperature) ||
          !finite_positive(species.prandtl) || species.conductivity != 0.0) {
        return false;
      }
    } else {
      return false;
    }
    if (!prove_cp_greater_than_r(species.nasa7_low,
                                 spec.minimum_temperature,
                                 species.temperature_switch) ||
        !prove_cp_greater_than_r(species.nasa7_high,
                                 species.temperature_switch,
                                 spec.maximum_temperature)) {
      return false;
    }
    const long double switch_temperature = species.temperature_switch;
    if (!continuous(dimensionless_cp(species.nasa7_low,
                                     switch_temperature),
                    dimensionless_cp(species.nasa7_high,
                                     switch_temperature)) ||
        !continuous(dimensionless_h(species.nasa7_low,
                                    switch_temperature),
                    dimensionless_h(species.nasa7_high,
                                    switch_temperature))) {
      return false;
    }
  }
  return true;
}

bool canonicalize_spec(ThermophysicalSpec& spec) noexcept {
  if (!valid_spec(spec)) {
    return false;
  }
  for (SpeciesThermophysicalSpec& species : spec.species) {
    const double constant = canonical_high_enthalpy_constant(species);
    if (!std::isfinite(constant)) {
      return false;
    }
    species.nasa7_high[5U] = constant;
  }
  return valid_spec(spec);
}

class SpecHash64 final {
 public:
  void integer(std::uint64_t value) noexcept {
    for (unsigned byte_index = 0U; byte_index < 8U; ++byte_index) {
      byte(static_cast<std::uint8_t>(value & 0xffU));
      value >>= 8U;
    }
  }

  void real(double value) noexcept {
    if (value == 0.0) {
      value = 0.0;
    }
    std::uint64_t bits = 0U;
    static_assert(sizeof(bits) == sizeof(value),
                  "double fingerprint requires 64-bit binary storage");
    std::memcpy(&bits, &value, sizeof(bits));
    integer(bits);
  }

  void text(std::string_view value) noexcept {
    integer(value.size());
    for (const unsigned char character : value) {
      byte(character);
    }
  }

  PlanFingerprint finish() const noexcept {
    return value_ == 0U ? 1U : value_;
  }

 private:
  void byte(std::uint8_t value) noexcept {
    value_ ^= value;
    value_ *= kFnvPrime;
  }

  std::uint64_t value_{kFnvOffset};
};

PlanFingerprint fingerprint_spec(const ThermophysicalSpec& spec) noexcept {
  if (spec.data_file.empty() || !valid_spec(spec)) {
    return 0U;
  }
  try {
    ThermophysicalSpec canonical = spec;
    if (!canonicalize_spec(canonical)) {
      return 0U;
    }
    const std::string data_file = canonical.data_file.generic_string();
    if (data_file.empty()) {
      return 0U;
    }
    SpecHash64 hash;
    hash.text("HUNDUN-FLOW-v0.4-thermophysical-spec-v1");
    hash.text(data_file);
    hash.real(canonical.minimum_temperature);
    hash.real(canonical.maximum_temperature);
    hash.real(canonical.temperature_relative_tolerance);
    hash.integer(canonical.maximum_temperature_iterations);
    hash.real(canonical.closed_mass_relative_tolerance);
    hash.integer(canonical.maximum_closed_mass_iterations);
    hash.real(canonical.maximum_closed_mass_relative_step);
    hash.integer(canonical.species.size());
    for (const SpeciesThermophysicalSpec& species : canonical.species) {
      hash.text(species.stable_name);
      hash.real(species.molecular_weight);
      hash.real(species.temperature_switch);
      for (const double coefficient : species.nasa7_low) {
        hash.real(coefficient);
      }
      for (const double coefficient : species.nasa7_high) {
        hash.real(coefficient);
      }
      hash.integer(static_cast<std::uint64_t>(species.transport_law));
      hash.real(species.viscosity_reference);
      hash.real(species.transport_reference_temperature);
      hash.real(species.sutherland_temperature);
      hash.real(species.prandtl);
      hash.real(species.conductivity);
    }
    return hash.finish();
  } catch (...) {
    return 0U;
  }
}

PlanFingerprint fingerprint_scalar_catalog(
    Span<const TransportedScalarSpec> catalog) noexcept {
  if ((catalog.size != 0U && catalog.data == nullptr) ||
      catalog.size > 64U) {
    return 0U;
  }
  try {
    SpecHash64 hash;
    hash.text("HUNDUN-FLOW-v0.4-transported-scalar-catalog-v1");
    hash.integer(catalog.size);
    for (std::size_t index = 0U; index < catalog.size; ++index) {
      const TransportedScalarSpec& scalar = catalog.data[index];
      if (!valid_name(scalar.stable_name) ||
          static_cast<std::uint8_t>(scalar.role) >
              static_cast<std::uint8_t>(
                  TransportedScalarRole::passive_scalar)) {
        return 0U;
      }
      hash.text(scalar.stable_name);
      hash.integer(static_cast<std::uint64_t>(scalar.role));
    }
    return hash.finish();
  } catch (...) {
    return 0U;
  }
}

Status parse_text(std::string_view text, ThermophysicalSpec& out) noexcept {
  if (text.empty() || text.size() > kMaximumTextBytes) {
    return invalid(kInput);
  }
  try {
    std::vector<std::string> token_storage;
    const Status token_status = tokenize(text, token_storage);
    if (!token_status) {
      return token_status;
    }
    Tokens tokens(token_storage);
    ThermophysicalSpec candidate;
    std::uint32_t species_count = 0U;
    if (!tokens.exact("HUNDUN_THERMOPHYSICS_V1") ||
        !tokens.exact("temperature_bounds") ||
        !tokens.real(candidate.minimum_temperature) ||
        !tokens.real(candidate.maximum_temperature) ||
        !tokens.exact("temperature_inversion") ||
        !tokens.real(candidate.temperature_relative_tolerance) ||
        !tokens.u32(candidate.maximum_temperature_iterations) ||
        !tokens.exact("closed_mass_newton") ||
        !tokens.real(candidate.closed_mass_relative_tolerance) ||
        !tokens.u32(candidate.maximum_closed_mass_iterations) ||
        !tokens.real(candidate.maximum_closed_mass_relative_step) ||
        !tokens.exact("species_count") || !tokens.u32(species_count) ||
        species_count == 0U || species_count > kMaximumSpecies) {
      return invalid(kSyntax);
    }
    candidate.species.reserve(species_count);
    for (std::uint32_t species_index = 0U; species_index < species_count;
         ++species_index) {
      SpeciesThermophysicalSpec species;
      if (!tokens.exact("species") || !tokens.text(species.stable_name) ||
          !tokens.exact("molecular_weight") ||
          !tokens.real(species.molecular_weight) ||
          !tokens.exact("temperature_switch") ||
          !tokens.real(species.temperature_switch) ||
          !tokens.exact("nasa7_low")) {
        return invalid(kSyntax);
      }
      for (double& coefficient : species.nasa7_low) {
        if (!tokens.real(coefficient)) {
          return invalid(kSyntax);
        }
      }
      if (!tokens.exact("nasa7_high")) {
        return invalid(kSyntax);
      }
      for (double& coefficient : species.nasa7_high) {
        if (!tokens.real(coefficient)) {
          return invalid(kSyntax);
        }
      }
      if (tokens.exact("transport_constant")) {
        species.transport_law = TransportLaw::constant;
        if (!tokens.real(species.viscosity_reference) ||
            !tokens.real(species.conductivity)) {
          return invalid(kSyntax);
        }
      } else if (tokens.exact("transport_sutherland")) {
        species.transport_law = TransportLaw::sutherland;
        if (!tokens.real(species.viscosity_reference) ||
            !tokens.real(species.transport_reference_temperature) ||
            !tokens.real(species.sutherland_temperature) ||
            !tokens.real(species.prandtl)) {
          return invalid(kSyntax);
        }
      } else {
        return invalid(kSyntax);
      }
      if (!tokens.exact("end_species")) {
        return invalid(kSyntax);
      }
      candidate.species.push_back(std::move(species));
    }
    if (!tokens.exact("end") || !tokens.finished() ||
        !canonicalize_spec(candidate)) {
      return invalid(kValue);
    }
    out = std::move(candidate);
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, kInput};
  } catch (...) {
    return invalid(kSyntax);
  }
}

Status collective_status(MPI_Comm communicator, Status local, int rank,
                         int size, int& lowest) noexcept {
  const int candidate = local ? size : rank;
  lowest = size;
  if (MPI_Allreduce(&candidate, &lowest, 1, MPI_INT, MPI_MIN, communicator) !=
      MPI_SUCCESS) {
    lowest = -1;
    return {StatusCode::mpi_failure, kCollective};
  }
  if (lowest == size) {
    lowest = -1;
    return {};
  }
  std::array<std::uint64_t, 2U> selected{};
  if (rank == lowest) {
    selected[0] = static_cast<std::uint64_t>(local.code);
    selected[1] = local.detail;
  }
  if (MPI_Bcast(selected.data(), static_cast<int>(selected.size()),
                MPI_UINT64_T, lowest, communicator) != MPI_SUCCESS ||
      selected[0] > static_cast<std::uint64_t>(StatusCode::io_failure) ||
      selected[1] > UINT32_MAX) {
    return {StatusCode::mpi_failure, kCollective};
  }
  return {static_cast<StatusCode>(selected[0]),
          static_cast<std::uint32_t>(selected[1])};
}

}  // namespace

namespace detail {

bool valid_thermophysical_spec(const ThermophysicalSpec& spec) noexcept {
  return valid_spec(spec);
}

Status canonicalize_thermophysical_spec(
    ThermophysicalSpec& spec) noexcept {
  return canonicalize_spec(spec) ? Status{} : invalid(kValue);
}

PlanFingerprint thermophysical_spec_fingerprint(
    const ThermophysicalSpec& spec) noexcept {
  return fingerprint_spec(spec);
}

Status parse_thermophysical_text(std::string_view text,
                                 ThermophysicalSpec& out) noexcept {
  return parse_text(text, out);
}

}  // namespace detail

Status ThermophysicalCompiler::load_and_compile(
    MPI_Comm communicator, const ValidatedModel& model,
    ThermophysicalSpec& spec,
    ThermodynamicsPlan& thermodynamics, TransportPlan& transport,
    ThermophysicalCompileDiagnostics* diagnostics) noexcept {
  if (diagnostics != nullptr) {
    diagnostics->lowest_failing_rank = -1;
  }
  if (communicator == MPI_COMM_NULL) {
    return {StatusCode::invalid_plan, kCollective};
  }
  int rank = 0;
  int size = 0;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &size) != MPI_SUCCESS || size <= 0) {
    return {StatusCode::mpi_failure, kCollective};
  }

  // ValidatedModel is intentionally a plain public value type.  Do not trust
  // its case fingerprint alone after publication: include the complete typed
  // thermophysical payload in the collective identity so a rank-local
  // mutation cannot silently compile a divergent plan.
  const std::array<PlanFingerprint, 3U> local_identity{
      model.fingerprint,
      detail::thermophysical_spec_fingerprint(model.thermophysics),
      fingerprint_scalar_catalog(
          {model.transported_scalars.data(),
           model.transported_scalars.size()})};
  std::array<PlanFingerprint, 3U> authoritative_identity = local_identity;
  if (MPI_Bcast(authoritative_identity.data(),
                static_cast<int>(authoritative_identity.size()),
                MPI_UINT64_T, 0, communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kCollective};
  }
  Status local{};
  if (local_identity[0U] == 0U || local_identity[1U] == 0U ||
      local_identity[2U] == 0U ||
      local_identity != authoritative_identity) {
    local = {StatusCode::invalid_plan, kInput};
  }
  int lowest = -1;
  ThermophysicalSpec candidate_spec;
  ThermodynamicsPlan candidate_thermodynamics;
  TransportPlan candidate_transport;
  if (local) {
    try {
      candidate_spec = model.thermophysics;
    } catch (const std::bad_alloc&) {
      local = {StatusCode::allocation_failure, kInput};
    } catch (...) {
      local = {StatusCode::invalid_plan, kInput};
    }
  }
  if (local) {
    local = detail::canonicalize_thermophysical_spec(candidate_spec);
    if (!local) {
      local = {StatusCode::invalid_plan, kInput};
    }
  }
  if (local) {
    local = ThermodynamicsPlan::compile(
        candidate_spec,
        Span<const TransportedScalarSpec>{model.transported_scalars.data(),
                                          model.transported_scalars.size()},
        candidate_thermodynamics);
  }
  if (local) {
    local = TransportPlan::compile(candidate_spec, candidate_thermodynamics,
                                   candidate_transport);
  }
  const std::array<PlanFingerprint, 2U> local_plans{
      local ? candidate_thermodynamics.fingerprint() : 0U,
      local ? candidate_transport.fingerprint() : 0U};
  std::array<PlanFingerprint, 2U> authoritative_plans = local_plans;
  if (MPI_Bcast(authoritative_plans.data(),
                static_cast<int>(authoritative_plans.size()), MPI_UINT64_T,
                0, communicator) != MPI_SUCCESS) {
    return {StatusCode::mpi_failure, kCollective};
  }
  if (local && (local_plans != authoritative_plans ||
                local_plans[0U] == 0U || local_plans[1U] == 0U)) {
    local = {StatusCode::invalid_plan, kCollective};
  }
  Status consensus = collective_status(communicator, local, rank, size, lowest);
  if (!consensus) {
    if (diagnostics != nullptr) {
      diagnostics->lowest_failing_rank = lowest;
    }
    return consensus;
  }
  spec = std::move(candidate_spec);
  thermodynamics = std::move(candidate_thermodynamics);
  transport = std::move(candidate_transport);
  return {};
}

}  // namespace hundun::v04

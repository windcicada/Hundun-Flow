// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_cantera.hpp"
#include "physics_nasa_detail.hpp"
#include "models_orders_detail.hpp"
#include "models_roundoff_detail.hpp"
#include "models_kerosene_cantera_detail.hpp"

#include "cantera/base/Solution.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/NasaPoly2.h"
#include "cantera/thermo/Species.h"
#include "cantera/transport/Transport.h"
#include "cantera/zeroD/IdealGasConstPressureReactor.h"
#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/ReactorNet.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace hundun::v04::chemistry {
namespace {

struct Sha256 final {
  std::array<std::uint32_t, 8> state{
      UINT32_C(0x6a09e667), UINT32_C(0xbb67ae85), UINT32_C(0x3c6ef372),
      UINT32_C(0xa54ff53a), UINT32_C(0x510e527f), UINT32_C(0x9b05688c),
      UINT32_C(0x1f83d9ab), UINT32_C(0x5be0cd19)};
  std::array<std::uint8_t, 64> block{};
  std::uint64_t bytes{};
  std::size_t used{};
};

constexpr std::array<std::uint32_t, 64> kSha256{
    0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U, 0x3956c25bU,
    0x59f111f1U, 0x923f82a4U, 0xab1c5ed5U, 0xd807aa98U, 0x12835b01U,
    0x243185beU, 0x550c7dc3U, 0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U,
    0xc19bf174U, 0xe49b69c1U, 0xefbe4786U, 0x0fc19dc6U, 0x240ca1ccU,
    0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU, 0x983e5152U,
    0xa831c66dU, 0xb00327c8U, 0xbf597fc7U, 0xc6e00bf3U, 0xd5a79147U,
    0x06ca6351U, 0x14292967U, 0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU,
    0x53380d13U, 0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U,
    0xa2bfe8a1U, 0xa81a664bU, 0xc24b8b70U, 0xc76c51a3U, 0xd192e819U,
    0xd6990624U, 0xf40e3585U, 0x106aa070U, 0x19a4c116U, 0x1e376c08U,
    0x2748774cU, 0x34b0bcb5U, 0x391c0cb3U, 0x4ed8aa4aU, 0x5b9cca4fU,
    0x682e6ff3U, 0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U,
    0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U};

std::uint32_t rotate(std::uint32_t value, unsigned count) noexcept {
  return (value >> count) | (value << (32U - count));
}

void transform(Sha256 &hash) noexcept {
  std::array<std::uint32_t, 64> words{};
  for (std::size_t index = 0; index < 16U; ++index) {
    const std::size_t offset = 4U * index;
    words[index] = static_cast<std::uint32_t>(hash.block[offset]) << 24U |
                   static_cast<std::uint32_t>(hash.block[offset + 1U]) << 16U |
                   static_cast<std::uint32_t>(hash.block[offset + 2U]) << 8U |
                   static_cast<std::uint32_t>(hash.block[offset + 3U]);
  }
  for (std::size_t index = 16U; index < words.size(); ++index) {
    const std::uint32_t s0 = rotate(words[index - 15U], 7U) ^
                             rotate(words[index - 15U], 18U) ^
                             (words[index - 15U] >> 3U);
    const std::uint32_t s1 = rotate(words[index - 2U], 17U) ^
                             rotate(words[index - 2U], 19U) ^
                             (words[index - 2U] >> 10U);
    words[index] = words[index - 16U] + s0 + words[index - 7U] + s1;
  }
  auto a = hash.state[0];
  auto b = hash.state[1];
  auto c = hash.state[2];
  auto d = hash.state[3];
  auto e = hash.state[4];
  auto f = hash.state[5];
  auto g = hash.state[6];
  auto h = hash.state[7];
  for (std::size_t index = 0; index < words.size(); ++index) {
    const std::uint32_t s1 = rotate(e, 6U) ^ rotate(e, 11U) ^ rotate(e, 25U);
    const std::uint32_t choose = (e & f) ^ (~e & g);
    const std::uint32_t temporary1 =
        h + s1 + choose + kSha256[index] + words[index];
    const std::uint32_t s0 = rotate(a, 2U) ^ rotate(a, 13U) ^ rotate(a, 22U);
    const std::uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
    const std::uint32_t temporary2 = s0 + majority;
    h = g;
    g = f;
    f = e;
    e = d + temporary1;
    d = c;
    c = b;
    b = a;
    a = temporary1 + temporary2;
  }
  hash.state[0] += a;
  hash.state[1] += b;
  hash.state[2] += c;
  hash.state[3] += d;
  hash.state[4] += e;
  hash.state[5] += f;
  hash.state[6] += g;
  hash.state[7] += h;
}

void update(Sha256 &hash, const std::uint8_t *data, std::size_t size) noexcept {
  hash.bytes += static_cast<std::uint64_t>(size);
  for (std::size_t index = 0; index < size; ++index) {
    hash.block[hash.used++] = data[index];
    if (hash.used == hash.block.size()) {
      transform(hash);
      hash.used = 0U;
    }
  }
}

std::string finish(Sha256 hash) {
  const std::uint64_t bit_count = hash.bytes * 8U;
  hash.block[hash.used++] = 0x80U;
  if (hash.used > 56U) {
    std::fill(hash.block.begin() + static_cast<std::ptrdiff_t>(hash.used),
              hash.block.end(), 0U);
    transform(hash);
    hash.used = 0U;
  }
  std::fill(hash.block.begin() + static_cast<std::ptrdiff_t>(hash.used),
            hash.block.begin() + 56, 0U);
  for (std::size_t index = 0; index < 8U; ++index) {
    hash.block[63U - index] =
        static_cast<std::uint8_t>(bit_count >> (8U * index));
  }
  transform(hash);
  std::ostringstream output;
  output << std::hex << std::setfill('0');
  for (const auto word : hash.state) {
    output << std::setw(8) << word;
  }
  return output.str();
}

std::string file_sha256(const std::filesystem::path &path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::invalid_argument("could not open reacting mechanism");
  }
  Sha256 hash;
  std::array<std::uint8_t, 16384> buffer{};
  while (input) {
    input.read(reinterpret_cast<char *>(buffer.data()),
               static_cast<std::streamsize>(buffer.size()));
    const auto count = input.gcount();
    if (count > 0) {
      update(hash, buffer.data(), static_cast<std::size_t>(count));
    }
  }
  if (!input.eof()) {
    throw std::invalid_argument("could not read complete reacting mechanism");
  }
  return finish(hash);
}

CompositionIdentity composition(const Cantera::ThermoPhase &thermo) {
  CompositionIdentity result;
  result.element_names.reserve(thermo.nElements());
  for (std::size_t element = 0; element < thermo.nElements(); ++element) {
    result.element_names.push_back(thermo.elementName(element));
  }
  const auto weights = thermo.molecularWeights();
  result.species.reserve(thermo.nSpecies());
  for (std::size_t species = 0; species < thermo.nSpecies(); ++species) {
    SpeciesIdentity identity;
    identity.name = thermo.speciesName(species);
    identity.molecular_weight_kg_per_kmol = weights[species];
    identity.element_counts.reserve(thermo.nElements());
    for (std::size_t element = 0; element < thermo.nElements(); ++element) {
      const double count = thermo.nAtoms(species, element);
      const double rounded = std::round(count);
      if (count != rounded ||
          rounded >
              static_cast<double>(std::numeric_limits<std::int32_t>::max())) {
        throw std::invalid_argument(
            "mechanism contains unsupported elemental composition");
      }
      identity.element_counts.push_back(static_cast<std::int32_t>(rounded));
    }
    result.species.push_back(std::move(identity));
  }
  result.fingerprint = composition_identity_fingerprint(result);
  validate_composition_identity(result);
  return result;
}

struct Workspace final {
  std::shared_ptr<Cantera::Solution> solution;
  std::shared_ptr<Cantera::ThermoPhase> thermo;
  std::shared_ptr<Cantera::Kinetics> kinetics;
  std::shared_ptr<Cantera::Transport> transport;
  std::shared_ptr<Cantera::Reactor> reactor;
  std::unique_ptr<Cantera::ReactorNet> network;
  std::unique_ptr<detail::KeroseneMaterial> kerosene;
  std::vector<double> query_diffusion, query_enthalpies, query_rates;
  std::vector<double> advance_fractions, advance_delta;
  bool reference_tolerances_initialized{};
  std::vector<double> reference_tolerances;
};

void continue_enthalpy(Cantera::ThermoPhase& thermo) {
  for (std::size_t i = 0; i < thermo.nSpecies(); ++i) {
    auto species = thermo.species(i);
    const auto nasa = std::dynamic_pointer_cast<Cantera::NasaPoly2>(species->thermo);
    if (!nasa) continue;
    std::array<double, 15> coefficients{};
    std::size_t index{};
    int type{};
    double low_temperature{}, high_temperature{}, pressure{};
    nasa->reportParameters(index, type, low_temperature, high_temperature,
                           pressure, coefficients.data());
    std::array<double, 7> low{}, high{};
    std::copy_n(coefficients.data() + 1, 7, high.data());
    std::copy_n(coefficients.data() + 8, 7, low.data());
    coefficients[6] = hundun::v04::detail::nasa_high_enthalpy_constant(
        coefficients[0], low, high);
    species->thermo = std::make_shared<Cantera::NasaPoly2>(
        low_temperature, high_temperature, pressure, coefficients.data());
    thermo.modifySpecies(i, species);
  }
}

// Cantera's relative HP stopping criterion can leave a residual above the
// public absolute-h gate near the enthalpy reference zero. Finish against
// that same gate with fixed-composition Newton corrections.
void solve_pressure_enthalpy(Cantera::ThermoPhase& thermo,double target,double pressure) {
  thermo.setState_HP(target,pressure,1e-12);
  const double tolerance=1e-10*std::max(1.,std::abs(target));
  for(unsigned correction=0;correction<4;++correction) {
    const double residual=target-thermo.enthalpy_mass();
    if(std::isfinite(residual) && std::abs(residual)<=tolerance)return;
    const double cp=thermo.cp_mass();
    const double temperature=thermo.temperature()+residual/cp;
    if(!std::isfinite(cp) || cp<=0 || !std::isfinite(temperature) || temperature<=0)
      throw std::runtime_error("invalid PH correction");
    thermo.setState_TP(temperature,pressure);
  }
  if(std::abs(target-thermo.enthalpy_mass())>tolerance)
    throw std::runtime_error("PH correction residual exceeds query tolerance");
}

Workspace make_workspace(const std::filesystem::path &mechanism,
                         const std::string &phase, bool continuous_enthalpy) {
  Workspace result;
  result.solution =
      Cantera::newSolution(mechanism.string(), phase, "mixture-averaged");
  result.thermo = result.solution->thermo();
  if (continuous_enthalpy) continue_enthalpy(*result.thermo);
  detail::install_order_rates(result.solution);
  result.kinetics = result.solution->kinetics();
  result.kerosene = detail::kerosene_material(result.solution);
  result.transport = result.solution->transport();
  result.reactor = std::make_shared<Cantera::IdealGasConstPressureReactor>(
      result.solution, false);
  result.network = std::make_unique<Cantera::ReactorNet>(result.reactor);
  result.query_diffusion.resize(result.thermo->nSpecies());
  result.query_enthalpies.resize(result.thermo->nSpecies());
  result.query_rates.resize(result.thermo->nSpecies());
  result.advance_fractions.resize(result.thermo->nSpecies());
  result.advance_delta.resize(result.thermo->nSpecies());
  return result;
}

} // namespace

struct CanteraBackendRuntime::Impl final {
  std::filesystem::path mechanism;
  std::string sha256;
  std::string phase;
  std::string_view reaction_model;
  bool continuous_enthalpy{};
  double minimum_temperature{}, maximum_temperature{};
  CompositionIdentity composition;
  portable::GasIdentity gas_identity;
  combustion::ChemistryIdentity closure_identity;
};

struct CanteraWorkspacePool::Impl final {
  std::shared_ptr<int> lifetime{std::make_shared<int>(0)};
  std::weak_ptr<const CanteraBackendRuntime> runtime;
  std::vector<Workspace> workspaces;
  std::size_t next_lane{};
};

struct CanteraBackend::Impl final {
  std::shared_ptr<const CanteraBackendRuntime> runtime;
  CanteraWorkspacePool::Impl *pool{};
  std::weak_ptr<int> pool_lifetime;
  std::size_t lane{};
  ChemistrySolverConfig controls;
};

CanteraBackendRuntime::CanteraBackendRuntime(const CanteraBackendConfig &config)
    : impl_(std::make_unique<Impl>()) {
  impl_->mechanism = config.mechanism.file;
  impl_->sha256 = file_sha256(impl_->mechanism);
  if (impl_->sha256 != config.mechanism.sha256) {
    throw std::invalid_argument("reacting mechanism SHA-256 mismatch");
  }
  impl_->phase = config.mechanism.phase;
  impl_->continuous_enthalpy = config.continuous_enthalpy;
  const bool source_interval = config.minimum_temperature == 0 && config.maximum_temperature == 0;
  if (!source_interval && (!std::isfinite(config.minimum_temperature) ||
      !std::isfinite(config.maximum_temperature) || config.minimum_temperature<=0 ||
      config.maximum_temperature<=config.minimum_temperature))
    throw std::invalid_argument("invalid thermodynamic admission interval");
  impl_->minimum_temperature=config.minimum_temperature;
  impl_->maximum_temperature=config.maximum_temperature;
  auto probe = Cantera::newSolution(impl_->mechanism.string(), impl_->phase,
                                    "mixture-averaged");
  if (impl_->continuous_enthalpy) continue_enthalpy(*probe->thermo());
  detail::install_order_rates(probe);
  const auto kerosene = detail::kerosene_material(probe);
  impl_->reaction_model = kerosene ? detail::KeroseneMaterial::model : "mechanism_defined_v1";
  impl_->composition = chemistry::composition(*probe->thermo());
  if (probe->thermo()->type() != "ideal-gas")
    throw std::invalid_argument(
        "portable constant-pressure backend requires ideal gas");
  auto &gas = impl_->gas_identity;
  gas.mechanism_sha256 = impl_->sha256;
  gas.phase = impl_->phase;
  gas.enthalpy_reference = "absolute-standard-formation-298.15K-v1";
  gas.thermodynamic_model = impl_->continuous_enthalpy
      ? "continuous-nasa7-v1" : "source-thermo-v1";
  if (!source_interval) {
    std::uint64_t lower{}, upper{};
    std::memcpy(&lower,&config.minimum_temperature,sizeof(lower));
    std::memcpy(&upper,&config.maximum_temperature,sizeof(upper));
    gas.thermodynamic_model += "/admission-v1:"+std::to_string(lower)+":"+std::to_string(upper);
  }
  gas.element_names = impl_->composition.element_names;
  gas.composition_fingerprint = impl_->composition.fingerprint;
  auto &closure = impl_->closure_identity;
  closure.element_count = gas.element_names.size();
  probe->thermo()->setState_TP(298.15, Cantera::OneAtm);
  std::vector<double> formation(probe->thermo()->nSpecies());
  probe->thermo()->getPartialMolarEnthalpies(formation.data());
  for (std::size_t i = 0; i < impl_->composition.species.size(); ++i) {
    const auto &species = impl_->composition.species[i];
    gas.species_names.push_back(species.name);
    gas.molecular_weights_kg_per_kmol.push_back(
        species.molecular_weight_kg_per_kmol);
    combustion::ChemicalSpeciesIdentity mapped;
    mapped.molecular_weight_kg_per_kmol = species.molecular_weight_kg_per_kmol;
    mapped.formation_enthalpy_j_per_kg =
        formation[i] / species.molecular_weight_kg_per_kmol;
    if (!std::isfinite(mapped.formation_enthalpy_j_per_kg))
      throw std::invalid_argument("nonfinite reference enthalpy");
    for (auto count : species.element_counts) {
      if (count < 0 || count > std::numeric_limits<std::uint16_t>::max())
        throw std::invalid_argument("closure element count not representable");
      gas.element_counts.push_back(static_cast<std::uint32_t>(count));
      mapped.element_counts.push_back(static_cast<std::uint16_t>(count));
    }
    closure.species.push_back(std::move(mapped));
  }
  closure.fingerprint = combustion::chemistry_identity_fingerprint(closure);
  gas.closure_fingerprint = closure.fingerprint;
  if (!matches(config)) {
    throw std::invalid_argument("reacting mechanism composition mismatch");
  }
}

CanteraBackendRuntime::~CanteraBackendRuntime() noexcept = default;

const CompositionIdentity &CanteraBackendRuntime::composition() const noexcept {
  return impl_->composition;
}

std::string_view CanteraBackendRuntime::mechanism_sha256() const noexcept {
  return impl_->sha256;
}

std::string_view CanteraBackendRuntime::mechanism_phase() const noexcept {
  return impl_->phase;
}

std::string_view CanteraBackendRuntime::reaction_model() const noexcept {
  return impl_->reaction_model;
}

bool CanteraBackendRuntime::matches(const CanteraBackendConfig &config) const {
  if (config.mechanism.file != impl_->mechanism ||
      config.mechanism.sha256 != impl_->sha256 ||
      config.mechanism.phase != impl_->phase ||
      config.continuous_enthalpy != impl_->continuous_enthalpy ||
      config.minimum_temperature != impl_->minimum_temperature ||
      config.maximum_temperature != impl_->maximum_temperature ||
      config.species_names.size() != impl_->composition.species.size()) {
    return false;
  }
  for (std::size_t index = 0; index < config.species_names.size(); ++index) {
    if (config.species_names[index] != impl_->composition.species[index].name) {
      return false;
    }
  }
  return true;
}

CanteraWorkspacePool::CanteraWorkspacePool(
    std::shared_ptr<const CanteraBackendRuntime> runtime,
    std::size_t workspace_count)
    : impl_(std::make_unique<Impl>()) {
  if (!runtime || workspace_count == 0U) {
    throw std::invalid_argument(
        "Cantera workspace pool requires runtime and positive size");
  }
  impl_->runtime = runtime;
  impl_->workspaces.reserve(workspace_count);
  for (std::size_t index = 0; index < workspace_count; ++index) {
    impl_->workspaces.push_back(
        make_workspace(runtime->impl_->mechanism, runtime->impl_->phase,
                       runtime->impl_->continuous_enthalpy));
  }
}

CanteraWorkspacePool::~CanteraWorkspacePool() noexcept = default;

std::size_t CanteraWorkspacePool::workspace_count() const noexcept {
  return impl_->workspaces.size();
}

bool CanteraWorkspacePool::workspaces_are_distinct() const noexcept {
  for (std::size_t left = 0; left < impl_->workspaces.size(); ++left) {
    for (std::size_t right = left + 1U; right < impl_->workspaces.size();
         ++right) {
      const auto &a = impl_->workspaces[left];
      const auto &b = impl_->workspaces[right];
      if (a.solution.get() == b.solution.get() ||
          a.thermo.get() == b.thermo.get() ||
          a.kinetics.get() == b.kinetics.get() ||
          a.transport.get() == b.transport.get() ||
          a.reactor.get() == b.reactor.get() ||
          a.network.get() == b.network.get() ||
          (a.kerosene && a.kerosene.get() == b.kerosene.get())) {
        return false;
      }
    }
  }
  return true;
}

void validate_point(const ThermochemicalPoint &point,
                    std::size_t species_count) {
  if (!std::isfinite(point.p0_pa) || point.p0_pa <= 0.0 ||
      !std::isfinite(point.h_tc_j_per_kg) ||
      point.mass_fractions.size() != species_count) {
    throw std::invalid_argument("invalid Cantera thermochemical point");
  }
  double sum = 0.0;
  for (const double fraction : point.mass_fractions) {
    if (!std::isfinite(fraction) || fraction < 0.0) {
      throw std::invalid_argument("invalid Cantera mass fraction");
    }
    sum += fraction;
  }
  const double tolerance = 1.0e-12 * static_cast<double>(species_count);
  if (!std::isfinite(sum) || std::abs(sum - 1.0) > tolerance) {
    throw std::invalid_argument("Cantera mass fractions must sum to one");
  }
}

CanteraBackend::CanteraBackend(std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}

CanteraBackend::~CanteraBackend() = default;

std::size_t CanteraBackend::lane_index() const noexcept { return impl_->lane; }

const CompositionIdentity &CanteraBackend::composition() const noexcept {
  return impl_->runtime->composition();
}

const portable::GasIdentity &CanteraBackend::gas_identity() const noexcept {
  return impl_->runtime->impl_->gas_identity;
}
const combustion::ChemistryIdentity &
CanteraBackend::closure_identity() const noexcept {
  return impl_->runtime->impl_->closure_identity;
}
portable::Status
CanteraBackend::query_gas(const portable::GasQuery &q,
                          portable::GasQueryOutput &out) noexcept {
  out.sample = {};
  if (impl_->pool_lifetime.expired())
    return portable::Status::unavailable;
  const auto n = composition().species.size();
  if (q.composition_fingerprint != composition().fingerprint)
    return portable::Status::identity_mismatch;
  if (out.capacity < n)
    return portable::Status::capacity_exceeded;
  if (!q.mass_fractions || q.species_count != n ||
      !out.diffusivities_m2_per_s || !out.species_enthalpies_j_per_kg ||
      !out.net_mass_rates_kg_per_m3_s || q.revision.algorithm_version != 1 ||
      !std::isfinite(q.pressure_pa) || q.pressure_pa <= 0)
    return portable::Status::invalid_input;
  double sum = 0;
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(q.mass_fractions[i]) || q.mass_fractions[i] < 0 ||
        q.mass_fractions[i] > 1)
      return portable::Status::invalid_input;
    sum += q.mass_fractions[i];
  }
  if (std::abs(sum - 1) > 1e-12 * n)
    return portable::Status::invalid_input;
  auto &w = impl_->pool->workspaces[impl_->lane];
  auto &thermo = *w.thermo;
  const auto& config=*impl_->runtime->impl_;
  const double minimum=config.minimum_temperature>0?config.minimum_temperature:thermo.minTemp();
  const double maximum=config.maximum_temperature>0?config.maximum_temperature:thermo.maxTemp();
  try {
    thermo.setMassFractions_NoNorm(q.mass_fractions);
    if (q.coordinates == portable::GasStateCoordinates::pressure_temperature) {
      if (!std::isfinite(q.temperature_k) ||
          q.temperature_k < minimum ||
          q.temperature_k > maximum)
        return portable::Status::invalid_input;
      thermo.setState_TP(q.temperature_k, q.pressure_pa);
    } else if (q.coordinates ==
               portable::GasStateCoordinates::pressure_enthalpy) {
      if (!std::isfinite(q.enthalpy_j_per_kg))
        return portable::Status::invalid_input;
      // A cell query starts from its own hint (or a fixed reference), rather
      // than the previous cell's reactor endpoint. PH remains authoritative.
      const double seed = std::isfinite(q.temperature_k) &&
          q.temperature_k >= minimum && q.temperature_k <= maximum
          ? q.temperature_k : std::clamp(298.15, minimum, maximum);
      thermo.setState_TP(seed, q.pressure_pa);
      // Resolve PH more tightly than the public query conservation gate.
      solve_pressure_enthalpy(thermo, q.enthalpy_j_per_kg, q.pressure_pa);
    } else
      return portable::Status::invalid_input;
    portable::GasSample sample{q.revision,
                               composition().fingerprint,
                               thermo.pressure(),
                               thermo.temperature(),
                               thermo.density(),
                               thermo.enthalpy_mass(),
                               thermo.cp_mass(),
                               w.transport->viscosity(),
                               w.transport->thermalConductivity()};
    if (!std::isfinite(sample.temperature_k) ||
        sample.temperature_k < minimum ||
        sample.temperature_k > maximum)
      return portable::Status::invalid_input;
    if (!std::isfinite(sample.density_kg_per_m3) ||
        sample.density_kg_per_m3 <= 0 ||
        !std::isfinite(sample.enthalpy_j_per_kg) ||
        !std::isfinite(sample.cp_j_per_kg_k) || sample.cp_j_per_kg_k <= 0 ||
        !std::isfinite(sample.viscosity_pa_s) || sample.viscosity_pa_s <= 0 ||
        !std::isfinite(sample.conductivity_w_per_m_k) ||
        sample.conductivity_w_per_m_k <= 0)
      return portable::Status::provider_failure;
    if (std::abs(sample.pressure_pa - q.pressure_pa) >
            1e-12 * std::max(1., q.pressure_pa) ||
        (q.coordinates == portable::GasStateCoordinates::pressure_enthalpy &&
         std::abs(sample.enthalpy_j_per_kg - q.enthalpy_j_per_kg) >
             1e-10 * std::max(1., std::abs(q.enthalpy_j_per_kg))))
      return portable::Status::provider_failure;
    w.transport->getMixDiffCoeffs(w.query_diffusion.data());
    thermo.getPartialMolarEnthalpies(w.query_enthalpies.data());
    if (w.kerosene) w.kerosene->molar_rates(thermo, w.query_rates.data());
    else w.kinetics->getNetProductionRates(w.query_rates.data());
    for (std::size_t i = 0; i < n; ++i) {
      const auto mw = composition().species[i].molecular_weight_kg_per_kmol;
      // Cantera molar units are kmol: J/kmol -> J/kg and kmol/m3/s -> kg/m3/s.
      w.query_enthalpies[i] /= mw;
      w.query_rates[i] *= mw;
      // A pure species has no mixture diffusion; Cantera may report zero.
      // Preserve that value. The evaporation film requires its own positive D.
      if (!std::isfinite(w.query_diffusion[i]) || w.query_diffusion[i] < 0 ||
          !std::isfinite(w.query_enthalpies[i]) ||
          !std::isfinite(w.query_rates[i]))
        return portable::Status::provider_failure;
    }
    std::copy(w.query_diffusion.begin(), w.query_diffusion.end(),
              out.diffusivities_m2_per_s);
    std::copy(w.query_enthalpies.begin(), w.query_enthalpies.end(),
              out.species_enthalpies_j_per_kg);
    std::copy(w.query_rates.begin(), w.query_rates.end(),
              out.net_mass_rates_kg_per_m3_s);
    out.sample = sample;
    return portable::Status::success;
  } catch (const std::bad_alloc &) {
    return portable::Status::capacity_exceeded;
  } catch (...) {
    return portable::Status::provider_failure;
  }
}

ThermodynamicProperties
CanteraBackend::evaluate(const ThermochemicalPoint &point) const {
  if (impl_->pool_lifetime.expired())
    throw std::runtime_error("Cantera workspace pool expired");
  validate_point(point, composition().species.size());
  auto &thermo = *impl_->pool->workspaces[impl_->lane].thermo;
  const auto& config=*impl_->runtime->impl_;
  const double minimum=config.minimum_temperature>0?config.minimum_temperature:thermo.minTemp();
  const double maximum=config.maximum_temperature>0?config.maximum_temperature:thermo.maxTemp();
  try {
    thermo.setMassFractions_NoNorm(point.mass_fractions.data());
    thermo.setState_TP(std::clamp(298.15, minimum, maximum), point.p0_pa);
    solve_pressure_enthalpy(thermo, point.h_tc_j_per_kg, point.p0_pa);
  } catch (const std::exception &) {
    throw std::runtime_error("Cantera thermodynamic state inversion failed");
  }
  ThermodynamicProperties result{thermo.temperature(), thermo.density(),
                                 thermo.cp_mass(),
                                 thermo.meanMolecularWeight()};
  if (!std::isfinite(result.temperature_k) || result.temperature_k < minimum || result.temperature_k > maximum ||
      !std::isfinite(result.density_kg_per_m3) ||
      result.density_kg_per_m3 <= 0.0 || !std::isfinite(result.cp_j_per_kg_k) ||
      result.cp_j_per_kg_k <= 0.0 ||
      !std::isfinite(result.mixture_molecular_weight_kg_per_kmol) ||
      result.mixture_molecular_weight_kg_per_kmol <= 0.0) {
    throw std::runtime_error(
        "Cantera thermodynamic state inversion produced invalid output");
  }
  const double pressure_tolerance =
      1.0e-12 * std::max(1.0, std::abs(point.p0_pa));
  const double enthalpy_tolerance =
      1.0e-10 * std::max(1.0, std::abs(point.h_tc_j_per_kg));
  if (std::abs(thermo.pressure() - point.p0_pa) > pressure_tolerance ||
      std::abs(thermo.enthalpy_mass() - point.h_tc_j_per_kg) >
          enthalpy_tolerance) {
    throw std::runtime_error(
        "Cantera thermodynamic state inversion changed input authority");
  }
  return result;
}

portable::Status
CanteraBackend::advance_gas(const portable::GasAdvanceQuery &r,
                            portable::GasAdvanceOutput &out) noexcept {
  out.final_sample = {};
  out.completed_duration_s = 0;
  out.internal_step_count = 0;
  out.integrated_heat_release_j_per_m3 = 0;
  if (impl_->pool_lifetime.expired())
    return portable::Status::unavailable;
  const auto n = composition().species.size();
  if (out.capacity < n)
    return portable::Status::capacity_exceeded;
  if (!out.final_mass_fractions ||
      !out.integrated_species_density_delta_kg_per_m3 ||
      r.state.coordinates != portable::GasStateCoordinates::pressure_enthalpy ||
      !std::isfinite(r.start_time_s) || r.start_time_s < 0 ||
      !std::isfinite(r.duration_s) || r.duration_s < 0 ||
      !std::isfinite(r.start_time_s + r.duration_s))
    return portable::Status::invalid_input;
  auto &w = impl_->pool->workspaces[impl_->lane];
  // Each retry starts from the same PH/Y authority. Solver tolerances are
  // upper error limits; conservation determines whether an interval publishes.
  unsigned total_steps = 0;
  const auto attempt = [&](double relative_tolerance,
                           double absolute_tolerance) noexcept {
    bool integration_started = false, steps_recorded = false;
    portable::GasQueryOutput sample{{},
                                    w.query_diffusion.data(),
                                    w.query_enthalpies.data(),
                                    w.query_rates.data(),
                                    n};
    auto status = query_gas(r.state, sample);
    if (status != portable::Status::success)
      return status;
    const double initial_rho = sample.sample.density_kg_per_m3;
    std::uint32_t steps = 0;
    try {
      if (r.duration_s > 0) {
        long count = 0;
        if (w.kerosene) {
          integration_started = true;
          w.kerosene->integrate(*w.thermo, r.duration_s, relative_tolerance,
              absolute_tolerance, impl_->controls.maximum_internal_steps - total_steps,
              w.advance_fractions,impl_->controls.molar_reference_controls);
          count = w.kerosene->steps();
        } else {
          // Independent cell intervals share storage, not reactor volume history.
          w.reactor->setInitialVolume(1.0);
          w.reactor->setEnergyEnabled(true);
          w.reactor->syncState();
          w.network->setInitialTime(0.0);
          if(!impl_->controls.molar_reference_controls)
            w.network->setTolerances(relative_tolerance, absolute_tolerance);
          w.network->setMaxSteps(impl_->controls.maximum_internal_steps - total_steps);
          w.network->reinitialize();
          if(impl_->controls.molar_reference_controls && !w.reference_tolerances_initialized) {
            // COAST integrates Y/W with atol=1e-10 kmol/kg, rtol=0.
            // ReactorNet integrates Y: the equivalent component limit is W*atol.
            // Auxiliary mass/temperature tolerances keep numerical Jacobian
            // perturbations in physical units; PH closes the reported temperature.
            w.reference_tolerances.assign(w.network->neq(),absolute_tolerance);
            w.reference_tolerances[w.reactor->componentIndex("mass")]=absolute_tolerance*initial_rho;
            w.reference_tolerances[w.reactor->componentIndex("temperature")]=1e-3;
            for(std::size_t i=0;i<n;++i)
              w.reference_tolerances[w.reactor->componentIndex(w.thermo->speciesName(i))]=
                  w.thermo->molecularWeight(i)*absolute_tolerance;
            w.network->integrator().setTolerances(0.,w.reference_tolerances.size(),w.reference_tolerances.data());
            w.network->integrator().initialize(0.,*w.network);
            w.reference_tolerances_initialized=true;
          }
          integration_started = true;
          w.network->advance(r.duration_s);
          count = w.network->solverStats()["steps"].asInt();
          auto thermo = w.reactor->phase()->thermo();
          thermo->getMassFractions(w.advance_fractions.data());
        }
        if (count <= 0 || static_cast<unsigned long>(count) >
                              std::numeric_limits<std::uint32_t>::max()) {
          if (count < 0) total_steps = impl_->controls.maximum_internal_steps;
          return portable::Status::provider_failure;
        }
        steps = static_cast<std::uint32_t>(count);
        total_steps += steps;
        steps_recorded = true;
        if (!impl_->controls.molar_reference_controls && !w.kerosene && std::abs(w.thermo->enthalpy_mass() - r.state.enthalpy_j_per_kg) >
            1e-8 * std::max(1., std::abs(r.state.enthalpy_j_per_kg))) {
          return portable::Status::conservation_failure;
        }
        if(impl_->controls.molar_reference_controls && w.thermo->speciesIndex("N2")!=Cantera::npos) {
          const auto dependent=w.thermo->speciesIndex("N2");
          long double sum=0.;
          for(std::size_t i=0;i<n;++i)if(i!=dependent) {
            if(!std::isfinite(w.advance_fractions[i]))return portable::Status::provider_failure;
            w.advance_fractions[i]=std::max(0.,w.advance_fractions[i]);
            sum+=w.advance_fractions[i];
          }
          w.advance_fractions[dependent]=static_cast<double>(1.-sum);
        } else if(!detail::bound_chemistry_roundoff(w.advance_fractions,impl_->controls.absolute_tolerance))
          return portable::Status::provider_failure;
      } else
        std::copy(r.state.mass_fractions, r.state.mass_fractions + n,
                  w.advance_fractions.begin());
      // Publish complete compositions at the same closure tolerance used by
      // ESF. Reintegrate an inaccurate endpoint; retain every elemental
      // inventory.
      const double fraction_sum = std::accumulate(
          w.advance_fractions.begin(), w.advance_fractions.end(), 0.);
      if (!std::isfinite(fraction_sum) || std::abs(fraction_sum - 1) > 2e-12)
        return portable::Status::conservation_failure;
      double mass = 0, scale = 0, heat = 0;
      for (std::size_t i = 0; i < n; ++i) {
        const double y = w.advance_fractions[i];
        if (!std::isfinite(y) || y < 0 || y > 1)
          return portable::Status::provider_failure;
        const double delta = initial_rho * (y - r.state.mass_fractions[i]);
        w.advance_delta[i] = delta;
        mass += delta;
        scale += std::abs(delta);
        heat -=
            closure_identity().species[i].formation_enthalpy_j_per_kg * delta;
      }
      if (!std::isfinite(heat) || std::abs(mass) > 1e-12 + 1e-9 * scale)
        return portable::Status::conservation_failure;
      for (std::size_t e = 0; e < composition().element_names.size(); ++e) {
        double residual = 0, element_scale = 0;
        for (std::size_t i = 0; i < n; ++i) {
          const auto &species = composition().species[i];
          const double term = w.advance_delta[i] * species.element_counts[e] /
                              species.molecular_weight_kg_per_kmol;
          residual += term;
          element_scale += std::abs(term);
        }
        // Reference checkmass changes elemental inventory by its clipping.
        // The common source ledger reports that change; it does not reintegrate.
        if (!std::isfinite(residual) || (!impl_->controls.molar_reference_controls &&
            std::abs(residual) > 1e-12 + 1e-9 * element_scale))
          return portable::Status::conservation_failure;
      }
      auto final_query = r.state;
      final_query.mass_fractions = w.advance_fractions.data();
      status = query_gas(final_query, sample);
      if (status != portable::Status::success)
        return status;
      std::copy(w.advance_fractions.begin(), w.advance_fractions.end(),
                out.final_mass_fractions);
      std::copy(w.advance_delta.begin(), w.advance_delta.end(),
                out.integrated_species_density_delta_kg_per_m3);
      out.final_sample = sample.sample;
      out.completed_duration_s = r.duration_s;
      out.final_sample.enthalpy_j_per_kg = r.state.enthalpy_j_per_kg;
      out.internal_step_count = steps;
      out.integrated_heat_release_j_per_m3 = heat;
      return portable::Status::success;
    } catch (const std::bad_alloc &) {
      return portable::Status::capacity_exceeded;
    } catch (const std::exception& error) {
      std::fprintf(stderr,"chemistry_reference_exception %s\n",error.what());
      if (integration_started && !steps_recorded) {
        try {
          const auto count = w.kerosene ? w.kerosene->steps() :
              w.network->solverStats()["steps"].asInt();
          if (count < 0) total_steps = impl_->controls.maximum_internal_steps;
          else total_steps += static_cast<unsigned>(count);
        } catch (...) {
          total_steps =
              static_cast<unsigned>(impl_->controls.maximum_internal_steps);
        }
      }
      return portable::Status::provider_failure;
    }
  };
  if(impl_->controls.molar_reference_controls) {
    const auto status=attempt(0.,impl_->controls.absolute_tolerance);
    out.internal_step_count=total_steps;
    if(status!=portable::Status::success)std::fprintf(stderr,"chemistry_reference_failure status=%u steps=%u\n",unsigned(status),total_steps);
    return status;
  }
  double relative = impl_->controls.relative_tolerance;
  double absolute = impl_->controls.absolute_tolerance;
  // Keep four bounded refinement opportunities when a caller already uses
  // a tighter absolute tolerance than the legacy 1e-18 refinement floor.
  // Admission retains the original caller tolerance throughout the retries.
  const double absolute_floor = absolute < 1e-18
      ? std::max(std::numeric_limits<double>::min(), absolute * 1e-4) : 1e-18;
  auto result = portable::Status::provider_failure;
  for (unsigned refinement = 0; refinement < 5; ++refinement) {
    if (total_steps >=
        static_cast<unsigned>(impl_->controls.maximum_internal_steps))
      break;
    result = attempt(relative, absolute);
    if (result == portable::Status::success) {
      out.internal_step_count = total_steps;
      return result;
    }
    if (result != portable::Status::conservation_failure &&
        result != portable::Status::provider_failure)
      return result;
    const double next_relative =
        std::min(relative, std::max(1e-13, relative * .1));
    const double next_absolute =
        std::min(absolute, std::max(absolute_floor, absolute * .1));
    if (next_relative == relative && next_absolute == absolute)
      break;
    relative = next_relative;
    absolute = next_absolute;
  }
  return result;
}

TransportProperties
CanteraBackend::evaluate(const ThermochemicalPoint &point,
                         const ThermodynamicProperties &thermodynamics) const {
  const ThermodynamicProperties current = evaluate(point);
  const auto agrees = [](double left, double right) {
    return std::isfinite(left) && std::isfinite(right) &&
           std::abs(left - right) <=
               1.0e-12 * std::max({1.0, std::abs(left), std::abs(right)});
  };
  if (!agrees(current.temperature_k, thermodynamics.temperature_k) ||
      !agrees(current.density_kg_per_m3, thermodynamics.density_kg_per_m3) ||
      !agrees(current.cp_j_per_kg_k, thermodynamics.cp_j_per_kg_k) ||
      !agrees(current.mixture_molecular_weight_kg_per_kmol,
              thermodynamics.mixture_molecular_weight_kg_per_kmol)) {
    throw std::invalid_argument(
        "Cantera transport thermodynamic state mismatch");
  }
  auto &workspace = impl_->pool->workspaces[impl_->lane];
  TransportProperties result;
  try {
    result.viscosity_pa_s = workspace.transport->viscosity();
    result.conductivity_w_per_m_k = workspace.transport->thermalConductivity();
    result.mixture_diffusivity_m2_per_s.resize(workspace.thermo->nSpecies());
    workspace.transport->getMixDiffCoeffs(
        result.mixture_diffusivity_m2_per_s.data());
  } catch (const std::exception &) {
    throw std::runtime_error("Cantera transport evaluation failed");
  }
  const bool invalid_scalar = !std::isfinite(result.viscosity_pa_s) ||
                              result.viscosity_pa_s <= 0.0 ||
                              !std::isfinite(result.conductivity_w_per_m_k) ||
                              result.conductivity_w_per_m_k <= 0.0;
  const bool invalid_diffusivity =
      std::any_of(result.mixture_diffusivity_m2_per_s.begin(),
                  result.mixture_diffusivity_m2_per_s.end(), [](double value) {
                    return !std::isfinite(value) || value < 0.0;
                  });
  if (invalid_scalar || invalid_diffusivity ||
      result.mixture_diffusivity_m2_per_s.size() !=
          composition().species.size()) {
    throw std::runtime_error("Cantera transport produced invalid output");
  }
  return result;
}

ChemistryIntervalReport
CanteraBackend::integrate(const ChemistryIntervalRequest &request) {
  const auto failure = [&](ChemistryStatus status) {
    ChemistryIntervalReport report;
    report.final_state = request.state;
    report.integrated_rho_y_delta_kg_per_m3.assign(
        request.state.mass_fractions.size(), 0.0);
    report.status = status;
    return report;
  };
  if (impl_->pool_lifetime.expired())
    return failure(ChemistryStatus::workspace_failure);
  if (!std::isfinite(request.start_time_s) || request.start_time_s < 0.0 ||
      !std::isfinite(request.duration_s) || request.duration_s < 0.0) {
    return failure(ChemistryStatus::invalid_input);
  }
  try {
    validate_point(request.state, composition().species.size());
  } catch (const std::invalid_argument &) {
    return failure(ChemistryStatus::invalid_input);
  }
  if (request.duration_s == 0.0) {
    auto result = failure(ChemistryStatus::success);
    return result;
  }

  auto &workspace = impl_->pool->workspaces[impl_->lane];
  if (workspace.kerosene) {
    // Both public chemistry interfaces share the PH/Y retry and conservation
    // transaction. Legacy callers keep their existing request/report types.
    auto result = failure(ChemistryStatus::integration_failure);
    portable::GasAdvanceQuery query;
    query.state.revision.algorithm_version = 1;
    query.state.composition_fingerprint = composition().fingerprint;
    query.state.pressure_pa = request.state.p0_pa;
    query.state.enthalpy_j_per_kg = request.state.h_tc_j_per_kg;
    query.state.mass_fractions = request.state.mass_fractions.data();
    query.state.species_count = request.state.mass_fractions.size();
    query.start_time_s = request.start_time_s;
    query.duration_s = request.duration_s;
    portable::GasAdvanceOutput output{{}, result.final_state.mass_fractions.data(),
        result.integrated_rho_y_delta_kg_per_m3.data(), query.state.species_count};
    const auto status = advance_gas(query, output);
    if (status != portable::Status::success) {
      if (status == portable::Status::conservation_failure)
        return failure(ChemistryStatus::conservation_failure);
      if (status == portable::Status::invalid_input)
        return failure(ChemistryStatus::invalid_input);
      return failure(ChemistryStatus::integration_failure);
    }
    result.status = ChemistryStatus::success;
    result.completed_duration_s = output.completed_duration_s;
    result.internal_step_count = output.internal_step_count;
    return result;
  }
  ThermodynamicProperties initial;
  try {
    initial = evaluate(request.state);
  } catch (const std::exception &) {
    return failure(ChemistryStatus::state_inversion_failure);
  }
  try {
    workspace.reactor->setInitialVolume(1.0);
    workspace.reactor->setEnergyEnabled(true);
    workspace.reactor->syncState();
    workspace.network->setInitialTime(0.0);
    workspace.network->setTolerances(impl_->controls.relative_tolerance,
                                     impl_->controls.absolute_tolerance);
    workspace.network->setMaxSteps(impl_->controls.maximum_internal_steps);
    workspace.network->reinitialize();
    workspace.network->advance(request.duration_s);
  } catch (const std::exception &) {
    return failure(ChemistryStatus::integration_failure);
  }

  const auto final_thermo = workspace.reactor->phase()->thermo();
  std::vector<double> final_fractions(final_thermo->nSpecies());
  final_thermo->getMassFractions(final_fractions.data());
  if (!detail::bound_chemistry_roundoff(final_fractions,
                                       impl_->controls.absolute_tolerance))
    return failure(ChemistryStatus::non_finite_output);
  const bool finite = std::all_of(
      final_fractions.begin(), final_fractions.end(),
      [](double value) { return std::isfinite(value) && value >= 0.0; });
  if (!finite || !std::isfinite(final_thermo->temperature()) ||
      final_thermo->temperature() <= 0.0) {
    return failure(ChemistryStatus::non_finite_output);
  }
  const double fraction_sum =
      std::accumulate(final_fractions.begin(), final_fractions.end(), 0.0);
  const double enthalpy_error =
      std::abs(final_thermo->enthalpy_mass() - request.state.h_tc_j_per_kg);
  const double enthalpy_tolerance =
      1.0e-8 * std::max(1.0, std::abs(request.state.h_tc_j_per_kg));
  if (std::abs(fraction_sum - 1.0) >
          1.0e-10 * static_cast<double>(final_fractions.size()) ||
      enthalpy_error > enthalpy_tolerance) {
    return failure(ChemistryStatus::conservation_failure);
  }

  ChemistryIntervalReport result;
  result.final_state = request.state;
  result.final_state.mass_fractions = final_fractions;
  result.integrated_rho_y_delta_kg_per_m3.resize(final_fractions.size());
  for (std::size_t index = 0; index < final_fractions.size(); ++index) {
    result.integrated_rho_y_delta_kg_per_m3[index] =
        initial.density_kg_per_m3 *
        (final_fractions[index] - request.state.mass_fractions[index]);
  }
  result.status = ChemistryStatus::success;
  result.completed_duration_s = request.duration_s;
  try {
    const long int steps = workspace.network->solverStats()["steps"].asInt();
    if (steps > 0 && static_cast<unsigned long>(steps) <=
                         std::numeric_limits<std::uint32_t>::max()) {
      result.internal_step_count = static_cast<std::uint32_t>(steps);
    }
  } catch (const std::exception &) {
    result.internal_step_count = 0U;
  }
  if (result.internal_step_count == 0U) {
    return failure(ChemistryStatus::integration_failure);
  }
  return result;
}

std::unique_ptr<CanteraBackend>
make_cantera_backend(const CanteraBackendConfig &config,
                     CanteraWorkspacePool &pool) {
  // The donor obtained these checks from its case schema. This independent
  // backend has no schema caller, so reject controls before claiming a lane.
  if (!std::isfinite(config.chemistry.relative_tolerance) ||
      (config.chemistry.relative_tolerance < 0.0 ||
       (config.chemistry.relative_tolerance == 0.0 && !config.chemistry.molar_reference_controls)) ||
      !std::isfinite(config.chemistry.absolute_tolerance) ||
      config.chemistry.absolute_tolerance <= 0.0 ||
      config.chemistry.maximum_internal_steps <= 0) {
    throw std::invalid_argument("invalid chemistry solver controls");
  }
  const auto runtime = pool.impl_->runtime.lock();
  if (!runtime) {
    throw std::runtime_error("Cantera backend runtime has expired");
  }
  if (!runtime->matches(config)) {
    throw std::invalid_argument("Cantera backend case identity mismatch");
  }
  if (pool.impl_->next_lane >= pool.impl_->workspaces.size()) {
    throw std::runtime_error("no unclaimed Cantera workspace is available");
  }
  auto impl = std::make_unique<CanteraBackend::Impl>();
  impl->runtime = runtime;
  impl->pool = pool.impl_.get();
  impl->pool_lifetime = pool.impl_->lifetime;
  impl->lane = pool.impl_->next_lane++;
  impl->controls = config.chemistry;
  return std::unique_ptr<CanteraBackend>(new CanteraBackend(std::move(impl)));
}

} // namespace hundun::v04::chemistry

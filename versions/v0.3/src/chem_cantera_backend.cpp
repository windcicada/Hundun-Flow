// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "chem_workspace_detail.hpp"

#include "cantera/base/Solution.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
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

namespace hundun::chemistry {
namespace {

struct Sha256 final {
  std::array<std::uint32_t, 8> state{
      UINT32_C(0x6a09e667), UINT32_C(0xbb67ae85),
      UINT32_C(0x3c6ef372), UINT32_C(0xa54ff53a),
      UINT32_C(0x510e527f), UINT32_C(0x9b05688c),
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
    const std::uint32_t s0 =
        rotate(words[index - 15U], 7U) ^
        rotate(words[index - 15U], 18U) ^ (words[index - 15U] >> 3U);
    const std::uint32_t s1 =
        rotate(words[index - 2U], 17U) ^
        rotate(words[index - 2U], 19U) ^ (words[index - 2U] >> 10U);
    words[index] =
        words[index - 16U] + s0 + words[index - 7U] + s1;
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
    const std::uint32_t s1 =
        rotate(e, 6U) ^ rotate(e, 11U) ^ rotate(e, 25U);
    const std::uint32_t choose = (e & f) ^ (~e & g);
    const std::uint32_t temporary1 =
        h + s1 + choose + kSha256[index] + words[index];
    const std::uint32_t s0 =
        rotate(a, 2U) ^ rotate(a, 13U) ^ rotate(a, 22U);
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
          rounded > static_cast<double>(
                        std::numeric_limits<std::int32_t>::max())) {
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
};

Workspace make_workspace(const std::filesystem::path &mechanism,
                         const std::string &phase) {
  Workspace result;
  result.solution =
      Cantera::newSolution(mechanism.string(), phase, "mixture-averaged");
  result.thermo = result.solution->thermo();
  result.kinetics = result.solution->kinetics();
  result.transport = result.solution->transport();
  result.reactor = std::make_shared<Cantera::IdealGasConstPressureReactor>(
      result.solution, false);
  result.network = std::make_unique<Cantera::ReactorNet>(result.reactor);
  return result;
}

} // namespace

struct CanteraBackendRuntime::Impl final {
  std::filesystem::path mechanism;
  std::string sha256;
  std::string phase;
  CompositionIdentity composition;
};

struct CanteraWorkspacePool::Impl final {
  std::weak_ptr<const CanteraBackendRuntime> runtime;
  std::vector<Workspace> workspaces;
  std::size_t next_lane{};
};

struct CanteraBackend::Impl final {
  std::shared_ptr<const CanteraBackendRuntime> runtime;
  CanteraWorkspacePool::Impl *pool{};
  std::size_t lane{};
  config::ChemistrySolverConfig controls;
};

CanteraBackendRuntime::CanteraBackendRuntime(
    const config::ResolvedReactingCaseV4 &config)
    : impl_(std::make_unique<Impl>()) {
  impl_->mechanism = config.mechanism.file;
  impl_->sha256 = file_sha256(impl_->mechanism);
  if (impl_->sha256 != config.mechanism.sha256) {
    throw std::invalid_argument("reacting mechanism SHA-256 mismatch");
  }
  impl_->phase = config.mechanism.phase;
  auto probe = Cantera::newSolution(impl_->mechanism.string(), impl_->phase,
                                    "mixture-averaged");
  impl_->composition = chemistry::composition(*probe->thermo());
  if (!matches(config)) {
    throw std::invalid_argument("reacting mechanism composition mismatch");
  }
}

CanteraBackendRuntime::~CanteraBackendRuntime() noexcept = default;

const CompositionIdentity &
CanteraBackendRuntime::composition() const noexcept {
  return impl_->composition;
}

std::string_view CanteraBackendRuntime::mechanism_sha256() const noexcept {
  return impl_->sha256;
}

std::string_view CanteraBackendRuntime::mechanism_phase() const noexcept {
  return impl_->phase;
}

bool CanteraBackendRuntime::matches(
    const config::ResolvedReactingCaseV4 &config) const {
  if (config.mechanism.file != impl_->mechanism ||
      config.mechanism.sha256 != impl_->sha256 ||
      config.mechanism.phase != impl_->phase ||
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
        make_workspace(runtime->impl_->mechanism, runtime->impl_->phase));
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
          a.network.get() == b.network.get()) {
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

std::size_t CanteraBackend::lane_index() const noexcept {
  return impl_->lane;
}

const CompositionIdentity &CanteraBackend::composition() const noexcept {
  return impl_->runtime->composition();
}

ThermodynamicProperties
CanteraBackend::evaluate(const ThermochemicalPoint &point) const {
  validate_point(point, composition().species.size());
  auto &thermo = *impl_->pool->workspaces[impl_->lane].thermo;
  try {
    thermo.setMassFractions_NoNorm(point.mass_fractions.data());
    thermo.setState_HP(point.h_tc_j_per_kg, point.p0_pa);
  } catch (const std::exception &) {
    throw std::runtime_error("Cantera thermodynamic state inversion failed");
  }
  ThermodynamicProperties result{
      thermo.temperature(), thermo.density(), thermo.cp_mass(),
      thermo.meanMolecularWeight()};
  if (!std::isfinite(result.temperature_k) || result.temperature_k <= 0.0 ||
      !std::isfinite(result.density_kg_per_m3) ||
      result.density_kg_per_m3 <= 0.0 ||
      !std::isfinite(result.cp_j_per_kg_k) ||
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

TransportProperties CanteraBackend::evaluate(
    const ThermochemicalPoint &point,
    const ThermodynamicProperties &thermodynamics) const {
  const ThermodynamicProperties current = evaluate(point);
  const auto agrees = [](double left, double right) {
    return std::isfinite(left) && std::isfinite(right) &&
           std::abs(left - right) <=
               1.0e-12 *
                   std::max({1.0, std::abs(left), std::abs(right)});
  };
  if (!agrees(current.temperature_k, thermodynamics.temperature_k) ||
      !agrees(current.density_kg_per_m3,
              thermodynamics.density_kg_per_m3) ||
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
    result.conductivity_w_per_m_k =
        workspace.transport->thermalConductivity();
    result.mixture_diffusivity_m2_per_s.resize(
        workspace.thermo->nSpecies());
    workspace.transport->getMixDiffCoeffs(
        result.mixture_diffusivity_m2_per_s.data());
  } catch (const std::exception &) {
    throw std::runtime_error("Cantera transport evaluation failed");
  }
  const bool invalid_scalar =
      !std::isfinite(result.viscosity_pa_s) ||
      result.viscosity_pa_s <= 0.0 ||
      !std::isfinite(result.conductivity_w_per_m_k) ||
      result.conductivity_w_per_m_k <= 0.0;
  const bool invalid_diffusivity = std::any_of(
      result.mixture_diffusivity_m2_per_s.begin(),
      result.mixture_diffusivity_m2_per_s.end(),
      [](double value) { return !std::isfinite(value) || value <= 0.0; });
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
  ThermodynamicProperties initial;
  try {
    initial = evaluate(request.state);
  } catch (const std::exception &) {
    return failure(ChemistryStatus::state_inversion_failure);
  }
  try {
    workspace.reactor->setEnergyEnabled(true);
    workspace.reactor->syncState();
    workspace.network->setInitialTime(request.start_time_s);
    workspace.network->setTolerances(
        impl_->controls.relative_tolerance,
        impl_->controls.absolute_tolerance);
    workspace.network->setMaxSteps(
        impl_->controls.maximum_internal_steps);
    workspace.network->reinitialize();
    workspace.network->advance(request.start_time_s + request.duration_s);
  } catch (const std::exception &) {
    return failure(ChemistryStatus::integration_failure);
  }

  const auto final_thermo = workspace.reactor->phase()->thermo();
  std::vector<double> final_fractions(final_thermo->nSpecies());
  final_thermo->getMassFractions(final_fractions.data());
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
      std::abs(final_thermo->enthalpy_mass() -
               request.state.h_tc_j_per_kg);
  const double enthalpy_tolerance =
      1.0e-8 *
      std::max(1.0, std::abs(request.state.h_tc_j_per_kg));
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
    if (steps > 0 &&
        static_cast<unsigned long>(steps) <=
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
make_cantera_backend(const config::ResolvedReactingCaseV4 &config,
                     CanteraWorkspacePool &pool) {
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
  impl->lane = pool.impl_->next_lane++;
  impl->controls = config.chemistry;
  return std::unique_ptr<CanteraBackend>(
      new CanteraBackend(std::move(impl)));
}

} // namespace hundun::chemistry

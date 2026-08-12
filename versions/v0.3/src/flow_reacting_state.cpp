// SPDX-License-Identifier: Apache-2.0

#include "hundun/flow_reacting_state.hpp"

#include "chem_thermodynamics_service_detail.hpp"

#include "hundun/rt_field_storage.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace hundun::flow {
namespace {

runtime::FieldDescriptor descriptor(std::string name, std::string unit) {
  return {std::move(name),
          std::move(unit),
          "reacting-flow",
          runtime::FunctionSpace::cell_average,
          runtime::ScalarType::float64,
          1U,
          2,
          true,
          runtime::RestartPolicy::persistent,
          runtime::OutputPolicy::selected};
}

std::size_t cell_count(runtime::Int3 extent) {
  if (extent.x <= 0 || extent.y <= 0 || extent.z <= 0) {
    throw std::invalid_argument("reacting state extent must be positive");
  }
  return static_cast<std::size_t>(extent.x) *
         static_cast<std::size_t>(extent.y) *
         static_cast<std::size_t>(extent.z);
}

std::array<int, 3> coordinates(std::size_t flat, runtime::Int3 extent) {
  const std::size_t nx = static_cast<std::size_t>(extent.x);
  const std::size_t ny = static_cast<std::size_t>(extent.y);
  return {static_cast<int>(flat % nx),
          static_cast<int>((flat / nx) % ny),
          static_cast<int>(flat / (nx * ny))};
}

double find_enthalpy(const chemistry::ThermodynamicsService &service,
                     double p0_pa, double target_temperature,
                     const std::vector<double> &mass_fractions) {
  double enthalpy = target_temperature * 1000.0;
  chemistry::ThermodynamicProperties current;
  bool valid = false;
  for (std::size_t attempt = 0; attempt < 32U && !valid; ++attempt) {
    try {
      current = service.evaluate({p0_pa, enthalpy, mass_fractions});
      valid = true;
    } catch (const std::exception &) {
      enthalpy += std::ldexp(1.0e5, static_cast<int>(attempt));
    }
  }
  if (!valid) {
    throw std::invalid_argument("could not bracket reacting initial enthalpy");
  }
  for (std::size_t iteration = 0; iteration < 32U; ++iteration) {
    const double error = current.temperature_k - target_temperature;
    if (std::abs(error) <=
        1.0e-11 * std::max(1.0, target_temperature)) {
      return enthalpy;
    }
    const double increment =
        1.0e-6 * std::max(1.0, std::abs(enthalpy));
    const auto shifted =
        service.evaluate({p0_pa, enthalpy + increment, mass_fractions});
    const double derivative =
        (shifted.temperature_k - current.temperature_k) / increment;
    if (!std::isfinite(derivative) || derivative <= 0.0) {
      throw std::invalid_argument("reacting enthalpy inversion is not monotone");
    }
    enthalpy -= error / derivative;
    current = service.evaluate({p0_pa, enthalpy, mass_fractions});
  }
  throw std::invalid_argument("reacting initial enthalpy did not converge");
}

} // namespace

ReactingFieldIds
declare_reacting_fields(runtime::FieldRegistry &registry,
                        const chemistry::CompositionIdentity &identity) {
  chemistry::validate_composition_identity(identity);
  ReactingFieldIds result;
  result.species_density.reserve(identity.species.size());
  for (const auto &species : identity.species) {
    result.species_density.push_back(registry.declare_field(
        descriptor("reacting.rhoY." + species.name, "kg/m^3")));
  }
  result.total_thermochemical_enthalpy = registry.declare_field(
      descriptor("reacting.rho_h_tc", "J/m^3"));
  return result;
}

struct ReactingFlowState::Impl final {
  const runtime::FieldRegistry *registry{};
  runtime::Int3 extent{};
  ReactingFieldIds fields;
  chemistry::CompositionIdentity composition;
  std::array<std::unique_ptr<runtime::FieldStorage>, 3> layers;
  bool attempt_active{};
  double p0_pa{};
  std::uint64_t writer_epoch{1U};
  mutable std::uint64_t cache_epoch{};
  mutable ReactingLayer cached_layer{ReactingLayer::committed};
  mutable std::size_t cached_cell{};
  mutable std::optional<chemistry::ThermodynamicProperties> cache;
};

ReactingFlowState ReactingFlowState::create(
    const runtime::FieldRegistry &registry, runtime::Int3 extent,
    ReactingFieldIds fields, chemistry::CompositionIdentity composition) {
  chemistry::validate_composition_identity(composition);
  if (!registry.frozen() ||
      fields.species_density.size() != composition.species.size()) {
    throw std::invalid_argument("invalid reacting state registry contract");
  }
  static_cast<void>(cell_count(extent));
  auto impl = std::make_unique<Impl>();
  impl->registry = &registry;
  impl->extent = extent;
  impl->fields = std::move(fields);
  impl->composition = std::move(composition);
  for (auto &layer : impl->layers) {
    layer = std::make_unique<runtime::FieldStorage>(registry, extent);
  }
  return ReactingFlowState(std::move(impl));
}

ReactingFlowState::ReactingFlowState(std::unique_ptr<Impl> impl) noexcept
    : impl_(std::move(impl)) {}
ReactingFlowState::~ReactingFlowState() noexcept = default;
ReactingFlowState::ReactingFlowState(ReactingFlowState &&) noexcept = default;
ReactingFlowState &
ReactingFlowState::operator=(ReactingFlowState &&) noexcept = default;

void ReactingFlowState::initialize_uniform(
    const chemistry::ThermodynamicsService &service, double p0_pa,
    double temperature_k, const std::vector<double> &mass_fractions) {
  if (!std::isfinite(p0_pa) || p0_pa <= 0.0 ||
      !std::isfinite(temperature_k) || temperature_k <= 0.0 ||
      mass_fractions.size() != impl_->fields.species_density.size()) {
    throw std::invalid_argument("invalid reacting initial state");
  }
  const double enthalpy =
      find_enthalpy(service, p0_pa, temperature_k, mass_fractions);
  const auto properties =
      service.evaluate({p0_pa, enthalpy, mass_fractions});
  if (!std::isfinite(properties.density_kg_per_m3) ||
      properties.density_kg_per_m3 <= 0.0) {
    throw std::invalid_argument("invalid reacting initial density");
  }
  impl_->p0_pa = p0_pa;
  const std::size_t cells = cell_count(impl_->extent);
  for (auto &layer : impl_->layers) {
    for (std::size_t species = 0; species < mass_fractions.size(); ++species) {
      auto view =
          layer->view<double>(impl_->fields.species_density[species]);
      for (std::size_t cell = 0; cell < cells; ++cell) {
        const auto c = coordinates(cell, impl_->extent);
        view(c[0], c[1], c[2], 0) =
            properties.density_kg_per_m3 * mass_fractions[species];
      }
    }
    auto enthalpy_view = layer->view<double>(
        impl_->fields.total_thermochemical_enthalpy);
    for (std::size_t cell = 0; cell < cells; ++cell) {
      const auto c = coordinates(cell, impl_->extent);
      enthalpy_view(c[0], c[1], c[2], 0) =
          properties.density_kg_per_m3 * enthalpy;
    }
  }
  ++impl_->writer_epoch;
  impl_->cache.reset();
}

ReactingLayerValues ReactingFlowState::snapshot(ReactingLayer selected) const {
  const auto &storage =
      *impl_->layers[static_cast<std::size_t>(selected)];
  const std::size_t cells = cell_count(impl_->extent);
  ReactingLayerValues result;
  result.rho_y_kg_per_m3.resize(impl_->fields.species_density.size());
  for (std::size_t species = 0; species < result.rho_y_kg_per_m3.size();
       ++species) {
    const auto view =
        storage.view<double>(impl_->fields.species_density[species]);
    result.rho_y_kg_per_m3[species].resize(cells);
    for (std::size_t cell = 0; cell < cells; ++cell) {
      const auto c = coordinates(cell, impl_->extent);
      result.rho_y_kg_per_m3[species][cell] =
          view(c[0], c[1], c[2], 0);
    }
  }
  const auto enthalpy =
      storage.view<double>(impl_->fields.total_thermochemical_enthalpy);
  result.rho_h_tc_j_per_m3.resize(cells);
  for (std::size_t cell = 0; cell < cells; ++cell) {
    const auto c = coordinates(cell, impl_->extent);
    result.rho_h_tc_j_per_m3[cell] = enthalpy(c[0], c[1], c[2], 0);
  }
  return result;
}

void ReactingFlowState::begin_attempt() {
  if (impl_->attempt_active) {
    throw std::logic_error("reacting attempt is already active");
  }
  const auto committed = snapshot(ReactingLayer::committed);
  impl_->attempt_active = true;
  for (std::size_t cell = 0; cell < committed.rho_h_tc_j_per_m3.size();
       ++cell) {
    std::vector<double> values(committed.rho_y_kg_per_m3.size());
    for (std::size_t species = 0; species < values.size(); ++species) {
      values[species] = committed.rho_y_kg_per_m3[species][cell];
    }
    write_trial_cell(cell, values, committed.rho_h_tc_j_per_m3[cell]);
  }
}

void ReactingFlowState::write_trial_cell(
    std::size_t cell, const std::vector<double> &rho_y,
    double rho_h_tc) {
  if (!impl_->attempt_active || cell >= cell_count(impl_->extent) ||
      rho_y.size() != impl_->fields.species_density.size() ||
      !std::isfinite(rho_h_tc)) {
    throw std::invalid_argument("invalid reacting trial write");
  }
  double density = 0.0;
  for (double value : rho_y) {
    if (!std::isfinite(value) || value < 0.0) {
      throw std::invalid_argument("invalid reacting species density");
    }
    density += value;
  }
  if (!std::isfinite(density) || density <= 0.0) {
    throw std::invalid_argument("reacting density must be positive");
  }
  const auto c = coordinates(cell, impl_->extent);
  auto &trial = *impl_->layers[static_cast<std::size_t>(ReactingLayer::trial)];
  for (std::size_t species = 0; species < rho_y.size(); ++species) {
    trial.view<double>(impl_->fields.species_density[species])(
        c[0], c[1], c[2], 0) = rho_y[species];
  }
  trial.view<double>(impl_->fields.total_thermochemical_enthalpy)(
      c[0], c[1], c[2], 0) = rho_h_tc;
  ++impl_->writer_epoch;
  impl_->cache.reset();
}

void ReactingFlowState::rollback_attempt() {
  if (!impl_->attempt_active) {
    throw std::logic_error("reacting attempt is not active");
  }
  impl_->attempt_active = false;
  ++impl_->writer_epoch;
  impl_->cache.reset();
}

void ReactingFlowState::commit_attempt() {
  if (!impl_->attempt_active) {
    throw std::logic_error("reacting attempt is not active");
  }
  const auto committed = snapshot(ReactingLayer::committed);
  const auto trial = snapshot(ReactingLayer::trial);
  auto publish = [&](ReactingLayer layer, const ReactingLayerValues &values) {
    impl_->attempt_active = true;
    auto old_trial = std::move(impl_->layers[2]);
    impl_->layers[2] =
        std::move(impl_->layers[static_cast<std::size_t>(layer)]);
    for (std::size_t cell = 0; cell < values.rho_h_tc_j_per_m3.size(); ++cell) {
      std::vector<double> rho_y(values.rho_y_kg_per_m3.size());
      for (std::size_t species = 0; species < rho_y.size(); ++species) {
        rho_y[species] = values.rho_y_kg_per_m3[species][cell];
      }
      write_trial_cell(cell, rho_y, values.rho_h_tc_j_per_m3[cell]);
    }
    impl_->layers[static_cast<std::size_t>(layer)] =
        std::move(impl_->layers[2]);
    impl_->layers[2] = std::move(old_trial);
  };
  publish(ReactingLayer::history, committed);
  publish(ReactingLayer::committed, trial);
  impl_->attempt_active = false;
  ++impl_->writer_epoch;
  impl_->cache.reset();
}

chemistry::ThermodynamicProperties ReactingFlowState::thermodynamics(
    ReactingLayer layer, std::size_t cell,
    const chemistry::ThermodynamicsService &service) const {
  if (cell >= cell_count(impl_->extent)) {
    throw std::out_of_range("reacting cell is outside state");
  }
  if (impl_->cache && impl_->cache_epoch == impl_->writer_epoch &&
      impl_->cached_layer == layer && impl_->cached_cell == cell) {
    return *impl_->cache;
  }
  const auto values = snapshot(layer);
  double density = 0.0;
  std::vector<double> fractions(values.rho_y_kg_per_m3.size());
  for (std::size_t species = 0; species < fractions.size(); ++species) {
    density += values.rho_y_kg_per_m3[species][cell];
  }
  if (!std::isfinite(density) || density <= 0.0) {
    throw std::runtime_error("reacting cached density is invalid");
  }
  for (std::size_t species = 0; species < fractions.size(); ++species) {
    fractions[species] = values.rho_y_kg_per_m3[species][cell] / density;
  }
  impl_->cache = service.evaluate(
      {impl_->p0_pa, values.rho_h_tc_j_per_m3[cell] / density, fractions});
  impl_->cache_epoch = impl_->writer_epoch;
  impl_->cached_layer = layer;
  impl_->cached_cell = cell;
  return *impl_->cache;
}

std::uint64_t ReactingFlowState::writer_epoch() const noexcept {
  return impl_->writer_epoch;
}
std::uint64_t ReactingFlowState::cache_epoch() const noexcept {
  return impl_->cache_epoch;
}

} // namespace hundun::flow

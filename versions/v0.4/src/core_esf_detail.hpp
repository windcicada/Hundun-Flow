// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "core_reaction_detail.hpp"
#include "core_tcr_history_detail.hpp"
#include "models_esf_detail.hpp"

namespace hundun::v04::detail {
// Persistent fields live in native StateLayers; typed TCR histories are staged
// alongside them.
// the native attempt transaction is the sole acceptance authority.
class ProductEsf {
public:
  Status configure(const ValidatedModel &model,
                   const ProductReactionSources &gas, Int3 cells) {
    if (!model.reaction.esf)
      return {};
    spec_ = *model.reaction.esf;
    cells_ = cells;
    ns_ = gas.gas_identity().species_names.size();
    stride_ = ns_ + 1;
    count_ = std::size_t(cells.x) * cells.y * cells.z;
    c_z_ = model.reaction.mixing_c_z;
    sc_t_ = model.reaction.turbulent_schmidt;
    species_scheme_ = model.schemes.species;
    enthalpy_scheme_ = model.schemes.enthalpy;
    if (!gas.gas_advance() || !gas.gas_query() || ns_ < 2 || ns_ >= UINT8_MAX ||
        model.time.scheme != TimeScheme::backward_euler ||
        model.immersed_boundary || spec_.tcr.mode == TcrMode::validated)
      return invalid();
    for (const auto &boundary : model.boundaries)
      if (boundary.flow_kind != BoundaryKind::periodic)
        return invalid();
    if (!spec_.initial_species_offsets.empty() &&
        spec_.initial_species_offsets.size() != spec_.fields * (ns_ - 1))
      return invalid();
    for (std::size_t s = 0; s + 1 < ns_; ++s) {
      double sum = 0;
      for (std::size_t f = 0; f < spec_.fields; ++f)
        sum += offset(f, s);
      if (std::abs(sum) > 1e-14)
        return invalid();
    }
    const bool tcr = spec_.tcr.mode != TcrMode::off;
    if (tcr) {
      if (spec_.tcr.progress_weights.size() != ns_)
        return invalid();
      for (const auto &name : spec_.tcr.reactants) {
        const auto &names = gas.gas_identity().species_names;
        const auto found = std::find(names.begin(), names.end(), name);
        if (found == names.end())
          return invalid();
        const auto index = std::size_t(found - names.begin());
        if (std::find(reactants_.begin(), reactants_.end(), index) !=
            reactants_.end())
          return invalid();
        reactants_.push_back(index);
      }
    }
    const std::size_t per_cell =
        sizeof(double) * (4 * spec_.fields * stride_ + 4) +
        (tcr ? 2 * (sizeof(tcr::detail::History) +
                    ProductTcrHistory::record_bytes)
             : 0);

    if (count_ > SIZE_MAX / per_cell ||
        count_ * per_cell > model.mesh.limits.max_memory_bytes_per_rank)
      return {StatusCode::allocation_failure, 10215};
    if (tcr)
      tcr_history.configure(gas.fingerprint(), count_,
                            spec_.tcr.initialization_sign);
    workspace_ = std::make_unique<esf::detail::Workspace>(ns_);
    rates_.resize(count_ * spec_.fields * stride_);
    gradients_.resize(3 * rates_.size());
    scratch_.resize(3 * count_);
    mass_divergence_.resize(count_);
    tuple_.resize(spec_.fields * stride_);
    means_.resize(stride_);
    independent_.resize(ns_ - 1);
    diffusion_.resize(ns_);
    enthalpies_.resize(ns_);
    query_rates_.resize(ns_);
    return {};
  }
  Status validate_restart_cell(const std::vector<RestartImageField> &fields,
                               std::size_t start, std::size_t cell,
                               Span<const double> mean_species, double mean_h,
                               double pressure,
                               const ProductReactionSources &gas,
                               const ThermodynamicsPlan &thermo,
                               portable::Revision revision) noexcept {
    if (fields.size() != start + spec_.fields + 1 ||
        mean_species.size != ns_ - 1)
      return invalid();
    std::fill(means_.begin(), means_.end(), 0.0);
    for (std::size_t f = 0; f < spec_.fields; ++f) {
      const auto &field = fields[start + f];
      if (field.role != RestartFieldRole::stochastic_field ||
          field.components != stride_ ||
          field.values.size() != count_ * stride_)
        return invalid();
      const double *row = field.values.data() + cell * stride_;
      portable::GasSample sample;
      auto status = query(gas, thermo, row, pressure, revision, sample);
      if (!status)
        return status;
      for (std::size_t c = 0; c < stride_; ++c)
        means_[c] += row[c] / spec_.fields;
    }
    const auto mapping = gas.species_indices();
    for (std::size_t s = 0; s < mapping.size; ++s)
      if (std::abs(means_[mapping.data[s]] - mean_species.data[s]) > 2e-12)
        return numerical();
    if (std::abs(means_[ns_] - mean_h) > 2e-12 * std::max(1., std::abs(mean_h)))
      return numerical();
    const auto &cache = fields.back();
    if (cache.role != RestartFieldRole::stochastic_transport ||
        cache.components != 2 || cache.values.size() != 2 * count_ ||
        !std::isfinite(cache.values[2 * cell]) ||
        !(cache.values[2 * cell] > 0) ||
        !std::isfinite(cache.values[2 * cell + 1]) ||
        cache.values[2 * cell + 1] < 0)
      return numerical();
    return {};
  }
  bool enabled() const noexcept { return bool(workspace_); }
  ProductTcrHistory tcr_history;
  double offset(std::size_t f, std::size_t s) const noexcept {
    return spec_.initial_species_offsets.empty()
               ? 0
               : spec_.initial_species_offsets[f * (ns_ - 1) + s];
  }
  Status initialize(Span<FieldView> fields,
                    Span<const ConstFieldView> mean_species,
                    ConstFieldView mean_h, const ProductReactionSources &gas,
                    const ThermodynamicsPlan &thermo,
                    double pressure) noexcept {
    if (fields.size != spec_.fields || mean_species.size != ns_ - 1)
      return invalid();
    const auto mapping = gas.species_indices();
    for (int z = 0; z < cells_.z; ++z)
      for (int y = 0; y < cells_.y; ++y)
        for (int x = 0; x < cells_.x; ++x) {
          Int3 cell{x, y, z};
          for (std::size_t f = 0; f < fields.size; ++f) {
            double sum = 0;
            auto *row = tuple_.data() + f * stride_;
            for (std::size_t s = 0; s < mapping.size; ++s) {
              row[mapping.data[s]] =
                  mean_species.data[s].unchecked(cell, 0) + offset(f, s);
              sum += row[mapping.data[s]];
            }
            row[gas.dependent_index()] = 1 - sum;
            row[ns_] = mean_h.unchecked(cell, 0);
            portable::GasSample sample;
            auto status = query(gas, thermo, row, pressure, {0, 1, 1}, sample);
            if (!status)
              return status;
            for (std::size_t c = 0; c < stride_; ++c)
              fields.data[f].unchecked(cell, c) = row[c];
          }
        }
    return {};
  }
  template <class ScalarView>
  Status refresh_transport(const TransportPlan &transport,
                           ConstFieldView temperature, Span<ScalarView> species,
                           ConstFieldView molecular, ConstFieldView effective,
                           ConstFieldView cp, FieldView conductivity,
                           FieldView gamma) noexcept {
    if (species.size != independent_.size())
      return invalid();
    for (int z = 0; z < cells_.z; ++z)
      for (int y = 0; y < cells_.y; ++y)
        for (int x = 0; x < cells_.x; ++x) {
          const Int3 cell{x, y, z};
          for (std::size_t s = 0; s < species.size; ++s)
            independent_[s] = species.data[s].unchecked(cell, 0);
          MolecularTransportState intrinsic;
          auto status = transport.evaluate(
              temperature.unchecked(cell, 0),
              {independent_.data(), independent_.size()}, intrinsic);
          if (!status)
            return status;
          const double turbulent =
              effective.unchecked(cell, 0) - molecular.unchecked(cell, 0);
          const double heat_capacity = cp.unchecked(cell, 0);
          const double coefficient =
              intrinsic.conductivity / heat_capacity + turbulent / sc_t_;
          if (!(heat_capacity > 0) || !std::isfinite(turbulent) ||
              turbulent < 0 || !std::isfinite(coefficient) ||
              !(coefficient > 0))
            return numerical();
          gamma.unchecked(cell, 0) = coefficient;
          conductivity.unchecked(cell, 0) = coefficient * heat_capacity;
        }
    return {};
  }
  Status cache_transport(FieldView cache, ConstFieldView rho, ConstFieldView mu,
                         ConstFieldView mu_eff, ConstFieldView lambda,
                         ConstFieldView cp) noexcept {
    for (int z = 0; z < cells_.z; ++z)
      for (int y = 0; y < cells_.y; ++y)
        for (int x = 0; x < cells_.x; ++x) {
          Int3 cell{x, y, z};
          const double density = rho.unchecked(cell, 0),
                       heat_capacity = cp.unchecked(cell, 0);
          const double turbulent =
              mu_eff.unchecked(cell, 0) - mu.unchecked(cell, 0);
          const double gamma = lambda.unchecked(cell, 0) / heat_capacity;
          if (!(density > 0) || !(heat_capacity > 0) || !std::isfinite(gamma) ||
              !(gamma > 0) || !std::isfinite(turbulent) || turbulent < 0)
            return numerical();
          cache.unchecked(cell, 0) = gamma;
          cache.unchecked(cell, 1) = turbulent / (sc_t_ * density);
        }
    return {};
  }
  Status prepare(const CartesianKernelPlan &kernels,
                 const ProductReactionSources &gas,
                 const ThermodynamicsPlan &thermo,
                 Span<const ConstFieldView> accepted, Span<FieldView> trial,
                 ConstFieldView cache, ConstFieldView rho, ConstFieldView pi,
                 Span<const ConstFieldView> mean_species, ConstFieldView mean_h,
                 double pressure_reference, ConstFaceFluxView flux, double time,
                 double dt, std::uint64_t step, RevisionToken generation,
                 Span<FieldView> sources) noexcept {
    if (accepted.size != spec_.fields || trial.size != spec_.fields ||
        sources.size != ns_ - 1 || mean_species.size != ns_ - 1 || !(dt > 0) ||
        !std::isfinite(dt))
      return invalid();
    auto scratch = scratch_view(rho, 1);
    KernelInvocation call{{}, {&scratch, 1}, {{0, 0, 0}, cells_}, 0, 0,
                          1,  flux.revision};
    auto status = cartesian_face_divergence(kernels, flux, call);
    if (!status)
      return status;
    std::copy(scratch_.begin(), scratch_.begin() + count_,
              mass_divergence_.begin());
    ConstFieldView gamma = cache;
    gamma.components = 1;
    for (std::size_t f = 0; f < spec_.fields; ++f) {
      const auto input = accepted.data[f];
      for (std::uint8_t c = 0; c < stride_; ++c) {
        scratch = scratch_view(rho, 1);
        call.reads = {&input, 1};
        call.read_component_begin = c;
        call.required_face_flux_revision = flux.revision;
        status = cartesian_convection(
            kernels, c == ns_ ? enthalpy_scheme_ : species_scheme_, flux, call);
        if (!status)
          return status;
        for (std::size_t i = 0; i < count_; ++i)
          rates_[slot(i, f, c)] = -scratch_[i];
        call.required_face_flux_revision = 0;
        status = cartesian_diffusion(kernels, gamma, call);
        if (!status)
          return status;
        for (std::size_t i = 0; i < count_; ++i)
          rates_[slot(i, f, c)] += scratch_[i];
        scratch = scratch_view(rho, 3);
        status = cartesian_gradient(kernels, call);
        if (!status)
          return status;
        for (std::size_t i = 0; i < count_; ++i)
          for (std::size_t d = 0; d < 3; ++d)
            gradients_[3 * slot(i, f, c) + d] = scratch_[i + d * count_];
      }
    }
    const auto mapping = gas.species_indices();
    for (int z = 0; z < cells_.z; ++z)
      for (int y = 0; y < cells_.y; ++y)
        for (int x = 0; x < cells_.x; ++x) {
          Int3 cell{x, y, z};
          const auto i = std::size_t(x) + std::size_t(cells_.x) *
                                              (y + std::size_t(cells_.y) * z);
          const double density = rho.unchecked(cell, 0);
          if (!(density > 0) || !std::isfinite(density))
            return numerical();
          for (std::size_t f = 0; f < spec_.fields; ++f)
            for (std::size_t c = 0; c < stride_; ++c) {
              const double value = accepted.data[f].unchecked(cell, c);
              tuple_[f * stride_ + c] = value;
              // Specific transport uses the same committed mass flux, including
              // its continuity correction, rather than a second density solver.
              rates_[slot(i, f, c)] =
                  (rates_[slot(i, f, c)] + value * mass_divergence_[i]) /
                  density;
            }
          const portable::Revision revision{step, generation, 1};
          esf::detail::Request request;
          request.accepted = {revision, gas.chemistry_identity().fingerprint,
                              spec_.fields, ns_, tuple_.data()};
          request.expected_revision = revision;
          request.identity = &gas.chemistry_identity();
          request.dt_s = dt;
          const double delta = std::cbrt(cell_volume(kernels, cell));
          request.mixing_time_s =
              c_z_ * delta * delta / (2 * cache.unchecked(cell, 0) / density);
          request.turbulent_diffusivity_m2_s = cache.unchecked(cell, 1);
          request.deterministic_rates = rates_.data() + slot(i, 0, 0);
          request.gradients = gradients_.data() + 3 * slot(i, 0, 0);
          request.random = {spec_.seed, step, 1, 0, 0, 1};
          if (tcr_history.enabled()) {
            double sum = 0;
            for (std::size_t a = 0; a < mapping.size; ++a) {
              means_[mapping.data[a]] = mean_species.data[a].unchecked(cell, 0);
              sum += means_[mapping.data[a]];
            }
            means_[gas.dependent_index()] = 1 - sum;
            means_[ns_] = mean_h.unchecked(cell, 0);
            const double pressure = pressure_reference + pi.unchecked(cell, 0);
            std::array<double, 4> rates{};
            double psr_rate{};
            status = progress_rate(gas, thermo, means_.data(), pressure,
                                   revision, psr_rate);
            for (std::size_t f = 0; f < spec_.fields && status; ++f)
              status = progress_rate(gas, thermo, tuple_.data() + f * stride_,
                                     pressure, revision, rates[f]);
            if (!status)
              return status;
            const auto &identity = gas.gas_identity();
            const auto mapped =
                tcr::detail::ideal_gas_reactant_mole_fraction_v1(
                    {revision, revision, identity.composition_fingerprint,
                     identity.composition_fingerprint, means_.data(),
                     identity.molecular_weights_kg_per_kmol.data(), ns_,
                     reactants_.data(), reactants_.size(), rates.data(),
                     spec_.fields, psr_rate, spec_.tcr.weak_rate_threshold});
            const auto &history = tcr_history.accepted(i);
            tcr::detail::TrialRequest tcr_request;
            tcr_request.expected_revision = history.revision;
            tcr_request.mapping = mapped;
            tcr_request.mode = spec_.tcr.mode == TcrMode::shadow
                                   ? tcr::detail::Mode::shadow
                                   : tcr::detail::Mode::experimental;
            tcr_request.initialization_sign =
                history.initialized ? 0 : spec_.tcr.initialization_sign;
            const auto candidate = tcr::detail::prepare(history, tcr_request);
            if (!candidate.available)
              return {StatusCode::numerical_failure,
                      10220U + static_cast<std::uint32_t>(candidate.status)};
            if (!tcr_history.stage(i, candidate, step))
              return numerical();
            request.tcr_control = candidate.mixer_control;
          }
          auto moved = workspace_->advance(request);
          if (moved.status != portable::Status::success)
            return numerical();
          std::array<double, 4> pressures{}, densities{};
          for (std::size_t f = 0; f < spec_.fields; ++f) {
            portable::GasSample sample;
            pressures[f] = pressure_reference + pi.unchecked(cell, 0);
            status = query(gas, thermo, moved.candidate.values + f * stride_,
                           pressures[f], revision, sample);
            if (!status)
              return status;
            densities[f] = sample.density_kg_per_m3;
          }
          esf::detail::ReactionRequest reaction{moved.candidate,
                                                revision,
                                                &gas.chemistry_identity(),
                                                &gas.gas_identity(),
                                                pressures.data(),
                                                densities.data(),
                                                time,
                                                dt};
          const auto reacted = workspace_->react(reaction, *gas.gas_advance());
          if (reacted.status != portable::Status::success ||
              !reacted.mean_integrated_species_density_delta_kg_per_m3)
            return numerical();
          for (std::size_t s = 0; s < mapping.size; ++s)
            sources.data[s].unchecked(cell, 0) =
                reacted.mean_integrated_species_density_delta_kg_per_m3
                    [mapping.data[s]] /
                dt;
          for (std::size_t f = 0; f < spec_.fields; ++f)
            for (std::size_t c = 0; c < stride_; ++c)
              trial.data[f].unchecked(cell, c) =
                  reacted.candidate.values[f * stride_ + c];
        }
    tcr_history.seal();
    return {};
  }
  Status reconcile(Span<FieldView> fields, Span<const ConstFieldView> species,
                   ConstFieldView h, ConstFieldView pi, double pressure,
                   const ProductReactionSources &gas,
                   const ThermodynamicsPlan &thermo,
                   portable::Revision revision) noexcept {
    const auto mapping = gas.species_indices();
    for (int z = 0; z < cells_.z; ++z)
      for (int y = 0; y < cells_.y; ++y)
        for (int x = 0; x < cells_.x; ++x) {
          Int3 cell{x, y, z};
          std::fill(means_.begin(), means_.end(), 0);
          for (std::size_t c = 0; c < stride_; ++c)
            for (std::size_t f = 0; f < spec_.fields; ++f)
              means_[c] += fields.data[f].unchecked(cell, c) / spec_.fields;
          double sum = 0;
          for (std::size_t s = 0; s < mapping.size; ++s) {
            const double target = species.data[s].unchecked(cell, 0);
            sum += target;
            means_[mapping.data[s]] -= target;
          }
          means_[gas.dependent_index()] -= 1 - sum;
          means_[ns_] -= h.unchecked(cell, 0);
          for (std::size_t f = 0; f < spec_.fields; ++f) {
            auto *row = tuple_.data() + f * stride_;
            for (std::size_t c = 0; c < stride_; ++c)
              row[c] = fields.data[f].unchecked(cell, c) - means_[c];
            portable::GasSample sample;
            auto status =
                query(gas, thermo, row, pressure + pi.unchecked(cell, 0),
                      revision, sample);
            if (!status)
              return status; // Reject whole attempt; never clip fields.
            for (std::size_t c = 0; c < stride_; ++c)
              fields.data[f].unchecked(cell, c) = row[c];
          }
        }
    return {};
  }

private:
  static Status invalid() noexcept { return {StatusCode::invalid_plan, 10215}; }
  static Status numerical() noexcept {
    return {StatusCode::numerical_failure, 10216};
  }
  std::size_t slot(std::size_t i, std::size_t f, std::size_t c) const noexcept {
    return (i * spec_.fields + f) * stride_ + c;
  }
  FieldView scratch_view(ConstFieldView authority,
                         std::uint8_t components) noexcept {
    return {scratch_.data(),
            cells_,
            {},
            components,
            std::size_t(cells_.x),
            std::size_t(cells_.x) * cells_.y,
            count_,
            0,
            authority.field,
            authority.revision,
            reinterpret_cast<std::uintptr_t>(scratch_.data()),
            authority.revision_domain};
  }
  Status query(const ProductReactionSources &gas,
               const ThermodynamicsPlan &thermo, const double *row,
               double pressure, portable::Revision revision,
               portable::GasSample &sample) noexcept {
    double sum = 0;
    for (std::size_t s = 0; s < ns_; ++s) {
      if (!std::isfinite(row[s]) || row[s] < 0 || row[s] > 1)
        return numerical();
      sum += row[s];
    }
    if (std::abs(sum - 1) > 2e-12 || !std::isfinite(row[ns_]))
      return numerical();
    const auto mapping = gas.species_indices();
    for (std::size_t s = 0; s < mapping.size; ++s)
      independent_[s] = row[mapping.data[s]];
    ThermoState native;
    if (!thermo.evaluate(pressure, row[ns_],
                         {independent_.data(), independent_.size()}, {},
                         native))
      return numerical();
    portable::GasQuery request{revision,
                               gas.gas_identity().composition_fingerprint,
                               portable::GasStateCoordinates::pressure_enthalpy,
                               pressure,
                               row[ns_],
                               0,
                               row,
                               ns_};
    portable::GasQueryOutput output{
        {}, diffusion_.data(), enthalpies_.data(), query_rates_.data(), ns_};
    if (!portable::same_gas_identity(gas.gas_identity(),
                                     gas.gas_query()->gas_identity()) ||
        gas.gas_query()->query_gas(request, output) !=
            portable::Status::success ||
        output.diffusivities_m2_per_s != diffusion_.data() ||
        output.species_enthalpies_j_per_kg != enthalpies_.data() ||
        output.net_mass_rates_kg_per_m3_s != query_rates_.data() ||
        output.capacity != ns_)
      return numerical();
    const auto close = [](double a, double b) {
      return std::isfinite(a) && std::isfinite(b) &&
             std::abs(a - b) <= 1e-8 * std::max({1., std::abs(a), std::abs(b)});
    };
    const auto &v = output.sample;
    if (v.revision != revision ||
        v.composition_fingerprint != request.composition_fingerprint ||
        !close(v.pressure_pa, pressure) ||
        !close(v.enthalpy_j_per_kg, row[ns_]) ||
        !close(v.temperature_k, native.temperature) ||
        !close(v.density_kg_per_m3, native.rho) ||
        !close(v.cp_j_per_kg_k, native.cp) || !(v.cp_j_per_kg_k > 0))
      return numerical();
    sample = v;
    return {};
  }
  Status progress_rate(const ProductReactionSources &gas,
                       const ThermodynamicsPlan &thermo, const double *row,
                       double pressure, portable::Revision revision,
                       double &rate) noexcept {
    portable::GasSample sample;
    auto status = query(gas, thermo, row, pressure, revision, sample);
    if (!status)
      return status;
    double total = 0, magnitude = 0;
    rate = 0;
    for (std::size_t a = 0; a < ns_; ++a) {
      const double value = query_rates_[a];
      if (!std::isfinite(value))
        return numerical();
      total += value;
      magnitude += std::abs(value);
      rate += spec_.tcr.progress_weights[a] * value / sample.density_kg_per_m3;
    }
    if (!std::isfinite(rate) || !std::isfinite(magnitude) ||
        std::abs(total) > 1e-12 + 1e-10 * magnitude)
      return numerical();
    const auto &identity = gas.gas_identity();
    for (std::size_t e = 0; e < identity.element_names.size(); ++e) {
      total = magnitude = 0;
      for (std::size_t a = 0; a < ns_; ++a) {
        const double value =
            query_rates_[a] *
            identity.element_counts[a * identity.element_names.size() + e] /
            identity.molecular_weights_kg_per_kmol[a];
        total += value;
        magnitude += std::abs(value);
      }
      if (!std::isfinite(magnitude) ||
          std::abs(total) > 1e-12 + 1e-10 * magnitude)
        return numerical();
    }
    return {};
  }
  std::vector<std::size_t> reactants_;
  EsfSpec spec_;
  Int3 cells_{};
  std::size_t ns_{}, stride_{}, count_{};
  double c_z_{}, sc_t_{};
  ConvectionScheme species_scheme_{}, enthalpy_scheme_{};
  std::unique_ptr<esf::detail::Workspace> workspace_;
  std::vector<double> rates_, gradients_, scratch_, mass_divergence_, tuple_,
      means_;
  std::vector<double> independent_, diffusion_, enthalpies_, query_rates_;
};
} // namespace hundun::v04::detail

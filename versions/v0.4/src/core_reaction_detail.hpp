// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "app_reaction_detail.hpp"
#include "hundun/v04_product.hpp"
#include "models_chemistry_adapter_detail.hpp"
#include "solver_cartesian_detail.hpp"
#if defined(HUNDUN_V04_REACTING_CANTERA)
#include "hundun/v04_cantera.hpp"
#endif
#include <algorithm>
#include <cmath>
#include <cstring>
#include <memory>
#include <vector>

namespace hundun::v04::detail {
// The finite-rate derivative joins the existing conservative BDF/EX2 rate
// history. It is not a second flow driver or an implicit chemistry integrator.
class ProductReactionSources {
public:
  Status configure(const ValidatedModel &model,
                   ProductCouplingBindings bindings,
                   const std::filesystem::path &case_root) {
    if (!valid_reaction_spec(model.reaction))
      return invalid();
    if (model.reaction.mode == ReactionMode::none)
      return bindings.gas_query == nullptr && bindings.gas_advance == nullptr &&
                     bindings.chemistry_identity == nullptr
                 ? Status{}
                 : invalid();
    if (model.reaction.mode != ReactionMode::finite_rate_mean &&
        model.reaction.mode != ReactionMode::pasr_algebraic_v1 &&
        model.reaction.mode != ReactionMode::esf_tpdf)
      return invalid();
    mode_ = model.reaction.mode;
    mixing_c_z_ = model.reaction.mixing_c_z;
    turbulent_schmidt_ = model.reaction.turbulent_schmidt;
    if (bindings.gas_query == nullptr) {
      const auto &r = model.reaction;
      if (r.representation == ReactionSpec::Representation::analytic_isomer) {
        const bool reversed = model.thermophysics.species.size() == 2 &&
                              model.thermophysics.species[0].stable_name == "B";
        auto backend =
            std::make_unique<chemistry::detail::AnalyticIsomerBackend>(
                r.analytic_rate_s, reversed, r.analytic_cp_j_per_kg_k);
        advance_provider_ = backend.get();
        closure_ = backend->closure_identity();
        owned_provider_ = std::move(backend);
      } else if (r.representation ==
                 ReactionSpec::Representation::direct_cantera) {
#if defined(HUNDUN_V04_REACTING_CANTERA)
        chemistry::CanteraBackendConfig config;
        config.mechanism = {case_root / r.mechanism_file, r.mechanism_sha256,
                            r.phase};
        config.chemistry = {r.relative_tolerance, r.absolute_tolerance,
                            int(r.maximum_internal_steps)};
        for (const auto &species : model.thermophysics.species)
          config.species_names.push_back(species.stable_name);
        cantera_runtime_ =
            std::make_shared<chemistry::CanteraBackendRuntime>(config);
        cantera_pool_ = std::make_unique<chemistry::CanteraWorkspacePool>(
            cantera_runtime_, 1U);
        auto backend = chemistry::make_cantera_backend(config, *cantera_pool_);
        advance_provider_ = backend.get();
        closure_ = backend->closure_identity();
        owned_provider_ = std::move(backend);
#else
        (void)case_root;
        return invalid();
#endif
      }
      bindings.gas_query = owned_provider_.get();
    }
    if (bindings.gas_query == nullptr)
      return invalid();
    const auto &gas = bindings.gas_query->gas_identity();
    const auto ns = model.thermophysics.species.size();
    if (ns < 2 || gas.species_names.size() != ns ||
        gas.molecular_weights_kg_per_kmol.size() != ns ||
        gas.element_names.empty() ||
        gas.element_counts.size() != ns * gas.element_names.size() ||
        gas.mechanism_sha256.size() != 64 ||
        gas.mechanism_sha256 != model.reaction.mechanism_sha256 ||
        gas.phase.empty() || gas.phase != model.reaction.phase ||
        gas.enthalpy_reference.empty() || gas.composition_fingerprint == 0 ||
        gas.closure_fingerprint == 0)
      return invalid();
    for (char c : gas.mechanism_sha256)
      if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f')))
        return invalid();
    for (std::size_t i = 0; i < ns; ++i) {
      const auto &spec = model.thermophysics.species[i];
      if (gas.species_names[i] != spec.stable_name ||
          !std::isfinite(gas.molecular_weights_kg_per_kmol[i]) ||
          gas.molecular_weights_kg_per_kmol[i] <= 0 ||
          std::abs(gas.molecular_weights_kg_per_kmol[i] -
                   spec.molecular_weight) > 1e-12 * spec.molecular_weight)
        return invalid();
    }
    if (bindings.gas_advance != nullptr ||
        bindings.chemistry_identity != nullptr) {
      if (bindings.gas_advance == nullptr ||
          bindings.chemistry_identity == nullptr ||
          !portable::same_gas_identity(gas,
                                       bindings.gas_advance->gas_identity()) ||
          bindings.chemistry_identity->fingerprint != gas.closure_fingerprint ||
          combustion::chemistry_identity_fingerprint(
              *bindings.chemistry_identity) != gas.closure_fingerprint)
        return invalid();
      advance_provider_ = bindings.gas_advance;
      closure_ = *bindings.chemistry_identity;
    }
    identity_ = gas;
    for (const auto &s : model.transported_scalars)
      if (s.role == TransportedScalarRole::species) {
        auto it = std::find(gas.species_names.begin(), gas.species_names.end(),
                            s.stable_name);
        if (it == gas.species_names.end())
          return invalid();
        species_.push_back(
            static_cast<std::size_t>(it - gas.species_names.begin()));
      }
    if (species_.size() + 1 != ns)
      return invalid();
    for (std::size_t i = 0; i < ns; ++i)
      if (std::find(species_.begin(), species_.end(), i) == species_.end())
        dependent_ = i;
    independent_.resize(ns - 1);
    y_.resize(ns);
    diffusion_.resize(ns);
    enthalpies_.resize(ns);
    rates_.resize(ns);
    outputs_.resize(ns - 1);
    views_.resize(ns - 1);
    std::uint64_t h = UINT64_C(1469598103934665603);
    const auto integer = [&](std::uint64_t x) {
      h ^= x;
      h *= UINT64_C(1099511628211);
    };
    const auto string = [&](const std::string &s) {
      integer(s.size());
      for (unsigned char c : s)
        integer(c);
    };
    string(gas.mechanism_sha256);
    string(gas.phase);
    string(gas.enthalpy_reference);
    integer(static_cast<unsigned>(model.reaction.mode));
    integer(gas.composition_fingerprint);
    integer(gas.closure_fingerprint);
    for (const auto &n : gas.species_names)
      string(n);
    for (const auto &n : gas.element_names)
      string(n);
    for (auto v : gas.element_counts)
      integer(v);
    for (double v : gas.molecular_weights_kg_per_kmol) {
      std::uint64_t bits;
      std::memcpy(&bits, &v, 8);
      integer(bits);
    }
    if (mode_ == ReactionMode::pasr_algebraic_v1) {
      string("les-scalar-dissipation-v1;reactant-depletion-l1-v1");
      for (double v : {mixing_c_z_, turbulent_schmidt_}) {
        std::uint64_t bits;
        std::memcpy(&bits, &v, 8);
        integer(bits);
      }
    }
    fingerprint_ = h ? h : 1;
    provider_ = bindings.gas_query;
    return {};
  }
  bool enabled() const noexcept { return provider_ != nullptr; }
  bool esf_enabled() const noexcept { return mode_ == ReactionMode::esf_tpdf; }
  portable::GasQueryProvider *gas_query() const noexcept { return provider_; }
  portable::GasAdvanceProvider *gas_advance() const noexcept {
    return advance_provider_;
  }
  const portable::GasIdentity &gas_identity() const noexcept {
    return identity_;
  }
  const combustion::ChemistryIdentity &chemistry_identity() const noexcept {
    return closure_;
  }
  Span<const std::size_t> species_indices() const noexcept {
    return {species_.data(), species_.size()};
  }
  std::size_t dependent_index() const noexcept { return dependent_; }

  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  Status bind(Span<const FieldId> conserved,
              Span<const FieldId> source) noexcept {
    if (!enabled())
      return conserved.size == 0 && source.size == 0 ? Status{} : invalid();
    if (conserved.size != views_.size() || source.size != views_.size())
      return invalid();
    for (std::size_t i = 0; i < views_.size(); ++i) {
      auto &v = views_[i];
      v.conserved_quantity = conserved.data[i];
      v.explicit_source_field = source.data[i];
      v.stage = esf_enabled() ? 2U : 1U;
      v.units.si_exponents = {1, -3, -1, 0, 0, 0, 0};
      v.capability = esf_enabled() ? ContributionCapability::reacting
                                   : ContributionCapability::chemistry;
      v.source_identity = fingerprint_;
    }
    return {};
  }
  Span<const EquationContributionView> contributions() const noexcept {
    return esf_enabled() ? Span<const EquationContributionView>{}
                         : Span<const EquationContributionView>{views_.data(),
                                                                views_.size()};
  }
  Status prepare(const EquationStateView &state,
                 const ThermodynamicsPlan &thermodynamics,
                 const EquationMaterialView &material,
                 const CartesianKernelPlan &kernels, StateLayers &layers,
                 Int3 cells, Span<const std::uint8_t> activity,
                 std::uint64_t step) noexcept {
    if (!enabled() || esf_enabled())
      return {};
    if (!portable::same_gas_identity(identity_, provider_->gas_identity()) ||
        state.independent_species.size != species_.size() ||
        (activity.size != 0 &&
         activity.size != std::size_t(cells.x) * cells.y * cells.z))
      return invalid();
    for (std::size_t i = 0; i < views_.size(); ++i) {
      auto status = layers.revise_runtime(FieldLifetime::persistent_workspace,
                                          views_[i].explicit_source_field);
      if (status)
        status =
            layers.runtime_view(FieldLifetime::persistent_workspace,
                                views_[i].explicit_source_field, outputs_[i]);
      if (!status)
        return status;
      views_[i].explicit_source_density = as_const(outputs_[i]);
    }
    for (int z = 0; z < cells.z; ++z)
      for (int y = 0; y < cells.y; ++y)
        for (int x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          const auto flat = std::size_t(x) + std::size_t(cells.x) *
                                                 (y + std::size_t(cells.y) * z);
          if (activity.size && activity.data[flat] == 0) {
            for (auto &f : outputs_)
              f.unchecked(cell, 0) = 0;
            continue;
          }
          double sum = 0;
          for (std::size_t s = 0; s < species_.size(); ++s) {
            double f =
                state.independent_species.data[s].trial.unchecked(cell, 0);
            if (!std::isfinite(f) || f < 0 || f > 1)
              return numerical();
            independent_[s] = f;
            y_[species_[s]] = f;
            sum += f;
          }
          if (sum > 1)
            return numerical();
          y_[dependent_] = 1 - sum;
          portable::GasQuery q;
          q.revision = {
              step, std::max(RevisionToken{1}, state.enthalpy.trial.revision),
              1};
          q.composition_fingerprint = identity_.composition_fingerprint;
          q.coordinates = portable::GasStateCoordinates::pressure_enthalpy;
          q.pressure_pa = state.pressure_reference +
                          state.pressure_perturbation.trial.unchecked(cell, 0);
          q.enthalpy_j_per_kg = state.enthalpy.trial.unchecked(cell, 0);
          q.mass_fractions = y_.data();
          q.species_count = y_.size();
          portable::GasQueryOutput out{{},
                                       diffusion_.data(),
                                       enthalpies_.data(),
                                       rates_.data(),
                                       rates_.size()};
          if (provider_->query_gas(q, out) != portable::Status::success ||
              out.diffusivities_m2_per_s != diffusion_.data() ||
              out.species_enthalpies_j_per_kg != enthalpies_.data() ||
              out.net_mass_rates_kg_per_m3_s != rates_.data() ||
              out.capacity != rates_.size())
            return numerical();
          const auto close = [](double a, double b) {
            return std::isfinite(a) && std::isfinite(b) &&
                   std::abs(a - b) <=
                       1e-8 * std::max({1., std::abs(a), std::abs(b)});
          };
          double native_h = 0, native_cp = 0, native_r = 0;
          if (!thermodynamics.mixture_enthalpy(
                  state.temperature.trial.unchecked(cell, 0),
                  {independent_.data(), independent_.size()}, native_h,
                  native_cp, native_r))
            return numerical();
          if (out.sample.revision != q.revision ||
              out.sample.composition_fingerprint != q.composition_fingerprint ||
              !close(out.sample.pressure_pa, q.pressure_pa) ||
              !close(out.sample.enthalpy_j_per_kg, q.enthalpy_j_per_kg) ||
              !close(out.sample.temperature_k,
                     state.temperature.trial.unchecked(cell, 0)) ||
              !close(out.sample.density_kg_per_m3,
                     state.density.trial.unchecked(cell, 0)) ||
              !close(out.sample.cp_j_per_kg_k, native_cp) ||
              !close(native_h, q.enthalpy_j_per_kg) ||
              !(out.sample.density_kg_per_m3 > 0) ||
              !(out.sample.cp_j_per_kg_k > 0))
            return numerical();
          double mass = 0, scale = 0;
          for (double v : rates_) {
            if (!std::isfinite(v))
              return numerical();
            mass += v;
            scale += std::abs(v);
          }
          if (std::abs(mass) > 1e-12 + 1e-10 * scale)
            return numerical();
          for (std::size_t e = 0; e < identity_.element_names.size(); ++e) {
            double residual = 0, magnitude = 0;
            for (std::size_t s = 0; s < rates_.size(); ++s) {
              double v =
                  rates_[s] *
                  identity_
                      .element_counts[s * identity_.element_names.size() + e] /
                  identity_.molecular_weights_kg_per_kmol[s];
              residual += v;
              magnitude += std::abs(v);
            }
            if (std::abs(residual) > 1e-12 + 1e-10 * magnitude)
              return numerical();
          }
          double kappa = 1.0;
          if (mode_ == ReactionMode::pasr_algebraic_v1) {
            double diffusion = 0.0, consumed_density = 0.0, consumption = 0.0;
            for (std::size_t s = 0; s < rates_.size(); ++s) {
              if (!std::isfinite(diffusion_[s]) || diffusion_[s] < 0)
                return numerical();
              diffusion += y_[s] * diffusion_[s];
              if (rates_[s] < 0) {
                consumed_density += out.sample.density_kg_per_m3 * y_[s];
                consumption -= rates_[s];
              }
            }
            const double mu = material.molecular_viscosity.unchecked(cell, 0);
            const double mu_eff =
                material.effective_viscosity.unchecked(cell, 0);
            if (!std::isfinite(mu) || !std::isfinite(mu_eff) || mu <= 0 ||
                mu_eff < mu || consumption <= 0)
              return numerical();
            const auto mixing = combustion::evaluate_mixing_time(
                {std::cbrt(detail::cell_volume(kernels, cell)), diffusion,
                 (mu_eff - mu) / out.sample.density_kg_per_m3,
                 turbulent_schmidt_, mixing_c_z_});
            if (!mixing.succeeded())
              return numerical();
            const auto fraction = combustion::evaluate_pasr_reacting_fraction(
                mixing.tau_mix_s, consumed_density / consumption);
            if (fraction.status !=
                combustion::PasrReactingFractionStatus::success)
              return numerical();
            kappa = fraction.kappa;
          }
          for (std::size_t s = 0; s < species_.size(); ++s)
            outputs_[s].unchecked(cell, 0) = kappa * rates_[species_[s]];
        }
    return {};
  }

private:
  static Status invalid() noexcept { return {StatusCode::invalid_plan, 10214}; }
  static Status numerical() noexcept {
    return {StatusCode::numerical_failure, 10214};
  }
#if defined(HUNDUN_V04_REACTING_CANTERA)
  std::shared_ptr<const chemistry::CanteraBackendRuntime> cantera_runtime_;
  std::unique_ptr<chemistry::CanteraWorkspacePool> cantera_pool_;
#endif
  std::unique_ptr<portable::GasQueryProvider> owned_provider_;
  portable::GasQueryProvider *provider_{};
  portable::GasAdvanceProvider *advance_provider_{};
  combustion::ChemistryIdentity closure_;
  portable::GasIdentity identity_;
  PlanFingerprint fingerprint_{};
  ReactionMode mode_{ReactionMode::none};
  double mixing_c_z_{1.0}, turbulent_schmidt_{0.7};
  std::size_t dependent_{};
  std::vector<std::size_t> species_;
  std::vector<double> independent_, y_, diffusion_, enthalpies_, rates_;
  std::vector<FieldView> outputs_;
  std::vector<EquationContributionView> views_;
};
} // namespace hundun::v04::detail

// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "app_reaction_detail.hpp"
#include "core_response.hpp"
#include "core_striped_batch.hpp"
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
// Mean-rate sources enter the conservative species equations and accepted
// nonadvective history; the common driver owns their time-step transaction.
class ProductReactionSources {
public:
  Status configure(const ValidatedModel &model,
                   ProductCouplingBindings bindings,
                   const std::filesystem::path &case_root) {
    if (!valid_reaction_spec(model.reaction))
      return invalid();
    fold_provider_ = bindings.tcr_fold;
    fold_identity_ = fold_provider_ ? fold_provider_->fingerprint() : 0;
    if (fold_provider_ &&
        (model.reaction.mode != ReactionMode::esf_tpdf || !model.reaction.esf ||
         model.reaction.esf->tcr.mode == TcrMode::off || fold_identity_ == 0))
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
        config.continuous_enthalpy = true;
        config.minimum_temperature=model.thermophysics.minimum_temperature;
        config.maximum_temperature=model.thermophysics.maximum_temperature;
        config.mechanism = {case_root / r.mechanism_file, r.mechanism_sha256,
                            r.phase};
        config.chemistry = {r.relative_tolerance, r.absolute_tolerance,
                            int(r.maximum_internal_steps),r.mode==ReactionMode::esf_tpdf};
        if(config.chemistry.molar_reference_controls) {
          config.chemistry.frozen_material_interval=true;
          config.chemistry.relative_tolerance=0.;
          config.chemistry.absolute_tolerance=1e-10;
        }
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
    interval_enabled_ = model.time.scheme == TimeScheme::cn_be &&
        mode_ == ReactionMode::finite_rate_mean && advance_provider_ != nullptr;
    pasr_interval_enabled_ = model.time.scheme == TimeScheme::cn_be &&
        mode_ == ReactionMode::pasr_algebraic_v1;
    if (pasr_interval_enabled_ && advance_provider_ == nullptr) return invalid();
    response_enabled_ = interval_enabled_ && owned_provider_ &&
        model.reaction.representation == ReactionSpec::Representation::direct_cantera;
    // Reserve one tenth of the primitive composition coupling tolerance for
    // lagging an interval response. The gas integrator keeps its own controls.
    response_relative_ = 0.1 * std::min(1e-10, model.solver.cold_stopping
        ? model.solver.cold_stopping->species : 1e-10);
    response_reference_time_ = model.solver.cold_stopping
        ? model.solver.cold_stopping->reference_time : 0.0;
    response_absolute_ = std::min(model.reaction.absolute_tolerance, response_relative_);
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
    final_y_.resize(ns);
    integrated_delta_.resize(ns);
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
    if (!gas.thermodynamic_model.empty()) string(gas.thermodynamic_model);
    integer(static_cast<unsigned>(model.reaction.mode));
    if (interval_enabled_) string("mean-transport-reactor-cnbe-v1");
    if (pasr_interval_enabled_) string("pasr-accepted-interval-source-cnbe-v1");
    if (response_enabled_) {
      string("bounded-interval-response-reference-time-v2");
      for(double value : {response_relative_,response_absolute_,response_reference_time_}) {
        std::uint64_t bits{}; std::memcpy(&bits,&value,sizeof(bits)); integer(bits);
      }
    }
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
    if (mode_ == ReactionMode::pasr_algebraic_v1 ||
        mode_ == ReactionMode::esf_tpdf) {
      string(mode_ == ReactionMode::pasr_algebraic_v1
                 ? "les-scalar-dissipation-v1;reactant-depletion-l1-v1"
                 : "esf-shared-gamma-total-h-v1");
      for (double v : {mixing_c_z_, turbulent_schmidt_}) {
        std::uint64_t bits;
        std::memcpy(&bits, &v, 8);
        integer(bits);
      }
    }
    if (mode_ == ReactionMode::esf_tpdf) {
      const auto &spec = *model.reaction.esf;
      string("esf-mu-mut-iem-v1;transport-full-reactor-v1");
      string("esf-molecular-stochastic-tuple-bounds-v1");
      const bool dual_pressure = (spec.tcr.mode!=TcrMode::experimental || dynamic_tcr_model(spec.tcr.model)) &&
          effective_coupling(model.time.scheme,model.solver.coupling)==CouplingKind::outer_corrected;
      string(dual_pressure ? "esf-dual-physical-mean-field0-pressure-v1;whole-tuple-flux-correction-v1;realized-statistical-transport-ledger-v1;statistical-face-energy-ledger-v1;thermal-frozen-transport-v1"
                           : "esf-bounded-mean-recenter-v1");
      if(spec.tcr.model==TcrModel::cdphyso_dynamic_v1) {
        string("cdphyso-dynamic-v1;species-interval-rates;eta-0.3;cd-cadence-4;fluid-filter;favre-volume-power-v2");
        string(spec.tcr.fuel);
      }
      if(spec.tcr.model==TcrModel::dyn711_v1) {
        string("dyn711-v1;rate-window-8-evaluate-9;cphi-1-2-9;mass-filter;local-sgs-times-v1");
        string("resolved-sgs-k-limit-v1");
        string("physical-mean-psr-ph;single-interval-conservative-remix;terminal-statistics-v1;bilger-cn-v1");
        string(spec.tcr.fuel);string(spec.tcr.mixture_fraction);
        std::uint64_t bits;std::memcpy(&bits,&spec.tcr.oxidizer_oxygen_mass_fraction,8);integer(bits);
      }
      string("esf-global-transport-then-chemistry-v1");
      string("esf-ibm-complete-scalar-flux-v1");
      if ((spec.tcr.mode!=TcrMode::experimental || dynamic_tcr_model(spec.tcr.model))) {
        string("esf-mean-first-joint-implicit-iem-pressure-v2");
        string("esf-frozen-transport-mass-chemical-increment-v1");
      }
      if (model.schemes.species == ConvectionScheme::tvd2 ||
          model.schemes.species == ConvectionScheme::limited_central2)
        string("esf-common-mixture-face-v1;mean-predictor-v1");
      integer(spec.fields);
      integer(spec.seed);
      for (double value : spec.initial_species_offsets) {
        std::uint64_t bits;
        std::memcpy(&bits, &value, sizeof(bits));
        integer(bits);
      }
      integer(static_cast<unsigned>(spec.tcr.mode));
      integer(static_cast<unsigned>(spec.tcr.initialization_sign + 1));
      for (const auto &name : spec.tcr.reactants)
        string(name);
      for (double value : spec.tcr.progress_weights) {
        std::uint64_t bits;
        std::memcpy(&bits, &value, sizeof(bits));
        integer(bits);
      }
      std::uint64_t bits;
      std::memcpy(&bits, &spec.tcr.weak_rate_threshold, sizeof(bits));
      integer(bits);
    }
    if (fold_provider_) {
      string("native-tcr-explicit-fold-evidence-v1");
      integer(fold_identity_);
    }
    if(model.thermophysics.fixed_pressure_pa>0) {
      string("fixed-thermodynamic-pressure-v1");
      std::uint64_t bits{};
      std::memcpy(&bits,&model.thermophysics.fixed_pressure_pa,8);
      integer(bits);
    }
    fingerprint_ = h ? h : 1;
    provider_ = bindings.gas_query;
    return {};
  }
  bool enabled() const noexcept { return provider_ != nullptr; }
  bool interval_enabled() const noexcept { return interval_enabled_; }
  std::string_view reaction_model() const noexcept {
    if (!enabled()) return "none";
#if defined(HUNDUN_V04_REACTING_CANTERA)
    if (cantera_runtime_) return cantera_runtime_->reaction_model();
#endif
    return owned_provider_ ? "analytic_isomer" : "external";
  }
  std::size_t interval_workspace_bytes() const noexcept {
    return sizeof(double)*(interval_h_.capacity()+interval_p_.capacity()+interval_density_.capacity()) + response_.bytes() + batch_.bytes() + sizeof(double)*(batch_input_.capacity()+batch_output_.capacity());
  }
  Status prepare_distribution(MPI_Comm comm,std::size_t cells) noexcept {
    if(!response_enabled_)return {};
    int ranks{};if(MPI_Comm_size(comm,&ranks)!=MPI_SUCCESS)return {StatusCode::mpi_failure,10360};
    if(ranks==1)return {};
    auto status=batch_.prepare(comm,cells,3+y_.size(),10+2*y_.size());
    if(!status)return status;
    try {batch_input_.reserve(cells*(3+y_.size()));batch_output_.reserve(cells*(10+2*y_.size()));}
    catch(const std::bad_alloc&){status={StatusCode::allocation_failure,10360};}
    return StripedBatch::agree(comm,status);
  }
  double interval_response_tolerance() const noexcept { return response_enabled_ ? response_relative_ : 0.0; }
  double interval_response_absolute_tolerance() const noexcept { return response_enabled_ ? response_absolute_ : 0.0; }
  double interval_response_error() const noexcept { return response_error_; }
  std::uint64_t interval_response_reused_cells() const noexcept { return response_reused_; }
  bool mixing_enabled() const noexcept {
    return mode_ == ReactionMode::pasr_algebraic_v1;
  }
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

  const ProductTcrFoldProvider *fold_provider() const noexcept {
    return fold_provider_;
  }
  PlanFingerprint fold_identity() const noexcept { return fold_identity_; }
  PlanFingerprint fingerprint() const noexcept { return fingerprint_; }
  Status bind(Span<const FieldId> conserved,
              Span<const FieldId> source, std::size_t local_cells) noexcept {
    if (!enabled())
      return conserved.size == 0 && source.size == 0 ? Status{} : invalid();
    if (conserved.size != views_.size() || source.size != views_.size())
      return invalid();
    if (interval_enabled_) {
      const auto scalars = 3U + (response_enabled_ ? y_.size() : 0U);
      if(local_cells==0 || scalars>SIZE_MAX/sizeof(double)-1 ||
         local_cells>SIZE_MAX/(scalars*sizeof(double)+(response_enabled_ ? 1U : 0U))) return invalid();
      try {
        interval_h_.resize(local_cells);
        interval_p_.resize(local_cells);
        interval_density_.resize(local_cells);
        if(response_enabled_) {
          auto status=response_.prepare(local_cells,y_.size(),dependent_,response_relative_,response_absolute_);
          if(!status) return status;
        }
      } catch(const std::bad_alloc&) {
        return {StatusCode::allocation_failure,10240};
      }
    }
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
  // The mean BE species solve consumes the integrated stochastic chemistry
  // as a current-step source. Stored EX2 rates keep using contributions(),
  // which excludes this interval source from the transport history.
  Span<const EquationContributionView> coupling_contributions() const noexcept {
    return esf_enabled() && !esf_sources_ready_ ? Span<const EquationContributionView>{}
        : Span<const EquationContributionView>{views_.data(),views_.size()};
  }
  StageId coupling_source_stage() const noexcept {
    return esf_enabled() ? 2U : enabled() ? 1U : 177U;
  }
  Status publish_esf_sources(Span<const FieldView> sources) noexcept {
    if (!esf_enabled() || sources.size!=views_.size() || !sources.data) return invalid();
    for(std::size_t s=0;s<sources.size;++s) {
      const auto& field=sources.data[s];
      if(field.field!=views_[s].explicit_source_field || !field.base ||
          field.components!=1 || field.revision==0) return invalid();
    }
    for(std::size_t s=0;s<sources.size;++s)
      views_[s].explicit_source_density=as_const(sources.data[s]);
    esf_sources_ready_=true;
    return {};
  }
  Status clear_interval(StateLayers& layers, Int3 cells, bool begin_attempt = false) noexcept {
    if(begin_attempt) response_.reset();
    if (!interval_enabled_) return invalid();
    const auto count=std::size_t(cells.x)*cells.y*cells.z;
    if(interval_h_.size()!=count || interval_p_.size()!=count ||
       interval_density_.size()!=count) return invalid();
    for (std::size_t s = 0; s < outputs_.size(); ++s) {
      auto status = layers.revise_runtime(FieldLifetime::persistent_workspace,
                                          views_[s].explicit_source_field);
      if (status) status = layers.runtime_view(FieldLifetime::persistent_workspace,
          views_[s].explicit_source_field, outputs_[s]);
      if (!status) return status;
      views_[s].explicit_source_density = as_const(outputs_[s]);
      for (int z=0;z<cells.z;++z) for(int y=0;y<cells.y;++y) for(int x=0;x<cells.x;++x)
        outputs_[s].unchecked({x,y,z},0)=0;
    }
    return {};
  }

  // Integrate the retained pre-reaction transport state over this same
  // physical interval. Final species storage uses the reacted state, while
  // the final transport audit uses the retained transport state.
  Status advance_transport(const EquationStateView& state,
      const ThermodynamicsPlan& thermo, StateLayers& layers,
      Span<const FieldView> candidate, Int3 cells, Span<const std::uint8_t> activity,
      double start, double dt, RevisionToken step, RevisionToken generation,
      double& maximum_change, std::uint64_t& chemistry_steps, MPI_Comm comm) noexcept {
    maximum_change=0; chemistry_steps=0;
    response_error_=0; response_reused_=0;
    int ranks{};
    if(MPI_Comm_size(comm,&ranks)!=MPI_SUCCESS)return {StatusCode::mpi_failure,10360};
    // Caller-owned providers may have rank-local state; only independently
    // reset built-in Cantera intervals use the distributed route.
    const bool balanced=response_enabled_ && ranks>1;
    auto ready=[&]()->Status {
    if (!interval_enabled_ || candidate.size != species_.size() ||
        !candidate.data || !(dt>0) || !std::isfinite(dt) ||
        !std::isfinite(start) || start<0 || !(start+dt>start) ||
        !portable::same_gas_identity(identity_,advance_provider_->gas_identity()))
      return {StatusCode::invalid_plan,10240};
    for (std::size_t s=0;s<outputs_.size();++s) {
      auto status=layers.revise_runtime(FieldLifetime::persistent_workspace,
                                       views_[s].explicit_source_field);
      if (status) status=layers.runtime_view(FieldLifetime::persistent_workspace,
          views_[s].explicit_source_field,outputs_[s]);
      if (!status) return status;
      views_[s].explicit_source_density=as_const(outputs_[s]);
    }
      return {};
    }();
    if(balanced)ready=StripedBatch::agree(comm,ready);
    if(!ready)return ready;
    auto prepare_cell = [&](Int3 cell, portable::GasQuery& query, ThermoState& transported, double& storage_density)->Status {
      storage_density=state.density.trial.unchecked(cell,0);
      if (!(storage_density>0) || !std::isfinite(storage_density))
        return {StatusCode::numerical_failure,10241};
      // Use the same represented composition as ThermodynamicsPlan and the
      // species solve; an extended sum can reject a valid near-pure trace.
      double sum=0;
      for(std::size_t s=0;s<species_.size();++s) {
        const double value=candidate.data[s].unchecked(cell,0);
        const long double transported=value;
        if(!std::isfinite(transported) || transported<0 || transported>1)
          return {StatusCode::rejected_step,10242};
        independent_[s]=double(transported);
        y_[species_[s]]=independent_[s];
        sum+=independent_[s];
      }
      if(sum>1) return {StatusCode::rejected_step,10242};
      y_[dependent_]=double(1-sum);
      query={{step,generation,1},identity_.composition_fingerprint,
          portable::GasStateCoordinates::pressure_enthalpy,
          state.eos_pressure(state.pressure_reference,state.pressure_perturbation.trial.unchecked(cell,0)),
          state.enthalpy.trial.unchecked(cell,0),0,y_.data(),y_.size()};
      auto status=thermo.evaluate(query.pressure_pa,query.enthalpy_j_per_kg,
          {independent_.data(),independent_.size()},{},transported);
      if(!status) return status;
      query.temperature_k=transported.temperature;
      return {};
    };
    std::size_t response_index{};
    if(balanced) {
      Status ready;
      batch_input_.clear();
      try {
        for(int z=0;z<cells.z && ready;++z)for(int y=0;y<cells.y && ready;++y)for(int x=0;x<cells.x && ready;++x) {
          const auto flat=(std::size_t(z)*cells.y+y)*cells.x+x;
          if(activity.size && activity.data[flat]==0)continue;
          portable::GasQuery q;ThermoState t;double rho{};
          ready=prepare_cell({x,y,z},q,t,rho);if(!ready)break;
          batch_input_.insert(batch_input_.end(),{q.pressure_pa,q.enthalpy_j_per_kg,q.temperature_k});
          batch_input_.insert(batch_input_.end(),y_.begin(),y_.end());
        }
      }catch(const std::bad_alloc&){ready={StatusCode::allocation_failure,10360};}
      ready=StripedBatch::agree(comm,ready);if(!ready)return ready;
      ready=batch_.run(comm,batch_input_,3+y_.size(),10+2*y_.size(),batch_output_,
        [&](const double* in,double* out)->Status {
          portable::GasQuery q{{step,generation,1},identity_.composition_fingerprint,
            portable::GasStateCoordinates::pressure_enthalpy,in[0],in[1],in[2],in+3,y_.size()};
          portable::GasAdvanceOutput o{{},out+10,out+10+y_.size(),y_.size()};
          const auto result=advance_provider_->advance_gas({q,start,dt},o);
          if(result!=portable::Status::success)return {StatusCode::numerical_failure,10250U+std::uint32_t(result)};
          if(o.final_mass_fractions!=out+10 || o.integrated_species_density_delta_kg_per_m3!=out+10+y_.size() ||
             o.capacity!=y_.size() || o.final_sample.revision!=q.revision ||
             o.final_sample.composition_fingerprint!=q.composition_fingerprint)return {StatusCode::numerical_failure,10243};
          out[0]=o.completed_duration_s;out[1]=o.internal_step_count;out[2]=o.integrated_heat_release_j_per_m3;
          const auto& t=o.final_sample;
          out[3]=t.pressure_pa;out[4]=t.temperature_k;out[5]=t.density_kg_per_m3;out[6]=t.enthalpy_j_per_kg;
          out[7]=t.cp_j_per_kg_k;out[8]=t.viscosity_pa_s;out[9]=t.conductivity_w_per_m_k;
          return {};
        });
      if(!ready)return ready;
    }
    for(int z=0;z<cells.z;++z) for(int y=0;y<cells.y;++y) for(int x=0;x<cells.x;++x) {
      const Int3 cell{x,y,z};
      const auto flat=(std::size_t(z)*cells.y+y)*cells.x+x;
      if(activity.size && activity.data[flat]==0) continue;
      portable::GasQuery query;ThermoState transported;double storage_density{};
      auto status=prepare_cell(cell,query,transported,storage_density);if(!status)return status;
      interval_density_[flat]=storage_density;interval_h_[flat]=query.enthalpy_j_per_kg;interval_p_[flat]=query.pressure_pa;
      portable::GasAdvanceOutput output{{},final_y_.data(),integrated_delta_.data(),y_.size()};
      if(balanced) {
        const double* in=batch_output_.data()+response_index++*(10+2*y_.size());
        output.completed_duration_s=in[0];output.internal_step_count=static_cast<std::uint32_t>(in[1]);
        output.integrated_heat_release_j_per_m3=in[2];
        output.final_sample={query.revision,query.composition_fingerprint,in[3],in[4],in[5],in[6],in[7],in[8],in[9]};
        std::copy_n(in+10,y_.size(),final_y_.data());std::copy_n(in+10+y_.size(),y_.size(),integrated_delta_.data());
      } else {
      const auto result=advance_provider_->advance_gas({query,start,dt},output);
      if(result!=portable::Status::success)
        return {StatusCode::numerical_failure,10250U+std::uint32_t(result)};
      }
      const auto close=[](double a,double b) {
        return std::isfinite(a)&&std::isfinite(b)&&
            std::abs(a-b)<=1e-8*std::max({1.,std::abs(a),std::abs(b)});
      };
      if(output.final_mass_fractions!=final_y_.data() ||
          output.integrated_species_density_delta_kg_per_m3!=integrated_delta_.data() ||
          output.capacity!=y_.size() || output.completed_duration_s!=dt ||
          output.final_sample.revision!=query.revision ||
          output.final_sample.composition_fingerprint!=query.composition_fingerprint ||
          !close(output.final_sample.enthalpy_j_per_kg,query.enthalpy_j_per_kg) ||
          !close(output.final_sample.pressure_pa,query.pressure_pa))
        return {StatusCode::numerical_failure,10243};
      const auto check_response = [&](bool provider_delta) -> Status {
        double mass=0,scale=0;
        for(std::size_t j=0;j<y_.size();++j) {
          const double delta=transported.rho*(final_y_[j]-y_[j]);
          if(!std::isfinite(final_y_[j]) || final_y_[j]<0 || final_y_[j]>1 ||
              (provider_delta && !close(delta,integrated_delta_[j])))
            return {StatusCode::numerical_failure,10244};
          const double physical_delta=storage_density*(final_y_[j]-y_[j]);
          mass+=physical_delta; scale+=std::abs(physical_delta);
        }
        if(std::abs(mass)>1e-12+1e-10*scale)
          return {StatusCode::numerical_failure,10245};
        for(std::size_t e=0;e<identity_.element_names.size();++e) {
          double residual=0,magnitude=0;
          for(std::size_t j=0;j<y_.size();++j) {
            const double value=storage_density*(final_y_[j]-y_[j])*
                identity_.element_counts[j*identity_.element_names.size()+e]/
                identity_.molecular_weights_kg_per_kmol[j];
            residual+=value; magnitude+=std::abs(value);
          }
          if(std::abs(residual)>1e-12+1e-10*magnitude)
            return {StatusCode::numerical_failure,10246};
        }
        return {};
      };
      status=check_response(true);
      if(!status) return status;
      if(response_enabled_) {
        IntervalReactionResponse::Report response;
        // The cache stores a dimensionless interval increment. Reference
        // equation residuals measure rates over reference_time, so their
        // permissible increment scales with this physical interval's dt.
        const double interval_scale = response_reference_time_ > 0
            ? std::min(1.0, dt / response_reference_time_) : 1.0;
        status=response_.select(flat,{y_.data(),y_.size()},
            {final_y_.data(),final_y_.size()},response,interval_scale);
        if(!status) return status;
        if(response.reused) {
          status=check_response(false);
          if(!status) return status;
          ++response_reused_;
          response_error_=std::max(response_error_,response.error_ratio);
        }
      }
      for(std::size_t j=0;j<species_.size();++j) {
        maximum_change=std::max(maximum_change,std::abs(
            final_y_[species_[j]]-candidate.data[j].unchecked(cell,0)));
        outputs_[j].unchecked(cell,0)=storage_density*
            (final_y_[species_[j]]-y_[species_[j]])/dt;
        candidate.data[j].unchecked(cell,0)=final_y_[species_[j]];
      }
      chemistry_steps+=output.internal_step_count;
    }
    return {};
  }

  // Chemistry changes mass fractions within the cell's target mass. Pressure
  // corrections update that mass; apply the common density factor to every
  // species source while retaining the integrated specific increment.
  Status reweight_interval(const EquationStateView& state, StateLayers& layers,
      Int3 cells, Span<const std::uint8_t> activity,
      Span<const PrimitiveHistory> transported, Span<const FieldView> reacted,
      double dt) noexcept {
    if (transported.size!=outputs_.size() || reacted.size!=outputs_.size() ||
        !transported.data || !reacted.data || !(dt>0) || !std::isfinite(dt))
      return {StatusCode::invalid_plan,10247};
    for(std::size_t s=0;s<outputs_.size();++s) {
      auto status=layers.revise_runtime(FieldLifetime::persistent_workspace,
          views_[s].explicit_source_field);
      if(status) status=layers.runtime_view(FieldLifetime::persistent_workspace,
          views_[s].explicit_source_field,outputs_[s]);
      if(!status) return status;
      views_[s].explicit_source_density=as_const(outputs_[s]);
    }
    for(int z=0;z<cells.z;++z) for(int y=0;y<cells.y;++y) for(int x=0;x<cells.x;++x) {
      const auto i=(std::size_t(z)*cells.y+y)*cells.x+x;
      if(activity.size && activity.data[i]==0) continue;
      const Int3 c{x,y,z};
      const double density=state.density.trial.unchecked(c,0);
      const double ratio=density/interval_density_[i];
      if(!(ratio>0) || !std::isfinite(ratio))
        return {StatusCode::numerical_failure,10247};
      for(std::size_t j=0;j<outputs_.size();++j) {
        const double increment=reacted.data[j].unchecked(c,0)-
            transported.data[j].trial.unchecked(c,0);
        const double original=interval_density_[i]*increment/dt;
        if(outputs_[j].unchecked(c,0)!=original)
          return {StatusCode::numerical_failure,10248};
        outputs_[j].unchecked(c,0)=density*increment/dt;
      }
      interval_density_[i]=density;
    }
    return {};
  }

  double interval_input_residual(const EquationStateView& state,
      ConstFieldView cp, Int3 cells, Span<const std::uint8_t> activity) const noexcept {
    double maximum=0;
    for(int z=0;z<cells.z;++z) for(int y=0;y<cells.y;++y) for(int x=0;x<cells.x;++x) {
      const auto i=(std::size_t(z)*cells.y+y)*cells.x+x;
      if(activity.size && activity.data[i]==0) continue;
      const Int3 c{x,y,z};
      const double pressure=state.eos_pressure(state.pressure_reference,
          state.pressure_perturbation.trial.unchecked(c,0));
      const double h_scale=std::max(1.,cp.unchecked(c,0)*state.temperature.trial.unchecked(c,0));
      const double error=std::max(
          std::abs(state.enthalpy.trial.unchecked(c,0)-interval_h_[i])/h_scale,
          std::abs(pressure-interval_p_[i])/std::max(1.,pressure));
      if(!std::isfinite(error)) return std::numeric_limits<double>::infinity();
      maximum=std::max(maximum,error);
    }
    return maximum;
  }

  Status prepare(const EquationStateView &state,
                 const ThermodynamicsPlan &thermodynamics,
                 const EquationMaterialView &material,
                 const CartesianKernelPlan &kernels, StateLayers &layers,
                 Int3 cells, Span<const std::uint8_t> activity,
                 std::uint64_t step, double start=0., double dt=0.) noexcept {
    if (!enabled() || esf_enabled())
      return {};
    if (!std::isfinite(dt) || dt<0 || !std::isfinite(start) || start<0 ||
        (dt>0 && !(start+dt>start))) return invalid();
    if (pasr_interval_enabled_ && dt>0 &&
        !portable::same_gas_identity(identity_,advance_provider_->gas_identity()))
      return invalid();
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
          q.pressure_pa = state.eos_pressure(state.pressure_reference,
                          state.pressure_perturbation.trial.unchecked(cell, 0));
          q.enthalpy_j_per_kg = state.enthalpy.trial.unchecked(cell, 0);
          q.temperature_k = state.temperature.trial.unchecked(cell, 0);
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
                mu_eff < mu)
              return numerical();
            // Zero net chemistry preserves its state, as in REFERENCE's inactive
            // progress branch. A chemical timescale is needed for active rates.
            if (scale == 0.0) {
              for (auto &output : outputs_) output.unchecked(cell, 0) = 0.0;
              continue;
            }
            if (consumption <= 0) return numerical();
            const auto fraction = combustion::evaluate_pasr_reacting_fraction(
                {std::cbrt(detail::cell_volume(kernels, cell)), diffusion,
                 (mu_eff - mu) / out.sample.density_kg_per_m3,
                 turbulent_schmidt_, mixing_c_z_}, consumed_density / consumption);
            if (fraction.status !=
                combustion::PasrReactingFractionStatus::success)
              return numerical();
            kappa = fraction.kappa;
          }
          if (pasr_interval_enabled_ && dt>0 && kappa>0) {
            // Integrate once from the accepted state, before fluid coupling.
            // A convex PaSR increment cannot consume more than that state's
            // inventory. Freezing an instantaneous stiff sink has no such
            // property, even when kappa is a valid timescale fraction.
            portable::GasAdvanceOutput advanced{{},final_y_.data(),
                integrated_delta_.data(),y_.size()};
            const auto result=advance_provider_->advance_gas({q,start,dt},advanced);
            if(result!=portable::Status::success)
              return {StatusCode::numerical_failure,10250U+std::uint32_t(result)};
            if(advanced.final_mass_fractions!=final_y_.data() ||
               advanced.integrated_species_density_delta_kg_per_m3!=integrated_delta_.data() ||
               advanced.capacity!=y_.size() || advanced.completed_duration_s!=dt ||
               advanced.final_sample.revision!=q.revision ||
               advanced.final_sample.composition_fingerprint!=q.composition_fingerprint ||
               !close(advanced.final_sample.pressure_pa,q.pressure_pa) ||
               !close(advanced.final_sample.enthalpy_j_per_kg,q.enthalpy_j_per_kg) ||
               !std::isfinite(advanced.integrated_heat_release_j_per_m3))
              return {StatusCode::numerical_failure,10243};
            double total=0., magnitude=0.;
            for(std::size_t s=0;s<y_.size();++s) {
              if(!std::isfinite(final_y_[s]) || final_y_[s]<0 || final_y_[s]>1)
                return {StatusCode::numerical_failure,10244};
              const double change=final_y_[s]-y_[s];
              if(!close(out.sample.density_kg_per_m3*change,integrated_delta_[s]))
                return {StatusCode::numerical_failure,10244};
              const double delta=state.density.trial.unchecked(cell,0)*change;
              rates_[s]=delta/dt;
              if(!std::isfinite(rates_[s])) return numerical();
              total+=delta; magnitude+=std::abs(delta);
            }
            if(std::abs(total)>1e-12+1e-10*magnitude)
              return {StatusCode::numerical_failure,10245};
            for(std::size_t e=0;e<identity_.element_names.size();++e) {
              double defect=0., magnitude=0.;
              for(std::size_t s=0;s<y_.size();++s) {
                const double amount=(final_y_[s]-y_[s])*
                    identity_.element_counts[s*identity_.element_names.size()+e]/
                    identity_.molecular_weights_kg_per_kmol[s];
                defect+=amount; magnitude+=std::abs(amount);
              }
              if(std::abs(defect)>1e-12+1e-10*magnitude)
                return {StatusCode::numerical_failure,10246};
            }
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
  const ProductTcrFoldProvider *fold_provider_{};
  PlanFingerprint fold_identity_{};
  ReactionMode mode_{ReactionMode::none};
  bool interval_enabled_{};
  bool pasr_interval_enabled_{};
  bool esf_sources_ready_{};
  bool response_enabled_{};
  double response_reference_time_{};
  double response_relative_{}, response_absolute_{}, response_error_{};
  std::uint64_t response_reused_{};
  IntervalReactionResponse response_;
  double mixing_c_z_{1.0}, turbulent_schmidt_{0.7};
  std::size_t dependent_{};
  std::vector<std::size_t> species_;
  std::vector<double> independent_, y_, diffusion_, enthalpies_, rates_;
  std::vector<double> final_y_, integrated_delta_, interval_h_, interval_p_, interval_density_;
  std::vector<FieldView> outputs_;
  std::vector<EquationContributionView> views_;
  StripedBatch batch_;
  std::vector<double> batch_input_,batch_output_;
};
} // namespace hundun::v04::detail

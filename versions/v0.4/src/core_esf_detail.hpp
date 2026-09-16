// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "core_reaction_detail.hpp"
#include "core_tcr_history_detail.hpp"
#include "core_tcr_dynamic_detail.hpp"
#include "hundun/v04_ibm.hpp"
#include "models_esf_detail.hpp"
#include "models_exchange_batch_detail.hpp"
#include "solver_heat_boundary_detail.hpp"
#include "solver_ibm_scalar_transport_detail.hpp"
#include "solver_cold.hpp"
#include "solver_mass_source_detail.hpp"
#include "field_view_interval_detail.hpp"

namespace hundun::v04::detail {
// Persistent fields live in native StateLayers; typed TCR histories are staged
// alongside them.
// the native attempt transaction is the sole acceptance authority.
class ProductEsf {
public:
  // Gamma, turbulent scalar diffusivity, molecular viscosity, eddy viscosity.
  static constexpr std::uint8_t transport_components=4;
  static constexpr StageId transport_source_stage=178U;
  Status configure(const ValidatedModel &model,
                   const ProductReactionSources &gas, MeshPatch patch,
                   Int3 global_cells) {
    if (!model.reaction.esf)
      return {};
    const auto cells = patch.cells;
    begin_ = patch.begin;
    global_cells_ = global_cells;
    for(unsigned face=0;face<6;++face)
      physical_boundary_[face]=model.boundaries[face].flow_kind!=BoundaryKind::periodic;
    spec_ = *model.reaction.esf;
    cells_ = cells;
    ns_ = gas.gas_identity().species_names.size();
    gas_fingerprint_=gas.fingerprint();
    stride_ = ns_ + 1;
    count_ = std::size_t(cells.x) * cells.y * cells.z;
    c_z_ = model.reaction.mixing_c_z;
    sc_t_ = model.reaction.turbulent_schmidt;
    for (const auto &scalar : model.transported_scalars)
      if (scalar.role == TransportedScalarRole::passive_scalar)
        passive_schmidt_.push_back(
            {scalar.molecular_schmidt, scalar.turbulent_schmidt});
    species_scheme_ = model.schemes.species;
    enthalpy_scheme_ = model.schemes.enthalpy;
    common_transport_ = species_scheme_ == ConvectionScheme::tvd2 ||
                        species_scheme_ == ConvectionScheme::limited_central2;
    if (!esf::valid_field_count(spec_.fields) || !gas.gas_advance() || !gas.gas_query() || ns_ < 2 || ns_ >= UINT8_MAX ||
        (model.time.scheme != TimeScheme::backward_euler &&
         !(model.time.scheme == TimeScheme::cn_be &&
           (spec_.tcr.mode == TcrMode::off || spec_.tcr.mode == TcrMode::shadow ||
            spec_.tcr.model == TcrModel::cdphyso_dynamic_v1))) ||
        spec_.tcr.mode == TcrMode::validated)
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
    if (tcr && !dynamic_tcr()) {
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
        sizeof(double) * (4 * spec_.fields * stride_ + spec_.fields + 5) +
        (spec_.tcr.mode!=TcrMode::experimental || dynamic_tcr() ? sizeof(ColdPressureRow)+sizeof(double) : 0) +
        (dynamic_tcr() ? 32 * (5*ns_+3) + 48 :
         tcr ? 2 * (sizeof(tcr::detail::History) +
                    ProductTcrHistory::record_bytes)
             : 0);

    const std::size_t thermal_count = common_transport_ ?
        std::size_t(cells.x + 4) * (cells.y + 4) * (cells.z + 4) : 0;
    const std::size_t face_count = common_transport_ ?
        std::size_t(cells.x + 1) * cells.y * cells.z +
        std::size_t(cells.y + 1) * cells.x * cells.z +
        std::size_t(cells.z + 1) * cells.x * cells.y : 0;
    const auto extra_bytes = (thermal_count + face_count +
        2*spec_.fields*stride_ + 6*spec_.fields + 1 + stride_) * sizeof(double) +
        (common_transport_ ? (ns_ - 1) * sizeof(ConstFieldView) : 0);
    const auto maximum_bytes = model.mesh.limits.max_memory_bytes_per_rank;
    if (count_ > SIZE_MAX / per_cell || extra_bytes > maximum_bytes ||
        count_ * per_cell > maximum_bytes - extra_bytes)
      return {StatusCode::allocation_failure, 10215};
    if (common_transport_) {
      const auto allocated = FaceFluxStorage::allocate_workspace(cells, 1, mixture_storage_);
      if (!allocated) return allocated;
      thermal_coordinate_.resize(thermal_count);
      mixture_species_.resize(ns_ - 1);
    }
    if (dynamic_tcr()) {
      const auto &names=gas.gas_identity().species_names;
      const auto find=[&](std::string_view name) {
        return std::size_t(std::find(names.begin(),names.end(),name)-names.begin());
      };
      fuel_=find(spec_.tcr.fuel);
      product_species_=find("H2O");
      if(product_species_==ns_) product_species_=find("CO2");
      dynamic_groups_={product_species_,find("O2"),find("OH")};
      if(fuel_==ns_ || product_species_==ns_ || dynamic_groups_[1]==ns_ ||
          ns_>120 || 2/c_z_<1 || 2/c_z_>16) return invalid();
      mixing_groups_.assign(ns_+1,tcr::detail::MixingGroup::radical);
      for(std::size_t s=0;s<ns_;++s) {
        if(names[s]=="H2O" || names[s]=="CO2")mixing_groups_[s]=tcr::detail::MixingGroup::product;
        else if(names[s]=="O2" || s==fuel_)mixing_groups_[s]=tcr::detail::MixingGroup::reactant;
      }
      mixing_groups_[ns_]=tcr::detail::MixingGroup::product;
      tcr_history.configure_dynamic(gas.fingerprint(),count_,ns_,2/c_z_);
    } else if (tcr)
      tcr_history.configure(gas.fingerprint(), count_,
                            spec_.tcr.initialization_sign);
    workspace_ = std::make_unique<esf::detail::Workspace>(ns_, spec_.fields);
    if (implicit_transport()) {
      transport_rows_.resize(count_);
      transport_pc_=std::make_unique<ColdPressureDilu>(transport_rows_,cells_,LinearIdentity{});
    }
    rates_.resize(count_ * spec_.fields * stride_);
    reactor_density_.resize(count_ * spec_.fields);
    gradients_.resize(3 * rates_.size());
    scratch_.resize(3 * count_);
    mass_divergence_.resize(count_);
    transport_carrier_density_.resize(count_);
    tuple_.resize(spec_.fields * stride_);
    noise_bounds_.resize(2*tuple_.size()+1);
    noise_.resize(stride_);
    wiener_.resize(spec_.fields);
    field_rates_.resize(spec_.fields);
    field_pressures_.resize(spec_.fields);
    field_densities_.resize(spec_.fields);
    means_.resize(stride_);
    independent_.resize(ns_ - 1);
    diffusion_.resize(ns_);
    enthalpies_.resize(ns_);
    query_rates_.resize(ns_);
    return {};
  }
  Status configure_immersed_exchange(MPI_Comm comm,
                                     const CartesianGeometryPlan &geometry,
                                     MeshPatch patch,
                                     const QuadraticStencilPlan &reconstruction,
                                     Span<const RemoteDonorFieldSpec> fields,
                                     std::uint64_t maximum_bytes) {
    auto status = RemoteDonorExchangePlan::analyze(
        comm, geometry.global_cells(), patch, reconstruction, fields, 10,
        immersed_halo_);
    if (!status)
      return status;
    const auto stats = immersed_halo_.stats();
    const auto old_bytes = owned_bytes();
    std::uint64_t bad = 0;
    if (stats.received_cells > UINT64_MAX - stats.supplied_cells ||
        stats.received_cells + stats.supplied_cells > UINT64_MAX / 256 ||
        stats.bytes_per_exchange > UINT64_MAX / 2)
      bad = 1;
    const auto metadata =
        bad ? 0 : (stats.received_cells + stats.supplied_cells) * 256;
    const auto buffers = bad ? 0 : 2 * stats.bytes_per_exchange;
    if (old_bytes > maximum_bytes || metadata > maximum_bytes - old_bytes ||
        buffers > maximum_bytes - old_bytes - metadata)
      bad = 1;
    if (MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_UINT64_T, MPI_MAX, comm) !=
        MPI_SUCCESS)
      return {StatusCode::mpi_failure, 10215};
    if (bad)
      return {StatusCode::allocation_failure, 10215};
    immersed_bytes_ = metadata + buffers;
    status = immersed_halo_.bind(comm);
    if (status)
      immersed_exchange_ = true;
    return status;
  }
  RemoteDonorExchangeStats halo_stats() const noexcept {
    auto result = immersed_exchange_ ? immersed_halo_.stats()
                                    : RemoteDonorExchangeStats{};
    if(dynamic_plan_) {
      const auto s=dynamic_plan_->stats();
      result.received_cells+=s.received_cells;
      result.supplied_cells+=s.supplied_cells;
      result.bytes_per_exchange+=5*s.bytes_per_exchange;
      result.peer_messages+=5*s.peer_messages;
    }
    return result;
  }
  RemoteDonorExchangeCounters immersed_counters() const noexcept {
    return immersed_exchange_ ? immersed_halo_.runtime_counters()
                              : RemoteDonorExchangeCounters{};
  }
  Status
  preflight_immersed_exchange(Span<const FieldView> fields) const noexcept {
    return immersed_exchange_ ? immersed_halo_.preflight_exchange(10, fields)
                              : Status{};
  }
  Status exchange_immersed(Span<FieldView> fields) noexcept {
    return immersed_exchange_ ? immersed_halo_.exchange(10, fields) : Status{};
  }
  void bind_immersed(const IbmEquationInterfacePlan &plan) noexcept {
    immersed_ = &plan;
  }
  Status validate_restart_cell(const std::vector<RestartImageField> &fields,
                               std::size_t start, std::size_t cell,
                               Span<const double> mean_species, double mean_h,
                               double pressure,
                               const ProductReactionSources &gas,
                               const ThermodynamicsPlan &thermo,
                               portable::Revision revision) noexcept {
    if (fields.size() != start + spec_.fields + 2 ||
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
    const auto &cache = fields[start + spec_.fields];
    if (cache.role != RestartFieldRole::stochastic_transport ||
        cache.components != 4 || cache.values.size() != 4 * count_ ||
        !std::isfinite(cache.values[4 * cell]) ||
        !(cache.values[4 * cell] > 0) ||
        !std::isfinite(cache.values[4 * cell + 1]) ||
        cache.values[4 * cell + 1] < 0 ||
        !std::isfinite(cache.values[4 * cell + 2]) ||
        !(cache.values[4 * cell + 2] > 0) ||
        !std::isfinite(cache.values[4 * cell + 3]) ||
        cache.values[4 * cell + 3] < 0)
      return numerical();
    const auto& auxiliary=fields.back();
    if (auxiliary.role!=RestartFieldRole::stochastic_auxiliary ||
        auxiliary.components!=stride_ || auxiliary.values.size()!=count_*stride_)
      return invalid();
    return query_auxiliary(gas,thermo,auxiliary.values.data()+cell*stride_,pressure,revision);
  }
  std::uint64_t owned_bytes() const noexcept {
    if (!enabled())
      return 0;
    std::uint64_t bytes =
        sizeof(*this) + immersed_bytes_ + workspace_->owned_bytes() +
        mixture_storage_.counters().aligned_payload_bytes +
        mixture_species_.capacity() * sizeof(ConstFieldView) +
        tcr_history.owned_bytes() +
        (dynamic_plan_ ? dynamic_plan_->owned_bytes() : 0) +
        mixing_groups_.capacity()*sizeof(tcr::detail::MixingGroup) +
        spec_.tcr.fuel.capacity() + 1 +
        transport_rows_.capacity()*sizeof(ColdPressureRow) +
        (transport_pc_ ? sizeof(ColdPressureDilu)+transport_pc_->owned_payload_bytes() : 0) +
        reactants_.capacity() * sizeof(std::size_t) +
        passive_schmidt_.capacity() * sizeof(std::array<double, 2>) +
        (spec_.initial_species_offsets.capacity() +
         spec_.tcr.progress_weights.capacity()) *
            sizeof(double) +
        spec_.tcr.reactants.capacity() * sizeof(std::string);
    for (const auto &name : spec_.tcr.reactants)
      bytes += name.capacity() + 1;
    for (const auto *v :
         {&rates_, &gradients_, &scratch_, &mass_divergence_, &tuple_, &means_,
          &noise_bounds_, &noise_, &reactor_density_, &transport_carrier_density_,
          &independent_, &diffusion_, &enthalpies_, &query_rates_, &thermal_coordinate_,
          &field_rates_, &field_pressures_, &field_densities_})
      bytes += v->capacity() * sizeof(double);
    return bytes + wiener_.capacity()*sizeof(std::array<double,3>);
  }
  bool enabled() const noexcept { return bool(workspace_); }
  bool dynamic_tcr() const noexcept { return spec_.tcr.model==TcrModel::cdphyso_dynamic_v1; }
  bool implicit_transport() const noexcept { return enabled() && (spec_.tcr.mode!=TcrMode::experimental || dynamic_tcr()); }
  std::vector<ColdPressureRow>& transport_rows() noexcept { return transport_rows_; }
  ColdPressureDilu& transport_preconditioner() noexcept { return *transport_pc_; }
  double mixing_cd() const noexcept { return 2/c_z_; }
  double mixing_cd(std::size_t cell,std::size_t component) const noexcept {
    if(!dynamic_tcr() || spec_.tcr.mode==TcrMode::shadow)return mixing_cd();
    const auto *p=tcr_history.dynamic()->accepted(cell)+5*ns_;
    return tcr::detail::dynamic_group_cd({p[0],p[1],p[2]},mixing_groups_[component]);
  }
  double mixing_control(std::size_t cell,std::size_t component) const noexcept {
    if(!dynamic_tcr() || spec_.tcr.mode==TcrMode::shadow)return 1.;
    const auto species=component==ns_ ? product_species_ : component;
    return tcr_history.dynamic()->accepted(cell)[5*species+3];
  }
  Status configure_dynamic_exchange(MPI_Comm comm,const CartesianGeometryPlan &geometry,
      MeshPatch patch,std::uint64_t maximum_bytes) {
    if(!dynamic_tcr())return {};
    dynamic_plan_.emplace();
    return dynamic_plan_->configure(comm,geometry,patch,
        {!physical_boundary_[0],!physical_boundary_[2],!physical_boundary_[4]},
        ns_,dynamic_groups_,maximum_bytes);
  }
  Status finish_dynamic(Span<const FieldView> fields,ConstFieldView auxiliary,
      ConstFieldView density,std::uint64_t step) noexcept {
    return dynamic_tcr() && dynamic_plan_ ? dynamic_plan_->finish(
        *tcr_history.dynamic(),fields,auxiliary,density,
        immersed_ ? immersed_->cell_activity() : Span<const std::uint8_t>{},step) :
        dynamic_tcr() ? invalid() : Status{};
  }
  double frozen_noise_source(std::size_t cell,std::size_t field,std::size_t component) const noexcept {
    return rates_[slot(cell,field,component)];
  }
  double pressure_source(std::size_t cell) const noexcept { return mass_divergence_[cell]; }
  // PDF transport uses the last accepted pressure increment, divided by its
  // own time interval, plus the refreshed accepted U dot grad(p). Reuse the
  // mass-divergence workspace after prepare_transport has frozen all noise.
  Status prepare_pressure_work(const CartesianKernelPlan& kernels,
      PrimitiveHistory pressure, ConstFieldView velocity, double reference,
      double previous_reference, double previous_dt) noexcept {
    if (!implicit_transport() || !transport_pending_ || !std::isfinite(previous_dt) ||
        previous_dt<0 || !std::isfinite(reference) || !std::isfinite(previous_reference))
      return invalid();
    auto gradient=scratch_view(pressure.accepted,3);
    const auto p=pressure.accepted;
    auto status=cartesian_gradient(kernels,{{&p,1},{&gradient,1},
        {{0,0,0},cells_},0,0,1,0,nullptr});
    if(status && immersed_)status=immersed_->correct_pressure_gradient(p,gradient);
    if(!status)return status;
    std::size_t i=0;
    for(int z=0;z<cells_.z;++z)for(int y=0;y<cells_.y;++y)for(int x=0;x<cells_.x;++x,++i) {
      mass_divergence_[i]=0;
      if(!active(i))continue;
      const Int3 cell{x,y,z};
      double work=previous_dt>0 ? ((reference-previous_reference)+
          (p.unchecked(cell,0)-pressure.previous.unchecked(cell,0)))/previous_dt : 0;
      for(unsigned a=0;a<3;++a)work+=velocity.unchecked(cell,a)*gradient.unchecked(cell,a);
      if(!std::isfinite(work))return numerical();
      mass_divergence_[i]=work;
    }
    return {};
  }
  bool common_transport() const noexcept { return enabled() && common_transport_; }
  Status prepare_mixture_faces(const CartesianKernelPlan &kernels,
                               const ThermodynamicsPlan &thermo,
                               Span<const ConstFieldView> species,
                               ConstFieldView enthalpy, ConstFieldView gamma,
                               ConstFaceFluxView flux, RevisionToken generation,
                               MixtureTransportFaces &mixture) noexcept {
    if (!common_transport() || species.size != ns_ - 1 || !species.data)
      return invalid();
    const std::size_t nx = std::size_t(cells_.x) + 4;
    const std::size_t ny = std::size_t(cells_.y) + 4;
    auto thermal = scratch_view(enthalpy, 1);
    thermal.base = thermal_coordinate_.data() + 2 + 2 * nx + 2 * nx * ny;
    thermal.ghosts = {2, 2, 2};
    thermal.stride_y = nx; thermal.stride_z = nx * ny;
    thermal.component_stride = thermal_coordinate_.size();
    thermal.storage_identity = reinterpret_cast<std::uintptr_t>(thermal_coordinate_.data());
    // Each stochastic field and the mean use their own h/Y coordinates,
    // with the same compiled reference and Cartesian face stencil.
    for (int z = -2; z < cells_.z + 2; ++z)
      for (int y = -2; y < cells_.y + 2; ++y)
        for (int x = -2; x < cells_.x + 2; ++x) {
          if (int(x < 0 || x >= cells_.x) + int(y < 0 || y >= cells_.y) +
              int(z < 0 || z >= cells_.z) > 1) continue;
          const Int3 c{x, y, z};
          for (std::size_t s = 0; s < species.size; ++s)
            independent_[s] = species.data[s].unchecked(c, 0);
          const auto status = thermo.transport_enthalpy_coordinate_from_h(
              enthalpy.unchecked(c, 0), {independent_.data(), independent_.size()},
              thermal.unchecked(c, 0));
          if (!status) return status;
        }
    FaceFluxView workspace;
    auto status = mixture_storage_.workspace_view(0, generation, workspace);
    if (status) status = prepare_cartesian_mixture_transport(kernels,
        species, as_const(thermal), gamma, flux, workspace, generation,
        mixture, immersed_, implicit_transport() ? MixtureFlatStencilPolicy::upwind_constraint
                                                : MixtureFlatStencilPolicy::ignore_roundoff);
    return status;
  }
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
  // Explicit transport carries its own statistical auxiliary mean. The
  // implicit COAST schedule replaces this with its independent field0 solve.
  Status capture_auxiliary(Span<const FieldView> fields, FieldView auxiliary) const noexcept {
    if(fields.size!=spec_.fields || auxiliary.components!=stride_ ||
       auxiliary.interior.x!=cells_.x || auxiliary.interior.y!=cells_.y || auxiliary.interior.z!=cells_.z)
      return invalid();
    for(std::size_t f=0;f<fields.size;++f)
      if(fields.data[f].components!=stride_ || fields.data[f].interior.x!=cells_.x ||
         fields.data[f].interior.y!=cells_.y || fields.data[f].interior.z!=cells_.z) return invalid();
    for(int z=0;z<cells_.z;++z)for(int y=0;y<cells_.y;++y)for(int x=0;x<cells_.x;++x)
      for(std::size_t c=0;c<stride_;++c) {
        long double mean=0;
        for(std::size_t f=0;f<fields.size;++f)mean+=static_cast<long double>(fields.data[f].unchecked({x,y,z},c))/fields.size;
        auxiliary.unchecked({x,y,z},c)=static_cast<double>(mean);
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
    if (cache.components!=4) return invalid();
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
              !(gamma > 0) || !std::isfinite(turbulent) || turbulent < 0 ||
              !std::isfinite(mu.unchecked(cell,0)) || mu.unchecked(cell,0)<=0)
            return numerical();
          cache.unchecked(cell, 0) = gamma;
          cache.unchecked(cell, 1) = turbulent / (sc_t_ * density);
          // Frozen material belongs to the accepted step, including retries.
          // IEM uses dynamic viscosities; diffusion keeps its own coefficient.
          cache.unchecked(cell, 2) = mu.unchecked(cell,0);
          cache.unchecked(cell, 3) = turbulent;
        }
    return {};
  }
  // Apply one common conservative exchange to every stochastic field before
  // querying TCR or forming any transport gradient. These are workspace views;
  // accepted native fields and the gas transaction remain untouched.
  Status deposit_sources(const portable::ExchangeBatchReport &exchange,
                         const ProductReactionSources &gas,
                         const ThermodynamicsPlan &thermo,
                         const TransportPlan &transport,
                         portable::Revision revision, ConstFieldView density,
                         ConstFieldView pi, double reference,
                         Span<const ConstFieldView> fields,
                         ConstFieldView accepted_cache, FieldView post_density,
                         Span<FieldView> post_fields,
                         Span<FieldView> post_species, FieldView post_enthalpy,
                         FieldView post_cache) noexcept {
    if (!enabled() || !exchange.available || exchange.cell_count != count_ ||
        !exchange.cells || fields.size != spec_.fields ||
        post_fields.size != spec_.fields || post_species.size != ns_ - 1 ||
        accepted_cache.components!=4 || post_cache.components!=4)
      return invalid();
    const auto mapping = gas.species_indices();
    std::size_t i = 0;
    for (int z = 0; z < cells_.z; ++z)
      for (int y = 0; y < cells_.y; ++y)
        for (int x = 0; x < cells_.x; ++x, ++i) {
          const Int3 cell{x, y, z};
          const auto &row = exchange.cells[i];
          const double mass = density.unchecked(cell, 0) * row.volume_m3;
          const double next = row.gas.gas_mass_candidate_kg;
          if (!(mass > 0) || !(next > 0) || !(row.volume_m3 > 0) ||
              std::abs(next - mass - row.gas.gas_mass_delta_kg) >
                  1e-12 * std::max(mass, next))
            return numerical();
          const double new_density = row.gas.gas_mass_delta_kg == 0
                                         ? density.unchecked(cell, 0)
                                         : next / row.volume_m3;
          post_density.unchecked(cell, 0) = new_density;
          std::fill(means_.begin(), means_.end(), 0.);
          for (std::size_t f = 0; f < spec_.fields; ++f) {
            auto *tuple = tuple_.data() + f * stride_;
            for (std::size_t c = 0; c < stride_; ++c) {
              const double delta =
                  c == ns_ ? row.gas.gas_thermochemical_enthalpy_delta_j
                           : row.gas_species_mass_delta_kg[c];
              tuple[c] =
                  row.gas.gas_mass_delta_kg == 0 && delta == 0
                      ? fields.data[f].unchecked(cell, c)
                      : (mass * fields.data[f].unchecked(cell, c) + delta) /
                            next;
              means_[c] += tuple[c] / spec_.fields;
            }
            portable::GasSample sample;
            auto status =
                query(gas, thermo, tuple, reference + pi.unchecked(cell, 0),
                      revision, sample);
            if (!status)
              return status;
            for (std::size_t c = 0; c < stride_; ++c)
              post_fields.data[f].unchecked(cell, c) = tuple[c];
          }
          for (std::size_t s = 0; s < mapping.size; ++s) {
            independent_[s] = means_[mapping.data[s]];
            post_species.data[s].unchecked(cell, 0) = independent_[s];
          }
          post_enthalpy.unchecked(cell, 0) = means_[ns_];
          portable::GasSample sample;
          auto status =
              query(gas, thermo, means_.data(),
                    reference + pi.unchecked(cell, 0), revision, sample);
          MolecularTransportState intrinsic;
          if (status)
            status = transport.evaluate(
                sample.temperature_k,
                {independent_.data(), independent_.size()}, intrinsic);
          if (!status)
            return status;
          // Read the accepted eddy viscosity from the shared material cache.
          const double turbulent = accepted_cache.unchecked(cell, 3);
          const double gamma =
              intrinsic.conductivity / sample.cp_j_per_kg_k + turbulent / sc_t_;
          if (!std::isfinite(gamma) || gamma <= 0 ||
              !std::isfinite(turbulent) || turbulent < 0)
            return numerical();
          post_cache.unchecked(cell, 0) = gamma;
          post_cache.unchecked(cell, 1) = turbulent / (sc_t_ * new_density);
          post_cache.unchecked(cell, 2) = intrinsic.viscosity;
          post_cache.unchecked(cell, 3) = turbulent;
        }
    return {};
  }
  Status source_diffusion(const CartesianKernelPlan &kernels,
                          const EnthalpyEquationPlan &enthalpy_plan,
                          ConstFieldView post_cache, ConstFieldView enthalpy,
                          Span<const ConstFieldView> species,
                          FieldView enthalpy_rate,
                          Span<FieldView> species_rates,
                          Span<const ConstFieldView> passives,
                          Span<FieldView> passive_rates,
                          FieldView diffusivity) noexcept {
    if (species.size != ns_ - 1 || species_rates.size != species.size ||
        passives.size != passive_schmidt_.size() ||
        passive_rates.size != passives.size)
      return invalid();
    auto gamma = post_cache;
    gamma.components = 1;
    const auto evaluate = [&](ConstFieldView field, FieldView rate) noexcept {
      KernelInvocation call{{&field, 1}, {&rate, 1}, {{0, 0, 0}, cells_}, 0, 0,
                            1,           0};
      auto result = cartesian_diffusion(kernels, gamma, call);
      if (result && immersed_)
        result =
            field.field == enthalpy.field
                ? immersed_->correct_zero_normal_diffusion(field, gamma, rate)
                : immersed_->correct_impermeable_scalar_diffusion(field, gamma,
                                                                  rate);
      return result;
    };
    auto status = evaluate(enthalpy, enthalpy_rate);
    if (status)
      status = apply_heat_flux_boundary(enthalpy_plan, kernels, enthalpy,
          gamma, {{0, 0, 0}, cells_}, enthalpy_rate,
          immersed_ ? immersed_->cell_activity() : Span<const std::uint8_t>{});
    for (std::size_t s = 0; s < species.size && status; ++s)
      status = evaluate(species.data[s], species_rates.data[s]);
    for (std::size_t s = 0; s < passives.size && status; ++s) {
      for (int z = -2; z < cells_.z + 2; ++z)
        for (int y = -2; y < cells_.y + 2; ++y)
          for (int x = -2; x < cells_.x + 2; ++x) {
            const Int3 c{x, y, z};
            diffusivity.unchecked(c, 0) =
                post_cache.unchecked(c, 2) / passive_schmidt_[s][0] +
                post_cache.unchecked(c, 3) / passive_schmidt_[s][1];
          }
      const auto field = passives.data[s];
      auto rate = passive_rates.data[s];
      KernelInvocation call{{&field, 1}, {&rate, 1}, {{0, 0, 0}, cells_}, 0, 0,
                            1,           0};
      status = cartesian_diffusion(kernels, as_const(diffusivity), call);
      if (status && immersed_)
        status = immersed_->correct_impermeable_scalar_diffusion(
            field, as_const(diffusivity), rate);
    }
    return status;
  }
  Status prepare_transport(MPI_Comm communicator, const CartesianKernelPlan &kernels,
                 const EnthalpyEquationPlan &enthalpy_plan,
                 const ProductReactionSources &gas,
                 const ThermodynamicsPlan &thermo,
                 Span<const ConstFieldView> accepted, Span<FieldView> trial,
                 ConstFieldView cache, ConstFieldView rho, ConstFieldView pi,
                 Span<const ConstFieldView> mean_species, ConstFieldView mean_h,
                 double pressure_reference, ConstFaceFluxView flux, double time,
                 double dt, std::uint64_t step, RevisionToken generation,
                 Span<FieldView> sources) noexcept {
    transport_ready_ = false;
    transport_pending_ = false;
    correction_seed_ready_ = false;
    unsigned invalid_input = accepted.size != spec_.fields || trial.size != spec_.fields ||
        sources.size != ns_ - 1 || mean_species.size != ns_ - 1 || !(dt > 0) ||
        !std::isfinite(dt) || !std::isfinite(1/dt) || cache.components!=4 ||
        !valid_cell_view(rho,cells_,0,1,0) || !valid_flux_view(flux,cells_,flux.revision);
    if(MPI_Allreduce(MPI_IN_PLACE,&invalid_input,1,MPI_UNSIGNED,MPI_MAX,
                     communicator)!=MPI_SUCCESS) return {StatusCode::mpi_failure,10215};
    if(invalid_input) return invalid();
    if(dynamic_tcr()) {
      const auto begun=tcr_history.dynamic()->begin(step);
      if(!begun)return begun;
    }
    // Freeze the transport's extensive carrier before any chemistry. It is
    // distinct from each reactor's PH/EOS density and from field0's density.
    std::size_t carrier_cell=0;
    for(int z=0;z<cells_.z;++z)for(int y=0;y<cells_.y;++y)
      for(int x=0;x<cells_.x;++x,++carrier_cell) {
        const Int3 c{x,y,z};
        const double accepted=rho.unchecked(c,0),volume=cell_volume(kernels,c);
        const double carrier=implicit_transport() && active(carrier_cell)
            ? frozen_density_carrier_mass(flux,c,accepted*volume,1/dt)/volume : accepted;
        if(!std::isfinite(accepted) || accepted<=0 || !std::isfinite(carrier) || carrier<=0)
          invalid_input=1;
        transport_carrier_density_[carrier_cell]=carrier;
      }
    if(MPI_Allreduce(MPI_IN_PLACE,&invalid_input,1,MPI_UNSIGNED,MPI_MAX,communicator)!=MPI_SUCCESS)
      return {StatusCode::mpi_failure,10215};
    if(invalid_input)return numerical();
    transport_density_authority_=rho;
    // One frozen global envelope per field/component for this attempt.
    // Signed maxima share the minimum reduction with the lower bounds.
    const auto size=tuple_.size();
    std::fill(noise_bounds_.begin(),noise_bounds_.end(),
              std::numeric_limits<double>::infinity());
    noise_bounds_.back()=0;
    std::size_t index{};
    for(int z=0;z<cells_.z;++z) for(int y=0;y<cells_.y;++y)
      for(int x=0;x<cells_.x;++x,++index) {
        if(!active(index)) continue;
        for(std::size_t f=0;f<spec_.fields;++f) for(std::size_t c=0;c<stride_;++c) {
          const double value=accepted.data[f].unchecked({x,y,z},c);
          const auto j=f*stride_+c;
          if(!std::isfinite(value)) noise_bounds_.back()=-1;
          noise_bounds_[j]=std::min(noise_bounds_[j],value);
          noise_bounds_[size+j]=std::min(noise_bounds_[size+j],-value);
        }
      }
    if(MPI_Allreduce(MPI_IN_PLACE,noise_bounds_.data(),int(noise_bounds_.size()),
                     MPI_DOUBLE,MPI_MIN,communicator)!=MPI_SUCCESS)
      return {StatusCode::mpi_failure,10215};
    if(noise_bounds_.back()<0) return numerical();
    for(std::size_t j=0;j<size;++j) {
      if(j%stride_<ns_) noise_bounds_[j]=std::max(0.,noise_bounds_[j]);
      noise_bounds_[size+j]=-noise_bounds_[size+j];
    }
    const auto wiener=esf::detail::balanced_wiener(spec_.fields,dt,
        {spec_.seed,step,1,0,0,1},wiener_.data(),wiener_.size());
    if(wiener!=portable::Status::success) return invalid();
    auto scratch = scratch_view(rho, 1);
    KernelInvocation call{{}, {&scratch, 1}, {{0, 0, 0}, cells_}, 0, 0,
                          1,  flux.revision};
    Status status;
    if (!implicit_transport()) {
      status = cartesian_face_divergence(kernels, flux, call);
      if (!status) return status;
      std::copy(scratch_.begin(), scratch_.begin() + count_, mass_divergence_.begin());
    }
    ConstFieldView gamma = cache;
    gamma.components = 1;
    const auto mapping = gas.species_indices();
    for (std::size_t f = 0; f < spec_.fields; ++f) {
      const auto input = accepted.data[f];
      MixtureTransportFaces mixture;
      if (common_transport_ && !implicit_transport()) {
        for (std::size_t s = 0; s < mapping.size; ++s) {
          auto scalar = input;
          scalar.base += mapping.data[s] * scalar.component_stride;
          scalar.components = 1;
          mixture_species_[s] = scalar;
        }
        auto h = input;
        h.base += ns_ * h.component_stride;
        h.components = 1;
        status = prepare_mixture_faces(kernels, thermo,
            {mixture_species_.data(), mixture_species_.size()}, h, gamma,
            flux, generation, mixture);
        if (!status) return status;
      }
      for (std::uint8_t c = 0; c < stride_; ++c) {
        scratch = scratch_view(rho, 1);
        call.reads = {&input, 1};
        call.read_component_begin = c;
        call.required_face_flux_revision = flux.revision;
        auto scalar = input;
        scalar.base += c * scalar.component_stride;
        scalar.components = 1;
        if (!implicit_transport()) {
        status = common_transport_
            ? cartesian_mixture_transport(kernels, mixture, gamma, flux, call)
            : cartesian_convection(
                  kernels, c == ns_ ? enthalpy_scheme_ : species_scheme_, flux, call);
        if (status && immersed_) {
          IbmScalarTransport::Field field{IbmScalarTransport::Quantity::enthalpy,0};
          if (c == gas.dependent_index())
            field.quantity = IbmScalarTransport::Quantity::dependent_species;
          else if (c != ns_) {
            field.quantity = IbmScalarTransport::Quantity::independent_species;
            const auto found = std::find(mapping.data,mapping.data+mapping.size,c);
            if (found == mapping.data+mapping.size) return invalid();
            field.component = std::size_t(found-mapping.data);
          }
          status = common_transport_
              ? IbmScalarTransport::transport(*immersed_,field,scalar,gamma,flux,
                    {{0,0,0},cells_},scratch,mixture)
              : IbmScalarTransport::convection(*immersed_,field,
                    c == ns_ ? enthalpy_scheme_ : species_scheme_,scalar,flux,
                    {{0,0,0},cells_},scratch);
        }
        if (!status)
          return status;
        for (std::size_t i = 0; i < count_; ++i)
          rates_[slot(i, f, c)] = -scratch_[i];
        call.required_face_flux_revision = 0;
        if (!common_transport_)
          status = cartesian_diffusion(kernels, gamma, call);
        else
          std::fill_n(scratch_.begin(), count_, 0.0);
        if (status && immersed_ && !common_transport_)
          status = IbmScalarTransport::diffusion(*immersed_,scalar,gamma,
              {{0,0,0},cells_},scratch);
        if (status && c == ns_)
          status = apply_heat_flux_boundary(enthalpy_plan, kernels, scalar,
              gamma, {{0, 0, 0}, cells_}, scratch,
              immersed_ ? immersed_->cell_activity() : Span<const std::uint8_t>{});
        if (!status)
          return status;
        for (std::size_t i = 0; i < count_; ++i)
          rates_[slot(i, f, c)] += scratch_[i];
        }
        call.required_face_flux_revision=0;
        scratch = scratch_view(rho, 3);
        status = cartesian_gradient(kernels, call);
        if (status && immersed_)
          // This native operator is the scalar homogeneous Neumann gradient;
          // it replaces each cut-link value using its sealed reconstruction.
          status = immersed_->correct_pressure_gradient(scalar, scratch);
        if (!status)
          return status;
        for (std::size_t i = 0; i < count_; ++i)
          for (std::size_t d = 0; d < 3; ++d)
            gradients_[3 * slot(i, f, c) + d] = scratch_[i + d * count_];
      }
    }
    for (int z = 0; z < cells_.z; ++z)
      for (int y = 0; y < cells_.y; ++y)
        for (int x = 0; x < cells_.x; ++x) {
          Int3 cell{x, y, z};
          const auto i = std::size_t(x) + std::size_t(cells_.x) *
                                              (y + std::size_t(cells_.y) * z);
          if (!active(i)) {
            for (std::size_t f = 0; f < spec_.fields; ++f)
              for (std::size_t c = 0; c < stride_; ++c)
                trial.data[f].unchecked(cell, c) =
                    accepted.data[f].unchecked(cell, c);
            for (std::size_t s = 0; s < sources.size; ++s)
              sources.data[s].unchecked(cell, 0) = 0;
            if (tcr_history.enabled() && !dynamic_tcr() && !tcr_history.stage_inactive(i, step))
              return numerical();
            continue;
          }
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
                  implicit_transport() ? 0. : (rates_[slot(i, f, c)] + value * mass_divergence_[i]) /
                  density;
            }
          const portable::Revision revision{step, generation, 1};
          if(noise_allowed(cell)) {
            for(std::size_t f=0;f<spec_.fields;++f) {
              const esf::detail::StochasticSourceRequest noise_request{
                  density,cache.unchecked(cell,2),cache.unchecked(cell,3),.7,
                  !implicit_transport() ? .7 : .5,dt,
                  wiener_[f],tuple_.data()+f*stride_,
                  gradients_.data()+3*slot(i,f,0),
                  noise_bounds_.data()+f*stride_,noise_bounds_.data()+size+f*stride_,stride_,
                  implicit_transport()};
              const auto result=esf::detail::stochastic_source(noise_request,noise_.data());
              if(result.status!=portable::Status::success) return numerical();
              for(std::size_t c=0;c<stride_;++c)
                rates_[slot(i,f,c)]+=implicit_transport() ? noise_[c] : noise_[c]/density;
            }
          }
          esf::detail::Request request;
          request.accepted = {revision, gas.chemistry_identity().fingerprint,
                              spec_.fields, ns_, tuple_.data()};
          request.expected_revision = revision;
          request.identity = &gas.chemistry_identity();
          request.dt_s = dt;
          esf::detail::IemSource mixing;
          // The existing c_z input maps to Cd=2/c_z; c_z=1 is GTMC Cd=2.
          // TCR control is applied once by Workspace::advance below.
          if (esf::detail::iem_source(cell_volume(kernels,cell),
                  cache.unchecked(cell,2),cache.unchecked(cell,3),
                  2/c_z_,1,0,mixing)!=portable::Status::success ||
              !(mixing.implicit_sink_density>0)) return numerical();
          request.mixing_time_s = density/(2*mixing.implicit_sink_density);
          // The frozen, tuple-limited stochastic RHS is included above.
          request.turbulent_diffusivity_m2_s = 0;
          request.deterministic_rates = rates_.data() + slot(i, 0, 0);
          request.gradients = nullptr;
          request.random = {spec_.seed, step, 1, 0, 0, 1};
          if (tcr_history.enabled() && !dynamic_tcr()) {
            double sum = 0;
            for (std::size_t a = 0; a < mapping.size; ++a) {
              means_[mapping.data[a]] = mean_species.data[a].unchecked(cell, 0);
              sum += means_[mapping.data[a]];
            }
            means_[gas.dependent_index()] = 1 - sum;
            means_[ns_] = mean_h.unchecked(cell, 0);
            const double pressure = pressure_reference + pi.unchecked(cell, 0);
            auto& rates = field_rates_;
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
            if (gas.fold_provider() &&
                mapped.status == tcr::detail::Status::success) {
              const auto *provider = gas.fold_provider();
              const auto global_cell =
                  std::uint64_t(x + begin_.x) +
                  std::uint64_t(global_cells_.x) *
                      (std::uint64_t(y + begin_.y) +
                       std::uint64_t(global_cells_.y) * (z + begin_.z));
              const ProductTcrFoldQuery query{
                  global_cell,        history.revision,
                  revision,           history.initialized,
                  history.input.eta,  history.input.rate_ratio,
                  mapped.input.eta,   mapped.input.rate_ratio,
                  history.branch_sign};
              ProductTcrFoldEvidence evidence;
              if (provider->fingerprint() != gas.fold_identity())
                return invalid();
              status = provider->query(query, evidence);
              if (!status)
                return status;
              if (provider->fingerprint() != gas.fold_identity())
                return invalid();
              if (evidence.supplied) {
                if (evidence.source_identity != gas.fold_identity() ||
                    evidence.global_cell != global_cell ||
                    evidence.history_revision != history.revision ||
                    evidence.input_revision != revision)
                  return invalid();
                tcr_request.fold = {true,
                                    history.revision,
                                    {evidence.eta, evidence.rate_ratio},
                                    evidence.departure_sign};
              }
            }
            const auto candidate = tcr::detail::prepare(history, tcr_request);
            if (!candidate.available)
              return {StatusCode::numerical_failure,
                      10220U + static_cast<std::uint32_t>(candidate.status)};
            if (!tcr_history.stage(i, candidate, step))
              return numerical();
            request.tcr_control = candidate.mixer_control;
          }
          if (implicit_transport()) {
            for (std::size_t f=0;f<spec_.fields;++f) for (std::size_t c=0;c<stride_;++c)
              trial.data[f].unchecked(cell,c)=tuple_[f*stride_+c];
            continue;
          }
          auto moved = workspace_->advance(request);
          if (moved.status != portable::Status::success)
            return numerical();
          for (std::size_t f = 0; f < spec_.fields; ++f) {
            portable::GasSample sample;
            status = query(gas, thermo, moved.candidate.values + f * stride_,
                           pressure_reference + pi.unchecked(cell, 0), revision, sample);
            if (!status)
              return status;
            reactor_density_[i * spec_.fields + f] = sample.density_kg_per_m3;
            for (std::size_t c = 0; c < stride_; ++c)
              trial.data[f].unchecked(cell, c) = moved.candidate.values[f * stride_ + c];
          }
        }
    transport_revision_ = {step, generation, 1};
    transport_time_ = time;
    transport_dt_ = dt;
    transport_ready_ = !implicit_transport();
    transport_pending_ = implicit_transport();
    return {};
  }
  Status finish_implicit_transport(const ProductReactionSources& gas,
      const ThermodynamicsPlan& thermo, Span<const FieldView> trial,
      ConstFieldView pi,double pressure,double time,double dt,
      std::uint64_t step,RevisionToken generation) noexcept {
    if (!implicit_transport() || !transport_pending_ || transport_ready_ || trial.size!=spec_.fields ||
        transport_revision_!=portable::Revision{step,generation,1} ||
        transport_time_!=time || transport_dt_!=dt) return invalid();
    transport_pending_=false;
    std::size_t i{};
    for(int z=0;z<cells_.z;++z) for(int y=0;y<cells_.y;++y) for(int x=0;x<cells_.x;++x,++i) {
      if(!active(i))continue;
      const Int3 cell{x,y,z};
      for(std::size_t f=0;f<spec_.fields;++f) {
        double sum=0;
        for(std::size_t c=0;c<stride_;++c) {
          const double q=trial.data[f].unchecked(cell,c);
          if(!std::isfinite(q) || (c<ns_ && (q<0 || q>1))) return numerical();
          tuple_[c]=q;if(c<ns_)sum+=q;
        }
        if(std::abs(sum-1)>2e-12)return numerical();
        portable::GasSample sample;
        const auto status=query(gas,thermo,tuple_.data(),pressure+pi.unchecked(cell,0),
            transport_revision_,sample);
        if(!status)return status;
        reactor_density_[i*spec_.fields+f]=sample.density_kg_per_m3;
      }
    }
    transport_ready_=true;return {};
  }
  // Enter only after all ranks have accepted the complete spatial/mixing
  // phase. Thermodynamic queries and TCR preparation precede this boundary;
  // the interval provider then reads the assembled field state.
  Status react(const ProductReactionSources &gas, const ThermodynamicsPlan& thermo,
               Span<FieldView> trial, FieldView auxiliary,
               ConstFieldView carrier_density, ConstFieldView pi, double pressure_reference, double time,
               double dt, std::uint64_t step, RevisionToken generation,
               Span<FieldView> sources) noexcept {
    const portable::Revision revision{step, generation, 1};
    const bool ready = transport_ready_;
    transport_ready_ = false;
    correction_seed_ready_ = false;
    if (!ready || transport_revision_ != revision || time != transport_time_ ||
        dt != transport_dt_ || trial.size != spec_.fields || sources.size != ns_ - 1 ||
        auxiliary.components!=stride_ ||
        carrier_density.base!=transport_density_authority_.base ||
        carrier_density.field!=transport_density_authority_.field ||
        carrier_density.revision!=transport_density_authority_.revision ||
        carrier_density.storage_identity!=transport_density_authority_.storage_identity ||
        carrier_density.revision_domain!=transport_density_authority_.revision_domain)
      return invalid();
    const auto mapping = gas.species_indices();
    std::size_t i = 0;
    for (int z = 0; z < cells_.z; ++z)
      for (int y = 0; y < cells_.y; ++y)
        for (int x = 0; x < cells_.x; ++x, ++i) {
          const Int3 cell{x, y, z};
          if (!active(i)) {
            for (std::size_t s = 0; s < sources.size; ++s)
              sources.data[s].unchecked(cell, 0) = 0;
            continue;
          }
          auto& pressures = field_pressures_;
          for (std::size_t f = 0; f < spec_.fields; ++f) {
            pressures[f] = pressure_reference + pi.unchecked(cell, 0);
            for (std::size_t c = 0; c < stride_; ++c)
              tuple_[f * stride_ + c] = trial.data[f].unchecked(cell, c);
          }
          const esf::detail::View transported{revision,
              gas.chemistry_identity().fingerprint, spec_.fields, ns_, tuple_.data()};
          const esf::detail::ReactionRequest reaction{transported, revision,
              &gas.chemistry_identity(), &gas.gas_identity(), pressures.data(),
              reactor_density_.data() + i * spec_.fields, time, dt,
              esf::detail::ReactionIntervals::full};
          const auto reacted = workspace_->react(reaction, *gas.gas_advance());
          const double rho=transport_carrier_density_[i];
          if (reacted.status != portable::Status::success || !(rho>0) || !std::isfinite(rho) ||
              !reacted.mean_integrated_mass_fraction_delta)
            return numerical();
          if(dynamic_tcr()) {
            std::array<double,UINT8_MAX> raw{},normalized{},final_y{},delta{};
            for(std::size_t c=0;c<stride_;++c)raw[c]=auxiliary.unchecked(cell,c);
            const auto coordinates=esf::detail::auxiliary_eos_coordinates(
                {revision,gas.chemistry_identity().fingerprint,1,ns_,raw.data()},
                revision,normalized.data(),ns_);
            if(coordinates.status!=portable::Status::success)return {StatusCode::numerical_failure,10238};
            normalized[ns_]=coordinates.query_enthalpy_j_per_kg;
            portable::GasSample psr;
            auto s=query(gas,thermo,normalized.data(),pressures[0],revision,psr);
            if(!s)return s;
            const portable::GasQuery input{revision,gas.gas_identity().composition_fingerprint,
                portable::GasStateCoordinates::pressure_enthalpy,pressures[0],
                normalized[ns_],psr.temperature_k,normalized.data(),ns_};
            portable::GasAdvanceOutput output{{},final_y.data(),delta.data(),ns_};
            const auto advanced=gas.gas_advance()->advance_gas({input,time,dt},output);
            if(advanced!=portable::Status::success ||
                output.final_mass_fractions!=final_y.data() ||
                output.integrated_species_density_delta_kg_per_m3!=delta.data() ||
                output.capacity!=ns_ || output.completed_duration_s!=dt ||
                output.final_sample.revision!=revision ||
                output.final_sample.composition_fingerprint!=input.composition_fingerprint ||
                !std::isfinite(output.final_sample.enthalpy_j_per_kg) ||
                std::abs(output.final_sample.enthalpy_j_per_kg-input.enthalpy_j_per_kg)>
                    1e-8*std::max(1.,std::abs(input.enthalpy_j_per_kg)) ||
                !std::isfinite(output.final_sample.pressure_pa) ||
                std::abs(output.final_sample.pressure_pa-input.pressure_pa)>
                    32*std::numeric_limits<double>::epsilon()*input.pressure_pa)
              return {StatusCode::numerical_failure,10239};
            double sum{};auto *history=tcr_history.dynamic()->candidate(i);
            for(std::size_t q=0;q<ns_;++q) {
              if(!std::isfinite(final_y[q]) || final_y[q]<0 || final_y[q]>1)return numerical();
              sum+=final_y[q];
              // Source rates are specific mole-number rates. Convert both
              // intervals by the same species molecular weight and time.
              const double mw=gas.gas_identity().molecular_weights_kg_per_kmol[q];
              const double pdf=reacted.mean_integrated_mass_fraction_delta[q]/(dt*mw);
              const double rate=coordinates.positive_weight_sum*(final_y[q]-normalized[q])/(dt*mw);
              const auto control=tcr::detail::cdphyso_species_control(.3,pdf,rate,spec_.tcr.weak_rate_threshold);
              if(!control.available)return {StatusCode::numerical_failure,10240};
              history[5*q]=pdf;history[5*q+1]=rate;history[5*q+2]=control.selected;
              history[5*q+3]=control.effective;history[5*q+4]=unsigned(control.state);
            }
            if(std::abs(sum-1)>2e-12)return numerical();
          }
          for (std::size_t s = 0; s < mapping.size; ++s)
            sources.data[s].unchecked(cell, 0) =
                rho * reacted.mean_integrated_mass_fraction_delta[mapping.data[s]] / dt;
          std::array<double,UINT8_MAX> auxiliary_row{};
          for(std::size_t s=0;s<ns_;++s)
            auxiliary_row[s]=auxiliary.unchecked(cell,s)+reacted.mean_integrated_mass_fraction_delta[s];
          long double delta_h=0;
          for(std::size_t f=0;f<spec_.fields;++f)
            delta_h+=(static_cast<long double>(reacted.candidate.values[f*stride_+ns_])-tuple_[f*stride_+ns_])/spec_.fields;
          auxiliary_row[ns_]=auxiliary.unchecked(cell,ns_)+static_cast<double>(delta_h);
          const auto auxiliary_status=query_auxiliary(gas,thermo,auxiliary_row.data(),
              pressure_reference+pi.unchecked(cell,0),revision);
          if(!auxiliary_status)return auxiliary_status;
          for(std::size_t c=0;c<stride_;++c)auxiliary.unchecked(cell,c)=auxiliary_row[c];
          for (std::size_t f = 0; f < spec_.fields; ++f)
            for (std::size_t c = 0; c < stride_; ++c)
              trial.data[f].unchecked(cell, c) = reacted.candidate.values[f * stride_ + c];
        }
    // Spatial gradients have completed their only use in this proposal.
    // Reuse their owned payload for immutable post-reactor tuples, in native
    // component-major field layout. Outer pressure corrections can restart
    // from these exact tuples without advancing chemistry or recentering Y.
    if (implicit_transport()) {
      std::size_t cell_index{};
      for (int z=0;z<cells_.z;++z) for (int y=0;y<cells_.y;++y)
        for (int x=0;x<cells_.x;++x,++cell_index) {
          const Int3 cell{x,y,z};
          for (std::size_t f=0;f<=spec_.fields;++f)
            for (std::size_t c=0;c<stride_;++c)
              gradients_[(f*stride_+c)*count_+cell_index]=f==spec_.fields
                  ? auxiliary.unchecked(cell,c) : trial.data[f].unchecked(cell,c);
        }
      correction_seed_ready_=true;
    }
    return dynamic_tcr() ? Status{} : tcr_history.seal();
  }
  // Random fields occupy [0,fields); the independent signed auxiliary tuple
  // occupies fields. Views are valid through this proposal's correction
  // phase; every subsequent prepare/react invalidates the phase identity.
  Status correction_seed(std::size_t field, FieldId semantic,
      portable::Revision revision, ConstFieldView& out) const noexcept {
    if (!implicit_transport() || !correction_seed_ready_ || field>spec_.fields ||
        revision!=transport_revision_) return invalid();
    ConstFieldView view;
    view.base=gradients_.data()+field*stride_*count_;
    view.interior=cells_;view.components=static_cast<std::uint8_t>(stride_);
    view.stride_y=cells_.x;view.stride_z=std::size_t(cells_.x)*cells_.y;
    view.component_stride=count_;view.field=semantic;view.revision=revision.input_revision;
    view.storage_identity=reinterpret_cast<StorageIdentity>(gradients_.data());
    view.revision_domain=reinterpret_cast<RevisionDomainIdentity>(this);
    out=view;
    return {};
  }
  // Restart images store an AoS tuple. Reuse the same raw field0 closure
  // used by the live pressure corrector without constructing a strided view.
  Status auxiliary_pressure(Span<const double> row, double pressure,
      const ProductReactionSources& gas, const ThermodynamicsPlan& thermo,
      portable::Revision revision, esf::detail::AuxiliaryPressureState& out) noexcept {
    if (!implicit_transport() || gas.fingerprint()!=gas_fingerprint_ ||
        !row.data || row.size!=stride_ || revision.algorithm_version!=1 ||
        revision.input_revision==0 || !(pressure>0) || !std::isfinite(pressure))
      return invalid();
    return query_auxiliary(gas,thermo,row.data,pressure,revision,&out);
  }
  struct DualState {
    // Physical temperature/material EOS remains explicitly distinct from
    // the density and tangents selected by pressure coupling.
    ThermoState physical{};
    esf::detail::AuxiliaryPressureState pressure{};
    double statistical_density{};
    double physical_enthalpy{};
    double auxiliary_enthalpy{};
    double enthalpy_increment{};
  };
  // The field/auxiliary tuples are the same frozen flux-correction anchor.
  // A common pressure-work increment shifts every h_f and h0 by delta_h;
  // their compositions and enthalpy offsets remain independent. Derivatives
  // in pressure are therefore also derivatives with respect to physical h.
  // This read-only closure also serves accepted/restarted states. The caller
  // supplies their revision; proposal-seed access has its separate lifetime
  // checks. Evaluate into scratch first: an EOS failure preserves outputs.
  Status evaluate_dual_state_cell(Int3 cell, Span<const FieldView> fields,
      ConstFieldView auxiliary, double pressure, double delta_h, Real3 velocity,
      const ProductReactionSources& gas, const ThermodynamicsPlan& thermo,
      portable::Revision revision, Span<double> physical_mean, DualState& out) noexcept {
    if (!implicit_transport() || gas.fingerprint()!=gas_fingerprint_ ||
        revision.algorithm_version!=1 || revision.input_revision==0 ||
        fields.size!=spec_.fields || !fields.data || physical_mean.size!=stride_ ||
        !physical_mean.data || !std::isfinite(pressure) || !(pressure>0) ||
        !std::isfinite(delta_h) || cell.x<0 || cell.x>=cells_.x ||
        cell.y<0 || cell.y>=cells_.y || cell.z<0 || cell.z>=cells_.z ||
        !valid_cell_view(auxiliary,cells_,0,stride_,0)) return invalid();
    for (std::size_t f=0;f<fields.size;++f)
      if (!valid_cell_view(as_const(fields.data[f]),cells_,0,stride_,0)) return invalid();
    if (field_view_overlaps_storage(auxiliary,physical_mean.data,stride_)) return invalid();
    for (std::size_t f=0;f<fields.size;++f)
      if (field_view_overlaps_storage(fields.data[f],physical_mean.data,stride_)) return invalid();
    std::array<double,UINT8_MAX> raw_auxiliary{};
    auto& densities = field_densities_;
    for (std::size_t c=0;c<stride_;++c)
      raw_auxiliary[c]=auxiliary.unchecked(cell,c)+(c==ns_ ? delta_h : 0.);
    esf::detail::AuxiliaryPressureState pressure_state;
    auto status=query_auxiliary(gas,thermo,raw_auxiliary.data(),pressure,revision,&pressure_state);
    if (!status) return status;
    for (std::size_t f=0;f<fields.size;++f) {
      for (std::size_t c=0;c<stride_;++c)
        tuple_[f*stride_+c]=fields.data[f].unchecked(cell,c)+(c==ns_ ? delta_h : 0.);
      portable::GasSample sample;
      status=query(gas,thermo,tuple_.data()+f*stride_,pressure,revision,sample);
      if (!status) return status;
      densities[f]=sample.density_kg_per_m3;
    }
    const auto moments=esf::detail::dual_state_moments(
        {{revision,gas.chemistry_identity().fingerprint,spec_.fields,ns_,tuple_.data()},
         {revision,gas.chemistry_identity().fingerprint,1,ns_,raw_auxiliary.data()},
         revision,densities.data(),pressure_state.density_kg_per_m3},
        {means_.data(),noise_.data(),stride_});
    if (moments.status!=portable::Status::success) return numerical();
    const auto mapping=gas.species_indices();
    for (std::size_t c=0;c<mapping.size;++c) independent_[c]=means_[mapping.data[c]];
    DualState candidate;
    status=thermo.evaluate(pressure,means_[ns_],{independent_.data(),independent_.size()},
                          velocity,candidate.physical);
    if (!status) return status;
    candidate.pressure=pressure_state;
    candidate.statistical_density=moments.statistical_density_kg_per_m3;
    candidate.physical_enthalpy=means_[ns_];
    candidate.auxiliary_enthalpy=raw_auxiliary[ns_];
    candidate.enthalpy_increment=delta_h;
    for (std::size_t c=0;c<stride_;++c)physical_mean.data[c]=means_[c];
    out=candidate;
    return {};
  }
  // Finalize one uncommitted cell after its pressure/energy correction.
  // Y is published solely from the random fields. A common h increment
  // preserves their spread and the independent field0 enthalpy offset.
  // Callers retain the immutable proposal seed for any outer-loop replay.
  Status publish_dual_state_cell(Int3 cell, Span<FieldView> fields,
      FieldView auxiliary, FieldView physical_mean, double pressure,
      double delta_h, Real3 velocity, const ProductReactionSources& gas,
      const ThermodynamicsPlan& thermo, portable::Revision revision,
      DualState& out) noexcept {
    if (!correction_seed_ready_ || revision!=transport_revision_ ||
        !valid_cell_view(as_const(physical_mean),cells_,0,stride_,0) ||
        !valid_cell_view(as_const(auxiliary),cells_,0,stride_,0) ||
        fields.size!=spec_.fields || !fields.data ||
        field_views_overlap(physical_mean,auxiliary)) return invalid();
    for (std::size_t f=0;f<fields.size;++f) {
      if (!valid_cell_view(as_const(fields.data[f]),cells_,0,stride_,0) ||
          field_views_overlap(fields.data[f],physical_mean) ||
          field_views_overlap(fields.data[f],auxiliary)) return invalid();
      for (std::size_t g=0;g<f;++g)
        if (field_views_overlap(fields.data[f],fields.data[g])) return invalid();
    }
    std::array<double,UINT8_MAX> mean{};
    DualState candidate;
    auto status=evaluate_dual_state_cell(cell,{fields.data,fields.size},as_const(auxiliary),
        pressure,delta_h,velocity,gas,thermo,revision,{mean.data(),stride_},candidate);
    if (!status) return status;
    for (std::size_t f=0;f<fields.size;++f)fields.data[f].unchecked(cell,ns_)+=delta_h;
    auxiliary.unchecked(cell,ns_)+=delta_h;
    for (std::size_t c=0;c<stride_;++c)physical_mean.unchecked(cell,c)=mean[c];
    out=candidate;
    return {};
  }
  void discard() noexcept {
    transport_ready_=transport_pending_=correction_seed_ready_=false;
    tcr_history.discard();
  }
  void commit() noexcept {
    transport_ready_=transport_pending_=correction_seed_ready_=false;
    tcr_history.commit();
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
          const auto i = std::size_t(x) + std::size_t(cells_.x) *
                                              (y + std::size_t(cells_.y) * z);
          if (!active(i))
            continue;
          double sum = 0;
          for (std::size_t s = 0; s < mapping.size; ++s) {
            const double target = species.data[s].unchecked(cell, 0);
            sum += target;
            means_[mapping.data[s]] = target;
          }
          means_[gas.dependent_index()] = 1 - sum;
          means_[ns_] = h.unchecked(cell, 0);
          for(std::size_t f=0;f<spec_.fields;++f)
            for(std::size_t c=0;c<stride_;++c)
              tuple_[f*stride_+c]=fields.data[f].unchecked(cell,c);
          const auto centered=workspace_->recenter(
              {revision,gas.chemistry_identity().fingerprint,spec_.fields,ns_,tuple_.data()},
              means_.data());
          if(centered.status!=portable::Status::success) return numerical();
          for (std::size_t f = 0; f < spec_.fields; ++f) {
            const auto *row = centered.candidate.values + f * stride_;
            portable::GasSample sample;
            auto status =
                query(gas, thermo, row, pressure + pi.unchecked(cell, 0),
                      revision, sample);
            if (!status)
              return status;
            for (std::size_t c = 0; c < stride_; ++c)
              fields.data[f].unchecked(cell, c) = row[c];
          }
        }
    return {};
  }

private:
  bool noise_allowed(Int3 cell) const noexcept {
    const int position[]{cell.x+begin_.x,cell.y+begin_.y,cell.z+begin_.z};
    const int dimensions[]{global_cells_.x,global_cells_.y,global_cells_.z};
    for(unsigned a=0;a<3;++a)
      if((position[a]==0 && physical_boundary_[2*a]) ||
         (position[a]==dimensions[a]-1 && physical_boundary_[2*a+1])) return false;
    return true;
  }
  bool active(std::size_t cell) const noexcept {
    return !immersed_ || immersed_->cell_activity().data[cell] !=
                             static_cast<std::uint8_t>(RegionFlag::solid);
  }
  RemoteDonorExchangePlan immersed_halo_;
  const IbmEquationInterfacePlan *immersed_{};
  std::uint64_t immersed_bytes_{};
  bool immersed_exchange_{};
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
  Status query_auxiliary(const ProductReactionSources& gas,
      const ThermodynamicsPlan& thermo,const double* row,double pressure,
      portable::Revision revision, esf::detail::AuxiliaryPressureState* output=nullptr) noexcept {
    std::array<double,UINT8_MAX> normalized{};
    const auto coordinates=esf::detail::auxiliary_eos_coordinates(
        {revision,gas.chemistry_identity().fingerprint,1,ns_,row},revision,normalized.data(),ns_);
    if(coordinates.status!=portable::Status::success)return numerical();
    normalized[ns_]=coordinates.query_enthalpy_j_per_kg;
    portable::GasSample sample;
    const auto status=query(gas,thermo,normalized.data(),pressure,revision,sample);
    if(!status)return status;
    esf::detail::AuxiliaryPressureState pressure_state;
    if (esf::detail::auxiliary_pressure_state(coordinates,sample,revision,
        gas.gas_identity().composition_fingerprint,pressure_state)!=portable::Status::success)
      return numerical();
    if (output) *output=pressure_state;
    return {};
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
  std::optional<DynamicTcrPlan> dynamic_plan_;
  std::array<std::size_t,3> dynamic_groups_{};
  std::vector<tcr::detail::MixingGroup> mixing_groups_;
  std::size_t fuel_{}, product_species_{};
  EsfSpec spec_;
  Int3 cells_{}, begin_{}, global_cells_{};
  std::size_t ns_{}, stride_{}, count_{};
  PlanFingerprint gas_fingerprint_{};
  double c_z_{}, sc_t_{};
  std::vector<std::array<double, 2>> passive_schmidt_;
  ConvectionScheme species_scheme_{}, enthalpy_scheme_{};
  bool common_transport_{};
  bool transport_ready_{},transport_pending_{},correction_seed_ready_{};
  portable::Revision transport_revision_{};
  double transport_time_{}, transport_dt_{};
  std::array<bool,6> physical_boundary_{};
  FaceFluxStorage mixture_storage_;
  std::vector<ConstFieldView> mixture_species_;
  std::vector<double> thermal_coordinate_;
  std::vector<std::array<double, 3>> wiener_;
  std::vector<double> field_rates_, field_pressures_, field_densities_;
  std::unique_ptr<esf::detail::Workspace> workspace_;
  std::vector<ColdPressureRow> transport_rows_;
  std::unique_ptr<ColdPressureDilu> transport_pc_;
  std::vector<double> rates_, gradients_, scratch_, mass_divergence_, tuple_,
      means_, noise_bounds_, noise_, reactor_density_, transport_carrier_density_;
  ConstFieldView transport_density_authority_{};
  std::vector<double> independent_, diffusion_, enthalpies_, query_rates_;
};
} // namespace hundun::v04::detail

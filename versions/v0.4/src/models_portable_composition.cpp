// SPDX-License-Identifier: Apache-2.0
#include "models_portable_composition_detail.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
namespace hundun::v04::portable {
namespace {
void add(ExchangeDelta &a, const ExchangeDelta &b) noexcept {
  a.mass_kg += b.mass_kg;
  a.thermochemical_enthalpy_j += b.thermochemical_enthalpy_j;
  a.kinetic_energy_j += b.kinetic_energy_j;
  for (int d = 0; d < 3; ++d)
    a.momentum_kg_m_per_s[d] += b.momentum_kg_m_per_s[d];
}
bool finite_exchange(const ExchangeDelta &d) noexcept {
  if (!std::isfinite(d.mass_kg) ||
      !std::isfinite(d.thermochemical_enthalpy_j) ||
      !std::isfinite(d.kinetic_energy_j))
    return false;
  for (double p : d.momentum_kg_m_per_s)
    if (!std::isfinite(p))
      return false;
  return true;
}
bool fractions(const double *y, std::size_t ns) noexcept {
  double sum = 0;
  for (std::size_t i = 0; i < ns; ++i) {
    if (!std::isfinite(y[i]) || y[i] < 0 || y[i] > 1)
      return false;
    sum += y[i];
  }
  return std::abs(sum - 1) <= 2e-12;
}
bool same_bytes(MPI_Comm comm, const void *data, std::size_t bytes) noexcept {
  std::uint64_t size = bytes, lo = 0, hi = 0;
  if (MPI_Allreduce(&size, &lo, 1, MPI_UINT64_T, MPI_MIN, comm) !=
          MPI_SUCCESS ||
      MPI_Allreduce(&size, &hi, 1, MPI_UINT64_T, MPI_MAX, comm) !=
          MPI_SUCCESS ||
      lo != hi)
    return false;
  const auto *source = static_cast<const unsigned char *>(data);
  bool equal = true;
  for (std::size_t offset = 0; offset < bytes; offset += 256) {
    unsigned char reference[256];
    const int n = static_cast<int>(std::min<std::size_t>(256, bytes - offset));
    std::memcpy(reference, source + offset, n);
    if (MPI_Bcast(reference, n, MPI_BYTE, 0, comm) != MPI_SUCCESS)
      return false;
    equal &= std::memcmp(reference, source + offset, n) == 0;
  }
  int local = equal, all = 0;
  return MPI_Allreduce(&local, &all, 1, MPI_INT, MPI_MIN, comm) ==
             MPI_SUCCESS &&
         all;
}
bool same_identity(MPI_Comm comm, const GasIdentity &id) noexcept {
  auto string = [&](const std::string &s) {
    return same_bytes(comm, s.data(), s.size());
  };
  auto strings = [&](const std::vector<std::string> &ss) {
    const std::uint64_t n = ss.size();
    if (!same_bytes(comm, &n, sizeof n))
      return false;
    for (const auto &s : ss)
      if (!string(s))
        return false;
    return true;
  };
  return string(id.mechanism_sha256) && string(id.phase) &&
         string(id.enthalpy_reference) && strings(id.species_names) &&
         strings(id.element_names) &&
         same_bytes(comm, id.element_counts.data(),
                    id.element_counts.size() * sizeof(std::uint32_t)) &&
         same_bytes(comm, id.molecular_weights_kg_per_kmol.data(),
                    id.molecular_weights_kg_per_kmol.size() * sizeof(double)) &&
         same_bytes(comm, &id.composition_fingerprint,
                    sizeof id.composition_fingerprint) &&
         same_bytes(comm, &id.closure_fingerprint,
                    sizeof id.closure_fingerprint);
}
} // namespace
bool valid_vapor_mapping(const spray::SprayParcelState &p,
                         const spray::detail::LiquidAsset &a,
                         const GasIdentity &g, std::size_t index) noexcept {
  return a.content_fingerprint != 0 &&
         a.pack.material_fingerprint == a.content_fingerprint &&
         p.liquid_material_fingerprint == a.pack.material_fingerprint &&
         index == a.vapor_species_index && index < g.species_names.size() &&
         index < g.molecular_weights_kg_per_kmol.size() &&
         a.vapor_species_name == g.species_names[index] &&
         a.vapor_molecular_weight_kg_per_kmol ==
             g.molecular_weights_kg_per_kmol[index] &&
         same_gas_identity(a.gas_identity, g);
}
CompositionWorkspace::CompositionWorkspace(std::size_t nc, std::size_t np,
                                           std::size_t segments, std::size_t ns)
    : species_(ns), parcel_capacity_(np), segment_capacity_(segments),
      exchange_(nc, segments, ns), routing_(nc, segments), attempt_ids_(np),
      esf_(ns), cell_inputs_(nc), segments_(segments), cells_(nc), parcels_(np),
      jobs_(np), vapor_indices_(np) {
  if (!ns || ns > std::numeric_limits<std::size_t>::max() / 8 - 1 ||
      nc > std::numeric_limits<std::size_t>::max() / (4 * (ns + 1)))
    throw std::invalid_argument("invalid portable composition capacity");
  fields_.resize(nc * 4 * (ns + 1));
  scratch_.resize(4 * (ns + 1));
}
bool CompositionWorkspace::current(const CompositionReport &r) const noexcept {
  return r.available && r.generation == generation_ && r.cells == cells_.data();
}
CompositionReport
CompositionWorkspace::prepare_parcels(const CompositionInput &q) noexcept {
  auto fail = [&](Status s, std::uint32_t module, std::size_t index) {
    CompositionReport r;
    r.status = s;
    r.failure_module = module;
    r.failure_index = index;
    r.revision = q.revision;
    r.generation = generation_;
    return r;
  };
  const auto ns = species_, n = q.field_count, stride = ns + 1;
  if (q.cell_count > cells_.size() || q.parcel_count > parcel_capacity_)
    return fail(Status::capacity_exceeded, 1, 0);
  if (q.revision.algorithm_version != 1 || (n != 2 && n != 4) || !q.identity ||
      !q.gas || (q.reaction_enabled && !q.chemistry) ||
      (q.cell_count && !q.cells) ||
      (q.injector_count && !q.candidate_injectors) ||
      q.injector_count > parcel_capacity_ ||
      (q.parcel_count &&
       (!q.parcels || !q.vapor_species_indices || !q.parcel_materials ||
        !q.migration || !q.location)) ||
      !std::isfinite(q.start_time_s) || q.start_time_s < 0 ||
      !std::isfinite(q.duration_s) || q.duration_s <= 0 ||
      !std::isfinite(q.start_time_s + q.duration_s))
    return fail(Status::invalid_input, 1, 0);
  const auto &gas_id = q.gas->gas_identity();
  if (q.reaction_enabled &&
      !same_gas_identity(gas_id, q.chemistry->gas_identity()))
    return fail(Status::identity_mismatch, 1, 0);
  for (std::size_t i = 0; i < q.injector_count; ++i) {
    spray::detail::DeterministicInjector validator;
    const auto &v = q.candidate_injectors[i];
    if (!validator.configure(v.spec, v.accepted))
      return fail(Status::invalid_input, 1, i);
    for (std::size_t j = 0; j < i; ++j)
      if (v.spec.injector_id == q.candidate_injectors[j].spec.injector_id)
        return fail(Status::invalid_input, 1, i);
  }
  if (gas_id.closure_fingerprint != q.identity->fingerprint ||
      gas_id.species_names.size() != ns || q.identity->species.size() != ns ||
      gas_id.molecular_weights_kg_per_kmol.size() != ns)
    return fail(Status::identity_mismatch, 1, 0);
  for (std::size_t c = 0; c < q.cell_count; ++c) {
    const auto &in = q.cells[c];
    if (!in.fields || !std::isfinite(in.pressure_pa) || in.pressure_pa <= 0 ||
        in.tcr_history.revision != q.revision ||
        in.tcr_trial.expected_revision != q.revision)
      return fail(Status::invalid_input, 1, c);
    cell_inputs_[c] = in.inventory;
    for (std::size_t f = 0; f < n; ++f)
      if (!fractions(in.fields + f * stride, ns) ||
          !std::isfinite(in.fields[f * stride + ns]))
        return fail(Status::invalid_input, 1, c);
  }
  CompositionReport result;
  result.revision = q.revision;
  result.generation = generation_;
  std::size_t job_count = q.parcel_count, segment_count = 0, final_count = 0;
  id_count_ = q.parcel_count;
  for (std::size_t i = 0; i < job_count; ++i) {
    if (q.parcels[i].revision != q.revision ||
        q.parcels[i].duration_s != q.duration_s ||
        q.parcels[i].elapsed_offset_s != 0 || q.vapor_species_indices[i] >= ns)
      return fail(Status::invalid_input, 1, i);
    if (!q.parcel_materials[i] ||
        !valid_vapor_mapping(q.parcels[i].accepted_parcel,
                             *q.parcel_materials[i], gas_id,
                             q.vapor_species_indices[i]))
      return fail(Status::identity_mismatch, 1, i);
    for (std::size_t j = 0; j < i; ++j)
      if (q.parcels[j].accepted_parcel.id == q.parcels[i].accepted_parcel.id)
        return fail(Status::invalid_input, 1, i);
    jobs_[i] = q.parcels[i];
    vapor_indices_[i] = q.vapor_species_indices[i];
    attempt_ids_[i] = q.parcels[i].accepted_parcel.id;
  }
  for (std::size_t i = 0; i < job_count; ++i) {
    const auto parcel = spray::detail::integrate_parcel_events(jobs_[i]);
    if (!parcel.available)
      return fail(Status::provider_failure, 2, i);
    if (parcel.revision != q.revision)
      return fail(Status::stale_revision, 2, i);
    if (parcel.segment_count > segment_capacity_ - segment_count)
      return fail(Status::capacity_exceeded, 2, i);
    for (std::size_t j = 0; j < parcel.segment_count; ++j) {
      const auto &s = parcel.segments[j];
      const auto &d = s.exchange;
      ExchangeSegment e;
      e.revision = q.revision;
      e.global_cell = s.global_cell;
      e.parcel_id = jobs_[i].accepted_parcel.id;
      e.segment_ordinal = s.ordinal;
      e.vapor_species_index = vapor_indices_[i];
      e.delta = {d.parcel_liquid_mass_delta_kg,
                 d.parcel_momentum_delta_kg_m_per_s,
                 d.parcel_thermochemical_enthalpy_delta_j,
                 d.parcel_kinetic_energy_delta_j, d.thermal_exchange_to_gas_j};
      segments_[segment_count++] = e;
    }
    add(result.outlet, parcel.outlet_inventory);
    add(result.wall, parcel.wall_exchange);
    if (parcel.breakup_requested) {
      const auto &b = parcel.breakup_budget;
      const double dk = b.bulk_kinetic_energy_residual_j;
      const double residual =
          b.surface_energy_increase_j + dk - b.supplied_deformation_energy_j;
      if (!std::isfinite(residual) ||
          !std::isfinite(b.surface_energy_increase_j) || !std::isfinite(dk) ||
          !std::isfinite(b.supplied_deformation_energy_j) ||
          std::abs(residual) >
              jobs_[i].energy_absolute_tolerance_j +
                  jobs_[i].relative_tolerance *
                      (std::abs(b.surface_energy_increase_j) + std::abs(dk) +
                       std::abs(b.supplied_deformation_energy_j)))
        return fail(Status::conservation_failure, 2, i);
      result.breakup_energy.deformation_consumed_j +=
          b.supplied_deformation_energy_j;
      result.breakup_energy.surface_increase_j += b.surface_energy_increase_j;
      result.breakup_energy.bulk_kinetic_increase_j += dk;
      result.breakup_energy.residual_j += residual;
      ++result.breakup_energy.event_count;
      if (!parcel.children.available ||
          parcel.children.child_count > parcel_capacity_ - job_count)
        return fail(Status::capacity_exceeded, 2, i);
      for (std::size_t j = 0; j < parcel.children.child_count; ++j) {
        if (id_count_ == attempt_ids_.size())
          return fail(Status::capacity_exceeded, 2, i);
        attempt_ids_[id_count_++] = parcel.children.children[j].id;
        auto next = jobs_[i];
        next.accepted_parcel = parcel.children.children[j];
        next.accepted_auxiliary = {};
        if (parcel.children_remaining_duration_s > 0) {
          next.elapsed_offset_s += parcel.advanced_duration_s;
          next.duration_s = parcel.children_remaining_duration_s;
          jobs_[job_count] = next;
          vapor_indices_[job_count++] = vapor_indices_[i];
        } else {
          if (final_count >= parcels_.size())
            return fail(Status::capacity_exceeded, 2, i);
          parcels_[final_count++] = {next.accepted_parcel, 0, 0};
        }
      }
    } else if (!parcel.parent_removed) {
      if (final_count >= parcels_.size())
        return fail(Status::capacity_exceeded, 2, i);
      parcels_[final_count++] = {parcel.parcel,
                                 parcel.auxiliary.tab_deformation,
                                 parcel.auxiliary.tab_deformation_rate_per_s,
                                 parcel.auxiliary.breakup_ordinal};
    }
  }
  segment_count_ = segment_count;
  if (!finite_exchange(result.outlet) || !finite_exchange(result.wall))
    return fail(Status::conservation_failure, 2, 0);
  const auto &b = result.breakup_energy;
  if (!std::isfinite(b.deformation_consumed_j) ||
      !std::isfinite(b.surface_increase_j) ||
      !std::isfinite(b.bulk_kinetic_increase_j) || !std::isfinite(b.residual_j))
    return fail(Status::conservation_failure, 2, 0);
  result.available = true;
  result.status = Status::success;
  result.parcel_count = final_count;
  return result;
}
CompositionReport
CompositionWorkspace::prepare_cells(const CompositionInput &q,
                                    CompositionReport result) noexcept {
  auto fail = [&](Status s, std::uint32_t module, std::size_t index) {
    CompositionReport r;
    r.status = s;
    r.failure_module = module;
    r.failure_index = index;
    r.revision = q.revision;
    r.generation = generation_;
    return r;
  };
  const auto ns = species_, n = q.field_count, stride = ns + 1;
  const auto &gas_id = q.gas->gas_identity();
  const auto exchange = exchange_.evaluate(
      q.revision, cell_inputs_.data(), q.cell_count, segments_.data(),
      segment_count_, q.thermal_absolute_tolerance_j,
      q.thermal_relative_tolerance);
  if (!exchange.available)
    return fail(exchange.status, 3, exchange.failure_index);
  GasQueryOutput gas_out{
      {}, scratch_.data(), scratch_.data() + ns, scratch_.data() + 2 * ns, ns};
  auto query_gas = [&](const GasQuery &query) noexcept {
    gas_out = {{},
               scratch_.data(),
               scratch_.data() + ns,
               scratch_.data() + 2 * ns,
               ns};
    if (q.gas->query_gas(query, gas_out) != Status::success ||
        gas_out.diffusivities_m2_per_s != scratch_.data() ||
        gas_out.species_enthalpies_j_per_kg != scratch_.data() + ns ||
        gas_out.net_mass_rates_kg_per_m3_s != scratch_.data() + 2 * ns ||
        gas_out.capacity != ns)
      return false;
    const auto &s = gas_out.sample;
    if (s.revision != q.revision ||
        s.composition_fingerprint != gas_id.composition_fingerprint ||
        s.pressure_pa != query.pressure_pa ||
        !std::isfinite(s.enthalpy_j_per_kg) ||
        std::abs(s.enthalpy_j_per_kg - query.enthalpy_j_per_kg) >
            1e-8 + 1e-12 * std::abs(query.enthalpy_j_per_kg) ||
        !std::isfinite(s.density_kg_per_m3) || s.density_kg_per_m3 <= 0 ||
        !std::isfinite(s.temperature_k) || s.temperature_k <= 0 ||
        !std::isfinite(s.cp_j_per_kg_k) || s.cp_j_per_kg_k <= 0 ||
        !std::isfinite(s.viscosity_pa_s) || s.viscosity_pa_s <= 0 ||
        !std::isfinite(s.conductivity_w_per_m_k) ||
        s.conductivity_w_per_m_k <= 0)
      return false;
    for (std::size_t i = 0; i < ns; ++i)
      if (!std::isfinite(scratch_[i]) || scratch_[i] < 0 ||
          !std::isfinite(scratch_[ns + i]) ||
          !std::isfinite(scratch_[2 * ns + i]))
        return false;
    return true;
  };
  for (std::size_t c = 0; c < q.cell_count; ++c) {
    const auto &in = q.cells[c];
    const auto &e = exchange.cells[c];
    auto *field = fields_.data() + c * 4 * stride;
    const double m = in.inventory.gas_mass_kg,
                 mnew = e.gas.gas_mass_candidate_kg;
    for (std::size_t f = 0; f < n; ++f) {
      for (std::size_t s = 0; s < ns; ++s)
        field[f * stride + s] =
            (m * in.fields[f * stride + s] + e.gas_species_mass_delta_kg[s]) /
            mnew;
      field[f * stride + ns] = (m * in.fields[f * stride + ns] +
                                e.gas.gas_thermochemical_enthalpy_delta_j) /
                               mnew;
      if (!fractions(field + f * stride, ns) ||
          !std::isfinite(field[f * stride + ns]))
        return fail(Status::conservation_failure, 3, c);
    }
    auto trial = in.tcr_trial;
    if (trial.mode != tcr::detail::Mode::off && in.progress_species_weights) {
      double rates[4]{}, psr_rate = 0;
      // means occupy the fourth scratch region, disjoint from provider outputs.
      auto *mean = scratch_.data() + 3 * ns;
      std::fill(mean, mean + stride, 0.0);
      for (std::size_t f = 0; f < n; ++f)
        for (std::size_t j = 0; j < stride; ++j)
          mean[j] += field[f * stride + j] / n;
      for (std::size_t f = 0; f <= n; ++f) {
        const double *state = f < n ? field + f * stride : mean;
        GasQuery query{q.revision,
                       gas_id.composition_fingerprint,
                       GasStateCoordinates::pressure_enthalpy,
                       in.pressure_pa,
                       state[ns],
                       0,
                       state,
                       ns};
        if (!query_gas(query))
          return fail(Status::provider_failure, 4, c);
        double rate = 0;
        for (std::size_t s = 0; s < ns; ++s)
          rate += in.progress_species_weights[s] *
                  gas_out.net_mass_rates_kg_per_m3_s[s] /
                  gas_out.sample.density_kg_per_m3;
        if (f < n)
          rates[f] = rate;
        else
          psr_rate = rate;
      }
      trial.mapping = tcr::detail::ideal_gas_reactant_mole_fraction_v1(
          {q.revision, q.revision, gas_id.composition_fingerprint,
           gas_id.composition_fingerprint, mean,
           gas_id.molecular_weights_kg_per_kmol.data(), ns, in.reactant_indices,
           in.reactant_count, rates, n, psr_rate, in.weak_progress_rate_per_s});
    }
    auto mixing = in.transport;
    mixing.accepted = {q.revision, q.identity->fingerprint, n, ns, field};
    mixing.expected_revision = q.revision;
    mixing.identity = q.identity;
    mixing.dt_s = q.duration_s;
    mixing.random.accepted_step = q.revision.accepted_step;
    const auto mixed = tcr::detail::mix(in.tcr_history, trial, mixing, esf_);
    if (!mixed.available)
      return fail(Status::provider_failure, 5, c);
    std::copy(mixed.esf.candidate.values,
              mixed.esf.candidate.values + n * stride, field);
    double pressures[4]{}, densities[4]{};
    for (std::size_t f = 0; f < n; ++f) {
      pressures[f] = in.pressure_pa;
      GasQuery query{q.revision,
                     gas_id.composition_fingerprint,
                     GasStateCoordinates::pressure_enthalpy,
                     in.pressure_pa,
                     field[f * stride + ns],
                     0,
                     field + f * stride,
                     ns};
      if (!query_gas(query))
        return fail(Status::provider_failure, 4, c);
      densities[f] = gas_out.sample.density_kg_per_m3;
    }
    esf::detail::Report reacted;
    if (q.reaction_enabled) {
      reacted = esf_.react({{q.revision, q.identity->fingerprint, n, ns, field},
                            q.revision,
                            q.identity,
                            &gas_id,
                            pressures,
                            densities,
                            q.start_time_s,
                            q.duration_s},
                           *q.chemistry);
    } else {
      reacted.status = Status::success;
      reacted.candidate.values = field;
      std::copy(densities, densities + n,
                reacted.final_densities_kg_per_m3.begin());
    }
    if (reacted.status != Status::success)
      return fail(reacted.status, 6, c);
    if (q.reaction_enabled)
      std::copy(reacted.candidate.values, reacted.candidate.values + n * stride,
                field);
    CompositionCellCandidate candidate;
    candidate.inventory = {in.inventory.global_cell, in.inventory.volume_m3,
                           mnew, e.gas.gas_momentum_candidate_kg_m_per_s};
    candidate.fields = field;
    candidate.eos_densities_kg_per_m3 = reacted.final_densities_kg_per_m3;
    candidate.tcr = mixed.tcr;
    candidate.heat_release_report_j_per_m3 =
        reacted.ensemble_heat_release_j_per_m3;
    cells_[c] = candidate;
    result.chemistry_calls += reacted.chemistry_call_count;
  }
  result.available = true;
  result.status = Status::success;
  result.reaction_enabled = q.reaction_enabled;
  result.cells = cells_.data();
  result.cell_count = q.cell_count;
  result.field_count = n;
  result.species_count = ns;
  result.parcels = parcels_.data();
  return result;
}
CompositionReport
CompositionWorkspace::prepare(MPI_Comm comm,
                              const CompositionInput &q) noexcept {
  ++generation_;
  if (comm == MPI_COMM_NULL)
    return {};
  if (q.migration)
    q.migration->discard();
  std::uint64_t clock[]{q.revision.accepted_step, q.revision.input_revision,
                        q.revision.algorithm_version, q.field_count, species_};
  std::uint64_t low[5]{}, high[5]{};
  if (MPI_Allreduce(clock, low, 5, MPI_UINT64_T, MPI_MIN, comm) !=
          MPI_SUCCESS ||
      MPI_Allreduce(clock, high, 5, MPI_UINT64_T, MPI_MAX, comm) != MPI_SUCCESS)
    return {};
  if (!std::equal(low, low + 5, high)) {
    CompositionReport r;
    r.status = Status::stale_revision;
    r.failure_module = 1;
    return r;
  }
  const double interval[]{q.start_time_s, q.duration_s};
  if (!same_bytes(comm, interval, sizeof interval) ||
      !same_bytes(comm, &q.reaction_enabled, sizeof q.reaction_enabled)) {
    CompositionReport r;
    r.status = Status::stale_revision;
    r.failure_module = 1;
    return r;
  }
  auto result = prepare_parcels(q);
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  auto consensus = [&]() noexcept {
    int failed = result.available ? std::numeric_limits<int>::max() : rank,
        global = 0;
    if (MPI_Allreduce(&failed, &global, 1, MPI_INT, MPI_MIN, comm) !=
        MPI_SUCCESS) {
      result = {};
      result.status = Status::provider_failure;
      return false;
    }
    if (global != std::numeric_limits<int>::max()) {
      std::uint64_t payload[]{static_cast<std::uint64_t>(result.status),
                              result.failure_module, result.failure_index};
      if (MPI_Bcast(payload, 3, MPI_UINT64_T, global, comm) != MPI_SUCCESS) {
        result = {};
        return false;
      }
      result = {};
      result.status = static_cast<Status>(payload[0]);
      result.failure_module = static_cast<std::uint32_t>(payload[1]);
      result.failure_index = payload[2];
      result.lowest_failing_rank = global;
      result.revision = q.revision;
      result.generation = generation_;
      return false;
    }
    return true;
  };
  if (!consensus())
    return result;
  if (!same_identity(comm, q.gas->gas_identity())) {
    result = {};
    result.status = Status::identity_mismatch;
    result.failure_module = 1;
    return result;
  }
  const auto ids = routing_.audit_ids(comm, attempt_ids_.data(), id_count_);
  if (ids != Status::success) {
    result = {};
    result.status = ids;
    result.failure_module = 2;
    return result;
  }
  const auto routed = routing_.prepare(comm, cell_inputs_.data(), q.cell_count,
                                       segments_.data(), segment_count_);
  if (!routed.available) {
    result = {};
    result.status = routed.status;
    result.failure_module = 3;
    return result;
  }
  if (routed.count)
    std::copy(routed.segments, routed.segments + routed.count,
              segments_.begin());
  segment_count_ = routed.count;
  result = prepare_cells(q, result);
  if (!consensus())
    return result;
  int use = q.migration ? 1 : 0, min_use = 0, max_use = 0;
  MPI_Allreduce(&use, &min_use, 1, MPI_INT, MPI_MIN, comm);
  MPI_Allreduce(&use, &max_use, 1, MPI_INT, MPI_MAX, comm);
  int missing = use && !q.location, any_missing = 0;
  MPI_Allreduce(&missing, &any_missing, 1, MPI_INT, MPI_MAX, comm);
  if (min_use != max_use || any_missing) {
    result = {};
    result.status = Status::invalid_input;
    result.failure_module = 7;
    return result;
  }
  if (use) {
    spray::detail::ParcelMigrationReport migration;
    const auto status =
        q.migration->prepare_checked({parcels_.data(), result.parcel_count},
                                     *q.location, q.revision, migration);
    if (!status || !migration.available) {
      result = {};
      result.status = Status::provider_failure;
      result.failure_module = 7;
      result.lowest_failing_rank = migration.lowest_failing_rank;
      return result;
    }
    auto values = q.migration->candidates();
    int short_capacity = values.size > parcels_.size(), any_short = 0;
    MPI_Allreduce(&short_capacity, &any_short, 1, MPI_INT, MPI_MAX, comm);
    if (any_short) {
      q.migration->discard();
      result = {};
      result.status = Status::capacity_exceeded;
      return result;
    }
    if (values.size)
      std::copy(values.data, values.data + values.size, parcels_.begin());
    result.parcel_count = values.size;
    q.migration->discard();
  }
  return result;
}
CompositionValueReport restore_composition_values(
    const CompositionSnapshot &s, const GasIdentity &id, Revision expected,
    std::size_t max_cells, std::size_t max_parcels,
    const spray::detail::ParcelLocationProvider *location) noexcept {
  CompositionValueReport out;
  if (s.cells.size() > max_cells || s.parcels.size() > max_parcels) {
    out.status = Status::capacity_exceeded;
    return out;
  }
  if (s.version != 1 || s.accepted_revision != expected ||
      expected.algorithm_version != 1 || s.parcel_rng_algorithm_version != 1 ||
      s.injectors.size() > max_parcels ||
      (s.field_count != 2 && s.field_count != 4) ||
      !std::isfinite(s.accepted_time_s) || s.accepted_time_s < 0 ||
      (!s.parcels.empty() && !location))
    return out;
  const auto &a = s.gas_identity;
  const auto ns = id.species_names.size();
  if (ns == 0 || a.mechanism_sha256 != id.mechanism_sha256 ||
      a.phase != id.phase || a.enthalpy_reference != id.enthalpy_reference ||
      a.species_names != id.species_names ||
      a.element_names != id.element_names ||
      a.element_counts != id.element_counts ||
      a.molecular_weights_kg_per_kmol != id.molecular_weights_kg_per_kmol ||
      a.composition_fingerprint != id.composition_fingerprint ||
      a.closure_fingerprint != id.closure_fingerprint) {
    out.status = Status::identity_mismatch;
    return out;
  }
  for (std::size_t i = 0; i < s.cells.size(); ++i) {
    const auto &c = s.cells[i];
    if ((i &&
         s.cells[i - 1].inventory.global_cell >= c.inventory.global_cell) ||
        !std::isfinite(c.inventory.gas_mass_kg) ||
        c.inventory.gas_mass_kg <= 0 || !std::isfinite(c.inventory.volume_m3) ||
        c.inventory.volume_m3 <= 0 || !std::isfinite(c.pressure_pa) ||
        c.pressure_pa <= 0 || c.fields.size() != s.field_count * (ns + 1) ||
        c.tcr_history.revision != expected ||
        c.random.accepted_step != expected.accepted_step ||
        !tcr::detail::valid_history(c.tcr_history))
      return out;
    for (double p : c.inventory.gas_momentum_kg_m_per_s)
      if (!std::isfinite(p))
        return out;
    for (std::size_t f = 0; f < s.field_count; ++f)
      if (!fractions(c.fields.data() + f * (ns + 1), ns) ||
          !std::isfinite(c.fields[f * (ns + 1) + ns]))
        return out;
  }
  for (std::size_t i = 0; i < s.parcels.size(); ++i) {
    const auto &p = s.parcels[i];
    if (spray::validate_parcel_state(p.parcel) !=
            spray::ParcelStateStatus::success ||
        !std::isfinite(p.tab_deformation) ||
        !std::isfinite(p.tab_deformation_rate_per_s))
      return out;
    for (std::size_t j = 0; j < i; ++j)
      if (p.parcel.id == s.parcels[j].parcel.id)
        return out;
    spray::detail::ParcelLocation located;
    if (!location->locate(p.parcel.position_m, expected, located) ||
        located.revision != expected || located.owner_rank < 0 ||
        located.global_cell != p.parcel.owner_global_cell)
      return out;
  }
  for (std::size_t i = 0; i < s.injectors.size(); ++i) {
    spray::detail::DeterministicInjector validator;
    if (!validator.configure(s.injectors[i].spec, s.injectors[i].accepted))
      return out;
    for (std::size_t j = 0; j < i; ++j)
      if (s.injectors[i].spec.injector_id == s.injectors[j].spec.injector_id)
        return out;
  }
  try {
    out.candidate = s;
  } catch (...) {
    out = {};
    out.status = Status::capacity_exceeded;
    return out;
  }
  out.status = Status::success;
  out.available = true;
  return out;
}
CompositionValueReport
prepare_accepted_values(MPI_Comm comm, const CompositionWorkspace &workspace,
                        const CompositionInput &in, const CompositionReport &r,
                        Revision next) noexcept {
  if (comm == MPI_COMM_NULL)
    return {};
  std::uint64_t clocks[]{next.accepted_step, next.input_revision,
                         next.algorithm_version};
  std::uint64_t low[3]{}, high[3]{};
  if (MPI_Allreduce(clocks, low, 3, MPI_UINT64_T, MPI_MIN, comm) !=
          MPI_SUCCESS ||
      MPI_Allreduce(clocks, high, 3, MPI_UINT64_T, MPI_MAX, comm) !=
          MPI_SUCCESS ||
      !std::equal(low, low + 3, high)) {
    CompositionValueReport out;
    out.status = Status::stale_revision;
    return out;
  }
  auto build = [&]() noexcept -> CompositionValueReport {
    CompositionValueReport out;
    if (!workspace.current(r) || !in.gas || r.revision != in.revision ||
        r.cell_count != in.cell_count || r.field_count != in.field_count ||
        next.algorithm_version != 1 ||
        in.revision.accepted_step ==
            std::numeric_limits<std::uint64_t>::max() ||
        next.accepted_step != in.revision.accepted_step + 1 ||
        next.input_revision <= in.revision.input_revision)
      return out;
    try {
      CompositionSnapshot candidate;
      candidate.accepted_revision = next;
      candidate.gas_identity = in.gas->gas_identity();
      candidate.field_count = r.field_count;
      candidate.accepted_time_s = in.start_time_s + in.duration_s;
      candidate.reaction_enabled = in.reaction_enabled;
      candidate.parcel_rng_seed = in.parcel_rng_seed;
      if (in.injector_count)
        candidate.injectors.assign(in.candidate_injectors,
                                   in.candidate_injectors + in.injector_count);
      candidate.cells.reserve(r.cell_count);
      for (std::size_t c = 0; c < r.cell_count; ++c) {
        const auto accepted =
            tcr::detail::accept(in.cells[c].tcr_history, r.cells[c].tcr, next);
        if (!accepted.available)
          return out;
        CompositionCellValue v;
        v.inventory = r.cells[c].inventory;
        v.pressure_pa = in.cells[c].pressure_pa;
        v.tcr_history = accepted.candidate;
        v.random = in.cells[c].transport.random;
        v.random.accepted_step = next.accepted_step;
        v.fields.assign(r.cells[c].fields,
                        r.cells[c].fields +
                            r.field_count * (r.species_count + 1));
        candidate.cells.push_back(std::move(v));
      }
      if (r.parcel_count)
        candidate.parcels.assign(r.parcels, r.parcels + r.parcel_count);
      out.candidate = std::move(candidate);
    } catch (...) {
      out = {};
      out.status = Status::capacity_exceeded;
      return out;
    }
    out.status = Status::success;
    out.available = true;
    return out;
  };
  auto out = build();
  int ok = out.available ? 1 : 0, all_ok = 0;
  if (MPI_Allreduce(&ok, &all_ok, 1, MPI_INT, MPI_MIN, comm) != MPI_SUCCESS ||
      !all_ok) {
    out = {};
    out.status = Status::provider_failure;
  }
  return out;
}
} // namespace hundun::v04::portable

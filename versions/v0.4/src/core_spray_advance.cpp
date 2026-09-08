// SPDX-License-Identifier: Apache-2.0
#include "core_spray_advance_detail.hpp"
#include <climits>
#include <cmath>
#include <cstring>
#include <limits>

namespace hundun::v04::detail {
namespace {
Status invalid(std::uint32_t detail = 10260) noexcept {
  return {StatusCode::invalid_plan, detail};
}
Status capacity() noexcept { return {StatusCode::allocation_failure, 10261}; }
Status unavailable() noexcept { return {StatusCode::rejected_step, 10262}; }
std::uint64_t cell_id(Int3 p, Int3 n) noexcept {
  return std::uint64_t(p.x) +
         std::uint64_t(n.x) *
             (std::uint64_t(p.y) + std::uint64_t(n.y) * std::uint64_t(p.z));
}
} // namespace
Status ProductParcelAdvance::agree(Status local) noexcept {
  int failed = local ? INT_MAX : rank_, first = 0;
  if (MPI_Allreduce(&failed, &first, 1, MPI_INT, MPI_MIN, comm_) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, 10263};
  if (first == INT_MAX)
    return {};
  std::uint32_t value[]{static_cast<std::uint32_t>(local.code), local.detail};
  if (MPI_Bcast(value, 2, MPI_UINT32_T, first, comm_) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, 10263};
  return {static_cast<StatusCode>(value[0]), value[1]};
}
Status ProductParcelAdvance::configure(MPI_Comm comm,
                                       const CartesianGeometryPlan &geometry,
                                       MeshPatch patch, const SpraySpec &spec,
                                       std::uint64_t maximum_bytes) noexcept {
  if (comm == MPI_COMM_NULL)
    return invalid();
  comm_ = comm;
  if (MPI_Comm_rank(comm, &rank_) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, 10263};
  Status local;
  if (geometry_ || spec.liquid_fingerprint != asset_.content_fingerprint ||
      !geometry.fingerprint() || !spec.maximum_local_parcels ||
      !spec.maximum_local_segments || !std::isfinite(spec.maximum_substep_s) ||
      !std::isfinite(spec.minimum_substep_s) || spec.minimum_substep_s <= 0 ||
      spec.maximum_substep_s < spec.minimum_substep_s ||
      !std::isfinite(spec.relative_tolerance) || spec.relative_tolerance <= 0 ||
      spec.relative_tolerance >= 1 || patch.cells.x <= 0 ||
      patch.cells.y <= 0 || patch.cells.z <= 0)
    local = invalid();
  auto status = agree(local);
  if (!status)
    return status;
  std::uint64_t controls[4]{}, minima[4]{}, maxima[4]{};
  std::memcpy(controls, &spec.maximum_substep_s, 8);
  std::memcpy(controls + 1, &spec.minimum_substep_s, 8);
  std::memcpy(controls + 2, &spec.relative_tolerance, 8);
  controls[3] = spec.tab_breakup;
  if (MPI_Allreduce(controls, minima, 4, MPI_UINT64_T, MPI_MIN, comm_) !=
          MPI_SUCCESS ||
      MPI_Allreduce(controls, maxima, 4, MPI_UINT64_T, MPI_MAX, comm_) !=
          MPI_SUCCESS)
    return {StatusCode::mpi_failure, 10263};
  if (!std::equal(minima, minima + 4, maxima))
    return invalid();
  const std::uint64_t nc =
      std::uint64_t(patch.cells.x) * patch.cells.y * patch.cells.z;
  const std::uint64_t np = spec.maximum_local_parcels,
                      ns = spec.maximum_local_segments;
  const std::uint64_t species = asset_.gas_identity.species_names.size();
  if (!species || species > 255 ||
      nc > SIZE_MAX / sizeof(portable::ExchangeCell))
    local = invalid();
  // Includes this lane, batch buffers, and conservative migration buffers.
  // The gas/film adapters and injectors are borrowed and budgeted by their
  // owner.
  int ranks = 0;
  if (MPI_Comm_size(comm, &ranks) != MPI_SUCCESS)
    local = {StatusCode::mpi_failure, 10263};
  const std::uint64_t per_cell = sizeof(portable::ExchangeCell) +
                                 sizeof(portable::CellExchange) +
                                 species * sizeof(double);
  const std::uint64_t fixed =
      sizeof(*this) + sizeof(portable::ExchangeWorkspace) +
      np * (2 * sizeof(Parcel) + sizeof(Job)) +
      ns * (2 * sizeof(portable::ExchangeSegment) + sizeof(std::size_t)) +
      std::uint64_t(ranks) * 5 * sizeof(int) +
      np * (sizeof(int) + 36 * sizeof(std::uint64_t) + sizeof(Parcel)) +
      ns * (2 * sizeof(std::uint64_t) + sizeof(spray::ParcelId)) + 24;
  if (fixed > maximum_bytes || nc > (maximum_bytes - fixed) / per_cell)
    local = capacity();
  status = agree(local);
  if (!status)
    return status;
  const std::uint64_t own = fixed + nc * per_cell;
  status = route_.configure(comm, geometry.global_cells(), patch, ns, ns,
                            species, maximum_bytes - own);
  if (!status)
    return status;
  status = migration_.configure(comm, geometry.global_cells(), patch, np, ns);
  if (!status)
    return status;
  try {
    current_.reserve(np);
    next_.reserve(np);
    jobs_.reserve(np);
    segments_.reserve(ns);
    ids_.reserve(ns);
    cells_.resize(nc);
    batch_ = std::make_unique<portable::ExchangeWorkspace>(nc, ns, species);
  } catch (...) {
    local = capacity();
  }
  status = agree(local);
  if (!status)
    return status;
  geometry_ = &geometry;
  patch_ = patch;
  parcel_capacity_ = np;
  segment_capacity_ = ns;
  maximum_step_ = spec.maximum_substep_s;
  minimum_step_ = spec.minimum_substep_s;
  relative_tolerance_ = spec.relative_tolerance;
  tab_enabled_ = spec.tab_breakup;
  minimum_width_ = std::numeric_limits<double>::max();
  for (unsigned d = 0; d < 3; ++d) {
    const auto widths = geometry.axis(static_cast<CartesianAxis>(d)).widths();
    for (std::size_t i = 0; i < widths.size; ++i)
      minimum_width_ = std::min(minimum_width_, widths.data[i]);
  }
  local_anchor_ = cell_id(patch.begin, geometry.global_cells());
  owned_bytes_ = own + route_.owned_bytes();
  return {};
}
Status
ProductParcelAdvance::append(const portable::ExchangeSegment &s) noexcept {
  if (segments_.size() == segment_capacity_)
    return capacity();
  segments_.push_back(s);
  return {};
}
Status ProductParcelAdvance::audit_id(spray::ParcelId id,
                                      portable::Revision revision) noexcept {
  if (ids_.size() == segment_capacity_)
    return capacity();
  portable::ExchangeSegment row;
  row.revision = revision;
  row.global_cell = local_anchor_;
  row.parcel_id = id;
  row.channel = portable::ExchangeChannel::external;
  row.segment_ordinal = UINT64_MAX;
  ids_.push_back(row);
  return {};
}
Status ProductParcelAdvance::ledger(const portable::ExchangeDelta &delta,
                                    const Parcel &parcel,
                                    portable::Revision revision,
                                    portable::ExchangeChannel channel,
                                    std::uint64_t ordinal) noexcept {
  portable::ExchangeSegment row;
  row.revision = revision;
  row.global_cell = local_anchor_;
  row.parcel_id = parcel.parcel.id;
  row.segment_ordinal = ordinal;
  row.channel = channel;
  row.delta = {delta.mass_kg,
               {delta.momentum_kg_m_per_s[0], delta.momentum_kg_m_per_s[1],
                delta.momentum_kg_m_per_s[2]},
               delta.thermochemical_enthalpy_j,
               delta.kinetic_energy_j,
               0};
  return append(row);
}
Status ProductParcelAdvance::inventory(const spray::SprayParcelState &p,
                                       portable::Revision revision) noexcept {
  const auto h =
      spray::detail::evaluate_liquid_enthalpy(asset_, p.temperature_k);
  if (!h.available)
    return invalid();
  portable::ExchangeDelta delta;
  delta.mass_kg = p.droplet_mass_kg * p.multiplicity;
  delta.thermochemical_enthalpy_j = delta.mass_kg * h.liquid_enthalpy_j_per_kg;
  for (unsigned d = 0; d < 3; ++d) {
    delta.momentum_kg_m_per_s[d] = delta.mass_kg * p.velocity_m_per_s[d];
    delta.kinetic_energy_j +=
        .5 * delta.mass_kg * p.velocity_m_per_s[d] * p.velocity_m_per_s[d];
  }
  return ledger(delta, {p, 0, 0, 0}, revision,
                portable::ExchangeChannel::external, UINT64_MAX - 1);
}
Status ProductParcelAdvance::move() noexcept {
  Status status;
  for (auto &value : next_) {
    spray::detail::ParcelLocation located;
    status = gas_.locate(value.parcel.position_m, revision_, located);
    if (!status)
      break;
    value.parcel.owner_global_cell = located.global_cell;
  }
  status = agree(status);
  if (!status)
    return status;
  spray::detail::ParcelMigrationReport report;
  status = migration_.prepare_checked({next_.data(), next_.size()}, gas_,
                                      revision_, report);
  if (!status)
    return status;
  const auto candidates = migration_.candidates();
  current_.clear();
  if (candidates.size)
    current_.assign(candidates.data, candidates.data + candidates.size);
  return {};
}
Status ProductParcelAdvance::wave(portable::Revision revision, double start,
                                  double duration,
                                  std::uint32_t index) noexcept {
  using namespace spray::detail;
  jobs_.clear();
  next_.clear();
  for (const auto &p : current_)
    jobs_.push_back({p, start, duration});
  std::size_t processed = 0;
  while (!jobs_.empty()) {
    const auto job = jobs_.back();
    jobs_.pop_back();
    if (++processed > segment_capacity_)
      return capacity();
    ParcelEventsInput in;
    in.accepted_parcel = job.value.parcel;
    in.revision = revision;
    in.accepted_auxiliary = {job.value.tab_deformation,
                             job.value.tab_deformation_rate_per_s,
                             job.value.breakup_ordinal};
    in.interval = &interval_;
    in.geometry = &events_;
    in.tab_evolution = tab_enabled_ ? &tab_ : nullptr;
    in.breakup = tab_enabled_ ? &breakup_ : nullptr;
    in.elapsed_offset_s = job.start;
    in.duration_s = job.duration;
    in.initial_substep_s = std::min(maximum_step_, duration);
    in.minimum_substep_s = minimum_step_;
    in.relative_tolerance = relative_tolerance_;
    const auto result = integrate_parcel_events(in);
    if (!result.available)
      return result.status == ParcelEventsStatus::provider_failure
                 ? unavailable()
                 : invalid(10270 + unsigned(result.status));
    const std::uint64_t prefix = std::uint64_t(index) << 32;
    for (std::size_t i = 0; i < result.segment_count; ++i) {
      const auto &segment = result.segments[i];
      spray::Vector3 midpoint{};
      for (unsigned d = 0; d < 3; ++d)
        midpoint[d] =
            segment.begin.position_m[d] +
            .5 * (segment.end.position_m[d] - segment.begin.position_m[d]);
      const auto stencil = gas_.stencil(midpoint);
      if (!stencil.succeeded())
        return invalid();
      portable::ExchangeSegment row;
      row.revision = revision;
      row.parcel_id = job.value.parcel.id;
      row.segment_ordinal = prefix | segment.ordinal;
      row.vapor_species_index = asset_.vapor_species_index;
      const auto &e = segment.exchange;
      row.delta = {
          e.parcel_liquid_mass_delta_kg, e.parcel_momentum_delta_kg_m_per_s,
          e.parcel_thermochemical_enthalpy_delta_j,
          e.parcel_kinetic_energy_delta_j, e.thermal_exchange_to_gas_j};
      for (std::size_t j = 0; j < stencil.entry_count; ++j) {
        if (!stencil.entries[j].weight)
          continue;
        row.global_cell = stencil.entries[j].global_cell;
        row.deposition_weight = stencil.entries[j].weight;
        const auto status = append(row);
        if (!status)
          return status;
      }
    }
    if (result.physical_outlet) {
      const auto status = ledger(result.outlet_inventory, job.value, revision,
                                 portable::ExchangeChannel::outlet,
                                 prefix | UINT32_C(0xfffffffd));
      if (!status)
        return status;
    }
    if ((result.wall_exchange.momentum_kg_m_per_s[0] != 0 ||
         result.wall_exchange.momentum_kg_m_per_s[1] != 0 ||
         result.wall_exchange.momentum_kg_m_per_s[2] != 0) ||
        result.wall_exchange.kinetic_energy_j != 0) {
      const auto status = ledger(result.wall_exchange, job.value, revision,
                                 portable::ExchangeChannel::wall,
                                 prefix | UINT32_C(0xfffffffe));
      if (!status)
        return status;
    }
    if (result.breakup_requested) {
      const auto &b = result.breakup_budget;
      const double residual = b.surface_energy_increase_j +
                              b.bulk_kinetic_energy_residual_j -
                              b.supplied_deformation_energy_j;
      if (!std::isfinite(residual) ||
          std::abs(residual) >
              in.energy_absolute_tolerance_j +
                  relative_tolerance_ *
                      (std::abs(b.surface_energy_increase_j) +
                       std::abs(b.bulk_kinetic_energy_residual_j) +
                       std::abs(b.supplied_deformation_energy_j)))
        return invalid();
      breakup_energy_.deformation_consumed_j += b.supplied_deformation_energy_j;
      breakup_energy_.surface_increase_j += b.surface_energy_increase_j;
      breakup_energy_.bulk_kinetic_increase_j +=
          b.bulk_kinetic_energy_residual_j;
      breakup_energy_.residual_j += residual;
      ++breakup_energy_.event_count;
      if (!std::isfinite(breakup_energy_.deformation_consumed_j) ||
          !std::isfinite(breakup_energy_.surface_increase_j) ||
          !std::isfinite(breakup_energy_.bulk_kinetic_increase_j) ||
          !std::isfinite(breakup_energy_.residual_j))
        return invalid();
      if (!result.children.available ||
          result.children.child_count >
              parcel_capacity_ - jobs_.size() - next_.size())
        return capacity();
      for (std::size_t child = 0; child < result.children.child_count;
           ++child) {
        Parcel value{result.children.children[child], 0, 0, 0};
        auto status = audit_id(value.parcel.id, revision);
        if (!status)
          return status;
        if (result.children_remaining_duration_s > 0)
          jobs_.push_back({value, job.start + result.advanced_duration_s,
                           result.children_remaining_duration_s});
        else
          next_.push_back(value);
      }
    } else if (!result.parent_removed) {
      if (next_.size() + jobs_.size() >= parcel_capacity_)
        return capacity();
      next_.push_back({result.parcel, result.auxiliary.tab_deformation,
                       result.auxiliary.tab_deformation_rate_per_s,
                       result.auxiliary.breakup_ordinal});
    }
  }
  return {};
}
Status ProductParcelAdvance::prepare(Span<const Parcel> accepted,
                                     Span<Injector *const> injectors,
                                     portable::Revision revision,
                                     double duration, ConstFieldView rho,
                                     ConstFieldView velocity) noexcept {
  discard();
  segments_.clear();
  ids_.clear();
  current_.clear();
  next_.clear();
  waves_ = 0;
  if (comm_ == MPI_COMM_NULL)
    return invalid();
  std::size_t started_injectors = 0;
  breakup_energy_ = {};
  auto fail = [&](Status s) noexcept {
    discard();
    if (injectors.data)
      for (std::size_t i = 0; i < started_injectors; ++i)
        if (injectors.data[i] && injectors.data[i]->trial_active())
          (void)injectors.data[i]->rollback_trial();
    return s;
  };
  Status local;
  if (!geometry_ || accepted.size > parcel_capacity_ ||
      (accepted.size && !accepted.data) ||
      (injectors.size && !injectors.data) || !std::isfinite(duration) ||
      duration <= 0 || revision.algorithm_version != 1 ||
      !revision.input_revision || !gas_.bound_to(revision, duration) ||
      !valid_cell_view(rho, patch_.cells, 0, 1, 0) ||
      !valid_cell_view(velocity, patch_.cells, 0, 3, 0))
    local = invalid();
  auto status = agree(local);
  if (!status)
    return fail(status);
  std::uint64_t clock[]{revision.accepted_step, revision.input_revision, 0},
      low[3]{}, high[3]{};
  std::memcpy(clock + 2, &duration, 8);
  if (MPI_Allreduce(clock, low, 3, MPI_UINT64_T, MPI_MIN, comm_) !=
          MPI_SUCCESS ||
      MPI_Allreduce(clock, high, 3, MPI_UINT64_T, MPI_MAX, comm_) !=
          MPI_SUCCESS)
    return fail({StatusCode::mpi_failure, 10263});
  if (!std::equal(low, low + 3, high))
    return fail(invalid());
  revision_ = revision;
  local = events_.bind_revision(revision);
  if (film_.bind_revision(revision) != portable::Status::success)
    local = invalid();
  const auto n = geometry_->global_cells();
  double speed = 0;
  std::size_t cell = 0;
  for (int z = 0; z < patch_.cells.z; ++z)
    for (int y = 0; y < patch_.cells.y; ++y)
      for (int x = 0; x < patch_.cells.x; ++x, ++cell) {
        const Int3 p{x, y, z},
            g{x + patch_.begin.x, y + patch_.begin.y, z + patch_.begin.z};
        auto &c = cells_[cell];
        c.global_cell = cell_id(g, n);
        c.volume_m3 = geometry_->x().widths().data[g.x] *
                      geometry_->y().widths().data[g.y] *
                      geometry_->z().widths().data[g.z];
        c.gas_mass_kg = rho.unchecked(p, 0) * c.volume_m3;
        if (!std::isfinite(c.gas_mass_kg) || c.gas_mass_kg <= 0)
          local = invalid();
        for (unsigned d = 0; d < 3; ++d) {
          const double u = velocity.unchecked(p, d);
          c.gas_momentum_kg_m_per_s[d] = c.gas_mass_kg * u;
          if (!std::isfinite(u))
            local = invalid();
          speed = std::max(speed, std::abs(u));
        }
      }
  for (std::size_t i = 0; i < accepted.size && local; ++i) {
    if (accepted.data[i].parcel.liquid_material_fingerprint !=
        asset_.pack.material_fingerprint) {
      local = invalid();
      break;
    }
    next_.push_back(accepted.data[i]);
    local = audit_id(accepted.data[i].parcel.id, revision);
  }
  for (std::size_t i = 0; i < injectors.size && local; ++i) {
    auto *injector = injectors.data[i];
    if (!injector || !injector->configured() || injector->trial_active() ||
        injector->configured_spec().liquid_material_fingerprint !=
            asset_.pack.material_fingerprint) {
      local = invalid();
      break;
    }
    started_injectors = i + 1;
    const auto report = injector->begin_trial(revision.accepted_step, duration);
    if (!report.succeeded() ||
        report.parcel_count > parcel_capacity_ - next_.size()) {
      local = capacity();
      break;
    }
    for (std::size_t j = 0; j < report.parcel_count && local; ++j) {
      Parcel value;
      if (!injector->candidate_at(j, value.parcel)) {
        local = invalid();
        break;
      }
      next_.push_back(value);
      local = audit_id(value.parcel.id, revision);
      if (local)
        local = inventory(value.parcel, revision);
    }
  }
  status = agree(local);
  if (!status)
    return fail(status);
  // This initial checked migration audits IDs before any trajectory work.
  status = move();
  if (!status)
    return fail(status);
  for (const auto &p : current_)
    for (double u : p.parcel.velocity_m_per_s)
      speed = std::max(speed, std::abs(u));
  double global_speed = 0;
  if (MPI_Allreduce(&speed, &global_speed, 1, MPI_DOUBLE, MPI_MAX, comm_) !=
      MPI_SUCCESS)
    return fail({StatusCode::mpi_failure, 10263});
  const double initial_wave =
      global_speed > 0
          ? std::min(maximum_step_, .25 * minimum_width_ / global_speed)
          : maximum_step_;
  double elapsed = 0;
  while (elapsed < duration) {
    if (waves_ == 2048)
      return fail(capacity());
    double dt = std::min(initial_wave, duration - elapsed);
    const auto nsegments = segments_.size(), nids = ids_.size();
    const auto old_breakup = breakup_energy_;
    for (;;) {
      if (elapsed + dt <= elapsed)
        return fail(unavailable());
      segments_.resize(nsegments);
      ids_.resize(nids);
      breakup_energy_ = old_breakup;
      status = agree(wave(revision, elapsed, dt, waves_));
      if (status)
        break;
      if (status.detail != 10262 || dt / 2 < minimum_step_)
        return fail(status);
      dt /= 2;
    }
    status = move();
    if (!status)
      return fail(status);
    elapsed += dt;
    ++waves_;
  }
  auto routed = route_.prepare(revision, ids_.data(), ids_.size());
  if (routed != portable::Status::success)
    return fail(invalid(10300 + unsigned(routed)));
  routed = route_.prepare(revision, segments_.data(), segments_.size());
  if (routed != portable::Status::success)
    return fail(invalid(10300 + unsigned(routed)));
  exchange_ =
      batch_->evaluate_routed(revision, cells_.data(), cells_.size(),
                              route_.candidates(), 1e-16, relative_tolerance_);
  status =
      agree(exchange_.available ? Status{}
                                : invalid(10320 + unsigned(exchange_.status)));
  if (!status)
    return fail(status);
  available_ = true;
  return {};
}
} // namespace hundun::v04::detail

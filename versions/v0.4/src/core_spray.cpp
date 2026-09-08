// SPDX-License-Identifier: Apache-2.0
#include "core_spray_detail.hpp"
#include <climits>

namespace hundun::v04::detail {
Status ProductSpray::agree(Status local) const noexcept {
  int failing = local ? INT_MAX : rank_, first = 0;
  if (MPI_Allreduce(&failing, &first, 1, MPI_INT, MPI_MIN, comm_) !=
      MPI_SUCCESS)
    return {StatusCode::mpi_failure, 10341};
  if (first == INT_MAX)
    return {};
  std::uint32_t words[]{static_cast<std::uint32_t>(local.code), local.detail};
  if (MPI_Bcast(words, 2, MPI_UINT32_T, first, comm_) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, 10341};
  return {static_cast<StatusCode>(words[0]), words[1]};
}
Status ProductSpray::configure_local(const ValidatedModel &model,
                                     const std::filesystem::path &root,
                                     const ProductReactionSources &reaction,
                                     const CartesianGeometryPlan &geometry,
                                     MeshPatch patch, int rank,
                                     RestartCellRecordsView tcr) {
  if (!model.spray)
    return {};
  if (enabled() || !reaction.gas_query() ||
      model.time.scheme != TimeScheme::backward_euler)
    return invalid();
  std::array<bool, 6> walls{};
  for (unsigned d = 0; d < 3; ++d) {
    periodic_[d] = model.boundaries[2 * d].flow_kind == BoundaryKind::periodic;
    if (periodic_[d] !=
        (model.boundaries[2 * d + 1].flow_kind == BoundaryKind::periodic))
      return invalid();
  }
  for (unsigned face = 0; face < 6; ++face) {
    switch (model.boundaries[face].flow_kind) {
    case BoundaryKind::no_slip_wall:
    case BoundaryKind::slip:
    case BoundaryKind::symmetry:
      walls[face] = true;
      break;
    case BoundaryKind::moving_wall:
      // The current event law has a stationary wall energy/impulse ledger.
      return invalid();
    default:
      break;
    }
  }
  if (model.immersed_boundary)
    return invalid();
  spec_ = *model.spray;
  maximum_bytes_ = model.mesh.limits.max_memory_bytes_per_rank;
  const auto ns = reaction.gas_identity().species_names.size();
  if (ns < 2 || ns > 252 || spec_.injectors.empty() ||
      spec_.injectors.size() > 64)
    return invalid();
  // A bounded asset contains <=65536 bytes. Reserve an additional MiB for its
  // parsed identity/strings and the cold control vectors; no solver arena is
  // hidden in this allowance. Injector SoA and ID scratch are charged per row.
  local_bytes_ =
      sizeof(*this) + UINT64_C(1048576) +
      spec_.injectors.size() *
          (sizeof(ProductSprayHistory::Injector) +
           std::uint64_t(spec_.maximum_local_parcels) *
               (sizeof(spray::SprayParcelState) + sizeof(spray::ParcelId))) +
      ns * (sizeof(std::size_t) + 2 * sizeof(ConstFieldView) +
            8 * sizeof(double));
  if (local_bytes_ > maximum_bytes_)
    return {StatusCode::allocation_failure, 10342};
  auto loaded = spray::detail::load_liquid_asset(root / spec_.liquid_file,
                                                 spec_.liquid_fingerprint,
                                                 reaction.gas_identity());
  if (!loaded.available)
    return invalid();
  asset_ = std::move(loaded.asset);
  geometry_ = &geometry;
  patch_ = patch;
  rank_ = rank;
  auto status = gas_.configure(geometry, patch, periodic_,
                               reaction.gas_identity().composition_fingerprint,
                               reaction.species_indices(),
                               reaction.dependent_index(), true);
  if (status)
    status = events_.configure(geometry, periodic_, walls);
  if (!status)
    return status;
  spray::detail::LiquidPropertyService liquid(&asset_.pack, 1);
  for (const auto &input : spec_.injectors) {
    spray::detail::ParcelLocation location;
    const spray::Vector3 position{input.origin_m.x, input.origin_m.y,
                                  input.origin_m.z};
    status = gas_.locate(position, {0, 1, 1}, location);
    if (!status)
      return status;
    if (location.owner_rank != rank)
      continue;
    const auto properties = liquid.evaluate(
        {asset_.pack.material_fingerprint, input.temperature_k});
    if (!properties.succeeded())
      return invalid();
    auto injector = std::make_unique<ProductSprayHistory::Injector>();
    spray::detail::InjectorSpec injection;
    injection.seed = spec_.seed;
    injection.injector_id = input.id;
    injection.shape = input.cone_half_angle_rad == 0
                          ? spray::detail::InjectorShape::point
                          : spray::detail::InjectorShape::cone;
    injection.origin_m = position;
    injection.axis = {input.axis.x, input.axis.y, input.axis.z};
    injection.cone_half_angle_rad = input.cone_half_angle_rad;
    injection.injection_speed_m_per_s = input.speed_m_per_s;
    injection.mass_flow_rate_kg_per_s = input.mass_flow_rate_kg_per_s;
    injection.represented_mass_per_parcel_kg =
        input.represented_mass_per_parcel_kg;
    injection.droplet_diameter_m = input.droplet_diameter_m;
    injection.droplet_mass_kg = properties.properties.density_kg_per_m3 *
                                std::acos(-1.) / 6 * input.droplet_diameter_m *
                                input.droplet_diameter_m *
                                input.droplet_diameter_m;
    injection.temperature_k = input.temperature_k;
    injection.liquid_material_fingerprint = asset_.pack.material_fingerprint;
    injection.owner_global_cell = location.global_cell;
    status = injector->reserve(spec_.maximum_local_parcels);
    if (status)
      status = injector->configure(injection);
    if (!status)
      return status;
    injector_views_.push_back(injector.get());
    injectors_.push_back(std::move(injector));
  }
  PlanFingerprint identity =
      (model.fingerprint ^ UINT64_C(0x5350524159000001)) *
      UINT64_C(1099511628211);
  if (!identity)
    identity = 1;
  status = history.configure(identity, patch, geometry.global_cells(),
                             spec_.maximum_local_parcels,
                             asset_.pack.material_fingerprint,
                             {injector_views_.data(), injector_views_.size()},
                             tcr, maximum_bytes_ - local_bytes_);
  if (!status)
    return status;
  local_bytes_ += history.owned_bytes();
  film_ = std::make_unique<spray::detail::FilmEnvironmentBridge>(
      asset_, *reaction.gas_query(), gas_, portable::Revision{0, 1, 1});
  advance_ =
      std::make_unique<ProductParcelAdvance>(gas_, events_, asset_, *film_);
  halo_views_.resize(ns + 2);
  species_views_.resize(ns - 1);
  const auto n = geometry.global_cells();
  // Gather the entire accepted sampling fringe, including U/p faces which
  // are not part of the thermophysical predictor's primary field Halo.
  for (int z = -2; z < patch.cells.z + 2; ++z)
    for (int y = -2; y < patch.cells.y + 2; ++y)
      for (int x = -2; x < patch.cells.x + 2; ++x) {
        if (x >= 0 && x < patch.cells.x && y >= 0 && y < patch.cells.y &&
            z >= 0 && z < patch.cells.z)
          continue;
        auto wrap = [](int v, int count) {
          const int r = v % count;
          return r < 0 ? r + count : r;
        };
        int index[]{x + patch.begin.x, y + patch.begin.y, z + patch.begin.z};
        const int extent[]{n.x, n.y, n.z};
        bool physical_ghost = false;
        for (unsigned d = 0; d < 3; ++d) {
          if (periodic_[d])
            index[d] = wrap(index[d], extent[d]);
          else if (index[d] < 0 || index[d] >= extent[d])
            physical_ghost = true;
        }
        // Sampling/deposition use one-sided interior cell values at physical
        // boundaries. Preserve the native gas solver's physical ghost closure.
        if (physical_ghost)
          continue;
        const Int3 g{index[0], index[1], index[2]};
        halo_cells_.push_back(std::uint64_t(g.x) +
                              std::uint64_t(n.x) * (std::uint64_t(g.y) +
                                                    std::uint64_t(n.y) * g.z));
        halo_indices_.push_back({x, y, z});
      }
  identity_ = identity;
  return {};
}
Status
ProductSpray::configure_collective(MPI_Comm comm,
                                   Span<const RemoteDonorFieldSpec> fields) {
  if (!enabled())
    return {};
  comm_ = comm;
  Status status =
      agree(fields.size == halo_views_.size() ? Status{} : invalid());
  if (!status)
    return status;
  const RemoteDonorTargets targets{
      geometry_->fingerprint(),
      2,
      periodic_,
      {halo_cells_.data(), halo_cells_.size()},
      {halo_indices_.data(), halo_indices_.size()}};
  status = RemoteDonorExchangePlan::analyze_cells(
      comm, geometry_->global_cells(), patch_, targets, fields, 10, halo_);
  if (!status)
    return status;
  const auto stats = halo_.stats();
  // Conservative bound for donor metadata, requests and scalar send/receive
  // buffers, including target aliases retained by the prepared plan.
  const std::uint64_t fringe =
      (stats.received_cells + stats.supplied_cells) * 256 +
      2 * stats.bytes_per_exchange +
      halo_cells_.capacity() * sizeof(GlobalCellId) +
      halo_indices_.capacity() * sizeof(Int3);
  const auto np = spec_.maximum_local_parcels,
             ns = spec_.maximum_local_segments;
  int ranks = 0;
  if (MPI_Comm_size(comm, &ranks) != MPI_SUCCESS)
    return {StatusCode::mpi_failure, 10341};
  const std::uint64_t restore =
      std::uint64_t(ranks) * 5 * sizeof(int) +
      std::uint64_t(np) * (sizeof(int) + 36 * sizeof(std::uint64_t) +
                           sizeof(ProductSprayHistory::Parcel)) +
      std::uint64_t(ns) *
          (2 * sizeof(std::uint64_t) + sizeof(spray::ParcelId)) +
      24;
  status = agree(local_bytes_ > maximum_bytes_ ||
                         fringe > maximum_bytes_ - local_bytes_ ||
                         restore > maximum_bytes_ - local_bytes_ - fringe
                     ? Status{StatusCode::allocation_failure, 10342}
                     : Status{});
  if (!status)
    return status;
  local_bytes_ += fringe + restore;
  status = halo_.bind(comm);
  if (status)
    status = restore_migration_.configure(comm, geometry_->global_cells(),
                                          patch_, np, ns);
  if (status)
    status = advance_->configure(comm, *geometry_, patch_, spec_,
                                 maximum_bytes_ - local_bytes_);
  if (status)
    owned_bytes_ = local_bytes_ + advance_->owned_bytes();
  return status;
}
Status ProductSpray::configure_source_transport(
    Span<const RemoteDonorFieldSpec> fields) {
  if (!enabled() || source_halo_bound_ || fields.size == 0)
    return invalid();
  const RemoteDonorTargets targets{
      geometry_->fingerprint(),
      2,
      periodic_,
      {halo_cells_.data(), halo_cells_.size()},
      {halo_indices_.data(), halo_indices_.size()}};
  auto status = RemoteDonorExchangePlan::analyze_cells(
      comm_, geometry_->global_cells(), patch_, targets, fields, 10,
      source_halo_);
  if (!status)
    return status;
  const auto stats = source_halo_.stats();
  status = agree(stats.peer_messages > UINT32_MAX - halo_.stats().peer_messages
                     ? invalid() : Status{});
  if (!status)
    return status;
  const auto bytes = (stats.received_cells + stats.supplied_cells) * 256 +
                     2 * stats.bytes_per_exchange;
  status = agree(owned_bytes_ > maximum_bytes_ ||
                         bytes > maximum_bytes_ - owned_bytes_
                     ? Status{StatusCode::allocation_failure, 10342}
                     : Status{});
  if (status)
    status = source_halo_.bind(comm_);
  if (status) {
    owned_bytes_ += bytes;
    source_halo_bound_ = true;
  }
  return status;
}
Status ProductSpray::exchange_source_transport(
    Span<FieldView> fields, RevisionToken boundary,
    Span<ThermophysicalGhostAuthority> ghosts) noexcept {
  auto status = agree(source_halo_bound_ && ghosts.size == fields.size &&
                              ghosts.data && boundary
                          ? Status{}
                          : invalid());
  if (status)
    status =
        agree(source_halo_.preflight_exchange(10, {fields.data, fields.size}));
  if (status)
    status = source_halo_.exchange(10, fields);
  if (!status)
    return status;
  for (std::size_t i = 0; i < fields.size; ++i) {
    const auto f = fields.data[i];
    ghosts.data[i] = {reinterpret_cast<std::uintptr_t>(&source_halo_),
                      f.field,
                      f.revision,
                      f.storage_identity,
                      f.revision_domain,
                      geometry_->topology_revision(),
                      boundary,
                      2};
  }
  return {};
}
Status ProductSpray::prepare(portable::Revision revision, double duration,
                             double pressure_reference, ConstFieldView density,
                             FieldView pressure, FieldView enthalpy,
                             FieldView velocity,
                             Span<const FieldView> independent) noexcept {
  if (!enabled())
    return {};
  discard();
  auto status =
      agree(independent.size == species_views_.size() && independent.data
                ? Status{}
                : invalid());
  if (!status)
    return status;
  halo_views_[0] = pressure;
  halo_views_[1] = enthalpy;
  halo_views_[2] = velocity;
  for (std::size_t i = 0; i < independent.size; ++i)
    halo_views_[i + 3] = independent.data[i];
  status = agree(
      halo_.preflight_exchange(10, {halo_views_.data(), halo_views_.size()}));
  if (status)
    status = halo_.exchange(10, {halo_views_.data(), halo_views_.size()});
  if (!status)
    return status;
  for (std::size_t i = 0; i < independent.size; ++i)
    species_views_[i] = as_const(halo_views_[i + 3]);
  status = agree(gas_.bind(revision, duration, pressure_reference,
                           as_const(halo_views_[0]), as_const(halo_views_[1]),
                           as_const(halo_views_[2]),
                           {species_views_.data(), species_views_.size()}));
  if (status)
    status = advance_->prepare(history.accepted_parcels(),
                               {injector_views_.data(), injector_views_.size()},
                               revision, duration, density,
                               as_const(halo_views_[2]));
  if (status)
    revision_ = revision;
  return status;
}
Status ProductSpray::stage(RestartCellRecordsView tcr) noexcept {
  if (!enabled())
    return {};
  return agree(
      history.stage_next(advance_->parcels(), revision_.accepted_step, tcr));
}
Status ProductSpray::stage_restore(const RestartImage &image) noexcept {
  if (!enabled())
    return invalid();
  discard();
  auto status = agree(history.stage_restore(image));
  if (status) {
    revision_ = {image.step, 1, 1};
    status = agree(gas_.bind_location_revision(revision_));
  }
  return status;
}
Status ProductSpray::validate_restored_positions() noexcept {
  auto values = history.prepared_parcels();
  Status local;
  spray::detail::LiquidPropertyService liquid(&asset_.pack, 1);
  for (std::size_t i = 0; i < values.size; ++i) {
    const auto &p = values.data[i].parcel;
    const auto properties =
        liquid.evaluate({p.liquid_material_fingerprint, p.temperature_k});
    if (p.liquid_material_fingerprint != asset_.pack.material_fingerprint ||
        !properties.succeeded()) {
      local = invalid();
      break;
    }
    const double expected_mass = properties.properties.density_kg_per_m3 *
                                 std::acos(-1.) / 6 * p.droplet_diameter_m *
                                 p.droplet_diameter_m * p.droplet_diameter_m;
    if (!std::isfinite(expected_mass) || expected_mass <= 0 ||
        std::abs(expected_mass - p.droplet_mass_kg) >
            1e-10 * std::max(expected_mass, p.droplet_mass_kg)) {
      local = invalid();
      break;
    }
  }
  auto status = agree(local);
  if (!status)
    return status;
  spray::detail::ParcelMigrationReport report;
  return restore_migration_.prepare_checked(values, gas_, revision_, report);
}
} // namespace hundun::v04::detail

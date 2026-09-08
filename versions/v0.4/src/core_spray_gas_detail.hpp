// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "models_spray_film_bridge_detail.hpp"
#include "models_spray_mechanics_detail.hpp"
#include "models_spray_migration_detail.hpp"
#include "solver_cartesian_detail.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace hundun::v04::detail {
// Borrows one accepted native MeanState/Halo snapshot. Every query uses the
// same global stencil that the exchange depositor consumes. No hot allocation
// or per-parcel MPI call; a trajectory beyond available Halo is unavailable.
class ProductParcelGas final : public spray::detail::ParcelGasStateProvider,
                               public spray::detail::ParcelLocationProvider {
public:
  Status configure(const CartesianGeometryPlan &geometry, MeshPatch patch,
                   std::array<bool, 3> periodic, PlanFingerprint composition,
                   Span<const std::size_t> independent, std::size_t dependent,
                   bool trial_boundary_continuation = false) {
    const auto ns = independent.size + 1;
    if (!geometry.fingerprint() || !composition || !independent.data ||
        independent.size == 0 || ns > 255 || dependent >= ns)
      return invalid();
    std::vector<bool> seen(ns);
    seen[dependent] = true;
    for (std::size_t i = 0; i < independent.size; ++i) {
      if (independent.data[i] >= ns || seen[independent.data[i]])
        return invalid();
      seen[independent.data[i]] = true;
    }
    surface_ = nullptr;
    fluid_mask_ = {};
    geometry_ = &geometry;
    minimum_width_ = std::numeric_limits<double>::infinity();
    for (unsigned d = 0; d < 3; ++d) {
      const auto faces = geometry.axis(static_cast<CartesianAxis>(d)).faces();
      for (std::size_t i = 1; i < faces.size; ++i)
        minimum_width_ =
            std::min(minimum_width_, faces.data[i] - faces.data[i - 1]);
    }
    patch_ = patch;
    periodic_ = periodic;
    trial_boundary_continuation_ = trial_boundary_continuation;
    composition_ = composition;
    dependent_ = dependent;
    indices_.assign(independent.data, independent.data + independent.size);
    species_.resize(independent.size);
    scratch_.resize(ns);
    return {};
  }
  // Borrow the immutable native topology mask after its complete donor
  // fringe has been exchanged. Its 0/1 values describe cell centres, while
  // visibility against the actual surface constrains physical query points.
  Status bind_immersed(const ImmersedSurfacePlan &surface, ConstFieldView mask,
                       PlanFingerprint geometry) noexcept {
    if (!geometry_ || geometry != geometry_->fingerprint() ||
        !surface.fingerprint() || !valid_cell_view(mask, patch_.cells, 0, 1, 0))
      return invalid();
    for (int z = -mask.ghosts.z; z < patch_.cells.z + mask.ghosts.z; ++z)
      for (int y = -mask.ghosts.y; y < patch_.cells.y + mask.ghosts.y; ++y)
        for (int x = -mask.ghosts.x; x < patch_.cells.x + mask.ghosts.x; ++x) {
          const double value = mask.unchecked({x, y, z}, 0);
          if (value != 0 && value != 1)
            return invalid();
        }
    surface_ = &surface;
    fluid_mask_ = mask;
    return {};
  }
  Status bind(portable::Revision revision, double duration,
              double pressure_reference, ConstFieldView pressure,
              ConstFieldView enthalpy, ConstFieldView velocity,
              Span<const ConstFieldView> independent) noexcept {
    available_ = false;
    if (!geometry_ || revision.algorithm_version != 1 ||
        !revision.input_revision || !std::isfinite(duration) || duration <= 0 ||
        !std::isfinite(pressure_reference) || pressure_reference <= 0 ||
        independent.size != indices_.size() || !independent.data ||
        !valid_cell_view(pressure, patch_.cells, 0, 1, 0) ||
        !valid_cell_view(enthalpy, patch_.cells, 0, 1, 0) ||
        !valid_cell_view(velocity, patch_.cells, 0, 3, 0))
      return invalid();
    for (std::size_t i = 0; i < independent.size; ++i)
      if (!valid_cell_view(independent.data[i], patch_.cells, 0, 1, 0))
        return invalid();
    revision_ = revision;
    duration_ = duration;
    pressure_reference_ = pressure_reference;
    pressure_ = pressure;
    enthalpy_ = enthalpy;
    velocity_ = velocity;
    std::copy(independent.data, independent.data + independent.size,
              species_.begin());
    available_ = true;
    return {};
  }
  Status bind_location_revision(portable::Revision revision) noexcept {
    if (revision.algorithm_version != 1 || !revision.input_revision)
      return invalid();
    available_ = false;
    revision_ = revision;
    return {};
  }
  bool bound_to(portable::Revision revision, double duration) const noexcept {
    return available_ && revision_ == revision && duration_ == duration;
  }
  spray::detail::ParcelGridCouplingStencil
  stencil(spray::Vector3 position) const noexcept {
    auto result = geometry_ ? spray::detail::build_parcel_grid_coupling_stencil(
                                  *geometry_, patch_, position, periodic_)
                            : spray::detail::ParcelGridCouplingStencil{};
    if (!surface_ || !result.succeeded())
      return result;
    canonical_position(position);
    unsigned count = 0, largest = 0;
    double sum = 0;
    for (unsigned i = 0; i < result.entry_count; ++i) {
      auto entry = result.entries[i];
      if (entry.weight == 0)
        continue;
      Int3 local{};
      if (!local_index(fluid_mask_, entry.global_index, local))
        return {};
      if (fluid_mask_.unchecked(local, 0) == 0)
        continue;
      const int index[]{entry.global_index.x, entry.global_index.y,
                        entry.global_index.z};
      spray::Vector3 centre{};
      for (unsigned d = 0; d < 3; ++d) {
        const auto &axis = geometry_->axis(static_cast<CartesianAxis>(d));
        centre[d] = axis.centres().data[index[d]];
        if (periodic_[d]) {
          const auto faces = axis.faces();
          const double length = faces.data[faces.size - 1] - faces.data[0];
          centre[d] += std::round((position[d] - centre[d]) / length) * length;
        }
      }
      SurfaceSegmentIntersection hit;
      if (centre != position &&
          !surface_->first_segment_intersection(
              {centre[0], centre[1], centre[2]},
              {position[0], position[1], position[2]}, hit))
        return {};
      // A query on the surface is permitted. A donor across the wall is not.
      if (hit.triangle != kInvalidSurfaceTriangle &&
          hit.segment_fraction <
              1 - 4096 * std::numeric_limits<double>::epsilon())
        continue;
      result.entries[count] = entry;
      if (count == 0 || entry.weight > result.entries[largest].weight)
        largest = count;
      sum += entry.weight;
      ++count;
    }
    if (count == 0 || !(sum > 0) || !std::isfinite(sum))
      return {};
    result.entry_count = count;
    double normalized = 0;
    for (unsigned i = 0; i < count; ++i) {
      result.entries[i].weight /= sum;
      normalized += result.entries[i].weight;
    }
    result.entries[largest].weight += 1 - normalized;
    return result;
  }
  spray::detail::ParcelGasSample
  sample(const spray::SprayParcelState &parcel, double elapsed,
         spray::detail::ParcelPass pass, portable::Revision revision,
         double *full_y, std::size_t capacity) const noexcept override {
    spray::detail::ParcelGasSample out;
    if (revision != revision_) {
      out.status = portable::Status::stale_revision;
      return out;
    }
    if (!available_ || !full_y || capacity < scratch_.size() ||
        !std::isfinite(elapsed) || elapsed < 0 || elapsed > duration_ ||
        (pass != spray::detail::ParcelPass::predictor &&
         pass != spray::detail::ParcelPass::corrector))
      return out;
    // Event-aware integration probes its unconstrained endpoint before it can
    // truncate the interval at a wall/outlet. Continue only that PH query by
    // the existing one-sided end-cell value, at most one boundary cell wide.
    // Public deposition stencils and owner location remain domain constrained.
    auto position = parcel.position_m;
    if (trial_boundary_continuation_)
      for (unsigned d = 0; d < 3; ++d) {
        if (periodic_[d])
          continue;
        const auto faces =
            geometry_->axis(static_cast<CartesianAxis>(d)).faces();
        const double low = faces.data[0], high = faces.data[faces.size - 1];
        if (!std::isfinite(position[d]) ||
            position[d] < low - (faces.data[1] - low) ||
            position[d] > high + (high - faces.data[faces.size - 2]))
          return out;
        position[d] = std::clamp(position[d], low, high);
      }
    auto weights = stencil(position);
    if (!weights.succeeded() && surface_ && trial_boundary_continuation_) {
      canonical_position(position);
      ClosestSurfacePoint point;
      if (!surface_->closest_point({position[0], position[1], position[2]},
                                   point) ||
          point.triangle == kInvalidSurfaceTriangle ||
          point.squared_distance > minimum_width_ * minimum_width_)
        return out;
      weights = stencil({point.point.x, point.point.y, point.point.z});
    }
    if (!weights.succeeded())
      return out;
    std::fill(scratch_.begin(), scratch_.end(), 0);
    for (std::size_t i = 0; i < weights.entry_count; ++i) {
      const auto &entry = weights.entries[i];
      if (entry.weight == 0)
        continue;
      Int3 p{}, h{}, u{};
      if (!local_index(pressure_, entry.global_index, p) ||
          !local_index(enthalpy_, entry.global_index, h) ||
          !local_index(velocity_, entry.global_index, u))
        return out;
      out.pressure_pa += entry.weight * pressure_.unchecked(p, 0);
      out.enthalpy_j_per_kg += entry.weight * enthalpy_.unchecked(h, 0);
      for (unsigned c = 0; c < 3; ++c)
        out.velocity_m_per_s[c] += entry.weight * velocity_.unchecked(u, c);
      for (std::size_t j = 0; j < species_.size(); ++j) {
        Int3 cell{};
        if (!local_index(species_[j], entry.global_index, cell))
          return out;
        scratch_[indices_[j]] += entry.weight * species_[j].unchecked(cell, 0);
      }
    }
    double independent_sum = 0;
    for (double y : scratch_)
      independent_sum += y;
    scratch_[dependent_] = 1 - independent_sum;
    out.pressure_pa += pressure_reference_;
    if (!std::isfinite(out.pressure_pa) || out.pressure_pa <= 0 ||
        !std::isfinite(out.enthalpy_j_per_kg))
      return out;
    for (double u : out.velocity_m_per_s)
      if (!std::isfinite(u))
        return out;
    for (double y : scratch_)
      if (!std::isfinite(y) || y < 0 || y > 1)
        return out;
    out.revision = revision;
    out.composition_fingerprint = composition_;
    out.species_count = scratch_.size();
    out.status = portable::Status::success;
    std::copy(scratch_.begin(), scratch_.end(), full_y);
    return out;
  }
  Status locate(const spray::Vector3 &position, portable::Revision revision,
                spray::detail::ParcelLocation &out) const noexcept override {
    if (!geometry_ || revision != revision_)
      return {StatusCode::invalid_plan,
              static_cast<std::uint32_t>(
                  spray::detail::MigrationDetail::stale_revision)};
    if (!stencil(position).succeeded())
      return invalid();
    // A fluid physical point can occupy a cell whose centre is solid. Owner
    // routing follows physical cell faces, not the filtered donor list.
    const auto weights = spray::detail::build_parcel_grid_coupling_stencil(
        *geometry_, patch_, position, periodic_);
    std::array<int, 3> cell{};
    for (unsigned d = 0; d < 3; ++d) {
      const auto faces = geometry_->axis(static_cast<CartesianAxis>(d)).faces();
      const double low = faces.data[0], high = faces.data[faces.size - 1];
      double x = position[d];
      if (periodic_[d] && (x < low || x >= high)) {
        double offset = std::fmod(x - low, high - low);
        if (offset < 0)
          offset += high - low;
        x = low + offset;
        if (x >= high)
          x = low;
      }
      // A rebound can end exactly on the upper physical wall; its owner is
      // the adjacent interior cell. Outlet parcels are removed before locate.
      if (!periodic_[d] && x == high)
        x = std::nextafter(high, low);
      if (!std::isfinite(x) || x < low || x >= high)
        return invalid();
      cell[d] = int(std::upper_bound(faces.data, faces.data + faces.size, x) -
                    faces.data) -
                1;
    }
    for (std::size_t i = 0; i < weights.entry_count; ++i) {
      const auto &entry = weights.entries[i];
      if (entry.global_index.x == cell[0] && entry.global_index.y == cell[1] &&
          entry.global_index.z == cell[2]) {
        out = {revision, entry.global_cell, entry.owner_rank};
        return {};
      }
    }
    return invalid();
  }

private:
  void canonical_position(spray::Vector3 &position) const noexcept {
    for (unsigned d = 0; d < 3; ++d) {
      if (!periodic_[d])
        continue;
      const auto faces = geometry_->axis(static_cast<CartesianAxis>(d)).faces();
      const double low = faces.data[0], high = faces.data[faces.size - 1];
      if (position[d] < low || position[d] >= high) {
        double offset = std::fmod(position[d] - low, high - low);
        if (offset < 0)
          offset += high - low;
        position[d] = low + offset;
        if (position[d] >= high)
          position[d] = low;
      }
    }
  }
  bool local_index(ConstFieldView view, Int3 global, Int3 &out) const noexcept {
    const auto n = geometry_->global_cells();
    int index[]{global.x - patch_.begin.x, global.y - patch_.begin.y,
                global.z - patch_.begin.z};
    const int extent[]{patch_.cells.x, patch_.cells.y, patch_.cells.z};
    const int ghosts[]{view.ghosts.x, view.ghosts.y, view.ghosts.z};
    const int global_extent[]{n.x, n.y, n.z};
    for (unsigned d = 0; d < 3; ++d) {
      if (periodic_[d]) {
        if (index[d] < -ghosts[d])
          index[d] += global_extent[d];
        else if (index[d] >= extent[d] + ghosts[d])
          index[d] -= global_extent[d];
      }
      if (index[d] < -ghosts[d] || index[d] >= extent[d] + ghosts[d])
        return false;
    }
    out = {index[0], index[1], index[2]};
    return true;
  }
  static Status invalid() noexcept { return {StatusCode::invalid_plan, 10240}; }
  const CartesianGeometryPlan *geometry_{};
  MeshPatch patch_{};
  const ImmersedSurfacePlan *surface_{};
  ConstFieldView fluid_mask_{};
  double minimum_width_{};
  std::array<bool, 3> periodic_{};
  PlanFingerprint composition_{};
  std::size_t dependent_{};
  portable::Revision revision_{0, 1, 1};
  double duration_{}, pressure_reference_{};
  ConstFieldView pressure_{}, enthalpy_{}, velocity_{};
  std::vector<std::size_t> indices_;
  std::vector<ConstFieldView> species_;
  mutable std::vector<double> scratch_;
  bool available_{}, trial_boundary_continuation_{};
};
} // namespace hundun::v04::detail

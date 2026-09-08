// SPDX-License-Identifier: Apache-2.0
#pragma once
#include "models_spray_events_detail.hpp"
#include "models_spray_mechanics_detail.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace hundun::v04::detail {
// Native fixed-geometry query on the adaptive integrator's trajectory segment.
// Positions stay unwrapped at periodic crossings. The gas sampler and final
// owner locator use the same period; no discontinuous jump enters exchange.
class ProductParcelGeometry final
    : public spray::detail::ParcelEventGeometryProvider {
public:
  Status
  configure(const CartesianGeometryPlan &geometry, std::array<bool, 3> periodic,
            std::array<bool, 6> wall_faces = {},
            const ImmersedSurfacePlan *surface = nullptr,
            ImmersedFluidSide side = ImmersedFluidSide::outside) noexcept {
    if (!geometry.fingerprint() || (surface && !surface->fingerprint()) ||
        (side != ImmersedFluidSide::outside &&
         side != ImmersedFluidSide::inside))
      return invalid();
    const auto n = geometry.global_cells();
    if (n.x <= 0 || n.y <= 0 || n.z <= 0)
      return invalid();
    std::uint64_t count = 1;
    for (int extent : {n.x, n.y, n.z}) {
      if (count > UINT64_MAX / 12 / std::uint64_t(extent))
        return invalid();
      count *= extent;
    }
    if (surface && surface->triangles().size >= (UINT64_C(1) << 63))
      return invalid();
    for (unsigned d = 0; d < 3; ++d)
      if (periodic[d] && (wall_faces[2 * d] || wall_faces[2 * d + 1]))
        return invalid();
    geometry_ = &geometry;
    periodic_ = periodic;
    walls_ = wall_faces;
    surface_ = surface;
    side_ = side;
    count_ = count;
    return {};
  }
  Status bind_revision(portable::Revision revision) noexcept {
    if (revision.algorithm_version != 1 || !revision.input_revision)
      return invalid();
    revision_ = revision;
    return {};
  }
  spray::detail::ParcelEventQueryReport
  query(const spray::SprayParcelState &begin,
        const spray::SprayParcelState &end, double t0, double t1,
        spray::detail::ParcelPass pass,
        portable::Revision revision) const noexcept override {
    using namespace spray::detail;
    ParcelEventQueryReport out;
    out.revision = revision_;
    if (!geometry_ || revision != revision_ || !std::isfinite(t0) ||
        !std::isfinite(t1) || t0 < 0 || t1 <= t0 ||
        begin.owner_global_cell >= count_ ||
        (pass != ParcelPass::predictor && pass != ParcelPass::corrector))
      return out;
    const auto n = geometry_->global_cells();
    const int extent[]{n.x, n.y, n.z};
    auto id = begin.owner_global_cell;
    int cell[3];
    for (unsigned d = 0; d < 3; ++d) {
      cell[d] = int(id % extent[d]);
      id /= extent[d];
    }
    for (unsigned d = 0; d < 3; ++d) {
      const double a = begin.position_m[d], b = end.position_m[d];
      if (!std::isfinite(a) || !std::isfinite(b) ||
          !std::isfinite(begin.velocity_m_per_s[d]))
        return out;
      const auto faces = geometry_->axis(static_cast<CartesianAxis>(d)).faces();
      double low = faces.data[cell[d]], high = faces.data[cell[d] + 1];
      if (periodic_[d]) {
        const double length = faces.data[faces.size - 1] - faces.data[0];
        const double shift =
            std::round((a - (low + .5 * (high - low))) / length) * length;
        low += shift;
        high += shift;
      }
      const double tolerance =
          64 * std::numeric_limits<double>::epsilon() *
          std::max({1., std::abs(a), std::abs(low), std::abs(high)});
      if (a < low - tolerance || a > high + tolerance)
        return out;
      const double dx = b - a;
      if (dx == 0)
        continue;
      const bool upper = dx > 0;
      const double face = upper ? high : low;
      if (upper ? b < face : b > face)
        continue;
      const double fraction = std::clamp((face - a) / dx, 0., 1.);
      ParcelEvent event;
      event.elapsed_time_s = t0 + fraction * (t1 - t0);
      event.identity = 6 * begin.owner_global_cell + 2 * d + unsigned(upper);
      event.has_contact_position = true;
      for (unsigned c = 0; c < 3; ++c)
        event.contact_position_m[c] =
            begin.position_m[c] +
            fraction * (end.position_m[c] - begin.position_m[c]);
      event.contact_position_m[d] = face;
      int next[3]{cell[0], cell[1], cell[2]};
      next[d] += upper ? 1 : -1;
      if (next[d] < 0 || next[d] >= extent[d]) {
        if (periodic_[d])
          next[d] = upper ? 0 : extent[d] - 1;
        else {
          event.kind = walls_[2 * d + unsigned(upper)]
                           ? ParcelEventKind::wall_collision
                           : ParcelEventKind::physical_outlet;
          event.wall_normal[d] = upper ? -1 : 1;
          out.events[out.count++] = event;
          continue;
        }
      }
      event.kind = ParcelEventKind::internal_cell_crossing;
      event.next_global_cell =
          std::uint64_t(next[0]) +
          std::uint64_t(n.x) * (std::uint64_t(next[1]) +
                                std::uint64_t(n.y) * std::uint64_t(next[2]));
      if (event.next_global_cell != begin.owner_global_cell)
        out.events[out.count++] = event;
    }
    if (surface_) {
      Real3 query_begin{begin.position_m[0], begin.position_m[1],
                        begin.position_m[2]};
      const Real3 query_end{end.position_m[0], end.position_m[1],
                            end.position_m[2]};
      // At most one departure retry. Only the immutable surface query is
      // shifted, never the physical trajectory or an extensive inventory.
      for (unsigned attempt = 0; attempt < 2; ++attempt) {
        SurfaceSegmentIntersection hit;
        if (!surface_->first_segment_intersection(query_begin, query_end, hit))
          return {};
        if (hit.triangle == kInvalidSurfaceTriangle)
          break;
        if (hit.triangle >= surface_->triangles().size ||
            !std::isfinite(hit.segment_fraction) || hit.segment_fraction < 0 ||
            hit.segment_fraction > 1)
          return {};
        const auto normal = surface_->triangles()
                                .data[std::size_t(hit.triangle)]
                                .geometric_outward_normal;
        spray::Vector3 fluid_normal{normal.x, normal.y, normal.z};
        if (side_ == ImmersedFluidSide::inside)
          for (auto &v : fluid_normal)
            v = -v;
        const spray::Vector3 point{hit.point.x, hit.point.y, hit.point.z};
        long double length2 = 0, projection = 0, incoming = 0;
        for (unsigned d = 0; d < 3; ++d) {
          const long double dx = end.position_m[d] - begin.position_m[d];
          length2 += dx * dx;
          projection += (point[d] - begin.position_m[d]) * dx;
        }
        if (length2 == 0)
          break;
        const double fraction =
            std::clamp(double(projection / length2), 0., 1.);
        for (unsigned d = 0; d < 3; ++d)
          incoming += (begin.velocity_m_per_s[d] +
                       fraction * (end.velocity_m_per_s[d] -
                                   begin.velocity_m_per_s[d])) *
                      fluid_normal[d];
        if (incoming >= 0) {
          if (fraction > 4096 * std::numeric_limits<double>::epsilon() ||
              attempt)
            return {};
          const double scale = std::max(
              {1., std::abs(point[0]), std::abs(point[1]), std::abs(point[2])});
          const double shift =
              4096 * std::numeric_limits<double>::epsilon() * scale;
          query_begin = {point[0] + shift * fluid_normal[0],
                         point[1] + shift * fluid_normal[1],
                         point[2] + shift * fluid_normal[2]};
          continue;
        }
        out.events[out.count++] = {ParcelEventKind::wall_collision,
                                   t0 + fraction * (t1 - t0),
                                   (UINT64_C(1) << 63) | hit.triangle,
                                   0,
                                   fluid_normal,
                                   1,
                                   false,
                                   point,
                                   true};
        break;
      }
    }
    out.available = true;
    return out;
  }

private:
  static Status invalid() noexcept { return {StatusCode::invalid_plan, 10250}; }
  const CartesianGeometryPlan *geometry_{};
  const ImmersedSurfacePlan *surface_{};
  std::array<bool, 3> periodic_{};
  std::array<bool, 6> walls_{};
  ImmersedFluidSide side_{ImmersedFluidSide::outside};
  std::uint64_t count_{};
  portable::Revision revision_{0, 1, 1};
};
} // namespace hundun::v04::detail

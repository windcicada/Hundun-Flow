// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "models_spray_mechanics_detail.hpp"

#include "mesh_focus_detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace hundun::v04::spray::detail {
namespace {

constexpr double kRoundoff = 256.0 * std::numeric_limits<double>::epsilon();

bool finite(Vector3 value) noexcept {
  return std::isfinite(value[0U]) && std::isfinite(value[1U]) &&
         std::isfinite(value[2U]);
}

Vector3 add(Vector3 left, Vector3 right) noexcept {
  return {left[0U] + right[0U], left[1U] + right[1U],
          left[2U] + right[2U]};
}

Vector3 subtract(Vector3 left, Vector3 right) noexcept {
  return {left[0U] - right[0U], left[1U] - right[1U],
          left[2U] - right[2U]};
}

Vector3 multiply(double scalar, Vector3 value) noexcept {
  return {scalar * value[0U], scalar * value[1U], scalar * value[2U]};
}

double dot(Vector3 left, Vector3 right) noexcept {
  return left[0U] * right[0U] + left[1U] * right[1U] +
         left[2U] * right[2U];
}

double norm(Vector3 value) noexcept {
  return std::sqrt(dot(value, value));
}

Vector3 interpolate(Vector3 begin, Vector3 end, double fraction) noexcept {
  return add(begin, multiply(fraction, subtract(end, begin)));
}

Vector3 to_vector(Real3 value) noexcept {
  return {value.x, value.y, value.z};
}

bool valid_kind(TrajectoryCandidateKind kind) noexcept {
  return kind == TrajectoryCandidateKind::predictor ||
         kind == TrajectoryCandidateKind::corrector;
}

std::int32_t component(Int3 value, std::size_t axis) noexcept {
  return axis == 0U ? value.x : (axis == 1U ? value.y : value.z);
}

bool valid_patch(Int3 global, const MeshPatch& patch) noexcept {
  if (global.x <= 0 || global.y <= 0 || global.z <= 0 ||
      patch.process_grid.x <= 0 || patch.process_grid.y <= 0 ||
      patch.process_grid.z <= 0 ||
      patch.process_grid.x > global.x || patch.process_grid.y > global.y ||
      patch.process_grid.z > global.z || patch.process_coord.x < 0 ||
      patch.process_coord.y < 0 || patch.process_coord.z < 0 ||
      patch.process_coord.x >= patch.process_grid.x ||
      patch.process_coord.y >= patch.process_grid.y ||
      patch.process_coord.z >= patch.process_grid.z) {
    return false;
  }
  const std::uint64_t px =
      static_cast<std::uint64_t>(patch.process_grid.x);
  const std::uint64_t py =
      static_cast<std::uint64_t>(patch.process_grid.y);
  const std::uint64_t pz =
      static_cast<std::uint64_t>(patch.process_grid.z);
  if (px > std::numeric_limits<std::uint64_t>::max() / py ||
      px * py > std::numeric_limits<std::uint64_t>::max() / pz) {
    return false;
  }
  const std::uint64_t rank_count = px * py * pz;
  if (rank_count == 0U ||
      rank_count > static_cast<std::uint64_t>(
                       std::numeric_limits<std::int32_t>::max())) {
    return false;
  }
  const std::int64_t rank =
      patch.process_coord.x +
      static_cast<std::int64_t>(patch.process_grid.x) *
          (patch.process_coord.y +
           static_cast<std::int64_t>(patch.process_grid.y) *
               patch.process_coord.z);
  if (rank < 0 || static_cast<std::uint64_t>(rank) >= rank_count) {
    return false;
  }
  MeshPatch authoritative{};
  if (!::hundun::v04::detail::make_mesh_patch(
          static_cast<std::int32_t>(rank),
          static_cast<std::int32_t>(rank_count), global, authoritative)) {
    return false;
  }
  return patch.begin.x == authoritative.begin.x &&
         patch.begin.y == authoritative.begin.y &&
         patch.begin.z == authoritative.begin.z &&
         patch.cells.x == authoritative.cells.x &&
         patch.cells.y == authoritative.cells.y &&
         patch.cells.z == authoritative.cells.z &&
         patch.process_grid.x == authoritative.process_grid.x &&
         patch.process_grid.y == authoritative.process_grid.y &&
         patch.process_grid.z == authoritative.process_grid.z &&
         patch.process_coord.x == authoritative.process_coord.x &&
         patch.process_coord.y == authoritative.process_coord.y &&
         patch.process_coord.z == authoritative.process_coord.z;
}

std::int32_t owner_coordinate(std::int32_t global_index,
                              std::int32_t global_cells,
                              std::int32_t partitions) noexcept {
  const std::int64_t quotient = global_cells / partitions;
  const std::int64_t remainder = global_cells % partitions;
  const std::int64_t wider_end = (quotient + 1) * remainder;
  if (global_index < wider_end) {
    return static_cast<std::int32_t>(global_index / (quotient + 1));
  }
  return static_cast<std::int32_t>(
      remainder + (global_index - wider_end) / quotient);
}

bool global_cell_id(Int3 index, Int3 cells, std::uint64_t& out) noexcept {
  if (index.x < 0 || index.y < 0 || index.z < 0 || index.x >= cells.x ||
      index.y >= cells.y || index.z >= cells.z) {
    return false;
  }
  const std::uint64_t x = static_cast<std::uint64_t>(index.x);
  const std::uint64_t y = static_cast<std::uint64_t>(index.y);
  const std::uint64_t z = static_cast<std::uint64_t>(index.z);
  const std::uint64_t nx = static_cast<std::uint64_t>(cells.x);
  const std::uint64_t ny = static_cast<std::uint64_t>(cells.y);
  if (y > (std::numeric_limits<std::uint64_t>::max() - x) / nx) {
    return false;
  }
  const std::uint64_t xy = x + nx * y;
  if (z > (std::numeric_limits<std::uint64_t>::max() - xy) / nx / ny) {
    return false;
  }
  out = xy + nx * ny * z;
  return true;
}

struct AxisStencil {
  std::array<std::int32_t, 2U> indices{};
  std::array<double, 2U> weights{};
  std::uint8_t count{};
  bool clamped{};
};

ParcelGridCouplingStencil stencil_failure(
    ParcelGridStencilStatus status) noexcept {
  ParcelGridCouplingStencil result;
  result.status = status;
  return result;
}

bool build_axis_stencil(const AxisMetrics &metrics, double coordinate,
                        AxisStencil &out, bool periodic) noexcept {
  const Span<const double> centres = metrics.centres();
  if (!std::isfinite(coordinate) || centres.data == nullptr ||
      centres.size == 0U) {
    return false;
  }
  if (periodic) {
    const auto faces = metrics.faces();
    if (!faces.data || faces.size != centres.size + 1)
      return false;
    const double length = faces.data[faces.size - 1] - faces.data[0];
    if (!std::isfinite(length) || length <= 0)
      return false;
    if (centres.size == 1) {
      out.indices[0] = 0;
      out.weights[0] = 1;
      out.count = 1;
      return true;
    }
    const auto last = centres.size - 1;
    if (coordinate < centres.data[0] || coordinate > centres.data[last]) {
      // Translate the last centre to the previous periodic image. Its gap
      // to the first centre is determined by the actual two end-cell widths.
      if (coordinate > centres.data[last])
        coordinate -= length;
      const double left = centres.data[last] - length;
      const double spacing = centres.data[0] - left;
      const double high_weight = (coordinate - left) / spacing;
      if (!std::isfinite(spacing) || spacing <= 0 ||
          !std::isfinite(high_weight) || high_weight < 0 || high_weight > 1)
        return false;
      out.indices = {static_cast<std::int32_t>(last), 0};
      out.weights = {1 - high_weight, high_weight};
      out.count = 2;
      return true;
    }
  }
  if (centres.size == 1U || coordinate <= centres.data[0U]) {
    out.indices[0U] = 0;
    out.weights[0U] = 1.0;
    out.count = 1U;
    out.clamped = coordinate != centres.data[0U];
    return std::isfinite(centres.data[0U]);
  }
  const std::size_t last = centres.size - 1U;
  if (coordinate >= centres.data[last]) {
    out.indices[0U] = static_cast<std::int32_t>(last);
    out.weights[0U] = 1.0;
    out.count = 1U;
    out.clamped = coordinate != centres.data[last];
    return std::isfinite(centres.data[last]);
  }
  const double* right_ptr =
      std::lower_bound(centres.data, centres.data + centres.size, coordinate);
  const std::size_t right = static_cast<std::size_t>(right_ptr - centres.data);
  if (*right_ptr == coordinate) {
    out.indices[0U] = static_cast<std::int32_t>(right);
    out.weights[0U] = 1.0;
    out.count = 1U;
    return true;
  }
  const std::size_t left = right - 1U;
  const double spacing = centres.data[right] - centres.data[left];
  const double upper_weight = (coordinate - centres.data[left]) / spacing;
  const double lower_weight = 1.0 - upper_weight;
  if (!(spacing > 0.0) || !std::isfinite(spacing) ||
      !std::isfinite(lower_weight) || !std::isfinite(upper_weight) ||
      lower_weight < 0.0 || upper_weight < 0.0) {
    return false;
  }
  out.indices = {static_cast<std::int32_t>(left),
                 static_cast<std::int32_t>(right)};
  out.weights = {lower_weight, upper_weight};
  out.count = 2U;
  return true;
}

ParcelTrajectoryCandidate trajectory_failure(
    const ParcelTrajectoryInput& input, ParcelTrajectoryStatus status) noexcept {
  ParcelTrajectoryCandidate result;
  result.kind = input.kind;
  result.requested_duration_s = input.duration_s;
  result.status = status;
  return result;
}

StaticIbmReboundCandidate rebound_failure(
    TrajectoryCandidateKind kind, StaticIbmReboundStatus status) noexcept {
  StaticIbmReboundCandidate result;
  result.kind = kind;
  result.status = status;
  return result;
}

PhysicalDomainOutletCandidate outlet_failure(
    PhysicalDomainOutletStatus status) noexcept {
  PhysicalDomainOutletCandidate result;
  result.status = status;
  return result;
}

bool valid_segment(const ParcelTrajectorySegment& segment) noexcept {
  return std::isfinite(segment.begin_time_s) &&
         std::isfinite(segment.end_time_s) &&
         segment.end_time_s > segment.begin_time_s &&
         finite(segment.begin_position_m) && finite(segment.end_position_m) &&
         finite(segment.begin_velocity_m_per_s) &&
         finite(segment.end_velocity_m_per_s);
}

bool same(Vector3 left, Vector3 right) noexcept {
  return left[0U] == right[0U] && left[1U] == right[1U] &&
         left[2U] == right[2U];
}

bool valid_trajectory_candidate(
    const ParcelTrajectoryCandidate& trajectory) noexcept {
  if (!trajectory.succeeded() || !valid_kind(trajectory.kind) ||
      validate_parcel_state(trajectory.parcel) != ParcelStateStatus::success ||
      !std::isfinite(trajectory.requested_duration_s) ||
      !std::isfinite(trajectory.advanced_duration_s) ||
      trajectory.requested_duration_s < 0.0 ||
      trajectory.advanced_duration_s < 0.0 ||
      trajectory.advanced_duration_s > trajectory.requested_duration_s ||
      trajectory.segment_count > trajectory.segments.size()) {
    return false;
  }
  if (trajectory.segment_count == 0U) {
    return trajectory.advanced_duration_s == 0.0;
  }
  for (std::uint32_t ordinal = 0U; ordinal < trajectory.segment_count;
       ++ordinal) {
    const ParcelTrajectorySegment& segment = trajectory.segments[ordinal];
    if (!valid_segment(segment) || segment.ordinal != ordinal ||
        (ordinal == 0U && segment.begin_time_s != 0.0)) {
      return false;
    }
    if (ordinal != 0U) {
      const ParcelTrajectorySegment& previous =
          trajectory.segments[ordinal - 1U];
      if (segment.begin_time_s != previous.end_time_s ||
          !same(segment.begin_position_m, previous.end_position_m) ||
          !same(segment.begin_velocity_m_per_s,
                previous.end_velocity_m_per_s)) {
        return false;
      }
    }
  }
  const ParcelTrajectorySegment& last =
      trajectory.segments[trajectory.segment_count - 1U];
  return last.end_time_s == trajectory.advanced_duration_s &&
         same(last.end_position_m, trajectory.parcel.position_m) &&
         same(last.end_velocity_m_per_s,
              trajectory.parcel.velocity_m_per_s);
}

bool append_segment(StaticIbmReboundCandidate& candidate,
                    double begin_time, double end_time, Vector3 begin_position,
                    Vector3 end_position, Vector3 begin_velocity,
                    Vector3 end_velocity) noexcept {
  if (!(end_time > begin_time)) {
    return true;
  }
  if (candidate.segment_count >= candidate.segments.size()) {
    return false;
  }
  const std::uint32_t ordinal = candidate.segment_count;
  candidate.segments[ordinal] = {
      ordinal, begin_time, end_time, begin_position, end_position,
      begin_velocity, end_velocity};
  ++candidate.segment_count;
  return true;
}

Vector3 reflect(Vector3 value, Vector3 fluid_normal,
                double normal_restitution,
                double tangential_restitution) noexcept {
  const double normal_component = dot(value, fluid_normal);
  const Vector3 normal_part = multiply(normal_component, fluid_normal);
  const Vector3 tangential_part = subtract(value, normal_part);
  return add(multiply(tangential_restitution, tangential_part),
             multiply(-normal_restitution, normal_part));
}

Vector3 apply_rebound_history(
    Vector3 value, const StaticIbmReboundCandidate& candidate,
    const StaticIbmReboundConfig& config) noexcept {
  for (std::uint32_t event = 0U; event < candidate.event_count; ++event) {
    value = reflect(value, candidate.events[event].fluid_side_normal,
                    config.normal_restitution,
                    config.tangential_restitution);
  }
  return value;
}

bool inside_domain(Vector3 point, Real3 lower, Real3 upper) noexcept {
  return finite(point) && point[0U] >= lower.x && point[0U] <= upper.x &&
         point[1U] >= lower.y && point[1U] <= upper.y &&
         point[2U] >= lower.z && point[2U] <= upper.z;
}

}  // namespace

ParcelGridCouplingStencil
build_parcel_grid_coupling_stencil(const CartesianGeometryPlan &geometry,
                                   const MeshPatch &patch, Vector3 position_m,
                                   std::array<bool, 3U> periodic) noexcept {
  ParcelGridCouplingStencil result;
  const Int3 cells = geometry.global_cells();
  const Real3 lower = geometry.lower();
  const Real3 upper = geometry.upper();
  if (geometry.fingerprint() == 0U || !valid_patch(cells, patch) ||
      !finite(position_m)) {
    return stencil_failure(ParcelGridStencilStatus::invalid_input);
  }
  for (std::size_t axis = 0; axis < 3; ++axis) {
    const double low = axis == 0 ? lower.x : axis == 1 ? lower.y : lower.z;
    const double high = axis == 0 ? upper.x : axis == 1 ? upper.y : upper.z;
    if (periodic[axis]) {
      const double length = high - low;
      if (!std::isfinite(length) || length <= 0)
        return stencil_failure(ParcelGridStencilStatus::invalid_input);
      if (position_m[axis] < low || position_m[axis] >= high) {
        double offset = std::fmod(position_m[axis] - low, length);
        if (!std::isfinite(offset))
          return stencil_failure(ParcelGridStencilStatus::invalid_input);
        if (offset < 0)
          offset += length;
        position_m[axis] = low + offset;
        if (position_m[axis] >= high)
          position_m[axis] = low;
      }
    } else if (position_m[axis] < low || position_m[axis] > high) {
      return stencil_failure(ParcelGridStencilStatus::outside_domain);
    }
  }
  const std::array<const AxisMetrics*, 3U> metrics{
      &geometry.x(), &geometry.y(), &geometry.z()};
  std::array<AxisStencil, 3U> axes{};
  for (std::size_t axis = 0U; axis < axes.size(); ++axis) {
    if (metrics[axis]->centres().size !=
            static_cast<std::size_t>(component(cells, axis)) ||
        !build_axis_stencil(*metrics[axis], position_m[axis], axes[axis],
                            periodic[axis])) {
      return stencil_failure(ParcelGridStencilStatus::invalid_input);
    }
    result.boundary_clamped |= axes[axis].clamped;
  }

  for (std::uint8_t z = 0U; z < axes[2U].count; ++z) {
    for (std::uint8_t y = 0U; y < axes[1U].count; ++y) {
      for (std::uint8_t x = 0U; x < axes[0U].count; ++x) {
        const Int3 index{axes[0U].indices[x], axes[1U].indices[y],
                         axes[2U].indices[z]};
        std::uint64_t id = 0U;
        if (!global_cell_id(index, cells, id) ||
            result.entry_count >= result.entries.size()) {
          return stencil_failure(ParcelGridStencilStatus::invalid_input);
        }
        const std::int32_t owner_x = owner_coordinate(
            index.x, cells.x, patch.process_grid.x);
        const std::int32_t owner_y = owner_coordinate(
            index.y, cells.y, patch.process_grid.y);
        const std::int32_t owner_z = owner_coordinate(
            index.z, cells.z, patch.process_grid.z);
        const std::int64_t owner =
            owner_x + static_cast<std::int64_t>(patch.process_grid.x) *
                          (owner_y +
                           static_cast<std::int64_t>(patch.process_grid.y) *
                               owner_z);
        const double weight = axes[0U].weights[x] * axes[1U].weights[y] *
                              axes[2U].weights[z];
        if (owner < 0 || owner > std::numeric_limits<std::int32_t>::max() ||
            !std::isfinite(weight) || weight < 0.0) {
          return stencil_failure(ParcelGridStencilStatus::non_finite_output);
        }
        result.entries[result.entry_count++] = {
            id, index, static_cast<std::int32_t>(owner), weight};
      }
    }
  }
  std::sort(result.entries.begin(),
            result.entries.begin() + result.entry_count,
            [](const ParcelGridStencilEntry& left,
               const ParcelGridStencilEntry& right) noexcept {
              return left.global_cell < right.global_cell;
            });
  double sum = 0.0;
  for (std::size_t entry = 0U; entry < result.entry_count; ++entry) {
    sum += result.entries[entry].weight;
  }
  if (!(sum > 0.0) || !std::isfinite(sum)) {
    return stencil_failure(ParcelGridStencilStatus::non_finite_output);
  }
  std::size_t largest = 0U;
  double normalized_sum = 0.0;
  for (std::size_t entry = 0U; entry < result.entry_count; ++entry) {
    result.entries[entry].weight /= sum;
    normalized_sum += result.entries[entry].weight;
    if (result.entries[entry].weight > result.entries[largest].weight) {
      largest = entry;
    }
  }
  result.entries[largest].weight += 1.0 - normalized_sum;
  if (!std::isfinite(result.entries[largest].weight) ||
      result.entries[largest].weight < 0.0) {
    return stencil_failure(ParcelGridStencilStatus::non_finite_output);
  }
  result.available = true;
  result.status = ParcelGridStencilStatus::success;
  return result;
}

ParcelTrajectoryCandidate make_parcel_trajectory_candidate(
    const ParcelTrajectoryInput& input) noexcept {
  if (validate_parcel_state(input.parcel) != ParcelStateStatus::success ||
      !valid_kind(input.kind) || !std::isfinite(input.duration_s) ||
      input.duration_s < 0.0 ||
      !std::isfinite(input.event_stop_time_s) ||
      (input.event_stop_time_s >= 0.0 &&
       input.event_stop_time_s > input.duration_s)) {
    return trajectory_failure(input, ParcelTrajectoryStatus::invalid_input);
  }
  const double duration = input.event_stop_time_s >= 0.0
                              ? input.event_stop_time_s
                              : input.duration_s;
  if (!std::isfinite(input.parcel.age_s + duration)) {
    return trajectory_failure(input,
                              ParcelTrajectoryStatus::non_finite_output);
  }

  ParcelTrajectoryCandidate result;
  result.kind = input.kind;
  result.requested_duration_s = input.duration_s;
  result.parcel = input.parcel;
  if (duration == 0.0) {
    result.available = true;
    result.status = ParcelTrajectoryStatus::success;
    return result;
  }
  if (!std::isfinite(input.characteristic_length_m) ||
      !(input.characteristic_length_m > 0.0) ||
      !std::isfinite(input.maximum_particle_cfl) ||
      !(input.maximum_particle_cfl > 0.0) ||
      input.maximum_substeps == 0U ||
      input.maximum_substeps > kMaximumTrajectorySegments ||
      input.acceleration == nullptr) {
    return trajectory_failure(input, ParcelTrajectoryStatus::invalid_input);
  }

  const double distance_limit =
      input.characteristic_length_m * input.maximum_particle_cfl;
  if (!(distance_limit > 0.0) || !std::isfinite(distance_limit)) {
    return trajectory_failure(input, ParcelTrajectoryStatus::invalid_input);
  }
  Vector3 position = input.parcel.position_m;
  Vector3 velocity = input.parcel.velocity_m_per_s;
  double elapsed = 0.0;
  while (elapsed < duration) {
    if (result.segment_count >= input.maximum_substeps) {
      return trajectory_failure(input, ParcelTrajectoryStatus::substep_limit);
    }
    const ParcelKinematicSample begin_sample{position, velocity, elapsed};
    const ParcelAccelerationSample begin_acceleration =
        input.acceleration(begin_sample,
                           input.immutable_acceleration_context);
    if (!begin_acceleration.available) {
      return trajectory_failure(
          input, ParcelTrajectoryStatus::acceleration_unavailable);
    }
    if (!finite(begin_acceleration.acceleration_m_per_s2)) {
      return trajectory_failure(input,
                                ParcelTrajectoryStatus::non_finite_output);
    }
    const double speed = norm(velocity);
    const double acceleration_norm =
        norm(begin_acceleration.acceleration_m_per_s2);
    if (!std::isfinite(speed) || !std::isfinite(acceleration_norm)) {
      return trajectory_failure(input,
                                ParcelTrajectoryStatus::non_finite_output);
    }
    double step = duration - elapsed;
    if (speed > 0.0 || acceleration_norm > 0.0) {
      const double root = std::sqrt(speed * speed +
                                    2.0 * acceleration_norm * distance_limit);
      const double bounded =
          acceleration_norm > 0.0
              ? 2.0 * distance_limit / (speed + root)
              : distance_limit / speed;
      step = std::min(step, bounded);
    }
    if (!(step > 0.0) || !std::isfinite(step)) {
      return trajectory_failure(input, ParcelTrajectoryStatus::substep_limit);
    }

    Vector3 end_position{};
    Vector3 end_velocity{};
    bool accepted_step = false;
    for (std::uint32_t retry = 0U; retry < 64U; ++retry) {
      const Vector3 midpoint_position =
          add(add(position, multiply(0.5 * step, velocity)),
              multiply(0.125 * step * step,
                       begin_acceleration.acceleration_m_per_s2));
      const Vector3 midpoint_velocity =
          add(velocity, multiply(0.5 * step,
                                 begin_acceleration.acceleration_m_per_s2));
      const ParcelKinematicSample midpoint_sample{
          midpoint_position, midpoint_velocity, elapsed + 0.5 * step};
      const ParcelAccelerationSample midpoint_acceleration =
          input.acceleration(midpoint_sample,
                             input.immutable_acceleration_context);
      if (!midpoint_acceleration.available) {
        return trajectory_failure(
            input, ParcelTrajectoryStatus::acceleration_unavailable);
      }
      if (!finite(midpoint_acceleration.acceleration_m_per_s2)) {
        return trajectory_failure(input,
                                  ParcelTrajectoryStatus::non_finite_output);
      }
      end_position =
          add(add(position, multiply(step, velocity)),
              multiply(0.5 * step * step,
                       midpoint_acceleration.acceleration_m_per_s2));
      end_velocity =
          add(velocity, multiply(step,
                                 midpoint_acceleration.acceleration_m_per_s2));
      const double travelled = norm(subtract(end_position, position));
      if (!finite(end_position) || !finite(end_velocity) ||
          !std::isfinite(travelled)) {
        return trajectory_failure(input,
                                  ParcelTrajectoryStatus::non_finite_output);
      }
      if (travelled <= distance_limit * (1.0 + kRoundoff)) {
        accepted_step = true;
        break;
      }
      step *= 0.5;
      if (!(step > 0.0) || elapsed + step == elapsed) {
        break;
      }
    }
    if (!accepted_step) {
      return trajectory_failure(input, ParcelTrajectoryStatus::substep_limit);
    }
    double end_time = elapsed + step;
    if (duration - end_time <=
        kRoundoff * std::max(1.0, std::abs(duration))) {
      end_time = duration;
      step = end_time - elapsed;
    }
    const std::uint32_t ordinal = result.segment_count;
    result.segments[ordinal] = {ordinal,
                                elapsed,
                                end_time,
                                position,
                                end_position,
                                velocity,
                                end_velocity};
    ++result.segment_count;
    elapsed = end_time;
    position = end_position;
    velocity = end_velocity;
  }
  result.parcel.position_m = position;
  result.parcel.velocity_m_per_s = velocity;
  result.parcel.age_s += duration;
  if (validate_parcel_state(result.parcel) != ParcelStateStatus::success) {
    return trajectory_failure(input,
                              ParcelTrajectoryStatus::non_finite_output);
  }
  result.advanced_duration_s = duration;
  result.available = true;
  result.status = ParcelTrajectoryStatus::success;
  return result;
}

StaticIbmReboundCandidate resolve_static_ibm_rebounds(
    const ParcelTrajectoryCandidate& trajectory,
    const ImmersedSurfacePlan& surface,
    StaticIbmReboundConfig config) noexcept {
  StaticIbmReboundCandidate result;
  result.kind = trajectory.kind;
  if (!valid_trajectory_candidate(trajectory) ||
      surface.fingerprint() == 0U || surface.triangles().data == nullptr ||
      surface.triangles().size == 0U ||
      (config.fluid_side != ImmersedFluidSide::outside &&
       config.fluid_side != ImmersedFluidSide::inside) ||
      !std::isfinite(config.normal_restitution) ||
      config.normal_restitution < 0.0 || config.normal_restitution > 1.0 ||
      !std::isfinite(config.tangential_restitution) ||
      config.tangential_restitution < 0.0 ||
      config.tangential_restitution > 1.0 ||
      config.maximum_events > kMaximumReboundEvents) {
    return rebound_failure(trajectory.kind,
                           StaticIbmReboundStatus::invalid_input);
  }
  result.parcel = trajectory.parcel;
  Vector3 final_position = trajectory.parcel.position_m;
  Vector3 final_velocity = trajectory.parcel.velocity_m_per_s;
  for (std::uint32_t source = 0U; source < trajectory.segment_count;
       ++source) {
    const ParcelTrajectorySegment& segment = trajectory.segments[source];
    if (!valid_segment(segment) || segment.ordinal != source) {
      return rebound_failure(trajectory.kind,
                             StaticIbmReboundStatus::invalid_input);
    }
    Vector3 current_position =
        source == 0U ? segment.begin_position_m : final_position;
    Vector3 current_velocity =
        apply_rebound_history(segment.begin_velocity_m_per_s, result, config);
    const Vector3 source_displacement =
        subtract(segment.end_position_m, segment.begin_position_m);
    Vector3 target_position =
        add(current_position,
            apply_rebound_history(source_displacement, result, config));
    Vector3 target_velocity =
        apply_rebound_history(segment.end_velocity_m_per_s, result, config);
    double current_time = segment.begin_time_s;
    const double target_time = segment.end_time_s;
    Vector3 previous_normal{};
    bool just_rebounded = false;
    while (current_time < target_time) {
      const Vector3 displacement = subtract(target_position, current_position);
      if (!(norm(displacement) > 0.0)) {
        if (!append_segment(result, current_time, target_time,
                            current_position, target_position,
                            current_velocity, target_velocity)) {
          return rebound_failure(trajectory.kind,
                                 StaticIbmReboundStatus::segment_limit);
        }
        current_position = target_position;
        current_velocity = target_velocity;
        current_time = target_time;
        break;
      }

      Vector3 query_begin = current_position;
      if (just_rebounded) {
        const double scene_scale = std::max(
            {1.0, std::abs(current_position[0U]),
             std::abs(current_position[1U]),
             std::abs(current_position[2U]), norm(displacement)});
        query_begin = add(current_position,
                          multiply(4096.0 *
                                       std::numeric_limits<double>::epsilon() *
                                       scene_scale,
                                   previous_normal));
      }
      SurfaceSegmentIntersection hit{};
      const Status queried = surface.first_segment_intersection(
          {query_begin[0U], query_begin[1U], query_begin[2U]},
          {target_position[0U], target_position[1U], target_position[2U]}, hit);
      if (!queried) {
        return rebound_failure(
            trajectory.kind, StaticIbmReboundStatus::surface_query_failure);
      }
      if (hit.triangle == kInvalidSurfaceTriangle) {
        if (!append_segment(result, current_time, target_time,
                            current_position, target_position,
                            current_velocity, target_velocity)) {
          return rebound_failure(trajectory.kind,
                                 StaticIbmReboundStatus::segment_limit);
        }
        current_position = target_position;
        current_velocity = target_velocity;
        current_time = target_time;
        break;
      }
      if (hit.triangle >= surface.triangles().size ||
          !std::isfinite(hit.segment_fraction) ||
          hit.segment_fraction < 0.0 || hit.segment_fraction > 1.0) {
        return rebound_failure(
            trajectory.kind, StaticIbmReboundStatus::surface_query_failure);
      }
      const SurfaceTriangle& triangle =
          surface.triangles().data[static_cast<std::size_t>(hit.triangle)];
      Vector3 normal = to_vector(triangle.geometric_outward_normal);
      if (config.fluid_side == ImmersedFluidSide::inside) {
        normal = multiply(-1.0, normal);
      }
      const Vector3 hit_position = to_vector(hit.point);
      const double query_distance = norm(subtract(target_position, query_begin));
      const double full_distance = norm(displacement);
      double path_fraction = hit.segment_fraction;
      if (query_distance > 0.0 && full_distance > 0.0) {
        path_fraction = norm(subtract(hit_position, current_position)) /
                        full_distance;
      }
      path_fraction = std::clamp(path_fraction, 0.0, 1.0);
      const double event_time =
          current_time + path_fraction * (target_time - current_time);
      const Vector3 incoming_velocity =
          interpolate(current_velocity, target_velocity, path_fraction);
      if (dot(incoming_velocity, normal) >= 0.0) {
        if (path_fraction > kRoundoff) {
          return rebound_failure(trajectory.kind,
                                 StaticIbmReboundStatus::invalid_input);
        }
        // A parcel already on the wall and moving toward the fluid must not
        // rebound at fraction zero.  Nudge only the next immutable query;
        // the reported path remains anchored at the exact wall point.
        previous_normal = normal;
        just_rebounded = true;
        continue;
      }
      if (result.event_count >= config.maximum_events) {
        return rebound_failure(trajectory.kind,
                               StaticIbmReboundStatus::event_limit);
      }
      if (!append_segment(result, current_time, event_time,
                          current_position, hit_position, current_velocity,
                          incoming_velocity)) {
        return rebound_failure(trajectory.kind,
                               StaticIbmReboundStatus::segment_limit);
      }
      const Vector3 outgoing_velocity =
          reflect(incoming_velocity, normal, config.normal_restitution,
                  config.tangential_restitution);
      const std::uint32_t event_ordinal = result.event_count;
      result.events[event_ordinal] = {
          event_ordinal, source, hit.triangle, event_time, hit_position,
          normal, incoming_velocity, outgoing_velocity};
      ++result.event_count;

      const Vector3 remaining_displacement =
          subtract(target_position, hit_position);
      target_position =
          add(hit_position,
              reflect(remaining_displacement, normal,
                      config.normal_restitution,
                      config.tangential_restitution));
      target_velocity =
          reflect(target_velocity, normal, config.normal_restitution,
                  config.tangential_restitution);
      current_position = hit_position;
      current_velocity = outgoing_velocity;
      current_time = event_time;
      previous_normal = normal;
      just_rebounded = true;
      if (target_time - current_time <=
          kRoundoff * std::max(1.0, std::abs(target_time))) {
        current_time = target_time;
      }
    }
    final_position = current_position;
    final_velocity = current_velocity;
  }
  result.parcel.position_m = final_position;
  result.parcel.velocity_m_per_s = final_velocity;
  if (validate_parcel_state(result.parcel) != ParcelStateStatus::success) {
    return rebound_failure(trajectory.kind,
                           StaticIbmReboundStatus::non_finite_output);
  }
  result.available = true;
  result.status = StaticIbmReboundStatus::success;
  return result;
}

PhysicalDomainOutletCandidate resolve_physical_domain_outlet(
    const ParcelTrajectoryCandidate& trajectory,
    const CartesianGeometryPlan& geometry) noexcept {
  PhysicalDomainOutletCandidate result;
  if (!valid_trajectory_candidate(trajectory) ||
      geometry.fingerprint() == 0U) {
    return outlet_failure(PhysicalDomainOutletStatus::invalid_input);
  }
  const Real3 lower = geometry.lower();
  const Real3 upper = geometry.upper();
  if (!(lower.x < upper.x) || !(lower.y < upper.y) ||
      !(lower.z < upper.z)) {
    return outlet_failure(PhysicalDomainOutletStatus::invalid_input);
  }
  result.parcel = trajectory.parcel;
  if (trajectory.segment_count == 0U) {
    if (!inside_domain(trajectory.parcel.position_m, lower, upper)) {
      return outlet_failure(PhysicalDomainOutletStatus::invalid_input);
    }
    result.available = true;
    result.status = PhysicalDomainOutletStatus::success;
    result.parcel_retained = true;
    return result;
  }
  if (!inside_domain(trajectory.segments[0U].begin_position_m, lower, upper)) {
    return outlet_failure(PhysicalDomainOutletStatus::invalid_input);
  }

  const std::array<double, 3U> minimum{lower.x, lower.y, lower.z};
  const std::array<double, 3U> maximum{upper.x, upper.y, upper.z};
  for (std::uint32_t ordinal = 0U; ordinal < trajectory.segment_count;
       ++ordinal) {
    const ParcelTrajectorySegment& segment = trajectory.segments[ordinal];
    if (!valid_segment(segment) || segment.ordinal != ordinal ||
        !inside_domain(segment.begin_position_m, lower, upper)) {
      return outlet_failure(PhysicalDomainOutletStatus::invalid_input);
    }
    double event_fraction = std::numeric_limits<double>::infinity();
    PhysicalDomainFace event_face = PhysicalDomainFace::x_min;
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      const double begin = segment.begin_position_m[axis];
      const double end = segment.end_position_m[axis];
      double candidate = std::numeric_limits<double>::infinity();
      PhysicalDomainFace face = static_cast<PhysicalDomainFace>(2U * axis);
      if (end < minimum[axis] && end < begin) {
        candidate = (minimum[axis] - begin) / (end - begin);
      } else if (end > maximum[axis] && end > begin) {
        candidate = (maximum[axis] - begin) / (end - begin);
        face = static_cast<PhysicalDomainFace>(2U * axis + 1U);
      }
      if (candidate >= 0.0 && candidate <= 1.0 &&
          (candidate < event_fraction ||
           (candidate == event_fraction &&
            static_cast<std::uint8_t>(face) <
                static_cast<std::uint8_t>(event_face)))) {
        event_fraction = candidate;
        event_face = face;
      }
    }
    if (!std::isfinite(event_fraction)) {
      continue;
    }
    Vector3 event_position =
        interpolate(segment.begin_position_m, segment.end_position_m,
                    event_fraction);
    const std::size_t face_index =
        static_cast<std::size_t>(event_face);
    const std::size_t event_axis = face_index / 2U;
    event_position[event_axis] =
        face_index % 2U == 0U ? minimum[event_axis] : maximum[event_axis];
    const Vector3 event_velocity =
        interpolate(segment.begin_velocity_m_per_s,
                    segment.end_velocity_m_per_s, event_fraction);
    const double event_time =
        segment.begin_time_s +
        event_fraction * (segment.end_time_s - segment.begin_time_s);
    const double represented_mass =
        trajectory.parcel.droplet_mass_kg * trajectory.parcel.multiplicity;
    const Vector3 represented_momentum =
        multiply(represented_mass, event_velocity);
    if (!finite(event_position) || !finite(event_velocity) ||
        !std::isfinite(event_time) || !(represented_mass > 0.0) ||
        !std::isfinite(represented_mass) || !finite(represented_momentum)) {
      return outlet_failure(PhysicalDomainOutletStatus::non_finite_output);
    }
    result.parcel.position_m = event_position;
    result.parcel.velocity_m_per_s = event_velocity;
    const double remaining_time =
        trajectory.advanced_duration_s - event_time;
    result.parcel.age_s -= remaining_time;
    if (validate_parcel_state(result.parcel) != ParcelStateStatus::success) {
      return outlet_failure(PhysicalDomainOutletStatus::non_finite_output);
    }
    result.ledger.recorded = true;
    result.ledger.requires_parcel_removal = true;
    result.ledger.permits_gas_source_deposition = false;
    result.ledger.parcel_id = trajectory.parcel.id;
    result.ledger.face = event_face;
    result.ledger.event_time_s = event_time;
    result.ledger.position_m = event_position;
    result.ledger.velocity_m_per_s = event_velocity;
    result.ledger.represented_liquid_mass_kg = represented_mass;
    result.ledger.represented_momentum_kg_m_per_s = represented_momentum;
    result.available = true;
    result.status = PhysicalDomainOutletStatus::success;
    result.parcel_retained = false;
    return result;
  }
  if (!inside_domain(trajectory.parcel.position_m, lower, upper)) {
    return outlet_failure(PhysicalDomainOutletStatus::invalid_input);
  }
  result.available = true;
  result.status = PhysicalDomainOutletStatus::success;
  result.parcel_retained = true;
  return result;
}

}  // namespace hundun::v04::spray::detail

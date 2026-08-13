// SPDX-License-Identifier: Apache-2.0

#include "mesh_focus_detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <new>
#include <queue>
#include <utility>
#include <vector>

namespace hundun::v04::detail {
namespace {

constexpr std::uint32_t kFocusInput = 101U;
constexpr std::uint32_t kFocusCellCount = 102U;
constexpr std::uint32_t kFocusSpacing = 103U;
constexpr std::uint32_t kFocusGrowth = 104U;
constexpr std::uint32_t kFocusMemory = 105U;
constexpr std::uint64_t kMaximumMetricSamples = 4U * 1024U * 1024U;
constexpr std::uint64_t kMaximumFaceProjectionVisits = 100000000U;
constexpr double kConstraintTolerance = 1.0e-12;

bool checked_add(std::uint64_t left, std::uint64_t right,
                 std::uint64_t& out) noexcept {
  if (right > std::numeric_limits<std::uint64_t>::max() - left) {
    return false;
  }
  out = left + right;
  return true;
}

bool checked_multiply(std::uint64_t left, std::uint64_t right,
                      std::uint64_t& out) noexcept {
  if (left != 0U && right > std::numeric_limits<std::uint64_t>::max() / left) {
    return false;
  }
  out = left * right;
  return true;
}

bool within_peak(std::uint64_t persistent_bytes, std::uint64_t scratch_bytes,
                 std::uint64_t limit) noexcept {
  std::uint64_t peak = 0U;
  return limit != 0U && checked_add(persistent_bytes, scratch_bytes, peak) &&
         peak <= limit;
}

double component(Real3 value, std::size_t axis) noexcept {
  return axis == 0U ? value.x : (axis == 1U ? value.y : value.z);
}

std::int32_t component(Int3 value, std::size_t axis) noexcept {
  return axis == 0U ? value.x : (axis == 1U ? value.y : value.z);
}

struct IntervalTarget {
  double lower{};
  double upper{};
  double spacing{};
};

struct FocusStartLess {
  bool operator()(const IntervalTarget& left,
                  const IntervalTarget& right) const noexcept {
    if (left.lower != right.lower) {
      return left.lower < right.lower;
    }
    if (left.upper != right.upper) {
      return left.upper < right.upper;
    }
    return left.spacing < right.spacing;
  }
};

struct ActiveFocus {
  double spacing{};
  double upper{};
  std::size_t id{};
};

struct ActiveFocusGreater {
  bool operator()(const ActiveFocus& left,
                  const ActiveFocus& right) const noexcept {
    if (left.spacing != right.spacing) {
      return left.spacing > right.spacing;
    }
    if (left.upper != right.upper) {
      return left.upper > right.upper;
    }
    return left.id > right.id;
  }
};

bool fill_spacing_envelope(const std::vector<double>& sample_centres,
                           double base_spacing, double minimum_spacing,
                           double growth_ratio,
                           std::vector<IntervalTarget> focuses,
                           std::vector<double>& spacing) {
  std::sort(focuses.begin(), focuses.end(), FocusStartLess{});
  spacing.assign(sample_centres.size(), base_spacing);
  std::priority_queue<ActiveFocus, std::vector<ActiveFocus>,
                      ActiveFocusGreater>
      active;
  std::size_t next = 0U;
  for (std::size_t sample = 0U; sample < sample_centres.size(); ++sample) {
    const double coordinate = sample_centres[sample];
    while (next < focuses.size() && focuses[next].lower <= coordinate) {
      active.push({focuses[next].spacing, focuses[next].upper, next});
      ++next;
    }
    while (!active.empty() && active.top().upper < coordinate) {
      active.pop();
    }
    if (!active.empty()) {
      spacing[sample] = std::min(base_spacing, active.top().spacing);
    }
  }

  if (growth_ratio > 1.0 && sample_centres.size() > 1U) {
    // Apply an O(samples) conservative physical-space envelope.  A spacing h
    // may grow by at most (r-1)*distance when measured in local h-sized cells;
    // the later exact cell pass still checks the geometric adjacent ratio.
    const double slope = growth_ratio - 1.0;
    for (std::size_t index = 1U; index < spacing.size(); ++index) {
      const double distance = sample_centres[index] - sample_centres[index - 1U];
      spacing[index] =
          std::min(spacing[index], spacing[index - 1U] + slope * distance);
    }
    for (std::size_t index = spacing.size() - 1U; index > 0U; --index) {
      const double distance = sample_centres[index] - sample_centres[index - 1U];
      spacing[index - 1U] =
          std::min(spacing[index - 1U], spacing[index] + slope * distance);
    }
  }
  for (double& value : spacing) {
    value = std::max(minimum_spacing, std::min(base_spacing, value));
  }
  return true;
}

bool make_metric_samples(double lower, double upper, double base_spacing,
                         double minimum_spacing, double growth_ratio,
                         const std::vector<IntervalTarget>& focuses,
                         std::vector<double>& coordinates,
                         std::vector<double>& cumulative) {
  std::vector<double> knots;
  knots.reserve(2U + focuses.size() * 2U);
  knots.push_back(lower);
  knots.push_back(upper);
  for (const IntervalTarget& focus : focuses) {
    knots.push_back(focus.lower);
    knots.push_back(focus.upper);
  }
  std::sort(knots.begin(), knots.end());
  knots.erase(std::unique(knots.begin(), knots.end()), knots.end());

  coordinates.clear();
  cumulative.clear();
  std::uint64_t total_pieces = 0U;
  for (std::size_t interval = 0U; interval + 1U < knots.size(); ++interval) {
    const double a = knots[interval];
    const double b = knots[interval + 1U];
    if (!(a < b)) {
      continue;
    }
    const double scale =
        std::max(minimum_spacing, std::min(base_spacing, b - a));
    std::uint64_t pieces = static_cast<std::uint64_t>(
        std::ceil((b - a) /
                  std::max(scale * 0.25, minimum_spacing * 0.25)));
    pieces = std::max<std::uint64_t>(pieces, 8U);
    if (pieces > kMaximumMetricSamples - total_pieces) {
      return false;
    }
    total_pieces += pieces;
  }
  coordinates.reserve(static_cast<std::size_t>(total_pieces) + 1U);
  cumulative.reserve(static_cast<std::size_t>(total_pieces) + 1U);
  coordinates.push_back(lower);
  for (std::size_t interval = 0U; interval + 1U < knots.size(); ++interval) {
    const double a = knots[interval];
    const double b = knots[interval + 1U];
    if (!(a < b)) {
      continue;
    }
    const double scale = std::max(minimum_spacing,
                                  std::min(base_spacing, b - a));
    std::uint64_t pieces = static_cast<std::uint64_t>(
        std::ceil((b - a) / std::max(scale * 0.25, minimum_spacing * 0.25)));
    pieces = std::max<std::uint64_t>(pieces, 8U);
    for (std::uint64_t item = 1U; item <= pieces; ++item) {
      const double x = item == pieces
                           ? b
                           : a + (b - a) * static_cast<double>(item) /
                                     static_cast<double>(pieces);
      coordinates.push_back(x);
    }
  }
  if (coordinates.size() < 2U || coordinates.back() != upper) {
    return false;
  }
  std::vector<double> sample_centres(total_pieces);
  for (std::size_t index = 0U; index < sample_centres.size(); ++index) {
    sample_centres[index] =
        coordinates[index] + 0.5 * (coordinates[index + 1U] - coordinates[index]);
  }
  std::vector<double> spacing;
  if (!fill_spacing_envelope(sample_centres, base_spacing, minimum_spacing,
                             growth_ratio, focuses, spacing)) {
    return false;
  }
  cumulative.resize(coordinates.size());
  cumulative[0] = 0.0;
  for (std::size_t index = 0U; index < spacing.size(); ++index) {
    const double integral = cumulative[index] +
                            (coordinates[index + 1U] - coordinates[index]) /
                                spacing[index];
    if (!std::isfinite(integral) || !(integral > cumulative[index])) {
      return false;
    }
    cumulative[index + 1U] = integral;
  }
  return true;
}

double inverse_cumulative(const std::vector<double>& coordinates,
                          const std::vector<double>& cumulative,
                          double target) noexcept {
  const auto found = std::lower_bound(cumulative.begin(), cumulative.end(),
                                      target);
  if (found == cumulative.begin()) {
    return coordinates.front();
  }
  if (found == cumulative.end()) {
    return coordinates.back();
  }
  const std::size_t right = static_cast<std::size_t>(found - cumulative.begin());
  const std::size_t left = right - 1U;
  const double weight = (target - cumulative[left]) /
                        (cumulative[right] - cumulative[left]);
  return coordinates[left] + weight * (coordinates[right] - coordinates[left]);
}

enum class FaceBuildResult : std::uint8_t {
  feasible,
  too_coarse,
  too_fine
};

FaceBuildResult build_faces_for_count(
    double lower, double upper, double base_spacing, double minimum_spacing,
    double growth_ratio, const std::vector<IntervalTarget>& focuses,
    const std::vector<double>& coordinates,
    const std::vector<double>& cumulative, std::uint64_t cell_count,
    std::vector<double>& faces) {
  faces.resize(static_cast<std::size_t>(cell_count) + 1U);
  faces.front() = lower;
  faces.back() = upper;
  const double integral_step = cumulative.back() / static_cast<double>(cell_count);
  for (std::uint64_t index = 1U; index < cell_count; ++index) {
    faces[static_cast<std::size_t>(index)] =
        inverse_cumulative(coordinates, cumulative,
                           integral_step * static_cast<double>(index));
  }

  // Equidistribution through a sampled monitor can leave a small ratio
  // overshoot at a sharp focus transition.  Project face positions onto the
  // The two sweeps are deterministic and allocation-free.
  constexpr double projection_tolerance = 1.0e-13;
  for (std::size_t iteration = 0U; iteration < 1024U; ++iteration) {
    bool changed = false;
    for (std::size_t index = 1U; index + 1U < faces.size(); ++index) {
      const double left = faces[index] - faces[index - 1U];
      const double right = faces[index + 1U] - faces[index];
      if (left > growth_ratio * right * (1.0 + projection_tolerance)) {
        const double replacement =
            (faces[index - 1U] + growth_ratio * faces[index + 1U]) /
            (1.0 + growth_ratio);
        if (replacement < faces[index]) {
          faces[index] = replacement;
          changed = true;
        }
      }
    }
    for (std::size_t index = faces.size() - 2U; index > 0U; --index) {
      const double left = faces[index] - faces[index - 1U];
      const double right = faces[index + 1U] - faces[index];
      if (right > growth_ratio * left * (1.0 + projection_tolerance)) {
        const double replacement =
            (growth_ratio * faces[index - 1U] + faces[index + 1U]) /
            (1.0 + growth_ratio);
        if (replacement > faces[index]) {
          faces[index] = replacement;
          changed = true;
        }
      }
    }
    if (!changed) {
      break;
    }
  }

  const double tolerance = 128.0 * std::numeric_limits<double>::epsilon();
  for (std::size_t index = 0U; index + 1U < faces.size(); ++index) {
    const double width = faces[index + 1U] - faces[index];
    if (!(width > 0.0) || !std::isfinite(width) ||
        width + tolerance * std::max(1.0, std::abs(width)) < minimum_spacing) {
      return FaceBuildResult::too_fine;
    }
    if (width > base_spacing * (1.0 + tolerance)) {
      return FaceBuildResult::too_coarse;
    }
  }
  for (std::size_t index = 1U; index + 1U < faces.size(); ++index) {
    const double left = faces[index] - faces[index - 1U];
    const double right = faces[index + 1U] - faces[index];
    const double ratio = std::max(left / right, right / left);
    if (ratio > growth_ratio * (1.0 + kConstraintTolerance)) {
      return FaceBuildResult::too_coarse;
    }
  }

  // Check all overlapping focus targets with the same start-event/min-heap
  // sweep used to construct the monitor.  This is O((N+F) log F), not F*N.
  std::vector<IntervalTarget> ordered_focuses = focuses;
  std::sort(ordered_focuses.begin(), ordered_focuses.end(), FocusStartLess{});
  std::priority_queue<ActiveFocus, std::vector<ActiveFocus>,
                      ActiveFocusGreater>
      active;
  std::size_t next_focus = 0U;
  for (std::size_t index = 0U; index + 1U < faces.size(); ++index) {
    const double centre = 0.5 * (faces[index] + faces[index + 1U]);
    while (next_focus < ordered_focuses.size() &&
           ordered_focuses[next_focus].lower <= centre) {
      active.push({ordered_focuses[next_focus].spacing,
                   ordered_focuses[next_focus].upper, next_focus});
      ++next_focus;
    }
    while (!active.empty() && active.top().upper < centre) {
      active.pop();
    }
    if (!active.empty()) {
      const double width = faces[index + 1U] - faces[index];
      if (width > active.top().spacing * (1.0 + 16.0 * tolerance)) {
        return FaceBuildResult::too_coarse;
      }
    }
  }
  return FaceBuildResult::feasible;
}

bool create_stretched_axis(double lower, double upper, double base_spacing,
                           double minimum_spacing, double growth_ratio,
                           const std::vector<IntervalTarget>& focuses,
                           bool has_exact_cells, std::int32_t exact_cells,
                           std::uint64_t maximum_count,
                           std::uint64_t persistent_face_bytes,
                           std::uint64_t max_peak_payload_bytes,
                           std::vector<double>& faces) {
  // Before metric sampling, bound its actual worst-case payload: four dense
  // double arrays plus knots/focus copies and heap storage.  The fixed sample
  // hard cap remains independent from the case memory budget.
  const long double raw_samples =
      std::ceil((static_cast<long double>(upper) - lower) /
                std::max(static_cast<long double>(base_spacing) * 0.25L,
                         static_cast<long double>(minimum_spacing) * 0.25L));
  const long double interval_count =
      static_cast<long double>(focuses.size() * 2U + 1U);
  const long double conservative_samples = raw_samples + 8.0L * interval_count;
  const std::uint64_t estimated_samples =
      conservative_samples >= static_cast<long double>(kMaximumMetricSamples)
          ? kMaximumMetricSamples
          : static_cast<std::uint64_t>(std::ceil(conservative_samples));
  std::uint64_t sample_bytes = 0U;
  std::uint64_t focus_bytes = 0U;
  if (!checked_multiply(estimated_samples + 1U, 4U * sizeof(double),
                        sample_bytes) ||
      // Simultaneous focus-related storage includes caller focus, sorted
      // envelope copy, target-check copy, heap backing, and 2F+2 knots.
      !checked_multiply(static_cast<std::uint64_t>(focuses.size()),
                        8U * sizeof(IntervalTarget), focus_bytes) ||
      !checked_add(sample_bytes, focus_bytes, sample_bytes) ||
      !within_peak(persistent_face_bytes, sample_bytes,
                   max_peak_payload_bytes)) {
    return false;
  }
  std::vector<double> coordinates;
  std::vector<double> cumulative;
  if (!make_metric_samples(lower, upper, base_spacing, minimum_spacing,
                           growth_ratio, focuses, coordinates, cumulative)) {
    return false;
  }

  const std::uint64_t monitor_minimum =
      static_cast<std::uint64_t>(std::ceil(cumulative.back()));
  std::uint64_t cell_count = has_exact_cells
                                 ? static_cast<std::uint64_t>(exact_cells)
                                 : monitor_minimum;
  if (cell_count == 0U || cell_count > maximum_count ||
      (has_exact_cells && cell_count < monitor_minimum)) {
    return false;
  }
  std::uint64_t candidate_bytes = 0U;
  if (!checked_multiply(cell_count + 1U, sizeof(double), candidate_bytes) ||
      !checked_add(sample_bytes, candidate_bytes, candidate_bytes) ||
      !within_peak(persistent_face_bytes, candidate_bytes,
                   max_peak_payload_bytes)) {
    return false;
  }
  if (has_exact_cells) {
    const std::uint64_t visits_per_iteration =
        cell_count > 1U ? 2U * (cell_count - 1U) : 0U;
    if (visits_per_iteration >
        kMaximumFaceProjectionVisits / 1024U) {
      return false;
    }
    return build_faces_for_count(lower, upper, base_spacing, minimum_spacing,
                                 growth_ratio, focuses, coordinates, cumulative,
                                 cell_count, faces) ==
           FaceBuildResult::feasible;
  }

  // Do not binary-search this predicate: the minimum-spacing and focus/growth
  // constraints together do not prove monotonic feasibility in N.  Starting
  // at the monitor lower bound, visit every candidate until the first feasible
  // plan, a physically too-fine plan, or the hard count bound.
  std::uint64_t remaining_visits = kMaximumFaceProjectionVisits;
  for (;;) {
    if (!checked_multiply(cell_count + 1U, sizeof(double), candidate_bytes) ||
        !checked_add(sample_bytes, candidate_bytes, candidate_bytes) ||
        !within_peak(persistent_face_bytes, candidate_bytes,
                     max_peak_payload_bytes)) {
      return false;
    }
    const std::uint64_t visits_per_iteration =
        cell_count > 1U ? 2U * (cell_count - 1U) : 0U;
    if (visits_per_iteration >
        remaining_visits / 1024U) {
      return false;
    }
    remaining_visits -= visits_per_iteration * 1024U;
    const FaceBuildResult result = build_faces_for_count(
        lower, upper, base_spacing, minimum_spacing, growth_ratio, focuses,
        coordinates, cumulative, cell_count, faces);
    if (result == FaceBuildResult::feasible) {
      return true;
    }
    if (cell_count == maximum_count) {
      return false;
    }
    ++cell_count;
  }
}

std::uint64_t maximum_axis_cells(double lower, double upper,
                                 double minimum_spacing,
                                 std::uint64_t global_limit) noexcept {
  if (global_limit == 0U) {
    return 0U;
  }
  const long double span = static_cast<long double>(upper) -
                           static_cast<long double>(lower);
  const long double spacing = static_cast<long double>(minimum_spacing);
  if (!(span > 0.0L) || !(spacing > 0.0L)) {
    return 0U;
  }
  const long double ratio = span / spacing;
  const long double safe_ratio =
      ratio * (1.0L + static_cast<long double>(kConstraintTolerance));
  const std::uint64_t integer_limit = static_cast<std::uint64_t>(
      std::numeric_limits<std::int32_t>::max());
  std::uint64_t spacing_limit = integer_limit;
  if (safe_ratio < static_cast<long double>(integer_limit)) {
    spacing_limit = safe_ratio >= 1.0L
                        ? static_cast<std::uint64_t>(std::floor(safe_ratio))
                        : 0U;
  }
  return std::min(global_limit, spacing_limit);
}

}  // namespace

Status generate_cartesian_faces(
    const CartesianMeshSpec& mesh, std::uint64_t max_peak_payload_bytes,
    std::array<std::vector<double>, 3U>& faces) noexcept {
  try {
    std::uint64_t persistent_face_bytes = 0U;
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      const double lower = component(mesh.lower, axis);
      const double upper = component(mesh.upper, axis);
      const double minimum_spacing = component(mesh.minimum_spacing, axis);
      if (!std::isfinite(lower) || !std::isfinite(upper) ||
          !std::isfinite(minimum_spacing) || !(lower < upper) ||
          !(minimum_spacing > 0.0)) {
        return {StatusCode::invalid_plan, kFocusInput};
      }
      if (mesh.kind == GeometryKind::uniform) {
        const std::int32_t cells = component(mesh.exact_cells, axis);
        if (!mesh.has_exact_cells || cells <= 0 || mesh.has_base_spacing ||
            !mesh.focus_regions.empty() || mesh.max_growth_ratio != 1.0) {
          return {StatusCode::invalid_plan, kFocusInput};
        }
        const double width = (upper - lower) / static_cast<double>(cells);
        if (!(width >= minimum_spacing) || !std::isfinite(width)) {
          return {StatusCode::invalid_plan, kFocusSpacing};
        }
        std::uint64_t axis_bytes = 0U;
        std::uint64_t next_persistent = 0U;
        if (!checked_multiply(static_cast<std::uint64_t>(cells) + 1U,
                              sizeof(double), axis_bytes) ||
            !checked_add(persistent_face_bytes, axis_bytes, next_persistent) ||
            !within_peak(next_persistent, 0U, max_peak_payload_bytes)) {
          return {StatusCode::invalid_plan, kFocusMemory};
        }
        faces[axis].resize(static_cast<std::size_t>(cells) + 1U);
        persistent_face_bytes = next_persistent;
        for (std::int32_t index = 0; index <= cells; ++index) {
          faces[axis][static_cast<std::size_t>(index)] =
              index == cells
                  ? upper
                  : lower + width * static_cast<double>(index);
        }
        continue;
      }

      const double base_spacing = component(mesh.base_spacing, axis);
      const std::int32_t exact_cells = component(mesh.exact_cells, axis);
      const std::uint64_t maximum_count = maximum_axis_cells(
          lower, upper, minimum_spacing, mesh.limits.max_global_cells);
      if (!mesh.has_base_spacing || !std::isfinite(base_spacing) ||
          base_spacing < minimum_spacing ||
          !std::isfinite(mesh.max_growth_ratio) ||
          mesh.max_growth_ratio < 1.0 ||
          maximum_count == 0U ||
          (mesh.has_exact_cells &&
           (exact_cells <= 0 ||
            static_cast<std::uint64_t>(exact_cells) > maximum_count))) {
        return {StatusCode::invalid_plan, kFocusInput};
      }
      std::vector<IntervalTarget> focuses;
      focuses.reserve(mesh.focus_regions.size());
      for (const FocusRegionSpec& region : mesh.focus_regions) {
        const double focus_lower = component(region.lower, axis);
        const double focus_upper = component(region.upper, axis);
        const double focus_spacing = component(region.target_spacing, axis);
        if (!std::isfinite(focus_lower) || !std::isfinite(focus_upper) ||
            !std::isfinite(focus_spacing) || focus_lower < lower ||
            focus_upper > upper || !(focus_lower < focus_upper) ||
            focus_spacing < minimum_spacing || focus_spacing > base_spacing) {
          return {StatusCode::invalid_plan, kFocusInput};
        }
        focuses.push_back({focus_lower, focus_upper, focus_spacing});
      }
      if (!create_stretched_axis(lower, upper, base_spacing, minimum_spacing,
                                 mesh.max_growth_ratio, focuses,
                                 mesh.has_exact_cells, exact_cells, maximum_count,
                                 persistent_face_bytes,
                                 max_peak_payload_bytes,
                                 faces[axis])) {
        return {StatusCode::invalid_plan,
                mesh.has_exact_cells ? kFocusCellCount : kFocusGrowth};
      }
      std::uint64_t axis_bytes = 0U;
      if (!checked_multiply(static_cast<std::uint64_t>(faces[axis].size()),
                            sizeof(double), axis_bytes) ||
          !checked_add(persistent_face_bytes, axis_bytes,
                       persistent_face_bytes)) {
        return {StatusCode::invalid_plan, kFocusMemory};
      }
    }
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, 0U};
  } catch (...) {
    return {StatusCode::invalid_plan, kFocusInput};
  }
}

Status finish_axis_metrics(std::vector<double> faces,
                           AxisPayload& out) noexcept {
  try {
    if (faces.size() < 2U) {
      return {StatusCode::invalid_plan, kFocusCellCount};
    }
    AxisPayload candidate;
    candidate.faces = std::move(faces);
    const std::size_t cells = candidate.faces.size() - 1U;
    candidate.centres.resize(cells);
    candidate.widths.resize(cells);
    candidate.inverse_widths.resize(cells);
    candidate.inverse_centre_distances.resize(cells > 0U ? cells - 1U : 0U);
    for (std::size_t index = 0U; index < cells; ++index) {
      const double left = candidate.faces[index];
      const double right = candidate.faces[index + 1U];
      const double width = right - left;
      if (!std::isfinite(left) || !std::isfinite(right) || !(width > 0.0)) {
        return {StatusCode::invalid_plan, kFocusSpacing};
      }
      candidate.centres[index] = left + 0.5 * width;
      candidate.widths[index] = width;
      candidate.inverse_widths[index] = 1.0 / width;
    }
    for (std::size_t index = 0U;
         index < candidate.inverse_centre_distances.size(); ++index) {
      candidate.inverse_centre_distances[index] =
          1.0 / (candidate.centres[index + 1U] - candidate.centres[index]);
    }
    candidate.uniform = true;
    candidate.uniform_width = candidate.widths.front();
    candidate.uniform_inverse_width = candidate.inverse_widths.front();
    const double tolerance =
        64.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, std::abs(candidate.uniform_width));
    for (const double width : candidate.widths) {
      if (std::abs(width - candidate.uniform_width) > tolerance) {
        candidate.uniform = false;
        candidate.uniform_width = 0.0;
        candidate.uniform_inverse_width = 0.0;
        break;
      }
    }
    out = std::move(candidate);
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, 0U};
  } catch (...) {
    return {StatusCode::invalid_plan, kFocusSpacing};
  }
}

}  // namespace hundun::v04::detail

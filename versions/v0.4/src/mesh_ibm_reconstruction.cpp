// SPDX-License-Identifier: Apache-2.0

#include "hundun/v04_ibm.hpp"

#include "mesh_ibm_qr_detail.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <type_traits>
#include <utility>
#include <vector>

namespace hundun::v04 {
namespace {

constexpr std::uint32_t kQuadraticInput = 1301U;
constexpr std::uint32_t kQuadraticDonors = 1302U;
constexpr std::uint32_t kQuadraticGeometry = 1303U;
constexpr std::uint32_t kQuadraticCoverage = 1304U;
constexpr std::uint32_t kQuadraticRank = 1305U;
constexpr std::uint32_t kQuadraticCondition = 1306U;
constexpr std::uint32_t kQuadraticFunctional = 1307U;
constexpr std::uint32_t kQuadraticOverflow = 1308U;
constexpr std::uint32_t kQuadraticEvaluation = 1309U;
constexpr std::uint64_t kFnvOffset = UINT64_C(14695981039346656037);
constexpr std::uint64_t kFnvPrime = UINT64_C(1099511628211);
constexpr std::size_t kBasisSize = 10U;
constexpr std::uint8_t kHardMinimumDonors = 14U;
constexpr std::uint8_t kHardMaximumDonors = 32U;
constexpr std::uint8_t kHardMaximumReach = 4U;
constexpr std::uint8_t kHardMinimumNormalBands = 3U;
constexpr double kHardConditionLimit = 1.0e8;

static_assert(sizeof(double) == sizeof(std::uint64_t) &&
                  std::numeric_limits<double>::is_iec559,
              "quadratic fingerprint requires IEEE-754 binary64");
static_assert(std::is_nothrow_move_assignable_v<QuadraticStencilPlan>,
              "quadratic publication must not throw");

class Hash64 {
 public:
  void integer(std::uint64_t value) noexcept {
    for (std::size_t byte = 0U; byte < sizeof(value); ++byte) {
      value_ ^= (value >> (byte * 8U)) & UINT64_C(0xff);
      value_ *= kFnvPrime;
    }
  }

  void signed_integer(std::int32_t value) noexcept {
    integer(static_cast<std::uint32_t>(value));
  }

  void real(double value) noexcept {
    std::uint64_t bits = 0U;
    std::memcpy(&bits, &value, sizeof(bits));
    integer(bits);
  }

  PlanFingerprint finish() const noexcept { return value_ == 0U ? 1U : value_; }

 private:
  std::uint64_t value_{kFnvOffset};
};

bool finite(Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

bool same(Int3 left, Int3 right) noexcept {
  return left.x == right.x && left.y == right.y && left.z == right.z;
}

double dot(Real3 left, Real3 right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

Real3 subtract(Real3 left, Real3 right) noexcept {
  return {left.x - right.x, left.y - right.y, left.z - right.z};
}

double norm(Real3 value) noexcept {
  return std::sqrt(dot(value, value));
}

Real3 cross(Real3 left, Real3 right) noexcept {
  return {left.y * right.z - left.z * right.y,
          left.z * right.x - left.x * right.z,
          left.x * right.y - left.y * right.x};
}

bool valid_frame(const QuadraticFrame& frame) noexcept {
  if (!finite(frame.origin) || !finite(frame.normal) ||
      !finite(frame.tangent1) || !finite(frame.tangent2) ||
      !std::isfinite(frame.scale) || !(frame.scale > 0.0)) {
    return false;
  }
  constexpr double kFrameTolerance = 512.0 *
                                     std::numeric_limits<double>::epsilon();
  return std::abs(norm(frame.normal) - 1.0) <= kFrameTolerance &&
         std::abs(norm(frame.tangent1) - 1.0) <= kFrameTolerance &&
         std::abs(norm(frame.tangent2) - 1.0) <= kFrameTolerance &&
         std::abs(dot(frame.normal, frame.tangent1)) <= kFrameTolerance &&
         std::abs(dot(frame.normal, frame.tangent2)) <= kFrameTolerance &&
         std::abs(dot(frame.tangent1, frame.tangent2)) <= kFrameTolerance &&
         std::abs(dot(cross(frame.normal, frame.tangent1), frame.tangent2) -
                  1.0) <= kFrameTolerance;
}

std::uint8_t logical_reach(const QuadraticFrame& frame,
                           const QuadraticDonorCell& donor,
                           bool& valid) noexcept {
  const auto difference = [](std::int32_t left, std::int32_t right,
                             std::uint64_t& out) noexcept {
    const std::int64_t delta =
        static_cast<std::int64_t>(left) - static_cast<std::int64_t>(right);
    out = static_cast<std::uint64_t>(delta < 0 ? -delta : delta);
  };
  std::uint64_t x = 0U;
  std::uint64_t y = 0U;
  std::uint64_t z = 0U;
  difference(donor.global_index.x, frame.anchor_global_cell.x, x);
  difference(donor.global_index.y, frame.anchor_global_cell.y, y);
  difference(donor.global_index.z, frame.anchor_global_cell.z, z);
  const std::uint64_t maximum = std::max({x, y, z});
  valid = maximum <= std::numeric_limits<std::uint8_t>::max();
  return valid ? static_cast<std::uint8_t>(maximum) : UINT8_MAX;
}

struct LocalCoordinates {
  double n{};
  double t1{};
  double t2{};
};

LocalCoordinates coordinates(const QuadraticFrame& frame,
                             Real3 point) noexcept {
  const Real3 offset = subtract(point, frame.origin);
  const double inverse_scale = 1.0 / frame.scale;
  return {dot(offset, frame.normal) * inverse_scale,
          dot(offset, frame.tangent1) * inverse_scale,
          dot(offset, frame.tangent2) * inverse_scale};
}

std::array<double, kBasisSize> donor_basis(
    const QuadraticFrame& frame,
    const QuadraticDonorCell& donor) noexcept {
  const LocalCoordinates centre = coordinates(frame, donor.centre);
  const double inverse_scale_squared = 1.0 / (frame.scale * frame.scale);
  const double variance_x = donor.widths.x * donor.widths.x / 12.0;
  const double variance_y = donor.widths.y * donor.widths.y / 12.0;
  const double variance_z = donor.widths.z * donor.widths.z / 12.0;
  const auto covariance = [&](Real3 left, Real3 right) noexcept {
    return (left.x * right.x * variance_x +
            left.y * right.y * variance_y +
            left.z * right.z * variance_z) *
           inverse_scale_squared;
  };
  return {1.0,
          centre.n,
          centre.t1,
          centre.t2,
          centre.n * centre.n + covariance(frame.normal, frame.normal),
          centre.n * centre.t1 + covariance(frame.normal, frame.tangent1),
          centre.n * centre.t2 + covariance(frame.normal, frame.tangent2),
          centre.t1 * centre.t1 +
              covariance(frame.tangent1, frame.tangent1),
          centre.t1 * centre.t2 +
              covariance(frame.tangent1, frame.tangent2),
          centre.t2 * centre.t2 +
              covariance(frame.tangent2, frame.tangent2)};
}

std::array<double, kBasisSize> value_functional(
    const QuadraticFrame& frame, Real3 point) noexcept {
  const LocalCoordinates value = coordinates(frame, point);
  return {1.0, value.n, value.t1, value.t2, value.n * value.n,
          value.n * value.t1, value.n * value.t2,
          value.t1 * value.t1, value.t1 * value.t2,
          value.t2 * value.t2};
}

std::array<double, kBasisSize> derivative_functional(
    const QuadraticFrame& frame, Real3 point, Real3 direction) noexcept {
  const LocalCoordinates value = coordinates(frame, point);
  const double inverse_scale = 1.0 / frame.scale;
  const LocalCoordinates derivative{
      dot(direction, frame.normal) * inverse_scale,
      dot(direction, frame.tangent1) * inverse_scale,
      dot(direction, frame.tangent2) * inverse_scale};
  return {0.0,
          derivative.n,
          derivative.t1,
          derivative.t2,
          2.0 * value.n * derivative.n,
          value.n * derivative.t1 + value.t1 * derivative.n,
          value.n * derivative.t2 + value.t2 * derivative.n,
          2.0 * value.t1 * derivative.t1,
          value.t1 * derivative.t2 + value.t2 * derivative.t1,
          2.0 * value.t2 * derivative.t2};
}

bool finite_basis(const std::array<double, kBasisSize>& values) noexcept {
  return std::all_of(values.begin(), values.end(),
                     [](double value) { return std::isfinite(value); });
}

void subtract_constraint(std::array<double, kBasisSize>& functional,
                         QuadraticConstraint constraint,
                         double& wall_value_weight,
                         double& wall_gradient_weight_m) noexcept {
  wall_value_weight = 0.0;
  wall_gradient_weight_m = 0.0;
  if (constraint == QuadraticConstraint::origin_value) {
    wall_value_weight = functional[0U];
    functional[0U] = 0.0;
  } else if (constraint == QuadraticConstraint::origin_normal_gradient) {
    // The n coefficient is prescribed in normalized coordinates.  Donor
    // cell averages therefore require a separate constrained least-squares
    // solve; merely removing f_n from the full solve would not eliminate the
    // donor n moment from A^T w.
  }
}

std::uint8_t quadrant(double t1, double t2) noexcept {
  const std::uint8_t high_t1 = t1 >= 0.0 ? 1U : 0U;
  const std::uint8_t high_t2 = t2 >= 0.0 ? 1U : 0U;
  return static_cast<std::uint8_t>(1U << (high_t1 * 2U + high_t2));
}

std::uint8_t count_normal_bands(
    const std::array<double, detail::ibm_qr::kMaximumRows>& values,
    std::size_t count) noexcept {
  std::array<double, detail::ibm_qr::kMaximumRows> sorted = values;
  std::sort(sorted.begin(), sorted.begin() + count);
  std::uint8_t bands = 0U;
  double previous_band = 0.0;
  constexpr double kMultiplier = 512.0;
  const double epsilon = std::numeric_limits<double>::epsilon();
  for (std::size_t index = 0U; index < count; ++index) {
    const double tolerance =
        kMultiplier * epsilon * std::max(1.0, std::abs(sorted[index]));
    if (bands == 0U || sorted[index] - previous_band > tolerance) {
      ++bands;
      previous_band = sorted[index];
    }
  }
  return bands;
}

bool valid_limits(QuadraticStencilLimits limits) noexcept {
  return limits.minimum_donors >= kHardMinimumDonors &&
         limits.minimum_donors <= limits.maximum_donors &&
         limits.maximum_donors <= kHardMaximumDonors &&
         limits.maximum_reach > 0U &&
         limits.maximum_reach <= kHardMaximumReach &&
         limits.minimum_normal_bands >= kHardMinimumNormalBands &&
         std::isfinite(limits.condition_limit) &&
         limits.condition_limit >= 1.0 &&
         limits.condition_limit <= kHardConditionLimit;
}

void hash_int3(Hash64& hash, Int3 value) noexcept {
  hash.signed_integer(value.x);
  hash.signed_integer(value.y);
  hash.signed_integer(value.z);
}

void hash_real3(Hash64& hash, Real3 value) noexcept {
  hash.real(value.x);
  hash.real(value.y);
  hash.real(value.z);
}

PlanFingerprint pivot_fingerprint(
    const detail::ibm_qr::Factorization& factor) noexcept {
  Hash64 hash;
  hash.integer(factor.columns);
  for (std::size_t index = 0U; index < factor.columns; ++index) {
    hash.integer(factor.permutation[index]);
  }
  return hash.finish();
}

bool checked_add(std::size_t left, std::size_t right,
                 std::size_t& out) noexcept {
  if (right > std::numeric_limits<std::size_t>::max() - left) {
    return false;
  }
  out = left + right;
  return true;
}

bool checked_multiply(std::size_t left, std::size_t right,
                      std::size_t& out) noexcept {
  if (left != 0U && right > std::numeric_limits<std::size_t>::max() / left) {
    return false;
  }
  out = left * right;
  return true;
}

Status validate_request(const QuadraticStencilRequest& request,
                        QuadraticStencilLimits limits) noexcept {
  if (!valid_frame(request.frame) || request.donors.data == nullptr ||
      request.functionals.data == nullptr || request.functionals.size == 0U ||
      request.functionals.size > UINT8_MAX) {
    return {StatusCode::invalid_plan, kQuadraticInput};
  }
  if (request.donors.size < limits.minimum_donors ||
      request.donors.size > limits.maximum_donors) {
    return {StatusCode::invalid_plan, kQuadraticDonors};
  }
  return {};
}

struct BuildPayload {
  std::vector<QuadraticStencilGroup> groups;
  std::vector<QuadraticAffineRow> rows;
  std::vector<GlobalCellId> donor_global_cells;
  std::vector<Int3> donor_local_indices;
  std::vector<double> weights;
  std::uint8_t maximum_halo_reach{};

  void clear() noexcept {
    groups.clear();
    rows.clear();
    donor_global_cells.clear();
    donor_local_indices.clear();
    weights.clear();
    maximum_halo_reach = 0U;
  }
};

struct DonorPool {
  std::array<std::size_t, detail::ibm_qr::kMaximumRows> order{};
  std::array<LocalCoordinates, detail::ibm_qr::kMaximumRows> local{};
  std::array<std::array<double, kBasisSize>,
             detail::ibm_qr::kMaximumRows>
      basis{};
  std::array<double, detail::ibm_qr::kMaximumRows> distance_squared{};
  std::array<std::uint8_t, detail::ibm_qr::kMaximumRows> quadrant_bits{};
  std::array<std::uint8_t, detail::ibm_qr::kMaximumRows> normal_band{};
  std::array<std::uint8_t, detail::ibm_qr::kMaximumRows> reach{};
  std::size_t size{};
};

struct CompileWorkspace {
  BuildPayload trial;
  BuildPayload best;
};

bool canonical_donor_less(const QuadraticDonorCell& left,
                          const QuadraticDonorCell& right) noexcept {
  if (left.global_cell != right.global_cell) {
    return left.global_cell < right.global_cell;
  }
  if (left.global_index.x != right.global_index.x) {
    return left.global_index.x < right.global_index.x;
  }
  if (left.global_index.y != right.global_index.y) {
    return left.global_index.y < right.global_index.y;
  }
  return left.global_index.z < right.global_index.z;
}

Status prepare_donor_pool(const QuadraticStencilRequest& request,
                          QuadraticStencilLimits limits,
                          DonorPool& pool) noexcept {
  const Status preflight = validate_request(request, limits);
  if (!preflight) {
    return preflight;
  }
  pool = {};
  pool.size = request.donors.size;
  std::array<std::size_t, detail::ibm_qr::kMaximumRows> canonical{};
  for (std::size_t index = 0U; index < request.donors.size; ++index) {
    canonical[index] = index;
  }
  std::sort(canonical.begin(), canonical.begin() + request.donors.size,
            [&](std::size_t left, std::size_t right) {
              const auto& a = request.donors.data[left];
              const auto& b = request.donors.data[right];
              return canonical_donor_less(a, b);
            });
  GlobalCellId previous_global = 0U;
  bool have_previous = false;
  for (std::size_t row = 0U; row < request.donors.size; ++row) {
    const std::size_t source = canonical[row];
    const QuadraticDonorCell& donor = request.donors.data[source];
    if (have_previous && donor.global_cell == previous_global) {
      return {StatusCode::invalid_plan, kQuadraticDonors};
    }
    have_previous = true;
    previous_global = donor.global_cell;
    for (std::size_t previous = 0U; previous < row; ++previous) {
      const QuadraticDonorCell& prior =
          request.donors.data[canonical[previous]];
      if (same(donor.global_index, prior.global_index) ||
          same(donor.local_index, prior.local_index)) {
        return {StatusCode::invalid_plan, kQuadraticDonors};
      }
    }
    if (!finite(donor.centre) || !finite(donor.widths) ||
        !(donor.widths.x > 0.0) || !(donor.widths.y > 0.0) ||
        !(donor.widths.z > 0.0)) {
      return {StatusCode::invalid_plan, kQuadraticGeometry};
    }
    bool reach_valid = false;
    const std::uint8_t reach =
        logical_reach(request.frame, donor, reach_valid);
    if (!reach_valid || reach > limits.maximum_reach) {
      return {StatusCode::invalid_plan, kQuadraticCoverage};
    }
    const LocalCoordinates local = coordinates(request.frame, donor.centre);
    if (!std::isfinite(local.n) || !std::isfinite(local.t1) ||
        !std::isfinite(local.t2) || !(local.n > 0.0)) {
      return {StatusCode::invalid_plan, kQuadraticCoverage};
    }
    const auto basis = donor_basis(request.frame, donor);
    if (!finite_basis(basis)) {
      return {StatusCode::invalid_plan, kQuadraticGeometry};
    }
    const double distance_squared =
        local.n * local.n + local.t1 * local.t1 + local.t2 * local.t2;
    if (!std::isfinite(distance_squared)) {
      return {StatusCode::invalid_plan, kQuadraticGeometry};
    }
    pool.local[source] = local;
    pool.basis[source] = basis;
    pool.distance_squared[source] = distance_squared;
    pool.quadrant_bits[source] = quadrant(local.t1, local.t2);
    pool.reach[source] = reach;
  }

  std::array<std::size_t, detail::ibm_qr::kMaximumRows> by_normal = canonical;
  std::sort(by_normal.begin(), by_normal.begin() + request.donors.size,
            [&](std::size_t left, std::size_t right) {
              if (pool.local[left].n != pool.local[right].n) {
                return pool.local[left].n < pool.local[right].n;
              }
              return canonical_donor_less(request.donors.data[left],
                                          request.donors.data[right]);
            });
  std::uint8_t band = 0U;
  double previous_band = 0.0;
  constexpr double kBandMultiplier = 512.0;
  const double epsilon = std::numeric_limits<double>::epsilon();
  for (std::size_t index = 0U; index < request.donors.size; ++index) {
    const std::size_t source = by_normal[index];
    const double value = pool.local[source].n;
    const double tolerance =
        kBandMultiplier * epsilon * std::max(1.0, std::abs(value));
    if (index == 0U || value - previous_band > tolerance) {
      previous_band = value;
      band = static_cast<std::uint8_t>(index == 0U ? 0U : band + 1U);
    }
    pool.normal_band[source] = band;
  }

  // Build one canonical nested family.  Coverage gain is resolved before
  // distance, so a short prefix cannot accidentally omit a tangential
  // quadrant or a distinct positive-normal layer.  Every remaining tie has
  // a total physical/global key and is therefore independent of caller
  // permutation.
  std::array<bool, detail::ibm_qr::kMaximumRows> used{};
  std::uint32_t selected_bands = 0U;
  std::uint8_t selected_quadrants = 0U;
  for (std::size_t position = 0U; position < request.donors.size; ++position) {
    std::size_t selected = request.donors.size;
    int selected_gain = -1;
    bool selected_quadrant_gain = false;
    bool selected_band_gain = false;
    for (std::size_t source = 0U; source < request.donors.size; ++source) {
      if (used[source]) {
        continue;
      }
      const bool quadrant_gain =
          (selected_quadrants & pool.quadrant_bits[source]) == 0U;
      const std::uint32_t band_bit = UINT32_C(1)
                                     << pool.normal_band[source];
      const bool band_gain = (selected_bands & band_bit) == 0U;
      const int gain = static_cast<int>(quadrant_gain) +
                       static_cast<int>(band_gain);
      bool better = selected == request.donors.size || gain > selected_gain;
      if (!better && gain == selected_gain &&
          quadrant_gain != selected_quadrant_gain) {
        better = quadrant_gain;
      }
      if (!better && gain == selected_gain &&
          quadrant_gain == selected_quadrant_gain &&
          band_gain == selected_band_gain) {
        if (pool.distance_squared[source] !=
            pool.distance_squared[selected]) {
          better = pool.distance_squared[source] <
                   pool.distance_squared[selected];
        } else {
          better = canonical_donor_less(request.donors.data[source],
                                        request.donors.data[selected]);
        }
      }
      if (better) {
        selected = source;
        selected_gain = gain;
        selected_quadrant_gain = quadrant_gain;
        selected_band_gain = band_gain;
      }
    }
    pool.order[position] = selected;
    used[selected] = true;
    selected_quadrants = static_cast<std::uint8_t>(
        selected_quadrants | pool.quadrant_bits[selected]);
    selected_bands |= UINT32_C(1) << pool.normal_band[selected];
  }
  return {};
}

Status validate_functionals(const QuadraticStencilRequest& request) noexcept {
  for (std::size_t row = 0U; row < request.functionals.size; ++row) {
    const QuadraticFunctionalRequest& functional =
        request.functionals.data[row];
    if (!finite(functional.point) ||
        functional.kind > QuadraticFunctionalKind::directional_derivative ||
        functional.constraint >
            QuadraticConstraint::origin_normal_gradient ||
        (functional.kind == QuadraticFunctionalKind::directional_derivative &&
         (!finite(functional.direction) ||
          !(norm(functional.direction) > 0.0)))) {
      return {StatusCode::invalid_plan, kQuadraticFunctional};
    }
  }
  return {};
}

Status compile_prefix(const QuadraticStencilRequest& request,
                      const DonorPool& pool, std::size_t donor_count,
                      QuadraticStencilLimits limits,
                      BuildPayload& candidate) {
  candidate.clear();
  std::array<std::size_t, detail::ibm_qr::kMaximumRows> selected_order{};
  std::copy_n(pool.order.begin(), donor_count, selected_order.begin());
  std::sort(selected_order.begin(), selected_order.begin() + donor_count,
            [&](std::size_t left, std::size_t right) {
              return canonical_donor_less(request.donors.data[left],
                                          request.donors.data[right]);
            });
  std::array<double, detail::ibm_qr::kMaximumRows * kBasisSize> matrix{};
  std::array<double, detail::ibm_qr::kMaximumRows> normal_coordinates{};
  std::uint8_t quadrant_mask = 0U;
  std::uint8_t maximum_reach = 0U;
  for (std::size_t row = 0U; row < donor_count; ++row) {
    const std::size_t source = selected_order[row];
    const QuadraticDonorCell& donor = request.donors.data[source];
    normal_coordinates[row] = pool.local[source].n;
    quadrant_mask = static_cast<std::uint8_t>(
        quadrant_mask | pool.quadrant_bits[source]);
    maximum_reach = std::max(maximum_reach, pool.reach[source]);
    std::copy(pool.basis[source].begin(), pool.basis[source].end(),
              matrix.begin() +
                  static_cast<std::ptrdiff_t>(row * kBasisSize));
    candidate.donor_global_cells.push_back(donor.global_cell);
    candidate.donor_local_indices.push_back(donor.local_index);
  }
  const std::uint8_t normal_bands =
      count_normal_bands(normal_coordinates, donor_count);
  if (normal_bands < limits.minimum_normal_bands || quadrant_mask != 0x0fU) {
    return {StatusCode::invalid_plan, kQuadraticCoverage};
  }
  detail::ibm_qr::Factorization factor;
  if (!detail::ibm_qr::factorize(matrix.data(), donor_count,
                                 kBasisSize, limits.condition_limit, factor)) {
    return {StatusCode::invalid_plan,
            factor.rank < kBasisSize ? kQuadraticRank : kQuadraticCondition};
  }

  const PlanFingerprint qr_fingerprint = pivot_fingerprint(factor);
  constexpr std::array<std::size_t, kBasisSize - 1U> value_free_basis{
      1U, 2U, 3U, 4U, 5U, 6U, 7U, 8U, 9U};
  constexpr std::array<std::size_t, kBasisSize - 1U> gradient_free_basis{
      0U, 2U, 3U, 4U, 5U, 6U, 7U, 8U, 9U};
  std::array<double, detail::ibm_qr::kMaximumRows *
                         (kBasisSize - 1U)>
      value_constrained_matrix{};
  std::array<double, detail::ibm_qr::kMaximumRows *
                         (kBasisSize - 1U)>
      gradient_constrained_matrix{};
  for (std::size_t donor = 0U; donor < donor_count; ++donor) {
    for (std::size_t column = 0U; column < kBasisSize - 1U; ++column) {
      value_constrained_matrix[donor * (kBasisSize - 1U) + column] =
          matrix[donor * kBasisSize + value_free_basis[column]];
      gradient_constrained_matrix[donor * (kBasisSize - 1U) + column] =
          matrix[donor * kBasisSize + gradient_free_basis[column]];
    }
  }
  detail::ibm_qr::Factorization value_constrained_factor;
  detail::ibm_qr::Factorization gradient_constrained_factor;
  if (!detail::ibm_qr::factorize(value_constrained_matrix.data(),
                                 donor_count, kBasisSize - 1U,
                                 limits.condition_limit,
                                 value_constrained_factor) ||
      !detail::ibm_qr::factorize(gradient_constrained_matrix.data(),
                                 donor_count, kBasisSize - 1U,
                                 limits.condition_limit,
                                 gradient_constrained_factor)) {
    return {StatusCode::invalid_plan,
            value_constrained_factor.rank < kBasisSize - 1U ||
                    gradient_constrained_factor.rank < kBasisSize - 1U
                ? kQuadraticRank
                : kQuadraticCondition};
  }
  const double condition_estimate =
      std::max({factor.condition_estimate,
                value_constrained_factor.condition_estimate,
                gradient_constrained_factor.condition_estimate});
  double maximum_functional_l1 = 0.0;
  Hash64 group_hash;
  group_hash.integer(donor_count);
  group_hash.integer(request.functionals.size);
  hash_real3(group_hash, request.frame.origin);
  hash_real3(group_hash, request.frame.normal);
  hash_real3(group_hash, request.frame.tangent1);
  hash_real3(group_hash, request.frame.tangent2);
  group_hash.real(request.frame.scale);
  hash_int3(group_hash, request.frame.anchor_global_cell);
  group_hash.integer(qr_fingerprint);
  group_hash.real(condition_estimate);

  for (std::size_t row_index = 0U;
       row_index < request.functionals.size; ++row_index) {
    const QuadraticFunctionalRequest& request_row =
        request.functionals.data[row_index];
    std::array<double, kBasisSize> functional =
        request_row.kind == QuadraticFunctionalKind::value
            ? value_functional(request.frame, request_row.point)
            : derivative_functional(request.frame, request_row.point,
                                    request_row.direction);
    double wall_value_weight = 0.0;
    double wall_gradient_weight = 0.0;
    subtract_constraint(functional, request_row.constraint, wall_value_weight,
                        wall_gradient_weight);
    if (!finite_basis(functional)) {
      return {StatusCode::invalid_plan, kQuadraticFunctional};
    }
    std::array<double, detail::ibm_qr::kMaximumRows> weights{};
    bool solved = false;
    if (request_row.constraint == QuadraticConstraint::none) {
      solved = detail::ibm_qr::solve_transpose_minimum_norm(
          factor, functional.data(), weights.data());
    } else if (request_row.constraint == QuadraticConstraint::origin_value) {
      std::array<double, kBasisSize - 1U> reduced{};
      for (std::size_t column = 0U; column < reduced.size(); ++column) {
        reduced[column] = functional[value_free_basis[column]];
      }
      solved = detail::ibm_qr::solve_transpose_minimum_norm(
          value_constrained_factor, reduced.data(), weights.data());
      if (solved) {
        double constrained_moment = 0.0;
        for (std::size_t donor = 0U; donor < donor_count; ++donor) {
          constrained_moment += weights[donor] *
                                matrix[donor * kBasisSize + 0U];
        }
        wall_value_weight -= constrained_moment;
      }
    } else {
      std::array<double, kBasisSize - 1U> reduced{};
      for (std::size_t column = 0U; column < reduced.size(); ++column) {
        reduced[column] = functional[gradient_free_basis[column]];
      }
      solved = detail::ibm_qr::solve_transpose_minimum_norm(
          gradient_constrained_factor, reduced.data(), weights.data());
      if (solved) {
        double constrained_moment = 0.0;
        for (std::size_t donor = 0U; donor < donor_count; ++donor) {
          constrained_moment += weights[donor] *
                                matrix[donor * kBasisSize + 1U];
        }
        wall_gradient_weight =
            request.frame.scale * (functional[1U] - constrained_moment);
      }
    }
    if (!solved || !std::isfinite(wall_value_weight) ||
        !std::isfinite(wall_gradient_weight)) {
      return {StatusCode::numerical_failure, kQuadraticFunctional};
    }
    const std::size_t weight_begin = candidate.weights.size();
    if (weight_begin > UINT32_MAX) {
      return {StatusCode::invalid_plan, kQuadraticOverflow};
    }
    double l1 = std::abs(wall_value_weight) +
                std::abs(wall_gradient_weight);
    for (std::size_t donor = 0U; donor < donor_count; ++donor) {
      candidate.weights.push_back(weights[donor]);
      l1 += std::abs(weights[donor]);
      group_hash.real(weights[donor]);
    }
    if (!std::isfinite(l1)) {
      return {StatusCode::numerical_failure, kQuadraticFunctional};
    }
    maximum_functional_l1 = std::max(maximum_functional_l1, l1);
    candidate.rows.push_back(
        {static_cast<std::uint32_t>(candidate.groups.size()),
         static_cast<std::uint32_t>(weight_begin), wall_value_weight,
         wall_gradient_weight});
    group_hash.integer(static_cast<std::uint8_t>(request_row.kind));
    group_hash.integer(static_cast<std::uint8_t>(request_row.constraint));
    hash_real3(group_hash, request_row.point);
    hash_real3(group_hash, request_row.direction);
    group_hash.real(wall_value_weight);
    group_hash.real(wall_gradient_weight);
  }

  for (std::size_t donor = 0U;
       donor < candidate.donor_global_cells.size(); ++donor) {
    group_hash.integer(candidate.donor_global_cells[donor]);
    hash_int3(group_hash, candidate.donor_local_indices[donor]);
  }
  QuadraticStencilGroup group;
  group.donor_begin = 0U;
  group.row_begin = 0U;
  group.row_count = static_cast<std::uint8_t>(request.functionals.size);
  group.quality = {
      static_cast<std::uint8_t>(donor_count), normal_bands,
      quadrant_mask, maximum_reach, factor.rank, condition_estimate,
      maximum_functional_l1, qr_fingerprint};
  group.fingerprint = group_hash.finish();
  candidate.groups.push_back(group);
  candidate.maximum_halo_reach =
      std::max(candidate.maximum_halo_reach, maximum_reach);
  return {};
}

bool better_candidate(const BuildPayload& left,
                      const BuildPayload& right) noexcept {
  if (right.groups.empty()) {
    return true;
  }
  const QuadraticStencilGroup& a = left.groups.front();
  const QuadraticStencilGroup& b = right.groups.front();
  if (a.quality.functional_l1 != b.quality.functional_l1) {
    return a.quality.functional_l1 < b.quality.functional_l1;
  }
  if (a.quality.condition_estimate != b.quality.condition_estimate) {
    return a.quality.condition_estimate < b.quality.condition_estimate;
  }
  if (a.quality.donor_count != b.quality.donor_count) {
    return a.quality.donor_count < b.quality.donor_count;
  }
  return a.fingerprint < b.fingerprint;
}

Status append_group(const BuildPayload& source,
                    BuildPayload& destination) {
  if (source.groups.size() != 1U) {
    return {StatusCode::invalid_plan, kQuadraticInput};
  }
  const std::size_t group_index = destination.groups.size();
  const std::size_t donor_begin = destination.donor_global_cells.size();
  const std::size_t row_begin = destination.rows.size();
  const std::size_t weight_begin = destination.weights.size();
  if (group_index > UINT32_MAX || donor_begin > UINT32_MAX ||
      row_begin > UINT32_MAX || weight_begin > UINT32_MAX) {
    return {StatusCode::invalid_plan, kQuadraticOverflow};
  }
  destination.donor_global_cells.insert(
      destination.donor_global_cells.end(),
      source.donor_global_cells.begin(), source.donor_global_cells.end());
  destination.donor_local_indices.insert(
      destination.donor_local_indices.end(),
      source.donor_local_indices.begin(), source.donor_local_indices.end());
  destination.weights.insert(destination.weights.end(), source.weights.begin(),
                             source.weights.end());
  for (const QuadraticAffineRow& source_row : source.rows) {
    if (static_cast<std::size_t>(source_row.weight_begin) + weight_begin >
        UINT32_MAX) {
      return {StatusCode::invalid_plan, kQuadraticOverflow};
    }
    QuadraticAffineRow row = source_row;
    row.group = static_cast<std::uint32_t>(group_index);
    row.weight_begin = static_cast<std::uint32_t>(
        static_cast<std::size_t>(source_row.weight_begin) + weight_begin);
    destination.rows.push_back(row);
  }
  QuadraticStencilGroup group = source.groups.front();
  group.donor_begin = static_cast<std::uint32_t>(donor_begin);
  group.row_begin = static_cast<std::uint32_t>(row_begin);
  destination.groups.push_back(group);
  destination.maximum_halo_reach = std::max(
      destination.maximum_halo_reach, source.maximum_halo_reach);
  return {};
}

Status compile_group(const QuadraticStencilRequest& request,
                     QuadraticStencilLimits limits, BuildPayload& candidate,
                     CompileWorkspace& workspace) {
  DonorPool pool;
  const Status pool_status = prepare_donor_pool(request, limits, pool);
  if (!pool_status) {
    return pool_status;
  }
  const Status functional_status = validate_functionals(request);
  if (!functional_status) {
    return functional_status;
  }
  const std::size_t weights = request.donors.size * request.functionals.size;
  for (BuildPayload* payload : {&workspace.trial, &workspace.best}) {
    payload->groups.reserve(1U);
    payload->rows.reserve(request.functionals.size);
    payload->donor_global_cells.reserve(request.donors.size);
    payload->donor_local_indices.reserve(request.donors.size);
    payload->weights.reserve(weights);
    payload->clear();
  }
  Status last_failure{StatusCode::invalid_plan, kQuadraticRank};
  for (std::size_t donor_count = limits.minimum_donors;
       donor_count <= request.donors.size; ++donor_count) {
    const Status status = compile_prefix(request, pool, donor_count, limits,
                                         workspace.trial);
    if (!status) {
      last_failure = status;
      continue;
    }
    if (better_candidate(workspace.trial, workspace.best)) {
      workspace.best = workspace.trial;
    }
  }
  if (workspace.best.groups.empty()) {
    return last_failure;
  }
  return append_group(workspace.best, candidate);
}

bool index_in_view(Int3 index, ConstFieldView field) noexcept {
  const auto in_axis = [](std::int32_t value, std::int32_t interior,
                          std::int32_t ghosts) noexcept {
    return interior > 0 && ghosts >= 0 && value >= -ghosts &&
           static_cast<std::int64_t>(value) <
               static_cast<std::int64_t>(interior) + ghosts;
  };
  return in_axis(index.x, field.interior.x, field.ghosts.x) &&
         in_axis(index.y, field.interior.y, field.ghosts.y) &&
         in_axis(index.z, field.interior.z, field.ghosts.z);
}

}  // namespace

Status QuadraticStencilCompiler::compile(
    Span<const QuadraticStencilRequest> requests,
    QuadraticStencilLimits limits, QuadraticStencilPlan& out) noexcept {
  if (requests.data == nullptr || requests.size == 0U ||
      !valid_limits(limits)) {
    return {StatusCode::invalid_plan, kQuadraticInput};
  }
  try {
    std::size_t donor_capacity = 0U;
    std::size_t row_capacity = 0U;
    std::size_t weight_capacity = 0U;
    for (std::size_t index = 0U; index < requests.size; ++index) {
      std::size_t request_weights = 0U;
      if (!checked_add(donor_capacity, requests.data[index].donors.size,
                       donor_capacity) ||
          !checked_add(row_capacity, requests.data[index].functionals.size,
                       row_capacity) ||
          !checked_multiply(requests.data[index].donors.size,
                            requests.data[index].functionals.size,
                            request_weights) ||
          !checked_add(weight_capacity, request_weights, weight_capacity) ||
          donor_capacity > UINT32_MAX || row_capacity > UINT32_MAX ||
          weight_capacity > UINT32_MAX) {
        return {StatusCode::invalid_plan, kQuadraticOverflow};
      }
    }
    BuildPayload payload;
    payload.groups.reserve(requests.size);
    payload.rows.reserve(row_capacity);
    payload.donor_global_cells.reserve(donor_capacity);
    payload.donor_local_indices.reserve(donor_capacity);
    payload.weights.reserve(weight_capacity);
    CompileWorkspace workspace;
    for (std::size_t index = 0U; index < requests.size; ++index) {
      const Status status = compile_group(requests.data[index], limits,
                                          payload, workspace);
      if (!status) {
        return status;
      }
    }
    QuadraticStencilPlan candidate;
    candidate.groups_ = std::move(payload.groups);
    candidate.rows_ = std::move(payload.rows);
    candidate.donor_global_cells_ = std::move(payload.donor_global_cells);
    candidate.donor_local_indices_ = std::move(payload.donor_local_indices);
    candidate.weights_ = std::move(payload.weights);
    candidate.maximum_halo_reach_ = payload.maximum_halo_reach;
    candidate.refresh_fingerprint();
    out = std::move(candidate);
    return {};
  } catch (const std::bad_alloc&) {
    return {StatusCode::allocation_failure, 0U};
  } catch (...) {
    return {StatusCode::invalid_plan, kQuadraticInput};
  }
}

Status evaluate_quadratic_row(const QuadraticStencilPlan& plan,
                              std::uint32_t row_index,
                              ConstFieldView field,
                              std::uint8_t component, double wall_value,
                              double wall_normal_gradient,
                              double& out) noexcept {
  const Span<const QuadraticAffineRow> rows = plan.rows();
  const Span<const QuadraticStencilGroup> groups = plan.groups();
  const Span<const Int3> donor_indices = plan.donor_local_indices();
  const Span<const double> weights = plan.weights();
  if (plan.fingerprint() == 0U || row_index >= rows.size ||
      field.base == nullptr || field.components == 0U ||
      component >= field.components || field.stride_y == 0U ||
      field.stride_z == 0U || field.component_stride == 0U ||
      !std::isfinite(wall_value) ||
      !std::isfinite(wall_normal_gradient)) {
    return {StatusCode::invalid_plan, kQuadraticEvaluation};
  }
  const QuadraticAffineRow& row = rows.data[row_index];
  if (row.group >= groups.size) {
    return {StatusCode::invalid_plan, kQuadraticEvaluation};
  }
  const QuadraticStencilGroup& group = groups.data[row.group];
  const std::size_t donor_begin = group.donor_begin;
  const std::size_t donor_end = donor_begin + group.quality.donor_count;
  const std::size_t weight_end =
      static_cast<std::size_t>(row.weight_begin) + group.quality.donor_count;
  if (donor_end > donor_indices.size || weight_end > weights.size) {
    return {StatusCode::invalid_plan, kQuadraticEvaluation};
  }

  double candidate = row.wall_value_weight * wall_value +
                     row.wall_normal_gradient_weight_m *
                         wall_normal_gradient;
  for (std::size_t donor = 0U; donor < group.quality.donor_count; ++donor) {
    const Int3 index = donor_indices.data[donor_begin + donor];
    if (!index_in_view(index, field)) {
      return {StatusCode::invalid_plan, kQuadraticEvaluation};
    }
    const double value = field.unchecked(index, component);
    const double weight = weights.data[row.weight_begin + donor];
    if (!std::isfinite(value) || !std::isfinite(weight)) {
      return {StatusCode::numerical_failure, kQuadraticEvaluation};
    }
    candidate += weight * value;
  }
  if (!std::isfinite(candidate)) {
    return {StatusCode::numerical_failure, kQuadraticEvaluation};
  }
  out = candidate;
  return {};
}

}  // namespace hundun::v04

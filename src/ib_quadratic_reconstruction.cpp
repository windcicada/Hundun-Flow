// SPDX-License-Identifier: Apache-2.0

#include "hundun/ib_quadratic_reconstruction.hpp"

#include "ib_deterministic_qr_detail.hpp"
#include "ib_quadratic_reconstruction_detail.hpp"

#include "hundun/mesh_geometry.hpp"
#include "hundun/rt_error.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace hundun::immersed {
namespace {

constexpr std::size_t kBasisSize = 10U;
constexpr std::size_t kMinimumDonors = 14U;
constexpr std::size_t kMaximumDonors = 32U;
constexpr std::uint32_t kMaximumHaloReach = 4U;
constexpr double kConditionLimit = 1.0e8;

thread_local const mesh::MeshTopology *validated_topology = nullptr;
thread_local const mesh::MeshGeometry *validated_geometry = nullptr;
thread_local bool relaxed_boundary_authority_coverage = false;

bool finite(runtime::Real3 value) noexcept {
  return std::isfinite(value.x) && std::isfinite(value.y) &&
         std::isfinite(value.z);
}

runtime::Real3 add(runtime::Real3 lhs, runtime::Real3 rhs) noexcept {
  return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

runtime::Real3 subtract(runtime::Real3 lhs, runtime::Real3 rhs) noexcept {
  return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

runtime::Real3 multiply(runtime::Real3 value, double factor) noexcept {
  return {value.x * factor, value.y * factor, value.z * factor};
}

double dot(runtime::Real3 lhs, runtime::Real3 rhs) noexcept {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

runtime::Real3 cross(runtime::Real3 lhs, runtime::Real3 rhs) noexcept {
  return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z,
          lhs.x * rhs.y - lhs.y * rhs.x};
}

double norm(runtime::Real3 value) noexcept {
  return std::hypot(value.x, value.y, value.z);
}

double determinant(runtime::Real3 a, runtime::Real3 b,
                   runtime::Real3 c) noexcept {
  return dot(a, cross(b, c));
}

void hash_u64(std::uint64_t &hash, std::uint64_t value) noexcept {
  constexpr std::uint64_t kPrime = 1099511628211ULL;
  for (int byte = 0; byte < 8; ++byte) {
    hash ^= (value >> static_cast<unsigned>(8 * byte)) & 0xffU;
    hash *= kPrime;
  }
}

std::array<double, kBasisSize>
basis_at(runtime::Real3 point_m, runtime::Real3 origin_m, runtime::Real3 normal,
         runtime::Real3 tangent1, runtime::Real3 tangent2, double scale_m) {
  const auto offset = subtract(point_m, origin_m);
  const double n = dot(offset, normal) / scale_m;
  const double t1 = dot(offset, tangent1) / scale_m;
  const double t2 = dot(offset, tangent2) / scale_m;
  return {1.0, n, t1, t2, n * n, n * t1, n * t2, t1 * t1, t1 * t2, t2 * t2};
}

std::array<runtime::Real3, 8>
cell_vertices(runtime::Int3 cell, const mesh::MeshGeometry &geometry) {
  const std::array<runtime::Int3, 8> offsets{
      runtime::Int3{0, 0, 0}, runtime::Int3{1, 0, 0}, runtime::Int3{1, 1, 0},
      runtime::Int3{0, 1, 0}, runtime::Int3{0, 0, 1}, runtime::Int3{1, 0, 1},
      runtime::Int3{1, 1, 1}, runtime::Int3{0, 1, 1}};
  std::array<runtime::Real3, 8> vertices{};
  for (std::size_t index = 0; index < vertices.size(); ++index) {
    vertices[index] = geometry.vertex_position_m({cell.x + offsets[index].x,
                                                  cell.y + offsets[index].y,
                                                  cell.z + offsets[index].z});
    if (!finite(vertices[index])) {
      throw runtime::Error(
          "quadratic reconstruction donor vertex must be finite");
    }
  }
  return vertices;
}

struct CellMoments final {
  std::array<double, kBasisSize> average{};
  runtime::Real3 centroid_m{};
  double volume_m3{};
};

CellMoments cell_moments(runtime::Int3 cell, runtime::Real3 origin_m,
                         runtime::Real3 normal, runtime::Real3 tangent1,
                         runtime::Real3 tangent2, double scale_m,
                         const mesh::MeshTopology &topology,
                         const mesh::MeshGeometry &geometry) {
  const auto vertices = cell_vertices(cell, geometry);
  runtime::Real3 reference{};
  for (const auto vertex : vertices) {
    reference = add(reference, multiply(vertex, 0.125));
  }

  constexpr std::array<std::array<std::size_t, 3>, 12> triangles{{
      {{0, 4, 7}},
      {{0, 7, 3}},
      {{1, 2, 6}},
      {{1, 6, 5}},
      {{0, 1, 5}},
      {{0, 5, 4}},
      {{3, 7, 6}},
      {{3, 6, 2}},
      {{0, 3, 2}},
      {{0, 2, 1}},
      {{4, 5, 6}},
      {{4, 6, 7}},
  }};

  std::array<double, kBasisSize> integrals{};
  runtime::Real3 centroid_integral{};
  double volume = 0.0;
  const auto reference_basis =
      basis_at(reference, origin_m, normal, tangent1, tangent2, scale_m);
  for (const auto triangle : triangles) {
    const auto a = vertices[triangle[0]];
    const auto b = vertices[triangle[1]];
    const auto c = vertices[triangle[2]];
    const double tetra_volume =
        determinant(subtract(a, reference), subtract(b, reference),
                    subtract(c, reference)) /
        6.0;
    if (!std::isfinite(tetra_volume) || tetra_volume <= 0.0) {
      throw runtime::Error(
          "quadratic reconstruction tetra volume must be finite and positive");
    }
    const std::array<std::array<double, kBasisSize>, 4> samples{{
        reference_basis,
        basis_at(a, origin_m, normal, tangent1, tangent2, scale_m),
        basis_at(b, origin_m, normal, tangent1, tangent2, scale_m),
        basis_at(c, origin_m, normal, tangent1, tangent2, scale_m),
    }};

    std::array<double, 3> sums{};
    std::array<double, 3> square_sums{};
    const auto coordinate = [](const std::array<double, kBasisSize> &basis,
                               std::size_t index) { return basis[index]; };
    for (std::size_t axis = 0; axis < 3; ++axis) {
      for (const auto &sample : samples) {
        const double value = coordinate(sample, axis + 1U);
        sums[axis] += value;
        square_sums[axis] += value * value;
      }
    }
    const auto product_average = [&](std::size_t lhs, std::size_t rhs) {
      double diagonal_sum = 0.0;
      for (const auto &sample : samples) {
        diagonal_sum +=
            coordinate(sample, lhs + 1U) * coordinate(sample, rhs + 1U);
      }
      return (sums[lhs] * sums[rhs] + diagonal_sum) / 20.0;
    };
    const std::array<double, kBasisSize> tetra_average{
        1.0,
        sums[0] / 4.0,
        sums[1] / 4.0,
        sums[2] / 4.0,
        product_average(0, 0),
        product_average(0, 1),
        product_average(0, 2),
        product_average(1, 1),
        product_average(1, 2),
        product_average(2, 2),
    };
    for (std::size_t basis = 0; basis < kBasisSize; ++basis) {
      integrals[basis] += tetra_volume * tetra_average[basis];
    }
    const auto tetra_centroid =
        multiply(add(add(reference, a), add(b, c)), 0.25);
    centroid_integral =
        add(centroid_integral, multiply(tetra_centroid, tetra_volume));
    volume += tetra_volume;
  }
  if (!std::isfinite(volume) || volume <= 0.0) {
    throw runtime::Error(
        "quadratic reconstruction cell volume must be finite and positive");
  }

  const auto global_id = topology.global_cell_id(cell);
  if (const auto local = topology.find_local_cell(global_id);
      local.has_value()) {
    const double expected = geometry.cell_volume_m3(*local);
    const double scale = std::max({1.0, std::abs(expected), std::abs(volume)});
    const double tolerance =
        256.0 * std::numeric_limits<double>::epsilon() * scale;
    if (!std::isfinite(expected) || std::abs(volume - expected) > tolerance) {
      throw runtime::Error(
          "quadratic reconstruction cell volume does not match mesh geometry");
    }
  }

  CellMoments result;
  result.volume_m3 = volume;
  result.centroid_m = multiply(centroid_integral, 1.0 / volume);
  for (std::size_t basis = 0; basis < kBasisSize; ++basis) {
    result.average[basis] = integrals[basis] / volume;
    if (!std::isfinite(result.average[basis])) {
      throw runtime::Error(
          "quadratic reconstruction cell moment must be finite");
    }
  }
  const double moment_tolerance =
      256.0 * std::numeric_limits<double>::epsilon();
  if (std::abs(result.average[0] - 1.0) > moment_tolerance) {
    throw runtime::Error(
        "quadratic reconstruction constant moment is inconsistent");
  }
  return result;
}

void validate_frame(runtime::Real3 origin, runtime::Real3 normal,
                    runtime::Real3 tangent1, runtime::Real3 tangent2,
                    double scale) {
  if (!finite(origin) || !finite(normal) || !finite(tangent1) ||
      !finite(tangent2) || !std::isfinite(scale) || scale <= 0.0) {
    throw runtime::Error(
        "quadratic reconstruction frame must be finite with positive scale");
  }
  const double tolerance = 512.0 * std::numeric_limits<double>::epsilon();
  if (std::abs(norm(normal) - 1.0) > tolerance ||
      std::abs(norm(tangent1) - 1.0) > tolerance ||
      std::abs(norm(tangent2) - 1.0) > tolerance ||
      std::abs(dot(normal, tangent1)) > tolerance ||
      std::abs(dot(normal, tangent2)) > tolerance ||
      std::abs(dot(tangent1, tangent2)) > tolerance ||
      std::abs(dot(cross(normal, tangent1), tangent2) - 1.0) > tolerance) {
    throw runtime::Error(
        "quadratic reconstruction frame must be right-handed orthonormal");
  }
}

std::uint32_t logical_reach(runtime::Int3 anchor, runtime::Int3 donor) {
  const auto magnitude = [](int lhs, int rhs) -> std::uint32_t {
    const auto wide =
        static_cast<std::int64_t>(lhs) - static_cast<std::int64_t>(rhs);
    return static_cast<std::uint32_t>(wide < 0 ? -wide : wide);
  };
  return std::max({magnitude(donor.x, anchor.x), magnitude(donor.y, anchor.y),
                   magnitude(donor.z, anchor.z)});
}

void validate_coverage(const std::vector<CellMoments> &moments,
                       runtime::Real3 origin, runtime::Real3 normal,
                       runtime::Real3 tangent1, runtime::Real3 tangent2,
                       double scale) {
  std::vector<double> normal_coordinates;
  normal_coordinates.reserve(moments.size());
  std::array<bool, 4> quadrants{};
  for (const auto &moment : moments) {
    const auto offset = subtract(moment.centroid_m, origin);
    const double n = dot(offset, normal) / scale;
    const double t1 = dot(offset, tangent1) / scale;
    const double t2 = dot(offset, tangent2) / scale;
    if (!std::isfinite(n) || !std::isfinite(t1) || !std::isfinite(t2) ||
        (!relaxed_boundary_authority_coverage && n <= 0.0)) {
      throw runtime::Error(
          "quadratic reconstruction donors must lie on the positive-normal "
          "fluid side");
    }
    normal_coordinates.push_back(n);
    const std::size_t quadrant = (t1 < 0.0 ? 2U : 0U) + (t2 < 0.0 ? 1U : 0U);
    quadrants[quadrant] = true;
  }
  if (relaxed_boundary_authority_coverage)
    return;
  std::sort(normal_coordinates.begin(), normal_coordinates.end());
  std::size_t bands = 0U;
  double previous = 0.0;
  for (const double coordinate : normal_coordinates) {
    const double tolerance = 512.0 * std::numeric_limits<double>::epsilon() *
                             std::max(1.0, std::abs(coordinate));
    if (bands == 0U || coordinate - previous > tolerance) {
      ++bands;
      previous = coordinate;
    }
  }
  if (bands < 3U) {
    throw runtime::Error(
        "quadratic reconstruction requires three positive-normal donor bands");
  }
  if (!std::all_of(quadrants.begin(), quadrants.end(),
                   [](bool present) { return present; })) {
    throw runtime::Error(
        "quadratic reconstruction requires four tangential donor quadrants");
  }
}

} // namespace

namespace detail {

std::array<double, kQuadraticBasisSize>
quadratic_basis_at(runtime::Real3 point, const QuadraticFrame &frame) {
  validate_frame(frame.origin_m, frame.normal, frame.tangent1, frame.tangent2,
                 frame.scale_m);
  return basis_at(point, frame.origin_m, frame.normal, frame.tangent1,
                  frame.tangent2, frame.scale_m);
}

std::array<double, kQuadraticBasisSize>
quadratic_directional_derivative_basis(runtime::Real3 point,
                                       runtime::Real3 direction,
                                       const QuadraticFrame &frame) {
  validate_frame(frame.origin_m, frame.normal, frame.tangent1, frame.tangent2,
                 frame.scale_m);
  if (!finite(point) || !finite(direction))
    throw runtime::Error(
        "quadratic directional derivative input must be finite");
  const auto value = basis_at(point, frame.origin_m, frame.normal,
                              frame.tangent1, frame.tangent2, frame.scale_m);
  const double dn = dot(direction, frame.normal) / frame.scale_m;
  const double dt1 = dot(direction, frame.tangent1) / frame.scale_m;
  const double dt2 = dot(direction, frame.tangent2) / frame.scale_m;
  return {0.0,
          dn,
          dt1,
          dt2,
          2.0 * value[1U] * dn,
          value[2U] * dn + value[1U] * dt1,
          value[3U] * dn + value[1U] * dt2,
          2.0 * value[2U] * dt1,
          value[3U] * dt1 + value[2U] * dt2,
          2.0 * value[3U] * dt2};
}

std::array<double, kQuadraticBasisSize>
quadratic_cell_average_basis(runtime::Int3 cell, const QuadraticFrame &frame,
                             const mesh::MeshTopology &topology,
                             const mesh::MeshGeometry &geometry) {
  validate_frame(frame.origin_m, frame.normal, frame.tangent1, frame.tangent2,
                 frame.scale_m);
  return cell_moments(cell, frame.origin_m, frame.normal, frame.tangent1,
                      frame.tangent2, frame.scale_m, topology, geometry)
      .average;
}

ValidatedGeometryScope::ValidatedGeometryScope(
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry) {
  if (validated_topology != nullptr || validated_geometry != nullptr) {
    throw runtime::Error(
        "quadratic reconstruction validated geometry scope is nested");
  }
  validated_topology = &topology;
  validated_geometry = &geometry;
}

ValidatedGeometryScope::~ValidatedGeometryScope() noexcept {
  validated_geometry = nullptr;
  validated_topology = nullptr;
}

bool geometry_validation_is_scoped(
    const mesh::MeshTopology &topology,
    const mesh::MeshGeometry &geometry) noexcept {
  return validated_topology == &topology && validated_geometry == &geometry;
}

BoundaryAuthorityCoverageScope::BoundaryAuthorityCoverageScope() {
  if (relaxed_boundary_authority_coverage)
    throw runtime::Error(
        "quadratic reconstruction boundary-authority scope is nested");
  relaxed_boundary_authority_coverage = true;
}

BoundaryAuthorityCoverageScope::~BoundaryAuthorityCoverageScope() noexcept {
  relaxed_boundary_authority_coverage = false;
}

bool boundary_authority_coverage_is_scoped() noexcept {
  return relaxed_boundary_authority_coverage;
}

struct QuadraticReconstructionStorage final {
  runtime::Real3 origin_m{};
  runtime::Real3 normal{};
  runtime::Real3 tangent1{};
  runtime::Real3 tangent2{};
  double scale_m{};
  runtime::Box3 owned_global_box{};
  std::vector<runtime::Int3> donor_global_cells;
  std::vector<mesh::GlobalCellId> donor_global_ids;
  std::vector<double> donor_normal_moments;
  std::vector<double> coefficient_weights;
  std::vector<double> origin_constrained_coefficient_weights;
  std::vector<double> origin_normal_gradient_constrained_coefficient_weights;
  ReconstructionQuality quality{};
  std::uint64_t boundary_authority_link{
      std::numeric_limits<std::uint64_t>::max()};
  int boundary_authority_owner_rank{-1};
  std::shared_ptr<const QuadraticReconstructionStorage>
      boundary_authority_storage;
  std::shared_ptr<const std::vector<BoundaryAuthorityOwner>>
      boundary_authority_catalog;
};

std::vector<double> DeterministicQr::functional_weights(
    const std::vector<double> &functional) const {
  if (functional.size() != columns || rank != columns ||
      thin_q.size() != rows * columns || upper_r.size() != columns * columns ||
      pivots.size() != columns) {
    throw runtime::Error("deterministic QR functional layout is inconsistent");
  }
  std::vector<double> y(columns, 0.0);
  for (std::size_t row = 0; row < columns; ++row) {
    double value = functional[pivots[row]];
    for (std::size_t column = 0; column < row; ++column) {
      value -= upper_r[column * columns + row] * y[column];
    }
    const double diagonal = upper_r[row * columns + row];
    if (!std::isfinite(value) || !std::isfinite(diagonal) || diagonal == 0.0) {
      throw runtime::Error(
          "deterministic QR functional solve encountered breakdown");
    }
    y[row] = value / diagonal;
  }
  std::vector<double> weights(rows, 0.0);
  for (std::size_t row = 0; row < rows; ++row) {
    for (std::size_t column = 0; column < columns; ++column) {
      weights[row] += thin_q[row * columns + column] * y[column];
    }
    if (!std::isfinite(weights[row])) {
      throw runtime::Error("deterministic QR functional weight must be finite");
    }
  }
  return weights;
}

DeterministicQr factorize_design_matrix(const std::vector<double> &matrix,
                                        std::size_t rows, std::size_t columns) {
  if (rows < columns || columns == 0U ||
      rows > std::numeric_limits<std::size_t>::max() / columns ||
      matrix.size() != rows * columns) {
    throw runtime::Error("deterministic QR matrix size is inconsistent");
  }
  if (columns > kBasisSize)
    throw runtime::Error(
        "deterministic QR column count exceeds the quadratic basis");
  for (const double value : matrix) {
    if (!std::isfinite(value)) {
      throw runtime::Error("deterministic QR matrix must be finite");
    }
  }

  std::vector<double> work = matrix;
  std::vector<std::uint32_t> pivots(columns);
  std::iota(pivots.begin(), pivots.end(), 0U);
  std::vector<std::vector<double>> reflectors(columns,
                                              std::vector<double>(rows, 0.0));
  std::vector<double> betas(columns, 0.0);
  std::vector<double> upper(columns * columns, 0.0);
  double maximum_initial_column_norm = 0.0;
  for (std::size_t column = 0; column < columns; ++column) {
    double squared = 0.0;
    for (std::size_t row = 0; row < rows; ++row) {
      squared += work[row * columns + column] * work[row * columns + column];
    }
    maximum_initial_column_norm =
        std::max(maximum_initial_column_norm, std::sqrt(squared));
  }
  const double rank_tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
                                std::max(1.0, maximum_initial_column_norm);

  std::uint32_t rank = 0U;
  for (std::size_t step = 0; step < columns; ++step) {
    std::size_t selected = step;
    double selected_norm = -1.0;
    for (std::size_t column = step; column < columns; ++column) {
      double squared = 0.0;
      for (std::size_t row = step; row < rows; ++row) {
        const double value = work[row * columns + column];
        squared += value * value;
      }
      const double candidate_norm = std::sqrt(squared);
      const double tie_tolerance =
          64.0 * std::numeric_limits<double>::epsilon() *
          std::max(1.0, std::max(selected_norm, candidate_norm));
      if (candidate_norm > selected_norm + tie_tolerance ||
          (std::abs(candidate_norm - selected_norm) <= tie_tolerance &&
           pivots[column] < pivots[selected])) {
        selected = column;
        selected_norm = candidate_norm;
      }
    }
    if (!std::isfinite(selected_norm) || selected_norm <= rank_tolerance) {
      break;
    }
    if (selected != step) {
      for (std::size_t row = 0; row < rows; ++row) {
        std::swap(work[row * columns + step], work[row * columns + selected]);
      }
      std::swap(pivots[step], pivots[selected]);
    }

    auto &reflector = reflectors[step];
    double norm_x = 0.0;
    for (std::size_t row = step; row < rows; ++row) {
      reflector[row] = work[row * columns + step];
      norm_x = std::hypot(norm_x, reflector[row]);
    }
    const double alpha = reflector[step] >= 0.0 ? -norm_x : norm_x;
    reflector[step] -= alpha;
    double squared_v = 0.0;
    for (std::size_t row = step; row < rows; ++row) {
      squared_v += reflector[row] * reflector[row];
    }
    if (!std::isfinite(squared_v) || squared_v == 0.0) {
      break;
    }
    const double beta = 2.0 / squared_v;
    betas[step] = beta;
    for (std::size_t column = step; column < columns; ++column) {
      double projection = 0.0;
      for (std::size_t row = step; row < rows; ++row) {
        projection += reflector[row] * work[row * columns + column];
      }
      projection *= beta;
      for (std::size_t row = step; row < rows; ++row) {
        work[row * columns + column] -= projection * reflector[row];
      }
    }
    work[step * columns + step] = alpha;
    ++rank;
  }
  if (rank != columns) {
    throw runtime::Error(
        "deterministic QR design matrix rank is below its column count");
  }

  double minimum_diagonal = std::numeric_limits<double>::infinity();
  double maximum_diagonal = 0.0;
  for (std::size_t row = 0; row < columns; ++row) {
    for (std::size_t column = row; column < columns; ++column) {
      upper[row * columns + column] = work[row * columns + column];
    }
    const double diagonal = std::abs(upper[row * columns + row]);
    minimum_diagonal = std::min(minimum_diagonal, diagonal);
    maximum_diagonal = std::max(maximum_diagonal, diagonal);
  }
  const double condition = maximum_diagonal / minimum_diagonal;
  if (!std::isfinite(condition) || condition > kConditionLimit) {
    throw runtime::Error(
        "deterministic QR quadratic design matrix exceeds condition limit: " +
        std::to_string(condition));
  }

  std::vector<double> thin_q(rows * columns, 0.0);
  for (std::size_t column = 0; column < columns; ++column) {
    std::vector<double> vector(rows, 0.0);
    vector[column] = 1.0;
    for (std::size_t reverse = columns; reverse-- > 0U;) {
      double projection = 0.0;
      for (std::size_t row = reverse; row < rows; ++row) {
        projection += reflectors[reverse][row] * vector[row];
      }
      projection *= betas[reverse];
      for (std::size_t row = reverse; row < rows; ++row) {
        vector[row] -= projection * reflectors[reverse][row];
      }
    }
    for (std::size_t row = 0; row < rows; ++row) {
      thin_q[row * columns + column] = vector[row];
    }
  }

  std::uint64_t fingerprint = 1469598103934665603ULL;
  for (const auto pivot : pivots) {
    hash_u64(fingerprint, pivot);
  }
  return {rows,
          columns,
          rank,
          condition,
          std::move(pivots),
          fingerprint,
          std::move(thin_q),
          std::move(upper)};
}

namespace {

std::vector<WeightedDonor>
combine_functional_weights(const QuadraticReconstructionStorage &storage,
                           const std::array<double, kBasisSize> &functional) {
  std::vector<WeightedDonor> result(storage.donor_global_ids.size());
  for (std::size_t donor = 0; donor < result.size(); ++donor) {
    double weight = 0.0;
    for (std::size_t coefficient = 0; coefficient < kBasisSize; ++coefficient) {
      weight +=
          functional[coefficient] *
          storage.coefficient_weights[coefficient * result.size() + donor];
    }
    if (!std::isfinite(weight)) {
      throw runtime::Error(
          "quadratic reconstruction functional weight must be finite");
    }
    result[donor] = {storage.donor_global_ids[donor], weight};
  }
  return result;
}

QuadraticReconstructionWeights::AffineBoundaryFunctional
combine_origin_normal_gradient_constrained_functional_weights(
    const QuadraticReconstructionStorage &storage,
    const std::array<double, kBasisSize> &functional) {
  constexpr std::array<std::size_t, kBasisSize - 1U> free_basis{
      0U, 2U, 3U, 4U, 5U, 6U, 7U, 8U, 9U};
  const std::size_t donor_count = storage.donor_global_ids.size();
  if (storage.origin_normal_gradient_constrained_coefficient_weights.size() !=
          free_basis.size() * donor_count ||
      storage.donor_normal_moments.size() != donor_count)
    throw runtime::Error(
        "quadratic reconstruction normal-gradient affine weights are "
        "inconsistent");

  QuadraticReconstructionWeights::AffineBoundaryFunctional result{};
  result.donors.resize(donor_count);
  double constrained_moment = 0.0;
  for (std::size_t donor = 0U; donor < donor_count; ++donor) {
    double weight = 0.0;
    for (std::size_t coefficient = 0U; coefficient < free_basis.size();
         ++coefficient)
      weight +=
          functional[free_basis[coefficient]] *
          storage.origin_normal_gradient_constrained_coefficient_weights
              [coefficient * donor_count + donor];
    if (!std::isfinite(weight))
      throw runtime::Error(
          "quadratic reconstruction normal-gradient affine donor weight is "
          "non-finite");
    result.donors[donor] = {storage.donor_global_ids[donor], weight};
    constrained_moment += weight * storage.donor_normal_moments[donor];
  }
  result.boundary_coefficient =
      storage.scale_m * (functional[1U] - constrained_moment);
  if (!std::isfinite(result.boundary_coefficient))
    throw runtime::Error(
        "quadratic reconstruction normal-gradient affine boundary weight is "
        "non-finite");
  return result;
}

} // namespace

std::vector<WeightedDonor> QuadraticReconstructionWeights::value_weights(
    const QuadraticReconstruction &reconstruction, runtime::Real3 point_m) {
  if (!finite(point_m)) {
    throw runtime::Error(
        "quadratic reconstruction evaluation point must be finite");
  }
  const auto &storage = *reconstruction.storage_;
  return combine_functional_weights(
      storage, basis_at(point_m, storage.origin_m, storage.normal,
                        storage.tangent1, storage.tangent2, storage.scale_m));
}

std::vector<WeightedDonor> QuadraticReconstructionWeights::cell_average_weights(
    const QuadraticReconstruction &reconstruction, runtime::Int3 global_cell,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry) {
  const auto &storage = *reconstruction.storage_;
  const auto moments = cell_moments(
      global_cell, storage.origin_m, storage.normal, storage.tangent1,
      storage.tangent2, storage.scale_m, topology, geometry);
  return combine_functional_weights(storage, moments.average);
}

std::vector<WeightedDonor>
QuadraticReconstructionWeights::directional_gradient_weights(
    const QuadraticReconstruction &reconstruction, runtime::Real3 point_m,
    runtime::Real3 direction) {
  if (!finite(point_m) || !finite(direction)) {
    throw runtime::Error(
        "quadratic reconstruction gradient functional must be finite");
  }
  const auto &storage = *reconstruction.storage_;
  const auto local =
      basis_at(point_m, storage.origin_m, storage.normal, storage.tangent1,
               storage.tangent2, storage.scale_m);
  const double n = local[1];
  const double t1 = local[2];
  const double t2 = local[3];
  const double derivative_n = dot(direction, storage.normal) / storage.scale_m;
  const double derivative_t1 =
      dot(direction, storage.tangent1) / storage.scale_m;
  const double derivative_t2 =
      dot(direction, storage.tangent2) / storage.scale_m;
  const std::array<double, kBasisSize> functional{
      0.0,
      derivative_n,
      derivative_t1,
      derivative_t2,
      2.0 * n * derivative_n,
      t1 * derivative_n + n * derivative_t1,
      t2 * derivative_n + n * derivative_t2,
      2.0 * t1 * derivative_t1,
      t2 * derivative_t1 + t1 * derivative_t2,
      2.0 * t2 * derivative_t2,
  };
  return combine_functional_weights(storage, functional);
}

std::vector<WeightedDonor>
QuadraticReconstructionWeights::origin_constrained_directional_gradient_weights(
    const QuadraticReconstruction &reconstruction, runtime::Real3 point_m,
    runtime::Real3 direction) {
  if (!finite(point_m) || !finite(direction)) {
    throw runtime::Error(
        "quadratic reconstruction constrained gradient functional must be "
        "finite");
  }
  const auto &storage = *reconstruction.storage_;
  const auto local =
      basis_at(point_m, storage.origin_m, storage.normal, storage.tangent1,
               storage.tangent2, storage.scale_m);
  const double n = local[1];
  const double t1 = local[2];
  const double t2 = local[3];
  const double derivative_n = dot(direction, storage.normal) / storage.scale_m;
  const double derivative_t1 =
      dot(direction, storage.tangent1) / storage.scale_m;
  const double derivative_t2 =
      dot(direction, storage.tangent2) / storage.scale_m;
  const std::array<double, kBasisSize> functional{
      0.0,
      derivative_n,
      derivative_t1,
      derivative_t2,
      2.0 * n * derivative_n,
      t1 * derivative_n + n * derivative_t1,
      t2 * derivative_n + n * derivative_t2,
      2.0 * t1 * derivative_t1,
      t2 * derivative_t1 + t1 * derivative_t2,
      2.0 * t2 * derivative_t2,
  };
  const std::size_t donor_count = storage.donor_global_ids.size();
  if (storage.origin_constrained_coefficient_weights.size() !=
      (kBasisSize - 1U) * donor_count) {
    throw runtime::Error(
        "quadratic reconstruction constrained functional weights are "
        "inconsistent");
  }
  std::vector<WeightedDonor> result(donor_count);
  for (std::size_t donor = 0U; donor < donor_count; ++donor) {
    double weight = 0.0;
    for (std::size_t coefficient = 1U; coefficient < kBasisSize;
         ++coefficient) {
      weight +=
          functional[coefficient] *
          storage.origin_constrained_coefficient_weights[(coefficient - 1U) *
                                                             donor_count +
                                                         donor];
    }
    if (!std::isfinite(weight)) {
      throw runtime::Error(
          "quadratic reconstruction constrained functional weight must be "
          "finite");
    }
    result[donor] = {storage.donor_global_ids[donor], weight};
  }
  return result;
}

QuadraticReconstructionWeights::AffineBoundaryFunctional
QuadraticReconstructionWeights::
    origin_normal_gradient_constrained_value_weights(
        const QuadraticReconstruction &reconstruction,
        runtime::Real3 point_m) {
  if (!finite(point_m))
    throw runtime::Error(
        "quadratic reconstruction evaluation point must be finite");
  const auto &storage = *reconstruction.storage_;
  return combine_origin_normal_gradient_constrained_functional_weights(
      storage, basis_at(point_m, storage.origin_m, storage.normal,
                        storage.tangent1, storage.tangent2, storage.scale_m));
}

QuadraticReconstructionWeights::AffineBoundaryFunctional
QuadraticReconstructionWeights::
    origin_normal_gradient_constrained_directional_gradient_weights(
        const QuadraticReconstruction &reconstruction,
        runtime::Real3 point_m, runtime::Real3 direction) {
  if (!finite(point_m) || !finite(direction))
    throw runtime::Error(
        "quadratic reconstruction constrained gradient functional must be "
        "finite");
  const auto &storage = *reconstruction.storage_;
  const auto local =
      basis_at(point_m, storage.origin_m, storage.normal, storage.tangent1,
               storage.tangent2, storage.scale_m);
  const double n = local[1];
  const double t1 = local[2];
  const double t2 = local[3];
  const double derivative_n = dot(direction, storage.normal) / storage.scale_m;
  const double derivative_t1 =
      dot(direction, storage.tangent1) / storage.scale_m;
  const double derivative_t2 =
      dot(direction, storage.tangent2) / storage.scale_m;
  const std::array<double, kBasisSize> functional{
      0.0,
      derivative_n,
      derivative_t1,
      derivative_t2,
      2.0 * n * derivative_n,
      t1 * derivative_n + n * derivative_t1,
      t2 * derivative_n + n * derivative_t2,
      2.0 * t1 * derivative_t1,
      t2 * derivative_t1 + t1 * derivative_t2,
      2.0 * t2 * derivative_t2,
  };
  return combine_origin_normal_gradient_constrained_functional_weights(
      storage, functional);
}

const std::vector<mesh::GlobalCellId> &
QuadraticReconstructionWeights::donor_global_ids(
    const QuadraticReconstruction &reconstruction) noexcept {
  return reconstruction.storage_->donor_global_ids;
}

} // namespace detail

QuadraticReconstruction QuadraticReconstruction::create(
    runtime::Real3 origin_m, runtime::Real3 normal, runtime::Real3 tangent1,
    runtime::Real3 tangent2, double scale_m,
    runtime::Int3 stencil_anchor_global_cell,
    const std::vector<runtime::Int3> &donor_global_cells,
    const mesh::MeshTopology &topology, const mesh::MeshGeometry &geometry) {
  if (!detail::geometry_validation_is_scoped(topology, geometry))
    geometry.require_compatible(topology);
  validate_frame(origin_m, normal, tangent1, tangent2, scale_m);
  if (donor_global_cells.size() < kMinimumDonors ||
      (donor_global_cells.size() > kMaximumDonors &&
       !detail::boundary_authority_coverage_is_scoped())) {
    throw runtime::Error(
        "quadratic reconstruction requires between fourteen and thirty-two "
        "donors");
  }
  static_cast<void>(topology.global_cell_id(stencil_anchor_global_cell));

  std::vector<std::pair<mesh::GlobalCellId, runtime::Int3>> ordered;
  ordered.reserve(donor_global_cells.size());
  for (const auto donor : donor_global_cells) {
    ordered.emplace_back(topology.global_cell_id(donor), donor);
  }
  std::sort(
      ordered.begin(), ordered.end(),
      [](const auto &lhs, const auto &rhs) { return lhs.first < rhs.first; });
  for (std::size_t index = 1; index < ordered.size(); ++index) {
    if (ordered[index - 1U].first == ordered[index].first) {
      throw runtime::Error(
          "quadratic reconstruction donor global IDs must be unique");
    }
  }

  std::vector<CellMoments> moments;
  moments.reserve(ordered.size());
  std::uint32_t maximum_reach = 0U;
  std::vector<runtime::Int3> canonical_donors;
  canonical_donors.reserve(ordered.size());
  std::vector<mesh::GlobalCellId> canonical_global_ids;
  canonical_global_ids.reserve(ordered.size());
  for (const auto &[global_id, donor] : ordered) {
    maximum_reach = std::max(maximum_reach,
                             logical_reach(stencil_anchor_global_cell, donor));
    canonical_donors.push_back(donor);
    canonical_global_ids.push_back(global_id);
    moments.push_back(cell_moments(donor, origin_m, normal, tangent1, tangent2,
                                   scale_m, topology, geometry));
  }
  if (maximum_reach > kMaximumHaloReach) {
    throw runtime::Error(
        "quadratic reconstruction donor exceeds maximum halo reach four");
  }
  validate_coverage(moments, origin_m, normal, tangent1, tangent2, scale_m);

  std::vector<double> design(ordered.size() * kBasisSize, 0.0);
  for (std::size_t row = 0; row < ordered.size(); ++row) {
    for (std::size_t basis = 0; basis < kBasisSize; ++basis) {
      design[row * kBasisSize + basis] = moments[row].average[basis];
    }
  }
  const auto qr =
      detail::factorize_design_matrix(design, ordered.size(), kBasisSize);
  std::vector<double> coefficient_weights(kBasisSize * ordered.size(), 0.0);
  for (std::size_t coefficient = 0; coefficient < kBasisSize; ++coefficient) {
    std::vector<double> functional(kBasisSize, 0.0);
    functional[coefficient] = 1.0;
    const auto weights = qr.functional_weights(functional);
    for (std::size_t donor = 0; donor < ordered.size(); ++donor) {
      coefficient_weights[coefficient * ordered.size() + donor] =
          weights[donor];
    }
  }
  constexpr std::size_t constrained_basis_size = kBasisSize - 1U;
  std::vector<double> constrained_design(
      ordered.size() * constrained_basis_size, 0.0);
  for (std::size_t row = 0; row < ordered.size(); ++row)
    for (std::size_t basis = 0; basis < constrained_basis_size; ++basis)
      constrained_design[row * constrained_basis_size + basis] =
          moments[row].average[basis + 1U];
  const auto constrained_qr = detail::factorize_design_matrix(
      constrained_design, ordered.size(), constrained_basis_size);
  std::vector<double> origin_constrained_coefficient_weights(
      constrained_basis_size * ordered.size(), 0.0);
  for (std::size_t coefficient = 0; coefficient < constrained_basis_size;
       ++coefficient) {
    std::vector<double> functional(constrained_basis_size, 0.0);
    functional[coefficient] = 1.0;
    const auto weights = constrained_qr.functional_weights(functional);
    for (std::size_t donor = 0; donor < ordered.size(); ++donor)
      origin_constrained_coefficient_weights[coefficient * ordered.size() +
                                             donor] = weights[donor];
  }
  constexpr std::array<std::size_t, kBasisSize - 1U> normal_gradient_free_basis{
      0U, 2U, 3U, 4U, 5U, 6U, 7U, 8U, 9U};
  std::vector<double> normal_gradient_constrained_design(
      ordered.size() * normal_gradient_free_basis.size(), 0.0);
  for (std::size_t row = 0; row < ordered.size(); ++row)
    for (std::size_t basis = 0; basis < normal_gradient_free_basis.size();
         ++basis)
      normal_gradient_constrained_design[row *
                                             normal_gradient_free_basis.size() +
                                         basis] =
          moments[row].average[normal_gradient_free_basis[basis]];
  const auto normal_gradient_constrained_qr = detail::factorize_design_matrix(
      normal_gradient_constrained_design, ordered.size(),
      normal_gradient_free_basis.size());
  std::vector<double> origin_normal_gradient_constrained_coefficient_weights(
      normal_gradient_free_basis.size() * ordered.size(), 0.0);
  for (std::size_t coefficient = 0;
       coefficient < normal_gradient_free_basis.size(); ++coefficient) {
    std::vector<double> functional(normal_gradient_free_basis.size(), 0.0);
    functional[coefficient] = 1.0;
    const auto weights =
        normal_gradient_constrained_qr.functional_weights(functional);
    for (std::size_t donor = 0; donor < ordered.size(); ++donor)
      origin_normal_gradient_constrained_coefficient_weights
          [coefficient * ordered.size() + donor] = weights[donor];
  }

  auto storage = std::make_shared<detail::QuadraticReconstructionStorage>();
  storage->origin_m = origin_m;
  storage->normal = normal;
  storage->tangent1 = tangent1;
  storage->tangent2 = tangent2;
  storage->scale_m = scale_m;
  storage->owned_global_box = topology.owned_global_box();
  storage->donor_global_cells = std::move(canonical_donors);
  storage->donor_global_ids = std::move(canonical_global_ids);
  storage->donor_normal_moments.reserve(moments.size());
  for (const auto &moment : moments)
    storage->donor_normal_moments.push_back(moment.average[1]);
  storage->coefficient_weights = std::move(coefficient_weights);
  storage->origin_constrained_coefficient_weights =
      std::move(origin_constrained_coefficient_weights);
  storage->origin_normal_gradient_constrained_coefficient_weights =
      std::move(origin_normal_gradient_constrained_coefficient_weights);
  storage->quality = {qr.rank, qr.condition_estimate, maximum_reach,
                      qr.pivot_fingerprint};
  return QuadraticReconstruction(std::move(storage));
}

namespace {

std::array<double, kBasisSize>
coefficients(const detail::QuadraticReconstructionStorage &storage,
             const runtime::FieldView<const double> &field,
             std::size_t component) {
  if (component >= field.components()) {
    throw runtime::Error(
        "quadratic reconstruction field component is out of bounds");
  }
  const auto extent = field.interior_extent();
  const auto box = storage.owned_global_box;
  if (extent.x != box.end.x - box.begin.x ||
      extent.y != box.end.y - box.begin.y ||
      extent.z != box.end.z - box.begin.z) {
    throw runtime::Error(
        "quadratic reconstruction field layout does not match owned box");
  }
  std::array<double, kBasisSize> result{};
  const auto component_index = static_cast<int>(component);
  for (std::size_t donor = 0; donor < storage.donor_global_cells.size();
       ++donor) {
    const auto global = storage.donor_global_cells[donor];
    const int i = global.x - box.begin.x;
    const int j = global.y - box.begin.y;
    const int k = global.z - box.begin.z;
    const double value = field(i, j, k, component_index);
    if (!std::isfinite(value)) {
      throw runtime::Error(
          "quadratic reconstruction donor field value must be finite");
    }
    for (std::size_t coefficient = 0; coefficient < kBasisSize; ++coefficient) {
      result[coefficient] +=
          storage.coefficient_weights[coefficient *
                                          storage.donor_global_cells.size() +
                                      donor] *
          value;
    }
  }
  if (!std::all_of(result.begin(), result.end(),
                   [](double value) { return std::isfinite(value); })) {
    throw runtime::Error(
        "quadratic reconstruction coefficients must be finite");
  }
  return result;
}

std::array<double, kBasisSize> origin_constrained_coefficients(
    const detail::QuadraticReconstructionStorage &storage,
    const runtime::FieldView<const double> &field, std::size_t component,
    double value_at_origin) {
  if (component >= field.components() || !std::isfinite(value_at_origin))
    throw runtime::Error(
        "quadratic reconstruction origin constraint is invalid");
  const auto extent = field.interior_extent();
  const auto box = storage.owned_global_box;
  if (extent.x != box.end.x - box.begin.x ||
      extent.y != box.end.y - box.begin.y ||
      extent.z != box.end.z - box.begin.z)
    throw runtime::Error(
        "quadratic reconstruction field layout does not match owned box");
  const std::size_t donor_count = storage.donor_global_cells.size();
  if (storage.origin_constrained_coefficient_weights.size() !=
      (kBasisSize - 1U) * donor_count)
    throw runtime::Error(
        "quadratic reconstruction constrained weights are inconsistent");
  std::array<double, kBasisSize> result{};
  result[0] = value_at_origin;
  const auto component_index = static_cast<int>(component);
  for (std::size_t donor = 0; donor < donor_count; ++donor) {
    const auto global = storage.donor_global_cells[donor];
    const int i = global.x - box.begin.x;
    const int j = global.y - box.begin.y;
    const int k = global.z - box.begin.z;
    const double donor_value = field(i, j, k, component_index);
    if (!std::isfinite(donor_value))
      throw runtime::Error(
          "quadratic reconstruction donor field value must be finite");
    const double shifted = donor_value - value_at_origin;
    for (std::size_t coefficient = 1U; coefficient < kBasisSize; ++coefficient)
      result[coefficient] +=
          storage.origin_constrained_coefficient_weights[(coefficient - 1U) *
                                                             donor_count +
                                                         donor] *
          shifted;
  }
  if (!std::all_of(result.begin(), result.end(),
                   [](double value) { return std::isfinite(value); }))
    throw runtime::Error(
        "quadratic reconstruction constrained coefficients must be finite");
  return result;
}

std::array<double, kBasisSize> origin_normal_gradient_constrained_coefficients(
    const detail::QuadraticReconstructionStorage &storage,
    const runtime::FieldView<const double> &field, std::size_t component,
    double normal_gradient_at_origin) {
  if (component >= field.components() ||
      !std::isfinite(normal_gradient_at_origin))
    throw runtime::Error(
        "quadratic reconstruction origin normal-gradient constraint is "
        "invalid");
  const auto extent = field.interior_extent();
  const auto box = storage.owned_global_box;
  if (extent.x != box.end.x - box.begin.x ||
      extent.y != box.end.y - box.begin.y ||
      extent.z != box.end.z - box.begin.z)
    throw runtime::Error(
        "quadratic reconstruction field layout does not match owned box");
  constexpr std::array<std::size_t, kBasisSize - 1U> free_basis{
      0U, 2U, 3U, 4U, 5U, 6U, 7U, 8U, 9U};
  const std::size_t donor_count = storage.donor_global_cells.size();
  if (storage.origin_normal_gradient_constrained_coefficient_weights.size() !=
          free_basis.size() * donor_count ||
      storage.donor_normal_moments.size() != donor_count)
    throw runtime::Error(
        "quadratic reconstruction normal-gradient constrained weights are "
        "inconsistent");
  std::array<double, kBasisSize> result{};
  result[1] = storage.scale_m * normal_gradient_at_origin;
  const auto component_index = static_cast<int>(component);
  for (std::size_t donor = 0; donor < donor_count; ++donor) {
    const auto global = storage.donor_global_cells[donor];
    const int i = global.x - box.begin.x;
    const int j = global.y - box.begin.y;
    const int k = global.z - box.begin.z;
    const double donor_value = field(i, j, k, component_index);
    if (!std::isfinite(donor_value))
      throw runtime::Error(
          "quadratic reconstruction donor field value must be finite");
    const double shifted =
        donor_value - result[1] * storage.donor_normal_moments[donor];
    for (std::size_t coefficient = 0; coefficient < free_basis.size();
         ++coefficient)
      result[free_basis[coefficient]] +=
          storage.origin_normal_gradient_constrained_coefficient_weights
              [coefficient * donor_count + donor] *
          shifted;
  }
  if (!std::all_of(result.begin(), result.end(),
                   [](double value) { return std::isfinite(value); }))
    throw runtime::Error(
        "quadratic reconstruction normal-gradient constrained coefficients "
        "must be finite");
  return result;
}

double evaluate_value(const detail::QuadraticReconstructionStorage &storage,
                      runtime::Real3 point_m,
                      const std::array<double, kBasisSize> &recovered) {
  const auto basis =
      basis_at(point_m, storage.origin_m, storage.normal, storage.tangent1,
               storage.tangent2, storage.scale_m);
  double result = 0.0;
  for (std::size_t index = 0; index < kBasisSize; ++index)
    result += recovered[index] * basis[index];
  if (!std::isfinite(result))
    throw runtime::Error("quadratic reconstruction value must be finite");
  return result;
}

runtime::Real3
evaluate_gradient(const detail::QuadraticReconstructionStorage &storage,
                  runtime::Real3 point_m,
                  const std::array<double, kBasisSize> &recovered) {
  const auto basis =
      basis_at(point_m, storage.origin_m, storage.normal, storage.tangent1,
               storage.tangent2, storage.scale_m);
  const double n = basis[1];
  const double t1 = basis[2];
  const double t2 = basis[3];
  const double inverse_scale = 1.0 / storage.scale_m;
  const double derivative_n = (recovered[1] + 2.0 * recovered[4] * n +
                               recovered[5] * t1 + recovered[6] * t2) *
                              inverse_scale;
  const double derivative_t1 = (recovered[2] + recovered[5] * n +
                                2.0 * recovered[7] * t1 + recovered[8] * t2) *
                               inverse_scale;
  const double derivative_t2 = (recovered[3] + recovered[6] * n +
                                recovered[8] * t1 + 2.0 * recovered[9] * t2) *
                               inverse_scale;
  const auto result = add(add(multiply(storage.normal, derivative_n),
                              multiply(storage.tangent1, derivative_t1)),
                          multiply(storage.tangent2, derivative_t2));
  if (!finite(result))
    throw runtime::Error("quadratic reconstruction gradient must be finite");
  return result;
}

} // namespace

namespace detail {

double QuadraticReconstructionWeights::value_with_origin_constraint(
    const QuadraticReconstruction &reconstruction, runtime::Real3 point_m,
    const runtime::FieldView<const double> &field, std::size_t component,
    double value_at_origin) {
  if (!finite(point_m))
    throw runtime::Error(
        "quadratic reconstruction evaluation point must be finite");
  const auto recovered = origin_constrained_coefficients(
      *reconstruction.storage_, field, component, value_at_origin);
  return evaluate_value(*reconstruction.storage_, point_m, recovered);
}

runtime::Real3 QuadraticReconstructionWeights::gradient_with_origin_constraint(
    const QuadraticReconstruction &reconstruction, runtime::Real3 point_m,
    const runtime::FieldView<const double> &field, std::size_t component,
    double value_at_origin) {
  if (!finite(point_m))
    throw runtime::Error(
        "quadratic reconstruction evaluation point must be finite");
  const auto recovered = origin_constrained_coefficients(
      *reconstruction.storage_, field, component, value_at_origin);
  return evaluate_gradient(*reconstruction.storage_, point_m, recovered);
}

double QuadraticReconstructionWeights::value_with_origin_normal_gradient(
    const QuadraticReconstruction &reconstruction, runtime::Real3 point_m,
    const runtime::FieldView<const double> &field, std::size_t component,
    double normal_gradient_at_origin) {
  if (!finite(point_m))
    throw runtime::Error(
        "quadratic reconstruction evaluation point must be finite");
  const auto recovered = origin_normal_gradient_constrained_coefficients(
      *reconstruction.storage_, field, component, normal_gradient_at_origin);
  return evaluate_value(*reconstruction.storage_, point_m, recovered);
}

runtime::Real3
QuadraticReconstructionWeights::gradient_with_origin_normal_gradient(
    const QuadraticReconstruction &reconstruction, runtime::Real3 point_m,
    const runtime::FieldView<const double> &field, std::size_t component,
    double normal_gradient_at_origin) {
  if (!finite(point_m))
    throw runtime::Error(
        "quadratic reconstruction evaluation point must be finite");
  const auto recovered = origin_normal_gradient_constrained_coefficients(
      *reconstruction.storage_, field, component, normal_gradient_at_origin);
  return evaluate_gradient(*reconstruction.storage_, point_m, recovered);
}

QuadraticReconstruction QuadraticReconstructionWeights::with_boundary_authority(
    const QuadraticReconstruction &point_reconstruction, std::uint64_t link,
    int owner_rank, const QuadraticReconstruction &authority_reconstruction,
    std::shared_ptr<const std::vector<BoundaryAuthorityOwner>> catalog) {
  if (link == std::numeric_limits<std::uint64_t>::max() || owner_rank < 0 ||
      !catalog || catalog->empty())
    throw runtime::Error(
        "quadratic reconstruction boundary authority is invalid");
  auto storage = std::make_shared<QuadraticReconstructionStorage>(
      *point_reconstruction.storage_);
  storage->boundary_authority_link = link;
  storage->boundary_authority_owner_rank = owner_rank;
  storage->boundary_authority_storage = authority_reconstruction.storage_;
  storage->boundary_authority_catalog = std::move(catalog);
  return QuadraticReconstruction(std::move(storage));
}

std::uint64_t QuadraticReconstructionWeights::boundary_authority_link(
    const QuadraticReconstruction &reconstruction) {
  if (reconstruction.storage_->boundary_authority_link ==
          std::numeric_limits<std::uint64_t>::max() ||
      reconstruction.storage_->boundary_authority_owner_rank < 0 ||
      !reconstruction.storage_->boundary_authority_storage)
    throw runtime::Error(
        "quadratic reconstruction boundary authority is unavailable");
  return reconstruction.storage_->boundary_authority_link;
}

int QuadraticReconstructionWeights::boundary_authority_owner_rank(
    const QuadraticReconstruction &reconstruction) {
  static_cast<void>(boundary_authority_link(reconstruction));
  return reconstruction.storage_->boundary_authority_owner_rank;
}

QuadraticReconstruction
QuadraticReconstructionWeights::boundary_authority_reconstruction(
    const QuadraticReconstruction &reconstruction) {
  static_cast<void>(boundary_authority_link(reconstruction));
  return QuadraticReconstruction(
      reconstruction.storage_->boundary_authority_storage);
}

const std::vector<BoundaryAuthorityOwner> &
QuadraticReconstructionWeights::boundary_authority_catalog(
    const QuadraticReconstruction &reconstruction) {
  static_cast<void>(boundary_authority_link(reconstruction));
  if (!reconstruction.storage_->boundary_authority_catalog)
    throw runtime::Error(
        "quadratic reconstruction boundary authority catalog is unavailable");
  return *reconstruction.storage_->boundary_authority_catalog;
}

double
value_with_origin_constraint(const QuadraticReconstruction &reconstruction,
                             runtime::Real3 point_m,
                             const runtime::FieldView<const double> &field,
                             std::size_t component, double value_at_origin) {
  return QuadraticReconstructionWeights::value_with_origin_constraint(
      reconstruction, point_m, field, component, value_at_origin);
}

runtime::Real3
gradient_with_origin_constraint(const QuadraticReconstruction &reconstruction,
                                runtime::Real3 point_m,
                                const runtime::FieldView<const double> &field,
                                std::size_t component, double value_at_origin) {
  return QuadraticReconstructionWeights::gradient_with_origin_constraint(
      reconstruction, point_m, field, component, value_at_origin);
}

double value_with_origin_normal_gradient(
    const QuadraticReconstruction &reconstruction, runtime::Real3 point_m,
    const runtime::FieldView<const double> &field, std::size_t component,
    double normal_gradient_at_origin) {
  return QuadraticReconstructionWeights::value_with_origin_normal_gradient(
      reconstruction, point_m, field, component, normal_gradient_at_origin);
}

runtime::Real3 gradient_with_origin_normal_gradient(
    const QuadraticReconstruction &reconstruction, runtime::Real3 point_m,
    const runtime::FieldView<const double> &field, std::size_t component,
    double normal_gradient_at_origin) {
  return QuadraticReconstructionWeights::gradient_with_origin_normal_gradient(
      reconstruction, point_m, field, component, normal_gradient_at_origin);
}

QuadraticReconstruction with_boundary_authority(
    const QuadraticReconstruction &point_reconstruction, std::uint64_t link,
    int owner_rank, const QuadraticReconstruction &authority_reconstruction,
    std::shared_ptr<const std::vector<BoundaryAuthorityOwner>> catalog) {
  return QuadraticReconstructionWeights::with_boundary_authority(
      point_reconstruction, link, owner_rank, authority_reconstruction,
      std::move(catalog));
}

std::uint64_t
boundary_authority_link(const QuadraticReconstruction &reconstruction) {
  return QuadraticReconstructionWeights::boundary_authority_link(
      reconstruction);
}

int boundary_authority_owner_rank(
    const QuadraticReconstruction &reconstruction) {
  return QuadraticReconstructionWeights::boundary_authority_owner_rank(
      reconstruction);
}

QuadraticReconstruction boundary_authority_reconstruction(
    const QuadraticReconstruction &reconstruction) {
  return QuadraticReconstructionWeights::boundary_authority_reconstruction(
      reconstruction);
}

const std::vector<BoundaryAuthorityOwner> &
boundary_authority_catalog(const QuadraticReconstruction &reconstruction) {
  return QuadraticReconstructionWeights::boundary_authority_catalog(
      reconstruction);
}

} // namespace detail

double
QuadraticReconstruction::value(runtime::Real3 point_m,
                               const runtime::FieldView<const double> &field,
                               std::size_t component) const {
  if (!finite(point_m)) {
    throw runtime::Error(
        "quadratic reconstruction evaluation point must be finite");
  }
  const auto recovered = coefficients(*storage_, field, component);
  return evaluate_value(*storage_, point_m, recovered);
}

runtime::Real3
QuadraticReconstruction::gradient(runtime::Real3 point_m,
                                  const runtime::FieldView<const double> &field,
                                  std::size_t component) const {
  if (!finite(point_m)) {
    throw runtime::Error(
        "quadratic reconstruction evaluation point must be finite");
  }
  const auto recovered = coefficients(*storage_, field, component);
  return evaluate_gradient(*storage_, point_m, recovered);
}

const ReconstructionQuality &QuadraticReconstruction::quality() const noexcept {
  return storage_->quality;
}

} // namespace hundun::immersed

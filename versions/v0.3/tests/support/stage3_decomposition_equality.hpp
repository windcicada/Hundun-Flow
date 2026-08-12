// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/ib_wall_force.hpp"
#include "hundun/rt_error.hpp"
#include "hundun/rt_types.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace hundun::test::stage3 {

inline constexpr double kDecompositionRelativeTolerance = 5.0e-12;

struct DecompositionFieldComparison final {
  bool equal{};
  double tolerance{};
  double maximum_difference{};
  std::size_t maximum_index{};
};

inline DecompositionFieldComparison
compare_decomposition_field(const std::vector<double> &reference,
                            const std::vector<double> &candidate) noexcept {
  DecompositionFieldComparison result{};
  if (reference.size() != candidate.size())
    return result;
  double scale = 1.0;
  for (const double value : reference) {
    if (!std::isfinite(value))
      return result;
    scale = std::max(scale, std::abs(value));
  }
  result.tolerance = kDecompositionRelativeTolerance * scale;
  for (std::size_t index = 0U; index < reference.size(); ++index) {
    if (!std::isfinite(candidate[index]))
      return result;
    const double difference = std::abs(reference[index] - candidate[index]);
    if (difference > result.maximum_difference) {
      result.maximum_difference = difference;
      result.maximum_index = index;
    }
  }
  result.equal = result.maximum_difference <= result.tolerance;
  return result;
}

inline bool decomposition_vector_equal(runtime::Real3 reference,
                                       runtime::Real3 candidate) noexcept {
  const double values[] = {reference.x, reference.y, reference.z,
                           candidate.x, candidate.y, candidate.z};
  for (const double value : values)
    if (!std::isfinite(value))
      return false;
  const double scale = std::max({1.0, std::abs(reference.x),
                                 std::abs(reference.y), std::abs(reference.z)});
  const double tolerance = kDecompositionRelativeTolerance * scale;
  return std::abs(reference.x - candidate.x) <= tolerance &&
         std::abs(reference.y - candidate.y) <= tolerance &&
         std::abs(reference.z - candidate.z) <= tolerance;
}

inline bool
decomposition_force_equal(const immersed::ForceComponents &reference,
                          const immersed::ForceComponents &candidate) noexcept {
  return decomposition_vector_equal(reference.pressure_N,
                                    candidate.pressure_N) &&
         decomposition_vector_equal(reference.viscous_N, candidate.viscous_N) &&
         decomposition_vector_equal(reference.total_N, candidate.total_N);
}

inline void require_decomposition_field(const char *name,
                                        const std::vector<double> &reference,
                                        const std::vector<double> &candidate) {
  const auto result = compare_decomposition_field(reference, candidate);
  if (result.equal)
    return;
  std::ostringstream message;
  message.precision(17);
  message << name << " decomposition comparison failed";
  if (reference.size() != candidate.size()) {
    message << ": size " << candidate.size() << " differs from reference "
            << reference.size();
  } else {
    message << ": maximum difference " << result.maximum_difference
            << " exceeds " << result.tolerance << " at "
            << result.maximum_index;
    if (!reference.empty())
      message << " (reference " << reference[result.maximum_index]
              << ", candidate " << candidate[result.maximum_index] << ')';
  }
  throw runtime::Error(message.str());
}

inline void
require_decomposition_force(const char *name,
                            const immersed::ForceComponents &reference,
                            const immersed::ForceComponents &candidate) {
  if (decomposition_force_equal(reference, candidate))
    return;
  throw runtime::Error(std::string(name) +
                       " decomposition force comparison failed");
}

inline bool decomposition_equality_oracle_is_mutation_sensitive() {
  const std::vector<double> field{1.0, -2.0, 3.0};
  auto changed_field = field;
  changed_field[1] += 1.0e-8;
  auto changed_size = field;
  changed_size.push_back(4.0);
  auto non_finite_field = field;
  non_finite_field[0] = std::numeric_limits<double>::quiet_NaN();
  const runtime::Real3 vector{1.0, -2.0, 3.0};
  auto changed_vector = vector;
  changed_vector.y += 1.0e-8;
  auto non_finite_vector = vector;
  non_finite_vector.z = std::numeric_limits<double>::infinity();
  const immersed::ForceComponents force{
      {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {5.0, 7.0, 9.0}};
  auto changed_force = force;
  changed_force.viscous_N.z += 1.0e-8;
  return compare_decomposition_field(field, field).equal &&
         !compare_decomposition_field(field, changed_field).equal &&
         !compare_decomposition_field(field, changed_size).equal &&
         !compare_decomposition_field(field, non_finite_field).equal &&
         decomposition_vector_equal(vector, vector) &&
         !decomposition_vector_equal(vector, changed_vector) &&
         !decomposition_vector_equal(vector, non_finite_vector) &&
         decomposition_force_equal(force, force) &&
         !decomposition_force_equal(force, changed_force);
}

} // namespace hundun::test::stage3

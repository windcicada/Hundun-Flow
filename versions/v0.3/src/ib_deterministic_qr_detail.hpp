// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/ib_quadratic_reconstruction.hpp"
#include "ib_quadratic_reconstruction_detail.hpp"

#include <cstddef>
#include <cstdint>
#include <array>
#include <vector>

namespace hundun::immersed::detail {

inline constexpr std::size_t kQuadraticBasisSize = 10U;

struct QuadraticFrame final {
  runtime::Real3 origin_m{};
  runtime::Real3 normal{};
  runtime::Real3 tangent1{};
  runtime::Real3 tangent2{};
  double scale_m{};
};

std::array<double, kQuadraticBasisSize>
quadratic_basis_at(runtime::Real3, const QuadraticFrame &);
std::array<double, kQuadraticBasisSize>
quadratic_directional_derivative_basis(runtime::Real3, runtime::Real3,
                                       const QuadraticFrame &);
std::array<double, kQuadraticBasisSize>
quadratic_cell_average_basis(runtime::Int3, const QuadraticFrame &,
                             const mesh::MeshTopology &,
                             const mesh::MeshGeometry &);

struct DeterministicQr final {
  std::size_t rows{};
  std::size_t columns{};
  std::uint32_t rank{};
  double condition_estimate{};
  std::vector<std::uint32_t> pivots;
  std::uint64_t pivot_fingerprint{};
  std::vector<double> thin_q;
  std::vector<double> upper_r;

  std::vector<double>
  functional_weights(const std::vector<double> &functional) const;
};

DeterministicQr factorize_design_matrix(const std::vector<double> &matrix,
                                        std::size_t rows, std::size_t columns);

class ValidatedGeometryScope final {
public:
  ValidatedGeometryScope(const mesh::MeshTopology &,
                         const mesh::MeshGeometry &);
  ~ValidatedGeometryScope() noexcept;

  ValidatedGeometryScope(const ValidatedGeometryScope &) = delete;
  ValidatedGeometryScope &operator=(const ValidatedGeometryScope &) = delete;
  ValidatedGeometryScope(ValidatedGeometryScope &&) = delete;
  ValidatedGeometryScope &operator=(ValidatedGeometryScope &&) = delete;
};

bool geometry_validation_is_scoped(const mesh::MeshTopology &,
                                   const mesh::MeshGeometry &) noexcept;

class BoundaryAuthorityCoverageScope final {
public:
  BoundaryAuthorityCoverageScope();
  ~BoundaryAuthorityCoverageScope() noexcept;

  BoundaryAuthorityCoverageScope(const BoundaryAuthorityCoverageScope &) =
      delete;
  BoundaryAuthorityCoverageScope &
  operator=(const BoundaryAuthorityCoverageScope &) = delete;
};

bool boundary_authority_coverage_is_scoped() noexcept;

class QuadraticReconstructionWeights final {
public:
  struct AffineBoundaryFunctional final {
    std::vector<WeightedDonor> donors;
    double boundary_coefficient{};
  };

  static std::vector<WeightedDonor>
  value_weights(const QuadraticReconstruction &, runtime::Real3 point_m);
  static std::vector<WeightedDonor>
  cell_average_weights(const QuadraticReconstruction &,
                       runtime::Int3 global_cell, const mesh::MeshTopology &,
                       const mesh::MeshGeometry &);
  static std::vector<WeightedDonor>
  directional_gradient_weights(const QuadraticReconstruction &,
                               runtime::Real3 point_m,
                               runtime::Real3 direction);
  static std::vector<WeightedDonor>
  origin_constrained_directional_gradient_weights(
      const QuadraticReconstruction &, runtime::Real3 point_m,
      runtime::Real3 direction);
  static AffineBoundaryFunctional
  origin_normal_gradient_constrained_value_weights(
      const QuadraticReconstruction &, runtime::Real3 point_m);
  static AffineBoundaryFunctional
  origin_normal_gradient_constrained_directional_gradient_weights(
      const QuadraticReconstruction &, runtime::Real3 point_m,
      runtime::Real3 direction);
  static const std::vector<mesh::GlobalCellId> &
  donor_global_ids(const QuadraticReconstruction &) noexcept;
  static const std::vector<runtime::Int3> &
  donor_stencil_images(const QuadraticReconstruction &) noexcept;
  static double donor_extrema_bounded_value(
      const QuadraticReconstruction &, runtime::Real3 point_m,
      const runtime::FieldView<const double> &, std::size_t component);
  static double
  value_with_origin_constraint(const QuadraticReconstruction &, runtime::Real3,
                               const runtime::FieldView<const double> &,
                               std::size_t, double);
  static runtime::Real3 gradient_with_origin_constraint(
      const QuadraticReconstruction &, runtime::Real3,
      const runtime::FieldView<const double> &, std::size_t, double);
  static double value_with_origin_normal_gradient(
      const QuadraticReconstruction &, runtime::Real3,
      const runtime::FieldView<const double> &, std::size_t, double);
  static runtime::Real3 gradient_with_origin_normal_gradient(
      const QuadraticReconstruction &, runtime::Real3,
      const runtime::FieldView<const double> &, std::size_t, double);
  static QuadraticReconstruction with_boundary_authority(
      const QuadraticReconstruction &, std::uint64_t, int,
      const QuadraticReconstruction &,
      std::shared_ptr<const std::vector<BoundaryAuthorityOwner>>);
  static std::uint64_t boundary_authority_link(const QuadraticReconstruction &);
  static int boundary_authority_owner_rank(const QuadraticReconstruction &);
  static QuadraticReconstruction
  boundary_authority_reconstruction(const QuadraticReconstruction &);
  static const std::vector<BoundaryAuthorityOwner> &
  boundary_authority_catalog(const QuadraticReconstruction &);
};

} // namespace hundun::immersed::detail

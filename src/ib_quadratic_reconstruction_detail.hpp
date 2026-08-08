// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/ib_quadratic_reconstruction.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace hundun::immersed::detail {

struct BoundaryAuthorityOwner final {
  std::uint64_t link{};
  int owner_rank{};
};

double value_with_origin_constraint(const QuadraticReconstruction &,
                                    runtime::Real3 point_m,
                                    const runtime::FieldView<const double> &,
                                    std::size_t component,
                                    double value_at_origin);
runtime::Real3
gradient_with_origin_constraint(const QuadraticReconstruction &,
                                runtime::Real3 point_m,
                                const runtime::FieldView<const double> &,
                                std::size_t component, double value_at_origin);
double value_with_origin_normal_gradient(
    const QuadraticReconstruction &, runtime::Real3 point_m,
    const runtime::FieldView<const double> &, std::size_t component,
    double normal_gradient_at_origin);
runtime::Real3 gradient_with_origin_normal_gradient(
    const QuadraticReconstruction &, runtime::Real3 point_m,
    const runtime::FieldView<const double> &, std::size_t component,
    double normal_gradient_at_origin);

QuadraticReconstruction with_boundary_authority(
    const QuadraticReconstruction &point_reconstruction, std::uint64_t link,
    int owner_rank, const QuadraticReconstruction &authority_reconstruction,
    std::shared_ptr<const std::vector<BoundaryAuthorityOwner>> catalog);
std::uint64_t boundary_authority_link(const QuadraticReconstruction &);
int boundary_authority_owner_rank(const QuadraticReconstruction &);
QuadraticReconstruction
boundary_authority_reconstruction(const QuadraticReconstruction &);
const std::vector<BoundaryAuthorityOwner> &
boundary_authority_catalog(const QuadraticReconstruction &);

} // namespace hundun::immersed::detail

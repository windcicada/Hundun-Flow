// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "hundun/bc_basic_boundary.hpp"
#include "hundun/fvm_matrix_free_poisson.hpp"

namespace hundun::finite_volume {

PoissonBoundarySpec
make_poisson_boundary_spec(const boundary::BoundaryRegistry &registry);

} // namespace hundun::finite_volume

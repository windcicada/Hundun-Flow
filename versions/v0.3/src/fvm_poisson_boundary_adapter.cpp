// SPDX-License-Identifier: Apache-2.0

#include "hundun/fvm_poisson_boundary_adapter.hpp"

#include "hundun/rt_error.hpp"

namespace hundun::finite_volume {

PoissonBoundarySpec
make_poisson_boundary_spec(const boundary::BoundaryRegistry &registry) {
  if (!registry.open_domain()) {
    return {PressureConstraintMode::constant_nullspace, std::nullopt};
  }
  const auto outlet = registry.pressure_outlet_patch_id();
  if (!outlet.has_value()) {
    throw runtime::Error(
        "open boundary registry has no pressure reference patch");
  }
  return {PressureConstraintMode::pressure_reference_patch, outlet};
}

} // namespace hundun::finite_volume

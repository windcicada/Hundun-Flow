// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_flow.hpp"

namespace hundun::v04::detail {

inline MgHierarchyPolicy production_pressure_mg_policy() noexcept {
  MgHierarchyPolicy policy;
  policy.pre_sweeps = 1U;
  policy.post_sweeps = 2U;
  policy.point_smoother = MgPointSmootherKind::chebyshev_jacobi;
  policy.cycle = MgCycleKind::f_cycle;
  policy.chebyshev_lower_spectrum_fraction = 0.3;
  return policy;
}

#if defined(HUNDUN_V04_ENABLE_TEST_ACCESS)

// Arm on every participating rank.  The next refresh injects a local
// post-halo thermophysical revalidation failure only on failing_rank.
void arm_piso_post_halo_revalidation_failure_once_for_test(
    int failing_rank) noexcept;
void clear_piso_post_halo_revalidation_failure_for_test() noexcept;

#endif

struct PressureAssemblyBinding {
  const CartesianKernelPlan* kernels{};
  PisoCouplerWorkspace workspace{};
  Int3 cells{};
  PlanFingerprint coupler{};
  PlanFingerprint pressure_reference_plan{};
  PisoIntermediateCertificate current{};
  BoundaryThermophysicalGhostContext thermophysical_context{};
  PlanFingerprint geometry_fingerprint{};
};

Status assemble_pressure_system_impl(
    const PressureAssemblyBinding& binding,
    const PressureCorrectionInput& input,
    PressureCorrectionSystemView system,
    PressureCorrectionCertificate& certificate) noexcept;

}  // namespace hundun::v04::detail

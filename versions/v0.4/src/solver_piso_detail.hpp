// SPDX-License-Identifier: Apache-2.0

#pragma once

#include "hundun/v04_flow.hpp"

namespace hundun::v04::detail {

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

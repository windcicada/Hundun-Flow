// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#error "material-density transport test access is test-build only"
#endif

#include <cstddef>

namespace hundun::flow {
class FlowState;
class MaterialDensityTransport;

namespace test {

class MaterialDensityTransportTestAccess final {
public:
  static void reset() noexcept;
  static void set_density_residual(double value, int rank = -1) noexcept;
  static void set_transport_residual(std::size_t field, double value,
                                     int rank = -1) noexcept;
  static void set_mass_conservation_defect(double value,
                                           int rank = -1) noexcept;
  static void set_transport_conservation_defect(std::size_t field, double value,
                                                int rank = -1) noexcept;
  static void force_attempt_identity_wrap(FlowState &) noexcept;
  static void
  force_finalization_identity_wrap(MaterialDensityTransport &) noexcept;
};

} // namespace test
} // namespace hundun::flow

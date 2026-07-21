// SPDX-License-Identifier: Apache-2.0
#pragma once

#ifndef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#error "material-density transport test access is test-build only"
#endif

#include <cstddef>

namespace hundun::flow {
class FlowState;
class MaterialDensityTransport;
class MaterialDensityTransportReport;

namespace test {

struct MaterialConservationScaleInput final {
  double current_integral{};
  double next_integral{};
  double history_integral{};
  double history_coefficient{};
  double boundary_scale{};
  double cancellation_current{};
  double cancellation_next{};
  double cancellation_history{};
};

enum class MaterialReportCorruption : unsigned char {
  scalar,
  vector_size,
  vector_element,
  availability,
  seal
};

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
  static void set_accepted_transport_value(FlowState &, bool history,
                                           std::size_t field,
                                           std::size_t local_cell,
                                           double value);
  static void
  force_finalization_identity_wrap(MaterialDensityTransport &) noexcept;
  static double
  conservation_denominator(const MaterialConservationScaleInput &) noexcept;
  static void corrupt_report(MaterialDensityTransportReport &,
                             MaterialReportCorruption) noexcept;
};

} // namespace test
} // namespace hundun::flow

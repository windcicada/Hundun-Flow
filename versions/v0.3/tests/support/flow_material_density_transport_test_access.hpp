// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#ifndef HUNDUN_FLOW_ENABLE_TEST_ACCESS
#error "material-density transport test access is test-build only"
#endif

#include "flow_checkpoint_v2_detail.hpp"
#include "flow_density_closure_detail.hpp"
#include "hundun/rt_error.hpp"

#include <cmath>
#include <cstddef>
#include <limits>

namespace hundun::flow {
class FlowState;
class MaterialDensityTransport;
class MaterialDensityTransportReport;

namespace detail {
void material_transport_reset_raw() noexcept;
void material_transport_set_density_residual_raw(double, int) noexcept;
void material_transport_set_transport_residual_raw(std::size_t, double,
                                                   int) noexcept;
void material_transport_set_mass_conservation_defect_raw(double,
                                                         int) noexcept;
void material_transport_set_transport_conservation_defect_raw(
    std::size_t, double, int) noexcept;
double material_transport_conservation_denominator_raw(
    double, double, double, double, double, double, double,
    double) noexcept;
} // namespace detail

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
  static void reset() noexcept { detail::material_transport_reset_raw(); }
  static void set_density_residual(double value, int rank = -1) noexcept {
    detail::material_transport_set_density_residual_raw(value, rank);
  }
  static void set_transport_residual(std::size_t field, double value,
                                     int rank = -1) noexcept {
    detail::material_transport_set_transport_residual_raw(field, value, rank);
  }
  static void set_mass_conservation_defect(double value,
                                           int rank = -1) noexcept {
    detail::material_transport_set_mass_conservation_defect_raw(value, rank);
  }
  static void set_transport_conservation_defect(std::size_t field, double value,
                                                int rank = -1) noexcept {
    detail::material_transport_set_transport_conservation_defect_raw(
        field, value, rank);
  }
  static void force_attempt_identity_wrap(FlowState &) noexcept;
  static void set_accepted_transport_value(FlowState &, bool history,
                                           std::size_t field,
                                           std::size_t local_cell,
                                           double value);
  static void set_accepted_density_value(FlowState &, bool history,
                                         std::size_t local_cell,
                                         double value);
  static void
  force_finalization_identity_wrap(MaterialDensityTransport &) noexcept;
  static double conservation_denominator(
      const MaterialConservationScaleInput &input) noexcept {
    return detail::material_transport_conservation_denominator_raw(
        input.current_integral, input.next_integral, input.history_integral,
        input.history_coefficient, input.boundary_scale,
        input.cancellation_current, input.cancellation_next,
        input.cancellation_history);
  }
  static void corrupt_report(MaterialDensityTransportReport &,
                             MaterialReportCorruption) noexcept;
};

#ifdef HUNDUN_FLOW_ENABLE_TEST_ACCESS
inline void MaterialDensityTransportTestAccess::force_attempt_identity_wrap(
    FlowState &state) noexcept {
  state.impl_->attempt_identity = std::numeric_limits<std::uint64_t>::max();
}

inline void MaterialDensityTransportTestAccess::set_accepted_transport_value(
    FlowState &state, bool history, std::size_t field, std::size_t local_cell,
    double value) {
  if (field >= state.impl_->fields.transported_cell_fields.size())
    throw runtime::Error("material test transport index is invalid");
  const auto extent = state.impl_->layout.cell_interior_extent;
  const auto x = static_cast<std::size_t>(extent.x);
  const auto y = static_cast<std::size_t>(extent.y);
  const auto z = static_cast<std::size_t>(extent.z);
  if (y != 0U && x > std::numeric_limits<std::size_t>::max() / y)
    throw runtime::Error("flow-state cell count overflows");
  const auto xy = x * y;
  if (z != 0U && xy > std::numeric_limits<std::size_t>::max() / z)
    throw runtime::Error("flow-state cell count overflows");
  if (local_cell >= xy * z)
    throw runtime::Error("material test transport index is invalid");
  const int i =
      static_cast<int>(local_cell % static_cast<std::size_t>(extent.x));
  const std::size_t yz = local_cell / static_cast<std::size_t>(extent.x);
  const int j = static_cast<int>(yz % static_cast<std::size_t>(extent.y));
  const int k = static_cast<int>(yz / static_cast<std::size_t>(extent.y));
  auto &storage = history ? state.impl_->history : state.impl_->committed;
  auto view = storage.acquire_write<double>(
      state.impl_->access, runtime::PhaseId{1800U}, runtime::ActorId{1800U},
      state.impl_->fields.transported_cell_fields[field]);
  view(i, j, k, 0) = value;
}

inline void MaterialDensityTransportTestAccess::set_accepted_density_value(
    FlowState &state, bool history, std::size_t local_cell, double value) {
  const auto extent = state.impl_->layout.cell_interior_extent;
  const auto x = static_cast<std::size_t>(extent.x);
  const auto y = static_cast<std::size_t>(extent.y);
  const auto z = static_cast<std::size_t>(extent.z);
  if (y != 0U && x > std::numeric_limits<std::size_t>::max() / y)
    throw runtime::Error("flow-state cell count overflows");
  const auto xy = x * y;
  if (z != 0U && xy > std::numeric_limits<std::size_t>::max() / z)
    throw runtime::Error("flow-state cell count overflows");
  if (local_cell >= xy * z)
    throw runtime::Error("material test density index is invalid");
  const int i =
      static_cast<int>(local_cell % static_cast<std::size_t>(extent.x));
  const std::size_t yz = local_cell / static_cast<std::size_t>(extent.x);
  const int j = static_cast<int>(yz % static_cast<std::size_t>(extent.y));
  const int k = static_cast<int>(yz / static_cast<std::size_t>(extent.y));
  auto &storage = history ? state.impl_->history : state.impl_->committed;
  auto view = storage.acquire_write<double>(
      state.impl_->access, runtime::PhaseId{1800U}, runtime::ActorId{1800U},
      state.impl_->fields.density);
  view(i, j, k, 0) = value;
}

inline void
MaterialDensityTransportTestAccess::force_finalization_identity_wrap(
    MaterialDensityTransport &transport) noexcept {
  detail::DensityClosureBridge::
      force_transport_finalization_identity_wrap_raw(transport);
}

inline void MaterialDensityTransportTestAccess::corrupt_report(
    MaterialDensityTransportReport &report,
    MaterialReportCorruption corruption) noexcept {
  switch (corruption) {
  case MaterialReportCorruption::scalar:
    ++report.finalization_identity_;
    break;
  case MaterialReportCorruption::vector_size:
    report.transport_normalized_l2_.push_back(0.0);
    break;
  case MaterialReportCorruption::vector_element:
    if (!report.transport_normalized_l2_.empty())
      report.transport_normalized_l2_.front() =
          std::nextafter(report.transport_normalized_l2_.front(), 1.0);
    break;
  case MaterialReportCorruption::availability:
    report.density_residual_available_ = !report.density_residual_available_;
    break;
  case MaterialReportCorruption::seal:
    report.seal_ ^= 1U;
    break;
  }
}
#endif

} // namespace test
} // namespace hundun::flow

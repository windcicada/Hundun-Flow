// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/v04_flow.hpp"
#include "field_view_interval_detail.hpp"
#include "solver_cartesian_detail.hpp"
#include "solver_equation_detail.hpp"
#include <cstring>

namespace hundun::v04::detail {

// Temporal Schur diagonal of the conservative p/h system at fixed Y:
// Ep - Eh Cp/Ch = a0 V (-1 - rho rho_p/rho_h) = a0 V Cv/R.
// Bound EOS views are immutable for one Krylov solve. Spatially dominated
// systems retain their multigrid preconditioner; selection belongs to the
// caller's collective matrix policy. Inactive IBM rows return zero.
class PressureEnergyTemporalPreconditioner final : public LinearPreconditioner {
 public:
  Status bind(const CartesianKernelPlan& kernels, double a0,
              ConstFieldView density, ConstFieldView rho_p,
              ConstFieldView rho_h, PressureEnergyCellActivity activity,
              LinearPreconditionerCertificate lifecycle) noexcept {
    const Int3 shape = kernels.cells();
    const std::size_t count = static_cast<std::size_t>(shape.x) * shape.y * shape.z;
    if (kernels.fingerprint() == 0U || !std::isfinite(a0) || a0 <= 0.0 ||
        lifecycle.identity.fingerprint == 0U ||
        lifecycle.collective_fingerprint == 0U ||
        !valid_cell_view(density, shape, 0U, 1U, 0U) ||
        !valid_cell_view(rho_p, shape, 0U, 1U, 0U) ||
        !valid_cell_view(rho_h, shape, 0U, 1U, 0U) ||
        (activity.cells.size != 0U &&
         (activity.cells.size != count || activity.cells.data == nullptr ||
          activity.local_fingerprint == 0U || activity.collective_fingerprint == 0U)))
      return {StatusCode::invalid_plan, 1581U};
    PressureEnergyTemporalPreconditioner next;
    next.kernels_ = &kernels;
    next.a0_ = a0;
    next.density_ = density;
    next.rho_p_ = rho_p;
    next.rho_h_ = rho_h;
    next.activity_ = activity;
    bool valid = true;
    next.cells([&](Int3 c, std::size_t i) {
      if (!next.active(i)) return;
      const double d = next.diagonal(c);
      valid = valid && std::isfinite(density.unchecked(c, 0U)) &&
          density.unchecked(c, 0U) > 0.0 &&
          std::isfinite(rho_p.unchecked(c, 0U)) && rho_p.unchecked(c, 0U) > 0.0 &&
          std::isfinite(rho_h.unchecked(c, 0U)) && rho_h.unchecked(c, 0U) < 0.0 &&
          std::isfinite(d) && d > 0.0;
    });
    if (!valid) return {StatusCode::numerical_failure, 1581U};
    next.certificate_ = lifecycle;
    // Both supported nonsymmetric solvers accept the fixed-general contract.
    next.certificate_.preconditioner_class = LinearPreconditionerClass::fixed_general;
    next.certificate_.status_scope = LinearPreconditionerStatusScope::rank_local;
    next.certificate_.apply_lifecycle = LinearPreconditionerApplyLifecycle::per_call_checked;
    auto& fingerprint = next.certificate_.collective_fingerprint;
    const auto mix = [&](std::uint64_t value) { fingerprint ^= value; fingerprint *= UINT64_C(1099511628211); };
    std::uint64_t a0_bits{};
    std::memcpy(&a0_bits, &a0, sizeof(a0_bits));
    mix(UINT64_C(0x7363687572706374));
    mix(a0_bits);
    mix(activity.collective_fingerprint);
    for (const auto field : {density, rho_p, rho_h}) {
      mix(field.field);
      mix(field.revision);
    }
    if (fingerprint == 0U) fingerprint = 1U;
    *this = next;
    return {};
  }

  LinearPreconditionerCertificate certificate() const noexcept override {
    return certificate_;
  }
  std::uint64_t applications() const noexcept { return applications_; }

  Status apply(ConstFieldView input, FieldView output,
               std::uint32_t) noexcept override {
    if (kernels_ == nullptr || certificate_.collective_fingerprint == 0U ||
        !valid_cell_view(input, kernels_->cells(), 0U, 1U, 0U) ||
        !valid_cell_view(output, kernels_->cells(), 0U, 1U) ||
        field_views_overlap(input, as_const(output)))
      return {StatusCode::invalid_plan, 1582U};
    for (const auto field : {density_, rho_p_, rho_h_})
      if (field_views_overlap(input, field) ||
          field_views_overlap(as_const(output), field))
        return {StatusCode::invalid_plan, 1582U};
    bool finite = true;
    cells([&](Int3 c, std::size_t i) {
      const double value = active(i) ? input.unchecked(c, 0U) / diagonal(c) : 0.0;
      finite = finite && std::isfinite(value);
      output.unchecked(c, 0U) = value;
    });
    if (finite) ++applications_;
    return finite ? Status{} : Status{StatusCode::numerical_failure, 1582U};
  }

 private:
  template<class F> void cells(F f) const noexcept {
    const auto n = kernels_->cells();
    std::size_t i = 0U;
    for (int z = 0; z < n.z; ++z)
      for (int y = 0; y < n.y; ++y)
        for (int x = 0; x < n.x; ++x, ++i) f(Int3{x, y, z}, i);
  }
  bool active(std::size_t i) const noexcept {
    return activity_.cells.size == 0U || activity_.cells.data[i] != 0U;
  }
  double diagonal(Int3 c) const noexcept {
    return a0_ * cell_volume(*kernels_, c) *
        (-1.0 - density_.unchecked(c, 0U) * rho_p_.unchecked(c, 0U) /
                    rho_h_.unchecked(c, 0U));
  }
  const CartesianKernelPlan* kernels_{};
  double a0_{};
  ConstFieldView density_{}, rho_p_{}, rho_h_{};
  PressureEnergyCellActivity activity_{};
  LinearPreconditionerCertificate certificate_{};
  std::uint64_t applications_{};
};
} // namespace hundun::v04::detail

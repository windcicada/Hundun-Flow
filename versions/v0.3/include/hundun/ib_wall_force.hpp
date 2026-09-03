// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/ib_ghost_stencil_plan.hpp"
#include "hundun/rt_field_view.hpp"
#include "hundun/rt_types.hpp"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace hundun::runtime {
class MpiContext;
}

namespace hundun::immersed {

struct ForceComponents final {
  runtime::Real3 pressure_N{};
  runtime::Real3 viscous_N{};
  runtime::Real3 total_N{};
};

struct MomentComponents final {
  runtime::Real3 pressure_N_m{};
  runtime::Real3 viscous_N_m{};
  runtime::Real3 total_N_m{};
};

struct WallForceSample final {
  ForceComponents surface_traction;
  MomentComponents moment_about_global_origin;
  runtime::Real3 area_vector_closure_m2{};
  std::uint64_t quadrature_point_count{};
  int lowest_failing_rank{-1};
};

class WallForceIntegrator final {
public:
  static WallForceIntegrator create(const WallQuadraturePlan &,
                                    const runtime::MpiContext &);

  WallForceSample
  integrate(const runtime::FieldView<const double> &mechanical_pressure,
            const runtime::FieldView<const double> &velocity,
            const runtime::FieldView<const double> &velocity_gradient,
            const runtime::FieldView<const double> &mu_eff_by_cell) const;

private:
  WallForceIntegrator(const WallQuadraturePlan &plan,
                      const runtime::MpiContext &mpi,
                      std::vector<std::size_t> local_order) noexcept
      : plan_(&plan), mpi_(&mpi), local_order_(std::move(local_order)) {}

  const WallQuadraturePlan *plan_{};
  const runtime::MpiContext *mpi_{};
  std::vector<std::size_t> local_order_;
};

} // namespace hundun::immersed

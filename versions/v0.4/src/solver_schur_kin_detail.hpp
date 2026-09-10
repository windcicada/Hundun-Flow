// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09
#pragma once
#include "hundun/v04_execution.hpp"
#include "hundun/v04_flow.hpp"
#include "hundun/v04_ibm.hpp"
#include "solver_cartesian_detail.hpp"
#include <vector>
namespace hundun::v04::detail {
// Temporal kinetic-energy pressure response. The caller
// already exchanged Cartesian pressure ghosts and applied pressure BCs.
class SchurKinetic {
public:
  const CartesianKernelPlan *kernels{};
  const IbmEquationInterfacePlan *immersed{};
  RemoteDonorExchangePlan *donors{};
  StageId stage{};
  FieldId pressure_field{};
  ConstFieldView velocity, density, r_au;
  double a0{};
  PressureContinuityActivityView activity;
  mutable std::vector<double> memory;
  mutable FieldView gradient;
  mutable std::vector<double> pressure_memory;
  mutable FieldView pressure;
  mutable LinearOperatorFailureProvenance failure{};
  void reserve(Int3 n, int donor_reach, FieldId pressure_id,
               FieldId gradient_id) {
    pressure_field = pressure_id;
    std::size_t count = std::size_t(n.x) * n.y * n.z;
    memory.resize(3 * count);
    gradient.base = memory.data();
    gradient.interior = n;
    gradient.components = 3;
    gradient.stride_y = n.x;
    gradient.stride_z = std::size_t(n.x) * n.y;
    gradient.component_stride = count;
    gradient.field = gradient_id;
    gradient.revision = 1;
    gradient.storage_identity =
        reinterpret_cast<StorageIdentity>(memory.data());
    gradient.revision_domain = reinterpret_cast<RevisionDomainIdentity>(this);
    const int g = std::max(1, donor_reach);
    const std::size_t sy = n.x + 2 * g, sz = sy * (n.y + 2 * g),
                      nc = sz * (n.z + 2 * g);
    pressure_memory.resize(nc);
    pressure.base = pressure_memory.data() + g + g * sy + g * sz;
    pressure.interior = n;
    pressure.ghosts = {g, g, g};
    pressure.components = 1;
    pressure.stride_y = sy;
    pressure.stride_z = sz;
    pressure.component_stride = nc;
    pressure.field = pressure_field;
    pressure.revision = 1;
    pressure.storage_identity =
        reinterpret_cast<StorageIdentity>(pressure_memory.data());
    pressure.revision_domain = reinterpret_cast<RevisionDomainIdentity>(this);
  }

  Status validate() const noexcept {
    if (!kernels || !same_cells(pressure.interior, kernels->cells()) ||
        !std::isfinite(a0) || a0 <= 0.0 ||
        !valid_cell_view(velocity, kernels->cells(), 0U, 3U, 0U) ||
        !valid_cell_view(density, kernels->cells(), 0U, 1U, 0U) ||
        !valid_cell_view(r_au, kernels->cells(), 0U, 3U, 0U))
      return {StatusCode::invalid_plan, 1585U};
    const auto n = kernels->cells();
    const auto count = std::size_t(n.x) * n.y * n.z;
    if (activity.cells.size &&
        (!activity.cells.data || activity.cells.size != count))
      return {StatusCode::invalid_plan, 1585U};
    std::size_t i = 0;
    for (int z = 0; z < n.z; ++z)
      for (int y = 0; y < n.y; ++y)
        for (int x = 0; x < n.x; ++x, ++i) {
          if (activity.cells.size && !activity.cells.data[i])
            continue;
          const Int3 c{x, y, z};
          if (!std::isfinite(density.unchecked(c, 0)) ||
              density.unchecked(c, 0) <= 0.0)
            return {StatusCode::numerical_failure, 1585U};
          for (int a = 0; a < 3; ++a)
            if (!std::isfinite(velocity.unchecked(c, a)) ||
                !std::isfinite(r_au.unchecked(c, a)))
              return {StatusCode::numerical_failure, 1585U};
        }
    return {};
  }

  Status add(FieldView p, FieldView out) const {
    failure = {};
    auto shape = p.interior;
    for (int z = -1; z <= shape.z; ++z)
      for (int y = -1; y <= shape.y; ++y)
        for (int x = -1; x <= shape.x; ++x) {
          unsigned outside = (x < 0 || x >= shape.x) + (y < 0 || y >= shape.y) +
                             (z < 0 || z >= shape.z);
          if (outside <= 1)
            pressure.unchecked({x, y, z}, 0) = p.unchecked({x, y, z}, 0);
        }
    if (donors) {
      auto s = donors->exchange(stage, {&pressure, 1});
      if (!s) {
        failure = {s, LinearOperatorStatusScope::collective, -1};
        return s;
      }
    }
    const auto read = as_const(pressure);
    const KernelInvocation invocation{
        {&read, 1}, {&gradient, 1}, {{0, 0, 0}, p.interior}, 0, 0, 1,
        0,          nullptr};
    auto status = cartesian_gradient(*kernels, invocation);
    if (status && immersed)
      status = immersed->correct_pressure_gradient(read, gradient);
    if (!status)
      return status;
    auto n = p.interior;
    std::size_t i = 0;
    for (int z = 0; z < n.z; ++z)
      for (int y = 0; y < n.y; ++y)
        for (int x = 0; x < n.x; ++x, ++i) {
          if (activity.cells.size && !activity.cells.data[i])
            continue;
          Int3 c{x, y, z};
          double dk = 0;
          for (int a = 0; a < 3; ++a)
            dk -= velocity.unchecked(c, a) * r_au.unchecked(c, a) *
                  gradient.unchecked(c, a);
          out.unchecked(c, 0) +=
              a0 * cell_volume(*kernels, c) * density.unchecked(c, 0) * dk;
        }
    return {};
  }
};
} // namespace hundun::v04::detail

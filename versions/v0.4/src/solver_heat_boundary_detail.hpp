// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#pragma once
#include "hundun/v04_flow.hpp"
#include "solver_cartesian_detail.hpp"
#include <cmath>

namespace hundun::v04::detail {

// Replace the Cartesian ghost-derived flux by the prescribed outward heat.
// The same lift is used by the equation residual and accepted transport rates.
inline Status apply_heat_flux_boundary(
    const EnthalpyEquationPlan& plan, const CartesianKernelPlan& kernels,
    ConstFieldView coordinate, ConstFieldView coefficient, KernelBox box,
    FieldView rate, Span<const std::uint8_t> activity = {}) noexcept {
  if (!plan.has_prescribed_heat_flux()) return {};
  const Int3 cells = kernels.cells();
  const auto component = [](Int3 p, unsigned a) { return a == 0 ? p.x : (a == 1 ? p.y : p.z); };
  const auto set = [](Int3& p, unsigned a, int v) { (a == 0 ? p.x : (a == 1 ? p.y : p.z)) = v; };
  for (unsigned f = 0; f < 6; ++f) {
    double q{};
    if (!plan.prescribed_heat_flux(static_cast<CartesianFace>(f), q)) continue;
    const unsigned a = f / 2;
    const bool high = f % 2;
    const int owner = high ? component(cells, a) - 1 : 0;
    Int3 begin = box.begin;
    Int3 end{begin.x + box.cells.x, begin.y + box.cells.y, begin.z + box.cells.z};
    if (owner < component(begin, a) || owner >= component(end, a)) continue;
    set(begin, a, owner);
    set(end, a, owner + 1);
    for (int z = begin.z; z < end.z; ++z)
      for (int y = begin.y; y < end.y; ++y)
        for (int x = begin.x; x < end.x; ++x) {
          const Int3 c{x, y, z};
          const auto flat = static_cast<std::size_t>(x) + static_cast<std::size_t>(cells.x) *
              (static_cast<std::size_t>(y) + static_cast<std::size_t>(cells.y) * z);
          if (activity.size && activity.data[flat] == 0) continue;
          Int3 ghost = c, face = c;
          set(ghost, a, owner + (high ? 1 : -1));
          set(face, a, high ? owner + 1 : owner);
          const auto axis = static_cast<CartesianAxis>(a);
          const double original = positive_transmissibility(kernels, coefficient, axis, face) *
              (coordinate.unchecked(ghost, 0) - coordinate.unchecked(c, 0));
          const double correction = (-q * face_area(kernels, axis, face) - original) /
              cell_volume(kernels, c);
          const double value = rate.unchecked(c, 0) + correction;
          if (!std::isfinite(value)) return {StatusCode::numerical_failure, 17820};
          rate.unchecked(c, 0) = value;
        }
  }
  return {};
}
}  // namespace hundun::v04::detail

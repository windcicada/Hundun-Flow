// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_flow.hpp"
#include "hundun/v04_ibm.hpp"
#include "field_view_interval_detail.hpp"
#include "solver_equation_detail.hpp"

namespace hundun::v04::detail {

// One allocation-free face authority for the target residual, persisted
// spatial rates and independent boundary energy ledger. Fourier conduction
// remains separate. Positive face_flux is heat into the left cell, i.e.
// sum_s (h_s-h_dep) Gamma_s (Y_right-Y_left) A/d; J_dep=-sum_s J_s.
class MixtureEnthalpyDiffusion {
 public:
  static Status validate(const EnthalpyEquationPlan& plan,
                         const EquationStateView& state,
                         const EquationMaterialView& material) noexcept {
    if (plan.unity_lewis_total_enthalpy_) return {};
    const auto count = plan.species_specs_.size();
    if (plan.kernels_ == nullptr || plan.fingerprint_ == 0U ||
        state.independent_species.size != count ||
        (count != 0U && state.independent_species.data == nullptr) ||
        (material.scalar_mass_diffusivity.size != 0U &&
         material.scalar_mass_diffusivity.data == nullptr))
      return {StatusCode::invalid_plan, 1431U};
    if (count == 0U) return {};
    const Int3 cells = plan.cells_;
    if (!valid_cell_view(state.temperature.trial, cells, 0U, 1U, 1U) ||
        (material.scalar_mass_diffusivity.size == 0U &&
         (!valid_cell_view(material.molecular_viscosity, cells, 0U, 1U, 1U) ||
          !valid_cell_view(material.effective_viscosity, cells, 0U, 1U, 1U))) ||
        (material.scalar_mass_diffusivity.size != 0U &&
         (material.scalar_mass_diffusivity.data == nullptr ||
          material.scalar_mass_diffusivity.size < count)))
      return {StatusCode::invalid_plan, 1431U};
    for (std::size_t s = 0U; s < count; ++s) {
      const auto field = state.independent_species.data[s].trial;
      if (field.field != plan.species_specs_[s].field ||
          !valid_cell_view(field, cells, 0U, 1U, 1U) ||
          (material.scalar_mass_diffusivity.size != 0U &&
           !valid_cell_view(material.scalar_mass_diffusivity.data[s],
                            cells, 0U, 1U, 1U)))
        return {StatusCode::invalid_plan, 1431U};
    }
    return {};
  }

  // Inputs have passed validate(). No output is published for a bad face.
  static Status face_flux(const EnthalpyEquationPlan& plan,
                          const EquationStateView& state,
                          const EquationMaterialView& material,
                          CartesianAxis axis, Int3 face,
                          double& out) noexcept {
    out = 0.0;
    if (plan.unity_lewis_total_enthalpy_ || plan.species_specs_.empty()) return {};
    const auto& kernels = *plan.kernels_;
    Int3 left = face;
    const auto normal = axis == CartesianAxis::x ? face.x :
                        (axis == CartesianAxis::y ? face.y : face.z);
    if (axis == CartesianAxis::x) --left.x;
    else if (axis == CartesianAxis::y) --left.y;
    else --left.z;
    const double t = interpolate_face(
        kernels, axis, normal, state.temperature.trial.unchecked(left, 0U),
        state.temperature.trial.unchecked(face, 0U));
    const double location = face_coordinate(kernels, axis, normal);
    const double dl = location - centre_coordinate(kernels, axis, normal - 1);
    const double dr = centre_coordinate(kernels, axis, normal) - location;
    double flux = 0.0;
    for (std::size_t s = 0U; s < plan.species_specs_.size(); ++s) {
      const auto y = state.independent_species.data[s].trial;
      const double jump = y.unchecked(face, 0U) - y.unchecked(left, 0U);
      if (jump == 0.0) continue;
      const auto& spec = plan.species_specs_[s];
      const auto gamma = [&](Int3 c) {
        if (material.scalar_mass_diffusivity.size != 0U)
          return material.scalar_mass_diffusivity.data[s].unchecked(c, 0U);
        const double molecular = material.molecular_viscosity.unchecked(c, 0U);
        const double effective = material.effective_viscosity.unchecked(c, 0U);
        double turbulent = effective - molecular;
        if (turbulent < 0.0 &&
            turbulent > -1.0e-12 * std::max(1.0, molecular)) turbulent = 0.0;
        if (!std::isfinite(molecular) || molecular < 0.0 ||
            !std::isfinite(effective) || turbulent < 0.0)
          return std::numeric_limits<double>::quiet_NaN();
        return molecular / spec.molecular_schmidt +
               turbulent / spec.turbulent_schmidt;
      };
      const double gl = gamma(left), gr = gamma(face);
      if (!std::isfinite(gl) || gl <= 0.0 || !std::isfinite(gr) || gr <= 0.0 ||
          !(dl > 0.0) || !(dr > 0.0)) {
        return {StatusCode::numerical_failure, 1432U};
      }
      // Exactly the species kernel's harmonic transmissibility, including
      // the physical inlet-face material authority (G23).
      const double conductance = face_area(kernels, axis, face) /
          (kernels.physical_inlet_material(static_cast<std::size_t>(axis), normal)
               ? (dl + dr) / (normal == 0 ? gl : gr) : dl / gl + dr / gr);
      double difference = 0.0;
      const Status status = plan.thermodynamics_.independent_species_enthalpy_difference(
          s, t, difference);
      if (!status) return status;
      flux += difference * conductance * jump;
    }
    if (!std::isfinite(flux)) return {StatusCode::numerical_failure, 1432U};
    out = flux;
    return {};
  }

  static Status add_rate(const EnthalpyEquationPlan& plan,
                         const EquationStateView& state,
                         const EquationMaterialView& material,
                         const IbmEquationInterfacePlan* immersed,
                         KernelBox box, FieldView rate) noexcept {
    Status status = validate(plan, state, material);
    const Int3 cells = plan.cells_;
    if (!status) return status;
    if (!valid_kernel_box(box, cells) || !valid_cell_view(rate, cells, 0U, 1U) ||
        (immersed != nullptr &&
         (immersed->fingerprint() == 0U || immersed->topology_ == nullptr)))
      return {StatusCode::invalid_plan, 1431U};
    if (plan.unity_lewis_total_enthalpy_ || plan.species_specs_.empty()) return {};
    const auto overlaps = [&](ConstFieldView input) {
      return input.base != nullptr && field_views_overlap(input, as_const(rate));
    };
    if (overlaps(state.temperature.trial) || overlaps(material.molecular_viscosity) ||
        overlaps(material.effective_viscosity))
      return {StatusCode::invalid_plan, 1431U};
    for (std::size_t s = 0U; s < plan.species_specs_.size(); ++s)
      if (overlaps(state.independent_species.data[s].trial) ||
          (material.scalar_mass_diffusivity.size != 0U &&
           overlaps(material.scalar_mass_diffusivity.data[s])))
        return {StatusCode::invalid_plan, 1431U};
    const Int3 end{box.begin.x + box.cells.x, box.begin.y + box.cells.y,
                   box.begin.z + box.cells.z};
    const auto region = immersed == nullptr ? Span<const std::uint8_t>{} :
                                             immersed->cell_activity();
    for (int z = box.begin.z; z < end.z; ++z)
      for (int y = box.begin.y; y < end.y; ++y)
        for (int x = box.begin.x; x < end.x; ++x) {
          const Int3 cell{x, y, z};
          const auto flat = (static_cast<std::size_t>(z) * cells.y + y) * cells.x + x;
          if (region.size != 0U && region.data[flat] == 0U) continue;
          double divergence = 0.0;
          for (std::size_t a = 0U; a < 3U; ++a) {
            const auto axis = static_cast<CartesianAxis>(a);
            Int3 plus = cell;
            if (a == 0U) ++plus.x;
            else if (a == 1U) ++plus.y;
            else ++plus.z;
            double minus_flux = 0.0, plus_flux = 0.0;
            status = face_flux(plan, state, material, axis, cell, minus_flux);
            if (status) status = face_flux(plan, state, material, axis, plus, plus_flux);
            if (!status) return status;
            divergence += plus_flux - minus_flux;
          }
          const double value = rate.unchecked(cell, 0U) +
                               divergence / cell_volume(*plan.kernels_, cell);
          if (!std::isfinite(value)) return {StatusCode::numerical_failure, 1432U};
          rate.unchecked(cell, 0U) = value;
        }
    if (immersed != nullptr) {
      // Match correct_impermeable_scalar_diffusion: remove exactly the
      // Cartesian cut-face contribution. Prescribed internal inlets also
      // have zero species diffusion; their advective enthalpy is separate.
      const auto links = immersed->topology_->links();
      for (std::size_t index = 0U; index < links.size; ++index) {
        const auto& link = links.data[index];
        const auto c = link.fluid_local_index, s = link.solid_local_index;
        if (c.x < box.begin.x || c.x >= end.x || c.y < box.begin.y || c.y >= end.y ||
            c.z < box.begin.z || c.z >= end.z) continue;
        const auto axis = c.x != s.x ? CartesianAxis::x :
                          (c.y != s.y ? CartesianAxis::y : CartesianAxis::z);
        const int delta = axis == CartesianAxis::x ? s.x - c.x :
                          (axis == CartesianAxis::y ? s.y - c.y : s.z - c.z);
        double flux = 0.0;
        status = face_flux(plan, state, material, axis, delta > 0 ? s : c, flux);
        if (!status) return status;
        rate.unchecked(c, 0U) -= (delta > 0 ? flux : -flux) /
                                 cell_volume(*plan.kernels_, c);
      }
    }
    return {};
  }
};

}  // namespace hundun::v04::detail

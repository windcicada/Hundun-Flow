// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#pragma once

#include "hundun/v04_app.hpp"
#include "solver_equation_detail.hpp"
#include "solver_mixture_enthalpy_diffusion_detail.hpp"
#include "solver_mixture_enthalpy_convection_detail.hpp"
#include "solver_viscous_detail.hpp"

namespace hundun::v04::detail {

struct ProductBoundaryBalanceHistory {
  bool valid{};
  std::uint64_t epoch_start_step{};
  long double mass{}, previous_mass{}, energy{}, previous_energy{};
};

// All fields are already certified by the final equation assembly. Only
// rank-local scalar accumulators are owned here; no field or state is changed.
inline Status collect_terminal_equations(
    const CartesianKernelPlan& kernels, const EquationStateView& state,
    BdfCoefficients bdf, ConstFaceFluxView flux,
    ConstFieldView momentum_residual, ConstFieldView energy_residual,
    Span<const std::uint8_t> activity, ReductionEngine& reductions,
    DriverTerminalEquationReport& out, Int3 global_cells, MeshPatch patch,
    int rank, Span<const std::uint32_t> interface_cells) noexcept {
  out = {};
  const Int3 cells = kernels.cells();
  std::array<long double, 10U> sum{};
  std::array<double, 3U> maximum{};
  std::array<ReductionMaximumLocation, 3U> worst{};
  std::array<double, 6U> region_max{};
  std::array<long double, 8U> region_sum{};
  std::size_t interface_index = 0U;
  Status local;
  std::size_t flat = 0U;
  for (std::int32_t z = 0; z < cells.z && local; ++z)
    for (std::int32_t y = 0; y < cells.y && local; ++y)
      for (std::int32_t x = 0; x < cells.x; ++x, ++flat) {
        if (activity.size != 0U && activity.data[flat] == 0U) continue;
        const Int3 cell{x, y, z};
        const double volume = cell_volume(kernels, cell);
        const double rho = state.density.trial.unchecked(cell, 0U);
        // EBTopology freezes interface_cells in increasing local flat order.
        while (interface_index < interface_cells.size &&
               interface_cells.data[interface_index] < flat) ++interface_index;
        const std::size_t region = interface_index < interface_cells.size &&
            interface_cells.data[interface_index] == flat ? 0U : 1U;
        region_sum[6U + region] += 1.0;
        const std::uint64_t gid =
            (static_cast<std::uint64_t>(patch.begin.z + z) * global_cells.y +
             patch.begin.y + y) * global_cells.x + patch.begin.x + x;
        const double h = state.enthalpy.trial.unchecked(cell, 0U);
        const double p = state.pressure_reference +
                         state.pressure_perturbation.trial.unchecked(cell, 0U);
        const double temporal =
            volume * (bdf.a0 * rho +
                      bdf.a1 * state.density.accepted.unchecked(cell, 0U) +
                      (bdf.order == 2U
                           ? bdf.a2 * state.density.previous.unchecked(cell, 0U)
                           : 0.0));
        const double continuity = temporal +
            flux.x.unchecked({x + 1, y, z}) - flux.x.unchecked(cell) +
            flux.y.unchecked({x, y + 1, z}) - flux.y.unchecked(cell) +
            flux.z.unchecked({x, y, z + 1}) - flux.z.unchecked(cell);
        const double energy = energy_residual.unchecked(cell, 0U);
        double kinetic = 0.0;
        double momentum_work = 0.0;
        for (std::uint8_t component = 0U; component < 3U; ++component) {
          const double u = state.velocity.trial.unchecked(cell, component);
          const double residual = momentum_residual.unchecked(cell, component);
          if (!std::isfinite(u) || !std::isfinite(residual))
            local = {StatusCode::numerical_failure, 10212U};
          kinetic += 0.5 * u * u;
          momentum_work += u * residual;
          maximum[component] = std::max(maximum[component], std::abs(residual));
          sum[component] += std::abs(residual);
          const double correction = std::abs(residual) / (bdf.a0 * rho * volume);
          auto& selected = worst[component];
          if (!selected.valid || correction > selected.value ||
              (correction == selected.value && gid < selected.global_location))
            selected = {true, correction, gid, rank,
                {residual, rho * volume, static_cast<double>(patch.begin.x + x),
                 static_cast<double>(patch.begin.y + y), static_cast<double>(patch.begin.z + z)}};
          const auto bin = region * 3U + component;
          region_max[bin] = std::max(region_max[bin], correction);
          region_sum[bin] += static_cast<long double>(correction) * correction;
        }
        const double internal = (rho * h - p) * volume;
        const double kinetic_energy = rho * kinetic * volume;
        const double defect = energy + momentum_work - kinetic * continuity;
        if (!std::isfinite(continuity) || !std::isfinite(energy) ||
            !std::isfinite(internal) || !std::isfinite(kinetic_energy) ||
            !std::isfinite(defect)) {
          local = {StatusCode::numerical_failure, 10212U};
          break;
        }
        sum[3U] += continuity;
        sum[4U] += energy;
        sum[5U] += std::abs(energy);
        sum[6U] += defect;
        sum[7U] += rho * volume;
        sum[8U] += internal;
        sum[9U] += kinetic_energy;
      }
  std::array<double, 10U> local_values{}, global{};
  for (std::size_t i = 0U; i < sum.size(); ++i)
    local_values[i] = static_cast<double>(sum[i]);
  // ProductDriver guarantees at least eight reduction slots even for the
  // smallest Krylov configuration. Keep the audit inside that capacity.
  Status status = reductions.checked_sum(
      {local_values.data(), 8U}, {global.data(), 8U}, local);
  if (status)
    status = reductions.checked_sum(
        {local_values.data() + 8U, 2U}, {global.data() + 8U, 2U}, {});
  std::array<double, 3U> global_maximum{};
  if (status)
    status = reductions.checked_max(
        {maximum.data(), maximum.size()},
        {global_maximum.data(), global_maximum.size()}, {});
  if (!status) return status;
  DriverTerminalEquationReport report;
  report.valid = true;
  report.final_flux = flux.revision;
  report.momentum_linf = global_maximum;
  std::copy_n(global.begin(), 3U, report.momentum_l1.begin());
  report.continuity_signed = global[3U];
  report.energy_signed = global[4U];
  report.energy_l1 = global[5U];
  report.total_equation_defect = global[6U];
  report.mass = global[7U];
  report.internal_energy = global[8U];
  report.kinetic_energy = global[9U];
  std::array<double, 8U> local_regions{}, regions{};
  for (std::size_t i = 0U; i < region_sum.size(); ++i)
    local_regions[i] = static_cast<double>(region_sum[i]);
  std::array<double, 6U> region_maximum{};
  status = reductions.checked_sum({local_regions.data(), 8U}, {regions.data(), 8U}, {});
  if (status) status = reductions.checked_max({region_max.data(), 6U}, {region_maximum.data(), 6U}, {});
  if (status) status = reductions.checked_max_locations(
      {worst.data(), worst.size()}, {report.momentum_worst.data(), report.momentum_worst.size()}, {});
  if (!status) return status;
  const double reference = std::sqrt(2.0 * report.kinetic_energy / report.mass);
  report.momentum_reference_velocity = reference;
  report.momentum_normalization_valid = std::isfinite(reference) && reference > 0.0;
  for (std::size_t region = 0U; region < 2U; ++region) {
    report.momentum_region_cells[region] = static_cast<std::uint64_t>(regions[6U + region]);
    if (!report.momentum_normalization_valid) continue;
    for (std::size_t c = 0U; c < 3U; ++c) {
      report.momentum_normalized_linf[c] = report.momentum_worst[c].value / reference;
      report.momentum_region_normalized_linf[region][c] = region_maximum[region * 3U + c] / reference;
      report.momentum_region_normalized_rms[region][c] = regions[6U + region] > 0.0
          ? std::sqrt(regions[region * 3U + c] / regions[6U + region]) / reference : 0.0;
    }
  }
  out = report;
  return {};
}

inline Status collect_boundary_balance(
    const EnthalpyEquationPlan& enthalpy_plan,
    const CartesianKernelPlan& kernels, const SchemePlan& schemes,
    const BoundaryPlan& boundary, const EquationStateView& state,
    const EquationMaterialView& material, ConstFieldView gradient,
    BdfCoefficients bdf, ConstFaceFluxView flux,
    Span<const std::uint8_t> activity, FieldView kinetic, FieldView scratch,
    std::uint64_t accepted_step, const ProductBoundaryBalanceHistory& history,
    ReductionEngine& reductions, DriverConservationReport& out,
    ProductBoundaryBalanceHistory& pending,
    const IbmEquationInterfacePlan* immersed_interface = nullptr) noexcept {
  out = {};
  pending = {};
  const bool direct_h = enthalpy_plan.unity_lewis_total_enthalpy();
  const ConstFieldView thermal_coordinate = direct_h
      ? state.enthalpy.trial : state.temperature.trial;
  const ConstFieldView thermal_coefficient = direct_h
      ? material.enthalpy_diffusivity : material.thermal_conductivity;
  const Int3 cells = kernels.cells();
  const std::int32_t reach = kernels.reach();
  Status local;
  local = MixtureEnthalpyDiffusion::validate(enthalpy_plan, state, material);
  if (!valid_cell_view(as_const(kinetic), cells, 0U, 1U, kernels.reach()) ||
      !valid_cell_view(state.velocity.trial, cells, 0U, 3U, kernels.reach()))
    local = {StatusCode::invalid_plan, 10212U};
  const auto fluid = [&](Int3 cell) {
    const std::size_t index = (static_cast<std::size_t>(cell.z) * cells.y +
                               cell.y) * cells.x + cell.x;
    return activity.size == 0U || activity.data[index] != 0U;
  };
  // K is reconstructed as a transported scalar using the momentum scheme.
  // This defines a separate physical boundary energy flux; it is not assumed
  // to equal U times the discretised momentum convection (product-rule error).
  for (std::int32_t z = -reach; z < cells.z + reach && local; ++z)
    for (std::int32_t y = -reach; y < cells.y + reach; ++y)
      for (std::int32_t x = -reach; x < cells.x + reach; ++x) {
        const Int3 cell{x, y, z};
        double value = 0.0;
        for (std::uint8_t c = 0U; c < 3U; ++c) {
          const double u = state.velocity.trial.unchecked(cell, c);
          value += 0.5 * u * u;
        }
        kinetic.unchecked(cell, 0U) = value;
      }
  std::array<long double, 12U> sum{};
  const KernelBox box{{0, 0, 0}, cells};
  const auto convection = [&](ConstFieldView field, ConvectionScheme scheme,
                              IbmInterfaceInletFieldKind inlet_field) {
    const std::array<ConstFieldView, 1U> reads{field};
    const std::array<FieldView, 1U> writes{scratch};
    Status status = cartesian_convection(
        kernels, scheme, flux,
        {{reads.data(), reads.size()}, {writes.data(), writes.size()},
         box, 0U, 0U, 1U, flux.revision, nullptr});
    if (status && immersed_interface != nullptr)
      status = immersed_interface->add_source_convection_correction(
          {inlet_field, 0U}, scheme, field, 1.0, scratch, box);
    return status;
  };
  if (local) local = convection(state.enthalpy.trial, schemes.enthalpy(),
                                IbmInterfaceInletFieldKind::enthalpy);
  if (local) local = MixtureEnthalpyConvection::add_correction(
      enthalpy_plan,state,flux,immersed_interface,box,scratch);
  if (local)
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          if (fluid(cell))
            sum[7U] += cell_volume(kernels, cell) * scratch.unchecked(cell, 0U);
        }
  if (local) local = convection(as_const(kinetic), schemes.momentum(),
                                IbmInterfaceInletFieldKind::kinetic_energy);
  const std::array<ConstFieldView, 3U> density{
      state.density.trial, state.density.accepted, state.density.previous};
  const std::array<ConstFieldView, 3U> enthalpy{
      state.enthalpy.trial, state.enthalpy.accepted, state.enthalpy.previous};
  const std::array<ConstFieldView, 3U> pressure{
      state.pressure_perturbation.trial, state.pressure_perturbation.accepted,
      state.pressure_perturbation.previous};
  const std::array<ConstFieldView, 3U> velocity{
      state.velocity.trial, state.velocity.accepted, state.velocity.previous};
  const std::array<double, 3U> reference{
      state.pressure_reference, state.accepted_pressure_reference,
      state.previous_pressure_reference};
  if (local)
    for (std::int32_t z = 0; z < cells.z; ++z)
      for (std::int32_t y = 0; y < cells.y; ++y)
        for (std::int32_t x = 0; x < cells.x; ++x) {
          const Int3 cell{x, y, z};
          if (!fluid(cell)) continue;
          const double volume = cell_volume(kernels, cell);
          for (std::size_t layer = 0U; layer < 3U; ++layer) {
            const std::size_t i = layer == 2U && bdf.order == 1U ? 1U : layer;
            const double rho = density[i].unchecked(cell, 0U);
            double k = 0.0;
            for (std::uint8_t c = 0U; c < 3U; ++c) {
              const double u = velocity[i].unchecked(cell, c);
              k += 0.5 * u * u;
            }
            sum[layer] += volume * rho;
            sum[layer + 3U] += volume *
                (rho * (enthalpy[i].unchecked(cell, 0U) + k) -
                 (reference[i] + pressure[i].unchecked(cell, 0U)));
          }
          sum[6U] += flux.x.unchecked({x + 1, y, z}) - flux.x.unchecked(cell) +
                     flux.y.unchecked({x, y + 1, z}) - flux.y.unchecked(cell) +
                     flux.z.unchecked({x, y, z + 1}) - flux.z.unchecked(cell);
          sum[8U] += volume * scratch.unchecked(cell, 0U);
        }
  if (local && immersed_interface != nullptr) {
    double inlet_work = 0.0;
    local = immersed_interface->inlet_viscous_work_input(
        state.velocity.trial, material.effective_viscosity, inlet_work);
    if (local) sum[10U] += inlet_work;
  }
  // Below, only true external faces contribute. Periodic and partition faces cancel
  // in the conservative transport sums; stationary adiabatic IBM walls have
  // zero physical work/heat, irrespective of Cartesian ghost truncation error.
  for (std::uint8_t axis_index = 0U; axis_index < 3U && local; ++axis_index) {
    const auto axis = static_cast<CartesianAxis>(axis_index);
    const auto set_axis = [=](Int3& cell, std::int32_t value) {
      if (axis_index == 0U) cell.x = value;
      else if (axis_index == 1U) cell.y = value;
      else cell.z = value;
    };
    const std::int32_t count = axis_index == 0U ? cells.x :
                               (axis_index == 1U ? cells.y : cells.z);
    constexpr std::array<CartesianFace, 6U> external_faces{
        CartesianFace::x_min, CartesianFace::x_max, CartesianFace::y_min,
        CartesianFace::y_max, CartesianFace::z_min, CartesianFace::z_max};
    for (bool high : {false, true}) {
      const BoundaryFacePlan* face_plan = nullptr;
      local = boundary.face(external_faces[2U * axis_index + (high ? 1U : 0U)],
                            face_plan);
      if (!local) break;
      if (!face_plan->local_owner || face_plan->periodic) continue;
      Int3 end = cells;
      set_axis(end, 1);
      for (std::int32_t z = 0; z < end.z; ++z)
        for (std::int32_t y = 0; y < end.y; ++y)
          for (std::int32_t x = 0; x < end.x; ++x) {
            Int3 face{x, y, z}, owner{x, y, z};
            set_axis(face, high ? count : 0);
            set_axis(owner, high ? count - 1 : 0);
            if (!fluid(owner)) continue;
            Int3 left = face;
            const std::int32_t normal = high ? count : 0;
            set_axis(left, normal - 1);
            const double sign = high ? 1.0 : -1.0;
            const auto interpolate = [&](ConstFieldView f, std::uint8_t c) {
              return interpolate_face(kernels, axis, normal,
                                      f.unchecked(left, c), f.unchecked(face, c));
            };
            sum[9U] += sign * positive_transmissibility(
                kernels, thermal_coefficient, axis, face) *
                (thermal_coordinate.unchecked(face, 0U) -
                 thermal_coordinate.unchecked(left, 0U));
            double species_heat = 0.0;
            if (local)
              local = MixtureEnthalpyDiffusion::face_flux(
                  enthalpy_plan, state, material, axis, face, species_heat);
            sum[11U] += sign * species_heat;
            for (std::uint8_t c = 0U; c < 3U; ++c) {
              const double traction_area = viscous_face_traction_area(
                  kernels, state.velocity.trial, gradient,
                  material.effective_viscosity, axis, face, c);
              sum[10U] += sign * interpolate(state.velocity.trial, c) *
                          traction_area;
            }
          }
    }
  }
  std::array<double, 12U> values{}, global{};
  for (std::size_t i = 0U; i < sum.size(); ++i)
    values[i] = static_cast<double>(sum[i]);
  Status status = reductions.checked_sum(
      {values.data(), 8U}, {global.data(), 8U}, local);
  if (status)
    status = reductions.checked_sum(
        {values.data() + 8U, 4U}, {global.data() + 8U, 4U}, {});
  if (!status) return status;
  const long double mass = history.valid ? history.mass : global[1U];
  const long double previous_mass = history.valid ? history.previous_mass : global[2U];
  const long double energy = history.valid ? history.energy : global[4U];
  const long double previous_energy = history.valid ? history.previous_energy : global[5U];
  const long double energy_out = static_cast<long double>(global[7U]) +
                                global[8U] - global[9U] - global[10U] - global[11U];
  ProductBoundaryBalanceHistory next;
  next.valid = true;
  next.epoch_start_step = history.valid ? history.epoch_start_step : accepted_step;
  next.previous_mass = mass;
  next.previous_energy = energy;
  next.mass = (-bdf.a1 * mass - bdf.a2 * previous_mass - global[6U]) / bdf.a0;
  next.energy = (-bdf.a1 * energy - bdf.a2 * previous_energy - energy_out) / bdf.a0;
  DriverConservationReport report;
  report.epoch_start_step = next.epoch_start_step;
  report.mass_outflow = global[6U];
  report.enthalpy_outflow = global[7U];
  report.kinetic_energy_outflow = global[8U];
  report.conductive_heat_input = global[9U];
  report.species_enthalpy_diffusion_input = global[11U];
  report.viscous_work_input = global[10U];
  report.mass_bdf_rate = static_cast<double>(
      static_cast<long double>(bdf.a0) * global[0U] +
      static_cast<long double>(bdf.a1) * global[1U] +
      static_cast<long double>(bdf.a2) * global[2U]);
  report.total_energy_bdf_rate = static_cast<double>(
      static_cast<long double>(bdf.a0) * global[3U] +
      static_cast<long double>(bdf.a1) * global[4U] +
      static_cast<long double>(bdf.a2) * global[5U]);
  report.mass_balance_defect = report.mass_bdf_rate + global[6U];
  report.total_energy_balance_defect =
      report.total_energy_bdf_rate + static_cast<double>(energy_out);
  report.cumulative_mass_defect = static_cast<double>(global[0U] - next.mass);
  report.cumulative_energy_defect = static_cast<double>(global[3U] - next.energy);
  if (!std::isfinite(report.mass_balance_defect) ||
      !std::isfinite(report.total_energy_balance_defect) ||
      !std::isfinite(report.cumulative_mass_defect) ||
      !std::isfinite(report.cumulative_energy_defect))
    return {StatusCode::numerical_failure, 10212U};
  report.valid = true;
  out = report;
  pending = next;
  return {};
}

}  // namespace hundun::v04::detail

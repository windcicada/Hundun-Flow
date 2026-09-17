// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#pragma once
#include "solver_cold.hpp"
#include "solver_equation_detail.hpp"

namespace hundun::v04::detail {

// Complete correction matrix for an assembled scalar equation with frozen
// material, flux and common central/VLS face coefficients. Independent
// limited scalars use an upwind correction with the full reconstruction
// retained in the residual. The reference
// supplies integrated storage, sources, diffusion diagonal and true residual.
// Boundary values enter that residual; their homogeneous relations close the
// correction matrix. Rows are returned per volume for native Krylov solvers.
// Scratch and rows are attempt-local; callers agree status before solving.
inline Status close_mixture_scalar_rows(
    const CartesianKernelPlan& kernels, const BoundaryPlan& boundary,
    BoundaryStage boundary_stage, FieldId field, ConvectionScheme scheme,
    ConstFieldView boundary_velocity, const EquationAssemblyContext& context,
    const EquationSystemView& reference,
    const EquationAssemblyCertificate& certificate,
    FieldView boundary_variation, Span<ColdPressureRow> rows) noexcept {
  constexpr std::uint32_t invalid=17860, numerical=17861;
  const auto cells=kernels.cells();
  const auto count=std::size_t(cells.x)*cells.y*cells.z;
  const ConstFaceFluxView diffusion{as_const(reference.x_coefficient),
      as_const(reference.y_coefficient),as_const(reference.z_coefficient),
      context.face_flux,{}};
  if (!certificate.valid() || certificate.time!=context.time ||
      certificate.geometry!=context.geometry || certificate.face_flux!=context.face_flux ||
      certificate.scope!=context.scope || certificate.dt!=context.dt ||
      boundary.revision()!=context.boundary || rows.size!=count || !rows.data ||
      !valid_cell_view(reference.diagonal,cells,0,1) ||
      !valid_cell_view(reference.residual,cells,0,1) ||
      !valid_cell_view(as_const(boundary_variation),cells,0,1,1) ||
      !valid_equation_faces(reference,cells) ||
      !valid_flux_view(context.mass_flux,cells,context.face_flux) ||
      field_views_overlap(as_const(boundary_variation),reference.diagonal) ||
      field_views_overlap(as_const(boundary_variation),reference.residual) ||
      field_views_overlap(as_const(boundary_variation),reference.rhs) ||
      cell_input_aliases_faces(as_const(boundary_variation),reference,true) ||
      cell_face_views_overlap(as_const(boundary_variation),context.mass_flux.x) ||
      cell_face_views_overlap(as_const(boundary_variation),context.mass_flux.y) ||
      cell_face_views_overlap(as_const(boundary_variation),context.mass_flux.z))
    return {StatusCode::invalid_plan,invalid};
  if (context.mixture_transport &&
      (context.mixture_transport->linearization!=context.time ||
       context.mixture_transport->face_flux!=context.face_flux))
    return {StatusCode::invalid_plan,invalid};
  const auto activity=context.immersed_interface
      ? context.immersed_interface->cell_activity() : Span<const std::uint8_t>{};
  if (activity.size && (activity.size!=count || !activity.data))
    return {StatusCode::invalid_plan,invalid};
  for (int z=0;z<cells.z;++z) for (int y=0;y<cells.y;++y) for (int x=0;x<cells.x;++x)
    boundary_variation.unchecked({x,y,z},0)=1.;
  auto status=apply_homogeneous_scalar_boundary_ghosts(boundary_stage,boundary,
      field,boundary_variation,1,boundary_velocity);
  if (!status) return status;
  const std::array<ConstFaceFieldView,3> D{diffusion.x,diffusion.y,diffusion.z};
  const std::array<ConstFaceFieldView,3> F{
      context.mass_flux.x,context.mass_flux.y,context.mass_flux.z};
  std::size_t index{};
  for (int z=0;z<cells.z;++z) for (int y=0;y<cells.y;++y) for (int x=0;x<cells.x;++x,++index) {
    auto& row=rows.data[index]; row={};
    if (activity.size && activity.data[index]==0) {
      row.diagonal=1.; continue;
    }
    const Int3 cell{x,y,z};
    const double volume=cell_volume(kernels,cell);
    row.diagonal=reference.diagonal.unchecked(cell,0)/volume;
    row.rhs=-reference.residual.unchecked(cell,0)/volume;
    const int coordinate[3]{x,y,z},extent[3]{cells.x,cells.y,cells.z};
    for (unsigned f=0;f<6;++f) {
      const unsigned a=f/2;
      const bool high=f%2;
      const auto axis=static_cast<CartesianAxis>(a);
      auto neighbor=cell;
      (a==0 ? neighbor.x : a==1 ? neighbor.y : neighbor.z)+=high ? 1 : -1;
      const auto face=high ? neighbor : cell;
      double prescribed{};
      if (context.immersed_interface &&
          context.immersed_interface->prescribed_face_flux(axis,face,prescribed))
        continue; // The inlet composition is fixed in the physical residual.
      const int normal=coordinate[a]+(high ? 1 : 0);
      const double lower=interpolate_face(kernels,axis,normal,1.,0.);
      const double mass=F[a].unchecked(face)*(high ? 1. : -1.);
      // Limited scalar reconstruction stays in the complete residual. Use
      // its monotone upwind part for the deferred-correction matrix.
      const double owner=context.mixture_transport || scheme==ConvectionScheme::central2
          ? (high ? lower : 1.-lower) : (mass>=0. ? 1. : 0.);
      const double conductance=D[a].unchecked(face);
      if (!std::isfinite(conductance) || conductance<0 || !std::isfinite(mass))
        return {StatusCode::numerical_failure,numerical};
      row.diagonal+=mass*owner/volume;
      row.neighbour[f]=(conductance-mass*(1.-owner))/volume;
      if (coordinate[a]==(high ? extent[a]-1 : 0)) {
        const BoundaryFacePlan* face_plan{};
        status=boundary.face(static_cast<CartesianFace>(f),face_plan);
        if (!status || !face_plan) return {StatusCode::invalid_plan,invalid};
        if (face_plan->local_owner && !face_plan->periodic) {
          row.diagonal-=row.neighbour[f]*boundary_variation.unchecked(neighbor,0);
          row.neighbour[f]=0;
        }
      }
      if (!std::isfinite(row.neighbour[f])) return {StatusCode::numerical_failure,numerical};
    }
    if (!std::isfinite(row.diagonal) || row.diagonal<=0 || !std::isfinite(row.rhs))
      return {StatusCode::numerical_failure,numerical};
  }
  return {};
}

// Frozen-density stochastic transport uses rho*(q-q_old)/dt and the
// advective operator div(F*q)-q*div(F). COAST condif forms its diagonal
// from the six neighbour coefficients, followed by step(rho) and mixer.
// The ordinary conservative equation remains the reference assembly;
// this closure applies the explicit change of equation and its derivative.
inline Status close_frozen_density_scalar_rows(
    const CartesianKernelPlan& kernels, const BoundaryPlan& boundary,
    BoundaryStage stage, FieldId field, ConvectionScheme scheme,
    ConstFieldView velocity, const EquationAssemblyContext& context,
    const PrimitiveHistory& density, ConstFieldView scalar,
    const EquationSystemView& reference,
    const EquationAssemblyCertificate& certificate,
    FieldView variation, Span<ColdPressureRow> rows) noexcept {
  constexpr std::uint32_t invalid=17862, numerical=17863;
  const auto cells=kernels.cells();
  if (!valid_cell_view(scalar,cells,0,1,0) || scalar.field!=field ||
      !valid_cell_view(density.trial,cells,0,1,0) ||
      !valid_cell_view(density.accepted,cells,0,1,0) ||
      !valid_cell_view(density.previous,cells,0,1,0) ||
      context.reaction_endpoint.base || context.bdf.a2!=0 ||
      context.bdf.a1!=-context.bdf.a0 ||
      !std::isfinite(context.dt*context.bdf.a0) ||
      std::abs(context.dt*context.bdf.a0-1.)>4*std::numeric_limits<double>::epsilon() ||
      field_views_overlap(scalar,as_const(variation)) ||
      field_views_overlap(density.trial,as_const(variation)) ||
      field_views_overlap(density.accepted,as_const(variation)) ||
      field_views_overlap(density.previous,as_const(variation)))
    return {StatusCode::invalid_plan,invalid};
  for (int z=0;z<cells.z;++z) for (int y=0;y<cells.y;++y) for (int x=0;x<cells.x;++x) {
    const Int3 c{x,y,z};
    const double rho=density.trial.unchecked(c,0);
    if (!std::isfinite(rho) || rho<=0 ||
        rho!=density.accepted.unchecked(c,0) || rho!=density.previous.unchecked(c,0) ||
        !std::isfinite(scalar.unchecked(c,0)))
      return {StatusCode::invalid_plan,invalid};
  }
  auto status=close_mixture_scalar_rows(kernels,boundary,stage,field,scheme,
      velocity,context,reference,certificate,variation,rows);
  if (!status) return status;
  const auto activity=context.immersed_interface
      ? context.immersed_interface->cell_activity() : Span<const std::uint8_t>{};
  const auto flux=context.mass_flux;
  std::size_t index{};
  for (int z=0;z<cells.z;++z) for (int y=0;y<cells.y;++y) for (int x=0;x<cells.x;++x,++index) {
    if (activity.size && activity.data[index]==0) continue;
    const Int3 c{x,y,z};
    const double divergence=(flux.x.unchecked({x+1,y,z})-flux.x.unchecked(c)+
        flux.y.unchecked({x,y+1,z})-flux.y.unchecked(c)+
        flux.z.unchecked({x,y,z+1})-flux.z.unchecked(c))/cell_volume(kernels,c);
    auto& row=rows.data[index];
    row.diagonal-=divergence;
    row.rhs+=scalar.unchecked(c,0)*divergence;
    if (!std::isfinite(row.diagonal) || row.diagonal<=0 || !std::isfinite(row.rhs))
      return {StatusCode::numerical_failure,numerical};
  }
  return {};
}
} // namespace hundun::v04::detail

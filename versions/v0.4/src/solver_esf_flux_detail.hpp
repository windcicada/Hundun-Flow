// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "solver_mass_source_detail.hpp"
#include "field_view_interval_detail.hpp"
#include "solver_scalar_boundary_detail.hpp"
#include "hundun/v04_ibm.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace hundun::v04::detail {

struct StatisticalFluxReport {
  unsigned iterations{};
  double residual{};
  double convergence_residual{};
  double mass_pairing_residual{};
  double composition_closure{};
  // Rank-local failure evidence, retained across the consensus boundary.
  unsigned failure_reason{}; // 1 carrier mass, 2 solid face, 3 tuple, 4 closure, 5 inlet authority
  Int3 failure_cell{};
  std::array<double,3> failure_values{};
};

// A complete [Y_1,...,Y_N,h] tuple uses one pressure-flux correction matrix:
//   M q + sum_faces dm q_upwind = M_star q_star,
//   M_star = M_n - dt div(F_star), M = M_n - dt div(F),
//   dm = dt (F - F_star).
// q_star is the immutable post-transport/post-chemistry state. M here is the
// continuity mass, distinct from the provisional pressure EOS density. Their
// agreement remains a final coupled-flow acceptance condition. This map
// commutes with ensemble averaging and h -> h + sum_s c_s Y_s, including a
// signed auxiliary tuple. It contributes a flux, never a chemical source.
//
// The caller owns the candidate/scratch and the attempt transaction. A failed
// call leaves that candidate uncommitted; seed, accepted density and both
// fluxes are read-only. close_ghosts closes scalar AND enthalpy boundaries for
// this realization, after its MPI/periodic halo. Physical ghosts represent
// face values through the same mirror convention as the public equations.
template<class CloseGhosts>
Status correct_statistical_flux(const CartesianKernelPlan& kernels,
    ConstFieldView accepted_density, ConstFieldView seed,
    ConstFaceFluxView frozen_flux, ConstFaceFluxView final_flux, double dt,
    Span<const std::uint8_t> activity, bool signed_auxiliary,
    const BoundaryPlan& boundary, ConstFieldView velocity,
    FieldView iterate, FieldView next, HaloEngine& halo, StageId stage,
    ReductionEngine& reductions, CloseGhosts&& close_ghosts,
    StatisticalFluxReport& report, ConstFieldView candidate={},
    bool audit_only=false, const IbmEquationInterfacePlan* immersed=nullptr) noexcept {
  constexpr std::uint32_t invalid = 10231U, unconverged = 10232U;
  constexpr double tolerance = 32 * std::numeric_limits<double>::epsilon();
  report = {};
  const auto cells = seed.interior;
  const auto nc = seed.components;
  const std::size_t count = std::size_t(cells.x) * cells.y * cells.z;
  const double a0 = 1 / dt;
  Status local;
  if (kernels.fingerprint()==0 || kernels.cells().x!=cells.x ||
      kernels.cells().y!=cells.y || kernels.cells().z!=cells.z ||
      !valid_cells(cells) || nc < 3 || !std::isfinite(dt) || !(dt > 0) ||
      !std::isfinite(a0) ||
      !valid_cell_view(seed, cells, 0, nc, 0) ||
      !valid_cell_view(accepted_density, cells, 0, 1, 0) ||
      !valid_cell_view(velocity, cells, 0, 3, 0) ||
      !valid_cell_view(as_const(iterate), cells, 0, nc, 1) ||
      !valid_cell_view(as_const(next), cells, 0, nc, 0) ||
      field_views_overlap(seed,iterate) || field_views_overlap(seed,next) ||
      field_views_overlap(iterate,next) ||
      field_views_overlap(accepted_density,iterate) || field_views_overlap(accepted_density,next) ||
      field_views_overlap(velocity,iterate) || field_views_overlap(velocity,next) ||
      (candidate.base && (!valid_cell_view(candidate,cells,0,nc,0) ||
          field_views_overlap(candidate,iterate) || field_views_overlap(candidate,next))) ||
      (audit_only && !candidate.base) ||
      !valid_flux_view(frozen_flux, cells, frozen_flux.revision) ||
      !valid_flux_view(final_flux, cells, final_flux.revision) ||
      (activity.size && (activity.size != count || !activity.data)))
    local = {StatusCode::invalid_plan, invalid};
  const auto immersed_activity=immersed ? immersed->cell_activity() : Span<const std::uint8_t>{};
  if(local && immersed) {
    if(immersed_activity.size!=count || !immersed_activity.data || activity.size!=count)
      local={StatusCode::invalid_plan,invalid};
    if(local)local=immersed->validate_interface_flux(frozen_flux);
    if(local)local=immersed->validate_interface_flux(final_flux);
    if(!local)report.failure_reason=5;
  }
  const auto active = [&](std::size_t i) { return !activity.size || activity.data[i]; };
  const std::array<ConstFaceFieldView,3> base{frozen_flux.x,frozen_flux.y,frozen_flux.z};
  const std::array<ConstFaceFieldView,3> final{final_flux.x,final_flux.y,final_flux.z};
  // Arena identity is shared by disjoint fields. Actual storage intervals
  // establish aliasing, including immutable face and density inputs.
  for (auto flux : {frozen_flux, final_flux})
    for (auto face : {flux.x,flux.y,flux.z})
      if (cell_face_views_overlap(iterate,face) || cell_face_views_overlap(next,face))
        local={StatusCode::invalid_plan,invalid};
  // Complete local preflight precedes candidate publication and any halo.
  for (auto flux : {frozen_flux, final_flux})
    for (auto face : {flux.x,flux.y,flux.z})
      for (int z=0;z<face.extents.z && local;++z)
        for (int y=0;y<face.extents.y && local;++y)
          for (int x=0;x<face.extents.x;++x)
            if (!std::isfinite(face.unchecked({x,y,z}))) {
              local={StatusCode::rejected_step,invalid}; break;
            }
  std::size_t i{};
  for (int z=0;z<cells.z && local;++z) for (int y=0;y<cells.y && local;++y)
    for (int x=0;x<cells.x;++x,++i) {
      const Int3 c{x,y,z};
      if(immersed && bool(activity.data[i])!=bool(immersed_activity.data[i])) {
        local={StatusCode::invalid_plan,invalid};break;
      }
      const double mn=accepted_density.unchecked(c,0)*cell_volume(kernels,c);
      const double ms=active(i) ? frozen_density_carrier_mass(frozen_flux,c,mn,a0) : mn;
      const double mf=active(i) ? frozen_density_carrier_mass(final_flux,c,mn,a0) : mn;
      if (!(mn>0) || !(ms>0) || !(mf>0) || !std::isfinite(mn+ms+mf)) {
        report.failure_reason=1;report.failure_cell=c;report.failure_values={mn,ms,mf};
        local={StatusCode::rejected_step,invalid}; break;
      }
      if (!active(i)) {
        for (unsigned f=0;f<6;++f) {
          const unsigned axis=f/2; Int3 face=c;
          if (f%2) ++(axis==0 ? face.x : axis==1 ? face.y : face.z);
          if (base[axis].unchecked(face)!=0 || final[axis].unchecked(face)!=0) {
            double prescribed{};
            // A fixed IBM inlet carries the same already-accounted flux in
            // both histories. Its correction dm is exactly zero; every
            // other solid face retains the impermeable contract.
            if(immersed && immersed->prescribed_face_flux(static_cast<CartesianAxis>(axis),face,prescribed) &&
                base[axis].unchecked(face)==prescribed && final[axis].unchecked(face)==prescribed) continue;
            report.failure_reason=2;report.failure_cell=c;
            report.failure_values={double(f),base[axis].unchecked(face),final[axis].unchecked(face)};
            local={StatusCode::rejected_step,invalid};
          }
        }
      }
      long double sum{};
      for (unsigned s=0;s<nc;++s) {
        const double q=seed.unchecked(c,s);
        if (!std::isfinite(q) || !std::isfinite(ms*q) ||
            (!signed_auxiliary && active(i) && s+1<nc && (q<0 || q>1))) {
          report.failure_reason=3;report.failure_cell=c;report.failure_values={double(s),q,ms*q};
          local={StatusCode::rejected_step,invalid}; break;
        }
        if (s+1<nc) sum+=q;
      }
      if (local && active(i) && std::abs(sum-1)>2e-12L) {
        report.failure_reason=4;report.failure_cell=c;report.failure_values={double(sum),0,0};
        local={StatusCode::rejected_step,invalid};
      }
    }
  auto status=reductions.consensus(local);
  if (!status) return status;
  for (int z=0;z<cells.z;++z) for (int y=0;y<cells.y;++y)
    for (int x=0;x<cells.x;++x) for (unsigned s=0;s<nc;++s)
      iterate.unchecked({x,y,z},s)=(candidate.base ? candidate : seed).unchecked({x,y,z},s);
  for (unsigned iteration=0;iteration<128;++iteration) {
    ++iterate.revision;
    HaloTicket ticket;
    status=halo.begin(stage,{&iterate,1},{},ticket);
    if (status) status=halo.finish(ticket,{&iterate,1});
    if (status) status=close_ghosts(iterate);
    status=reductions.consensus(status);
    if (!status) return status;
    double maximum[4]{};
    i=0;
    for (int z=0;z<cells.z && local;++z) for (int y=0;y<cells.y && local;++y)
      for (int x=0;x<cells.x;++x,++i) {
        const Int3 c{x,y,z};
        if (!active(i)) {
          for (unsigned s=0;s<nc;++s) next.unchecked(c,s)=seed.unchecked(c,s);
          continue;
        }
        const double mn=accepted_density.unchecked(c,0)*cell_volume(kernels,c);
        const double ms=frozen_density_carrier_mass(frozen_flux,c,mn,a0);
        const double mf=frozen_density_carrier_mass(final_flux,c,mn,a0);
        std::array<double,6> dm{};
        std::array<Int3,6> neighbours{};
        long double diagonal=mf, net{};
        for (unsigned f=0;f<6;++f) {
          const unsigned a=f/2;
          const int sign=f%2 ? 1 : -1;
          Int3 face=c, neighbour=c;
          auto axis=[&](Int3& p)->int& {return a==0 ? p.x : a==1 ? p.y : p.z;};
          if (sign>0) ++axis(face);
          axis(neighbour)+=sign;
          neighbours[f]=neighbour;
          dm[f]=sign*(final[a].unchecked(face)-base[a].unchecked(face))/a0;
          diagonal+=std::max(dm[f],0.);
          net+=dm[f];
        }
        maximum[2]=std::max(maximum[2],double(std::abs(mf-ms+net)/(mf+ms)));
        long double sum{};
        for (unsigned s=0;s<nc;++s) {
          auto scalar=as_const(iterate);
          scalar.base+=s*scalar.component_stride; scalar.components=1;
          long double rhs=static_cast<long double>(ms)*seed.unchecked(c,s);
          for (unsigned f=0;f<6;++f) if (dm[f]<0)
            rhs-=static_cast<long double>(dm[f])*
                scalar_upwind_donor(boundary,scalar,neighbours[f],velocity);
          const long double lhs=diagonal*iterate.unchecked(c,s);
          const long double residual=std::abs(lhs-rhs), scale=std::abs(lhs)+std::abs(rhs);
          const long double quantization=.5L*diagonal*std::numeric_limits<double>::denorm_min();
          const double q=double(rhs/diagonal);
          if (!std::isfinite(q) || !std::isfinite(double(scale)) ||
              (!signed_auxiliary && s+1<nc && (q<0 || q>1))) {
            local={StatusCode::rejected_step,invalid}; break;
          }
          next.unchecked(c,s)=q;
          maximum[0]=std::max(maximum[0],scale==0 ? 0. : double(residual/scale));
          maximum[1]=std::max(maximum[1],scale==0 ? 0. :
              double(std::max(0.L,residual-quantization)/scale));
          if (s+1<nc) sum+=iterate.unchecked(c,s);
        }
        maximum[3]=std::max(maximum[3],double(std::abs(sum-1)));
      }
    double global[4]{};
    status=reductions.checked_max({maximum,4},{global,4},local);
    if (!status) return status;
    report={iteration,global[0],global[1],global[2],global[3]};
    if (audit_only) return {};
    if (report.convergence_residual<=tolerance) {
      if (report.composition_closure>2e-12 || report.mass_pairing_residual>tolerance)
        return {StatusCode::rejected_step,invalid};
      return {};
    }
    for (int z=0;z<cells.z;++z) for (int y=0;y<cells.y;++y)
      for (int x=0;x<cells.x;++x) for (unsigned s=0;s<nc;++s)
        iterate.unchecked({x,y,z},s)=next.unchecked({x,y,z},s);
  }
  return {StatusCode::rejected_step,unconverged};
}

} // namespace hundun::v04::detail

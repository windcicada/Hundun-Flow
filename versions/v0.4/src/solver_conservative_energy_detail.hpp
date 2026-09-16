// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once
#include "hundun/v04_flow.hpp"
#include "solver_cartesian_detail.hpp"

namespace hundun::v04::detail {
// Extra pressure work in rho*h-p0+rho*K when EOS pressure is fixed.
// The oriented, authoritative mass flux also defines this volume flux.
inline double mechanical_pressure_face_work(const CartesianKernelPlan& kernels,
    const EquationStateView& state, CartesianAxis axis, Int3 face,
    double mass_rate) noexcept {
  if (state.fixed_thermodynamic_pressure <= 0 || mass_rate == 0) return 0;
  Int3 lower=face;
  (axis==CartesianAxis::x ? lower.x : axis==CartesianAxis::y ? lower.y : lower.z)--;
  const auto index=[&](Int3 p){return axis==CartesianAxis::x ? p.x : axis==CartesianAxis::y ? p.y : p.z;};
  const double a=centre_coordinate(kernels,axis,index(lower)),
               b=centre_coordinate(kernels,axis,index(face)),
               w=(b-face_coordinate(kernels,axis,index(face)))/(b-a);
  const double rho=w*state.density.trial.unchecked(lower,0)+
      (1-w)*state.density.trial.unchecked(face,0);
  const double pi=state.pressure_reference-state.fixed_thermodynamic_pressure+
      w*state.pressure_perturbation.trial.unchecked(lower,0)+
      (1-w)*state.pressure_perturbation.trial.unchecked(face,0);
  return mass_rate*pi/rho;
}
inline double kinetic_energy(ConstFieldView velocity, Int3 cell) noexcept {
  double value=0.0;
  for(std::uint8_t c=0U;c<3U;++c) {
    const double u=velocity.unchecked(cell,c);
    value+=0.5*u*u;
  }
  return value;
}

// Uses the same reconstruction as ordinary Cartesian scalar convection,
// sampling K(U) on demand. Caller validates velocity reach and flux once.
double kinetic_convection_face(const CartesianKernelPlan& kernels,
    ConvectionScheme scheme, ConstFieldView velocity, CartesianAxis axis,
    Int3 face, double mass_rate) noexcept;

// Four axial samples at face offsets -2,-1,0,+1, reconstructed by the exact
// production scalar kernel. The caller owns validated stencil/flux authority.
double sampled_convection_face(const CartesianKernelPlan& kernels,
    ConvectionScheme scheme, const std::array<double,4U>& samples,
    CartesianAxis axis, Int3 face, double mass_rate) noexcept;

// Directional action of the same sampled face reconstruction. Uses the
// production semismooth zero-slope choice at limiter kinks; samples and
// variations share offsets -2,-1,0,+1. Caller validates stencil authority.
double sampled_convection_direction(const CartesianKernelPlan& kernels,
    ConvectionScheme scheme, const std::array<double,4U>& samples,
    const std::array<double,4U>& variation, CartesianAxis axis,
    Int3 face, double mass_rate) noexcept;

// Same scalar face leaf, preserving fractional subnormal values until the
// caller accumulates physical face products. Validated stencil/flux required.
long double precise_scalar_convection_face(const CartesianKernelPlan& kernels,
    ConvectionScheme scheme, ConstFieldView scalar, CartesianAxis axis,
    Int3 face, double mass_rate) noexcept;

inline double kinetic_temporal_density(const EquationStateView& state,
                                        BdfCoefficients bdf, Int3 c) noexcept {
  double rate=bdf.a0*state.density.trial.unchecked(c,0U)*
      kinetic_energy(state.velocity.trial,c)+
      bdf.a1*state.density.accepted.unchecked(c,0U)*
      kinetic_energy(state.velocity.accepted,c);
  if(bdf.order==2U)
    rate+=bdf.a2*state.density.previous.unchecked(c,0U)*
        kinetic_energy(state.velocity.previous,c);
  return rate;
}
} // namespace hundun::v04::detail

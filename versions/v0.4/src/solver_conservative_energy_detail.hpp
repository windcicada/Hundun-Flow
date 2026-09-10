// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once
#include "hundun/v04_flow.hpp"

namespace hundun::v04::detail {
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

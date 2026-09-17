// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once
#include "hundun/v04_ibm.hpp"
#include "hundun/v04_flow.hpp"

namespace hundun::v04::detail {
// Private complete Cartesian-rate boundary replacement. The physical flux
// contract differs from the public additive-correction compatibility API.
class IbmScalarTransport {
 public:
  enum class Quantity { independent_species, dependent_species, enthalpy, passive_scalar };
  struct Field { Quantity quantity; std::size_t component{}; };
  static Status convection(const IbmEquationInterfacePlan& plan,Field field,
      ConvectionScheme scheme,ConstFieldView q,ConstFaceFluxView flux,
      KernelBox box,FieldView rate,
      const MixtureTransportFaces* mixture=nullptr) noexcept;
  static Status transport(const IbmEquationInterfacePlan& plan,Field field,
      ConstFieldView q,ConstFieldView gamma,ConstFaceFluxView flux,
      KernelBox box,FieldView rate,const MixtureTransportFaces& mixture) noexcept;
  static Status convection(const IbmEquationInterfacePlan& plan,std::size_t species,
      ConvectionScheme scheme,ConstFieldView q,ConstFaceFluxView flux,
      KernelBox box,FieldView rate,
      const MixtureTransportFaces* mixture=nullptr) noexcept;
  static Status transport(const IbmEquationInterfacePlan& plan,std::size_t species,
      ConstFieldView q,ConstFieldView gamma,ConstFaceFluxView flux,
      KernelBox box,FieldView rate,const MixtureTransportFaces& mixture) noexcept;
  // Apply the same physical cut-face authority to an exported Cartesian
  // transport ledger. Sealed links carry zero; prescribed links carry phi*q_in.
  static Status constrain_flux(const IbmEquationInterfacePlan& plan,Field field,
      ConstFaceFluxView mass_flux,std::array<FaceFieldView,3> output) noexcept;
  static Status diffusion(const IbmEquationInterfacePlan& plan,ConstFieldView q,
      ConstFieldView gamma,KernelBox box,FieldView rate) noexcept;
  static double diffusion_diagonal(const IbmEquationInterfacePlan& plan,
      ConstFieldView gamma,Int3 cell) noexcept;
  // The assembler supplies its physical residual and masked diffusion
  // diagonal. Freeze solid corrections and remove their stored face edges.
  static Status constrain_rows(const IbmEquationInterfacePlan& plan,
      ConstFieldView q,KernelBox box,EquationSystemView system) noexcept;
 private:
  static Status convection_impl(const IbmEquationInterfacePlan& plan,Field field,
      ConvectionScheme scheme,ConstFieldView q,ConstFaceFluxView flux,
      KernelBox box,FieldView rate,const MixtureTransportFaces* mixture,
      ConstFieldView physical) noexcept;
};
} // namespace hundun::v04::detail

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
  static Status convection(const IbmEquationInterfacePlan& plan,std::size_t species,
      ConvectionScheme scheme,ConstFieldView q,ConstFaceFluxView flux,
      KernelBox box,FieldView rate) noexcept;
  static Status diffusion(const IbmEquationInterfacePlan& plan,ConstFieldView q,
      ConstFieldView gamma,KernelBox box,FieldView rate) noexcept;
};
} // namespace hundun::v04::detail

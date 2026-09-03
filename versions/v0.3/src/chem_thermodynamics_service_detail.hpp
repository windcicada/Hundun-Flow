// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09
#pragma once

#include "hundun/chem_composition.hpp"

namespace hundun::chemistry {

class ThermodynamicsService {
public:
  virtual ~ThermodynamicsService() = default;
  virtual const CompositionIdentity &composition() const noexcept = 0;
  virtual ThermodynamicProperties
  evaluate(const ThermochemicalPoint &) const = 0;
};

class TransportPropertyService {
public:
  virtual ~TransportPropertyService() = default;
  virtual TransportProperties
  evaluate(const ThermochemicalPoint &,
           const ThermodynamicProperties &) const = 0;
};

} // namespace hundun::chemistry

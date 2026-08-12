// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "src/chem_backend_detail.hpp"
#include "src/chem_thermodynamics_service_detail.hpp"

#include <memory>

namespace hundun::test {

class AnalyticReactingBackend final
    : public chemistry::ThermodynamicsService,
      public chemistry::TransportPropertyService,
      public chemistry::ChemistryBackend {
public:
  ~AnalyticReactingBackend() override;

  const chemistry::CompositionIdentity &composition() const noexcept override;
  chemistry::ThermodynamicProperties
  evaluate(const chemistry::ThermochemicalPoint &point) const override;
  chemistry::TransportProperties evaluate(
      const chemistry::ThermochemicalPoint &point,
      const chemistry::ThermodynamicProperties &thermodynamics) const override;
  chemistry::ChemistryIntervalReport
  integrate(const chemistry::ChemistryIntervalRequest &request) override;

private:
  friend std::unique_ptr<AnalyticReactingBackend>
  make_analytic_reacting_backend_for_tests();
  AnalyticReactingBackend();
  chemistry::CompositionIdentity composition_;
};

std::unique_ptr<AnalyticReactingBackend>
make_analytic_reacting_backend_for_tests();

} // namespace hundun::test

// SPDX-License-Identifier: Apache-2.0

#include "../support/turbulence_fixture.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string_view>

namespace {

using namespace hundun::v04;
using namespace hundun::v04::test;

bool expect(bool condition, std::string_view description) {
  if (!condition) {
    std::cerr << "FAIL: " << description << '\n';
  }
  return condition;
}

bool test_wale_invariants_and_wall_order() {
  bool passed = true;
  VelocityGradient shear;
  shear.value[1U] = 2.0;
  double viscosity = -1.0;
  passed &= expect(wale_kinematic_viscosity(shear, 0.1, 0.325,
                                            viscosity) &&
                       viscosity == 0.0,
                   "WALE is exactly zero for one-component pure shear");

  const auto near_wall = [](double distance, double& out) {
    VelocityGradient gradient;
    gradient.value[1U] = 1.0;
    gradient.value[3U] = distance;
    return wale_kinematic_viscosity(gradient, 0.1, 0.325, out);
  };
  double coarse = 0.0;
  double medium = 0.0;
  double fine = 0.0;
  passed &= expect(near_wall(0.008, coarse) && near_wall(0.004, medium) &&
                       near_wall(0.002, fine) && coarse > medium &&
                       medium > fine,
                   "WALE analytic near-wall family is positive and monotone");
  const double order_0 = std::log(coarse / medium) / std::log(2.0);
  const double order_1 = std::log(medium / fine) / std::log(2.0);
  passed &= expect(order_0 >= 2.8 && order_1 >= 2.8,
                   "WALE analytic near-wall viscosity has cubic scaling");

  VelocityGradient invalid = shear;
  invalid.value[4U] = std::numeric_limits<double>::quiet_NaN();
  const double marker = viscosity;
  passed &= expect(wale_kinematic_viscosity(invalid, 0.1, 0.325,
                                            viscosity)
                           .code == StatusCode::numerical_failure &&
                       viscosity == marker,
                   "WALE non-finite input rejects without publishing output");
  return passed;
}

bool test_wale_plan_lifecycle(GeometryKind kind) {
  TurbulencePlanSpec spec;
  spec.kind = TurbulenceKind::wale;
  TurbulenceFixture fixture;
  bool passed = expect(fixture.initialize(spec, kind),
                       "WALE plan compiles with one mu_eff authority");
  if (!passed) {
    return false;
  }
  for (VelocityGradient& gradient : fixture.gradients) {
    gradient.value[0U] = 0.4;
    gradient.value[4U] = -0.2;
    gradient.value[8U] = -0.2;
    gradient.value[1U] = 0.3;
  }
  TurbulenceCertificate certificate;
  passed &= expect(fixture.plan.update(fixture.input(), fixture.effective.view,
                                       certificate) &&
                       certificate.valid() && fixture.plan.update_count() == 1U,
                   "WALE publishes one certified effective viscosity field");
  for (double value : fixture.effective.storage) {
    passed &= expect(std::isfinite(value) && value >= 1.8e-5,
                     "WALE effective viscosity is finite and nonnegative");
  }
  const TurbulenceCertificate first = certificate;
  passed &= expect(fixture.plan.update(fixture.input(), fixture.effective.view,
                                       certificate) &&
                       certificate.state == first.state &&
                       fixture.plan.update_count() == 1U,
                   "equal revision tuple reuses the shared gradient result");
  TurbulenceUpdateInput revised = fixture.input(105U);
  passed &= expect(fixture.plan.update(revised, fixture.effective.view,
                                       certificate) &&
                       certificate.state != first.state &&
                       fixture.plan.update_count() == 2U,
                   "gradient revision invalidates only the turbulence numeric state");

  const std::vector<double> before = fixture.effective.storage;
  const TurbulenceCertificate before_certificate = certificate;
  fixture.gradients[1U].value[2U] =
      std::numeric_limits<double>::quiet_NaN();
  revised.gradient_revision = 106U;
  passed &= expect(fixture.plan
                           .update(revised, fixture.effective.view, certificate)
                           .code == StatusCode::numerical_failure &&
                       fixture.effective.storage == before &&
                       certificate.state == before_certificate.state,
                   "failed WALE update is atomic");
  return passed;
}

}  // namespace

int main(int argc, char** argv) {
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    return 2;
  }
  bool passed = test_wale_invariants_and_wall_order();
  passed &= test_wale_plan_lifecycle(GeometryKind::uniform);
  passed &= test_wale_plan_lifecycle(GeometryKind::tensor_stretched);
  MPI_Finalize();
  return passed ? 0 : 1;
}

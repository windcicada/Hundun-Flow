// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat: windcicada | Year.M: 2026.09

#include "src/chem_backend_detail.hpp"
#include "src/chem_thermodynamics_service_detail.hpp"
#include "tests/support/chem_analytic_backend.hpp"
#include "tests/support/test_main.hpp"

#include <cmath>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

static_assert(std::has_virtual_destructor_v<
              hundun::chemistry::ThermodynamicsService>);
static_assert(std::has_virtual_destructor_v<
              hundun::chemistry::TransportPropertyService>);
static_assert(
    std::has_virtual_destructor_v<hundun::chemistry::ChemistryBackend>);

namespace {

using hundun::chemistry::ChemistryBackend;
using hundun::chemistry::ChemistryIntervalReport;
using hundun::chemistry::ChemistryIntervalRequest;
using hundun::chemistry::ChemistryStatus;
using hundun::chemistry::CompositionIdentity;
using hundun::chemistry::ThermochemicalPoint;
using hundun::chemistry::ThermodynamicProperties;
using hundun::chemistry::ThermodynamicsService;
using hundun::chemistry::TransportProperties;
using hundun::chemistry::TransportPropertyService;

ThermochemicalPoint point(double y0) {
  return {101325.0, 300000.0, {y0, 1.0 - y0}};
}

class IndependentBackend final : public ThermodynamicsService,
                                 public TransportPropertyService,
                                 public ChemistryBackend {
public:
  IndependentBackend() {
    composition_.element_names = {"B"};
    composition_.species = {{"B", 1.0, {1}}, {"B2", 2.0, {2}}};
    composition_.fingerprint =
        hundun::chemistry::composition_identity_fingerprint(composition_);
  }

  const CompositionIdentity &composition() const noexcept override {
    return composition_;
  }
  ThermodynamicProperties evaluate(const ThermochemicalPoint &) const override {
    return {500.0, 3.0, 900.0, 1.5};
  }
  TransportProperties evaluate(const ThermochemicalPoint &,
                               const ThermodynamicProperties &) const override {
    return {2.0e-5, 0.03, {3.0e-5, 4.0e-5}};
  }
  ChemistryIntervalReport
  integrate(const ChemistryIntervalRequest &request) override {
    ChemistryIntervalReport report;
    report.final_state = request.state;
    report.integrated_rho_y_delta_kg_per_m3 = {0.0, 0.0};
    report.status = ChemistryStatus::success;
    report.completed_duration_s = request.duration_s;
    return report;
  }

private:
  CompositionIdentity composition_;
};

void test_backend_neutral_services_and_temperature_inversion() {
  auto backend = hundun::test::make_analytic_reacting_backend_for_tests();
  const ThermodynamicsService &thermo = *backend;
  const TransportPropertyService &transport = *backend;
  const ChemistryBackend &chemistry = *backend;

  HUNDUN_CHECK(thermo.composition().fingerprint ==
               chemistry.composition().fingerprint);
  const auto thermodynamics = thermo.evaluate(point(0.5));
  HUNDUN_CHECK_NEAR(thermodynamics.temperature_k, 300.19, 1.0e-12);
  HUNDUN_CHECK_NEAR(thermodynamics.density_kg_per_m3, 2.0, 0.0);
  HUNDUN_CHECK_NEAR(thermodynamics.cp_j_per_kg_k, 1000.0, 0.0);
  const auto properties = transport.evaluate(point(0.5), thermodynamics);
  HUNDUN_CHECK(properties.mixture_diffusivity_m2_per_s.size() == 2U);
  HUNDUN_CHECK_NEAR(properties.mixture_diffusivity_m2_per_s[1], 2.0e-5,
                    0.0);

  IndependentBackend independent;
  const ThermodynamicsService &independent_thermo = independent;
  HUNDUN_CHECK_NEAR(independent_thermo.evaluate(point(0.5)).temperature_k,
                    500.0, 0.0);
}

void test_integrated_delta_and_call_order_independence() {
  auto first = hundun::test::make_analytic_reacting_backend_for_tests();
  auto second = hundun::test::make_analytic_reacting_backend_for_tests();
  const ChemistryIntervalRequest request_a{point(0.75), 1.0, 0.5};
  const ChemistryIntervalRequest request_b{point(0.5), 2.0, 0.25};

  const auto a_then = first->integrate(request_a);
  const auto b_then = first->integrate(request_b);
  const auto b_first = second->integrate(request_b);
  const auto a_second = second->integrate(request_a);
  HUNDUN_CHECK(a_then.succeeded());
  HUNDUN_CHECK(b_then.succeeded());
  HUNDUN_CHECK_NEAR(a_then.final_state.mass_fractions[0],
                    0.25 + 0.5 * std::exp(-1.0), 1.0e-14);
  HUNDUN_CHECK_NEAR(a_then.integrated_rho_y_delta_kg_per_m3[0],
                    2.0 * (a_then.final_state.mass_fractions[0] - 0.75),
                    1.0e-14);
  HUNDUN_CHECK_NEAR(
      a_then.integrated_rho_y_delta_kg_per_m3[0] +
          a_then.integrated_rho_y_delta_kg_per_m3[1],
      0.0, 1.0e-14);
  HUNDUN_CHECK(a_then.final_state.p0_pa == request_a.state.p0_pa);
  HUNDUN_CHECK(a_then.final_state.h_tc_j_per_kg ==
               request_a.state.h_tc_j_per_kg);
  HUNDUN_CHECK(a_then.final_state.mass_fractions ==
               a_second.final_state.mass_fractions);
  HUNDUN_CHECK(b_then.final_state.mass_fractions ==
               b_first.final_state.mass_fractions);
}

void test_composition_mismatch_and_failure_report() {
  auto analytic = hundun::test::make_analytic_reacting_backend_for_tests();
  IndependentBackend independent;
  bool mismatch_rejected = false;
  try {
    hundun::chemistry::validate_backend_composition(
        analytic->composition(), independent);
  } catch (const std::invalid_argument &) {
    mismatch_rejected = true;
  }
  HUNDUN_CHECK(mismatch_rejected);

  auto invalid = point(0.5);
  invalid.p0_pa = 0.0;
  const auto report = analytic->integrate({invalid, 0.0, 0.5});
  HUNDUN_CHECK(report.status == ChemistryStatus::invalid_input);
  HUNDUN_CHECK(report.final_state.mass_fractions == invalid.mass_fractions);
  HUNDUN_CHECK(report.integrated_rho_y_delta_kg_per_m3 ==
               std::vector<double>({0.0, 0.0}));
  HUNDUN_CHECK(report.completed_duration_s == 0.0);
  HUNDUN_CHECK(report.internal_step_count == 0U);
}

} // namespace

int main() {
  return hundun::test::run([] {
    test_backend_neutral_services_and_temperature_inversion();
    test_integrated_delta_and_call_order_independence();
    test_composition_mismatch_and_failure_report();
  });
}

// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Github/Wechat:
// windcicada | Year.M: 2026.09
#include "models_chemistry_adapter_detail.hpp"
#include "models_esf_detail.hpp"
#include <cmath>
#include <iostream>
using namespace hundun::v04;
bool near(double actual, double expected, double absolute = 1e-14,
          double relative = 1e-13) {
  return std::isfinite(actual) && std::isfinite(expected) &&
         std::abs(actual - expected) <=
             absolute + relative * std::abs(expected);
}
int main() {
  for (bool reversed : {false, true}) {
    chemistry::detail::AnalyticIsomerBackend backend(2, reversed);
    const auto &gas = backend.gas_identity();
    const auto &id = backend.closure_identity();
    double values[12]{};
    for (std::size_t f = 0; f < 4; ++f) {
      values[f * 3 + (reversed ? 1 : 0)] = 1;
      values[f * 3 + 2] = 101850;
    }
    double p[]{101325, 101325, 101325, 101325}, rho[4];
    double d[2], h[2], w[2];
    portable::GasQuery query{{0, 1, 1},
                             gas.composition_fingerprint,
                             portable::GasStateCoordinates::pressure_enthalpy,
                             101325,
                             101850,
                             0,
                             values,
                             2};
    portable::GasQueryOutput sample{{}, d, h, w, 2};
    if (backend.query_gas(query, sample) != portable::Status::success)
      return 1;
    for (auto &r : rho)
      r = sample.sample.density_kg_per_m3;
    esf::detail::Workspace workspace(2);
    for (std::size_t fields : {2U, 4U}) {
      esf::detail::ReactionRequest q{
          {query.revision, id.fingerprint, fields, 2, values},
          query.revision,
          &id,
          &gas,
          p,
          rho,
          0,
          0.5};
      auto result = workspace.react(q, backend);
      if (result.status != portable::Status::success) {
        std::cerr << "status " << int(result.status) << " field "
                  << result.failure_field << '\n';
        return 2;
      }
      if (result.chemistry_call_count != 2 * fields ||
          !near(result.candidate.values[reversed ? 1 : 0],
                0.36787944117144233) ||
          result.candidate.values[2] != 101850) {
        std::cerr.precision(17);
        std::cerr << result.chemistry_call_count << ' '
                  << result.candidate.values[reversed ? 1 : 0] << ' '
                  << result.candidate.values[2] << '\n';
        return 3;
      }
      // Final T = 300 + 100(1-exp(-1)); constant pressure gives inverse-T
      // density.
      if (!near(result.final_densities_kg_per_m3[0] / rho[0],
                0.825963772790506964766474928233575))
        return 4;
    }
  }
  std::cout << "P1 P2 analytic ensemble PASS\n";
}

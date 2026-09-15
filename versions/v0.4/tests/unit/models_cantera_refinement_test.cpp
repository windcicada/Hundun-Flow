// SPDX-License-Identifier: Apache-2.0
#include "hundun/v04_cantera.hpp"
#include "../support/chemistry_test_support.hpp"
#include <array>
#include <fstream>

using namespace hundun::v04;
int main(int argc, char** argv) {
  HUNDUN_CHECK(argc == 3);
  std::ifstream input(argv[2]);
  double p{}, h{}, start{}, dt{};
  std::array<double,7> y{}, final{}, delta{};
  HUNDUN_CHECK(input >> p >> h >> start >> dt);
  for (auto& value : y) HUNDUN_CHECK(input >> value);
  chemistry::CanteraBackendConfig config;
  config.mechanism = {argv[1],
      "a1ca61b0b97847c59ea59e2c0373451c0d1c4cb6a46d3ddde21f7a1d64c05fa1", "jl4"};
  config.species_names = {"CH4","O2","CO2","CO","H2O","H2","N2"};
  config.chemistry = {1e-13,1e-24,100000};
  config.continuous_enthalpy = true;
  config.minimum_temperature = 273.15;
  config.maximum_temperature = 3500;
  const auto run = [&](chemistry::CanteraBackendConfig controls,
                       std::array<double,7>& endpoint) {
    auto runtime = std::make_shared<chemistry::CanteraBackendRuntime>(controls);
    chemistry::CanteraWorkspacePool pool(runtime,1);
    auto gas = chemistry::make_cantera_backend(controls,pool);
    const portable::GasAdvanceQuery query{
        {{1,32006,1},gas->gas_identity().composition_fingerprint,
         portable::GasStateCoordinates::pressure_enthalpy,p,h,0,y.data(),7},start,dt};
    portable::GasAdvanceOutput out{{},endpoint.data(),delta.data(),7};
    const auto status = gas->advance_gas(query,out);
    if (status == portable::Status::success) {
      HUNDUN_CHECK(out.completed_duration_s == dt);
      HUNDUN_CHECK(out.internal_step_count > 0 &&
                   out.internal_step_count <= unsigned(controls.chemistry.maximum_internal_steps));
      double sum{};
      for (auto value : endpoint) {
        HUNDUN_CHECK(std::isfinite(value) && value >= 0 && value <= 1);
        sum += value;
      }
      HUNDUN_CHECK_NEAR(sum,1.,2e-12);
    }
    return status;
  };
  // Captured GTMC endpoint at the original dt exceeds the 8*atol roundoff
  // budget. Internal refinement retains that budget and the same interval.
  HUNDUN_CHECK(run(config,final) == portable::Status::success);
  auto finer = config;
  finer.chemistry.absolute_tolerance = 1e-25;
  std::array<double,7> reference{}, repeated{};
  HUNDUN_CHECK(run(finer,reference) == portable::Status::success);
  for (unsigned i=0;i<7;++i) HUNDUN_CHECK_NEAR(final[i],reference[i],1e-12);
  HUNDUN_CHECK(run(config,repeated) == portable::Status::success && final == repeated);
  auto limited = config;
  limited.chemistry.maximum_internal_steps = 1;
  repeated.fill(-7.);
  HUNDUN_CHECK(run(limited,repeated) == portable::Status::provider_failure);
  for (auto value : repeated) HUNDUN_CHECK(value == -7.);
  return 0;
}

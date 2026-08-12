// SPDX-License-Identifier: Apache-2.0

#include "cantera/base/Solution.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/transport/Transport.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

std::uint64_t bits(double value) noexcept {
  std::uint64_t result{};
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

std::uint64_t hash(std::uint64_t state, double value) noexcept {
  state ^= bits(value);
  state *= UINT64_C(1099511628211);
  return state;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    throw std::runtime_error("usage: cantera_cpp_link_spike MECHANISM");
  }
  auto solution =
      Cantera::newSolution(argv[1], "synthetic-gas", "mixture-averaged");
  auto thermo = solution->thermo();
  thermo->setState_TPX(1000.0, Cantera::OneAtm, "A:0.9, B:0.1");
  std::uint64_t digest = UINT64_C(14695981039346656037);
  digest = hash(digest, thermo->temperature());
  digest = hash(digest, thermo->pressure());
  digest = hash(digest, thermo->density());
  digest = hash(digest, thermo->cp_mass());
  digest = hash(digest, solution->transport()->viscosity());
  std::vector<double> rates(thermo->nSpecies());
  solution->kinetics()->getNetProductionRates(rates.data());
  for (double rate : rates) {
    if (!std::isfinite(rate)) {
      throw std::runtime_error("non-finite Cantera production rate");
    }
    digest = hash(digest, rate);
  }
  std::cout << std::hex << digest << '\n';
  Cantera::appdelete();
}

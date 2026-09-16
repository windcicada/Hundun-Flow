// SPDX-License-Identifier: Apache-2.0
#include "../../src/solver_cold.hpp"
#include <iomanip>
#include <iostream>

int main() {
  using namespace hundun::v04::detail;
  double storage{}, dt{}, old{};
  while (std::cin >> storage >> dt >> old) {
    std::array<double, 6> neighbours{};
    ColdPressureRow spatial, row;
    for (auto& value : neighbours) std::cin >> value;
    std::cin >> spatial.diagonal;
    for (auto& value : spatial.neighbour) std::cin >> value;
    std::cin >> spatial.rhs;
    if (!std::cin || !time_centre_cold_row(spatial, old, neighbours, storage, dt, row))
      return 1;
    std::cout << std::setprecision(17) << row.diagonal;
    for (auto value : row.neighbour) std::cout << ' ' << value;
    std::cout << ' ' << row.rhs << '\n';
  }
}

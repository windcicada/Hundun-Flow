// SPDX-License-Identifier: Apache-2.0
#include "../../src/solver_cartesian_detail.hpp"
#include <iostream>
#include <iomanip>
int main() {
  double rho,volume,dt;std::array<double,6> flux;
  while(std::cin>>rho>>volume>>dt) {
    for(auto& f:flux)std::cin>>f;
    hundun::v04::detail::CellConvectiveCflResult result;
    if(hundun::v04::detail::evaluate_cell_convective_cfl(rho,volume,flux,{{1,1,1,1,1,1}},dt,result)
        !=hundun::v04::detail::CellConvectiveCflStatus::success)return 1;
    std::cout<<std::setprecision(17)<<result.directional_max<<' '<<result.out<<' '<<result.absolute<<'\n';
  }
}

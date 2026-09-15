// SPDX-License-Identifier: Apache-2.0
// Developed by WANG YUDONG | Email: wangyudong@buaa.edu.cn | Year.M: 2026.09
#include "core_esf_energy_detail.hpp"
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
using namespace hundun::v04;
namespace { bool observe{}; std::size_t calls{}, allocated{}; }
void* operator new(std::size_t bytes, std::align_val_t alignment) {
  if (observe) { ++calls; allocated += bytes; }
  void* memory{};
  if (posix_memalign(&memory,static_cast<std::size_t>(alignment),bytes ? bytes : 1))
    throw std::bad_alloc{};
  return memory;
}
void operator delete(void* memory,std::align_val_t) noexcept {std::free(memory);}
void operator delete(void* memory,std::size_t,std::align_val_t) noexcept {std::free(memory);}
int main() {
  EnthalpyEquationPlan equation;
  bool passed=true;
  for (const Int3 cells : {Int3{3,2,1},Int3{8,7,5},Int3{1,1,1}}) {
    detail::StatisticalEnergyLedger ledger;
    std::size_t predicted{};
    passed &= bool(FaceFluxStorage::workspace_bytes(cells,5,predicted));
    observe=true; calls=allocated=0;
    const auto first=ledger.initialize(cells,2,41,equation);
    observe=false;
    passed &= bool(first) && calls==1 && allocated==predicted &&
              ledger.owned_payload_bytes()==predicted;
    observe=true; calls=allocated=0;
    const auto retry=ledger.initialize(cells,2,42,equation);
    const auto rejected=ledger.initialize(cells,2,0,equation);
    const auto changed_shape=ledger.initialize({cells.x+1,cells.y,cells.z},2,43,equation);
    const auto recovered=ledger.initialize(cells,4,44,equation);
    observe=false;
    passed &= bool(retry) && !rejected && !changed_shape && bool(recovered) &&
              calls==0 && allocated==0 && !ledger.start_correction() &&
              ledger.balance({}).generation==0 && ledger.owned_payload_bytes()==predicted;
    ledger.discard_attempt();
    passed &= !ledger.start_correction();
    std::cout << "shape=" << cells.x << ',' << cells.y << ',' << cells.z
              << " bytes=" << predicted << " hot_allocations=" << calls << '\n';
  }
  for (const Int3 cells : {Int3{0,1,1},Int3{-1,1,1},
                           Int3{std::numeric_limits<int>::max(),1,1},
                           Int3{1000000000,1000000000,1000000000}}) {
    std::size_t bytes=17;
    passed &= !FaceFluxStorage::workspace_bytes(cells,5,bytes) && bytes==0;
  }
  std::size_t bytes=17;
  passed &= !FaceFluxStorage::workspace_bytes({2,2,2},0,bytes) && bytes==0;
  passed &= !FaceFluxStorage::workspace_bytes({2,2,2},std::numeric_limits<std::size_t>::max(),bytes) && bytes==0;
  std::cout << "statistical_energy_workspace passed=" << passed << '\n';
  return passed ? 0 : 1;
}
